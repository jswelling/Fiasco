/************************************************************
 *                                                          *
 *  smoother.c                                              *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1998 Department of Statistics,         *
 *                        Carnegie Mellon University        *
 *                                                          *
 *  This program is distributed in the hope that it will    *
 *  be useful, but WITHOUT ANY WARRANTY; without even the   *
 *  implied warranty of MERCHANTABILITY or FITNESS FOR A    *
 *  PARTICULAR PURPOSE.  Neither Carnegie Mellon University *
 *  nor any of the authors assume any liability for         *
 *  damages, incidental or otherwise, caused by the         *
 *  installation or use of this software.                   *
 *                                                          *
 *  CLINICAL APPLICATIONS ARE NOT RECOMMENDED, AND THIS     *
 *  SOFTWARE HAS NOT BEEN EVALUATED BY THE UNITED STATES    *
 *  FDA FOR ANY CLINICAL USE.                               *
 *                                                          *
 *                                                          *
 *  Original programming by Joel Welling, 8/98              *
 ************************************************************/
/* This package implements smoothing methods. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <assert.h>
#include "mri.h"
#include "fmri.h" /* Includes smoother.h */
#include "misc.h"
#include "stdcrg.h"

/* Notes:
  -At this moment, smoother->k is unused.
*/

static char rcsid[] = "$Id: smoother.c,v 1.17 2007/03/21 23:50:20 welling Exp $";

static sm_type default_type= SM_GAUSSIAN;
static double default_bandwidth= 45.0;
static double default_k= 60.0;
static double default_threshold= 0.0;
static SmootherThreshTest default_thresh_test= NULL;
static int default_sm_dim= 't';

/* Initialize the package */
void sm_init()
{
  /* Nothing to do */
}

void sm_parse_cl_opts()
{
  char buf[512];

  /* We assume that things like cl_scan() have already been called;
   * this parsing is for smoother-related args only.
   */
  if (cl_get("smoother_type","%option %s",buf)) {
    if (!strncasecmp(buf,"gauss",5)) default_type= SM_GAUSSIAN;
    else if (!strncasecmp(buf,"tri",3)) default_type= SM_TRIANGULAR;
    else if (!strncasecmp(buf,"pow",3)) default_type= SM_POWER;
    else if (!strncasecmp(buf,"none",4)) default_type= SM_NONE;
    else if (!strncasecmp(buf,"ddx",3)) default_type= SM_DDX;
    else if (!strncasecmp(buf,"shift",5)) default_type= SM_SHIFT;
    else if (!strncasecmp(buf,"linterp",7)) default_type= SM_LINTERP;
    else if (!strncasecmp(buf,"runningsum",10)) default_type= SM_RUNNINGSUM;
    else if (!strncasecmp(buf,"median",10)) default_type= SM_MEDIAN;
    else Abort("sm_parse_cl_opts: unknown smoother type <%s>\n",buf);
  }
  (void)cl_get("smoother_bandwidth","%option %lf",&default_bandwidth);
  (void)cl_get("smoother_k","%option %lf",&default_k);
  (void)cl_get("smoother_threshold", "%option %lf",&default_threshold);
}

static void trim_bounds( int here, int* low, int* high,
			 float* data, double thresh )
{
  int i;

  for (i=here-1; i>=*low; i--) {
    if (fabs(data[i]-data[here]) > thresh) {
      *low= i+1;
      break;
    }
  }
  for (i=here+1; i<=*high; i++) {
    if (fabs(data[i]-data[here]) > thresh) {
      *high= i-1;
      break;
    }
  }
}

static void trim_bounds_group( Smoother* sm, int here, int* low, int* high,
			       float** data_tbl, int ndata )
{
  int i;
  double thresh= sm->threshold;
  SmootherThreshTest test= sm->thresh_test;

  for (i=here-1; i>=*low; i--) {
    if (thresh < (*test)( sm, data_tbl, ndata, i, here )) {
      *low= i+1;
      break;
    }
  }
  for (i=here+1; i<=*high; i++) {
    if (thresh < (*test)( sm, data_tbl, ndata, i, here )) {
      *high= i-1;
      break;
    }
  }
}

static void check_kernel( Smoother* sm, int n )
{
  double* kernel;
  double tmp;
  int i;

  if (sm->data && n==sm->n) return;

  if (sm->data) free(sm->data);
  if (!(kernel=(double*)malloc(n*sizeof(double))))
    Abort("smoother:check_kernel: unable to allocate %d doubles!\n",n);
  sm->n= n;

  sm->data= (void*)kernel;
  switch (sm->type) {
  case SM_GAUSSIAN:
    for (i=0; i<n; i++) {
      double dist= (double)i/sm->bandwidth;
      kernel[i]= exp(-dist*dist);
    }
    sm->last= n-1;
    break;
  case SM_TRIANGULAR:
    for (i=0; i<n; i++) {
      double dist= (double)i/sm->bandwidth;
      if (dist<=0.5) kernel[i]= 1.0-dist;
      else kernel[i]= 0.0;
    }
    sm->last= ((int)sm->bandwidth)/2;
    break;
  case SM_POWER:
    tmp= pow(0.5,2.0/sm->bandwidth); /* band = full width at half max */
    for (i=0; i<n; i++) {
      kernel[i]= pow( tmp, (double)i );
    }
    sm->last= n-1;
    break;
  default: Abort("smoother: check_kernel: unknown kernel type %d!\n",
		 sm->type);
  }
}

static void kernel_smooth( Smoother* sm, 
			   float* data_in, float* data_out, int n,
			   unsigned char** missing, int z)
{
  int i;
  int k;
  double weight, total_weight;
  int band_low;
  int band_high;
  double val;
  double* kernel;

  check_kernel(sm, n);
  kernel= (double*)sm->data;

  if (!missing || (sm->smDim != 't' && sm->smDim != 'z')) {
    for (i=0; i<n; i++) {
      total_weight = 0.0;
      val = 0.0;
      band_low= 0;
      band_high= n-1;
      if (sm->threshold>0.0) trim_bounds( i, &band_low, &band_high, 
					  data_in, sm->threshold );
      if (i-band_low > sm->last) band_low= i - sm->last;
      for (k=band_low; k<i; k++) {
	val += kernel[i-k]*data_in[k];
	total_weight += kernel[i-k];
      }
      for (k=i; k<=band_high; k++) {
	val += kernel[k-i]*data_in[k];
	total_weight += kernel[k-i];
	if ((k-i)>=sm->last) break;
      } 
      if (total_weight != 0.0)
	data_out[i]= val/total_weight;
    }
  }
  else if (sm->smDim == 't') {
    assert( missing != NULL );
    for (i=0; i<n; i++) {
      total_weight = 0.0;
      val = 0.0;
      band_low= 0;
      band_high= n-1;
      if (sm->threshold>0.0) trim_bounds( i, &band_low, &band_high, 
					  data_in, sm->threshold );
      if (i-band_low > sm->last) band_low= i - sm->last;
      for (k=band_low; k<i; k++) {
	if (!(missing[k][z])) {
	  val += kernel[i-k]*data_in[k];
	  total_weight += kernel[i-k];
	}
      }
      for (k=i; k<=band_high; k++) {
	if (!(missing[k][z])) {
	  val += kernel[k-i]*data_in[k];
	  total_weight += kernel[k-i];
	} 
	if ((k-i)>=sm->last) break;
      }
      if (total_weight != 0.0)
	data_out[i]= val/total_weight;
    }
  }
  else {
    assert( missing != NULL );
    assert( sm->smDim == 'z' );
    for (i=0; i<n; i++) {
      total_weight = 0.0;
      val = 0.0;
      band_low= 0;
      band_high= n-1;
      if (sm->threshold>0.0) trim_bounds( i, &band_low, &band_high, 
					  data_in, sm->threshold );
      if (i-band_low > sm->last) band_low= i - sm->last;
      for (k=band_low; k<i; k++) {
	if (!(missing[z][k])) {
	  val += kernel[i-k]*data_in[k];
	  total_weight += kernel[i-k];
	}
      }
      for (k=i; k<=band_high; k++) {
	if (!(missing[z][k])) {
	  val += kernel[k-i]*data_in[k];
	  total_weight += kernel[k-i];
	} 
	if ((k-i)>=sm->last) break;
      }
      if (total_weight != 0.0)
	data_out[i]= val/total_weight;
    }
  }

}

static void kernel_smooth_group( Smoother* sm,
				 float** dtbl_in, float** dtbl_out, 
				 int ndata, int n,
				 unsigned char** missing, int z )
{
  int i;
  int j;
  int k;
  double weight, total_weight;
  int band_low;
  int band_high;
  double* kernel;

  check_kernel(sm, n);
  kernel= (double*)sm->data;

  if (!missing || (sm->smDim != 't' && sm->smDim != 'z')) {
    for (i=0; i<n; i++) {
      total_weight = 0.0;
      band_low= 0;
      band_high= n-1;
      if (sm->threshold>0.0) 
	trim_bounds_group( sm, i, &band_low, &band_high, dtbl_in, ndata );
      if (i-band_low > sm->last) band_low= i - sm->last;
      for (j=0; j<ndata; j++) dtbl_out[j][i]= 0.0;
      for (k=band_low; k<i; k++) {
	for (j=0; j<ndata; j++) dtbl_out[j][i] += kernel[i-k]*dtbl_in[j][k];
	total_weight += kernel[i-k];
      }
      for (k=i; k<=band_high; k++) {
	for (j=0; j<ndata; j++) dtbl_out[j][i] += kernel[k-i]*dtbl_in[j][k];
	total_weight += kernel[k-i];
	if ((k-i)>=sm->last) break;
      } 
      if (total_weight != 0.0)
	for (j=0; j<ndata; j++) dtbl_out[j][i]*= (1.0/total_weight);
    }
  }
  else if (sm->smDim == 't') {
    assert( missing != NULL );
    for (i=0; i<n; i++) {
      total_weight = 0.0;
      band_low= 0;
      band_high= n-1;
      if (sm->threshold>0.0) 
	trim_bounds_group( sm, i, &band_low, &band_high, dtbl_in, ndata );
      if (i-band_low > sm->last) band_low= i - sm->last;
      for (j=0; j<ndata; j++) dtbl_out[j][i]= 0.0;
      for (k=band_low; k<i; k++) {
	if (!(missing[k][z])) {
	  for (j=0; j<ndata; j++) dtbl_out[j][i] += kernel[i-k]*dtbl_in[j][k];
	  total_weight += kernel[i-k];
	}
      }
      for (k=i; k<=band_high; k++) {
	if (!(missing[k][z])) {
	  for (j=0; j<ndata; j++) dtbl_out[j][i] += kernel[k-i]*dtbl_in[j][k];
	  total_weight += kernel[k-i];
	} 
	if ((k-i)>=sm->last) break;
      }
      if (total_weight != 0.0)
	for (j=0; j<ndata; j++) dtbl_out[j][i]*= (1.0/total_weight);
    }
  }
  else {
    assert( missing != NULL );
    assert( sm->smDim == 'z' );
    for (i=0; i<n; i++) {
      total_weight = 0.0;
      band_low= 0;
      band_high= n-1;
      if (sm->threshold>0.0) 
	trim_bounds_group( sm, i, &band_low, &band_high, dtbl_in, ndata );
      if (i-band_low > sm->last) band_low= i - sm->last;
      for (j=0; j<ndata; j++) dtbl_out[j][i]= 0.0;
      for (k=band_low; k<i; k++) {
	if (!(missing[z][k])) {
	  for (j=0; j<ndata; j++) dtbl_out[j][i] += kernel[i-k]*dtbl_in[j][k];
	  total_weight += kernel[i-k];
	}
      }
      for (k=i; k<=band_high; k++) {
	if (!(missing[z][k])) {
	  for (j=0; j<ndata; j++) dtbl_out[j][i] += kernel[k-i]*dtbl_in[j][k];
	  total_weight += kernel[k-i];
	}
	if ((k-i)>=sm->last) break;
      }
      if (total_weight != 0.0)
	for (j=0; j<ndata; j++) dtbl_out[j][i]*= (1.0/total_weight);
    }
  }

}

static void ddx_smooth( Smoother* sm, 
			float* data_in, float* data_out, int n,
			unsigned char** missing, int z)
{
  int i;

  data_out[0]= 0.0;
  for (i=1; i<n; i++) {
    data_out[i]= data_in[i] - data_in[i-1];
  }
}

static void ddx_smooth_group( Smoother* sm,
			      float** dtbl_in, float** dtbl_out, 
			      int ndata, int n,
			      unsigned char** missing, int z )
{
  int i;
  int j;

  for (j=0; j<ndata; j++) {
    dtbl_out[j][0]= 0.0;
    for (i=1; i<n; i++)
      dtbl_out[j][i]= dtbl_in[j][i] - dtbl_in[j][i-1];
  }
}

static void linterp_smooth( Smoother* sm, 
			    float* data_in, float* data_out, int n,
			    unsigned char** missing, int z)
{
  int i;

  if (n<3) Abort("smoother: can't linearly interpolate for n<3!\n");
  for (i=1; i<n-1; i++) {
    data_out[i]= 0.5*(data_in[i-1] + data_in[i+1]);
  }
  data_out[0]= 2.0*data_in[1] - data_in[2];
  data_out[n-1]= 2.0*data_in[n-2] - data_in[n-3];
}

static void linterp_smooth_group( Smoother* sm,
				  float** dtbl_in, float** dtbl_out, 
				  int ndata, int n,
				  unsigned char** missing, int z )
{
  int i;
  int j;

  if (n<3) Abort("smoother: can't linearly interpolate for n<3!\n");
  for (j=0; j<ndata; j++) {
    for (i=1; i<n-1; i++) {
      dtbl_out[j][i]= 0.5*(dtbl_in[j][i-1] + dtbl_in[j][i+1]);
    }
    dtbl_out[j][0]= 2.0*dtbl_in[j][1] - dtbl_in[j][2];
    dtbl_out[j][n-1]= 2.0*dtbl_in[j][n-2] - dtbl_in[j][n-3];
  }
}

static void shift_smooth( Smoother* sm, 
			  float* data_in, float* data_out, int n,
			  unsigned char** missing, int z)
{
  int i;
  int shift= (int)rint(sm->bandwidth);
  int low= (shift>=0) ? 0 : -shift;
  int high= (shift>=0) ? n-(shift+1) : n-1;
  if (abs(shift)>n) 
    Abort("smoother: shift_smooth: shift is larger than range!");
  for (i=low; i<=high; i++) data_out[i]= data_in[i+shift];
  for (i=0; i<low; i++) data_out[i]= 0.0;
  for (i=high+1; i<n; i++) data_out[i]= 0.0;
}

static void shift_smooth_group( Smoother* sm,
				float** dtbl_in, float** dtbl_out, 
				int ndata, int n,
				unsigned char** missing, int z )
{
  int i;
  int j;
  int shift= (int)rint(sm->bandwidth);
  int low= (shift>=0) ? 0 : -shift;
  int high= (shift>=0) ? n-(shift+1) : n-1;
  if (abs(shift)>n) 
    Abort("smoother: shift_smooth: shift is larger than range!");

  for (j=0; j<ndata; j++) {
    for (i=low; i<=high; i++) dtbl_out[j][i]= dtbl_in[j][i+shift];
    for (i=0; i<low; i++) dtbl_out[j][i]= 0.0;
    for (i=high+1; i<n; i++) dtbl_out[j][i]= 0.0;
  }
}

static void runningsum_smooth( Smoother* sm, 
			       float* data_in, float* data_out, int n,
			       unsigned char** missing, int z)
{
  int i;
  data_out[0]= data_in[0];
  for (i=1; i<n; i++) {
    data_out[i]= data_out[i-1]+data_in[i];
  }
}

static void runningsum_smooth_group( Smoother* sm,
				     float** dtbl_in, float** dtbl_out, 
				     int ndata, int n,
				     unsigned char** missing, int z )
{
  int i;
  int j;
  for (j=0; j<ndata; j++) {
    dtbl_out[j][0]= dtbl_in[j][0];
    for (i=1; i<n; i++) {
      dtbl_out[j][i]= dtbl_out[j][i-1]+dtbl_in[j][i];
    }
  }
}

static void check_sortbuf( Smoother* sm, int n )
{
  float* sortbuf;
  int i;

  if (sm->data && n==sm->n) return;

  if (sm->data) free(sm->data);
  if (!(sortbuf=(float*)malloc(n*sizeof(float))))
    Abort("smoother:check_sortbuf: unable to allocate %d floats!\n",n);
  sm->n= n;
  sm->data= (void*)sortbuf;
}

static int compare_floats( const void *p1, const void *p2 )
{
  float* f1= (float*)p1;
  float* f2= (float*)p2;
  if (*f1<*f2) return -1;
  else if (*f1>*f2) return 1;
  else return 0;
}

static void median_smooth( Smoother* sm, 
			   float* data_in, float* data_out, int n,
			   unsigned char** missing, int z)
{
  int i;
  int half_band= (int)(rint(0.5*sm->bandwidth));
  float* sortbuf= NULL;

  check_sortbuf(sm, n);
  sortbuf= (float*)sm->data;

  for (i=0; i<n; i++) {
    int low= i-half_band;
    int high= i+(half_band-1);
    int nSort= 0;
    if (low<0) low= 0;
    if (high>=n) high=n-1;
    if (sm->threshold>0.0)
      trim_bounds( i, &low, &high, data_in, sm->threshold );
    if (!missing || ((sm->smDim != 't') && (sm->smDim != 'z'))) {
      nSort= high+1-low;
      bcopy(data_in+low, sortbuf, nSort*sizeof(float));
    }
    else {
      assert(missing!=NULL);
      if (sm->smDim=='t') {
	int j;
	for (j=low; j<=high; j++) 
	  if (!missing[j][z]) {
	    sortbuf[nSort++]= data_in[j];
	  }
      }
      else {
	int j;
	assert(sm->smDim=='z');
	for (j=low; j<=high; j++) 
	  if (!missing[z][j]) {
	    sortbuf[nSort++]= data_in[j];
	  }
      }
    }
    qsort(sortbuf, nSort, sizeof(float), compare_floats);
    data_out[i]= sortbuf[nSort/2];
  }
}

static void median_smooth_group( Smoother* sm,
				 float** dtbl_in, float** dtbl_out, 
				 int ndata, int n,
				 unsigned char** missing, int z )
{
  int i;
  int j;
  int half_band= (int)(rint(0.5*sm->bandwidth));
  float* sortbuf= NULL;

  check_sortbuf(sm, n);
  sortbuf= (float*)sm->data;

  for (j=0; j<ndata; j++) 
    for (i=0; i<n; i++) {
      int low= i-half_band;
      int high= i+(half_band-1);
      int nSort= 0;
      if (low<0) low= 0;
      if (high>=n) high=n-1;
      if (sm->threshold>0.0) 
	trim_bounds_group( sm, i, &low, &high, dtbl_in, ndata );
      if (!missing || ((sm->smDim != 't') && (sm->smDim != 'z'))) {
	nSort= high+1-low;
	bcopy(&(dtbl_in[j][low]), sortbuf, nSort*sizeof(float));
      }
      else {
	int k;
	assert(missing != NULL);
	if (sm->smDim=='t') {
	  for (k=low; k<=high; k++)
	    if (!missing[k][z]) {
	      sortbuf[nSort++]= dtbl_in[j][k];
	    }
	}
	else {
	  assert(sm->smDim=='z');
	  for (k=low; k<=high; k++)
	    if (!missing[z][k]) {
	      sortbuf[nSort++]= dtbl_in[j][k];
	    }
	}
      }
      qsort(sortbuf, nSort, sizeof(float), compare_floats);
      dtbl_out[j][i]= sortbuf[nSort/2];
    }
}

static void none_smooth( Smoother* sm, 
			 float* data_in, float* data_out, int n,
			 unsigned char** missing, int z)
{
  int i;

  for (i=0; i<n; i++) data_out[i]= data_in[i];
}

static void none_smooth_group( Smoother* sm,
			       float** dtbl_in, float** dtbl_out, 
			       int ndata, int n,
			       unsigned char** missing, int z )
{
  int i;
  int j;

  for (j=0; j<ndata; j++) 
    for (i=0; i<n; i++)
      dtbl_out[j][i]= dtbl_in[j][i];
}

Smoother* sm_create_smoother()
{
  Smoother* result;

  if ( !(result= (Smoother*)malloc(sizeof(Smoother))) )
    Abort("sm_create_smoother: unable to allocate %d bytes!\n",
	  sizeof(Smoother));

  result->type= default_type;
  result->bandwidth= default_bandwidth;
  result->k= default_k;
  result->threshold= default_threshold;
  result->thresh_test= default_thresh_test;
  result->n= 0;
  result->data= NULL;
  result->last= 0;
  result->smDim= default_sm_dim;

  switch (default_type) {
  case SM_GAUSSIAN: 
    {
      result->smooth= kernel_smooth; 
      result->smooth_group= kernel_smooth_group;
    }
  break;
  case SM_TRIANGULAR: 
    {
      result->smooth= kernel_smooth; 
      result->smooth_group= kernel_smooth_group;
    }
  break;
  case SM_POWER: 
    {
      result->smooth= kernel_smooth; 
      result->smooth_group= kernel_smooth_group;
    }
  break;
  case SM_NONE:
    {
      result->smooth= none_smooth;
      result->smooth_group= none_smooth_group;
    }
  break;
  case SM_DDX:
    {
      result->smooth= ddx_smooth;
      result->smooth_group= ddx_smooth_group;
    }
  break;
  case SM_SHIFT:
    {
      result->smooth= shift_smooth;
      result->smooth_group= shift_smooth_group;
    }
  break;
  case SM_LINTERP:
    {
      result->smooth= shift_smooth;
      result->smooth_group= shift_smooth_group;
    }
  break;
  case SM_RUNNINGSUM:
    {
      result->smooth= runningsum_smooth;
      result->smooth_group= runningsum_smooth_group;
    }
  break;
  case SM_MEDIAN:
    {
      result->smooth= median_smooth;
      result->smooth_group= median_smooth_group;
    }
  break;
  }

  return result;
}

void sm_destroy( Smoother* smoother )
{
  free(smoother->data);
  free(smoother);
}

void sm_set_params( sm_type type_in, float bandwidth_in, float k_in,
		    float thresh_in, SmootherThreshTest test_in )
{
  default_type= type_in;
  default_bandwidth= bandwidth_in;
  default_k= k_in;
  default_threshold= thresh_in;
  default_thresh_test= test_in;
}

void sm_get_params( sm_type* type_out, float* bandwidth_out, float* k_out,
		    float* thresh_out, SmootherThreshTest* test_out )
{
  if (type_out) *type_out= default_type;
  if (bandwidth_out) *bandwidth_out= default_bandwidth;
  if (k_out) *k_out= default_k;
  if (thresh_out) *thresh_out= default_threshold;
  if (test_out) *test_out= default_thresh_test;
}

void sm_set_direction( Smoother* sm, int dir )
{
  sm->smDim= dir;
}

int sm_get_direction( Smoother* sm )
{
  return sm->smDim;
}
