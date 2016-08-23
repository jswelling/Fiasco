/************************************************************
 *                                                          *
 *  mri_subsample.c                                            *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1997 Department of Statistics,         *
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
 *  Original programming by Joel Welling, 5/98              *
 ************************************************************/
/*************************************************************

  DESCRIPTION OF MRI_SUBSAMPLE

  mri_subsample takes a pgh MRI dataset of any type, and outputs
  a dataset of the same type containing subsampled data.
  The subsampled data spans the range of one of the dimensions
  of the input data, but has a smaller extent than that of
  the input dimension.  The output data is calculated by dividing
  the input data up into zones, one zone per step in the output
  extent.  Each zone thus contains one or more input data values.
  The output value for the zone is calculated by one of several
  methods;  see below.  The subsampling operation acts on all 
  chunks in the dataset.

**************************************************************/

#include <stdlib.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "misc.h"
#include "stdcrg.h"

#define KEYBUF_SIZE 512

static char rcsid[] = "$Id: mri_subsample.c,v 1.21 2007/07/06 18:45:53 welling Exp $";

typedef enum { 
  SMPL_MAX, SMPL_MIN, SMPL_MEAN, SMPL_SUM, SMPL_COUNT, SMPL_CLOSEST,
  SMPL_MEDIAN, SMPL_Q1, SMPL_Q3, SMPL_IQR
} MethodType;

typedef unsigned char uchar;

static MRI_Dataset *Input = NULL, *Output = NULL;
static char selected_dim[512]= "";
static int new_extent;
static float fwindow; /* Scaled window from command line */
static float foffset; /* Scaled offset from command line */
static int data_changed= 0;
static char* progname;
static int verbose_flag= 0;
static int reallyverbose_flag= 0;
static MethodType smpl_method= SMPL_MAX;
static unsigned char** missing= NULL; /* non-NULL if missing info present */
static int missing_dz= 0;             /* dim of missing info z */
static int missing_dt= 0;             /* dim of missing info t */
static int missing_z_stride= 0;       /* set for specific chunks */
static int missing_t_stride= 0;       /* set for specific chunks */
static int missing_relevant= 0;       /* set for specific chunks */
static void* sortbuf= NULL;

static void safe_copy(char* str1, char* str2) {
  strncpy(str1, str2, KEYBUF_SIZE);
  str1[KEYBUF_SIZE-1]= '\0';
}

static void safe_concat(char* str1, char* str2) {
  strncat(str1, str2, (KEYBUF_SIZE-strlen(str1))-1);
}

static int get_chunk_type(MRI_Dataset* ds, char* chunk)
{
  char key_buf[KEYBUF_SIZE];

  safe_copy(key_buf,chunk);
  safe_concat(key_buf,".datatype");
  if (mri_has(ds, key_buf)) {
    char* type_name= mri_get_string(ds, key_buf);
    if (!strcmp(type_name,"uint8")) return MRI_UNSIGNED_CHAR;
    else if (!strcmp(type_name,"int16")) return MRI_SHORT;
    else if (!strcmp(type_name,"int32")) return MRI_INT;
    else if (!strcmp(type_name,"int64")) return MRI_LONGLONG;
    else if (!strcmp(type_name,"float32")) return MRI_FLOAT;
    else if (!strcmp(type_name,"float64")) return MRI_DOUBLE;
    else Abort("%s: unknown data type for key %s!\n",progname,key_buf);
  }
  else Abort("%s: missing tag %s!\n",progname,key_buf);

  return 0; /* not reached */
}

static int safe_get_extent(MRI_Dataset* ds, char* chunk, char* dim)
{
  char key_buf[KEYBUF_SIZE];
  char dim_buf[4];
  dim_buf[0]= *dim;
  dim_buf[1]= '\0';
  safe_copy(key_buf,chunk);
  safe_concat(key_buf,".extent.");
  safe_concat(key_buf,dim_buf);
  if (mri_has(Input,key_buf)) return mri_get_int(Input,key_buf);
  else Abort("%s: input missing tag %s!\n",progname,key_buf);
  return 0; /* not reached */
}

static void calc_sizes(char* this_chunk, char* dimstr, char* selected,
		       int* fast_blocksize_out, int* slow_blocksize_out ) {
  /* This routine will fail if *selected is not in dimstr! */
  int fast_blocksize;
  int slow_blocksize;
  char* this_dim;

  fast_blocksize= 1;

  this_dim= dimstr;
  while (*this_dim != *selected) 
    fast_blocksize *= safe_get_extent(Input,this_chunk,this_dim++);

  this_dim++; /* step over selected dim */

  slow_blocksize= 1;
  while (*this_dim)
    slow_blocksize *= safe_get_extent(Input, this_chunk, this_dim++);

  *fast_blocksize_out= fast_blocksize;
  *slow_blocksize_out= slow_blocksize;
}

#define DEF_COMPARATOR( type ) \
static int compare_ ## type ( const void* p1, const void* p2 ) \
{ type v1= *(type*)p1; type v2= *(type*)p2; \
  if (v1<v2) return -1; else if (v1==v2) return 0; else return 1; }

DEF_COMPARATOR( uchar )
DEF_COMPARATOR( short )
DEF_COMPARATOR( int )
DEF_COMPARATOR( float )
DEF_COMPARATOR( double )

#undef DEF_COMPARATOR

#define MAKEIQR( type, compar ) { \
type* in= ibuf; type* out= obuf; unsigned char* m=mbuf; type* s; \
in += obase*stride; m += obase*stride; \
for (j=0; j<stride; j++) { \
  offset= j; s= sortbuf; \
  for (i=0; i<n; i++) { \
    if (!m[offset]) { *s++= in[offset]; } \
    offset += stride; \
  } \
  howmany= (s-(type*)sortbuf); \
  if (howmany>1) { \
    qsort( sortbuf, howmany, sizeof(type), compar ); \
    out[j]= ((type*)sortbuf)[(3*howmany)/4] - ((type*)sortbuf)[(howmany)/4]; \
  } \
  else { \
    out[j]= 0; \
  } \
}}

static void subsample_iqr( void* obuf, void* ibuf, unsigned char* mbuf,
			      int stride, int obase, int n, int type )
{
  int i;
  int j;
  int offset;
  int howmany;

  switch (type) {
  case MRI_UNSIGNED_CHAR: MAKEIQR(unsigned char, compare_uchar); break;
  case MRI_SHORT: MAKEIQR(short, compare_short); break;
  case MRI_INT: MAKEIQR(int, compare_int); break;
  case MRI_LONGLONG: MAKEIQR(long long, compare_int); break;
  case MRI_FLOAT: MAKEIQR(float, compare_float); break;
  case MRI_DOUBLE: MAKEIQR(double, compare_double); break;
  default:
    Abort("%s: switching on type: unrecognized type %d!\n",progname,type);
  }
}

#undef MAKEIQR

#define PICKSORTED( type, compar, numerator, denominator ) { \
type* in= ibuf; type* out= obuf; unsigned char* m=mbuf; type* s; \
in += obase*stride; m += obase*stride; \
for (j=0; j<stride; j++) { \
  offset= j; s= sortbuf; \
  for (i=0; i<n; i++) { \
    if (!m[offset]) { *s++= in[offset]; } \
    offset += stride; \
  } \
  howmany= (s-(type*)sortbuf); \
  if (howmany>0) { \
    qsort( sortbuf, howmany, sizeof(type), compar ); \
    out[j]= ((type*)sortbuf)[(numerator*howmany)/denominator]; \
  } \
  else { \
    out[j]= in[j]; \
  } \
}}

static void subsample_median( void* obuf, void* ibuf, unsigned char* mbuf,
			      int stride, int obase, int n, int type )
{
  int i;
  int j;
  int offset;
  int howmany;

  switch (type) {
  case MRI_UNSIGNED_CHAR: PICKSORTED(unsigned char, compare_uchar, 1, 2); break;
  case MRI_SHORT: PICKSORTED(short, compare_short, 1, 2); break;
  case MRI_INT: PICKSORTED(int, compare_int, 1, 2); break;
  case MRI_LONGLONG: PICKSORTED(long long, compare_int, 1, 2); break;
  case MRI_FLOAT: PICKSORTED(float, compare_float, 1, 2); break;
  case MRI_DOUBLE: PICKSORTED(double, compare_double, 1, 2); break;
  default:
    Abort("%s: switching on type: unrecognized type %d!\n",progname,type);
  }
}

static void subsample_q1( void* obuf, void* ibuf, unsigned char* mbuf,
			  int stride, int obase, int n, int type )
{
  int i;
  int j;
  int offset;
  int howmany;

  switch (type) {
  case MRI_UNSIGNED_CHAR: PICKSORTED(unsigned char, compare_uchar, 1, 4); break;
  case MRI_SHORT: PICKSORTED(short, compare_short, 1, 4); break;
  case MRI_INT: PICKSORTED(int, compare_int, 1, 4); break;
  case MRI_LONGLONG: PICKSORTED(long long, compare_int, 1, 4); break;
  case MRI_FLOAT: PICKSORTED(float, compare_float, 1, 4); break;
  case MRI_DOUBLE: PICKSORTED(double, compare_double, 1, 4); break;
  default:
    Abort("%s: switching on type: unrecognized type %d!\n",progname,type);
  }
}

static void subsample_q3( void* obuf, void* ibuf, unsigned char* mbuf,
			  int stride, int obase, int n, int type )
{
  int i;
  int j;
  int offset;
  int howmany;

  switch (type) {
  case MRI_UNSIGNED_CHAR: PICKSORTED(unsigned char, compare_uchar, 3, 4); break;
  case MRI_SHORT: PICKSORTED(short, compare_short, 3, 4); break;
  case MRI_INT: PICKSORTED(int, compare_int, 3, 4); break;
  case MRI_LONGLONG: PICKSORTED(long long, compare_int, 3, 4); break;
  case MRI_FLOAT: PICKSORTED(float, compare_float, 3, 4); break;
  case MRI_DOUBLE: PICKSORTED(double, compare_double, 3, 4); break;
  default:
    Abort("%s: switching on type: unrecognized type %d!\n",progname,type);
  }
}

#undef PICKSORTED

#define MAXIMIZE( type, accumtype ) { \
type* in= ibuf; type* out= obuf; accumtype val; unsigned char* m=mbuf; \
int init; \
in += obase*stride; m += obase*stride; \
for (j=0; j<stride; j++) { \
  init=1; \
  offset= j; \
  for (i=0; i<n; i++) { \
    if (!m[offset]) { \
      if (init) { init=0; val=in[j]; } \
      else if (in[offset]>val) val= in[offset]; \
    } \
    offset += stride; \
  } \
  out[j]= (init ? in[j] : val); \
}}

static void subsample_max( void* obuf, void* ibuf, unsigned char* mbuf,
			   int stride, int obase, int n, int type )
{
  int i;
  int j;
  int offset;

  switch (type) {
  case MRI_UNSIGNED_CHAR: MAXIMIZE(unsigned char,int); break;
  case MRI_SHORT: MAXIMIZE(short,int); break;
  case MRI_INT: MAXIMIZE(int,int); break;
  case MRI_LONGLONG: MAXIMIZE(long long,long long); break;
  case MRI_FLOAT: MAXIMIZE(float,float); break;
  case MRI_DOUBLE: MAXIMIZE(double,double); break;
  default:
    Abort("%s: switching on type: unrecognized type %d!\n",progname,type);
  }
}

#undef MAXIMIZE

#define MINIMIZE( type, accumtype ) { \
type* in= ibuf; type* out= obuf; accumtype val; unsigned char* m=mbuf; \
int init; \
in += obase*stride; m += obase*stride; \
for (j=0; j<stride; j++) { \
  init=1; \
  offset= j; \
  for (i=0; i<n; i++) { \
    if (!m[offset]) { \
      if (init) { init=0; val=in[j]; } \
      else if (in[offset]<val) val= in[offset]; \
    } \
    offset += stride; \
  } \
  out[j]= (init ? in[j] : val); \
}}

static void subsample_min( void* obuf, void* ibuf, unsigned char* mbuf,
			   int stride, int obase, int n, int type )
{
  int i;
  int j;
  int offset;

  switch (type) {
  case MRI_UNSIGNED_CHAR: MINIMIZE(unsigned char,int); break;
  case MRI_SHORT: MINIMIZE(short,int); break;
  case MRI_INT: MINIMIZE(int,int); break;
  case MRI_LONGLONG: MINIMIZE(long long,long long); break;
  case MRI_FLOAT: MINIMIZE(float,float); break;
  case MRI_DOUBLE: MINIMIZE(double,double); break;
  default:
    Abort("%s: switching on type: unrecognized type %d!\n",progname,type);
  }
}

#undef MINIMIZE

#define SUM( type, accumtype ) { \
type* in= ibuf; type* out= obuf; accumtype val; unsigned char* m=mbuf; \
int init; \
in += obase*stride; m += obase*stride; \
for (j=0; j<stride; j++) { \
  init=1; \
  offset= j; \
  for (i=0; i<n; i++) { \
    if (!m[offset]) { \
      if (init) { init=0; val=in[j]; } \
      else val += in[offset]; \
    } \
    offset += stride; \
  } \
  out[j]= (init ? in[j] : val); \
}}

static void subsample_sum( void* obuf, void* ibuf, unsigned char* mbuf,
			   int stride, int obase, int n, int type )
{
  int i;
  int j;
  int offset;

  switch (type) {
  case MRI_UNSIGNED_CHAR: SUM(unsigned char,int); break;
  case MRI_SHORT: SUM(short,int); break;
  case MRI_INT: SUM(int,int); break;
  case MRI_LONGLONG: SUM(long long,long long); break;
  case MRI_FLOAT: SUM(float,double); break;
  case MRI_DOUBLE: SUM(double,double); break;
  default:
    Abort("%s: switching on type: unrecognized type %d!\n",progname,type);
  }
}

#undef SUM
#define MEAN( type, accumtype ) { \
type* in= ibuf; type* out= obuf; accumtype val; unsigned char* m=mbuf; \
int init; int count;\
in += obase*stride; m += obase*stride; \
for (j=0; j<stride; j++) { \
  init=1; \
  offset= j; \
  count=0; \
  for (i=0; i<n; i++) { \
    if (!m[offset]) { \
      count += 1; \
      if (init) { init=0; val=in[j]; } \
      else val += in[offset]; \
    } \
    offset += stride; \
  } \
  out[j]= (init ? in[j] : (type)((double)val/(double)count)); \
}}

static void subsample_mean( void* obuf, void* ibuf, unsigned char* mbuf,
			    int stride, int obase, int n, int type )
{
  int i;
  int j;
  int offset;

  switch (type) {
  case MRI_UNSIGNED_CHAR: MEAN(unsigned char,int); break;
  case MRI_SHORT: MEAN(short,int); break;
  case MRI_INT: MEAN(int,int); break;
  case MRI_LONGLONG: MEAN(long long,long long); break;
  case MRI_FLOAT: MEAN(float,double); break;
  case MRI_DOUBLE: MEAN(double,double); break;
  default:
    Abort("%s: switching on type: unrecognized type %d!\n",progname,type);
  }
}

#undef MEAN

#define COUNT( type, accumtype ) { \
type* in= ibuf; type* out= obuf; accumtype val; unsigned char* m=mbuf; \
int count;\
in += obase*stride; m += obase*stride; \
for (j=0; j<stride; j++) { \
  offset= j; \
  count=0; \
  for (i=0; i<n; i++) { \
    if (!m[offset]) { \
      count += 1; \
    } \
    offset += stride; \
  } \
  out[j]= (type)count; \
}}

static void subsample_count( void* obuf, void* ibuf, unsigned char* mbuf,
			     int stride, int obase, int n, int type )
{
  int i;
  int j;
  int offset;

  switch (type) {
  case MRI_UNSIGNED_CHAR: COUNT(unsigned char,int); break;
  case MRI_SHORT: COUNT(short,int); break;
  case MRI_INT: COUNT(int,int); break;
  case MRI_LONGLONG: COUNT(long long,long long); break;
  case MRI_FLOAT: COUNT(float,double); break;
  case MRI_DOUBLE: COUNT(double,double); break;
  default:
    Abort("%s: switching on type: unrecognized type %d!\n",progname,type);
  }
}

#undef COUNT

#define CLOSEST( type ) { \
type* in= ibuf; type* out= obuf; type val; \
for (j=0; j<stride; j++) { \
  out[j]= in[(i*stride)+j]; \
} }

static void subsample_closest( void* obuf, void* ibuf, unsigned char* mbuf,
			       int stride, int ibot, int n, float lower, 
			       float ratio, int type )
{
  int i;
  int j;

  i= (int)(lower + (0.5*ratio) + 0.5) - ibot; /* rounded center of region */

  switch (type) {
  case MRI_UNSIGNED_CHAR: CLOSEST(unsigned char); break;
  case MRI_SHORT: CLOSEST(short); break;
  case MRI_INT: CLOSEST(int); break;
  case MRI_LONGLONG: CLOSEST(long long); break;
  case MRI_FLOAT: CLOSEST(float); break;
  case MRI_DOUBLE: CLOSEST(double); break;
  default:
    Abort("%s: switching on type: unrecognized type %d!\n",progname,type);
  }
}

#undef CLOSEST

static void subsample_missing( void* obuf, void* ibuf, unsigned char* mbuf,
			       int stride, int ibot, int n, 
			       int type, float lower, float ratio,
			       int ibase, int iwindow )
{
  /* The missing chunk is 1 if anything is missing and 0 otherwise.
   * Given this, it's easy to construct its appropriate value for
   * any given operation.
   */
  switch (smpl_method) {
  case SMPL_MAX:
  case SMPL_MIN:
  case SMPL_MEAN:
  case SMPL_SUM:
  case SMPL_COUNT:
  case SMPL_MEDIAN:
  case SMPL_Q1:
  case SMPL_Q3:
  case SMPL_IQR:
    /* At least *some* data will probably be there to contribute */
    subsample_min( obuf, ibuf, mbuf, stride, ibase, iwindow, type );
    break;
  case SMPL_CLOSEST:
    /* Either the intended datum is present or it is not. */
    subsample_closest( obuf, ibuf, mbuf, stride, ibot, n, lower, ratio, type );
    break;
  }
}

static void subsample( void* obuf, void* ibuf, unsigned char* mbuf,
		       int stride, int ibot, int n, 
		       int type, float lower, float ratio,
		       int ibase, int iwindow )
{
  switch (smpl_method) {
  case SMPL_MAX: {
    subsample_max( obuf, ibuf, mbuf, stride, ibase, iwindow, type );
  }
  break;
  
  case SMPL_MIN: {
    subsample_min( obuf, ibuf, mbuf, stride, ibase, iwindow, type );
  }
  break;
  
  case SMPL_MEAN: {
    subsample_mean( obuf, ibuf, mbuf, stride, ibase, iwindow, type );
  }
  break;
  
  case SMPL_SUM: {
    subsample_sum( obuf, ibuf, mbuf, stride, ibase, iwindow, type );
  }
  break;
  
  case SMPL_COUNT: {
    subsample_count( obuf, ibuf, mbuf, stride, ibase, iwindow, type );
  }
  break;
  
  case SMPL_CLOSEST: {
    subsample_closest( obuf, ibuf, mbuf, stride, ibot, n, lower, ratio, type );
  }
  break;
  
  case SMPL_MEDIAN: {
    subsample_median( obuf, ibuf, mbuf, stride, ibase, iwindow, type );
  }
  break;
  
  case SMPL_Q1: {
    subsample_q1( obuf, ibuf, mbuf, stride, ibase, iwindow, type );
  }
  break;
  
  case SMPL_Q3: {
    subsample_q3( obuf, ibuf, mbuf, stride, ibase, iwindow, type );
  }
  break;

  case SMPL_IQR: {
    subsample_iqr( obuf, ibuf, mbuf, stride, ibase, iwindow, type );
  }

  break;
  
  }
}

static void missing_setup(char* this_chunk) {
  char* dimstr;
  char key_buf[KEYBUF_SIZE];
  char* chunk_t;
  char* chunk_z;

  if (!missing || !strcmp(this_chunk,"missing")
      || (selected_dim[0] != 't' && selected_dim[0] != 'z')) {
    missing_relevant= 0;
    return;
  }

  safe_copy(key_buf, this_chunk);
  safe_concat(key_buf, ".dimensions");
  dimstr= mri_get_string(Input,key_buf);
  chunk_t= strchr(dimstr,'t');
  chunk_z= strchr(dimstr,'z');

  if (chunk_t==NULL || chunk_z==NULL) {
    missing_relevant= 0;
    return;
  }

  if (strchr(dimstr,'t') && strchr(dimstr,'z')) {
    int junk;
    missing_relevant= 1;
    calc_sizes(this_chunk, dimstr, "z", &missing_z_stride, &junk);
    calc_sizes(this_chunk, dimstr, "t", &missing_t_stride, &junk);
  }
  else missing_relevant= 0;
}

static void missing_ready(unsigned char* mbuf, const int m_size,
			  int *slow_cycle, int *fast_cycle)
{
  int i,j;
  int zTmp, tTmp;
  int zReps, tReps;
  unsigned char* here;

  if (!missing_relevant) {
    for (i=0; i<m_size; i++) mbuf[i]= 0;
  }
  else {
    here= mbuf;
    if (selected_dim[0]=='t') {
      if (missing_z_stride>missing_t_stride) {
	/* Subsampling in t; z varies more slowly.  We have
	 * to go through many cycles of t before z increments.
	 */
	tReps= missing_t_stride;
	zReps= missing_z_stride/(missing_t_stride*missing_dt);
	zTmp= *slow_cycle;
	for (tTmp=0; tTmp<missing_dt; tTmp++) {
	  for (i=0; i<tReps; i++) *here++= missing[tTmp][zTmp];
	}
	*fast_cycle += 1;
	if (*fast_cycle>=zReps) {
	  *fast_cycle= 0;
	  *slow_cycle += 1;
	}
      }
      else {
	/* Subsampling in t; z varies faster.
	 * Hence both t and z start out 0 for each run.
	 */
	tReps= missing_t_stride/(missing_z_stride*missing_dz);
	zReps= missing_z_stride;
	for (tTmp=0; tTmp<missing_dt; tTmp++) {
	  for (i=0; i<tReps; i++) {
	    for (zTmp=0; zTmp<missing_dz; zTmp++) {
	      for (j=0; j<zReps; j++) *here++= missing[tTmp][zTmp];
	    }
	  }
	}

	*slow_cycle= *fast_cycle= 0;
      }
    }
    else if (selected_dim[0]=='z') {
      if (missing_t_stride>missing_z_stride) {
	/* Subsampling in z; t varies more slowly.  We have
	 * to go through many cycles of z before t increments.
	 */
	zReps= missing_z_stride;
	tReps= missing_t_stride/(missing_z_stride*missing_dz);
	tTmp= *slow_cycle;
	for (zTmp=0; zTmp<missing_dz; zTmp++) {
	  for (i=0; i<zReps; i++) *here++= missing[tTmp][zTmp];
	}
	*fast_cycle += 1;
	if (*fast_cycle>=tReps) {
	  *fast_cycle= 0;
	  *slow_cycle += 1;
	}

      }
      else {
	/* Subsampling in z; t varies faster.
	 * Hence both t and z start out 0 for each run.
	 */
	zReps= missing_z_stride/(missing_t_stride*missing_dt);
	tReps= missing_t_stride;
	for (zTmp=0; zTmp<missing_dz; zTmp++) {
	  for (i=0; i<zReps; i++) {
	    for (tTmp=0; tTmp<missing_dt; tTmp++) {
	      for (j=0; j<tReps; j++) *here++= missing[tTmp][zTmp];
	    }
	  }
	}

	*slow_cycle= *fast_cycle= 0;
      }
    }
  }
}

static void transfer_data(char* this_chunk,
			  int fast_blksize, int slow_blksize, 
			  int selected_extent, int new_extent)
{
  int in_offset= 0;
  int out_offset= 0;
  int in_framestart= 0;
  int ifast;
  int islow;
  int type;
  int typesize;
  void* obuf= NULL;
  void* ibuf= NULL;
  unsigned char* mbuf= NULL;
  int max_old_per_new;
  double ratio;
  double in_lower;
  double epsilon;
  int ibot;
  int itop;
  int irange;
  int iwindow;
  int ibase;
  int missing_offset;
  int missing_flag= 0;  /* True if this *is* the 'missing' chunk */
  int missing_slow_cycle= 0; /* some offset info for accessing missing */
  int missing_fast_cycle= 0; /* ditto */

  type= get_chunk_type(Input,this_chunk);
  typesize= get_typesize(Input,this_chunk);
  missing_setup(this_chunk);

  missing_flag= (!strcmp(this_chunk,"missing")); /* special treatment */

  epsilon= 1.0/selected_extent;

  if (!(obuf= (void*)malloc(fast_blksize*typesize)))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname, fast_blksize*typesize);
  if (!(mbuf= (unsigned char*)malloc(fast_blksize*selected_extent)))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname, fast_blksize*selected_extent);
  if ((smpl_method==SMPL_MEDIAN) || (smpl_method==SMPL_Q1)
      || (smpl_method==SMPL_Q3) || (smpl_method==SMPL_IQR)) {
    if (!(sortbuf= (void*)malloc(selected_extent*typesize)))
      Abort("%s: unable to allocate %d bytes!\n",selected_extent*typesize);
  }
  else sortbuf= NULL;

  ratio= ((double)selected_extent)/((double)new_extent);

  missing_slow_cycle= missing_fast_cycle= 0;
  for (islow=0; islow<slow_blksize; islow++) {
    in_offset= in_framestart;
    in_lower= 0.0-epsilon;
    missing_ready(mbuf, fast_blksize*selected_extent,
		  &missing_slow_cycle, &missing_fast_cycle);
    missing_offset= 0;
    for (ifast=0; ifast<new_extent; ifast++) {
      ibot= ceil(in_lower);
      itop= floor(in_lower+ratio);
      irange= itop - ibot + 1;
      iwindow= (irange+0.5)*fwindow;
      ibase= (irange+0.5)*foffset;
      if (reallyverbose_flag) {
	fprintf(stderr,
                "block %d %d: %f %f -> range %d to %d, wndw %d at %d:\n",
		islow, ifast, in_lower, ratio, ibot, itop, iwindow, ibase);
	fprintf(stderr,
		"	      %d at %d -> %d at %d\n",
		fast_blksize*irange, in_offset,
		fast_blksize, out_offset);
      }
      ibuf= mri_get_chunk(Input, this_chunk, fast_blksize*irange,
			    in_offset, type);
      in_offset += fast_blksize*irange;
      if (missing_flag) 
	subsample_missing( obuf, ibuf, mbuf+missing_offset, fast_blksize, 
			   ibot, irange, type, in_lower, ratio, 
			   ibase, iwindow );
      else 
	subsample( obuf, ibuf, mbuf+missing_offset, fast_blksize, 
		   ibot, irange, type, in_lower, ratio, ibase, iwindow );
      mri_set_chunk( Output, this_chunk, fast_blksize, out_offset, type,
		     obuf );
      out_offset += fast_blksize;
      missing_offset += fast_blksize*irange;
      in_lower += ratio;
    }
    in_framestart += fast_blksize*selected_extent;
  }

  free(obuf);
  free(mbuf);
  if (sortbuf) free(sortbuf);
}

static void subsample_chunk(char* this_chunk) {
  char key_buf[KEYBUF_SIZE];
  char* dimstr;
  int selected_extent;

  if (reallyverbose_flag) fprintf(stderr,"subsampling chunk <%s>\n",this_chunk);
  safe_copy(key_buf, this_chunk);
  safe_concat(key_buf, ".dimensions");
  if (mri_has(Input,key_buf)) {
    dimstr= mri_get_string(Input,key_buf);
    if (reallyverbose_flag) fprintf(stderr,"dimstr %s\n",dimstr);
    if (strchr(dimstr,*selected_dim)) {
      safe_copy(key_buf, this_chunk);
      safe_concat(key_buf, ".extent.");
      safe_concat(key_buf, selected_dim);
      if (mri_has(Input,key_buf)) 
	selected_extent= mri_get_int(Input,key_buf);
      else Abort("%s: input missing tag %s!\n",progname,key_buf);
      if (selected_extent < new_extent)
	Abort("%s: selected dim in chunk %s is too small to subsample!",
	      progname,this_chunk);
      if (((fwindow*selected_extent)/new_extent)<1.0) 
	Abort("%s: window is too small to use with chunk %s!\n",
	      progname,this_chunk);
      if (reallyverbose_flag) 
	fprintf(stderr,"selected dim extent on input %d\n",selected_extent);
      if ((new_extent != selected_extent) 
	  || (new_extent == selected_extent && smpl_method==SMPL_COUNT)) {
	int fast_blksize;
	int slow_blksize;

	data_changed= 1;
	mri_set_int(Output, key_buf, new_extent);
	calc_sizes(this_chunk, dimstr, selected_dim,
		   &fast_blksize, &slow_blksize);
	transfer_data(this_chunk, fast_blksize, slow_blksize, 
		      selected_extent, new_extent);
      }
    }
    else {
      /* Chunk copied correctly in initial dataset copy */
    }
  }
}

int main( int argc, char* argv[] ) 
{
  char infile[512], outfile[512];
  char* this_key;
  char this_chunk[KEYBUF_SIZE];
  int min_flag= 0;
  int max_flag= 0;
  int mean_flag= 0;
  int sum_flag= 0;
  int count_flag= 0;
  int closest_flag= 0;
  int median_flag= 0;
  int q1_flag= 0;
  int q3_flag= 0;
  int iqr_flag= 0;
  float window_in;
  float offset_in;
  float base_in;

  progname= argv[0];

  /* Print version number */
  if (verbose_flag) Message( "# %s\n", rcsid );

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Deprecate old options */

  if (cl_present( "e" ))
     Abort ("Option e(extent) has been replaced by length|len|l.  Please see help file.\n");

  if (cl_present( "o" ))
     Abort ("Option o(offset) has been replaced by shift|shi|s.  Please see help file.\n");

  if (cl_present( "w" ))
     Abort ("Option w has been expanded to window|win.  Please see help file.\n");

  if (cl_present( "q1" ))
     Abort ("Option q1(first quartile) has been replaced by 1st_quartile|1qr.  Please see help file.\n");

  if (cl_present( "q3" ))
     Abort ("Option q3(third quartile) has been replaced by 3rd_quartile|3qr.  Please see help file.\n");


  /* Get filenames */
  if (!cl_get("dimension|dim|d", "%option %s[t]",selected_dim)) {
    fprintf(stderr,"%s: dimension to be subsampled not given.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  if (strlen(selected_dim)>1) {
    fprintf(stderr,"%s: Selected dim name must be 1 char long.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  if (!cl_get("length|len|l", "%option %d",&new_extent)) {
    fprintf(stderr,"%s: new extent not given.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  cl_get("b","%option %f[1.0]",&base_in);
  if (!cl_get("window|win","%option %f",&window_in)) window_in= base_in;
  cl_get("shift|shi|s","%option %f[0.0]",&offset_in);
  fwindow= window_in/base_in;
  foffset= offset_in/base_in;
  if (fwindow<0.0 || fwindow>1.0) {
    fprintf(stderr,"%s: invalid window or base given.\n",argv[0]);
    Help("usage");
  }
  if (foffset<0.0 || foffset>(1.0-fwindow)) {
    fprintf(stderr,"%s: offset negative or too large for window.\n",argv[0]);
    Help("usage");
  }
  if (cl_present("min")) min_flag= 1;
  if (cl_present("max")) max_flag= 1;
  if (cl_present("mean")) mean_flag= 1;
  if (cl_present("sum")) sum_flag= 1;
  if (cl_present("count|cnt")) count_flag= 1;
  if (cl_present("closest|cls")) closest_flag= 1;
  if (cl_present("median|med")) median_flag= 1;
  if (cl_present("1st_quartile|1qr")) q1_flag= 1;
  if (cl_present("3rd_quartile|3qr")) q3_flag= 1;
  if (cl_present("inter_quartile|iqr")) iqr_flag= 1;
  if (min_flag+max_flag+mean_flag+sum_flag+count_flag+closest_flag
      +median_flag + q1_flag + q3_flag + iqr_flag > 1) {
    fprintf(stderr,"%s: mutually exclusive method flags given.\n",argv[0]);
    Help( "usage" );
    exit(-1);
  }
  else if (min_flag+max_flag+mean_flag+sum_flag+count_flag+closest_flag
	   + median_flag + q1_flag + q3_flag + iqr_flag< 1)
    max_flag= 1;
  if (closest_flag && ((fwindow != 1.0) || (foffset != 0.0))) {
    fprintf(stderr,
	    "%s: window, offset, and base are incompatible with -closest.\n",
	    argv[0]);
    Help( "usage" );
    exit(-1);
  }
    
  if (max_flag) smpl_method= SMPL_MAX;
  else if (min_flag) smpl_method= SMPL_MIN;
  else if (mean_flag) smpl_method= SMPL_MEAN;
  else if (sum_flag) smpl_method= SMPL_SUM;
  else if (count_flag) smpl_method= SMPL_COUNT;
  else if (closest_flag) smpl_method= SMPL_CLOSEST;
  else if (median_flag) smpl_method= SMPL_MEDIAN;
  else if (q1_flag) smpl_method= SMPL_Q1;
  else if (q3_flag) smpl_method= SMPL_Q3;
  else if (iqr_flag) smpl_method= SMPL_IQR;

  verbose_flag= cl_present("v|verbose");
  reallyverbose_flag= cl_present("V|debug");
  if (reallyverbose_flag) verbose_flag= 1;

  if (!cl_get("", "%s", infile)) {
    fprintf(stderr,"%s: Input file name not given.\n",argv[0]);
    Help( "usage" );
    exit(-1);
  }
  if (!cl_get("", "%s", outfile)) {
    fprintf(stderr,"%s: Output file name not given.\n",argv[0]);
    Help( "usage" );
    exit(-1);
  }

  if (cl_cleanup_check()) {
    int i;
    fprintf(stderr,"%s: invalid argument in command line:\n    ",argv[0]);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(-1);
  }

  /*** End command-line parsing ***/  

  /* Open input and output datasets */
  if( !strcmp( infile, outfile ) )
    Abort( "Input and output files must be distinct." );
  Input = mri_open_dataset( infile, MRI_READ );
  Output = mri_copy_dataset( outfile, Input );
  hist_add_cl( Output, argc, argv );

  /* Grab the missing info from the input dataset if present. */
  if( mri_has( Input, "missing" ) ) {
    if (mri_has( Input, "missing.size" ) &&
	mri_has( Input, "missing.dimensions" ) &&
	mri_has( Input, "missing.datatype" ) &&
	mri_has( Input, "missing.extent.t" ) &&
	mri_has( Input, "missing.extent.z" ) &&
	( ( mri_get_int( Input, "missing.extent.t" ) *
	    mri_get_int( Input, "missing.extent.z" ) ) ==
	  mri_get_int( Input, "missing.size" ) ) &&
	!strcmp( mri_get_string( Input, "missing.dimensions" ), "zt" ) &&
	!strcmp( mri_get_string( Input, "missing.datatype" ), "uint8" ) ) {
      missing_dz= mri_get_int(Input,"missing.extent.z");
      missing_dt= mri_get_int(Input,"missing.extent.t");
      missing= Matrix(missing_dt, missing_dz, unsigned char);
      memcpy( *missing, (unsigned char *)
	      mri_get_chunk( Input, "missing", 
			     (int) ( missing_dt * missing_dz ),
			     0, MRI_UNSIGNED_CHAR ),
	      (long) ( missing_dt * missing_dz ) );
    }
    else Abort("%s: missing info present but invalid!\n",argv[0]);
  }
  else missing= NULL;

  /* Walk through the input handling each chunk in turn */
  mri_iterate_over_keys(Input);
  while ((this_key= mri_next_key(Input)) != NULL) {
    if (strlen(this_key)>=KEYBUF_SIZE) 
      Abort("%s: key too long!\n",argv[0]);
    if (!strcmp(mri_get_string(Input,this_key),"[chunk]")) {
      strncpy(this_chunk, this_key, KEYBUF_SIZE);
      this_chunk[KEYBUF_SIZE-1]= '\0';
      subsample_chunk(this_chunk);
    }
  }

  /* Write and close data-sets */
  mri_close_dataset( Input );
  mri_close_dataset( Output );

  /* Clean up */
  if (missing) FreeMatrix(missing);
  
  if (verbose_flag) Message( "#      Subsampling complete.\n" );
  if (!data_changed) 
    Message("#      Warning: input and output datasets identical!\n");

  return 0;
}

