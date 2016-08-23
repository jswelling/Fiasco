/************************************************************
 *                                                          *
 *  algorithm.c                                             *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1995 Department of Statistics,         *
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
 *  Original programming by Joel Welling 3/00               *
 *      10/02: Parallelization, Jenn Bakal                  *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#ifdef DARWIN
#include <sys/time.h>
#endif
#include <time.h>
#include <math.h>
#include <limits.h>
#include <sys/resource.h>
#include <unistd.h>
#include "../fmri/lapack.h"
#include "mri.h"
#include "fmri.h"
#include "par.h"
#include "stdcrg.h"
#include "misc.h"
#include "array.h"
#include "acct.h"
#include "estireg_utils.h"
#include "algorithm.h"

static char rcsid[] = "$Id: algorithm.c,v 1.8 2007/03/21 23:47:26 welling Exp $";

/* Several mutual info histograms may be used; pick some sizes. */
static int mutualInfoNBins[N_MUTUAL_INFO_CONTEXTS]= { 241 };
#ifdef never
static int mutualInfoNBins[N_MUTUAL_INFO_CONTEXTS]= 
  { 211, 229, 241, 263, 277 };
#endif

int algNeedsMask(Algorithm* alg)
{
  if (alg->objective_method==OBJECTIVE_MI) {
    if (alg->weight_method != WEIGHT_CONST) return 1;
    else return 0;
  }
  else return 0;
}

void algSetDebugLevel( Algorithm* alg, const int lvl )
{
  int i;
  alg->debugLevel= lvl;
  alg->opt->setDebugLevel(alg->opt, lvl); 
  for (i=0; i<N_MUTUAL_INFO_CONTEXTS; i++)
    if (alg->mutualInfoContext[i] != NULL) {
      if (lvl) ent_setMIVerbose(alg->mutualInfoContext[i],(lvl!=0));
      if (lvl>1) ent_setMIDebug(alg->mutualInfoContext[i],(lvl!=0));
    }
}

int algGetDebugLevel( Algorithm* alg )
{
  return alg->debugLevel;
}

int algNeedsStdv( Algorithm* alg )
{
  switch (alg->weight_method) {
  case WEIGHT_CONST: return 0;
  case WEIGHT_INV_STDV: return 1;
  case WEIGHT_ALIGN: return 0;
  case WEIGHT_SMTHALIGN: return 0;
  case WEIGHT_INVALID: return 0; /* well, it doesn't! */
  default: return 0;
  }
}

void algParseOpts( Algorithm* alg )
{
  sm_parse_cl_opts();
}

static char* objectiveMethodName( ObjectiveMethod mthd )
{
  switch (mthd) {
  case OBJECTIVE_MSE: return "mse";
  case OBJECTIVE_MI:  return "mutualinfo";
  case OBJECTIVE_JE:  return "jointentropy";
  default: return NULL;
  }
}

static char* searchMethodName( SearchMethod mthd )
{
  switch (mthd) {
  case SEARCH_TRILIN: return "trilin"; 
  case SEARCH_SHEAR4_FFT: return "shear4"; 
  case SEARCH_SHEAR7_FFT: return "shear7"; 
  case SEARCH_SHEAR13_FFT: return "shear13"; 
  case SEARCH_INVALID: return "invalid interpolation!"; 
  default: return NULL;
  }
}

static char* weightMethodName( WeightMethod mthd )
{
  switch (mthd) {
  case WEIGHT_CONST: return "const";
  case WEIGHT_INV_STDV: return "inv-stdv";
  case WEIGHT_ALIGN: return "align";
  case WEIGHT_SMTHALIGN: return "smoothalign";
  case WEIGHT_INVALID: return "invalid weighting!";
  default: return NULL;
  }
}

static char* qualMethodName( int mthd )
{
  switch (mthd) {
  case FR3D_QUAL_COX: return "cox";
  case FR3D_QUAL_SUM_ABS: return "sabs";
  case FR3D_QUAL_SUM_SQR: return "ssqr";
  case FR3D_QUAL_UNIT_CELL: return "ucell";
  case FR3D_QUAL_UNSET: return "invalid quality measure!";
  default: return NULL;
  }
}

static ObjectiveMethod parseObjectiveMethod( char* s )
{
  if (!strcasecmp(s,"mse")) return OBJECTIVE_MSE;
  else if (!strcasecmp(s,"mutualinfo")) return OBJECTIVE_MI;
  else if (!strcasecmp(s,"jointentropy")) return OBJECTIVE_JE;
  else return OBJECTIVE_INVALID;
}

static Optimizer* parseOptMethod(char* s, const int ndim)
{
  /* Two-sided 95% T cutoffs for DF up to 20 */
  static float tcrit[]= 
    { 0.0, 12.706176098809, 4.30263617020413,
      3.18244545308107, 2.77644502408949, 2.57058182376418,
      2.44680441622004, 2.36457539633475, 2.30597941403506, 
      2.26214359559273, 2.2281309162355, 2.20098027417758,
      2.17880969240071, 2.16036657035061, 2.14478525936996, 
      2.13144854260425, 2.11990457994238, 2.10981505255424, 
      2.10092165066567, 2.09302376167972, 2.08596322490106
  };

  if (ndim> (sizeof(tcrit)/sizeof(float))-1)
    Abort("algorithm:parseOptMethod: tcrit for %d degrees of freedom is not compiled in!\n",
          ndim);

  if (!strcasecmp(s,"nelmin"))
    return createNelminOptimizer(0.02,1.0,12);
  else if (!strcasecmp(s,"nelmin_t"))
    return createNelminTOptimizer(0.02,1.0,tcrit[ndim],12);
  else if (!strcasecmp(s,"praxis"))
    return createPraxisOptimizer(0.000001,0.1);
  else if (!strcasecmp(s,"none")) return createNoneOptimizer();
  else return NULL;
}

static SearchMethod parseSearchMethod(char* s)
{
  if (!strcasecmp(s,"trilin")) return SEARCH_TRILIN;
  else if (!strcasecmp(s,"shear4")) return SEARCH_SHEAR4_FFT;
  else if (!strcasecmp(s,"shear7")) return SEARCH_SHEAR7_FFT;
  else if (!strcasecmp(s,"shear13")) return SEARCH_SHEAR13_FFT;
  else return SEARCH_INVALID;
}

static WeightMethod parseWeightMethod(char* s)
{
  if (!strcasecmp(s,"const")) return WEIGHT_CONST;
  else if (!strcasecmp(s,"inv-stdv")) return WEIGHT_INV_STDV;  
  else if (!strcasecmp(s,"align")) return WEIGHT_ALIGN;
  else if (!strcasecmp(s,"smoothalign")) return WEIGHT_SMTHALIGN;
  else return WEIGHT_INVALID;
}

static int parseQualMethod(char* s)
{
  if (!strcasecmp(s,"cox")) 
    fshrot3d_set(FR3D_QUAL_MEASURE, FR3D_QUAL_COX);
  else if (!strcasecmp(s,"sabs")) 
    fshrot3d_set(FR3D_QUAL_MEASURE, FR3D_QUAL_SUM_ABS);
  else if (!strcasecmp(s,"ssqr")) 
    fshrot3d_set(FR3D_QUAL_MEASURE, FR3D_QUAL_SUM_SQR);
  else if (!strcasecmp(s,"ucell")) 
    fshrot3d_set(FR3D_QUAL_MEASURE, FR3D_QUAL_UNIT_CELL);
  else return 0;

  return 1;
}

int algParseInfoString( Algorithm* alg, const char* string )
{
  char* work= strdup(string); /* make guaranteed-writable copy */
  char* tok= NULL;
  char* optString= NULL;
  char* tptr;
  double tol= 0.0;
  int tolSet= 0;
  double scale= 0.0;
  int scaleSet= 0;
  int i;

  work= strdup(string); 
  tok= strtok_r(work, " ,+&:;",&tptr);
  while (tok != NULL) {
    if (!strcasecmp(tok,"smooth")) {
      alg->smooth_flag= 1;
    }
    else if (!strcasecmp(tok,"nosmooth")) {
      alg->smooth_flag= 0;
    }
    else if (!strcasecmp(tok,"inplane")) {
      alg->inplane_flag= 1;
    }
    else if (!strcasecmp(tok,"noinplane")) {
      alg->inplane_flag= 0;
    }
    else if (!strcasecmp(tok,"rotonly")) {
      alg->rot_only_flag= 1;
    }
    else if (!strcasecmp(tok,"norotonly")) {
      alg->rot_only_flag= 0;
    }
    else if (!strcasecmp(tok,"transonly")) {
      alg->trans_only_flag= 1;
    }
    else if (!strcasecmp(tok,"notransonly")) {
      alg->trans_only_flag= 0;
    }
    else if (!strcasecmp(tok,"xonly")) {
      alg->x_only_flag= 1;
    }
    else if (!strcasecmp(tok,"noxonly")) {
      alg->x_only_flag= 0;
    }
    else if (!strncasecmp(tok,"opt=",4)) {
      optString= strdup(tok+4);
    }
    else if (!strncasecmp(tok,"weight=",7)) {
      if ((alg->weight_method=parseWeightMethod(tok+7))==WEIGHT_INVALID)
	{ free(work); return 0; }
    }
    else if (!strncasecmp(tok,"wtfloor=",8)) {
      alg->weight_floor=atof(tok+8);
    }
    else if (!strncasecmp(tok,"inner=",6)) {
      if ((alg->inner_search_method=parseSearchMethod(tok+6))
	  ==SEARCH_INVALID)
	{ free(work); return 0; }
    }
    else if (!strncasecmp(tok,"outer=",6)) {
      if ((alg->outer_search_method=parseSearchMethod(tok+6))
	  ==SEARCH_INVALID)
	{ free(work); return 0; }
    }
    else if (!strncasecmp(tok,"qual=",5)) {
      if (!parseQualMethod(tok+5)) { free(work); return 0; }
    }
    else if (!strncasecmp(tok,"obj=",4)) {
      if ((alg->objective_method=parseObjectiveMethod(tok+4))
	  ==OBJECTIVE_INVALID)
	{ free(work); return 0; }
    }
    else if (!strncasecmp(tok,"optol=",6)) {
      tol= atof(tok+6);
      tolSet= 1;
    }
    else if (!strncasecmp(tok,"opscale=",8)) {
      scale= atof(tok+8);
      scaleSet= 1;
    }
    else { free(work); return 0; }
    tok= strtok_r(NULL, " ,+&:;",&tptr);
  }

  /* We need to do this last, because the other elements of the
   * algorithm effect the number of degrees of freedom, which
   * in turn effects the optimizer.
   */
  if (optString != NULL) {
    if (alg->opt) alg->opt->destroySelf(alg->opt);
    alg->opt= parseOptMethod(optString,alg->getNDimCB(alg));
    if (!alg->opt) {
      Abort("%s: internal error: unrecognized optimizer type <%s>!\n",
	    alg->progname, optString);
    }
    alg->opt->setDebugLevel(alg->opt, alg->debugLevel);
    if (tolSet) alg->opt->setTol(alg->opt, tol);
    if (scaleSet) alg->opt->setScale(alg->opt, scale);
    free(optString);
  }
  
  /* This turns out to be a good place for a necessary reality check */
  if (alg->getNDimCB(alg)>MAX_DOF)
    Abort("%s: algorithm implies %d dimensions, but compiled limit is %d!\n",
	  alg->progname,alg->getNDimCB(alg),MAX_DOF);

  /* The algorithm may or may not need a mutualInfo context. */
  for (i=0; i<N_MUTUAL_INFO_CONTEXTS; i++) 
    if (alg->mutualInfoContext[i]) 
      ent_destroyMIContext(alg->mutualInfoContext[i]);
  switch (alg->objective_method) {
  case OBJECTIVE_MSE: break;
  case OBJECTIVE_MI: 
  case OBJECTIVE_JE: 
    for (i=0; i<N_MUTUAL_INFO_CONTEXTS; i++)
      alg->mutualInfoContext[i]= ent_createMIContext();
    break;
  default: Abort("%s: internal error: unknown objective method!\n",
		 alg->progname);
  }
  for (i=0; i<N_MUTUAL_INFO_CONTEXTS; i++)
    if (alg->mutualInfoContext[i] != NULL) {
      ent_setMINBins(alg->mutualInfoContext[i],mutualInfoNBins[i]);
      if (alg->debugLevel>1) ent_setMIDebug(alg->mutualInfoContext[i],1);
      if (alg->debugLevel>0) ent_setMIVerbose(alg->mutualInfoContext[i],1);
    }

  free(work);
  return 1;
}

int algInit( Algorithm* alg, const char* string, const char* progname_in,
	     int (*getNDimCallback_in)())
{
  int i;

  alg->progname= strdup(progname_in);
  alg->getNDimCB= getNDimCallback_in;
  alg->opt= NULL;
  alg->smoother= sm_create_smoother();
  for (i=0; i<N_MUTUAL_INFO_CONTEXTS; i++)
    alg->mutualInfoContext[i]= NULL;
  alg->debugLevel= 0;
  alg->inner_search_method= SEARCH_INVALID;
  alg->outer_search_method= SEARCH_INVALID;
  alg->weight_method= WEIGHT_INVALID;
  alg->smooth_flag= 0;
  alg->inplane_flag= 0;
  alg->rot_only_flag= 0;
  alg->trans_only_flag= 0;
  alg->x_only_flag= 0;

  algParseInfoString(alg, string);

  return(1);
}

const char* algGetInfoString(const Algorithm* alg)
{
  char scratch[128];
  static char result[512];
  int offset= 0;

  result[0]= '\0';

  sprintf(scratch,"qual=%s,",
	  qualMethodName(fshrot3d_get(FR3D_QUAL_MEASURE)));
  strncat(result,scratch,sizeof(result)-(offset+1));
  offset += strlen(scratch);
  if (offset>=sizeof(result)-1) return result;

  if (alg->smooth_flag) strcpy(scratch,"smooth,");
  else strcpy(scratch,"nosmooth,");
  strncat(result,scratch,sizeof(result)-(offset+1));
  offset += strlen(scratch);
  if (offset>=sizeof(result)-1) return result;

  if (alg->inplane_flag) strcpy(scratch,"inplane,");
  else strcpy(scratch,"noinplane,");
  strncat(result,scratch,sizeof(result)-(offset+1));
  offset += strlen(scratch);
  if (offset>=sizeof(result)-1) return result;

  if (alg->rot_only_flag) strcpy(scratch,"rotonly,");
  else strcpy(scratch,"norotonly,");
  strncat(result,scratch,sizeof(result)-(offset+1));
  offset += strlen(scratch);
  if (offset>=sizeof(result)-1) return result;

  if (alg->trans_only_flag) strcpy(scratch,"transonly,");
  else strcpy(scratch,"notransonly,");
  strncat(result,scratch,sizeof(result)-(offset+1));
  offset += strlen(scratch);
  if (offset>=sizeof(result)-1) return result;

  if (alg->x_only_flag) strcpy(scratch,"xonly,");
  else strcpy(scratch,"noxonly,");
  strncat(result,scratch,sizeof(result)-(offset+1));
  offset += strlen(scratch);
  if (offset>=sizeof(result)-1) return result;

  sprintf(scratch,"weight=%s,",weightMethodName(alg->weight_method));
  strncat(result,scratch,sizeof(result)-(offset+1));
  offset += strlen(scratch);
  if (offset>=sizeof(result)-1) return result;  

  sprintf(scratch,"wtfloor=%g,",alg->weight_floor);
  strncat(result,scratch,sizeof(result)-(offset+1));
  offset += strlen(scratch);
  if (offset>=sizeof(result)-1) return result;

  sprintf(scratch,"opt=%s,",alg->opt->getMethodName(alg->opt));
  strncat(result,scratch,sizeof(result)-(offset+1));
  offset += strlen(scratch);
  if (offset>=sizeof(result)-1) return result;  

  sprintf(scratch,"inner=%s,",searchMethodName(alg->inner_search_method));
  strncat(result,scratch,sizeof(result)-(offset+1));
  offset += strlen(scratch);
  if (offset>=sizeof(result)-1) return result;

  sprintf(scratch,"outer=%s,",searchMethodName(alg->outer_search_method));
  strncat(result,scratch,sizeof(result)-(offset+1));
  offset += strlen(scratch);
  if (offset>=sizeof(result)-1) return result;

  sprintf(scratch,"obj=%s,",objectiveMethodName(alg->objective_method));
  strncat(result,scratch,sizeof(result)-(offset+1));
  offset += strlen(scratch);
  if (offset>=sizeof(result)-1) return result;

  sprintf(scratch,"optol=%f,",alg->opt->getTol(alg->opt));
  strncat(result,scratch,sizeof(result)-(offset+1));
  offset += strlen(scratch);
  if (offset>=sizeof(result)-1) return result;

  sprintf(scratch,"opscale=%f",alg->opt->getScale(alg->opt));
  strncat(result,scratch,sizeof(result)-(offset+1));
  offset += strlen(scratch);
  if (offset>=sizeof(result)-1) return result;

  return result;
}

void algFinalize(Algorithm* alg)
{
  int i;
  if (alg->opt != NULL) {
    alg->opt->destroySelf(alg->opt);
    alg->opt= NULL;
  }
  for (i=0; i<N_MUTUAL_INFO_CONTEXTS; i++)
    if (alg->mutualInfoContext[i]!=NULL) { 
      ent_destroyMIContext(alg->mutualInfoContext[i]);
      alg->mutualInfoContext[i]= NULL;
    }
  if (alg->smoother != NULL) {
    sm_destroy( alg->smoother );
    alg->smoother= NULL;
  }
}

void algPack( const Algorithm* alg )
{
  const char* cs;
  char* s;

  cs= algGetInfoString( alg );
  par_pkint(strlen(cs));
  par_pkstr((char*)cs); /* don't free s; algGetInfoString uses a static buffer */

  s= alg->opt->getStringRep(alg->opt);
  par_pkint(strlen(s));
  par_pkstr(s);
  free(s);

  par_pkint((int)(alg->smoother->type));
  par_pkdouble(alg->smoother->bandwidth);
  par_pkdouble(alg->smoother->k);
  par_pkdouble(alg->smoother->threshold);
}

void algUnpack( Algorithm* alg, const char* progname_in,
		int (*getNDimCallback)(Algorithm* alg) )
{
  int l;
  char tbuf[256];
  sm_type smOldType, smType;
  float smOldBand, smBand;
  float smOldK, smK;
  float smOldThresh, smThresh;
  SmootherThreshTest smTest;

  /* The following block ends up building the optimizer twice, once
   * when it unpacks the algorithm and then again when it unpacks
   * the copy of the master's actual optimizer.  We just throw the
   * first one away.
   */
  l= par_upkint();
  if (l > (sizeof(tbuf)-1))
    Abort("%s: cannot unpack algorithm; string rep is too long! (%d vs %d)\n",
	  alg->progname, l, sizeof(tbuf)-1);
  par_upkstr(tbuf);
  if (!algInit(alg, tbuf, progname_in, getNDimCallback)) 
    Abort("%s: internal error: master sent invalid alg string <%s>!\n",
	  alg->progname,tbuf);
  l= par_upkint();
  if (l > (sizeof(tbuf)-1))
    Abort("%s: cannot unpack optimizer; string rep is too long! (%d vs %d)\n",
	  alg->progname, l, sizeof(tbuf)-1);
  par_upkstr(tbuf);
  if (alg->opt != NULL) alg->opt->destroySelf(alg->opt);
  alg->opt= optimizerFromStringRep(tbuf);
  if (!alg->opt)
      Abort("%s: internal error unpacking unrecognized optimizer <%s>!\n",
	    alg->progname, tbuf);
  
  sm_get_params( &smOldType, &smOldBand, &smOldK, &smOldThresh, &smTest );
  smType= (sm_type)par_upkint();
  smBand= par_upkdouble();
  smK= par_upkdouble();
  smThresh= par_upkdouble();
  sm_set_params( smType, smBand, smK, smThresh, smTest );
  if (alg->smoother != NULL) sm_destroy( alg->smoother );
  alg->smoother= sm_create_smoother();
  sm_set_params( smOldType, smOldBand, smOldK, smOldThresh, smTest );  
}

void algPrep()
{
  /* Initialize the smoother package and set default behaviors, in case
   * smoothing is turned on.
   */
  sm_init();
  sm_set_params( SM_GAUSSIAN, 1.0, 0.0, 0.0, NULL );
}

void algSmoothImage( Algorithm* alg, float* img, int dx, int dy, int dz,
		     unsigned char** missing, int t )
{
  /* Function to slightly blur images */
  float *scratch_in;
  float *scratch_out;
  long max_dim;
  long i,j,k;
  
  max_dim= dx;
  if (dy>max_dim) max_dim= dy;
  if (dz>max_dim) max_dim= dz;
  
  if (!(scratch_in= (float*)malloc(max_dim*sizeof(float))))
    Abort("%s: unable to allocate %d bytes!\n",max_dim*sizeof(float));
  if (!(scratch_out= (float*)malloc(max_dim*sizeof(float))))
    Abort("%s: unable to allocate %d bytes!\n",max_dim*sizeof(float));
  
  /* Smooth in X */
  sm_set_direction(alg->smoother,'x');
  for (j=0; j<dy; j++) {
    for (k=0; k<dz; k++) {
      for (i=0; i<dx; i++) scratch_in[i]= MEM(img,dx,dy,dz,i,j,k);
      SM_SMOOTH(alg->smoother, scratch_in, scratch_out, dx, missing, t);
      for (i=0; i<dx; i++) MEM(img,dx,dy,dz,i,j,k)= scratch_out[i];
    }
  }
  
  /* Smooth in Y */
  sm_set_direction(alg->smoother,'y');
  for (k=0; k<dz; k++) {
    for (i=0; i<dx; i++) {
      for (j=0; j<dy; j++) scratch_in[j]= MEM(img,dx,dy,dz,i,j,k);
      SM_SMOOTH(alg->smoother, scratch_in, scratch_out, dy, missing, t);
      for (j=0; j<dy; j++) MEM(img,dx,dy,dz,i,j,k)= scratch_out[j];
    }
  }
  
  /* Smooth in Z */
  sm_set_direction(alg->smoother,'z');
  for (i=0; i<dx; i++) {
    for (j=0; j<dy; j++) {
      for (k=0; k<dz; k++) scratch_in[k]= MEM(img,dx,dy,dz,i,j,k);
      SM_SMOOTH(alg->smoother, scratch_in, scratch_out, dz, missing, t);
      for (k=0; k<dz; k++) MEM(img,dx,dy,dz,i,j,k)= scratch_out[k];
    }
  }
  
  free(scratch_in);
  free(scratch_out);
}

void algMaybeSmoothImage( Algorithm* alg, float* img, int dx, int dy, int dz,
			  unsigned char** missing, int t )
{
  if (alg->smooth_flag) algSmoothImage(alg,img,dx,dy,dz,missing,t);
}

void algSmoothImageComplex( Algorithm* alg, FComplex* img,
			    int dx, int dy, int dz, 
			    unsigned char** missing, int t )
{
  float *scratch_in;
  float *scratch_out;
  long max_dim;
  long i,j,k;
  
  max_dim= dx;
  if (dy>max_dim) max_dim= dy;
  if (dz>max_dim) max_dim= dz;
  
  /* Note that we're only interested in smoothing the real part of this! */
  
  if (!(scratch_in= (float*)malloc(max_dim*sizeof(float))))
    Abort("%s: unable to allocate %d bytes!\n",max_dim*sizeof(float));
  if (!(scratch_out= (float*)malloc(max_dim*sizeof(float))))
    Abort("%s: unable to allocate %d bytes!\n",max_dim*sizeof(float));
  
  /* Smooth in X */
  sm_set_direction(alg->smoother,'x');
  for (j=0; j<dy; j++) {
    for (k=0; k<dz; k++) {
      for (i=0; i<dx; i++) scratch_in[i]= MEM(img,dx,dy,dz,i,j,k).real;
      SM_SMOOTH(alg->smoother, scratch_in, scratch_out, dx, NULL, t);
      for (i=0; i<dx; i++) MEM(img,dx,dy,dz,i,j,k).real= scratch_out[i];
    }
  }
  
  /* Smooth in Y */
  sm_set_direction(alg->smoother,'y');
  for (k=0; k<dz; k++) {
    for (i=0; i<dx; i++) {
      for (j=0; j<dy; j++) scratch_in[j]= MEM(img,dx,dy,dz,i,j,k).real;
      SM_SMOOTH(alg->smoother, scratch_in, scratch_out, dy, NULL, t);
      for (j=0; j<dy; j++) MEM(img,dx,dy,dz,i,j,k).real= scratch_out[j];
    }
  }
  
  /* Smooth in Z */
  sm_set_direction(alg->smoother,'z');
  for (i=0; i<dx; i++) {
    for (j=0; j<dy; j++) {
      for (k=0; k<dz; k++) scratch_in[k]= MEM(img,dx,dy,dz,i,j,k).real;
      SM_SMOOTH(alg->smoother, scratch_in, scratch_out, dz, NULL, t);
      for (k=0; k<dz; k++) MEM(img,dx,dy,dz,i,j,k).real= scratch_out[k];
    }
  }
  
  free(scratch_in);
  free(scratch_out);
}

void algMaybeSmoothImageComplex( Algorithm* alg, FComplex* img, 
				 int dx, int dy, int dz, 
				 unsigned char** missing, int t )
{
  if (alg->smooth_flag) algSmoothImageComplex(alg,img,dx,dy,dz,missing,t);
}

/* This routine calculates the current weight function based on the global
 * task info.  This value is ultimately returned by mse().
 */
double algCalcChiSqr(Algorithm* alg, FComplex* moved_image, 
		     float* align_image, float* weight_image, int* mask,
		     char* check, int dx, int dy, int dz)
{
  double sse;
  long x, y, z;
  long xMin, xMax, yMin, yMax, zMin, zMax;
  register double tmp;
  register double tmp1=0.0, tmp2=0.0;

  /* This sets the boundary over which we take SSE */
  if (alg->x_only_flag) {
    xMin= 0;
    xMax= dx;
    yMin= 0;
    yMax= dy;
    zMin= 0;
    zMax= dz;
  }
  else if (alg->inplane_flag) {
    xMin= 2;
    xMax= dx-2;
    yMin= 2;
    yMax= dy-2;
    zMin= 0;
    zMax= dz;
  }
  else {
    xMin= 2;
    xMax= dx-2;
    yMin= 2;
    yMax= dy-2;
    zMin= 2;
    zMax= dz-2;
  }


  if (xMax<=xMin)
    Abort("%s: x dimension of %d is not enough for alignment!\n",alg->progname,dx);
  if (yMax<=yMin)
    Abort("%s: y dimension of %d is not enough for alignment!\n",alg->progname,dy);
  if (zMax<=zMin)
    Abort("%s: %d slices is not enough for alignment!\n",alg->progname,dz);

#ifdef never
  xMin= dx/4;
  xMax= (3*dx)/4;
  yMin= dx/4;
  yMax= (3*dy)/4;
  if (alg->inplane_flag) {
    zMin= 0;
    zMax= dz;
  }
  else {
    zMin= 2;
    zMax= dz-2;
  }
#endif

  /* Find mean-squared difference */
  sse = 0.0;

  switch (alg->weight_method) {
  case WEIGHT_CONST:
    {
      register long count= 0;
      if (alg->inner_search_method==SEARCH_TRILIN) {
	for( x = xMin; x < xMax; x++ )
	  for( y = yMin; y < yMax; y++ )
	    for( z = zMin; z < zMax; z++ )
	      if( MEM(check,dx,dy,dz,x,y,z) ) {
		tmp1 = MEM(moved_image,dx,dy,dz,x,y,z).real;
		tmp2 = MEM(align_image,dx,dy,dz,x,y,z);
		tmp = tmp1-tmp2;
		sse += tmp*tmp;

#ifdef JENN
		tmp= MEM(moved_image,dx,dy,dz,x,y,z).real
		  - MEM(align_image,dx,dy,dz,x,y,z);
		sse += tmp*tmp;
#endif


		count++;
	      }

#ifdef never
	fprintf(tp, "1a:sse = %.14g, count = %.14g\n", sse, count);
	fprintf(tp, "1a:tmp1 = %.14g, tmp2 = %.14g\n",tmp1, tmp2); 
	fflush(tp);
#endif
      }
      else {
	for( x = xMin; x < xMax; x++ )
	  for( y = yMin; y < yMax; y++ )
	    for( z = zMin; z < zMax; z++ ) {
		tmp1 = MEM(moved_image,dx,dy,dz,x,y,z).real;
		tmp2 = MEM(align_image,dx,dy,dz,x,y,z);
		tmp = tmp1-tmp2;
		sse += tmp*tmp;



#ifdef JENN
	      tmp= MEM(moved_image,dx,dy,dz,x,y,z).real
		- MEM(align_image,dx,dy,dz,x,y,z);
	      sse += tmp*tmp;
#endif

	      count++;
	    }
#ifdef never
	fprintf(stderr, "1b:sse = %.14g, count = %.14g\n", sse, count);
	fprintf(stderr, "1b:tmp1 = %.14g, tmp2 = %.14g\n",tmp1, tmp2); 
	fflush(tp);
#endif
      }
      if (count != 0) sse /= (double)count;

#ifdef never
      fprintf(tp, "2:sse = %.14g \n", sse); 
      fflush(tp);
#endif
    }
    break;
  case WEIGHT_INV_STDV:
  case WEIGHT_ALIGN:
  case WEIGHT_SMTHALIGN:
    {
      /* put print statements in here jenn */

      register double weight= 0.0;
      register double tot_weight= 0.0;
      if (alg->inner_search_method==SEARCH_TRILIN) {
	for( x = xMin; x < xMax; x++ )
	  for( y = yMin; y < yMax; y++ )
	    for( z = zMin; z < zMax; z++ ) {
	      weight= MEM(weight_image,dx,dy,dz,x,y,z);
	      if (weight != 0.0 && MEM(check,dx,dy,dz,x,y,z)) {
		tmp1 = MEM(moved_image,dx,dy,dz,x,y,z).real;
		tmp2 = MEM(align_image,dx,dy,dz,x,y,z);
		tmp = tmp1-tmp2;
		sse += weight*tmp*tmp;
		tot_weight += weight;
	      }
	    }
#ifdef never
	fprintf(stderr, "1a:sse = %.14g, tot_weight = %.14g, weight = %.14g\n", sse, tot_weight, weight);
	fprintf(stderr, "1a:tmp1 = %.14g, tmp2 = %.14g\n",tmp1, tmp2); 
	fflush(tp);
#endif
      }
      else {
	for( x = xMin; x < xMax; x++ )
	  for( y = yMin; y < yMax; y++ )
	    for( z = zMin; z < zMax; z++ ) {
	      weight= MEM(weight_image,dx,dy,dz,x,y,z);
	      if (weight != 0.0) {
		tmp1 = MEM(moved_image,dx,dy,dz,x,y,z).real;
		tmp2 = MEM(align_image,dx,dy,dz,x,y,z);
		tmp = tmp1-tmp2;
#ifdef JENN
		tmp= MEM(moved_image,dx,dy,dz,x,y,z).real
		  - MEM(align_image,dx,dy,dz,x,y,z);
#endif
		sse += weight*tmp*tmp;
		tot_weight += weight;
	      }
	    }

#ifdef never
	fprintf(stderr, "1b:sse = %.14g, tot_weight = %.14g, weight = %.14g\n", sse, tot_weight, weight);
	fprintf(stderr, "1b:tmp1 = %.14g, tmp2 = %.14g\n",tmp1, tmp2); 
	fflush(tp);
#endif
      }

      
      if (tot_weight != 0.0) sse /= tot_weight;

#ifdef never
      fprintf(tp, "2:sse = %.14g \n", sse); 
      fflush(tp);
#endif
  }
  break;
 case WEIGHT_INVALID:
    Abort("%s: internal error: invalid weight method in algCalcChiSqr!\n",
	  alg->progname);
  }

#ifdef never  
      fprintf(tp, "3:sse = %.14g\n", sse); 
      fflush(tp);
#endif

  return sse;
}

/* This routine calculates the current weight function when mutual information
 * is used, based on the global task info.  This value is ultimately returned 
 * by mse().
 */
double algCalcMutualInfo(Algorithm* alg, FComplex* moved_image, 
			 float* align_image, float* weight_image, int* mask,
			 char* check, int dx, int dy, int dz)
{
  long zMin;
  long zMax;
  double sum;
  int iCtx;

  /* This sets the boundary over which we take mutual info */
  if (alg->inplane_flag || alg->x_only_flag) {
    zMin= 0;
    zMax= dz;
  }
  else {
    zMin= 2;
    zMax= dz-2;
  }
  if (zMax<=zMin)
    Abort("%s: %d slices is not enough for alignment!\n",alg->progname,dz);

  /* Calculate the mutual info.  If this is our first pass, set the bounds 
   * for the histograms first. */
  if (algNeedsMask(alg)) {
    if (!ent_getMIMax1Set(alg->mutualInfoContext[0])) {
      long i;
      long j;
      long k;
      int init= 0;
      double min1= 0.0;
      double max1= 0.0;
      double min2= 0.0;
      double max2= 0.0;
      double ave1= 0.0;
      double ave2= 0.0;
      for (k=zMin; k<zMax; k++)
	for (j=0; j<dy; j++)
	  for (i=0; i<dx; i++) 
	    if (MEM(mask,dx,dy,dz,i,j,k)) {
	      if (!init) {
		max1= min1= MEM(align_image,dx,dy,dz,i,j,k);
		max2= min2= MEM(moved_image,dx,dy,dz,i,j,k).real;
		init= 1;
	      }
	      else {
		if (MEM(align_image,dx,dy,dz,i,j,k)<min1) 
		  min1= MEM(align_image,dx,dy,dz,i,j,k);
		if (MEM(align_image,dx,dy,dz,i,j,k)>max1) 
		  max1= MEM(align_image,dx,dy,dz,i,j,k);
		if (MEM(moved_image,dx,dy,dz,i,j,k).real<min2) 
		  min2= MEM(moved_image,dx,dy,dz,i,j,k).real;
		if (MEM(moved_image,dx,dy,dz,i,j,k).real>max2) 
		  max2= MEM(moved_image,dx,dy,dz,i,j,k).real;
	      }
	    }
      /* Widen by 10% at both ends */
      ave1= 0.5*(min1+max1);
      min1= ave1 - 1.1*(ave1-min1);
      max1= ave1 + 1.1*(max1-ave1);
      ave2= 0.5*(min2+max2);
      min2= ave2 - 1.1*(ave2-min2);
      max2= ave2 + 1.1*(max2-ave2);
      for (iCtx=0; iCtx<N_MUTUAL_INFO_CONTEXTS; iCtx++) {
	ent_setMIMin1(alg->mutualInfoContext[iCtx],min1);
	ent_setMIMax1(alg->mutualInfoContext[iCtx],max1);
	ent_setMIMin2(alg->mutualInfoContext[iCtx],min2);
	ent_setMIMax2(alg->mutualInfoContext[iCtx],max2);
      }
    }
    sum= 0.0;
    for (iCtx=0; iCtx<N_MUTUAL_INFO_CONTEXTS; iCtx++)
      sum += ent_calcMaskedMutualInformationFloat(alg->mutualInfoContext[iCtx],
						  align_image + zMin*dx*dy, 
						  (float*)(moved_image
							   +zMin*dx*dy),
						  mask + zMin*dx*dy, 
						  dx, dy, (zMax-zMin), 
						  1, 2, 1);
    return -1.0*sum/(double)N_MUTUAL_INFO_CONTEXTS;
  }
  else {
    if (!ent_getMIMax1Set(alg->mutualInfoContext[0])) {
      long i;
      long j;
      long k;
      double min1= MEM(align_image,dx,dy,dz,0,0,zMin);
      double max1= min1;
      double min2= MEM(moved_image,dx,dy,dz,0,0,zMin).real;
      double max2= min2;
      double ave1;
      double ave2;
      for (k=zMin; k<zMax; k++)
	for (j=0; j<dy; j++)
	  for (i=0; i<dx; i++) {
	    if (MEM(align_image,dx,dy,dz,i,j,k)<min1) 
	      min1= MEM(align_image,dx,dy,dz,i,j,k);
	    if (MEM(align_image,dx,dy,dz,i,j,k)>max1) 
	      max1= MEM(align_image,dx,dy,dz,i,j,k);
	    if (MEM(moved_image,dx,dy,dz,i,j,k).real<min2) 
	      min2= MEM(moved_image,dx,dy,dz,i,j,k).real;
	    if (MEM(moved_image,dx,dy,dz,i,j,k).real>max2) 
	      max2= MEM(moved_image,dx,dy,dz,i,j,k).real;
	  }
      ave1= 0.5*(min1+max1);
      min1= ave1 - 1.1*(ave1-min1);
      max1= ave1 + 1.1*(max1-ave1);
      ave2= 0.5*(min2+max2);
      min2= ave2 - 1.1*(ave2-min2);
      max2= ave2 + 1.1*(max2-ave2);
      for (iCtx= 0; iCtx<N_MUTUAL_INFO_CONTEXTS; iCtx++) {
	ent_setMIMin1(alg->mutualInfoContext[iCtx],min1);
	ent_setMIMax1(alg->mutualInfoContext[iCtx],max1);
	ent_setMIMin2(alg->mutualInfoContext[iCtx],min2);
	ent_setMIMax2(alg->mutualInfoContext[iCtx],max2);
      }
    }
    sum= 0.0;
    for (iCtx=0; iCtx<N_MUTUAL_INFO_CONTEXTS; iCtx++) 
      sum += ent_calcMutualInformationFloat(alg->mutualInfoContext[iCtx],
					    align_image + zMin*dx*dy, 
					    (float*)(moved_image+zMin*dx*dy),
					    dx, dy, (zMax-zMin), 1, 2);
    return -1.0*sum/(double)N_MUTUAL_INFO_CONTEXTS;
  }
}

/* This routine calculates the current weight function when mutual information
 * is used, based on the global task info.  This value is ultimately returned 
 * by mse().
 */
double algCalcJointEntropy(Algorithm* alg, FComplex* moved_image, 
			   float* align_image, float* weight_image, int* mask,
			   char* check, int dx, int dy, int dz)
{
  long zMin;
  long zMax;
  int iCtx;
  double sum;

  /* This sets the boundary over which we take mutual info */
  if (alg->inplane_flag || alg->x_only_flag) {
    zMin= 0;
    zMax= dz;
  }
  else {
    zMin= 2;
    zMax= dz-2;
  }
  if (zMax<=zMin)
    Abort("%s: %d slices is not enough for alignment!\n",alg->progname,dz);

  /* Calculate the mutual info.  If this is our first pass, set the bounds 
   * for the histograms first. */
  if (algNeedsMask(alg)) {
    if (!ent_getMIMax1Set(alg->mutualInfoContext[0])) {
      long i;
      long j;
      long k;
      int init= 0;
      double min1= 0.0;
      double max1= 0.0;
      double min2= 0.0;
      double max2= 0.0;
      double ave1= 0.0;
      double ave2= 0.0;
      for (k=zMin; k<zMax; k++)
	for (j=0; j<dy; j++)
	  for (i=0; i<dx; i++) 
	    if (MEM(mask,dx,dy,dz,i,j,k)) {
	      if (!init) {
		max1= min1= MEM(align_image,dx,dy,dz,i,j,k);
		max2= min2= MEM(moved_image,dx,dy,dz,i,j,k).real;
		init= 1;
	      }
	      else {
		if (MEM(align_image,dx,dy,dz,i,j,k)<min1) 
		  min1= MEM(align_image,dx,dy,dz,i,j,k);
		if (MEM(align_image,dx,dy,dz,i,j,k)>max1) 
		  max1= MEM(align_image,dx,dy,dz,i,j,k);
		if (MEM(moved_image,dx,dy,dz,i,j,k).real<min2) 
		  min2= MEM(moved_image,dx,dy,dz,i,j,k).real;
		if (MEM(moved_image,dx,dy,dz,i,j,k).real>max2) 
		  max2= MEM(moved_image,dx,dy,dz,i,j,k).real;
	      }
	    }
      /* Widen by 10% at both ends */
      ave1= 0.5*(min1+max1);
      min1= ave1 - 1.1*(ave1-min1);
      max1= ave1 + 1.1*(max1-ave1);
      ave2= 0.5*(min2+max2);
      min2= ave2 - 1.1*(ave2-min2);
      max2= ave2 + 1.1*(max2-ave2);
      for (iCtx=0; iCtx<N_MUTUAL_INFO_CONTEXTS; iCtx++) {
	ent_setMIMin1(alg->mutualInfoContext[iCtx],min1);
	ent_setMIMax1(alg->mutualInfoContext[iCtx],max1);
	ent_setMIMin2(alg->mutualInfoContext[iCtx],min2);
	ent_setMIMax2(alg->mutualInfoContext[iCtx],max2);
      }
    }
    sum= 0.0;
    for (iCtx=0; iCtx<N_MUTUAL_INFO_CONTEXTS; iCtx++) 
      sum += ent_calcMaskedJointEntropyFloat(alg->mutualInfoContext[iCtx],
					     align_image + zMin*dx*dy, 
					     (float*)(moved_image+zMin*dx*dy),
					     mask + zMin*dx*dy, 
					     dx, dy, (zMax-zMin), 1, 2, 1);
    return sum/(double)N_MUTUAL_INFO_CONTEXTS;
  }
  else {
    if (!ent_getMIMax1Set(alg->mutualInfoContext[0])) {
      long i;
      long j;
      long k;
      double min1= MEM(align_image,dx,dy,dz,0,0,zMin);
      double max1= min1;
      double min2= MEM(moved_image,dx,dy,dz,0,0,zMin).real;
      double max2= min2;
      double ave1;
      double ave2;
      for (k=zMin; k<zMax; k++)
	for (j=0; j<dy; j++)
	  for (i=0; i<dx; i++) {
	    if (MEM(align_image,dx,dy,dz,i,j,k)<min1) 
	      min1= MEM(align_image,dx,dy,dz,i,j,k);
	    if (MEM(align_image,dx,dy,dz,i,j,k)>max1) 
	      max1= MEM(align_image,dx,dy,dz,i,j,k);
	    if (MEM(moved_image,dx,dy,dz,i,j,k).real<min2) 
	      min2= MEM(moved_image,dx,dy,dz,i,j,k).real;
	    if (MEM(moved_image,dx,dy,dz,i,j,k).real>max2) 
	      max2= MEM(moved_image,dx,dy,dz,i,j,k).real;
	  }
      ave1= 0.5*(min1+max1);
      min1= ave1 - 1.1*(ave1-min1);
      max1= ave1 + 1.1*(max1-ave1);
      ave2= 0.5*(min2+max2);
      min2= ave2 - 1.1*(ave2-min2);
      max2= ave2 + 1.1*(max2-ave2);
      for (iCtx=0; iCtx<N_MUTUAL_INFO_CONTEXTS; iCtx++) {
	ent_setMIMin1(alg->mutualInfoContext[iCtx],min1);
	ent_setMIMax1(alg->mutualInfoContext[iCtx],max1);
	ent_setMIMin2(alg->mutualInfoContext[iCtx],min2);
	ent_setMIMax2(alg->mutualInfoContext[iCtx],max2);
      }
    }
    sum= 0.0;
    for (iCtx=0; iCtx<N_MUTUAL_INFO_CONTEXTS; iCtx++) 
      sum += ent_calcJointEntropyFloat(alg->mutualInfoContext[iCtx],
				     align_image + zMin*dx*dy, 
				     (float*)(moved_image+zMin*dx*dy),
				     dx, dy, (zMax-zMin), 1, 2);
    return sum/(double)N_MUTUAL_INFO_CONTEXTS;
  }
}

void algMaybeBuildMask(Algorithm* alg, float* weight_image, int* mask, 
		       int dx, int dy, int dz)
{
  long x,y,z;

  if (algNeedsMask(alg)) {
    for (x=0; x<dx; x++)
      for (y=0; y<dy; y++)
	for (z=0; z<dz; z++) {
	  /* When master calculated the weight, all values below
	   * the appropriate threshold were set to zero.
	   */
	  MEM(mask,dx,dy,dz,x,y,z)= 
	    (MEM(weight_image,dx,dy,dz,x,y,z) != 0.0);
	}
  }
}

void algBuildWeight(Algorithm* alg, float* weight_image, float* align_image, 
		    MRI_Dataset* Align, MRI_Dataset* Stdv, 
		    int dx, int dy, int dz, 
		    unsigned char** alignMissing, unsigned char** stdvMissing)
{
  switch (alg->weight_method) {
  case WEIGHT_CONST:
    {
      long i;
      for (i=0; i<dx*dy*dz; i++) weight_image[i]= 1.0;
    }
    break;
  case WEIGHT_INV_STDV:
    {
      long i;
      loadImage( weight_image, Stdv, 0, dx, dy, dz );
      algMaybeSmoothImage(alg,weight_image,dx,dy,dz,stdvMissing,0);
      invertStdvImage(weight_image, dx, dy, dz);
      for (i=0; i<dx*dy*dz; i++) 
	weight_image[i]= weight_image[i]*weight_image[i];
    }
    break;
  case WEIGHT_ALIGN:
    {
      long i;
      double maxWeight= 0.0;

      /* By definition, we want smoothed version if smooth option set */
      copyImage( weight_image, align_image, dx, dy, dz );
      for (i=0; i<dx*dy*dz; i++) 
	if (fabs(weight_image[i])>maxWeight) maxWeight= fabs(weight_image[i]);
      for (i=0; i<dx*dy*dz; i++)
	if (fabs(weight_image[i])<alg->weight_floor*maxWeight)
	  weight_image[i]= 0.0;
    }
    break;
  case WEIGHT_SMTHALIGN:
    {
      long i;
      double maxWeight= 0.0;

      /* We want to avoid possibly smoothing twice, so we reload */
      loadImage( weight_image, Align, 0, dx, dy, dz );
      algSmoothImage( alg, weight_image, dx, dy, dz, alignMissing,0 );
      for (i=0; i<dx*dy*dz; i++) 
	if (fabs(weight_image[i])>maxWeight) maxWeight= fabs(weight_image[i]);
      for (i=0; i<dx*dy*dz; i++)
	if (fabs(weight_image[i])<alg->weight_floor*maxWeight)
	  weight_image[i]= 0.0;
    }
    break;
  }
}
