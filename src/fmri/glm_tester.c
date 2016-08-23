/************************************************************
 *                                                          *
 *  glm_tester.c                                            *
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
 *  Original programming by Joel Welling, 6/98              *
 ************************************************************/
/* This package tests glm */

#include <stdio.h>
#include <stdlib.h>
#include "glm.h"

#ifdef never
#define REGRESSORTYPE GLM_TYPE_LLSQ
#define COMPLEX_FLAG 0
#define NFACTORS 4
#define NOBS 10

double* counts= NULL;
double obs[NOBS]= 
/*{ 1.0, 2.0, 9.0, 28.0, 65.0, 126.0, 217.0, 344.0, 513.0, 730.0 };*/
/*{ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };*/
/*{ 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };*/
/*{ 1.0, 2.0, 5.0, 10.0, 17.0, 26.0, 37.0, 50.0, 65.0, 82.0 };*/
{0.508828445016,
   0.761312108596,
   0.940834004308,
   0.624964552941,
   0.337026517753,
   0.771881697977,
   0.603018572652,
   0.275402021434,
   0.563790022863,
   0.528420860147};

/*  { 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0 }; */

/*
double factors[NFACTORS*NOBS]= {
  1.0, 0.0, 0.0, 0.0,
  1.0, 1.0, 1.0, 1.0,
  1.0, 2.0, 4.0, 8.0, 
  1.0, 3.0, 9.0, 27.0,
  1.0, 4.0, 16.0, 64.0,
  1.0, 5.0, 25.0, 125.0,
  1.0, 6.0, 36.0, 216.0,
  1.0, 7.0, 49.0, 343.0,
  1.0, 8.0, 64.0, 512.0,
  1.0, 9.0, 81.0, 729.0
};
*/

/* Orthogonal unit vectors going as the 0th through 3rd powers */
double factors[NFACTORS*NOBS]= {
  0.31622776601683794, -0.49543369430686218, 0.52223296786709272, -0.45342519294190498, 
  0.31622776601683794, -0.38533731779422614, 0.1740776559556973, 0.1511417309806373, 
  0.31622776601683794, -0.27524094128159005, -0.087038827977849564, 0.37785432745159075, 
  0.31622776601683794, -0.16514456476895406, -0.26111648393354708, 0.33467097574283711, 
  0.31622776601683794, -0.055048188256317979, -0.34815531191139593, 0.1295500551262605, 
  0.31622776601683794, 0.055048188256318034, -0.3481553119113957, -0.12955005512625861, 
  0.31622776601683794, 0.16514456476895409, -0.26111648393354692, -0.33467097574283478, 
  0.31622776601683794, 0.27524094128159016, -0.087038827977848954, -0.37785432745158876, 
  0.31622776601683794, 0.38533731779422625, 0.17407765595569771, -0.1511417309806343, 
  0.31622776601683794, 0.49543369430686218, 0.52223296786709394, 0.45342519294190781
};


/*
double factors[NFACTORS*NOBS]= {
  1.0, -4.5,
  1.0, -3.5,
  1.0, -2.5,
  1.0, -1.5,
  1.0, -0.5,
  1.0, 0.5,
  1.0, 1.5,
  1.0, 2.5,
  1.0, 3.5,
  1.0, 4.5,
};
*/

#endif

#ifdef never
#define REGRESSORTYPE GLM_TYPE_LLSQ
#define COMPLEX_FLAG 1
double* counts= NULL;
double obs[2*NOBS]= { 
1.0, 1.0, 
2.0, 2.0,
5.0, 5.0,
10.0, 10.0,
17.0, 17.0,
26.0, 26.0,
37.0, 37.0,
50.0, 50.0,
65.0, 65.0,
82.0, 82.0
};
#endif

/* G. Rodrigues' birth control regression dataset */
#define REGRESSORTYPE GLM_TYPE_LOGISTIC
#define COMPLEX_FLAG 0
#define NFACTORS 4
#define NOBS 16

double counts[]= {59, 14, 264, 60,
		  74, 29, 209, 92,
		  145, 157, 164, 146,
		  41, 94, 16, 43}; /* Total counts in each category */
double obs[]=    {6, 4, 52, 10,
		  14, 10, 54, 27,
		  33, 80, 46, 78,
		  6, 48, 8, 31}; /* Yes answers in each category */
/* Factors: first col is age group, then education, then desires kids */
double factors[]= {1,0,0,0,
		   1,0,0,1,
		   1,0,1,0,
		   1,0,1,1,
		   1,1,0,0,
		   1,1,0,1,
		   1,1,1,0,
		   1,1,1,1,
		   1,2,0,0,
		   1,2,0,1,
		   1,2,1,0,
		   1,2,1,1,
		   1,3,0,0,
		   1,3,0,1,
		   1,3,1,0,
		   1,3,1,1};

#ifdef never

/* Blindness-by-age example from 'Intro to R' */
#define REGRESSORTYPE GLM_TYPE_LOGISTIC
#define COMPLEX_FLAG 0
#define NFACTORS 2
#define NOBS 5

double counts[]= {50, 50, 50, 50, 50}; /* Total subjects in each group */
double obs[]=    {6,17,26,37,44}; /* Subjects who are blind */
double factors[]= { 1, 20, 
		    1, 35, 
		    1, 45, 
		    1, 55, 
		    1, 70 }; /* Constant, and Ages of each group */
#endif

#ifdef never
/* From R 'example(glm)' */
#define REGRESSORTYPE GLM_TYPE_POISSON
#define COMPLEX_FLAG 0
#define NFACTORS 5
#define NOBS 9

double obs[]= {18, 17, 15, 20, 10, 20, 25, 13, 12};
double *counts= NULL;
/* factors are treatment and outcome, treated as catagorical */
/* So the X matrix which follows is: constant, treatment2, treatment3,
 * outcome, outcome3.
 */
double factors[]= { 1, 0, 0, 0, 0,
		    1, 0, 0, 1, 0,
		    1, 0, 0, 0, 1,
		    1, 1, 0, 0, 0,
		    1, 1, 0, 1, 0,
		    1, 1, 0, 0, 1,
		    1, 0, 1, 0, 0,
		    1, 0, 1, 1, 0,
		    1, 0, 1, 0, 1};
#endif

#ifdef never

#define REGRESSORTYPE GLM_TYPE_POISSON
#define COMPLEX_FLAG 0
#define NFACTORS 3
#define NOBS 9

double obs[]= {18, 17, 15, 20, 10, 20, 25, 13, 12};
double *counts= NULL;
/* factors are treatment, treatment=treatment2, and treatment=treatment3*/
double factors[]= { 1, 0, 0,
		    1, 0, 0,
		    1, 0, 0,
		    1, 1, 0,
		    1, 1, 0,
		    1, 1, 0,
		    1, 0, 1,
		    1, 0, 1,
		    1, 0, 1};
#endif

static void dump_data(Regressor* r, int afterFlag,
		      double* obs, double* factors, double* params,
		      int nparams, int complex_flag)
{
  int i;
  int j;
  int iparam;
  double* param_off;

  if (afterFlag) {
    if (glm_get(r,GLM_RESIDUALS))
      fprintf(stderr,"residuals:\n");
    else 
      fprintf(stderr,"observations:\n");
    if (complex_flag) {
      for (i=0; i<NOBS; i++) fprintf(stderr,"(%f,%f) ",obs[2*i],obs[(2*i)+1]);
    }
    else {
      for (i=0; i<NOBS; i++) fprintf(stderr,"%f ",obs[i]);
    }
    fprintf(stderr,"\n");
  }
  else {
    fprintf(stderr,"observations:\n");
    if (complex_flag) {
      for (i=0; i<NOBS; i++) fprintf(stderr,"(%f,%f) ",obs[2*i],obs[(2*i)+1]);
    }
    else {
      for (i=0; i<NOBS; i++) fprintf(stderr,"%f ",obs[i]);
    }
    fprintf(stderr,"\n");
  }
    
  fprintf(stderr,"factors:\n");
  for (i=0; i<NOBS; i++) {
    for (j=0; j<NFACTORS; j++) fprintf(stderr,"%f ",factors[(i*NFACTORS)+j]);
    fprintf(stderr,"\n");
  }
    
  fprintf(stderr,"b estimates:\n");
  if (complex_flag) {
    for (i=0; i<NFACTORS; i++) fprintf(stderr,"(%f,%f) ",
				       params[2*i],params[(2*i)+1]);
    param_off= params+(2*NFACTORS);
  }
  else {
    for (i=0; i<NFACTORS; i++) fprintf(stderr,"%g ",params[i]);
    param_off= params+NFACTORS;
  }
  fprintf(stderr,"\n");

  if (glm_get(r,GLM_VARIANCES)) {
    fprintf(stderr,"variances:\n");
    for (i=0; i<NFACTORS; i++) fprintf(stderr,"%f ",param_off[i]);
    fprintf(stderr,"\n");
    param_off += NFACTORS;
  }

  if (glm_get(r,GLM_COVARIANCES)) {
    fprintf(stderr,"covariances:\n");
    if (complex_flag) {
      for (i=0; i<2*NFACTORS; i++) {
	for (j=0; j<2*NFACTORS; j++) fprintf(stderr,"%f ",
					     param_off[(i*NFACTORS)+j]);
	fprintf(stderr,"\n");
      }
      param_off += 2*2*NFACTORS*NFACTORS;
    }
    else {
      for (i=0; i<NFACTORS; i++) {
	for (j=0; j<NFACTORS; j++) fprintf(stderr,"%f ",
					   param_off[(i*NFACTORS)+j]);
	fprintf(stderr,"\n");
      }
      param_off += NFACTORS*NFACTORS;
    }
  }

  if (glm_get(r,GLM_SSQR)) {
    fprintf(stderr,"sum square terms\n");
    for (i=0; i<NFACTORS+2; i++) fprintf(stderr,"%f ",param_off[i]);
    fprintf(stderr,"\n");
    param_off += NFACTORS+2;
  }

  if (glm_get(r,GLM_ORTHO)) {
    fprintf(stderr,"ortho measure: %g\n",param_off[0]);
    param_off++;
  }

  if (glm_get(r,GLM_DEVIANCE)) {
    fprintf(stderr,"deviance: %g\n",param_off[0]);
    param_off++;
  }

  if (glm_get(r,GLM_DEBUG)) fprintf(stderr,"Debug flag is set\n");
  else fprintf(stderr,"Debug flag is NOT set\n");
}

int main()
{
  int nparams;
  int complex_flag= COMPLEX_FLAG;
  int retcode;
  int complex_fac;
  double* params;
  int i;
  int j;
  double tmp;
  Regressor* r= NULL;
  double* ssqr_terms= NULL;

  switch (REGRESSORTYPE) {
  case GLM_TYPE_LLSQ: 
    r=glm_create_llsq_regressor();
    break;
  case GLM_TYPE_LOGISTIC: 
    r=glm_create_logistic_regressor();
    break;
  case GLM_TYPE_POISSON: 
    r=glm_create_poisson_regressor();
    break;
  }

  if (!r) {
    fprintf(stderr,"Unable to create regressor!\n");
    exit(-1);
  }

  glm_set(r,GLM_COMPLEX,complex_flag);
  if (glm_error_msg()) {
    fprintf(stderr,"glm error: <%s>\n",glm_error_msg());
    glm_clear_error_msg();
  }
			       
  glm_set(r,GLM_VARIANCES,1);
  if (glm_error_msg()) {
    fprintf(stderr,"glm error: <%s>\n",glm_error_msg());
    glm_clear_error_msg();
  }

  glm_set(r,GLM_COVARIANCES,1);
  if (glm_error_msg()) {
    fprintf(stderr,"glm error: <%s>\n",glm_error_msg());
    glm_clear_error_msg();
  }

  glm_set(r,GLM_DEVIANCE,1);
  if (glm_error_msg()) {
    fprintf(stderr,"glm error: <%s>\n",glm_error_msg());
    glm_clear_error_msg();
  }

  glm_set(r,GLM_RESIDUALS,1);
  if (glm_error_msg()) {
    fprintf(stderr,"glm error: <%s>\n",glm_error_msg());
    glm_clear_error_msg();
  }

  glm_set(r,GLM_SSQR,1);
  if (glm_error_msg()) {
    fprintf(stderr,"glm error: <%s>\n",glm_error_msg());
    glm_clear_error_msg();
  }

  glm_set(r,GLM_ORTHO,1);
  if (glm_error_msg()) {
    fprintf(stderr,"glm error: <%s>\n",glm_error_msg());
    glm_clear_error_msg();
  }

  glm_set(r,GLM_DEBUG,1);
  if (glm_error_msg()) {
    fprintf(stderr,"glm error: <%s>\n",glm_error_msg());
    glm_clear_error_msg();
  }

  nparams= glm_n_params(r,NFACTORS);

  complex_fac= (complex_flag ? 2 : 1);
  params= (double*)malloc(nparams*complex_fac*sizeof(double));

  /* Some code to generate observations and factors on the fly */
#ifdef never
  obs= (double*)malloc(NOBS*complex_fac*sizeof(double));
  factors= (double*)malloc(NOBS*NFACTORS*complex_fac*sizeof(double));

  if (complex_flag) {
    exit(-1);
#ifdef never
    for (i=0; i<NOBS; i++) {
      obs[2*i]= (double)i;
      obs[(2*i)+1]= ((double)i + 0.5);
      factors[2*(i*NFACTORS)]= 1.0; /* for const component of fit */
      factors[(2*(i*NFACTORS))+1]= 1.0; /* for const component of fit */
      for (j=1; j<NFACTORS; j++) {
	factors[2*((i*NFACTORS) + j)]= 100.0*i+j;
	factors[(2*((i*NFACTORS) + j))+1]= (100.0*i+j) + 0.5;
      }
    }
#endif
  }
  else {
    for (i=0; i<NOBS; i++) {
      obs[i]= 1.23*factors[i*NFACTORS]
	+ 4.56*factors[i*NFACTORS+1]
	+ 7.89*factors[i*NFACTORS+2]
	+ 1.35*factors[i*NFACTORS+3];
#ifdef never
      factors[(i*NFACTORS)]= 1.0; /* for const component of fit */
      for (j=1; j<NFACTORS; j++) {
	factors[(i*NFACTORS) + j]= 100.0*i+j;
      }
#endif
    }
  }
#endif

  fprintf(stderr,"*****Before:******\n");
  dump_data(r, 0, obs, factors, params, nparams, complex_flag);
  fprintf(stderr,"*****Fitting:*****\n");

  retcode= glm_fit(r,obs, factors, counts, params, NOBS, NFACTORS);
  if (retcode) fprintf(stderr,"glm_fit error: <%s>\n",glm_error_msg());

  fprintf(stderr,"*****After:******\n");
  fprintf(stderr,"retcode %d\n",retcode);
  dump_data(r, 1, obs, factors, params, nparams, complex_flag);

  if (glm_get(r,GLM_SSQR)) {
    if (complex_flag) {
      ssqr_terms= params+2*NFACTORS;
      if (glm_get(r,GLM_VARIANCES)) ssqr_terms += 2*NFACTORS;
      if (glm_get(r,GLM_COVARIANCES)) ssqr_terms += 2*2*NFACTORS*NFACTORS;
    }
    else {
      ssqr_terms= params+NFACTORS;
      if (glm_get(r,GLM_VARIANCES)) ssqr_terms += NFACTORS;
      if (glm_get(r,GLM_COVARIANCES)) ssqr_terms += NFACTORS*NFACTORS;
    }
  }

  if (glm_get(r,GLM_TYPE)==GLM_TYPE_LLSQ && glm_get(r,GLM_RESIDUALS)) {
    fprintf(stderr,"Sum Squares Check:\n");
    tmp= 0.0;
    if (complex_flag) {
      for (i=0; i<NOBS; i++) tmp += (obs[2*i]*obs[2*i]) 
			       + (obs[(2*i)+1]*obs[(2*i)+1]);
    }
    else {
      for (i=0; i<NOBS; i++) tmp += (obs[i]*obs[i]);
    }
    fprintf(stderr,"  SSE: %f\n",tmp);
    fprintf(stderr,"  SSTO: %f\n",ssqr_terms[0]);
    tmp= 0.0;
    for (i=0; i<NFACTORS; i++) tmp += ssqr_terms[2+i];
    fprintf(stderr,"  SSR: %f\n",tmp-ssqr_terms[1]);
  }
  glm_destroy(r);
  return 0;
}
