/************************************************************
 *                                                          *
 *  glm_irls.c                                          *
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
/* This package fits a general linear model to data */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <string.h>
#include "lapack.h"
#include <glm.h>

#ifdef DARWIN
#define finite( foo ) isfinite( foo )
#endif

/* Coding notes-
   -In fact, I suspect the notion of complex data is nonsense
    here and I should just rip it out.  At least for logistic
    regression.
   -Each regressor should really have its own error message buffer.
 */

/* Usage Notes-
 * -be sure to use glm_set to request particular outputs *before* using
 *  glm_n_params to find number of parameters returned.
 * -options:
 *    GLM_RESIDUALS:   calculate deviance residuals.
 *    GLM_COVARIANCES: returns the covariance matrix of the estimate beta,
 *                 given by (X'X)^-1 of the inner llsq regressor
 *    GLM_DEVIANCE: returns the deviance statistic for the fit.
 * -parameters returned: 
 *  if complex_flag is set, params contains:
 *   2*nfactors entries for the b value estimates
 *   2*2*nfactors*nfactors entries for the covariances of the b values
 *   1 scalar giving the deviance statistic
 *  if complex_flag is not set, params contains:
 *   1*nfactors entries for the b value estimates
 *   1*nfactors*nfactors entries for covariances
 *   1 scalar giving the deviance statistic
 *
 * Implementation Notes-
 * -subscripts in the range i,j,k,etc. denote values 0<i<nobs, or
 *  values the range of which is intended to be general (like for matrix
 *  math).  Subscripts in the range a,b,c,etc. denote values 0<a<nfactors.
 */

#define MALLOC_FAILURE(num,type)\
  { fprintf(stderr,"%s:%d: %s: unable to allocate %ld %ss\n",\
            __FILE__,__LINE__,__func__,(long)num,#type); exit(-1); }

#define CHECK_DEBUG(r) (r->get(r,GLM_DEBUG))

#define MAX_ITERATIONS 25

typedef void (*CALCINITIALWEIGHTSFUNC)(double* weights, int nobs);
typedef void (*CALCINITIALZFUNC)(double* z, const double* y, const double* n,
				 int nobs);
typedef void (*CALCZFUNC)(double* z, const double* n, const double* y,
			  const double* eta, int nobs);
typedef void (*CALCWFUNC)(double* w, const double* n, const double* eta, 
			  int nobs);
typedef double (*CALCDEVIANCEFUNC)( const double* obs, const double* counts,
				    const double* eta, int nobs );
typedef void (*CALCDEVIANCERESIDUALSFUNC)( double* obs, 
					   const double* counts,
					   const double* eta, int nobs );

typedef struct regressor_irls_data_struct {
  Regressor* innerRegressor;
  CALCINITIALWEIGHTSFUNC calcInitialWeights;
  CALCINITIALZFUNC calcInitialZ;
  CALCZFUNC calcZ;
  CALCWFUNC calcW;
  CALCDEVIANCEFUNC calcDeviance;
  CALCDEVIANCEFUNC calcDevianceComplex;
  CALCDEVIANCERESIDUALSFUNC calcDevianceResiduals;
  CALCDEVIANCERESIDUALSFUNC calcDevianceResidualsComplex;
  double* prevBeta;
  double* weights;
  double* sqrtWeights;
  double* z;
  double* eta;
  double* scaledZ;
  double* scaledXMatrix;
  int nfactorsAllocated;
  int nobsAllocated;
  int complexFlagAllocated;
  int lastFitValid;
} IRLSData;

static char rcsid[] = "$Id: glm_irls.c,v 1.5 2008/01/29 00:57:21 welling Exp $";

static char err_buf[512]; /* not used for all messages */


static void irls_destroy_self(Regressor* r)
{
  if (r->hook) {
    IRLSData* data= (IRLSData*)(r->hook);
    if (data->innerRegressor) glm_destroy(data->innerRegressor);
    if (data->prevBeta) free(data->prevBeta);
    if (data->weights) free(data->weights);
    if (data->sqrtWeights) free(data->sqrtWeights);
    if (data->z) free(data->z);
    if (data->eta) free(data->eta);
    if (data->scaledZ) free(data->scaledZ);
    if (data->scaledXMatrix) free(data->scaledXMatrix);
  }
  glm_base_destroy_self(r);
}

static int irls_n_params(Regressor* r, int nfactors)
{
  int sum;
  sum= (r->get(r,GLM_COMPLEX)) ? 2*nfactors : nfactors;
  if (r->get(r,GLM_COVARIANCES)) {
    if (r->get(r,GLM_COMPLEX)) sum += 2*2*nfactors*nfactors;
    else sum += nfactors*nfactors;
  }
  if (r->get(r,GLM_DEVIANCE)) {
    sum += 1;
  }
  return sum;
}

static int irls_is_settable( Regressor* r, glm_feature feature )
{
  switch ((int)feature) {
  case GLM_COVARIANCES:
  case GLM_DEVIANCE:
  case GLM_RESIDUALS:
    return 1;
  default:
    return glm_base_is_settable(r, feature);
  }
}

static int irls_context_valid( Regressor* r, int nobs, int nfactors )
{
  IRLSData* data= (IRLSData*)(r->hook);

  if (!data->lastFitValid) {
    _glm_err_txt= "No previous fit, or previous fit failed";
    return 0;
  }
  if (!data->nobsAllocated || !data->nfactorsAllocated) {
    _glm_err_txt= "No previous call to glm_fit";
    return 0;
  }
  if (nobs != data->nobsAllocated) {
    sprintf(err_buf,
	    "Observation dimension %d does not match that of last fit (%d)!",
	    nobs, data->nobsAllocated);
    _glm_err_txt= err_buf;
    return 0;
  }
  if (nfactors != data->nfactorsAllocated) {
    sprintf(err_buf,
	    "Factor dimension %d does not match that of last fit (%d)!",
	    nfactors, data->nfactorsAllocated);
    _glm_err_txt= err_buf;
    return 0;
  }
  if (r->get(r,GLM_COMPLEX) != data->complexFlagAllocated) {
    /* This can happen if glm_set has been used since the last fit. */
    _glm_err_txt= "Complex flag has been changed since last fit!";
    return 0;
  }

  return ( data->innerRegressor->context_valid(data->innerRegressor,
					       nobs, nfactors) 
	   && glm_base_context_valid(r, nobs, nfactors) );
}

static void irls_check_memory(Regressor* r, int nobs, int nfactors) 
{
  IRLSData* data= (IRLSData*)(r->hook);

  if ((nobs>data->nobsAllocated) || (nfactors>data->nfactorsAllocated) 
      || (r->get(r,GLM_COMPLEX) != data->complexFlagAllocated)) {
    int complex_fac= (r->get(r,GLM_COMPLEX) ? 2 : 1);

    if (data->prevBeta) free(data->prevBeta);
    if (!(data->prevBeta=(double*)malloc(nfactors*sizeof(double))))
      MALLOC_FAILURE(nfactors,double);
    
    if (data->weights) free(data->weights);
    if (!(data->weights=(double*)malloc(nobs*sizeof(double))))
      MALLOC_FAILURE(nobs,double);
    
    if (data->sqrtWeights) free(data->sqrtWeights);
    if (!(data->sqrtWeights=(double*)malloc(nobs*sizeof(double))))
      MALLOC_FAILURE(nobs,double);
    
    if (data->z) free(data->z);
    if (!(data->z=(double*)malloc(nobs*sizeof(double))))
      MALLOC_FAILURE(nobs,double);
    
    if (data->eta) free(data->eta);
    if (!(data->eta=(double*)malloc(nobs*sizeof(double))))
      MALLOC_FAILURE(nobs,double);
    
    if (data->scaledZ) free(data->scaledZ);
    if (!(data->scaledZ=(double*)malloc(nobs*sizeof(double))))
      MALLOC_FAILURE(nobs,double);
    
    if (data->scaledXMatrix) free(data->scaledXMatrix);
    if (!(data->scaledXMatrix=(double*)malloc(nobs*nfactors*sizeof(double))))
      MALLOC_FAILURE(nobs*nfactors,double);
    
    data->nobsAllocated= nobs;
    data->nfactorsAllocated= nfactors;
    data->complexFlagAllocated= r->get(r,GLM_COMPLEX);
  }
}

static void logit_calcInitialWeights(double* weights, int nobs)
{
  int i;
  for (i=0; i<nobs; i++) weights[i]= 1.0;
}

static void log_calcInitialWeights(double* weights, int nobs)
{
  int i;
  for (i=0; i<nobs; i++) weights[i]= 1.0;
}

static void logit_calcInitialZ(double* z, const double* y, const double* n,
			       int nobs)
{
  int i;
  for (i=0; i<nobs; i++) z[i]= log( (y[i]+0.5)/((n[i]-y[i])+0.5) );
}

static void log_calcInitialZ(double* z, const double* y, const double* n,
			       int nobs)
{
  int i;
  for (i=0; i<nobs; i++) z[i]= log( y[i]+0.5 );
#ifdef never
  fprintf(stderr,"Initial z: ");
  for (i=0; i<nobs; i++) fprintf(stderr,"%g, ",z[i]);
  fprintf(stderr,"\n");
#endif
}

static void copy_and_apply_scale( double* data_out, const double* data, 
				  const double* scale, 
				  int complex_flag, int n, int ncols )
{
  int icol;
  int i;
  int offset;
  
  for (icol=0; icol<ncols; icol++) {
    if (complex_flag) {
      offset= 2*icol;
      for (i=0; i<n; i++) {
	data_out[offset]= data[offset]*scale[i];
	data_out[offset+1]= data[offset+1]*scale[i];
	offset += 2*ncols;
      }
    }
    else {
      offset= icol;
      for (i=0; i<n; i++) {
	data_out[offset]= data[offset] * scale[i];
	offset += ncols;
      }
    }
  }
}

static int calcBeta(double* beta, 
		    const double* xMat, double* scaledXMat,
		    const double* weights, double* sqrtWeights, 
		    const double* z, double* scaledZ,
		    int nfactors, int nobs, Regressor* llsqRegressor)
{
  int i;
  int a;

  for (i=0; i<nobs; i++) {
    sqrtWeights[i]= sqrt(weights[i]);
    scaledZ[i]= sqrtWeights[i]*z[i];
  }

  copy_and_apply_scale(scaledXMat, xMat, sqrtWeights, 0, nobs, nfactors);
  glm_clear_error_msg();
  if (llsqRegressor->fit(llsqRegressor, scaledZ, scaledXMat, NULL,
			 beta, nobs, nfactors)) {
    /* Regression failed- not a good sign */
    snprintf(err_buf,sizeof(err_buf),"An iterative regression failed: %s",
	     glm_error_msg());
    _glm_err_txt= err_buf;
    return 1;
  }
  /* Check for non-finite values */
  for (a=0; a<nfactors; a++)
    if (!finite(beta[a])) {
      if (glm_error_msg()) {
	snprintf(err_buf,sizeof(err_buf),"An iterative regression failed: %s",
	     glm_error_msg());
      }
      else {
	snprintf(err_buf,sizeof(err_buf),
		 "Iterative regression seems to have diverged");
      }
      _glm_err_txt= err_buf;
      return 1;
    }
  return 0;
}

static void calcEta( double* eta, const double* xMatrix, const double* beta, 
		     int nfactors, int nobs )
{
  int i;
  int a;
  double one= 1.0;
  double zero= 0.0;
  int int_one= 1;

  DGEMV("t", &nfactors, &nobs, &one, (double*)xMatrix, &nfactors, 
	(double*)beta, &int_one,
	&zero, eta, &int_one);
#ifdef never
  fprintf(stderr,"Eta: ");
  for (i=0; i<nobs; i++) fprintf(stderr,"%g, ",eta[i]);
  fprintf(stderr,"\n");
#endif
}

static void logit_calcZ( double* z, const double* n, const double* y, 
			 const double* eta, int nobs )
{
  int i;

  for (i=0; i<nobs; i++) {
    double pi_i= exp(eta[i])/(1.0+exp(eta[i])); /* inverse logit */
    z[i]= eta[i] + ((y[i]-(n[i]*pi_i))/(n[i]*pi_i*(1.0-pi_i)));
  }
#ifdef never
  fprintf(stderr,"Z: ");
  for (i=0; i<nobs; i++) fprintf(stderr,"%g, ",z[i]);
  fprintf(stderr,"\n");
#endif
}

static void log_calcZ( double* z, const double* n, const double* y, 
		       const double* eta, int nobs )
{
  int i;

  for (i=0; i<nobs; i++) {
    double mu_i= exp(eta[i]); /* inverse log-linear link */
    z[i]= eta[i] + ((y[i]-mu_i)/mu_i);
  }
#ifdef never
  fprintf(stderr,"Z: ");
  for (i=0; i<nobs; i++) fprintf(stderr,"%g, ",z[i]);
  fprintf(stderr,"\n");
#endif
}

static void logit_calcW( double* w, const double* n, const double* eta, 
			 int nobs )
{
  int i;
  for (i=0; i<nobs; i++) {
    double pi_i= exp(eta[i])/(1.0+exp(eta[i])); /* inverse logit */
    w[i]= n[i]*pi_i*(1.0-pi_i);
  }
#ifdef never
  fprintf(stderr,"W: ");
  for (i=0; i<nobs; i++) fprintf(stderr,"%g, ",w[i]);
  fprintf(stderr,"\n");
#endif
}

static void log_calcW( double* w, const double* n, const double* eta, 
			 int nobs )
{
  int i;
  for (i=0; i<nobs; i++) {
    double mu_i= exp(eta[i]); /* inverse log-linear link */
    w[i]= mu_i;
  }
#ifdef never
  fprintf(stderr,"W: ");
  for (i=0; i<nobs; i++) fprintf(stderr,"%g, ",w[i]);
  fprintf(stderr,"\n");
#endif
}

static int converged(const double* beta, const double* oldBeta, int nfactors)
{
  double sum_squares= 0.0;
  double sum_diff_squares= 0.0;
  int a;

  for (a=0; a<nfactors; a++) {
    double diff= beta[a]-oldBeta[a];
    sum_squares += beta[a]*beta[a];
    sum_diff_squares= diff*diff;
  }
  /* We should do some clever log-based convergence test, 
   * but for the moment...
   */
  return (sum_diff_squares/sum_squares <= DLAMCH("e"));
}

static int calcCovariances(Regressor* r, double* buf_out, int nobs,
			   int nfactors)
{
  IRLSData* data= (IRLSData*)(r->hook);

  if (!(data->innerRegressor->getXtXInv(data->innerRegressor, buf_out,
					nobs, nfactors))) {
    snprintf(err_buf,sizeof(err_buf),
	     "Error getting covariances: inner regressor says: %s",
	     glm_error_msg());
    _glm_err_txt= err_buf;
    return 0;
  }
  return 1;
}

static double logit_calcDeviance( const double* obs, const double* counts, 
				  const double* eta, int nobs)
{
  double sum= 0.0;
  int i;

  for (i=0; i<nobs; i++) {
    double pi_i= exp(eta[i])/(1.0+exp(eta[i])); /* inverse logit */
    double mu_i= counts[i]*pi_i;
    double term= obs[i]*log(obs[i]/mu_i) 
      + (counts[i]-obs[i])*log((counts[i]-obs[i])/(counts[i]-mu_i));
    sum += term;
  }

  return 2.0*sum;
}

static double log_calcDeviance( const double* obs, const double* counts, 
				  const double* eta, int nobs)
{
  double sum= 0.0;
  int i;

  for (i=0; i<nobs; i++) {
    double mu_i= exp(eta[i]); /* inverse log-linear link */
    double term= obs[i]*log(obs[i]/mu_i) - (obs[i]-mu_i);
    sum += term;
  }

  return 2.0*sum;
}

static double unimp_calcDevianceComplex( const double* obs, 
					 const double* counts, 
					 const double* eta, int nobs )
{
  fprintf(stderr,"calcDevianceComplex: not implemented!\n");
  exit(-1);
  return 0.0;
}

static void logit_calcDevianceResiduals( double* obs, const double* counts, 
					 const double* eta, int nobs)
{
  /* This overwrites obs! */
  int i;

  for (i=0; i<nobs; i++) {
    double pi_i= exp(eta[i])/(1.0+exp(eta[i])); /* inverse logit */
    double mu_i= counts[i]*pi_i;
    double term= obs[i]*log(obs[i]/mu_i) 
      + (counts[i]-obs[i])*log((counts[i]-obs[i])/(counts[i]-mu_i));
    if (obs[i]>mu_i)
      obs[i]= sqrt(2.0*term);
    else
      obs[i]= -sqrt(2.0*term);
  }
}

static void log_calcDevianceResiduals( double* obs, const double* counts, 
				       const double* eta, int nobs)
{
  /* This overwrites obs! */
  int i;

  for (i=0; i<nobs; i++) {
    double mu_i= exp(eta[i]); /* inverse log-linear link */
    double term= obs[i]*log(obs[i]/mu_i) - (obs[i]-mu_i);
    if (obs[i]>mu_i)
      obs[i]= sqrt(2.0*term);
    else
      obs[i]= -sqrt(2.0*term);
  }
}

static void unimp_calcDevianceResidualsComplex( double* obs, 
						const double* counts, 
						const double* eta,
						int nobs )
{
  fprintf(stderr,"calcDevianceResidualsComplex: not implemented!\n");
  exit(-1);
}

static int irls_fit(Regressor* r, double* obs, const double* factors,
			const double* counts, double* param_out, 
			int nobs, int nfactors)
{
  int nparam;
  double* param_offset;
  double *bvals= NULL;;
  double *covariances= NULL;
  double *deviance= NULL;
  IRLSData* data= (IRLSData*)(r->hook);
  int iteration= 0;

  if (nobs<nfactors) {
    _glm_err_txt= "Number of factors is greater than number of observations!";
    return 1;
  }

  /* Breakdown of parameter vector */
  bvals= param_out;
  if (r->get(r,GLM_COMPLEX)) param_offset= param_out + (2*nfactors);
  else param_offset= param_out + nfactors;
  if (r->get(r,GLM_COVARIANCES)) {
    covariances= param_offset;
    if (r->get(r,GLM_COMPLEX)) param_offset += 2*2*nfactors*nfactors;
    else param_offset += nfactors*nfactors;
  }
  if (r->get(r,GLM_DEVIANCE)) {
    deviance= param_offset;
    param_offset += 1;
  }
  nparam= param_offset - param_out;

  irls_check_memory(r, nobs, nfactors);

  data->calcInitialWeights(data->weights, nobs);
  data->calcInitialZ(data->z, obs, counts, nobs);

  if (calcBeta(bvals, factors, data->scaledXMatrix, 
	       data->weights, data->sqrtWeights, 
	       data->z, data->scaledZ,
	       nfactors, nobs, data->innerRegressor)) {
    /* Something went wrong in the inner iterative regression.
     *  Error message should already have been set by the inner regressor.
     */
    return 1;
  }
  if (CHECK_DEBUG(r)){
    int a;
    fprintf(stderr,"Iteraion %d: Beta: ", iteration);
    for (a=0; a<nfactors; a++) fprintf(stderr,"%g, ",bvals[a]);
    fprintf(stderr,"\n");
  }

  while (iteration==0 || !converged(bvals,data->prevBeta,nfactors)) {
    int a;
    iteration++;
    for (a=0; a<nfactors; a++) data->prevBeta[a]= bvals[a];
    calcEta(data->eta, factors, bvals, nfactors, nobs);
    data->calcZ(data->z, counts, obs, data->eta, nobs);
    data->calcW(data->weights, counts, data->eta, nobs);
    if (calcBeta(bvals, factors, data->scaledXMatrix, 
		 data->weights, data->sqrtWeights, 
		 data->z, data->scaledZ,
		 nfactors, nobs, data->innerRegressor)) {
      /* Something went wrong in the inner iterative regression */
      return 1;
    }
    if (CHECK_DEBUG(r)) {
      int a;
      fprintf(stderr,"Iteraion %d: Beta: ", iteration);
      for (a=0; a<nfactors; a++) fprintf(stderr,"%g, ",bvals[a]);
      fprintf(stderr,"\n");
    }
    if (iteration>=MAX_ITERATIONS) {
      if (CHECK_DEBUG(r)) fprintf(stderr,"Convergence failed!\n");
      snprintf(err_buf,sizeof(err_buf),
	       "Convergence of inner iteration failed!");
      _glm_err_txt= err_buf;
      return 1;
    }
  }
  if (CHECK_DEBUG(r)) {
    int i;
    fprintf(stderr,"Converged in %d iterations!\n",iteration);
    fprintf(stderr,"Final weights: ");
    for (i=0; i<nobs; i++) fprintf(stderr,"%g, ",data->weights[i]);
    fprintf(stderr,"\n");
  }
  data->lastFitValid= 1;

  if (r->get(r,GLM_COVARIANCES)) {
    if (!calcCovariances(r,covariances,nobs,nfactors)) 
      return 1;
  }

  if (r->get(r,GLM_DEVIANCE)) {
    if (r->get(r,GLM_COMPLEX)) 
      *deviance= data->calcDevianceComplex(obs, counts, data->eta, nobs);
    else
      *deviance= data->calcDeviance(obs, counts, data->eta, nobs);
  }

  /* This has to come last because it overwrites obs! */
  if (r->get(r,GLM_RESIDUALS)) {
    if (r->get(r,GLM_COMPLEX)) 
      data->calcDevianceResidualsComplex(obs, counts, data->eta, nobs);
    else
      data->calcDevianceResiduals(obs, counts, data->eta, nobs);
  }

  return 0;
}

static Regressor* create_base_irls_regressor()
{
  Regressor* result= glm_create_base_regressor();
  IRLSData* data= NULL;
 
  result->is_settable= irls_is_settable;
  result->n_params= irls_n_params;
  result->destroy_self= irls_destroy_self;
  result->context_valid= irls_context_valid;
  result->fit= irls_fit;
  /* Leave project and normproject with base methods, 
   * which fail, until we get around to implementing them.
   */

  if (!(data=(IRLSData*)malloc(sizeof(IRLSData))))
    MALLOC_FAILURE(sizeof(IRLSData),byte);
  result->hook= data;

  data->nfactorsAllocated= 0;
  data->nobsAllocated= 0;
  data->complexFlagAllocated= 0;
  data->lastFitValid= 0;
  data->innerRegressor= glm_create_llsq_regressor();
  data->calcInitialWeights= NULL;
  data->calcInitialZ= NULL;
  data->calcZ= NULL;
  data->calcW= NULL;
  data->calcDeviance= NULL;
  data->calcDevianceComplex= NULL;
  data->calcDevianceResiduals= NULL;
  data->calcDevianceResidualsComplex= NULL;
  data->prevBeta= NULL;
  data->weights= NULL;
  data->sqrtWeights= NULL;
  data->z= NULL;
  data->eta= NULL;
  data->scaledZ= NULL;
  data->scaledXMatrix= NULL;

  return result;
}

Regressor* glm_create_logistic_regressor()
{
  Regressor* result= create_base_irls_regressor();
  IRLSData* data= (IRLSData*)(result->hook);
 
  result->feature_array[GLM_TYPE]= GLM_TYPE_LOGISTIC;

  data->calcInitialWeights= logit_calcInitialWeights;
  data->calcInitialZ= logit_calcInitialZ;
  data->calcZ= logit_calcZ;
  data->calcW= logit_calcW;
  data->calcDeviance= logit_calcDeviance;
  data->calcDevianceComplex= unimp_calcDevianceComplex;
  data->calcDevianceResiduals= logit_calcDevianceResiduals;
  data->calcDevianceResidualsComplex= unimp_calcDevianceResidualsComplex;

  return result;
}

Regressor* glm_create_poisson_regressor()
{
  Regressor* result= create_base_irls_regressor();
  IRLSData* data= (IRLSData*)(result->hook);
 
  result->feature_array[GLM_TYPE]= GLM_TYPE_POISSON;

  data->calcInitialWeights= log_calcInitialWeights;
  data->calcInitialZ= log_calcInitialZ;
  data->calcZ= log_calcZ;
  data->calcW= log_calcW;
  data->calcDeviance= log_calcDeviance;
  data->calcDevianceComplex= unimp_calcDevianceComplex;
  data->calcDevianceResiduals= log_calcDevianceResiduals;
  data->calcDevianceResidualsComplex= unimp_calcDevianceResidualsComplex;

  return result;
}

