/************************************************************
 *                                                          *
 *  glm.c                                                   *
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
#include "lapack.h"
#include <glm.h>

#ifdef DARWIN
#define finite( foo ) isfinite( foo )
#endif

/* Coding notes-
 * -It's really unfortunate that the return code convention (whether 1 or 0
 *  is the success code) seems to differ between the fitting routines and
 *  the parameter set/test routines.
 */

/* Usage Notes-
 * -To fit a constant parameter, the user must include a column of 1's
 *  in the factor matrix, increasing the number of factors given by 1.
 * -be sure to use glm_set to request particular outputs *before* using
 *  glm_n_params to find number of parameters returned.
 * -options:
 *    GLM_COMPLEX:     data and parameters (but not factors) are complex.
 *    GLM_RESIDUALS:   calculate residuals.
 *    GLM_VARIANCES:   calculate variances.
 *    GLM_COVARIANCES: calculate covariances.
 *    GLM_SSQR:        calculate sums of squares.
 *    GLM_ORTHO:       calculate factor orthogonality measure.
 * -parameters returned: 
 *  if complex_flag is set, params contains:
 *   2*nfactors entries for the b value estimates
 *   2*nfactors entries for the variances of the b values, if requested
 *   2*2*nfactors*nfactors entries for the covariances of the b values,
 *     if requested
 *   1 entry for SSTO, if requested
 *   1 entry for (1/nobs)(sum(obs)*sum(obs)) (the average term) if requested
 *   nfactors entries for the SSR components associated w/ each factor if req.
 *  if complex_flag is not set, params contains:
 *   1*nfactors entries for the b value estimates
 *   1*nfactors entries for the variances of the b values if requested
 *   1*nfactors*nfactors entries for the covariances of the b values 
 *     if requested
 *   1 entry for SSTO if requested
 *   1 entry for (1/nobs)(sum(obs)*sum(obs)) (the average term) if requested
 *   nfactors entries for the SSR components associated w/ each factor if req.
 *   one entry for the ratio of the determinant of (X transpose)X where X is
 *     the factor matrix to the product of squared diagonal elements of the
 *     factor matrix, as a measure of orthogonality of the factors, if 
 *     requested.
 *
 * Implementation Notes-
 * -subscripts in the range i,j,k,etc. denote values 0<i<nobs, or
 *  values the range of which is intended to be general (like for matrix
 *  math).  Subscripts in the range a,b,c,etc. denote values 0<a<nfactors.
 */

#define MALLOC_FAILURE(num,type)\
  { fprintf(stderr,"%s:%d: %s: unable to allocate %d %ss\n",\
            __FILE__,__LINE__,__func__,num,#type); exit(-1); }

#define CHECK_DEBUG(r) (r->get(r,GLM_DEBUG))

typedef struct regressor_llsq_data_struct {
  double* factors_or_u_ftn;
  double* v_ftn;
  double* singular;
  double* inv_singular;
  double* residuals;
  double* covariance_ftn;
  double* rwork;
  int lwork;
  int nfactorsAllocated;
  int nobsAllocated;
  int complexFlagAllocated;
  int lastFitValid;
} LLSqData;

static char rcsid[] = "$Id: glm.c,v 1.29.2.6 2008/01/29 00:57:21 welling Exp $";

char* _glm_err_txt= NULL;
static char err_buf[512]; /* not used for all messages */

/* Must correspond to the glm_feature enum! */
static char* featureNames[]= {
  "complex",
  "residuals",
  "variances",
  "covariances",
  "ssqr",
  "ortho",
  "debug",
  "deviance",
  "type" /* must be last */
};
#define N_FEATURE_NAMES (sizeof(featureNames)/sizeof(char*))

/* Must correspond to glm_type enum! */
static char* regressorTypeNames[]= {
  "llsq",
  "logistic",
  "poisson",
  "base",
};

static int mach_prec_initialized= 0;
static double mach_precision;

const char* glm_get_feature_name(glm_feature feature) 
{
  if ((int)feature>=0 && (int)feature<=(int)GLM_TYPE)
    return featureNames[(int)feature];
  else {
    snprintf(err_buf,sizeof(err_buf),"Unknown regressor feature id %d!",
	     (int)feature);
    _glm_err_txt= err_buf;
    return NULL;

  }
  return NULL; /* not reached */
}

const char* glm_get_regressor_type_name(glm_type type)
{
  if ((int)type>=0 && (int)type<=(int)GLM_TYPE_BASE)
    return regressorTypeNames[(int)type];
  else {
    snprintf(err_buf,sizeof(err_buf),"Unknown regressor type id %d!",
	     (int)type);
    _glm_err_txt= err_buf;
    return NULL;
  }
}

Regressor* glm_create_regressor_by_type(glm_type type)
{
  Regressor* result= NULL;
  switch (type) {
  case GLM_TYPE_BASE: return glm_create_base_regressor();
  case GLM_TYPE_LLSQ: return glm_create_llsq_regressor();
  case GLM_TYPE_LOGISTIC: return glm_create_logistic_regressor();
  case GLM_TYPE_POISSON: return glm_create_poisson_regressor();
  default:
    snprintf(err_buf,sizeof(err_buf),"Unknown regressor type id %d!",
	     (int)type);
    _glm_err_txt= err_buf;
    return NULL;
  }
  return NULL; /* not reached */
}

char* glm_error_msg()
{
  return _glm_err_txt;
}

void glm_clear_error_msg()
{
  _glm_err_txt= NULL;
}

int glm_fit(Regressor* r, double* obs, const double* factors, 
	    const double* counts, double* param_out,
	    int nobs, int nfactors)
{
  return( (r->fit)(r, obs, factors, counts, param_out, nobs, nfactors) );
}

int glm_project(Regressor* r, const double* obs_vec_in, 
		double* factor_vec_out, int nobs, int nfactors)
{
  return( (r->project)(r, obs_vec_in, factor_vec_out, nobs, nfactors) );
}

int glm_normproject(Regressor* r, const double* obs_vec_in,
		    double* factor_vec_out, int nobs, int nfactors)
{
  return( (r->normproject)(r, obs_vec_in, factor_vec_out, nobs, nfactors) );
}

int glm_context_valid( Regressor* r, int nobs, int nfactors )
{
  return( (r->context_valid)(r, nobs, nfactors) );
}

int glm_is_settable( Regressor* r, glm_feature feature )
{
  return( (r->is_settable)(r, feature) );
}

void glm_set(Regressor* r, glm_feature feature, int flag) 
{
  (r->set)(r,feature,flag);
}

int glm_get(Regressor* r, glm_feature feature)
{
  return( (r->get)(r,feature) );
}

int glm_n_params(Regressor* r, int nfactors)
{
  return (r->n_params)(r,nfactors);
}

void glm_destroy(Regressor* target)
{
  target->destroy_self(target);
}

int glm_base_context_valid( Regressor* r, int nobs, int nfactors )
{
  return 1; 
}

void glm_base_destroy_self(Regressor* r)
{
  if (r->hook) free(r->hook);
  free(r);
}

int glm_base_is_settable( Regressor* r, glm_feature feature )
{
  if (feature==GLM_DEBUG) return 1;
  else return 0;
}

static void base_set(Regressor* r, glm_feature feature, int flag)
{
  if (r->is_settable(r, feature))
    r->feature_array[feature]= flag;
  else {
    if (feature==GLM_TYPE) {
      fprintf(stderr,"Cannot change the type of an existing Regressor!\n");
      exit(-1);
    }
    else {
      if (feature<0) {
	fprintf(stderr,"glm.c: invalid feature ID %d!\n",feature);
	exit(-1);
      }
      else if (feature>=GLM_TYPE) {
	/* GLM_TYPE must be last in table! */
	snprintf(err_buf,sizeof(err_buf),"unknown feature ID %d",feature);
	_glm_err_txt= err_buf;
      }
      else {
	assert(N_FEATURE_NAMES==GLM_N_FEATURES);
	snprintf(err_buf,sizeof(err_buf),
		 "feature %s is not settable for regressor type %s",
		 featureNames[feature],
		 regressorTypeNames[r->feature_array[GLM_TYPE]]);
	_glm_err_txt= err_buf;
      }
    }
  }
}

static int base_get(Regressor* r, glm_feature feature)
{
  if (feature<0) {
    fprintf(stderr,"glm.c: invalid feature ID %d!\n",feature);
    exit(-1);
  }
  else if (feature>GLM_TYPE) { /* GLM_TYPE must be last */
    snprintf(err_buf,sizeof(err_buf),"unknown feature ID %d",feature);
    _glm_err_txt= err_buf;
    return 0;
  }
  else return r->feature_array[feature];
}

static int base_n_params(Regressor* r, int nfactors)
{
  int sum;
  sum= (r->get(r,GLM_COMPLEX)) ? 2*nfactors : nfactors;
  return sum;
}

static int base_fit(Regressor* r, double* obs, const double* factors,
		    const double* counts, double* param_out, 
		    int nobs, int nfactors)
{
  fprintf(stderr,"glm.c: called 'fit' method on base Regressor class!\n");
  exit(-1);
  return 0;
}

static int base_project(Regressor* r, const double* obs_vec_in, 
			double* factor_vec_out, int nobs, int nfactors)
{
  fprintf(stderr,"glm.c: called 'project' method on base Regressor class!\n");
  exit(-1);
  return 0;
}

static int base_normproject(Regressor* r, const double* obs_vec_in, 
			    double* factor_vec_out, int nobs, int nfactors)
{
  fprintf(stderr,
	  "glm.c: called 'normproject' method on base Regressor class!\n");
  exit(-1);
  return 0;
}

static int base_getXtXInv(Regressor* r, double* buf_out, 
			  int nobs, int nfactors)
{
  fprintf(stderr,
	  "glm.c: called 'getXtXInv' method on base Regressor class!\n");
  exit(-1);
  return 0;
}

Regressor* glm_create_base_regressor(void)
{
  Regressor* result= NULL;
  int i;

  if (!(result=(Regressor*)malloc(sizeof(Regressor)))) 
    MALLOC_FAILURE(sizeof(Regressor),byte);
  result->destroy_self= glm_base_destroy_self;
  result->context_valid= glm_base_context_valid;
  result->is_settable= glm_base_is_settable;
  result->set= base_set;
  result->get= base_get;
  result->n_params= base_n_params;
  result->fit= base_fit;
  result->project= base_project;
  result->normproject= base_normproject;
  result->getXtXInv= base_getXtXInv;
  result->feature_array[GLM_TYPE]= GLM_TYPE_BASE;
  for (i=0; i<GLM_TYPE; i++) result->feature_array[i]= 0;
  result->hook= NULL;
  return result;
}

static void flip_to_ftn( const double* in, double* out, int xdim, int ydim )
{
  int i;
  int j;

  for (i=0; i<xdim; i++) 
    for (j=0; j<ydim; j++) out[(j*xdim)+i]= in[(i*ydim)+j];
}

static void flip_to_c( double* in, double* out, int xdim, int ydim )
{
  flip_to_ftn(in, out, ydim, xdim );
}

static void llsq_check_memory(Regressor* r, int nobs, int nfactors) 
{
  LLSqData* data= (LLSqData*)(r->hook);

  if ((nobs>data->nobsAllocated) || (nfactors>data->nfactorsAllocated) 
      || (r->get(r,GLM_COMPLEX) != data->complexFlagAllocated)) {
    int complex_fac= (r->get(r,GLM_COMPLEX) ? 2 : 1);
    
    if (data->factors_or_u_ftn) free(data->factors_or_u_ftn);
    if (!(data->factors_or_u_ftn=
	  (double*)malloc(nobs*nfactors*sizeof(double))))
      MALLOC_FAILURE(nobs*nfactors,double);
      
    if (data->v_ftn) free(data->v_ftn);
    if (!(data->v_ftn=
	  (double*)malloc(nfactors*nfactors*sizeof(double))))
      MALLOC_FAILURE(nfactors*nfactors,double);

    if (data->singular) free(data->singular);
    if (!(data->singular= (double*)malloc(nfactors*sizeof(double))))
      MALLOC_FAILURE(nfactors,double);

    if (data->inv_singular) free(data->inv_singular);
    if (!(data->inv_singular= (double*)malloc(nfactors*sizeof(double))))
      MALLOC_FAILURE(nfactors,double);

    if (data->residuals) free(data->residuals);
    if (!(data->residuals= (double*)malloc(nobs*complex_fac*sizeof(double))))
      MALLOC_FAILURE(nobs*complex_fac,double)

    if (data->covariance_ftn) free(data->covariance_ftn);
    if (!(data->covariance_ftn=
	  (double*)malloc(complex_fac*complex_fac*nfactors*nfactors*sizeof(double)))) 
      MALLOC_FAILURE(complex_fac*complex_fac*nfactors*nfactors,double);

    data->lwork= 5*nfactors+nobs; /* upper limit of what everyone needs */
    if (data->rwork) free(data->rwork);
    if (!(data->rwork= (double*)malloc(data->lwork*sizeof(double))))
      MALLOC_FAILURE(data->lwork,double);

    data->nobsAllocated= nobs;
    data->nfactorsAllocated= nfactors;
    data->complexFlagAllocated= r->get(r,GLM_COMPLEX);
  }
}

static void zero_singular_coeffs( const double* data_in, double* data_out,
				  int n, int nsamples )
{
  double dlim;
  int i;

  /* Numerical Recipes suggests zeroing out SVD singular coeffs with
   * value ratios less than sqrt(nsamples)*machine precision.
   */

  if (!mach_prec_initialized) {
    mach_precision= DLAMCH("p");
    mach_prec_initialized= 1;
  }

  dlim= data_in[0];
  for (i=1; i<n; i++) if (data_in[i]>dlim) dlim= data_in[i];
  dlim *= sqrt((double)nsamples) * mach_precision;
  for (i=0; i<n; i++) {
    if (!finite(data_in[i]) || data_in[i]<dlim) data_out[i]= 0.0;
    else data_out[i]= 1.0/data_in[i];
  }
}

static void llsq_calc_b(double* bvals, double* v_ftn, double* inv_singular, 
			double* factors_or_u_ftn, const double* obs, 
			int nfactors, int nobs,
			int complex_flag)
{
  int a;
  int b;
  int i;
  double t1;
  double t1_i;
  double t2;
  double t2_i;

  if (complex_flag) {
    for (a=0; a<nfactors; a++) {
      t1= t1_i= 0.0;
      for (b=0; b<nfactors; b++) {
	t2= t2_i= 0.0;
	for (i=0; i<nobs; i++) {
	  t2 += factors_or_u_ftn[(b*nobs)+i]*obs[2*i];
	  t2_i += factors_or_u_ftn[(b*nobs)+i]*obs[(2*i)+1];
	}
	t1 += v_ftn[a*nfactors+b] * inv_singular[b] * t2;
	t1_i += v_ftn[a*nfactors+b] * inv_singular[b] * t2_i;
      }
      bvals[2*a]= t1;
      bvals[(2*a)+1]= t1_i;
    }
  }
  else {
    for (a=0; a<nfactors; a++) {
      t1= 0.0;
      for (b=0; b<nfactors; b++) {
	t2= 0.0;
	for (i=0; i<nobs; i++) t2 += factors_or_u_ftn[(b*nobs)+i]*obs[i];
	t1+= v_ftn[a*nfactors+b] * inv_singular[b] * t2;
      }
      bvals[a]= t1;
    }
  }
}

static void llsq_calc_residuals(double* res, 
				const double* obs, const double* factors, 
				double* bvals, 
				int nfactors, int nobs, int complex_flag) 
{
  int i;
  int j;
  register double t;
  register double t_i;

  if (complex_flag) {
    for (i=0; i<nobs; i++) {
      t= obs[2*i];
      t_i= obs[(2*i)+1];
      for (j=0; j<nfactors; j++) {
	t -= factors[(i*nfactors)+j]*bvals[2*j];
	t_i -= factors[(i*nfactors)+j]*bvals[(2*j)+1];
      }
      res[2*i]= t;
      res[(2*i)+1]= t_i;
    }
  }
  else {
    for (i=0; i<nobs; i++) {
      t= obs[i];
      for (j=0; j<nfactors; j++) 
	t -= factors[(i*nfactors)+j]*bvals[j];
      res[i]= t;
    }
  }

}

static void calc_covariances(double* covariance_ftn, double* residuals, 
			   double* v_ftn, double* inv_singular,
			   int nfactors, int nobs, int complex_flag)
{
  int i;
  int a;
  int b;
  int c;

  if (nobs>nfactors) {
    if (complex_flag) {
      double sigmasqr_rr;
      double sigmasqr_ii;
      double sigmasqr_ri;
      sigmasqr_rr= 0.0;
      sigmasqr_ii= 0.0;
      sigmasqr_ri= 0.0;
      for (i=0; i<nobs; i++) {
	sigmasqr_rr += (residuals[2*i]*residuals[2*i]);
	sigmasqr_ii += (residuals[(2*i)+1]*residuals[(2*i)+1]);
	sigmasqr_ri += (residuals[(2*i)]*residuals[(2*i)+1]);
      }
      sigmasqr_rr /= (double)(nobs-nfactors);
      sigmasqr_ii /= (double)(nobs-nfactors);
      sigmasqr_ri /= (double)(nobs-nfactors);
      for (a=0; a<nfactors; a++)
	for (b=0; b<nfactors; b++) {
	  double sum= 0.0;
	  for (c=0; c<nfactors; c++) 
	    sum += v_ftn[(a*nfactors)+c] * v_ftn[(b*nfactors)+c] 
	      * inv_singular[c] * inv_singular[c];
	  covariance_ftn[2*( 2*((b*nfactors)+a)   )   ]= sum*sigmasqr_rr;
	  covariance_ftn[2*( 2*((b*nfactors)+a)   ) +1]= sum*sigmasqr_ri;
	  covariance_ftn[2*( 2*((b*nfactors)+a)+1 )   ]= sum*sigmasqr_ri;
	  covariance_ftn[2*( 2*((b*nfactors)+a)+1 ) +1]= sum*sigmasqr_ii;
	}
    }
    else {
      double sigmasqr;
      sigmasqr= 0.0;
      for (i=0; i<nobs; i++) sigmasqr += residuals[i]*residuals[i];
      sigmasqr /= (double)(nobs-nfactors);
      for (a=0; a<nfactors; a++)
	for (b=0; b<nfactors; b++) {
	  double sum= 0.0;
	  for (c=0; c<nfactors; c++) 
	    sum += v_ftn[(a*nfactors)+c] * v_ftn[(b*nfactors)+c] 
	      * inv_singular[c] * inv_singular[c];
	  covariance_ftn[(b*nfactors)+a]= sum*sigmasqr;
	}
    }
    
  }
  else {
    /* set them all to zero */
    if (complex_flag) {
      for (a=0; a<nfactors; a++)
	for (b=0; b<nfactors; b++) {
	  covariance_ftn[2*( 2*((b*nfactors)+a)   )   ]= 0.0;
	  covariance_ftn[2*( 2*((b*nfactors)+a)   ) +1]= 0.0;
	  covariance_ftn[2*( 2*((b*nfactors)+a)+1 )   ]= 0.0;
	  covariance_ftn[2*( 2*((b*nfactors)+a)+1 ) +1]= 0.0;
	}
    }
    else {
      for (a=0; a<nfactors; a++)
	for (b=0; b<nfactors; b++)
	  covariance_ftn[(b*nfactors)+a]= 0.0;
    }
  }
}

static void calc_sum_squares(double* ssto, double* y_sqr_by_n, 
			     double* ssr_terms, const double* obs, 
			     double* bvals, double* singular,
			     double* v_ftn, int nobs, int nfactors, 
			     int complex_flag)
{
  int i;
  int j;
  int a;
  int b;
  int c;
  int d;
  register double t1;
  register double t1_i;
  register double t2;
  register double t2_i;
  register double t3;
  register double sum;
  register double sum_i;

  /* y squared by n */
  sum= 0.0;
  if (complex_flag) {
    sum_i= 0.0;
    for (i=0; i<nobs; i++) { 
      sum += obs[2*i];
      sum_i += obs[(2*i)+1];
    }
    *y_sqr_by_n= ((sum*sum)+(sum_i*sum_i))/((double)nobs);
  }
  else {
    for (i=0; i<nobs; i++) sum += obs[i];
    *y_sqr_by_n= (sum*sum)/((double)nobs);
  }
  
  /* SSTO */
  sum= 0.0;
  if (complex_flag) {
    for (i=0; i<nobs; i++) 
      sum += (obs[2*i]*obs[2*i]) + (obs[(2*i)+1]*obs[(2*i)+1]);
  }
  else {
    for (i=0; i<nobs; i++) sum += obs[i]*obs[i];
  }
  *ssto= sum - (*y_sqr_by_n);

  /* Terms of SSR, one for each factor */
  if (complex_flag) {
    for (a=0; a<nfactors; a++) {
      t3= 0.0;
      for (b=0; b<nfactors; b++) {
	t2= 0.0;
	t2_i= 0.0;
	for (c=0; c<=a; c++) {
	  t1= singular[b]*v_ftn[(c*nfactors)+b]*bvals[2*c];
	  t1_i= singular[b]*v_ftn[(c*nfactors)+b]*bvals[(2*c)+1];
	  t2 += t1;
	  t2_i += t1_i;
	}
	t3 += t2*t2 + t2_i*t2_i;
      }
      for (d=0; d<a; d++) t3 -= ssr_terms[d];
      ssr_terms[a]= t3;
      if (ssr_terms[a]<0.0) ssr_terms[a]= 0.0; /* worry about rounding error */
    }
  }
  else {
    for (a=0; a<nfactors; a++) {
      t3= 0.0;
      for (b=0; b<nfactors; b++) {
	t2= 0.0;
	for (c=0; c<=a; c++) {
	  t1= singular[b]*v_ftn[(c*nfactors)+b]*bvals[c];
	  t2 += t1;
	}
	t3 += t2*t2;
      }
      for (d=0; d<a; d++) t3 -= ssr_terms[d];
      ssr_terms[a]= t3;
      if (ssr_terms[a]<0.0) ssr_terms[a]= 0.0; /* worry about rounding error */
    }
  }

}

static double calc_ortho_measure( double* v_ftn, double* singular, int nfactors )
{
  int a;
  int b;
  double eigen_prod= 1.0;
  double trace_prod= 1.0;
  double sum;

  for (a=0; a<nfactors; a++) {
    eigen_prod *= singular[a]*singular[a];
    sum= 0.0;
    for (b=0; b<nfactors; b++) 
      sum += v_ftn[(b*nfactors)+a]*v_ftn[(b*nfactors)+a]
	*singular[b]*singular[b];
    trace_prod *= sum;
  }

  return(eigen_prod/trace_prod);
}

static int llsq_context_valid( Regressor* r, int nobs, int nfactors )
{
  LLSqData* data= (LLSqData*)(r->hook);

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

  return glm_base_context_valid(r, nobs, nfactors);
}

static int llsq_project(Regressor* r, const double* obs_vec_in, 
			double* factor_vec_out, int nobs, int nfactors)
{
  LLSqData* data= (LLSqData*)(r->hook);

  if (!glm_context_valid(r,nobs,nfactors)) return 1;

  /* Decomposition matrices should still be valid, so we can use the same
   * matrix to project the input values in "observation space" down to
   * "factor space".  This is the same operation performed by llsq_calc_b().
   */
  llsq_calc_b(factor_vec_out, data->v_ftn, data->inv_singular, 
	      data->factors_or_u_ftn, 
	      obs_vec_in, nfactors, nobs, r->get(r,GLM_COMPLEX));

  return 0;
}

int llsq_normproject(Regressor* r, const double* obs_vec_in, 
		     double* factor_vec_out, int nobs, int nfactors)
{
  LLSqData* data= (LLSqData*)(r->hook);
  int a;
  int b;
  int c;
  int d;
  int i;
  double t1;
  double t1_i;
  double t2;
  double t2_i;
  double t3;
  double t3_i;

  if (!glm_context_valid(r,nobs,nfactors)) return 1;

  /* Decomposition matrices should still be valid, so we can use the same
   * matrix to project the input values in "observation space" down to
   * "factor space".
   */
  if (r->get(r,GLM_COMPLEX)) {
    for (a=0; a<nfactors; a++) {
      c= a; /* diagonal elements only */
      t1= t1_i= 0.0;
      for (b=0; b<nfactors; b++) {
	t2= t2_i= 0.0;
	for (d=0; d<nfactors; d++) {
	  t3= t3_i= 0.0;
	  for (i=0; i<nobs; i++) {
	    t3 += data->factors_or_u_ftn[(b*nobs)+i]*obs_vec_in[2*i]
	      *obs_vec_in[2*i]*data->factors_or_u_ftn[(d*nobs)+i];
	    t3_i += data->factors_or_u_ftn[(b*nobs)+i]*obs_vec_in[(2*i)+1]
	      *obs_vec_in[(2*i)+1]*data->factors_or_u_ftn[(d*nobs)+i];
	  }
	  t2 += data->v_ftn[c*nfactors+d] * data->inv_singular[d] * t3;
	  t2_i += data->v_ftn[c*nfactors+d] * data->inv_singular[d] * t3_i;
	}
	t1 += data->v_ftn[a*nfactors+b] * data->inv_singular[b] * t2;
	t1_i += data->v_ftn[a*nfactors+b] * data->inv_singular[b] * t2_i;
      }
      factor_vec_out[2*a]= sqrt(t1);
      factor_vec_out[(2*a)+1]= sqrt(t1_i);
    }
  }
  else {
    for (a=0; a<nfactors; a++) {
      c= a; /* diagonal elements only */
      t1= 0.0;
      for (b=0; b<nfactors; b++) {
	t2= 0.0;
	for (d=0; d<nfactors; d++) {
	  t3= 0.0;
	  for (i=0; i<nobs; i++) {
	    t3 += data->factors_or_u_ftn[(b*nobs)+i]*obs_vec_in[i]
	      *obs_vec_in[i]*data->factors_or_u_ftn[(d*nobs)+i];
	  }
	  t2 += data->v_ftn[c*nfactors+d] * data->inv_singular[d] * t3;
	}
	t1 += data->v_ftn[a*nfactors+b] * data->inv_singular[b] * t2;
      }
      factor_vec_out[a]= sqrt(t1);
    }
  }

  return 0;
}

int llsq_fit(Regressor* r, double* obs, const double* factors, 
	     const double* counts, double* param_out, int nobs, int nfactors)
{
  int nparam;
  int i;
  int a,b;
  int lapack_retcode;
  double sing_max;
  double *bvals;
  double *variances;
  double *covariances;
  double *ssto;
  double *y_sqr_by_n;
  double *ssr_terms;
  double *param_offset;
  double *ortho_measure;
  LLSqData* data= (LLSqData*)(r->hook);

  if (nobs<nfactors) {
    _glm_err_txt= "Number of factors is greater than number of observations!";
    return 1;
  }

  /* Breakdown of parameter vector */
  bvals= param_out;
  if (r->get(r,GLM_COMPLEX)) param_offset= param_out + (2*nfactors);
  else param_offset= param_out + nfactors;
  if (r->get(r,GLM_VARIANCES)) {
    variances= param_offset;
    if (r->get(r,GLM_COMPLEX)) param_offset += 2*nfactors;
    else param_offset += nfactors;
  }
  if (r->get(r,GLM_COVARIANCES)) {
    covariances= param_offset;
    if (r->get(r,GLM_COMPLEX)) param_offset += 2*2*nfactors*nfactors;
    else param_offset += nfactors*nfactors;
  }
  if (r->get(r,GLM_SSQR)) {
    ssto= param_offset;
    y_sqr_by_n= ssto+1;
    ssr_terms= y_sqr_by_n + 1;
    param_offset += (nfactors+2);
  }
  if (r->get(r,GLM_ORTHO)) {
    ortho_measure= param_offset++;
  }
  nparam= param_offset - param_out;

  llsq_check_memory(r, nobs, nfactors);

  /* Perform a check for non-finite values in the factors.  For
   * example, this can happen if the user fails to detect an underflow.
   * We'll just initialize everything to the same non-finite value 
   * if we find any.
   */
  for (i=0; i<nobs*nfactors; i++) {
    if (!finite(factors[i])) {
      int j;
      for (j=0; j<nparam; j++) param_out[j]= factors[i];
      data->lastFitValid= 0;
      _glm_err_txt= "An input factor value was not finite";
      return 1;
    }
  }
  
  flip_to_ftn(factors, data->factors_or_u_ftn, nobs, nfactors);

  if (CHECK_DEBUG(r)) {
    /* Some handy diagnostics */
    if (r->get(r,GLM_COMPLEX)) {
      double tmp= 0.0, tmp_i= 0.0;
      for (i=0; i<nobs; i++) {
	tmp += obs[2*i];
	tmp += obs[2*i+1];
	fprintf(stderr,"scaled complex obs %d: %f %f (running sum %f %f)\n",
		i,obs[2*i],obs[2*i+1],tmp,tmp_i);
      }
    }
    else {
      double tmp= 0.0;
      for (i=0; i<nobs; i++) {
	tmp += obs[i];
	fprintf(stderr,"scaled obs %d: %f (running sum %f)\n",i,obs[i],tmp);
      }
    }
  }

  /* Do SVD.  This destroys the factors_or_u_ftn matrix. */
  (void)DGESVD("O", "A", &nobs, &nfactors, data->factors_or_u_ftn,
	       &nobs, data->singular, NULL, &nobs, data->v_ftn, &nfactors,
	       data->rwork, &(data->lwork), &lapack_retcode);
  if (lapack_retcode<0) {
    sprintf(err_buf,
	    "DGESVD argument %d had an illegal value",
	    -lapack_retcode);
    _glm_err_txt= err_buf;
    data->lastFitValid= 0;
    return 1;
  }
  else if (lapack_retcode>0) {
    sprintf(err_buf,
	    "DBDSQR did not converge in SGESVD; %d superdiagonals failed",
	    lapack_retcode);
    _glm_err_txt= err_buf;
    data->lastFitValid= 0;
    return 1;
  }

  /* Zero out appropriate parts of the diagonal matrix */
  zero_singular_coeffs(data->singular, data->inv_singular, nfactors, nobs);

  /* We can now have faith in v_ftn, u_ftn (in factors_or_u_ftn), 
   * and inv_singular 
   */
  data->lastFitValid= 1;

  /* Calculate best fit parameters */
  llsq_calc_b(bvals, data->v_ftn, data->inv_singular, data->factors_or_u_ftn, 
	      obs, nfactors, nobs, r->get(r,GLM_COMPLEX));

  /* Calculate residues.  Have to use factors
   * rather than factors_or_u_ftn because the SVD operation overwrote
   * that version of factors. 
   */
  if (r->get(r,GLM_RESIDUALS) || r->get(r,GLM_VARIANCES) || r->get(r,GLM_COVARIANCES))
    llsq_calc_residuals(data->residuals, obs, factors, bvals, nfactors, nobs, 
			r->get(r,GLM_COMPLEX));

  /* Calculate covariances */
  if (r->get(r,GLM_VARIANCES) || r->get(r,GLM_COVARIANCES)) {
    calc_covariances(data->covariance_ftn, data->residuals, data->v_ftn, 
		     data->inv_singular, nfactors, nobs, r->get(r,GLM_COMPLEX));
    if (r->get(r,GLM_VARIANCES)) {
      if (r->get(r,GLM_COMPLEX)) {
	for (a=0; a<nfactors; a++) {
	  variances[2*a]= data->covariance_ftn[2*(2*((a*nfactors)+a))];
	  variances[(2*a)+1]= data->covariance_ftn[(2*((a*nfactors)+a)+1)+1];
	}
      }
      else {
	for (a=0; a<nfactors; a++) 
	  variances[a]= data->covariance_ftn[(a*nfactors)+a];
      }
    }
    if (r->get(r,GLM_COVARIANCES)) {
      if (r->get(r,GLM_COMPLEX)) {
	for (a=0; a<nfactors; a++) 
	  for (b=0; b<nfactors; b++) {
	    covariances[(2*a)*nfactors+(2*b)]=
	      data->covariance_ftn[(2*b)*nfactors+(2*a)];
	    covariances[(2*a+1)*nfactors+(2*b)]=
	      data->covariance_ftn[(2*b)*nfactors+(2*a+1)];
	    covariances[(2*a)*nfactors+(2*b+1)]=
	      data->covariance_ftn[(2*b+1)*nfactors+(2*a)];
	    covariances[(2*a+1)*nfactors+(2*b+1)]=
	      data->covariance_ftn[(2*b+1)*nfactors+(2*a+1)];
	  }
      }
      else {
	for (a=0; a<nfactors; a++) 
	  for (b=0; b<nfactors; b++)
	    covariances[a*nfactors+b]= data->covariance_ftn[(b*nfactors)+a];
      }
    }
  }

  /* Calculate sums of squares */
  if (r->get(r,GLM_SSQR)) {
    calc_sum_squares(ssto, y_sqr_by_n, ssr_terms, obs, bvals, data->singular,
		     data->v_ftn, nobs, nfactors, r->get(r,GLM_COMPLEX));
  }

  /* Calculate factor orthogonality measure if requested */
  if (r->get(r,GLM_ORTHO)) 
    *ortho_measure= calc_ortho_measure(data->v_ftn, data->singular, nfactors);

  /* Copy in residuals */
  if (r->get(r,GLM_RESIDUALS)) {
    if (r->get(r,GLM_COMPLEX)) {
      for (i=0; i<2*nobs; i++) obs[i]= data->residuals[i];
    }
    else {
      for (i=0; i<nobs; i++) obs[i]= data->residuals[i];
    }
  }

  if (CHECK_DEBUG(r)) {
    /* Some handly diagnostics */
    double tmp= 0.0, tmp_i= 0.0;
    if (r->get(r,GLM_RESIDUALS)) {
      if (r->get(r,GLM_COMPLEX)) {
	for (i=0; i<nobs; i++) 
	  tmp += data->residuals[2*i]*data->residuals[2*i]
	    + data->residuals[2*i+1]*data->residuals[2*i+1];
      }
      else {
	for (i=0; i<nobs; i++) tmp += data->residuals[i]*data->residuals[i];
      }
      fprintf(stderr,"  SSE: %f (by direct sum of residuals)\n",tmp);
    }
    if (r->get(r,GLM_SSQR)) {
      fprintf(stderr,"  SSTO: %f\n",*ssto);
      fprintf(stderr,"  y_sqr_by_n: %f\n",*y_sqr_by_n);
      tmp= 0.0;
      for (i=0; i<nfactors; i++) {
	fprintf(stderr,"   %d: %f\n",i,ssr_terms[i]);
	tmp += ssr_terms[i];
      }
      fprintf(stderr,"  SSR: %f\n",tmp-(*y_sqr_by_n));
    }
    if (r->get(r,GLM_ORTHO))
      fprintf(stderr,"  ortho measure: %f\n",*ortho_measure);
    
    if (r->get(r,GLM_COMPLEX)) {
      tmp= 0.0;
      tmp_i= 0.0;
      for (i=0; i<nobs; i++) {
	double tmp2= 0.0, tmp2_i= 0.0;
	fprintf(stderr,"obs %d: ",i);
	for (b=0; b<nfactors; b++) {
	  if (b != 0) fprintf(stderr," + ");
	  tmp2 += factors[(i*nfactors)+b]*bvals[2*b];
	  tmp2_i += factors[(i*nfactors)+b]*bvals[2*b+1];
	  fprintf(stderr,"%f*( %f %f )",
		  factors[(i*nfactors)+b],bvals[2*b],bvals[2*b+1]);
	  if (!((b+1)%2)) fprintf(stderr,"\n   ");
	}
	if (r->get(r,GLM_RESIDUALS)) {
	  fprintf(stderr,"   = %f %f (residual %f %f)\n",
		  tmp2,tmp2_i,data->residuals[2*i],data->residuals[2*i+1]);
	  tmp += data->residuals[2*i];
	  tmp_i += data->residuals[2*i+1];
	}
	else
	  fprintf(stderr,"   = %f %f\n",tmp2,tmp2_i);
      }
      if (r->get(r,GLM_RESIDUALS))
	fprintf(stderr,"Sum of residuals is %f %f\n",tmp,tmp_i);
    }
    else {
      tmp= 0.0;
      for (i=0; i<nobs; i++) {
	double tmp2= 0.0;
	fprintf(stderr,"obs %d: ",i);
	for (b=0; b<nfactors; b++) {
	  if (b != 0) fprintf(stderr," + ");
	  tmp2 += factors[(i*nfactors)+b]*bvals[b];
	  fprintf(stderr,"%f*%f",factors[(i*nfactors)+b],bvals[b]);
	  if (!((b+1)%3)) fprintf(stderr,"\n   ");
	}
	if (r->get(r,GLM_RESIDUALS)) {
	  tmp += data->residuals[i];
	  fprintf(stderr,"   = %f (residual %f, running sum %f)\n",
		  tmp2,data->residuals[i],tmp);
	}
	else
	  fprintf(stderr,"   = %f\n",tmp2);
      }
      if (r->get(r,GLM_RESIDUALS))
	fprintf(stderr,"Sum of residuals is %f\n",tmp);
    }
  }

  return 0;
}

static void llsq_destroy_self(Regressor* r)
{
  if (r->hook) {
    LLSqData* data= (LLSqData*)(r->hook);
    if (data->factors_or_u_ftn) free(data->factors_or_u_ftn);
    if (data->v_ftn) free(data->v_ftn);
    if (data->singular) free(data->singular);
    if (data->inv_singular) free(data->inv_singular);
    if (data->residuals) free(data->residuals);
    if (data->covariance_ftn) free(data->covariance_ftn);
    if (data->rwork) free(data->rwork);
  }
  glm_base_destroy_self(r);
}

static int llsq_n_params(Regressor* r, int nfactors)
{
  int sum;
  sum= (r->get(r,GLM_COMPLEX)) ? 2*nfactors : nfactors;
  if (r->get(r,GLM_VARIANCES)) {
    if (r->get(r,GLM_COMPLEX)) sum += 2*nfactors;
    else sum += nfactors;
  }
  if (r->get(r,GLM_COVARIANCES)) {
    if (r->get(r,GLM_COMPLEX)) sum += 2*2*nfactors*nfactors;
    else sum += nfactors*nfactors;
  }
  if (r->get(r,GLM_SSQR)) sum += nfactors+2;
  if (r->get(r,GLM_ORTHO)) sum += 1;
  return sum;
}

static int llsq_is_settable( Regressor* r, glm_feature feature )
{
  switch ((int)feature) {
  case GLM_COMPLEX:
  case GLM_RESIDUALS:
  case GLM_VARIANCES:
  case GLM_COVARIANCES:
  case GLM_SSQR:
  case GLM_ORTHO:
    return 1;
  default:
    return glm_base_is_settable(r, feature);
  }
}

static int llsq_getXtXInv(Regressor* r, double* buf_out, 
			  int nobs, int nfactors)
{
  LLSqData* data= (LLSqData*)(r->hook);
  int i;
  int a;
  int b;
  int c;

  if (!(r->context_valid)(r,nobs,nfactors)) {
    snprintf(err_buf,sizeof(err_buf),
	     "Requested (X'X)^-1, but no previous fit set the context!");
    _glm_err_txt= err_buf;
    return 0;
  }
  
  if (nobs>nfactors) {
    for (a=0; a<nfactors; a++)
      for (b=0; b<nfactors; b++) {
	double sum= 0.0;
	for (c=0; c<nfactors; c++) 
	  sum += data->v_ftn[(a*nfactors)+c] * data->v_ftn[(b*nfactors)+c] 
	    * data->inv_singular[c] * data->inv_singular[c];
	buf_out[(a*nfactors)+b]= sum;
      }
  }
  else {
    /* set them all to zero */
    for (a=0; a<nfactors; a++)
      for (b=0; b<nfactors; b++)
	buf_out[(a*nfactors)+b]= 0.0;
  }

  return 1;
}

Regressor* glm_create_llsq_regressor()
{
  Regressor* result= glm_create_base_regressor();
  LLSqData* data= NULL;
 
  result->feature_array[GLM_TYPE]= GLM_TYPE_LLSQ;
  result->is_settable= llsq_is_settable;
  result->n_params= llsq_n_params;
  result->destroy_self= llsq_destroy_self;
  result->context_valid= llsq_context_valid;
  result->fit= llsq_fit;
  result->project= llsq_project;
  result->normproject= llsq_normproject;
  result->getXtXInv= llsq_getXtXInv;
  if (!(data=(LLSqData*)malloc(sizeof(LLSqData))))
    MALLOC_FAILURE(sizeof(LLSqData),byte);
  result->hook= data;
  data->factors_or_u_ftn= NULL;
  data->v_ftn= NULL;
  data->singular= NULL;
  data->inv_singular= NULL;
  data->residuals= NULL;
  data->covariance_ftn= NULL;
  data->rwork= NULL;
  data->lwork= 0;
  data->nfactorsAllocated= 0;
  data->nobsAllocated= 0;
  data->complexFlagAllocated= 0;
  data->lastFitValid= 0;

  return result;
}

#ifdef never
int main()
{
  Regressor* r= glm_create_llsq_regressor();
  int nfactors= 17;
  glm_set(r,GLM_ORTHO,1);
  printf("%d %d %d %d %d %d %d %d on %d factors -> n_params= %d\n",
	 glm_get(r,GLM_COMPLEX),
	 glm_get(r,GLM_RESIDUALS),
	 glm_get(r,GLM_VARIANCES),
	 glm_get(r,GLM_COVARIANCES),
	 glm_get(r,GLM_SSQR),
	 glm_get(r,GLM_ORTHO),
	 glm_get(r,GLM_DEBUG),
	 glm_get(r,GLM_TYPE),
	 nfactors,
	 glm_n_params(r,nfactors));
  glm_destroy(r);
}
#endif
