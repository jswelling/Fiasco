/************************************************************
 *                                                          *
 *  glm.h                                                   *
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
/* Header file for glm.c */

typedef enum { 
  GLM_COMPLEX, 
  GLM_RESIDUALS, 
  GLM_VARIANCES, 
  GLM_COVARIANCES,
  GLM_SSQR,
  GLM_ORTHO,
  GLM_DEBUG,
  GLM_DEVIANCE,
  GLM_TYPE /* must be last! */
} glm_feature;

#define GLM_N_FEATURES ((int)GLM_TYPE)+1

typedef enum {
  GLM_TYPE_LLSQ,      /* linear least squares */
  GLM_TYPE_LOGISTIC,  /* logistic regression */
  GLM_TYPE_POISSON,    /* Poisson regresson (link function is ln()) */
  GLM_TYPE_BASE,      /* does nothing; must be last! */
} glm_type;

/* The base type doesn't count... */
#define GLM_N_TYPES ((int)GLM_TYPE_BASE)

typedef struct regressor_struct {
  void (*destroy_self)(struct regressor_struct* self);
  int (*context_valid)(struct regressor_struct* self, 
		       int nobs, int nfactors);
  int (*is_settable)(struct regressor_struct* self, glm_feature feature);
  void (*set)(struct regressor_struct* self, 
	      glm_feature feature, int flag);
  int (*get)(struct regressor_struct* self, glm_feature feature);
  int (*n_params)(struct regressor_struct* self, int nfactors);
  int (*fit)(struct regressor_struct* self, 
	     double* obs, const double* factors, const double* counts, 
	     double* param_out, int nobs, int nfactors);
  int (*project)(struct regressor_struct* self, 
		 const double* obs_vec_in, 
		 double* factor_vec_out,
		 int nobs, int nfactors);
  int (*normproject)(struct regressor_struct* self, 
		     const double* obs_vec_in, 
		     double* factor_vec_out, int nobs, int nfactors);
  /* The following is of use to iterative regressors using sub-regressors */
  int (*getXtXInv)(struct regressor_struct* self, double* buf_out,
		   int nobs, int nfactors); /* return (X'X)^-1 */
  /* type field has been replaced by feature_array[GLM_TYPE] */
  int feature_array[GLM_N_FEATURES];
  void* hook;
} Regressor;

Regressor* glm_create_base_regressor(void);
void glm_base_destroy_self( Regressor*  self );
int glm_base_is_settable(Regressor* r, glm_feature feature);
int glm_base_context_valid(Regressor* r, int nobs, int nfactors);

Regressor* glm_create_llsq_regressor(void);
Regressor* glm_create_logistic_regressor(void);
Regressor* glm_create_poisson_regressor(void);

int glm_is_settable(Regressor* r, glm_feature feature);
void glm_set(Regressor* r, glm_feature feature, int flag);
int glm_get(Regressor* r, glm_feature feature);
int glm_n_params(Regressor* r, int nfactors);
void glm_destroy( Regressor* target );

int glm_context_valid(Regressor* r, int nobs, int nfactors);

int glm_fit(Regressor* r, double* obs, const double* factors, 
	    const double* counts, double* param_out, 
	    int nobs, int nfactors);

int glm_project( Regressor* r, const double* obs_vec_in, 
		 double* factor_vec_out, int nobs, int nfactors );

int glm_normproject( Regressor* r, const double* obs_vec_in, 
		     double* factor_vec_out, int nobs, int nfactors );

char* glm_error_msg(void);
void glm_clear_error_msg(void);
char* _glm_err_txt;

const char* glm_get_feature_name(glm_feature feature);
const char* glm_get_regressor_type_name(glm_type type);
Regressor* glm_create_regressor_by_type(glm_type type);

