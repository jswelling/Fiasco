/************************************************************
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *                                                          *
 *     Copyright (c) 1999 Carnegie Mellon University        *
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
 ***********************************************************/

#include "tmperror.h"  /* ATTN: To be replaced by ErrorStream.h.  Only loaded once */

/*** Constants and Macros ***/

static const double BigBound = 1.0e100;
static const double NegBigBound = -1.0e100;

enum   qnb_int_configs      /* Integer Options */
{
    QNB_CONFIG_BOUND_TYPE, QNB_CONFIG_ITER_MAX,  QNB_CONFIG_EVAL_MAX,
    QNB_CONFIG_GRAD_MAX,   QNB_CONFIG_HESS_TYPE, QNB_CONFIG_NUM_CORRECTIONS,
    NUM_QNB_INT_CONFIGS
};
enum   qnb_double_configs   /* Real Options    */
{
    QNB_CONFIG_GRAD_TOL,     QNB_CONFIG_STEP_TOL, QNB_CONFIG_REL_FUNC_TOL,
    QNB_CONFIG_ABS_FUNC_TOL, QNB_FALSE_CONV_TOL,  QNB_CONFIG_MAX_STEP,
    QNB_CONFIG_TRUST_SIZE,
    NUM_QNB_DOUBLE_CONFIGS
};

    /* Error and Warning Messages */

enum   qnb_errors
{
    QNB_ERR_NONE,       QNB_ERR_INFEASIBLE_START, QNB_ERR_INFEASIBLE_BOUNDS,
    QNB_ERR_ITER_LIMIT, QNB_ERR_EVAL_LIMIT,       QNB_ERR_NOT_CONVERGED,
    QNB_ERR_INIT_VALUE_NOT_IMPROVED,
    NUM_QNB_ERRORS
};


/*** Public Function Prototypes ***/

void qnb_minimize( int n, int nbig, double* p, double* plb, double* pub, 
		   double* psc, double fsc, double* work, int lwork, 
		   int* iwork, int liwork,
		   double* val, int* piters, int* pevals, double* grad_tol, 
		   void (*func)(int*, double*, double*), 
		   void (*dfunc)(int*, double*, double*),
		   ErrorStream es );
void            qnb_initialize( int n, int nbig );
void            qnb_cleanup();

int             qnb_config_int( const char*, int );
int             qnb_config_double( const char*, double );


void   minimize_nmb( int n, int nbig, double *p, double *val, int *itmax, 
		     double tol,
                     double *plb, double *pub, double *psc, double fsc, 
		     double *work, int lwork,
                     void (*func)(int *, double *, double *), 
		     void (*dfunc)(int *, double *, double *),
                     int *info );

void   min_config_nmb_int( const char *, int val );
void   min_config_nmb_double( const char *, double val );
void   initialize_min_nmb( int n, int nbig );
void   cleanup_min_nmb( int n, int nbig );
char  *min_errors_nmb( int );

void     find_drift_df( int D, double *P, double df, double *pu, 
			double *pulow, double *puhigh,
			double *eig, double *work, int lwork, double tol );
