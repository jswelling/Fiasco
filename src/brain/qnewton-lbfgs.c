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

static char rcsid[] = "$Id: qnewton-lbfgs.c,v 1.6 2007/03/21 23:45:49 welling Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include  "stdcrg.h"
#include  <math.h>

#include  "tmperror.h"  /* ATTN: TEMPORARY. To be replaced by ErrorStream.h in version 1*/
#include  "minqnb.h"

/*** Constants and Macros ***/

    /* L-BFGS-B Interface */

enum lbfgs_indices
{
    TOTAL_ITERS=29, TOTAL_EVALS=33,    FREE_VARIABLES=37,    ACTIVE_CONSTRAINTS=38,
    LAST_VALUE=1,   MACH_PRECISION=4,  GRAD_INF_NORM=12
};

    /* Error and Warning Messages */

static const char* qnb_error_msgs[] =
                   {
                       "",
                       "Starting value for optimization infeasible",
                       "Variable bounds for optimization infeasible",
                       "Iteration limit exceeded for optimization",
                       "Evaluation limit exceeded for optimization",
                       "Optimization did not converge, best approximation given",
                       "Optimization did not improve starting value"
                   };

/*** Private Function Prototypes ***/

static void     set_default_params( int*, double* );


    /* L-BFGS-B Routines */

void            setulb( int*, int*, double*, double*, double*, int*, double*, double*, double*, double*,
                        double*, int*, char*, int*, char*, long*, long*, double* );


/*** Static and Global Variables ***/

    /* Algorithm Parameters */

static int        Algorithm_Initialized = 0;
static int        Default_Params_Set = 0;
static int        IParam[NUM_QNB_INT_CONFIGS];
static double     RParam[NUM_QNB_DOUBLE_CONFIGS];
static int        opt_counter; /* convenient for debugging */


    /* Once Allocated Working Arrays */

static int*       nbd;
static double*    grad;


/***  Public Source  ***/

/*
 * qnb_minimize() applies the lbfgs code to minimize the function func
 * with derivative function dfunc.  Only the following arguments are
 * referenced:
 *
 *           ...
 * 
 * The rest are included for compatibility with other versions of qnb_minimize
 *
 * All errors encountered are written to the given error stream.
 *
 */

void            qnb_minimize( int n, int nbig, 
			      double* p, double* plb, double* pub, 
                              double* psc, double fsc, 
			      double* work, int lwork, int* iwork, int liwork,
                              double* val, int* piters, 
			      int* pevals, double* grad_tol, 
                              void (*func)(int*, double*, double*), 
			      void (*dfunc)(int*, double*, double*),
                              ErrorStream es )
{

    /* Algorithm Parameters */

    int          i;      
    int     iprint = -1; /* set -1 for normal operation */
    int     corr;
    int     itmax;
    int     evmax;

    char    task[61];
    char    csave[61];

    long    lsave[4];
    long    isave[44];

    double  factr;
    double  pgtol;
    
    double  dsave[29];
    
        
    /* Other local variables */

    int     its;

    /* Initialization --- Keep this minimal since this will be done every time */

    if( !Algorithm_Initialized )
        qnb_initialize( n, nbig );

    strcpy( task, "START" );

        /* ATTN: These vars can be eliminated and pointers into Param vecs supplied */
        /* I put them in temporarily here for clarity but we can remove them later  */

    corr  = IParam[QNB_CONFIG_NUM_CORRECTIONS];    
    itmax = IParam[QNB_CONFIG_ITER_MAX];          
    evmax = IParam[QNB_CONFIG_EVAL_MAX];

    factr = RParam[QNB_CONFIG_REL_FUNC_TOL];
    pgtol = RParam[QNB_CONFIG_GRAD_TOL];

    /* Check work area sizes */
    if (lwork < (2*corr*n +  4*n + 11*corr*corr + 8*corr)) {
      es_write_mesg( es, ES_FATAL, 
		     "Minimizer double work area too small: %d vs. %d",
		     lwork, (2*corr*n +  4*n + 11*corr*corr + 8*corr) );
    }
    if (liwork < (3*n)) {
      es_write_mesg( es, ES_FATAL, 
		     "Minimizer double work area too small: %d vs. %d",
		     liwork, (3*n) );
    }

    /* Initialize bounds flag array, allocated in qnb_initialize */
    for( i = 0; i < nbig; i++ )
      {
	if( plb[i] > NegBigBound )
	  {
	    if( pub[i] < BigBound )
	      nbd[i] = 2;          /* Both upper and lower bounds */
	    else
	      nbd[i] = 1;          /* Only lower bound            */
	  }
	else if( pub[i] < BigBound )
	  nbd[i] = 3;              /* Only upper bound            */
	else
	  nbd[i] = 0;              /* No bounds on variable       */
      }

    /* Main Loop */

    its = 0;

    for( ;; )
    {
        setulb( &n, &corr, p, plb, pub, nbd, val, grad, &factr, &pgtol,
                work, iwork, task, &iprint, csave, lsave, isave, dsave );

        if( task[0] == 'F' && task[1] == 'G' )  
        {
            /* Need new evaluation of function and gradient */

            func( &n, p, val );
            dfunc( &n, p, grad );

            if( isave[TOTAL_EVALS] + 2 > evmax )
            {
	        strcpy( task, "STOP" );

                es_write_error( es, ES_WARN_D, QNB_ERR_EVAL_LIMIT, qnb_error_msgs[QNB_ERR_EVAL_LIMIT] );
            }
        }
        else if( task[0] == 'N' && task[1] == 'E' && task[2] == 'W' && task[3] == '_' && task[4] == 'X' )
        {
            /* One iteration completed...decide whether or not to continue */

            if( ++its > itmax )
            {
	        strcpy( task, "STOP" );

                es_write_error( es, ES_WARN_D, QNB_ERR_ITER_LIMIT, qnb_error_msgs[QNB_ERR_ITER_LIMIT] );
            }
            
        }
        else if( (task[0] == 'S' && task[1] == 'T' && task[2] == 'O' && task[3] == 'P') ||
                 (task[0] == 'C' && task[1] == 'O' && task[2] == 'N' && task[3] == 'V') ) {
            break;  /* Done */
	}
        else
        {
            if( task[0] == 'A' )
                es_write_error( es, ES_WARN_D, QNB_ERR_NOT_CONVERGED,
                                "%s (mesg: %s) on optimization %d",
                                qnb_error_msgs[QNB_ERR_NOT_CONVERGED], task,
				opt_counter);
            else
                es_write_error( es, ES_WARN_D, QNB_ERR_INFEASIBLE_START,
                                "%s (mesg: %s) on optimization %d",
                                qnb_error_msgs[QNB_ERR_INFEASIBLE_START], 
				task, opt_counter );
            break;

            /* Insert other stopping rules here if desired */

        }
    }

    *piters = its;
    *pevals = isave[TOTAL_EVALS];
    *grad_tol = dsave[GRAD_INF_NORM];

    opt_counter++;
}

void   qnb_initialize( int n, int nbig )
{
    if( !Algorithm_Initialized )
    {
        if( !Default_Params_Set )
            set_default_params( IParam, RParam );

        grad = Calloc( nbig, double );
        nbd = Malloc( nbig, int );
    }

    Algorithm_Initialized = 1;
    opt_counter= 0;
}

void   qnb_cleanup()
{
    Algorithm_Initialized = 0;

    Free( nbd );
    Free( grad );
}

/*
 * qnb_config_int() and qnb_config_double() set the
 * integer- and real-valued algorithm parameters for
 * the qnb.  They take a string and set the corresponding
 * option to the given value.
 *
 * Both return 1 if the given string matches an option name
 * of the appropriate type and 0 otherwise.
 *
 */

int    qnb_config_int( const char* str, int val )
{
    static const char   *QNB_Int_Options[] =
              {
                  "bounds",       "max iterations", "max evals", "max grad evals",
                  "hessian type", "corrections"
              };

    int                  i;

    if( !Default_Params_Set )
	set_default_params( IParam, RParam );

    for( i = 0; i < NUM_QNB_INT_CONFIGS; i++ )
        if( !strcasecmp( str, QNB_Int_Options[i] ) )
        {
            IParam[i] = val;

            return( 1 );
        }

    return( 0 );
}


int    qnb_config_double( const char* str, double val )
{
    static const char   *QNB_Real_Options[] =
        {
            "gradient tolerance", "step tolerance", "relative function tolerance",
            "absolute function tolerance", "false convergence tolerance",
            "max step size", "trust region size" 
        };

    int            i;

    if( !Default_Params_Set )
	set_default_params( IParam, RParam );

    for( i = 0; i < NUM_QNB_DOUBLE_CONFIGS; i++ )
        if( !strcasecmp( str, QNB_Real_Options[i] ) )
        {
            RParam[i] = val;

            return( 1 );
        }

    return( 0 );
}


   /***  Private Source  ***/

static
void            set_default_params( int* iparams, double* rparams )
{
    /* Set conservative defaults of those parameters used by lbfgs */

    IParam[QNB_CONFIG_ITER_MAX]        = 10000;
    IParam[QNB_CONFIG_EVAL_MAX]        = 100000;
    IParam[QNB_CONFIG_NUM_CORRECTIONS] = 5;
    
    RParam[QNB_CONFIG_GRAD_TOL]        = 1.0e-4;
    RParam[QNB_CONFIG_REL_FUNC_TOL]    = 1.0e8;

    Default_Params_Set = 1;
}



/*
 * Local Variables:        
 * mode: c
 * eval: (make-local-variable 'compile-command)
 * compile-command: "make qnewton-lbfgs.o"
 * End:
 *
 */
