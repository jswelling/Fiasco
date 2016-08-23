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

static char rcsid[] = "$Id: amoeba.c,v 1.5 2007/03/21 23:45:49 welling Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include  "stdcrg.h"
#include  <math.h>


/*** Constants and Macros ***/

static double const Alpha        = 1.0;  /* Reflection               */
static double const Beta         = 0.5;  /* One-variable contraction */
static double const Gamma        = 2.0;  /* Reflection and Expansion */

#define set_error(i)           (errc=(i))
#define set_warning(i)         (warnc=(i))


/*** Function Prototypes ***/

    /* Public  */

void   minimize_nmb( int n, int nbig, double *p, double *val, int *itmax, double tol,
		     double *plb, double *pub, double *psc, double fsc, double *work, int lwork,
		     void (*func)(int *, double *, double *), void (*dfunc)(int *, double *, double *),
		     int *info );

void   min_config_nmb_int( const char *, int val );
void   min_config_nmb_double( const char *, double val );
void   initialize_min_nmb( int n, int nbig );
void   cleanup_min_nmb( int n, int nbig );
char  *min_errors_nmb( int );

    /* Private */

static void   init_simplex( int, int , double *, double *, double *, double *, double,
			    void (*)(int *,double *,double *) );

static void   amoeba( double **, double y[], double *, double *, int ,
		      void (*)(int *, double *, double *), int *, int * );

static double amotry( double **, double *, double *, double *, double *, int,
		      void (*func)(int *, double *, double *), int ihi, int *nfunc, double fac);

        /* Utility Routines */

static int    *ivector( int, int );
static double *dvector( int, int );
static double **dmatrix( int, int, int, int );
static void free_dmatrix( double **m, int nrl, int nrh, int ncl, int nch );
static void free_dvector( double *v, int nl, int nh );

/*** Static and Global Variables ***/

    /* Private Values */

static int      init_alloc = 0;

static double   Tolerance = 1.0e-9;

static double   *simpvals;
static double   *ptry;
static double   *psum;

static double  **simplex;

    /* Configurable Parameters */

enum   configP { ITER_MAX, RESTARTS, MAX_TRIES, AUTO_CLEAN, NUM_CONFIG_P };
enum   configQ { LARGE_VALUE, INIT_CHANGE, SEARCH_UP, SEARCH_DOWN, SCALE, BASE, THRESH, NUM_CONFIG_Q };

static int    AmoebaP[] = { 5000, 1, 25, 1 };
static double AmoebaQ[] = { 1.0e256, 0.05, 1.618034, 0.618034, 0.5, 1.0, 2.0 };

    /* Error and Warning Messages */

enum     err_codes { NO_ERR, INFEASIBLE_BOUNDS, ITERATION_LIMIT, NUM_ERR_CODES };
enum     wrn_codes { NO_WRN, UNMATCHED_CONFIG, INFEASIBLE_START, INIT_VALUE_NOT_IMPROVED, NUM_WRN_CODES };

static   int       errc = 0;
static   int       warnc = 0;

static   char     *errors[] = {
                                "no error",
                                "Error in minimize [method=amoeba]:  Given variable bounds infeasible",
				"Error in minimize [method=amoeba]:  Iteration limit exceeded"
			      };

static   char     *warnings[] = {
                                  "no warning",
                                  "Warning in min_config: Given config tag not matched",
				  "Warning in minimize_nmb: Starting point infeasible",
				  "Warning in minimize_nmb: Initial objective value varies slowly in at least one direction"
			        };

   /***  Public Source  ***/

void     minimize_nmb( int n, int nbig, double *p, double *val, int *itmax, double tol,
		       double *plb, double *pub, double *psc, double fsc, double *work, int lwork,
		       void (*func)(int *, double *, double *),
		       void (*dfunc)(int *, double *, double *), int *info )
{
    int      i, ilo, iter = 0, iter_save = 0, restarts, nfunc;

    /* Since the NRIC functions are 1-indexed, ptrs needs to be adjusted */
    /* to reflect this.  Note that we are assuming that func and dfunc   */
    /* take a zero-indexed vector, so pointers that are passed to func   */
    /* do need to be re-adjusted back.                                   */

    p--;
    plb--;
    pub--;
    psc--;

    /* Basic Settings */

    errc = 0;             /* Clear errors and warnings */
    warnc = 0;

    Tolerance = tol;

    if( *itmax > 0 )
	AmoebaP[ITER_MAX] = *itmax;

    if( !init_alloc )
	initialize_min_nmb( n, nbig );

    /* Main Processing */

    for( restarts = 0; restarts <= AmoebaP[RESTARTS]; restarts++ )
    {
	init_simplex( n, nbig, p, plb, pub, psc, fsc, func );

	if( errc )
	{
	    *info = errc;
	    return;
	}
	else if( warnc )
	    *info = -warnc;

	amoeba( simplex, simpvals, plb, pub, n, func, &nfunc, &iter );

	if( errc )
	{
	    *info = errc;
	    return;
	}
	else if( warnc )
	    *info = -warnc;

	for( ilo = 1, i = 2; i <= n + 1; i++ )
	    if( simpvals[i] < simpvals[ilo] )
		ilo = i;

	for( i = 1; i <= nbig; i++ )
	    p[i] = simplex[ilo][i];

	*val = simpvals[ilo];

	iter_save += iter;
    }
    
    *itmax = iter_save;    /* Save Iteration Count  */
    work[0] = nfunc;       /* Save # Function Calls */

    if( AmoebaP[AUTO_CLEAN] )
	cleanup_min_nmb( n, nbig );
}

void   initialize_min_nmb( int n, int nbig )
{
    if( !init_alloc )
    {
	simplex = dmatrix( 1, n+1, 1, nbig );

	simpvals = dvector( 1, n + 1 );
	psum = dvector( 1, nbig );
	ptry = dvector( 1, nbig );

	init_alloc = 1;
    }
}

void   cleanup_min_nmb( int n, int nbig )
{
    if( init_alloc )
    {
	free_dvector( simpvals, 1, n + 1 );
	free_dvector( psum, 1, nbig );
	free_dvector( ptry, 1, nbig );

	free_dmatrix( simplex, 1, n + 1, 1, nbig );

	init_alloc = 0;
    }
}

void   min_config_nmb_int( const char *str, int val )
{
    static const char   *StringsP[] = { "iterations", "restarts", "max tries", "auto clean" };

    int                  i;

    for( i = 0; i < NUM_CONFIG_P; i++ )
        if( !strcasecmp( str, StringsP[i] ) )
        {
            AmoebaP[i] = val;
            break;
        }

    if( i == NUM_CONFIG_P )
        set_warning( UNMATCHED_CONFIG );
}

void   min_config_nmb_double( const char *str, double val )
{
    static const char   *StringsQ[] =
        {
            "large value", "init change", "search up", "search down",
            "default scale", "base", "threshold"
        };

    int                  i;

    for( i = 0; i < NUM_CONFIG_Q; i++ )
        if( !strcasecmp( str, StringsQ[i] ) )
        {
            AmoebaQ[i] = val;
            break;
        }

    if( i == NUM_CONFIG_Q )
        set_warning( UNMATCHED_CONFIG );
}

char  *min_errors_nmb( int code )
{
    if( !code || (code > NUM_ERR_CODES) || (code < -NUM_WRN_CODES) )
	return( NULL );

    if( code > 0 )   /* Error */
	return( errors[code] );

    return( warnings[-code] );
}


   /***  Private Source  ***/

static void   init_simplex( int n, int nbig, double *pinit, double *plb, double *pub, double *psc,
			    double fsc, void (*func)(int *, double *, double *) )
{
    int        i, j;
    double     u, w, x, z;
    double     v, vlower, vupper;

    /* Since this NRIC code is 1-indexed, pinit needs to be adjusted  */
    /* but we are assuming that this is called by minimize which will */
    /* have already done the adjustment.  Nonetheless, pointers that  */
    /* are passed to (*func) do need to be re-adjusted back.          */
    
    for( i = 1; i <= nbig; i++ )
    {
	if( pub[i] <= plb[i] )  /* If a variable is desired fixed, it can be given index n < i <= nbig */
	{
	    set_error( INFEASIBLE_BOUNDS );
	    return;
	}

	if( (pinit[i] > pub[i]) || (pinit[i] < plb[i]) )
	    set_warning( INFEASIBLE_START );

	pinit[i] = Min( Max( plb[i], pinit[i] ), pub[i] );   /* Make starting point feasible */
	ptry[i] = pinit[i];
    }

    for( i = 1; i <= n + 1; i++ )
        for( j = 1; j <= nbig; j++ )
            simplex[i][j] = pinit[j];

    (*func)( &n, simplex[1] + 1, &simpvals[1] );
    
    if( fsc == 0.0 )
	x = AmoebaQ[INIT_CHANGE] * fabs( simpvals[1] );
    else
	x = fabs( fsc );
	    
    for( j = 1; j <= n; j++ )
    {
	if( psc[j] == 0.0 )
	    v = (AmoebaQ[SCALE] * fabs(pinit[j]) + AmoebaQ[BASE]);
	else
	    v = fabs( psc[j] );

	if( pinit[j] == pub[j] )   /* Against upper bound */
	    v = -Min( v, 0.5 * (pub[j] - plb[j]) );
	else
	    v = Min( v, 0.5 * (pub[j] - pinit[j]) );
	
	/* Find decreasing direction, feasible first jump */

	ptry[j] = pinit[j] + v;      /* This is designed to be feasible */

	(*func)( &n, ptry + 1, &u );

	if( (u > simpvals[1]) && (pinit[j] - v < pub[j]) && (pinit[j] - v > plb[j]) )
	{
	    ptry[j] = pinit[j] - v;
	    (*func)( &n, ptry + 1, &w );

	    if( w < simpvals[1] )
	    {
		u = w;
		v *= -1;
	    }
	}
	
	z = fabs( u - simpvals[1] );

	if( z < x || ( (z >= AmoebaQ[THRESH] * x) && (u >= simpvals[1]) ) )
	{
	    if( z >= x )
	    {
		vupper = v;
		vlower = 0.0;
		v *= AmoebaQ[SEARCH_DOWN];
	    }
	    else
	    {
		vlower = v;
		vupper = (v > 0.0) ? (pub[j] - pinit[j]) : (plb[j] - pinit[j]);
			
		v = Minabs( vupper, v * AmoebaQ[SEARCH_UP] );
	    }

	    /* Improve jump to within threshold */

	    for( i = 0; i < AmoebaP[MAX_TRIES]; i++ )
	    {
		ptry[j] = pinit[j] + v;

		if( (ptry[j] > pub[j]) || (ptry[j] < plb[j]) )  /* Outside bounds */
		    u = AmoebaQ[LARGE_VALUE];
		else
		    (*func)( &n, ptry + 1, &u );
	    
		z = fabs( u - simpvals[1] );

		if( z >= x && ( (z < AmoebaQ[THRESH] * x) || (u < simpvals[1]) ) )
		    break;
		else if( z >= x )
		{
		    vupper = v;
		    v = vlower + AmoebaQ[SEARCH_DOWN] * (vupper - vlower);
		}
		else
		{
		    vlower = v;
		    v = Minabs( v * AmoebaQ[SEARCH_UP], vlower + AmoebaQ[SEARCH_DOWN]*(vupper - vlower) );
		}
	    }

	    if( i == AmoebaP[MAX_TRIES] )
		set_warning( INIT_VALUE_NOT_IMPROVED );
	}

	if( u < simpvals[1] )  /* Take one more step, trying to improve things further */
	{
	    vlower = v;
	    vupper = (v > 0.0) ? (pub[j] - pinit[j]) : (plb[j] - pinit[j]);
			
	    v = Minabs( vupper, v * AmoebaQ[SEARCH_UP] );

	    ptry[j] = pinit[j] + v;
	    
	    if( (ptry[j] >= pub[j]) || (ptry[j] <= plb[j]) )  /* At or Outside bounds */
		w = AmoebaQ[LARGE_VALUE];
	    else
		(*func)( &n, ptry + 1, &w );

	    if( w <= u )
		u = w;
	    else
		ptry[j] = pinit[j] + vlower;
	}
	
	
        simplex[j+1][j] = ptry[j];
	simpvals[j+1] = u;

	ptry[j] = pinit[j];
    }

    for( i = 1; i <= nbig; i++ )
        psum[i] = ptry[i] = pinit[i];
}

static void amoeba( double **p, double y[], double *plb, double *pub, int ndim,
		    void (*func)(int *, double *, double *), int *nfunc, int *piter )
{
    int     i, j;
    int     ilo, ihi, inhi;
    int     mpts = ndim + 1;
    int     iter;
    int     outside_bounds;

    double  ytry, ysave, sum, rtol;
    
    for( j = 1; j <= ndim; j++ )
    {
	for( i = 1, sum = 0.0; i <= mpts; i++ )
	    sum += p[i][j];

	psum[j]=sum;
    }    

    for( iter = 0; ; iter++ ) 
    {
        ilo=1;
        ihi = y[1]>y[2] ? (inhi=2,1) : (inhi=1,2);

        for (i=1;i<=mpts;i++) 
        {
            if (y[i] < y[ilo])
		ilo=i;

            if (y[i] > y[ihi])
	    {
                inhi = ihi;
                ihi = i;
            }
	    else if( (y[i] > y[inhi]) && (i != ihi) )
                inhi=i;
        }

        rtol = 2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo]));

#       if DIAGNOSTIC > 1 
        Message( "Amoeba Search:  Lo = %16.8lg,  Hi = %16.8lg,  Tol = %16.8lg\n", y[ilo], y[ihi], rtol );
#       endif        

        if( rtol < Tolerance )
	{
	    *piter = iter;
            break;
	}

        if( iter >= AmoebaP[ITER_MAX] )
	{
            set_error( ITERATION_LIMIT );
	    *piter = iter;
	    return;
	}
        
        ytry = amotry( p, y, psum, plb, pub, ndim, func, ihi, nfunc, -Alpha );
        
        if( ytry <= y[ilo] )
            ytry = amotry( p, y, psum, plb, pub, ndim, func, ihi, nfunc, Gamma );
        else if( ytry >= y[inhi] )
        {
            ysave = y[ihi];
            ytry = amotry( p, y, psum, plb, pub, ndim, func, ihi, nfunc, Beta );

            if (ytry >= ysave)
            {
                for( i = 1; i <= mpts; i++ )
                {
                    if( i != ilo )
                    {
                        for( j = 1, outside_bounds = 0; j <= ndim; j++ ) 
                        {
                            psum[j]=0.5*(p[i][j]+p[ilo][j]);
                            p[i][j]=psum[j];
			    outside_bounds += ((psum[j] > pub[j]) || (psum[j] < plb[j]));
                        }

			if( outside_bounds )
			    y[i] = AmoebaQ[LARGE_VALUE];
			else
			{
			    (*func)( &ndim, psum + 1, &y[i] );
			    (*nfunc)++;
			}
                    }
                }

		for( j = 1; j <= ndim; j++ )
		{
		    for( i = 1, sum = 0.0; i <= mpts; i++ )
			sum += p[i][j];

		    psum[j] = sum;
		}    
            }
        }
    }
}

static double amotry( double **p, double *y, double *psum, double *plb, double *pub, int ndim,
		      void (*func)(int *, double *, double *), int ihi, int *nfunc, double fac )
{
    int        j;
    int        outside_bounds = 0;
    double     fac1, fac2, ytry;

    fac1 = (1.0 - fac)/ndim;
    fac2 = fac1 - fac;

    for( j = 1, outside_bounds = 0; j <= ndim; j++ )
    {
	ptry[j] = psum[j]*fac1 - p[ihi][j]*fac2;
	outside_bounds += ((ptry[j] > pub[j]) || (ptry[j] < plb[j]));
    }

    if( outside_bounds )
	ytry = AmoebaQ[LARGE_VALUE];
    else
    {
	(*func)( &ndim, ptry + 1, &ytry );
	++(*nfunc);
    }

    if( ytry < y[ihi] ) 
    {
	y[ihi] = ytry;
	for( j = 1; j <= ndim; j++ )
	{
	    psum[j] += ptry[j] - p[ihi][j];
	    p[ihi][j] = ptry[j];
	}
    }

    return( ytry );
}

/*** NRIC utility routines     ***/

static int    *ivector( int nl, int nh)
{
    int *v;

    v = Malloc( (unsigned)(nh-nl+1), int );

    return( v - nl );
}

static double *dvector( int nl, int nh)
{
    double *v;

    v = Malloc( (unsigned)(nh-nl+1), double );

    return( v - nl );
}

static double **dmatrix( int nrl, int nrh, int ncl, int nch )
{
    int i;
    double **m;

    m = Malloc( (unsigned)(nrh-nrl+1), double * );
    m -= nrl;

    for( i = nrl; i <= nrh; i++ )
    {
	m[i] = Malloc( (unsigned)(nch-ncl+1), double );
	m[i] -= ncl;
    }

    return( m );
}


static void free_dvector( double *v, int nl, int nh )
{
    Free( (void *)(v+nl) );
}


static void free_dmatrix( double **m, int nrl, int nrh, int ncl, int nch )
{
    int i;

    for( i = nrh; i >= nrl; i-- )
	Free( (void *)(m[i]+ncl) );

    Free( (void *)(m+nrl) );
}


/*
 * Local Variables:        
 * mode: c
 * eval: (make-local-variable 'compile-command)
 * compile-command: "make amoeba.o"
 * End:
 *
 */

