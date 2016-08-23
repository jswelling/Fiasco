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

static char rcsid[] = "$Id: info.c,v 1.4 2000/10/06 05:32:04 welling Exp $";

/*********************************************************************
 *  info.c                                                           *
 *                                                                   *
 *  The routine hessian_scd computes a second-central difference     *
 *  approximation of negative Hessian of a function at the given     *
 *  point.  If the function is a log-likelihood, this computes       *
 *  the observed Fisher information.                                 *
 *                                                                   *
 *  The routine hessian_fcd computes a first-central difference      *
 *  approximation of negative Hessian of a function at the given     *
 *  point, using a user-supplied gradient.  If the function is a     *
 *  log-likelihood, this computes the observed Fisher information.   *
 *                                                                   *
 *  Two approaches are provided for both cases:                      *
 *                                                                   *
 *  (1) Direct Evaluation at  epsilon^p x_c                          *
 *  (2) Richardson Extrapolation                                     *
 *                                                                   *
 *  Central difference are used throughout, unless a parameter is    *
 *  at its boundary in which case directional differences are used.  *
 *                                                                   *
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "stdcrg.h"
#include "lapack.h"
#include <math.h>

/* Constants and Macros */

#define   REQUIRED_SPACE(n)      (7*(n))

enum methods  { DIRECT, RICHARDSON, NUM_METHODS };
enum errors   { NOERR, ERROR_WORK_SPACE, ERROR_BAD_METHOD, ERROR_BAD_H, NUM_ERRORS };
enum warnings { NOWRN, WARNING_NO_RICH_INIT, WARNING_SMALL_DENOM, WARNING_SMALL_NUMER,
		WARNING_ERR_EXCEEDS_TOL, WARNING_CRAZY_VALUE, NUM_WARNS=16 };

static double const Scale = 1.0e-1;
static double const Under = 1.0e-34;

static double const Rich_Err_Scale  = 16.0;

#define  set_error(i)           (err_code = (i))
#define  set_warning(i,k)       (warn_code = -((i) + (k)*NUM_WARNS))


/* Function Prototypes  */

void         hessian_scd( int n, int nbig, double *p, double *plb, double *pub, double *I, int method,
			  double *hinit, int *skip, double maxmin, double errf, double tol,
			  double *work, int worksize, void (*func)(int *, double *, double *), int *info );

void         hessian_fcd( int n, int nbig, double *p, double *plb, double *pub, double *I, int method,
			  double *hinit, int *skip, double maxmin, double errf, double tol,
			  double *work, int worksize, void (*func)(int *, double *, double *), int *info );

static void  dummy_func( double * );

static void  scd( int n, int nbig, double *x, double h, double *dx, double *v, double *err,
		  double *w, void (*func)(int *, double *, double *), double maxval );
static void  fcd( int n, int nbig, double *x, double h, double *dx, double *v, double *err,
		  double *w, void (*func)(int *, double *, double *), double maxval );

static void  rich_extrap_scd( int n, int nbig, double *x, double h, double *dx, double *v, double *err,
			      double *w, void (*func)(int *, double *, double *), double maxval );

static void  rich_extrap_fcd( int n, int nbig, double *x, double h, double *dx, double *v, double *err,
			      double *w, void (*func)(int *, double *, double *), double maxval );

static void  sfd( int n, int nbig, double *x, double h, double *dx, double *v, double *err,
		  double *w, void (*func)(int *, double *, double *), double maxval );
static void  ffd( int n, int nbig, double *x, double h, double *dx, double *v, double *err,
		  double *w, void (*func)(int *, double *, double *), double maxval );

static void  rich_extrap_sfd( int n, int nbig, double *x, double h, double *dx, double *v, double *err,
			      double *w, void (*func)(int *, double *, double *), double maxval );

static void  rich_extrap_ffd( int n, int nbig, double *x, double h, double *dx, double *v, double *err,
			      double *w, void (*func)(int *, double *, double *), double maxval );

double      dlamch( char * );

/* Global Variables     */

static  int err_code  = 0;       /* Error Code   */
static  int warn_code = 0;       /* Warning Code */

static  void  (*SD[2*NUM_METHODS])( int, int, double *, double, double *, double *, double *,
				    double *, void (*)(int *, double *, double *), double )
                 = { scd, sfd, rich_extrap_scd, rich_extrap_sfd };




/*
 * HESSIAN_SCD computes an approximate negative hessian by second-central
 * differences of the given function at the given point.  The algorithm uses
 * one of two methods:  (1) direct evaluation and (2) Richardson extrapolation,
 * depending on the method argument. 
 *
 */

void     hessian_scd( int n, int nbig, double *p, double *plb, double *pub, double *I, int method,
		      double *hinit, int *skip, double maxmin, double errf, double tol,
		      double *work, int worksize, void (*func)(int *, double *, double *), int *info )
{
    char       ech = 'E';
    
    int        i, j, ns;
    int        indi, indj;
    int        which;

    double     temp, maxval;
    double     h, deriv = 0.0;
    double     err, errf_for, eps;
    double     ui, uj;

    double     *x, *dx, *se;
    double     *cap, *dir;


    /* Clear Error Indicators */

    set_error( 0 );
    set_warning( 0, 0 );
    *info = 0;

    /* Check that we have enough space */

    if( worksize < REQUIRED_SPACE(nbig) )
    {
	*info = set_error( ERROR_WORK_SPACE );
	return;
    }

    /* Set-Up */
    
    x = work;
    ns = n;
    work += nbig;

    dx = work;

    for( i = 0; i < nbig; i++ )
	dx[i] = 0.0;

    work += nbig;

    se = work;
    work += nbig;

    dir = work;
    work += nbig;

    cap = work;
    work += nbig;
    

    if( method == DIRECT )
    {
	if( errf == 0.0 )       
	{
	    eps = dlamch( &ech );
	    errf = pow( eps, 0.25 );
	    errf_for = pow( eps, 0.33 );
	}
	else
	    errf_for = pow( errf, 1.33 );
    }
    
    
    /*
     * Convert step size so that p+h - p is an exactly
     * machine representable number.  This is the role
     * of the odd assignment to temp; the call to
     * dummy_func() is intended to fool optimizers that
     * remove this code.
     *
     */

    for( i = 0; i < n; i++ )
    {
	x[i] = p[i];

	if( skip[i] )
	{
	    ns--;
	    continue;
	}
	else
	{
	    cap[i] = Min( p[i] - plb[i], pub[i] - p[i] );
	    dir[i] = (p[i] - plb[i] <= pub[i] - p[i]) ? 1.0 : -1.0;

	    if( method == DIRECT )
	    {
		if( hinit[i] == 0.0 )
		    hinit[i] = dMax( 1.0, fabs(p[i]) );

		if( errf * hinit[i] > cap[i] )
		    temp = p[i] + dir[i] * errf_for * hinit[i];
		else
		    temp = p[i] + errf * hinit[i];
	    }
	    else if( method == RICHARDSON )
	    {
		if( hinit[i] == 0.0 )
		{
		    hinit[i] = Scale;       /* Not recommended */
		    set_warning( WARNING_NO_RICH_INIT, i*(n+1) );
		}

		if( fabs(hinit[i]) > cap[i] )
		    temp = p[i] + dir[i] * fabs(hinit[i]);
		else
		    temp = p[i] + hinit[i];
	    }
	    else
	    {
		*info = set_error( ERROR_BAD_METHOD );
		return;
	    }
	}
	

	dummy_func( &temp );          /* fool the optimizing compiler */
	hinit[i] = temp - p[i];
    }

    for( ; i < nbig; i++ )
    {
	x[i] = p[i];
	hinit[i] = 0.0;
    }

    /* Record extreme value */

    (*func)( &n, p, &maxval );


    /* Start with unmixed second derivatives */

    for( i = 0, indi = 0; i < n; i++ )
    {
	if( skip[i] )
	    continue;
	
	h = hinit[i];
	dx[i] = 1.0;
	which = (fabs(h) > cap[i]) + method*NUM_METHODS;   /* Which difference function to use */

	(*SD[which])( n, nbig, x, h, dx, &deriv, &err, work, func, maxval );
	deriv *= maxmin;

	if( err > tol )
	{
	    h /= Rich_Err_Scale;
	    which = (fabs(h) > cap[i]) + method*NUM_METHODS;  /* Check again */
		
	    (*SD[which])( n, nbig, x, h, dx, &deriv, &err, work, func, maxval );
	    deriv *= maxmin;

	    if( err > tol )
		set_warning( WARNING_ERR_EXCEEDS_TOL, i*(n+1) );
	}

 	if( deriv <= 0.0 )     /* Completely inaccurate numerical derivative --- skip it */
	{
	    set_warning( WARNING_CRAZY_VALUE, i*(n+1) );
	    se[i] = 0.0;
	    skip[i] = NUM_WARNS;   /* This value must be anything > 1 in absolute value */
	}
	else                  /* Value at least reasonable -- store it as usual         */
	{
	    I[ indi + (indi * (indi + 1))/2 ] = deriv;
	    indi++;
	    se[i] = sqrt( deriv );
	}

	x[i] = p[i];                       /* Restore copy of parameter vector */
	dx[i] = 0.0;
    }

    /* Next mixed second derivatives using unmixed results */

    for( i = 0, indi = 0; i < n; i++ )
    {
	if( skip[i] )
	    continue;

	ui = 1.0/se[i];
	dx[i] = (1 - 2*(hinit[i] < 0.0)) * ui;

	for( j = 0, indj = 0; j < i; j++ )
	{
	    if( skip[j] )
		continue;
	   
	    h = fabs(hinit[i]) + fabs(hinit[j]);

	    uj = 1.0/se[j];
	    dx[j] = (1 - 2*(hinit[i] < 0.0)) * uj;

	    /* Which difference function to use? */

            which = ((h*ui > cap[i]) | (h*uj > cap[j])) + method*NUM_METHODS;

	    (*SD[which])( n, nbig, x, h, dx, &deriv, &err, work, func, maxval );
	    deriv *= maxmin;

	    if( err > tol )
	    {
		h /= Rich_Err_Scale;
		
		which = ((h*ui > cap[i]) | (h*uj > cap[j])) + method*NUM_METHODS;   /* Check again */

		(*SD[which])( n, nbig, x, h, dx, &deriv, &err, work, func, maxval );
		deriv *= maxmin;

		if( err > tol )
		    set_warning( WARNING_ERR_EXCEEDS_TOL, i + j*n );
	    }

	    I[indj + (indi*(indi+1))/2] = 0.5 * (deriv - 2.0) * se[i] * se[j];

	    indj++;

	    x[j] = p[j];
	    dx[j] = 0.0;
	}

	indi++;
	
	x[i] = p[i];
	dx[i] = 0.0;
    }

    *info = (err_code) ? err_code : warn_code;
}


void     hessian_fcd( int n, int nbig, double *p, double *plb, double *pub, double *I, int method,
		      double *hinit, int *skip, double maxmin, double errf, double tol,
		      double *work, int worksize, void (*func)(int *, double *, double *), int *info )
{
}


/* Central Difference Computation */

static void  scd( int n, int nbig, double *x, double h, double *dx, double *v, double *err,
		  double *w, void (*func)(int *, double *, double *), double maxval )
{
    int       i;
    double    vp, vn, numer, denom;

    *err = 0.0;                      /* No error estimate available, so don't flag problems */
    
    for( i = 0; i < nbig; i++ )
	w[i] = x[i] + h*dx[i];
    
    (*func)( &n, w, &vp );

    for( i = 0; i < nbig; i++ )
	w[i] = x[i] - h*dx[i];

    (*func)( &n, w, &vn );
	    
    numer = 2.0 * maxval - vp - vn;
    denom = h * h;

    if( denom < Under )
	set_warning( WARNING_SMALL_DENOM, 0 );
    else if( fabs(numer) < Under )
	set_warning( WARNING_SMALL_NUMER, 0 );

    *v = numer/denom;
}


static void  fcd( int n, int nbig, double *x, double h, double *dx, double *v, double *err,
		  double *w, void (*func)(int *, double *, double *), double maxval )
{
}


/* Forward Difference Computations */


static void  sfd( int n, int nbig, double *x, double h, double *dx, double *v, double *err,
		  double *w, void (*func)(int *, double *, double *), double maxval )
{
    int       i;
    double    vp, vn, numer, denom;

    
    *err = 0.0;                      /* No error estimate available, so don't flag problems */
    
    for( i = 0; i < nbig; i++ )
	w[i] = x[i] + 2.0*h*dx[i];
    
    (*func)( &n, w, &vp );

    for( i = 0; i < nbig; i++ )
	w[i] = x[i] + h*dx[i];

    (*func)( &n, w, &vn );
	    
    numer = 2.0 * vn - maxval - vp;    /* -( f(x + 2h) - 2 f(x + h) + f(x) ) */
    denom = h * h;

    if( denom < Under )
	set_warning( WARNING_SMALL_DENOM, 0 );
    else if( fabs(numer) < Under )
	set_warning( WARNING_SMALL_NUMER, 0 );

    *v = numer/denom;
}

static void  ffd( int n, int nbig, double *x, double h, double *dx, double *v, double *err,
		  double *w, void (*func)(int *, double *, double *), double maxval )
{
}


/*
 * Richardson Extrapolation algorithm for evaluation of F(0) given 
 * specified values of
 *
 *     F(h) = a_0 + a_1 h^p_1 + a_2 h^p_2 + ...
 *
 * for h small.  Uses deferred approach to the limit with an algorithm
 * that makes use of a Neville tableau for polynomial interpolation.
 *
 * This code is adapted from that originally presented in Numerical Recipes in C,
 * second edition, Press et al., page 188ff.  Configurable parameters below
 * are set as in the book.  See Numerical Methods by A. Bjork (translated by N. Anderson)
 * for a nice description of the method.
 *
 */

#define RICH_NTAB     20

static double const   Large_Value = 1.0e30;
static double const   CON         = 1.4;      /* Stepsize increment */
static double const   CONSQ       = 1.96;     /* CON squared        */
static double const   SAFE        = 2.0;      /* Error Threshold    */


static void   rich_extrap_scd( int n, int nbig, double *x, double h, double *dx, double *v, double *err,
			       double *w, void (*func)(int *, double *, double *), double maxval )
{
    int       i, j;

    double    errt, errsd, fac, hh, ans;
    double    tab[RICH_NTAB][RICH_NTAB];

    if( h == 0.0 )
    {
	set_error( ERROR_BAD_H );
	return;
    }

    hh = h;
    *err = Large_Value;

    scd( n, nbig, x, hh, dx, &tab[0][0], &errsd, w, func, maxval );

    for( i = 1; i < RICH_NTAB; i++ )
    {
	/* Successive columns in the tableau will go to smaller */
        /* stepsizes and higher orders of extrapolation         */

	hh /= CON;
	fac = CONSQ;

	scd( n, nbig, x, hh, dx, &tab[0][i], &errsd, w, func, maxval );

	for( j = 1; j <= i; j++ )
	{
	    /* Compute extrapolations of various orders */
	    /* No new function evaluations are needed   */

	    tab[j][i] = (fac * tab[j-1][i] - tab[j-1][i-1])/(fac - 1.0);
	    fac *= CONSQ;
	    
	        /* Compare each new extrapolation to one order lower, */
	        /* both at the current and previous stepsizes.        */

	    errt = Max( fabs(tab[j][i] - tab[j-1][i]), fabs(tab[j][i] - tab[j-1][i-1]) );

	    if( errt <= *err )
	    {
		*err = errt;
		ans = tab[j][i];
	    }
	}

	/* If higher order is worse by a significant fctor SAFE, then quit early */

	if( fabs( tab[i][i] - tab[i-1][i-1]) >= SAFE * (*err) )
	{
	    *v = ans;
	    return;
	}
    }

    *v = ans;
    return;
}


static void   rich_extrap_sfd( int n, int nbig, double *x, double h, double *dx, double *v, double *err,
			       double *w, void (*func)(int *, double *, double *), double maxval )
{
    int       i, j;

    double    errt, errfd, fac, hh, ans;
    double    tab[RICH_NTAB][RICH_NTAB];

    if( h == 0.0 )
    {
	set_error( ERROR_BAD_H );
	return;
    }

    hh = h;
    *err = Large_Value;

    sfd( n, nbig, x, hh, dx, &tab[0][0], &errfd, w, func, maxval );

    for( i = 1; i < RICH_NTAB; i++ )
    {
	/* Successive columns in the tableau will go to smaller */
        /* stepsizes and higher orders of extrapolation         */

	hh /= CON;
	fac = CON;

	sfd( n, nbig, x, hh, dx, &tab[0][i], &errfd, w, func, maxval );

	for( j = 1; j <= i; j++ )
	{
	    /* Compute extrapolations of various orders */
	    /* No new function evaluations are needed   */

	    tab[j][i] = (fac * tab[j-1][i] - tab[j-1][i-1])/(fac - 1.0);
	    fac *= CON;
	    
	        /* Compare each new extrapolation to one order lower, */
	        /* both at the current and previous stepsizes.        */

	    errt = Max( fabs(tab[j][i] - tab[j-1][i]), fabs(tab[j][i] - tab[j-1][i-1]) );

	    if( errt <= *err )
	    {
		*err = errt;
		ans = tab[j][i];
	    }
	}

	/* If higher order is worse by a significant fctor SAFE, then quit early */

	if( fabs( tab[i][i] - tab[i-1][i-1]) >= SAFE * (*err) )
	{
	    *v = ans;
	    return;
	}
    }

    *v = ans;
    return;
}

static void dummy_func( double *t )
{
}


/*
 * Local Variables:        
 * mode: c
 * eval: (make-local-variable 'compile-command)
 * compile-command: "make -k info.o"
 * End:
 *
 */


