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

static char rcsid[] = "$Id: drift.c,v 1.8 2000/11/02 20:04:47 welling Exp $";

/*
 * drift.c
 *
 * Routines for dealing with smooth drift terms in mri model.
 * The drift is assumed modeled by a low-degree polynomial
 * or a spline with a small number of knots and low degree.
 *
 * Routines herein build appropriate bases, regress on these
 * bases, and compute various derivative bases as well.
 *
 * Several of these routines depend on LAPACK library routines.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include  "fmri.h"
#include  "stdcrg.h"
#include  "lapack.h"
#include  "minqnb.h"
#include  <math.h>

static  int        iOne = 1;
static  int        iZero = 0;
static  int        imOne = -1;

static  float      fOne = 1.0;
static  float      fZero = 0.0;
static  float      fmOne = -1.0;

static  double     One = 1.0;
static  double     Zero = 0.0;
static  double     mOne = -1.0;

/** Public Function Prototypes **/

void     make_drift_basis( int, int, int, double **, double *, double *, double *, int );
void     make_drift_qforms( int, int, int, double *, double *, double *, double *,
			    double *, double *, int );

void     find_drift_df( int D, double *P, double df, double *pu, 
			double *pulow, double *puhigh,
			double *eig, double *work, int lwork, double tol );

void     init_detrend( int, int, int, double *, int );
void     detrend( int, int, int, double *, double **, double *, double *, int );

void     init_binom_coefs( void );


#if !defined(USE_B_SPLINES)

/** Private Function Prototypes **/

double   pb_integral( double a, double b, int d, int D );

/** Macros and Constants **/

static double const Tolerance = 1.0e-8;

enum basis_limits { DEGREE_MAX = 4, KNOT_MAX = 4096, BASIS_MAX = 4100 };

static   int      dt_initialized = 0;

static   double  *dt_design_mat;
static   double   dt_tau[BASIS_MAX + 3];





/*
 * Construct drift basis matrix Q and transformation matrix R
 * from the power basis on the given grid with given knots.
 *
 */

void     make_drift_basis( int nt, int deg, int dk, double **basis, double *knots,
			   double *R, double *work, int lwork )
{
    int      i, j, t, nc, info = 0;
    double   u, v;

    nc = 1 + deg + dk;

    if( nt < nc )
	Abort( "Basis too large for given data vector length. (%d > %d)", nc, nt );

    /* Start with power basis */

    for( t = 0; t < nt; t++ )
    {
	u = (double)t/(nt - 1);
	    
	for( i = 0, v = 1.0; i <= deg; i++, v *= u )
	    R[t + i*nt] = v;

	for( i = 0; i < dk; i++ )
	    R[t + (deg + i + 1)*nt] = (u > knots[i]) ? pow( u - knots[i], (double)deg ) : 0.0;
    }

    /* Orthonormalize with respect to data */

    dgeqrf( &nt, &nc, R, &nt, dt_tau, work, &lwork, &info );

    if( info )
	Abort( "Could not compute QR decomposition in dgeqrf, info = %d.", info );
    else if( work[0] > lwork )
	Warning( DIAGNOSTIC, "Workspace insufficient for optimal block size in QR (%d < %d).",
		 lwork, (int)work[0] );

    for( i = 0, j = nc*nt; i < j; i++ )  /* Set to block identity */
	basis[0][i] = 0.0;
    
    for( i = 0; i < nc; i++ )
	basis[i][i] = 1.0;

    dormqr( &Left, &NoTrans, &nt, &nc, &nc, R, &nt, dt_tau, basis[0],
	    &nt, work, &lwork, &info );

    if( info )
	Abort( "Could not form columns of Q in dormqr, info = %d.", info );

    /* Restore signs to keep consistent signs on diagonal of R but keep basis orthogonal */
    /* Adjust both R and basis according to this sign convention                         */

    for( i = 0; i < nc; i++ )
    {
	u = (R[i + i*nt] > 0.0) ? 1.0 : -1.0;     /* sign(diag(R)) */

	for( t = 0; t < nt; t++ )
	    basis[i][t] *= u;

	for( j = i; j < nc; j++ )
	    R[i + j*nt] *= u; 
    }
}


void     make_drift_qforms( int nt, int D, int dk, double *knots, double *Pn, double *Pc,
			    double *XtX, double *sknots, double *work, int lwork )
{
    int       i, j;
    int       m, ddm1, iim1;
    int       indi, indj, diag;

    double    u, v, w;
    double    hm3, h, ntm1;

    m = 1 + D + dk;  /* Include intercept term */

    /* First: Norm Penalty Matrix */

        /* Pure Polynomial Piece */

    for( i = 0; i <= D; i++ )
    {
	for( j = 0; j < i; j++ )
	    Pn[i + j*m] = Pn[j + i*m] = 1.0/(1.0 + i + j);    /* 1/(1 + d + d')  0 <= d,d' <= D */

	Pn[i + i*m] = 1.0/(1.0 + 2.0*i);
    }

        /* Polynomial . Spline Piece */

    for( i = 0; i <= D; i++ )
    {
	for( j = D + 1; j < m; j++ )
	    Pn[i + j*m] = Pn[j + i*m] = pb_integral( knots[j-D-1], 1.0, i, D );
    }
    
        /* Spline . Spline Piece */

    for( i = D + 1; i < m; i++ )
    {
	for( j = D + 1; j < i; j++ )
	    Pn[i + j*m] = Pn[j + i*m] =
	     pb_integral( fabs(knots[i-D-1] - knots[j-D-1]), 1.0 - Min(knots[i-D-1],knots[j-D-1]), D, D );

	Pn[i + i*m] = pb_integral( 0.0, 1.0 - knots[i-D-1], D, D );
    }

    /* Second: Curvature Penalty Matrix */

    if( D == 1 )       /* Use Second Difference Matrix */
    {
	/* Pure Polynomial Piece */

	Pc[0] = Pc[1] = Pc[m] = Pc[m+1] = 0.0;

        /* Polynomial . Spline Piece */

	for( j = 2; j < m; j++ )
	    Pc[j*m] = Pc[j] = 0.0;
    
        /* Spline . Spline Piece */

	hm3 = Cube(nt-1);
	ntm1 = nt - 1.0;
	h = 1.0/(nt - 1.0);

	for( diag = 1, i = 1; i < dk; i++ )
	    diag *= (sknots[i] - sknots[i-1] > h);  /* Matrix diagonal if all knots separated */

	if( diag )
	{
	    for( i = D + 1; i < m; i++ )
	    {
		u = ceil( ntm1 * knots[i-D-1] )/ntm1 - knots[i-D-1];
		
		for( j = D + 1; j < i; j++ )
		    Pc[i + j*m] = Pc[j + i*m] = 0.0;

		Pc[i + i*m] = hm3 * ( (h - u) * (h - u) + u * u );
	    }
	}
	else
	{
	    for( i = D + 1; i < m; i++ )
	    {
		w = ceil( ntm1 * knots[i-D-1] );
		u = w/ntm1 - knots[i-D-1];
		indi = (int)w;
		
		for( j = D + 1; j < i; j++ )
		{
		    if( fabs(knots[i-D-1] - knots[j-D-1]) < h )
		    {
			w = ceil( ntm1 * knots[j-D-1] );
			v = w/ntm1 - knots[j-D-1];
			indj = (int)w;
		    
			if( indi < indj )
			    Pc[i + j*m] = Pc[j + i*m] = hm3 * (h - u) * v;
			else if( indi > indj )
			    Pc[i + j*m] = Pc[j + i*m] = hm3 * u * (h - v);
			else
			    Pc[i + j*m] = Pc[j + i*m] = hm3 * (u * v + (h - u) * (h - v));
		    }
		    else
			Pc[i + j*m] = Pc[j + i*m] = 0.0;
		}

		Pc[i + i*m] = hm3 * ( (h - u) * (h - u) + u * u );
	    }
	}
    }
    else if( D == 2 )
    {
	/* Pure Polynomial Piece */

	for( i = 0; i <= D; i++ )
	{
	    for( j = 0; j < i; j++ )
		Pc[i + j*m] = Pc[j + i*m] = 0.0;

	    Pc[i + i*m] = 0.0;
	}

	Pc[2 + 2*m] = 4.0;

        /* Polynomial . Spline Piece */

	for( j = D + 1; j < m; j++ )
	{
	    Pc[j]     = Pc[j*m]     = 0.0;    /* Constant Column and Row */
	    Pc[j + m] = Pc[1 + j*m] = 0.0;    /* Linear Column and Row   */
	}

	for( j = D + 1; j < m; j++ )
	    Pc[2 + j*m] = Pc[j + 2*m] = 4.0 * (1.0 - knots[j-D-1]);
    
        /* Spline . Spline Piece */

	for( i = D + 1; i < m; i++ )
	{
	    for( j = D + 1; j < i; j++ )
		Pc[i + j*m] = Pc[j + i*m] = 4.0 * (1.0 - Max( knots[i-D-1], knots[j-D-1] ) );

	    Pc[i + i*m] = 4.0 * (1.0 - knots[i-D-1]);
	}
    }
    else
    {
	ddm1 = D * (D - 1);
	
        /* Pure Polynomial Piece */

	for( i = 0; i < m; i++ )          /* Constant and Linear piece have zero second derivative */
	{
	    Pc[i]     = Pc[    i*m] = 0.0;
	    Pc[i + m] = Pc[1 + i*m] = 0.0;
	}
	    
	for( i = 2; i <= D; i++ )
	{
	    iim1 = i * (i - 1);

	    for( j = 2; j < i; j++ )
		Pc[i + j*m] = Pc[j + i*m] = iim1 * j * (j-1.0)/(i + j - 3.0);

	    Pc[i + i*m] = iim1 * (double)iim1/(2.0*i - 3.0);
	}

        /* Polynomial . Spline Piece */

	for( i = 2; i <= D; i++ )
	{
	    iim1 = i * (i - 1);
	    
	    for( j = D + 1; j < m; j++ )
		Pc[i + j*m] = Pc[j + i*m] = ddm1 * iim1 * pb_integral( knots[j-D-1], 1.0, i-2, D-2 );
	}
    
        /* Spline . Spline Piece */

	for( i = D + 1; i < m; i++ )
	{
	    for( j = D + 1; j < i; j++ )
		Pc[i + j*m] = Pc[j + i*m] = ddm1 * ddm1 *
		    pb_integral( fabs(knots[i-D-1]-knots[j-D-1]), 1.0 - Min(knots[i-D-1],knots[j-D-1]), D-2, D-2 );

	    Pc[i + i*m] = ddm1 * ddm1 * pb_integral( 0.0, 1.0 - knots[i-D-1], D-2, D-2 );
	}
    }
}


/*
 * FIND_DRIFT_DF finds the value of the normalized smoothing parameter $\sigma^2/\lambda$
 * that achieves a specified degree of freedom.   The effective degrees of freedom is
 * defined as the trace of the smoothing matrix.  In the orthogonal spline case,
 * this is of the form
 *
 *                 trace( (I + u P)^-1 )
 *
 * where u = \sigma^2/\lambda = 1/(r \lambda).  Since the identity is in this expression,
 * we can decompose P = U D U' to get
 *
 *   trace( U (I + u D)^-1 U' ) = trace( U U' (I + u D)^-1 ) = trace( (I + u D)^-1 )
 *
 * Hence, we need only find the eigenvalues eta_i of P and set
 *
 *    Sum[ i = 0, D - 1;   1/(1 + u \eta_i) ] =  c.
 *
 *
 * This function returns the corresponding value of u, the list of eigenvalues, and
 * the values ulow and uhigh that correspond to deviations of 1 degree of freedom in
 * either direction.
 *
 * The B-spline version of this function, given later in this file, is more complicated.
 *
 */


void     find_drift_df( int D, double *P, double df, double *pu, double *pulow, double *puhigh,
			double *eig, double *work, int lwork, double tol )
{
    int       i;
    int       info = 0;

    double    u, ulow, uhigh, ua, ub;
    double    v, vlow, vhigh, va, vb;
    double    sumeig, sumeigsq;
    double    deriv, tmp;

    int       j;
    int       IWork[8192];
    char      range = 'A';

    if( df >= D || df <= 0.0 )
	Abort( "Target degrees of freedom (%lg) out of range, need 0 < df < %d.", df, D );

    dsyev( &No, &Upper, &D, P, &D, eig, work, &lwork, &info );

    if( info )
	Abort( "Cannot find eigendecomposition of smoothing matrix (info = %d).", info );

    /*
     * Use Taylor approximation  Sum[ i; 1/(1 + u \eta_i) ] ~ Sum[ i; (1 - u\eta_i + u^2 \eta_i^2)].
     * to set initial value.
     *
     */

    for( i = 0, sumeig = sumeigsq = 0.0; i < D; i++ )
    {
	if( eig[i] < 0.0 )  /* Any negatives are numerical glitches, P is pos. def. */
	    eig[i] = 0.0;

	sumeig += eig[i];
	sumeigsq += eig[i] * eig[i];
    }

    if( sumeig == 0.0 )
	Abort( "Zero trace encountered in search for drift smoothing target. This shouldn't happen." );

    /* u = (sumeig - sqrt( sumeig*sumeig - 4.0*(D - df)*sumeigsq ))/(2.0*sumeigsq); */

    u = (D - df)/sumeig;

    ulow = uhigh = -1.0;
    ua = ub = -1.0;
    vlow = -1.0e10;
    vhigh = 1.0e10;

    while( ua < 0.0 || ub < 0.0 )
    {
	for( i = 0, v = -df, deriv = 0.0; i < D; i++ )
	    v += 1.0/(1.0 + u*eig[i]);

	if( v < 0.0 )
	{
	    ub = u;
	    vb = v;

	    if( fabs(v + 1.0) < fabs(vhigh + 1.0) ) /* Record bounds closest to 1 */
	    {
		uhigh = u;
		vhigh = v;
	    }

	    if( ua >= 0.0 )
	    {
		u = ua + (ub - ua) * va/(va - vb);  /* Inverse Linear Interpolation */
		/* u = 0.5 * (ua + ub); Simple Bisection */
	    }
	    else
		u *= 0.61803398874989484820;
	}
	else if( v > 0.0 )
	{
	    ua = u;
	    va = v;

	    if( fabs(v - 1.0) < fabs(vlow - 1.0) ) /* Record bounds closest to 1 */
	    {
		ulow = u;
		vlow = v;
	    }

	    if( ub >= 0.0 )
	    {
		u = ua + (ub - ua) * va/(va - vb);  /* Inverse Linear Interpolation */
		/* u = 0.5 * (ua + ub); Simple Bisection */
	    }
	    else
		u *= 1.61803398874989484820;
	}
	else
	{
	    ub = u * 1.61803398874989484820;
	    ua = u * 0.61803398874989484820;
	}
    }

    /* Search from here with Newton steps */

    for( i = 0, v = -df, deriv = 0.0; i < D; i++ )
    {
	tmp = 1.0/(1.0 + u*eig[i]);
	v += tmp;
	deriv -= eig[i] * (tmp * tmp);
    }

    /* vlast = vlast2 = 1.0e10; */

    while( fabs(v) >= tol )
    {
	/*
        vlast2 = vlast;
	vlast = v;

	if( fabs( v - vlast2 ) < tol/2.0 )
	    Restart at random location in (ua,ub) ; 
	*/

	
	if( u - v/deriv <= ua )
	{
	    tmp = 0.5 * (ua + ub);
	    ua = u;
	    u = tmp;
	}
	else if( u - v/deriv >= ub )
	{
	    tmp = 0.5 * (ua + ub);
	    ub = u;
	    u = tmp;
	}
	else
	    u -= v/deriv;

	for( i = 0, v = -df, deriv = 0.0; i < D; i++ )
	{
	    tmp = 1.0/(1.0 + u*eig[i]);
	    v += tmp;
	    deriv -= eig[i] * (tmp * tmp);
	}
    }

    /* Update the low and high points */

    if( pulow && puhigh )
    {
	if( ulow < 0.0 && uhigh < 0.0 )
	{
	    uhigh = u * 1.61803398874989484820;
	    ulow  = u * 0.61803398874989484820;

	    for( i = 0, vlow = vhigh = -df; i < D; i++ )
	    {
		vhigh += 1.0/(1.0 + uhigh*eig[i]);
		vlow  += 1.0/(1.0 +  ulow*eig[i]);
	    }
	}
    }
    
        /* FILL IN HERE */

    /* Record and Done */

    *pu = u;

    if( pulow && puhigh )
    {
	*puhigh = uhigh;
	*pulow = ulow;
    }
}


static int made_binom_coef = 0;
static double binom_coef[DEGREE_MAX + 1][DEGREE_MAX + 1];

double   pb_integral( double a, double b, int d, int D )
{
    int       i;

    double    u, v, w, x, y;

    /* assert( 0.0 <= a && a < b && b <= 1.0 ) */           /* Assume 0 < a < b < 1 */

    if( a < Tolerance && b < 1.0 )
	return( pow( b, (double)(1 + d + D) )/(1 + d + D) );
    else if( a < Tolerance )
	return( 1.0/(1 + d + D) );
    else
    {
	u = (1 - 2 * (D % 2)) * pow( a, (double)(1 + d + D) );  /* (-1)^D a^(1+d+D) */
	v = 0.0;
	w = b/a;
	x = pow( w, (double)(1+d) );
	y = 1.0;

	for( i = 0; i <= D; i++ )
	{
	    v += y * binom_coef[D][i] * (x - 1.0)/(double)(1 + d + i);

	    x *= w;
	    y *= -1.0;
	}

	return( u * v );
    }
}

void   init_binom_coefs( void )
{
    int           i, d;
    double        u;
    
    if( !made_binom_coef )    /* Set up-binomial coefs */
    {
	binom_coef[0][0] = 1.0;

	for( d = 1; d <= DEGREE_MAX; d++ )
	{
	    binom_coef[d][0] = binom_coef[d][d] = 1.0;

	    for( i = 1, u = 1.0; i <= d/2; i++ )
	    {
		u *= ((double)(d + 1 - i))/i;
		binom_coef[d][i] = binom_coef[d][d-i] = u;
	    }

	    for( i = d + 1; i <= DEGREE_MAX; i++ )
		binom_coef[d][i] = 0.0;
	}
	

	made_binom_coef = 1;
    }
}


/*** Detrending for finding initial approximate max ***/

void     init_detrend( int nt, int deg, int nknots, double *work, int lwork )
{
    /* Nothing to do here */
}

void     detrend( int nt, int deg, int nknots, double *y, double **basis, double *pars,
		  double *work, int lwork )
{
    int       i, j, l, nc, info = 0;

    if( !dt_initialized )
    {
	dt_design_mat = Malloc( (deg + nknots + 1) * nt, double );
	dt_initialized = 1;
    }

    /* Set-up Design Matrix.  */

    nc = 1 + deg + nknots;
    info = 0;
    
    for( j = 0; j < nc; j++ )
	for( i = 0; i < nt; i++ )
	    dt_design_mat[i + j*nt] = basis[j][i];

    dgeqrf( &nt, &nc, dt_design_mat, &nt, dt_tau, work, &lwork, &info );

    if( info )
	Abort( "Could not compute QR decomposition in dgeqrf, info = %d.", info );

    
    /* Least Squares Computation: Residual and Coefficients */

        /* Q'y */

    dormqr( &Left, &Trans, &nt, &iOne, &nc, dt_design_mat, &nt, dt_tau, y, &nt, work, &lwork, &info );

    if( info )
	Abort( "Could not form Q'y in dormqr, info = %d.", info );

        /* Save coefs and zero out first part of y for resid computation */

    for( l = 0; l < nc; l++ )
    {
	pars[l] = y[l];
	y[l] = 0.0;
    }

        /* Compute Coefs by solving triangular system: (Q'y) = R pars */

    dtrtrs( &Upper, &NoTrans, &NoUnit, &nc, &iOne, dt_design_mat, &nt, pars, &nc, &info );

    if( info )
	Abort( "Could not compute detrend coefs in dtrtrs, info = %d.", info );

    /* Compute residuals by applying Q */

    dormqr( &Left, &NoTrans, &nt, &iOne, &nc, dt_design_mat, &nt, dt_tau, y, &nt, work, &lwork, &info );
}


#else

#endif


