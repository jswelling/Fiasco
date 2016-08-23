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

static char rcsid[] = "$Id: addendum.c,v 1.3 2000/10/06 05:32:04 welling Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Prototypes */

void   lnpost_fixed( int *npar, double *p, double *val );
void   Dlnpost_fixed( int *npar, double *p, double *deriv );
void   lnpost_variable( int *npar, double *p, double *val );
void   Dlnpost_variable( int *npar, double *p, double *deriv );
void   profiles_aligned( double *p, double *dr_prof, double *ac_prof, double *work, int lwork );
void   profiles_unaligned( double *p, double *dr_prof, double *ac_prof, double *work, int lwork );
void   poly_bell4( long t, double *b, long *blen, double *p, double stimlen, double offset, double IAI );
void   poly_bell6( long t, double *b, long *blen, double *p, double stimlen, double offset, double IAI );
void   poly_bell8( long t, double *b, long *blen, double *p, double stimlen, double offset, double IAI );
void   shape_deriv_aligned4( double *p, double **dbells, double *stims, double offset, double IAI );
void   shape_deriv_aligned6( double *p, double **dbells, double *stims, double offset, double IAI );
void   shape_deriv_aligned8( double *p, double **dbells, double *stims, double offset, double IAI );


/* Source     */

/*
 * LNPOST_FIXED evaluates the NEGATIVE of the log-likelihood at the given point.
 * Assumes that the drift knots (and thus the drift basis) are fixed.
 *
 * Side Effects: Sets the values of R and SSR to the residual vector and
 *               residual sum of squares respectively.
 *
 */

void   lnpost_fixed( int *npar, double *p, double *val )
{
    char     noevecs = 'N';
    
    int      i, j, m;
    int      d, k, s, t;
    int      lwork, info = 0;

    double   c, u, v;
    double   lnr, ll;
    double   dscaling;

    double   knots[KNOT_MAX];

    /* Compute and store profiles  */

    if( Design_Aligned )
	profiles_aligned( p, Drift_prof, Active_prof, Work, WorkSize );
    else
	profiles_unaligned( p, Drift_prof, Active_prof, Work, WorkSize );


    /* Construct Residuals and RSS */

    for( t = 0, SSR = 0.0; t < T; t++ )
    {
        R[t] = u = Y[t] - p[Mu] - Drift_prof[t] - Active_prof[t];
        SSR += u*u;
    }

    SSR *= 0.5;


    /* Likelihood Term */

    lnr = Log(p[o_SigmaSq]);

    ll = SSR * p[o_SigmaSq] - 0.5 * T * lnr;
    
    /* Prior Terms */
                                                          /* Baseline */

    ll += Log( 1.0 + Mu_b * (p[Mu] - Mu_a)*(p[Mu] - Mu_a) );

                                                          /* Drift    */

    dscaling = (Drift_scaled) ? (1.0/(p[Mu] * p[Mu])) : 1.0;

    /* Compute -0.5 delta' R^-1' P R^-1 delta */

    Work[0] = 0.0;    /* Drift profile orthogonal to constants, so intercept zero */

    for( i = 0; i < D; i++ )
	Work[1+i] = p[Coefs[i]];

    m = D + 1;  /* Include intercept term because 0 int in Q basis => non-zero int in power basis */
    
    dtrtrs( &Upper, &NoTrans, &NoUnit, &m, &iOne, basisR, &T, Work, &m, &info );

    if( info )
	Warning( 1, "Could not solve triangular system required for drift prior (info = %d).", info );

    for( i = 0, u = 0.0; i < m; i++ )       
    {
	u += Work[i] * Work[i] * Pcomb[i + i*m];
	
	for( j = 0; j < i; j++ )
	    u += 2.0 * Work[i] * Work[j] * Pcomb[j + i*m];
    }

    ll += 0.5 * dscaling * u;
    

    for( k = 0; k < K; k++ )                              /* Responsiveness */
	ll += Resp_b * p[Resp[k]];
    
    for( s = 0; s < SHAPE_PARAMS; s++ )                   /* Shape */
    {
        ll += Shape_b[s]*p[Shape[s]] - (Shape_a[s] - 1.0) * Log(p[Shape[s]]);
    }

                                                          /* o_SigmaSq */

    ll += 0.5 * (lnr + SigmaSq_a)*(lnr + SigmaSq_a)*o_SigmaSq_bsq - lnr;

    *val = ll;
}

void   Dlnpost_fixed( int *npar, double *p, double *deriv )
{
    int      d, i, k, s, t;

    double   v, u, r;

    double   *df;

    /* Compute Profiles and Profile Derivatives */

    if( Design_Aligned )
    {
	profiles_aligned( p, Drift_prof, Active_prof, Work, WorkSize );
	shape_deriv_aligned( p, Shape_Params, DBells, Work, WorkSize );
    }
    else
    {
	profiles_unaligned( p, Drift_prof, Active_prof, Work, WorkSize );
	shape_deriv_unaligned( p, Shape_Params, DBells, Work, WorkSize );
    }

    /* Construct Residuals and RSS */

    for( t = 0, SSR = 0.0; t < T; t++ )
    {
        R[t] = u = Y[t] - p[Mu] - Drift_prof[t] - Active_prof[t];
        SSR += u*u;
    }

    SSR *= 0.5;

    /* Set derivative components one at a time */

    r = p[o_SigmaSq];

        /* Baseline */

    for( t = 0, u = 1.0/p[Mu], v = 0.0; t < T; t++ )
        v += R[t] * (1.0 + u * Active_prof[t]);

    u = p[Mu] - Mu_a;                     

    deriv[Mu] = -v * r  +                                 /* Likelihood Term */
                (2.0 * Mu_b * u/(1.0  +  Mu_b * u * u));  /* Prior Term      */


        /* Drift */

            /* Prepare Prior Term First:  Compute -(R^-1' P R^-1) delta */

    Work[0] = 0.0;    /* Drift profile orthogonal to constants, so intercept zero */

    for( i = 0; i < D; i++ )
	Work[1+i] = p[Coefs[i]];

    m = D + 1;  /* Include intercept term because 0 int in Q basis => non-zero int in power basis */
    
    dtrtrs( &Upper, &NoTrans, &NoUnit, &m, &iOne, basisR, &T, Work, &m, &info );

    if( info )
	Warning( 1, "Could not solve triangular system for drift prior deriv (info = %d).", info );

    for( i = 0; i < m; i++ )           /* Compute  P (R^-1 delta) */
    {
	for( j = 0, u = 0.0; j < m; j++ )
	    u += Work[j] * Pcomb[j + i*m];

	Work[m + i] = u;
    }

    dtrtrs( &Upper, &Trans, &NoUnit, &m, &iOne, basisR, &T, Work + m, &m, &info );

    if( info )
	Warning( 1, "Couldn't solve 2nd triangular system for drift prior deriv (info = %d).", info );

    dscaling = (Drift_scaled) ? (1.0/(p[Mu] * p[Mu])) : 1.0;

            /* Now combine with Likelihood Term */

    for( d = 0; d < D; d++ )
    {
	for( t = 0, df = Basis[d], v = 0.0; t < T; t++ )
	    v += R[t] * df[t];

	deriv[Coefs[d]] = -v * r  +  dscaling * Work[m + 1 + d];     /* Likelihood + Prior Term */

    }

        /* Responsiveness */

    for( k = 0, u = p[Mu]; k < Keff; k++ )
    {
	for( i = CondBlock[RespEff + k][0], j = 0, v = 0.0; j < i; j++ )
	{
	    s = CondBlock[RespEff + k][j];
	    d = 2 * (RespEff + k - 1);
	    
	    for( t = BlockStart[s]; t < BlockEnd[s]; t++ )
		v += R[t] * Bells[t + d*T];
	}

	deriv[RespEff + k] = -v * r * u + Resp_b;  /* Likelihood + Prior Term */  
    }


        /* Shape */    
    
    for( s = 0; s < Shape_Params; s++ )
	deriv[Shape[s]] = (Shape_b[s] - (Shape_a[s] - 1.0)/p[Shape[s]]);   /* Prior Term      */

    for( b = 0, v = 0.0; b < B; b++ )
    {
	k = BlockStart[b];
	j = WhichStim[b];
	hgt = p[Mu] * p[Resp[Cond[b]]];

	if( hgt > 0.0 )
	{
	    l = Min( Bell_Length[j], T - k - 1 );

	    for( s = 0, v = 0.0; s < Shape_Params; s++ )
	    {
		for( i = 0, df = DBells[s]; i <= l; i++ )
		    v += R[i + k] * df[i + j*T];

		deriv[Shape[s]] += -v * r;    /* Likelihood Term */
	    }
	}
    }
	    

        /* o_SigmaSq */

    deriv[o_SigmaSq] = (SSR * r * r - 0.5 * T * r) +                        /* Likelihood Term */
                       r * ( 1.0 - (Log(r) + SigmaSq_a) * o_SigmaSq_bsq );  /* Prior Term      */
}


void   lnpost_variable( int *npar, double *p, double *val )
{
}

void   Dlnpost_variable( int *npar, double *p, double *deriv )
{
}


void   profiles_aligned( double *p, double *dr_prof, double *ac_prof, double *work, int lwork )
{
    int       i, j, k, l;
    int       b;

    double    len, hgt;

    /*
     * Drift Profile
     *
     * Need to construct basis for current drift knots and second derivative (wrt time) profile
     *
     * Assumes that Mu and Drift[0] are contiguous.
     *
     */

    make_drift_profs( T, Dg, Dk, Basis, &p[Drift[0]], dr_prof, work, lwork );

    /* Generate Activation Profile */

      /* Set Bells for each unique stimulus length */
      /* Do not compute bell for fixed blocks      */

    for( i = 0; i < UniqueStims; i++ )
        poly_bell( T, &Bells[i*2*T], &Bell_Length[i], p, Stims[i] );
    
      /* Combine over blocks using smooth_paste() */

          /* Set first block then initialize profile */

    j = WhichStim[0];
    hgt = p[Mu] * p[Resp[Cond[b]]];
    l = Bell_Length[j];

    for( i = 0; i <= l; i++ )
	ac_prof[i] = hgt * Bells[i + j*T];

    for( ; i < T; i++ )
	ac_prof[i] = 0.0;

          /* Set remaining blocks */

    for( b = 1; b < B; b++ )
    {
        k = BlockStart[b];
        j = WhichStim[b];
        hgt = p[Mu] * p[Resp[Cond[b]]];

	if( hgt > 0.0 )
	{
	    l = Min( Bell_Length[j], T - k - 1 );

	    for( i = 0; i <= l; i++ )
		ac_prof[i + k] += hgt * Bells[i + j*T];
	}
    }
}

void   profiles_unaligned( double *p, double *dr_prof, double *ac_prof, double *work, int lwork )
{
    int       i, j, k, l;
    int       b;

    double    len, hgt;

    /*
     * Drift Profile
     *
     * Need to construct basis for current drift knots and second derivative (wrt time) profile
     *
     * Assumes that Mu and Drift[0] are contiguous.
     *
     */

    make_drift_profs( T, Dg, Dk, Basis, &p[Drift[0]], dr_prof, work, lwork );

    /* Generate Activation Profile */

      /* Set Bells for each unique stimulus length */
      /* Do not compute bell for fixed blocks      */

    for( i = 0; i < UniqueStims; i++ )
        poly_bell( T, &Bells[i*2*T], &Bell_Length[i], p, Stims[i] );
    
      /* Combine over blocks using smooth_paste() */

          /* Set first block then initialize profile */

    j = WhichStim[0];
    hgt = p[Mu] * p[Resp[Cond[b]]];
    l = Bell_Length[j];

    for( i = 0; i <= l; i++ )
	ac_prof[i] = hgt * Bells[i + j*T];

    for( ; i < T; i++ )
	ac_prof[i] = 0.0;

          /* Set remaining blocks */

    for( b = 1; b < B; b++ )
    {
        k = BlockStart[b];
        j = WhichStim[b];
        hgt = p[Mu] * p[Resp[Cond[b]]];

	if( hgt > 0.0 )
	{
	    l = Min( Bell_Length[j], T - k - 1 );

	    for( i = 0; i <= l; i++ )
		ac_prof[i + k] +=  hgt * Bells[i + j*T];
	}
    }
}

/*
 * POLY_BELL: Aligned version; 4, 6, and 8 parameter models
 *
 * Computes a bell function corresponding to a given stimulus length,
 * slight offset, and inter-acquisition interval.  Assumes that the
 * condition beginnings are aligned with the image acquisitions.
 *
 * The bells are pre-computed on a fine grid in [0,1] but must be interpolated
 * onto the image acquisition grid.  ...
 *
 * Assuming that the image acquisition grid is aligned with the stimulus presentation
 * times (i.e., every stimulus presentation lies on an image acquisition time).
 * When this assumption does not hold, unaligned profiles should be used
 * (see profiles_unaligned).  Given that the grids are aligned, we assume that
 * the fine grid is chosen to fit some multiple into a single inter-image space.
 * 
 *
 * The pre-computed functions are AttackRamp and DecayRamp.  Both are vectors
 * of length RampLen + 1; they sample an attack and decay ramp between 0 and 1
 * over a fine regular grid over [0,1].  (Actually a bit longer and padded with
 * 0's and 1's as appropriate to prevent overrun.)
 *
 * Eventually this will include a dip and varying rise/fall rates.
 *
 * WARNING:  Care must be taken that the rounding schemes used here are
 *           portable.  Use of floor() throughout would greatly increase
 *           the expense, but this approach should be handled carefully.
 *
 *
 */

void  poly_bell4( long t, double *b, long *blen, double *p, double stimlen, double offset, double IAI )
{
    int       i, ja, jd;
    int       len = 0;

    double    u, v;
    double    a1, a2, d1, d2;

    double    lag_on, lag_off;
    double    attack, decay;

    double    awgt, dwgt;    /* Interpolation Weights       */
    double    astart, aend;  /* Beginning and End of Attack */
    double    dstart, dend;  /* Beginning and End of Decay  */

    /* Store shape parameters for later use */

    lag_on = p[Shape[LAG_ON]];
    lag_off = p[Shape[LAG_OFF]];
    attack = p[Shape[ATTACK]];
    decay = p[Shape[DECAY]];

    /*
     * Compute b = AB * DB
     *
     * Compute AB and DB piecewise so that we loop over the vector only
     * once, multiplying them together as we go.
     *
     */

        /* Mark important boundaries */

    astart = Min( (lag_on)/IAI, t );
    aend   = Min( (lag_on + attack)/IAI, t );
    dstart = Min( (stimlen + lag_off)/IAI, t );
    dend   = Min( (stimlen + lag_off + decay)/IAI, t );

        /* Pre-compute values that will be needed */

    a1 = RampLen * (offset - lag_on)/attack;
    a2 = RampLen * IAI/attack;
    d1 = RampLen * (offset - stimlen - lag_off)/decay;
    d2 = RampLen * IAI/decay;

        /* Create the bell function */

    for( i = 0; i < astart; i++ )
	b[i] = 0.0;

    if( dstart < aend )
    {
	for( ; i < dstart; i++ )
	{
	    u = a1 + i * a2;
	    ja = u;           /* Want effect of floor(u) here */
	    awgt = u - ja;
	 
	    b[i] = (1.0 - awgt) * AttackRamp[ja]  +  awgt * AttackRamp[ja + 1];
	}

	if( aend < dend )
	    for( ; i <= aend; i++ )
	    {
		u = a1 + i * a2;
		ja = u;           /* Want effect of floor(u) here */
		awgt = u - ja;
	 
		v = d1 + i * d2;
		jd = v;           /* Want effect of floor(v) here */
		dwgt = v - jd;

		b[i] = ((1.0 - awgt) * AttackRamp[ja]  +  awgt * AttackRamp[ja + 1]) *
	   	       ((1.0 - dwgt) * DecayRamp[jd]   +  dwgt * DecayRamp[jd + 1]);
	    }
	else
	    for( ; i <= dend; i++ )
	    {
		u = a1 + i * a2;
		ja = u;           /* Want effect of floor(u) here */
		awgt = u - ja;
	 
		v = d1 + i * d2;
		jd = v;           /* Want effect of floor(v) here */
		dwgt = v - jd;

		b[i] = ((1.0 - awgt) * AttackRamp[ja]  +  awgt * AttackRamp[ja + 1]) *
	   	       ((1.0 - dwgt) * DecayRamp[jd]   +  dwgt * DecayRamp[jd + 1]);
	    }
    }
    else
    {
	for( ; i < aend; i++ )
	{
	    u = a1 + i * a2;
	    ja = u;           /* Want effect of floor(u) here */
	    awgt = u - ja;
	 
	    b[i] = (1.0 - awgt) * AttackRamp[ja]  +  awgt * AttackRamp[ja + 1];
	}

	for( ; i < dstart; i++ )
	    b[i] = 1.0;
    }
    

    for( ; i < dend; i++ )
    {
	v = d1 + i * d2;
	jd = v;           /* Want effect of floor(v) here */
	dwgt = v - jd;

	b[i] = (1.0 - dwgt) * DecayRamp[jd]  +  dwgt * DecayRamp[jd + 1];
    }

    b[i] = 0.0;
    *blen = i;  /* No need for Min( i, t ) since Min's done above for dend */
}


void  poly_bell6( long t, double *b, long *blen, double *p, double stimlen, double offset, double IAI )
{
}

void  poly_bell8( long t, double *b, long *blen, double *p, double stimlen, double offset, double IAI )
{
}


/*
 * SHAPE_DERIV:  Aligned version; 4, 6, or 8 parameter models
 *
 * Computes derivative bells corresponding to each shape parameter
 * for every unique stimulus length.
 *
 *
 */

void      shape_deriv_aligned4( double *p, double **dbells, double *stims, double offset, double IAI )
{
    int       i, ja, jd, k;
    int       len = 0;

    double    u, v;
    double    a1, a2, d1, d2;

    double    lag_on, lag_off;
    double    attack, decay;
    double    o_na, o_nasq, o_nd, o_ndsq;

    double    awgt, dwgt;    /* Interpolation Weights       */
    double    astart, aend;  /* Beginning and End of Attack */
    double    dstart, dend;  /* Beginning and End of Decay  */

    double    *b;

    /* Store shape parameters for later use */

    lag_on = p[Shape[LAG_ON]];
    lag_off = p[Shape[LAG_OFF]];
    attack = p[Shape[ATTACK]];
    decay = p[Shape[DECAY]];

    o_na = -1.0/attack;
    o_nasq = -IAI/(attack * attack);
    o_nd = -1.0/decay;
    o_ndsq = -IAI/(decay * decay);


    /*
     * Compute AB and DB and their derivatives piecewise so
     * that we loop over the vector only once, multiplying them together as we go.
     * The derivatives for each of the parameters is computed in the inner loop
     * to avoid repeating the overhead of re-computing the interpolation.
     *
     * See poly_bell for discussion of this algorithm.
     *
     */

    for( k = 0; k < UniqueStims; k++ )
    {
	stimlen = stims[k];
	ind = 2*k*T;

        /* Mark important boundaries */

	astart = Min( (lag_on)/IAI, t );
	aend   = Min( (lag_on + attack)/IAI, t );
	dstart = Min( (stimlen + lag_off)/IAI, t );
	dend   = Min( (stimlen + lag_off + decay)/IAI, t );

        /* Pre-compute values that will be needed */

	a1 = RampLen * (offset - lag_on)/attack;
	a2 = RampLen * IAI/attack;
	d1 = RampLen * (offset - stimlen - lag_off)/decay;
	d2 = RampLen * IAI/decay;

        /* Create the bell function */

	for( i = 0; i < astart; i++ )
	{
	    dbells[LAG_ON][i + ind] = 0.0;
	    dbells[ATTACK][i + ind] = 0.0;
	    dbells[LAG_OFF][i + ind] = 0.0;
	    dbells[DECAY][i + ind] = 0.0;
	}
	
	if( dstart < aend )
	{
	    for( ; i < dstart; i++ )
	    {
		u = a1 + i * a2;
		ja = u;           /* Want effect of floor(u) here */
		awgt = u - ja;

		dbells[LAG_ON][i + ind] = o_na * ((1.0 - awgt)*DAttackRamp[ja] + awgt*DAttackRamp[ja + 1]);
		dbells[ATTACK][i + ind] = o_nasq * (i - astart + offset) *
		                                 ((1.0 - awgt)*DAttackRamp[ja] + awgt*DAttackRamp[ja + 1]);
		dbells[LAG_OFF][i + ind] = 0.0;
		dbells[DECAY][i + ind] = 0.0;
	    }

	    if( aend < dend )
		for( ; i <= aend; i++ )
		{
		    u = a1 + i * a2;
		    ja = u;           /* Want effect of floor(u) here */
		    awgt = u - ja;
	 
		    v = d1 + i * d2;
		    jd = v;           /* Want effect of floor(v) here */
		    dwgt = v - jd;

		    dbells[LAG_ON][i + ind] = o_na *
			                         ((1.0 - awgt)*DAttackRamp[ja] + awgt*DAttackRamp[ja + 1]) *
                  			         ((1.0 - dwgt)*DecayRamp[jd]   + dwgt*DecayRamp[jd + 1]);
		    dbells[ATTACK][i + ind] = o_nasq * (i - astart + offset) *
 		                                 ((1.0 - awgt)*DAttackRamp[ja] + awgt*DAttackRamp[ja + 1]) *
                  			         ((1.0 - dwgt)*DecayRamp[jd]   + dwgt*DecayRamp[jd + 1]);
		    dbells[LAG_OFF][i + ind] = o_nd *
			                         ((1.0 - awgt)*AttackRamp[ja] + awgt*AttackRamp[ja + 1]) *
                  			         ((1.0 - dwgt)*DDecayRamp[jd]   + dwgt*DDecayRamp[jd + 1]);
		    dbells[DECAY][i + ind] = o_ndsq * (i - dstart + offset) *
 		                                 ((1.0 - awgt)*AttackRamp[ja] + awgt*AttackRamp[ja + 1]) *
                  			         ((1.0 - dwgt)*DDecayRamp[jd]   + dwgt*DDecayRamp[jd + 1]);
		}
	    else
		for( ; i <= dend; i++ )
		{
		    u = a1 + i * a2;
		    ja = u;           /* Want effect of floor(u) here */
		    awgt = u - ja;
	 
		    v = d1 + i * d2;
		    jd = v;           /* Want effect of floor(v) here */
		    dwgt = v - jd;

		    dbells[LAG_ON][i + ind] = o_na *
			                         ((1.0 - awgt)*DAttackRamp[ja] + awgt*DAttackRamp[ja + 1]) *
                  			         ((1.0 - dwgt)*DecayRamp[jd]   + dwgt*DecayRamp[jd + 1]);
		    dbells[ATTACK][i + ind] = o_nasq * (i - astart + offset) *
 		                                 ((1.0 - awgt)*DAttackRamp[ja] + awgt*DAttackRamp[ja + 1]) *
                  			         ((1.0 - dwgt)*DecayRamp[jd]   + dwgt*DecayRamp[jd + 1]);
		    dbells[LAG_OFF][i + ind] = o_nd *
			                         ((1.0 - awgt)*AttackRamp[ja] + awgt*AttackRamp[ja + 1]) *
                  			         ((1.0 - dwgt)*DDecayRamp[jd]   + dwgt*DDecayRamp[jd + 1]);
		    dbells[DECAY][i + ind] = o_ndsq * (i - dstart + offset) *
 		                                 ((1.0 - awgt)*AttackRamp[ja] + awgt*AttackRamp[ja + 1]) *
                  			         ((1.0 - dwgt)*DDecayRamp[jd]   + dwgt*DDecayRamp[jd + 1]);
		}
	}
	else
	{
	    for( ; i < aend; i++ )
	    {
		u = a1 + i * a2;
		ja = u;           /* Want effect of floor(u) here */
		awgt = u - ja;
	 
		dbells[LAG_ON][i + ind] = o_na * ((1.0 - awgt)*DAttackRamp[ja] + awgt*DAttackRamp[ja + 1]);
		dbells[ATTACK][i + ind] = o_nasq * (i - astart + offset) *
		                                 ((1.0 - awgt)*DAttackRamp[ja] + awgt*DAttackRamp[ja + 1]);
		dbells[LAG_OFF][i + ind] = 0.0;
		dbells[DECAY][i + ind] = 0.0;
	    }

	    for( ; i < dstart; i++ )
	    {
		dbells[LAG_ON][i + ind] = 0.0;
		dbells[ATTACK][i + ind] = 0.0;
		dbells[LAG_OFF][i + ind] = 0.0;
		dbells[DECAY][i + ind] = 0.0;
	    }
	}
    

	for( ; i < dend; i++ )
	{
	    v = d1 + i * d2;
	    jd = v;           /* Want effect of floor(v) here */
	    dwgt = v - jd;

	    dbells[LAG_ON][i + ind] = 0.0;
	    dbells[ATTACK][i + ind] = 0.0;
	    dbells[LAG_OFF][i + ind] = o_nd * ((1.0 - dwgt)*DDecayRamp[jd] + dwgt*DDecayRamp[jd + 1]);
	    dbells[DECAY][i + ind] = o_ndsq * (i - dstart + offset) *
		                                 ((1.0 - dwgt)*DDecayRamp[jd] + dwgt*DDecayRamp[jd + 1]);
	}

	b[i] = 0.0;
    }
}



/* Drift Profile Computations */


#if !defined(NO_INLINE_DRIFTPROF)

    /*
     * Assume that Basis has been just set to contain basis without constant term.
     *
     * Now, construct profile.  (Do not include constant term)
     *
     */

#define  make_drift_profs( nt, deg, dk, basis, coefs, prof, work, lwork ) \
{\
    int        nc;\
    nc = (deg) + (dk);\
    dgemv( &NoTrans, &(nt), &(nc), &One, (basis)[0], &(nt), (coefs), &iOne, &Zero, (prof), &iOne );\
}\

#else

void     make_drift_profs( int nt, int deg, int dk, double **basis, double *coefs,
			   double *prof, double *work, int lwork )
{
    int        nc;
    
    /*
     * Assume that Basis has been just set to contain basis without constant term.
     *
     * Now, construct profile.  (Do not include constant term)
     *
     */

    nc = deg + dk;

    dgemv( &NoTrans, &nt, &nc, &One, basis[0], &nt, coefs, &iOne, &Zero, prof, &iOne );
}

#endif
