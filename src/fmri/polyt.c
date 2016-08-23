/************************************************************
 *                                                          *
 *  polyt.c                                                 *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *                                                          *
 *  Copyright (c) 1999 Department of Statistics             *
 *                     Carnegie Mellon University           *
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
 ************************************************************/
/* polyt.f -- translated by f2c (version of 3 February 1990  3:36:42).
   You must link the resulting object file with the libraries:
	-lf77 -li77 -lm -lc   (in that order)
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"

static char rcsid[] = "$Id: polyt.c,v 1.8 2005/02/07 07:08:04 welling Exp $";

#define TTEST( val, mean, var, tcritsqr ) \
  (((val-mean)*(val-mean))>(tcritsqr*var))

/* Subroutine */ 
int nelmin_t(float (*fn)(float*,void*), void (*rst)(float*,void*), int* n, 
	     float* start, float* mmin, float* ynewl, float* reqmin, 
	     float* step, int* konvge, int* kcount, int* icount, 
	     int* numres, int* ifault, int* maxres, float* tcritsqr,
	     void* userHook)
{
    /* Initialized data */

    static float rcoeff = 1.;
    static float ecoeff = 2.;
    static float ccoeff = .5;

    /* System generated locals */
    int i_1, i_2;
    float d_1;

    /* Local variables */
    static float pbar[20], dsum;
    static int i, j, l;
    static float p[420];	/* was [20][21] */
    static float x;
    static float y[21], z, pstar[20], ystar, p2star[20], dn, y2star;
    static int nn;
    static float curmin, durmin;
    static int jcount;
    static float del;
    static int ihi;
    static float dnn;
    static int ilo;
    static float ylo, sum;
    static double prevylo;
    static float npar[21];
    static float sigma_dist; /* square root of variance reqmin, plus scale */

    /* Parameter adjustments */
    --start;
    --mmin;
    --step;

    /* Function Body */

/*        ALGORITHM AS 47 APPLIED STATISTICS (J.R.STATIST.SOC C), */
/*        (1971) VOL.20, NO.3 by R. O'Neill */

/*        AS MODIFIED IN REMARK AS R28 by I. D. Hill */

/*      THE NELDER-MEAD SIMPLEX MINIMISATION PROCEDURE */


/*        PURPOSE :: TO FIND THE MINIMUM VALUE OF A USER-SPECIFIED */
/*                   POSITIVE FUNCTION. */


/*  FORMAL PARAMETERS :: */
/*          N : INPUT : THE NUMBER OF VARIABLES OVER WHICH WE ARE */
/*                    : MINIMISING */
/*      START : INPUT : ARRAY; CONTAINS THE COORDINATES OF THE */
/*                    : STARTING POINT. */
/*             OUTPUT : VALUES MAY BE OVERWRITTEN */
/*      MMIN : OUTPUT : ARRAY; CONTAINS THE COORDINATES OF THE */
/*                    : MINIMUM. */
/*     YNEWL : OUTPUT : THE MINIMUM VALUE OF THE FUNCTION. */
/*     REQMIN : INPUT : THE TERMINATING LIMIT FOR THE VARIANCE OF */
/*                    : FUNCTION VALUES. */
/*       STEP : INPUT : ARRAY; DETERMINES THE SIZE AND SHAPE OF THE */
/*                    : INITIAL SIMPLEX. THE RELATIVE MAGNITUDES OF */
/*                    : ITS N ELEMENTS SHOULD REFLECT THE UNITS OF */
/*                    : THE N VARIABLES. */
/*     KONVGE : INPUT : THE CONVERGENCE CHECK IS CARRIED OUT EVERY */
/*                    : KONVGE ITERATIONS. */
/*     KCOUNT : INPUT : MAXIMUM NUMBER OF FUNCTION EVALUATIONS. */
/*                    : If kcount==0, conditions causing a restart */
/*                    : cause the function to return with ifault=0, */
/*                    : rather than ifault=2 as occurs when icount>=kcount */
/*    ICOUNT : OUTPUT : ACTUAL NUMBER OF FUNCTION EVALUATIONS */
/*    NUMRES : OUTPUT : NUMBER OF RESTARTS */
/*    IFAULT : OUTPUT : 0  NO ERROR */
/*                    : 1  REQMIN, N OR KONVGE HAS ILLEGAL VALUE */
/*                    : 2  TERMINATION BECAUSE KCOUNT WAS EXCEEDED */
/*                    :    BEFORE CONVERGENCE WAS ACHIEVED */
/*     MAXRES : INPUT : MAXIMUM NUMBER OF RESTARTS */
/*        ALL VARIABLES AND ARRAYS ARE TO BE DECLARED IN THE CALLING */
/*        PROGRAM AS FLOATS. */


/*        AUXILIARY ALGORITHM :: THE REAL FUNCTION */
/*        SUBPROGRAM FN(A) CALCULATES THE FUNCTION VALUE AT POINT A. */
/*        A IS FLOAT WITH N ELEMENTS. */


/*        REFERENCE :: NELDER,J.A. AND MEAD, R.(1965). A SIMPLEX METHOD */

/*        FOR FUNCTION MINIMIZATION. COMPUTER J.,VOL.7,308-313 */

/* ***************************************************** */




/*        REFLECTION,EXTENSION AND CONTRACTION COEFFICIENTS. */

/*        VALIDITY CHECKS ON INPUT PARAMETERS. */

    *ifault = 1;
    if (*reqmin <= 0. || *n < 1 || *n > 20 || *konvge < 1) {
	return 0;
    }
    *ifault = 2;
    *icount = 0;
    *numres = 0;

    jcount = *konvge;
    dn = (float) (*n);
    nn = *n + 1;
    dnn = (float) nn;
    del = (float)1.;

    /* Construct population std deviation for use in T tests */
    sigma_dist= (float)sqrt((*reqmin)/((float)(*n)));

/*        CONSTRUCTION OF INITIAL SIMPLEX */

L1001:
    (*rst)(&start[1],userHook);
    i_1 = *n;
    for (i = 1; i <= i_1; ++i) {
/* L1: */
	p[i + nn * 20 - 21] = start[i];
    }
    z = prevylo = (*fn)(&start[1],userHook);
    y[nn - 1] = z;
    i_1 = *n;
    for (j = 1; j <= i_1; ++j) {
	x = start[j];
	start[j] += step[j] * del;
	i_2 = *n;
	for (i = 1; i <= i_2; ++i) {
/* L3: */
	    p[i + j * 20 - 21] = start[i];
	}
	z = (*fn)(&start[1],userHook);
	y[j - 1] = z;
/* L2: */
	start[j] = x;
    }
    *icount += nn;

/*        SIMPLEX CONSTRUCTION COMPLETE */

/*        FIND HIGHEST AND LOWEST Y VALUES. YNEWL ( =Y(IHI) ) INDICATES */

/*        THE VERTEX OF THE SIMPLEX TO BE REPLACED. */

L1000:
    ylo = y[0];
    *ynewl = ylo;
    ilo = 1;
    ihi = 1;
    i_1 = nn;
    for (i = 2; i <= i_1; ++i) {
	if (y[i - 1] >= ylo) {
	    goto L4;
	}
	ylo = y[i - 1];
	ilo = i;
L4:
	if (y[i - 1] <= *ynewl) {
	    goto L5;
	}
	*ynewl = y[i - 1];
	ihi = i;
L5:
    ;}

    if(ylo<prevylo){
        for(j = 1; j <= i_1; ++j)
	    npar[j - 1] = p[j + ilo * 20 - 21];
        (*rst)(npar,userHook);
	prevylo = ylo;
    }

/*        CALCULATE PBAR, THE CENTROID OF THE SIMPLEX VERTICES */
/*        EXCEPTING THAT WITH Y VALUE YNEWL. */

    i_1 = *n;
    for (i = 1; i <= i_1; ++i) {
	z = (float)0.;
	i_2 = nn;
	for (j = 1; j <= i_2; ++j) {
/* L6: */
	    z += p[i + j * 20 - 21];
	}
	z -= p[i + ihi * 20 - 21];
/* L7: */
	pbar[i - 1] = z / dn;
    }

/*        REFLECTION THROUGH THE CENTROID. */

    i_1 = *n;
    for (i = 1; i <= i_1; ++i) {
/* L8: */
	pstar[i - 1] = (rcoeff + (float)1.) * pbar[i - 1] - rcoeff * p[i + 
		ihi * 20 - 21];
    }
    ystar = (*fn)(pstar,userHook);
    ++(*icount);
    if (ystar >= ylo) {
	goto L12;
    }

/*        SUCCESSFUL REFLECTION, SO EXTENSION */

    i_1 = *n;
    for (i = 1; i <= i_1; ++i) {
/* L9: */
	p2star[i - 1] = ecoeff * pstar[i - 1] + ((float)1. - ecoeff) * pbar[i 
		- 1];
    }
    y2star = (*fn)(p2star,userHook);
    ++(*icount);

/*        RETAIN EXTENSION OR CONTRACTION. */

    if (y2star >= ystar) {
	goto L19;
    }
L10:
    i_1 = *n;
    for (i = 1; i <= i_1; ++i) {
/* L11: */
	p[i + ihi * 20 - 21] = p2star[i - 1];
    }
    y[ihi - 1] = y2star;
    goto L901;
/*        NO EXTENSION. */
L12:
    l = 0;
    i_1 = nn;
    for (i = 1; i <= i_1; ++i) {
	if (y[i - 1] > ystar) {
	    ++l;
	}
/* L13: */
    }
    if (l > 1) {
	goto L19;
    }
    if (l == 0) {
	goto L15;
    }

/*        CONTRACTION ON THE REFLECTION SIDE OF THE CENTROID. */

    i_1 = *n;
    for (i = 1; i <= i_1; ++i) {
/* L14: */
	p[i + ihi * 20 - 21] = pstar[i - 1];
    }
    y[ihi - 1] = ystar;

/*        CONTRACTION ON THE Y(IHI) SIDE OF THE CENTROID. */

L15:
    i_1 = *n;
    for (i = 1; i <= i_1; ++i) {
/* L16: */
	p2star[i - 1] = ccoeff * p[i + ihi * 20 - 21] + ((float)1. - ccoeff) *
		 pbar[i - 1];
    }
    y2star = (*fn)(p2star,userHook);
    ++(*icount);
    if (y2star <= y[ihi - 1]) {
	goto L10;
    }

/*        CONTRACT WHOLE SIMPLEX */

    i_1 = nn;
    for (j = 1; j <= i_1; ++j) {
	i_2 = *n;
	for (i = 1; i <= i_2; ++i) {
	    p[i + j * 20 - 21] = (p[i + j * 20 - 21] + p[i + ilo * 20 - 21]) *
		     (float).5;
/* L17: */
	    mmin[i] = p[i + j * 20 - 21];
	}
/* L18: */
	y[j - 1] = (*fn)(&mmin[1],userHook);
    }
    *icount += nn;
    goto L901;
/*        RETAIN REFLECTION */
L19:
    i_1 = *n;
    for (i = 1; i <= i_1; ++i) {
/* L20: */
	p[i + ihi * 20 - 21] = pstar[i - 1];
    }
    y[ihi - 1] = ystar;
L901:
    --jcount;
    if (jcount != 0) {
	goto L1000;
    }

/*        CHECK TO SEE IF MINIMUM REACHED. */

    if (*icount > *kcount) {
      if (y[ihi - 1] > y[ilo - 1]) {
	ihi = ilo;
      }
      i_1 = *n;
      for (i = 1; i <= i_1; ++i) {
	mmin[i] = p[i + ihi * 20 - 21];
      }
      *ynewl = y[ihi - 1];
      return 0;
    }
    jcount = *konvge;
    sum = (float)0.;
    i_1 = nn;
    for (i = 1; i <= i_1; ++i) {
/* L902: */
	sum += y[i - 1];
    }
    dsum = sum / dnn;
    durmin = (float)0.;
    i_1 = nn;
    for (i = 1; i <= i_1; ++i) {
/* L903: */
/* Computing 2nd power */
	d_1 = y[i - 1] - dsum;
	durmin += d_1 * d_1;
    }
    curmin = durmin / dn;

/*        CURMIN IS THE VARIANCE OF THE N+1 FN VALUES AT THE VERTICES. */

    if (curmin >= *reqmin) {
	goto L1000;
    }

/*        FACTORIAL TESTS TO CHECK THAT YNEWL IS A LOCAL MINIMUM */

L22:
    if (y[ihi - 1] > y[ilo - 1]) {
	ihi = ilo;
    }
    i_1 = *n;
    for (i = 1; i <= i_1; ++i) {
/* L23: */
	mmin[i] = p[i + ihi * 20 - 21];
    }
    *ynewl = y[ihi - 1];
    if (*icount > *kcount) {
	return 0;
    }
    i_1 = *n;
    for (i = 1; i <= i_1; ++i) {
	del = step[i] * (float).001;
	mmin[i] += del;
	z = (*fn)(&mmin[1],userHook);
	++(*icount);
	if ((z<*ynewl) && TTEST(z,dsum,sigma_dist,*tcritsqr)) {
	    goto L25;
	}
	mmin[i] = mmin[i] - del - del;
	z = (*fn)(&mmin[1],userHook);
	++(*icount);
	if ((z<*ynewl) && TTEST(z,dsum,sigma_dist,*tcritsqr)) {
	    goto L25;
	}
/* L24: */
	mmin[i] += del;
    }
    *ifault = 0;
    return 0;

/*        RESTART PROCEDURE */

L25:
    if (*numres >= *maxres) {
      if (*maxres==0) {
	*ifault= 0; /* The user wants to quit */
	*ynewl= (*fn)(&mmin[1],userHook);
      }
      return 0;
    }
    i_1 = *n;
    for (i = 1; i <= i_1; ++i) {
/* L26: */
	start[i] = mmin[i];
    }
    del = (float).001;
    ++(*numres);
    goto L1001;
} /* nelmin_ */

/* Subroutine */ 
int nelmin(float (*fn)(float*,void*), void (*rst)(float*,void*), int* n, 
	   float* start, float* mmin, float* ynewl, float* reqmin, 
	   float* step, int* konvge, int* kcount, int* icount, int* numres, 
	   int* ifault, int* maxres, void* userHook)
{
    /* Initialized data */

    static float rcoeff = 1.;
    static float ecoeff = 2.;
    static float ccoeff = .5;

    /* System generated locals */
    int i_1, i_2;
    float d_1;

    /* Local variables */
    static float pbar[20], dsum;
    static int i, j, l;
    static float p[420];	/* was [20][21] */
    static float x;
    static float y[21], z, pstar[20], ystar, p2star[20], dn, y2star;
    static int nn;
    static float curmin, durmin;
    static int jcount;
    static float del;
    static int ihi;
    static float dnn;
    static int ilo;
    static float ylo, sum;
    static double prevylo;
    static float npar[21];

    /* Parameter adjustments */
    --start;
    --mmin;
    --step;

    /* Function Body */

/*        ALGORITHM AS 47 APPLIED STATISTICS (J.R.STATIST.SOC C), */
/*        (1971) VOL.20, NO.3 by R. O'Neill */

/*        AS MODIFIED IN REMARK AS R28 by I. D. Hill */

/*      THE NELDER-MEAD SIMPLEX MINIMISATION PROCEDURE */


/*        PURPOSE :: TO FIND THE MINIMUM VALUE OF A USER-SPECIFIED */
/*                   POSITIVE FUNCTION. */


/*  FORMAL PARAMETERS :: */
/*          N : INPUT : THE NUMBER OF VARIABLES OVER WHICH WE ARE */
/*                    : MINIMISING */
/*      START : INPUT : ARRAY; CONTAINS THE COORDINATES OF THE */
/*                    : STARTING POINT. */
/*             OUTPUT : VALUES MAY BE OVERWRITTEN */
/*      MMIN : OUTPUT : ARRAY; CONTAINS THE COORDINATES OF THE */
/*                    : MINIMUM. */
/*     YNEWL : OUTPUT : THE MINIMUM VALUE OF THE FUNCTION. */
/*     REQMIN : INPUT : THE TERMINATING LIMIT FOR THE VARIANCE OF */
/*                    : FUNCTION VALUES. */
/*       STEP : INPUT : ARRAY; DETERMINES THE SIZE AND SHAPE OF THE */
/*                    : INITIAL SIMPLEX. THE RELATIVE MAGNITUDES OF */
/*                    : ITS N ELEMENTS SHOULD REFLECT THE UNITS OF */
/*                    : THE N VARIABLES. */
/*     KONVGE : INPUT : THE CONVERGENCE CHECK IS CARRIED OUT EVERY */
/*                    : KONVGE ITERATIONS. */
/*     KCOUNT : INPUT : MAXIMUM NUMBER OF FUNCTION EVALUATIONS. */
/*                    : If kcount==0, conditions causing a restart */
/*                    : cause the function to return with ifault=0, */
/*                    : rather than ifault=2 as occurs when icount>=kcount */
/*    ICOUNT : OUTPUT : ACTUAL NUMBER OF FUNCTION EVALUATIONS */
/*    NUMRES : OUTPUT : NUMBER OF RESTARTS */
/*    IFAULT : OUTPUT : 0  NO ERROR */
/*                    : 1  REQMIN, N OR KONVGE HAS ILLEGAL VALUE */
/*                    : 2  TERMINATION BECAUSE KCOUNT WAS EXCEEDED */
/*                    :    BEFORE CONVERGENCE WAS ACHIEVED */
/*     MAXRES : INPUT : MAXIMUM NUMBER OF RESTARTS */
/*        ALL VARIABLES AND ARRAYS ARE TO BE DECLARED IN THE CALLING */
/*        PROGRAM AS FLOATS. */


/*        AUXILIARY ALGORITHM :: THE REAL FUNCTION */
/*        SUBPROGRAM FN(A) CALCULATES THE FUNCTION VALUE AT POINT A. */
/*        A IS FLOAT WITH N ELEMENTS. */


/*        REFERENCE :: NELDER,J.A. AND MEAD, R.(1965). A SIMPLEX METHOD */

/*        FOR FUNCTION MINIMIZATION. COMPUTER J.,VOL.7,308-313 */

/* ***************************************************** */




/*        REFLECTION,EXTENSION AND CONTRACTION COEFFICIENTS. */

/*        VALIDITY CHECKS ON INPUT PARAMETERS. */

    *ifault = 1;
    if (*reqmin <= 0. || *n < 1 || *n > 20 || *konvge < 1) {
	return 0;
    }
    *ifault = 2;
    *icount = 0;
    *numres = 0;

    jcount = *konvge;
    dn = (float) (*n);
    nn = *n + 1;
    dnn = (float) nn;
    del = (float)1.;

/*        CONSTRUCTION OF INITIAL SIMPLEX */

L1001:
    (*rst)(&start[1],userHook);
    i_1 = *n;
    for (i = 1; i <= i_1; ++i) {
/* L1: */
	p[i + nn * 20 - 21] = start[i];
    }
    z = prevylo = (*fn)(&start[1],userHook);
    y[nn - 1] = z;
    i_1 = *n;
    for (j = 1; j <= i_1; ++j) {
	x = start[j];
	start[j] += step[j] * del;
	i_2 = *n;
	for (i = 1; i <= i_2; ++i) {
/* L3: */
	    p[i + j * 20 - 21] = start[i];
	}
	z = (*fn)(&start[1],userHook);
	y[j - 1] = z;
/* L2: */
	start[j] = x;
    }
    *icount += nn;

/*        SIMPLEX CONSTRUCTION COMPLETE */

/*        FIND HIGHEST AND LOWEST Y VALUES. YNEWL ( =Y(IHI) ) INDICATES */

/*        THE VERTEX OF THE SIMPLEX TO BE REPLACED. */

L1000:
    ylo = y[0];
    *ynewl = ylo;
    ilo = 1;
    ihi = 1;
    i_1 = nn;
    for (i = 2; i <= i_1; ++i) {
	if (y[i - 1] >= ylo) {
	    goto L4;
	}
	ylo = y[i - 1];
	ilo = i;
L4:
	if (y[i - 1] <= *ynewl) {
	    goto L5;
	}
	*ynewl = y[i - 1];
	ihi = i;
L5:
    ;}

    if(ylo<prevylo){
        for(j = 1; j <= i_1; ++j)
	    npar[j - 1] = p[j + ilo * 20 - 21];
        (*rst)(npar,userHook);
	prevylo = ylo;
    }

/*        CALCULATE PBAR, THE CENTROID OF THE SIMPLEX VERTICES */
/*        EXCEPTING THAT WITH Y VALUE YNEWL. */

    i_1 = *n;
    for (i = 1; i <= i_1; ++i) {
	z = (float)0.;
	i_2 = nn;
	for (j = 1; j <= i_2; ++j) {
/* L6: */
	    z += p[i + j * 20 - 21];
	}
	z -= p[i + ihi * 20 - 21];
/* L7: */
	pbar[i - 1] = z / dn;
    }

/*        REFLECTION THROUGH THE CENTROID. */

    i_1 = *n;
    for (i = 1; i <= i_1; ++i) {
/* L8: */
	pstar[i - 1] = (rcoeff + (float)1.) * pbar[i - 1] - rcoeff * p[i + 
		ihi * 20 - 21];
    }
    ystar = (*fn)(pstar,userHook);
    ++(*icount);
    if (ystar >= ylo) {
	goto L12;
    }

/*        SUCCESSFUL REFLECTION, SO EXTENSION */

    i_1 = *n;
    for (i = 1; i <= i_1; ++i) {
/* L9: */
	p2star[i - 1] = ecoeff * pstar[i - 1] + ((float)1. - ecoeff) * pbar[i 
		- 1];
    }
    y2star = (*fn)(p2star,userHook);
    ++(*icount);

/*        RETAIN EXTENSION OR CONTRACTION. */

    if (y2star >= ystar) {
	goto L19;
    }
L10:
    i_1 = *n;
    for (i = 1; i <= i_1; ++i) {
/* L11: */
	p[i + ihi * 20 - 21] = p2star[i - 1];
    }
    y[ihi - 1] = y2star;
    goto L901;
/*        NO EXTENSION. */
L12:
    l = 0;
    i_1 = nn;
    for (i = 1; i <= i_1; ++i) {
	if (y[i - 1] > ystar) {
	    ++l;
	}
/* L13: */
    }
    if (l > 1) {
	goto L19;
    }
    if (l == 0) {
	goto L15;
    }

/*        CONTRACTION ON THE REFLECTION SIDE OF THE CENTROID. */

    i_1 = *n;
    for (i = 1; i <= i_1; ++i) {
/* L14: */
	p[i + ihi * 20 - 21] = pstar[i - 1];
    }
    y[ihi - 1] = ystar;

/*        CONTRACTION ON THE Y(IHI) SIDE OF THE CENTROID. */

L15:
    i_1 = *n;
    for (i = 1; i <= i_1; ++i) {
/* L16: */
	p2star[i - 1] = ccoeff * p[i + ihi * 20 - 21] + ((float)1. - ccoeff) *
		 pbar[i - 1];
    }
    y2star = (*fn)(p2star,userHook);
    ++(*icount);
    if (y2star <= y[ihi - 1]) {
	goto L10;
    }

/*        CONTRACT WHOLE SIMPLEX */

    i_1 = nn;
    for (j = 1; j <= i_1; ++j) {
	i_2 = *n;
	for (i = 1; i <= i_2; ++i) {
	    p[i + j * 20 - 21] = (p[i + j * 20 - 21] + p[i + ilo * 20 - 21]) *
		     (float).5;
/* L17: */
	    mmin[i] = p[i + j * 20 - 21];
	}
/* L18: */
	y[j - 1] = (*fn)(&mmin[1],userHook);
    }
    *icount += nn;
    goto L901;
/*        RETAIN REFLECTION */
L19:
    i_1 = *n;
    for (i = 1; i <= i_1; ++i) {
/* L20: */
	p[i + ihi * 20 - 21] = pstar[i - 1];
    }
    y[ihi - 1] = ystar;
L901:
    --jcount;
    if (jcount != 0) {
	goto L1000;
    }

/*        CHECK TO SEE IF MINIMUM REACHED. */

    if (*icount > *kcount) {
	goto L22;
    }
    jcount = *konvge;
    sum = (float)0.;
    i_1 = nn;
    for (i = 1; i <= i_1; ++i) {
/* L902: */
	sum += y[i - 1];
    }
    dsum = sum / dnn;
    durmin = (float)0.;
    i_1 = nn;
    for (i = 1; i <= i_1; ++i) {
/* L903: */
/* Computing 2nd power */
	d_1 = y[i - 1] - dsum;
	durmin += d_1 * d_1;
    }
    curmin = durmin / dn;

/*        CURMIN IS THE VARIANCE OF THE N+1 FN VALUES AT THE VERTICES. */

    if (curmin >= *reqmin) {
	goto L1000;
    }

/*        FACTORIAL TESTS TO CHECK THAT YNEWL IS A LOCAL MINIMUM */

L22:
    if (y[ihi - 1] > y[ilo - 1]) {
	ihi = ilo;
    }
    i_1 = *n;
    for (i = 1; i <= i_1; ++i) {
/* L23: */
	mmin[i] = p[i + ihi * 20 - 21];
    }
    *ynewl = y[ihi - 1];
    if (*icount > *kcount) {
	return 0;
    }
    i_1 = *n;
    for (i = 1; i <= i_1; ++i) {
	del = step[i] * (float).001;
	mmin[i] += del;
	z = (*fn)(&mmin[1],userHook);
	++(*icount);
	if (z < *ynewl) {
	    goto L25;
	}
	mmin[i] = mmin[i] - del - del;
	z = (*fn)(&mmin[1],userHook);
	++(*icount);
	if (z < *ynewl) {
	    goto L25;
	}
/* L24: */
	mmin[i] += del;
    }
    *ifault = 0;
    return 0;

/*        RESTART PROCEDURE */

L25:
    if (*numres >= *maxres) {
      if (*maxres==0) {
	*ifault= 0; /* The user wants to quit */
	*ynewl= (*fn)(&mmin[1],userHook);
      }
      return 0;
    }
    i_1 = *n;
    for (i = 1; i <= i_1; ++i) {
/* L26: */
	start[i] = mmin[i];
    }
    del = (float).001;
    ++(*numres);
    goto L1001;
} /* nelmin_ */

