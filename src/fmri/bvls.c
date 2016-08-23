/* bvls.f -- translated by f2c (version 19970805).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/* Table of constant values */

static long c__9 = 9;
static long c__1 = 1;
static double c_b38 = 1.;

#define MIN( a, b ) ( (a) < (b) ? (a) : (b) )
#define MAX( a, b ) ( (a) > (b) ? (a) : (b) )

static int qr(long* m, long* n, double* a, double* b, 
	       double* x, double* resq);

/* Summary: */
/*BVLS solves linear least-squares problems with upper and lower bounds on the
*/
/* variables, using an active set strategy.  It is documented in the J. of */
/* Computational Statistics, and can be used iteratively to solve minimum */
/* l-1, l-2 and l-infinity fitting problems. */
/* Statement: */
/*This algorithm may be used freely for non-commercial purposes, and may be*/
/*freely distributed for non-commercial purposes.  The authors do not warrant
*/
/* the software in any way: use it at your own risk. */
/* ======================================================================= */
/* Subroutine */ int bvls(long* key, long* m, long* n, 
			   double* a, double* b, double* bl, double* bu, 
			   double* x, double* w, double* act, double* zz, 
			   long* istate, long* loopa)
{
    /* Initialized data */

    double eps = 1e-11;

    /* System generated locals */
    long a_dim1, a_offset, act_dim1, act_offset, i__1, i__2, i__3;
    double d__1, d__2;

    /* Local variables */
    long iact, nact;
    double resq;
    long i__;
    double bdiff;
    long j, k;
    double alpha;
    long noldb;
    double bound, bnorm;
    long k1;
    double worst;
    long ifrom5, jj, kk, mm;
    double ri;
    long it, ks;
    double sj;
    long nbound, mm1;
    double bad, alf, obj, bsq;


/* =======================================================================
 */

/* $$$$ calls qr */
/* --------------------Bounded Variable Least Squares---------------------
 */

/*        Robert L. Parker and Philip B. Stark    Version 3/19/90 */

/*  Robert L. Parker                           Philip B. Stark */
/*  Scripps Institution of Oceanography        Department of Statistics */
/*  University of California, San Diego        University of California */
/*  La Jolla CA 92093                          Berkeley CA 94720-3860 */
/*  rlparker@ucsd.edu                          stark@stat.berkeley.edu */

/*  Copyright of this software is reserved by the authors; however, this 
*/
/*  algorithm and subroutine may be used freely for non-commercial */
/*  purposes, and may be distributed freely for non-commercial purposes. 
*/

/*  The authors do not warrant this software in any way: use it at your */
/*  own risk. */


/*  See the article ``Bounded Variable Least Squares:  An Algorithm and */
/*  Applications'' by P.B. Stark and R.L. Parker, in the journal */
/*  Computational Statistics, in press (1995) for further description */
/*  and applications to minimum l-1, l-2 and l-infinity fitting problems, 
*/
/*  as well as finding bounds on linear functionals subject to bounds on 
*/
/*  variables and fitting linear data within l-1, l-2 or l-infinity */
/*  measures of misfit. */

/*  BVLS solves the problem: */

/*          min  || a.x - b ||     such that   bl <= x <= bu */
/*                            2 */
/*    where */
/*               x  is an unknown n-vector */
/*               a  is a given m by n matrix */
/*               b  is a given  m-vector */
/*               bl is a given n-vector of lower bounds on the */
/*                                components of x. */
/*               bu is a given n-vector of upper bounds on the */
/*                                components of x. */


/* -----------------------------------------------------------------------
 */
/*    Input parameters: */

/*  m, n, a, b, bl, bu   see above.   Let mm=min(m,n). */

/*  If key = 0, the subroutine solves the problem from scratch. */

/*  If key > 0 the routine initializes using the user's guess about */
/*   which components of  x  are `active', i.e. are stricly within their 
*/
/*   bounds, which are at their lower bounds, and which are at their */
/*   upper bounds.  This information is supplied through the array */
/*   istate.  istate(n+1) should contain the total number of components */
/*   at their bounds (the `bound variables').  The absolute values of the 
*/
/*   first nbound=istate(n+1) entries of  istate  are the indices */
/*   of these `bound' components of  x.  The sign of istate(j), j=1,..., 
*/
/*   nbound, indicates whether  x(|istate(j)|) is at its upper or lower */
/*   bound.  istate(j) is positive if the component is at its upper */
/*   bound, negative if the component is at its lower bound. */
/*   istate(j), j=nbound+1,...,n  contain the indices of the components */
/*   of  x  that are active (i.e. are expected to lie strictly within */
/*   their bounds).  When key > 0, the routine initially sets the active 
*/
/*   components to the averages of their upper and lower bounds: */
/*   x(j)=(bl(j)+bu(j))/2, for j in the active set. */

/* -----------------------------------------------------------------------
 */
/*    Output parameters: */

/*  x       the solution vector. */

/*  w(1)    the minimum 2-norm || a.x-b ||. */

/*  istate  vector indicating which components of  x  are active and */
/*          which are at their bounds (see the previous paragraph). */
/*          istate can be supplied to the routine to give it a good */
/*          starting guess for the solution. */

/*  loopA   number of iterations taken in the main loop, Loop A. */

/* -----------------------------------------------------------------------
 */
/*    Working  arrays: */

/*  w      dimension n.               act      dimension m*(mm+2). */
/*  zz     dimension m.               istate   dimension n+1. */

/* -----------------------------------------------------------------------
 */
/*  Method: active variable method along the general plan of NNLS by */
/*  Lawson & Hanson, "Solving Least Squares Problems," 1974.  See */
/*  Algorithm 23.10.  Step numbers in comment statements refer to their */
/*  scheme. */
/*  For more details and further uses, see the article */
/*  "Bounded Variable Least Squares:  An Algorithm and Applications" */
/*  by Stark and Parker in 1995 Computational Statistics. */

/* -----------------------------------------------------------------------
 */
/*  A number of measures are taken to enhance numerical reliability: */

/* 1. As noted by Lawson and Hanson, roundoff errors in the computation */
/*   of the gradient of the misfit may cause a component on the bounds */
/*   to appear to want to become active, yet when the component is added 
*/
/*   to the active set, it moves away from the feasible region.  In this 
*/
/*   case the component is not made active, the gradient of the misfit */
/*   with respect to a change in that component is set to zero, and the */
/*   program returns to the Kuhn-Tucker test.  Flag  ifrom5  is used in */
/*   this test, which occurs at the end of Step 6. */


/* 2. When the least-squares minimizer after Step 6 is infeasible, it */
/*   is used in a convex interpolation with the previous solution to */
/*   obtain a feasible vector.  The constant in this interpolation is */
/*   supposed to put at least one component of  x   on a bound. There can 
*/
/*   be difficulties: */

/* 2a. Sometimes, due to roundoff, no interpolated component ends up on */
/*   a bound.  The code in Step 11 uses the flag  jj, computed in Step 8, 
*/
/*   to ensure that at least the component that determined the */
/*   interpolation constant  alpha  is moved to the appropriate bound. */
/*   This guarantees that what Lawson and Hanson call `Loop B' is finite. 
*/

/* 2b. The code in Step 11 also incorporates Lawson and Hanson's feature 
*/
/*   that any components remaining infeasible at this stage (which must */
/*   be due to roundoff) are moved to their nearer bound. */


/* 3. If the columns of  a  passed to qr are linearly dependent, the new 
*/
/*   potentially active component is not introduced: the gradient of the 
*/
/*   misfit with respect to that component is set to zero, and control */
/*   returns to the Kuhn-Tucker test. */


/* 4. When some of the columns of  a  are approximately linearly */
/*   dependent, we have observed cycling of active components: a */
/*   component just moved to a bound desires immediately to become */
/*   active again; qr allows it to become active and a different */
/*   component is moved to its bound.   This component immediately wants 
*/
/*   to become active, which qr allows, and the original component is */
/*   moved back to its bound.  We have taken two steps to avoid this */
/*   problem: */

/* 4a. First, the column of the matrix  a  corresponding to the new */
/*   potentially active component is passed to qr as the last column of */
/*   its matrix.  This ordering tends to make a component recently moved 
*/
/*   to a bound fail the test mentioned in (1), above. */

/* 4b. Second, we have incorporated a test that prohibits short cycles. */
/*   If the most recent successful change to the active set was to move */
/*   the component x(jj) to a bound, x(jj) is not permitted to reenter */
/*   the solution at this stage.  This test occurs just after checking */
/*   the Kuhn-Tucker conditions, and uses the flag  jj, set in Step 8. */
/*   The flag  jj  is reset after Step 6 if Step 6 was entered from */
/*   Step 5 indicating that a new component has successfully entered the 
*/
/*   active set. The test for resetting  jj  uses the flag  ifrom5, */
/*   which will not equal zero in case Step 6 was entered from Step 5. */


/*     dimension w(n), act(m,min(m,n)+2), zz(m), istate(n+1) */

    /* Parameter adjustments */
    --zz;
    act_dim1 = *m;
    act_offset = act_dim1 + 1;
    act -= act_offset;
    --b;
    --istate;
    --w;
    --x;
    --bu;
    --bl;
    a_dim1 = *m;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */

/* ----------------------First Executable Statement-----------------------
 */

/*  Step 1.  Initialize everything--active and bound sets, initial */
/*   values, etc. */

/*  Initialize flags, etc. */
    mm = MIN(*m,*n);
    mm1 = mm + 1;
    jj = 0;
    ifrom5 = 0;
/*  Check consistency of given bounds  bl, bu. */
    bdiff = (float)0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	d__1 = bdiff, d__2 = bu[j] - bl[j];
	bdiff = MAX(d__1,d__2);
	if (bl[j] > bu[j]) {
	  fprintf(stderr,"Inconsistent bounds in BVLS.\n");
	  exit(-1);
	}
/* L1005: */
    }
    if (bdiff == (float)0.) {
      fprintf(stderr,"No free variables in BVLS--check input bounds.\n");
      exit(-1);
    }

/*  In a fresh initialization (key = 0) bind all variables at their lower 
*/
/*   bounds.  If (key != 0), use the supplied  istate  vector to */
/*   initialize the variables.  istate(n+1) contains the number of */
/*   bound variables.  The absolute values of the first */
/*   nbound=istate(n+1) entries of  istate  are the indices of the bound 
*/
/*   variables.  The sign of each entry determines whether the indicated 
*/
/*   variable is at its upper (positive) or lower (negative) bound. */
    if (*key == 0) {
	nbound = *n;
	nact = 0;
	i__1 = nbound;
	for (j = 1; j <= i__1; ++j) {
	    istate[j] = -j;
/* L1010: */
	}
    } else {
	nbound = istate[*n + 1];
    }
    nact = *n - nbound;
    if (nact > mm) {
      fprintf(stderr,"Too many active variables in BVLS starting solution!\n");
      exit(-1);
    }
    i__1 = nbound;
    for (k = 1; k <= i__1; ++k) {
	j = (i__2 = istate[k], abs(i__2));
	if (istate[k] < 0) {
	    x[j] = bl[j];
	}
	if (istate[k] > 0) {
	    x[j] = bu[j];
	}
/* L1100: */
    }

/*  In a warm start (key != 0) initialize the active variables to */
/*   (bl+bu)/2.  This is needed in case the initial qr results in */
/*   active variables out-of-bounds and Steps 8-11 get executed the */
/*   first time through. */
    i__1 = *n;
    for (k = nbound + 1; k <= i__1; ++k) {
	kk = istate[k];
	x[kk] = (bu[kk] + bl[kk]) / 2;
/* L1150: */
    }

/*  Compute bnorm, the norm of the data vector b, for reference. */
    bsq = (float)0.;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = b[i__];
	bsq += d__1 * d__1;
/* L1200: */
    }
    bnorm = sqrt(bsq);

/* -----------------------------Main Loop---------------------------------
 */

/*  Initialization complete.  Begin major loop (Loop A). */
    i__1 = *n * 3;
    for (*loopa = 1; *loopa <= i__1; ++(*loopa)) {

/*  Step 2. */
/*  Initialize the negative gradient vector w(*). */
/* L2000: */
	obj = (float)0.;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    w[j] = (float)0.;
/* L2050: */
	}

/*  Compute the residual vector b-a.x , the negative gradient vector 
*/
/*   w(*), and the current objective value obj = || a.x - b ||. */
/*   The residual vector is stored in the mm+1'st column of act(*,*). 
*/
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ri = b[i__];
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		ri -= a[i__ + j * a_dim1] * x[j];
/* L2100: */
	    }
/* Computing 2nd power */
	    d__1 = ri;
	    obj += d__1 * d__1;
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		w[j] += a[i__ + j * a_dim1] * ri;
/* L2200: */
	    }
	    act[i__ + mm1 * act_dim1] = ri;
/* L2300: */
	}

/*  Converged?  Stop if the misfit << || b ||, or if all components ar
e */
/*   active (unless this is the first iteration from a warm start). */
	if (sqrt(obj) <= bnorm * eps || *loopa > 1 && nbound == 0) {
	    istate[*n + 1] = nbound;
	    w[1] = sqrt(obj);
	    return 0;
	}

/*  Add the contribution of the active components back into the residu
al. */
	i__2 = *n;
	for (k = nbound + 1; k <= i__2; ++k) {
	    j = istate[k];
	    i__3 = *m;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		act[i__ + mm1 * act_dim1] += a[i__ + j * a_dim1] * x[j];
/* L2400: */
	    }
/* L2500: */
	}

/*  The first iteration in a warm start requires immediate qr. */
	if (*loopa == 1 && *key != 0) {
	    goto L6000;
	}

/*  Steps 3, 4. */
/*  Find the bound element that most wants to be active. */
L3000:
	worst = (float)0.;
	it = 1;
	i__2 = nbound;
	for (j = 1; j <= i__2; ++j) {
	    ks = (i__3 = istate[j], abs(i__3));
	    /* bad = w[ks] * i_sign(&c__1, &istate[j]); */
	    if (istate[j]<0) bad = -w[ks];
	    else bad = w[ks];
	    if (bad < worst) {
		it = j;
		worst = bad;
		iact = ks;
	    }
/* L3100: */
	}

/*  Test whether the Kuhn-Tucker condition is met. */
	if (worst >= (float)0.) {
	    istate[*n + 1] = nbound;
	    w[1] = sqrt(obj);
	    return 0;
	}

/*  The component  x(iact)  is the one that most wants to become activ
e. */
/*   If the last successful change in the active set was to move x(iac
t) */
/*   to a bound, don't let x(iact) in now: set the derivative of the 
*/
/*   misfit with respect to x(iact) to zero and return to the Kuhn-Tuc
ker */
/*   test. */
	if (iact == jj) {
	    w[jj] = (float)0.;
	    goto L3000;
	}

/*  Step 5. */
/*  Undo the effect of the new (potentially) active variable on the */
/*   residual vector. */
	if (istate[it] > 0) {
	    bound = bu[iact];
	}
	if (istate[it] < 0) {
	    bound = bl[iact];
	}
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    act[i__ + mm1 * act_dim1] += bound * a[i__ + iact * a_dim1];
/* L5100: */
	}

/*  Set flag ifrom5, indicating that Step 6 was entered from Step 5. 
*/
/*   This forms the basis of a test for instability: the gradient */
/*   calculation shows that x(iact) wants to join the active set; if 
*/
/*   qr puts x(iact) beyond the bound from which it came, the gradient
 */
/*   calculation was in error and the variable should not have been */
/*   introduced. */
	ifrom5 = istate[it];

/*  Swap the indices (in istate) of the new active variable and the */
/*   rightmost bound variable; `unbind' that location by decrementing 
*/
/*   nbound. */
	istate[it] = istate[nbound];
	--nbound;
	++nact;
	istate[nbound + 1] = iact;

	if (mm < nact) {
	  fprintf(stderr,"Too many free variables in BVLS!\n");
	  exit(-1);
	}

/*  Step 6. */
/*  Load array  act  with the appropriate columns of  a  for qr.  For 
*/
/*   added stability, reverse the column ordering so that the most */
/*   recent addition to the active set is in the last column.  Also */
/*   copy the residual vector from act(., mm1) into act(., mm1+1). */
L6000:
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    act[i__ + (mm1 + 1) * act_dim1] = act[i__ + mm1 * act_dim1];
	    i__3 = *n;
	    for (k = nbound + 1; k <= i__3; ++k) {
		j = istate[k];
		act[i__ + (nact + 1 - k + nbound) * act_dim1] = a[i__ + j * 
			a_dim1];
/* L6100: */
	    }
/* L6200: */
	}

	qr(m, &nact, &act[act_offset], &act[(mm1 + 1) * act_dim1 + 1], &zz[1]
		, &resq);

/*  Test for linear dependence in qr, and for an instability that move
s */
/*   the variable just introduced away from the feasible region */
/*   (rather than into the region or all the way through it). */
/*   In either case, remove the latest vector introduced from the */
/*   active set and adjust the residual vector accordingly. */
/*   Set the gradient component (w(iact)) to zero and return to */
/*   the Kuhn-Tucker test. */
	if (resq < (float)0. || ifrom5 > 0 && zz[nact] > bu[iact] || ifrom5 < 
		0 && zz[nact] < bl[iact]) {
	    ++nbound;
	    d__1 = x[iact] - bu[iact];
	    /*istate[nbound] = (long) (istate[nbound] * d_sign(&c_b38, &d__1));*/
	    istate[nbound] = (long) (istate[nbound] * copysign(c_b38, d__1)
		    );
	    --nact;
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		act[i__ + mm1 * act_dim1] -= x[iact] * a[i__ + iact * a_dim1];
/* L6500: */
	    }
	    ifrom5 = 0;
	    w[iact] = (float)0.;
	    goto L3000;
	}

/*  If Step 6 was entered from Step 5 and we are here, a new variable 
*/
/*   has been successfully introduced into the active set; the last */
/*   variable that was fixed at a bound is again permitted to become 
*/
/*   active. */
	if (ifrom5 != 0) {
	    jj = 0;
	}
	ifrom5 = 0;

/*   Step 7.  Check for strict feasibility of the new qr solution. */
	i__2 = nact;
	for (k = 1; k <= i__2; ++k) {
	    k1 = k;
	    j = istate[k + nbound];
	    if (zz[nact + 1 - k] < bl[j] || zz[nact + 1 - k] > bu[j]) {
		goto L8000;
	    }
/* L7100: */
	}
	i__2 = nact;
	for (k = 1; k <= i__2; ++k) {
	    j = istate[k + nbound];
	    x[j] = zz[nact + 1 - k];
/* L7200: */
	}
/*  New iterate is feasible; back to the top. */
	goto L15000;

/*  Steps 8, 9. */
L8000:
	alpha = (float)2.;
	alf = alpha;
	i__2 = nact;
	for (k = k1; k <= i__2; ++k) {
	    j = istate[k + nbound];
	    if (zz[nact + 1 - k] > bu[j]) {
		alf = (bu[j] - x[j]) / (zz[nact + 1 - k] - x[j]);
	    }
	    if (zz[nact + 1 - k] < bl[j]) {
		alf = (bl[j] - x[j]) / (zz[nact + 1 - k] - x[j]);
	    }
	    if (alf < alpha) {
		alpha = alf;
		jj = j;
		d__1 = zz[nact + 1 - k] - bl[j];
		/*sj = d_sign(&c_b38, &d__1);*/
		sj = copysign(c_b38, d__1);
	    }
/* L8200: */
	}

/*  Step 10 */
	i__2 = nact;
	for (k = 1; k <= i__2; ++k) {
	    j = istate[k + nbound];
	    x[j] += alpha * (zz[nact + 1 - k] - x[j]);
/* L10000: */
	}

/*  Step 11. */
/*  Move the variable that determined alpha to the appropriate bound. 
*/
/*   (jj is its index; sj is + if zz(jj)> bu(jj), - if zz(jj)<bl(jj) )
. */
/*   If any other component of  x  is infeasible at this stage, it mus
t */
/*   be due to roundoff.  Bind every infeasible component and every */
/*   component at a bound to the appropriate bound.  Correct the */
/*   residual vector for any variables moved to bounds.  Since at leas
t */
/*   one variable is removed from the active set in this step, Loop B 
*/
/*   (Steps 6-11) terminates after at most  nact  steps. */
	noldb = nbound;
	i__2 = nact;
	for (k = 1; k <= i__2; ++k) {
	    j = istate[k + noldb];
	    if (bu[j] - x[j] <= (float)0. || j == jj && sj > (float)0.) {
/*  Move x(j) to its upper bound. */
		x[j] = bu[j];
		istate[k + noldb] = istate[nbound + 1];
		istate[nbound + 1] = j;
		++nbound;
		i__3 = *m;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    act[i__ + mm1 * act_dim1] -= bu[j] * a[i__ + j * a_dim1];
/* L11100: */
		}
	    } else if (x[j] - bl[j] <= (float)0. || j == jj && sj < (float)0.)
		     {
/*  Move x(j) to its lower bound. */
		x[j] = bl[j];
		istate[k + noldb] = istate[nbound + 1];
		istate[nbound + 1] = -j;
		++nbound;
		i__3 = *m;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    act[i__ + mm1 * act_dim1] -= bl[j] * a[i__ + j * a_dim1];
/* L11150: */
		}
	    }
/* L11200: */
	}
	nact = *n - nbound;

/*  If there are still active variables left repeat the qr; if not, */
/*    go back to the top. */
	if (nact > 0) {
	    goto L6000;
	}

L15000:
	;
    }

    fprintf(stderr,"BVLS fails to converge!\n");
    exit(-1);
    return 1;
} /* bvls */

/* ====================================================================== */
/* Subroutine */ static int qr(long* m, long* n, double* a, double* b, 
				double* x, double* resq)
{
    /* System generated locals */
    long a_dim1, a_offset, i__1, i__2, i__3;
    double d__1;

    /* Local variables */
    long i__, j;
    double const__;
    long j1;
    double u1;
    long ii, jj;
    double sq, qv1, dot, sum;

/* ====================================================================== 
*/
/* $$$$ calls no other routines */
/*  Relies on FORTRAN77 do-loop conventions! */
/*  Solves over-determined least-squares problem  ax ~ b */
/*  where  a  is an  m by n  matrix,  b  is an m-vector . */
/*  resq  is the sum of squared residuals of optimal solution.  Also used 
*/
/*  to signal error conditions - if -2 , system is underdetermined,  if */
/*  -1,  system is singular. */
/*  Method - successive Householder rotations.  See Lawson & Hanson - */
/*  Solving Least Squares Problems (1974). */
/*  Routine will also work when m=n. */
/* *****   CAUTION -  a and b  are overwritten by this routine. */

    /* Parameter adjustments */
    --b;
    --x;
    a_dim1 = *m;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    *resq = (float)-2.;
    if (*m < *n) {
	return 0;
    }
    *resq = (float)-1.;
/*   Loop ending on 1800 rotates  a  into upper triangular form. */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*   Find constants for rotation and diagonal entry. */
	sq = (float)0.;
	i__2 = *m;
	for (i__ = j; i__ <= i__2; ++i__) {
/* Computing 2nd power */
	    d__1 = a[i__ + j * a_dim1];
	    sq = d__1 * d__1 + sq;
/* L1100: */
	}
	if (sq == (float)0.) {
	    return 0;
	}
	d__1 = sqrt(sq);
	/*qv1 = -d_sign(&d__1, &a[j + j * a_dim1]);*/
	qv1 = -copysign(d__1, a[j + j * a_dim1]);
	u1 = a[j + j * a_dim1] - qv1;
	a[j + j * a_dim1] = qv1;
	j1 = j + 1;
/*  Rotate remaining columns of sub-matrix. */
	i__2 = *n;
	for (jj = j1; jj <= i__2; ++jj) {
	    dot = u1 * a[j + jj * a_dim1];
	    i__3 = *m;
	    for (i__ = j1; i__ <= i__3; ++i__) {
		dot = a[i__ + jj * a_dim1] * a[i__ + j * a_dim1] + dot;
/* L1200: */
	    }
	    const__ = dot / (d__1 = qv1 * u1, abs(d__1));
	    i__3 = *m;
	    for (i__ = j1; i__ <= i__3; ++i__) {
		a[i__ + jj * a_dim1] -= const__ * a[i__ + j * a_dim1];
/* L1300: */
	    }
	    a[j + jj * a_dim1] -= const__ * u1;
/* L1400: */
	}
/*  Rotate  b  vector. */
	dot = u1 * b[j];
	i__2 = *m;
	for (i__ = j1; i__ <= i__2; ++i__) {
	    dot = b[i__] * a[i__ + j * a_dim1] + dot;
/* L1600: */
	}
	const__ = dot / (d__1 = qv1 * u1, abs(d__1));
	b[j] -= const__ * u1;
	i__2 = *m;
	for (i__ = j1; i__ <= i__2; ++i__) {
	    b[i__] -= const__ * a[i__ + j * a_dim1];
/* L1700: */
	}
/* L1800: */
    }
/*  Solve triangular system by back-substitution. */
    i__1 = *n;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = *n - ii + 1;
	sum = b[i__];
	i__2 = *n;
	for (j = i__ + 1; j <= i__2; ++j) {
	    sum -= a[i__ + j * a_dim1] * x[j];
/* L2100: */
	}
	if (a[i__ + i__ * a_dim1] == (float)0.) {
	    return 0;
	}
	x[i__] = sum / a[i__ + i__ * a_dim1];
/* L2200: */
    }
/*  Find residual in overdetermined case. */
    *resq = (float)0.;
    i__1 = *m;
    for (i__ = *n + 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = b[i__];
	*resq = d__1 * d__1 + *resq;
/* L2300: */
    }
    return 0;
} /* qr */

