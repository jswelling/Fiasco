#include <math.h>
#include <stdio.h>
#include "fmri.h"

/* The squared inverse of the golden ratio */
#define INV_GOLD ((3. - sqrt(5.)) * .5)

double fmin1D(double ax, double bx, double (*func)(double), double tol, double eps, int debug)
{
  /* Local variables */
  double a, b, d__, e, p, q, r__, u, v, w, x;
  double t2, fu, fv, fw, fx, xm, tol1, tol3;
  int niter= 0;

/*  an approximation  x  to the point where  f  attains a minimum  on
 *  the interval  (ax,bx)  is determined.
 *
 *  input..
 *
 *  ax    left endpoint of initial interval
 *  bx    right endpoint of initial interval
 *  f     function subprogram which evaluates  f(x)  for any  x
 *        in the interval  (ax,bx)
 *  tol   desired length of the interval of uncertainty of the final
 *        result (.ge.0.)
 *  eps   Tolerance; typically sqrt of machine precision
 *  debug requests debugging output
 *
 *  output..
 *
 *  fmin  abcissa approximating the point where  f  attains a
 *        minimum
 *
 *      the method used is a combination of  golden  section  search  and
 *  successive parabolic interpolation.  convergence is never much slower
 *  than  that  for  a  fibonacci search.  if  f  has a continuous second
 *  derivative which is positive at the minimum (which is not  at  ax  or
 *  bx),  then  convergence  is  superlinear, and usually of the order of
 *  about  1.324....
 *      the function  f  is never evaluated at two points closer together
 *  than  eps*fabs(fmin)+(tol/3), where eps is  approximately  the  square
 *  root  of  the  relative  machine  precision.   if   f   is a unimodal
 *  function and the computed values of   f   are  always  unimodal  when
 *  separated  by  at least  eps*fabs(x)+(tol/3), then  fmin  approximates
 *  the abcissa of the global minimum of  f  on the interval  ax,bx  with
 *  an error less than  3*eps*fabs(fmin)+tol.  if   f   is  not  unimodal,
 *  then fmin may approximate a local, but perhaps non-global, minimum to
 *  the same accuracy.
 *      this function subprogram is a slightly modified  version  of  the
 *  algol  60 procedure  localmin  given in richard brent, algorithms for
 *  minimization without derivatives, prentice-hall, inc. (1973).
 */

  tol1 = (eps*eps) + 1.0;
  a = ax;
  b = bx;
  v = a + INV_GOLD * (b - a);
  w = v;
  x = v;
  e = 0.;
  fx = func(x);
  fv = fx;
  fw = fx;
  tol3 = tol / 3.;
  
  /*  main loop starts here */
  
  while (1) {
    xm = (a + b) * .5;
    tol1 = eps * fabs(x) + tol3;
    t2 = tol1 * 2.;
    
    /*  check stopping criterion */
    if (fabs(x - xm) <= t2 - (b - a) * .5) break;
    niter++;

    p = 0.;
    q = 0.;
    r__ = 0.;

    if (fabs(e) > tol1) {
      /*  fit parabola */
      r__ = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r__;
      q = (q - r__) * 2.;
      if (q <= 0.) q = -q;
      else p = -p;
      r__ = e;
      e = d__;
    }

    if (fabs(p) >= fabs(q * .5 * r__) || p <= q * (a - x) || p >= q * (b - x)) {
      /*  a golden-section step */
      if (x >= xm) e = a - x;
      else e = b - x;
      d__ = INV_GOLD * e;
      if (debug) fprintf(stderr,"Golden section! %g %g\n",x,d__);
    }
    else {    
      /*  a parabolic-interpolation step */
      d__ = p / q;
      u = x + d__;
      
      /*  f must not be evaluated too close to ax or bx */
      if (!(u - a >= t2 && b - u >= t2)) {
	d__ = tol1;
	if (x >= xm) d__ = -d__;
      }
      if (debug) fprintf(stderr,"Parabolic step! %g %g\n",x,d__);
    }
    
    /*  f must not be evaluated too close to x */
    if (fabs(d__) < tol1) {
      if (d__ <= 0.) u = x - tol1;
      else u = x + tol1;
    }
    else u = x + d__;

    fu = func(u);
    
    /*  update  a, b, v, w, and x */
    
    if (fx <= fu) {
      if (u >= x) b = u;
      else a = u;
    }
    if (fu > fx) {
      if (fu > fw && w != x) {
	if (!(fu > fv && v != x && v != w)) {
	  v = u;
	  fv = fu;
	}
      }
      else {
	v = w;
	fv = fw;
	w = u;
	fw = fu;
      }
    }
    else {
      if (u >= x) a = x;
      else b = x;
      v = w;
      fv = fw;
      w = x;
      fw = fx;
      x = u;
      fx = fu;
    }
  }  /*  end of main loop */
  
  if (debug) fprintf(stderr,"fmin1D returning %f after %d iterations\n",x,niter);
  return x;
}

