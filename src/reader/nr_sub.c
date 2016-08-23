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
/* These routines substitute for missing Numerical Recipes utilities */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fmri.h"
#include "misc.h"
#include "nr_sub.h"

float **nr_matrix(int nrl, int nrh, int ncl, int nch)
{
  int i;
  float** m;

  if (!(m=(float**)malloc((unsigned)(nrh-nrl+1)*sizeof(float*))))
    Abort("Unable to allocate %d bytes!\n",
	  (unsigned)(nrh-nrl+1)*sizeof(float*));
  m -= nrl;

  for (i=nrl; i<=nrh; i++) {
    if (!(m[i]= (float*)malloc((unsigned)(nch-ncl+1)*sizeof(float))))
      Abort("Unable to allocate %d bytes!\n",
	    (unsigned)(nch-ncl+1)*sizeof(float));
    m[i] -= ncl;
  }

  return m;
}

void nr_free_matrix(float**m, int nrl, int nrh, int ncl, int nch)
{
  int i;

  for (i=nrh; i>=nrl;i--) free((void*)(m[i]+ncl));
  free((void*)(m+nrl));
}

fcomplex nr_Csub(fcomplex a, fcomplex b)
{
  fcomplex c;
  c.r= a.r - b.r;
  c.i= a.i - b.i;
  return c;
}

fcomplex nr_Cmul(fcomplex a, fcomplex b)
{
  fcomplex c;
  c.r= (a.r*b.r) - (a.i*b.i);
  c.i= (a.r*b.i) + (a.i*b.r);
  return c;
}

fcomplex nr_Complex(float re, float im)
{
  fcomplex c;
  c.r= re;
  c.i= im;
  return c;
}

fcomplex nr_Conjg(fcomplex z)
{
  fcomplex c;
  c.r= z.r;
  c.i= -z.i;
  return c;
}

float nr_Cabs(fcomplex z)
{
  float result;

  result= (float)sqrt((z.r*z.r) + (z.i*z.i));
  return result;
}

fcomplex nr_Csqrt(fcomplex z)
{
  fcomplex c;
  float x,y,w,r;
  if ((z.r==0.0) && (z.i==0.0)) {
    c.r= c.i= 0.0;
    return c;
  }
  else {
    x= fabs(z.r);
    y= fabs(z.i);
    if (x>=y) {
      r=y/x;
      w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
    }
    else {
      r= x/y;
      w= sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
    }
    if (z.r>=0) {
      c.r= w;
      c.i= z.i/(2.0*w);
    }
    else {
      c.i= (z.i >= 0.0) ? w : -w;
      c.r= z.i/(2.0*c.i);
    }
    return c;
  }
}

