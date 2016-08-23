/****************************************************************************
 * spline.c
 * Author Chris Rodriguez, Joel Welling
 * Copyright 1999, Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Permission use, copy, and modify this software and its documentation
 * without fee for personal use or use within your organization is hereby
 * granted, provided that the above copyright notice is preserved in all
 * copies and that that copyright and this permission notice appear in
 * supporting documentation.  Permission to redistribute this software to
 * other organizations or individuals is not granted;  that must be
 * negotiated with the PSC.  Neither the PSC nor Carnegie Mellon
 * University make any representations about the suitability of this
 * software for any purpose.  It is provided "as is" without express or
 * implied warranty.
 *****************************************************************************/
/*
 * This module provides a simple spline class.  It is derived from the starter
 * code for CMU 15-462 (Intro to Computer Graphics) Assignment 2, which
 * contained the following information:
 *
 *
 * spline.c
 * Starter code for 15-462, Computer Graphics 1,
 * Assignment 2: Generalized Cylinders
 *
 * Chris Rodriguez
 *
 * 4 Feb 1999
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mri.h"
#include "fmri.h"
#include "spline.h"

/* We need a constant to tell us when to actually worry about a spline
 * calc with a value not between 0.0 and 1.0 .  Some small outside values
 * can arise by rounding errors.
 */
#define LOC_TOL 0.0001

static void* safe_malloc( int nbytes ) 
{
  void* result;
  if (!(result= malloc(nbytes))) 
    Abort("spline:safe_malloc: Unable to allocate %d bytes!\n",nbytes);
  return result;
}

Spline* spl_create(long stride_in, long extent_in, SplineType type_in, 
		   double* pts_in)
{
  Spline *temp = NULL;

  if (type_in==SPL_BEZIER && (extent_in % 3) != 1) {
    Abort("spl_create: Bezier splines must have 3n+1 points; found %d!\n",
	  extent_in);
  }
  
  if (extent_in < 4) {
    Abort("spl_create: need at least 4 samples to spline; found %d!\n",
	  extent_in);
  }  

  temp= (Spline *)safe_malloc(sizeof(Spline));
  temp->stride= stride_in;
  temp->extent= extent_in;
  temp->type= type_in;
  temp->aux= 0.5;
  temp->dataPtr= pts_in;
  
  return (temp);
}

void spl_reset( Spline* cp, long stride_in, long extent_in, double* pts_in )
{
  if (cp->stride != stride_in)
    Abort("spl_reset: tried to change stride from %ld to %ld!\n",
	  cp->stride, stride_in);
  if (cp->extent != extent_in)
    Abort("spl_reset: tried to change extent from %ld to %ld!\n",
	  cp->extent, extent_in);
  cp->dataPtr= pts_in;
}

void spl_destroy(Spline *cp) 
{
  free(cp);
}

void spl_set_tension(Spline *sp, double val) 
{
  sp->aux= val;
}

double spl_get_tension(const Spline *sp) 
{
  return sp->aux;
}

void spl_dump(FILE* ofile, const Spline *cp) 
{
  int i;
  switch(cp->type) {
    case SPL_BEZIER:
      fprintf(ofile,"Bezier Spline: ");
      break;
    case SPL_CATMULLROM:
      fprintf(ofile,"Catmull Rom Spline; tension %f: ",cp->aux);
      break;
    case SPL_BSPLINE:
      fprintf(ofile,"B-spline: ");
      break;
    default:
      fprintf(ofile,"Unknown Spline; tension %f: ",cp->aux);
      break;
  }
  fprintf(ofile,
	  "stride %d; extent %d, data loc 0x%lx\n",
	  cp->stride,cp->extent,
	  (long)cp->dataPtr);
}

Spline* spl_copy(Spline *in) 
{
  Spline* out= (Spline*)safe_malloc(sizeof(Spline));
  out->type= in->type;
  out->aux= in->aux;
  out->stride= in->stride;
  out->extent= in->extent;
  out->dataPtr= in->dataPtr;
  return out;
}

void spl_calc( double* out, Spline* ctl, long n, long offset,
	       double loc)
{
  int i;

  if (loc<0.0 || loc>1.0) {
    if (loc<-1.0*LOC_TOL || loc>1.0+LOC_TOL)
      Abort("spl_calc: out-of-range loc %f\n",loc);
    if (loc<0.0) loc= 0.0;
    if (loc>1.0) loc= 1.0;
  }

  switch (ctl->type) {
  case SPL_BEZIER:
    {
      double u;
      long nseg;
      long seg;
      double A,B,C,D;
      double* c1;
      double* c2;
      double* c3;
      double* c4;

      nseg= (ctl->extent - 1)/3;
      u= nseg*loc;
      seg= (int)u;
      if (seg==nseg) seg--;
      u -= seg;

      A= (1.0-u)*(1.0-u)*(1.0-u);
      B= 3*u*(1.0-u)*(1.0-u);
      C= 3*u*u*(1.0-u);
      D= u*u*u;

      c1= ctl->dataPtr + 3*ctl->stride*seg + offset;
      c2= c1+ctl->stride;
      c3= c2+ctl->stride;
      c4= c3+ctl->stride;
      for (i=0; i<n; i++)
	*out++= A*(*c1++) + B*(*c2++) + C*(*c3++) + D*(*c4++);
    }
    break;
  case SPL_CATMULLROM:
    {
      double u;
      double aux;
      long seg;
      double A,B,C,D;
      double* c1;
      double* c2;
      double* c3;
      double* c4;

      aux= ctl->aux;
      u= ((ctl->extent - 1)*loc);
      seg= (int)u;
      if (seg==ctl->extent-1) seg--;
      u -= seg;
#ifdef never
      A= ((-aux*u*u*u)+(2.0*aux*u*u)+(-aux*u)+0.0);
      B= (((2-aux)*u*u*u)+((aux-3)*u*u)+1.0);
      C= (((aux-2)*u*u*u)+((3.0-2.0*aux)*u*u)+(aux*u));
      D= ((aux*u*u*u)-(aux*u*u));
#endif
      A= ((-aux*u+2.0*aux)*u - aux)*u+0.0;
      B= (((2-aux)*u+(aux-3))*u*u)+1.0;
      C= (((aux-2)*u+(3.0-2.0*aux))*u + aux)*u;
      D= aux*(u-1)*u*u;
      c1= ctl->dataPtr + ctl->stride*(seg-1) + offset;
      c2= c1+ctl->stride;
      c3= c2+ctl->stride;
      c4= c3+ctl->stride;
      if (seg==0) c1=c2; /* at left boundary, replicate 1st point */
      if (seg==ctl->extent-2) c4=c3; /* at right boundary, replicate last pt */
      for (i=0; i<n; i++)
	*out++= A*(*c1++) + B*(*c2++) + C*(*c3++) + D*(*c4++);
    }
    break;
  case SPL_BSPLINE:
    {
      double u;
      double aux;
      long seg;
      double A,B,C,D;
      double* c1;
      double* c2;
      double* c3;
      double* c4;

      aux= ctl->aux;
      u= ((ctl->extent - 1)*loc);
      seg= (int)u;
      if (seg==ctl->extent-1) seg--;
      u -= seg;

      A= ((1.0-u)*(1.0-u)*(1.0-u))/6.0;
      B= (3.0*u*u*u - 6.0*u*u + 4.0)/6.0;
      C= (-3.0*u*u*u + 3.0*u*u + 3.0*u +1)/6.0;
      D= (u*u*u)/6.0;
#ifdef never
      A= ((1.0-u)*(1.0-u)*(1.0-u))/6.0;
      B= ((3.0*u-6.0)*u*u + 4.0)/6.0;
      C= (((-3.0*u+3.0)*u+3.0)*u +1)/6.0;
      D= (u*u*u)/6.0;
#endif
      c1= ctl->dataPtr + ctl->stride*(seg-1) + offset;
      c2= c1+ctl->stride;
      c3= c2+ctl->stride;
      c4= c3+ctl->stride;
      if (seg==0) c1=c2; /* at left boundary, replicate 1st point */
      if (seg==ctl->extent-2) c4=c3; /* at right boundary, replicate last pt */
      for (i=0; i<n; i++)
	*out++= A*(*c1++) + B*(*c2++) + C*(*c3++) + D*(*c4++);
    }
    break;
  default: 
    {
      Abort("spl_calc: Unimplemented spline type %d\n",ctl->type);
    }
  }


}

void spl_grad( double* out, Spline* ctl, long n, long offset,
	       double loc) 
{
  int i;

  if (loc<0.0 || loc>1.0) {
    if (loc<-1.0*LOC_TOL || loc>1.0+LOC_TOL)
      Abort("spl_grad: out-of-range loc %f\n",loc);
    if (loc<0.0) loc= 0.0;
    if (loc>1.0) loc= 1.0;
  }

  switch (ctl->type) {
  case SPL_BEZIER:
    {
      double u;
      long nseg;
      long seg;
      double A,B,C,D;
      double* c1;
      double* c2;
      double* c3;
      double* c4;

      nseg= (ctl->extent - 1)/3;
      u= nseg*loc;
      seg= (int)u;
      if (seg==nseg) seg--;
      u -= seg;

      A= -3.0*(1.0-u)*(1.0-u);
      B= (3.0-9.0*u)*(1.0-u);
      C= (6.0-9.0*u)*u;
      D= 3.0*u*u;

      c1= ctl->dataPtr + 3*ctl->stride*seg + offset;
      c2= c1+ctl->stride;
      c3= c2+ctl->stride;
      c4= c3+ctl->stride;
      for (i=0; i<n; i++)
	*out++= A*(*c1++) + B*(*c2++) + C*(*c3++) + D*(*c4++);
    }
    break;
  case SPL_CATMULLROM:
    {
      double u;
      double aux;
      long seg;
      double A,B,C,D;
      double* c1;
      double* c2;
      double* c3;
      double* c4;

      aux= ctl->aux;
      u= ((ctl->extent - 1)*loc);
      seg= (int)u;
      if (seg==ctl->extent-1) seg--;
      u -= seg;

#ifdef never
      A= ((-3.0*aux*u*u)+(4.0*aux*u)-aux);
      B= ((3.0*(2-aux)*u*u)+(2.0*(aux-3)*u));
      C= ((3.0*(aux-2)*u*u)+(2.0*(3.0-2.0*aux)*u)+aux);
      D= ((3.0*aux*u*u)-(2.0aux*u));
#endif
      A= ((-3.0*aux*u)+(4.0*aux))*u - aux;
      B= ((3.0*(2-aux)*u)+(2.0*(aux-3)))*u;
      C= ((3.0*(aux-2)*u)+(2.0*(3.0-2.0*aux)))*u + aux;
      D= ((3.0*aux*u)-(2.0*aux))*u;

      c1= ctl->dataPtr + ctl->stride*(seg-1) + offset;
      c2= c1+ctl->stride;
      c3= c2+ctl->stride;
      c4= c3+ctl->stride;
      if (seg==0) c1=c2; /* at left boundary, replicate 1st point */
      if (seg==ctl->extent-2) c4=c3; /* at right boundary, replicate last pt */
      for (i=0; i<n; i++)
	*out++= A*(*c1++) + B*(*c2++) + C*(*c3++) + D*(*c4++);
    }
    break;
  case SPL_BSPLINE:
    {
      double u;
      double aux;
      long seg;
      double A,B,C,D;
      double* c1;
      double* c2;
      double* c3;
      double* c4;

      aux= ctl->aux;
      u= ((ctl->extent - 1)*loc);
      seg= (int)u;
      if (seg==ctl->extent-1) seg--;
      u -= seg;

#ifdef never
      A= -3.0*((1.0-u)*(1.0-u))/6.0;
      B= (9.0*u*u - 12.0*u)/6.0;
      C= (-9.0*u*u + 6.0*u + 3.0)/6.0;
      D= 3.0*(u*u)/6.0;
#endif
      A= -((1.0-u)*(1.0-u))/2.0;
      B= ((3.0*u - 4.0)*u)/2.0;
      C= ((-3.0*u + 2.0)*u + 1.0)/2.0;
      D= (u*u)/2.0;
      c1= ctl->dataPtr + ctl->stride*(seg-1) + offset;
      c2= c1+ctl->stride;
      c3= c2+ctl->stride;
      c4= c3+ctl->stride;
      if (seg==0) c1=c2; /* at left boundary, replicate 1st point */
      if (seg==ctl->extent-2) c4=c3; /* at right boundary, replicate last pt */
      for (i=0; i<n; i++)
	*out++= A*(*c1++) + B*(*c2++) + C*(*c3++) + D*(*c4++);
    }
    break;
  default: 
    {
      Abort("spl_grad: Unimplemented spline type %d\n",ctl->type);
    }
  }


}
void spl_gradsqr( double* out, Spline* ctl, long n, long offset,
		  double loc)
{
  int i;

  if (loc<0.0 || loc>1.0) {
    if (loc<-1.0*LOC_TOL || loc>1.0+LOC_TOL)
      Abort("spl_gradsqr: out-of-range loc %f\n",loc);
    if (loc<0.0) loc= 0.0;
    if (loc>1.0) loc= 1.0;
  }

  switch (ctl->type) {
  case SPL_BEZIER:
    {
      double u;
      long nseg;
      long seg;
      double A,B,C,D;
      double* c1;
      double* c2;
      double* c3;
      double* c4;

      nseg= (ctl->extent - 1)/3;
      u= nseg*loc;
      seg= (int)u;
      if (seg==nseg) seg--;
      u -= seg;

      A= 6.0*(1.0-u);
      B= (-12.0+18.0*u);
      C= (6.0-18.0*u);
      D= 6.0*u;

      c1= ctl->dataPtr + 3*ctl->stride*seg + offset;
      c2= c1+ctl->stride;
      c3= c2+ctl->stride;
      c4= c3+ctl->stride;
      for (i=0; i<n; i++)
	*out++= A*(*c1++) + B*(*c2++) + C*(*c3++) + D*(*c4++);
    }
    break;
  case SPL_CATMULLROM:
    {
      double u;
      double aux;
      long seg;
      double A,B,C,D;
      double* c1;
      double* c2;
      double* c3;
      double* c4;

      aux= ctl->aux;
      u= ((ctl->extent - 1)*loc);
      seg= (int)u;
      if (seg==ctl->extent-1) seg--;
      u -= seg;

      A= (-6.0*aux*u)+(4.0*aux);
      B= (6.0*(2-aux)*u)+(2.0*(aux-3.0));
      C= (6.0*(aux-2)*u)+(2.0*(3.0-2.0*aux));
      D= (6.0*aux*u)-(2.0*aux);

      c1= ctl->dataPtr + ctl->stride*(seg-1) + offset;
      c2= c1+ctl->stride;
      c3= c2+ctl->stride;
      c4= c3+ctl->stride;
      if (seg==0) c1=c2; /* at left boundary, replicate 1st point */
      if (seg==ctl->extent-2) c4=c3; /* at right boundary, replicate last pt */
      for (i=0; i<n; i++)
	*out++= A*(*c1++) + B*(*c2++) + C*(*c3++) + D*(*c4++);
    }
    break;
  case SPL_BSPLINE:
    {
      double u;
      double aux;
      long seg;
      double A,B,C,D;
      double* c1;
      double* c2;
      double* c3;
      double* c4;

      aux= ctl->aux;
      u= ((ctl->extent - 1)*loc);
      seg= (int)u;
      if (seg==ctl->extent-1) seg--;
      u -= seg;

      A= 1.0-u;
      B= 3.0*u - 2.0;
      C= -3.0*u + 1.0;
      D= u;
      c1= ctl->dataPtr + ctl->stride*(seg-1) + offset;
      c2= c1+ctl->stride;
      c3= c2+ctl->stride;
      c4= c3+ctl->stride;
      if (seg==0) c1=c2; /* at left boundary, replicate 1st point */
      if (seg==ctl->extent-2) c4=c3; /* at right boundary, replicate last pt */
      for (i=0; i<n; i++)
	*out++= A*(*c1++) + B*(*c2++) + C*(*c3++) + D*(*c4++);
    }
    break;
  default: 
    {
      Abort("spl_gradsqr: Unimplemented spline type %d\n",ctl->type);
    }
  }

}

