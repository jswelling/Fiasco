/************************************************************
 *                                                          *
 *  vpolygon.c                                               *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1997 Department of Statistics,         *
 *                        Carnegie Mellon University        *
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
 *  Original programming by Joel Welling, 2/03              *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "par.h"
#include "bio.h"
#include "array.h"
#include "acct.h"
#include "mri.h"
#include "fmri.h"
#include "misc.h"
#include "stdcrg.h"
#include "vpolygon.h"

static char rcsid[] = "$Id: vpolygon.c,v 1.7 2005/07/07 20:02:27 welling Exp $";

static FILE* voronoiPostscript= NULL;

static int debug= 0;

typedef double Vec2[2];
#define V2_COPY( vout, vin ) { vout[0]=vin[0]; vout[1]=vin[1]; }
#define V2_MINUS( v1, v2 ) { v1[0] -= v2[0]; v1[1] -= v2[1]; }
#define V2_PLUS( v1, v2 ) { v1[0] += v2[0]; v1[1] += v2[1]; }
#define V2_DOT( v1, v2 ) ( v1[0]*v2[0] + v1[1]*v2[1] )
#define V2_CROSSMAG( v1, v2 ) ( v1[0]*v2[1] - v2[0]*v1[1] )
#define V2_ISLEFT( pvec, edge ) ( V2_CROSSMAG( pvec, edge ) >= 0.0 )
#define V2_MULT( v, factor ) { v[0] *= factor; v[1] *= factor; }

VPoly* vply_create(double ctrX, double ctrY, int initialN)
{
  VPoly* result= NULL;
  VWing* wings= NULL;

  if (!(result=(VPoly*)malloc(sizeof(VPoly))))
    Abort("%s: unable to allocate %d bytes!\n",sizeof(VPoly));
  if (!(result->wings=(VWing*)malloc(initialN*sizeof(VWing))))
    Abort("%s: unable to allocate %d bytes!\n",sizeof(VWing));
  result->wingBufSize= initialN;
  result->nWings= 0;
  result->ctr[0]= ctrX;
  result->ctr[1]= ctrY;
  result->maxDist= result->maxDistSqr= 0.0;
  return result;
}

void vply_destroy(VPoly* vp)
{
  free(vp->wings);
  free(vp);
}

void vply_setDebug(const int i)
{
  debug= i;
}

void vply_dump(VPoly* vp, FILE* ofile)
{
  int i;
  fprintf(ofile,"VPoly at ctr (%f %f) has %d wings of %d allocated:\n",
	  vp->ctr[0],vp->ctr[1], vp->nWings, vp->wingBufSize);
  for (i=0; i<vp->nWings; i++) {
    VWing w= vp->wings[i];
    fprintf(ofile,"   (%f %f); theta= %f\n",w.vtx[0],w.vtx[1],w.theta);
  }
}

void vply_add(VPoly* vp, double x, double y)
{
  VWing* thisWing;
  double distSqr;

  if (vp->nWings >= vp->wingBufSize) {
    int newBufSize= 2*vp->wingBufSize;
    if (debug) fprintf(stderr,"Reallocating vp %lx from %d to %d\n",
		       (long)vp,vp->wingBufSize,newBufSize);
    if (!(vp->wings=(VWing*)realloc(vp->wings, newBufSize*sizeof(VWing))))
      Abort("%s: unable to reallocate up to %d bytes!\n",
	    newBufSize*sizeof(VWing));
    vp->wingBufSize= newBufSize;
  }
  thisWing= vp->wings + vp->nWings++;
  thisWing->vtx[0]= x;
  thisWing->vtx[1]= y;
  thisWing->theta= 0.0;
  distSqr= (x - vp->ctr[0])*(x - vp->ctr[0])
    + (y - vp->ctr[1])*(y - vp->ctr[1]);
  if (distSqr>vp->maxDistSqr) {
    vp->maxDistSqr= distSqr;
    vp->maxDist= sqrt(distSqr);
  }

}

static int wing_compare( const void* p1, const void* p2 )
{
  VWing* wing1= (VWing*)p1;
  VWing* wing2= (VWing*)p2;
  if (wing1->theta<wing2->theta) return -1;
  else if (wing1->theta>wing2->theta) return 1;
  else return 0;
}

void vply_sort(VPoly* vp)
{
  int i;
  VWing* wings= vp->wings;
  
  for (i=0; i<vp->nWings; i++) {
    wings[i].theta= atan2( wings[i].vtx[1] - vp->ctr[1],
			   wings[i].vtx[0] - vp->ctr[0] );
  }
  qsort(wings, vp->nWings, sizeof(VWing), wing_compare);

#ifdef never
  fprintf(stderr,"Sorted:\n");
  for (i=0; i<vp->nWings; i++) {
    fprintf(stderr,"vtx %d: %f %f; theta= %f\n",i,
	    wings[i].vtx[0],wings[i].vtx[1],wings[i].theta);
  }
#endif
}

double vply_calcArea( VPoly* vp )
{
  VWing* wings= vp->wings;
  VWing* w1;
  VWing* w2;
  Vec2 ctr;
  double result= 0.0;
  double dResult= 0.0;
  double dx1;
  double dx2;
  double dy1;
  double dy2;
  int i;

  /* We need to check to make sure this isn't a clipped polygon
   * consisting of only a single edge.
   */
  if (vp->nWings<3) return 0.0;

  /* Find an internal point of the convex polygon */
  V2_COPY(ctr, vp->wings[0].vtx);
  for (i=1; i<vp->nWings; i++)
    V2_PLUS(ctr, vp->wings[i].vtx);
  V2_MULT(ctr, 1.0/(double)(vp->nWings));

  w2= wings+(vp->nWings-1);;
  dx2= w2->vtx[0] - ctr[0];
  dy2= w2->vtx[1] - ctr[1];
  for (i=0; i<vp->nWings; i++) {
    w1= wings+i;
    dx1= w1->vtx[0] - ctr[0];
    dy1= w1->vtx[1] - ctr[1];
    dResult= 0.5*(dx2*dy1 - dx1*dy2);  /* half the cross product */
#ifdef never
    fprintf(stderr,"Delta: %f\n",dResult);
#endif
    if (dResult<0.0) {
      vply_dump(vp,stderr);
      Abort("vply_calcArea: internal error: got an incorrectly oriented wing pair!\n");
    }
    result += dResult;
    w2= w1;
    dx2= dx1;
    dy2= dy1;
  }

  return result;
}

void vply_writePS(VPoly* vp)
{
  VWing* wings= vp->wings;
  FILE* f= voronoiPostscript;
  int i;

#define SCALE 9.0

  if (vp->nWings>1) {
    for (i=0; i<vp->nWings; i++) {
      if (i==0) {
	fprintf(f,"newpath %f %f moveto\n",
		SCALE*wings[0].vtx[0]+300,SCALE*wings[0].vtx[1]+350);
      }
      else {
	fprintf(f,"%f %f lineto\n",
		SCALE*wings[i].vtx[0]+300,SCALE*wings[i].vtx[1]+350);
      }
    }
    fprintf(f,"closepath stroke \n\n");
  }
  fprintf(f,"%f %f marksample\n",
	  SCALE*vp->ctr[0]+300,SCALE*vp->ctr[1]+350);
}

void vply_prepPS(char* fname)
{
  if (voronoiPostscript) {
    /* Close off the old one */
    vply_finishPS();
  }

  voronoiPostscript= fopen("polys.ps","w");
  fprintf(voronoiPostscript,
	  "%%!PS-Adobe-1.0\n%%%%BoundingBox: 0 0 612 792\n");
  fprintf(voronoiPostscript,
	  "/marksample { moveto 0 1 rlineto 1 0 rlineto -1 -1 rlineto closepath stroke} def\n");
}

void vply_finishPS()
{	  
  if (!voronoiPostscript)
    Abort("vpolygon: vply_finishPS called without vply_prepPS!\n");
  fprintf(voronoiPostscript,"showpage\n");
  fclose(voronoiPostscript);
  voronoiPostscript= NULL;
}

static void calcIntercept( const Vec2 p1, const Vec2 left, 
			   double* oldPt, const Vec2 newPt )
{
  Vec2 delta;
  Vec2 tmp;
  double lamda;

  V2_COPY( delta, newPt );
  V2_MINUS( delta, oldPt );

  V2_COPY( tmp, p1 );
  V2_MINUS( tmp, oldPt );

  lamda= V2_DOT( tmp, left ) / V2_DOT( delta, left );
  V2_MULT( delta, lamda );
  V2_PLUS( oldPt, delta );
}

static VPoly* vply_clipAgainstHalfplane( VPoly* vp, Vec2 p1, Vec2 p2 )
{
  Vec2 edge; /* vector from p1 to p2 */
  Vec2 left; /* vector perp to edge on left side */
  Vec2 tmp;
  Vec2 lastInside;
  Vec2 lastOutside;
  int i;
  VPoly* newPoly= vply_create(vp->ctr[0], vp->ctr[1], vp->nWings);
  int changed= 0;
  typedef enum {IN, OUT, STARTING} State;
  State state, firstState;
  
  V2_COPY( edge, p2 );
  V2_MINUS( edge, p1 );
  left[0]= -edge[1];
  left[1]= edge[0];

  state= STARTING;
  firstState= STARTING; /* convince complier it's initialized */
  for (i=0; i<vp->nWings; i++) {
    V2_COPY( tmp, vp->wings[i].vtx );
    V2_MINUS( tmp, p1 );
    if (V2_ISLEFT(tmp, edge)) {
      /* Point is inside */
      if (debug) fprintf(stderr,"i= %d; pt= (%f %f), state %d, inside\n",
			 i, tmp[0]+p1[0], tmp[1]+p1[1], (int)state);
      switch (state) {
      case IN:
	/* Still in */
	V2_COPY(lastInside, vp->wings[i].vtx);
	vply_add(newPoly, lastInside[0], lastInside[1]);
	break;
      case OUT:
	/* Transition to in; add transition point and new point */
	calcIntercept( p1, left, lastOutside, vp->wings[i].vtx );
	vply_add(newPoly, lastOutside[0], lastOutside[1]);
	V2_COPY(lastInside, vp->wings[i].vtx);
	vply_add(newPoly, lastInside[0], lastInside[1]);
	break;
      case STARTING:
	V2_COPY(lastInside, vp->wings[i].vtx);
	vply_add(newPoly, lastInside[0], lastInside[1]);
	firstState= IN;
	break;
      }
      state= IN;
    }
    else {
      /* point is outside */
      changed= 1;
      if (debug) fprintf(stderr,"i= %d; pt= (%f %f), state %d, outside\n",
			 i, tmp[0]+p1[0], tmp[1]+p1[1], (int)state);
      switch (state) {
      case IN: 
	/* Transition to out; add transition point */
	calcIntercept( p1, left, lastInside, vp->wings[i].vtx );
	vply_add(newPoly, lastInside[0], lastInside[1]);
	break;
      case OUT:
	/* Still out */
	break;
      case STARTING:
	/* Starting up outside; do nothing */
	firstState= OUT;
	break;
      }
      state= OUT;
      V2_COPY(lastOutside, vp->wings[i].vtx);
    }
  }

  if (state != firstState) {
    if (state==IN) {
      calcIntercept( p1, left, lastInside, vp->wings[0].vtx );
      vply_add(newPoly, lastInside[0], lastInside[1]);
    }
    else {
      calcIntercept( p1, left, lastOutside, vp->wings[0].vtx );
      vply_add(newPoly, lastOutside[0], lastOutside[1]);
    }
  }

  if (debug) {
    fprintf(stderr,"Clipping complete; original:\n");
    vply_dump(vp,stderr);
    fprintf(stderr,"          clipped version:\n");
    vply_dump(newPoly,stderr);
  }

  vply_destroy(vp);
  return newPoly;
}

static void calcCircleIntercept( const Vec2 ctr, const double radSqr, 
				 double* oldPt, const Vec2 newPt )
{
  Vec2 delta;
  Vec2 tmp;
  double lamda;
  double a;
  double b;
  double c;
  double disc;
  double discRoot;

  V2_COPY( delta, newPt );
  V2_MINUS( delta, oldPt );

  V2_COPY( tmp, oldPt );
  V2_MINUS( tmp, ctr );

  a= V2_DOT( delta, delta );
  b= 2.0*V2_DOT( tmp, delta );
  c= V2_DOT( tmp, tmp ) - radSqr;
  disc= b*b - 4*a*c;
  if (disc<0.0) 
    Abort("vply_calcCircleIntercept: internal error; < 0.0! a=%f, b= %f, c= %f\n",
	  a, b, c);
  discRoot= sqrt(disc);
  if ( -b - discRoot >= 0.0 ) lamda= (-b - discRoot)/(2.0*a);
  else lamda= (-b + discRoot)/(2.0*a);

  V2_MULT( delta, lamda );
  V2_PLUS( oldPt, delta );
}

VPoly* vply_clip(VPoly* vp, const VPoly* clipAgainst)
{
  VPoly* result= vp;
  int i;

  if (clipAgainst->nWings<1) {
    /* Return an empty poly, since clip region has 0 area */
    result= vply_create(vp->ctr[0],vp->ctr[1],10);
    vply_destroy(vp);
  }
  else {
    for (i=1; i<clipAgainst->nWings; i++)
      result= vply_clipAgainstHalfplane(result, 
					clipAgainst->wings[i].vtx,
					clipAgainst->wings[i-1].vtx);
  }

  return result;
}

VPoly* vply_clipToCircle(VPoly* vp, const double xCtr, const double yCtr, 
			 const double radius)
{
  Vec2 tmp;
  Vec2 ctr;
  Vec2 lastInside;
  Vec2 lastOutside;
  int i;
  VPoly* newPoly= vply_create(vp->ctr[0], vp->ctr[1], vp->nWings);
  int changed= 0;
  typedef enum {IN, OUT, STARTING} State;
  State state, firstState;
  double distSqr;
  
  ctr[0]= xCtr;
  ctr[1]= yCtr;
  state= STARTING;
  firstState= STARTING; /* convince complier it's initialized */
  for (i=0; i<vp->nWings; i++) {
    V2_COPY( tmp, vp->wings[i].vtx );
    V2_MINUS( tmp, ctr );
    distSqr= V2_DOT( tmp, tmp );
    if (distSqr <= radius*radius) {
      /* Point is inside */
      if (debug) fprintf(stderr,"i= %d; pt= (%f %f), state %d, inside\n",
			 i, tmp[0]+ctr[0], tmp[1]+ctr[1], (int)state);
      switch (state) {
      case IN:
	/* Still in */
	V2_COPY(lastInside, vp->wings[i].vtx);
	vply_add(newPoly, lastInside[0], lastInside[1]);
	break;
      case OUT:
	/* Transition to in; add transition point and new point */
	calcCircleIntercept( ctr, radius*radius, 
			     lastOutside, vp->wings[i].vtx );
	vply_add(newPoly, lastOutside[0], lastOutside[1]);
	V2_COPY(lastInside, vp->wings[i].vtx);
	vply_add(newPoly, lastInside[0], lastInside[1]);
	break;
      case STARTING:
	V2_COPY(lastInside, vp->wings[i].vtx);
	vply_add(newPoly, lastInside[0], lastInside[1]);
	firstState= IN;
	break;
      }
      state= IN;
    }
    else {
      /* point is outside */
      changed= 1;
      if (debug) fprintf(stderr,"i= %d; pt= (%f %f), state %d, outside\n",
			 i, tmp[0]+ctr[0], tmp[1]+ctr[1], (int)state);
      switch (state) {
      case IN: 
	/* Transition to out; add transition point */
	calcCircleIntercept( ctr, radius*radius, 
			     lastInside, vp->wings[i].vtx );
	vply_add(newPoly, lastInside[0], lastInside[1]);
	break;
      case OUT:
	/* Still out */
	break;
      case STARTING:
	/* Starting up outside; do nothing */
	firstState= OUT;
	break;
      }
      state= OUT;
      V2_COPY(lastOutside, vp->wings[i].vtx);
    }
  }

  if (state != firstState) {
    if (state==IN) {
      calcCircleIntercept( ctr, radius*radius, 
			 lastInside, vp->wings[0].vtx );
      vply_add(newPoly, lastInside[0], lastInside[1]);
    }
    else {
      calcCircleIntercept( ctr, radius*radius,
			 lastOutside, vp->wings[0].vtx );
      vply_add(newPoly, lastOutside[0], lastOutside[1]);
    }
  }

  if (debug) {
    fprintf(stderr,"Circle Clipping complete; original:\n");
    vply_dump(vp,stderr);
    fprintf(stderr,"          clipped version:\n");
    vply_dump(newPoly,stderr);
  }

  vply_destroy(vp);
  return newPoly;
}

double vply_getMaxDistance( VPoly* vp )
{
  return vp->maxDist;
}
