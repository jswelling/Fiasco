/************************************************************
 *                                                          *
 *  fshrot3d.c                                              *
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
 *  Original programming by Joel Welling 2/99               *
 ************************************************************/

/*************************************************************
 * fourier_shift_rot performs a translation and rotation     *
 *   of a complex image using a decomposition of the         *
 *   rotation and translation into a series of axis-         *
 *   aligned shears.                                         *
 *                                                           *
 *   q is the input quaternion specifying the rotation       *
 *   dx, dy, and dz (input) specify the translation          *
 *   orig_image and moved_image are the input and            *
 *     output images respectively.  Both are assumed         *
 *     to have been allocated by the caller.                 *
 *   nx, ny, and nz are the image dimensions                 *
 *   length_x, length_y, length_z are voxel edge lengths     *
 *************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"

/* Notes-
 * -Note that the routine expects data presented in z-fastest order!
 */

static char rcsid[] = "$Id: fshrot3d.c,v 1.18 2004/03/19 18:27:18 welling Exp $";

/* How bad does cancellation have to get before we consider a quat bad? */
#define CANCELLATION_TOL 0.01

/* How long to look for a good quat */
#define REPAIR_QUAT_MAX_STEPS 10   

/* How long to iterate looking for Euler angle decompositions of quats */
#define EULER_DECOMP_MAX_ITER 100

/* How big of steps to take looking for a good quat */
#define REPAIR_QUAT_STEPSIZE 0.01

/* Accessor for 3D arrays */
#define MEM(matrix,nx,ny,nz,x,y,z) matrix[((((x)*ny)+(y))*nz)+(z)]

/* Angle conversion */
#define DEG2RAD (M_PI/180.0)
#define RAD2DEG (180.0/M_PI)

/* Decompositions of interest */
#define SHEAR_NONE 0 /* for housekeeping with failed factorizations */
#define SHEAR_YZXY 1
#define SHEAR_YXZY -1
#define SHEAR_ZXYZ 2
#define SHEAR_ZYXZ -2
#define SHEAR_XYZX 3
#define SHEAR_XZYX -3
#define SHEAR_YZXYXZY 4
#define SHEAR_ZXYZYXZ 5
#define SHEAR_XYZXZYX 6

typedef struct shear_param_struct {
  int decomposition;
  double a, b, c, d, e, f, g, h;
  double abar, bbar, cbar, dbar, ebar, fbar, gbar, hbar; /*long decomps only*/
  double quality; /* low values are good */
  int qualType;
} ShearParams;

typedef struct trans_param_struct {
  double dx, dy, dz;
} TransParams;

static int debug= 0;
static int shear_pattern= FR3D_SHEAR_4;
static int qual_measure_type= FR3D_QUAL_COX;
static int count_shear_x= 0;
static int count_shear_y= 0;
static int count_shear_z= 0;
static int count_set_phase= 0;
static int count_calls= 0;

void fshrot3d_set( int flag, int value )
{
  switch (flag) {
  case FR3D_DEBUG: {
    debug= value; 
  }
  break;
  case FR3D_SHEAR_PATTERN: {
    shear_pattern= value;
  }
  break;
  case FR3D_QUAL_MEASURE: {
    qual_measure_type= value;
  }
  break;
  }
}

int fshrot3d_get( int flag )
{
  switch (flag) {
  case FR3D_DEBUG: {
    return debug; 
  }
    /* break; -NOTREACHED- */
  case FR3D_SHEAR_PATTERN: {
    return shear_pattern;
  }
    /* break; -NOTREACHED- */
  case FR3D_QUAL_MEASURE: {
    return qual_measure_type;
  }
    /* break; -NOTREACHED- */
  default:
    Abort("fshrot3d_get: unknown flag %d!\n",flag);
  }
  return 0; /* not reached */
}

void fshrot3d_set_debug( int flag )
{
  fshrot3d_set( FR3D_DEBUG, flag );
}

void fshrot3d_clear_shear_counts()
{
  count_shear_x= count_shear_y= count_shear_z= count_calls= 
    count_set_phase= 0;
}

void fshrot3d_get_shear_counts( int* xcount, int* ycount, int* zcount,
				int* phase_count, 
				double *mean_shears_per_call)
{
  *xcount= count_shear_x;
  *ycount= count_shear_y;
  *zcount= count_shear_z;
  *phase_count= count_set_phase;
  *mean_shears_per_call= 
    (double)(count_shear_x + count_shear_y + count_shear_z + count_set_phase)
    / (double)count_calls;
}

static char* decompName(int dtype) {
  switch (dtype) {
  case SHEAR_YZXY: return "YZXY";
  case SHEAR_YXZY: return "YXZY";
  case SHEAR_ZXYZ: return "ZXYZ";
  case SHEAR_ZYXZ: return "ZYXZ";
  case SHEAR_XYZX: return "XYZX";
  case SHEAR_XZYX: return "XZYX";
  case SHEAR_YZXYXZY: return "YZXYXZY";
  case SHEAR_ZXYZYXZ: return "ZXYZYXZ";
  case SHEAR_XYZXZYX: return "ZXYZYXZ";
  case SHEAR_NONE: return "NONE";
  default: return "*UNKNOWN*";
  }
}

static int decompIsLong(int dtype) {
  return ((dtype == SHEAR_YZXYXZY)
	  || (dtype == SHEAR_ZXYZYXZ)
	  || (dtype == SHEAR_XYZXZYX));
}

static char* qualMeasureName(int qtype) {
  switch (qtype) {
  case FR3D_QUAL_UNSET: return "unset";
  case FR3D_QUAL_COX: return "Cox";
  case FR3D_QUAL_SUM_ABS: return "sum of absolute values";
  case FR3D_QUAL_SUM_SQR: return "sum of squares";
  case FR3D_QUAL_UNIT_CELL: return "unit cell";
  default: return "*UNKNOWN*";
  }
}

static void writeShearTrans(ShearParams *s, TransParams* t) 
{
  fprintf(stderr,
	  "decomp %7s : %f %f %f %f\n                 %f %f %f %f\n",
	  decompName(s->decomposition),
	  s->a,s->b,s->c,s->d,
	  s->e,s->f,s->g,s->h);
  if (decompIsLong(s->decomposition))
    fprintf(stderr,
	    "                 %f %f %f %f\n                 %f %f %f %f\n",
	  s->abar,s->bbar,s->cbar,s->dbar,
	  s->ebar,s->fbar,s->gbar,s->hbar);
  fprintf(stderr,"   quality %f (%s)\n",
	  s->quality,qualMeasureName(s->qualType));
  fprintf(stderr,"Adjusted trans params %f %f %f\n", t->dx, t->dy, t->dz);
}

/*
 * This routine finds Euler angles xout, yout, and zout such that the series of
 * axis-aligned rotations RzRyRxRyRz for the given angles (y and z being used
 * twice) is equal to the input quaternion rotation q.  Returns 1 on success,
 * zero on failure.
 */
static int euler5_decompose( Quat* q, double* xout, double* yout, double* zout ) {
  double theta_x=0.0, theta_y=0.0, theta_z=0.0;
  Quat Rx, Ry, Rz, temp, residual;
  double x,y,z,phi;
  int i;
  int iter= 0;

  if (debug) fprintf(stderr,"decomposing (%g %g %g %g)\n",q->x,q->y,q->z,q->w);

  do {
    (void)quat_from_axis_angle(&Rx,1.0,0.0,0.0,theta_x);
    (void)quat_from_axis_angle(&Ry,0.0,1.0,0.0,theta_y);
    (void)quat_from_axis_angle(&Rz,0.0,0.0,1.0,theta_z);

    /* Calculate RzRyRxRxRyRz */
    (void)quat_copy(&temp,q);
    (void)quat_conjugate(&temp);
    (void)quat_mult_right(&temp, &Rz);
    (void)quat_mult_right(&temp, &Ry);
    (void)quat_mult_right(&temp, &Rx);
    (void)quat_mult_right(&temp, &Rx);
    (void)quat_mult_right(&temp, &Ry);
    (void)quat_mult_right(&temp, &Rz);

    /* Make a new guess */
    quat_to_axis_angle(&temp,&x,&y,&z,&phi);
    /*
      fprintf(stderr,"new guess %g %g %g -> %g about (%g,%g,%g)\n",
      RAD2DEG*theta_x,RAD2DEG*theta_y,RAD2DEG*theta_z,
      RAD2DEG*phi,x,y,z);
      */
    theta_x -= 0.5*x*phi;
    theta_y -= 0.5*y*phi;
    theta_z -= 0.5*z*phi;
  } while (phi!=0.0 && iter++ < EULER_DECOMP_MAX_ITER);

  if (iter>= EULER_DECOMP_MAX_ITER) return 0; /* convergence failure */

  if (debug) {
    fprintf(stderr,"decomp results: %g %g %g after %d iterations\n", 
	    theta_x, theta_y, theta_z, iter);
    (void)quat_from_axis_angle(&Rx,1.0,0.0,0.0,theta_x);
    (void)quat_from_axis_angle(&Ry,0.0,1.0,0.0,theta_y);
    (void)quat_from_axis_angle(&Rz,0.0,0.0,1.0,theta_z);
    (void)quat_identity(&temp);
    (void)quat_mult_right(&temp, &Rz);
    (void)quat_mult_right(&temp, &Ry);
    (void)quat_mult_right(&temp, &Rx);
    (void)quat_mult_right(&temp, &Rx);
    (void)quat_mult_right(&temp, &Ry);
    (void)quat_mult_right(&temp, &Rz);
    fprintf(stderr,"goal: (%g %g %g %g) vs. result (%g %g %g %g)\n",
	    q->x,q->y,q->z,q->w,temp.x,temp.y,temp.z,temp.w);
  }

  *xout= 2*theta_x; /* to include both X rotations */
  *yout= theta_y;
  *zout= theta_z;
  return 1;
}

/* Returns 1 on success */
static int quat_to_shears_yzxy( Quat* q, ShearParams* s)
{
  double X= q->x;
  double Y= q->y;
  double Z= q->z;
  double W= q->w;
  double R= X*Y - Z*W;
  double S= Y*Z - X*W;

  if (debug) fprintf(stderr,"Decomposing quat (%f %f %f %f) as YZXY\n",
		     q->x,q->y,q->z,q->w);

  s->decomposition= SHEAR_YZXY;
  s->qualType= FR3D_QUAL_UNSET;

  if (X==0.0 && Y==0.0 && Z==0.0) {
    /* Special case the identity quaternion */
    s->a= s->b= s->c= s->d= s->e= s->f= s->g= s->h= 0.0;
    if (debug) {
      fprintf(stderr,
	      "Result of trivial decomp to YZXY: a=%f, b=%f, c=%f, d= %f\n",
	      s->a,s->b,s->c,s->d);
      fprintf(stderr,
	      "                                  e= %f, f=%f, g=%f, h=%f\n",
	      s->e,s->f,s->g,s->h); 
    }
    return 1;
  }
  else if (fabs(W)<CANCELLATION_TOL) {
    if (fabs(X)<CANCELLATION_TOL 
	|| fabs(Y)<CANCELLATION_TOL
	|| fabs(Z)<CANCELLATION_TOL) {
      /* Singular decomposition */
      if (debug)
	fprintf(stderr,"decomp failed- rot by 2pi in coordinate plane.\n");
      s->decomposition= SHEAR_NONE;
    }
    else {
      if (debug)
	fprintf(stderr,"small W expansion.\n");
      s->a= (Y*Y+X*X)*(1.0/(Y*Z) + (X*W)/(Y*Y*Z*Z) + (X*X*W*W)/(Y*Y*Y*Z*Z*Z));
      s->b= (-1.0-X*X)/(X*Y) 
	-W*(1.0-Y*Y+X*X*X*X+X*X*Y*Y)/(X*X*Y*Y*Z)
	-(1.0-X*X-2.0*Y*Y+Y*Y*Y*Y+X*X*X*X*Y*Y
	  +X*X*X*X*X*X+X*X*X*X+X*X*Y*Y)*W*W
	/(X*X*X*Y*Y*Y*Z*Z);
      s->c= (2.0*Z/X)+ (2.0*Z*Z*W)/(X*X*Y) + (2.0*Z*Z*Z*W*W)/(X*X*X*Y*Y);
      s->d= 2.0*( X*W - Y*Z );
      s->e= 2.0*( X*Y - Z*W );
      s->f= -(2.0*X/Z) - (2.0*X*X*W)/(Z*Z*Y) - (2.0*X*X*X*W*W)/(Z*Z*Z*Y*Y);
      s->g= (2.0-(X*X+Y*Y))/(Y*Z)
	+ ((2.0-2.0*(X*X+Y*Y)+X*X*(Y*Y+X*X))*W)/(X*Y*Y*Z*Z)
	-((-2.0+4.0*X*X+4.0*Y*Y+2.0*X*X*Y*Y*Y*Y-2.0*Y*Y*Y*Y
	   +3.0*X*X*X*X*Y*Y+X*X*X*X*X*X-4.0*X*X*X*X-6.0*X*X*Y*Y)*W*W)
	/(Z*Z*Z*X*X*Y*Y*Y) - (2.0*W*W)/(Y*Z);
      s->h= -(Y*Y+Z*Z)*(1.0/(X*Y) + (Z*W)/(X*X*Y*Y) + (Z*Z*W*W)/(X*X*X*Y*Y*Y));
    }
  }
  else if (Y==0.0) {
    if (debug) fprintf(stderr,"Special case of Y==0\n");
    s->a= -X/W;
    s->b= Z/W;
    s->c= 0.0;
    s->d= 2.0*X*W;
    s->e= -2.0*Z*W;
    s->f= 0.0;
    s->g= -X/W;
    s->h= Z/W;
  }
  else {
    /* We now do a long, painful dance through a minefield of real and
     * fake singularities.  Results are to third or fourth order in 
     * the small parameters, depending on lazyness on my part.
     */

    s->d= -2.0*S;
    
    s->e= 2.0*R;
      
    if (fabs(S)>=CANCELLATION_TOL && fabs(R)>=CANCELLATION_TOL) {
      /* Nothing is small */
      if (debug) fprintf(stderr,"Solving in general case\n");
      s->a= (X*X + Y*Y)/S;
      
      s->f= -2.0*X*Y/S;
      
      s->b= (-1.0+W*W-X*X)/R - 2.0*X*W*(X*X + Y*Y)/(R*S);
      
      s->c= 2.0*Y*Z/R;
      
      s->g= (X*X - Y*Y)/S + 2.0*X*Y*(Y*Y + Z*Z)/(R*S);
      
      s->h= -(Y*Y + Z*Z)/R;
      if (debug)
	fprintf(stderr,"Results: %f %f %f %f\n         %f %f %f %f \n",
		s->a, s->b, s->c, s->d, s->e, s->f, s->g, s->h);
    }
    else if (fabs(Y)<CANCELLATION_TOL) {
      if ((fabs(R)>=CANCELLATION_TOL || fabs(R/Y)>=CANCELLATION_TOL)
	  && (fabs(S)>=CANCELLATION_TOL || fabs(S/Y)>=CANCELLATION_TOL)) {
	/* Y is small enough to use as expansion parameter; 
	 * Expansion in Y, in terms of R and S 
	 */
	if (debug) 
	  fprintf(stderr,"Small Y expansion (Y= %f, W= %f, R= %f, S= %f)\n",
		  Y,W,R,S);
	s->a= S/(W*W) + (2.0*Y*R)/(W*W*W) + (2.0*Y*Y*S)/(W*W*W*W) 
	  + Y*(Y/S) + (Y*(Y/S)*R*R)/(W*W*W*W);

	s->b= -R/(W*W) + ((Y/R)*Y*S*S)/(W*W*W*W) + (2.0*Y*Y*R)/(W*W*W*W)
	  + Y*(Y/R) + (2.0*Y*Y*(Y/S)*R*R)/(W*W*W*W*W) 
	  + (2.0*Y*Y*(Y/S))/W;
	
	s->c= -2.0*Y/W - (2.0*Y*(Y/R)*S)/(W*W) - (2.0*Y*Y*Y)/(W*W*W);
	
	s->f= 2.0*Y/W + (2.0*Y*(Y/S)*R)/(W*W) + (2.0*Y*Y*Y)/(W*W*W);

	s->g= S/(W*W) - (2.0*Y*Y*S)/(W*W*W*W) - (Y*(Y/S)*R*R)/(W*W*W*W)
	  - Y*(Y/S) - (6.0*Y*Y*Y*R)/(W*W*W*W*W) 
	  - 2.0*Y*Y*(Y/R)*S*S/(W*W*W*W*W) - (2.0*Y*Y*(Y/R))/W;
	
	s->h= -R/(W*W) - (2.0*Y*S)/(W*W*W) - (2.0*Y*Y*R)/(W*W*W*W)
	  -(Y*(Y/R)*S*S)/(W*W*W*W) - Y*(Y/R);
	if (debug)
	  fprintf(stderr,"Results: %f %f %f %f\n         %f %f %f %f \n",
		  s->a, s->b, s->c, s->d, s->e, s->f, s->g, s->h);
      }
      else {
	/* Can't expand in Y; solution is singular */
	if (debug) fprintf(stderr,"R/Y (%f) or S/Y (%f) too small to expand in Y\n",
			   R/Y, S/Y);
	s->decomposition= SHEAR_NONE;
      }
    }
    else {
      /* No expansion; solution is singular */
      if (debug) 
	fprintf(stderr,"Y large, R (%f) or S (%f) small, so expansion fails.\n",
		R,S);
      s->decomposition= SHEAR_NONE;
    }
  }

  if (s->decomposition != SHEAR_YZXY) {
    if (debug) {
      fprintf(stderr,
	      "quat_to_shears_yzxy: decomposition of (%f %f %f %f) failed!\n",
		X,Y,Z,W);
    }
    return 0;
  }
  else return 1;
}

/* This is assumed to be symmetric in a...h and to use only absolute vals,
 * and is to return values greater than 0.0.  Small is good. 
 */
static double qualityMeasure( ShearParams* s ) 
{
  double result;
  if (s->decomposition==SHEAR_NONE) result= HUGE_VAL;
  else {
    switch (qual_measure_type) {
    case FR3D_QUAL_COX:
      {
	double tmp;
	result= fabs(s->a);
	if ((tmp=fabs(s->b))>result) result= tmp;
	if ((tmp=fabs(s->c))>result) result= tmp;
	if ((tmp=fabs(s->d))>result) result= tmp;
	if ((tmp=fabs(s->e))>result) result= tmp;
	if ((tmp=fabs(s->f))>result) result= tmp;
	if ((tmp=fabs(s->g))>result) result= tmp;
	if ((tmp=fabs(s->h))>result) result= tmp;
	if (decompIsLong(s->decomposition)) {
	  if ((tmp=fabs(s->abar))>result) result= tmp;
	  if ((tmp=fabs(s->bbar))>result) result= tmp;
	  if ((tmp=fabs(s->cbar))>result) result= tmp;
	  if ((tmp=fabs(s->dbar))>result) result= tmp;
	  if ((tmp=fabs(s->ebar))>result) result= tmp;
	  if ((tmp=fabs(s->fbar))>result) result= tmp;
	  if ((tmp=fabs(s->gbar))>result) result= tmp;
	  if ((tmp=fabs(s->hbar))>result) result= tmp;
	}
      }
      break;
    case FR3D_QUAL_SUM_ABS:
      {
	result= fabs(s->a);
	result += fabs(s->b);
	result += fabs(s->c);
	result += fabs(s->d);
	result += fabs(s->e);
	result += fabs(s->f);
	result += fabs(s->g);
	result += fabs(s->h);
	if (decompIsLong(s->decomposition)) {
	  result += fabs(s->abar);
	  result += fabs(s->bbar);
	  result += fabs(s->cbar);
	  result += fabs(s->dbar);
	  result += fabs(s->ebar);
	  result += fabs(s->fbar);
	  result += fabs(s->gbar);
	  result += fabs(s->hbar);
	}
      }
      break;
    case FR3D_QUAL_SUM_SQR:
      {
	result= s->a*s->a;
	result += s->b*s->b;
	result += s->c*s->c;
	result += s->d*s->d;
	result += s->e*s->e;
	result += s->f*s->f;
	result += s->g*s->g;
	result += s->h*s->h;
	if (decompIsLong(s->decomposition)) {
	  result += s->abar*s->abar;
	  result += s->bbar*s->bbar;
	  result += s->cbar*s->cbar;
	  result += s->dbar*s->dbar;
	  result += s->ebar*s->ebar;
	  result += s->fbar*s->fbar;
	  result += s->gbar*s->gbar;
	  result += s->hbar*s->hbar;
	}
      }
      break;
    case FR3D_QUAL_UNIT_CELL:
      {
	if (decompIsLong(s->decomposition)) {
	  result= fabs(s->a + s->g + s->abar + s->gbar);
	  result += fabs(s->b + s->h + s->bbar + s->hbar);
	  result += fabs(s->c + s->ebar);
	  result += fabs(s->d + s->fbar);
	  result += fabs(s->e + s->cbar);
	  result += fabs(s->f + s->dbar);
	}
	else {
	  result= fabs(s->a + s->g);
	  result += fabs(s->b + s->h);
	  result += fabs(s->c);
	  result += fabs(s->d);
	  result += fabs(s->e);
	  result += fabs(s->f);
	}
      }
      break;
    default:
      {
	fprintf(stderr,"Unknown quality measure type %d!\n",
		qual_measure_type);
	result= HUGE_VAL;
      }
    }
  }
  if (debug) fprintf(stderr,"qualityMeasure (%s): %f\n",
		     qualMeasureName(qual_measure_type), result);
  s->quality= result;
  s->qualType= qual_measure_type;
  return result;
}

/* non-zero on success; tries for a specific decomposition */
static int quat_to_shears( int decomp, Quat *q, ShearParams *s )
{
  Quat test_q;
  ShearParams test_s;

  s->decomposition= decomp;
  s->qualType= FR3D_QUAL_UNSET;

  switch (decomp) {
  case SHEAR_YZXY:
    {
      return(quat_to_shears_yzxy(q,s));
    }
    /* break; -NOTREACHED- */

  case SHEAR_YXZY:
    {
      test_q.x= q->x; test_q.y= q->y; test_q.z= q->z; test_q.w= -q->w;
      if (quat_to_shears_yzxy(&test_q,&test_s)) {
	s->a= -test_s.g; s->b= -test_s.h;
	s->c= -test_s.e; s->d= -test_s.f;
	s->e= -test_s.c; s->f= -test_s.d;
	s->g= -test_s.a; s->h= -test_s.b;
	return 1;
      }
      else {
	s->decomposition= SHEAR_NONE;
	return 0;
      }
    }
    /* break; -NOTREACHED- */

  case SHEAR_ZXYZ:
    {
      test_q.x= q->y; test_q.y= q->z; test_q.z= q->x; test_q.w= q->w;
      if (quat_to_shears_yzxy(&test_q,&test_s)) {
	s->a= test_s.a; s->b= test_s.b;
	s->c= test_s.c; s->d= test_s.d;
	s->e= test_s.e; s->f= test_s.f;
	s->g= test_s.g; s->h= test_s.h;
	return 1;
      }
      else {
	s->decomposition= SHEAR_NONE;
	return 0;
      }
    }
    /* break; -NOTREACHED- */

  case SHEAR_ZYXZ:
    {
      test_q.x= q->y; test_q.y= q->z; test_q.z= q->x; test_q.w= -q->w;
      if (quat_to_shears_yzxy(&test_q,&test_s)) {
	s->a= -test_s.g; s->b= -test_s.h;
	s->c= -test_s.e; s->d= -test_s.f;
	s->e= -test_s.c; s->f= -test_s.d;
	s->g= -test_s.a; s->h= -test_s.b;
	return 1;
      }
      else {
	s->decomposition= SHEAR_NONE;
	return 0;
      }
    }
    /* break; -NOTREACHED- */

  case SHEAR_XYZX:
    {
      test_q.x= q->z; test_q.y= q->x; test_q.z= q->y; test_q.w= q->w;
      if (quat_to_shears_yzxy(&test_q,&test_s)) {
	s->a= test_s.a; s->b= test_s.b;
	s->c= test_s.c; s->d= test_s.d;
	s->e= test_s.e; s->f= test_s.f;
	s->g= test_s.g; s->h= test_s.h;
	return 1;
      }
      else {
	s->decomposition= SHEAR_NONE;
	return 0;
      }
    }
    /* break; -NOTREACHED- */

  case SHEAR_XZYX:
    {
      test_q.x= q->z; test_q.y= q->x; test_q.z= q->y; test_q.w= -q->w;
      if (quat_to_shears_yzxy(&test_q,&test_s)) {
	s->a= -test_s.g; s->b= -test_s.h;
	s->c= -test_s.e; s->d= -test_s.f;
	s->e= -test_s.c; s->f= -test_s.d;
	s->g= -test_s.a; s->h= -test_s.b;
	return 1;
      }
      else {
	s->decomposition= SHEAR_NONE;
	return 0;
      }
    }
    /* break; -NOTREACHED- */

  case SHEAR_YZXYXZY:
    {
      Quat qlcl;
      Quat half;
      ShearParams fwd, bkwd;

      qlcl= *q;
      if (qlcl.w < 0.0) {
	/* We want to deal with negative of this quat, to minimize angles */
	qlcl.x *= -1.0;
	qlcl.y *= -1.0;
	qlcl.z *= -1.0;
	qlcl.w *= -1.0;
      }

      quat_nrm_sqrt(&half,&qlcl);
      if (quat_to_shears(SHEAR_YZXY,&half,&fwd)
	  && quat_to_shears(SHEAR_YXZY,&half,&bkwd)) {
	s->a= fwd.a; s->b= fwd.b; s->c= fwd.c; s->d= fwd.d; 
	s->e= fwd.e; s->f= fwd.f; s->g= fwd.g; s->h= fwd.h; 
	s->abar= bkwd.a; s->bbar= bkwd.b; s->cbar= bkwd.c; s->dbar= bkwd.d; 
	s->ebar= bkwd.e; s->fbar= bkwd.f; s->gbar= bkwd.g; s->hbar= bkwd.h; 
	return 1;
      }
      else {
	s->decomposition= SHEAR_NONE;
	return 0;
      }
    }
    /* break; -NOTREACHED- */

  case SHEAR_ZXYZYXZ:
    {
      Quat qlcl;
      Quat half;
      ShearParams fwd, bkwd;

      qlcl= *q;
      if (qlcl.w < 0.0) {
	/* We want to deal with negative of this quat, to minimize angles */
	qlcl.x *= -1.0;
	qlcl.y *= -1.0;
	qlcl.z *= -1.0;
	qlcl.w *= -1.0;
      }

      quat_nrm_sqrt(&half,&qlcl);
      if (quat_to_shears(SHEAR_ZXYZ,&half,&fwd)
	  && quat_to_shears(SHEAR_ZYXZ,&half,&bkwd)) {
	s->a= fwd.a; s->b= fwd.b; s->c= fwd.c; s->d= fwd.d; 
	s->e= fwd.e; s->f= fwd.f; s->g= fwd.g; s->h= fwd.h; 
	s->abar= bkwd.a; s->bbar= bkwd.b; s->cbar= bkwd.c; s->dbar= bkwd.d; 
	s->ebar= bkwd.e; s->fbar= bkwd.f; s->gbar= bkwd.g; s->hbar= bkwd.h; 
	return 1;
      }
      else {
	s->decomposition= SHEAR_NONE;
	return 0;
      }
    }
    /* break; -NOTREACHED- */

  case SHEAR_XYZXZYX:
    {
      Quat qlcl;
      Quat half;
      ShearParams fwd, bkwd;

      qlcl= *q;
      if (qlcl.w < 0.0) {
	/* We want to deal with negative of this quat, to minimize angles */
	qlcl.x *= -1.0;
	qlcl.y *= -1.0;
	qlcl.z *= -1.0;
	qlcl.w *= -1.0;
      }

      quat_nrm_sqrt(&half,&qlcl);
      if (quat_to_shears(SHEAR_XYZX,&half,&fwd)
	  && quat_to_shears(SHEAR_XZYX,&half,&bkwd)) {
	s->a= fwd.a; s->b= fwd.b; s->c= fwd.c; s->d= fwd.d; 
	s->e= fwd.e; s->f= fwd.f; s->g= fwd.g; s->h= fwd.h; 
	s->abar= bkwd.a; s->bbar= bkwd.b; s->cbar= bkwd.c; s->dbar= bkwd.d; 
	s->ebar= bkwd.e; s->fbar= bkwd.f; s->gbar= bkwd.g; s->hbar= bkwd.h; 
	return 1;
      }
      else {
	s->decomposition= SHEAR_NONE;
	return 0;
      }
    }
    /* break; -NOTREACHED- */

  default:
  case SHEAR_NONE:
    {
      fprintf(stderr,"Unknown shear decomposition in quat_to_shears!\n");
      return 0;
    }
    /* break; -NOTREACHED- */
  }
#if ( ! ( SGI5 || SGI6 || SGI64 || SGIMP64 ) )
  return 0; /* not reached */
#endif
}

/* non-zero on success */
static int quat_to_shears_any( Quat *q, ShearParams *s ) 
{
  Quat test_q;
  ShearParams test_s;

  (void)quat_to_shears( SHEAR_YZXY, q, &test_s );
  qualityMeasure( &test_s );
  *s= test_s;

  if (quat_to_shears(SHEAR_YXZY,q,&test_s)) {
    if (qualityMeasure(&test_s)<s->quality) *s= test_s;
  }

  if (quat_to_shears(SHEAR_ZXYZ,q,&test_s)) {
    if (qualityMeasure(&test_s)<s->quality) *s= test_s;
  }

  if (quat_to_shears(SHEAR_ZYXZ,q,&test_s)) {
    if (qualityMeasure(&test_s)<s->quality) *s= test_s;
  }

  if (quat_to_shears(SHEAR_XYZX,q,&test_s)) {
    if (qualityMeasure(&test_s)<s->quality) *s= test_s;
  }

  if (quat_to_shears(SHEAR_XZYX,q,&test_s)) {
    if (qualityMeasure(&test_s)<s->quality) *s= test_s;
  }

  return (s->decomposition != SHEAR_NONE);
}

/* non-zero on success */
static int quat_to_shears_long_any( Quat *q, ShearParams *s )
{
  Quat test_q;
  ShearParams test_s;

  (void)quat_to_shears( SHEAR_YZXYXZY, q, &test_s );
  qualityMeasure( &test_s );
  *s= test_s;

  if (quat_to_shears(SHEAR_ZXYZYXZ,q,&test_s)) {
    if (qualityMeasure(&test_s)<s->quality) *s= test_s;
  }

  if (quat_to_shears(SHEAR_XYZXZYX,q,&test_s)) {
    if (qualityMeasure(&test_s)<s->quality) *s= test_s;
  }

  return (s->decomposition != SHEAR_NONE);
}

static int safe_quat( Quat* q ) {
  /* Returns true if we are well away from singularities in the
   * demoninators of the shear coefficient equations.
   */
  ShearParams junk;
  return quat_to_shears_any( q, &junk );
}

static int safe_quat_long( Quat* q ) {
  /* Returns true if we are well away from singularities in the
   * demoninators of the shear coefficient equations.
   */
  ShearParams junk;
  return quat_to_shears_long_any( q, &junk );
}

static void trans_shear_adjust(ShearParams* s, TransParams* t_in,
			       TransParams* t_out)
{
  switch (s->decomposition) 
    {
    case SHEAR_YZXYXZY:
    case SHEAR_YZXY: {
      t_out->dx= t_in->dx;
      t_out->dy= t_in->dy - s->b*t_in->dx - s->a*t_in->dz;
      t_out->dz= t_in->dz - s->c*t_in->dx;
    }
    break;
    case SHEAR_YXZY: {
      t_out->dx= t_in->dx - s->d*t_in->dz;
      t_out->dy= t_in->dy - s->b*t_in->dx - s->a*t_in->dz;
      t_out->dz= t_in->dz;
    }
    break;

    case SHEAR_ZXYZYXZ:
    case SHEAR_ZXYZ: {
      t_out->dy= t_in->dy;
      t_out->dz= t_in->dz - s->b*t_in->dy - s->a*t_in->dx;
      t_out->dx= t_in->dx - s->c*t_in->dy;
    }
    break;
    case SHEAR_ZYXZ: {
      t_out->dy= t_in->dy - s->d*t_in->dx;
      t_out->dz= t_in->dz - s->b*t_in->dy - s->a*t_in->dx;
      t_out->dx= t_in->dx;
    }
    break;

    case SHEAR_XYZXZYX:
    case SHEAR_XYZX: {
      t_out->dz= t_in->dz;
      t_out->dx= t_in->dx - s->b*t_in->dz - s->a*t_in->dy;
      t_out->dy= t_in->dy - s->c*t_in->dz;
    }
    break;
    case SHEAR_XZYX: {
      t_out->dz= t_in->dz - s->d*t_in->dy;
      t_out->dx= t_in->dx - s->b*t_in->dz - s->a*t_in->dy;
      t_out->dy= t_in->dy;
    }
    break;

    case SHEAR_NONE:
    default:
      {
	Abort("fshrot3d internal error: unknown decomposition (%d) in trans_shear_adjust!\n",
	      (int)s->decomposition);
      }
    }
}

static void step_repairing_quat( Quat* result, Quat* q, int step ) {
  Quat qstep;
  Quat tmp;

  if (debug) fprintf(stderr,"Repair step %d\n",step);

  switch (step % 3) {
  case 0: quat_create(&qstep, REPAIR_QUAT_STEPSIZE, 0.0, 0.0, 1); break;
  case 1: quat_create(&qstep, 0.0, REPAIR_QUAT_STEPSIZE, 0.0, 1); break;
  case 2: quat_create(&qstep, 0.0, 0.0, REPAIR_QUAT_STEPSIZE, 1); break;
  }

  quat_mult_right(result, &qstep);
  if (debug) fprintf(stderr,"Trying (%f %f %f %f) with (%f %f %f %f)\n",
		     result->x, result->y, result->z, result->w,
		     q->x, q->y, q->z, q->w);
  quat_copy(&tmp, result);
  quat_mult_right(&tmp, q);
  if (safe_quat(&tmp)) return;
  else if (step<REPAIR_QUAT_MAX_STEPS) 
    step_repairing_quat( result, q, step+1 );
  else 
    Abort(
  "step_repairing_quat: could not find good soln for quaternion %f %f %f %f!\n",
	  q->x,q->y,q->z,q->w);
}

static void find_repairing_quat( Quat* result, Quat* q ) {
  quat_identity(result);
  step_repairing_quat( result, q, 0 );
}

static void step_repairing_quat_long( Quat* result, Quat* q, int step ) {
  Quat qstep;
  Quat tmp;

  if (debug) fprintf(stderr,"Repair step %d\n",step);

  switch (step % 3) {
  case 0: quat_create(&qstep, REPAIR_QUAT_STEPSIZE, 0.0, 0.0, 1); break;
  case 1: quat_create(&qstep, 0.0, REPAIR_QUAT_STEPSIZE, 0.0, 1); break;
  case 2: quat_create(&qstep, 0.0, 0.0, REPAIR_QUAT_STEPSIZE, 1); break;
  }

  quat_mult_right(result, &qstep);
  if (debug) fprintf(stderr,"Trying (%f %f %f %f) with (%f %f %f %f)\n",
		     result->x, result->y, result->z, result->w,
		     q->x, q->y, q->z, q->w);
  quat_copy(&tmp, result);
  quat_mult_right(&tmp, q);
  if (safe_quat_long(&tmp)) return;
  else if (step<REPAIR_QUAT_MAX_STEPS) 
    step_repairing_quat( result, q, step+1 );
  else 
    Abort(
  "step_repairing_quat_long: could not find good soln for quaternion %f %f %f %f!",
	  q->x,q->y,q->z,q->w);
}

static void find_repairing_quat_long( Quat* result, Quat* q ) {
  quat_identity(result);
  step_repairing_quat_long( result, q, 0 );
}

static void shear_x( FComplex* image, double a, double b, double delta,
		     long nx, long ny, long nz,
		     double length_x, double length_y, double length_z,
		     int real_flag )
{
  FComplex* here;
  long halfx= nx/2;
  long halfy= ny/2;
  long halfz= nz/2;
  int i;
  int j;
  int k;
  int lowerbound;
  double theta;
  double x_phase_scale;
  double y_phase_scale;
  double z_phase_scale;
  long nx_mod;
  long ny_mod;
  long nz_mod;

  /* Return if there is nothing to do */
  if (a==0.0 && b==0.0 && delta==0.0) return;

  if (debug) fprintf(stderr,"shear_x %f %f %f\n",a,b,delta);

  /* Adjustments for odd data dimensions */
  if (nx>1 && (nx % 2)) nx_mod= nx-1;
  else nx_mod= nx;
  if (ny>1 && (ny % 2)) ny_mod= ny-1;
  else ny_mod= ny;
  if (nz>1 && (nz % 2)) nz_mod= nz-1;
  else nz_mod= nz;

  /* Step counter */
  count_shear_x++;

  /* FFT in x */
  fft3d( image, nx, ny, nz, +1, "x" );

  /* delta is in voxels, a and b in fractional shears 
   * (typical range -1 to 1).
   */
  x_phase_scale= 2.0*M_PI*delta/((double)nx_mod);
  y_phase_scale= 2.0*M_PI*a*length_y/(length_x*ny_mod);
  z_phase_scale= 2.0*M_PI*b*length_z/(length_x*nz_mod);

  /* Apply phase changes */
  if (real_flag && !(nx%2)) lowerbound= 1; /* no phase shift at Nyquist freq */
  else lowerbound= 0;
  for (i=lowerbound; i<nx; i++) {
    for (j=0; j<ny; j++) {
      for (k=0; k<nz; k++) {
	FComplex t;
	float c;
	float s;
	here= &(MEM(image,nx,ny,nz,i,j,k));
	t.real= here->real;
	t.imag= here->imag;
	/* Signs of terms are determined by relationship between data
	 * coordinate system and geometric coordinate system.
	 */
	theta= (i-halfx)*x_phase_scale 
	  - (j-halfy)*(i-halfx)*y_phase_scale
	  + (k-halfz)*(i-halfx)*z_phase_scale;
	c= cos(theta);
	s= sin(theta);
	here->real= c*t.real - s*t.imag;
	here->imag= c*t.imag + s*t.real;
      }
    }
  }

  /* FFT back in x */
  fft3d( image, nx, ny, nz, -1, "x" );

}

static void shear_y( FComplex* image, double a, double b, double delta,
		     long nx, long ny, long nz,
		     double length_x, double length_y, double length_z,
		     int real_flag )
{
  FComplex* here;
  long halfx= nx/2;
  long halfy= ny/2;
  long halfz= nz/2;
  int i;
  int j;
  int k;
  int lowerbound;
  double theta;
  double x_phase_scale;
  double y_phase_scale;
  double z_phase_scale;
  long nx_mod;
  long ny_mod;
  long nz_mod;

  /* Return if there is nothing to do */
  if (a==0.0 && b==0.0 && delta==0.0) return;

  if (debug) fprintf(stderr,"shear_y %f %f %f\n",a,b,delta);

  /* Adjustments for odd data dimensions */
  if (nx>1 && (nx % 2)) nx_mod= nx-1;
  else nx_mod= nx;
  if (ny>1 && (ny % 2)) ny_mod= ny-1;
  else ny_mod= ny;
  if (nz>1 && (nz % 2)) nz_mod= nz-1;
  else nz_mod= nz;

  /* Step counter */
  count_shear_y++;

  /* FFT in y */
  fft3d( image, nx, ny, nz, +1, "y" );

  /* delta is in voxels, a and b in fractional shears 
   * (typical range -1 to 1).
   */
  y_phase_scale= 2.0*M_PI*delta/((double)ny_mod);
  z_phase_scale= 2.0*M_PI*a*length_z/(length_y*nz_mod);
  x_phase_scale= 2.0*M_PI*b*length_x/(length_y*nx_mod);

  /* Apply phase changes */
  if (real_flag && !(ny%2)) lowerbound= 1; /* no phase shift at Nyquist freq */
  else lowerbound= 0;
  for (i=0; i<nx; i++) {
    for (j=lowerbound; j<ny; j++) {
      for (k=0; k<nz; k++) {
	FComplex t;
	float c;
	float s;
	here= &(MEM(image,nx,ny,nz,i,j,k));
	t.real= here->real;
	t.imag= here->imag;
	/* Signs of terms are determined by relationship between data
	 * coordinate system and geometric coordinate system.
	 */
	theta= -(j-halfy)*y_phase_scale 
	  - (k-halfz)*(j-halfy)*z_phase_scale
	  - (i-halfx)*(j-halfy)*x_phase_scale;
	c= cos(theta);
	s= sin(theta);
	here->real= c*t.real - s*t.imag;
	here->imag= c*t.imag + s*t.real;
      }
    }
  }

  /* FFT back in y */
  fft3d( image, nx, ny, nz, -1, "y" );
}

static void shear_z( FComplex* image, double a, double b, double delta,
		     long nx, long ny, long nz,
		     double length_x, double length_y, double length_z,
		     int real_flag )
{
  FComplex* here;
  long halfx= nx/2;
  long halfy= ny/2;
  long halfz= nz/2;
  int i;
  int j;
  int k;
  int lowerbound;
  double theta;
  double x_phase_scale;
  double y_phase_scale;
  double z_phase_scale;
  long nx_mod;
  long ny_mod;
  long nz_mod;

  /* Return if there is nothing to do */
  if (a==0.0 && b==0.0 && delta==0.0) return;

  if (debug) fprintf(stderr,"shear_z %f %f %f\n",a,b,delta);

  /* Adjustments for odd data dimensions */
  if (nx>1 && (nx % 2)) nx_mod= nx-1;
  else nx_mod= nx;
  if (ny>1 && (ny % 2)) ny_mod= ny-1;
  else ny_mod= ny;
  if (nz>1 && (nz % 2)) nz_mod= nz-1;
  else nz_mod= nz;

  /* Step counter */
  count_shear_z++;

  /* FFT in z */
  fft3d( image, nx, ny, nz, +1, "z" );

  /* delta is in voxels, a and b in fractional shears 
   * (typical range -1 to 1).
   */
  z_phase_scale= 2.0*M_PI*delta/((double)nz_mod);
  x_phase_scale= 2.0*M_PI*a*length_x/(length_z*nx_mod);
  y_phase_scale= 2.0*M_PI*b*length_y/(length_z*ny_mod);

  /* Apply phase changes */
  if (real_flag && !(nz%2)) lowerbound= 1; /* no phase shift at Nyquist freq */
  else lowerbound= 0;
  for (i=0; i<nx; i++) {
    for (j=0; j<ny; j++) {
      for (k=lowerbound; k<nz; k++) {
	FComplex t;
	float c;
	float s;
	here= &(MEM(image,nx,ny,nz,i,j,k));
	t.real= here->real;
	t.imag= here->imag;
	/* Signs of terms are determined by relationship between data
	 * coordinate system and geometric coordinate system.
	 */
	theta= (k-halfz)*z_phase_scale 
	  + (i-halfx)*(k-halfz)*x_phase_scale
	  - (j-halfy)*(k-halfz)*y_phase_scale;
	c= cos(theta);
	s= sin(theta);
	here->real= c*t.real - s*t.imag;
	here->imag= c*t.imag + s*t.real;
      }
    }
  }

  /* FFT back in z */
  fft3d( image, nx, ny, nz, -1, "z" );
}

static void shift_only( TransParams* t, FComplex* moved_image, 
			long nx, long ny, long nz,
			double length_x, double length_y, double length_z,
			int real_flag )
{
  /* Easy as pie, since pure shifts commute */
  shear_x(moved_image, 0.0, 0.0, t->dx,
	  nx, ny, nz, length_x, length_y, length_z, real_flag);
  shear_y(moved_image, 0.0, 0.0, t->dy,
	  nx, ny, nz, length_x, length_y, length_z, real_flag);
  shear_z(moved_image, 0.0, 0.0, t->dz,
	  nx, ny, nz, length_x, length_y, length_z, real_flag);
}

static void rot_x( Quat* q, TransParams* t, FComplex* moved_image, 
		   long nx, long ny, long nz, 
		   double length_x, double length_y, double length_z,
		   int real_flag) 
{
  TransParams adjusted_trans;
  double alpha;
  double delta;

  alpha= -(q->x/q->w);
  delta= 2.0*q->x*q->w;
  adjusted_trans.dx= t->dx;
  adjusted_trans.dy= t->dy - alpha*t->dz;
  adjusted_trans.dz= t->dz;

  shear_y(moved_image, alpha, 0.0, 0.0, nx, ny, nz, 
	  length_x, length_y, length_z, real_flag);
  shear_z(moved_image, 0.0, delta, adjusted_trans.dz, nx, ny, nz, 
	  length_x, length_y, length_z, real_flag);
  shear_y(moved_image, alpha, 0.0, adjusted_trans.dy, nx, ny, nz, 
	  length_x, length_y, length_z, real_flag);
  shear_x(moved_image, 0.0, 0.0, adjusted_trans.dx, nx, ny, nz, 
	  length_x, length_y, length_z, real_flag);
}

static void rot_y( Quat* q, TransParams* t, FComplex* moved_image, 
		   long nx, long ny, long nz, 
		   double length_x, double length_y, double length_z,
		   int real_flag ) 
{
  TransParams adjusted_trans;
  double alpha;
  double delta;

  alpha= -(q->y/q->w);
  delta= 2.0*q->y*q->w;
  adjusted_trans.dx= t->dx;
  adjusted_trans.dy= t->dy;
  adjusted_trans.dz= t->dz - alpha*t->dx;

  shear_z(moved_image, alpha, 0.0, 0.0, nx, ny, nz, 
	  length_x, length_y, length_z, real_flag);
  shear_x(moved_image, 0.0, delta, adjusted_trans.dx, nx, ny, nz, 
	  length_x, length_y, length_z, real_flag);
  shear_z(moved_image, alpha, 0.0, adjusted_trans.dz, nx, ny, nz, 
	  length_x, length_y, length_z, real_flag);
  shear_y(moved_image, 0.0, 0.0, adjusted_trans.dy, nx, ny, nz, 
	  length_x, length_y, length_z, real_flag);
}

static void rot_z( Quat* q, TransParams* t, FComplex* moved_image, 
		   long nx, long ny, long nz, 
		   double length_x, double length_y, double length_z,
		   int real_flag ) 
{
  TransParams adjusted_trans;
  double alpha;
  double delta;

  alpha= -(q->z/q->w);
  delta= 2.0*q->z*q->w;
  adjusted_trans.dx= t->dx - alpha*t->dy;
  adjusted_trans.dy= t->dy;
  adjusted_trans.dz= t->dz;

  shear_x(moved_image, alpha, 0.0, 0.0, nx, ny, nz, 
	  length_x, length_y, length_z, real_flag);
  shear_y(moved_image, 0.0, delta, adjusted_trans.dy, nx, ny, nz, 
	  length_x, length_y, length_z, real_flag);
  shear_x(moved_image, alpha, 0.0, adjusted_trans.dx, nx, ny, nz, 
	  length_x, length_y, length_z, real_flag);
  shear_z(moved_image, 0.0, 0.0, adjusted_trans.dz, nx, ny, nz, 
	  length_x, length_y, length_z, real_flag);
}

static void apply_shear_decomp( FComplex* image,
				ShearParams* s, TransParams* t,
				long nx, long ny, long nz,
				double length_x, double length_y, 
				double length_z, int real_flag )
{
  TransParams adjusted_trans;

  /* Adjust the translation for the shears */
  trans_shear_adjust(s, t, &adjusted_trans);

  /* Tell any intrepid debuggers the plan */
  if (debug) writeShearTrans(s, &adjusted_trans);

  /* Do the shears in place */
  switch (s->decomposition) {
  case SHEAR_YZXY: {
    shear_y(image, s->g, s->h, 0.0, nx, ny, nz,
	    length_x, length_y, length_z, real_flag);
    shear_x(image, s->e, s->f, adjusted_trans.dx,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
    shear_z(image, s->c, s->d, adjusted_trans.dz,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
    shear_y(image, s->a, s->b, adjusted_trans.dy,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
  }
  break;
  case SHEAR_YXZY: {
    shear_y(image, s->g, s->h, 0.0,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
    shear_z(image, s->e, s->f, adjusted_trans.dz,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
    shear_x(image, s->c, s->d, adjusted_trans.dx,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
    shear_y(image, s->a, s->b, adjusted_trans.dy,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
  }
  break;
  
  case SHEAR_ZXYZ: {
    shear_z(image, s->g, s->h, 0.0, nx, ny, nz,
	    length_x, length_y, length_z, real_flag);
    shear_y(image, s->e, s->f, adjusted_trans.dy,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
    shear_x(image, s->c, s->d, adjusted_trans.dx,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
    shear_z(image, s->a, s->b, adjusted_trans.dz,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
  }
  break;
  case SHEAR_ZYXZ: {
    shear_z(image, s->g, s->h, 0.0,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
    shear_x(image, s->e, s->f, adjusted_trans.dx,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
    shear_y(image, s->c, s->d, adjusted_trans.dy,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
    shear_z(image, s->a, s->b, adjusted_trans.dz, 
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
  }
  break;
  
  case SHEAR_XYZX: {
    shear_x(image, s->g, s->h, 0.0, nx, ny, nz,
	    length_x, length_y, length_z, real_flag);
    shear_z(image, s->e, s->f, adjusted_trans.dz,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
    shear_y(image, s->c, s->d, adjusted_trans.dy,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
    shear_x(image, s->a, s->b, adjusted_trans.dx,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
  }
  break;
  case SHEAR_XZYX: {
    shear_x(image, s->g, s->h, 0.0,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
    shear_y(image, s->e, s->f, adjusted_trans.dy,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
    shear_z(image, s->c, s->d, adjusted_trans.dz,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
    shear_x(image, s->a, s->b, adjusted_trans.dx, 
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
  }
  break;
  
  case SHEAR_YZXYXZY: {
    shear_y(image, s->gbar, s->hbar, 0.0,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
    shear_z(image, s->ebar, s->fbar, 0.0,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
    shear_x(image, s->cbar, s->dbar, 0.0,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
    shear_y(image, s->g+s->abar, 
	    s->h+s->bbar, 0.0, nx, ny, nz,
	    length_x, length_y, length_z, real_flag);
    shear_x(image, s->e, s->f, adjusted_trans.dx,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
    shear_z(image, s->c, s->d, adjusted_trans.dz,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
    shear_y(image, s->a, s->b, adjusted_trans.dy,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
  }
  break;
  
  case SHEAR_ZXYZYXZ: {
    shear_z(image, s->gbar, s->hbar, 0.0,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
    shear_x(image, s->ebar, s->fbar, 0.0,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
    shear_y(image, s->cbar, s->dbar, 0.0,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
    shear_z(image, s->g+s->abar, 
	    s->h+s->bbar, 0.0, nx, ny, nz,
	    length_x, length_y, length_z, real_flag);
    shear_y(image, s->e, s->f, adjusted_trans.dy,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
    shear_x(image, s->c, s->d, adjusted_trans.dx,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
    shear_z(image, s->a, s->b, adjusted_trans.dz,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
  }
  break;
  
  case SHEAR_XYZXZYX: {
    shear_x(image, s->gbar, s->hbar, 0.0,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
    shear_y(image, s->ebar, s->fbar, 0.0,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
    shear_z(image, s->cbar, s->dbar, 0.0,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
    shear_x(image, s->g+s->abar, 
	    s->h+s->bbar, 0.0, nx, ny, nz,
	    length_x, length_y, length_z, real_flag);
    shear_z(image, s->e, s->f, adjusted_trans.dz,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
    shear_y(image, s->c, s->d, adjusted_trans.dy,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
    shear_x(image, s->a, s->b, adjusted_trans.dx,
	    nx, ny, nz, length_x, length_y, length_z, real_flag);
  }
  break;
  
  default:
    Abort("fshrot3d: apply_shear_decomp: invalid internal decomposition %d!\n",
	  s->decomposition);
  }
}

static void rot_13_shears( Quat* q, double dx, double dy, double dz,
			   FComplex* image,
			   long nx, long ny, long nz,
			   double length_x, double length_y, double length_z,
			   int real_flag )
{
  double theta_x;
  double theta_y;
  double theta_z;
  double xa, xb;
  double ya, yb;
  double za, zb;
  
  if (dx != 0.0 || dy != 0.0 || dz != 0.0)
    Abort("fshrot3d: translations in 13-shear mode not yet implemented!\n");
  
  if (!euler5_decompose(q, &theta_x, &theta_y, &theta_z))
    Abort("fshrot3d: unable to decompose (%g, %g, %g, %g) into 5 Euler angles!\n",
	  q->x,q->y,q->z,q->w);
  if (debug) fprintf(stderr,
		     "Decomposed (%g, %g, %g, %g) into 5 rots: %g %g %g\n",
		     q->x,q->y,q->z,q->w,theta_x,theta_y,theta_z);

  xa= -tan( 0.5*theta_x );
  xb= sin( theta_x );
  ya= -tan( 0.5*theta_y );
  yb= sin( theta_y );
  za= -tan( 0.5*theta_z );
  zb= sin( theta_z );

  shear_x(image, za, 0.0, 0.0, nx, ny, nz, 
	  length_x, length_y, length_z, real_flag);
  shear_y(image, 0.0, zb, 0.0, nx, ny, nz, 
	  length_x, length_y, length_z, real_flag);
  shear_x(image, za, 0.0, 0.0, nx, ny, nz, 
	  length_x, length_y, length_z, real_flag);
  shear_z(image, ya, 0.0, 0.0, nx, ny, nz, 
	  length_x, length_y, length_z, real_flag);
  shear_x(image, 0.0, yb, 0.0, nx, ny, nz, 
	  length_x, length_y, length_z, real_flag);
  shear_z(image, ya, -xa, 0.0, nx, ny, nz, 
	  length_x, length_y, length_z, real_flag);
  shear_y(image, -xb, 0.0, 0.0, nx, ny, nz, 
	  length_x, length_y, length_z, real_flag);
  shear_z(image, ya, -xa, 0.0, nx, ny, nz, 
	  length_x, length_y, length_z, real_flag);
  shear_x(image, 0.0, yb, 0.0, nx, ny, nz, 
	  length_x, length_y, length_z, real_flag);
  shear_z(image, ya, 0.0, 0.0, nx, ny, nz, 
	  length_x, length_y, length_z, real_flag);
  shear_x(image, za, 0.0, 0.0, nx, ny, nz, 
	  length_x, length_y, length_z, real_flag);
  shear_y(image, 0.0, zb, 0.0, nx, ny, nz, 
	  length_x, length_y, length_z, real_flag);
  shear_x(image, za, 0.0, 0.0, nx, ny, nz, 
	  length_x, length_y, length_z, real_flag);
}

static void rot_7_shears( Quat* q, double dx, double dy, double dz,
			  FComplex* image,
			  long nx, long ny, long nz,
			  double length_x, double length_y, double length_z,
			  int real_flag )
{
  /* This routine implements the general case rotation with 7 shears */
  ShearParams q_params;
  TransParams trans;
  
  /* Find all the right parameters */
  if (quat_to_shears_long_any(q, &q_params)) {
    trans.dx= dx;
    trans.dy= dy;
    trans.dz= dz;
    apply_shear_decomp( image, &q_params, &trans, nx, ny, nz,
			length_x, length_y, length_z, real_flag );
  }
  else {
    /*
     * Found trouble, so generate a repairing quaternion and its inverse, 
     * and frame the shift and rotate with them.  This results in 8 shears.
     */
    Quat repairing_quat;
    Quat repaired_q;
    Quat repairing_conj;
    
    if (debug) fprintf(stderr,"Repair is needed; becoming recursive\n");
    
    /* Find the repairing quaternion, and mix it in with q */
    find_repairing_quat_long(&repairing_quat, q);
    quat_copy(&repaired_q,&repairing_quat);
    quat_mult_right(&repaired_q,q);
    quat_copy(&repairing_conj, &repairing_quat);
    quat_conjugate(&repairing_conj);
    
    rot_7_shears(&repaired_q, 0.0, 0.0, 0.0, image, nx, ny, nz,
		 length_x, length_y, length_z, real_flag);
    rot_7_shears(&repairing_conj, dx, dy, dz, image, nx, ny, nz,
		 length_x, length_y, length_z, real_flag);
  }
}

/* Input and output data are assumed to be ordered such that z is fastest
 * in memory.  (Note that other expansions are used for special
 * cases).  The convention is that the rotation is applied *before*
 * the shift; that is, R = T*r where R is the whole transformation,
 * T is the translation (shift), and r is the pure rotation given by
 * the input quaternion.
 */
static void rot_4_shears( Quat* q, double dx, double dy, double dz,
			  FComplex* image,
			  long nx, long ny, long nz,
			  double length_x, double length_y, double length_z,
			  int real_flag )
{
  /* This routine implements the general case rotation with 4 shears */
    ShearParams q_params;
    TransParams trans;
    TransParams adjusted_trans;

    /* Find all the right parameters */
    if (quat_to_shears_any(q, &q_params)) {
      trans.dx= dx;
      trans.dy= dy;
      trans.dz= dz;
      apply_shear_decomp( image, &q_params, &trans, nx, ny, nz,
			  length_x, length_y, length_z, real_flag );
    }
    else {
      /*
       * Found trouble, so generate a repairing quaternion and its inverse, 
       * and frame the shift and rotate with them.  This results in 8 shears.
       */
      Quat repairing_quat;
      Quat repaired_q;
      Quat repairing_conj;
      
      if (debug) fprintf(stderr,"Repair is needed; becoming recursive\n");
      
      /* Find the repairing quaternion, and mix it in with q */
      find_repairing_quat(&repairing_quat, q);
      quat_copy(&repaired_q,&repairing_quat);
      quat_mult_right(&repaired_q,q);
      quat_copy(&repairing_conj, &repairing_quat);
      quat_conjugate(&repairing_conj);

      rot_4_shears(&repaired_q, 0.0, 0.0, 0.0, image, nx, ny, nz,
		   length_x, length_y, length_z, real_flag);
      rot_4_shears(&repairing_conj, dx, dy, dz, image, nx, ny, nz,
		   length_x, length_y, length_z, real_flag);
    }
}

/* This routine applies appropriate phase shifts to the given
 * image to cause it to be shifted by the given distances the
 * next time it undergoes a full 3D FFT.  If kspace_flag is
 * set, the coming FFT is expected to be an inverse FFT; 
 * otherwise a forward FFT is assumed.
 */
void fshrot3d_set_shift_phases( double dx, double dy, double dz,
				FComplex* image,
				long nx, long ny, long nz, double length_x, 
				double length_y, double length_z,
				int kspace_flag, int real_flag )
{
  FComplex* here;
  long halfx= nx/2;
  long halfy= ny/2;
  long halfz= nz/2;
  int i;
  int j;
  int k;
  int lowerx, lowery, lowerz;
  double x_phase_scale;
  double y_phase_scale;
  double z_phase_scale;
  double theta;
  double c,s ;
  FComplex tc;
  long nx_mod;
  long ny_mod;
  long nz_mod;

  /* Return if there is nothing to do */
  if (dx==0.0 && dy==0.0 && dz==0.0) return;

  if (debug) fprintf(stderr,"fshrot3d_set_shift_phases %f %f %f\n",
		     dx, dy, dz);

  /* Adjustments for odd data dimensions */
  if (nx>1 && (nx % 2)) nx_mod= nx-1;
  else nx_mod= nx;
  if (ny>1 && (ny % 2)) ny_mod= ny-1;
  else ny_mod= ny;
  if (nz>1 && (nz % 2)) nz_mod= nz-1;
  else nz_mod= nz;

  /* Step counter */
  count_set_phase++;

  /* No need for FFT here; we are assuming k-space data. */

  /* dx, dy, dz are in voxels.  Note that we need to reverse
   * phases if we are doing the unexpected case (the future
   * FFT is in the positive direction).
   */
  x_phase_scale= 2.0*M_PI*dx/((double)nx_mod);
  y_phase_scale= 2.0*M_PI*dy/((double)ny_mod);
  z_phase_scale= 2.0*M_PI*dz/((double)nz_mod);
  if (!kspace_flag) {
    x_phase_scale= -x_phase_scale;
    y_phase_scale= -y_phase_scale;
    z_phase_scale= -z_phase_scale;
  }

  /*** Apply phase changes ***/
  /* Signs of terms are determined by relationship between data
   * coordinate system and geometric coordinate system.
   */

  if (real_flag && !(nx%2)) lowerx= 1;
  else lowerx= 0;

  if (real_flag && !(ny%2)) lowery= 1;
  else lowery= 0;

  if (real_flag && !(nz%2)) lowerz= 1;
  else lowerz= 0;

  for (i=lowerx; i<nx; i++) {
    theta= (i-halfx)*x_phase_scale;
    c= cos(theta);
    s= sin(theta);
    for (j=0; j<ny; j++)
      for (k=0; k<nz; k++) {
	here= &(MEM(image,nx,ny,nz,i,j,k));
	tc= *here;
	here->real= c*tc.real - s*tc.imag;
	here->imag= c*tc.imag + s*tc.real;
      }
  }

  for (i=0; i<nx; i++)
    for (j=lowery; j<ny; j++) {
      theta= -(j-halfy)*y_phase_scale;
      c= cos(theta);
      s= sin(theta);
      for (k=0; k<nz; k++) {
	here= &(MEM(image,nx,ny,nz,i,j,k));
	tc= *here;
	here->real= c*tc.real - s*tc.imag;
	here->imag= c*tc.imag + s*tc.real;
      }
    }

  for (i=0; i<nx; i++)
    for (k=lowerz; k<nz; k++) {
      theta= (k-halfz)*z_phase_scale;
      c= cos(theta);
      s= sin(theta);
      for (j=0; j<ny; j++) {
	here= &(MEM(image,nx,ny,nz,i,j,k));
	tc= *here;
	here->real= c*tc.real - s*tc.imag;
	here->imag= c*tc.imag + s*tc.real;
      }
    }

  /* No need for FFT here; we are assuming k-space data. */

}

void fourier_shift_rot3d( Quat* q, double dx, double dy, double dz,
			  FComplex* orig_image,
			  FComplex* moved_image,
			  long nx, long ny, long nz,
			  double length_x, double length_y, double length_z,
			  int real_flag )
{
  TransParams trans;
  int i;

  trans.dx= dx;
  trans.dy= dy;
  trans.dz= dz;

  /* Step counter */
  count_calls++;

  /* Copy into output buffer */
  for (i=0; i<nx*ny*nz; i++) {
    moved_image[i].real= orig_image[i].real;
    moved_image[i].imag= orig_image[i].imag;
  }

  if (q->x==0.0 && q->y==0.0 && q->z==0.0) {
    shift_only( &trans, moved_image, nx, ny, nz,
		length_x, length_y, length_z, real_flag );
  }
  else if (q->x==0.0 && q->y==0.0 && fabs(q->w)>=CANCELLATION_TOL) {
    rot_z( q, &trans, moved_image, nx, ny, nz, 
	   length_x, length_y, length_z, real_flag );
  }
  else if (q->y==0.0 && q->z==0.0 && fabs(q->w)>=CANCELLATION_TOL) {
    rot_x( q, &trans, moved_image, nx, ny, nz, 
	   length_x, length_y, length_z, real_flag );
  }
  else if (q->z==0.0 && q->x==0.0 && fabs(q->w)>=CANCELLATION_TOL) {
    rot_y( q, &trans, moved_image, nx, ny, nz, 
	   length_x, length_y, length_z, real_flag );
  }
  else if (shear_pattern==FR3D_SHEAR_13) {
    /* Slow, fully symmetric decomp of rotation requested.  
     * (Axis-aligned rotations are inherently fully symmetrical.
     */
    rot_13_shears(q, dx, dy, dz, moved_image, nx, ny, nz, 
		 length_x, length_y, length_z, real_flag);
  }
  else if (shear_pattern==FR3D_SHEAR_7) {
    /* 7-shear rotation requested.  This decomp is symmetric for
     * small rotation angles.
     */
    rot_7_shears(q, dx, dy, dz, moved_image, nx, ny, nz, 
		 length_x, length_y, length_z, real_flag);
  }
  else if (shear_pattern==FR3D_SHEAR_4) {
    rot_4_shears(q, dx, dy, dz, moved_image, nx, ny, nz, 
		 length_x, length_y, length_z, real_flag);
  }
  else Abort("fshrot3d: unknown shear pattern %d requested!\n",
	     shear_pattern);
	     
}


