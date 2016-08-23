/************************************************************
 *                                                          *
 *  linrot3d.c                                              *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *                                                          *
 *  Copyright (c) 2000 Department of Statistics             *
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
 *  Original programming by Joel Welling 1/2000             *
 ************************************************************/
/* Implementation of 3D rotation and shift by linear interpolation */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"

static char rcsid[] = "$Id: linrot3d.c,v 1.10 2003/09/25 19:45:42 welling Exp $";

/* Notes-
 *
 */

/* Accessor for 3D arrays.  If this is changed, the striding math in
 * the main interpolation routine must also be changed. 
 */
#define MEM(matrix,nx,ny,nz,x,y,z) matrix[((((x)*ny)+(y))*nz)+(z)]

/* Angle conversion */
#define DEG2RAD (M_PI/180.0)
#define RAD2DEG (180.0/M_PI)

/* Tolerance for coord calculations (allow for numerical noise) */
#define EPSILON 0.000001

/* Linear interpolation, up to 3 steps.  Scale factors are:
 *
 * a is 0.0 on l side of l-r direction
 * b is 0.0 on d side of d-u direction
 * c is 0.0 on f side of f-b direction
 */


#define LIN( l, r, a )  ( ( (a)<=0.0 ) ? (l) : ( ((a)>=1.0) ? (r) : \
                                              ( (l)*(1.0-(a)) + (r)*(a) ) ) )
#define LIN2( ld, rd, lu, ru, a, b ) LIN( LIN(ld,rd,a), LIN(lu,ru,a), b )
#define LIN3( ldf, rdf, luf, ruf, ldb, rdb, lub, rub, a, b, c ) \
  LIN( LIN2( ldf, rdf, luf, ruf, a, b ), LIN2( ldb, rdb, lub, rub, a, b ), c )

static int debug= 0;
static int count_calls= 0;

void linrot3d_set( int flag, int value )
{
  switch (flag) {
  case LR3D_DEBUG: {
    debug= value; 
  }
  break;
  }
}

int linrot3d_get( int flag )
{
  switch (flag) {
  case LR3D_DEBUG: {
    return debug; 
  }
    /* break; -NOTREACHED- */
  default:
    Abort("linrot3d_get: unknown flag %d!\n",flag);
  }
  return 0; /* not reached */
}

void linrot3d_set_debug( int flag )
{
  linrot3d_set( LR3D_DEBUG, flag );
}

void linrot3d_clear_counts()
{
  count_calls= 0;
}

void linrot3d_get_counts( int* ncalls )
{
  *ncalls= count_calls;
}

static void v3_add( double val[3], double addme[3] )
{
  val[0] += addme[0];
  val[1] += addme[1];
  val[2] += addme[2];
}

static void v3_inv_transform( double out[3], double in[3], 
			      Quat* q, double dx, double dy, double dz )
{
  /* This routine performs the *inverse* of the transformation given
   * by the quaternion and displacements input.  Specifically, the
   * forward transform would be In = Q Out Q' + d (where ' is conjugate),
   * so this transform is: Out =  Q'(In - d)Q
   */
  Quat tqconj;
  Quat v;

  /* This code could be wildly more efficient if we unrolled these
   * quaternion operations.
   */
  quat_copy(&tqconj, q);
  quat_conjugate(&tqconj);
  v.x= in[0]-dx;
  v.y= in[1]-dy;
  v.z= in[2]-dz;
  v.w= 0.0;
  quat_mult_right(quat_mult_left(&tqconj,&v), q);
  out[0]= v.x;
  out[1]= v.y;
  out[2]= v.z;
}

/*
 * Quat* q and the dx, dy, dz values specify rotation and shift (in voxels)
 * orig_image is input
 * moved_image is output
 * check is output; non-zero where moved_image data is valid
 * nx, ny, nz are grid dimensions
 * length_x, length_y, length_z are voxel edge lengths
 * kspace_flag should be set if both real and imaginary parts should
 *   be transformed;  if it's zero only the real half of the input data
 *   will be used.
 *
 * Input and output data are assumed to be ordered such that z is 
 * fastest in memory.  The convention is that the rotation is applied 
 * *before* the shift; that is, R = T*r where R is the whole transformation,
 * T is the translation (shift), and r is the pure rotation given by
 * the input quaternion.
 */
void linear_shift_rot3d( Quat* q, double dx, double dy, double dz,
			 FComplex* orig_image,
			 FComplex* moved_image,
			 char* check,
			 long nx, long ny, long nz,
			 double length_x, double length_y, double length_z,
			 int kspace_flag)
{
  double pout[3]; /* point in output space (grid-aligned) */
  double p[3]; /* point in input space */
  double zoutstep[3]; /* one step in Z in output space */
  double zstep[3]; /* one step in Z in input space */
  int iout, jout, kout;
  int i, j, k;
  long halfx= nx/2;
  long halfy= ny/2;
  long halfz= nz/2;

  /* Step counter */
  count_calls++;

  /* Calculate rescaled Z-axis step direction.  Signs of terms are
   * determined by relationship between grid and 3D (radiological) coords.
   */
  zoutstep[0]= 0.0;
  zoutstep[1]= 0.0;
  zoutstep[2]= length_z/nz;
  v3_inv_transform( zstep, zoutstep, q, 0.0, 0.0, 0.0 );
  zstep[0] *= nx/length_x; /* rescale */
  zstep[1] *= -ny/length_y;
  zstep[2] *= nz/length_z;
  if (debug) fprintf(stderr,"stride along Z: %f %f %f <- %f %f %f\n",
		     zoutstep[0],zoutstep[1],zoutstep[2],
		     zstep[0],zstep[1],zstep[2]);

  for (iout=0; iout<nx; iout++) 
    for (jout=0; jout<ny; jout++) {

      /* Calculate input-space coords for beginning of this Z row.
       * Signs of terms are determined by relationship between grid
       * and 3D (radiological) coords.
       */
      pout[0]= (iout-halfx)*(length_x/nx);
      pout[1]= -(jout-halfy)*(length_y/ny);
      pout[2]= (0.0-halfz)*(length_z/nz);
      v3_inv_transform( p, pout, q, dx, dy, dz );
      p[0] *= nx/length_x; /* rescale */
      p[1] *= -ny/length_y;
      p[2] *= nz/length_z;
      p[0] += halfx;    /* shift to grid coords */
      p[1] += halfy;
      p[2] += halfz;

#ifdef never
      if (iout==31 && jout==31) {
	fprintf(stderr,"row start: %f %f %f <- %f %f %f\n",
		pout[0],pout[1],pout[2],p[0],p[1],p[2]);
	fprintf(stderr,"zstep is %f %f %f\n",zstep[0],zstep[1],zstep[2]);
      }
#endif

      /* Walk the row in the Z direction */
      for (kout=0; kout<nz; kout++) {
	if ((p[0]<(-EPSILON)) || (p[0]>(double)(nx-1))
	    || (p[1]<(-EPSILON)) || (p[1]>(double)(ny-1))
	    || (p[2]<(-EPSILON)) || (p[2]>(double)(nz-1))) {

	  /* This point maps outside the input array */
	  MEM(moved_image,nx,ny,nz,iout,jout,kout).real= 0.0;
	  if (kspace_flag)
	    MEM(moved_image,nx,ny,nz,iout,jout,kout).imag= 0.0;
	  MEM(check,nx,ny,nz,iout,jout,kout)= 0; /* mark out */

	}
	else {
	  /* Valid location */
	  float ax, ay, az;
	  FComplex* here;

	  i= (int)floor( p[0] );
	  if (i<0) i= 0;
	  ax= p[0]-i;
	  j= (int)floor( p[1] );
	  if (j<0) j= 0;
	  ay= p[1]-j;
	  k= (int)floor( p[2] );
	  if (k<0) k= 0;
	  az= p[2]-k;

	  here= (&(MEM(orig_image,nx,ny,nz,i,j,k)));

	  MEM(moved_image,nx,ny,nz,iout,jout,kout).real= 
	    LIN3( here->real, 
		  (here + ny*nz)->real,
		  (here + nz)->real,
		  (here + ny*nz + nz)->real,
		  (here + 1)->real, 
		  (here + ny*nz + 1)->real,
		  (here + nz + 1)->real,
		  (here + ny*nz + nz + 1)->real,
		  ax, ay, az );

	  if (kspace_flag) {
	    MEM(moved_image,nx,ny,nz,iout,jout,kout).imag= 
	      LIN3( here->imag, 
		    (here + ny*nz)->imag,
		    (here + nz)->imag,
		    (here + ny*nz + nz)->imag,
		    (here + 1)->imag, 
		    (here + ny*nz + 1)->imag,
		    (here + nz + 1)->imag,
		    (here + ny*nz + nz + 1)->imag,
		    ax, ay, az );
	  }

	  MEM(check,nx,ny,nz,iout,jout,kout)= 1; /* mark in */

	}	

	v3_add( p, zstep );
      }

    }

}


