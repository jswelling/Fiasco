/************************************************************
 *                                                          *
 *  closest_warp.c                                          *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *                                                          *
 *  Copyright (c) 2004 Department of Statistics             *
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
 *  Derived from linwarp by Joel Welling 11/2004            *
 ************************************************************/
/* Implementation of 3D rotation and shift by 'closest' interpolation */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"

static char rcsid[] = "$Id: closest_warp.c,v 1.3 2007/03/21 23:50:20 welling Exp $";

/* Notes-
 *
 */

/* Accessor for 3D arrays.  If this is changed, the striding math in
 * the main interpolation routine must also be changed. 
 */
#define MEM(matrix,nx,ny,nz,x,y,z) matrix[((((x)*ny)+(y))*nz)+(z)]

/* Tolerance for coord calculations (allow for numerical noise) */
#define EPSILON 0.000001

static int debug= 0;
static int count_calls= 0;

void clwarp_set( int flag, int value )
{
  switch (flag) {
  case LWARP_DEBUG: {
    debug= value; 
  }
  break;
  }
}

int clwarp_get( int flag )
{
  switch (flag) {
  case LWARP_DEBUG: {
    return debug; 
  }
    /* break; -NOTREACHED- */
  default:
    Abort("clwarp_get: unknown flag %d!\n",flag);
  }
  return 0; /* not reached */
}

void clwarp_set_debug( int flag )
{
  clwarp_set( LWARP_DEBUG, flag );
}

void clwarp_clear_counts()
{
  count_calls= 0;
}

void clwarp_get_counts( int* ncalls )
{
  *ncalls= count_calls;
}

/*
 * Transform t_in values specify the warp (in mm)
 * orig_image is input
 * moved_image is output
 * check is output; non-zero where moved_image data is valid
 * nx, ny, nz are grid dimensions
 * length_x, length_y, length_z are voxel edge lengths
 *
 * Input and output data are assumed to be ordered such that z is 
 * fastest in memory.
 */
void clwarp_warp( Transform t_in, 
		  FComplex* orig_image,
		  FComplex* moved_image,
		  char* check,
		  long nx, long ny, long nz,
		  double length_x, double length_y, double length_z )
{
  Vec4 pout; /* point in output space (grid-aligned) */
  Vec4 p; /* point in input space */
  long i, j, k;
  long iout, jout, kout;
  long halfx= nx/2;
  long halfy= ny/2;
  long halfz= nz/2;
  Transform t;
  Transform tTemp;

#ifdef never
  double xmin= 1000.0;
  double xmax= -1000.0;
  double ymin= 1000.0;
  double ymax= -1000.0;
  double zmin= 1000.0;
  double zmax= -1000.0;
#endif

  /* Step counter */
  count_calls++;

  /* Make a version of the transform in which the voxel scale
   * factors have been sucked into the matrix.  Signs of terms 
   * are determined by relationship between grid and 3D (radiological) 
   * coords.
   */
  trans_identity(tTemp);
  tTemp[0]= length_x/(double)nx;
  tTemp[5]= -length_y/(double)ny;
  tTemp[10]= length_z/(double)nz;
  trans_copy(t, t_in);
  trans_mult_right(t,tTemp);
  tTemp[0]= 1.0/tTemp[0];
  tTemp[5]= 1.0/tTemp[5];
  tTemp[10]= 1.0/tTemp[10];
  trans_mult_left(tTemp,t);

  /* Just do it */
  pout[3]= 1.0;
  for (iout=0; iout<nx; iout++) {
    pout[0]= (double)(iout-halfx);
    for (jout=0; jout<ny; jout++) {
      pout[1]= (double)(jout-halfy);
      for (kout=0; kout<nz; kout++) {
	pout[2]= (double)(kout-halfz);
	bcopy(pout,p,sizeof(pout));
	trans_vec_mult(t, p);
	p[0] += halfx;
	p[1] += halfy;
	p[2] += halfz;
#ifdef never
	if (p[0]<xmin) xmin= p[0];
	if (p[0]>xmax) xmax= p[0];
	if (p[1]<ymin) ymin= p[1];
	if (p[1]>ymax) ymax= p[1];
	if (p[2]<zmin) zmin= p[2];
	if (p[2]>zmax) zmax= p[2];
#endif
	if ((p[0]<(-EPSILON)) || (p[0]>(double)(nx-1)+EPSILON)
	    || (p[1]<(-EPSILON)) || (p[1]>(double)(ny-1)+EPSILON)
	    || (p[2]<(-EPSILON)) || (p[2]>(double)(nz-1))+EPSILON) {
	  
	  /* This point maps outside the input array */
	  MEM(moved_image,nx,ny,nz,iout,jout,kout).real= 0.0;
	  MEM(moved_image,nx,ny,nz,iout,jout,kout).imag= 0.0;
	  MEM(check,nx,ny,nz,iout,jout,kout)= 0; /* mark out */
	}
	else {
	  /* Valid location */
	  i= (int)rint( p[0] );
	  j= (int)rint( p[1] );
	  k= (int)rint( p[2] );
	  
	  MEM(moved_image,nx,ny,nz,iout,jout,kout)=
	    MEM(orig_image,nx,ny,nz,i,j,k);
	  
	  MEM(check,nx,ny,nz,iout,jout,kout)= 1; /* mark in */
	} 
      }
    }
  }
#ifdef never
  fprintf(stderr,"mapping %g<x<%g, %g<y<%g, %g<z<%g\n",
	  xmin,xmax,ymin,ymax,zmin,zmax);
#endif
}


