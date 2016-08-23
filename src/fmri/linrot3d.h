/************************************************************
 *                                                          *
 *  linrot3d.h                                              *
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
/* Header file for linrot3d.c */

/* Switches to be passed to setflag and getflag */
#define LR3D_DEBUG 0

/* Entry points for getting and setting linrot3d switches */
void linrot3d_set(int, int);
int linrot3d_get(int);

/* Turn debugging on or off (backward compatibility) */
void linrot3d_set_debug( int flag );

/*
  Quat* q and the dx, dy, dz values specify rotation and shift (in voxels)
  orig_image is input
  moved_image is output
  check is output; non-zero where moved_image data is valid
  nx, ny, nz are grid dimensions
  length_x, length_y, length_z are voxel edge lengths
  kspace_flag should be set if both real and imaginary parts should
    be transformed;  if it's zero only the real half of the input data
    will be used.
  Convention is that rotation is applied *before* shift
 */
void linear_shift_rot3d( Quat* q, double dx, double dy, double dz,
			 FComplex* orig_image,
			 FComplex* moved_image,
			 char* check,
			 long nx, long ny, long nz,
			 double length_x, double length_y, double length_z,
			 int kspace_flag );


/* Clear and get the counters for operations (for diagnostics) */
void linrot3d_clear_counts(void);
void linrot3d_get_counts( int* ncalls );

