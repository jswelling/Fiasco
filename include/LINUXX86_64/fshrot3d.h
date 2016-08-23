/************************************************************
 *                                                          *
 *  fshrot3d.h                                              *
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
/* Header file for fshrot3d.c */

/* Switches to be passed to setflag and getflag */
#define FR3D_DEBUG 0
#define FR3D_SHEAR_PATTERN 1
#define FR3D_QUAL_MEASURE 2

#define FR3D_SHEAR_4 0
#define FR3D_SHEAR_7 1
#define FR3D_SHEAR_13 2

#define FR3D_QUAL_UNSET 0
#define FR3D_QUAL_COX 1
#define FR3D_QUAL_SUM_ABS 2
#define FR3D_QUAL_SUM_SQR 3
#define FR3D_QUAL_UNIT_CELL 4

/* Entry points for getting and setting fshrot3d switches */
void fshrot3d_set(int, int);
int fshrot3d_get(int);

/* Turn debugging on or off (backward compatibility) */
void fshrot3d_set_debug( int flag );

/*
  Quat* q and the dx, dy, dz values specify rotation and shift (in voxels)
  orig_image is input
  moved_image is output
  nx, ny, nz are grid dimensions
  real_flag is zero if the data is complex (both real and
    imaginary parts of images are occupied), one if only real
    parts are occupied.
  length_x, length_y, length_z are voxel edge lengths
  Convention is that rotation is applied *before* shift
 */
void fourier_shift_rot3d( Quat* q, double dx, double dy, double dz,
			  FComplex* orig_image,
			  FComplex* moved_image,
			  long nx, long ny, long nz,
			  double length_x, double length_y, double length_z,
			  int kspace_flag);

/* This routine simply sets phases in an image in such a way that
 * a subsequent full 3D FFT causes a translation given by dx, dy,
 * and dz (which are in voxels).  The image is phase shifted in
 * place.  kspace_flag is nonzero if this is a k-space image, and note
 * that that means in all 3 directions, so that the expected coming 
 * 3D FFT is in the inverse direction.
 * real_flag has the following meaning:
 *   if (real_flag && !kspace_flag) the data is pure real (imaginary
 *       parts 0) and is in image space
 *   if (real_flag && kspace_flag) the data represents the FFT of 
 *       real data, so when it is inverse-fft'd the result should be
 *       pure real.
 *   if (!real_flag && !kspace_flag) the data is complex numbers in
 *       image space.
 *   if (!real_flag && kspace_flag) the data is complex numbers in kspace.
 */
void fshrot3d_set_shift_phases( double dx, double dy, double dz,
				FComplex* image,
				long nx, long ny, long nz, double length_x, 
				double length_y, double length_z,
				int kspace_flag, int real_flag );

/* Clear and get the counters for the shear operations (for diagnostics) */
void fshrot3d_clear_shear_counts(void);
void fshrot3d_get_shear_counts( int* xcount, int* ycount, int* zcount,
				int* phase_count, 
				double* mean_shears_per_call );

