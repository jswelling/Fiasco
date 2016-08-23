/************************************************************
 *                                                          *
 *  interpolator.h                                              *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *                                                          *
 *  Copyright (c) 2005 Department of Statistics             *
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
/* Header file for interpolator.c */

#ifndef INCL_INTERPOLATOR_H
#define INCL_INTERPOLATOR_H 1

#include "spline.h"

/* Switches to be passed to set and get options */
#define INTRP_OPT_DEBUG 0
#define INTRP_OPT_FASTBLK 1
#define INTRP_OPT_EXTENT 2
#define INTRP_OPT_TENSION 3

typedef enum { INTRP_CLOSEST, INTRP_LINEAR, INTRP_CATMULLROM, 
	       INTRP_BEZIER, INTRP_BSPLINE, INTRP_UNKNOWN } InterpolatorType;

extern InterpolatorType intrp_typeFromName( const char* name );
extern const char* intrp_nameFromType( InterpolatorType type );
extern void intrp_warpClearCounts(void);
extern void intrp_warpGetCounts(long* count);

typedef struct Interpolator_struct {
  void (*prep)( struct Interpolator_struct* self, double* data, long sz );
  void (*calc)( struct Interpolator_struct* self, 
		double* result, double* loc, long runLength, long offset );
  void (*destroySelf)( struct Interpolator_struct* self );
  void (*dumpSelf)( const struct Interpolator_struct* self, FILE* ofile );
  void (*setInt)( struct Interpolator_struct* self, int which, long val );
  void (*setDouble)( struct Interpolator_struct* self, int which, double val );
  long (*getInt)( const struct Interpolator_struct* self, int which );
  double (*getDouble)( const struct Interpolator_struct* self, int which );
  const char* typeName;
  double* dataField;
  long dataFieldLength;
  void* hook;
  int debug;
  long fast_blksize; /* usually corresponds to fast_blksize */
  long extent; /* dim of space in which to interp */
} Interpolator;

/*
 * Interpolator* interp specifies the interpolator to use (better be 3D!)
 * Transform t values specify the warp (in mm)
 * orig_image is input
 * moved_image is output
 * check is output; non-zero where moved_image data is valid
 * nx, ny, nz are grid dimensions; fast_blk is the blocksize to the
 * left of the X dimension (typically dv).
 * length_x, length_y, length_z are voxel edge lengths (mm)
 *
 * Input and output data are assumed to be ordered such that x is 
 * fastest in memory.  Note that this differs from closest_warp
 * and linwarp!
 */
void intrp_warpApply( Interpolator* interp, Transform t,
		      double* orig_image,
		      double* moved_image,
		      char* check,
		      long nx, long ny, long nz, long fast_blk,
		      double length_x, double length_y, double length_z );

/* fast_blk in the following factory methods is the dimensionality of the
 * field to be interpolated, not the dimensionality of the space that
 * the interpolation happens in.  For example, when interpolating a
 * rank-3 vector field along a line, use and Interpolator1D with fast_blk=3.
 * For 1D interpolation one can think of nx as selected_extent; for
 * 2D and 3D the ny and nz dimensions are the extents of the subsequent
 * (and slower-varying) indices.
 * The 'prep' method then provides selected_extent*fast_blksize doubles;
 * the sz parameter is used to check for consistency with the constructor.
 * The 'calc' method calculates 'runLength' interpolants, offset from
 * the beginning of fast_blksize by 'offset'.  Thus runLength+offset 
 * must be <= fast_blksize==fast_blk.  The location sampled is
 * loc of the way through the selected extent, so the full valid range
 * of loc is 0.0 <= loc <= (extent-1.0)
 */
extern Interpolator* 
intrp_createClosestInterpolator1D( long nx, long fast_blk );

extern Interpolator* 
intrp_createClosestInterpolator2D( long nx, long ny, long fast_blk );

extern Interpolator* 
intrp_createClosestInterpolator3D( long nx, long ny, long nz, long fast_blk );

extern Interpolator* 
intrp_createLinearInterpolator1D( long nx, long fast_blk );

extern Interpolator* 
intrp_createLinearInterpolator2D( long nx, long ny, long fast_blk );

extern Interpolator* 
intrp_createLinearInterpolator3D( long nx, long ny, long nz, long fast_blk );

extern Interpolator* 
intrp_createSplineInterpolator1D( SplineType type, long nx, long fast_blk );

extern Interpolator* 
intrp_createSplineInterpolator2D( SplineType type, long nx, long ny, 
				  long fast_blk );

extern Interpolator* 
intrp_createSplineInterpolator3D( SplineType type, long nX, long ny, long nz, 
				  long fast_blk );

extern Interpolator* 
intrp_createInterpolator1DByType( InterpolatorType type, 
				  long nx, long fast_blk );
extern Interpolator* 
intrp_createInterpolator2DByType( InterpolatorType type, 
				   long nx, long ny, long fast_blk );
extern Interpolator* 
intrp_createInterpolator3DByType( InterpolatorType type, 
				  long nx, long ny, long nz, long fast_blk );
						       
#endif
