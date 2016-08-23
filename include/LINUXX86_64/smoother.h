/************************************************************
 *                                                          *
 *  smoother.h                                              *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1998 Department of Statistics,         *
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
 *  Original programming by Joel Welling, 8/98              *
 ************************************************************/
/* This package implements smoothing methods. */

typedef enum { 
  SM_GAUSSIAN, SM_TRIANGULAR, SM_POWER, SM_NONE, 
  SM_DDX, SM_SHIFT, SM_LINTERP, SM_RUNNINGSUM,
  SM_MEDIAN
} sm_type;

struct smoother_struct;

typedef float (*SmootherThreshTest)( struct smoother_struct* sm,
				     float** dtbl_in, int ndata, 
				     int n1, int n2 );
  
typedef struct smoother_struct {
  void (*smooth)(struct smoother_struct*, float*, float*, int,
		 unsigned char**, int);
  void (*smooth_group)(struct smoother_struct*, float**, float**, int, int,
		       unsigned char**, int);
  sm_type type;
  double bandwidth;
  double k;
  double threshold;
  SmootherThreshTest thresh_test;
  int n;
  void* data;
  int last;
  int smDim;
} Smoother;

void sm_init(void); /* Initialize */

void sm_parse_cl_opts(void);

void sm_set_params( sm_type type_in, float bandwidth_in, float k_in, 
		    float threshold_in, SmootherThreshTest test_in );

void sm_get_params( sm_type* type_out, float* bandwidth_out, float* k_out,
		    float* threshold_out, SmootherThreshTest* test_out );

void sm_set_direction( Smoother* smoother, int dir );
int sm_get_direction( Smoother* smoother );

Smoother* sm_create_smoother(void); /* create based on currently set vals */

void sm_destroy( Smoother* smoother );

/************
 * Rules for the use of missing data by the smoothers:
 * 
 * If 'missing' data is not supplied, smoothing is done using all inputs.
 * If 'missing' data is supplied but sm_get_direction() is not t or z,
 *   missing information is ignored and smoothing is done using all inputs.
 * If 'missing' is given and sm_get_direction() is t, smoothing is
 *   done excluding those points for which missing[i][z] is non-zero.
 * If 'missing' is given and sm_get_direction() is z, smoothing is
 *   done excluding those points for which missing[z][i] is non-zero.
 *
 * Needless to say, the reasons for these conventions are lost in history.
 ***********/

#define SM_SMOOTH( smoother, data_in, data_out, n, missing, z ) \
  (*(smoother->smooth))(smoother, data_in, data_out, n, missing, z);

#define SM_SMOOTH_GROUP(smoother,datatbl_in,datatbl_out,ndata,n,miss,z) \
  (*(smoother->smooth_group))(smoother,datatbl_in,datatbl_out,ndata,n,miss,z)


