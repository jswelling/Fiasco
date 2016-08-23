/************************************************************
 *                                                          *
 *  fmri.h                                                   *
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
 ************************************************************/
/* This is the overall include file for things in the fmri directory */

#define MALLOC_FAILURE(num,type)\
  { fprintf(stderr,"%s:%d: %s: unable to allocate %d %ss\n",\
            __FILE__,__LINE__,__func__,num,#type); exit(-1); }

typedef struct FComplex
{
  float real;
  float imag;
} FComplex;

typedef struct SComplex
{
  short real;
  short imag;
} SComplex;

typedef struct DComplex
{
  double real;
  double imag;
} DComplex;

typedef struct RegPars
{
  float x_shift;           /* Unit is in pixels */
  float y_shift;           /* Unit is in pixels */
  float rotation;          /* Unit is in degrees */
  int valid;               /* if we need to know if other values are set */
} RegPars;

/* to allow inclusion of this w/o mri.h */
#ifndef MRI_MAX_DIMS
typedef struct dummy_mri_dataset_struct MRI_Dataset; 
#endif

#ifndef FreeMatrix
#define FreeMatrix(p) { free(*(p)); free((p)); }
#endif

/* Utility Functions (in fmri.c) */
extern unsigned char **get_missing(MRI_Dataset* d);
extern long get_typesize(MRI_Dataset* d, char* chunkname);
void realign_matrix( void** ptr, long n, size_t rowsize );
void efseek( FILE* fp, long offset, int whence, char* filename );
float Modulus( FComplex z );
float Phase( FComplex z );
void fqsrt( float* v, long left, long right );

/* Optimization functions (in polyt.c, praxis.c, bvls.c) */
extern int nelmin_t(float (*fn)(float*,void*), void (*rst)(float*,void*), 
		    int* n, float* start, float* mmin, float* ynewl, 
		    float* reqmin, float* step, int* konvge, int* kcount, 
		    int* icount, int* numres, int* ifault, int* maxres, 
		    float* tcritsqr, void* userHook);
extern int nelmin(float (*fn)(float*,void*), void (*rst)(float*,void*), 
		  int* n, float* start, float* mmin, float* ynewl, 
		  float* reqmin, float* step, int* konvge, int* kcount, 
		  int* icount, int* numres, int* ifault, int* maxres, 
		  void* userHook);

extern double praxis(double t0, double machep, double h0, int n, int prin, 
		     double* x, double (*f)(double*,int,void*), 
		     void (*reset)(double*,int,void*), double fmin, 
		     void* userHook);

extern int bvls(long* key, long* m, long* n, 
		double* a, double* b, double* bl, double* bu, 
		double* x, double* w, double* act, double* zz, 
		long* istate, long* loopa);

extern double fmin1D(double ax, double bx, double (*func)(double), 
		     double tol, double eps, int debug);


/* Other Functions (in fft2d.c, fft3d.c, and fshrot.c) */
void fft2d( FComplex** data, long nx, long ny, long sign, 
	    char rowcol, long from, long to );
void fft3d( FComplex* data, long nx, long ny, long nz, long sign,
	    char* rowcol );
void fourier_shift_rot( RegPars par, FComplex** orig_image, 
			FComplex** moved_image, long nx, long ny, 
			char domain );

/* FIAT version info */
#include "fiat.h"

/* Header file for glm.c */
#include "glm.h"

/* Header file for smoother.c */
#include "smoother.h"

/* Header for parsesplit.c */
#include "parsesplit.h"

/* Header for quaternion.c */
#include "quaternion.h"

/* Header file for optimizer.c */
#include "optimizer.h"

/* Header for fshrot3d.c */
#include "fshrot3d.h"

/* Header for linrot3d.c */
#include "linrot3d.h"

/* Header for linwarp.c */
#include "linwarp.h"

/* Header for closest_warp.c */
#include "closest_warp.h"

/* Header for history stuff */
#include "history.h"

/* Error and message handling routines */
#include "errors.h"

/* Header for file type identification routines */
#include "filetypes.h"

/* Header for generic key-value hash table */
#include "kvhash.h"

/* Header for RPN math engine */
#include "rpn_engine.h"

/* Header for image entropy and mutual information routines */
#include "entropy.h"

/* Header for spline.h */
#include "spline.h"

/* Header for interpolator.c */
#include "interpolator.h"

/* Header for slicepattern.c */
#include "slicepattern.h"

/* Header for MRI supplemental utilities */
#include "mriu.h"

/* Header for Kalman filter utilities */
#include "kalmanfilter.h"
