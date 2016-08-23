/************************************************************
 *                                                          *
 *  algorithm.h                                             *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1995 Department of Statistics,         *
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
 *  Original programming by Joel Welling 3/00               *
 *      10/02: Parallelization, Jenn Bakal                  *
 ************************************************************/


/* $Id: algorithm.h,v 1.6 2005/07/30 02:11:31 welling Exp $ */


/* Maximum free parameters */
#define MAX_DOF 20

/* Number of mutual info contexts */
#define N_MUTUAL_INFO_CONTEXTS 1

/* Access for 3D arrays */
#define MEM_BCK(matrix,nx,ny,nz,x,y,z) matrix[((((z)*ny)+(y))*nx)+(x)]
#define MEM(matrix,nx,ny,nz,x,y,z) matrix[((((x)*ny)+(y))*nz)+(z)]

typedef enum {
  OBJECTIVE_MSE,
  OBJECTIVE_MI,
  OBJECTIVE_JE,
  OBJECTIVE_INVALID
} ObjectiveMethod;

typedef enum {
  SEARCH_TRILIN,
  SEARCH_SHEAR4_FFT,
  SEARCH_SHEAR7_FFT,
  SEARCH_SHEAR13_FFT,
  SEARCH_INVALID
} SearchMethod;

typedef enum {
  WEIGHT_CONST,
  WEIGHT_INV_STDV,
  WEIGHT_ALIGN,
  WEIGHT_SMTHALIGN,
  WEIGHT_INVALID
} WeightMethod;

typedef enum {
  OPT_NELMIN,
  OPT_NELMIN_T,
  OPT_PRAXIS,
  OPT_NONE,
  OPT_INVALID
} OptMethod;

typedef struct algorithm_struct {
  char* progname;
  int debugLevel;
  SearchMethod inner_search_method;  /* this defines the algorithm to be */
  SearchMethod outer_search_method;  /* used.  Technically, the quality  */
  int smooth_flag;                   /* measure used by fshrot3d is also */
  int inplane_flag;                  /* part of the algorithm in effect  */
  int trans_only_flag;
  int rot_only_flag;
  int x_only_flag;
  WeightMethod weight_method;
  ObjectiveMethod objective_method;
  double weight_floor;
  Optimizer* opt;
  Smoother* smoother;
  MutualInfoContext* mutualInfoContext[N_MUTUAL_INFO_CONTEXTS];
  int (*getNDimCB)(struct algorithm_struct* alg);
} Algorithm;

int algInit( Algorithm* alg, const char* string, const char* progname_in,
	     int (*getNDimCallback)(Algorithm* alg) );
void algFinalize( Algorithm* alg );
void algParseOpts( Algorithm* alg );
const char* algGetInfoString( const Algorithm* alg );
int algParseInfoString( Algorithm* alg, const char* string );
void algPack( const Algorithm* alg );
void algUnpack( Algorithm* alg, const char* progname_in,
		int (*getNDimCallback)(Algorithm* alg));
void algSetDebugLevel( Algorithm* alg, const int lvl );
int algGetDebugLevel( Algorithm* alg );
int algNeedsMask( Algorithm* alg );
int algNeedsStdv( Algorithm* alg );
int algGetCurrentNDim( Algorithm* alg );
void algSmoothImage( Algorithm* alg, float* img, int dx, int dy, int dz,
		     unsigned char** missing, int t );
void algMaybeSmoothImage( Algorithm* alg, float* img, int dx, int dy, int dz,
			  unsigned char** missing, int t );
void algSmoothImageComplex( Algorithm* alg, FComplex* img, 
			    int dx, int dy, int dz, 
			    unsigned char** missing, int t );
void algMaybeSmoothImageComplex( Algorithm* alg, FComplex* img, 
				 int dx, int dy, int dz, 
				 unsigned char** missing, int t );
void algPrep(); /* called once before using the package */
double algCalcChiSqr(Algorithm* alg, FComplex* moved_image, 
		     float* align_image, float* weight_image, int* mask,
		     char* check, int dx, int dy, int dz);
double algCalcJointEntropy(Algorithm* alg, FComplex* moved_image, 
			   float* align_image, float* weight_image, int* mask,
			   char* check, int dx, int dy, int dz);
double algCalcMutualInfo(Algorithm* alg, FComplex* moved_image, 
			 float* align_image, float* weight_image, int* mask,
			 char* check, int dx, int dy, int dz);
void algMaybeBuildMask(Algorithm* alg, float* weight_image, int* mask, 
		       int dx, int dy, int dz);
void algBuildWeight(Algorithm* alg, float* weight_image, float* align_image, 
		    MRI_Dataset* Align, MRI_Dataset* Stdv, 
		    int dx, int dy, int dz, 
		    unsigned char** alignMissing, unsigned char** stdvMissing);

