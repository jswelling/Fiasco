/************************************************************
 *                                                          *
 *  kalmanfilter.h                                          *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *                                                          *
 *  Copyright (c) 2008 Pittsburgh Supercomputing Center     *
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
 *  Original programming by Joel Welling 1/2008             *
 ************************************************************/
/* Header file for kalmanfilter.c */

#ifndef INCL_KALMANFILTER_H
#define INCL_KALMANFILTER_H 1

/* Column major order, a.k.a. Fortran order */
#define DATA_ORDER_COLUMNMAJOR 0
/* Row major order, a.k.a. C order */
#define DATA_ORDER_ROWMAJOR 1

struct KalmanState_struct;

typedef struct KalmanProcess_struct {
  int L;
  int M;
  int debugFlag;
  int dataOrder; /* COLUMNMAJOR or ROWMAJOR */
  double* A; /* owned by the object; size M*M doubles */
  double* H; /* owned by the object; size L*M doubles */
  double* Q; /* owned by the object; size M*M doubles */
  double* R; /* owned by the object; size L*L doubles */
  double* scratch; /* owned by the object */
  long scratchSize;
  int* iScratch;   /* owned by the object */
  long iScratchSize;
  int validA, validH, validQ, validR;

  void (*destroySelf)( struct KalmanProcess_struct* self );
  void (*dumpSelf)( const struct KalmanProcess_struct* self, FILE* ofile );
  void (*setDebug)(struct KalmanProcess_struct* self, int val);
  int  (*getDebug)(struct KalmanProcess_struct* self);
  int  (*getDataOrder)(struct KalmanProcess_struct* self);

  /* sizeA must be >= M*M */
  void (*setA)( struct KalmanProcess_struct* self, 
		const double* AIn, int sizeA );
  void (*getA)( const struct KalmanProcess_struct* self, 
		double* AOut, int sizeA );
  /* sizeH must be >= L*M */
  void (*setH)( struct KalmanProcess_struct* self, 
		const double* HIn, int sizeH );
  void (*getH)( const struct KalmanProcess_struct* self, 
		double* HOut, int sizeH );
  /* sizeQ must be >= M*M */
  void (*setQ)( struct KalmanProcess_struct* self, 
		const double* QIn, int sizeQ );
  void (*getQ)( const struct KalmanProcess_struct* self, 
		double* QOut, int sizeQ );
  /* sizeR must be >= L*L */
  void (*setR)( struct KalmanProcess_struct* self, 
		const double* RIn, int sizeR );
  void (*getR)( const struct KalmanProcess_struct* self, 
		double* ROut, int sizeR );

  /* apply method returns deltaLogLikelihood */
  /* sizeZ must match or exceed L */
  double (*apply)( const struct KalmanProcess_struct* self,
		   struct KalmanState_struct* state, 
		   const double* z, int sizeZ, int missing, 
		   int calcLogLikelihoodDelta, int updatePandK );

} KalmanProcess;

typedef struct KalmanState_struct {
  int M;
  int dataOrder; /* COLUMNMAJOR or ROWMAJOR */
  double* x; /* hidden state, in caller's memory */
  double* P; /* Error covariance, in caller's memory */
  void (*destroySelf)( struct KalmanState_struct* self );
  void (*dumpSelf)( const struct KalmanState_struct* self, FILE* ofile );
  int  (*getDataOrder)(struct KalmanState_struct* self);
  
  /* sizeX must be >= M */
  void (*setX)( struct KalmanState_struct* self, double* x, int sizeX );
  /* sizeP must be >= M*M */
  void (*setP)( struct KalmanState_struct* self, double* P, int sizeP );
} KalmanState;

typedef struct KalmanFilter_struct {
  int L; /* learned from the input Process */
  int M; /* learned from the input Process */
  KalmanProcess* process;
  KalmanState* state;
  void (*destroySelf)( struct KalmanFilter_struct* self );
  void (*dumpSelf)( const struct KalmanFilter_struct* self, FILE* ofile );
  int  (*getDataOrder)(struct KalmanFilter_struct* self);
  /* Returns log likelihood */
  double (*run)( struct KalmanFilter_struct* self,
		 long tdim,
		 const double* zIn, long sizeZ,
		 double* xOut, long sizeXOut,
		 double* POut, long sizePOut,
		 int calcLogLikelihood );
} KalmanFilter;

extern KalmanProcess* klmn_createKalmanProcess(int L, int M);
extern KalmanProcess* klmn_createRowMajorKalmanProcess(int L, int M);
extern KalmanState* klmn_createKalmanState(int M);
extern KalmanState* klmn_createRowMajorKalmanState(int M);

/* The Process becomes the property of the new Filter, and is deleted
 * with it.
 */
extern KalmanFilter* klmn_createKalmanFilter(KalmanProcess* process);

#endif
