/************************************************************
 *                                                          *
 *  mri_kalman.c                                            *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1997 Department of Statistics,         *
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
 *  Original programming by Joel Welling, 6/03              *
 ************************************************************/
/*************************************************************

  DESCRIPTION OF MRI_KALMAN

  mri_kalman performs a Kalman filter on a Pgh MRI dataset.
**************************************************************/

/*
 * Notes
 * -add a reduction matrix, to save writing out all the estimates?
 */

#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "misc.h"
#include "stdcrg.h"
#include "../fmri/lapack.h"

#define KEYBUF_SIZE 512

static char rcsid[] = "$Id: mri_kalman.c,v 1.7 2008/02/12 01:13:43 welling Exp $";

static char* progname;
static int verbose_flg= 0;
static int debug_flg= 0;

#define MAYBE_MAKE_HASHMARK( i, n ) if (verbose_flg) makeHashMark(i,n)

static void makeHashMark( long i, long n )
{
  if (i==0) Message("      ");
  if ( (i==n-1) || (i+1)%60 == 0 ) {
    Message( "# %ld\n", i+1 );
    if ( i != n-1 ) Message( "      " );
  }
  else Message("#");
}

static void safe_copy(char* str1, const char* str2) {
  strncpy(str1, str2, KEYBUF_SIZE);
  str1[KEYBUF_SIZE-1]= '\0';
}

static void safe_concat(char* str1, const char* str2) {
  strncat(str1, str2, (KEYBUF_SIZE-strlen(str1))-1);
}

static int safe_get_extent(MRI_Dataset* ds, const char* chunk, 
			   const char dim)
{
  char key_buf[KEYBUF_SIZE];
  char dim_buf[4];
  dim_buf[0]= dim;
  dim_buf[1]= '\0';
  safe_copy(key_buf,chunk);
  safe_concat(key_buf,".extent.");
  safe_concat(key_buf,dim_buf);
  if (mri_has(ds,key_buf)) return mri_get_int(ds,key_buf);
  else Abort("%s: input missing tag %s!\n",progname,key_buf);
  return 0; /* not reached */
}

static char* safe_get_dims(MRI_Dataset* ds, const char* chunk)
{
  char key_buf[KEYBUF_SIZE];
  safe_copy(key_buf,chunk);
  safe_concat(key_buf,".dimensions");
  if (mri_has(ds,key_buf)) return mri_get_string(ds,key_buf);
  else Abort("%s: input missing tag %s!\n",progname,key_buf);
  return 0; /* not reached */
}

#ifdef never
static char* matrix_to_string(double* A, long a, long b)
{
  static char buf[512];
  char* here= buf;
  char* end= buf+(sizeof(buf)-1);
  int i;
  int j;
  for (i=0; i<a; i++) {
    for (j=0; j<b; j++) {
      snprintf(here,(int)(end-here),"%g ",A[j*a+i]);
      here += strlen(here);
    }
    snprintf(here,(int)(end-here),"\n");
    here += strlen(here);
  }
  return buf;
}

/* returns increment to log likelihood */
static double kalman_iterate_once( const double* A, const double* H,
				   const double* Q, const double* R,
				   const double* z, double* xHat, double* P,
				   long L, long M, int missingFlag,
				   int logLikelihoodFlag, int stepPandKFlag,
				   double* scratch, long scratchSize,
				   int* iScratch, long iScratchSize )
{
  double* xHatPrior= scratch;
  double* PPrior= scratch+M;
  double* K= PPrior + (M*M);
  double* tmpMM= K + (L*M);
  double* tmpLM= tmpMM + (M*M);
  double* tmpLL= tmpLM + (L*M); /* inverse of innovation covariance */
  double* tmpLL2= tmpLL + (L*L);
  double* innovation= tmpLL2 + (L*L);
  double* tmpL= innovation + L;
  double* determinantOfInnovCov= tmpL + L;
  double one= 1.0;
  double neg_one= -1.0;
  double zero= 0.0;
  double logLikelihoodIncrement= 0.0;
  int int_one= 1;
  int iL= (int)L;
  int iM= (int)M;
  int iLL= (int)(L*L);
  int info= 0;
  long i;

  if (scratchSize<2*L+M+2*(M*M)+2*(L*M)+2*(L*L)+1) 
    Abort("%s: internal error; scratch space too small!\n",progname);
  if (iScratchSize<L)
    Abort("%s: internal error; int scratch space too small!\n",progname);

  /* xHat and P are at time t-1; we mean to step them to time t */
  /* xHatPrior in scratch, size M; PPrior in scratch, size M */

  if (debug_flg) {
    fprintf(stderr,"A= %s\n",matrix_to_string((double*)A,M,M));
    fprintf(stderr,"H= %s\n",matrix_to_string((double*)H,L,M));
    fprintf(stderr,"Q= %s\n",matrix_to_string((double*)Q,M,M));
    fprintf(stderr,"R= %s\n",matrix_to_string((double*)R,L,L));
    fprintf(stderr,"z= %s\n",matrix_to_string((double*)z,L,1));
    fprintf(stderr,"xHat^T= %s\n",matrix_to_string(xHat,1,M));
    fprintf(stderr,"P= %s\n",matrix_to_string(P,M,M));
  }
  if (stepPandKFlag) {
    /* PPrior= A*P*transpose(A) + Q;  MxM * MxM * MxM */
    DGEMM("n","n",&iM,&iM,&iM,&one,(double*)A,&iM,P,&iM,&zero,tmpMM,&iM);
    bcopy(Q,PPrior,M*M*sizeof(double));
    DGEMM("n","t",&iM,&iM,&iM,&one,tmpMM,&iM,(double*)A,&iM,
	  &one,PPrior,&iM);
    if (debug_flg) {
      fprintf(stderr,"PPrior= %s\n",matrix_to_string(PPrior,M,M));
    }

    if (missingFlag) {
      for (i=0; i<L*M; i++) K[i]= 0.0;
      *determinantOfInnovCov= 1.0;
    }
    else {
      /* K= PPrior*transpose(H)*inverse(H*PPrior*transpose(H) + R); */
      /* The thing that gets inverted is the innovation covariance */
      /* Start by computing tmpLL= H*PPrior*transpose(H) + R */
      DGEMM("n","n",&iL,&iM,&iM,&one,(double*)H,&iL,PPrior,&iM,
	    &zero,tmpLM,&iL);
      bcopy(R,tmpLL,L*L*sizeof(double));
      DGEMM("n","t",&iL,&iL,&iM,&one,tmpLM,&iL,(double*)H,&iL,
	    &one,tmpLL,&iL);
      /* Compute inverse(tmpLL) in place */
      DGETRF( &iL, &iL, tmpLL, &iL, iScratch, &info );
      if (info!=0) 
	Abort("%s: matrix inversion failed computing K; err %d on DGETRF!\n",
	      progname,info);
      /* Snag the determinant while it's handy */
      *determinantOfInnovCov= 1.0;
      for (i=0; i<L*L; i+=(L+1))
	*determinantOfInnovCov *= tmpLL[i];
      DGETRI( &iL, tmpLL, &iL, iScratch, tmpLL2, &iLL, &info );
      if (info!=0) 
	Abort("%s: matrix inversion failed computing K; err %d on DGETRI!\n",
	      progname,info);
      /* Compute K= PPrior*transpose(H)*tmpLL */
      DGEMM("n","t",&iM, &iL, &iM, &one, PPrior, &iM, (double*)H, &iL,
	    &zero,tmpLM,&iM);
      DGEMM("n","n",&iM, &iL, &iL, &one, tmpLM, &iM, tmpLL, &iL,
	    &zero, K, &iM);
      /* !!!NOTE!!! we are keeping tmpLL around for possible use in
       * computing the log likelihood increment.  It contains the
       * inverse innovation covariance.
       */
    }
    if (debug_flg) fprintf(stderr,"K: %s\n",matrix_to_string(K,M,L));

    /* P= (1-K*H)*PPrior; */
    /* or more conveniently, P= -1*K*H*PPrior + 1*PPrior; */
    DGEMM("n","n",&iM, &iM, &iL, &one, K, &iM, (double*)H, &iL,
	  &zero, tmpMM, &iM);
    bcopy(PPrior, P, M*M*sizeof(double));
    DGEMM("n","n",&iM, &iM, &iM, &neg_one, tmpMM, &iM, PPrior, &iM,
	  &one, P, &iM);
  }

  /* xHatPrior= A*xHat;  MxM * M  (xHat has not yet been updated) */
  DGEMV("n", &iM, &iM, &one, (double*)A, &iM, xHat, &int_one, &zero, 
	xHatPrior, &int_one);
  if (debug_flg) {
    fprintf(stderr,"xHatPrior^T= %s\n",matrix_to_string(xHatPrior,1,M));
  }
  
  /* xHat= xHatPrior + K*(z - H*xHatPrior); */
  /* (z - H*xHatPrior) is the 'innovation' */
  bcopy(z, innovation, L*sizeof(double));
  DGEMV("n", &iL, &iM, &neg_one, (double*)H, &iL, xHatPrior, &int_one,
	&one, innovation, &int_one);
  bcopy(xHatPrior, xHat, M*sizeof(double));
  DGEMV("n", &iM, &iL, &one, K, &iM, innovation, &int_one,
	&one, xHat, &int_one);

  if (logLikelihoodFlag && !missingFlag) {
    /* log(f)= -(1/2) log[ (2pi)^L * det(H*PPrior*H^T + R) ]
     *
     *        -(1/2)(z-H*xHat)^T * (H*PPrior*H^T + R)^-1 * (z-H*xHat)
     */
    double dot= 0.0;
    DGEMV("n", &iL, &iL, &one, tmpLL, &iL, innovation, &iL, 
	  &zero, tmpL, &iL );
    for (i=0; i<L; i++) dot += innovation[i]*tmpL[i];
    logLikelihoodIncrement= 
      -0.5*(L*log(2.0*M_PI) + log(*determinantOfInnovCov) + dot);
    if (debug_flg) fprintf(stderr,"logLikelihoodIncrement= %g\n",
			   logLikelihoodIncrement);
  }
  else {
    /* save some pointless computation */
    logLikelihoodIncrement= 0.0;
  }
  /* And we're done. */
  return logLikelihoodIncrement;
}

/* returns increment to log likelihood */
/* Special case for the common L==1 */
static double kalman_iterate_once_Leq1( const double* A, const double* H,
					const double* Q, const double* R,
					const double* z, double* xHat, 
					double* P, long M, int missingFlag,
					int logLikelihoodFlag, 
					int stepPandKFlag,
					double* scratch, long scratchSize,
					int* iScratch, long iScratchSize )
{
  double* xHatPrior= scratch;
  double* PPrior= scratch+M;
  double* K= PPrior + (M*M);
  double* tmpMM= K + M;
  double* tmpLM= tmpMM + (M*M);
  double* tmpLL= tmpLM + M; /* inverse of innovation covariance */
  double* tmpLL2= tmpLL + 1;
  double* innovation= tmpLL2 + 1;
  double* tmpL= innovation + 1;
  double* determinantOfInnovCov= tmpL + 1;
  double one= 1.0;
  double neg_one= -1.0;
  double zero= 0.0;
  double logLikelihoodIncrement= 0.0;
  int int_one= 1;
  int iM= (int)M;
  int info= 0;
  long i;

  if (scratchSize<2+M+2*(M*M)+2*M+2+1) 
    Abort("%s: internal error; scratch space too small!\n",progname);
  if (iScratchSize<1)
    Abort("%s: internal error; int scratch space too small!\n",progname);

  /* xHat and P are at time t-1; we mean to step them to time t */
  /* xHatPrior in scratch, size M; PPrior in scratch, size M */

  if (debug_flg) {
    fprintf(stderr,"A= %s\n",matrix_to_string((double*)A,M,M));
    fprintf(stderr,"H= %s\n",matrix_to_string((double*)H,1,M));
    fprintf(stderr,"Q= %s\n",matrix_to_string((double*)Q,M,M));
    fprintf(stderr,"R= %g\n",R[0]);
    fprintf(stderr,"z= %g\n",z[0]);
    fprintf(stderr,"xHat^T= %s\n",matrix_to_string(xHat,1,M));
    fprintf(stderr,"P= %s\n",matrix_to_string(P,M,M));
  }
  if (stepPandKFlag) {
    /* PPrior= A*P*transpose(A) + Q;  MxM * MxM * MxM */
    DGEMM("n","n",&iM,&iM,&iM,&one,(double*)A,&iM,P,&iM,&zero,tmpMM,&iM);
    bcopy(Q,PPrior,M*M*sizeof(double));
    DGEMM("n","t",&iM,&iM,&iM,&one,tmpMM,&iM,(double*)A,&iM,
	  &one,PPrior,&iM);
    if (debug_flg) {
      fprintf(stderr,"PPrior= %s\n",matrix_to_string(PPrior,M,M));
    }

    if (missingFlag) {
      for (i=0; i<M; i++) K[i]= 0.0;
      *determinantOfInnovCov= 1.0;
    }
    else {
      /* K= PPrior*transpose(H)*inverse(H*PPrior*transpose(H) + R); */
      /* The thing that gets inverted is the innovation covariance */
      /* Start by computing tmpLL= H*PPrior*transpose(H) + R */
      DGEMV("n", &iM, &iM, &one, PPrior, &iM, (double*)H, &int_one, &zero,
	    tmpLM, &int_one);
      tmpLL[0]= R[0];
      for (i=0; i<M; i++) tmpLL[0] += H[i]*tmpLM[i];
      *determinantOfInnovCov= tmpLL[0];
      /* Compute inverse(tmpLL) in place */
      tmpLL[0]= 1.0/tmpLL[0];
      /* Compute K= PPrior*transpose(H)*tmpLL */
      DGEMV("n", &iM, &iM, tmpLL, PPrior, &iM, (double*)H, &int_one, &zero, K, 
	    &int_one);
      /* !!!NOTE!!! we are keeping tmpLL around for possible use in
       * computing the log likelihood increment.  It contains the
       * inverse innovation covariance.
       */
    }
    if (debug_flg) fprintf(stderr,"K: %s\n",matrix_to_string(K,M,1));

    /* P= (1-K*H)*PPrior; */
    /* or more conveniently, P= -1*K*H*PPrior + 1*PPrior; */
    DGEMM("n","n",&iM, &iM, &int_one, &one, K, &iM, (double*)H, &int_one,
	  &zero, tmpMM, &iM);
    bcopy(PPrior, P, M*M*sizeof(double));
    DGEMM("n","n",&iM, &iM, &iM, &neg_one, tmpMM, &iM, PPrior, &iM,
	  &one, P, &iM);
  }

  /* xHatPrior= A*xHat;  MxM * M  (xHat has not yet been updated) */
  DGEMV("n", &iM, &iM, &one, (double*)A, &iM, xHat, &int_one, &zero, 
	xHatPrior, &int_one);
  if (debug_flg) {
    fprintf(stderr,"xHatPrior^T= %s\n",matrix_to_string(xHatPrior,1,M));
  }
  
  /* xHat= xHatPrior + K*(z - H*xHatPrior); */
  /* (z - H*xHatPrior) is the 'innovation' */
  innovation[0]= z[0];
  for (i=0; i<M; i++) innovation[0] -= H[i]*xHatPrior[i];
  for (i=0; i<M; i++) xHat[i]= xHatPrior[i] + K[i]*innovation[0];

  if (logLikelihoodFlag && !missingFlag) {
    /* log(f)= -(1/2) log[ (2pi)^L * det(H*PPrior*H^T + R) ]
     *
     *        -(1/2)(z-H*xHat)^T * (H*PPrior*H^T + R)^-1 * (z-H*xHat)
     */
    double dot= innovation[0]*tmpLL[0]*innovation[0];
    logLikelihoodIncrement= 
      -0.5*(log(2.0*M_PI) + log(*determinantOfInnovCov) + dot);
    if (debug_flg) fprintf(stderr,"logLikelihoodIncrement= %g\n",
			   logLikelihoodIncrement);
  }
  else {
    /* save some pointless computation */
    logLikelihoodIncrement= 0.0;
  }
  /* And we're done. */
  return logLikelihoodIncrement;
}
#endif

static void apply_kalman( MRI_Dataset* input, 
			  MRI_Dataset* transMatrix, 
			  MRI_Dataset* obsMatrix, 
			  MRI_Dataset* procCovMatrix, 
			  MRI_Dataset* obsCovMatrix,
			  MRI_Dataset* output, 
			  MRI_Dataset* errCovOut,
			  MRI_Dataset* logLikelihoodOut,
			  const char* chunk,
			  const long L, const long M, 
			  const long dt, const long dz,
			  const long samplesPerSlice,
			  const long obsMatrixStride,
			  const long transMatrixStride )
{
  long samp= 0;
  long z= 0;
  long t= 0;
  long s= 0;
  long i= 0;
  double* procCovBuf= NULL;
  double* obsCovBuf= NULL;
  long inputStride= L*samplesPerSlice;
  long long inputOffset= 0;
  long outputStride= M*samplesPerSlice;
  long errCovStride= M*M;
  long logLikelihoodStride= samplesPerSlice;
  long long outputOffset= 0;
  long long errCovOffset= 0;
  long long logLikelihoodOffset= 0;
  long long transMatrixOffset= 0;
  long long obsMatrixOffset= 0;
  long scratchSize= 2*L+M+2*(M*M)+2*(L*M)+2*(L*L)+1;
  long iScratchSize= L;
  double* inputFrame= NULL;
  double* outputFrame= NULL;
  double* errCovFrame= NULL;
  double* logLikelihoodFrame= NULL;
  double* transMatrixFrame= NULL;
  double* obsMatrixFrame= NULL;
  double* scratch= NULL;
  int* iScratch= NULL;
  unsigned char** missing= get_missing(input);
  KalmanProcess* proc= NULL;
  KalmanState* state= NULL;

  if (!(proc= klmn_createKalmanProcess(L,M)))
    Abort("%s unable to create KalmanProcess instance!\n",
	  progname);
  if (!(state= klmn_createKalmanState(M)))
    Abort("%s unable to create KalmanState instance!\n",
	  progname);

  if (!(procCovBuf=(double*)malloc(M*M*sizeof(double))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,M*M*sizeof(double));
  if (!(obsCovBuf=(double*)malloc(L*L*sizeof(double))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,L*L*sizeof(double));

  if (!(inputFrame=(double*)malloc(inputStride*sizeof(double))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,inputStride*sizeof(double));
  if (!(outputFrame=(double*)malloc(outputStride*sizeof(double))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,outputStride*sizeof(double));
  if (!(errCovFrame=(double*)malloc(errCovStride*sizeof(double))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,errCovStride*sizeof(double));
  if (!(logLikelihoodFrame=
	(double*)malloc(logLikelihoodStride*sizeof(double))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,logLikelihoodStride*sizeof(double));

  if (!(transMatrixFrame=(double*)malloc(M*M*sizeof(double))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,M*M*sizeof(double));
  if (!(obsMatrixFrame=(double*)malloc(L*M*sizeof(double))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,L*M*sizeof(double));

#ifdef never
  if (!(scratch=(double*)malloc(scratchSize*sizeof(double))))
    Abort("%s: unable to allocate %d bytes!\n",
	  scratchSize*sizeof(double));
  if (!(iScratch=(int*)malloc(iScratchSize*sizeof(int))))
    Abort("%s: unable to allocate %d bytes!\n",
	  iScratchSize*sizeof(int));
#endif

  mri_read_chunk( procCovMatrix, chunk, M*M, 0, MRI_DOUBLE, procCovBuf );
  proc->setQ(proc,procCovBuf, (int)(M*M));
  mri_read_chunk( obsCovMatrix, chunk, L*L, 0, MRI_DOUBLE, obsCovBuf );
  proc->setR(proc,obsCovBuf, L*L);

  for (z=0; z<dz; z++) {
    double* inputSample= NULL;
    double* outputSample= NULL;
    double* errCovSample= NULL; 
    double* logLikelihoodSample= NULL;

    if (debug_flg) fprintf(stderr,"Starting z= %ld\n",z);
    inputOffset= z*inputStride;
    outputOffset= z*outputStride;
    errCovOffset= z*errCovStride;
    logLikelihoodOffset= z*logLikelihoodStride;
    obsMatrixOffset= z*obsMatrixStride;
    transMatrixOffset= z*transMatrixStride;

    for (t=0; t<dt; t++) {

      mri_read_chunk( input, chunk, inputStride, inputOffset,
		      MRI_DOUBLE, inputFrame );
      if (t==0 || obsMatrixStride!=0) {
	mri_read_chunk( obsMatrix, chunk, L*M, obsMatrixOffset,
			MRI_DOUBLE, obsMatrixFrame );
	proc->setH(proc,obsMatrixFrame,L*M);
      }
      if (t==0 || transMatrixStride!=0) {
	mri_read_chunk( transMatrix, chunk, M*M, transMatrixOffset,
			MRI_DOUBLE, transMatrixFrame );
	proc->setA(proc,transMatrixFrame, M*M);
      }
      
      inputSample= inputFrame;
      outputSample= outputFrame;
      errCovSample= errCovFrame;
      logLikelihoodSample= logLikelihoodFrame;
      for (s=0; s<samplesPerSlice; s++) { 
	double val;

	if (t==0) {
	  /* outputSample gets initialized to all 1's. */
	  for (i=0; i<M; i++) outputSample[i]= 1.0;

	  /* errCovSample gets initialized to procCovMatrix */
	  bcopy(procCovBuf,errCovSample,M*M*sizeof(double));
	  state->setP(state, errCovSample, M*M);

	  /* logLikelihoodSample gets initialized to 0.0 */
	  *logLikelihoodSample= 0.0;
	}

	/* Solve one iteration, and propogate the results forward */
	if (debug_flg)
	  fprintf(stderr,"-------------\n samp %ld t=%ld\n-------------\n",
		  s,t);
	
	/* outputSample and errCovSample are updated through their
	 * references in the KalmanState. 
	 */
	state->setX( state, outputSample, M );
	val= proc->apply( proc, state, inputSample, L, missing[t][z],
			  (logLikelihoodOut != NULL), (s==0) );

	if (logLikelihoodOut)
	  *logLikelihoodSample += val;

	inputSample += L;
	outputSample += M;
	logLikelihoodSample += 1;
      }

      if (output)
	mri_write_chunk( output, chunk, outputStride, outputOffset, 
			 MRI_DOUBLE, outputFrame );
      if (errCovOut) 
	mri_write_chunk( errCovOut, chunk, errCovStride, errCovOffset, 
			 MRI_DOUBLE, errCovFrame );

      if (verbose_flg) makeHashMark( z*dt + t, dz*dt );

      inputOffset += dz*inputStride;
      outputOffset += dz*outputStride;
      errCovOffset += dz*errCovStride;
      obsMatrixOffset += dz*obsMatrixStride;
      transMatrixOffset += dz*transMatrixStride;
    }
    if (logLikelihoodOut) 
      mri_write_chunk( logLikelihoodOut, chunk, logLikelihoodStride, 
		       logLikelihoodOffset, 
		       MRI_DOUBLE, logLikelihoodFrame );
  }
  
#ifdef never
  free(iScratch);
  free(scratch);
#endif
  free(obsMatrixFrame);
  free(transMatrixFrame);
  free(logLikelihoodFrame);
  free(errCovFrame);
  free(outputFrame);
  free(inputFrame);
  free(obsCovBuf);
  free(procCovBuf);
  state->destroySelf(state);
  proc->destroySelf(proc);
  
}

static int chunk_check( MRI_Dataset* ds, const char* chunk )
{
  return( mri_has(ds, chunk) 
	  && !strcmp(mri_get_string(ds,chunk),"[chunk]") );
}

static int structure_check_obsmatrix( MRI_Dataset* ds, 
				      const char* chunk, 
				      long L, long* pM,
				      long dz, long dt,
				      long* pObsMatrixStride )
{
  char* dimstr= safe_get_dims(ds,chunk);
  int ndims= strlen(dimstr);

  /*
  The dimension string for obsMatrix must be either two characters or
  four characters; if the string is four characters long the last two
  must be zt in that order.  The second dimension must match the first
  dimension of Infile in name and extent, so for the current example
  it must be a and have extent L.  The first dimension may be any
  unique character.  For the sake of this example, assume it is b and
  that it has extent M.  M will be the dimensionality of the hidden
  state of the Kalman process.
  */

  if (ndims==2) {
    if (safe_get_extent(ds,chunk,dimstr[1]) != L) return 0;
    *pM= safe_get_extent(ds,chunk,dimstr[0]);
    *pObsMatrixStride= 0;
  }
  else if (ndims==4) {
    if (safe_get_extent(ds,chunk,dimstr[1]) != L) return 0;
    if (dimstr[2] != 'z' || safe_get_extent(ds,chunk,'z') != dz)
      return 0;
    if (dimstr[3] != 't' || safe_get_extent(ds,chunk,'t') != dt)
      return 0;
    *pM= safe_get_extent(ds,chunk,dimstr[0]);
    *pObsMatrixStride= *pM * L;
  }
  else return 0;

  return 1;
}

static int structure_check_transmatrix( MRI_Dataset* ds, 
					const char* chunk, 
					long M, long dz, long dt,
					long* pTransMatrixStride )
{
  char* dimstr= safe_get_dims(ds,chunk);
  int ndims= strlen(dimstr);

  /*
  The dimension string of transMatrix must be either two characters or
  four characters; if the string is four characters long the last two
  must be zt in that order.  The first two dimensions must have the
  same extent, and that extent must be M, the extent of the second
  dimension of obsMatrix.  Remember that the leftmost index in a Pgh
  MRI datafile varies fastest, and that this is contrary to the
  convention for writing down algebraic matrices.  Thus the sequence
  of values in transMatrix is the same sequence (reading left-right
  then up-down) as that in the transpose of the transition matrix as
  it is usually written algebraically.
  */

  if (ndims==2) {
    if (safe_get_extent(ds,chunk,dimstr[0]) != M) return 0;
    if (safe_get_extent(ds,chunk,dimstr[1]) != M) return 0;
    *pTransMatrixStride= 0;
  }
  else if (ndims==4) {
    if (safe_get_extent(ds,chunk,dimstr[0]) != M) return 0;
    if (safe_get_extent(ds,chunk,dimstr[1]) != M) return 0;
    if (dimstr[2] != 'z' || safe_get_extent(ds,chunk,'z') != dz)
      return 0;
    if (dimstr[3] != 't' || safe_get_extent(ds,chunk,'t') != dt)
      return 0;
    *pTransMatrixStride= M*M;
  }
  else return 0;

  return 1;
}

static int structure_check_proccovmatrix( MRI_Dataset* ds, 
					  const char* chunk, 
					  long M )
{
  char* dimstr= safe_get_dims(ds,chunk);
  int ndims= strlen(dimstr);

  /*
  The procCovMatrix and obsCovMatrix inputs must have two-character
  dimension strings, and the extent of both dimensions must match the
  corresponding values of transMatrix and obsMatrix respectively.
  Thus for the current example, the extent of both dimensions for
  procCovMatrix must be M and the extent of both dimensions for
  obsCovMatrix must be L.  Varying these quantities by time (with
  additional dimensions zt) is not supported.
  */

  if (ndims != 2 
      || safe_get_extent(ds,chunk,dimstr[0]) != M
      || safe_get_extent(ds,chunk,dimstr[1]) != M)
    return 0;

  return 1;
}

static int structure_check_obscovmatrix( MRI_Dataset* ds, 
					 const char* chunk, 
					 long L )
{
  char* dimstr= safe_get_dims(ds,chunk);
  int ndims= strlen(dimstr);

  /*
  The procCovMatrix and obsCovMatrix inputs must have two-character
  dimension strings, and the extent of both dimensions must match the
  corresponding values of transMatrix and obsMatrix respectively.
  Thus for the current example, the extent of both dimensions for
  procCovMatrix must be M and the extent of both dimensions for
  obsCovMatrix must be L.  Varying these quantities by time (with
  additional dimensions zt) is not supported.
  */
  if (ndims != 2 
      || safe_get_extent(ds,chunk,dimstr[0]) != L
      || safe_get_extent(ds,chunk,dimstr[1]) != L)
    return 0;

  return 1;
}

static int structure_check_input( MRI_Dataset* ds, const char* chunk, 
				  long* pL, long* pdz, long* pdt,
				  long* pSamplesPerSlice )
{
  char* dimstr= safe_get_dims(ds,chunk);
  int ndims= strlen(dimstr);
  int i;
  long samples= 1;

  /*
  Assume that Infile has dimensions aXzt, where X is any combination
  of letters, and that the extent of a is L.  The dimensions X and z
  enumerate the pixels in the input dataset; each pixel is filtered
  separately.  The program will iterate the Kalman filter a number of
  times equal to the extent of t.

  Missing information in the input file is applied if present; the
  Kalman weight is set to zero for missing input data.  Because of the
  definition of the missing data chunk, missing information must
  change only by slice and by time.
  */

  if (ndims<3) return 0;
  *pL= safe_get_extent(ds,chunk,dimstr[0]);
  if (dimstr[ndims-2]=='z') 
    *pdz= safe_get_extent(ds,chunk,dimstr[ndims-2]);
  else return 0;
  if (dimstr[ndims-1]=='t') 
    *pdt= safe_get_extent(ds,chunk,dimstr[ndims-1]);
  else return 0;
  samples= 1;
  for (i=1; i<ndims-2; i++) {
    samples *= safe_get_extent(ds,chunk,dimstr[i]);
  }
  *pSamplesPerSlice= samples;

  return 1;
}

static void removeMissingChunk( MRI_Dataset* ds )
{
  /* If a missing chunk exists in the dataset, remove it and associated
   * keys since the notion of 'missing' data is not well defined for 
   * extent.t==1.
   */
  if (mri_has(ds, "missing") 
      && !strcmp(mri_get_string(ds,"missing"),"[chunk]")) {
    char* key= NULL;
    int len= strlen("missing");
    mri_iterate_over_keys(ds);
    while ((key = mri_next_key(ds)) != NULL)
      if (strncmp(key, "missing", len) == 0 &&
	  (key[len] == '\0' || key[len] == '.'))
	mri_remove(ds, key);
  }
}

static void setChunkDataType( MRI_Dataset* ds, const char* chunk )
{
  char key_buf[KEYBUF_SIZE];
  /*
   * We'll use float32 as the datatype unless
   * it's already double precision.
   */
  snprintf(key_buf,sizeof(key_buf),"%s.datatype",chunk);
  key_buf[sizeof(key_buf)-1]= '\0';
  if (strcmp(mri_get_string(ds,key_buf),"float64"))
    mri_set_string(ds,key_buf,"float32");
}

static MRI_Dataset* create_output_ds( const char* fname, 
				      MRI_Dataset* inDS, 
				      const char* chunk,
				      long M )
{
  MRI_Dataset* result= NULL;
  char key_buf[KEYBUF_SIZE];
  const char* dimstr= safe_get_dims(inDS,chunk);

  /*
  The output dataset Outfile will contain a chunk with dimensions
  matching Infile, except that the extent of the first dimension will
  be M.
  */

  result= mri_copy_dataset(fname,inDS);

  snprintf(key_buf,sizeof(key_buf),"%s.extent.%c",chunk,dimstr[0]);
  key_buf[sizeof(key_buf)-1]= '\0';
  mri_set_int( result,key_buf,M);

  setChunkDataType(result, chunk);
  removeMissingChunk(result);

  return result;
}

static MRI_Dataset* create_errcov_ds( const char* fname, 
				      MRI_Dataset* inDS, 
				      const char* chunk, long M )
{
  MRI_Dataset* result= NULL;
  char key_buf[KEYBUF_SIZE];
  char scratch_buf[256];
  const char* dimstr= safe_get_dims(inDS,chunk);
  char c;
  const char* here;

  /*
  The error covariance output errCovOut will only be produced if it is
  specified on the command line.  If so, it will have dimensions
  matching those of xOut, with an additional leading dimension of
  extent M and all in-plane dimensions set to extent 1.  (Since the
  error covariance does not depend on the sample data it does not vary
  by sample).
  */

  result= mri_copy_dataset(fname,inDS);

  snprintf(key_buf,sizeof(key_buf),"%s.extent.%c",chunk,dimstr[0]);
  key_buf[sizeof(key_buf)-1]= '\0';
  mri_set_int( result,key_buf,M);

  /* Set extents of all in-slice dimensions to 1 */
  for (here=dimstr+1; *here != 'z'; here++) {
    snprintf(key_buf,sizeof(key_buf),"%s.extent.%c",
	     chunk, *here);
    mri_set_int(result, key_buf, 1);
  }
  
  /* Find a unique lower-case letter.  This is assisted by the
   * fact that the ASCII codes for the letters are contiguous.
   */
  for (c='a'; c<='z'; c++) 
    if (!strchr(dimstr,c)) break;
  if (c>'z') 
    Abort("%s: couldn't find a unique dimension index for err covariance!\n",
	  progname);

  snprintf(key_buf,sizeof(key_buf),"%s.dimensions",chunk);
  key_buf[sizeof(key_buf)-1]= '\0';
  snprintf(scratch_buf,sizeof(scratch_buf),"%c%s",c,dimstr);
  scratch_buf[sizeof(scratch_buf)-1]= '\0';
  mri_set_string(result,key_buf,scratch_buf);

  snprintf(key_buf,sizeof(key_buf),"%s.extent.%c",chunk,c);
  key_buf[sizeof(key_buf)-1]= '\0';
  mri_set_int(result,key_buf,M);

  setChunkDataType(result, chunk);
  removeMissingChunk(result);

  return result;
}

static MRI_Dataset* create_loglike_ds( const char* fname, 
				       MRI_Dataset* inDS, 
				       const char* chunk )
{
  MRI_Dataset* result= NULL;
  char key_buf[KEYBUF_SIZE];
  char scratch_buf[256];
  const char* dimstr= safe_get_dims(inDS,chunk);
  char c;

  /*
  The log likelihood output logLikelihoodOut will only be produced if
  it is specified on the command line.  If so, it will have dimensions
  matching those of Outfile, but with the leading dimension and the t
  dimension both having extent 1.
  */

  result= mri_copy_dataset(fname,inDS);

  /* Set leading dimension extent to 1 */
  snprintf(key_buf,sizeof(key_buf),"%s.extent.%c",chunk,dimstr[0]);
  key_buf[sizeof(key_buf)-1]= '\0';
  mri_set_int(result,key_buf,1);

  /* Set t dimension extent to 1 */
  snprintf(key_buf,sizeof(key_buf),"%s.extent.t",chunk);
  key_buf[sizeof(key_buf)-1]= '\0';
  mri_set_int(result,key_buf,1);

  setChunkDataType(result, chunk);
  removeMissingChunk(result);

  return result;
}

int main( int argc, char* argv[] ) 
{
  char transMatrixName[512], obsMatrixName[512];
  char procCovMatrixName[512], obsCovMatrixName[512];
  char errCovOutName[512], logLikelihoodName[512];
  char inName[512], outName[512];
  MRI_Dataset *input= NULL, *output= NULL;
  MRI_Dataset *transMatrix= NULL, *obsMatrix= NULL;
  MRI_Dataset *procCovMatrix= NULL, *obsCovMatrix= NULL;
  MRI_Dataset *errCovOut= NULL, *logLikelihoodOut= NULL;
  char chunk[KEYBUF_SIZE];
  char key_buf[KEYBUF_SIZE];
  int errCovFlag= 0;
  int logLikelihoodFlag= 0;
  int outFlag= 0;
  long L= 0; /* dims in input, the observed quantity */
  long M= 0; /* dims in output, the hidden quantity */
  long dt= 0; /* time extent */
  long dz= 0; /* z extent */
  long samplesPerSlice= 0; /* typically pixels in a slice */
  long obsMatrixStride= 0; /* zero if obsMatrix doesn't vary */
  long transMatrixStride= 0; /* zero if transMatrix doesn't vary */

  progname= argv[0];

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /*** Parse command line ***/

  cl_scan( argc, argv );

  verbose_flg= cl_present("verbose|ver|v");
  debug_flg= cl_present("debug|deb");
  cl_get("chunk|chu|c", "%option %s[%]","images",chunk);
  if (!cl_get("transition|trn","%option %s",transMatrixName)) {
    fprintf(stderr,"%s: Transition matrix file name not given.\n",progname);
    Help( "usage" );
    exit(-1);
  }
  if (!cl_get("observation|obs","%option %s",obsMatrixName)) {
    fprintf(stderr,"%s: Observation matrix file name not given.\n",progname);
    Help( "usage" );
    exit(-1);
  }
  if (!cl_get("proccov|pcv","%option %s",procCovMatrixName)) {
    fprintf(stderr,"%s: Process covariance matrix file name not given.\n",
	    progname);
    Help( "usage" );
    exit(-1);
  }
  if (!cl_get("obscov|ocv","%option %s",obsCovMatrixName)) {
    fprintf(stderr,"%s: Observation covariance matrix file name not given.\n",
	    progname);
    Help( "usage" );
    exit(-1);
  }
  errCovFlag= cl_get("errcov|ecv","%option %s",errCovOutName);
  if (!errCovFlag) errCovOutName[0]= '\0';
  logLikelihoodFlag= cl_get("loglikelihood|lgl","%option %s",
			    logLikelihoodName);
  if (!logLikelihoodFlag) logLikelihoodName[0]= '\0';
  outFlag= cl_get("out","%option %s",outName);
  if (!outFlag) outName[0]= '\0';
  if (!cl_get("", "%s", inName)) {
    fprintf(stderr,"%s: Input file name not given.\n",progname);
    Help( "usage" );
    exit(-1);
  }
  if (cl_cleanup_check()) {
    int i;
    fprintf(stderr,"%s: invalid argument in command line:\n    ",progname);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(-1);
  }

  /*** End command-line parsing ***/
  
  /* Print version number */
  if (verbose_flg) Message( "# %s\n", rcsid );

  /* Check name consistency and open input datasets */
  if ( !strcmp( transMatrixName, obsMatrixName)
       || !strcmp( transMatrixName, procCovMatrixName )
       || !strcmp( transMatrixName, obsCovMatrixName )
       || !strcmp( transMatrixName, errCovOutName )
       || !strcmp( transMatrixName, logLikelihoodName )
       || !strcmp( transMatrixName, inName )
       || !strcmp( transMatrixName, outName )
       || !strcmp( obsMatrixName, procCovMatrixName )
       || !strcmp( obsMatrixName, obsCovMatrixName )
       || !strcmp( obsMatrixName, errCovOutName )
       || !strcmp( obsMatrixName, logLikelihoodName )
       || !strcmp( obsMatrixName, inName )
       || !strcmp( obsMatrixName, outName )
       || !strcmp( procCovMatrixName, obsCovMatrixName )
       || !strcmp( procCovMatrixName, errCovOutName )
       || !strcmp( procCovMatrixName, logLikelihoodName )
       || !strcmp( procCovMatrixName, inName )
       || !strcmp( procCovMatrixName, outName )
       || !strcmp( obsCovMatrixName, errCovOutName )
       || !strcmp( obsCovMatrixName, logLikelihoodName )
       || !strcmp( obsCovMatrixName, inName )
       || !strcmp( obsCovMatrixName, outName )
       || (errCovFlag && !strcmp( errCovOutName, logLikelihoodName ))
       || (errCovFlag && !strcmp( errCovOutName, inName ))
       || (errCovFlag && !strcmp( errCovOutName, outName ))
       || (logLikelihoodFlag && !strcmp( logLikelihoodName, inName ))
       || (logLikelihoodFlag && !strcmp( logLikelihoodName, outName ))
       || (outFlag && !strcmp( inName, outName )) ) 
    Abort( "%s: two or more of the specified filenames are identical.\n",
	   progname );

  /* Open input files */
  input = mri_open_dataset( inName, MRI_READ );
  transMatrix = mri_open_dataset( transMatrixName, MRI_READ );
  obsMatrix = mri_open_dataset( obsMatrixName, MRI_READ );
  procCovMatrix = mri_open_dataset( procCovMatrixName, MRI_READ );
  obsCovMatrix = mri_open_dataset( obsCovMatrixName, MRI_READ );

  /* Will they work with the program? */
#define CHUNK_CHECK( ds, fname ) \
  if (!chunk_check(ds, chunk)) \
    Abort("%s: %s has no chunk %s!\n",progname,fname,chunk)
  CHUNK_CHECK(input, inName);
  CHUNK_CHECK(transMatrix, transMatrixName);
  CHUNK_CHECK(obsMatrix, obsMatrixName);
  CHUNK_CHECK(procCovMatrix, procCovMatrixName);
  CHUNK_CHECK(obsCovMatrix, obsCovMatrixName);
#undef CHUNK_CHECK
  
  if (!structure_check_input(input, chunk, &L, &dz, &dt, 
			     &samplesPerSlice))
    Abort("%s: Dataset %s has an inappropriate dimensional structure.\n",
	  progname,inName);
  if (!structure_check_obsmatrix(obsMatrix, chunk, L, &M, dz, dt,
				 &obsMatrixStride))
    Abort("%s: Dataset %s has an inappropriate dimensional structure.\n",
	  progname,obsMatrixName);
  if (!structure_check_transmatrix(transMatrix, chunk, M, dz, dt,
				   &transMatrixStride))
    Abort("%s: Dataset %s has an inappropriate dimensional structure.\n",
	  progname,transMatrixName);
  if (!structure_check_proccovmatrix(procCovMatrix, chunk, M))
    Abort("%s: Dataset %s has an inappropriate dimensional structure.\n",
	  progname,procCovMatrixName);
  if (!structure_check_obscovmatrix(obsCovMatrix, chunk, L))
    Abort("%s: Dataset %s has an inappropriate dimensional structure.\n",
	  progname,obsCovMatrixName);
  if (verbose_flg) {
    Message("# Observation vector length %d, process vector length %d\n",
	    L,M);
    Message("# Filtering %d samples per slice, %d slices, %d times\n",
	    samplesPerSlice,dz,dt);
  }

  /* Open output datasets. */
  if (outFlag) {
    output = create_output_ds( outName, input, chunk, M );
    hist_add_cl( output, argc, argv );
  }
  else output= NULL;
  if (errCovFlag) {
    errCovOut = create_errcov_ds( errCovOutName, input, chunk, M );
    hist_add_cl( errCovOut, argc, argv );
  }
  else errCovOut= NULL;
  if (logLikelihoodFlag) {
    logLikelihoodOut= create_loglike_ds( logLikelihoodName, input, chunk );
    hist_add_cl( logLikelihoodOut, argc, argv );
  }
  else logLikelihoodOut= NULL;

  /* Do the eigenvalue solution */
  apply_kalman( input, transMatrix, obsMatrix, 
		procCovMatrix, obsCovMatrix,
		output, errCovOut, logLikelihoodOut,
		chunk, L, M, dt, dz,
		samplesPerSlice, obsMatrixStride, transMatrixStride );

  /* Write and close data-sets */
  mri_close_dataset( input );
  mri_close_dataset( transMatrix );
  mri_close_dataset( obsMatrix );
  mri_close_dataset( procCovMatrix );
  mri_close_dataset( obsCovMatrix );
  if (outFlag) mri_close_dataset( output );
  if (errCovFlag) mri_close_dataset( errCovOut );
  if (logLikelihoodFlag) mri_close_dataset( logLikelihoodOut );
  
  if (verbose_flg) Message( "#      Kalman filtering complete.\n" );

  return 0;
}

