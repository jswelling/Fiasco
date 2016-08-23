/************************************************************
 *                                                          *
 *  kalmanfilter.c                                          *
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

#include <math.h>
#include <stdio.h>
#include <strings.h>
#include <assert.h>
#include <stdlib.h>
#include "mri.h"
#include "fmri.h"
#include "lapack.h"
#include "kalmanfilter.h"

static char rcsid[] = "$Id: kalmanfilter.c,v 1.5 2008/03/04 19:13:52 welling Exp $";

/**********************
 * Notes-
 * -We need a doc!
 *********************/

static void dumpMatrix( const char* label, 
			const double* m, int rows, int cols, 
			int dataOrder, FILE* ofile )
{
  int i,j;

  /* kalmanfilter uses classic BLAS routines, which access data in
   * column-major (Fortran-type) order.
   */
  switch (dataOrder) {
  case DATA_ORDER_COLUMNMAJOR:
    {
      fprintf(ofile,"%s (column major data order)\n",label);
      for (i=0; i<rows; i++) {
	for (j=0; j<cols; j++) fprintf(ofile,"%7.3g ",m[j*rows + i]);
	fprintf(ofile,"\n");
      }
    }
    break;
  case DATA_ORDER_ROWMAJOR:
    {
      fprintf(ofile,"%s (row major data order)\n",label);
      for (i=0; i<rows; i++) {
	for (j=0; j<cols; j++) fprintf(ofile,"%7.3g ",m[j + i*cols]);
	fprintf(ofile,"\n");
      }
    }
    break;
  default:
    {
      fprintf(stderr,"kalmanfilter:dumpMatrix: unknown data order code %d!\n",
	      dataOrder);
    }
  }
}

static void destroyKalmanProcess( KalmanProcess* target )
{
  if (target->A) free(target->A);
  if (target->H) free(target->H);
  if (target->Q) free(target->Q);
  if (target->R) free(target->R);
  if (target->scratch) free(target->scratch);
  if (target->iScratch) free(target->iScratch);
  free(target);
}

static void dumpKalmanProcess( const KalmanProcess* p, FILE* ofile )
{
  int L= p->L;
  int M= p->M;
  fprintf(ofile,"Kalman process: L= %d, M= %d; debugFlag= %d, data order %s\n",
	  L,M,p->debugFlag,
	  ((p->dataOrder==DATA_ORDER_ROWMAJOR)?
	   "DATA_ORDER_ROWMAJOR":"DATA_ORDER_COLUMNMAJOR"));
  if (p->validA) dumpMatrix("process matrix A",
			    p->A,M,M,p->dataOrder,ofile);
  else fprintf(ofile,"process matrix A uninitialized\n");
  if (p->validH) dumpMatrix("observation matrix H",
			    p->H,L,M,p->dataOrder,ofile);
  else fprintf(ofile,"observation matrix H uninitialized\n");
  if (p->validQ) dumpMatrix("process covariance Q",
			    p->Q,M,M,p->dataOrder,ofile);
  else fprintf(ofile,"process covariance Q uninitialized\n");
  if (p->validR) dumpMatrix("observation covariance R",
			    p->R,L,L,p->dataOrder,ofile);
  else fprintf(ofile,"observation covariance R uninitialized\n");
}

static void setA( KalmanProcess* self, const double* AIn, int sizeA )
{
  if (sizeA<self->M*self->M) {
    fprintf(stderr,"KalmanProcess:setA: input matrix is too small!\n");
    exit(-1);
  }
  bcopy(AIn, self->A, self->M*self->M*sizeof(double));
  self->validA= 1;
}
static void getA( const KalmanProcess* self, double* AOut, int sizeA ) {
  if (sizeA<self->M*self->M) {
    fprintf(stderr,"KalmanProcess:getA: output matrix is too small!\n");
    exit(-1);
  }
  bcopy(self->A, AOut, self->M*self->M*sizeof(double));
}

static void setH( KalmanProcess* self, const double* HIn, int sizeH )
{
  if (sizeH<self->L*self->M) {
    fprintf(stderr,"KalmanProcess:setH: input matrix is too small!\n");
    exit(-1);
  }
  bcopy(HIn, self->H, self->L*self->M*sizeof(double));
  self->validH= 1;
}
static void getH( const KalmanProcess* self, double* HOut, int sizeH ) {
  if (sizeH<self->L*self->M) {
    fprintf(stderr,"KalmanProcess:getH: output matrix is too small!\n");
    exit(-1);
  }
  bcopy(self->H, HOut, self->L*self->M*sizeof(double));
}

static void setQ( KalmanProcess* self, const double* QIn, int sizeQ )
{
  if (sizeQ<self->M*self->M) {
    fprintf(stderr,"KalmanProcess:setQ: input matrix is too small!\n");
    exit(-1);
  }
  bcopy(QIn, self->Q, self->M*self->M*sizeof(double));
  self->validQ= 1;
}
static void getQ( const KalmanProcess* self, double* QOut, int sizeQ ) {
  if (sizeQ<self->M*self->M) {
    fprintf(stderr,"KalmanProcess:getQ: output matrix is too small!\n");
    exit(-1);
  }
  bcopy(self->Q, QOut, self->M*self->M*sizeof(double));
}

static void setR( KalmanProcess* self, const double* RIn, int sizeR )
{
  if (sizeR<self->L*self->L) {
    fprintf(stderr,"KalmanProcess:setR: input matrix is too small!\n");
    exit(-1);
  }
  bcopy(RIn, self->R, self->L*self->L*sizeof(double));
  self->validR= 1;
}
static void getR( const KalmanProcess* self, double* ROut, int sizeR ) {
  if (sizeR<self->L*self->L) {
    fprintf(stderr,"KalmanProcess:getR: output matrix is too small!\n");
    exit(-1);
  }
  bcopy(self->R, ROut, self->L*self->L*sizeof(double));
}

/* returns increment to log likelihood */
static double kalman_iterate_once( const double* A, const double* H,
				   const double* Q, const double* R,
				   const double* z, double* xHat, double* P,
				   long L, long M, 
				   int missingFlag, int logLikelihoodFlag, 
				   int stepPandKFlag, int debugFlag,
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

  assert(scratchSize>=2*L+M+2*(M*M)+2*(L*M)+2*(L*L)+1);
  assert(iScratchSize>=L);

  /* xHat and P are at time t-1; we mean to step them to time t */
  /* xHatPrior in scratch, size M; PPrior in scratch, size M*M */

  if (debugFlag) {
    dumpMatrix("### A: ###",A,M,M,DATA_ORDER_COLUMNMAJOR,stderr);
    dumpMatrix("### H: ###",H,L,M,DATA_ORDER_COLUMNMAJOR,stderr);
    dumpMatrix("### Q: ###",Q,M,M,DATA_ORDER_COLUMNMAJOR,stderr);
    dumpMatrix("### R: ###",R,L,L,DATA_ORDER_COLUMNMAJOR,stderr);
    dumpMatrix("### z: ###",z,L,1,DATA_ORDER_COLUMNMAJOR,stderr);
    dumpMatrix("### xHat^T: ###",xHat,1,M,DATA_ORDER_COLUMNMAJOR,stderr);
    dumpMatrix("### P :###",P,M,M,DATA_ORDER_COLUMNMAJOR,stderr);
  }
  if (stepPandKFlag) {
    /* PPrior= A*P*transpose(A) + Q;  MxM * MxM * MxM */
    DGEMM("n","n",&iM,&iM,&iM,&one,(double*)A,&iM,P,&iM,&zero,tmpMM,&iM);
    bcopy(Q,PPrior,M*M*sizeof(double));
    DGEMM("n","t",&iM,&iM,&iM,&one,tmpMM,&iM,(double*)A,&iM,
	  &one,PPrior,&iM);
    if (debugFlag) 
      dumpMatrix("### PPrior: ###",PPrior,M,M,DATA_ORDER_COLUMNMAJOR,stderr);

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
      if (debugFlag) 
	dumpMatrix("### InnovCov: ###",tmpLL,L,L,DATA_ORDER_COLUMNMAJOR,
		   stderr);
      /* Compute inverse(tmpLL) in place */
      DGETRF( &iL, &iL, tmpLL, &iL, iScratch, &info );
      if (info!=0) {
	fprintf(stderr,
		"KalmanProcess::apply: matrix inversion failed computing K; err %d on DGETRF!\n",
		info);
	exit(-1);
      }
#ifdef never
      if (debugFlag) dumpMatrix("### LU decomp of InnovCov: ###",
				tmpLL,L,L,DATA_ORDER_COLUMNMAJOR,stderr);
#endif
      /* Snag the determinant while it's handy */
      *determinantOfInnovCov= 1.0;
      for (i=0; i<L*L; i+=(L+1))
	*determinantOfInnovCov *= tmpLL[i];
      DGETRI( &iL, tmpLL, &iL, iScratch, tmpLL2, &iLL, &info );
      if (info!=0) {
	fprintf(stderr,
		"KalmanProcess::apply: matrix inversion failed computing K; err %d on DGETRI!\n",
		info);
	exit(-1);
      }
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
    if (debugFlag) 
      dumpMatrix("### K: ###",K,M,L,DATA_ORDER_COLUMNMAJOR,stderr);

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
  if (debugFlag) 
    dumpMatrix("### xHatPrior^T ###",xHatPrior,1,M,DATA_ORDER_COLUMNMAJOR,
	       stderr);
  
  /* xHat= xHatPrior + K*(z - H*xHatPrior); */
  /* (z - H*xHatPrior) is the 'innovation' */
  bcopy(z, innovation, L*sizeof(double));
  DGEMV("n", &iL, &iM, &neg_one, (double*)H, &iL, xHatPrior, &int_one,
	&one, innovation, &int_one);
  if (debugFlag)
    dumpMatrix("### innovation ###",innovation,L,1,DATA_ORDER_COLUMNMAJOR,
	       stderr);
  bcopy(xHatPrior, xHat, M*sizeof(double));
  DGEMV("n", &iM, &iL, &one, K, &iM, innovation, &int_one,
	&one, xHat, &int_one);

  if (logLikelihoodFlag && !missingFlag) {
    /* I think the xHat's in this expression are really xHatPrior's! */
    /* log(f)= -(1/2) log[ (2pi)^L * det(H*PPrior*H^T + R) ]
     *
     *        -(1/2)(z-H*xHat)^T * (H*PPrior*H^T + R)^-1 * (z-H*xHat)
     */
    double dot= 0.0;
    DGEMV("n", &iL, &iL, &one, tmpLL, &iL, innovation, &int_one, 
	   &zero, tmpL, &int_one );
    for (i=0; i<L; i++) dot += innovation[i]*tmpL[i];
    logLikelihoodIncrement= 
      -0.5*(L*log(2.0*M_PI) + log(*determinantOfInnovCov) + dot);
    if (debugFlag) fprintf(stderr,"logLikelihoodIncrement= %g\n",
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
/* row major (C-style) data order */
static double kalman_iterate_once_rwmjr( const double* A, const double* H,
					 const double* Q, const double* R,
					 const double* z, double* xHat, 
					 double* P,
					 long L, long M, 
					 int missingFlag, 
					 int logLikelihoodFlag, 
					 int stepPandKFlag, 
					 int debugFlag,
					 double* scratch, 
					 long scratchSize,
					 int* iScratch, 
					 long iScratchSize )
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

  assert(scratchSize>=2*L+M+2*(M*M)+2*(L*M)+2*(L*L)+1);
  assert(iScratchSize>=L);

  /* xHat and P are at time t-1; we mean to step them to time t */
  /* xHatPrior in scratch, size M; PPrior in scratch, size M*M */

  if (debugFlag) {
    dumpMatrix("### A: ###",A,M,M,DATA_ORDER_ROWMAJOR,stderr);
    dumpMatrix("### H: ###",H,L,M,DATA_ORDER_ROWMAJOR,stderr);
    dumpMatrix("### Q: ###",Q,M,M,DATA_ORDER_ROWMAJOR,stderr);
    dumpMatrix("### R: ###",R,L,L,DATA_ORDER_ROWMAJOR,stderr);
    dumpMatrix("### z: ###",z,L,1,DATA_ORDER_ROWMAJOR,stderr);
    dumpMatrix("### xHat^T: ###",xHat,1,M,DATA_ORDER_ROWMAJOR,stderr);
    dumpMatrix("### P :###",P,M,M,DATA_ORDER_ROWMAJOR,stderr);
  }
  if (stepPandKFlag) {
    /* PPrior= A*P*transpose(A) + Q;  MxM * MxM * MxM + MxM*/
    /* The things stored in the arrays are the transposes (by virtue
     * of the data order).  We want all intermediate results as
     * transposes as well, for compatibility.  To request a transpose,
     * we use "n" because the data is *already* the transpose.
     */
    /* Transpose: PPrior'= A*P'*A' + Q';  MxM * MxM * MxM + MxM*/
    DGEMM("t","n",&iM,&iM,&iM,&one,(double*)A,&iM,P,&iM,&zero,tmpMM,&iM);
    bcopy(Q,PPrior,M*M*sizeof(double));
    DGEMM("n","n",&iM,&iM,&iM,&one,tmpMM,&iM,(double*)A,&iM,
	  &one,PPrior,&iM);
    if (debugFlag) 
      dumpMatrix("### PPrior: ###",PPrior,M,M,DATA_ORDER_ROWMAJOR,stderr);

    if (missingFlag) {
      for (i=0; i<L*M; i++) K[i]= 0.0;
      *determinantOfInnovCov= 1.0;
    }
    else {
      /* K= PPrior*transpose(H)*inverse(H*PPrior*transpose(H) + R); */
      /* The thing that gets inverted is the innovation covariance */
      /* Start by computing tmpLL= H*PPrior*transpose(H) + R */
      /* Transpose: tmpLL'= H*PPrior'*H' + R' */
      DGEMM("t","n",&iL,&iM,&iM,&one,(double*)H,&iM,PPrior,&iM,
	    &zero,tmpLM,&iL);
      bcopy(R,tmpLL,L*L*sizeof(double));
      DGEMM("n","n",&iL,&iL,&iM,&one,tmpLM,&iL,(double*)H,&iM,
	    &one,tmpLL,&iL);
      if (debugFlag) 
	dumpMatrix("### InnovCov: ###",tmpLL,L,L,DATA_ORDER_ROWMAJOR,
		   stderr);
      /* Compute inverse(tmpLL) in place */
      DGETRF( &iL, &iL, tmpLL, &iL, iScratch, &info );
      if (info!=0) {
	fprintf(stderr,
		"KalmanProcess::apply: matrix inversion failed computing K; err %d on DGETRF!\n",
		info);
	exit(-1);
      }
#ifdef never
      if (debugFlag) dumpMatrix("### LU decomp of InnovCov: ###",
				tmpLL,L,L,DATA_ORDER_ROWMAJOR,stderr);
#endif
      /* Snag the determinant while it's handy */
      *determinantOfInnovCov= 1.0;
      for (i=0; i<L*L; i+=(L+1))
	*determinantOfInnovCov *= tmpLL[i];
      DGETRI( &iL, tmpLL, &iL, iScratch, tmpLL2, &iLL, &info );
      if (info!=0) {
	fprintf(stderr,
		"KalmanProcess::apply: matrix inversion failed computing K; err %d on DGETRI!\n",
		info);
	exit(-1);
      }
      /* Compute K= PPrior*transpose(H)*tmpLL */
      /* Transpose: K'= tmpLL'*H*PPrior', LxL * LxM * MxM*/
      DGEMM("n","t",&iL, &iM, &iL, &one, tmpLL, &iL, (double*)H, &iM,
	    &zero, tmpLM, &iL);
      DGEMM("n","n",&iL, &iM, &iM, &one, tmpLM, &iL, (double*)PPrior, &iM,
	    &zero,K,&iL);
      /* !!!NOTE!!! we are keeping tmpLL around for possible use in
       * computing the log likelihood increment.  It contains the
       * inverse innovation covariance.
       */
    }
    if (debugFlag) 
      dumpMatrix("### K: ###",K,M,L,DATA_ORDER_ROWMAJOR,stderr);

    /* P= (1-K*H)*PPrior; */
    /* or more conveniently, P= -1*K*H*PPrior + 1*PPrior; */
    /* transpose: P'= -1*PPrior'*H'*K' + 1*PPrior', MxM * MxL * LxM */
    DGEMM("n","n",&iM, &iM, &iL, &one, (double*)H, &iM, (double*)K, &iL,
	  &zero, tmpMM, &iM);
    bcopy(PPrior, P, M*M*sizeof(double));
    DGEMM("n","n",&iM, &iM, &iM, &neg_one, PPrior, &iM, tmpMM, &iM,
	  &one, P, &iM);
  }

  /* xHatPrior= A*xHat;  MxM * M  (xHat has not yet been updated) */
  /* transpose: xHatPrior'= (A*xHat)', MxM * M */
  DGEMV("t", &iM, &iM, &one, (double*)A, &iM, xHat, &int_one, &zero, 
	xHatPrior, &int_one);
  if (debugFlag) 
    dumpMatrix("### xHatPrior^T ###",xHatPrior,1,M,DATA_ORDER_ROWMAJOR,
	       stderr);
  
  /* xHat= xHatPrior + K*(z - H*xHatPrior); */
  /* (z - H*xHatPrior) is the 'innovation' */
  /* transpose: xHat'= xHatPrior' + (z - H*xHatPrior)'*K' */
  /*                 = xHatPrior' + (z' - (H*xHatPrior)')*K' (LxM * M)'*LxM */
  /*                 = xHatPrior' + (K*(z-(H*xHatPrior))' MxL * (LxM * M) */
  bcopy(z, innovation, L*sizeof(double));
  DGEMV("t", &iM, &iL, &neg_one, (double*)H, &iM, xHatPrior, &int_one,
	&one, innovation, &int_one);
  if (debugFlag)
    dumpMatrix("### innovation ###",innovation,L,1,DATA_ORDER_ROWMAJOR,
	       stderr);
  bcopy(xHatPrior, xHat, M*sizeof(double));
  DGEMV("t", &iL, &iM, &one, K, &iL, innovation, &int_one,
	&one, xHat, &int_one);

  if (logLikelihoodFlag && !missingFlag) {
    /* I think the xHat's in this expression are really xHatPrior's! */
    /* log(f)= -(1/2) log[ (2pi)^L * det(H*PPrior*H^T + R) ]
     *
     *        -(1/2)(z-H*xHat)^T * (H*PPrior*H^T + R)^-1 * (z-H*xHat)
     */
    double dot= 0.0;
    DGEMV("n", &iL, &iL, &one, tmpLL, &iL, innovation, &int_one, 
	   &zero, tmpL, &int_one );
    for (i=0; i<L; i++) dot += innovation[i]*tmpL[i];
    logLikelihoodIncrement= 
      -0.5*(L*log(2.0*M_PI) + log(*determinantOfInnovCov) + dot);
    if (debugFlag) fprintf(stderr,"logLikelihoodIncrement= %g\n",
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
					double* P, long M, 
					int missingFlag, int logLikelihoodFlag, 
					int stepPandKFlag, int debugFlag,
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

  assert(scratchSize>=2+M+2*(M*M)+2*M+2+1);
  assert(iScratchSize>=1);

  /* xHat and P are at time t-1; we mean to step them to time t */
  /* xHatPrior in scratch, size M; PPrior in scratch, size M*M */

  if (debugFlag) {
    dumpMatrix("### A: ###",A,M,M,DATA_ORDER_COLUMNMAJOR,stderr);
    dumpMatrix("### H: ###",H,1,M,DATA_ORDER_COLUMNMAJOR,stderr);
    dumpMatrix("### Q: ###",Q,M,M,DATA_ORDER_COLUMNMAJOR,stderr);
    fprintf(stderr,"### R= %g\n",R[0]);
    fprintf(stderr,"### z= %g\n",z[0]);
    dumpMatrix("### xHat^T: ###",xHat,1,M,DATA_ORDER_COLUMNMAJOR,stderr);
    dumpMatrix("### P :###",P,M,M,DATA_ORDER_COLUMNMAJOR,stderr);
  }
  if (stepPandKFlag) {
    /* PPrior= A*P*transpose(A) + Q;  MxM * MxM * MxM */
    DGEMM("n","n",&iM,&iM,&iM,&one,(double*)A,&iM,P,&iM,&zero,tmpMM,&iM);
    bcopy(Q,PPrior,M*M*sizeof(double));
    DGEMM("n","t",&iM,&iM,&iM,&one,tmpMM,&iM,(double*)A,&iM,
	  &one,PPrior,&iM);
    if (debugFlag) 
      dumpMatrix("### PPrior: ###",PPrior,M,M,DATA_ORDER_COLUMNMAJOR,stderr);

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
    if (debugFlag) 
      dumpMatrix("### K: ###",K,M,1,DATA_ORDER_COLUMNMAJOR,stderr);

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
  if (debugFlag) 
    dumpMatrix("### xHatPrior^T ###",xHatPrior,1,M,DATA_ORDER_COLUMNMAJOR,
	       stderr);
  
  /* xHat= xHatPrior + K*(z - H*xHatPrior); */
  /* (z - H*xHatPrior) is the 'innovation' */
  innovation[0]= z[0];
  for (i=0; i<M; i++) innovation[0] -= H[i]*xHatPrior[i];
  for (i=0; i<M; i++) xHat[i]= xHatPrior[i] + K[i]*innovation[0];

  if (logLikelihoodFlag && !missingFlag) {
    /* I think the xHat's in this expression are really xHatPrior's! */
    /* log(f)= -(1/2) log[ (2pi)^L * det(H*PPrior*H^T + R) ]
     *
     *        -(1/2)(z-H*xHat)^T * (H*PPrior*H^T + R)^-1 * (z-H*xHat)
     */
    double dot= innovation[0]*tmpLL[0]*innovation[0];
    logLikelihoodIncrement= 
      -0.5*(log(2.0*M_PI) + log(*determinantOfInnovCov) + dot);
    if (debugFlag) fprintf(stderr,"logLikelihoodIncrement= %g\n",
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
/* Special case for the common L==1, row major (C-style) data order */
static double kalman_iterate_once_Leq1_rwmjr( const double* A, const double* H,
					      const double* Q, const double* R,
					      const double* z, double* xHat, 
					      double* P, long M, 
					      int missingFlag, 
					      int logLikelihoodFlag, 
					      int stepPandKFlag, 
					      int debugFlag,
					      double* scratch, 
					      long scratchSize,
					      int* iScratch, 
					      long iScratchSize )
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

  assert(scratchSize>=2+M+2*(M*M)+2*M+2+1);
  assert(iScratchSize>=1);

  /* xHat and P are at time t-1; we mean to step them to time t */
  /* xHatPrior in scratch, size M; PPrior in scratch, size M*M */

  if (debugFlag) {
    dumpMatrix("### A: ###",A,M,M,DATA_ORDER_ROWMAJOR,stderr);
    dumpMatrix("### H: ###",H,1,M,DATA_ORDER_ROWMAJOR,stderr);
    dumpMatrix("### Q: ###",Q,M,M,DATA_ORDER_ROWMAJOR,stderr);
    fprintf(stderr,"### R= %g\n",R[0]);
    fprintf(stderr,"### z= %g\n",z[0]);
    dumpMatrix("### xHat^T: ###",xHat,1,M,DATA_ORDER_ROWMAJOR,stderr);
    dumpMatrix("### P :###",P,M,M,DATA_ORDER_ROWMAJOR,stderr);
  }
  if (stepPandKFlag) {
    /* PPrior= A*P*transpose(A) + Q;  MxM * MxM * MxM */
    /* The things stored in the arrays are the transposes (by virtue
     * of the data order).  We want all intermediate results as
     * transposes as well, for compatibility.  To request a transpose,
     * we use "n" because the data is *already* the transpose.
     */
    /* Transpose: PPrior'= A*P'*A' + Q';  MxM * MxM * MxM + MxM*/
    DGEMM("t","n",&iM,&iM,&iM,&one,(double*)A,&iM,P,&iM,&zero,tmpMM,&iM);
    bcopy(Q,PPrior,M*M*sizeof(double));
    DGEMM("n","n",&iM,&iM,&iM,&one,tmpMM,&iM,(double*)A,&iM,
	  &one,PPrior,&iM);
    if (debugFlag) 
      dumpMatrix("### PPrior: ###",PPrior,M,M,DATA_ORDER_ROWMAJOR,stderr);

    if (missingFlag) {
      for (i=0; i<M; i++) K[i]= 0.0;
      *determinantOfInnovCov= 1.0;
    }
    else {
      /* K= PPrior*transpose(H)*inverse(H*PPrior*transpose(H) + R); */
      /* The thing that gets inverted is the innovation covariance */
      /* Start by computing tmpLL= H*PPrior*transpose(H) + R */
      /* Transpose: tmpLL'= H*PPrior'*H' + R' */
      DGEMV("n", &iM, &iM, &one, PPrior, &iM, (double*)H, &int_one, &zero,
	    tmpLM, &int_one);
      tmpLL[0]= R[0];
      for (i=0; i<M; i++) tmpLL[0] += H[i]*tmpLM[i];
      *determinantOfInnovCov= tmpLL[0];
      /* Compute inverse(tmpLL) in place */
      tmpLL[0]= 1.0/tmpLL[0];
      /* Compute K= PPrior*transpose(H)*tmpLL; tmpLL is a scalar */
      /* transpose: K= (tmpLL)*H*PPrior'; 1x1 * 1xM * MxM */
      /* K is Mx1, so no need to deal with transpose */
      DGEMV("t", &iM, &iM, tmpLL, PPrior, &iM, (double*)H, &int_one, &zero, K, 
	    &int_one);
      /* !!!NOTE!!! we are keeping tmpLL around for possible use in
       * computing the log likelihood increment.  It contains the
       * inverse innovation covariance.
       */
    }
    if (debugFlag) 
      dumpMatrix("### K: ###",K,M,1,DATA_ORDER_ROWMAJOR,stderr);

    /* P= (1-K*H)*PPrior; */
    /* or more conveniently, P= -1*K*H*PPrior + 1*PPrior; */
    /* transpose: P'= -1*PPrior'*H'*K' + 1*PPrior';  MxM * MxL * LxM */
    DGEMV("n", &iM, &iM, &neg_one, PPrior, &iM, (double*)H, &int_one,
	  &zero, tmpLM, &int_one);
    bcopy(PPrior, P, M*M*sizeof(double));
    DGEMM("n","n",&iM, &iM, &int_one, &one, tmpLM, &iM, K, &int_one, 
	  &one, P, &iM);
  }

  /* xHatPrior= A*xHat;  MxM * M  (xHat has not yet been updated) */
  /* xHatPrior is Mx1, so we don't have to worry about transpose */
  DGEMV("t", &iM, &iM, &one, (double*)A, &iM, xHat, &int_one, &zero, 
	xHatPrior, &int_one);
  if (debugFlag) 
    dumpMatrix("### xHatPrior^T ###",xHatPrior,1,M,DATA_ORDER_ROWMAJOR,
	       stderr);
  
  /* xHat= xHatPrior + K*(z - H*xHatPrior); */
  /* (z - H*xHatPrior) is the 'innovation' */
  /* All of these are vectors, so we don't have to worry about transpose */
  innovation[0]= z[0];
  for (i=0; i<M; i++) innovation[0] -= H[i]*xHatPrior[i];
  for (i=0; i<M; i++) xHat[i]= xHatPrior[i] + K[i]*innovation[0];

  if (logLikelihoodFlag && !missingFlag) {
    /* I think the xHat's in this expression are really xHatPrior's! */
    /* log(f)= -(1/2) log[ (2pi)^L * det(H*PPrior*H^T + R) ]
     *
     *        -(1/2)(z-H*xHat)^T * (H*PPrior*H^T + R)^-1 * (z-H*xHat)
     */
    /* All scalars */
    double dot= innovation[0]*tmpLL[0]*innovation[0];
    logLikelihoodIncrement= 
      -0.5*(log(2.0*M_PI) + log(*determinantOfInnovCov) + dot);
    if (debugFlag) fprintf(stderr,"logLikelihoodIncrement= %g\n",
			   logLikelihoodIncrement);
  }
  else {
    /* save some pointless computation */
    logLikelihoodIncrement= 0.0;
  }
  /* And we're done. */
  return logLikelihoodIncrement;
}

double apply( const KalmanProcess* self, KalmanState* state,
	      const double* z, int sizeZ, int missing, 
	      int calcLogLikelihoodDelta, int updatePandK )
{
  if (!self->validA || !self->validH || !self->validQ || !self->validR) {
    fprintf(stderr,"KalmanProcess:apply: A, H, Q or R not set!\n");
    exit(-1);
  }

  if (self->M != state->M) {
    fprintf(stderr,
	    "KalmanProcess:apply: my M does not match the state provided!\n");
    exit(-1);
  }

  if (sizeZ<self->L) {
    fprintf(stderr,
	    "KalmanProcess:apply: not enough values in observation!\n");
    exit(-1);
  }

  if (self->dataOrder!=state->dataOrder) {
    fprintf(stderr,
	    "KalmanProcess:apply: my data order does not match state's!\n");
    exit(-1);
  }

  if (self->dataOrder==DATA_ORDER_COLUMNMAJOR)
    return kalman_iterate_once( self->A, self->H, self->Q, self->R,
				z, state->x, state->P, self->L, self->M, 
				missing, calcLogLikelihoodDelta, 
				updatePandK, self->debugFlag,
				self->scratch, self->scratchSize, 
				self->iScratch, self->iScratchSize );
  else if (self->dataOrder==DATA_ORDER_ROWMAJOR)
    return kalman_iterate_once_rwmjr( self->A, self->H, self->Q, self->R,
				      z, state->x, state->P, self->L, self->M, 
				      missing, calcLogLikelihoodDelta, 
				      updatePandK, self->debugFlag,
				      self->scratch, self->scratchSize, 
				      self->iScratch, self->iScratchSize );
  else
    fprintf(stderr,"KalmanProcess:apply: unknown data order code %d!\n",
	    self->dataOrder);
    exit(-1);

}

static double apply_Leq1( const KalmanProcess* self, KalmanState* state,
			  const double* z, int sizeZ, int missing, 
			  int calcLogLikelihoodDelta, int updatePandK )
{
  assert(self->L==1);
  if (!self->validA || !self->validH || !self->validQ || !self->validR) {
    fprintf(stderr,"KalmanProcess:apply_Leq1: A, H, Q or R not set!\n");
    exit(-1);
  }

  if (!state->x || !state->P) {
    fprintf(stderr,"KalmanProcess:apply: x or P in state not set!\n");
    exit(-1);
  }

  if (self->M != state->M) {
    fprintf(stderr,
	    "KalmanProcess:apply: my M does not match the state provided!\n");
    exit(-1);
  }

  if (sizeZ<self->L) {
    fprintf(stderr,
	    "KalmanProcess:apply: not enough values in observation!\n");
    exit(-1);
  }

  if (self->dataOrder!=state->dataOrder) {
    fprintf(stderr,
	    "KalmanProcess:apply: my data order does not match state's!\n");
    exit(-1);
  }

  if (self->dataOrder==DATA_ORDER_COLUMNMAJOR) 
    return kalman_iterate_once_Leq1( self->A, self->H, self->Q, self->R,
				     z, state->x, state->P, self->M, 
				     missing, calcLogLikelihoodDelta, 
				     updatePandK, self->debugFlag,
				     self->scratch, self->scratchSize, 
				     self->iScratch, self->iScratchSize );
  else if (self->dataOrder==DATA_ORDER_ROWMAJOR) 
    return kalman_iterate_once_Leq1_rwmjr(self->A, self->H, self->Q, self->R,
					  z, state->x, state->P, self->M, 
					  missing, calcLogLikelihoodDelta, 
					  updatePandK, self->debugFlag,
					  self->scratch, self->scratchSize, 
					  self->iScratch, self->iScratchSize);
  else {
    fprintf(stderr,"KalmanProcess:apply_Leq1: unknown data order code %d!\n",
	    self->dataOrder);
    exit(-1);
  }
}

static void processSetDebug( KalmanProcess* self, int val )
{
  self->debugFlag= val;
}

static int processGetDebug( KalmanProcess* self )
{
  return self->debugFlag;
}

static int processGetDataOrder( KalmanProcess* self )
{
  return self->dataOrder;
}

KalmanProcess* klmn_createKalmanProcess(int L, int M)
{
  KalmanProcess* result= (KalmanProcess*)malloc(sizeof(KalmanProcess));

  if (!result) MALLOC_FAILURE(sizeof(KalmanProcess),byte);
  result->L= L;
  result->M= M;
  result->debugFlag= 0;
  result->dataOrder= DATA_ORDER_COLUMNMAJOR;

  if (!(result->A=(double*)malloc(M*M*sizeof(double))))
    MALLOC_FAILURE(M*M,double);
  if (!(result->H=(double*)malloc(L*M*sizeof(double))))
    MALLOC_FAILURE(L*M,double);
  if (!(result->Q=(double*)malloc(M*M*sizeof(double))))
    MALLOC_FAILURE(M*M,double);
  if (!(result->R=(double*)malloc(L*L*sizeof(double))))
    MALLOC_FAILURE(L*L,double);
  result->validA= result->validH= result->validQ= result->validR= 0;

  result->scratchSize= 2*L+M+2*(M*M)+2*(L*M)+2*(L*L)+1;
  if (!(result->scratch=(double*)malloc(result->scratchSize*sizeof(double))))
    MALLOC_FAILURE(result->scratchSize,double);
  result->iScratchSize= L;
  if (!(result->iScratch=(int*)malloc(result->iScratchSize*sizeof(int))))
    MALLOC_FAILURE(result->iScratchSize,int);

  result->destroySelf= destroyKalmanProcess;
  result->dumpSelf= dumpKalmanProcess;
  result->setDebug= processSetDebug;
  result->getDebug= processGetDebug;
  result->getDataOrder= processGetDataOrder;
  result->setA= setA;
  result->getA= getA;
  result->setH= setH;
  result->getH= getH;
  result->setQ= setQ;
  result->getQ= getQ;
  result->setR= setR;
  result->getR= getR;
  if (result->L==1) result->apply= apply_Leq1;
  else result->apply= apply;

  return result;
}

KalmanProcess* klmn_createRowMajorKalmanProcess(int L, int M)
{
  KalmanProcess* result= (KalmanProcess*)malloc(sizeof(KalmanProcess));

  if (!result) MALLOC_FAILURE(sizeof(KalmanProcess),byte);
  result->L= L;
  result->M= M;
  result->debugFlag= 0;
  result->dataOrder= DATA_ORDER_ROWMAJOR;

  if (!(result->A=(double*)malloc(M*M*sizeof(double))))
    MALLOC_FAILURE(M*M,double);
  if (!(result->H=(double*)malloc(L*M*sizeof(double))))
    MALLOC_FAILURE(L*M,double);
  if (!(result->Q=(double*)malloc(M*M*sizeof(double))))
    MALLOC_FAILURE(M*M,double);
  if (!(result->R=(double*)malloc(L*L*sizeof(double))))
    MALLOC_FAILURE(L*L,double);
  result->validA= result->validH= result->validQ= result->validR= 0;

  result->scratchSize= 2*L+M+2*(M*M)+2*(L*M)+2*(L*L)+1;
  if (!(result->scratch=(double*)malloc(result->scratchSize*sizeof(double))))
    MALLOC_FAILURE(result->scratchSize,double);
  result->iScratchSize= L;
  if (!(result->iScratch=(int*)malloc(result->iScratchSize*sizeof(int))))
    MALLOC_FAILURE(result->iScratchSize,int);

  result->destroySelf= destroyKalmanProcess;
  result->dumpSelf= dumpKalmanProcess;
  result->setDebug= processSetDebug;
  result->getDebug= processGetDebug;
  result->getDataOrder= processGetDataOrder;
  result->setA= setA;
  result->getA= getA;
  result->setH= setH;
  result->getH= getH;
  result->setQ= setQ;
  result->getQ= getQ;
  result->setR= setR;
  result->getR= getR;
  if (result->L==1) result->apply= apply_Leq1;
  else result->apply= apply;

  return result;
}


static void destroyKalmanState( KalmanState* target )
{
  /* The State doesn't own its memory, so this is really simple */
  free(target);
}

static void dumpKalmanState( const KalmanState* s, FILE* ofile )
{
  fprintf(ofile,"Kalman state: M= %d, data order %s\n",
	  s->M,
	  (s->dataOrder==DATA_ORDER_ROWMAJOR)?
	   "DATA_ORDER_ROWMAJOR":"DATA_ORDER_COLUMNMAJOR");
  if (s->x != NULL)
    dumpMatrix("Hidden state x^T",
	       s->x,1,s->M,s->dataOrder,ofile);
  else fprintf(ofile,"Hidden state x not initialized\n");
  if (s->P != NULL)
    dumpMatrix("Error covariance P",s->P,s->M,s->M,s->dataOrder,ofile);
  else fprintf(ofile,"Error covariance P not initialized\n");
}

static int stateGetDataOrder( KalmanState* self )
{
  return self->dataOrder;
}

static void setX( KalmanState* self, double* x, int sizeX )
{
  if (sizeX<self->M) {
    fprintf(stderr,"KalmanState:setX: state vector is too small!\n");
    exit(-1);
  }
  self->x= x;
  
}

static void setP( KalmanState* self, double* P, int sizeP )
{
  if (sizeP<self->M*self->M) {
    fprintf(stderr,"KalmanState:setP: error covariance array is too small!\n");
    exit(-1);
  }
  self->P= P;
  
}

KalmanState* klmn_createKalmanState(int M)
{
  KalmanState* result= (KalmanState*)malloc(sizeof(KalmanState));

  if (!result) MALLOC_FAILURE(sizeof(KalmanState),byte);
  result->M= M;
  result->dataOrder= DATA_ORDER_COLUMNMAJOR;
  result->x= NULL;
  result->P= NULL;

  result->destroySelf= destroyKalmanState;
  result->dumpSelf= dumpKalmanState;
  result->getDataOrder= stateGetDataOrder;
  result->setX= setX;
  result->setP= setP;

  return result;
}

KalmanState* klmn_createRowMajorKalmanState(int M)
{
  KalmanState* result= (KalmanState*)malloc(sizeof(KalmanState));

  if (!result) MALLOC_FAILURE(sizeof(KalmanState),byte);
  result->M= M;
  result->dataOrder= DATA_ORDER_ROWMAJOR;
  result->x= NULL;
  result->P= NULL;

  result->destroySelf= destroyKalmanState;
  result->dumpSelf= dumpKalmanState;
  result->getDataOrder= stateGetDataOrder;
  result->setX= setX;
  result->setP= setP;

  return result;
}

static void destroyKalmanFilter( KalmanFilter* target )
{
  target->state->destroySelf(target->state);
  free(target);
}

static void dumpKalmanFilter( const KalmanFilter* f, FILE* ofile )
{
  fprintf(ofile,"###########################\n");
  fprintf(ofile,"Kalman filter: L= %d, M= %d, consisting of:\n",
	  f->L,f->M);
  f->process->dumpSelf(f->process,ofile);
  f->state->dumpSelf(f->state,ofile);
  fprintf(ofile,"###########################\n");
}

static int filterGetDataOrder( KalmanFilter* self )
{
  /* Our data order is the process' data order */
  return self->process->getDataOrder(self->process);
}

static double runKalmanFilter( KalmanFilter* self, long tdim,
			       const double* zIn, long sizeZ,
			       double* xOut, long sizeXOut,
			       double* POut, long sizePOut,
			       int calcLogLikelihood )
{
  long t= 0;
  const double* zFrame= NULL;
  double* xFrame= NULL;
  double* errCovFrame= NULL;
  double* PBuffer= NULL;
  double* xBuffer= NULL;
  double logLikelihood= 0.0;
  long i= 0;

  /* Verify that various buffers provided by the caller are large enough */
  if (sizeZ<tdim*self->L) {
    fprintf(stderr,
	    "KalmanFilter:run: not enough input data provided! (%ld vs %ld)\n",
	    sizeZ,tdim*self->L);
    exit(-1);
  }
  if (sizeXOut<tdim*self->M) {
    fprintf(stderr,
	    "KalmanFilter:run: not enough output data space provided! (%ld vs %ld)\n",
	    sizeXOut,tdim*self->M);
    exit(-1);
  }
  if (POut!=NULL && sizePOut<tdim*self->M*self->M) {
    fprintf(stderr,
	    "KalmanFilter:run: output space for P is too small! (%ld vs %ld)\n",
	    sizePOut,tdim*self->M*self->M);
    exit(-1);
  }
    
  zFrame= zIn;
  xFrame= xOut;
  errCovFrame= POut;

  if (!(xBuffer=(double*)malloc(self->M*sizeof(double))))
    MALLOC_FAILURE(self->M,double);
  if (!(PBuffer=(double*)malloc(self->M*self->M*sizeof(double))))
    MALLOC_FAILURE(self->M*self->M,double);

  /* x gets initialized to all 1's. */
  for (i=0; i<self->M; i++) xBuffer[i]= 1.0;

  /* errCovFrame (aka P) gets initialized to procCovMatrix */
  bcopy(self->process->Q,PBuffer,self->M*self->M*sizeof(double));

  /* logLikelihoodSample gets initialized to 0.0 */
  logLikelihood= 0.0;

  self->state->setX(self->state, xBuffer, self->M);
  self->state->setP(self->state, PBuffer, self->M*self->M);
  for (t=0; t<tdim; t++) {

    /* Solve one iteration, and propogate the results forward */
    if (self->process->debugFlag)
      fprintf(stderr,"-------------\n t=%ld\n-------------\n",t);
    
    logLikelihood += self->process->apply(self->process, self->state,
					  zFrame, self->L, 
					  0, calcLogLikelihood, 1);
    
    zFrame += self->L;

    bcopy(xBuffer,xFrame,self->M*sizeof(double));
    xFrame += self->M;

    if (POut!=NULL) {
      bcopy(PBuffer,errCovFrame,self->M*self->M*sizeof(double));
      errCovFrame+=self->M*self->M;
    }
  }

  free(xBuffer);
  free(PBuffer);
  return logLikelihood;
}

KalmanFilter* klmn_createKalmanFilter(KalmanProcess* process)
{
  KalmanFilter* result= (KalmanFilter*)malloc(sizeof(KalmanFilter));

  if (!result) MALLOC_FAILURE(sizeof(KalmanFilter),byte);

  result->L= process->L;
  result->M= process->M;
  result->process= process;
  if (process->getDataOrder(process)==DATA_ORDER_COLUMNMAJOR)
    result->state= klmn_createKalmanState(result->M);
  else
    result->state= klmn_createRowMajorKalmanState(result->M);
  result->dumpSelf= dumpKalmanFilter;
  result->getDataOrder= filterGetDataOrder;
  result->run= runKalmanFilter;

  return result;
}
