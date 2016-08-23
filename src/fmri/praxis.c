/************************************************************
 *                                                          *
 *  praxis.c                                                *
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
 ************************************************************/
/* This is the routine PRAXIS from netlib, by Richard Brent. 
 * version from Stanford Linear Accelerator Center, dated 3/1/73.
 * It was manually transcoded from Fortran to C by Joel Welling during 3/00.
 * This C version should theoretically be reentrant.
 */

#include <stdio.h>
#include <math.h>
#include <fmri.h>
#include <stdlib.h>
#include "lapack.h"

static char rcsid[] = "$Id: praxis.c,v 1.6 2004/03/01 05:41:42 welling Exp $";

/* Maximum dimensionality of problem */
#define WORK_DIM 20

/* Praxis convergence paramters.  See comments in praxis body. */
#define SCBD 1.0
#define ILLC 0
#define KTM 1

/* prx_minfit convergence parameters. */
#define MINFIT_KTM 30 

#define MAX( a, b ) (((a)>(b)) ? (a) : (b))
#define MIN( a, b ) (((a)<(b)) ? (a) : (b))

typedef struct gbl_state_struct {
  double fx;
  double ldt;
  double dmin;
  int nf;
  int nl;
  double v[WORK_DIM][WORK_DIM];
  double q0[WORK_DIM];
  double q1[WORK_DIM];
  double qa, qb, qc, qd0, qd1, qf1;
  void* userHook;
} GlobalState;

static int prx_minfit(int m, int n, double machep, double tol, 
		      double ab[][WORK_DIM], double* q);

static void prx_min(int n, int j, int nits, double* d2, double* x1,
		    double* f1, int fk, 
		    double (*f)(double*,int,void* userHook),
		    double* x, double t, double machep, double h, 
		    GlobalState* gs);

static double prx_flin(int n, int j, double l, 
		       double (*f)(double*,int,void* userHook),
		       double* x, GlobalState* gs);

static void prx_sort(int m, int n,double* d, double v[][WORK_DIM]);

static void prx_quad(int n, double (*f)(double*, int, void* userHook), 
		     double* x, double t, double machep, double h, 
		     GlobalState* gs);

static void prx_vcprnt(int option, double* v, int n);

static void prx_print(int n, double* x, int prin, double fmin, 
		      GlobalState* gs);

static void prx_maprnt(int option, double v[][WORK_DIM], int m, int n);

double praxis(double t0, double machep, double h0, int n, int prin, 
	      double* x, double (*f)(double*,int,void* userHook), 
	      void (reset)(double*,int,void* userHook),
	      double fmin, void* userHook)
{
  /*
    REAL*8 FUNCTION PRAXIS(T0,MACHEP,H0,N,PRIN,X,F,FMIN)
    C                             LAST MODIFIED 3/1/73
    C
    C     PRAXIS RETURNS THE MINIMUM OF THE FUNCTION F(X,N) OF N VARIABLES
    C     USING THE PRINCIPAL AXIS METHOD.  THE GRADIENT OF THE FUNCTION IS
    C     NOT REQUIRED.
    C
    C     FOR A DESCRIPTION OF THE ALGORITHM, SEE CHAPTER SEVEN OF
    C     "ALGORITHMS FOR FINDING ZEROS AND EXTREMA OF FUNCTIONS WITHOUT
    C     CALCULATING DERIVATIVES" BY RICHARD P BRENT.
    C
    C     THE PARAMETERS ARE:
    C     T0       IS A TOLERANCE.  PRAXIS ATTEMPTS TO RETURN PRAXIS=F(X)
    C              SUCH THAT IF X0 IS THE TRUE LOCAL MINIMUM NEAR X, THEN
    C              NORM(X-X0) < T0 + SQUAREROOT(MACHEP)*NORM(X).
    C     MACHEP   IS THE MACHINE PRECISION, THE SMALLEST NUMBER SUCH THAT
    C              1 + MACHEP > 1.  MACHEP SHOULD BE 16.**-13 (ABOUT
    C              2.22D-16) FOR REAL*8 ARITHMETIC ON THE IBM 360.
    C     H0       IS THE MAXIMUM STEP SIZE.  H0 SHOULD BE SET TO ABOUT THE
    C              MAXIMUM DISTANCE FROM THE INITIAL GUESS TO THE MINIMUM.
    C              (IF H0 IS SET TOO LARGE OR TOO SMALL, THE INITIAL RATE OF
    C              CONVERGENCE MAY BE SLOW.)
    C     N        (AT LEAST TWO) IS THE NUMBER OF VARIABLES UPON WHICH
    C              THE FUNCTION DEPENDS.
    C     PRIN     CONTROLS THE PRINTING OF INTERMEDIATE RESULTS.
    C              IF PRIN=0, NOTHING IS PRINTED.
    C              IF PRIN=1, F IS PRINTED AFTER EVERY N+1 OR N+2 LINEAR
    C              MINIMIZATIONS.  FINAL X IS PRINTED, BUT INTERMEDIATE X IS
    C              PRINTED ONLY IF N IS AT MOST 4.
    C              IF PRIN=2, THE SCALE FACTORS AND THE PRINCIPAL VALUES OF
    C              THE APPROXIMATING QUADRATIC FORM ARE ALSO PRINTED.
    C              IF PRIN=3, X IS ALSO PRINTED AFTER EVERY FEW LINEAR
    C              MINIMIZATIONS.
    C              IF PRIN=4, THE PRINCIPAL VECTORS OF THE APPROXIMATING
    C              QUADRATIC FORM ARE ALSO PRINTED.
    C     X        IS AN ARRAY CONTAINING ON ENTRY A GUESS OF THE POINT OF
    C              MINIMUM, ON RETURN THE ESTIMATED POINT OF MINIMUM.
    C     F(X,N,userHook)   IS THE FUNCTION TO BE MINIMIZED.  F SHOULD BE 
    C              A REAL*8 FUNCTION DECLARED EXTERNAL IN THE CALLING PROGRAM.
    C     reset(X,N,userHook) is a reset function associated with F
    C     FMIN     IS AN ESTIMATE OF THE MINIMUM, USED ONLY IN PRINTING
    C              INTERMEDIATE RESULTS.
    C     userHook provided to pass additional data to F and reset
    C     THE APPROXIMATING QUADRATIC FORM IS
    C              Q(X') = F(X,N) + (1/2) * (X'-X)-TRANSPOSE * A * (X'-X)
    C     WHERE X IS THE BEST ESTIMATE OF THE MINIMUM AND A IS
    C              INVERSE(V-TRANSPOSE) * D * INVERSE(V)
    C     (V(*,*) IS THE MATRIX OF SEARCH DIRECTIONS; D(*) IS THE ARRAY
    C     OF SECOND DIFFERENCES).  IF F HAS CONTINUOUS SECOND DERIVATIVES
    C     NEAR X0, A WILL TEND TO THE HESSIAN OF F AT X0 AS X APPROACHES X0.
    C
    C     IT IS ASSUMED THAT ON FLOATING-POINT UNDERFLOW THE RESULT IS SET
    C     TO ZERO.
    C     THE USER SHOULD OBSERVE THE COMMENT ON HEURISTIC NUMBERS AFTER
    C     THE INITIALIZATION OF MACHINE DEPENDENT NUMBERS.
    C
  */
  int i;
  int ii;
  int j;
  int k;
  int km1;
  int k2;
  int illc;
  int kl, kt, ktm, klmk;
  double s,sl,dn,f1,lds,t,h,sf,df;
  double m2,m4,small,vsmall,large,vlarge,scbd,ldfac,t2,dni,value;
  GlobalState gs;
  double d[WORK_DIM],y[WORK_DIM],z[WORK_DIM];
  /*
    C
    C.....IF N>20 OR IF N<20 AND YOU NEED MORE SPACE, CHANGE '20' TO THE
    C     LARGEST VALUE OF N IN THE NEXT CARD, IN THE CARD 'IDIM=20', AND
    C     IN THE DIMENSION STATEMENTS IN SUBROUTINES MINFIT,MIN,FLIN,QUAD.
    C
  */
  /*
    C
    C.....INITIALIZATION.....
    C     MACHINE DEPENDENT NUMBERS:
    C
  */
  if (n>WORK_DIM) {
    fprintf(stderr,
	    "praxis: called with dim %d greater than %d (max allowed)!\n",
	    n,WORK_DIM);
    exit(-1);
  }

  small= machep * machep;
  vsmall=small*small;
  large=1.0/small;
  vlarge=1.0/vsmall;
  m2=sqrt(machep);
  m4=sqrt(m2);
  /*
    C
    C     HEURISTIC NUMBERS:
    C     IF THE AXES MAY BE BADLY SCALED (WHICH IS TO BE AVOIDED IF
    C     POSSIBLE), THEN SET SCBD=10.  OTHERWISE SET SCBD=1.
    C     IF THE PROBLEM IS KNOWN TO BE ILL-CONDITIONED, SET ILLC=TRUE.
    C     OTHERWISE SET ILLC=FALSE.
    C     KTM IS THE NUMBER OF ITERATIONS WITHOUT IMPROVEMENT BEFORE THE
    C     ALGORITHM TERMINATES.  KTM=4 IS VERY CAUTIOUS; USUALLY KTM=1
    C     IS SATISFACTORY.
    C
  */
  scbd=SCBD;
  illc=ILLC;
  ktm=KTM;
  /*
    C
  */
  ldfac=0.01;
  if (illc) ldfac=0.1;
  kt=0;
 L30: /* Jump here to recover from bad problems */
  gs.userHook= userHook;
  gs.nl=0;
  gs.nf=1;
  t=small+fabs(t0);
  t2=t;
  gs.dmin=small;
  h=h0;
  if (h<100*t) h=100*t;
  gs.ldt=h;
  gs.fx=f(x,n,gs.userHook);
  /*
    C.....THE FIRST SET OF SEARCH DIRECTIONS V IS THE IDENTITY MATRIX.....
  */
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      gs.v[j][i]=0.0;
    }
    gs.v[i][i]=1.0;
  }
  d[0]=0.0;
  gs.qd0=0.0;
  for (i=0; i<n; i++) {
    gs.q0[i]=x[i];
    gs.q1[i]=x[i];
  }
  if (prin>0) prx_print(n,x,prin,fmin,&gs);
  /*
    C
    C.....THE MAIN LOOP STARTS HERE.....
  */
 L40:
  sf=d[1];
  d[1]=0.0;
  s=0.0;
  /*
    C
    C.....MINIMIZE ALONG THE FIRST DIRECTION V(*,1).
    C     FX MUST BE PASSED TO MIN BY VALUE.
  */
  (*reset)(x,n,gs.userHook);
  gs.fx=f(x,n,gs.userHook);
  gs.qf1=gs.fx;
  value=gs.fx;
  prx_min(n,0,2,&(d[0]),&s,&value,0,f,x,t,machep,h,&gs);
#ifdef never
  prx_print(n,x,prin,fmin,&gs);
  prx_maprnt(1,gs.v,WORK_DIM,n);
#endif
  if (s<=0.0) {
    for (i=0; i<n; i++)
      gs.v[0][i] *= -1;
  }
      
  if (!((sf>0.9*d[0]) && (0.9*sf<d[0]))) {
    for (i=1; i<n; i++)
      d[i]= 0.0;
  }

  /*
    C
    C.....THE INNER LOOP STARTS HERE.....
  */
  for (k=1; k<n; k++) {
    for (i=0; i<n; i++) 
      y[i]=x[i];
    sf=gs.fx;
    if (kt>0) illc=1;
  L80:
    kl=k;
    df=0.0;
    /*
      C
      C.....A RANDOM STEP FOLLOWS (TO AVOID RESOLUTION VALLEYS).
      C     PRAXIS ASSUMES THAT RANDOM RETURNS A RANDOM NUMBER UNIFORMLY
      C     DISTRIBUTED IN (0,1).
      C
    */
    if (illc) {
      for (i=0; i<n; i++) {
	s=(0.1*gs.ldt+t2*pow(10,kt))*(drand48()-0.5);
	z[i]=s;
	for (j=0; j<n; j++)
	  x[j]=x[j]+s*gs.v[i][j];
      }
      gs.fx=f(x,n,gs.userHook);
      gs.nf=gs.nf+1;
#ifdef never
      fprintf(stderr,"Random step! ldt= %g, t2=%g, kt= %d\n",gs.ldt,t2,kt);
      fprintf(stderr,"x is now:\n"); prx_vcprnt(17,x,n);
#endif
    }
    /*
      C
      C.....MINIMIZE ALONG THE "NON-CONJUGATE" DIRECTIONS V(*,K),...,V(*,N)
      C
    */
    for (k2=k; k2<n; k2++) {
      sl=gs.fx;
      s=0.0;
      value=gs.fx;
      prx_min(n,k2,2,&(d[k2]),&s,&value,0,f,x,t,machep,h,&gs);
#ifdef never
      fprintf(stderr,"non-conjugate search along %d; x is now:\n",k2);
      prx_vcprnt(17,x,n);
#endif
      if (!illc) {
	s=sl-gs.fx;
      }
      else {
	s=d[k2]*((s+z[k2])*(s+z[k2]));
      }
      if (!(df>s)) {
	df=s;
	kl=k2;
      }
    }
    if (!(illc || (df>fabs((100*machep)*gs.fx)))) {
      /*
	C
	C.....IF THERE WAS NOT MUCH IMPROVEMENT ON THE FIRST TRY, SET
	C     ILLC=TRUE AND START THE INNER LOOP AGAIN.....
	C
      */
      illc= 1;
      goto L80;
    }
    if (k==2 && prin>1) prx_vcprnt(1,d,n);
    /*
      C
      C.....MINIMIZE ALONG THE "CONJUGATE" DIRECTIONS V(*,1),...,V(*,K-1)
      C
    */
    km1=k-1;
    for (k2=0; k2<=km1; k2++) {
      s= 0;
      value= gs.fx;
      prx_min(n,k2,2,&(d[k2]),&s,&value,0,f,x,t,machep,h,&gs);
#ifdef never
      fprintf(stderr,"conjugate search along %d; x is now:\n",k2);
      prx_vcprnt(17,x,n);
#endif
    }
    f1=gs.fx;
    gs.fx=sf;
    lds=0;
    for (i=0; i<n; i++) {
      sl=x[i];
      x[i]=y[i];
      sl=sl-y[i];
      y[i]=sl;
      lds=lds+sl*sl;
    }
    lds=sqrt(lds);
    if (lds>small) {
      /*
	C
	C.....DISCARD DIRECTION V(*,KL).
	C     IF NO RANDOM STEP WAS TAKEN, V(*,KL) IS THE "NON-CONJUGATE"
	C     DIRECTION ALONG WHICH THE GREATEST IMPROVEMENT WAS MADE.....
	C
      */
      klmk=kl-k;
      if (klmk>=1) {
	for (ii=1; ii<=klmk; ii++) {
	  i=kl-ii;
	  for (j=0;j<n;j++) {
	    gs.v[i+1][j]=gs.v[i][j];
	  }
	  d[i+1]=d[i];
	}
      }
      d[k]=0;
      for (i=0;i<n;i++) 
	gs.v[k][i]=y[i]/lds;
      /*
	C
	C.....MINIMIZE ALONG THE NEW "CONJUGATE" DIRECTION V(*,K), WHICH IS
	C     THE NORMALIZED VECTOR:  (NEW X) - (0LD X).....
	C
      */
      value=f1;
      prx_min(n,k,4,&(d[k]),&lds,&value,1,f,x,t,machep,h,&gs);
#ifdef never
      fprintf(stderr,"new dir conjugate search along %d; x is now:\n",k2);
      prx_vcprnt(17,x,n);
#endif
      if (lds<=0.0) {
	lds=-lds;
	for (i=0; i<n; i++)
	  gs.v[k][i]=-gs.v[k][i];
      }
    }
    gs.ldt=ldfac*gs.ldt;
    if (gs.ldt<lds) gs.ldt=lds;
    if (prin>0) prx_print(n,x,prin,fmin,&gs);
    t2=0.0;
    for (i=0; i<n; i++) 
      t2=t2+x[i]*x[i];
    t2=m2*sqrt(t2)+t;
    /*
      C
      C.....SEE WHETHER THE LENGTH OF THE STEP TAKEN SINCE STARTING THE
      C     INNER LOOP EXCEEDS HALF THE TOLERANCE.....
      C
    */
    if (gs.ldt>(0.5*t2)) kt=-1;
    kt=kt+1;
    if (kt>ktm) goto L400;
  }
  /*
    C.....THE INNER LOOP ENDS HERE.
    C
    C     TRY QUADRATIC EXTRAPOLATION IN CASE WE ARE IN A CURVED VALLEY.
    C
  */
  prx_quad(n,f,x,t,machep,h,&gs);
  dn=0.0;
  for (i=0; i<n; i++) {
    d[i]=1.0/sqrt(d[i]);
    if (dn<d[i]) dn=d[i];
  }
  if (prin>3) prx_maprnt(1,gs.v,WORK_DIM,n);
  for (j=0; j<n; j++) {
    s=d[j]/dn;
    for (i=0; i<n; i++) 
      gs.v[j][i]=s*gs.v[j][i];
  }

  /*
    C
    C.....SCALE THE AXES TO TRY TO REDUCE THE CONDITION NUMBER.....
    C
  */
  if (!(scbd<=1.0)) {
    s=vlarge;
    for (i=0; i<n; i++) {
      sl=0.0;
      for (j=0; j<n; j++)
	sl=sl+gs.v[j][i]*gs.v[j][i];
      z[i]=sqrt(sl);
      if (z[i]<m4) z[i]=m4;
      if (s>z[i]) s=z[i];
    }
    for (i=0; i<n; i++) {
      sl=s/z[i];
      z[i]=1.0/sl;
      if (z[i]>scbd) {
	sl=1.0/scbd;
	z[i]=scbd;
      }
      for (j=0; j<n; j++) 
	gs.v[j][i]=sl*gs.v[j][i];
    }
  }
  /*
    C
    C.....CALCULATE A NEW SET OF ORTHOGONAL DIRECTIONS BEFORE REPEATING
    C     THE MAIN LOOP.
    C     FIRST TRANSPOSE V FOR MINFIT:
    C
  */
  for (i=1; i<n; i++) {
    int im1=i-1;
    for (j=0; j<=im1; j++) {
      s=gs.v[j][i];
      gs.v[j][i]=gs.v[i][j];
      gs.v[i][j]=s;
    }
  }
  /*
    C
    C.....CALL MINFIT TO FIND THE SINGULAR VALUE DECOMPOSITION OF V.
    C     THIS GIVES THE PRINCIPAL VALUES AND PRINCIPAL DIRECTIONS OF THE
    C     APPROXIMATING QUADRATIC FORM WITHOUT SQUARING THE CONDITION
    C     NUMBER.....
    C
    C     This has been modified to restart if prx_minfit doesn't converge
  */
  if (!prx_minfit(WORK_DIM,n,machep,vsmall,gs.v,d)) {
    fprintf(stderr,"praxis: Bailing out of bad SVD convergence!\n");
    if (gs.ldt>(0.5*t2)) kt=-1;
    kt += 1;
    if (kt>ktm) goto L400;
    goto L30; 
  }
  /*
    C
    C.....UNSCALE THE AXES.....
    C
  */
  if (scbd>1.0) {
    for (i=0; i<n; i++) {
      s=z[i];
      for (j=0; j<n; j++)
	gs.v[j][i]=s*gs.v[j][i];
    }
    for (i=0; i<n; i++) {
      s=0.0;
      for (j=0; j<n; j++) 
	s=s+gs.v[i][j]*gs.v[i][j];
      s=sqrt(s);
      d[i]=s*d[i];
      s=1/s;
      for (j=0; j<n; j++)
	gs.v[i][j]=s*gs.v[i][j];
    }
  }
  /*
    C
  */
  for (i=0; i<n; i++) {
    dni=dn*d[i];
    if (dni<=large) {
      if (dni>=small) {
	d[i]=1/(dni*dni);
      }
      else {
          d[i]=vlarge;
      }
    }
    else {
      d[i]=vsmall;
    }
  }
  /*
    C
    C.....SORT THE EIGENVALUES AND EIGENVECTORS.....
    C
  */
  prx_sort(WORK_DIM,n,d,gs.v);
  gs.dmin=d[n-1];
  if (gs.dmin<small) gs.dmin=small;
  illc=0;
  if (m2*d[0]>gs.dmin) illc=1;
  if (prin>1 && scbd>1.0) prx_vcprnt(2,z,n);
  if (prin>1) prx_vcprnt(3,d,n);
  if (prin>3) prx_maprnt(2,gs.v,WORK_DIM,n);
  /*
    C.....THE MAIN LOOP ENDS HERE.....
    C
  */
  goto L40;
  /*
    C
    C.....RETURN.....
    C
  */
  L400:
  if (prin>0) prx_vcprnt(4,x,n);
  return gs.fx;
}

static int prx_minfit(int m, int n, double machep, double tol, 
		       double ab[][WORK_DIM], double* q)
{
  /*
    C...AN IMPROVED VERSION OF MINFIT (SEE GOLUB AND REINSCH, 1969)
    C   RESTRICTED TO M=N,P=0.
    C   THE SINGULAR VALUES OF THE ARRAY AB ARE RETURNED IN Q AND AB IS
    C   OVERWRITTEN WITH THE ORTHOGONAL MATRIX V SUCH THAT U.DIAG(Q) = AB.V,
    C   WHERE U IS ANOTHER ORTHOGONAL MATRIX.
  */
  int work_dim= WORK_DIM;
  int one= 1;
  double scratch[10*WORK_DIM];
  double junk[1];
  int scratch_dim= 10*WORK_DIM;
  double t;
  int i;
  int j;
  int info;
  if (n==1) {
    q[0] = ab[0][0];
    ab[0][0] = 1.0;
  }
  else {
    DGESVD("N","O",&n,&n,(double*)(&ab[0][0]),&work_dim,q,junk,&one,
	   junk,&one,scratch,&scratch_dim,&info);
    if (info<0) {
      fprintf(stderr,
	      "praxis: internal error: DGESVD argument %d illegal value!\n",
	      -info);
      exit(-1);
    }
    else if (info>0) {
      fprintf(stderr,
	      "praxis: DBDSQR did not converge; %d superdiagonals failed.\n",
	      info);
      /* V has been corrupted, and we don't have another copy.  We'll return
       * something generic and hope that the next reset avoids this
       * singularity.
       */
      for (i=0; i<n; i++)
	for (j=0; j<n; j++) ab[i][j]= 0.0;
      for (i=0; i<n; i++) ab[i][i]= 1.0;
      return 0; /* signalling error */
      
    }
    /* Transpose the V matrix */
    for (i=0; i<n; i++)
      for (j=i+1; j<n; j++) {
	t= ab[j][i];
	ab[j][i]= ab[i][j];
	ab[i][j]= t;
      }
  }
  return 1;
}

static void prx_min(int n, int j, int nits, double* d2, double* x1,
		    double* f1, int fk, 
		    double (*f)(double*,int,void* userHook),
		    double* x, double t, double machep, double h, 
		    GlobalState* gs)
{
  /*
    C...THE SUBROUTINE MIN MINIMIZES F FROM X IN THE DIRECTION V(*,J) UNLESS
    C   J IS LESS THAN 1 (**0**), WHEN A QUADRATIC SEARCH IS MADE IN THE PLANE
    C   DEFINED BY Q0,Q1,X.
    C   D2 IS EITHER ZERO OR AN APPROXIMATION TO HALF F".
    C   ON ENTRY, X1 IS AN ESTIMATE OF THE DISTANCE FROM X TO THE MINIMUM
    C   ALONG V(*,J) (OR, IF J=0, A CURVE).  ON RETURN, X1 IS THE DISTANCE
    C   FOUND.
    C   IF FK=.TRUE., THEN F1 IS FLIN(X1).  OTHERWISE X1 AND F1 ARE IGNORED
    C   ON ENTRY UNLESS FINAL FX IS GREATER THAN F1.
    C   NITS CONTROLS THE NUMBER OF TIMES AN ATTEMPT WILL BE MADE TO HALVE
    C   THE INTERVAL.
  */
  int dz;
  double m2, m4;
  double sf1, sx1;
  int k;
  int i;
  double xm;
  double fm;
  double f0;
  double s;
  double temp;
  double t2;
  double x2;
  double f2;
  double d1;
  double small = machep*machep;
  m2 = sqrt(machep);
  m4 = sqrt(m2);
  sf1 = *f1;
  sx1 = *x1;
  k = 0;
  xm = 0.0;
  fm = gs->fx;
  f0 = gs->fx;
  dz= (*d2<machep);
  /*
    C...FIND THE STEP SIZE...
  */
  s = 0.0;
  for (i=0; i<n; i++)
    s = s + x[i]*x[i];
  s = sqrt(s);
  temp = *d2;
  if (dz) temp = gs->dmin;
  t2 = m4*sqrt(fabs(gs->fx)/temp + s*gs->ldt) + m2*gs->ldt;
  s = m4*s + t;
  if (dz && (t2>s)) t2 = s;
  t2 = MAX(t2,small);
  t2 = MIN(t2,.010*h);
  if (!((!fk) || (*f1>fm))) {
    xm = *x1;
    fm = *f1;
  }
  if (!(fk && (fabs(*x1)>=t2))) {
    temp=1.0;
    if (*x1<0.0) temp=-1.0;
    *x1=temp*t2;
    *f1 = prx_flin(n,j,*x1,f,x,gs);
  }
  if (*f1<=fm) {
    xm = *x1;
    fm = *f1;
  }
 L4:
  if (dz) {
    /*
      C...EVALUATE FLIN AT ANOTHER POINT AND ESTIMATE THE SECOND DERIVATIVE...
    */
    x2 = -*x1;
    if (f0>=*f1) x2 = 2.0*(*x1);
    f2 = prx_flin(n,j,x2,f,x,gs);
    if (f2<=fm) {
      xm = x2;
      fm = f2;
    }
    *d2 = (x2*(*f1 - f0)-(*x1)*(f2 - f0))/(((*x1)*x2)*((*x1) - x2));
  }
  /*
    C...ESTIMATE THE FIRST DERIVATIVE AT 0...
  */
  d1 = (*f1 - f0)/(*x1) - (*x1)*(*d2);
  dz = 1;
  /*
    C...PREDICT THE MINIMUM...
  */
  if (*d2<=small) {
    x2 = h;
    if (d1>=0.0) x2 = -x2;
  }
  else {
    x2 = (-0.50*d1)/(*d2);
  }
  if (fabs(x2)>h) {
    if (x2<=0.0) {
      x2 = -h;
    }
    else {
      x2 = h;
    }
  }
 L11:
  /*
    C...EVALUATE F AT THE PREDICTED MINIMUM...
  */
  f2 = prx_flin(n,j,x2,f,x,gs);
  if (!(k>=nits || f2<=f0)) {
    /*
      C...NO SUCCESS, SO TRY AGAIN...
    */
    k = k + 1;
    if (f0<*f1 && ((*x1)*x2)>0.0) goto L4;
    x2 = 0.5*x2;
    goto L11;
  }
  /*
    C...INCREMENT THE ONE-DIMENSIONAL SEARCH COUNTER...
  */
  gs->nl = gs->nl + 1;
  if (f2>fm) {
    x2 = xm;
  }
  else {
    fm = f2;
  }
  /*
    C...GET A NEW ESTIMATE OF THE SECOND DERIVATIVE...
  */
  if (fabs(x2*(x2 - (*x1)))>small) {
    (*d2) = (x2*(*f1-f0) - (*x1)*(fm-f0))/(((*x1)*x2)*((*x1) - x2));
  }
  else {
    if (k>0) (*d2) = 0.0;
  }
  if ((*d2)<=small) *d2 = small;
  *x1 = x2;
  gs->fx = fm;
  if (sf1<gs->fx) {
    gs->fx = sf1;
    *x1 = sx1;
  }
  /*
    C...UPDATE X FOR LINEAR BUT NOT PARABOLIC SEARCH...
  */
  if (j>=0) {
    for (i=0; i<n; i++) 
      x[i] = x[i] + (*x1)*gs->v[j][i];
  }
}

static double prx_flin(int n, int j, double l, 
		       double (*f)(double*,int,void* userHook),
		       double* x, GlobalState* gs)
{
  /*
    C...FLIN IS THE FUNCTION OF ONE REAL VARIABLE L THAT IS MINIMIZED
    C   BY THE SUBROUTINE MIN...
  */
  double t[WORK_DIM];
  int i;
  if (j>=0) {
    /*
      C...THE SEARCH IS LINEAR...
    */
    for (i=0; i<n; i++) 
      t[i] = x[i] + l*gs->v[j][i];
  }
  else {
    /*
      C...THE SEARCH IS ALONG A PARABOLIC SPACE CURVE...
    */
    gs->qa = (l*(l - gs->qd1))/(gs->qd0*(gs->qd0 + gs->qd1));
    gs->qb = ((l + gs->qd0)*(gs->qd1 - l))/(gs->qd0*gs->qd1);
    gs->qc = (l*(l + gs->qd0))/(gs->qd1*(gs->qd0 + gs->qd1));
    for (i=0; i<n; i++) 
      t[i] = (gs->qa*gs->q0[i] + gs->qb*x[i]) + gs->qc*gs->q1[i];
  }
  /*
    C...THE FUNCTION EVALUATION COUNTER NF IS INCREMENTED...
  */
  gs->nf += 1;
  return (*f)(t,n,gs->userHook); 
}

static void prx_sort(int m, int n,double* d, double v[][WORK_DIM])
{
  /*
    c...SORTS THE ELEMENTS OF D(N) INTO DESCENDING ORDER AND MOVES THE
    C   CORRESPONDING COLUMNS OF V(N,N).
    C   M IS THE ROW DIMENSION OF V AS DECLARED IN THE CALLING PROGRAM.
  */
  double s;
  int nm1;
  int ip1;
  int i;
  int j;
  int k;
  if (n==1) return;
  nm1 = n - 1;
  for (i=0; i<nm1; i++) {
    k=i;
    s = d[i];
    ip1 = i + 1;
    for (j=ip1; j<n; j++) {
      if (d[j]>s) {
	k = j;
	s = d[j];
      }
    }
    if (k>i) {
      d[k] = d[i];
      d[i] = s;
      for (j=0; j<n; j++) {
	s = v[i][j];
	v[i][j] = v[k][j];
	v[k][j] = s;
	}
    }
  }
}

static void prx_quad(int n, double (*f)(double*,int,void* userHook), 
		     double* x, double t, double machep, double h, 
		     GlobalState* gs)
{
  /*
    C...QUAD LOOKS FOR THE MINIMUM OF F ALONG A CURVE DEFINED BY Q0,Q1,X...
  */
  double l;
  int i;
  double s;
  double value;
  s = gs->fx;
  gs->fx = gs->qf1;
  gs->qf1 = s;
  gs->qd1 = 0.0;
  for (i=0; i<n; i++) {
    s = x[i];
    l = gs->q1[i];
    x[i] = l;
    gs->q1[i] = s;
    gs->qd1 = gs->qd1 + (s-l)*(s-l);
  }
  gs->qd1 = sqrt(gs->qd1);
  l = gs->qd1;
  s = 0.0;
  if (!(gs->qd0<=0.0 || gs->qd1<=0.0 || gs->nl < 3*n*n)) {
    value=gs->qf1;
    prx_min(n,0,2,&s,&l,&value,1,f,x,t,machep,h,gs);
    gs->qa = (l*(l-gs->qd1))/(gs->qd0*(gs->qd0+gs->qd1));
    gs->qb = ((l+gs->qd0)*(gs->qd1-l))/(gs->qd0*gs->qd1);
    gs->qc = (l*(l+gs->qd0))/(gs->qd1*(gs->qd0+gs->qd1));
  }
  else {
    gs->fx = gs->qf1;
    gs->qa = 0.0;
    gs->qb = gs->qa;
    gs->qc = 1.0;
  }
  gs->qd0 = gs->qd1;
  for (i=0; i<n; i++) {
    s = gs->q0[i];
    gs->q0[i] = x[i];
    x[i] = (gs->qa*s + gs->qb*x[i]) + gs->qc*gs->q1[i];
  }
}

static void prx_vcprnt(int option, double* v, int n)
{
  int i;
  switch (option) {
  case 1:
    fprintf(stderr,"praxis: the second difference array d(*) is:\n");
    break;
  case 2:
    fprintf(stderr,"praxis: the scale factors are:\n");
    break;
  case 3:
    fprintf(stderr,
      "praxis: the approximating quadratic form has the principal values:\n");
    break;
  case 4:
    fprintf(stderr,"praxis: x is:\n");
    break;
  }
  for (i=0; i<n; i++) {
    fprintf(stderr,"   %g",v[i]);
    if (!((i+1)%6)) fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");
}

static void prx_print(int n, double* x, int prin, double fmin, 
		      GlobalState* gs)
{
  double ln;
  int i;
  fprintf(stderr,
"praxis: after %d linear searches, the function has been evaluated %d times.\n",
	  gs->nl, gs->nf);
  fprintf(stderr,"praxis: the smallest value found is f(x)= %g\n",gs->fx);
  if (gs->fx>fmin) {
    ln= log10(gs->fx-fmin);
    fprintf(stderr,"praxis: log (f(x)-%g) = %g\n",fmin,ln);
  }
  else {
    fprintf(stderr,"praxis: log (f(x)-%g) is undefined\n",fmin);
  }
  if (n>4 && prin<=2) return;
  fprintf(stderr,"praxis: x is:\n");
  for (i=0; i<n; i++) {
    fprintf(stderr,"   %g",x[i]);
    if (!((i+1)%5)) fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");
}

static void prx_maprnt(int option, double v[][WORK_DIM], int m, int n)
{
  int low;
  int upp;
  int i;
  int j;
  /*
    C...THE SUBROUTINE MAPRNT PRINTS THE COLUMNS OF THE NXN MATRIX V
    C   WITH A HEADING AS SPECIFIED BY OPTION.
    C   M IS THE ROW DIMENSION OF V AS DECLARED IN THE CALLING PROGRAM...
  */
  low = 0;
  upp = 6;
  switch (option) {
  case 1:
    fprintf(stderr,"praxis: the new directions are:\n");
    break;
  case 2:
    fprintf(stderr,"praxis: and the principal axes:\n");
    break;
  }
 L3:
  if (n<upp) upp = n;
  for (i=0; i<n; i++) {
    for (j=low; j<upp; j++)
      fprintf(stderr,"  %g",v[j][i]);
    fprintf(stderr,"\n");
  }
  low = upp;
  if (n<=low) {
    fprintf(stderr,"\n");
    return;
  }
  upp = upp + 6;
  goto L3;
}

#ifdef never
      REAL*8 FUNCTION RANDOM(NAUGHT)
      REAL*8 RAN1,RAN3(127),HALF
      INTEGER RAN2,Q,R
      LOGICAL INIT
      DATA INIT/.FALSE./
      IF (INIT) GO TO 3
      R = MOD(NAUGHT,8190) + 1
      RAN2 = 128
      DO 2 I=1,127
         RAN2 = RAN2 - 1
         RAN1 = -2.D0**55
         DO 1 J=1,7
            R = MOD(1756*R,8191)
            Q = R/32
1           RAN1 = (RAN1 + Q)*(1.0D0/256)
2        RAN3(RAN2) = RAN1
      INIT = .TRUE.
3     IF (RAN2.EQ.1) RAN2 = 128
      RAN2 = RAN2 - 1
      RAN1 = RAN1 + RAN3(RAN2)
      HALF = .5D0
      IF (RAN1.GE.0.D0) HALF = -HALF
      RAN1 = RAN1 + HALF
      RAN3(RAN2) = RAN1
      RANDOM = RAN1 + .5D0
      RETURN
      END
#endif
