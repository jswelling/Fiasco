/************************************************************
 *                                                          *
 *  lapack.h                                                *
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
 *  Original programming by Joel Welling, 6/98              *
 ************************************************************/
/* This module includes prototypes for LAPACK entry points, and
 * provides name redefinition for portability.
 * various systems. 
 */

#if defined(FORTRAN_ADD_UNDERSCORE)

#define SGESVD sgesvd_
#define CGESVD cgesvd_
#define DGESVD dgesvd_
#define ZGESVD zgesvd_
#define SLAMCH slamch_
#define DLAMCH dlamch_

#define  DASUM   dasum_
#define  DDOT    ddot_
#define  DGEMM   dgemm_
#define  ZGEMM   zgemm_
#define  DGEMV   dgemv_
#define  ZGEMV   zgemv_
#define  DNRM2   dnrm2_
#define  DSPMV   dspmv_
#define  DTPMV   dtpmv_
#define  DTRMM   dtrmm_
#define  DTRMV   dtrmv_

#define  DGETRF  dgetrf_
#define  DGETRI  dgetri_
#define  DGEQRF  dgeqrf_
#define  DGEQPF  dgeqpf_
#define  DGELS   dgels_
#define  DGELSX  dgelsx_
#define  DLANGE  dlange_

#define  DLAIC1  dlaic1_
#define  ILAENV  ilaenv_
#define  DORMQR  dormqr_
#define  DORGQR  dorgqr_

#define  DPOTRF  dpotrf_
#define  DPPTRF  dpptrf_
#define  DPPTRF  dpptrf_

#define  DSPCON  dspcon_
#define  DSPEV   dspev_
#define  DSPEVX  dspevx_
#define  DSPSV   dspsv_
#define  DSPTRI  dsptri_
#define  DSPTRF  dsptrf_
#define  DLANSP  dlansp_

#define  DSYEV   dsyev_
#define  DSYTRF  dsytrf_
#define  DSYTRI  dsytri_
#define  DSYTRS  dsytrs_

#define  DTRTRI  dtrtri_
#define  DTRTRS  dtrtrs_
#define  DSYEVX  dsyevx_

#else

#define SGESVD sgesvd
#define CGESVD cgesvd
#define DGESVD dgesvd
#define ZGESVD zgesvd
#define SLAMCH slamch
#define DLAMCH dlamch

#define  DASUM   dasum
#define  DDOT    ddot
#define  DGEMM   dgemm
#define  ZGEMM   zgemm
#define  DGEMV   dgemv
#define  ZGEMV   zgemv
#define  DNRM2   dnrm2
#define  DSPMV   dspmv
#define  DTPMV   dtpmv
#define  DTRMM   dtrmm
#define  DTRMV   dtrmv

#define  DGETRF  dgetrf
#define  DGETRI  dgetri
#define  DGEQRF  dgeqrf
#define  DGEQPF  dgeqpf
#define  DGELS   dgels
#define  DGELSX  dgelsx
#define  DLANGE  dlange

#define  DLAIC1  dlaic1
#define  ILAENV  ilaenv
#define  DORMQR  dormqr
#define  DORGQR  dorgqr

#define  DPOTRF  dpotrf
#define  DPPTRF  dpptrf
#define  DPPTRF  dpptrf

#define  DSPCON  dspcon
#define  DSPEV   dspev
#define  DSPEVX  dspevx
#define  DSPSV   dspsv
#define  DSPTRI  dsptri
#define  DSPTRF  dsptrf
#define  DLANSP  dlansp

#define  DSYEV   dsyev
#define  DSYTRF  dsytrf
#define  DSYTRI  dsytri
#define  DSYTRS  dsytrs

#define  DTRTRI  dtrtri
#define  DTRTRS  dtrtrs
#define  DSYEVX  dsyevx

#endif

typedef struct complex_struct {
  float r;
  float i;
} complex;

int SGESVD(char *jobu, char *jobvt, int *m, int *n, float
           *a, int * lda, float *s, float *u, int *ldu, float *vt, int
           *ldvt, float *work, int *lwork, int *info);

int CGESVD(char *jobu, char *jobvt, int *m, int *n,
           complex *a, int *lda, float *s, complex *u, int *ldu, complex
           * vt, int *ldvt, complex *work, int *lwork, float *rwork,
           int *info);

int DGESVD(char *jobu, char *jobvt, int *m, int *n, double
           *a, int * lda, double *s, double *u, int *ldu, double *vt, int
           *ldvt, double *work, int *lwork, int *info);

int ZGESVD(char *jobu, char *jobvt, int *m, int *n, double
           *a, int * lda, double *s, double *u, int *ldu, double *vt, int
           *ldvt, double *work, int* lwork, double* rwork, int *info);

int DSYEVX(char *jobz, char *range, char *uplo, int *n, double *a, 
	   int *lda, double *vl, double *vu, int *il, int *iu, 
	   double *abstol, int *m, double *w, double *z, int *ldz, 
	   double *work, int *lwork, int *iwork, int *ifail, int *info);

float SLAMCH(char *cmach);
double DLAMCH(char *cmach);

    /* BLAS Routines   */

double  DASUM( int *, double *, int * );
double  DDOT( int *, double *, int *, double *, int * );
void    DGEMM( char*, char*, int*, int*, int*, double*, double*, int*, double*,  int*, double*, double*, int* );
void    ZGEMM( char*, char*, int*, int*, int*, double*, double*, int*, double*,  int*, double*, double*, int* );
void    DGEMV( char *, int *, int *, double *, double *, int *, double *, int *, double *, double*, int * );
void    ZGEMV( char *, int *, int *, double *, double *, int *, double *, int *, double *, double*, int * );
double  DNRM2( int *, double *, int * );
double  DSPMV( char *, int *, double *, double *, double *, int *,  double *, double *, int * );
double  DTPMV( char *, char *, char *, int *, double *, double *, int *);
void    DTRMM( char *, char *, char *, char *, int *, int *, double *, double *, int *, double *, int * );
void    DTRMV( char *, char *, char *, int *, double *, double *, double * );

    /* LAPACK Routines */

void    DGETRF( int *, int *, double *, int *, int *, int * );
void    DGETRI( int *, double *, int *, int *, double *, int *, int * );
void    DGEQRF( int *, int *, double *, int *, double *, double *, int *, int * );
void    DGEQPF( int *, int *, double *, int *, int *, double *, double *, int * );
void    DGELS( char *, int *, int *, int *, double *, int *, double *, int *, double *, int *, int * );
void    DGELSX( int *, int *, int *, double *, int *, double *, int *, int *, double *, int *, double *, int * );
double  DLANGE( char *, int *, int *, double *, int *, double * );

void    DLAIC1( int *, int *, double *, double *, double *, double *, double *, double *, double * );
int     ILAENV( int *, char *, char *, int *, int *, int *, int * );

void    DORMQR( char *, char *, int *, int *, int *, double *, int *, double *,	double *, int *, double *, int *, int * );
void    DORGQR( int *, int *, int *, double *, int *, double *, double *, int *, int * );

void    DPOTRF( char *, int *, double *, int *, int * );
void    DPOTRS( char *, int *, int *, double *, int *, double *, int *, int * );
void    DPPTRF( char *, int *, double *, int * );
void    DPPTRF( char *, int *, double *, int * );

void    DSPCON( char *, int *, double *, int *, double *, double *, double *, int *, int * );
void    DSPEV( char *, char *, int *, double *, double *, double *, int *, double *, int * );
void    DSPEVX( char *, char *, char *, int *, double *, double *, double *, int *, int *, double *, int *, double *, double *, int *, double *, int *, int *, int * );
void    DSPSV( char *, int *, int *, double *, int *, double *, int *, int * );
void    DSPTRI( char *, int *, double *, int *, double *, int * );
void    DSPTRF( char *, int *, double *, int *, int * );
double  DLANSP( char *, char *, int *, double *, double * );

void    DSYEV( char *, char *, int *, double *, int *, double *, double *, int *, int * );
void    DSYTRF( char *, int *, double *, int *, int *, double *, int *, int * );
void    DSYTRI( char *, int *, double *, int *, int *, double *, int * );
void    DSYTRS( char *, int *, int *, double *, int *, int *, double *, int *, int * );

void    DTRTRI( char *, char *, int *, double *, int *, int * );
void    DTRTRS( char *, char *, char *, int *, int *, double *, int *, double *, int *, int * );

