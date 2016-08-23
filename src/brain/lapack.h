/************************************************************
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *                                                          *
 *     Copyright (c) 1999 Carnegie Mellon University        *
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
 ***********************************************************/

/*
 * C prototypes and definitions to support LAPACK and BLAS calls
 *
 */

/*** Character flags for LAPACK Routines ***/

static char   Left = 'L';
static char   Right = 'R';
static char   Trans = 'T';
static char   NoTrans = 'N';
static char   Upper = 'U';
static char   Lower = 'L';
static char   Unit = 'U';
static char   NoUnit = 'N';
static char   No = 'N';


/*** Name Translation ***/

    /*
     * Some linkers append an underscore to the names of fortran functions
     * so this needs to be added in the call.  Here, we automate the
     * process.  If the preprocessor macro FORTRAN_ADD_UNDERSCORE is defined
     * before this file is included, the functions prototyped below will
     * be translated as follows.
     *
     * To update the list when prototypes are added, copy the prototypes
     * below to here, eliminate the translations, and then apply replace-regexp
     * with the regular-expressions lapack-underscore-from-regexp and
     * lapack-underscore-to-regexp defined below in the local variables list.
     * The lines in the prototype list should not be broken.
     *
     */

#if defined(FORTRAN_ADD_UNDERSCORE)

#define  dasum   dasum_
#define  ddot    ddot_
#define  dgemm   dgemm_
#define  dgemv   dgemv_
#define  dnrm2   dnrm2_
#define  dspmv   dspmv_
#define  dtpmv   dtpmv_
#define  dtrmm   dtrmm_
#define  dtrmv   dtrmv_

#define  dgetrf  dgetrf_
#define  dgetri  dgetri_
#define  dgeqrf  dgeqrf_
#define  dgeqpf  dgeqpf_
#define  dgesvd  dgesvd_
#define  dgels   dgels_
#define  dgelsx  dgelsx_
#define  dlange  dlange_

#define  dlaic1  dlaic1_
#define  dlamch  dlamch_
#define  ilaenv  ilaenv_

#define  dormqr  dormqr_
#define  dorgqr  dorgqr_

#define  dpotrf  dpotrf_
#define  dpptrf  dpptrf_
#define  dpptrf  dpptrf_

#define  dspcon  dspcon_
#define  dspev   dspev_
#define  dspevx  dspevx_
#define  dspsv   dspsv_
#define  dsptri  dsptri_
#define  dsptrf  dsptrf_
#define  dlansp  dlansp_

#define  dsyev   dsyev_
#define  dsytrf  dsytrf_
#define  dsytri  dsytri_
#define  dsytrs  dsytrs_

#define  dtrtri  dtrtri_
#define  dtrtrs  dtrtrs_

#endif



/*** Function Prototypes ***/

    /* BLAS Routines   */

double  dasum( int *, double *, int * );
double  ddot( int *, double *, int *, double *, int * );
void    dgemm( char*, char*, int*, int*, int*, double*, double*, int*, double*,  int*, double*, double*, int* );
void    dgemv( char *, int *, int *, double *, double *, int *, double *, int *, double *, double*, int * );
double  dnrm2( int *, double *, int * );
double  dspmv( char *, int *, double *, double *, double *, int *,  double *, double *, int * );
double  dtpmv( char *, char *, char *, int *, double *, double *, int *);
void    dtrmm( char *, char *, char *, char *, int *, int *, double *, double *, int *, double *, int * );
void    dtrmv( char *, char *, char *, int *, double *, double *, double * );

    /* LAPACK Routines */

void    dgetrf( int *, int *, double *, int *, int *, int * );
void    dgetri( int *, double *, int *, int *, double *, int *, int * );
void    dgeqrf( int *, int *, double *, int *, double *, double *, int *, int * );
void    dgeqpf( int *, int *, double *, int *, int *, double *, double *, int * );
void    dgesvd( char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);
void    dgels( char *, int *, int *, int *, double *, int *, double *, int *, double *, int *, int * );
void    dgelsx( int *, int *, int *, double *, int *, double *, int *, int *, double *, int *, double *, int * );
double  dlange( char *, int *, int *, double *, int *, double * );

void    dlaic1( int *, int *, double *, double *, double *, double *, double *, double *, double * );
double  dlamch( char * );
int     ilaenv( int *, char *, char *, int *, int *, int *, int * );

void    dormqr( char *, char *, int *, int *, int *, double *, int *, double *,	double *, int *, double *, int *, int * );
void    dorgqr( int *, int *, int *, double *, int *, double *, double *, int *, int * );

void    dpotrf( char *, int *, double *, int *, int * );
void    dpotrs( char *, int *, int *, double *, int *, double *, int *, int * );
void    dpptrf( char *, int *, double *, int * );
void    dpptrf( char *, int *, double *, int * );

void    dspcon( char *, int *, double *, int *, double *, double *, double *, int *, int * );
void    dspev( char *, char *, int *, double *, double *, double *, int *, double *, int * );
void    dspevx( char *, char *, char *, int *, double *, double *, double *, int *, int *, double *, int *, double *, double *, int *, double *, int *, int *, int * );
void    dspsv( char *, int *, int *, double *, int *, double *, int *, int * );
void    dsptri( char *, int *, double *, int *, double *, int * );
void    dsptrf( char *, int *, double *, int *, int * );
double  dlansp( char *, char *, int *, double *, double * );

void    dsyev( char *, char *, int *, double *, int *, double *, double *, int *, int * );
void    dsytrf( char *, int *, double *, int *, int *, double *, int *, int * );
void    dsytri( char *, int *, double *, int *, int *, double *, int * );
void    dsytrs( char *, int *, int *, double *, int *, int *, double *, int *, int * );

void    dtrtri( char *, char *, int *, double *, int *, int * );
void    dtrtrs( char *, char *, char *, int *, int *, double *, int *, double *, int *, int * );



/*
 * Local Variables:        
 * mode: c
 * eval: (make-local-variable 'lapack-underscore-from-regexp)
 * eval: (make-local-variable 'lapack-underscore-to-regexp)
 * lapack-underscore-from-regexp: "\\sw+\\s +\\(\\sw+\\)(.*$"
 * lapack-underscore-to-regexp: "#define  \\1  \\1_"
 * End:
 *
 */
