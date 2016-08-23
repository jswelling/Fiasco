/************************************************************
 *                                                          *
 *  mri_esa.c                                               *
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

  DESCRIPTION OF MRI_ESA

  mri_esa performs eigen-systems analysis on Pgh MRI files.
**************************************************************/

#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "misc.h"
#include "stdcrg.h"
#include "../fmri/lapack.h"

#define KEYBUF_SIZE 512
#define MAX_AT_ONCE (64*1024*1024)

static char rcsid[] = "$Id: mri_esa.c,v 1.6 2007/07/06 18:45:53 welling Exp $";

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

static int get_chunk_type(MRI_Dataset* ds, char* chunk)
{
  char key_buf[KEYBUF_SIZE];

  safe_copy(key_buf,chunk);
  safe_concat(key_buf,".datatype");
  if (mri_has(ds, key_buf)) {
    char* type_name= mri_get_string(ds, key_buf);
    if (!strcmp(type_name,"uint8")) return MRI_UNSIGNED_CHAR;
    else if (!strcmp(type_name,"int16")) return MRI_SHORT;
    else if (!strcmp(type_name,"int32")) return MRI_INT;
    else if (!strcmp(type_name,"float32")) return MRI_FLOAT;
    else if (!strcmp(type_name,"float64")) return MRI_DOUBLE;
    else if (!strcmp(type_name,"int64")) return MRI_LONGLONG;
    else Abort("%s: unknown data type for key %s!\n",progname,key_buf);
  }
  else Abort("%s: missing tag %s!\n",progname,key_buf);

  return 0; /* not reached */
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

static void calc_sizes(MRI_Dataset* ds, const char* this_chunk, 
		       const char* dimstr, const char selected_dim,
		       long* fast_blocksize_out, long* slow_blocksize_out ) {
  /* This routine will fail if selected_dim is not in dimstr! */
  long fast_blocksize;
  long slow_blocksize;
  const char* this_dim;

  fast_blocksize= 1;

  this_dim= dimstr;
  while (*this_dim != selected_dim) 
    fast_blocksize *= safe_get_extent(ds,this_chunk,*this_dim++);

  this_dim++; /* step over selected dim */

  slow_blocksize= 1;
  while (*this_dim)
    slow_blocksize *= safe_get_extent(ds, this_chunk, *this_dim++);

  *fast_blocksize_out= fast_blocksize;
  *slow_blocksize_out= slow_blocksize;
}

static void eigensolve_once(const int in_rank, const int nEvals,
			    const int upperFlag, double* in_buf, 
			    double* eval_buf, double* evec_buf)
{
  char* jobz;
  char* range;
  char* uplo;
  int il= (in_rank+1)-nEvals;
  int iu= in_rank;
  static double abstol= 0.0;
  int m_out;
  int info;
  double dzero= 0.0;
  int in_rank_lcl= in_rank;
  int nEvals_lcl= nEvals;
  static int lwork= 0;
  static double* work= NULL;
  static int in_rank_old= 0;
  static int nEvals_old= 0;
  static int* iwork= NULL;
  static int* ifail= NULL;

  if (abstol==0.0) {
    abstol= 2.0*DLAMCH("S");
  }

  if (evec_buf) jobz= "V";
  else jobz= "N";

  if (nEvals==in_rank) range= "A";
  else range= "I";

  if (upperFlag) uplo= "U";
  else uplo= "L";

  if (in_rank != in_rank_old || nEvals != nEvals_old) {
    /* Work space size needs to be recalculated */
    if (work) free(work);
    if (iwork) free(iwork);
    if (ifail) free(ifail);
    lwork= 0;
    in_rank_old= in_rank;
    nEvals_old= nEvals;
  }
  if (lwork==0) {
    int neg_one= -1;
    double f_opt_worksize;

    if (!(iwork=(int*)malloc(5*in_rank*sizeof(int))))
      Abort("%s: unable to allocate %d bytes!\n",
	    progname, 5*in_rank*sizeof(int));
    if (!(ifail=(int*)malloc(5*in_rank*sizeof(int))))
      Abort("%s: unable to allocate %d bytes!\n",
	    progname, 5*in_rank*sizeof(int));

    /* Ask for optimum work space size, and allocate it */
    (void)DSYEVX(jobz, range, uplo, &in_rank_lcl, in_buf, &in_rank_lcl, 
		 &dzero, &dzero, &il, &iu, &abstol, &m_out, 
		 eval_buf, evec_buf, &in_rank_lcl, 
		 &f_opt_worksize, &neg_one, iwork, ifail, &info);
    lwork= (int)f_opt_worksize;
    if (debug_flg) 
      fprintf(stderr,"DSYEVX requests work space of %d doubles.\n",lwork);
    if (!(work= (double*)malloc(lwork*sizeof(double))))
      Abort("%s: unable to allocate %d bytes!\n",lwork*sizeof(double));
  }

  (void)DSYEVX(jobz, range, uplo, &in_rank_lcl, in_buf, &in_rank_lcl, 
	       &dzero, &dzero, &il, &iu, &abstol, &m_out, 
	       eval_buf, evec_buf, &in_rank_lcl, 
	       work, &lwork, iwork, ifail, &info);

  if (info != 0) {
    if (info<0) Abort("%s: DSYEVX error; argument %d has an illegal value\n",
		      progname, -info);
    else Error("%s: DSYEVX warning: %d eigenvectors failed to converge.\n",
		 progname, -info);
  }
  if (m_out != nEvals) 
    Error("%s: found only %d of requested %d eigenvalues\n",
	    progname, m_out, nEvals);
}

static void eigensolve_chunk(MRI_Dataset* in, MRI_Dataset* evals, 
			     MRI_Dataset* evecs, const char* chunk, 
			     const int nEvals, const int upperFlag,
			     const int descendFlag)
{
  char* dimstr= safe_get_dims(in,chunk);
  long in_rank;
  long in_fast_blksize= 0;
  long in_slow_blksize= 0;
  long long in_base_offset= 0;
  long long eval_base_offset= 0;
  long long evec_base_offset= 0;
  long in_slow;
  long in_offset= 0;
  double* eval_buf= NULL;
  double* eval_reverse_buf= NULL;
  double* evec_buf= NULL;
  double* evec_reverse_buf= NULL;
  double* in_buf= NULL;

  in_rank= safe_get_extent(in, chunk, dimstr[1]);
  calc_sizes(in, chunk, dimstr, dimstr[1], 
	     &in_fast_blksize, &in_slow_blksize);
  if (in_rank != in_fast_blksize)
    Abort("%s: internal error: matrix isn't square!\n",progname);

  if (!(eval_buf= (double*)malloc(nEvals*sizeof(double))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,nEvals*sizeof(double));
  if (descendFlag) {
    if (!(eval_reverse_buf= (double*)malloc(nEvals*sizeof(double))))
      Abort("%s: unable to allocate %d bytes!\n",
	    progname,nEvals*sizeof(double));
  }

  if (evecs) {
    if (!(evec_buf= (double*)malloc(nEvals*in_rank*sizeof(double))))
      Abort("%s: unable to allocate %d bytes!\n",
	    progname,nEvals*in_rank*sizeof(double));
    if (descendFlag) {
      if (!(evec_reverse_buf= (double*)malloc(nEvals*in_rank*sizeof(double))))
	Abort("%s: unable to allocate %d bytes!\n",
	      progname,nEvals*in_rank*sizeof(double));
    }
  }
  else evec_buf= NULL;
  
  if (debug_flg) 
    fprintf(stderr,
	    "input rank %d, solving for %d eigenvalue for each of %d blocks\n",
	    in_rank, nEvals, in_slow_blksize);
  if (verbose_flg)
    Message("# %d matrices to solve\n",in_slow_blksize);
  for (in_slow=0; in_slow<in_slow_blksize; in_slow++) {
    if (debug_flg)
      fprintf(stderr,"Reading %d from input at %lld\n",
	      in_fast_blksize*in_rank, in_base_offset);
    in_buf= mri_get_chunk(in, chunk, 
			  in_fast_blksize*in_rank,
			  in_base_offset, MRI_DOUBLE);
    eigensolve_once(in_rank, nEvals, upperFlag, in_buf, eval_buf, evec_buf);
    if (descendFlag) {
      int i;
      int j;
      for (i=0; i<nEvals; i++) eval_reverse_buf[i]= eval_buf[(nEvals-1)-i];
      mri_set_chunk(evals, chunk,
		    nEvals, eval_base_offset, MRI_DOUBLE, eval_reverse_buf);
      if (evec_buf) {
	for (i=0; i<nEvals; i++)
	  memcpy(evec_reverse_buf+(i*in_rank), evec_buf+((nEvals-1)-i)*in_rank,
		 in_rank*sizeof(double));
	mri_set_chunk(evecs, chunk,
		      in_rank*nEvals, evec_base_offset, MRI_DOUBLE, 
		      evec_reverse_buf);
      }
    }
    else {
      mri_set_chunk(evals, chunk,
		    nEvals, eval_base_offset, MRI_DOUBLE, eval_buf);
      if (evec_buf) {
	mri_set_chunk(evecs, chunk,
		      in_rank*nEvals, evec_base_offset, MRI_DOUBLE, evec_buf);
      }
    }
    in_base_offset += in_fast_blksize*in_rank;
    eval_base_offset += nEvals;
    evec_base_offset += in_rank*nEvals;
  }

  free(eval_buf);
  if (evec_buf) free(evec_buf);
  if (eval_reverse_buf) free(eval_reverse_buf);
  if (evec_reverse_buf) free(evec_reverse_buf);
}

static int chunk_check( MRI_Dataset* ds, const char* chunk )
{
  return( mri_has(ds, chunk) && !strcmp(mri_get_string(ds,chunk),"[chunk]") );
}

static int structure_check( MRI_Dataset* ds, const char* chunk, 
			    const int nEvals )
{
  char* dimstr= safe_get_dims(ds,chunk);
  int dim0= safe_get_extent(ds,chunk,dimstr[0]);
  int dim1= 0;
  char* here;

  /* We'll solve anything with at least two dimensions, as long
   * as the extents of the first two dimensions are the same.  The
   * user is responsible for making sure the input is symmetrical.
   */
  if (strlen(dimstr)<2) {
    Error("%s: input dataset is not a matrix!\n",progname);
    return 0;
  }
  else if (dim0 < 2) {
    Error("%s: first dimension of input has an extent of 1!\n",
	  progname);
    return 0;
  }
  else if (dim0 != (dim1=safe_get_extent(ds,chunk,dimstr[1]))) {
    Error("%s: first two dims of input dataset have different extents!\n",
	  progname);
    return 0;
  }
  else if (dim1 < 2) {
    Error("%s: second dimension of input has an extent of 1!\n",
	  progname);
    return 0;    
  }
  else if (nEvals > dim0) {
    Error("%s: requested %d eigenvalues, but rank is only %d!\n",
	  progname, nEvals, dim0);
    return 0;
  }
  else {
    if (verbose_flg) 
      Message("# Solving dimensions %c%c; matrices are %d by %d\n",
	      dimstr[0],dimstr[1],dim0,dim1);
    return 1;
  }
}

static void restructure_eval_dims( MRI_Dataset* evals, const char* chunk,
			      MRI_Dataset* in, const int nEvals )
{
  /* This routine changes the dimensions of the output to be those
   * appropriate for the eigenvalues of the input.  structure_check()
   * has already verified that the input dataset have the necessary
   * structure for these steps to work.
   */
  char key_buf[KEYBUF_SIZE];
  const char* dimstr= safe_get_dims(in,chunk);
  const char lostDim= dimstr[0];
  const char* here= dimstr+1;

  key_buf[KEYBUF_SIZE-1]= '\0';

  /* Truncate dimensions */
  if (snprintf(key_buf, KEYBUF_SIZE-1, "%s.dimensions",chunk)
      != strlen(chunk)+strlen(".dimensions"))
    Abort("%s: unreasonable long chunk name!\n",progname);
  mri_set_string(evals, key_buf, here);

  /* Delete the lost dimension's extent */
  if (snprintf(key_buf, KEYBUF_SIZE-1, "%s.extent.%c",chunk,lostDim)
      != strlen(chunk)+strlen(".extent.")+1)
    Abort("%s: unreasonable long chunk name!\n",progname);
  mri_remove(evals, key_buf);

  /* Set the output dimension */
  if (snprintf(key_buf, KEYBUF_SIZE-1, "%s.extent.%c",chunk,*here)
      != strlen(chunk)+strlen(".extent.")+1)
    Abort("%s: unreasonable long chunk name!\n",progname);
  mri_set_int(evals,key_buf,nEvals);
}

static void restructure_evec_dims( MRI_Dataset* evecs, const char* chunk,
				   MRI_Dataset* in, const int nEvals )
{
  /* This routine changes the dimensions of the output to be those
   * appropriate for the eigenvalues of the input.  structure_check()
   * has already verified that the input dataset have the necessary
   * structure for these steps to work.
   */
  char key_buf[KEYBUF_SIZE];
  const char* dimstr= safe_get_dims(in,chunk);
  int i;

  key_buf[KEYBUF_SIZE-1]= '\0';

  /* Set the size of the leading output dimension */
    if (snprintf(key_buf, KEYBUF_SIZE-1, "%s.extent.%c",chunk,dimstr[1])
	!= strlen(chunk)+strlen(".extent.")+1)
      Abort("%s: unreasonable long chunk name!\n",progname);
    mri_set_int(evecs,key_buf,nEvals);
}

int main( int argc, char* argv[] ) 
{
  char inName[512], evalName[512], evecName[512];
  int summed_dim;
  int nEvals= 0;
  int evecFlag= 0;
  int upperFlag= 0;
  int descendFlag= 0;
  MRI_Dataset *input = NULL, *eigenvals = NULL, *eigenvecs = NULL;
  char chunk[KEYBUF_SIZE];
  char key_buf[KEYBUF_SIZE];

  progname= argv[0];

  /* Print version number */
  if (verbose_flg) Message( "# %s\n", rcsid );

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Deprecate old options */

  if (cl_present( "n" ))
     Abort ("Option n(number) has been expanded to number|num.  Please see help file.\n");
  if (cl_present( "evec" ))
     Abort ("Option evec has been replaced by evc.  Please see help file.\n");


  verbose_flg= cl_present("verbose|ver|v");
  debug_flg= cl_present("debug|deb");
  cl_get("chunk|chu|c", "%option %s[%]","images",chunk);
  cl_get("number|num", "%option %d[%]",0,&nEvals);
  evecFlag= cl_get("eigenvectors|evc", "%option %s", evecName);
  upperFlag= cl_present("upper|upp");
  descendFlag= cl_present("descend");
  if (!cl_get("", "%s", inName)) {
    fprintf(stderr,"%s: Input file name not given.\n",progname);
    Help( "usage" );
    exit(-1);
  }
  if (!cl_get("", "%s", evalName)) {
    fprintf(stderr,"%s: Eigenvalue file name not given.\n",progname);
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
  
  /* Check name consistency and open input datasets */
  if( !strcmp( inName, evalName ) )
    Abort( "%s: input and eigenvalue files must be distinct.", progname );
  if (evecFlag) {
    if( !strcmp( inName, evecName ) )
      Abort( "%s: input and eigenvector files must be distinct.", progname );
    if( !strcmp( evecName, evalName ) )
      Abort( "%s: eigenvalue and eigenvector files must be distinct.", 
	     progname );
  }
  input = mri_open_dataset( inName, MRI_READ );

  /* Will they work with the program? */
  if (!chunk_check(input, chunk))
    Abort("%s: %s has no chunk %s!\n",progname,inName,chunk);
  if (!nEvals) {
    /* User didn't specify; calculate them all */
    const char* dimstr= safe_get_dims(input, chunk);
    nEvals= safe_get_extent(input,chunk,dimstr[0]);
  }
  if (!structure_check(input, chunk, nEvals))
    Abort("%s: Dataset %s cannot be solved.\n",progname,inName);

  /* Open output dataset.  We'll use float32 as the datatype unless
   * it's already double precision.
   */
  safe_copy(key_buf,chunk);
  safe_concat(key_buf,".datatype");
  eigenvals = mri_copy_dataset( evalName, input );
  hist_add_cl( eigenvals, argc, argv );
  if (strcmp(mri_get_string(input,key_buf),"float64"))
    mri_set_string(eigenvals,key_buf,"float32");
  restructure_eval_dims( eigenvals, chunk, input, nEvals );
  if (evecFlag) {
    eigenvecs = mri_copy_dataset( evecName, input );
    hist_add_cl( eigenvecs, argc, argv );
    if (strcmp(mri_get_string(input,key_buf),"float64"))
      mri_set_string(eigenvecs,key_buf,"float32");
    restructure_evec_dims( eigenvecs, chunk, input, nEvals );
  }
  else eigenvecs= NULL;

  /* Do the eigenvalue solution */
  eigensolve_chunk(input, eigenvals, eigenvecs, chunk, nEvals, upperFlag,
		   descendFlag);

  /* Write and close data-sets */
  mri_close_dataset( input );
  mri_close_dataset( eigenvals );
  if (eigenvecs) mri_close_dataset( eigenvecs );
  
  if (verbose_flg) Message( "#      Eigen-analysis complete.\n" );

  return 0;
}

