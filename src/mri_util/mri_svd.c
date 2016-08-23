/************************************************************
 *                                                          *
 *  mri_svd.c                                               *
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

  DESCRIPTION OF MRI_SVD

  mri_svd performs singular value decomposition on Pgh MRI files.
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

static char rcsid[] = "$Id: mri_svd.c,v 1.3 2007/07/06 18:45:53 welling Exp $";

static char* progname;
static int verbose_flg= 0;
static int debug_flg= 0;

#define MAYBE_MAKE_HASHMARK( i, n ) if (verbose_flg) makeHashMark(i,n)

typedef struct context_struct {
  int uFlag;
  int vtFlag;
  int wFlag;
  int complexFlag;
  int dimsSetFlag;
  MRI_Dataset *input, *uDSet, *vDSet, *wDSet;
  char dims[3]; /* extra '\0' at the end for convenience */
  long extents[2];
  char* chunk;
} Context;

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
    else if (!strcmp(type_name,"int64")) return MRI_LONGLONG;
    else if (!strcmp(type_name,"float32")) return MRI_FLOAT;
    else if (!strcmp(type_name,"float64")) return MRI_DOUBLE;
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

static void calc_sizes(Context* ctx,
		       long* fast_blocksize_out, 
		       long long* slow_blocksize_out ) {
  /* This routine will fail if selected_dim is not in dimstr! */
  long fast_blocksize;
  long long slow_blocksize;
  const char* this_dim;
  const char* dimstr= safe_get_dims(ctx->input, ctx->chunk);

  fast_blocksize= 1;

  this_dim= dimstr;
  while (*this_dim != ctx->dims[0]) 
    fast_blocksize *= safe_get_extent(ctx->input,ctx->chunk,*this_dim++);

  this_dim++; /* step over selected dims */
  this_dim++; /* step over selected dims */

  slow_blocksize= 1;
  while (*this_dim)
    slow_blocksize *= safe_get_extent(ctx->input, ctx->chunk, *this_dim++);

  *fast_blocksize_out= fast_blocksize;
  *slow_blocksize_out= slow_blocksize;
}

static void svdsolve_once(Context* ctx, 
			  long fast_blksize, long in_matrix_size, 
			  long u_matrix_size, long vt_matrix_size, 
			  long w_matrix_size, long work_buf_size, 
			  double* in_buf, double* u_out_buf, 
			  double* vt_out_buf, double* w_out_buf, 
			  double* a_buf, double* u_buf, double* vt_buf, 
			  double* w_buf, double* work_buf)
{
  long  fast_loop;
  int d0= (int)ctx->extents[0];
  int d1= (int)ctx->extents[1];
  int worksize= (int)work_buf_size;
  
  for (fast_loop=0; fast_loop<fast_blksize; fast_loop++) {
    long i;
    int lapack_retcode= 0;

    /* Copy a into place */
    for (i=0; i<in_matrix_size; i++) 
      a_buf[i]= in_buf[i*fast_blksize + fast_loop];

    /* Do SVD.  This destroys the factors_ftn matrix. */
    if (u_buf != NULL) {
      if (vt_buf != NULL) {
	(void)DGESVD("A", "A", &d0, &d1, a_buf, &d0, w_buf, u_buf, &d0,
		     vt_buf, &d1, work_buf, &worksize, &lapack_retcode);
      }
      else {
	(void)DGESVD("A", "N", &d0, &d1, a_buf, &d0, w_buf, u_buf, &d0,
		     NULL, &d1, work_buf, &worksize, &lapack_retcode);
      }
    }
    else {
      if (vt_buf != NULL) {
	(void)DGESVD("N", "A", &d0, &d1, a_buf, &d0, w_buf, NULL, &d0,
		     vt_buf, &d1, work_buf, &worksize, &lapack_retcode);
      }
      else {
	(void)DGESVD("N", "N", &d0, &d1, a_buf, &d0, w_buf, NULL, &d0,
		     NULL, &d1, work_buf, &worksize, &lapack_retcode);
      }
    }
    if (lapack_retcode<0) {
      Abort("%s: DGESVD argument %d had an illegal value",
	    progname, -lapack_retcode);
    }
    else if (lapack_retcode>0) {
      Abort("%s: DBDSQR did not converge in DGESVD; %d superdiagonals failed",
	    progname, lapack_retcode);
    }
    for (i=0; i<w_matrix_size; i++)
      w_out_buf[i*fast_blksize + fast_loop]= w_buf[i];
    if (u_buf != NULL && u_out_buf != NULL)
      for (i=0; i<u_matrix_size; i++)
	u_out_buf[i*fast_blksize + fast_loop]= u_buf[i];
    if (vt_buf != NULL && vt_out_buf != NULL)
      for (i=0; i<vt_matrix_size; i++)
	vt_out_buf[i*fast_blksize + fast_loop]= vt_buf[i];
  }
}

static void svdsolve_once_complex(Context* ctx, 
				  long fast_blksize, long in_matrix_size, 
				  long u_matrix_size, long vt_matrix_size, 
				  long w_matrix_size, long work_buf_size, 
				  double* in_buf, double* u_out_buf, 
				  double* vt_out_buf, double* w_out_buf, 
				  double* a_buf, double* u_buf, 
				  double* vt_buf, double* w_buf,
				  double* work_buf, double* rwork_buf)
{
  long  fast_loop;
  int d0= (int)ctx->extents[0];
  int d1= (int)ctx->extents[1];
  int worksize= (int)work_buf_size;
  
  for (fast_loop=0; fast_loop<fast_blksize; fast_loop++) {
    long i;
    int lapack_retcode= 0;

    /* Copy a into place */
    for (i=0; i<in_matrix_size; i+=2) {
      a_buf[i]= in_buf[i*fast_blksize + fast_loop];
      a_buf[i+1]= in_buf[i*fast_blksize + fast_loop + 1];
    }

    /* Do SVD.  This destroys the factors_ftn matrix. */
    if (u_buf != NULL) {
      if (vt_buf != NULL) {
	(void)ZGESVD("A", "A", &d0, &d1, a_buf, &d0, w_buf, u_buf, &d0,
		     vt_buf, &d1, work_buf, &worksize, rwork_buf,
		     &lapack_retcode);
      }
      else {
	(void)ZGESVD("A", "N", &d0, &d1, a_buf, &d0, w_buf, u_buf, &d0,
		     NULL, &d1, work_buf, &worksize, rwork_buf,
		     &lapack_retcode);
      }
    }
    else {
      if (vt_buf != NULL) {
	(void)ZGESVD("N", "A", &d0, &d1, a_buf, &d0, w_buf, NULL, &d0,
		     vt_buf, &d1, work_buf, &worksize, rwork_buf,
		     &lapack_retcode);
      }
      else {
	(void)ZGESVD("N", "N", &d0, &d1, a_buf, &d0, w_buf, NULL, &d0,
		     NULL, &d1, work_buf, &worksize, rwork_buf,
		     &lapack_retcode);
      }
    }
    if (lapack_retcode<0) {
      Abort("%s: ZGESVD argument %d had an illegal value",
	    progname, -lapack_retcode);
    }
    else if (lapack_retcode>0) {
      Abort("%s: DBDSQR did not converge in ZGESVD; %d superdiagonals failed",
	    progname, lapack_retcode);
    }
    for (i=0; i<w_matrix_size; i+=2) {
      w_out_buf[i*fast_blksize + fast_loop]= w_buf[i];
      w_out_buf[i*fast_blksize + fast_loop + 1]= w_buf[i+1];
    }
    if (u_buf != NULL && u_out_buf != NULL)
      for (i=0; i<u_matrix_size; i+=2) {
	u_out_buf[i*fast_blksize + fast_loop]= u_buf[i];
	u_out_buf[i*fast_blksize + fast_loop + 1]= u_buf[i+1];
      }
    if (vt_buf != NULL && vt_out_buf != NULL)
      for (i=0; i<vt_matrix_size; i+=2) {
	vt_out_buf[i*fast_blksize + fast_loop]= vt_buf[i];
	vt_out_buf[i*fast_blksize + fast_loop + 1]= vt_buf[i+1];
      }
  }
}

static double* safeAllocDoubles( long n )
{
  double* result;
  if (!(result=(double*)malloc(n*sizeof(double))))
    Abort("%s: unable to allocate %ld bytes!\n",n*sizeof(double));
  return result;
}

static void svdsolve_chunk(Context* ctx)
{
  long fast_blksize= 0;
  long long slow_blksize= 0;
  long long in_base_offset= 0;
  long long u_base_offset= 0;
  long long vt_base_offset= 0;
  long long w_base_offset= 0;
  long long slow_loop;
  long long in_offset= 0;
  double* a_buf= NULL; /* allocated here and used in solver */
  double* u_buf= NULL; /* allocated here and used in solver */
  double* vt_buf= NULL;  /* allocated here and used in solver */
  double* w_buf= NULL;  /* allocated here and used in solver */
  double* work_buf= NULL;  /* allocated here and used in solver */
  double* rwork_buf= NULL;
  long work_buf_size= 0;
  long rwork_buf_size= 0;
  long in_matrix_size;
  long u_matrix_size;
  long vt_matrix_size;
  long w_matrix_size;
  double* in_buf= NULL;
  double* u_out_buf= NULL;
  double* vt_out_buf= NULL;
  double* w_out_buf= NULL;
  long w_extent= (ctx->extents[0]>ctx->extents[1] ? 
		  ctx->extents[1]:ctx->extents[0]);

  in_matrix_size= ctx->extents[0]*ctx->extents[1];
  u_matrix_size= ctx->extents[0]*ctx->extents[0];
  vt_matrix_size= ctx->extents[1]*ctx->extents[1];
  w_matrix_size= w_extent;
  calc_sizes(ctx,&fast_blksize, &slow_blksize);
  if (ctx->complexFlag) {
    in_matrix_size *= 2;
    u_matrix_size *= 2;
    vt_matrix_size *= 2;
    w_matrix_size *= 2;
    fast_blksize /= 2;
  }

  /* Let the appropriate routine figure out the work array size */
  if (ctx->complexFlag) {
    int lapack_retcode;
    int lwork= -1;
    int d0= (int)ctx->extents[0];
    int d1= (int)ctx->extents[1];
    double work;

    (void)ZGESVD("A", "A", &d0, &d1, NULL, &d0, NULL, NULL,
		 &d0, NULL, &d1, &work, &lwork, NULL, &lapack_retcode);
    if (lapack_retcode != 0)
      Abort("%s: ZGESVD died trying to calculate workspace size!\n",progname);
    work_buf_size= 2*(long)work;
    if (debug_flg) fprintf(stderr,"ZGESVD says workspace size is %ld\n",
			   work_buf_size);
  }
  else {
    int lapack_retcode;
    int lwork= -1;
    int d0= (int)ctx->extents[0];
    int d1= (int)ctx->extents[1];
    double work;

    (void)DGESVD("A", "A", &d0, &d1, NULL, &d0, NULL, NULL,
		 &d0, NULL, &d1, &work, &lwork, &lapack_retcode);
    if (lapack_retcode != 0)
      Abort("%s: DGESVD died trying to calculate workspace size!\n",progname);
    work_buf_size= (long)work;
    if (debug_flg) fprintf(stderr,"DGESVD says workspace size is %ld\n",
			   work_buf_size);
  }

  a_buf= safeAllocDoubles( in_matrix_size );
  w_buf= safeAllocDoubles( w_matrix_size );
  work_buf= safeAllocDoubles( work_buf_size );
  if (ctx->uFlag) {
    u_buf= safeAllocDoubles( u_matrix_size );
    u_out_buf= safeAllocDoubles( u_matrix_size*fast_blksize );
  }
  if (ctx->vtFlag) {
    vt_buf= safeAllocDoubles( vt_matrix_size );
    vt_out_buf= safeAllocDoubles( vt_matrix_size*fast_blksize );
  }
  if (ctx->wFlag) {
    w_out_buf= safeAllocDoubles( w_matrix_size*fast_blksize );
  }
  if (ctx->complexFlag) {
    rwork_buf_size= 5*w_extent;
    rwork_buf= safeAllocDoubles( rwork_buf_size );
  }
  
  if (debug_flg) 
    fprintf(stderr,
	    "fast block size %ld, slow block size %lld, matrix size %ld, %s\n",
	    fast_blksize, slow_blksize,in_matrix_size,
	    (ctx->complexFlag?"comlex":"scalar"));
  if (verbose_flg)
    Message("# %lld matrices to solve, each %ld x %ld %s\n",
	    slow_blksize*fast_blksize,ctx->extents[0],ctx->extents[1],
	    (ctx->complexFlag?"comlex":"scalar"));

  in_base_offset= u_base_offset= vt_base_offset= w_base_offset= 0;
  for (slow_loop=0; slow_loop<slow_blksize; slow_loop++) {
    if (debug_flg)
      fprintf(stderr,"Reading %ld from input at %lld\n",
	      fast_blksize*in_matrix_size, in_base_offset);
    in_buf= mri_get_chunk(ctx->input, ctx->chunk, 
			  fast_blksize*in_matrix_size,
			  in_base_offset, MRI_DOUBLE);
    if (ctx->complexFlag) {
      svdsolve_once_complex(ctx, fast_blksize, in_matrix_size, u_matrix_size,
		    vt_matrix_size, w_matrix_size, work_buf_size, 
		    in_buf, u_out_buf, vt_out_buf, w_out_buf, 
		    a_buf, u_buf, vt_buf, w_buf, work_buf, rwork_buf);
    }
    else {
      svdsolve_once(ctx, fast_blksize, in_matrix_size, u_matrix_size,
		    vt_matrix_size, w_matrix_size, work_buf_size, 
		    in_buf, u_out_buf, vt_out_buf, w_out_buf, 
		    a_buf, u_buf, vt_buf, w_buf, work_buf);
    }
    if (ctx->uFlag) {
      if (debug_flg) fprintf(stderr,"Writing %ld to U output at %lld\n",
			     u_matrix_size*fast_blksize, u_base_offset);
      mri_set_chunk(ctx->uDSet, ctx->chunk, u_matrix_size*fast_blksize, 
		    u_base_offset, MRI_DOUBLE, u_out_buf);
      u_base_offset += u_matrix_size*fast_blksize;
    }
    if (ctx->vtFlag) {
      if (debug_flg) fprintf(stderr,"Writing %ld to V output at %lld\n",
			     vt_matrix_size*fast_blksize, vt_base_offset);
      mri_set_chunk(ctx->vDSet, ctx->chunk, vt_matrix_size*fast_blksize, 
		    vt_base_offset, MRI_DOUBLE, vt_out_buf);
      vt_base_offset += vt_matrix_size*fast_blksize;
    }
    if (ctx->wFlag) {
      if (debug_flg) fprintf(stderr,"Writing %ld to W output at %lld\n",
			     w_matrix_size*fast_blksize, w_base_offset);
      mri_set_chunk(ctx->wDSet, ctx->chunk, w_matrix_size*fast_blksize, 
		    w_base_offset, MRI_DOUBLE, w_out_buf);
      w_base_offset += w_matrix_size*fast_blksize;
    }
    in_base_offset += in_matrix_size*fast_blksize;
  }

  if (a_buf) free(a_buf);
  if (u_buf) free(u_buf);
  if (vt_buf) free(vt_buf);
  if (w_buf) free(w_buf);
  if (work_buf) free(work_buf);
  if (u_out_buf) free(u_out_buf);
  if (vt_out_buf) free(vt_out_buf);
  if (w_out_buf) free(w_out_buf);
  if (rwork_buf) free(rwork_buf);
  
}

static int chunk_check( MRI_Dataset* ds, const char* chunk )
{
  return( mri_has(ds, chunk) && !strcmp(mri_get_string(ds,chunk),"[chunk]") );
}

static int structure_check( Context* ctx )
{
  char* dimstr= safe_get_dims(ctx->input,ctx->chunk);
  int dim0= safe_get_extent(ctx->input,ctx->chunk,dimstr[0]);
  int dim1= 0;
  char* here;

  if (ctx->complexFlag) {
    ctx->dims[2]= '\0';
    if (ctx->dimsSetFlag && strchr(ctx->dims,'v')) {
      Error("%s: can't solve on dimension v in a complex dataset!\n",progname);
      return 0;
    }
    if ( dimstr[0] != 'v' || dim0 != 2 ) {
      Error("%s: input does not have v dim first, or is not complex!\n",
	    progname);
      return 0;
    }
    dimstr++;
    dim0= safe_get_extent(ctx->input,ctx->chunk,dimstr[0]);
  }

  if (ctx->dimsSetFlag) {
    if (!strstr(dimstr,ctx->dims)) {
      Error("%s: requested dimensions do not appear in order in the input!\n",
	    progname);
      return 0;
    }
  }
  else {
    if (strlen(dimstr)<2) {
      Error("%s: input chunk must have at least %d dimensions!\n",
	    progname,(ctx->complexFlag?2:3));
      return 0;
    }
    strncpy(ctx->dims,dimstr,2);
    ctx->dimsSetFlag= 1;
  }

  if (verbose_flg) {
    if (ctx->complexFlag)
      Message("# Solving complex dimensions %c%c; matrices are %d by %d\n",
	      ctx->dims[0],ctx->dims[1],
	      safe_get_extent(ctx->input,ctx->chunk,ctx->dims[0]),
	      safe_get_extent(ctx->input,ctx->chunk,ctx->dims[1]));
    else
      Message("# Solving dimensions %c%c; matrices are %d by %d\n",
	      ctx->dims[0],ctx->dims[1],
	      safe_get_extent(ctx->input,ctx->chunk,ctx->dims[0]),
	      safe_get_extent(ctx->input,ctx->chunk,ctx->dims[1]));
    
  }

  return 1;
}

static void restructure_u_dims( Context* ctx )

{
  /*
   * Dimensions of U are XabY where X is fast block, Y is slow block,
   * and ab are the live dims (with both having extent of a)
  */
  MRI_Dataset* ds= ctx->uDSet;
  char key_buf[KEYBUF_SIZE];

  key_buf[KEYBUF_SIZE-1]= '\0';
  snprintf(key_buf,KEYBUF_SIZE-1,"%s.extent.%c",ctx->chunk,ctx->dims[1]);
  mri_set_int(ds,key_buf,ctx->extents[0]);
}

static void restructure_v_dims( Context* ctx )
{
  /*
   * Dimensions of V are XabY where X is fast block, Y is slow block,
   * and ab are the live dims (with both having extent of b)
  */
  MRI_Dataset* ds= ctx->vDSet;
  char key_buf[KEYBUF_SIZE];

  key_buf[KEYBUF_SIZE-1]= '\0';
  snprintf(key_buf,KEYBUF_SIZE-1,"%s.extent.%c",ctx->chunk,ctx->dims[0]);
  mri_set_int(ds,key_buf,ctx->extents[1]);
}

static void restructure_w_dims(Context* ctx)
{
  /*
   * Dimensions of W are XbY where X is fast block, Y is slow block,
   * and b is the second live dim (with extent the smaller of a and b)
  */
  MRI_Dataset* ds= ctx->wDSet;
  char key_buf[KEYBUF_SIZE];
  char* dimstr= strdup(safe_get_dims(ds,ctx->chunk));
  char* here= strchr(dimstr,ctx->dims[0]);
  long extent= ( ctx->extents[0]>ctx->extents[1] ?
		 ctx->extents[1]:ctx->extents[0] );

  while (*(here+1)) {
    *here= *(here+1);
    here++;
  }
  *here= '\0';

  key_buf[KEYBUF_SIZE-1]= '\0';
  snprintf(key_buf,KEYBUF_SIZE-1,"%s.dimensions",ctx->chunk);
  mri_set_string(ds,key_buf,dimstr);
  snprintf(key_buf,KEYBUF_SIZE-1,"%s.extent.%c",ctx->chunk,ctx->dims[1]);
  mri_set_int(ds,key_buf,extent);
  snprintf(key_buf,KEYBUF_SIZE-1,"%s.extent.%c",ctx->chunk,ctx->dims[0]);
  mri_remove(ds,key_buf);

  free(dimstr);
}

int main( int argc, char* argv[] ) 
{
  char inName[512], uName[512], vName[512], wName[512], dims_in[512];
  Context ctx;
  char chunk_in[KEYBUF_SIZE];
  char key_buf[KEYBUF_SIZE];

  progname= argv[0];

  /* Initialize the context */
  ctx.uFlag= 0;
  ctx.vtFlag= 0;
  ctx.wFlag= 0;
  ctx.complexFlag= 0;
  ctx.dimsSetFlag= 0;
  ctx.input = ctx.uDSet= ctx.vDSet= ctx.wDSet= NULL;
  ctx.dims[0]= ctx.dims[1]= ctx.dims[2]= '\0';
  ctx.extents[0]= ctx.extents[1]= 0;
  ctx.chunk= NULL;

  /* Print version number */
  if (verbose_flg) Message( "# %s\n", rcsid );

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Deprecate old options */

  if (cl_present( "evec" ))
     Abort ("Option evec has been replaced by evc.  Please see help file.\n");

  verbose_flg= cl_present("verbose|ver|v");
  debug_flg= cl_present("debug|deb");
  ctx.complexFlag= cl_present("complex|cpx");
  
  cl_get("chunk|chu|c", "%option %s[%]","images",chunk_in);
  ctx.uFlag= cl_get("umatrix|umt", "%option %s", uName);
  ctx.vtFlag= cl_get("vmatrix|vmt", "%option %s", vName);
  ctx.wFlag= cl_get("wvector|wvc", "%option %s", wName);
  ctx.dimsSetFlag= cl_get("dimensions|dim|d","%option %s",dims_in);
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
  
  if ( !ctx.uFlag && !ctx.vtFlag && !ctx.wFlag )
    Abort("%s: at least one of umatrix, vmatrix, or wvector must be given.\n",
	  progname);

  if ( ctx.dimsSetFlag ) {
    if (strlen(dims_in)!=2 )
      Abort("%s: dimension string must have length 2!\n",progname);
    ctx.dims[0]= dims_in[0];
    ctx.dims[1]= dims_in[1];
  }

  /* Check name consistency and open input dataset */
  if (ctx.wFlag) {
    if (!strcmp( inName, wName ))
      Abort("%s: input and W vector files must be distinct.",progname);
    if (ctx.uFlag && !strcmp( uName, wName ))
      Abort("%s: U matrix and W vector files must be distinct.",progname);
    if (ctx.vtFlag && !strcmp( vName, wName ))
      Abort("%s: V matrix and W vector files must be distinct.",progname);
  }
  if (ctx.uFlag) {
    if (!strcmp( inName, uName ))
      Abort("%s: input and U matrix files must be distinct.",progname);
    if (ctx.vtFlag && !strcmp( vName, uName ))
      Abort("%s: U matrix and V matrix files must be distinct.",progname);
  }
  if (ctx.vtFlag) {
    if (!strcmp( inName, vName ))
      Abort("%s: input and V matrix files must be distinct.",progname);
  }
  ctx.input = mri_open_dataset( inName, MRI_READ );
  ctx.chunk= strdup(chunk_in);

  /* Will it work with the program? The structure check routine will
   * also set the dimensions to use if they are not set from the
   * command line.
   */
  if (!chunk_check(ctx.input, ctx.chunk))
    Abort("%s: %s has no chunk %s!\n",progname,inName,ctx.chunk);
  if (!structure_check(&ctx))
    Abort("%s: Dataset %s cannot be solved.\n",progname,inName);
  key_buf[KEYBUF_SIZE-1]= '\0';
  snprintf(key_buf,KEYBUF_SIZE-1,"%s.extent.%c",ctx.chunk,ctx.dims[0]);
  ctx.extents[0]= mri_get_int(ctx.input,key_buf);
  snprintf(key_buf,KEYBUF_SIZE-1,"%s.extent.%c",ctx.chunk,ctx.dims[1]);
  ctx.extents[1]= mri_get_int(ctx.input,key_buf);

  /* Open output datasets.  We'll use float32 as the datatype unless
   * it's already double precision.
   */
  if (ctx.uFlag) {
    ctx.uDSet= mri_copy_dataset( uName, ctx.input );
    hist_add_cl( ctx.uDSet, argc, argv );
    restructure_u_dims( &ctx );
  }
  if (ctx.vtFlag) {
    ctx.vDSet= mri_copy_dataset( vName, ctx.input );
    hist_add_cl( ctx.vDSet, argc, argv );
    restructure_v_dims( &ctx );
  }
  if (ctx.wFlag) {
    ctx.wDSet= mri_copy_dataset( wName, ctx.input );
    hist_add_cl( ctx.wDSet, argc, argv );
    restructure_w_dims( &ctx );
  }
  safe_copy(key_buf,ctx.chunk);
  safe_concat(key_buf,".datatype");
  if (strcmp(mri_get_string(ctx.input,key_buf),"float64")) {
    if (ctx.uFlag) mri_set_string(ctx.uDSet,key_buf,"float32");
    if (ctx.vtFlag) mri_set_string(ctx.vDSet,key_buf,"float32");
    if (ctx.wFlag) mri_set_string(ctx.wDSet,key_buf,"float32");
  }

  /* Do the eigenvalue solution */
  svdsolve_chunk(&ctx);

  /* Write and close data-sets */
  mri_close_dataset( ctx.input );
  if (ctx.uFlag) mri_close_dataset( ctx.uDSet );
  if (ctx.vtFlag) mri_close_dataset( ctx.vDSet );
  if (ctx.wFlag) mri_close_dataset( ctx.wDSet );

  if (verbose_flg) Message( "#      SVD deomposition complete.\n" );

  return 0;
}

