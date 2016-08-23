/************************************************************
 *                                                          *
 *  mri_matmult.c                                           *
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
 *  Original programming by Joel Welling, 9/02              *
 ************************************************************/
/*************************************************************

  DESCRIPTION OF MRI_MATMULT

  mri_matmult multiplies two Pgh MRI files as if they were matrices.
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

static char rcsid[] = "$Id: mri_matmult.c,v 1.17 2007/07/06 18:45:53 welling Exp $";

static char* progname;
static int verbose_flg= 0;
static int debug_flg= 0;

#define MAYBE_MAKE_HASHMARK(  n ) if (verbose_flg) makeHashMark(n)

static void makeHashMark( long n )
{
  static long i= 0;
  if (i==0) Message("      ");
  if ( (i==n-1) || (i+1)%60 == 0 ) {
    Message( "# %ld\n", i+1 );
    if ( i != n-1 ) Message( "      " );
  }
  else Message("#");
  i++;
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

static const char* safe_get_dims(MRI_Dataset* ds, const char* chunk)
{
  char key_buf[KEYBUF_SIZE];
  safe_copy(key_buf,chunk);
  safe_concat(key_buf,".dimensions");
  if (mri_has(ds,key_buf)) return mri_get_string(ds,key_buf);
  else Abort("%s: input missing tag %s!\n",progname,key_buf);
  return 0; /* not reached */
}

static void calc_sizes(MRI_Dataset* ds, char* this_chunk, char* dimstr, 
		       const char selected_dim,
		       long* fast_blocksize_out, long* slow_blocksize_out ) {
  /* This routine will fail if selected_dim is not in dimstr! */
  long fast_blocksize;
  long slow_blocksize;
  char* this_dim;

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

static void mult_general(MRI_Dataset* Left, MRI_Dataset* Right,
			 MRI_Dataset* Out, char* chunk,
			 long left_fast_blksize, long left_slow_blksize,
			 long right_slow_blksize, long summed_extent,
			 long long left_foreach_offset, 
			 long long right_foreach_offset, 
			 long long out_foreach_offset,
			 long long foreach_blksize)
{
  double* left_fast_buf= NULL;
  double* accum_fast_buf= NULL;
  double* summed_buf= NULL;
  double* summed_buf_all= NULL;
  long long left_base_offset;
  long long left_sum_offset;
  long long right_offset;
  long long out_offset;
  int left_fast;
  int right_slow;
  int left_slow;
  int summed;
  int preread= 0;

  if (debug_flg) fprintf(stderr,"General method used.\n");

  if (right_slow_blksize*summed_extent < MAX_AT_ONCE) {
    if (debug_flg) fprintf(stderr,"Loading entire right matrix!\n");
    preread= 1;
    summed_buf_all= mri_get_chunk(Right, chunk, 
				  summed_extent*right_slow_blksize, 
				  right_foreach_offset, MRI_DOUBLE);
    mri_retain_buffer(Right, summed_buf_all);
  }

  if (!(accum_fast_buf= (double*)malloc(left_fast_blksize*sizeof(double))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,left_fast_blksize*sizeof(double));

  /* OK, here we go. */
  left_base_offset= left_foreach_offset;
  out_offset= out_foreach_offset;
  for (left_slow=0; left_slow<left_slow_blksize; left_slow++) {
    right_offset= right_foreach_offset;
    for (right_slow=0; right_slow<right_slow_blksize; right_slow++) { 
      if (preread) summed_buf= summed_buf_all+right_slow*summed_extent;
      else {
	if (debug_flg)
	  fprintf(stderr,"Reading %d from right at %lld\n", 
		  summed_extent,right_offset);
	summed_buf= mri_get_chunk(Right, chunk, summed_extent, 
				  right_offset, MRI_DOUBLE);
	right_offset += summed_extent;
      }
      for (left_fast=0; left_fast<left_fast_blksize; left_fast++)
	accum_fast_buf[left_fast]= 0.0;
      left_sum_offset= 0;
      for (summed=0; summed<summed_extent; summed++) {
	if (debug_flg)
	  fprintf(stderr,"Read %d from left at %lld\n",left_fast_blksize,
		  left_base_offset+left_sum_offset);
	left_fast_buf= mri_get_chunk(Left, chunk, left_fast_blksize,
				     left_base_offset+left_sum_offset, 
				     MRI_DOUBLE);
	left_sum_offset += left_fast_blksize;
	for (left_fast=0; left_fast<left_fast_blksize; left_fast++)
	  accum_fast_buf[left_fast] += 
	    summed_buf[summed]*left_fast_buf[left_fast];
      }
      if (debug_flg)
	fprintf(stderr,"Writing %d at %lld\n",left_fast_blksize,out_offset);
      mri_set_chunk(Out, chunk, left_fast_blksize, 
		    out_offset, MRI_DOUBLE, accum_fast_buf);
      out_offset += left_fast_blksize;
    }
    left_base_offset += left_fast_blksize*summed_extent;
    MAYBE_MAKE_HASHMARK( left_slow_blksize*foreach_blksize );
  }

  free(accum_fast_buf);
  if (preread) mri_discard_buffer(Right, summed_buf_all);
}

static void mult_general_complex(MRI_Dataset* Left, MRI_Dataset* Right,
				 MRI_Dataset* Out, char* chunk,
				 long left_fast_blksize, 
				 long left_slow_blksize,
				 long right_slow_blksize, 
				 long summed_extent,
				 long long left_foreach_offset, 
				 long long right_foreach_offset, 
				 long long out_foreach_offset,
				 long long foreach_blksize)

{
  double* left_fast_buf= NULL;
  double* accum_fast_buf= NULL;
  double* summed_buf= NULL;
  double* summed_buf_all= NULL;
  long long left_base_offset;
  long long left_sum_offset;
  long long right_offset;
  long long out_offset;
  int left_fast;
  int right_slow;
  int left_slow;
  int summed;
  int preread= 0;

  if (debug_flg) fprintf(stderr,"General complex method used.\n");

  if (2*right_slow_blksize*summed_extent < MAX_AT_ONCE) {
    if (debug_flg) fprintf(stderr,"Loading entire right matrix!\n");
    preread= 1;
    summed_buf_all= mri_get_chunk(Right, chunk, 
				  2*summed_extent*right_slow_blksize, 
				  right_foreach_offset, MRI_DOUBLE);
    mri_retain_buffer(Right, summed_buf_all);
  }

  if (!(accum_fast_buf= (double*)malloc(left_fast_blksize*sizeof(double))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,left_fast_blksize*sizeof(double));

  /* OK, here we go. */
  left_base_offset= left_foreach_offset;
  out_offset= out_foreach_offset;  
  for (left_slow=0; left_slow<left_slow_blksize; left_slow++) {
    right_offset= right_foreach_offset;
    for (right_slow=0; right_slow<right_slow_blksize; right_slow++) { 
      if (preread) summed_buf= summed_buf_all+2*right_slow*summed_extent;
      else {
	if (debug_flg)
	  fprintf(stderr,"Reading %d from right at %lld\n", 
		  summed_extent,right_offset);
	summed_buf= mri_get_chunk(Right, chunk, 2*summed_extent, 
				  right_offset, MRI_DOUBLE);
	right_offset += 2*summed_extent;
      }
      for (left_fast=0; left_fast<left_fast_blksize; left_fast += 2)
	accum_fast_buf[left_fast+1]= accum_fast_buf[left_fast]= 0.0;
      left_sum_offset= 0;
      for (summed=0; summed<summed_extent; summed++) {
	double rR= summed_buf[2*summed];
	double rI= summed_buf[2*summed+1];
	if (debug_flg)
	  fprintf(stderr,"Read %d from left at %lld\n",left_fast_blksize,
		  left_base_offset+left_sum_offset);
	left_fast_buf= mri_get_chunk(Left, chunk, left_fast_blksize,
				     left_base_offset+left_sum_offset, 
				     MRI_DOUBLE);
	left_sum_offset += left_fast_blksize;
	for (left_fast=0; left_fast<left_fast_blksize; left_fast+=2) {
	  double lR= left_fast_buf[left_fast];
	  double lI= left_fast_buf[left_fast+1];
	  accum_fast_buf[left_fast] += rR*lR - rI*lI;
	  accum_fast_buf[left_fast+1] += rR*lI + rI*lR;
	}
      }
      if (debug_flg)
	fprintf(stderr,"Writing %d at %lld\n",left_fast_blksize,out_offset);
      mri_set_chunk(Out, chunk, left_fast_blksize, out_offset,
		    MRI_DOUBLE, accum_fast_buf);
      out_offset += left_fast_blksize;
    }
    left_base_offset += left_fast_blksize*summed_extent;
    MAYBE_MAKE_HASHMARK( left_slow_blksize*foreach_blksize );
  }

  free(accum_fast_buf);
  if (preread) mri_discard_buffer(Right, summed_buf_all);
}

static void mult_leftsmall(MRI_Dataset* Left, MRI_Dataset* Right,
			   MRI_Dataset* Out, char* chunk,
			   long left_fast_blksize, long left_slow_blksize,
			   long right_slow_blksize, long summed_extent,
			   long long left_foreach_offset, 
			   long long right_foreach_offset, 
			   long long out_foreach_offset,
			   long long foreach_blksize)
{
  double* left_fast_buf= NULL;
  double* accum_fast_buf= NULL;
  double* summed_buf= NULL;
  double* summed_buf_all= NULL;
  long long left_base_offset;
  long long right_offset;
  long long out_offset;
  int left_fast;
  int right_slow;
  int left_slow;
  int summed;
  int preread= 0;

  if (debug_flg) fprintf(stderr,"Leftsmall method used.\n");

  if (right_slow_blksize*summed_extent < MAX_AT_ONCE) {
    if (debug_flg) fprintf(stderr,"Loading entire right matrix!\n");
    preread= 1;
    summed_buf_all= mri_get_chunk(Right, chunk, 
				  summed_extent*right_slow_blksize, 
				  right_foreach_offset, MRI_DOUBLE);
    mri_retain_buffer(Right, summed_buf_all);
  }

  if (!(accum_fast_buf= 
	(double*)malloc(left_fast_blksize*sizeof(double))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,left_fast_blksize*sizeof(double));

  /* OK, here we go. */
  left_base_offset= left_foreach_offset;
  out_offset= out_foreach_offset;
  for (left_slow=0; left_slow<left_slow_blksize; left_slow++) {
    if (debug_flg)
      fprintf(stderr,"Reading %d from left at %lld\n",
	      left_fast_blksize*summed_extent, left_base_offset);
    left_fast_buf= mri_get_chunk(Left, chunk, 
				 left_fast_blksize*summed_extent,
				 left_base_offset, MRI_DOUBLE);
    right_offset= right_foreach_offset;
    for (right_slow=0; right_slow<right_slow_blksize; right_slow++) { 
      if (preread) summed_buf= summed_buf_all+right_slow*summed_extent;
      else {
	if (debug_flg)
	  fprintf(stderr,"Reading %d from right at %lld\n", 
		  summed_extent,right_offset);
	summed_buf= mri_get_chunk(Right, chunk, summed_extent, 
				  right_offset, MRI_DOUBLE);
	right_offset += summed_extent;
      }
      for (left_fast=0; left_fast<left_fast_blksize; left_fast++)
	accum_fast_buf[left_fast]= 0.0;
      for (summed=0; summed<summed_extent; summed++)
	for (left_fast=0; left_fast<left_fast_blksize; left_fast++)
	  accum_fast_buf[left_fast] += 
	    summed_buf[summed]
	    *left_fast_buf[left_fast_blksize*summed+left_fast];
      if (debug_flg)
	fprintf(stderr,"Writing %d at %lld\n",left_fast_blksize,out_offset);
      mri_set_chunk(Out, chunk, left_fast_blksize, out_offset,
		    MRI_DOUBLE, accum_fast_buf);
      out_offset += left_fast_blksize;
    }
    left_base_offset += left_fast_blksize*summed_extent;
    MAYBE_MAKE_HASHMARK( left_slow_blksize*foreach_blksize );
  }

  free(accum_fast_buf);
  if (preread) mri_discard_buffer(Right, summed_buf_all);
}

static void mult_leftsmall_complex(MRI_Dataset* Left, MRI_Dataset* Right,
				   MRI_Dataset* Out, char* chunk,
				   long left_fast_blksize, 
				   long left_slow_blksize,
				   long right_slow_blksize, 
				   long summed_extent,
				   long long left_foreach_offset, 
				   long long right_foreach_offset, 
				   long long out_foreach_offset,
				   long long foreach_blksize)
{
  double* left_fast_buf= NULL;
  double* accum_fast_buf= NULL;
  double* summed_buf= NULL;
  double* summed_buf_all= NULL;
  long long left_base_offset;
  long long right_offset;
  long long out_offset;
  int left_fast;
  int right_slow;
  int left_slow;
  int summed;
  int preread= 0;

  if (debug_flg) fprintf(stderr,"Leftsmall complex method used.\n");

  if (2*right_slow_blksize*summed_extent < MAX_AT_ONCE) {
    if (debug_flg) fprintf(stderr,"Loading entire right matrix!\n");
    preread= 1;
    summed_buf_all= mri_get_chunk(Right, chunk, 
				  2*summed_extent*right_slow_blksize, 
				  right_foreach_offset, MRI_DOUBLE);
    mri_retain_buffer(Right, summed_buf_all);
  }

  if (!(accum_fast_buf= 
	(double*)malloc(left_fast_blksize*sizeof(double))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,left_fast_blksize*sizeof(double));

  /* OK, here we go. */
  left_base_offset= left_foreach_offset;
  out_offset= out_foreach_offset;
  for (left_slow=0; left_slow<left_slow_blksize; left_slow++) {
    if (debug_flg)
      fprintf(stderr,"Reading %d from left at %lld\n",
	      left_fast_blksize*summed_extent, left_base_offset);
    left_fast_buf= mri_get_chunk(Left, chunk, 
				 left_fast_blksize*summed_extent,
				 left_base_offset, MRI_DOUBLE);
    right_offset= right_foreach_offset;
    for (right_slow=0; right_slow<right_slow_blksize; right_slow++) { 
      if (preread) summed_buf= summed_buf_all+2*right_slow*summed_extent;
      else {
	if (debug_flg)
	  fprintf(stderr,"Reading %d from right at %lld\n", 
		  summed_extent,right_offset);
	summed_buf= mri_get_chunk(Right, chunk, 2*summed_extent, 
				  right_offset, MRI_DOUBLE);
	right_offset += 2*summed_extent;
      }
      for (left_fast=0; left_fast<left_fast_blksize; left_fast+=2)
	accum_fast_buf[left_fast+1]= accum_fast_buf[left_fast]= 0.0;
      for (summed=0; summed<summed_extent; summed++) {
	double rR= summed_buf[2*summed];
	double rI= summed_buf[2*summed+1];
	for (left_fast=0; left_fast<left_fast_blksize; left_fast+=2) {
	  double lR= left_fast_buf[left_fast_blksize*summed + left_fast];
	  double lI= left_fast_buf[left_fast_blksize*summed + left_fast+1];
	  accum_fast_buf[left_fast] += rR*lR - rI*lI;
	  accum_fast_buf[left_fast+1] += rR*lI + rI*lR;
	}
      }
      if (debug_flg)
	fprintf(stderr,"Writing %d at %lld\n",left_fast_blksize,out_offset);
      mri_set_chunk(Out, chunk, left_fast_blksize, out_offset,
		    MRI_DOUBLE, accum_fast_buf);
      out_offset += left_fast_blksize;
    }
    left_base_offset += left_fast_blksize*summed_extent;
    MAYBE_MAKE_HASHMARK( left_slow_blksize*foreach_blksize );
  }

  free(accum_fast_buf);
  if (preread) mri_discard_buffer(Right, summed_buf_all);
}

static void mult_leftsmall_accumsmall(MRI_Dataset* Left, MRI_Dataset* Right,
				      MRI_Dataset* Out, char* chunk,
				      long left_fast_blksize, 
				      long left_slow_blksize,
				      long right_slow_blksize, 
				      long summed_extent,
				      long long left_foreach_offset, 
				      long long right_foreach_offset, 
				      long long out_foreach_offset,
				      long long foreach_blksize)
{
  double* left_fast_buf= NULL;
  double* accum_fast_buf= NULL;
  double* summed_buf= NULL;
  double* summed_buf_all= NULL;
  double one= 1.0;
  double zero= 0.0;
  int int_one= 1;
  long long left_base_offset;
  long long right_offset;
  long long out_offset;
  long long accum_offset;
  int left_fast;
  int right_slow;
  int left_slow;
  int summed;
  int preread= 0;
  int rsb= (int)right_slow_blksize;
  int lfb= (int)left_fast_blksize;
  int se= (int)summed_extent;

  if (debug_flg) fprintf(stderr,"Leftsmall_Accumsmall method used.\n");
  
  if (right_slow_blksize*summed_extent < MAX_AT_ONCE) {
    if (debug_flg) fprintf(stderr,"Loading entire right matrix!\n");
    preread= 1;
    summed_buf_all= mri_get_chunk(Right, chunk, 
				  summed_extent*right_slow_blksize, 
				  right_foreach_offset, MRI_DOUBLE);
    mri_retain_buffer(Right, summed_buf_all);
  }
  else preread= 0;

  if (!(accum_fast_buf= 
	(double*)malloc(left_fast_blksize*right_slow_blksize*sizeof(double))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,left_fast_blksize*right_slow_blksize*sizeof(double));

  /* OK, here we go. */
  left_base_offset= left_foreach_offset;
  out_offset= out_foreach_offset;
  for (left_slow=0; left_slow<left_slow_blksize; left_slow++) {
    if (debug_flg)
      fprintf(stderr,"Reading %d from left at %lld\n",
	      left_fast_blksize*summed_extent, left_base_offset);
    left_fast_buf= mri_get_chunk(Left, chunk, 
				 left_fast_blksize*summed_extent,
				 left_base_offset, MRI_DOUBLE);
    right_offset= right_foreach_offset;
    accum_offset= 0;
    if (preread) {
      DGEMM( "n", "n", 
	     &lfb, &rsb, &se, 
	     &one, 
	     left_fast_buf, &lfb, 
	     summed_buf_all, &se, 
	     &zero, accum_fast_buf, &lfb );
    }
    else {
      for (right_slow=0; right_slow<right_slow_blksize; right_slow++) { 
	if (debug_flg)
	  fprintf(stderr,"Reading %d from right at %lld\n", 
		  summed_extent,right_offset);
	summed_buf= mri_get_chunk(Right, chunk, summed_extent, 
				  right_offset, MRI_DOUBLE);
	DGEMV( "n", &lfb, &se, 
	       &one, 
	       left_fast_buf, &lfb, 
	       summed_buf, &int_one,
	       &zero,
	       accum_fast_buf+accum_offset, &int_one );
	accum_offset += left_fast_blksize;
	right_offset += summed_extent;
      }
    }
    if (debug_flg) 
      fprintf(stderr,"Writing block %ld to %lld (size %ld)\n",
	      left_slow,out_offset,
	      left_fast_blksize*right_slow_blksize);
    mri_set_chunk(Out, chunk, left_fast_blksize*right_slow_blksize, 
		  out_offset, MRI_DOUBLE, accum_fast_buf);
    left_base_offset += left_fast_blksize*summed_extent;
    out_offset += left_fast_blksize*right_slow_blksize;
    MAYBE_MAKE_HASHMARK( left_slow_blksize*foreach_blksize );
  }

  free(accum_fast_buf);
  if (preread) mri_discard_buffer(Right, summed_buf_all);
}

static void mult_leftsmall_accumsmall_complex(MRI_Dataset* Left, 
					      MRI_Dataset* Right,
					      MRI_Dataset* Out, char* chunk,
					      long left_fast_blksize, 
					      long left_slow_blksize,
					      long right_slow_blksize, 
					      long summed_extent,
					      long long left_foreach_offset, 
					      long long right_foreach_offset, 
					      long long out_foreach_offset,
					      long long foreach_blksize)
{
  double* left_fast_buf= NULL;
  double* accum_fast_buf= NULL;
  double* summed_buf= NULL;
  double* summed_buf_all= NULL;
  double c_one[2]= {1.0,0.0};
  double c_zero[2]= {0.0,0.0};
  int int_one= 1;
  long long left_base_offset;
  long long right_offset;
  long long out_offset;
  long long accum_offset;
  int left_fast;
  int right_slow;
  int left_slow;
  int summed;
  int preread= 0;
  int rsb= (int)right_slow_blksize;
  int hlfb= (int)(left_fast_blksize/2); /* allow for complex */
  int se= (int)summed_extent;

  if (debug_flg) 
    fprintf(stderr,"Leftsmall_Accumsmall complex method used.\n");

  if (2*right_slow_blksize*summed_extent < MAX_AT_ONCE) {
    if (debug_flg) fprintf(stderr,"Loading entire right matrix!\n");
    preread= 1;
    summed_buf_all= mri_get_chunk(Right, chunk, 
				  2*summed_extent*right_slow_blksize, 
				  right_foreach_offset, MRI_DOUBLE);
    mri_retain_buffer(Right, summed_buf_all);
  }
  else preread= 0;

  if (!(accum_fast_buf= 
	(double*)malloc(left_fast_blksize*right_slow_blksize*sizeof(double))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,left_fast_blksize*right_slow_blksize*sizeof(double));

  /* OK, here we go. */
  left_base_offset= left_foreach_offset;
  out_offset= out_foreach_offset;
  for (left_slow=0; left_slow<left_slow_blksize; left_slow++) {
    if (debug_flg)
      fprintf(stderr,"Reading %d from left at %lld\n",
	      left_fast_blksize*summed_extent, left_base_offset);
    left_fast_buf= mri_get_chunk(Left, chunk, 
				 left_fast_blksize*summed_extent,
				 left_base_offset, MRI_DOUBLE);
    right_offset= right_foreach_offset;
    accum_offset= 0;
    if (preread) {
      ZGEMM( "n", "n", 
	     &hlfb, &rsb, &se, 
	     c_one, 
	     left_fast_buf, &hlfb, 
	     summed_buf_all, &se, 
	     c_zero, accum_fast_buf, &hlfb );
    }
    else {
      for (right_slow=0; right_slow<right_slow_blksize; right_slow++) { 
	if (debug_flg)
	  fprintf(stderr,"Reading %d from right at %lld\n", 
		  2*summed_extent,right_offset);
	summed_buf= mri_get_chunk(Right, chunk, 2*summed_extent, 
				  right_offset, MRI_DOUBLE);
	ZGEMV( "n", &hlfb, &se, 
	       c_one, 
	       left_fast_buf, &hlfb, 
	       summed_buf, &int_one,
	       c_zero,
	       accum_fast_buf+accum_offset, &int_one );
	accum_offset += left_fast_blksize;
	right_offset += 2*summed_extent;
      }
    }
    if (debug_flg) 
      fprintf(stderr,"Writing block %ld to %lld (size %ld)\n",
	      left_slow,out_offset,
	      left_fast_blksize*right_slow_blksize);
    mri_set_chunk(Out, chunk, left_fast_blksize*right_slow_blksize, 
		  out_offset, MRI_DOUBLE, accum_fast_buf);
    left_base_offset += left_fast_blksize*summed_extent;
    out_offset += left_fast_blksize*right_slow_blksize;
    MAYBE_MAKE_HASHMARK( left_slow_blksize*foreach_blksize );
  }

  free(accum_fast_buf);
  if (preread) mri_discard_buffer(Right, summed_buf_all);
}

static void mult_chunk(MRI_Dataset* Left, MRI_Dataset* Right, 
		       MRI_Dataset* Out, char* chunk, const int summed_dim,
		       const char* foreach_dims)
{
  char* dimstr1_orig= strdup(safe_get_dims(Left,chunk));
  char* dimstr1= dimstr1_orig;
  char* dimstr2_orig= strdup(safe_get_dims(Right,chunk));
  char* dimstr2= dimstr2_orig;
  long summed_extent;
  long left_fast_blksize;
  long left_slow_blksize;
  long right_fast_blksize;
  long right_slow_blksize;
  long long foreach_blksize= 1;
  long long left_foreach_offset= 0;
  long long right_foreach_offset= 0;
  long long out_foreach_offset= 0;
  long long foreach_loop;
  
  if (*foreach_dims) {
    char* here= strstr(dimstr1,foreach_dims);
    *here= '\0';
    here= strstr(dimstr2,foreach_dims);
    *here= '\0';
    here= (char*)foreach_dims;
    while (*here) foreach_blksize *= safe_get_extent(Left, chunk, *here++);
  }

  summed_extent= safe_get_extent(Left, chunk, summed_dim);
  calc_sizes(Left, chunk, dimstr1, summed_dim, 
	     &left_fast_blksize, &left_slow_blksize);
  calc_sizes(Right, chunk, dimstr2, summed_dim, 
	     &right_fast_blksize, &right_slow_blksize);
  if (debug_flg) {
    fprintf(stderr,"Left string %s, block sizes: %ld %ld %ld\n",
	    dimstr1,left_fast_blksize, summed_extent, left_slow_blksize);
    fprintf(stderr,"Right string %s, block sizes: %ld %ld %ld\n",
	    dimstr2,right_fast_blksize, summed_extent, right_slow_blksize);
    fprintf(stderr,"Foreach string <%s>, block size %lld\n", 
	    foreach_dims, foreach_blksize );
    fprintf(stderr,"Output dims %s\n",mri_get_string(Out,"images.dimensions"));
  }
  if (right_fast_blksize != 1)
    Abort("%s: internal error: right fast block is not 1!\n",progname);
  if (verbose_flg) 
    Message("Counting out %lld blocks:\n",left_slow_blksize*foreach_blksize);
  
  for (foreach_loop=0; foreach_loop<foreach_blksize; foreach_loop++) {
    if (left_fast_blksize*summed_extent <= MAX_AT_ONCE) {
      if (left_fast_blksize*right_slow_blksize <= MAX_AT_ONCE) {
	mult_leftsmall_accumsmall(Left, Right, Out, chunk,
				  left_fast_blksize, left_slow_blksize,
				  right_slow_blksize, summed_extent,
				  left_foreach_offset, right_foreach_offset, 
				  out_foreach_offset, foreach_blksize);
      }
      else {
	mult_leftsmall(Left, Right, Out, chunk,
		       left_fast_blksize, left_slow_blksize,
		       right_slow_blksize, summed_extent,
		       left_foreach_offset, right_foreach_offset, 
		       out_foreach_offset, foreach_blksize);
      }
    }
    else {
      mult_general(Left, Right, Out, chunk,
		   left_fast_blksize, left_slow_blksize,
		   right_slow_blksize, summed_extent,
		   left_foreach_offset, right_foreach_offset, 
		   out_foreach_offset, foreach_blksize);
    }
    left_foreach_offset += 
      left_fast_blksize*summed_extent*left_slow_blksize;
    right_foreach_offset += 
      right_fast_blksize*summed_extent*right_slow_blksize;
    out_foreach_offset += 
      left_fast_blksize*right_slow_blksize*left_slow_blksize;
  }

  free(dimstr1_orig);
  free(dimstr2_orig);
}

static void mult_chunk_complex(MRI_Dataset* Left, MRI_Dataset* Right, 
			       MRI_Dataset* Out, char* chunk, 
			       const int summed_dim, const char* foreach_dims)
{
  char* dimstr1_orig= strdup(safe_get_dims(Left,chunk));
  char* dimstr1= dimstr1_orig;
  char* dimstr2_orig= strdup(safe_get_dims(Right,chunk));
  char* dimstr2= dimstr2_orig;
  long summed_extent;
  long left_fast_blksize;
  long left_slow_blksize;
  long right_fast_blksize;
  long right_slow_blksize;
  long long foreach_blksize= 1;
  long long left_foreach_offset= 0;
  long long right_foreach_offset= 0;
  long long out_foreach_offset= 0;
  long long foreach_loop;
  
  if (*foreach_dims) {
    char* here= strstr(dimstr1,foreach_dims);
    *here= '\0';
    here= strstr(dimstr2,foreach_dims);
    *here= '\0';
    here= (char*)foreach_dims;
    while (*here) foreach_blksize *= safe_get_extent(Left, chunk, *here++);
  }
  
  summed_extent= safe_get_extent(Left, chunk, summed_dim);
  calc_sizes(Left, chunk, dimstr1, summed_dim, 
	     &left_fast_blksize, &left_slow_blksize);
  calc_sizes(Right, chunk, dimstr2, summed_dim, 
	     &right_fast_blksize, &right_slow_blksize);
  if (debug_flg) {
    fprintf(stderr,"Left string %s, block sizes: %ld %ld %ld\n",
	    dimstr1,left_fast_blksize, summed_extent, left_slow_blksize);
    fprintf(stderr,"Right string %s, block sizes: %ld %ld %ld\n",
	    dimstr2,right_fast_blksize, summed_extent, right_slow_blksize);
    fprintf(stderr,"Foreach string %s, block size %lld\n", 
	    foreach_dims,foreach_blksize );
    fprintf(stderr,"Output dims %s\n",mri_get_string(Out,"images.dimensions"));
  }
  if (right_fast_blksize != 2)
    Abort("%s: internal error: right fast block is not 2!\n",progname);
  if (verbose_flg) 
    Message("Counting out %d blocks:\n",left_slow_blksize*foreach_blksize);
  
  for (foreach_loop=0; foreach_loop<foreach_blksize; foreach_loop++){  
    if (left_fast_blksize*summed_extent <= MAX_AT_ONCE) {
      if (left_fast_blksize*right_slow_blksize <= MAX_AT_ONCE) {
	mult_leftsmall_accumsmall_complex(Left, Right, Out, chunk,
					  left_fast_blksize, left_slow_blksize,
					  right_slow_blksize, summed_extent,
					  left_foreach_offset, 
					  right_foreach_offset, 
					  out_foreach_offset, 
					  foreach_blksize);

      }
      else {
	mult_leftsmall_complex(Left, Right, Out, chunk,
			       left_fast_blksize, left_slow_blksize,
			       right_slow_blksize, summed_extent,
			       left_foreach_offset, right_foreach_offset, 
			       out_foreach_offset, foreach_blksize);
      }
    }
    else {
      mult_general_complex(Left, Right, Out, chunk,
			   left_fast_blksize, left_slow_blksize,
			   right_slow_blksize, summed_extent,
			   left_foreach_offset, right_foreach_offset,
			   out_foreach_offset, foreach_blksize);
    }
    left_foreach_offset += 
      left_fast_blksize*summed_extent*left_slow_blksize;
    right_foreach_offset += 
      right_fast_blksize*summed_extent*right_slow_blksize;
    out_foreach_offset += 
      left_fast_blksize*right_slow_blksize*left_slow_blksize;
  }

  free(dimstr1_orig);
  free(dimstr2_orig);
}

static int chunk_check( MRI_Dataset* ds, const char* chunk )
{
  return( mri_has(ds, chunk) && !strcmp(mri_get_string(ds,chunk),"[chunk]") );
}

static int structure_check( MRI_Dataset* Left, MRI_Dataset* Right,
			    const char* chunk, const int complex,
			    int* summed_dim_ptr, char** foreach_dims_ptr)
{
  char* dimstr1_orig= strdup(safe_get_dims(Left,chunk));
  char* dimstr1= dimstr1_orig;
  char* dimstr2_orig= strdup(safe_get_dims(Right,chunk));
  char* dimstr2= dimstr2_orig;
  int summed_dim= '\0';
  char* here;
  char* there;
  char* foreach_dims;

  /* If this is to be a complex matrix multiplication, both
   * inputs should have first input dimensions v with extent 2.
   */
  if (complex) {
    if (dimstr1[0] != 'v' || safe_get_extent(Left,chunk,'v')!=2 )
      Abort("%s: complex multiplication requested, but first input is not complex!\n",
	    progname);
    if (dimstr2[0] != 'v' || safe_get_extent(Right,chunk,'v')!=2 )
      Abort("%s: complex multiplication requested, but second input is not complex!\n",
	    progname);
    dimstr1++;
    dimstr2++;
  }

  /* Find any looped-over dimensions and clip them from the dim strings */
  here= dimstr1 + strlen(dimstr1) - 1;
  there= dimstr2 + strlen(dimstr2) - 1;
  while (*here==*there && here>dimstr1 && there>dimstr2) {
    if (safe_get_extent(Left,chunk,*here)
	!= safe_get_extent(Right,chunk,*here))
      Abort("%s: extents of looped-over dimension %c do not match!\n",
	    progname, *here);
    here--;
    there--;
  }
  foreach_dims= strdup(here+1);
  *(here+1)= *(there+1)= '\0';

  /* In remaining dim string, there must be exactly one dimension in common */
  for (here=dimstr1; *here; here++) {
    if ((there= strchr(dimstr2, *here)) != NULL) {
      if (summed_dim== '\0') summed_dim= *here;
      else Abort("%s: more than one duplicated dimension in input chunks!\n",
		 progname);
    }
  }

  /* It must be the left-most dimension of the right-hand dataset */
  if (dimstr2[0] != summed_dim)
    Abort("%s: Summed dimension %c must be left-most in input dataset #2!\n",
	  progname, summed_dim);

  /* This dimension must have the same extent in both factors */
  if (safe_get_extent(Left,chunk,summed_dim) 
      != safe_get_extent(Right,chunk,summed_dim))
    Abort("%s: The inputs have different extents for summed dimension %c!\n",
	  progname, summed_dim);

  *summed_dim_ptr= summed_dim;
  *foreach_dims_ptr= foreach_dims;
  free(dimstr1_orig);
  free(dimstr2_orig);
  return 1;
}

static void restructure_dims( MRI_Dataset* Output, const char* chunk,
			      MRI_Dataset* Factor1, MRI_Dataset* Factor2,
			      const int complex, const char* foreach_dims )
{
  /* This routine changes the dimensions of the output to be those
   * appropriate for the product of the inputs.  structure_check()
   * has already verified that the input dataset have the necessary
   * structure for these steps to work.
   */
  char key_buf[KEYBUF_SIZE];
  char dim_buf[KEYBUF_SIZE];
  char* dimstr1_orig= strdup(safe_get_dims(Factor1,chunk));
  char* dimstr1= dimstr1_orig;
  char* dimstr2_orig= strdup(safe_get_dims(Factor2,chunk));
  char* dimstr2= dimstr2_orig;
  int summed_dim; 
  char* work1= strdup(dimstr1);
  char* here= NULL;
  char* there= NULL;
  char* end;

  if (complex) { /* skip the v dimension */
    dimstr1++;
    dimstr2++;
  }

  summed_dim= dimstr2[0];
  here= strchr(work1,summed_dim);
  if (foreach_dims[0]) {
    there= strstr(dimstr1+1,foreach_dims);
    *there= '\0';
    there= strstr(dimstr2+1,foreach_dims);
    *there= '\0';
  }
  /* here now points to the live part of dimstr2, excluding the summed dim,
   * the foreach dims if any, and the v dim if any.
   */

  if (verbose_flg) {
    if (complex) {
      Message("Multiplying complex input chunk <%s>, dims v%s%s and v%s%s.\n",
	      chunk,dimstr1,foreach_dims,dimstr2,foreach_dims);
    }
    else {
      Message("Multiplying input chunk <%s>, dims %s%s and %s%s.\n",
	      chunk,dimstr1,foreach_dims,dimstr2,foreach_dims);
    }
    Message("Summing over dimension %c, looping over dimensions <%s>\n",
	    summed_dim, foreach_dims);
  }

  /* Break the left dimension string */
  *here= '\0';
  here++; /* now points to left slow dimensions + foreach dims */

  /* Define the new dimension string */
  safe_copy(dim_buf,work1);
  safe_concat(dim_buf,dimstr2 + 1);
  safe_concat(dim_buf,here);
  safe_copy(key_buf,chunk);
  safe_concat(key_buf,".dimensions");
  mri_set_string(Output, key_buf,dim_buf);
  if (verbose_flg) 
    Message("Output dimension string is <%s>\n",dim_buf);
  
  /* Delete the summed-over extent */
  safe_copy(key_buf,chunk);
  safe_concat(key_buf,".extent.");
  if (strlen(key_buf)>=KEYBUF_SIZE-2) 
    Abort("%s: this chunk name is ridiculously long!\n",progname);
  end= key_buf + strlen(key_buf);
  sprintf(end,"%c",summed_dim);
  mri_remove(Output,key_buf);

  /* Define the undefined extents */
  for (here= dimstr2+1; *here; here++) {
    sprintf(end,"%c",*here);
    mri_set_int(Output,key_buf,safe_get_extent(Factor2,chunk,*here));
  }
  
  free(work1);
  free(dimstr1_orig);
  free(dimstr2_orig);
}

int main( int argc, char* argv[] ) 
{
  char fname1[512], fname2[512], outfile[512];
  int summed_dim;
  char* foreach_dims= NULL;
  MRI_Dataset *Factor1 = NULL, *Factor2 = NULL, *Output = NULL;
  char chunk[KEYBUF_SIZE];
  char key_buf[KEYBUF_SIZE];
  int complex_flg= 0;

  progname= argv[0];

  /* Print version number */
  if (verbose_flg) Message( "# %s\n", rcsid );

  /* Check to see if help was requested */
  if( ( argc > 1 ) && !strcmp( argv[1], "-help" ) )
    {
      if( argc == 2 )
	Help( "selecttopic" );
      else
	Help( argv[2] );
    }

  /*** Parse command line ***/

  cl_scan( argc, argv );

  verbose_flg= cl_present("verbose|ver|v");
  debug_flg= cl_present("debug|deb");
  complex_flg= cl_present("complex|cpx");
  cl_get("chunk|chu|c", "%option %s[%]","images",chunk);
  if (!cl_get("outfile|out", "%option %s", outfile)) {
    fprintf(stderr,"%s: Output file name not given.\n",progname);
    Help( "usage" );
    exit(-1);
  }
  if (!cl_get("", "%s", fname1)) {
    fprintf(stderr,"%s: Factor1 file name not given.\n",progname);
    Help( "usage" );
    exit(-1);
  }
  if (!cl_get("", "%s", fname2)) {
    fprintf(stderr,"%s: Factor2 file name not given.\n",progname);
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
  if( !strcmp( fname1, outfile ) )
    Abort( "%s: Factor1 and output files must be distinct.", progname );
  if( !strcmp( fname2, outfile ) )
    Abort( "%s: Factor2 and output files must be distinct.", progname );
  if( !strcmp( fname1, fname2 ) )
    Abort( "%s: Factor1 and Factor2 files must be distinct.", progname );
  Factor1 = mri_open_dataset( fname1, MRI_READ );
  Factor2 = mri_open_dataset( fname2, MRI_READ );

  /* Will they work with the program? */
  if (!chunk_check(Factor1, chunk))
    Abort("%s: %s has no chunk %s!\n",progname,fname1,chunk);
  if (!chunk_check(Factor2, chunk))
    Abort("%s: %s has no chunk %s!\n",progname,fname2,chunk);
  if (!structure_check(Factor1, Factor2, chunk, complex_flg, &summed_dim,
		       &foreach_dims)
      || summed_dim=='\0' || foreach_dims==NULL)
    Abort("%s: Datasets %s and %s cannot be multiplied.\n",
	  progname,fname1,fname2);

  /* Open output dataset.  We'll use float32 as the datatype unless
   * it's already double precision.
   */
  Output = mri_copy_dataset( outfile, Factor1 );
  hist_add_cl( Output, argc, argv );
  safe_copy(key_buf,chunk);
  safe_concat(key_buf,".datatype");
  if (strcmp(mri_get_string(Output,key_buf),"float64"))
    mri_set_string(Output,key_buf,"float32");
  restructure_dims( Output, chunk, Factor1, Factor2, complex_flg, 
		    foreach_dims );

  /* Do the multiplication */
  if (complex_flg)
    mult_chunk_complex(Factor1, Factor2, Output, chunk, summed_dim,
		       foreach_dims);
  else
    mult_chunk(Factor1, Factor2, Output, chunk, summed_dim, foreach_dims);

  /* Write and close data-sets */
  mri_close_dataset( Factor1 );
  mri_close_dataset( Factor2 );
  mri_close_dataset( Output );
  
  if (verbose_flg) Message( "#      Multiplication complete.\n" );

  return 0;
}

