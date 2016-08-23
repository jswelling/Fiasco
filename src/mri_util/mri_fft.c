/************************************************************
 *                                                          *
 *  mri_fft.c                                     *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1999 Department of Statistics,         *
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
 *  Original programming by Joel Welling, 6/99              *
 ************************************************************/
/*************************************************************

  DESCRIPTION OF MRI_FFT

  mri_fft takes a pgh MRI dataset of any type, and outputs
  a dataset of the same type containing a ffted version of 
  the data.

**************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "misc.h"
#include "stdcrg.h"

#define KEYBUF_SIZE 512

static char rcsid[] = "$Id: mri_fft.c,v 1.15 2004/09/10 00:14:03 welling Exp $";

typedef enum { 
  RSLT_COMPLEX, RSLT_MODULUS, RSLT_PHASE, RSLT_SQMOD, RSLT_REAL, RSLT_IMAG } ResultType;

static MRI_Dataset *Input = NULL, *Output = NULL;
static char selected_dim[512]= "t";
static int data_changed= 0;
static char* progname;
static int verbose_flg= 0;
static long fft_sign= 1; /* implies reverse transform */
static ResultType resultType= RSLT_MODULUS;

static void safe_copy(char* str1, char* str2) {
  strncpy(str1, str2, KEYBUF_SIZE);
  str1[KEYBUF_SIZE-1]= '\0';
}

static void safe_concat(char* str1, char* str2) {
  strncat(str1, str2, (KEYBUF_SIZE-strlen(str1))-1);
}

static int safe_get_extent(MRI_Dataset* ds, char* chunk, char* dim)
{
  char key_buf[KEYBUF_SIZE];
  char dim_buf[4];
  dim_buf[0]= *dim;
  dim_buf[1]= '\0';
  safe_copy(key_buf,chunk);
  safe_concat(key_buf,".extent.");
  safe_concat(key_buf,dim_buf);
  if (mri_has(Input,key_buf)) return mri_get_int(Input,key_buf);
  else Abort("%s: input missing tag %s!\n",progname,key_buf);
  return 0; /* not reached */
}

static void calc_sizes(char* this_chunk, char* dimstr, 
		       long long* fast_blocksize_out, 
		       long long* slow_blocksize_out ) {
  /* This routine will fail if selected_dim is not in dimstr! */
  long long fast_blocksize;
  long long slow_blocksize;
  char* this_dim;

  fast_blocksize= 1;

  this_dim= dimstr;
  while (*this_dim != *selected_dim) 
    fast_blocksize *= safe_get_extent(Input,this_chunk,this_dim++);

  this_dim += strlen(selected_dim); /* step over selected dim */

  slow_blocksize= 1;
  while (*this_dim)
    slow_blocksize *= safe_get_extent(Input, this_chunk, this_dim++);

  *fast_blocksize_out= fast_blocksize;
  *slow_blocksize_out= slow_blocksize;
}

static void scalarize_chunk(MRI_Dataset* ds, char* chunk)
{
  char key_buf[KEYBUF_SIZE];

  /* We assume the chunk to be valid and complex */
  safe_copy(key_buf,chunk);
  safe_concat(key_buf,".extent.v");
  mri_set_int(ds,key_buf,1);
}

static void complexify_chunk(MRI_Dataset* ds, char* chunk)
{
  char key_buf[KEYBUF_SIZE];
  char new_dimstr[KEYBUF_SIZE];
  char* dimstr;

  /* We assume the chunk to be valid */
  safe_copy(key_buf, chunk);
  safe_concat(key_buf, ".dimensions");
  dimstr= mri_get_string(ds,key_buf);
  if (*dimstr != 'v') {
    safe_copy(new_dimstr,"v");
    safe_concat(new_dimstr,dimstr);
    mri_set_string(ds,key_buf,new_dimstr);
  }
  safe_copy(key_buf, chunk);
  safe_concat(key_buf, ".extent.v");
  mri_set_int(ds,key_buf,2);
  
}

static void transfer_data_scalar(char* this_chunk, 
				 long long fast_blksize, 
				 long long slow_blksize, 
				 long selected_extent_fast,
				 long selected_extent_slow)
{
  long long in_offset= 0;
  long long out_offset= 0;
  long long islow;
  long i;
  long collective_blksize;
  long oblock_size;
  FComplex* buf= NULL;
  float* oblock= NULL;

  collective_blksize= 
    (long)(fast_blksize * selected_extent_fast * selected_extent_slow);
  oblock_size= collective_blksize;
  if (resultType==RSLT_COMPLEX) oblock_size *= 2;

  if (!(buf=(FComplex*)malloc(collective_blksize*sizeof(FComplex))))
    Abort("%s: unable to allocate %ld bytes!\n",
	  progname, collective_blksize*sizeof(FComplex));
  if (!(oblock=(float*)malloc(oblock_size*sizeof(float))))
    Abort("%s: unable to allocate %ld bytes!\n",
	  progname, oblock_size*sizeof(float));

  for (islow=0; islow<slow_blksize; islow++) {
    float* block= mri_get_chunk(Input, this_chunk, collective_blksize,
				in_offset, MRI_FLOAT);
    for (i=0; i<collective_blksize; i++) {
      buf[i].real= block[i];
      buf[i].imag= 0.0;
    }
    if (selected_extent_slow==1)
      fft3d( buf, 1, selected_extent_fast, fast_blksize, fft_sign, "y" );
    else
      fft3d( buf, selected_extent_slow, selected_extent_fast, fast_blksize, 
	     fft_sign, "xy" );
    switch (resultType) {
    case RSLT_COMPLEX:
      for (i=0; i<collective_blksize; i++) {
	oblock[2*i]= buf[i].real;
	oblock[(2*i)+1]= buf[i].imag;
      }
      break;
    case RSLT_MODULUS:
      for (i=0; i<collective_blksize; i++) oblock[i]= Modulus( buf[i] );
      break;
    case RSLT_PHASE:
      for (i=0; i<collective_blksize; i++) oblock[i]= Phase( buf[i] );
      break;
    case RSLT_SQMOD:
      for (i=0; i<collective_blksize; i++) 
	oblock[i]= buf[i].real*buf[i].real + buf[i].imag*buf[i].imag;
      break;
    case RSLT_REAL:
      for (i=0; i<collective_blksize; i++) oblock[i]= buf[i].real;
      break;
    case RSLT_IMAG:
      for (i=0; i<collective_blksize; i++) oblock[i]= buf[i].imag;
      break;
    }
    mri_set_chunk( Output, this_chunk, oblock_size, out_offset,
		   MRI_FLOAT, oblock );
    if (verbose_flg) 
      fprintf(stderr,"block: %ld at %lld -> %ld at %lld\n",
	      collective_blksize, in_offset,
	      oblock_size, out_offset);
    in_offset += collective_blksize;
    out_offset += oblock_size;
  }

  free(buf);
  free(oblock);  
}

static void transfer_data_complex(char* this_chunk, 
				  long long fast_blksize, 
				  long long slow_blksize, 
				  long selected_extent_fast,
				  long selected_extent_slow)
{
  long long in_offset= 0;
  long long out_offset= 0;
  long long islow;
  long i;
  long collective_blksize;
  long oblock_size;
  FComplex* buf= NULL;
  float* oblock= NULL;

  collective_blksize= 
    (long)(fast_blksize * selected_extent_fast * selected_extent_slow);
  oblock_size= collective_blksize;
  if (resultType==RSLT_COMPLEX) oblock_size *= 2;

  if (!(buf=(FComplex*)malloc(collective_blksize*sizeof(FComplex))))
    Abort("%s: unable to allocate %ld bytes!\n",
	  progname, collective_blksize*sizeof(FComplex));
  if (!(oblock=(float*)malloc(oblock_size*sizeof(float))))
    Abort("%s: unable to allocate %ld bytes!\n",
	  progname, oblock_size*sizeof(float));

  for (islow=0; islow<slow_blksize; islow++) {
    float* block= mri_get_chunk(Input, this_chunk, 2*collective_blksize,
				in_offset, MRI_FLOAT);
    for (i=0; i<collective_blksize; i++) {
      buf[i].real= block[2*i];
      buf[i].imag= block[(2*i) + 1];
    }
    if (selected_extent_slow==1)
      fft3d( buf, 1, selected_extent_fast, fast_blksize, fft_sign, "y" );
    else
      fft3d( buf, selected_extent_slow, selected_extent_fast, fast_blksize, 
	     fft_sign, "xy" );
    switch (resultType) {
    case RSLT_COMPLEX:
      for (i=0; i<collective_blksize; i++) {
	oblock[2*i]= buf[i].real;
	oblock[(2*i)+1]= buf[i].imag;
      }
      break;
    case RSLT_MODULUS:
      for (i=0; i<collective_blksize; i++) oblock[i]= Modulus( buf[i] );
      break;
    case RSLT_PHASE:
      for (i=0; i<collective_blksize; i++) oblock[i]= Phase( buf[i] );
      break;
    case RSLT_SQMOD:
      for (i=0; i<collective_blksize; i++) 
	oblock[i]= buf[i].real*buf[i].real + buf[i].imag*buf[i].imag;
      break;
    case RSLT_REAL:
      for (i=0; i<collective_blksize; i++) oblock[i]= buf[i].real;
      break;
    case RSLT_IMAG:
      for (i=0; i<collective_blksize; i++) oblock[i]= buf[i].imag;
      break;
    }
    mri_set_chunk( Output, this_chunk, oblock_size, out_offset,
		   MRI_FLOAT, oblock );
    if (verbose_flg) 
      fprintf(stderr,"block: %ld at %lld -> %ld at %lld\n",
	      2*collective_blksize, in_offset,
	      oblock_size, out_offset);
    in_offset += 2*collective_blksize;
    out_offset += oblock_size;
  }

  free(buf);
  free(oblock);
}

static void transfer_data_3d_scalar(char* this_chunk, 
				    long long fast_blksize, 
				    long long slow_blksize, 
				    long selected_extent_fast,
				    long selected_extent_med,
				    long selected_extent_slow)
{
  long long in_offset= 0;
  long long out_offset= 0;
  long ifast;
  long long islow;
  long i;
  long collective_blksize;
  long fft_blksize;
  long oblock_size;
  FComplex* buf= NULL;
  float* oblock= NULL;

  fft_blksize= 
    selected_extent_fast * selected_extent_med * selected_extent_slow;
  collective_blksize= (long)(fast_blksize * fft_blksize);
  oblock_size= collective_blksize;
  if (resultType==RSLT_COMPLEX) oblock_size *= 2;

  if (!(buf=(FComplex*)malloc(fft_blksize*sizeof(FComplex))))
    Abort("%s: unable to allocate %ld bytes!\n",
	  progname, fft_blksize*sizeof(FComplex));
  if (!(oblock=(float*)malloc(oblock_size*sizeof(float))))
    Abort("%s: unable to allocate %ld bytes!\n",
	  progname, oblock_size*sizeof(float));

  for (islow=0; islow<slow_blksize; islow++) {
    float* block= mri_get_chunk(Input, this_chunk, collective_blksize,
				in_offset, MRI_FLOAT);
    for (ifast=0; ifast<fast_blksize; ifast++) {
      float* runner;
      FComplex* crunner= buf;
      for (runner= block+ifast; runner<block+collective_blksize;
	   runner += fast_blksize) {
	crunner->real= *runner;
	crunner->imag= 0.0;
	crunner++;
      }
      fft3d( buf, selected_extent_slow, selected_extent_med,
	     selected_extent_fast, fft_sign, "xyz" );
      switch (resultType) {
      case RSLT_COMPLEX:
	{
	  FComplex* orunner= ((FComplex*)oblock)+ifast;
	  for (crunner= buf; crunner<buf + fft_blksize; crunner++) {
	    orunner->real= crunner->real;
	    orunner->imag= crunner->imag;
	    orunner += fast_blksize;
	  }
	}
	break;
      case RSLT_MODULUS:
	{
	  float* orunner= oblock+ifast;
	  for (crunner= buf; crunner<buf + fft_blksize; crunner++) {
	    *orunner= Modulus( *crunner );
	    orunner += fast_blksize;
	  }
	}
	break;
      case RSLT_PHASE:
	{
	  float* orunner= oblock+ifast;
	  for (crunner= buf; crunner<buf + fft_blksize; crunner++) {
	    *orunner= Phase( *crunner );
	    orunner += fast_blksize;
	  }
	}
	break;
      case RSLT_SQMOD:
	{
	  float* orunner= oblock+ifast;
	  for (crunner= buf; crunner<buf + fft_blksize; crunner++) {
	    *orunner= 
	      crunner->real*crunner->real + crunner->imag*crunner->imag;
	    orunner += fast_blksize;
	  }
	}
	break;
      case RSLT_REAL:
	{
	  float* orunner= oblock+ifast;
	  for (crunner= buf; crunner<buf + fft_blksize; crunner++) {
	    *orunner= crunner->real;
	    orunner += fast_blksize;
	  }
	}
	break;
      case RSLT_IMAG:
	{
	  float* orunner= oblock+ifast;
	  for (crunner= buf; crunner<buf + fft_blksize; crunner++) {
	    *orunner= crunner->imag;
	    orunner += fast_blksize;
	  }
	}
	break;
      }
    }
    mri_set_chunk( Output, this_chunk, oblock_size, out_offset,
		   MRI_FLOAT, oblock );
    if (verbose_flg) 
      fprintf(stderr,"block: %ld at %lld -> %ld at %lld\n",
	      collective_blksize, in_offset,
	      oblock_size, out_offset);
    in_offset += collective_blksize;
    out_offset += oblock_size;
  }

  free(buf);
  free(oblock);  
}

static void transfer_data_3d_complex(char* this_chunk, 
				     long long fast_blksize, 
				     long long slow_blksize, 
				     long selected_extent_fast,
				     long selected_extent_med,
				     long selected_extent_slow)
{
  long long in_offset= 0;
  long long out_offset= 0;
  long ifast;
  long long islow;
  long i;
  long collective_blksize;
  long fft_blksize;
  long oblock_size;
  FComplex* buf= NULL;
  float* oblock= NULL;
  
  fft_blksize= 
    selected_extent_fast * selected_extent_med * selected_extent_slow;
  collective_blksize= (long)(fast_blksize * fft_blksize);
  oblock_size= collective_blksize;
  if (resultType==RSLT_COMPLEX) oblock_size *= 2;
  
  if (!(buf=(FComplex*)malloc(fft_blksize*sizeof(FComplex))))
    Abort("%s: unable to allocate %ld bytes!\n",
	  progname, fft_blksize*sizeof(FComplex));
  if (!(oblock=(float*)malloc(oblock_size*sizeof(float))))
    Abort("%s: unable to allocate %ld bytes!\n",
	  progname, oblock_size*sizeof(float));
  
  for (islow=0; islow<slow_blksize; islow++) {
    float* block= mri_get_chunk(Input, this_chunk, 2*collective_blksize,
				in_offset, MRI_FLOAT);
    for (ifast=0; ifast<fast_blksize; ifast++) {
      FComplex* irunner;
      FComplex* crunner= buf;
      for (irunner= ((FComplex*)block)+ifast;
	   irunner<((FComplex*)block)+collective_blksize;
	   irunner += fast_blksize) {
	crunner->real= irunner->real;
	crunner->imag= irunner->imag;
	crunner++;
      }
      fft3d( buf, selected_extent_slow, selected_extent_med,
	     selected_extent_fast, fft_sign, "xyz" );
      switch (resultType) {
      case RSLT_COMPLEX:
	{
	  FComplex* orunner= ((FComplex*)oblock)+ifast;
	  for (crunner= buf; crunner<buf + fft_blksize; crunner++) {
	    orunner->real= crunner->real;
	    orunner->imag= crunner->imag;
	    orunner += fast_blksize;
	  }
	}
	break;
      case RSLT_MODULUS:
	{
	  float* orunner= oblock+ifast;
	  for (crunner= buf; crunner<buf + fft_blksize; crunner++) {
	    *orunner= Modulus( *crunner );
	    orunner += fast_blksize;
	  }
	}
	break;
      case RSLT_PHASE:
	{
	  float* orunner= oblock+ifast;
	  for (crunner= buf; crunner<buf + fft_blksize; crunner++) {
	    *orunner= Phase( *crunner );
	    orunner += fast_blksize;
	  }
	}
	break;
      case RSLT_SQMOD:
	{
	  float* orunner= oblock+ifast;
	  for (crunner= buf; crunner<buf + fft_blksize; crunner++) {
	    *orunner= 
	      crunner->real*crunner->real + crunner->imag*crunner->imag;
	    orunner += fast_blksize;
	  }
	}
	break;
      case RSLT_REAL:
	{
	  float* orunner= oblock+ifast;
	  for (crunner= buf; crunner<buf + fft_blksize; crunner++) {
	    *orunner= crunner->real;
	    orunner += fast_blksize;
	  }
	}
	break;
      case RSLT_IMAG:
	{
	  float* orunner= oblock+ifast;
	  for (crunner= buf; crunner<buf + fft_blksize; crunner++) {
	    *orunner= crunner->imag;
	    orunner += fast_blksize;
	  }
	}
	break;
      }
    }
    mri_set_chunk( Output, this_chunk, oblock_size, out_offset,
		   MRI_FLOAT, oblock );
    if (verbose_flg) 
      fprintf(stderr,"block: %ld at %lld -> %ld at %lld\n",
	      2*collective_blksize, in_offset,
	      oblock_size, out_offset);
    in_offset += 2*collective_blksize;
    out_offset += oblock_size;
  }
  
  free(buf);
  free(oblock);
}

static void fft_chunk_3d(char* this_chunk) {
  char key_buf[KEYBUF_SIZE];
  char* dimstr;
  long selected_extent_fast, selected_extent_med, selected_extent_slow;
  int chunk_complex= 0; 
  
  if (strlen(selected_dim) != 3)
    Abort("%s: fft_chunk_3d called incorrectly!\n");
  if (verbose_flg) fprintf(stderr,"ffting chunk <%s> in 3D\n",this_chunk);
  safe_copy(key_buf, this_chunk);
  safe_concat(key_buf, ".dimensions");
  if (mri_has(Input,key_buf)) {
    dimstr= mri_get_string(Input,key_buf);
    if (verbose_flg) fprintf(stderr,"dimstr %s\n",dimstr);
    if (strstr(dimstr,selected_dim)) {
      selected_extent_fast= safe_get_extent(Input, this_chunk, selected_dim);
      selected_extent_med= safe_get_extent(Input, this_chunk, 
					   selected_dim+1);
      selected_extent_slow= safe_get_extent(Input, this_chunk, 
					    selected_dim+2);
      if (*dimstr == 'v') {
	if (safe_get_extent(Input, this_chunk, "v") == 2) chunk_complex= 1;
	else chunk_complex= 0;
      }
      
      if (verbose_flg) {
	if (chunk_complex)
	  fprintf(stderr,
		  "selected dim extent(s) on input %ld %ld %ld, chunk complex\n",
		  selected_extent_fast, selected_extent_med,
		  selected_extent_slow);
	else
	  fprintf(stderr,
		  "selected dim extent(s) on input %ld %ld %ld, chunk scalar\n",
		  selected_extent_fast, selected_extent_med, 
		  selected_extent_slow);
      }
      
      if (*selected_dim == 'v' && chunk_complex)
	Abort("%s: attempted to FFT the complex v dimension!\n",progname);
      
      if (selected_extent_fast*selected_extent_med*selected_extent_slow > 1) {
	long long fast_blksize;
	long long slow_blksize;
	
	data_changed= 1;
	
	calc_sizes(this_chunk, dimstr, &fast_blksize, &slow_blksize);
	if (chunk_complex) fast_blksize /= 2; /* count complex, not float */
	safe_copy(key_buf, this_chunk);
	safe_concat(key_buf, ".datatype");
	mri_set_string(Output, key_buf, "float32");
	if (chunk_complex) {
	  if ((resultType==RSLT_MODULUS) || (resultType==RSLT_PHASE)
	      || (resultType==RSLT_SQMOD) || (resultType==RSLT_REAL)
	      || (resultType==RSLT_IMAG))
	    scalarize_chunk(Output,this_chunk);
	  transfer_data_3d_complex(this_chunk, fast_blksize, slow_blksize, 
				   selected_extent_fast, selected_extent_med,
				   selected_extent_slow);
	}
	else {
	  if (resultType==RSLT_COMPLEX)
	    complexify_chunk(Output,this_chunk);
	  transfer_data_3d_scalar(this_chunk, fast_blksize, slow_blksize, 
				  selected_extent_fast, selected_extent_med,
				  selected_extent_slow);
	}
      }
    }
    else {
      /* Chunk copied correctly in initial dataset copy */
      if (verbose_flg) fprintf(stderr,"This chunk requires no fft\n");
    }
  }
}

static void fft_chunk(char* this_chunk) {
  char key_buf[KEYBUF_SIZE];
  char* dimstr;
  long selected_extent_fast, selected_extent_slow;
  int chunk_complex= 0; 
  
  if (verbose_flg) fprintf(stderr,"ffting chunk <%s>\n",this_chunk);
  safe_copy(key_buf, this_chunk);
  safe_concat(key_buf, ".dimensions");
  if (mri_has(Input,key_buf)) {
    dimstr= mri_get_string(Input,key_buf);
    if (verbose_flg) fprintf(stderr,"dimstr %s\n",dimstr);
    if (strstr(dimstr,selected_dim)) {
      selected_extent_fast= safe_get_extent(Input, this_chunk, selected_dim);
      if (strlen(selected_dim)==2) 
	selected_extent_slow= safe_get_extent(Input, this_chunk, 
					      selected_dim+1);
      else selected_extent_slow= 1;
      if (*dimstr == 'v') {
	if (safe_get_extent(Input, this_chunk, "v") == 2) chunk_complex= 1;
	else chunk_complex= 0;
      }
      
      if (verbose_flg) {
	if (chunk_complex)
	  fprintf(stderr,
		  "selected dim extent(s) on input %ld %ld, chunk complex\n",
		  selected_extent_fast, selected_extent_slow);
	else
	  fprintf(stderr,
		  "selected dim extent(s) on input %ld %ld, chunk scalar\n",
		  selected_extent_fast, selected_extent_slow);
      }
      
      if (*selected_dim == 'v' && chunk_complex)
	Abort("%s: attempted to FFT the complex v dimension!\n",progname);
      
      if (selected_extent_fast*selected_extent_slow > 1) {
	long long fast_blksize;
	long long slow_blksize;
	
	data_changed= 1;
	
	calc_sizes(this_chunk, dimstr, &fast_blksize, &slow_blksize);
	if (chunk_complex) fast_blksize /= 2; /* count complex, not float */
	safe_copy(key_buf, this_chunk);
	safe_concat(key_buf, ".datatype");
	mri_set_string(Output, key_buf, "float32");
	if (chunk_complex) {
	  if ((resultType==RSLT_MODULUS) || (resultType==RSLT_PHASE)
	      || (resultType==RSLT_SQMOD) || (resultType==RSLT_REAL)
	      || (resultType==RSLT_IMAG))
	    scalarize_chunk(Output,this_chunk);
	  transfer_data_complex(this_chunk, fast_blksize, slow_blksize, 
				selected_extent_fast, selected_extent_slow);
	}
	else {
	  if (resultType==RSLT_COMPLEX)
	    complexify_chunk(Output,this_chunk);
	  transfer_data_scalar(this_chunk, fast_blksize, slow_blksize, 
			       selected_extent_fast, selected_extent_slow);
	}
      }
    }
    else {
      /* Chunk copied correctly in initial dataset copy */
      if (verbose_flg) fprintf(stderr,"This chunk requires no fft\n");
    }
  }
}

static void flip_description_string( MRI_Dataset* ds, const char* chunk, 
				     const char c )
{
  char buf[KEYBUF_SIZE];

  buf[KEYBUF_SIZE-1]= '\0';
  snprintf(buf,KEYBUF_SIZE-1,"%s.description.%c",chunk,c);
  if (mri_has(ds,buf)) {
    const char* oldDesc= mri_get_string(ds, buf);
    if (strstr(oldDesc,"ungridded"))
      Warning(1,"%s: Warning: %c dimension is believed to be ungridded!",
	      progname,c);
    if (strstr(oldDesc,"discrete"))
      Warning(1,"%s: Warning: %c dimension is believed to be discrete!",
	      progname,c);
    if (strstr(oldDesc,"k-space"))
      mri_set_string(ds,buf,"gridded image-space");
    else if (strstr(oldDesc,"image-space"))
      mri_set_string(ds,buf,"gridded k-space");	      
  }
}

int main( int argc, char* argv[] ) 
{
  char infile[512], outfile[512];
  char* this_key;
  char this_chunk[KEYBUF_SIZE];
  unsigned char **missing = NULL;
  int count_set_type= 0;

  progname= argv[0];

  /* Print version number */
  Message( "# %s\n", rcsid );

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Deprecate old options */

  if (cl_present( "f" ))
     Abort ("Option f has been expanded to forward|fwd.  Please see help file.\n");
  if (cl_present( "i" ))
     Abort ("Option i has been expanded to inverse|inv.  Please see help file.\n");
  if (cl_present( "c" ))
     Abort ("Option c has been expanded to complex|cpx.  Please see help file.\n");
  if (cl_present( "m" ))
     Abort ("Option m has been expanded to modulus|mod.  Please see help file.\n");
  if (cl_present( "p" ))
     Abort ("Option p has been expanded to phase|pha.  Please see help file.\n");
  if (cl_present( "s" ))
     Abort ("Option s has been expanded to squared_modulus|sqr.  Please see help file.\n");


  /* Get filenames */
  cl_get("dimension|dim|d", "%option %s[t]",selected_dim);
  if (strlen(selected_dim)>3) {
    fprintf(stderr,
	    "%s: Selected dim name must be 3 or fewer chars long.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  *selected_dim= tolower(*selected_dim);

  verbose_flg= cl_present("verbose|ver|v");

  if (cl_present("forward|fwd")) {
    if (cl_present("inverse|inv"))
      Abort("%s: both -forward|fwd and -inverse|inv are present!\n",argv[0]);
    fft_sign= 1;
  }
  if (cl_present("inverse|inv")) fft_sign= -1;

  count_set_type= 0;
  if (cl_present("complex|cpx")) {
    count_set_type += 1;
    resultType= RSLT_COMPLEX;
  }
  if (cl_present("modulus|mod")) {
    count_set_type += 1;
    resultType= RSLT_MODULUS;
  }
  if (cl_present("phase|pha")) {
    count_set_type += 1;
    resultType= RSLT_PHASE;
  }
  if (cl_present("squared_modulus|sqr")) {
    count_set_type += 1;
    resultType= RSLT_SQMOD;
  }
  if (cl_present("real|rea")) {
    count_set_type += 1;
    resultType= RSLT_REAL;
  }
  if (cl_present("imaginary|ima")) {
    count_set_type += 1;
    resultType= RSLT_IMAG;
  }
  if (count_set_type>1)
      Abort("%s: at least two of -complex|cpx, -modulus|mod, -squared_modulus|sqr -phase|pha -real|rea -imaginary|ima are present!\n",argv[0]);

  if (!cl_get("", "%s", infile)) {
    fprintf(stderr,"%s: Input file name not given.\n",argv[0]);
    Help( "usage" );
    exit(-1);
  }
  if (!cl_get("", "%s", outfile)) {
    fprintf(stderr,"%s: Output file name not given.\n",argv[0]);
    Help( "usage" );
    exit(-1);
  }

  if (cl_cleanup_check()) {
    int i;
    fprintf(stderr,"%s: invalid argument in command line:\n    ",argv[0]);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(-1);
  }
  
  /*** End command-line parsing ***/
  
  /* Open input and output datasets */
  if( !strcmp( infile, outfile ) )
    Abort( "Input and output files must be distinct." );
  Input = mri_open_dataset( infile, MRI_READ );
  Output = mri_copy_dataset( outfile, Input );
  hist_add_cl( Output, argc, argv );
  
  /* Walk through the input handling each chunk in turn */
  mri_iterate_over_keys(Input);
  while ((this_key= mri_next_key(Input)) != NULL) {
    if (strlen(this_key)>=KEYBUF_SIZE) 
      Abort("%s: key too long!\n",argv[0]);
    if (!strcmp(mri_get_string(Input,this_key),"[chunk]")) {
      strncpy(this_chunk, this_key, KEYBUF_SIZE);
      this_chunk[KEYBUF_SIZE-1]= '\0';
      if (strcmp(this_chunk,"missing")) { /* don't fft missing chunk! */
	const char* runner;
	for (runner=selected_dim; *runner; runner++) {
	  flip_description_string( Output, this_chunk, *runner );
	}
  
	if (strlen(selected_dim)==3)
	  fft_chunk_3d(this_chunk);
	else
	  fft_chunk(this_chunk);
      }
    }
  }

  /* Write and close data-sets */
  mri_close_dataset( Input );
  mri_close_dataset( Output );
  
  Message( "#      FFT complete.\n" );
  if (!data_changed) 
    Message("#      Warning: input and output datasets identical!\n");
  
  return 0;
}

