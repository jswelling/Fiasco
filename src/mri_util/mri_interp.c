/************************************************************
 *                                                          *
 *  mri_interp.c                                            *
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
 *  Original programming by Joel Welling, 5/98              *
 ************************************************************/
/*************************************************************

  DESCRIPTION OF MRI_INTERP

  mri_interp takes a pgh MRI dataset of any type, and outputs
  a dataset of the same type containing interpolated data.
  The interpolated data spans the range of one of the dimensions
  of the input data, but has a different extent than that of
  the input dimension.  The interpolation operation acts on all 
  chunks in the dataset.

  Chunks not containing the selected dimension are copied verbatim.
  All chunks containing the given dimension will be interpolated.
  Note that this can give confusing results if the extents of the
  dimension vary between chunks!

**************************************************************/

#include <stdlib.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "misc.h"
#include "stdcrg.h"

#define KEYBUF_SIZE 512

static char rcsid[] = "$Id: mri_interp.c,v 1.14 2007/07/06 18:45:53 welling Exp $";

static MRI_Dataset *Input = NULL, *Output = NULL;
static char selected_dim[512]= "";
static int new_extent;
static int data_changed= 0;
static char* progname;
static int const_flag= 0;
static int verbose_flag= 0;
static int reallyverbose_flag= 0;

static void safe_copy(char* str1, char* str2) {
  strncpy(str1, str2, KEYBUF_SIZE);
  str1[KEYBUF_SIZE-1]= '\0';
}

static void safe_concat(char* str1, char* str2) {
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
		       int* fast_blocksize_out, int* slow_blocksize_out ) {
  /* This routine will fail if selected_dim is not in dimstr! */
  int fast_blocksize;
  int slow_blocksize;
  char* this_dim;

  fast_blocksize= 1;

  this_dim= dimstr;
  while (*this_dim != *selected_dim) 
    fast_blocksize *= safe_get_extent(Input,this_chunk,this_dim++);

  this_dim++; /* step over selected dim */

  slow_blocksize= 1;
  while (*this_dim)
    slow_blocksize *= safe_get_extent(Input, this_chunk, this_dim++);

  *fast_blocksize_out= fast_blocksize;
  *slow_blocksize_out= slow_blocksize;
}

static void missing_interp( void* out_out, void* low_in, void* high_in, 
			    float shift, int n, int type )
{
  int i;

  if (type==MRI_UNSIGNED_CHAR) {
    unsigned char* low= low_in;
    unsigned char* high= high_in;
    unsigned char* out= out_out;
    for (i=0; i<n; i++) 
      out[i]= (unsigned char)(high[i] || low[i]);
  }
  else Abort("%s: missing_interp: missing chunk is not unsigned char!\n",
	     progname);
}

static void type_interp( void* out_out, void* low_in, void* high_in, 
			 float shift, int n, int type )
{
  int i;

  switch (type) {

  case MRI_UNSIGNED_CHAR: {
    unsigned char* low= low_in;
    unsigned char* high= high_in;
    unsigned char* out= out_out;
    for (i=0; i<n; i++) 
      out[i]= (unsigned char)(shift*high[i] + (1.0-shift)*low[i]);
  }
  break;

  case MRI_SHORT: {
    short* low= low_in;
    short* high= high_in;
    short* out= out_out;
    for (i=0; i<n; i++) 
      out[i]= (short)(shift*high[i] + (1.0-shift)*low[i]);
  }
  break;

  case MRI_INT: {
    int* low= low_in;
    int* high= high_in;
    int* out= out_out;
    for (i=0; i<n; i++) 
      out[i]= (int)(shift*high[i] + (1.0-shift)*low[i]);
  }
  break;

  case MRI_LONGLONG: {
    long long* low= low_in;
    long long* high= high_in;
    long long* out= out_out;
    for (i=0; i<n; i++) 
      out[i]= (long long)(shift*high[i] + (1.0-shift)*low[i]);
  }
  break;

  case MRI_FLOAT: {
    float* low= low_in;
    float* high= high_in;
    float* out= out_out;
    for (i=0; i<n; i++) 
      out[i]= (shift*high[i] + (1.0-shift)*low[i]);
  }
  break;

  case MRI_DOUBLE: {
    double* low= low_in;
    double* high= high_in;
    double* out= out_out;
    for (i=0; i<n; i++) 
      out[i]= ((double)shift*high[i] + (1.0-shift)*low[i]);
  }
  break;

  default:
    Abort("%s: type_interp: unrecognized type %d!\n",progname,type);
  }
}

static void transfer_data(char* this_chunk, 
			  int fast_blksize, int slow_blksize, 
			  int selected_extent, int new_extent)
{
  int in_offset= 0;
  int out_offset= 0;
  int in_framestart= 0;
  int ifast;
  int islow;
  int type;
  int typesize;
  void* obuf= NULL;
  int missing_chunk= 0;

  type= get_chunk_type(Input,this_chunk);
  typesize= get_typesize(Input,this_chunk);

  missing_chunk= !strcmp(this_chunk,"missing");

  for (islow=0; islow<slow_blksize; islow++) {
    in_offset= in_framestart;
    if (const_flag) {
      /* const interpolation- this handles missing correctly */
      double ratio= ((double)selected_extent)/((double)new_extent);
      double row= 0.0;
      int old_row= -1;
      int i;
      void* in_chunk;
      if (reallyverbose_flag) 
	fprintf(stderr,"const block, ratio %f, %d * %d at %d -> %d at %d\n",
		ratio, fast_blksize, selected_extent, in_offset,
		new_extent, out_offset);
      for (i=0; i<new_extent; i++) {
	if ((int)row != old_row) {
	  /* load a new row */
	  in_chunk= mri_get_chunk(Input, this_chunk, fast_blksize, 
				  in_offset, type);
	  old_row= (int)row;
	  in_offset += fast_blksize;
	}
	mri_set_chunk(Output, this_chunk, fast_blksize, out_offset, type,
		      in_chunk);
	out_offset += fast_blksize;
	row += ratio;
      }
    }
    else {
      /* linear interpolation */
      double ratio= ((double)(selected_extent-1))/((double)(new_extent-1));
      double shift= 0.0;
      int i;
      void* low_chunk;
      void* high_chunk;
      if (!obuf) {
	if (!(obuf= (void*)malloc(fast_blksize*typesize)))
	  Abort("%s: unable to allocate %d bytes!\n",progname,
		fast_blksize*typesize);
      }
      if (reallyverbose_flag) 
	fprintf(stderr,"linear block, ratio %f, %d * %d at %d -> %d at %d\n",
		ratio, fast_blksize, selected_extent, in_offset,
		new_extent, out_offset);
      low_chunk= mri_get_chunk(Input, this_chunk, fast_blksize, 
			       in_offset, type);
      high_chunk= mri_get_chunk(Input, this_chunk, fast_blksize,
				in_offset + fast_blksize, type);
      for (i=0; i<new_extent; i++) {
	if (shift>=1.0) {
	  /* Load a new row */
	  shift -= 1.0;
	  in_offset += fast_blksize;
	  low_chunk= high_chunk;
	  if (i!=(new_extent-1)) /* very last row is on boundary */
	    high_chunk= mri_get_chunk(Input, this_chunk, fast_blksize,
				      in_offset + fast_blksize, type);
	}
	if (missing_chunk)
	  missing_interp( obuf, low_chunk, high_chunk, shift, fast_blksize,
			  type );
	else type_interp( obuf, low_chunk, high_chunk, shift, fast_blksize, 
			  type );
	mri_set_chunk(Output, this_chunk, fast_blksize, out_offset, type,
		      obuf);
	out_offset += fast_blksize;
	shift += ratio;
      }
    }
    in_framestart += fast_blksize*selected_extent;
  }
  free(obuf);
}

static void interp_chunk(char* this_chunk) {
  char key_buf[KEYBUF_SIZE];
  char* dimstr;
  int selected_extent;

  if (reallyverbose_flag) fprintf(stderr,"interpolating chunk <%s>\n",this_chunk);
  safe_copy(key_buf, this_chunk);
  safe_concat(key_buf, ".dimensions");
  if (mri_has(Input,key_buf)) {
    dimstr= mri_get_string(Input,key_buf);
    if (reallyverbose_flag) fprintf(stderr,"dimstr %s\n",dimstr);
    if (strchr(dimstr,*selected_dim)) {
      safe_copy(key_buf, this_chunk);
      safe_concat(key_buf, ".extent.");
      safe_concat(key_buf, selected_dim);
      if (mri_has(Input,key_buf)) 
	selected_extent= mri_get_int(Input,key_buf);
      else Abort("%s: input missing tag %s!\n",progname,key_buf);
      if ((!const_flag) && (selected_extent<2)) {
	Abort("%s: selected dim in chunk %s is too small!  Do you want -c flag?",
	      progname,this_chunk);
      }
      if (const_flag && (selected_extent<1)) {
	Abort("%s: selected dim in chunk %s is too small to interpolate!",
	      progname,this_chunk);
      }
      if (reallyverbose_flag) 
	fprintf(stderr,"selected dim extent on input %d\n",selected_extent);
      if (new_extent != selected_extent) {
	int fast_blksize;
	int slow_blksize;

	data_changed= 1;
	mri_set_int(Output, key_buf, new_extent);
	calc_sizes(this_chunk, dimstr, &fast_blksize, &slow_blksize);
	transfer_data(this_chunk, fast_blksize, slow_blksize, 
		      selected_extent, new_extent);
      }
    }
    else {
      /* Chunk copied correctly in initial dataset copy */
    }
  }
}

int main( int argc, char* argv[] ) 
{

  char infile[512], outfile[512];
  char* this_key;
  char this_chunk[KEYBUF_SIZE];

  progname= argv[0];

  /* Print version number */
  if (verbose_flag) Message( "# %s\n", rcsid );

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Deprecate old options */

  if (cl_present( "l" ))
     Abort ("Option l used to designate linear interpolation and has been expanded to linear|lin.  If you mean length, please use length|len.  Please see help file.\n");

  if (cl_present( "c" ))
     Abort ("Option c(constant) has been expanded to constant|con.  Please see help file.\n");

  if (cl_present( "e" ))
     Abort ("Option e(extent) has been replaced by length|len|l.  Please see help file.\n");


  /* Get filenames */
  if (!cl_get("dimension|dim|d", "%option %s[t]",selected_dim)) {
    fprintf(stderr,"%s: dimension to be interpolated not given.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  if (strlen(selected_dim)>1) {
    fprintf(stderr,"%s: Selected dim name must be 1 char long.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  *selected_dim= tolower(*selected_dim);
  if (!cl_get("length|len|l", "%option %d",&new_extent)) {
    fprintf(stderr,"%s: new extent not given.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
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
  if (cl_present("linear|lin")) {
    if (!cl_present("constant|con")) const_flag= 0;
    else {
      fprintf(stderr,"%s: -linear and -constant flags are mutually exclusive.\n",argv[0]);
      Help( "usage" );
      exit(-1);
    }
  }
  else {
    if (cl_present("constant|con")) const_flag= 1;
    else const_flag= 0;
  }
  verbose_flag= cl_present("verbose|ver|v");
  reallyverbose_flag= cl_present("V");
  if (reallyverbose_flag) verbose_flag= 1;
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
    Abort( "%s: Input and output files must be distinct.",argv[0] );
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
      interp_chunk(this_chunk);
    }
  }

  /* Write and close data-sets */
  mri_close_dataset( Input );
  mri_close_dataset( Output );
  
  if (verbose_flag) Message( "#      Interpolation complete.\n" );
  if (!data_changed) 
    Message("#      Warning: input and output datasets identical!\n");

  return 0;
}

