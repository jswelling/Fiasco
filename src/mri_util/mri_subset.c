/************************************************************
 *                                                          *
 *  mri_subset.c                                            *
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

  DESCRIPTION OF MRI_SUBSET

  mri_subset takes a pgh MRI dataset of any type, and outputs
  a dataset of the same type containing a subset of the data.
  The subset represents a subrange of one of the dimensions of
  the data.  The subset operation acts on all chunks in the
  dataset.

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

#define KEYBUF_SIZE 512

/* Never transfer more than 5M at a time */
#define MAX_BLOCKSIZE 5242880

static char rcsid[] = "$Id: mri_subset.c,v 1.20 2007/07/27 17:29:43 welling Exp $";

static MRI_Dataset *Input = NULL, *Output = NULL;
static char selected_dim[512]= "t";
static long offset;
static long range;
static int data_changed= 0;
static char* progname;
static int reallyverbose_flg= 0;
static int verbose_flg= 0;

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

static long safe_get_extent(MRI_Dataset* ds, char* chunk, char* dim)
{
  char key_buf[KEYBUF_SIZE];
  char dim_buf[4];
  dim_buf[0]= *dim;
  dim_buf[1]= '\0';
  safe_copy(key_buf,chunk);
  safe_concat(key_buf,".extent.");
  safe_concat(key_buf,dim_buf);
  if (mri_has(ds,key_buf)) return (long)mri_get_int(ds,key_buf);
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

  this_dim++; /* step over selected dim */

  slow_blocksize= 1;
  while (*this_dim)
    slow_blocksize *= safe_get_extent(Input, this_chunk, this_dim++);

  *fast_blocksize_out= fast_blocksize;
  *slow_blocksize_out= slow_blocksize;
}

static void transfer_data(char* this_chunk, 
			  long long fast_blksize, long long slow_blksize, 
			  long selected_extent, long truncated_range)
{
  long long in_offset= 0;
  long long out_offset= 0;
  long long ifast;
  long long islow;
  long long collective_blksize;
  int type;
  long num_this_block= 0;
  long long num_transferred;

  type= get_chunk_type(Input,this_chunk);

  collective_blksize= fast_blksize * truncated_range;
  for (islow=0; islow<slow_blksize; islow++) {
    num_transferred= 0;
    for (num_transferred=0; num_transferred<collective_blksize;
	 num_transferred += num_this_block) {
      if (collective_blksize-num_transferred > MAX_BLOCKSIZE) 
	num_this_block= MAX_BLOCKSIZE;
      else num_this_block= collective_blksize-num_transferred;
      mri_set_chunk( Output, this_chunk, num_this_block, 
		     out_offset+num_transferred, type,
		     mri_get_chunk(Input, this_chunk, num_this_block,
				   in_offset + (fast_blksize*offset) 
				   + num_transferred, 
				   type) );
    }
    if (reallyverbose_flg) 
      fprintf(stderr,"block: type %d, %lld at %lld -> %lld at %lld\n",
	      type, collective_blksize, in_offset + (fast_blksize*offset),
	      collective_blksize, out_offset);
    in_offset += fast_blksize * selected_extent;
    out_offset += collective_blksize;
  }
}

static int subsetIndexRemap(int oldIndex, void* hook)
{
  int newIndex= oldIndex-offset;
  if ((newIndex<0) || (newIndex>=range)) return -1;
  else return newIndex;
}

static void subset_chunk(char* this_chunk) {
  char key_buf[KEYBUF_SIZE];
  char* dimstr;
  long selected_extent;

  if (reallyverbose_flg) fprintf(stderr,"subsetting chunk <%s>\n",this_chunk);
  safe_copy(key_buf, this_chunk);
  safe_concat(key_buf, ".dimensions");
  if (mri_has(Input,key_buf)) {
    dimstr= mri_get_string(Input,key_buf);
    if (reallyverbose_flg) fprintf(stderr,"dimstr %s\n",dimstr);
    if (strchr(dimstr,*selected_dim)) {
      safe_copy(key_buf, this_chunk);
      safe_concat(key_buf, ".extent.");
      safe_concat(key_buf, selected_dim);
      if (mri_has(Input,key_buf)) 
	selected_extent= mri_get_int(Input,key_buf);
      else Abort("%s: input missing tag %s!\n",progname,key_buf);
      if (reallyverbose_flg) 
	fprintf(stderr,"selected dim extent on input %ld\n",selected_extent);
      if ((offset>0) || (selected_extent>range)) {
	long long fast_blksize;
	long long slow_blksize;
	long truncated_range;

	data_changed= 1;
	calc_sizes(this_chunk, dimstr, &fast_blksize, &slow_blksize);
	truncated_range= 
	  (range+offset <= selected_extent) ? range : (selected_extent-offset);
	if (truncated_range <= 0)
	  Abort("%s: chunk <%s> has a range too small for this subset!\n",
		progname, this_chunk);
	if (truncated_range != range)
	  Warning(1,"%s: range truncated to %ld!\n",
		  progname, truncated_range);
	mriu_updateLabels(Output, this_chunk, selected_dim[0], 
			  0, selected_extent,
			  subsetIndexRemap, NULL);
	mri_set_int(Output, key_buf, truncated_range);
	transfer_data(this_chunk, fast_blksize, slow_blksize, 
		      selected_extent, truncated_range);
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
  if (verbose_flg) Message( "# %s\n", rcsid );

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Deprecate old options */

  if (cl_present("r"))
    Abort ("Option r(range) has been replaced by length|len|l.  Please see help file. \n");

  if (cl_present("o"))
    Abort ("Option o(offset) has been replaced by shift|shi|s.  Please see help file. \n");

  /* Get filenames */
  cl_get("dimension|dim|d", "%option %s[t]",selected_dim);
  if (strlen(selected_dim)>1) {
    fprintf(stderr,"%s: Selected dim name must be 1 char long.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  *selected_dim= tolower(*selected_dim);

  if (!cl_get("shift|shi|s", "%option %ld",&offset)) {
    fprintf(stderr,"%s: Index shift not given.\n",argv[0]);
    Help("usage");
    exit(-1);
  }

  cl_get("length|len|l", "%option %ld[1]",&range);

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
  verbose_flg= cl_present("verbose|ver|v");
  reallyverbose_flg= cl_present("V|debug");
  if (reallyverbose_flg) verbose_flg= 1;
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
      subset_chunk(this_chunk);
    }
  }

  /* Write and close data-sets */
  mri_close_dataset( Input );
  mri_close_dataset( Output );
  
  if (verbose_flg) Message( "#      Subset extraction complete.\n" );
  if (!data_changed) 
    Message("#      Warning: input and output datasets identical!\n");

  return 0;
}
