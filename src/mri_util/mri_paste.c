/************************************************************
 *                                                          *
 *  mri_paste.c                                             *
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "misc.h"
#include "stdcrg.h"

#define KEYBUF_SIZE 512
#define DEFAULT_OUT_NAME "mri_paste_out"
#define MAX_INPUT_FILES 64
#define MAX_BUF 1000000

static char rcsid[] = "$Id: mri_paste.c,v 1.27 2007/07/06 18:45:53 welling Exp $";

typedef struct mrifile_struct {
  MRI_Dataset* ds;
  char* fname;
  char* chunk;
  long long offset;
  long long length;
} MRIFile;

static MRI_Dataset *Input[MAX_INPUT_FILES], *Output = NULL;
static char selected_dim[512]= "t";
static int range;
static long long offset;
static int data_changed= 0;
static float fillval= 0.0;
static char* progname;
static int verbose_flg= 0;
static int reallyverbose_flg= 0;
static int n_input_files= 0;
static int selected_extent[MAX_INPUT_FILES];
static int first_dim= 0; /* special case: paste first non-trivial dim */
static int last_dim=1;   /* special case: paste last dim */

static long get_max_buf( const int typesize )
{
  /* If the environment provides a memory size hint, use it.
   * Otherwise use the compiled-in value.  Note that the compiled-
   * in value is in "units" while the memory hint is in bytes.
   */
  char* here;
  long long default_memlimit;
  /* Allow the user to use more memory */
  if ((here=getenv("F_MEMSIZE_HINT")) != NULL) {
    long default_memlimit= atol(here);
    if (default_memlimit==0)
      Abort("%s: environment variable F_MEMSIZE_HINT is not a long integer!\n",
            progname);
    return default_memlimit/typesize;
  }
  else return MAX_BUF;
}

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
    else if (!strcmp(type_name,"int64")) return MRI_LONGLONG;
    else if (!strcmp(type_name,"float32")) return MRI_FLOAT;
    else if (!strcmp(type_name,"float64")) return MRI_DOUBLE;
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
  if (mri_has(ds,key_buf)) return mri_get_int(ds,key_buf);
  else Abort("%s: input file missing tag %s!\n",progname,key_buf);
  return 0; /* not reached */
}

static void check_order(char* this_chunk, char* dimstr)
{
  char key_buf[KEYBUF_SIZE];
  int i=1;
  
  safe_copy(key_buf, this_chunk);
  safe_concat(key_buf, ".dimensions");

  while (i<n_input_files) {
    if (strcmp(dimstr, mri_get_string(Input[i], key_buf)))
	Abort("%s: dimension order in input %d doesn't match order in first input\n",
	      progname, i+1);
    i++;
  }
}

static void calc_sizes(char* this_chunk, char* dimstr, 
		       long long* fast_blocksize_out, 
		       long long* slow_blksize_out) {
  /* This routine will fail if selected_dim is not in dimstr! */
  long long fast_blocksize;
  long long slow_blksize;
  int dim_extent;
  int j;
  char* this_dim;

  fast_blocksize= 1;

  this_dim= dimstr;
  while (*this_dim != *selected_dim) {
    dim_extent = safe_get_extent(Input[0],this_chunk,this_dim);
    for (j=1; j<n_input_files; j++)
      if (safe_get_extent(Input[j],this_chunk,this_dim)!=dim_extent)
	Abort("%s: dimension %c extent in input %d doesn't match extent in first input\n", progname, *this_dim, j+1);
    fast_blocksize *= dim_extent;
    this_dim++;
  }

  if (fast_blocksize == 1) first_dim= 1;
  else first_dim= 0;
  
  this_dim++; /* step over selected dim */

  slow_blksize= 1;
  last_dim= 1; 
  while (*this_dim){
    dim_extent = safe_get_extent(Input[0], this_chunk, this_dim);
    for (j=1; j<n_input_files; j++)
      if (safe_get_extent(Input[j],this_chunk,this_dim)!=dim_extent)
	Abort("%s: dimension %c extent doesn't match extent in input file %d\n", progname, *this_dim, j+1);
    slow_blksize *= dim_extent;
    this_dim++;
    last_dim=0;
  }

  *fast_blocksize_out= fast_blocksize;
  *slow_blksize_out= slow_blksize;
}

static void transfer_data(char* this_chunk, 
			  long long fast_blksize, long long slow_blksize)
{
  long long in_offset[MAX_INPUT_FILES];
  long long out_offset=0;
  long long islow;
  long long collective_blksize;
  long long out_blksize;
  int type;
  int typesize;
  void* obuf;
  void* ibuf;
  long long i,j;

  for (i=0; i<n_input_files; i++) in_offset[i]=0;

  type= get_chunk_type(Input[0],this_chunk);
  typesize= get_typesize(Input[0],this_chunk);

  for (islow=0; islow<slow_blksize; islow++) {

    out_blksize = fast_blksize * range;

    if (!(obuf=malloc(out_blksize * typesize))) {
      Abort("%s: unable to allocate %d bytes!\n",progname,
	  out_blksize*typesize);
    }


    for (j=0; j<n_input_files; j++) {

      collective_blksize= fast_blksize * selected_extent[j];
      
      if (type!=get_chunk_type(Input[j],this_chunk)) {
	Abort("%s: chunk type in input %d doesn't match type in first input\n",
	      progname, j+1);
      }
      
      ibuf= mri_get_chunk(Input[j], this_chunk, collective_blksize, 
			  in_offset[j], type);
      switch (type) {
      case MRI_UNSIGNED_CHAR:
	{
	  char* tobuf= (char*)obuf;
	  char* tibuf= (char*)ibuf;
	  for (i=0; i<collective_blksize; i++) *tobuf++= *tibuf++;
	}
	break;
      case MRI_SHORT:
	{
	  short* tobuf= (short*)obuf;
	  short* tibuf= (short*)ibuf;
	  for (i=0; i<collective_blksize; i++) *tobuf++= *tibuf++;
	}
	break;
      case MRI_INT:
	{
	  int* tobuf= (int*)obuf;
	  int* tibuf= (int*)ibuf;
	  for (i=0; i<collective_blksize; i++) *tobuf++= *tibuf++;
	}
	break;
      case MRI_LONGLONG:
	{
	  long long* tobuf= (long long*)obuf;
	  long long* tibuf= (long long*)ibuf;
	  for (i=0; i<collective_blksize; i++) *tobuf++= *tibuf++;
	}
	break;
      case MRI_FLOAT:
	{
	  float* tobuf= (float*)obuf;
	  float* tibuf= (float*)ibuf;
	  for (i=0; i<collective_blksize; i++) *tobuf++= *tibuf++;
	}
	break;
      case MRI_DOUBLE:
	{
	  double* tobuf= (double*)obuf;
	  double* tibuf= (double*)ibuf;
	  for (i=0; i<collective_blksize; i++) *tobuf++= *tibuf++;
	}
	break;
      }
      
      mri_set_chunk( Output, this_chunk, collective_blksize, out_offset, type,
		     obuf );

      if (reallyverbose_flg) 
	fprintf(stderr,"block: type %d, %lld at %lld -> %lld at %lld\n",
		type, collective_blksize, in_offset[j],
		collective_blksize, out_offset);
      
      in_offset[j] += collective_blksize;
      out_offset += collective_blksize; 
      
    }

    free(obuf);
  }

}

static void transfer_front_data(char* this_chunk, long long slow_blksize)
{
  /* Special case handling for pasting on first non-trivial dimension */

  long long in_offset[MAX_INPUT_FILES];
  long long out_offset=0;
  long long collective_blksize;
  int type;
  int typesize;
  long long blocks_per_step;
  long long blocks_this_step;
  long long blocks_to_go;
  long long buf_size;
  int block_shift;
  int i,j;
  long long iblock;
  void* obuf;
  void* ibuf;

  for (i=0; i<n_input_files; i++) in_offset[i]=0;

  type= get_chunk_type(Input[0],this_chunk);
  typesize= get_typesize(Input[0],this_chunk);

  blocks_per_step= (int)(get_max_buf(typesize)/range);
  if (blocks_per_step<1) blocks_per_step= 1;
  if (blocks_per_step>slow_blksize) blocks_per_step= slow_blksize;
  if (reallyverbose_flg) 
    fprintf(stderr,
	    "front paste special case; blocking factor %lld, slow blksize %lld\n",
	    blocks_per_step, slow_blksize);

  buf_size= range*blocks_per_step;
  if (!(obuf= (void*)malloc(buf_size*typesize)))
    Abort("%s: unable to allocate %d bytes!\n",buf_size*typesize);

  blocks_to_go= slow_blksize;
  while (blocks_to_go>0) {
    block_shift= 0;
    blocks_this_step= 
      (blocks_per_step>blocks_to_go) ? blocks_to_go : blocks_per_step;
    blocks_to_go -= blocks_this_step;
    for (i=0; i<n_input_files; i++) {
      void* orunner= (void*)((char*)obuf + typesize*block_shift);
      void* irunner;
      collective_blksize= blocks_this_step*selected_extent[i];
      irunner= ibuf= mri_get_chunk(Input[i], this_chunk, collective_blksize, 
				   in_offset[i], type);
      in_offset[i] += collective_blksize;
      for (iblock=0; iblock<blocks_this_step; iblock++) {
	memcpy(orunner, irunner, typesize*selected_extent[i]);
	irunner= (void*)((char*)irunner + typesize*selected_extent[i]);
	orunner= (void*)((char*)orunner + typesize*range);
      }
      block_shift += selected_extent[i];
    }
    mri_set_chunk( Output, this_chunk, range*blocks_this_step, 
		   out_offset, type, obuf );
    if (reallyverbose_flg) 
      fprintf(stderr,
	  "group: type %d, %lld blocks from %d inputs -> %lld at %lld (%lld to go)\n",
	      type, blocks_this_step, n_input_files, 
	      range*blocks_this_step, out_offset, blocks_to_go);
    out_offset += range*blocks_this_step;
  }

  free(obuf);
}

static void transfer_end_data(char* this_chunk, int fast_blksize)
{
  long long in_offset[MAX_INPUT_FILES];
  long long out_offset=0;
  long long islow;
  long long collective_blksize;
  long long transfer_total;
  int type;
  int typesize;
  void* obuf;
  void* ibuf;
  long long i;
  int j;

  for (j=0; j<n_input_files; j++) in_offset[j]=0;

  type= get_chunk_type(Input[0],this_chunk);
  typesize= get_typesize(Input[0],this_chunk);

  for (j=0; j<n_input_files; j++) {

    if (type!=get_chunk_type(Input[j],this_chunk)) {
      Abort("%s: chunk type in input %d doesn't match type in first input\n",
	  progname, j+1);
    }

    transfer_total= fast_blksize * selected_extent[j];

    while (transfer_total>0) {
      long max_buf= get_max_buf(typesize);
      if (transfer_total>max_buf)
	collective_blksize = max_buf;
      else
        collective_blksize = transfer_total;

      if (!(obuf=malloc(collective_blksize * typesize))) {
        Abort("%s: unable to allocate %d bytes!\n",progname,
	    collective_blksize*typesize);
      }

      ibuf= mri_get_chunk(Input[j], this_chunk, collective_blksize, 
	  		  in_offset[j], type);

      switch (type) {
      case MRI_UNSIGNED_CHAR:
        {
	  char* tobuf= (char*)obuf;
	  char* tibuf= (char*)ibuf;
	  for (i=0; i<collective_blksize; i++) *tobuf++= *tibuf++;
        }
        break;
      case MRI_SHORT:
        {
	  short* tobuf= (short*)obuf;
	  short* tibuf= (short*)ibuf;
	  for (i=0; i<collective_blksize; i++) *tobuf++= *tibuf++;
        }
        break;
      case MRI_INT:
        {
	  int* tobuf= (int*)obuf;
	  int* tibuf= (int*)ibuf;
	  for (i=0; i<collective_blksize; i++) *tobuf++= *tibuf++;
        }
        break;
      case MRI_LONGLONG:
        {
	  long long* tobuf= (long long*)obuf;
	  long long* tibuf= (long long*)ibuf;
	  for (i=0; i<collective_blksize; i++) *tobuf++= *tibuf++;
        }
        break;
      case MRI_FLOAT:
        {
	  float* tobuf= (float*)obuf;
	  float* tibuf= (float*)ibuf;
	  for (i=0; i<collective_blksize; i++) *tobuf++= *tibuf++;	
        }
        break;
      case MRI_DOUBLE:
        {
	  double* tobuf= (double*)obuf;
	  double* tibuf= (double*)ibuf;
	  for (i=0; i<collective_blksize; i++) *tobuf++= *tibuf++;
        }
        break;
      }

      mri_set_chunk( Output, this_chunk, collective_blksize, out_offset, type,
	 	     obuf );
      if (reallyverbose_flg) 
        fprintf(stderr,"block: type %d, %lld at %lld -> %lld at %lld\n",
	        type, collective_blksize, in_offset[j],
	        collective_blksize, out_offset);

      in_offset[j] += collective_blksize;
      out_offset += collective_blksize; 
      transfer_total -= collective_blksize; 

      free(obuf);
    }
  }
}

static void add_to_chunk(char* this_chunk) {
  char key_buf[KEYBUF_SIZE];
  char* dimstr;
  int i=0;
  long long fast_blksize;
  long long slow_blksize;

  range = 0;

  if (reallyverbose_flg) 
     fprintf(stderr,"adding chunk <%s> \n",this_chunk);
  safe_copy(key_buf, this_chunk);
  safe_concat(key_buf, ".dimensions");
  if (mri_has(Input[0],key_buf)) {
    dimstr= mri_get_string(Input[0],key_buf);
    if (reallyverbose_flg) fprintf(stderr,"dimstr %s\n",dimstr);
    if (!strchr(dimstr,*selected_dim))
      Abort("%s: chunk <%s> in first input has no dimension %c!\n",
	    progname,this_chunk,*selected_dim);
    check_order(this_chunk, dimstr);

    /* get range for this chunk */

    while ( i<n_input_files ) {

      if (strchr(dimstr,*selected_dim)) {
        safe_copy(key_buf, this_chunk);
        safe_concat(key_buf, ".extent.");
        safe_concat(key_buf, selected_dim);
        selected_extent[i]= mri_get_int(Input[i], key_buf);
        if (reallyverbose_flg) 
          fprintf(stderr,"selected dim extent on input %d is %d \n",i+1, selected_extent[i]);
        }
      range += selected_extent[i];
      i++;
    }

    data_changed= 1;
    mri_set_int(Output, key_buf, range);
    calc_sizes(this_chunk, dimstr, &fast_blksize, &slow_blksize);
    if (first_dim)
      transfer_front_data(this_chunk, slow_blksize);
    else if (last_dim)
      transfer_end_data(this_chunk, fast_blksize);
    else
      transfer_data(this_chunk, fast_blksize, slow_blksize);
  }
}

int
main( argc, argv ) 
     int argc;
     char **argv;
{

  char infile[512], outfile[512], junkfile[512];
  char* this_key;
  char this_chunk[KEYBUF_SIZE];
  char key_buf[KEYBUF_SIZE];
  int i,j;

  progname= argv[0];

  /* Print version number */
  if (verbose_flg) Message( "# %s\n", rcsid );

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Get filenames */

  if (cl_present("o"))
    Abort ("Option o(outfile) has been replaced by outfile|out.  Please see help file. \n");

  cl_get("dimension|dim|d", "%option %s[t]",selected_dim);
  if (strlen(selected_dim)>1) {
    fprintf(stderr,"%s: Selected dim name must be 1 char long.\n",progname);
    Help("usage");
    exit(-1);
  }
#ifdef never  
  /* I think removing this conversion allows the program to distinguish
   * between uppercase and lowercase dims.
   */
  *selected_dim= tolower(*selected_dim);
#endif

  cl_get("outfile|out", "%option %s[%]", DEFAULT_OUT_NAME, outfile);

  while (cl_get("", "%s", infile)) {
    if (n_input_files>=MAX_INPUT_FILES) {
      Abort("%s: too many input files; limit of %d compiled in!\n",
	    argv[0],MAX_INPUT_FILES);
    }
    if( !strcmp( infile, outfile )) {
      Abort( "Input and output files must be distinct." );}
    Input[n_input_files] = mri_open_dataset( infile, MRI_READ );
    n_input_files++;
  }

  if (n_input_files<=1) {
    fprintf(stderr,"%s: At least two input files required.\n",progname);
    Help( "usage" );
    exit(-1);
  }
  

  verbose_flg= cl_present("v|verbose");
  reallyverbose_flg= cl_present("V|VERBOSE");
  if (reallyverbose_flg) verbose_flg= 1;
  if (cl_cleanup_check()) {
    int i;
    fprintf(stderr,"%s: invalid argument in command line:\n    ",progname);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(-1);
  }

  /*** End command-line parsing ***/
  

  /* Open input and output datasets */

  Output = mri_copy_dataset( outfile, Input[0] );
  hist_add_cl( Output, argc, argv );

  /* Walk through the first input handling each chunk in turn */
  mri_iterate_over_keys(Input[0]);
  while ((this_key= mri_next_key(Input[0])) != NULL) {
    if (strlen(this_key)>=KEYBUF_SIZE) 
      Abort("%s: key too long!\n",progname);
    if (!strcmp(mri_get_string(Input[0],this_key),"[chunk]")) {
      safe_copy(this_chunk, this_key);
      safe_copy(key_buf, this_chunk);
      safe_concat(key_buf,".dimensions");
      if (!mri_has(Input[0],key_buf))
	Abort("%s: first input has no %s tag!\n",progname,key_buf);
      if (strchr(mri_get_string(Input[0],key_buf),*selected_dim)) {
	for (j=1; j<n_input_files; j++)
	  if (!mri_has(Input[j], this_chunk))
	    Abort("%s: input file %d does not have chunk <%s>\n", 
		  progname, j+1, this_chunk);
	add_to_chunk(this_chunk);
      }
      else {
	/* Ignore this chunk; it lacks the selected dimension.
	 * The copy from Input[0] has already been transcribed 
	 * to the output.
	 */
      }
    }
  }

  /* Write and close data-sets */

  for (i=0; i<n_input_files; i++) mri_close_dataset(Input[i]);
  if (reallyverbose_flg)
    fprintf(stderr, "%s: input files closed.\n", progname);
  mri_close_dataset( Output );
  
  if (verbose_flg) Message( "#      Pasting complete.\n" );
  if (!data_changed) 
    Message("#      Warning: first input and output datasets identical!\n");

  return 0;
}

