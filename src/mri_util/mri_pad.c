/************************************************************
 *                                                          *
 *  mri_pad.c                                               *
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

static char rcsid[] = "$Id: mri_pad.c,v 1.13 2007/07/27 17:42:07 welling Exp $";

static MRI_Dataset *Input = NULL, *Output = NULL;
static char selected_dim[512]= "t";
static int offset;
static int range;
static int data_changed= 0;
static float fillval= 0.0;
static char* progname;
static int verbose_flg= 0;
static int reallyverbose_flg= 0;

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

static void transfer_data(char* this_chunk, 
			  int fast_blksize, int slow_blksize, 
			  int selected_extent)
{
  int in_offset= 0;
  int out_offset= 0;
  int ifast;
  int islow;
  int collective_blksize;
  int out_blksize;
  int type;
  int typesize;
  void* obuf;
  void* ibuf;
  int i;

  type= get_chunk_type(Input,this_chunk);
  typesize= get_typesize(Input,this_chunk);

  collective_blksize= fast_blksize * selected_extent;
  out_blksize= fast_blksize * range;
  if (!(obuf= malloc(out_blksize*typesize))) {
    Abort("%s: unable to allocate %d bytes!\n",progname,
	  out_blksize*typesize);
  }
  
  for (islow=0; islow<slow_blksize; islow++) {

    ibuf= mri_get_chunk(Input, this_chunk, collective_blksize, 
			in_offset, type);

    switch (type) {
    case MRI_UNSIGNED_CHAR:
      {
	char* tobuf= (char*)obuf;
	char* tibuf= (char*)ibuf;
	int tail= range - (offset+selected_extent);
	for (i=0; i<offset*fast_blksize; i++) *tobuf++= (char)fillval;
	for (i=0; i<selected_extent*fast_blksize; i++) *tobuf++= *tibuf++;
	for (i=0; i<tail*fast_blksize; i++) *tobuf++= (char)fillval;
      }
      break;
    case MRI_SHORT:
      {
	short* tobuf= (short*)obuf;
	short* tibuf= (short*)ibuf;
	int tail= range - (offset+selected_extent);
	for (i=0; i<offset*fast_blksize; i++) *tobuf++= (short)fillval;
	for (i=0; i<selected_extent*fast_blksize; i++) *tobuf++= *tibuf++;
	for (i=0; i<tail*fast_blksize; i++) *tobuf++= (short)fillval;
      }
      break;
    case MRI_INT:
      {
	int* tobuf= (int*)obuf;
	int* tibuf= (int*)ibuf;
	int tail= range - (offset+selected_extent);
	for (i=0; i<offset*fast_blksize; i++) *tobuf++= (int)fillval;
	for (i=0; i<selected_extent*fast_blksize; i++) *tobuf++= *tibuf++;
	for (i=0; i<tail*fast_blksize; i++) *tobuf++= (int)fillval;
      }
      break;
    case MRI_LONGLONG:
      {
	long long* tobuf= (long long*)obuf;
	long long* tibuf= (long long*)ibuf;
	int tail= range - (offset+selected_extent);
	for (i=0; i<offset*fast_blksize; i++) *tobuf++= (long long)fillval;
	for (i=0; i<selected_extent*fast_blksize; i++) *tobuf++= *tibuf++;
	for (i=0; i<tail*fast_blksize; i++) *tobuf++= (long long)fillval;
      }
      break;
    case MRI_FLOAT:
      {
	float* tobuf= (float*)obuf;
	float* tibuf= (float*)ibuf;
	int tail= range - (offset+selected_extent);
	for (i=0; i<offset*fast_blksize; i++) *tobuf++= (float)fillval;
	for (i=0; i<selected_extent*fast_blksize; i++) *tobuf++= *tibuf++;
	for (i=0; i<tail*fast_blksize; i++) *tobuf++= (float)fillval;
      }
      break;
    case MRI_DOUBLE:
      {
	double* tobuf= (double*)obuf;
	double* tibuf= (double*)ibuf;
	int tail= range - (offset+selected_extent);
	for (i=0; i<offset*fast_blksize; i++) *tobuf++= (double)fillval;
	for (i=0; i<selected_extent*fast_blksize; i++) *tobuf++= *tibuf++;
	for (i=0; i<tail*fast_blksize; i++) *tobuf++= (double)fillval;
      }
      break;
    }

    mri_set_chunk( Output, this_chunk, out_blksize, out_offset, type,
		   obuf );

    if (reallyverbose_flg) 
      fprintf(stderr,"block: type %d, %d at %d -> %d at %d\n",
	      type, collective_blksize, in_offset,
	      out_blksize, out_offset);
    in_offset += collective_blksize;
    out_offset += out_blksize;
  }

  free(obuf);
}

static int padIndexRemap(int oldIndex, void* hook)
{
  int newIndex= oldIndex+offset;
  if ((newIndex<0) || (newIndex>=range)) return -1;
  else return newIndex;
}

static void pad_chunk(char* this_chunk) {
  char key_buf[KEYBUF_SIZE];
  char* dimstr;
  int selected_extent;

  if (reallyverbose_flg) fprintf(stderr,"padding chunk <%s>\n",this_chunk);
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
      if (selected_extent+offset>range)
	Abort("%s: chunk <%s> is too wide to pad!\n",progname,this_chunk);
      if (reallyverbose_flg) 
	fprintf(stderr,"selected dim extent on input %d\n",selected_extent);
      if (selected_extent<range) {
	int fast_blksize;
	int slow_blksize;

	data_changed= 1;
	mri_set_int(Output, key_buf, range);
	mriu_updateLabels(Output, this_chunk, selected_dim[0], 
			  0, selected_extent,
			  padIndexRemap, NULL);
	calc_sizes(this_chunk, dimstr, &fast_blksize, &slow_blksize);
	transfer_data(this_chunk, fast_blksize, slow_blksize, 
		      selected_extent);
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

  /* Get filenames */

  if (cl_present("r"))
    Abort ("Option r(range) has been replaced by length|len|l.  Please see help file.\n");

  if (cl_present("o"))
    Abort ("Option o(offset) has been replaced by shift|shi|s.  Please see help file.\n");
  
  if (cl_present( "f" ))
     Abort ("Option f has been expanded to fillvalue|fil .  Please see help file.\n");

  cl_get("dimension|dim|d", "%option %s[t]",selected_dim);
  if (strlen(selected_dim)>1) {
    fprintf(stderr,"%s: Selected dim name must be 1 char long.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  *selected_dim= tolower(*selected_dim);

  if (!cl_get("length|len|l", "%option %d",&range)) {
    fprintf(stderr,"%s: Index range not given.\n",argv[0]);
    Help("usage");
    exit(-1);
  }

  cl_get("shift|shi|s", "%option %d[0]",&offset);

  cl_get("fillvalue|fil", "%option %f[0.0]",&fillval);

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
  reallyverbose_flg= cl_present("V");
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
      pad_chunk(this_chunk);
    }
  }

  /* Write and close data-sets */
  mri_close_dataset( Input );
  mri_close_dataset( Output );
  
  if (verbose_flg) Message( "#      Padding complete.\n" );
  if (!data_changed) 
    Message("#      Warning: input and output datasets identical!\n");

  return 0;
}

