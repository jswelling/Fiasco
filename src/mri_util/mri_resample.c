/************************************************************
 *                                                          *
 *  mri_resample.c                                            *
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

static char rcsid[] = "$Id: mri_resample.c,v 1.4 2005/03/09 01:13:09 welling Exp $";

static MRI_Dataset *Input = NULL, *Output = NULL;
static char selected_dim[512]= "";
static int new_extent;
static int data_changed= 0;
static char* progname;
static int const_flag= 0;
static int verbose_flag= 0;
static int debug_flag= 0;
static InterpolatorType mode= INTRP_LINEAR;

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

  this_dim++; /* step over selected dim */

  slow_blocksize= 1;
  while (*this_dim)
    slow_blocksize *= safe_get_extent(Input, this_chunk, this_dim++);

  *fast_blocksize_out= fast_blocksize;
  *slow_blocksize_out= slow_blocksize;
}

static void transfer_data(char* this_chunk, 
			  long long fast_blksize, long long slow_blksize, 
			  int selected_extent, int new_extent,
			  double start, double end)
{
  long long in_offset= 0;
  long long out_offset= 0;
  long ifast;
  long long islow;
  double* ibuf= NULL;
  double* obuf= NULL;
  int missing_chunk= 0;
  Interpolator* interpolator= 
    intrp_createInterpolator1DByType(mode, selected_extent, fast_blksize);

  if (debug_flag) {
    interpolator->setInt(interpolator,INTRP_OPT_DEBUG,1);
    interpolator->dumpSelf(interpolator,stderr);
  }

  missing_chunk= !strcmp(this_chunk,"missing");
  if (missing_chunk) 
    Abort("%s: I can't handle the missing chunk!\n",progname);

  if (!(obuf= (double*)malloc(fast_blksize*new_extent*sizeof(double))))
    Abort("%s: unable to allocate %d bytes!\n",progname,
	  fast_blksize*new_extent*sizeof(double));

  for (islow=0; islow<slow_blksize; islow++) {
    double shift= start;
    double stride= (end-start)/(double)(new_extent-1);
    long frame_in= 0;
    long frame_out= 0;
    int i;
    if (debug_flag) 
      fprintf(stderr,
	      "linear block, shift %f, stride %f, %lld * %d at %lld -> %lld * %d at %lld\n",
	      shift, stride, fast_blksize, selected_extent, 
	      in_offset, fast_blksize, new_extent, out_offset);
    ibuf= mri_get_chunk(Input, this_chunk, fast_blksize*selected_extent, 
			     in_offset, MRI_DOUBLE);
    interpolator->prep(interpolator, ibuf, fast_blksize*selected_extent);
    for (i=0; i<new_extent; i++) {
      if (shift<=0.0) {
	for (ifast=0; ifast<fast_blksize; ifast++) {
	  obuf[frame_out+ifast]= ibuf[ifast];
        }
      }
      else if (shift>=(double)(selected_extent-1)) {
	for (ifast=0; ifast<fast_blksize; ifast++) {
	  obuf[frame_out+ifast]= 
	    ibuf[(fast_blksize*(selected_extent-1))+ifast];
          }
      }
      else {
	double offset= shift-floor(shift);
	frame_in= (long)floor(shift)*fast_blksize;
	interpolator->calc(interpolator,obuf+frame_out, &shift, 
			   fast_blksize,0);
#ifdef never
	for (ifast=0; ifast<fast_blksize; ifast++) {
	  obuf[frame_out+ifast]=
	    (1.0-offset)*ibuf[frame_in+ifast] 
	    + offset*ibuf[frame_in+fast_blksize+ifast];
	}
#endif
      }
      shift += stride;
      frame_in += fast_blksize;
      frame_out += fast_blksize;
    }
    mri_set_chunk(Output, this_chunk, fast_blksize*new_extent, 
		  out_offset, MRI_DOUBLE, obuf);
    in_offset += fast_blksize*selected_extent;
    out_offset += fast_blksize*new_extent;
  }
  free(obuf);
  interpolator->destroySelf(interpolator);
}

static void resample_chunk(char* this_chunk, double start, double end) {
  char key_buf[KEYBUF_SIZE];
  char* dimstr;
  int selected_extent;

  if (debug_flag) fprintf(stderr,"resampling chunk <%s>\n",this_chunk);
  safe_copy(key_buf, this_chunk);
  safe_concat(key_buf, ".dimensions");
  if (mri_has(Input,key_buf)) {
    dimstr= mri_get_string(Input,key_buf);
    if (debug_flag) fprintf(stderr,"dimstr %s\n",dimstr);
    if (strchr(dimstr,*selected_dim)) {
      safe_copy(key_buf, this_chunk);
      safe_concat(key_buf, ".extent.");
      safe_concat(key_buf, selected_dim);
      if (mri_has(Input,key_buf)) 
	selected_extent= mri_get_int(Input,key_buf);
      else Abort("%s: input missing tag %s!\n",progname,key_buf);
      if (debug_flag) 
	fprintf(stderr,"selected dim extent on input %d\n",selected_extent);
      if (selected_extent<2) {
	/* resampling a constant- chunk was copied correctly in initial
	 * dataset copy.
	 */
	return;
      }
      else {
	long long fast_blksize;
	long long slow_blksize;

	data_changed= 1;
	mri_set_int(Output, key_buf, new_extent);
	calc_sizes(this_chunk, dimstr, &fast_blksize, &slow_blksize);
	transfer_data(this_chunk, fast_blksize, slow_blksize, 
		      selected_extent, new_extent, start, end);
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
  char modeString[512];
  char* this_key;
  char this_chunk[KEYBUF_SIZE];
  char keybuf[KEYBUF_SIZE];
  double start;
  double end;

  progname= argv[0];

  /* Print version number */
  if (verbose_flag) Message( "# %s\n", rcsid );

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /*** Parse command line ***/

  cl_scan( argc, argv );

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

  if (!cl_get("start", "%option %lf",&start)) {
    fprintf(stderr,"%s: start point not given.\n",argv[0]);
    Help("usage");
    exit(-1);
  }

  if (!cl_get("end", "%option %lf",&end)) {
    fprintf(stderr,"%s: start point not given.\n",argv[0]);
    Help("usage");
    exit(-1);
  }

  if (cl_get("interp", "%option %s", modeString)) {
    mode= intrp_typeFromName(modeString);
    if (mode==INTRP_UNKNOWN) {
      fprintf(stderr,"%s: unrecognized interpolation mode <%s>.\n",
	      argv[0],modeString);
      Help("usage");
      exit(-1);
    }
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
  verbose_flag= cl_present("verbose|ver|v");
  debug_flag= cl_present("debug");
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

      safe_copy(keybuf,this_chunk);
      safe_concat(keybuf,".datatype");
      if (!mri_has(Output, keybuf)
	  || strncmp( mri_get_string( Output, keybuf ), "float", 5 )) {
	/* We will coerce the output datatype to 4-byte floats if it is
	 * not one of the floating types.
	 */
	mri_set_string( Output, keybuf, "float32" );
      }

      resample_chunk(this_chunk, start, end);
    }
  }

  /* Write and close data-sets */
  mri_close_dataset( Input );
  mri_close_dataset( Output );

  if (verbose_flag) Message( "#      Resampling complete.\n" );
  if (!data_changed) 
    Message("#      Warning: input and output datasets identical!\n");

  return 0;
}

