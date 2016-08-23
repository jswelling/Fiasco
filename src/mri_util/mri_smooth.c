/************************************************************
 *                                                          *
 *  mri_smooth.c                                     *
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

  DESCRIPTION OF MRI_SMOOTH

  mri_smooth takes a pgh MRI dataset of any type, and outputs
  a dataset of the same type containing a smoothed version of 
  the data.

**************************************************************/

#include <stdlib.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "misc.h"
#include "stdcrg.h"

#define KEYBUF_SIZE 512

static char rcsid[] = "$Id: mri_smooth.c,v 1.12 2007/03/21 23:55:39 welling Exp $";

typedef enum 
{ 
  MISSING_IGNORE, 
  MISSING_Z_SIMPLE, MISSING_Z_TFAST, MISSING_Z_TSLOW,
  MISSING_T_SIMPLE, MISSING_T_ZFAST, MISSING_T_ZSLOW 
} MissingCase;

static char* missingCaseNames[]= 
{
  "MISSING_IGNORE", 
  "MISSING_Z_SIMPLE", "MISSING_Z_TFAST", "MISSING_Z_TSLOW",
  "MISSING_T_SIMPLE", "MISSING_T_ZFAST", "MISSING_T_ZSLOW"
};

static MRI_Dataset *Input = NULL, *Output = NULL;
static char selected_dim[512]= "t";
static int data_changed= 0;
static char* progname;
static int verbose_flg= 0;
static int debug= 0;
static unsigned char **missing = NULL;

static void safe_copy(char* str1, const char* str2) {
  strncpy(str1, str2, KEYBUF_SIZE);
  str1[KEYBUF_SIZE-1]= '\0';
}

static void safe_concat(char* str1, char* str2) {
  strncat(str1, str2, (KEYBUF_SIZE-strlen(str1))-1);
}

static int safe_get_extent(MRI_Dataset* ds, const char* chunk, const char dim)
{
  char key_buf[KEYBUF_SIZE];
  char dim_buf[4];
  dim_buf[0]= dim;
  dim_buf[1]= '\0';
  safe_copy(key_buf,chunk);
  safe_concat(key_buf,".extent.");
  safe_concat(key_buf,dim_buf);
  if (mri_has(Input,key_buf)) return mri_get_int(Input,key_buf);
  else {
    Abort("%s: input missing tag %s!\n",progname,key_buf);
  }
  return 0; /* not reached */
}

static void calc_sizes(const char* this_chunk, const char* dimstr, const char liveDim,
		       long long* fast_blocksize_out, long long* slow_blocksize_out ) {
  /* This routine will fail if selected_dim is not in dimstr! */
  long long fast_blocksize;
  long long slow_blocksize;
  const char* this_dim;

  fast_blocksize= 1;

  this_dim= dimstr;
  while (*this_dim != liveDim) 
    fast_blocksize *= safe_get_extent(Input,this_chunk,*this_dim++);

  this_dim++; /* step over selected dim */

  slow_blocksize= 1;
  while (*this_dim)
    slow_blocksize *= safe_get_extent(Input, this_chunk, *this_dim++);

  *fast_blocksize_out= fast_blocksize;
  *slow_blocksize_out= slow_blocksize;
}

static void transfer_data(char* this_chunk, 
			  long long fast_blksize, long long slow_blksize, 
			  int selected_extent, Smoother* smoother,
			  MissingCase missingCase, long long missingScale,
			  long long missingMod)
{
  long long in_offset= 0;
  long long out_offset= 0;
  long long ifast;
  long long islow;
  long long i;
  long long collective_blksize;
  int type;
  float* ibuf= NULL;
  float* obuf= NULL;
  float* oblock= NULL;

  if (!(ibuf=(float*)malloc(selected_extent*sizeof(float))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname, selected_extent*sizeof(float));
  if (!(obuf=(float*)malloc(selected_extent*sizeof(float))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname, selected_extent*sizeof(float));
  if (!(oblock=(float*)malloc(fast_blksize*selected_extent*sizeof(float))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname, fast_blksize*selected_extent*sizeof(float));

  collective_blksize= fast_blksize * selected_extent;
  for (islow=0; islow<slow_blksize; islow++) {
    float* block= mri_get_chunk(Input, this_chunk, collective_blksize,
				in_offset, MRI_FLOAT);
    switch (missingCase) {
    case MISSING_T_SIMPLE: 
    case MISSING_Z_SIMPLE: 
      {
	for (ifast=0; ifast<fast_blksize; ifast++) {
	  for (i=0; i<selected_extent; i++) 
	    ibuf[i]= block[ifast + (i*fast_blksize)];
	  SM_SMOOTH(smoother, ibuf, obuf, selected_extent, 
		    missing, 0);
	  for (i=0; i<selected_extent; i++)
	    oblock[ifast + (i*fast_blksize)]= obuf[i];
	}
      }
      break;
    case MISSING_T_ZSLOW: 
    case MISSING_Z_TSLOW: 
      {
	for (ifast=0; ifast<fast_blksize; ifast++) {
	  for (i=0; i<selected_extent; i++) 
	    ibuf[i]= block[ifast + (i*fast_blksize)];
	  SM_SMOOTH(smoother, ibuf, obuf, selected_extent, 
		    missing, (int)((islow/missingScale)%missingMod));
	  for (i=0; i<selected_extent; i++)
	    oblock[ifast + (i*fast_blksize)]= obuf[i];
	}
      }
      break;
    case MISSING_T_ZFAST: 
    case MISSING_Z_TFAST: 
      {
	for (ifast=0; ifast<fast_blksize; ifast++) {
	  for (i=0; i<selected_extent; i++) 
	    ibuf[i]= block[ifast + (i*fast_blksize)];
	  SM_SMOOTH(smoother, ibuf, obuf, selected_extent, 
		    missing, (int)((ifast/missingScale)%missingMod));
	  for (i=0; i<selected_extent; i++)
	    oblock[ifast + (i*fast_blksize)]= obuf[i];
	}
      }
      break;
    case MISSING_IGNORE: 
      {
	for (ifast=0; ifast<fast_blksize; ifast++) {
	  for (i=0; i<selected_extent; i++) 
	    ibuf[i]= block[ifast + (i*fast_blksize)];
	  SM_SMOOTH(smoother, ibuf, obuf, selected_extent, NULL, 0);
	  for (i=0; i<selected_extent; i++)
	    oblock[ifast + (i*fast_blksize)]= obuf[i];
	}
      }
      break;
    default: 
      Abort("%s: internal error; unknown missingCase %d!\n",
	    progname,(int)missingCase);
    }
    mri_set_chunk( Output, this_chunk, collective_blksize, out_offset,
		   MRI_FLOAT, oblock );
    if (debug) 
      fprintf(stderr,"block: %lld at %lld -> %lld at %lld\n",
	      collective_blksize, in_offset,
	      collective_blksize, out_offset);
    in_offset += collective_blksize;
    out_offset += collective_blksize;
  }

  free(ibuf);
  free(obuf);
  free(oblock);
}

static void smooth_chunk(char* this_chunk, Smoother* smoother) {
  char key_buf[KEYBUF_SIZE];
  char* dimstr;
  int selected_extent;

  if (verbose_flg) fprintf(stderr,"smoothing chunk <%s>\n",this_chunk);
  safe_copy(key_buf, this_chunk);
  safe_concat(key_buf, ".dimensions");
  if (mri_has(Input,key_buf)) {
    dimstr= mri_get_string(Input,key_buf);
    if (verbose_flg) fprintf(stderr,"dimstr %s\n",dimstr);
    if (strchr(dimstr,*selected_dim)) {
      selected_extent= safe_get_extent(Input,this_chunk,*selected_dim);
      if (verbose_flg) 
	fprintf(stderr,"selected dim extent on input %d\n",selected_extent);
      if (selected_extent>1) {
	long long fast_blksize;
	long long slow_blksize;
	MissingCase missingCase;
	long long missingScale= 1;
	long long missingMod;
	long long dummy;

	data_changed= 1;
	calc_sizes(this_chunk, dimstr, *selected_dim, &fast_blksize, &slow_blksize);
	if (*selected_dim == 't') {
	  char* zptr= strchr(dimstr,'z');
	  char* tptr= strchr(dimstr,'t');
	  if (zptr) {
	    missingMod= safe_get_extent(Input,this_chunk,'z');
	    if (zptr < tptr) {
	      missingCase= MISSING_T_ZFAST;
	      calc_sizes(this_chunk, dimstr, 't', &missingScale, &dummy);
	    }
	    else {
	      missingCase= MISSING_T_ZSLOW;
	      calc_sizes(this_chunk, tptr+1, 'z', &missingScale, &dummy);
	    }
	  }
	  else {
	    missingCase= MISSING_T_SIMPLE;
	    missingScale= 1;
	    missingMod= 1;
	  }
	}
	else if (*selected_dim == 'z') {
	  char* zptr= strchr(dimstr,'z');
	  char* tptr= strchr(dimstr,'t');
	  if (tptr) {
	    missingMod= safe_get_extent(Input,this_chunk,'t');
	    if (tptr < zptr) {
	      missingCase= MISSING_Z_TFAST;
	      calc_sizes(this_chunk, dimstr, 't', &missingScale, &dummy);
	    }
	    else {
	      missingCase= MISSING_Z_TSLOW;
	      calc_sizes(this_chunk, zptr+1, 't', &missingScale, &dummy);
	    }
	  }
	  else {
	    missingCase= MISSING_Z_SIMPLE;
	    missingScale= 1;
	    missingMod= 1;
	  }
	}
	else {
	  missingCase= MISSING_IGNORE;
	  missingScale= 1;
	  missingMod= 1;
	}
	if (debug) { 
	  fprintf(stderr,"fast_blksize= %lld, slow_blksize= %lld\n",
		  fast_blksize, slow_blksize);
	  fprintf(stderr,"missing processing case %s, factor %lld, mod %lld\n",
		  missingCaseNames[(int)missingCase],missingScale,missingMod);
	}
	safe_copy(key_buf, this_chunk);
	safe_concat(key_buf, ".type");
	mri_set_string(Output, key_buf, "float32");
	transfer_data(this_chunk, fast_blksize, slow_blksize, 
		      selected_extent, smoother, missingCase, missingScale, missingMod);
      }
    }
    else {
      /* Chunk copied correctly in initial dataset copy */
    }
  }
}

int main( int argc, char* argv[] ) 
{
  char infile[512], outfile[512], kernel[512];
  char* this_key;
  char this_chunk[KEYBUF_SIZE];
  Smoother* smoother;
  sm_type smoother_type;
  float bandwidth;
  float threshold;

  progname= argv[0];

  sm_init();

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /*** Parse command line ***/

  cl_scan( argc, argv );

  if (cl_present( "b" ))
     Abort ("Option b(bandwidth) has been expanded to bandwidth|bdw.  Please see help file.\n");
  if (cl_present( "k" ))
     Abort ("Option k(kernel) has been expanded to kernel|ker.  Please see help file.\n");
  if (cl_present( "t" ))
     Abort ("Option t(threshold) has been expanded to threshold|thr.  Please see help file.\n");

  /* Get filenames */
  cl_get("dimension|dim|d", "%option %s[t]",selected_dim);
  if (strlen(selected_dim)>1) {
    fprintf(stderr,"%s: Selected dim name must be 1 char long.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  *selected_dim= tolower(*selected_dim);

  verbose_flg= cl_present("verbose|ver|v");
  debug= cl_present("debug");

  cl_get( "bandwidth|bdw", "%option %f[%]", 3.0, &bandwidth );
  cl_get( "kernel|ker", "%option %s[%]", "gaussian", kernel );
  cl_get( "threshold|thr", "%option %f[%]", 0.0, &threshold );

  if( !strcasecmp( kernel, "gaussian" ) || !strcasecmp( kernel, "gau" ) )
    {
      smoother_type= SM_GAUSSIAN;
    }
  else if( !strcasecmp( kernel, "triangular" ) || !strcasecmp( kernel, "tri" ) )
    {
      smoother_type= SM_TRIANGULAR;
    }
  else if( !strcasecmp( kernel, "exponential" ) || !strcasecmp( kernel, "exp" ) )
     {
       smoother_type= SM_POWER;
     }
  else if( !strcasecmp( kernel, "median" ) || !strcasecmp( kernel, "med" ) )
     {
       smoother_type= SM_MEDIAN;
     }
  else {
    Abort( "Kernel unrecognized (%s).\n", kernel );
  }

  sm_set_params( smoother_type, bandwidth, 0.0, threshold, NULL );
  sm_parse_cl_opts();
  sm_get_params( &smoother_type, &bandwidth, NULL, &threshold, NULL );

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
  
  if (verbose_flg) {
    /* Print version number */
    Message( "# %s\n", rcsid );

    Message("smoother type %d, bandwidth %f, threshold %f\n",
      smoother_type,bandwidth,threshold);
  }

  /* Open input and output datasets */
  if( !strcmp( infile, outfile ) )
    Abort( "Input and output files must be distinct." );
  Input = mri_open_dataset( infile, MRI_READ );
  Output = mri_copy_dataset( outfile, Input );
  hist_add_cl( Output, argc, argv );

  /* Create the smoother */
  smoother= sm_create_smoother();
  sm_set_direction(smoother,*selected_dim);

  /* Import missing information if appropriate */
  if (*selected_dim=='t' || *selected_dim=='z') {
    if (mri_has(Input,"missing") 
	&& !strcmp(mri_get_string(Input,"missing"), "[chunk]")) {
      if (debug) fprintf(stderr,"Respecting input file's missing chunk\n");
      missing= get_missing(Input);
    }
  }

  /* Walk through the input handling each chunk in turn */
  mri_iterate_over_keys(Input);
  while ((this_key= mri_next_key(Input)) != NULL) {
    if (strlen(this_key)>=KEYBUF_SIZE) 
      Abort("%s: key too long!\n",argv[0]);
    if (!strcmp(mri_get_string(Input,this_key),"[chunk]")) {
      strncpy(this_chunk, this_key, KEYBUF_SIZE);
      this_chunk[KEYBUF_SIZE-1]= '\0';
      if (strcmp(this_chunk,"missing")) /* don't smooth missing chunk! */
	smooth_chunk(this_chunk, smoother);
    }
  }

  /* Write and close data-sets */
  mri_close_dataset( Input );
  mri_close_dataset( Output );
  
  /* Clean up */
  sm_destroy(smoother);

  if (verbose_flg) Message( "#      Smoothing complete.\n" );
  if (!data_changed) 
    Message("#      Warning: input and output datasets identical!\n");

  return 0;
}

