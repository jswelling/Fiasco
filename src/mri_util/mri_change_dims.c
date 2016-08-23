/************************************************************
 *                                                          *
 *  mri_change_dims.c                                     *
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

  DESCRIPTION OF MRI_CHANGE_DIMS

  mri_change_dims changes the dimensions of one chunk in a pgh MRI
  file in a trivial way, adding or removing trivial dimensions.
  A dimension of extent 1 is added or removed from the dimensions
  for the given chunk.

  mri_change_dims [-c chunk] -d dimstr infile outfile

    -c chunk specifies the name of the chunk to change (default "images")
    -d dimstr specifies the dimension string that chunk is to have
       on output.
    infile specifies the input dataset
    outfile specifies the output dataset

**************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "misc.h"
#include "stdcrg.h"

#define KEYBUF_SIZE 512

static char rcsid[] = "$Id: mri_change_dims.c,v 1.6 2003/03/19 23:19:39 welling Exp $";

static MRI_Dataset *Input = NULL, *Output = NULL;
static char* progname;

static void safe_copy(char* str1, char* str2) {
  strncpy(str1, str2, KEYBUF_SIZE);
  str1[KEYBUF_SIZE-1]= '\0';
}

static void safe_concat(char* str1, char* str2) {
  strncat(str1, str2, (KEYBUF_SIZE-strlen(str1))-1);
}

static void convert_dimensions( MRI_Dataset* ds, char* chunk,
				char* dims_in, char* dims_out )
{
  char keybuf[KEYBUF_SIZE];
  char* in_runner= dims_in;
  char* out_runner= dims_out;
  char cbuf[2];
  char newdims[KEYBUF_SIZE];
  int new_length= 0;
  int i;

  if (strlen(dims_out)>KEYBUF_SIZE-1)
    Abort("%s: failed due to unreasonably long dim string!\n",progname);
  for (i=0; i<KEYBUF_SIZE; i++) newdims[i]= '\0';

  while (*in_runner || *out_runner) {
    if (*in_runner != *out_runner) {
      /* check for duplicate character */
      if (strchr(dims_out,*out_runner)<out_runner)
	Abort("%s: dimension %c is duplicated in output dim string!\n",
	      progname,*out_runner);
      /* if *out_runner exists farther along in *in_runner, its a deletion */
      if (strchr(in_runner,*out_runner)) {
	/* deletion */
	cbuf[0]= *in_runner;
	cbuf[1]= '\0';
	safe_copy(keybuf,chunk);
	safe_concat(keybuf,".extent.");
	safe_concat(keybuf,cbuf);
	if (!mri_has(ds,keybuf))
	  Abort("%s: input is missing %s tag!\n",progname,keybuf);
	if (mri_get_int(ds,keybuf) == 1) mri_remove(ds,keybuf);
	else {
	  Abort("%s: dimension %c has non-trivial extent; can't delete!\n",
		  progname,*in_runner);
	}
	in_runner++;
      }
      /* if *in_runner exists farther along in *out_runner, its an addition */
      else if (strchr(out_runner,*in_runner)) {
	/* addition */
	cbuf[0]= *out_runner;
	cbuf[1]= '\0';
	safe_copy(keybuf,chunk);
	safe_concat(keybuf,".extent.");
	safe_concat(keybuf,cbuf);
	mri_set_int(ds,keybuf,1);
	newdims[new_length++]= *out_runner;
	out_runner++;
      }
    }
    else {
      newdims[new_length++]= *out_runner;
      in_runner++;
      out_runner++;
    }
  }

  safe_copy(keybuf,chunk);
  safe_concat(keybuf,".dimensions");
  mri_set_string(ds,keybuf,newdims);
}

int main( int argc, char* argv[] ) 
{

  char infile[512], outfile[512];
  char keybuf[KEYBUF_SIZE];
  char chunkname[512];
  char dims_out[512];
  char* dims_in;

  progname= argv[0];

  /* Print version number */
  Message( "# %s\n", rcsid );
 
  
  /* Check to see if help was requested */
  if( ( argc > 1 ) && !strcmp( argv[1], "-help" ) )
    {
      if( argc == 2 )
	Help( "selecttopic" );
      else
	Help( argv[2] );
    }

 Abort("%s: mri_change_dims has been deprecated, please use mri_remap.\n", progname);

  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Get filenames */
  if (!cl_get("d", "%option %s",dims_out)) {
    fprintf(stderr,"%s: missing required flag with output dims.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  cl_get("c", "%options %s[%]", "images", chunkname);
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
  
  /* Open input and output datasets, checking for the requested chunk */
  if( !strcmp( infile, outfile ) )
    Abort( "Input and output files must be distinct." );
  Input = mri_open_dataset( infile, MRI_READ );
  if (!mri_has(Input,chunkname))
    Abort("%s: input has no chunk named %s!\n",argv[0],chunkname);
  Output = mri_copy_dataset( outfile, Input );
  hist_add_cl( Output, argc, argv );

  safe_copy(keybuf,chunkname);
  safe_concat(keybuf,".dimensions");
  if (mri_has(Output,keybuf)) {
    dims_in= mri_get_string(Output,keybuf);
    convert_dimensions(Output, chunkname, dims_in, dims_out);
  }
  else Abort("%s: Input file missing %s tag!\n",keybuf);

  /* Write and close data-sets */
  mri_close_dataset( Input );
  mri_close_dataset( Output );
  
  return 0;
}

