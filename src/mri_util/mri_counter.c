/************************************************************
 *                                                          *
 *  mri_counter.c                                     *
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
 *  Original programming by Joel Welling, 8/97              *
 ************************************************************/
/*************************************************************

  DESCRIPTION OF MRI_COUNTER

  mri_counter creates a pgh MRI file with one dimension, and stores
  a sequence of integers in that dimension.  The integers simply count
  upward from 0 to the maximum extent given.  This program is useful
  for a variety of tricks in combination with other utilities like
  mri_rpn_math and mri_subsample.


**************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "misc.h"
#include "stdcrg.h"

#define CHUNKSIZE 4096
#define KEYBUF_SIZE 512

static char rcsid[] = "$Id: mri_counter.c,v 1.7 2003/08/07 19:46:22 bakalj Exp $";

static float out_chunk[CHUNKSIZE];
static int verbose_flg= 0;

static void safe_copy(char* str1, char* str2) {
  strncpy(str1, str2, KEYBUF_SIZE);
  str1[KEYBUF_SIZE-1]= '\0';
}

static void safe_concat(char* str1, char* str2) {
  strncat(str1, str2, (KEYBUF_SIZE-strlen(str1))-1);
}

int main( int argc, char* argv[] ) 
{
  MRI_Dataset *Output = NULL;
  char dimstr[KEYBUF_SIZE];
  int extent;
  int i;
  char outfile[512];
  char keybuf[KEYBUF_SIZE];
  char chunk[KEYBUF_SIZE];
  int* data= NULL;

  /* Print version number */
  if (verbose_flg) Message( "# %s\n", rcsid );

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Get filenames */
  verbose_flg= cl_present("verbose|ver|v");
  cl_get("chunk|chu|c", "%option %s[%]","images",chunk);
  cl_get("dimension|dim|d", "%option %s[%]","t",dimstr);
  if (strlen(dimstr)>1) {
    fprintf(stderr,"%s: -d switch must be followed by a single character.",
	    argv[0]);
    Help("usage");
    exit(-1);
  }
  if (!cl_get("length|len|l", "%option %d",&extent)) {
    fprintf(stderr,"%s: required extent switch is missing.",argv[0]);
    Help( "usage" );
    exit(-1);
  }
  if (!cl_get("", "%s", outfile)) {
    fprintf(stderr,"%s: Output file name not given.",argv[0]);
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
  
  /* Open output dataset, and write the keywords */
  Output= mri_open_dataset( outfile, MRI_WRITE );
  hist_add_cl(Output,argc,argv);
  mri_create_chunk(Output, chunk);
  safe_copy(keybuf,chunk);
  safe_concat(keybuf,".datatype");
  mri_set_string(Output,keybuf,"int32");
  safe_copy(keybuf,chunk);
  safe_concat(keybuf,".dimensions");
  mri_set_string(Output,keybuf,dimstr);
  safe_copy(keybuf,chunk);
  safe_concat(keybuf,".extent.");
  safe_concat(keybuf,dimstr);
  mri_set_int(Output,keybuf,extent);
  safe_copy(keybuf,chunk);
  safe_concat(keybuf,".file");
  mri_set_string(Output,keybuf,".dat");

  /* Build and set the data */
  if ( !(data=(int*)malloc(extent*sizeof(int))) )
    Abort("%s: unable to allocate %d ints!\n",argv[0],extent);
  for (i=0; i<extent; i++) data[i]= i;
  mri_set_chunk(Output, chunk, extent, 0, MRI_INT, data);

  /* Close the resulting dataset;  we're done */
  mri_close_dataset( Output );
  
  return 0;

}

