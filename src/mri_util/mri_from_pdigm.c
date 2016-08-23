/************************************************************
 *                                                          *
 *  mri_from_pdigm.c                                            *
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

  DESCRIPTION OF MRI_FROM_PDIGM

  mri_from_pdigm takes a stream of ascii data as input.  That 
  stream consists of multiple lines, each of which starts with a 
  monotonically increasing float (a time stamp).  This is the
  output format produced by the "dmread" utility when it is used
  with the "-s" option.  mri_from_pdigm reads the stream and
  produces a Pittsburgh MRI dataset in "vt" format containing
  the data.

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
#define INITIAL_GUESS_EXTENT 1024

static char rcsid[] = "$Id: mri_from_pdigm.c,v 1.7 2003/08/07 20:22:07 bakalj Exp $";

static MRI_Dataset *Output = NULL;
static char chunk[512];
static float tmin= 0.0;
static float tmax= 0.0;
char* progname;

static void safe_copy(char* str1, char* str2) {
  strncpy(str1, str2, KEYBUF_SIZE);
  str1[KEYBUF_SIZE-1]= '\0';
}

static void safe_concat(char* str1, char* str2) {
  strncat(str1, str2, (KEYBUF_SIZE-strlen(str1))-1);
}

static void parse_tokens( float* databuf, char* inbuf, int max )
{
  int i;
  char* tok;

  for (i=0; i<max; i++) {
    tok= strtok( (i ? NULL : inbuf), " \t\n" );
    if (!tok) Abort("%s: found short line!\n",progname);
    databuf[i]= atof(tok);
  }
  if (strtok(NULL,inbuf))
    Abort("%s: found long line!\n",progname);
}

static void parse_and_save_stream()
{
  int guess_extent;
  int real_extent= 0;
  int n_components= 0;
  char inbuf[512];
  char inbuf2[512];
  char keybuf[KEYBUF_SIZE];
  float* databuf;
  int offset= 0;
  int in_range_flag= 0;

  guess_extent= INITIAL_GUESS_EXTENT;

  while (!n_components) { /* keep trying 'till we find some data */
    fgets(inbuf,512,stdin);
    strcpy(inbuf2,inbuf); /* keep the original */
    
    /* count tokens */
    if (strtok(inbuf2," \t\n")) n_components++;
    while (strtok(NULL," \t\n")) n_components++;
  }

  safe_copy(keybuf,chunk);
  safe_concat(keybuf,".extent.v");
  mri_set_int(Output,keybuf,n_components);
  safe_copy(keybuf,chunk);
  safe_concat(keybuf,".extent.t");
  mri_set_int(Output,keybuf,guess_extent);

  if (!(databuf= (float*)malloc(n_components*sizeof(float))))
    Abort("%s: unable to allocate %d floats!\n",progname,n_components);

  parse_tokens(databuf,inbuf,n_components);
  if (databuf[0]>=tmin) {
    if (databuf[0]>tmax) Abort("%s: first time stamp too high!\n",
			       progname);
    else {
      in_range_flag= 1;
      mri_set_chunk(Output,chunk,n_components,offset,MRI_FLOAT,databuf);
      offset += n_components;
      real_extent++;
    }
  }

  while (!feof(stdin) && !ferror(stdin)) {
    if (!fgets(inbuf,512,stdin)) break;
    parse_tokens(databuf,inbuf,n_components);
    if (in_range_flag) {
      if (databuf[0]>tmax) break;
    }
    else {
      if (databuf[0]>=tmin) in_range_flag= 1;
    }
    if (in_range_flag) {
      mri_set_chunk(Output,chunk,n_components,offset,MRI_FLOAT,databuf);
      offset += n_components;
      real_extent++;
      if (real_extent==guess_extent) {
	guess_extent= 2*guess_extent;
	safe_copy(keybuf,chunk);
	safe_concat(keybuf,".extent.t");
	mri_set_int(Output,keybuf,guess_extent);
      }
    }
  }

  /* Shrink back down to right size */
  safe_copy(keybuf,chunk);
  safe_concat(keybuf,".extent.t");
  mri_set_int(Output,keybuf,real_extent);
}

int main( int argc, char* argv[] ) 
{

  char outfile[512];
  char keybuf[KEYBUF_SIZE];

  progname= argv[0];

  /* Print version number */
  Message( "# %s\n", rcsid );

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Deprecate old options */

  if (cl_present( "l" ))
     Abort ("Option l(low) has been expanded to low.  Please see help file.\n");
  if (cl_present( "h" ))
     Abort ("Option h(high) has been expanded to high.  Please see help file.\n");


  /* Get parameters */
  cl_get("low","%option %f[0.0]",&tmin);
  if (!cl_get("high","%option %f",&tmax)) {
    fprintf(stderr,"%s: maximum time stamp not given.\n",argv[0]);
    Help( "usage" );
    exit(-1);
  }

  cl_get("chunk|chu|c","%option %s[%]","physio",&chunk);

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
  
  /* Reality check */
  if (tmin>=tmax)
    Abort("%s: invalid timestamp range!\n",argv[0]);

  /* Open input and output datasets */
  Output = mri_open_dataset( outfile, MRI_WRITE );
  hist_add_cl(Output,argc,argv);
  mri_create_chunk(Output,chunk);
  safe_copy(keybuf,chunk);
  safe_concat(keybuf,".datatype");
  mri_set_string( Output, keybuf, "float32" );
  safe_copy(keybuf,chunk);
  safe_concat(keybuf,".dimensions");
  mri_set_string( Output, keybuf, "vt" );
  safe_copy(keybuf,chunk);
  safe_concat(keybuf,".file");
  mri_set_string( Output, keybuf, ".dat" );
  
  parse_and_save_stream();

  /* Write and close data-sets */
  mri_close_dataset( Output );
  
  return 0;
}

