/*
 *	mri_copy_dataset.c - Copy a chunk from one dataset to another
 *
 *	This program copies one chunk from a given dataset to another.
 *
 *	Copyright (c) 1998 Pittsburgh Supercomputing Center
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
 *
 *	HISTORY:
 *		4/98 - Written by Greg Hood (PSC)
 *              10/01 - Adapted from mri_copy_chunk by Joel Welling
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mri.h"
#include "stdcrg.h"
#include "fmri.h"

static char rcsid[] = "$Id: mri_copy_dataset.c,v 1.3 2003/10/29 22:42:57 welling Exp $";

int main (argc, argv)
     int argc;
     char **argv;
{
  char infile[512];
  char outfile[512];
  MRI_Dataset *Input;
  MRI_Dataset *Output;
  int verbose_flg= 0;

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /*** Parse command line ***/

  cl_scan( argc, argv );

  verbose_flg= cl_present("verbose|ver|v");
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
  
  /* Print version number */
  if (verbose_flg) Message( "# %s\n", rcsid );

  /* Open input and copy to output */
  if( !strcmp( infile, outfile ) )
    Abort( "Input and output files must be distinct!\n" );
  Input = mri_open_dataset( infile, MRI_READ );
  if (Input == NULL)
    Abort("%s: cannot open input dataset %s!\n",argv[0],infile);
  Output = mri_copy_dataset( outfile, Input );
  if (Output == NULL)
    Abort("%s: cannot open output dataset %s!\n",argv[0],outfile);
  hist_add_cl( Output, argc, argv );

  /* And that's is. */
  mri_close_dataset(Input);
  mri_close_dataset(Output);
  if (verbose_flg) Message( "# Copy complete.\n" );
  return 0;
}
