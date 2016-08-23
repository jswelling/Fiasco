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
#include "slist.h"

static char rcsid[] = "$Id: mri_destroy_dataset.c,v 1.3 2006/07/20 21:52:00 welling Exp $";

int main (argc, argv)
     int argc;
     char **argv;
{
  char infile[512];
  MRI_Dataset *Input;
  int verbose_flg= 0;
  SList* filesToDelete= slist_create();

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /*** Parse command line ***/

  cl_scan( argc, argv );

  verbose_flg= cl_present("verbose|ver|v");
  while (cl_get("","%s",infile)) {
    slist_append(filesToDelete,strdup(infile));
  }
  if (slist_count(filesToDelete)==0) {
    fprintf(stderr,"%s: Name of dataset(s) to delete not given.\n",argv[0]);
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

  /* Open and delete the datasets.  First we open it to read to verify its
   * existence, then we open it to modify to verify we're allowed to 
   * destroy it.
   */
  mri_set_error_handling(MRI_IGNORE_ERRORS);
  while (!slist_empty(filesToDelete)) {
    char* fname= (char*)slist_pop(filesToDelete);
    Input = mri_open_dataset( fname, MRI_READ );
    if (Input == NULL)
      Abort("%s: dataset %s does not exist or is not readable!\n",
	    argv[0],fname);
    mri_close_dataset(Input);
    Input = mri_open_dataset( fname, MRI_MODIFY );
    if (Input == NULL)
      Abort("%s: cannot open dataset %s for writing!\n",
	    argv[0],fname);
    mri_set_error_handling(MRI_ABORT_ON_ERROR);
    mri_destroy_dataset( Input );
    if (verbose_flg) Message( "# Destroyed %s\n", fname );
    free(fname);
  }

  /* And that's is. */
  return 0;
}
