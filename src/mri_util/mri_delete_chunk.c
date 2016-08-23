/*
 *	mri_delete_chunk.c - Delete a chunk from a dataset
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
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mri.h"
#include "stdcrg.h"
#include "fmri.h"

static char rcsid[] = "$Id: mri_delete_chunk.c,v 1.5 2003/09/15 23:27:43 welling Exp $";

int main (argc, argv)
     int argc;
     char **argv;
{
  char file[512];
  char chunk[MRI_MAX_KEY_LENGTH+1];
  MRI_Dataset *ds;
  char *key;
  int len;

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /*** Parse command line and get input params ***/

  file[0] = '\0';
  chunk[0] = '\0';
  cl_scan( argc, argv );

  /* Deprecate old options */

  if (cl_present( "dataset|d" ))
     Abort ("Option dataset|d has been replaced by infile format.  Please see help file.\n");

  cl_get( "chunk|chu|c", "%option %s", chunk);

  if(!cl_get("", "%s", file)) {
    fprintf(stderr, "%s: Input file name not given.\n", argv[0]);
    exit(-1);
  }

  if (file[0] == '\0' || chunk[0] == '\0')
    {
      fprintf(stderr,"%s: Both dataset and chunk name must be specified\n",
	      argv[0]);
      Help( "usage" );
      exit(1);
    }
  if (cl_cleanup_check()) {
    int i;
    fprintf(stderr,"%s: invalid argument in command line:\n    ", argv[0]);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(1);
  }
  /*** End command-line parsing ***/

  ds = mri_open_dataset(file, MRI_MODIFY);
  if (ds == NULL)
    {
      fprintf(stderr, "mri_delete_chunk: cannot open dataset: %s\n", file);
      exit(1);
    }
  hist_add_cl(ds,argc,argv);

  /* check if the chunk exists in the input dataset */
  if (!mri_has(ds, chunk) ||
      strcmp(mri_get_string(ds, chunk), "[chunk]") != 0)
    {
      Message("mri_delete_chunk: Dataset %s has no chunk named %s\n",
	      file, chunk);
      mri_close_dataset(ds);
      exit(0);
    }

  /* remove all keys associated with the chunk */
  len = strlen(chunk);
  mri_iterate_over_keys(ds);
  while ((key = mri_next_key(ds)) != NULL)
    if (strncmp(key, chunk, len) == 0 &&
	(key[len] == '\0' || key[len] == '.'))
      mri_remove(ds, key);
  
  mri_close_dataset(ds);
  exit(0);
}
