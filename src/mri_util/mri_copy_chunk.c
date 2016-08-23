/*
 *	mri_copy_chunk.c - Copy a chunk from one dataset to another
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
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mri.h"
#include "stdcrg.h"
#include "fmri.h"

static char rcsid[] = "$Id: mri_copy_chunk.c,v 1.5 2003/08/06 20:51:01 bakalj Exp $";

#define BLOCK_SIZE	(16*65536)	/* maximum # of bytes to copy
					   at one time */

int main (argc, argv)
     int argc;
     char **argv;
{
  char in_file[512];
  char out_file[512];
  char in_chunk[MRI_MAX_KEY_LENGTH+1];
  char out_chunk[MRI_MAX_KEY_LENGTH+1];
  MRI_Dataset *ids;
  MRI_Dataset *ods;
  int replace;
  int in_size, out_size, size;
  void *buf;
  char *key;
  char new_key[MRI_MAX_KEY_LENGTH+1];
  int ilen, olen;
  int i;
  int offset;
  char out_chunk_file[MRI_MAX_FILENAME_LENGTH+1];

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /*** Parse command line and get input params ***/
  cl_scan( argc, argv );

  /* Deprecate old options */

  if (cl_present( "input|i" ))
    Abort ("Option input|i has been replaced by infile outfile format.  Please see help file.\n");
  if (cl_present( "output|o" ))
    Abort ("Option output|o has been replaced by infile outfile format.  Please see help file.\n");
  if (cl_present( "to|t" ))
    Abort ("Option to|t has been replaced by chunk_out|cho.  Please see help file.\n");
  if (cl_present( "file|f" ))
    Abort ("Option file|f has been replaced by infile outfile format.  Please see help file.\n");
  if (cl_present( "r" ))
    Abort ("Option r has been expanded to replace|rep.  Please see help file.\n");

  cl_get( "chunk|chu|c", "%option %s[%]", "images", in_chunk);
  cl_get( "chunk_out|cho", "%option %s[%]", in_chunk, out_chunk);

  /* Secret back door left in to make Joel happy; not included in help file */
  cl_get( "chunk_file|chf", "%option %s[%]", "", out_chunk_file);

  replace = cl_present("replace|rep");

  if(!cl_get("", "%s", in_file)) {
    fprintf(stderr, "%s: Input file name not given.\n", argv[0]);
    exit(-1);
  }
  if(!cl_get("", "%s", out_file)) {
    fprintf(stderr, "%s: Output file name not given.\n", argv[0]);
    exit(-1);
  }

  if (cl_cleanup_check()) {
    int i;
    fprintf(stderr,"%s: invalid argument in command line:\n    ", argv[0]);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(-1);
  }
  /*** End command-line parsing ***/

  /* check that flags make sense */
  if (argc <= 1)
    {
      Help( "usage" );
      exit(1);
    }
  if (strcmp(in_file, out_file) == 0 && strcmp(in_chunk, out_chunk) == 0)
    {
      fprintf(stderr, "mri_copy_chunk requires that you specify either a\n");
      fprintf(stderr, "  different output dataset or a different chunk name\n");
      exit(1);
    }

  /* open the input dataset */
  mri_set_error_handling(MRI_IGNORE_ERRORS);
  ids = mri_open_dataset(in_file, MRI_READ);
  if (ids == NULL)
    {
      fprintf(stderr, "%s: cannot open input dataset: %s\n", argv[0],
	      in_file);
      exit(1);
    }
  /* check if the chunk exists in the input dataset */
  if (!mri_has(ids, in_chunk) ||
      strcmp(mri_get_string(ids, in_chunk), "[chunk]") != 0)
    {
      fprintf(stderr, "%s: Input dataset %s has no chunk named %s\n", argv[0],
	      in_file, in_chunk);
      mri_close_dataset(ids);
      exit(1);
    }

  /* open the output dataset */
  if (strcmp(in_file, out_file) == 0)
    {
      mri_close_dataset(ids);
      ods = mri_open_dataset(out_file, MRI_MODIFY);
      if (ods == NULL)
	{
	  fprintf(stderr, "%s: cannot modify dataset: %s\n", argv[0],
		  out_file);
	  exit(1);
	}
      ids = ods;
    }
  else
    {
      /* check if output dataset already exists */
      ods = mri_open_dataset(out_file, MRI_MODIFY);
      if (ods == NULL)
	{
	  fprintf(stderr, "mri_copy_chunk: cannot open output dataset: %s\n",
		  out_file);
	  exit(1);
	}
    }
  mri_set_error_handling(MRI_ABORT_ON_ERROR);

  hist_add_cl(ods,argc,argv);

  /* clear out any conflicting chunk in the output dataset */
  olen = strlen(out_chunk);
  mri_iterate_over_keys(ods);
  while ((key = mri_next_key(ods)) != NULL)
    if (strncmp(key, out_chunk, olen) == 0 &&
	(key[olen] == '\0' || key[olen] == '.'))
      {
	if (!replace)
	  {
	    fprintf(stderr, "You must use the -replace flag if you want\n");
	    fprintf(stderr, "  to overwrite an existing chunk.\n");
	    mri_close_dataset(ids);
	    if (ods != ids)
	      mri_close_dataset(ods);
	    exit(1);
	  }
	mri_remove(ods, key);
      }
  
  /* create the chunk in the output dataset */
  /* first copy the associated keys */
  strcpy(new_key, out_chunk);
  ilen = strlen(in_chunk);
  mri_iterate_over_keys(ids);
  while ((key = mri_next_key(ids)) != NULL)
    if (strncmp(key, in_chunk, ilen) == 0 &&
	key[ilen] == '.')
      {
	strcpy(&new_key[olen], &key[ilen]);
	mri_set_string(ods, new_key, mri_get_string(ids, key));
      }
  /* fix up the file name */
  strcpy(&new_key[olen], ".file");
  if (out_chunk_file[0] != '\0') 
    {
      /* The output file name ".mri" means put it in the MRI file,
       * which is implemented by dropping the .file key.
       */
      if (!strcmp(out_chunk_file,".mri")) 
	mri_remove(ods, new_key);
      else
	mri_set_string(ods, new_key, out_chunk_file);
    }
  else
    {
      /* we prohibit absolute filenames from being copied since
         datasets shouldn't share files */
      if (mri_has(ods, new_key) &&
	  mri_get_string(ods, new_key)[0] != '.')
	mri_remove(ods, new_key);
    }
  /* remove any .order field that was copied */
  strcpy(&new_key[olen], ".order");
  if (mri_has(ods, new_key))
    mri_remove(ods, new_key);
  /* really create the chunk */
  mri_create_chunk(ods, out_chunk);

  in_size = mri_get_int(ids, mri_cat(in_chunk, ".size"));
  out_size = mri_get_int(ods, mri_cat(out_chunk, ".size"));
  if (in_size != out_size)
    {
      fprintf(stderr, "Inconsistent size field in input dataset!\n");
      mri_close_dataset(ids);
      if (ods != ids)
	mri_close_dataset(ods);
      exit(1);
    }
  size = in_size;
  if (size > BLOCK_SIZE)
    size = BLOCK_SIZE;

  for (offset = 0; offset < in_size; offset += BLOCK_SIZE)
    {
      size = in_size - offset;
      if (size > BLOCK_SIZE)
	size = BLOCK_SIZE;
      buf = mri_get_chunk(ids, in_chunk, size, offset, MRI_RAW);
      mri_set_chunk(ods, out_chunk, size, offset, MRI_RAW, buf);
    }

  mri_close_dataset(ids);
  if (ods != ids)
    mri_close_dataset(ods);
  exit(0);
}
