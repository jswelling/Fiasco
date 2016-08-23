/************************************************************
 *                                                          *
 *  mri_type_convert.c                                     *
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
 *  Original programming by Joel Welling, 7/98              *
 ************************************************************/
/*************************************************************

  DESCRIPTION OF MRI_TYPE_CONVERT

  mri_type_convert takes a pgh MRI dataset of any data type and
  converts a given chunk to another datatype.  

  The various type switches are mutually exclusive.  If a conversion 
  occurs which produces a value which is out of range for the new
  datatype, a warning message is given and conversion continues.

**************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "misc.h"
#include "stdcrg.h"

#define CHUNKSIZE 16384
#define KEYBUF_SIZE 512

static char rcsid[] = "$Id: mri_type_convert.c,v 1.7 2007/07/06 18:45:53 welling Exp $";

char* progname;
static float out_chunk[CHUNKSIZE];

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

int main( int argc, char* argv[] ) 
{

  MRI_Dataset *Input = NULL, *Output = NULL;
  int vec_length;
  long total_voxels;
  long voxels_moved;
  int voxels_this_chunk;
  int i;
  char infile[512], outfile[512];
  char keybuf[KEYBUF_SIZE];
  char chunk[KEYBUF_SIZE];
  char* dimstr= NULL;
  char* crunner;
  void* data_in= NULL;
  int type_in;
  int type= MRI_FLOAT;
  int type_flag_count= 0;
  int verbose_flg= 0;

  progname= argv[0];

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Deprecate old options */

  /* Get filenames */
  cl_get("chunk|chu|c", "%option %s[%]","images",chunk);
  verbose_flg= cl_present("v|ver|verbose");
  if (cl_present("char")) {
    type= MRI_UNSIGNED_CHAR;
    type_flag_count++;
  }
  if (cl_present("short|shr")) {
    type= MRI_SHORT;
    type_flag_count++;
  }
  if (cl_present("integer|int")) {
    type= MRI_INT;
    type_flag_count++;
  }
  if (cl_present("float|flt")) {
    type= MRI_FLOAT;
    type_flag_count++;
  }
  if (cl_present("double|dbl")) {
    type= MRI_DOUBLE;
    type_flag_count++;
  }
  if (type_flag_count>1) {
    fprintf(stderr,"%s: mutually exclusive type flags given.",argv[0]);
  }
  if (!cl_get("", "%s", infile)) {
    fprintf(stderr,"%s: Input file name not given.",argv[0]);
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

  /* Print version number */
  if (verbose_flg) Message( "# %s\n", rcsid );

  /* Open input dataset */
  if( !strcmp( infile, outfile ) )
    Abort( "Input and output files must be distinct." );
  Input = mri_open_dataset( infile, MRI_READ );

  /* Check that program will function on data-set */
  if( !mri_has( Input, chunk ) )
    Abort( "%s: requested chunk %s not present.\n", argv[0],chunk );
  safe_copy(keybuf,chunk);
  safe_concat(keybuf,".dimensions");
  if ( !mri_has( Input, keybuf ) ) 
    Abort( "%s needs %s info for input file.\n", argv[0], keybuf );
  dimstr= mri_get_string( Input, keybuf );

  /* Figure out how many instances of this vector */
  total_voxels= 1;
  for (crunner= dimstr; *crunner; crunner++) {
    char dimbuf[2];
    safe_copy(keybuf,chunk);
    safe_concat(keybuf,".extent.");
    dimbuf[0]= *crunner;
    dimbuf[1]= '\0';
    safe_concat(keybuf,dimbuf);
    if (!mri_has(Input, keybuf)) 
      Abort("%s: input has dimension %c but no %s tag.",
	    argv[0],*crunner,keybuf);
    else total_voxels *= mri_get_int( Input, keybuf );
  }
  
  /* Set up output parameters */
  Output = mri_copy_dataset( outfile, Input );
  hist_add_cl( Output, argc, argv );
  safe_copy(keybuf,chunk);
  safe_concat(keybuf,".datatype");
  switch (type) {
  case MRI_UNSIGNED_CHAR: mri_set_string( Output, keybuf, "uint8" ); break;
  case MRI_SHORT: mri_set_string( Output, keybuf, "int16" ); break;
  case MRI_INT: mri_set_string( Output, keybuf, "int32" ); break;
  case MRI_FLOAT: mri_set_string( Output, keybuf, "float32" ); break;
  case MRI_DOUBLE: mri_set_string( Output, keybuf, "float64" ); break;
  }

  /* Copy everything over, counting on libmri to do the actual 
   * translations.
   */
  type_in= get_chunk_type( Input, chunk );
  voxels_moved= 0;
  while (voxels_moved < total_voxels) {
    voxels_this_chunk= ((total_voxels - voxels_moved) > CHUNKSIZE) ?
      CHUNKSIZE : total_voxels - voxels_moved;
    if (!(data_in= mri_get_chunk(Input, chunk, voxels_this_chunk, 
				  voxels_moved, type_in)))
      Abort("%s: bad mri_get_chunk; length %ld, offset %ld\n",
	    argv[0], voxels_this_chunk, voxels_moved);
    mri_set_chunk(Output, chunk, voxels_this_chunk, voxels_moved, 
		  type_in, data_in);
    voxels_moved += voxels_this_chunk;
  }

  /* Write and close data-sets */
  mri_close_dataset( Input );
  mri_close_dataset( Output );
  
  if (verbose_flg) Message( "#      Type conversion complete.\n" );

  return 0;
}

