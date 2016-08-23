/*							    *
 *	mri_remap.c - remaps the dimensions of a chunk without  *
 *                moving any actual data around             *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *                                                          *
 *  Copyright (c) 1997 Pittsburgh Supercomputing Center     *
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
 *							    *
 *	HISTORY:                                            *
 *		2/97 - Written by Greg Hood                 */

/*************************************************************

  DESCRIPTION OF MRI_REMAP

  mri_remap is used to remap the specified chunk in the named dataset
  to have the new dimensions; it does not actually move chunk data
  around; it is good for adding new dimensions of length 1 to a chunk,
  or for subdividing one dimension into two others.  For example, we
  might have 36 images in the t dimension, and want to redimension the
  chunk to have a c ("condition") dimension of 4 and a t dimension of
  9.  The total number of elements in the input and output chunk must
  be the same.

**************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mri.h"
#include "stdcrg.h"
#include "fmri.h"

static char rcsid[] = "$Id: mri_remap.c,v 1.7 2005/07/01 21:06:06 welling Exp $";

static int verbose_flg= 0;
static int preserve_flg= 0;

static int legit_dim_string( const char* str )
{
  char* here= (char*)str;
  char* runner;

  while (*here) {
    runner= here+1;
    while (*runner) {
      if (*runner==*here) {
	/* characters at *here and *runner match */
	return 0;
      }
      runner++;
    }
    here++;
  }

  return 1;
}

int main( int argc, char* argv[] ) 
{

  MRI_Dataset *ds = NULL;
  char file[512];
  char chunk_name[256];
  char field[256];
  int old_size[256];
  char *old_dims;
  int old_total_size;
  int new_size[256];
  char new_dims[256];
  int new_total_size;
  char new_extents[256];
  int i;
  int j;
  char *str;
  char* here;
  char *p;
  int ns;
  int len;
  char *key;
  int trivial;

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Deprecate old options */

  if (cl_present( "extents|e" ))
     Abort ("Option extents|e has been replaced by length|len|l.  Please see help file.\n");

  if (cl_present( "file|f" ))
     Abort ("Option file|f has been replaced by infile format.  Please see help file.\n");


  /* Get filenames */
  cl_get( "chunk|chu|c", "%option %s[%]", "images", chunk_name);
  cl_get( "order|ord|o", "%option %s[%]", "vxyzt", new_dims);
  cl_get( "length|len|l", "%option %s[%]", "", new_extents);
  verbose_flg= cl_present( "verbose|ver|v" );
  preserve_flg= cl_present( "preserve" );
  if(!cl_get("", "%s", file)) {
    fprintf(stderr, "%s: Input file name not given.\n", argv[0]);
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

  /* Is chunk name a reasonable length? */
  if (strlen(chunk_name)>128)
    Abort("%s: chunk name <%s> is too long!\n",argv[0],chunk_name);

  /* Open dataset */
  ds = mri_open_dataset( file, MRI_MODIFY );
  hist_add_cl(ds,argc,argv);

  /* Check that program will function on data-set */
  if( !mri_has( ds, chunk_name ) )
    Abort( "%s: chunk %s is missing.\n", argv[0], chunk_name );
  sprintf(field, "%s.dimensions", chunk_name);
  if(!mri_has(ds, field) )
    Abort( "%s does not have the %s key.\n", file, field );

  /* Are new_dims legitimate? */
  if (!legit_dim_string(new_dims))
    Abort( "%s: dimension string <%s> is not valid!\n",argv[0],new_dims);



  old_dims = strdup(mri_get_string(ds, field));
  old_total_size = 1;
  for (i = 0; i < 256; ++i)
    {
      old_size[i] = 1;
      new_size[i] = 1;
    }
  for (i = 0; i < strlen(old_dims); ++i)
    {
      sprintf(field, "%s.extent.%c", chunk_name, old_dims[i]);
      if (mri_has(ds, field))
	old_size[(int) old_dims[i]] = mri_get_int(ds, field);
      else
	old_size[(int) old_dims[i]] = 1;
      if (old_size[(int) old_dims[i]] < 0)
	Abort("chunk has negative extent: %s = %d\n",
	      field, old_size[(int) old_dims[i]]);
      old_total_size *= old_size[(int) old_dims[i]];
    }

  new_total_size = 1;
  here= new_extents;
  for (i = 0; i < strlen(new_dims); ++i)
    {
      char* runner= here;
      while (*runner) {
	if (*runner==',' || *runner==':' || *runner==';' || *runner=='*' || 
	    *runner=='x' || *runner=='-' || *runner==' ' || *runner=='\t' ) 
	  break;
	runner++;
      }
      if (*runner) {
	*runner= '\0';
	str= here;
	here= runner+1;
      }
      else {
	if (runner != here) { str= here; here= runner; }
	else str= NULL;
      }
      
      if (str == NULL || str[0]=='\0')
	ns = old_size[(int) new_dims[i]];
      else
	if (sscanf(str, "%d", &ns) != 1)
	  Abort("Non-integral extent found (%s)\n", str);
      new_size[(int) new_dims[i]] = ns;
      new_total_size *= ns;
    }

  if (new_total_size != old_total_size)
    Abort("New chunk size (%d) does not match old chunk size (%d)\n",
	  new_total_size, old_total_size);

  /* Delete tags associated with dimensions which have vanished
   * or changed extent.
   */
  if (!preserve_flg) {
    i= 0;
    len = strlen(chunk_name);
    while (old_dims[i]) {
      if (!strchr(new_dims,old_dims[i]) 
	  || new_size[old_dims[i]] != old_size[old_dims[i]]) {
	mri_iterate_over_keys(ds);
	while ((key = mri_next_key(ds)) != NULL)
	  if (strncmp(key, chunk_name, len) == 0 &&
	      key[len] == '.' &&
	      (p = strchr(&key[len+1], '.')) != NULL &&
	      (p[1]==old_dims[i])&&
	      (p[2] == '.' || p[2] == '\0')) {
	    if (verbose_flg)
	      Message("Removing tag <%s>\n",key);
	    mri_remove(ds, key);
	  }
      }
      i++;
    }
  }

  /* add all the new keys */
  sprintf(field, "%s.dimensions", chunk_name);
  mri_set_string(ds, field, new_dims);
  for (i = 0; i < strlen(new_dims); ++i)
    {
      sprintf(field, "%s.extent.%c", chunk_name, new_dims[i]);
      mri_set_int(ds, field, new_size[(int) new_dims[i]]);
    }

  /* remove old extent tags */
  for (i=0; i<strlen(old_dims); i++) {
    if (!strchr(new_dims,old_dims[i])) {
      sprintf(field,"%s.extent.%c",chunk_name,old_dims[i]);
      if (mri_has(ds,field)) mri_remove(ds,field);
    }
  }

  /* Write and close dataset */
  mri_close_dataset(ds);
  
  if (verbose_flg) Message( "#      Remap complete.\n" );
  return 0;
}

