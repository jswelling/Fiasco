/************************************************************
 *                                                          *
 *  mri_sort.c                                               *
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

  DESCRIPTION OF MRI_SORT

  mri_sort sorts a selected chunk of a pgh MRI file.  The first
  dimension must be v; the data is sorted along its second
  dimension.

**************************************************************/

#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "misc.h"
#include "stdcrg.h"

#define KEYBUF_SIZE 512

static char rcsid[] = "$Id: mri_sort.c,v 1.8 2007/07/27 17:58:00 welling Exp $";

static MRI_Dataset *Input = NULL, *Output = NULL;
static char* progname;
static int verbose_flg= 0;
static int keyfield= 0;

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

static int safe_get_extent(MRI_Dataset* ds, char* chunk, const char* dim)
{
  char key_buf[KEYBUF_SIZE];
  char dim_buf[4];
  dim_buf[0]= *dim;
  dim_buf[1]= '\0';
  safe_copy(key_buf,chunk);
  safe_concat(key_buf,".extent.");
  safe_concat(key_buf,dim_buf);
  if (mri_has(ds,key_buf)) return mri_get_int(ds,key_buf);
  else Abort("%s: input missing tag %s!\n",progname,key_buf);
  return 0; /* not reached */
}

static void calc_sizes(char* this_chunk, char* dimstr, const char selected_dim,
		       int* fast_blocksize_out, int* slow_blocksize_out ) {
  /* This routine will fail if selected_dim is not in dimstr! */
  int fast_blocksize;
  int slow_blocksize;
  char* this_dim;

  fast_blocksize= 1;

  this_dim= dimstr;
  while (*this_dim != selected_dim) 
    fast_blocksize *= safe_get_extent(Input,this_chunk,this_dim++);

  this_dim++; /* step over selected dim */

  slow_blocksize= 1;
  while (*this_dim)
    slow_blocksize *= safe_get_extent(Input, this_chunk, this_dim++);

  *fast_blocksize_out= fast_blocksize;
  *slow_blocksize_out= slow_blocksize;
}

static int comp_char_a(const void* v1, const void* v2)
{
  char* cv1= (char*)v1;
  char* cv2= (char*)v2;
  if (cv1[keyfield]<cv2[keyfield]) return -1;
  else if (cv1[keyfield]>cv2[keyfield]) return +1;
  else return 0;
}

static int comp_char_d(const void* v1, const void* v2)
{
  char* cv1= (char*)v1;
  char* cv2= (char*)v2;
  if (cv1[keyfield]>cv2[keyfield]) return -1;
  else if (cv1[keyfield]<cv2[keyfield]) return +1;
  else return 0;
}

static int comp_short_a(const void* v1, const void* v2)
{
  short* cv1= (short*)v1;
  short* cv2= (short*)v2;
  if (cv1[keyfield]<cv2[keyfield]) return -1;
  else if (cv1[keyfield]>cv2[keyfield]) return +1;
  else return 0;
}

static int comp_short_d(const void* v1, const void* v2)
{
  short* cv1= (short*)v1;
  short* cv2= (short*)v2;
  if (cv1[keyfield]>cv2[keyfield]) return -1;
  else if (cv1[keyfield]<cv2[keyfield]) return +1;
  else return 0;
}

static int comp_int_a(const void* v1, const void* v2)
{
  int* cv1= (int*)v1;
  int* cv2= (int*)v2;
  if (cv1[keyfield]<cv2[keyfield]) return -1;
  else if (cv1[keyfield]>cv2[keyfield]) return +1;
  else return 0;
}

static int comp_int_d(const void* v1, const void* v2)
{
  int* cv1= (int*)v1;
  int* cv2= (int*)v2;
  if (cv1[keyfield]>cv2[keyfield]) return -1;
  else if (cv1[keyfield]<cv2[keyfield]) return +1;
  else return 0;
}

static int comp_longlong_a(const void* v1, const void* v2)
{
  long long* cv1= (long long*)v1;
  long long* cv2= (long long*)v2;
  if (cv1[keyfield]<cv2[keyfield]) return -1;
  else if (cv1[keyfield]>cv2[keyfield]) return +1;
  else return 0;
}

static int comp_longlong_d(const void* v1, const void* v2)
{
  long long* cv1= (long long*)v1;
  long long* cv2= (long long*)v2;
  if (cv1[keyfield]>cv2[keyfield]) return -1;
  else if (cv1[keyfield]<cv2[keyfield]) return +1;
  else return 0;
}

static int comp_float_a(const void* v1, const void* v2)
{
  float* cv1= (float*)v1;
  float* cv2= (float*)v2;
  if (isnan(cv1[keyfield])) return +1; /* Sort NaNs out the top */
  if (isnan(cv2[keyfield])) return -1;
  if (cv1[keyfield]<cv2[keyfield]) return -1;
  else if (cv1[keyfield]>cv2[keyfield]) return +1;
  else return 0;
}

static int comp_float_d(const void* v1, const void* v2)
{
  float* cv1= (float*)v1;
  float* cv2= (float*)v2;
  if (isnan(cv1[keyfield])) return -1; /* Sort NaNs out the bottom */
  if (isnan(cv2[keyfield])) return +1;
  if (cv1[keyfield]>cv2[keyfield]) return -1;
  else if (cv1[keyfield]<cv2[keyfield]) return +1;
  else return 0;
}

static int comp_double_a(const void* v1, const void* v2)
{
  double* cv1= (double*)v1;
  double* cv2= (double*)v2;
  if (isnan(cv1[keyfield])) return +1; /* Sort NaNs out the top */
  if (isnan(cv2[keyfield])) return -1;
  if (cv1[keyfield]<cv2[keyfield]) return -1;
  else if (cv1[keyfield]>cv2[keyfield]) return +1;
  else return 0;
}

static int comp_double_d(const void* v1, const void* v2)
{
  double* cv1= (double*)v1;
  double* cv2= (double*)v2;
  if (isnan(cv1[keyfield])) return -1; /* Sort NaNs out the bottom */
  if (isnan(cv2[keyfield])) return +1;
  if (cv1[keyfield]>cv2[keyfield]) return -1;
  else if (cv1[keyfield]<cv2[keyfield]) return +1;
  else return 0;
}

static void transfer_data(char* this_chunk, 
			  int v_extent, int slow_blksize, 
			  int selected_extent, int sort_ascending)
{
  int in_offset= 0;
  int ifast;
  int islow;
  int collective_blksize;
  int type;
  int typesize;
  void* ibuf;

  type= get_chunk_type(Input,this_chunk);
  typesize= get_typesize(Input,this_chunk);

  collective_blksize= v_extent * selected_extent;
  
  for (islow=0; islow<slow_blksize; islow++) {

    ibuf= mri_get_chunk(Input, this_chunk, collective_blksize, 
			in_offset, type);

    switch (type) {

    case MRI_UNSIGNED_CHAR:
      {
	if (sort_ascending)
	  qsort(ibuf,selected_extent,v_extent*typesize,comp_char_a);
	else
	  qsort(ibuf,selected_extent,v_extent*typesize,comp_char_d);
      }
      break;
    case MRI_SHORT:
      {
	if (sort_ascending)
	  qsort(ibuf,selected_extent,v_extent*typesize,comp_short_a);
	else
	  qsort(ibuf,selected_extent,v_extent*typesize,comp_short_d);
      }
      break;
    case MRI_INT:
      {
	if (sort_ascending)
	  qsort(ibuf,selected_extent,v_extent*typesize,comp_int_a);
	else
	  qsort(ibuf,selected_extent,v_extent*typesize,comp_int_d);
      }
      break;
    case MRI_LONGLONG:
      {
	if (sort_ascending)
	  qsort(ibuf,selected_extent,v_extent*typesize,comp_longlong_a);
	else
	  qsort(ibuf,selected_extent,v_extent*typesize,comp_longlong_d);
      }
      break;
    case MRI_FLOAT:
      {
	if (sort_ascending)
	  qsort(ibuf,selected_extent,v_extent*typesize,comp_float_a);
	else
	  qsort(ibuf,selected_extent,v_extent*typesize,comp_float_d);
      }
      break;
    case MRI_DOUBLE:
      {
	if (sort_ascending)
	  qsort(ibuf,selected_extent,v_extent*typesize,comp_double_a);
	else
	  qsort(ibuf,selected_extent,v_extent*typesize,comp_double_d);
      }
      break;
    }

    mri_set_chunk( Output, this_chunk, collective_blksize, in_offset, type,
		   ibuf );

    if (verbose_flg) 
      fprintf(stderr,"block: type %d, %d at %d\n",
	      type, collective_blksize, in_offset);
    in_offset += collective_blksize;
  }

}

static int sortIndexRemap(int oldIndex, void* hook)
{
  return -1; /* to cause label to be deleted */
}

static void sort_chunk(char* this_chunk, const char sort_dim, 
		       const int sort_ascending) {
  char key_buf[KEYBUF_SIZE];
  char* dimstr;
  int v_extent;
  int slow_blksize;
  int selected_extent;

  safe_copy(key_buf, this_chunk);
  safe_concat(key_buf, ".dimensions");
  dimstr= mri_get_string(Input,key_buf);
  selected_extent= safe_get_extent(Input,this_chunk,&sort_dim);

  /* Wipe away any labels, since different instances of the dimension
   * will sort into different places.
   */
  mriu_updateLabels(Output, this_chunk, sort_dim, 
		    0, selected_extent,
		    sortIndexRemap, NULL);
 
  calc_sizes(this_chunk, dimstr, sort_dim, &v_extent, &slow_blksize);
  transfer_data(this_chunk, v_extent, slow_blksize, 
		selected_extent, sort_ascending);
}

int main( int argc, char* argv[] ) 
{

  char infile[512], outfile[512];
  char this_chunk[KEYBUF_SIZE];
  char key_buf[KEYBUF_SIZE];
  char sort_dim;
  int sort_ascending;
  char* dimstr;

  progname= argv[0];

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Deprecate old options */

  if (cl_present( "a" ))
     Abort ("Option a(ascending) has been expanded to ascending|asc.  Please see help file.\n");

  if (cl_present( "k" ))
     Abort ("Option k(keyfield) has been expanded to keyfield|key.  Please see help file.\n");


  verbose_flg= cl_present("verbose|ver");
  cl_get("keyfield|key", "%option %d[0]",&keyfield);
  cl_get("chunk|chu|c", "%option %s[%]","images",this_chunk);
  sort_ascending= cl_present("ascending|asc");
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

  if (!strcmp(this_chunk,"missing")) 
    Abort("%s: it is not legal to sort the 'missing' chunk!\n",argv[0]);

  /* Open input dataset, and check consistency */
  if( !strcmp( infile, outfile ) )
    Abort( "%s: Input and output files must be distinct.", argv[0] );
  Input = mri_open_dataset( infile, MRI_READ );
  safe_copy(key_buf,this_chunk);
  safe_concat(key_buf,".dimensions");
  if (!mri_has(Input,key_buf))
    Abort("%s: input is missing %s tag!\n",argv[0],key_buf);
  dimstr= mri_get_string(Input,key_buf);
  if (dimstr[0] != 'v')
    Abort("%s: input chunk's first dimension must be 'v'!\n",argv[0]);
  if (safe_get_extent(Input,this_chunk,"v")<=keyfield)
    Abort("%s: key field index is too large for v extent!\n",argv[0]);
  if (strlen(dimstr)<2) 
    Abort("%s: input chunk has no dimension to sort!\n",argv[0]);
  sort_dim= dimstr[1];

  if (verbose_flg)
    fprintf(stderr,"Sorting input chunk <%s>, dims %s, on %c in %s order.\n",
	    this_chunk,dimstr,sort_dim,
	    (sort_ascending?"ascending":"descending"));

  /* Open output dataset */
  Output = mri_copy_dataset( outfile, Input );
  hist_add_cl( Output, argc, argv );

  /* Do the sort */
  sort_chunk(this_chunk, sort_dim, sort_ascending);

  /* Write and close data-sets */
  mri_close_dataset( Input );
  mri_close_dataset( Output );
  
  if (verbose_flg) Message( "#      Sorting complete.\n" );

  return 0;
}

