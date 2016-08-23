/************************************************************
 *                                                          *
 *  mri_from_ascii.c                                            *
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

  DESCRIPTION OF MRI_FROM_ASCII

  mri_from_ascii reads a stream of ascii data representing floats,
  and outputs a Pittsburgh MRI file containing the data.

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
#define BLOCK 2048

static char rcsid[] = "$Id: mri_from_ascii.c,v 1.10 2003/08/07 20:20:25 bakalj Exp $";

static MRI_Dataset *Output = NULL;
static char chunk[512];
char* progname;
static int verbose_flg= 0;

static void safe_copy(char* str1, char* str2) {
  strncpy(str1, str2, KEYBUF_SIZE);
  str1[KEYBUF_SIZE-1]= '\0';
}

static void safe_concat(char* str1, char* str2) {
  strncat(str1, str2, (KEYBUF_SIZE-strlen(str1))-1);
}

static void parse_and_save(float* buf, char* inbuf, int* ext_tbl,
			   int* stride_tbl, int ext_length, 
			   int index_length, int vals_length)
{
  float vals[64];
  int ind[64];
  int i;
  char* tok;
  int blank_flag= 0;
  int offset;

  if (*inbuf == '#') /* comment */
    return;

  if (index_length>64 || vals_length>64)
    Abort("%s: too many values per line!\n",progname);

  for (i=0; i<index_length; i++) {
    tok= strtok( (i ? NULL : inbuf), " \t\n" );
    if (!tok) {
      if (i==0) {
	blank_flag= 1;
	break;
      }
      else Abort("%s: short line found!\n",progname);
    }
    if (sscanf(tok,"%d",ind+i) != 1)
      Abort("%s: encountered invalid index <%s> in input!\n",
	    progname,tok);
    if (ind[i]<0 || ind[i]>=ext_tbl[i+(ext_length-index_length)])
      Abort("%s: encountered out-of-range index %d in column %d!\n",
	    progname,ind[i],i+1);
  }
  if (!blank_flag) {
    for (i=0; i<vals_length; i++) {
      tok= strtok( NULL, " \t\n" );
      if (!tok) Abort("%s: short line found!\n",progname);
      if (sscanf(tok,"%g",vals+i) != 1)
	Abort("%s: encountered invalid value <%s> in input!\n",
	      progname,tok);
    }
    if (strtok(NULL, " \t\n"))
      Abort("%s: long line found! <%s>\n",progname,tok);
  }

  /* Store the value */
  offset= 0;
  for (i=0; i<index_length; i++) 
    offset += ind[i]*stride_tbl[i+(ext_length-index_length)];
  for (i=0; i<vals_length; i++) buf[offset+i]= vals[i];
}

static void parse_indexed( int* ext_tbl, int* stride_tbl, int ext_length,
			   int index_length ) 
{
  int total_size;
  float* buf;
  int n_components;
  char inbuf[512];
  char inbuf2[512];
  int vals_length;
  int i;

  n_components= 0;
  while (!n_components) { /* keep trying 'till we find some data */
    fgets(inbuf,512,stdin);
    if (inbuf[0]=='#') continue; /* skip this comment */
    strcpy(inbuf2,inbuf); /* keep the original */
    
    /* count tokens */
    if (strtok(inbuf2," \t\n")) n_components++;
    while (strtok(NULL," \t\n")) n_components++;
  }

  /* check dimensionality */
  vals_length= stride_tbl[ext_length-index_length];
  if (n_components != index_length + vals_length)
    Abort("%s: wrong number of values per line (expected %d)\n",progname,
	  index_length + vals_length);

  /* In this case we read the whole lot into memory. */
  total_size= ext_tbl[ext_length-1]*stride_tbl[ext_length-1];
  if (!(buf=(float*)malloc(total_size*sizeof(float))))
    Abort("%s: unable to allocate %d bytes!\n",total_size*sizeof(float));
  for (i=0; i<total_size; i++) buf[i]= 0.0;

  parse_and_save(buf,inbuf, ext_tbl, stride_tbl, ext_length, 
		 index_length, vals_length);
  while (!feof(stdin) && !ferror(stdin)) {
    fgets(inbuf,512,stdin);
    if (feof(stdin) || ferror(stdin)) break;
    parse_and_save(buf, inbuf, ext_tbl, stride_tbl, ext_length, 
		   index_length, vals_length);
  }
  if (ferror(stdin)) {
    perror(progname);
    Abort("%s: error reading data!\n",progname);
  }
    
  /* Save the data and free the buffer */
  mri_set_chunk(Output,chunk,total_size,0,MRI_FLOAT,buf);
  free(buf);
}

static void parse_stream( int expected )
{
  int offset= 0;
  int i;
  float buf[BLOCK];

  while (offset<expected) {
    int count;
    int this_block= (BLOCK<expected-offset) ? BLOCK : expected-offset;
    for (i=0; i<this_block; i++) 
      if (scanf("%g",buf+i) != 1) break;
      
    if (i != this_block || feof(stdin) || ferror(stdin)) {
      if (ferror(stdin)) perror(progname);
      Abort("%s: out of data reading input!\n",progname);
    }

    mri_set_chunk(Output,chunk,this_block,offset,MRI_FLOAT,buf);
    offset += this_block;
  }
}

static void parse_extents(int* ext_tbl, int ext_length, char* extents)
{
  char* delimiters= ",:;*x- \t";
  int which= 0;
  char* str;
  int ns;
  
  while (*extents) {
    /* Check for initial delimiter */
    if (strchr(delimiters,*extents)) {
      ext_tbl[which++]= 1;
      extents++;
    }
    else {
      str= strtok(extents,delimiters);
      if (sscanf(str, "%d", &ns) != 1)
	Abort("%s: non-integral extent found (%s)", progname,str);
      ext_tbl[which++]= ns;
      extents += (strlen(str)+1);
    }
    if (which >= ext_length) break;
  }

  /* Handle trailing delimiter */
  if (*(extents-1) && strchr(delimiters,*(extents-1))) {
    ext_tbl[which++]= 1;
    extents++;
  }

  if (which != ext_length || *(extents-1) != '\0')
    Abort("%s: mismatch between dimension length and extents!\n",progname);
}

int main( int argc, char* argv[] ) 
{
  char outfile[512];
  char dimstr[512];
  char extents[512];
  char index_string[512];
  char keybuf[KEYBUF_SIZE];
  int total_size;
  char tbuf[2];
  int i;
  int ext_length;
  int* ext_tbl;
  int* stride_tbl;
  int index_mode= 0;

  progname= argv[0];

  /* Print version number */
  if (verbose_flg) Message( "# %s\n", rcsid );

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Deprecate old options */

  if (cl_present( "d" ))
     Abort ("Option d has been replaced by order|ord|o.  Please see help file.\n");

  if (cl_present( "e" ))
     Abort ("Option e has been replaced by length|len|l.  Please see help file.\n");

  if (cl_present( "i" ))
     Abort ("Option i has been expanded to index|ind.  Please see help file.\n");

  if (cl_present( "o" ))
     Abort ("Option o used to designate the outfile and has been replaced by outfile format.  If you mean order, please use order|ord.  Please see help file.\n");
 
  /* Get parameters */
  verbose_flg= cl_present("verbose|ver|v");
  cl_get("order|ord","%option %s[%]","t",dimstr);
  if (!cl_get("length|len|l","%option %s",extents)) {
    fprintf(stderr,"%s: required extent argument not given.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  cl_get("chunk|chu|c","%option %s[%]","images",&chunk);
  index_mode= cl_get("index|ind","%option %s",index_string);

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

  /* Reality-check arguments */
  if (index_mode) {
    if (strcmp(index_string, 
	       dimstr + (strlen(dimstr)-strlen(index_string)))) {
      fprintf(stderr,"%s: index string must be last segment of dim string!\n",
	      argv[0]);
      Help( "usage" );
      exit(-1);
    }
  }

  /* Open output dataset */
  Output = mri_open_dataset( outfile, MRI_WRITE );
  hist_add_cl(Output,argc,argv);
  mri_create_chunk(Output,chunk);
  safe_copy(keybuf,chunk);
  safe_concat(keybuf,".datatype");
  mri_set_string( Output, keybuf, "float32" );
  safe_copy(keybuf,chunk);
  safe_concat(keybuf,".dimensions");
  mri_set_string( Output, keybuf, dimstr );
  safe_copy(keybuf,chunk);
  safe_concat(keybuf,".file");
  mri_set_string( Output, keybuf, ".dat" );

  /* Parse and set extents */
  ext_length= strlen(dimstr);
  if (ext_length<=0)
    Abort("%s: internal error: invalid dimstring!\n");
  if (!(ext_tbl= (int*)malloc(ext_length*sizeof(int))))
    Abort("%s: Unable to allocate %d bytes!\n",ext_length*sizeof(int));
  if (!(stride_tbl= (int*)malloc(ext_length*sizeof(int))))
    Abort("%s: Unable to allocate %d bytes!\n",ext_length*sizeof(int));
  parse_extents(ext_tbl,ext_length,extents);
  total_size = 1;
  for (i = 0; i < strlen(dimstr); ++i)
    {
      safe_copy(keybuf,chunk);
      safe_concat(keybuf,".extent.");
      tbuf[0]= dimstr[i];
      tbuf[1]= '\0';
      safe_concat(keybuf, tbuf);
      mri_set_int( Output, keybuf, ext_tbl[i] );
      stride_tbl[i]= total_size;
      total_size *= ext_tbl[i];
    }
  
  if (index_mode) parse_indexed(ext_tbl,stride_tbl,ext_length,
				strlen(index_string));
  else parse_stream(total_size);

  /* Write and close data-sets */
  mri_close_dataset( Output );
  
  return 0;
}

