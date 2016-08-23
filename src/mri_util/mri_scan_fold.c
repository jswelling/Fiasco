/************************************************************
 *                                                          *
 *  mri_scan_fold.c                                            *
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

  DESCRIPTION OF MRI_SCAN_FOLD

  mri_scan_fold takes a Pgh MRI dataset as input, and folds the time
  dimension of the input dataset into a combination of slices (z
  values) and times.
**************************************************************/
/* Notes-
 * 
 */

#include <stdlib.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "misc.h"
#include "stdcrg.h"

#define DEFAULT_PATTERN "even/odd"
#define DEFAULT_DIRECTION 1
#define KEYBUF_SIZE 512

static char rcsid[] = "$Id: mri_scan_fold.c,v 1.12 2007/07/30 17:00:40 welling Exp $";

static MRI_Dataset *Input = NULL, *Output = NULL;
static char* progname= NULL;
static char cl_reorder_pattern[512];
static int cl_reorder_pattern_set= 0;
static long cl_z_extent= 0; /* zero means unset */
static int verbose_flag= 0;
static int debug_flag= 0;
static int force_flag= 0;
static int data_changed= 0;
/* Are we folding or un-folding?
 * direction==1 means acqisition order -> space order;
 * direction==0 is the reverse.
 */
static int cl_direction= 0;
static int cl_direction_set= 0;

static long get_chunk_type(MRI_Dataset* ds, const char* chunk)
{
  char key_buf[KEYBUF_SIZE];
  key_buf[KEYBUF_SIZE-1]= '\0';
  snprintf(key_buf,KEYBUF_SIZE-1,"%s.datatype",chunk);
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

static int safe_get_extent(MRI_Dataset* ds, const char* chunk, const char* dim)
{
  char key_buf[KEYBUF_SIZE];
  char dim_buf[4];
  dim_buf[0]= *dim;
  dim_buf[1]= '\0';
  key_buf[KEYBUF_SIZE-1]= '\0';
  snprintf(key_buf,KEYBUF_SIZE-1,"%s.extent.%c",chunk,*dim);
  if (mri_has(Input,key_buf)) return mri_get_int(Input,key_buf);
  else Abort("%s: input missing tag %s!\n",progname,key_buf);
  return 0; /* not reached */
}

static void calc_sizes(const char* this_chunk, const char* dimstr, 
		       long* fast_blocksize_out, 
		       long long* slow_blocksize_out ) {
  long fast_blocksize;
  long long slow_blocksize;
  const char* this_dim;

  fast_blocksize= 1;

  this_dim= dimstr;
  while (*this_dim != 'z' && *this_dim != 't') 
    fast_blocksize *= safe_get_extent(Input,this_chunk,this_dim++);

  /* step over selected dims, 'zt' or 't' */
  this_dim++; 
  if (*this_dim == 't') this_dim++;

  slow_blocksize= 1;
  while (*this_dim)
    slow_blocksize *= safe_get_extent(Input, this_chunk, this_dim++);

  *fast_blocksize_out= fast_blocksize;
  *slow_blocksize_out= slow_blocksize;
}

static const char* checkChunkForReorderPattern(const char* chunk)
{
  char key_buf[KEYBUF_SIZE];
  snprintf(key_buf,KEYBUF_SIZE,"%s.reorder_pattern",chunk);
  if (mri_has(Input,key_buf))
    return mri_get_string(Input,key_buf);
  else return NULL;
}

static int checkChunkForReorderDir(const char* chunk, int* val)
{
  char key_buf[KEYBUF_SIZE];
  snprintf(key_buf,KEYBUF_SIZE,"%s.reorder",chunk);
  if (mri_has(Input,key_buf)) {
    *val= mri_get_int(Input,key_buf);
    return 1;
  }
  else return 0;  
}

static void transfer_data(const char* this_chunk, long fast_blksize, 
			  long z_extent, long t_extent, long long slow_blksize,
			  int* indexMap)
{
  long long in_offset= 0;
  long long out_offset= 0;
  long type;
  void* in_chunk= NULL;
  long z;
  long t;
  long long islow;

  type= get_chunk_type(Input,this_chunk);

  for (islow=0; islow<slow_blksize; islow++) {
    for (t=0; t<t_extent; t++) {
      for (z=0; z<z_extent; z++) {
	in_chunk= mri_get_chunk(Input, this_chunk, fast_blksize,
				in_offset, type);
	in_offset += fast_blksize;
	mri_set_chunk(Output, this_chunk, fast_blksize,
		      out_offset+indexMap[z]*fast_blksize, type, in_chunk);
      }
      out_offset += z_extent*fast_blksize;
    }
  }
}

static int foldMapFunc(int oldIndex, void* hook)
{
  int* indexMap= (int*)hook;
  if (debug_flag)
    fprintf(stderr,"Mapping label %d -> %d\n",oldIndex, indexMap[oldIndex]);
  return indexMap[oldIndex];
}

static void restructure_tags(const char* this_chunk, long fast_blksize, 
			     long z_extent, long t_extent, 
			     long long slow_blksize, int direction,
			     const char* reorder_pattern, int* indexMap)
{
  char key_buf[KEYBUF_SIZE];
  char dim_buf[KEYBUF_SIZE]; /* wild overkill */
  const char* dimstr= NULL;
  const char* ztloc= NULL;
  const char* tloc= NULL;
  int i;

  key_buf[KEYBUF_SIZE-1]= '\0';
  dim_buf[KEYBUF_SIZE-1]= '\0';
  
  /* Update the chunk's embedded info about its reordering state */
  if (debug_flag) fprintf(stderr,"Setting <%s> dir %d pattern %s\n",
			  this_chunk,(direction?0:1),reorder_pattern);
  snprintf(key_buf,KEYBUF_SIZE,"%s.reorder",this_chunk);
  mri_set_int(Output,key_buf,(direction?0:1));
  snprintf(key_buf,KEYBUF_SIZE,"%s.reorder_pattern",this_chunk);
  mri_set_string(Output,key_buf,reorder_pattern);

  snprintf(key_buf,KEYBUF_SIZE-1,"%s.dimensions",this_chunk);
  dimstr= mri_get_string(Input,key_buf);
  /* Guaranteed to hit at least one of these from previous checks */
  ztloc= strstr(dimstr,"zt");
  tloc= strchr(dimstr,'t');
  if (ztloc) {
    /* Dimensions should be OK */
  }
  else if (tloc) {
    /* Need to paste in an appropriate z dimension, and resize t */
    for (i=0; (dimstr+i<tloc) && (i<KEYBUF_SIZE-2); i++) 
      dim_buf[i]= dimstr[i];
    if (i==KEYBUF_SIZE-2)
      Abort("%s: dimension string for chunk %s is too long!\n",
	    progname,this_chunk);
    dim_buf[i]= 'z';
    for (; (dimstr[i]!='\0') && (i<KEYBUF_SIZE-2); i++)
      dim_buf[i+1]= dimstr[i];
    dim_buf[i+1]= '\0';
    mri_set_string(Output,key_buf,dim_buf);
    snprintf(key_buf,KEYBUF_SIZE,"%s.extent.t",this_chunk);
    mri_set_int(Output,key_buf,t_extent);
    snprintf(key_buf,KEYBUF_SIZE,"%s.extent.z",this_chunk);
    mri_set_int(Output,key_buf,z_extent);
    snprintf(key_buf,KEYBUF_SIZE,"%s.description.z",this_chunk);
    mri_set_string(Output,key_buf,"gridded image-space");
  }
  else Abort("%s: internal error on chunk %s: can't find where to edit dimstr!\n",
	     progname,this_chunk);
  mriu_updateLabels(Output, this_chunk, 'z', 0, z_extent, 
		    foldMapFunc, indexMap);
}

static int fold_this(const char* this_chunk, long* fast_blksize, 
		     long* z_extent, long* t_extent, long long* slow_blksize)
{
  char key_buf[KEYBUF_SIZE];

  key_buf[KEYBUF_SIZE-1]= '\0';
  snprintf(key_buf,KEYBUF_SIZE-1,"%s.dimensions",this_chunk);
  if (mri_has(Input,key_buf)) {
    const char* dimstr= mri_get_string(Input,key_buf);
    const char* ztloc= strstr(dimstr,"zt");
    const char* tloc= strchr(dimstr,'t');
    
    if (debug_flag) fprintf(stderr,"chunk <%s> has dimstr %s\n",
			    this_chunk,dimstr);
    if (ztloc) {
      snprintf(key_buf,KEYBUF_SIZE-1,"%s.extent.z",this_chunk);
      *z_extent= mri_get_int(Input,key_buf);
      snprintf(key_buf,KEYBUF_SIZE-1,"%s.extent.t",this_chunk);
      *t_extent= mri_get_int(Input,key_buf);
      calc_sizes(this_chunk, dimstr, fast_blksize, slow_blksize);
      return 1;
    }
    else if (tloc) {
      if (cl_z_extent==0) {
	if (debug_flag) 
	  fprintf(stderr,"chunk <%s> has bare t dimension\n",
		  this_chunk);
	return 0;
      }
      else {
	long tdim;
	snprintf(key_buf,KEYBUF_SIZE-1,"%s.extent.t",this_chunk);
	tdim= mri_get_int(Input,key_buf);
	if (tdim<=0)
	  Abort("%s: chunk %s has invalid t extent %d!\n",
		progname,this_chunk,tdim);
	if (tdim % cl_z_extent) {
	  if (verbose_flag || debug_flag)
	    fprintf(stderr,
		    "chunk <%s> has bare t dim incommensurate with requested z dimension\n",
		    this_chunk);
	  return 0;
	}
	else {
	  *z_extent= cl_z_extent;
	  *t_extent= tdim/cl_z_extent;
	  calc_sizes(this_chunk, dimstr, fast_blksize, slow_blksize);
	  if (debug_flag) 
	    fprintf(stderr,"mapping %ld times -> %ld slices * %ld times\n",
		    tdim, *z_extent, *t_extent);
	  return 1;
	}
      }
    }
    else {
      return 0;
    }
  }
  else Abort("%s: chunk %s has no dimension key!\n",progname,this_chunk);

}

static int pickDirection(const char* chunk)
{
  int chunkDir;
  int chunkDirSet= checkChunkForReorderDir(chunk, &chunkDir);
  if (chunkDirSet) {
    if (cl_direction_set && (chunkDir != cl_direction)) {
      if (force_flag) {
	fprintf(stderr,"chunk <%s> direction %d overridden!\n",
		chunk, chunkDir);
	return cl_direction;
      }
      else
	Abort("%s: chunk <%s> wants an inconsistent reordering direction!\n",
	      progname,chunk);
    }
    return chunkDir;
  }
  else {
    if (cl_direction_set) return cl_direction;
    else return DEFAULT_DIRECTION;
  }
}

static int pickReorderPattern(const char* chunk, char* reorder_pattern, 
			      int size)
{
  const char* chunkPattern= checkChunkForReorderPattern(chunk);
  if (chunkPattern) {
    if (cl_reorder_pattern_set && strcmp(chunkPattern,cl_reorder_pattern)) {
      if (force_flag) {
	fprintf(stderr,"chunk <%s> pattern %s overridden!\n",
		chunk, chunkPattern);
	strncpy(reorder_pattern,cl_reorder_pattern,size);
	reorder_pattern[size-1]= '\0';
      }
      else
	Abort("%s: chunk <%s> wants an inconsistent reordering pattern!\n",
	      progname,chunk);
    }
    strncpy(reorder_pattern,chunkPattern,size);
    reorder_pattern[size-1]= '\0';
  }
  else {
    if (cl_reorder_pattern_set) 
      strncpy(reorder_pattern,cl_reorder_pattern,size);
    else
      strncpy(reorder_pattern,DEFAULT_PATTERN,size);
    reorder_pattern[size-1]= '\0';
  }
}

static void scan_fold_chunk(const char* this_chunk) {
  long fast_blksize;
  long long slow_blksize;
  long z_extent;
  long t_extent;

  if (fold_this(this_chunk, &fast_blksize, &z_extent, &t_extent,
		&slow_blksize)) {
    char reorder_pattern[KEYBUF_SIZE];
    int direction;
    int* indexMap= NULL;
    direction= pickDirection(this_chunk);
    pickReorderPattern(this_chunk, reorder_pattern, KEYBUF_SIZE);
    if (verbose_flag) 
      fprintf(stderr,"scan_foldolating chunk <%s> with %s %s\n",
	      this_chunk,(direction?"pattern":"inverse of pattern"),
	      reorder_pattern);
    if (direction)
      indexMap= slp_generateSlicePatternTable(z_extent, reorder_pattern);
    else
      indexMap= slp_generateInvertedSlicePatternTable(z_extent, 
						      reorder_pattern);

    restructure_tags(this_chunk, fast_blksize, z_extent, t_extent, 
		     slow_blksize,direction,reorder_pattern, indexMap);
    transfer_data(this_chunk, fast_blksize, z_extent, t_extent, slow_blksize,
		  indexMap);
    data_changed= 1;
    free(indexMap);
  }
  else {
    /* Chunk copied correctly in initial dataset copy */
    if (verbose_flag) 
      fprintf(stderr,"Not folding chunk <%s>\n",this_chunk);
  }
}

int main( int argc, char* argv[] ) 
{

  char infile[512], outfile[512];
  char* this_key;
  char creorder[512];

  progname= argv[0];

  /* Print version number */
  Message( "# %s\n", rcsid );

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Deprecate old options */

  if (cl_present( "r" ))
     Abort ("Option r(slice reordering) has been replaced by reorder|reo.  Please see help file.\n");
  if (cl_present( "s" ))
     Abort ("Option s(number of slices/z extent) has been replaced by zdm.  Please see help file.\n");

  /* Get params and filenames */
  if (cl_get("zdm", "%option %d",&cl_z_extent)) {
    if (cl_z_extent<=0) {
      fprintf(stderr,"%s: number of slices must be greater than zero.\n",argv[0]);
      Help("usage");
      exit(-1);
    }
  }
  else cl_z_extent= 0;

  cl_reorder_pattern_set= cl_get( "reorder|reo", "%option %s", creorder );
  cl_direction_set= cl_get("direction|dir", "%option %d", &cl_direction);

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
  force_flag= cl_present("force");
  verbose_flag= cl_present("verbose|ver|v");
  debug_flag= cl_present("debug|dbg");
  if (cl_cleanup_check()) {
    int i;
    fprintf(stderr,"%s: invalid argument in command line:\n    ",argv[0]);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(-1);
  }

  /*** End command-line parsing ***/
  
  /* Open input dataset */
  if( !strcmp( infile, outfile ) )
    Abort( "%s: Input and output files must be distinct.",argv[0] );
  Input = mri_open_dataset( infile, MRI_READ );

  /* Parse the reordering info from the command line */
  if (cl_reorder_pattern_set) {
    if (!strcmp(creorder,"0") || !strcmp(creorder,"F") 
	|| !strcmp(creorder,"FALSE") || !strcmp(creorder,"f")
	|| !strcmp(creorder,"false") || !strcmp(creorder,"False")) {
      strncpy(cl_reorder_pattern,"sequential",sizeof(cl_reorder_pattern));
      cl_reorder_pattern[sizeof(cl_reorder_pattern)-1]= '\0';
    }
    else {
      if (!strcmp(creorder,"1") || !strcmp(creorder,"T") 
	  || !strcmp(creorder,"TRUE") || !strcmp(creorder,"t")
	  || !strcmp(creorder,"true") || !strcmp(creorder,"True")) {
	strncpy(cl_reorder_pattern,DEFAULT_PATTERN,sizeof(cl_reorder_pattern));
	cl_reorder_pattern[sizeof(cl_reorder_pattern)-1]= '\0';
      }
      else {
	int* testTbl= slp_generateSlicePatternTable(2,creorder);
	if (!testTbl)
	  Abort("%s: unrecognized slice pattern name <%s>!\n",
		progname,creorder);
	free(testTbl);
	strncpy(cl_reorder_pattern,creorder,sizeof(cl_reorder_pattern));
	cl_reorder_pattern[sizeof(cl_reorder_pattern)-1]= '\0';
      }
    }
  }
  else {
    const char* pattern= NULL;
    if (pattern=checkChunkForReorderPattern("images")) {
      if (verbose_flag) 
	fprintf(stderr,"Reorder pattern set from <images> chunk\n");
    }
    else if (pattern=checkChunkForReorderPattern("samples")) {
      if (verbose_flag) 
	fprintf(stderr,"Reorder pattern set from <samples> chunk\n");
    }
    if (pattern) {
      strncpy(cl_reorder_pattern, pattern, sizeof(cl_reorder_pattern));
      cl_reorder_pattern[sizeof(cl_reorder_pattern)-1]= '\0';
      cl_reorder_pattern_set= 1;
    }
  }

  if (!cl_direction_set) {
    if (checkChunkForReorderDir("images",&cl_direction)) {
      if (verbose_flag) 
	fprintf(stderr,"Reorder direction set from <images> chunk\n");
      cl_direction_set= 1;
    }
    else if (checkChunkForReorderDir("samples",&cl_direction)) {
      if (verbose_flag) 
	fprintf(stderr,"Reorder direction set from <samples> chunk\n");
      cl_direction_set= 1;
    }
  }
  if (cl_direction_set) {
    if (cl_direction != 0 && cl_direction != 1)
      Abort("%s: direction must be 0 or 1; it is %d\n",progname,cl_direction);
  }

  /* Open output dataset */
  Output = mri_copy_dataset( outfile, Input );
  hist_add_cl( Output, argc, argv );

  /* Walk through the input handling each chunk in turn */
  mri_iterate_over_keys(Input);
  while ((this_key= mri_next_key(Input)) != NULL) {
    if (!strcmp(mri_get_string(Input,this_key),"[chunk]")) {
      scan_fold_chunk(this_key);
    }
  }

  /* Write and close data-sets */
  mri_close_dataset( Input );
  mri_close_dataset( Output );
  
  Message( "#      Scan folding complete.\n" );
  if (!data_changed) 
    Message("#      Warning: input and output datasets identical!\n");

  return 0;
}

