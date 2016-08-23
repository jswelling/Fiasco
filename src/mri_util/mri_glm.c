/************************************************************
 *                                                          *
 *  mri_glm.c                                               *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1998 Department of Statistics,         *
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
 *  Original programming by Joel Welling, 6/98              *
 ************************************************************/
/*************************************************************

  DESCRIPTION OF MRI_GLM

  mri_glm takes a pgh MRI dataset of type t... or vt... (for
  vector v of length 2) and performs multiple regression on
  each group of data, treating it as a separate time series.

**************************************************************/

/*****************
 * Notes-
 * -Can I still 'skip zero term' in the SSQR copy-out if it's not zero?
 *  Or does it need to contribute to SSE?
 * -checking return values for real input:
 * --b values (estimates) look OK
 * --variances look OK
 * --covariances: need to skip first row and column
 * --nobs
 * --SSTO: OK
 * --skip ave term
 * --synthesize SSE
 * --skip first SSR
 * --ortho measure looks OK
 * -checking return values for complex input:
 * --b values (estimates) look OK
 * --variances look OK
 * --covariances have extra first row and colum
 *
 *****************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "misc.h"
#include "stdcrg.h"
#include "slist.h"

#define KEYBUF_SIZE 512

static char rcsid[] = "$Id: mri_glm.c,v 1.42 2008/04/29 22:16:10 welling Exp $";

typedef struct mrifile_struct {
  MRI_Dataset *ds;
  char* fname;
  int rw;
  int refCount;
  struct mrifile_struct* next;
  struct mrifile_struct* prev;
} MRIFile;

typedef struct mrichunk_struct {
  MRIFile* file;
  char* chunk;
  char* dimstr;
  char* datatype;
  long long offset;
  long long length;
} MRIChunk;

static char* progname;

static MRIFile* openMRIFiles= NULL;

static Regressor* gbl_r= NULL;

/* Some command line flags.  In glm.h, GLM_TYPE is the last index used */
#define MRIGLM_COMPLEX   GLM_COMPLEX
#define MRIGLM_RESIDUALS GLM_RESIDUALS
#define MRIGLM_VAR       GLM_VARIANCES
#define MRIGLM_COVAR     GLM_COVARIANCES
#define MRIGLM_SSQR      GLM_SSQR
#define MRIGLM_ORTHO     GLM_ORTHO
#define MRIGLM_DEBUG     GLM_DEBUG
#define MRIGLM_DEVIANCE  GLM_DEVIANCE
#define MRIGLM_TYPE      GLM_TYPE
#define MRIGLM_ISTDV     (GLM_TYPE+1)
#define MRIGLM_VERBOSE   (GLM_TYPE+2)
#define MRIGLM_COUNTS    (GLM_TYPE+3)
#define MRIGLM_SCALE     (GLM_TYPE+4)
#define MRIGLM_CONSTFAC  (GLM_TYPE+5)

/* Must include all the MRIGLM_ values */
static int optionIndexArray[]= {
  MRIGLM_COMPLEX,
  MRIGLM_RESIDUALS,
  MRIGLM_VAR,
  MRIGLM_COVAR,
  MRIGLM_SSQR,
  MRIGLM_ORTHO,
  MRIGLM_DEBUG,
  MRIGLM_DEVIANCE,
  MRIGLM_TYPE,
  MRIGLM_ISTDV,
  MRIGLM_VERBOSE,
  MRIGLM_COUNTS,
  MRIGLM_SCALE,
  MRIGLM_CONSTFAC
};

/* This contains *only* the names of options beyond MRIGLM_TYPE! */
static char* optionIndexNames[]= {
  "istdv", "verbose", "counts", "scale", "include_const_factor"
};

#define NOPTIONS (sizeof(optionIndexArray)/sizeof(int))
#define NMRIGLMOPTIONS (NOPTIONS-((int)MRIGLM_TYPE+1))

static int gblOptionState[NOPTIONS];

#define OPT_NOT_ALLOWED 0
#define OPT_ALLOWED 1
#define OPT_REQUIRED 2
/* Columns are MRIGLM_ISTDV, MRIGLM_VERBOSE, MRIGLM_COUNTS, MRIGLM_SCALE,
 *             MRIGLM_CONSTFAC */
static int optionIsAllowedOrRequired[GLM_N_TYPES][NMRIGLMOPTIONS]= {
  {OPT_ALLOWED,     OPT_ALLOWED, OPT_NOT_ALLOWED, OPT_ALLOWED,
   OPT_ALLOWED},     /*LLSQ */
  {OPT_NOT_ALLOWED, OPT_ALLOWED, OPT_REQUIRED,   OPT_NOT_ALLOWED,
   OPT_ALLOWED}, /*LOGISTIC*/
  {OPT_NOT_ALLOWED, OPT_ALLOWED, OPT_NOT_ALLOWED,OPT_NOT_ALLOWED,
   OPT_ALLOWED}  /*POISSON*/
};

#define SETOPT(opt,val) { gblOptionState[opt]=val; }
#define GETOPT(opt) (gblOptionState[opt])

static void safe_copy(char* str1, char* str2) {
  strncpy(str1, str2, KEYBUF_SIZE);
  str1[KEYBUF_SIZE-1]= '\0';
}

static void safe_concat(char* str1, char* str2) {
  strncat(str1, str2, (KEYBUF_SIZE-strlen(str1))-1);
}

static int has_fname(char* fullname) {
  char* colon_offset;
  colon_offset= strchr(fullname,':');
  if (colon_offset==fullname) return 0;
  else return 1;
}

static int has_chunkname(char* fullname) {
  char* colon_offset;
  colon_offset= strchr(fullname,':');
  if (colon_offset && *(colon_offset+1)) return 1;
  else return 0;
}

static void parse_fname(char* fullname_in, 
			char* fname_default, char* fname_out,
			char* chunk_default, char* chunk_out) {
  char* colon_offset;
  colon_offset= strchr(fullname_in,':');
  if (!colon_offset) {
    /* filename only */
    if (!chunk_default)
      Abort("%s: chunk name required for %s!\n",progname,fullname_in);
    safe_copy(fname_out, fullname_in);
    safe_copy(chunk_out, chunk_default);
  }
  else if (colon_offset==fullname_in) {
    /* chunkname only */
    if (!fname_default)
      Abort("%s: filename required for %s!\n",progname,fullname_in);
    safe_copy(fname_out, fname_default);
    safe_copy(chunk_out, fullname_in+1); /* skip colon */
  }
  else {
    safe_copy(fname_out, fullname_in);
    fname_out[colon_offset-fullname_in]= '\0'; /* clip off the chunk */
    if (*(colon_offset+1)) safe_copy(chunk_out, colon_offset+1);
    else safe_copy(chunk_out, chunk_default);
  }
}

static void trim_off_mri(char* fname_trimmed)
{
  /* This just removes possible .mri suffices */
  if (strlen(fname_trimmed)>=4 
      && !strcmp(fname_trimmed+strlen(fname_trimmed)-4,".mri"))
    fname_trimmed[strlen(fname_trimmed)-4]= '\0';
}

static int fname_matches(MRIChunk* ch, char* fname2_raw, char* fname2_default)
{
  char fname_trimmed[512], chunk[512];
  char* fname1= ch->file->fname;
  parse_fname(fname2_raw, fname2_default, fname_trimmed,
	      "junkchunk", chunk);
  trim_off_mri(fname_trimmed);
  return (!strcmp(fname1,fname_trimmed));
}

static int fname_and_chunk_match(MRIChunk* ch, char* fname2_raw, 
				 char* fname2_default, char* chunk2_default)
{
  char fname_trimmed[512], chunk[512];
  char* fname1= ch->file->fname;
  char* chunk1= ch->chunk;
  parse_fname(fname2_raw, fname2_default, fname_trimmed,
	      chunk2_default, chunk);
  trim_off_mri(fname_trimmed);
  return (!strcmp(fname1,fname_trimmed) && !strcmp(chunk1,chunk));
}

static char* mriFileAccessString( int rw )
{
  switch (rw) {
  case MRI_READ: return "read";
  case MRI_WRITE: return "write";
  case MRI_MODIFY: return "modify";
  case MRI_MODIFY_DATA: return "modify_data";
  }
  return "***unknown access code!***";
}

static MRIFile* mriFile_copy( MRIFile* orig, char* newName )
{
  char fname_trimmed[512];
  MRIFile* thisMRIFile;
  char* here;

  safe_copy(fname_trimmed, newName);
  /* Worry about possible chunk name in newName */
  here= strchr(fname_trimmed,':');
  if (here != NULL) *here= '\0';
  trim_off_mri(fname_trimmed);
  thisMRIFile= openMRIFiles;
  while (thisMRIFile != NULL) {
    if (!strcmp(fname_trimmed,thisMRIFile->fname)) {
      Abort("%s: MRIFile_copy: can't copy <%s> as <%s>; target is already open!\n",
		 progname, orig->fname, fname_trimmed);
    }
    thisMRIFile= thisMRIFile->next;
  }
  /* OK, we ran off the end of the list of existing files.
   * This is a new file so we can create it. */
  if (!(thisMRIFile= (MRIFile*)malloc(sizeof(MRIFile))))
    Abort("%s: unable to allocate %d bytes!\n",sizeof(MRIFile));
  if (GETOPT(MRIGLM_DEBUG))
    fprintf(stderr,"Copying <%s> to <%s>\n",orig->fname,fname_trimmed);
  if (!(thisMRIFile->ds= mri_copy_dataset(fname_trimmed, orig->ds)))
    Abort("%s: cannot copy %s!\n",progname,fname_trimmed);
  thisMRIFile->fname= strdup(fname_trimmed);
  thisMRIFile->rw= MRI_WRITE;
  thisMRIFile->refCount= 0;
  thisMRIFile->next= openMRIFiles;
  if (openMRIFiles) openMRIFiles->prev= thisMRIFile;
  thisMRIFile->prev= NULL;
  openMRIFiles= thisMRIFile;
  return thisMRIFile;
}

static MRIFile* mriFile_open(char* fname, int rw)
{
  char fname_trimmed[512];
  MRIFile* thisMRIFile;

  safe_copy(fname_trimmed, fname);
  trim_off_mri(fname_trimmed);
  thisMRIFile= openMRIFiles;
  while (thisMRIFile != NULL) {
    if (!strcmp(fname_trimmed,thisMRIFile->fname)) {
      if (thisMRIFile->rw == rw) {
	return thisMRIFile;
      }
      else Abort("%s: MRIFile_open: can't open <%s> for %s; already open for %s!\n",
		 progname, fname_trimmed, mriFileAccessString(rw),
		 mriFileAccessString(thisMRIFile->rw));
    }
    thisMRIFile= thisMRIFile->next;
  }
  /* OK, we ran off the end of the list of existing files.
   * This is a new file so we create it. */
  if (!(thisMRIFile= (MRIFile*)malloc(sizeof(MRIFile))))
    Abort("%s: unable to allocate %d bytes!\n",sizeof(MRIFile));
  if (GETOPT(MRIGLM_DEBUG))
    fprintf(stderr,"Opening <%s>\n",fname_trimmed);
  if (!(thisMRIFile->ds= mri_open_dataset(fname_trimmed,rw)))
    Abort("%s: cannot open %s for %s!\n",progname,fname_trimmed,
	  mriFileAccessString(rw));
  thisMRIFile->fname= strdup(fname_trimmed);
  thisMRIFile->rw= rw;
  thisMRIFile->refCount= 0;
  thisMRIFile->next= openMRIFiles;
  if (openMRIFiles) openMRIFiles->prev= thisMRIFile;
  thisMRIFile->prev= NULL;
  openMRIFiles= thisMRIFile;
  return thisMRIFile;
}

static void mriFile_ref( MRIFile* file )
{
  file->refCount++;
  if (GETOPT(MRIGLM_DEBUG))
    fprintf(stderr,"Referencing %s; refCount now %d\n",
	    file->fname, file->refCount);
}

static void mriFile_unref( MRIFile* file )
{
  file->refCount--;
  if (GETOPT(MRIGLM_DEBUG))
    fprintf(stderr,"Unreferencing %s; refCount now %d\n",
	    file->fname,file->refCount);
  if (file->refCount<=0) {
    MRIFile* thisMRIFile= openMRIFiles;
    while (thisMRIFile != NULL) {
      if (thisMRIFile == file) break;
      thisMRIFile= thisMRIFile->next;
    }
    if (thisMRIFile != file) 
      Abort("%s: internal error: unknown MRI file!\n",progname);
    if (GETOPT(MRIGLM_DEBUG)) fprintf(stderr,"Closing <%s>\n",file->fname);
    if (file->next) file->next->prev= file->prev;
    if (file->prev) file->prev->next= file->next;
    mri_close_dataset(file->ds);
    free(file->fname);
    free(file);
  }
}

static MRIChunk* mriChunk_open_file(MRIFile *file, char* chunk_in, int rw) 
{
  char keybuf[KEYBUF_SIZE];
  char tbuf[64];
  char* crunner;
  long long l;
  MRIChunk* f;
  int created= 0;

  if (!(f=(MRIChunk*)malloc(sizeof(MRIChunk))))
    Abort("%s: unable to allocate %d bytes!\n",progname,sizeof(MRIChunk));

  f->file= file;
  mriFile_ref(f->file);
  f->offset= 0;
  f->chunk= strdup(chunk_in);

  if( !mri_has( f->file->ds, f->chunk ) ) {
    switch (rw) {
    case MRI_READ:
    case MRI_MODIFY_DATA:
      Abort( "%s did not find expected chunk %s in dataset <%s>!\n", progname, 
	     f->chunk, f->file->fname );
      break;
    case MRI_WRITE:
    case MRI_MODIFY:
      mri_create_chunk( f->file->ds, f->chunk );
      created= 1;
      break;
    default:
      Abort( "%s: mriChunk_open: internal error: unknown rw code %d!\n",
	     progname, rw );
    }
  }

  if (created) { 
    f->dimstr= NULL;
    f->datatype= NULL;
    f->length= 0;
  }
  else {
    safe_copy(keybuf, f->chunk);
    safe_concat(keybuf, ".datatype");
    if ( !mri_has( f->file->ds, keybuf ) ) 
      Abort( "%s needs %s info for input file <%s>!\n", progname, keybuf, 
	     f->file->fname );
    f->datatype= strdup(mri_get_string( f->file->ds, keybuf ));

    safe_copy(keybuf, f->chunk);
    safe_concat(keybuf, ".dimensions");
    if ( !mri_has( f->file->ds, keybuf ) ) 
      Abort( "%s needs %s info for input file <%s>!\n", progname, keybuf, 
	     f->file->fname );
    f->dimstr= strdup(mri_get_string( f->file->ds, keybuf ));
    
    /* Figure out how many instances of this vector */
    l= 1;
    for (crunner= f->dimstr; *crunner; crunner++) {
      safe_copy(keybuf,f->chunk);
      sprintf(tbuf,".extent.%c",*crunner);
      safe_concat(keybuf,tbuf);
      if (!mri_has(f->file->ds, keybuf)) 
	Abort("%s: file <%s> has dimension %c but no %s tag.",
	      progname,*crunner,keybuf);
      else {
	l *= mri_get_int( f->file->ds, keybuf );
      }
    }
    
    f->length= l;
  }

  if (GETOPT(MRIGLM_DEBUG))
    fprintf(stderr,"Opened file <%s>, chunk <%s>, dimstr %s, length %ld\n",
	    f->file->fname, f->chunk, f->dimstr, (long)f->length);

  return f;
}

static MRIChunk* mriChunk_open(char* fname_raw, 
			       char* default_fname_in, char* default_chunk_in, 
			       int rw) 
{
  char fname[512], chunkname[512];
  MRIChunk* f;

  parse_fname(fname_raw, default_fname_in, fname, default_chunk_in, chunkname);
  f= mriChunk_open_file( mriFile_open(fname,rw), chunkname, rw );
  if (GETOPT(MRIGLM_DEBUG)) fprintf(stderr,"Opened chunk %s:%s\n",
			     f->file->fname, f->chunk);
  return f;
}

static void mriChunk_rescan(MRIChunk* f) 
{
  char keybuf[KEYBUF_SIZE];
  char tbuf[64];
  char* dimstr;
  char* crunner;
  long long l;

  if (strlen(f->chunk)>KEYBUF_SIZE-64)
    Abort( "%s: chunk name too long!\n", progname );
  if( !mri_has( f->file->ds, f->chunk ) )
    Abort( "%s did not find expected chunk %s in dataset <%s>!\n", progname, 
	   f->chunk, f->file->fname );
  safe_copy(keybuf, f->chunk);
  safe_concat(keybuf, ".dimensions");
  if ( !mri_has( f->file->ds, keybuf ) ) 
    Abort( "%s needs %s info for input file <%s>!\n", progname, keybuf, 
	   f->file->fname );
  dimstr= mri_get_string( f->file->ds, keybuf );
  if (f->dimstr) free(f->dimstr);
  f->dimstr= strdup(dimstr);

  /* Figure out how many instances of this vector */
  l= 1;
  for (crunner= dimstr; *crunner; crunner++) {
    safe_copy(keybuf,f->chunk);
    sprintf(tbuf,".extent.%c",*crunner);
    safe_concat(keybuf,tbuf);
    if (!mri_has(f->file->ds, keybuf)) 
      Abort("%s: file <%s> has dimension %c but no %s tag.",
	    progname,*crunner,keybuf);
    else {
      l *= mri_get_int( f->file->ds, keybuf );
    }
  }
  
  f->length= l;

  if (GETOPT(MRIGLM_DEBUG)) 
    fprintf(stderr,"Rescanned chunk %s:%s; length now %ld\n",
	    f->file->fname, f->chunk, (long)f->length);
}

static void mriChunk_close(MRIChunk* f) 
{
  if (GETOPT(MRIGLM_DEBUG)) fprintf(stderr,"Closing chunk %s:%s\n",
			     f->file->fname, f->chunk);
  free( f->chunk );
  free( f->dimstr );
  mriFile_unref( f->file );
}

static void mriChunk_read(double* dest, MRIChunk* f, int n) 
{
  int nread= 0;
  int thisblock;
  long long offset;

  offset= f->offset;
  while (n>nread) {
    thisblock= (n-nread > f->length - offset) ? 
      (f->length - offset) : n-nread;
    if (!(mri_read_chunk(f->file->ds, f->chunk, thisblock, offset,
			 MRI_DOUBLE, dest+nread)))
      Abort("%s: bad mri_read_chunk on infile <%s>; length %d, offset %d",
	    f->file->fname,thisblock,offset);
    offset += thisblock;
    if (offset>=f->length) offset= 0;
    nread += thisblock;
  }
  if (GETOPT(MRIGLM_DEBUG)) fprintf(stderr,"Read %d from %s:%s at %ld\n",
			     n, f->file->fname, f->chunk, (long)f->offset);
}

static void mriChunk_write(double* source, MRIChunk* f, int n) 
{
  int nwritten= 0;
  int thisblock;
  long long offset;

  offset= f->offset;
  while (n>nwritten) {
    thisblock= (n-nwritten > f->length - offset) ? 
      (f->length - offset) : n-nwritten;
    mri_write_chunk(f->file->ds, f->chunk, thisblock, offset,
		    MRI_DOUBLE, source+nwritten);
    offset += thisblock;
    if (offset>=f->length) offset= 0;
    nwritten += thisblock;
  }
  if (GETOPT(MRIGLM_DEBUG)) fprintf(stderr,"Wrote %d to %s:%s at %ld\n",
			     n, f->file->fname, f->chunk, (long)f->offset);
}

static void mriChunk_advance(MRIChunk* f, long long n)
{
  f->offset= (f->offset + n) % f->length;
  if (GETOPT(MRIGLM_DEBUG)) fprintf(stderr,"Advance %s:%s; offset now %ld\n",
			     f->file->fname, f->chunk, (long)f->offset);
}

static int mriChunk_get_extent(MRIChunk* ch, char c)
{
  char keybuf[KEYBUF_SIZE];
  char tbuf[64];

  if (!strchr(ch->dimstr,c)) return 1; /* missing dim equiv to extent 1 */

  safe_copy(keybuf,ch->chunk);
  sprintf(tbuf,".extent.%c",c);
  safe_concat(keybuf,tbuf);

  /* file was checked for existence of this extent when file was opened */
  if (mri_has(ch->file->ds,keybuf)) 
    return mri_get_int(ch->file->ds, keybuf); 
  else return 1;
}

static int process_infile_keys(MRIChunk* Input,
			       int* is_complex,
			       int* tseries_length, long* tseries_per_slice,
			       int* zdim, char** err_msg)
{
  char keybuf[KEYBUF_SIZE];
  char* dim_runner;
  long n;
  int z_expected= 0;

  if (*(Input->dimstr)=='t') {
    *is_complex= 0;
    dim_runner= (Input->dimstr)+1;
  }
  else if (*(Input->dimstr)=='v' && *((Input->dimstr)+1)=='t') {
    int v_extent;
    v_extent= mriChunk_get_extent(Input,'v');
    dim_runner= (Input->dimstr)+2;
    if (v_extent==1) {
      /* no worries; ignore it */
      *is_complex= 0;
    }
    else if (v_extent==2) {
      *is_complex= 1;
    }
    else {
      *err_msg= "input chunk v extent is not 1 or 2, so not real or complex";
      return 0;
    }
  }
  else {
    *err_msg= "input chunk is not of type t... or vt...";
    return 0;
  }
  
  if (*((Input->dimstr) + strlen((Input->dimstr)) - 1) == 'z') {
    *zdim= mriChunk_get_extent(Input,'z');
    z_expected= 1;
  }
  else {
    z_expected= 0;
    *zdim= 1;
  }
  
  *tseries_length= mriChunk_get_extent(Input,'t');

  /* Count the number of time series. dim_runner set above to point
   * to first dim past t. 
   */
  n= 1;
  while (*dim_runner) {
    if (*dim_runner=='z') {
      if (!z_expected || (*(dim_runner+1))) {
	*err_msg= "input dimension string contains z, but it is not last";
	return 0;
      }
      dim_runner++;
    }
    else {
      n *= mriChunk_get_extent(Input,*dim_runner);
      dim_runner++;
    }
  }
  *tseries_per_slice= n;

  return 1;
}

static void create_param_chunk(MRIChunk* Params, MRIChunk* Input, 
			       int complex_flag, int nparams)
{
  char keybuf[KEYBUF_SIZE];
  char* dimstr;
  char dimbuf[KEYBUF_SIZE];

  safe_copy(keybuf,Params->chunk);
  safe_concat(keybuf,".datatype");
  if (!strcmp(Input->datatype,"float64"))
    mri_set_string( Params->file->ds, keybuf, "float64" );
  else
    mri_set_string( Params->file->ds, keybuf, "float32" );

  dimstr= Input->dimstr;
  dimstr= strchr(dimstr,'t') + 1; /* skip to past t */
  safe_copy(dimbuf,"v");
  safe_concat(dimbuf,dimstr);
  safe_copy(keybuf,Params->chunk);
  safe_concat(keybuf,".dimensions");
  mri_set_string(Params->file->ds, keybuf, dimbuf);
  safe_copy(keybuf,Params->chunk);
  safe_concat(keybuf,".extent.v");
  mri_set_int(Params->file->ds,keybuf, nparams );

  while (*dimstr) {
    int dim;

    dim= mriChunk_get_extent(Input,*dimstr);
    safe_copy(keybuf,Params->chunk);
    safe_concat(keybuf,".extent.");
    if (strlen(keybuf)<KEYBUF_SIZE-1) strncat(keybuf,dimstr,1);
    else Abort("%s: key buffer overflow!\n",progname);
    mri_set_int(Params->file->ds,keybuf,dim);
    
    dimstr++;
  }

  safe_copy(keybuf,Params->chunk);
  safe_concat(keybuf,".file");
  if (Params->file->refCount==1) 
    mri_set_string( Params->file->ds, keybuf, ".dat" );
  else mri_set_string( Params->file->ds, keybuf, ".glm_dat" );

  mriChunk_rescan(Params);
}

static int test_scale( MRIChunk* Scale,
		       int tseries_length, long* scale_tseries_per_slice,
		       int zdim, char** err_msg )
{
  const char* dimstr= Scale->dimstr;
  long tseries_per_slice= 1;

  /* Does the provided Scale dataset have the correct shape? */
  if (dimstr[0]=='v') {
    if (mriChunk_get_extent(Scale,'v')!=1) {
      *err_msg= "scale chunk v extent must be 1";
      return 0;
    }
    dimstr++; /* skip over unit dimension v */
  }

  if (dimstr[0]=='t') {
    if (mriChunk_get_extent(Scale,'t') != tseries_length) {
      *err_msg= "scale chunk t extent does not match input";
      return 0;
    }
  }
  else {
    *err_msg= "Dimension string for scale chunk is not (v)t...";
    return 0;
  }

  if (dimstr[strlen(Scale->dimstr)-1]=='z') {
    if (mriChunk_get_extent(Scale,'z')!=zdim) {
      *err_msg= "scale chunk z extent does not match input";
      return 0;
    }
  }
  else {
    if (zdim != 1) {
      *err_msg= "scale chunk lacks expected z dimension";
      return 0;
    }
  }

  dimstr++; /* skip over t */
  tseries_per_slice= 1;
  while (*dimstr && (*dimstr != 'z')) {
    tseries_per_slice *= mriChunk_get_extent(Scale,*dimstr);
    dimstr++;
  }
  
  *scale_tseries_per_slice= tseries_per_slice;
  return 1;
}

static int test_istdv( MRIChunk* Istdv,
		       int tseries_length, int zdim, char** err_msg )
{
  /* Does the provided Istdv dataset have the correct shape? */
  if (strcmp(Istdv->dimstr,"t") && strcmp(Istdv->dimstr,"vt")
      && strcmp(Istdv->dimstr,"vtz") && strcmp(Istdv->dimstr,"tz")) {
    *err_msg= "Dimension string for istdv chunk is not (v)t(z)";
    return 0;
  }

  if (strchr(Istdv->dimstr,'t')) {
    if (mriChunk_get_extent(Istdv,'t') != tseries_length) {
      *err_msg= "istdv chunk t extent does not match input";
      return 0;
    }
  }

  if (strchr(Istdv->dimstr,'v')) {
    if (mriChunk_get_extent(Istdv,'v') != (GETOPT(MRIGLM_COMPLEX)?2:1)) {
      *err_msg= "istdv chunk v extent must match input";
      return 0;
    }
  }

  if (strchr(Istdv->dimstr,'z')) {
    if (mriChunk_get_extent(Istdv,'z') != zdim) {
      *err_msg= "istdv chunk z extent does not match input";
      return 0;
    }
  }
  return 1;
}

static int test_counts( MRIChunk* Counts, MRIChunk* Input,
		       int tseries_length, char** err_msg )
{
  const char* dimstr= Counts->dimstr;
  const char* input_dimstr= Input->dimstr;

  /* Does the provided Counts dataset have the correct shape? */
  if (dimstr[0]=='v') {
    if (mriChunk_get_extent(Counts,'v')!=1) {
      *err_msg= "counts chunk v extent must be 1";
      return 0;
    }
    dimstr++; /* skip over unit dimension v */
  }

  if (input_dimstr[0]=='v') input_dimstr++;

  while (*dimstr || *input_dimstr) {
    if (*dimstr && *input_dimstr) {
      int ext= mriChunk_get_extent(Counts,*dimstr);
      int in_ext= mriChunk_get_extent(Input,*dimstr);
      if (*dimstr==*input_dimstr) {
	if (ext==in_ext) {
	  dimstr++;
	  input_dimstr++;
	}
	else if (ext==1) {
	  /* Trivial dim in Counts */
	  dimstr++;
	}
	else if (in_ext==1) {
	  /* Trivial dim in Infile */
	  input_dimstr++;
	}
      }
      else {
	if (ext==1) dimstr++;
	else if (in_ext==1) input_dimstr++;
	else {
	  *err_msg= "non-trivial mismatched dimensions";
	  return 0;
	}
      }
    }
    else {
      if (*dimstr) {
	/* Hopefully this is just a trailing dim of extent 1 */
	int ext= mriChunk_get_extent(Counts,*dimstr);
	if (ext==1) dimstr++;
	else {
	  *err_msg= "Counts input has trailing non-trivial dimensions";
	  return 0;
	}
      }
      else {
	int in_ext= mriChunk_get_extent(Input,*dimstr);
	if (in_ext==1) input_dimstr++;
	else {
	  *err_msg= "Counts input lacks trailing non-trivial dimensions";
	  return 0;
	}
      }
    }
  }
  return 1;
}

static MRIChunk* open_and_test_factor( char* fname_raw, 
				       char* default_fname, char* default_chunk,
				       int tseries_length, int zdim, int tdim,
				       char** err_msg)
{
  MRIChunk* factor= mriChunk_open( fname_raw, default_fname, default_chunk,
				   MRI_READ );

  if (strcmp(factor->dimstr,"t") && strcmp(factor->dimstr,"vt")
      && strcmp(factor->dimstr,"vtz") && strcmp(factor->dimstr,"tz")) {
    *err_msg= "Dimension string for a factor is not (v)t(z)";
    return 0;
  }
  
  if (strchr(factor->dimstr,'t')) {
    if (mriChunk_get_extent(factor,'t') != tseries_length) {
      *err_msg= "factor chunk t extent does not match input";
      return 0;
    }
  }

  if (strchr(factor->dimstr,'z')) {
    if (mriChunk_get_extent(factor,'z') != zdim) {
      *err_msg= "factor chunk z extent does not match input";
      return 0;
    }
  }
  return factor;
}

static void check_factor_scratch(int fac_dv, int fac_block, 
				 long* factor_scratch_size, 
				 double** factor_scratch)
{
  if (*factor_scratch_size<fac_dv*fac_block) {
    if (*factor_scratch) {
      free(*factor_scratch);
      *factor_scratch= NULL;
      *factor_scratch_size= 0;
    }
    if (!(*factor_scratch= (double*)malloc(fac_dv*fac_block*sizeof(double))))
      Abort("%s: unable to allocate %d doubles!\n",progname,fac_block);
    *factor_scratch_size= fac_dv*fac_block;
  }
}

static void load_factors(double* factor_raw,
			 int nfactors, int fac_block, SList* FactorList)
{
  int ifac;
  int ifile;
  int irun;
  int fac_stride;
  int fac_matrix_offset= 0;
  int fac_index;
  static double* factor_scratch= NULL;
  static long factor_scratch_size= 0;
  int fac_dv= 1;
  int v;

  fac_stride= nfactors;

  /* Synthesize the constant factor if necessary */
  if (GETOPT(MRIGLM_CONSTFAC)) {
    fac_index= fac_matrix_offset;
    for (irun=0; irun<fac_block; irun++) {
      factor_raw[fac_index]= 1.0;
      fac_index += fac_stride;
    }
    fac_matrix_offset += 1;
  }

  /* load and store factors */
  slist_totop(FactorList);
  while (!slist_atend(FactorList)) {
    MRIChunk* thisFactor= (MRIChunk*)slist_next(FactorList);
    fac_dv= mriChunk_get_extent(thisFactor,'v');
    check_factor_scratch(fac_dv, fac_block, 
			 &factor_scratch_size, &factor_scratch);
    mriChunk_read(factor_scratch, thisFactor, fac_dv*fac_block);
    mriChunk_advance(thisFactor, fac_dv*fac_block);
    fac_index= fac_matrix_offset;
    for (irun=0; irun<fac_block; irun++) {
      for (v=0; v<fac_dv; v++) {
	factor_raw[fac_index+v]= factor_scratch[irun*fac_dv + v];
      }
      fac_index += fac_stride;
    }
    fac_matrix_offset += fac_dv;
  }
}

static void copy_and_apply_scale( double* data_out, const double* data, 
				  const double* scale, 
				  int complex_flag, int n, int ncols )
{
  int icol;
  int i;
  int offset;
  
  for (icol=0; icol<ncols; icol++) {
    if (complex_flag) {
      offset= 2*icol;
      for (i=0; i<n; i++) {
	data_out[offset]= data[offset]*scale[i];
	data_out[offset+1]= data[offset+1]*scale[i];
	offset += 2*ncols;
      }
    }
    else {
      offset= icol;
      for (i=0; i<n; i++) {
	data_out[offset]= data[offset] * scale[i];
	offset += ncols;
      }
    }
  }
}

static void apply_scale( double* data, 
			 const double* scale, 
			 int complex_flag, int n, int ncols )
{
  int icol;
  int i;
  int offset;
  
  for (icol=0; icol<ncols; icol++) {
    if (complex_flag) {
      offset= 2*icol;
      for (i=0; i<n; i++) {
	data[offset] *= scale[i];
	data[offset+1] *= scale[i];
	offset += 2*ncols;
      }
    }
    else {
      offset= icol;
      for (i=0; i<n; i++) {
	data[offset] *= scale[i];
	offset += ncols;
      }
    }
  }
}

static void apply_unscale( double* data, double* scale, 
			   int complex_flag, int n, int ncols )
{
  int icol;
  int i;
  int offset;
  
  for (icol=0; icol<ncols; icol++) {
    if (complex_flag) {
      offset= 2*icol;
      for (i=0; i<n; i++) {
	if (scale[i]==0.0) {
	  data[offset]= 0.0;
	  data[offset+1]= 0.0;
	}
	else {
	  data[offset] /= scale[i];
	  data[offset+1] /= scale[i];
	}
	offset += 2*ncols;
      }
    }
    else {
      offset= icol;
      for (i=0; i<n; i++) {
	if (scale[i]==0.0) data[offset]= 0.0;
	else data[offset] /= scale[i];
	offset += ncols;
      }
    }
  }
}

static int pack_out_missing( double* data_in,
			     unsigned char** missing,
			     int tdim, int complex_flag, int nblocks,
			     int z, double* data_out )
{
  int iblock;
  int t;
  double* in;
  double* out;
  int n_valid;

  /* This program works on data for which the block index (corresponding
   * to the factor) varies faster than the tdim index (corresponding to
   * observation).
   */
  n_valid= 0;
  for (iblock=0; iblock<nblocks; iblock++) {
    if (complex_flag) {
      in= data_in + (2*iblock);
      out= data_out + (2*iblock);
      for (t=0; t<tdim; t++) {
	if (!missing[t][z]) {
	  *out= *in;
	  *(out+1)= *(in+1);
	  out += 2*nblocks;
	  n_valid++;
	}
	in += 2*nblocks;
      }
    }
    else {
      in= data_in + iblock;
      out= data_out + iblock;
      for (t=0; t<tdim; t++) {
	if (!missing[t][z]) {
	  *out= *in;
	  out += nblocks;
	  n_valid++;
	}
	in += nblocks;
      }
    }
  }

  return n_valid; /* note that this should be doubled for complex! */
}

static void pack_in_missing( double* data_in,
			     unsigned char** missing,
			     int tdim, int complex_flag, int nblocks,
			     int z, double fill_val, double* data_out )
{
  int iblock;
  int t;
  double* in;
  double* out;

  /* This program works on data for which the block index (corresponding
   * to the factor) varies faster than the tdim index (corresponding to
   * observation).
   */
  for (iblock=0; iblock<nblocks; iblock++) {
    if (complex_flag) {
      in= data_in + (2*iblock);
      out= data_out + (2*iblock);
      for (t=0; t<tdim; t++) {
	if (missing[t][z]) {
	  *out= fill_val;
	  *(out+1)= fill_val;
	}
	else {
	  *out= *in;
	  *(out+1)= *(in+1);
	  in += (2*nblocks);
	}
	out += (2*nblocks);
      }
    }
    else {
      in= data_in + iblock;
      out= data_out + iblock;
      for (t=0; t<tdim; t++) {
	if (missing[t][z]) *out= fill_val;
	else {
	  *out= *in;
	  in += nblocks;
	}
	out += nblocks;
      }
    }
  }
}

static int n_output_params( nfactors )
{
  int result= glm_n_params(gbl_r,nfactors);
  if (GETOPT(MRIGLM_COMPLEX)) {
    if (GETOPT(MRIGLM_COVAR) && GETOPT(MRIGLM_CONSTFAC)) {
      /* Compensate for presence of constant term.  glm expects 
       * (2*nfactors^2), but that nfactors includes const term.
       */
      result += ((nfactors-2)*(nfactors-2) - (nfactors*nfactors));
    }
    if (GETOPT(MRIGLM_SSQR)) {
      /* We are injecting a term for number of samples */
      result += 1;
      if (GETOPT(MRIGLM_CONSTFAC)) {
	/* We have one less SSE term */
	result -= 1;
      }
    }
  }
  else {
    if (GETOPT(MRIGLM_COVAR) && GETOPT(MRIGLM_CONSTFAC)) {
      /* Compensate for presence of constant term.  glm expects 
       * nfactors^2, but that nfactors includes const term.
       */
      result += ((nfactors-1)*(nfactors-1) - (nfactors*nfactors));
    }
    if (GETOPT(MRIGLM_SSQR)) {
      /* We are injecting a term for number of samples */
      result += 1;
      if (GETOPT(MRIGLM_CONSTFAC)) {
	/* We have one less SSE term */
	result -= 1;
      }
    }
  }
  return result;
}

static void format_output_params( double* params_out, double* params,
				  double* istdv_proj, int nsamples,
				  int nparams, int nparams_out, int nfactors )
{
  int i,j;
  int o_offset= 0;
  int i_offset= 0;

  /* reality check, against n_output_params() */
  if (nparams_out != n_output_params(nfactors)) {
    Abort("%s: internal error: failed reality check!\n",progname);
  }

  /* Estimates, including the constant term if present */
  for (i=0; i<nfactors; i++) params_out[o_offset++]= params[i_offset++];

  /* Variances, including the constant term if present */
  if (GETOPT(MRIGLM_VAR)) {
    if (istdv_proj) {
      for (i=0; i<nfactors; i++) 
	params_out[o_offset++]= 
	  params[i_offset++] + istdv_proj[i]*istdv_proj[i];
    }
    else {
      for (i=0; i<nfactors; i++) params_out[o_offset++]= params[i_offset++];
    }
  }

  /* Covariances.  We must exclude the constant term if it is present. */
  if (GETOPT(MRIGLM_COVAR)) {
    if (GETOPT(MRIGLM_CONSTFAC)) {
      if (istdv_proj) {
	for (j=0; j<nfactors; j++) i_offset++;
	for (i=1; i<nfactors; i++) {
	  i_offset++;
	  for (j=1; j<nfactors; j++) {
	    params_out[o_offset]= params[i_offset++];
	    if (i==j) params_out[o_offset] += istdv_proj[i-1]*istdv_proj[i-1];
	    o_offset++;
	  }
	}
      }
      else {
	for (j=0; j<nfactors; j++) i_offset++;
	for (i=1; i<nfactors; i++) {
	  i_offset++;
	  for (j=1; j<nfactors; j++) 
	    params_out[o_offset++]= params[i_offset++];
	}
      }
    }
    else {
      if (istdv_proj) {
	for (i=0; i<nfactors; i++) 
	  for (j=0; j<nfactors; j++) {
	    params_out[o_offset]= params[i_offset++];
	    if (i==j) params_out[o_offset] += istdv_proj[i]*istdv_proj[i];
	    o_offset++;
	  }
      }
      else {
	for (i=0; i<nfactors; i++) 
	  for (j=0; j<nfactors; j++) 
	    params_out[o_offset++]= params[i_offset++];
      }
    }
  }

  if (GETOPT(MRIGLM_SSQR)) {
    double ssto;
    double ssr;
    double sse;
    int sse_out_offset;

    params_out[o_offset++]= (double)nsamples;
    ssto= params_out[o_offset++]= params[i_offset++];
    i_offset++; /* skip average term */
    sse_out_offset= o_offset++; /* come back later for SSE */
    ssr= 0.0;
    if (GETOPT(MRIGLM_CONSTFAC)) {
      ssr += params[i_offset++]; /* skip constant term */
      for (i=1; i<nfactors; i++) {
	/* copy SSR terms */
	ssr += params[i_offset];
	params_out[o_offset++]= params[i_offset++];
      }      
    }
    else {
      for (i=0; i<nfactors; i++) {
	/* copy SSR terms */
	ssr += params[i_offset];
	params_out[o_offset++]= params[i_offset++];
      }
    }
    sse= ssto - ssr;
    if (sse<0.0) sse= 0.0; /* in case of numerical error */
    params_out[sse_out_offset]= sse;
  }

  if (GETOPT(MRIGLM_ORTHO)) params_out[o_offset++]= params[i_offset++];

  if (GETOPT(MRIGLM_DEVIANCE)) params_out[o_offset++]= params[i_offset++];
}

static void format_output_params_complex( double* params_out, double* params,
					  double* istdv_proj,
					  int nsamples, int nparams, 
					  int nparams_out, int nfactors )
{
  int i;
  int j;
  int v;
  int o_offset= 0;
  int i_offset= 0;

  /* reality check, against n_output_params() */
  if (nparams_out != n_output_params(nfactors)) {
    Abort("%s: internal error: failed reality check!\n",progname);
  }

  for (i=0; i<2*nfactors; i++) params_out[o_offset++]= params[i_offset++];
  if (GETOPT(MRIGLM_VAR)) {
    if (istdv_proj) {
      for (i=0; i<2*nfactors; i++) 
	params_out[o_offset++]= 
	  params[i_offset++] + istdv_proj[i]*istdv_proj[i];
    }
    else {
      for (i=0; i<2*nfactors; i++) params_out[o_offset++]= params[i_offset++];
    }
  }

  /* Covariances.  We must exclude the constant term if present. */
  if (GETOPT(MRIGLM_COVAR)) {
    if (GETOPT(MRIGLM_CONSTFAC)) {
      if (istdv_proj) {
	for (j=0; j<nfactors; j++) i_offset += 2;
	for (i=0; i<nfactors; i++) {
	  i_offset += 2;
	  for (j=1; j<nfactors; j++) 
	    for (v=0; v<2; v++) {
	      params_out[o_offset]= params[i_offset++];
	      if (i==j) 
		params_out[o_offset] += istdv_proj[2*i+v]*istdv_proj[2*i+v];
	      o_offset++;
	    }
	}
      }
      else {
	for (j=0; j<nfactors; j++) i_offset += 2;
	for (i=0; i<nfactors; i++) {
	  i_offset += 2;
	  for (j=1; j<nfactors; j++) 
	    for (v=0; v<2; v++) {
	      params_out[o_offset++]= params[i_offset++];
	    }
	}
      }
    }
    else {
      if (istdv_proj) {
	for (i=0; i<nfactors; i++) 
	  for (j=0; j<nfactors; j++) 
	    for (v=0; v<2; v++) {
	      params_out[o_offset]= params[i_offset++];
	      if (i==j) 
		params_out[o_offset] += istdv_proj[2*i+v]*istdv_proj[2*i+v];
	      o_offset++;
	    }
      }
      else {
	for (i=0; i<nfactors; i++) 
	  for (j=0; j<nfactors; j++) 
	    for (v=0; v<2; v++) {
	      params_out[o_offset++]= params[i_offset++];
	    }
      }
    }
  }
  
  if (GETOPT(MRIGLM_SSQR)) {
    double ssto;
    double ssr;
    int sse_out_offset;

    params_out[o_offset++]= (double)nsamples;
    ssto= params_out[o_offset++]= params[i_offset++];
    sse_out_offset= o_offset++; /* come back later for SSE */
    ssr= 0.0;
    if (GETOPT(MRIGLM_CONSTFAC)) {
      ssr += params[i_offset++]; /* skip constant term */
      for (i=1; i<nfactors; i++) {
	/* copy SSR terms */
	ssr += params[i_offset];
	params_out[o_offset++]= params[i_offset++];
      }
    }
    else {
      for (i=0; i<nfactors; i++) {
	/* copy SSR terms */
	ssr += params[i_offset];
	params_out[o_offset++]= params[i_offset++];
      }
    }
    params_out[sse_out_offset]= ssto - ssr;
  }

  if (GETOPT(MRIGLM_ORTHO)) params_out[o_offset++]= params[i_offset++];

  if (GETOPT(MRIGLM_DEVIANCE)) params_out[o_offset++]= params[i_offset++];
}

static void calc_residuals_direct( int complex_flag, 
				   double mean, double mean_i,
				   double* parameters, int nfactors, 
				   double* input_unpacked, 
				   double* factors_unpacked, 
				   double* tseries_unpacked, 
				   int tseries_length_unpacked )
{
  /* Given a set of parameter estimates, this routine calculates 
   * residuals in the original, unpacked representation.
   */
  int i; /* observations */
  int a; /* factors */

  if (complex_flag) {
    for (i=0; i<tseries_length_unpacked; i++) {
      tseries_unpacked[2*i]= input_unpacked[2*i] - mean;
      tseries_unpacked[2*i+1]= input_unpacked[2*i+1] - mean_i;
      for (a=0; a<nfactors; a++) {
	tseries_unpacked[2*i] -= 
	  parameters[2*a]*factors_unpacked[i*nfactors+a];
	tseries_unpacked[2*i+1] -= 
	  parameters[2*a+1]*factors_unpacked[i*nfactors+a];
      }
    }
  }
  else {
    for (i=0; i<tseries_length_unpacked; i++) {
      tseries_unpacked[i]= input_unpacked[i] - mean;
      for (a=0; a<nfactors; a++) {
	tseries_unpacked[i] -= 
	  parameters[a]*factors_unpacked[i*nfactors+a];
      }
    }
  }  
}

static void get_next_scale_block(MRIChunk* Scale, double* scale_unpacked,
				 double* scale, unsigned char** missing,
				 int z, int complex_flag,
				 int scale_block, int tseries_length)
{
  int scale_length_packed;
  double sum= 0.0;
  int i;
  mriChunk_read(scale_unpacked, Scale, scale_block);
  mriChunk_advance(Scale, scale_block);
  scale_length_packed= pack_out_missing(scale_unpacked, 
					missing, tseries_length, 
					complex_flag, 1, z, scale);
  /* Normalize the scale to an average weight of 1; 
   * the user may not have done so.
   */
  for (i=0; i<scale_length_packed; i++) sum += scale[i];
  if (sum!=0.0) {
    /* sum can be zero outside some region of interest, for example.
     * We'll just leave them as we found them in this case; the
     * user is not likely to be interested.
     */
    double n_by_sum= ((double)scale_length_packed)/sum;
    for (i=0; i<scale_length_packed; i++) scale[i] *= n_by_sum;
  }
}

static void fit_glm(MRIChunk* Input, 
		    unsigned char** missing, 
		    int nfactors,
		    SList* FactorList,
		    MRIChunk* Output,
		    MRIChunk* Params,
		    MRIChunk* Scale,
		    MRIChunk* Counts,
		    MRIChunk* Istdv,
		    int complex_flag, int tseries_length, 
		    long tseries_per_slice, long scale_tseries_per_slice,
		    int zdim)
{
  int z;
  int i;
  double* factors_unpacked;
  double* factors;
  double* factors_scaled;
  double* factor_means;
  double* tseries_unpacked;
  double* tseries;
  double* scale_unpacked;
  double* scale;
  double* counts_unpacked;
  double* counts;
  double* istdv_unpacked;
  double* istdv;
  double* istdv_proj;
  double* parameters;
  double* parameters_out;
  double* input_unpacked;
  int complex_fac;
  int in_block;
  int out_block;
  int par_block;
  int fac_block;
  int scale_block;
  int istdv_block;
  int nparams;
  int nparams_out;
  int retcode;
  int tseries_length_packed;
  int factors_length_packed;
  double mean= 0.0;
  double mean_i= 0.0;
  double mean_variance= 0.0;
  double mean_variance_i= 0.0;

  /* Here we walk through the data, actually fitting the GLM.  This
   * routine contains the main loop of the program.
   */

  complex_fac= (complex_flag ? 2 : 1); /* for convenience */
  nparams= glm_n_params(gbl_r,nfactors);
  nparams_out= n_output_params(nfactors);

  in_block= tseries_length*complex_fac;
  out_block= tseries_length*complex_fac;
  par_block= nparams_out;
  fac_block= tseries_length;
  scale_block= tseries_length;
  istdv_block= tseries_length*complex_fac;

  /* Allocate memory */
  if (!(input_unpacked= (double*)malloc(in_block*sizeof(double))))
    Abort("%s: unable to allocate %d doubles!\n",progname,in_block);
  if (!(factors_unpacked= (double*)malloc(tseries_length*nfactors*
					 sizeof(double))))
    Abort("%s: unable to allocate %d doubles!\n",
	  progname,tseries_length*nfactors);
  if (!(factors= (double*)malloc(tseries_length*nfactors*sizeof(double))))
    Abort("%s: unable to allocate %d doubles!\n",
	  progname,tseries_length*nfactors);
  if (!(factors_scaled= (double*)malloc(tseries_length*nfactors
					*sizeof(double))))
    Abort("%s: unable to allocate %d doubles!\n",
	  progname,tseries_length*nfactors);
  if (!(factor_means= (double*)malloc(nfactors*sizeof(double))))
    Abort("%s: unable to allocate %d doubles!\n", progname, nfactors);
  if (!(tseries_unpacked= (double*)malloc(tseries_length*complex_fac
					 *sizeof(double))))
    Abort("%s: unable to allocate %d doubles!\n",
	  progname,tseries_length*complex_fac);
  if (!(tseries= (double*)malloc(tseries_length*complex_fac*sizeof(double))))
    Abort("%s: unable to allocate %d doubles!\n",
	  progname,tseries_length*complex_fac);
  if (!(parameters= (double*)malloc(nparams*sizeof(double))))
    Abort("%s: unable to allocate %d doubles!\n",progname,nparams);
  if (!(parameters_out= (double*)malloc(nparams_out*sizeof(double))))
    Abort("%s: unable to allocate %d doubles!\n",progname,nparams_out);

  if (Scale) {
    if (!(scale_unpacked= (double*)malloc(scale_block*sizeof(double))))
      Abort("%s: unable to allocate %d doubles!\n",progname,scale_block);
    if (!(scale= (double*)malloc(scale_block*sizeof(double))))
      Abort("%s: unable to allocate %d doubles!\n",progname,scale_block);
  }
  else {
    scale_unpacked= NULL;
    scale= NULL;
  }
  if (Istdv) {
    if (!(istdv_unpacked= (double*)malloc(istdv_block*sizeof(double))))
      Abort("%s: unable to allocate %d doubles!\n",progname,istdv_block);
    if (!(istdv= (double*)malloc(istdv_block*sizeof(double))))
      Abort("%s: unable to allocate %d doubles!\n",progname,istdv_block);
    if (!(istdv_proj= (double*)malloc(nfactors*complex_fac*sizeof(double))))
      Abort("%s: unable to allocate %d doubles!\n",progname,
	    nfactors*complex_fac*sizeof(double));
  }
  else {
    istdv_unpacked= NULL;
    istdv= NULL;
    istdv_proj= NULL;
  }

  for (z=0; z<zdim; z++) {

    load_factors(factors_unpacked, nfactors, fac_block, FactorList);
    
    factors_length_packed= pack_out_missing(factors_unpacked,
					    missing, tseries_length, 0,
					    nfactors, z, factors);

    if (Scale) {
      if (scale_tseries_per_slice==1) {
	get_next_scale_block(Scale, scale_unpacked, scale, missing,
			     complex_flag, z, scale_block, tseries_length);
	copy_and_apply_scale( factors_scaled, factors, scale, complex_flag,
			      factors_length_packed/nfactors, nfactors );
      }
      else {
	/* Nothing- we'll set up scaled factors later */
      }
    }
    else bcopy( factors, factors_scaled, 
		nfactors*tseries_length*sizeof(double) );

    if (Istdv) {
      mriChunk_read(istdv_unpacked, Istdv, istdv_block);
      mriChunk_advance(Istdv, istdv_block);
      /* the packed length returned by this call to pack_out_missing
       * is guaranteed to match tseries_length_packed (below) because
       * it depends only on tseries_length and z.
       */
      (void)pack_out_missing(istdv_unpacked,
			     missing, tseries_length, 
			     complex_flag, 1, z, istdv);
    }

    for (i=0; i<tseries_per_slice; i++) {
      /* load input time series data */
      mriChunk_read(input_unpacked, Input, in_block);
      mriChunk_advance(Input, in_block);

      /* pack out missing data */
      tseries_length_packed= pack_out_missing(input_unpacked,
					      missing, tseries_length, 
					      complex_flag, 1, z, tseries);

      /* Apply scale if necessary */
      if (Scale) {
	if (scale_tseries_per_slice==1) {
	  /* Scale data for this slice has been loaded, and factors are
	   * already scaled.
	   */
	  apply_scale( tseries, scale, complex_flag, 
		       tseries_length_packed, 1 );
	}
	else {
	  /* Scale data for this voxel must be loaded, and both the time
	   * series data and the factor data must be scaled.
	   */
	  get_next_scale_block(Scale, scale_unpacked, scale, missing,
			       z, complex_flag, scale_block, tseries_length);
	  copy_and_apply_scale( factors_scaled, factors, scale, complex_flag,
				factors_length_packed/nfactors, nfactors );
	  apply_scale( tseries, scale, complex_flag, 
		       tseries_length_packed, 1 );
	}
      }

      /* Fit the data.  This may replace input with residuals. */
      if ((retcode= glm_fit(gbl_r,tseries, factors_scaled, counts, parameters, 
			    tseries_length_packed, nfactors)) != 0) {
	Warning(1,"%s: error in glm_fit: %s!\n",
		progname,glm_error_msg());
      }
      if (Istdv) {
	if ((retcode= glm_normproject(gbl_r,istdv, istdv_proj, 
				      tseries_length_packed, 
				      nfactors)) != 0) {
	  int i;
	  Warning(1,"%s: error in glm_normproject: %s!\n",
		  progname,glm_error_msg());
	  for (i=0; i<nfactors; i++) istdv_proj[i]= 0.0;
	}
      }

      if (Output) { /* residuals requested */
	/* Unscale if necessary */
	if (Scale) {
	  /* We can't simply divide out the scale because some scale
	   * values may have been zero.  We must wastefully and annoyingly 
	   * recalculate the residuals.  That's the only way I know of
	   * to get valid residuals for Y elements that had zero weight.
	   *
	   * Life is a sea of troubles.  Sigh.
	   */
	  calc_residuals_direct( complex_flag, mean, mean_i,
				 parameters, nfactors, 
				 input_unpacked,
				 factors_unpacked, 
				 tseries_unpacked, tseries_length );
	}
	else {
	  /* pack in 0's for missing data */
	  pack_in_missing(tseries, missing, tseries_length, 
			  complex_flag, 1, z, 0.0, tseries_unpacked);
	}

	/* write the time series */
	mriChunk_write(tseries_unpacked, Output, out_block);
	mriChunk_advance(Output, out_block);
      }
      if (complex_flag)
	format_output_params_complex(parameters_out, parameters, istdv_proj,
				     tseries_length_packed,
				     nparams, nparams_out, nfactors);
      else 
	format_output_params(parameters_out, parameters, istdv_proj,
			     tseries_length_packed, nparams, nparams_out, 
			     nfactors);
      mriChunk_write(parameters_out, Params, par_block);
      mriChunk_advance(Params, par_block);
    }
    if (GETOPT(MRIGLM_VERBOSE)) {
      if (zdim<=30) {
	Message("Finished slice %d\n",z);
      }
      else {
	/* Sometimes it's necessary to have very large zdim, because each
	 * voxel needs its own factor matrix.  We need a less verbose 
	 * form of output in that case.
	 */
	static int old_frac= 0;
#ifdef never
	int frac = ( 10*z )/zdim;
	if (frac != old_frac) {
	  Message("Finished %d%%\n",10*frac);
	  old_frac= frac;
	}
#endif
	int frac = ( 100*z )/zdim;
	if (frac != old_frac) {
	  Message("Finished %d%%\n",frac);
	  old_frac= frac;
	}
      }
    }
  }

  free( (void*)input_unpacked );
  free( (void*)factors_unpacked );
  free( (void*)factors );
  free( (void*)factors_scaled );
  free( (void*)factor_means );
  free( (void*)tseries_unpacked );
  free( (void*)tseries );
  free( (void*)parameters );
  free( (void*)parameters_out );
  if (scale_unpacked) free( (void*)scale_unpacked );
  if (scale) free( (void*)scale );
  if (istdv_unpacked) free( (void*)istdv_unpacked );
  if (istdv) free( (void*)istdv );
  if (istdv_proj) free( (void*)istdv_proj );
}

int main( int argc, char *argv[] ) 
{
  MRIChunk *Input = NULL, *Output = NULL, *Params = NULL;
  MRIChunk *Scale = NULL, *Counts= NULL, *Istdv = NULL;
  SList* FactorList= NULL;
  int i;
  char infile_raw[512], outfile_raw[512], paramfile_raw[512];
  char scalefile_raw[512], countsfile_raw[512], istdvfile_raw[512];
  SList* rawFactorNameList= NULL;
  char thisFactorRaw[512];
  int nfactors= 0;
  int nFactorFiles= 0;
  int tseries_length;
  long tseries_per_slice;
  long scale_tseries_per_slice;
  int zdim;
  char* err_msg;
  char keybuf[KEYBUF_SIZE];
  unsigned char** missing= NULL;

  progname= argv[0];

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Deprecate old options */

  if (cl_present( "input|i" ))
     Abort ("Option input|i has been replaced by multiple infiles format.  Please see help file.\n");
  if (cl_present( "o" ))
     Abort ("Option o has been expanded to output|out.  Please see help file.\n");
  if (cl_present( "R" ))
     Abort ("Option R has been replaced by variance|var.  Please see help file.\n");
  if (cl_present( "s" ))
     Abort ("Option s has been expanded to sum_of_squares|ssqr.  Please see help file.\n");
  if (cl_present( "S" ))
     Abort ("Option S has been replaced by scale|sca.  Please see help file.\n");
  if (cl_present( "param|p" ))
     Abort ("Option param|p has been replaced by estimates|est|e.  Please see help file.\n");
   if (cl_present( "factor|f" ))
     Abort ("Option factor|f has been replaced by multiple infiles format.  Please see help file.\n");
  if (cl_present( "d" ))
     Abort ("Option d has been expanded to debug|deb.  Please see help file.\n");
  if (cl_present( "istdv" ))
     Abort ("Option istdv has been replaced by stdv.  Please see help file.\n");

  /* Get flags */
  SETOPT(MRIGLM_VAR,cl_present("variance|var"));
  SETOPT(MRIGLM_COVAR,cl_present("covariance|cov"));
  SETOPT(MRIGLM_SSQR,cl_present("sum_of_squares|ssqr"));
  SETOPT(MRIGLM_ORTHO,cl_present("orthogonality|rth"));
  SETOPT(MRIGLM_VERBOSE,cl_present("verbose|ver|v"));
  SETOPT(MRIGLM_DEBUG,cl_present("debug|deb"));
  SETOPT(MRIGLM_TYPE,GLM_TYPE_LLSQ); /* This is the default */
  SETOPT(MRIGLM_CONSTFAC,!cl_present("noconst"));
  if (cl_present("logistic"))
    SETOPT(MRIGLM_TYPE,GLM_TYPE_LOGISTIC);
  if (cl_present("poisson"))
    SETOPT(MRIGLM_TYPE,GLM_TYPE_POISSON);

  /* Get filenames */
  SETOPT(MRIGLM_RESIDUALS,cl_get("output|out", "%option %s", outfile_raw)); 
  if (GETOPT(MRIGLM_RESIDUALS) && strchr(outfile_raw,':')) {
    fprintf(stderr,"%s: Output file chunk name not supported.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  
  if (!cl_get("estimates|est|e", "%option %s", paramfile_raw)) {
    if (GETOPT(MRIGLM_RESIDUALS)) {
      char* here;
      safe_copy(paramfile_raw, outfile_raw);
      here= strchr(paramfile_raw,':');
      if (here != NULL) *here= '\0';
      safe_concat(paramfile_raw,":glm");
    }
    else {
      fprintf(stderr,"%s: neither output file nor parameter file specified.\n",
	      argv[0]);
      Help("usage");
      exit(-1);
    }
  }

  SETOPT(MRIGLM_SCALE,cl_get("scale|sca", "%option %s", scalefile_raw));
  SETOPT(MRIGLM_ISTDV,cl_get("stdv", "%option %s", istdvfile_raw));
  SETOPT(MRIGLM_COUNTS,cl_get("counts", "%option %s", countsfile_raw));

  if(!cl_get("", "%s", infile_raw)) {
    fprintf(stderr, "%s: Input file name not given.\n", argv[0]);
    Help("usage");
    exit(-1);
  }
  if (!has_fname(infile_raw)) {
    fprintf(stderr,"%s: Input has no filename component.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  if (has_chunkname(infile_raw)) safe_concat(outfile_raw,strchr(infile_raw,':'));
  else safe_concat(outfile_raw,":images");

  rawFactorNameList= slist_create();
  while(cl_get("", "%s", thisFactorRaw))
    slist_append(rawFactorNameList, strdup(thisFactorRaw));
  nFactorFiles= slist_count(rawFactorNameList);
  if (!nFactorFiles) {
    fprintf(stderr,"%s: At least one factor must be given.\n",argv[0]);
    Help("usage");
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
  if (GETOPT(MRIGLM_VERBOSE)) Message( "# %s\n", rcsid );

  Input= mriChunk_open(infile_raw, NULL, "images", MRI_READ);
  missing= get_missing(Input->file->ds);
  {
    int flg_tmp;
    if (!process_infile_keys( Input, &flg_tmp,
			      &tseries_length, &tseries_per_slice, &zdim, 
			      &err_msg ))
      Abort("%s: %s.\n", argv[0], err_msg);
    SETOPT(MRIGLM_COMPLEX,flg_tmp);
  }

  /* Check for filename compatibility */
  if (GETOPT(MRIGLM_RESIDUALS) && fname_matches(Input, outfile_raw, NULL))
    Abort("%s: infile and outfile must be distinct!\n",argv[0]);

  if (fname_matches(Input,paramfile_raw,"glm"))
    Abort("%s: infile and paramfile must be distinct!\n",argv[0]);

  if (GETOPT(MRIGLM_SCALE) 
      && fname_and_chunk_match(Input,
			       scalefile_raw,NULL,"images"))
    Abort("%s: scale file and chunk are the same as input file and chunk!\n",
	  argv[0]);

  if (GETOPT(MRIGLM_COUNTS) 
      && fname_and_chunk_match(Input,
			       countsfile_raw,NULL,"images"))
    Abort("%s: scale file and chunk are the same as input file and chunk!\n",
	  argv[0]);

  if (GETOPT(MRIGLM_ISTDV)
      && fname_and_chunk_match(Input,
			       istdvfile_raw,NULL,"images"))
    Abort("%s: istdv file and chunk are the same as input file and chunk!\n",
	  argv[0]);

  if (GETOPT(MRIGLM_RESIDUALS)) {
    /* We need to worry about the case where parms and residuals go to
     * the same file.
     */
    MRIFile* outMRIFile= mriFile_copy( Input->file, outfile_raw );
    Output= mriChunk_open_file( outMRIFile, Input->chunk, MRI_WRITE );
    hist_add_cl( Output->file->ds, argc, argv );
    if (fname_and_chunk_match(Output,paramfile_raw,
			      Output->file->fname, "glm"))
      Abort("%s: output file and chunk are the same as parameter file and chunk!\n",
	    argv[0]);
    Params= mriChunk_open(paramfile_raw, Output->file->fname, "glm", 
			  MRI_WRITE);
    if (!mri_has(Params->file->ds,"history.1")) {
      /* avoid duplication if residuals and params go to same dataset */
      hist_add_cl( Params->file->ds, argc, argv );
    }
  }
  else {
    Params= mriChunk_open(paramfile_raw, NULL, "glm", MRI_WRITE);
    hist_add_cl( Params->file->ds, argc, argv );
    Output= NULL;
  }

  if (GETOPT(MRIGLM_SCALE)) {
    Scale= mriChunk_open( scalefile_raw, NULL, "images", MRI_READ );
    if (!test_scale(Scale, tseries_length, &scale_tseries_per_slice, 
		    zdim, &err_msg)) 
      Abort("%s: %s\n",argv[0], err_msg);
    if (scale_tseries_per_slice != 1 
	&& scale_tseries_per_slice != tseries_per_slice)
      Abort("%s: scale dataset is not commensurate with input\n",
	    argv[0]);
  }

  if (GETOPT(MRIGLM_ISTDV)) {
    Istdv= mriChunk_open( istdvfile_raw, NULL, "images", MRI_READ );
    if (!test_istdv(Istdv, tseries_length, zdim, &err_msg)) 
      Abort("%s: %s\n",argv[0], err_msg);

  }

  if (GETOPT(MRIGLM_COUNTS)) {
    Counts= mriChunk_open( countsfile_raw, NULL, "images", MRI_READ );
    if (!test_counts(Counts, Input, tseries_length, &err_msg)) 
      Abort("%s: %s\n",argv[0], err_msg);

  }

  /* Check for factor validity and open factors */
  nfactors= 0;
  if (GETOPT(MRIGLM_CONSTFAC)) nfactors += 1;
  slist_totop(rawFactorNameList);
  FactorList= slist_create();
  while (!slist_empty(rawFactorNameList)) {
    char* thisRawFactorName= (char*)slist_pop(rawFactorNameList);
    MRIChunk* thisFactor= NULL;

    if (GETOPT(MRIGLM_RESIDUALS) 
	&& fname_matches(Output, thisRawFactorName, NULL))
      Abort("%s: factor %s can't come from output file!\n",
	    argv[0], thisRawFactorName);
    if (fname_matches(Params, thisRawFactorName, NULL))
      Abort("%s: factor %s can't come from parameter file!\n",
	    argv[0], thisRawFactorName);
    if (GETOPT(MRIGLM_SCALE)
	&& fname_matches(Scale, thisRawFactorName, NULL))
      Abort("%s: factor %s can't come from scale file!\n",
	    argv[0], thisRawFactorName);
    if (GETOPT(MRIGLM_ISTDV)
	&& fname_matches(Istdv, thisRawFactorName, NULL))
      Abort("%s: factor %s can't come from istdv file!\n",
	    argv[0], thisRawFactorName);

    if (fname_and_chunk_match(Input, thisRawFactorName, NULL, "images")) 
      Abort("%s: factor %s comes from same file and chunk as input!\n",
	    argv[0], thisRawFactorName);

    if (GETOPT(MRIGLM_COUNTS)
	&& fname_and_chunk_match(Input, thisRawFactorName, NULL, "images")) 
      Abort("%s: factor %s comes from same file and chunk as counts!\n",
	    argv[0], thisRawFactorName);
    
    slist_totop(FactorList);
    while (!slist_atend(FactorList)) {
      MRIChunk* thatFactor= (MRIChunk*)slist_next(FactorList);
      if (fname_and_chunk_match(thatFactor, thisRawFactorName, NULL, "images"))
	Abort("%s: factor %s is repeated!\n",progname,thisRawFactorName);
    }

    if (!(thisFactor=
	  open_and_test_factor(thisRawFactorName, Input->file->fname, "images",
			       tseries_length, zdim, tseries_length,
			       &err_msg)))
      Abort("%s: factor %s: %s.\n", progname, thisRawFactorName, err_msg);
    slist_append(FactorList,thisFactor);
    nfactors += mriChunk_get_extent(thisFactor,'v');

    free(thisRawFactorName);
  }
  slist_destroy(rawFactorNameList,NULL);
  rawFactorNameList= NULL;
    
  /* Create the global regressor */
  if (!(gbl_r= glm_create_regressor_by_type(GETOPT(MRIGLM_TYPE))))
    Abort("%s: unable to create %s regressor!\n",progname,
	  glm_get_regressor_type_name(GETOPT(MRIGLM_TYPE)));

  /* Verify compatibility of options */
  for (i=0; i<NOPTIONS; i++) {
    int opt= optionIndexArray[i];
    if (GETOPT(opt)) {
      if (opt<MRIGLM_TYPE) {
	if (glm_is_settable(gbl_r, (glm_feature)opt))
	  glm_set(gbl_r, opt, 1);
	else {
	  Abort("%s: option %s is not compatible with %s regression!\n",
		progname,glm_get_feature_name(opt),
		glm_get_regressor_type_name(GETOPT(MRIGLM_TYPE)));
	}
      }
      else if (opt==MRIGLM_TYPE) { /* do nothing */ }
      else {
	int offset= opt-((int)MRIGLM_TYPE + 1);
	if (optionIsAllowedOrRequired[GETOPT(MRIGLM_TYPE)][offset]
	    == OPT_NOT_ALLOWED) {
	  Abort("%s: option %s is not compatible with %s regression!\n",
		progname,optionIndexNames[opt-((int)MRIGLM_TYPE+1)],
		glm_get_regressor_type_name(GETOPT(MRIGLM_TYPE)));
	}
      }
    }
    else {
      if (opt>MRIGLM_TYPE) {
	int offset= opt-((int)MRIGLM_TYPE + 1);
	if (optionIsAllowedOrRequired[GETOPT(MRIGLM_TYPE)][offset]
	    ==OPT_REQUIRED) {
	  Abort("%s: option %s is required with %s regression!\n",
		progname,optionIndexNames[opt-((int)MRIGLM_TYPE+1)],
		glm_get_regressor_type_name(GETOPT(MRIGLM_TYPE)));
	}
      }
    }
  }
  /* Output a report of option status */
  if (GETOPT(MRIGLM_VERBOSE)) {
    Message("Regressor type: %s\n",
	    glm_get_regressor_type_name(GETOPT(MRIGLM_TYPE)));
    Message("Option summary:\n");
    for (i=0; i<NOPTIONS; i++) {
      int opt= optionIndexArray[i];
      char* valstr= ( GETOPT(opt) ? "ON" : "OFF" );
      if (opt<MRIGLM_TYPE) {
	Message("  %s: %s\n",glm_get_feature_name(opt),valstr);
      }
      else if (opt==MRIGLM_TYPE) {
	/* Do nothing */
      }
      else {
	Message("  %s: %s\n",optionIndexNames[opt-((int)MRIGLM_TYPE+1)]
		,valstr);
      }
    }
  }

  /* Set up output datasets */
  create_param_chunk(Params, Input,
		     GETOPT(MRIGLM_COMPLEX), n_output_params(nfactors));

  /* Finally, fit glm */
  fit_glm(Input, missing, nfactors, FactorList, Output, Params, 
	  Scale, Counts, Istdv, GETOPT(MRIGLM_COMPLEX), 
	  tseries_length, tseries_per_slice, 
	  scale_tseries_per_slice, zdim);

  /* Write and close data-sets */
  mriChunk_close(Input);
  if (Scale) mriChunk_close( Scale );
  if (Output) mriChunk_close( Output );
  if (Counts) mriChunk_close( Counts );
  mriChunk_close( Params );
  slist_totop(FactorList);
  while (slist_atend(FactorList)) {
    MRIChunk* thisFactor= (MRIChunk*)slist_next(FactorList);
    mriChunk_close(thisFactor);
  }

  /* Clean up the regressor */
  glm_destroy(gbl_r);
  
  if (GETOPT(MRIGLM_VERBOSE)) Message( "#      Regression complete.\n" );
  exit(0);

}

