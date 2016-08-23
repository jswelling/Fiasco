/************************************************************
 *                                                          *
 *  mri_rpn_math.c                                          *
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
 *  Original programming by Joel Welling, 8/97              *
 ************************************************************/

/* Notes-
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include "mri.h"
#include "fmri.h"
#include "misc.h"
#include "stdcrg.h"

#define MAX_INPUT_FILES 20
#define DEFAULT_CHUNK_NAME "images"
#define MAX_SCRIPT_CHARS 512
#define KEYBUF_SIZE 512

static char rcsid[] = "$Id: mri_rpn_math.c,v 1.49 2005/03/10 01:39:18 welling Exp $";

typedef struct mrifile_struct {
  MRI_Dataset* ds;
  char* fname;
  char* chunk;
  long long offset;
  long long length;
  int is_complex;
} MRIFile;

static char* progname;
static int n_input_files= 0;
static MRIFile Input[MAX_INPUT_FILES];
static char chunkname[512];

static unsigned char** missing= NULL;

static int outfile_flg= 0; /* will we be producing output, or just printing? */
static int verbose_flg= 0;
static int debug_flg= 0;
static int complex_flg= 0; /* read complex from input files */

static void safe_copy(char* str1, char* str2) 
{
  strncpy(str1, str2, KEYBUF_SIZE);
  str1[KEYBUF_SIZE-1]= '\0';
}

static void safe_concat(char* str1, char* str2) 
{
  strncat(str1, str2, (KEYBUF_SIZE-strlen(str1))-1);
}

static void ds_open(MRIFile* f, char* fname, char* chunk_in, int rw) 
{
  char keybuf[KEYBUF_SIZE];
  char tbuf[64];
  char* dimstr;
  char* crunner;
  long long l;

  f->ds= mri_open_dataset(fname,rw);
  f->fname= strdup(fname);
  f->chunk= strdup(chunk_in);
  f->offset= 0;
  f->is_complex= 0;

  if (strlen(f->chunk)>KEYBUF_SIZE-64)
    Abort( "%s: chunk name too long!\n", progname );
  if( !mri_has( f->ds, f->chunk ) )
    Abort( "%s did not find expected chunk %s in dataset <%s>!\n", progname, 
	   f->chunk, fname );
  safe_copy(keybuf, chunkname);
  safe_concat(keybuf, ".dimensions");
  if ( !mri_has( f->ds, keybuf ) ) 
    Abort( "%s needs %s info for input file <%s>!\n", progname, keybuf, 
	   fname );
  dimstr= mri_get_string( f->ds, keybuf );

  /* Figure out how many instances of this vector */
  l= 1;
  for (crunner= dimstr; *crunner; crunner++) {
    safe_copy(keybuf,f->chunk);
    sprintf(tbuf,".extent.%c",*crunner);
    safe_concat(keybuf,tbuf);
    if (!mri_has(f->ds, keybuf)) 
      Abort("%s: file <%s> has dimension %c but no %s tag.",
	    progname,f->fname,*crunner,keybuf);
    else l *= mri_get_int( f->ds, keybuf );
  }
  
  f->length= l;

  /* Is this file complex? */
  if (dimstr[0]=='v') {
    keybuf[KEYBUF_SIZE-1]= '\0';
    snprintf(keybuf,KEYBUF_SIZE-1,"%s.extent.v",f->chunk);
    if (!mri_has(f->ds,keybuf)) 
      Abort("%s: file <%s> has dimension v but no %s tag.\n",
	    progname,f->fname,keybuf);
    if (mri_get_int(f->ds,keybuf)==2) f->is_complex= 1;
    else f->is_complex= 0;
  }
  else f->is_complex= 0;

  if (complex_flg && f->is_complex) f->length /= 2; /* count by pairs */
  if (debug_flg) 
    fprintf(stderr,"%s %s complex\n",fname,(f->is_complex ? "is":"is not"));
}

static void ds_close(MRIFile* f) 
{
  mri_close_dataset( f->ds );
  free( f->chunk );
  free( f->fname );
}

static void ds_read(double* dest, MRIFile* f, long n, long long offset) 
{
  long nread= 0;
  long thisblock;

  offset= offset % f->length;
  while (n>nread) {
    thisblock= (n-nread > f->length - offset) ? 
      (f->length - offset) : n-nread;
    if (!(mri_read_chunk(f->ds, f->chunk, thisblock, offset, MRI_DOUBLE,
			 dest+nread)))
      Abort("%s: bad mri_read_chunk on infile <%s>; length %d, offset %d",
	    f->fname,thisblock,offset);
    offset += thisblock;
    if (offset>=f->length) offset= 0;
    nread += thisblock;
  }
}

static void ds_read_complex(double* dest1, double* dest2, MRIFile* f, long n,
			    long long offset) 
{
  int nread= 0;
  int thisblock;
  int i;
  double* inbuf;

  offset= offset % f->length;
  if (f->is_complex) {
    while (n>nread) {
      thisblock= (n-nread) > (f->length - offset) ? 
	(f->length - offset) : n-nread;
      if (!(inbuf= mri_get_chunk(f->ds, f->chunk, 
				 2*thisblock, 2*offset, MRI_DOUBLE)))
	Abort("%s: bad mri_get_chunk on infile <%s>; length %d, offset %d",
	      f->fname,2*thisblock,2*offset);
      for (i=0; i<thisblock; i++) {
	*(dest1+nread+i)= *(inbuf+2*i);
	*(dest2+nread+i)= *(inbuf+2*i+1);
      }
      offset += thisblock;
      if (offset>=f->length) offset= 0;
      nread += thisblock;
    }
  }
  else {
    int i;
    ds_read(dest1, f, n, offset);
    for (i=0; i<n; i++) dest2[i]= 0.0;
  }
}

static const char* getDimensionsCB(const int which, void* usrHook)
{
  char buf[KEYBUF_SIZE];
  MRIFile* f= Input+which;
  snprintf(buf,sizeof(buf)-1,"%s.dimensions",f->chunk);
  buf[KEYBUF_SIZE-1]= '\0';
  if (!(mri_has(f->ds,buf)))
    Abort("%s: input %s has no %s tag!\n",progname,f->fname,buf);
  return mri_get_string(f->ds,buf);  
}

static const long getDimExtentCB(const int which, const char dim,
				 void* usrHook)
{
  char buf[KEYBUF_SIZE];
  MRIFile* f= Input+which;
  snprintf(buf,sizeof(buf)-1,"%s.extent.%c",f->chunk,dim);
  buf[KEYBUF_SIZE-1]= '\0';
  if (!(mri_has(f->ds,buf)))
    Abort("%s: input %s has no %s tag!\n",progname,f->fname,buf);
  return mri_get_int(f->ds,buf);
}

static void inputCB(const int which, const long n, 
		    const long long offset, double* buf,
		    void* usrHook)
{
  long long myOffset= offset;
  long myN= n;
  ds_read(buf, Input+which, myN, myOffset);
}

static void inputComplexCB(const int which, const long n,
			   const long long offset, double* buf1,
			   double* buf2, void* usrHook)
{
  long long myOffset= offset;
  long myN= n;
  ds_read_complex(buf1, buf2, Input+which, myN, myOffset);
}

static int missingCB(const long z, const long t, void* usrHook)
{
  if (missing!=NULL) return missing[t][z];
  else return 0;
}


int main( int argc, char* argv[] ) 
{
  MRI_Dataset *Output= NULL;
  int vec_length;
  long long voxels_moved;
  int voxels_this_chunk;
  int i;
  char infile[512], outfile[512], scriptfile[512];
  char script[MAX_SCRIPT_CHARS];
  RpnEngine* re= NULL;
  double* out_chunk;
  char keybuf[KEYBUF_SIZE];
  long rand_seed;
  int script_from_file= 0;
  const char* dimstr= NULL;

  progname= argv[0];

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /*** Parse command line ***/

  cl_scan( argc, argv );

  if (cl_present( "f" ))
     Abort ("Option f(expression file) has been replaced by expression|exp.  Please see help file.\n");
  if (cl_present( "r" ))
     Abort ("Option r(seed) has been replaced by seed.  Please see help file.\n");
  if (cl_present( "o" ))
     Abort ("Option o(outfile) has been expanded to outfile|out.  Please see help file.\n");


  /* Get options */
  verbose_flg= cl_present("verbose|ver|v");
  debug_flg= cl_present("debug");
  cl_get("chunk|chu|c", "%option %s[%]", DEFAULT_CHUNK_NAME, chunkname);
  outfile_flg= cl_get("out|outfile", "%option %s", outfile);
  if (!cl_get("seed", "%option %d",&rand_seed)) {
    /* Use milliseconds since the epoch */
    struct timeval tv;
    (void)gettimeofday(&tv, NULL);
    rand_seed= (long)(tv.tv_sec) ^ (long)(tv.tv_usec);
  }
  script_from_file= cl_get("expression|exp", "%option %s", scriptfile);
  complex_flg= cl_present("complex|cpx");

  if (!script_from_file) {
    if (!cl_get("", "%s", script)) {
      fprintf(stderr,"%s: RPN expression to execute not given.\n",argv[0]);
      Help( "usage" );
      exit(-1);
    }
  }

  /* Get filenames */
  if (!cl_get("", "%s", infile)) {
    fprintf(stderr,"%s: Required first input file name not given.\n",argv[0]);
    Help( "usage" );
    exit(-1);
  }
  if( outfile_flg && !strcmp( infile, outfile ) )
    Abort( "Input and output files must be distinct." );
  ds_open(Input+0, infile, chunkname, MRI_READ);

  n_input_files=1;
  while (cl_get("", "%s", infile)) {
    if (n_input_files>=MAX_INPUT_FILES) {
      Abort("%s: too many input files; limit of %d compiled in!\n",
	    argv[0],MAX_INPUT_FILES);
    }
    if( outfile_flg && !strcmp( infile, outfile ) )
      Abort( "Input and output files must be distinct." );
    ds_open(Input+n_input_files, infile, chunkname, MRI_READ);
    n_input_files++;
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

  if (verbose_flg && !script_from_file) Message("script: <%s>\n",script);

  dimstr= getDimensionsCB(0,NULL);
  if (strchr(dimstr,'z')!=NULL && strchr(dimstr,'t')!=NULL) {
    if (mri_has(Input->ds,"missing")) {
      missing= get_missing(Input->ds);
    }
    else {
      int dt= getDimExtentCB(0,'t',NULL);
      int dz= getDimExtentCB(0,'z',NULL);
      int t,z;
      missing= Matrix(dt,dz,unsigned char);
      for (t=0; t<dt; t++)
	for (z=0; z<dz; z++) 
	  missing[t][z]= 0;
    }
  }
  else {
    missing= NULL;
  }

  re= createRpnEngine(n_input_files, NULL, getDimensionsCB, getDimExtentCB,
		      inputCB, inputComplexCB, missingCB);
  rpnSetOutputFlag(re,outfile_flg);
  rpnSetVerbose(re,verbose_flg);
  rpnSetDebug(re,debug_flg);
  rpnSetComplex(re,complex_flg);

  /* Initialize random number generator */
  (void)srand48(rand_seed);
  
  /* Initialize the engine */
  if (!rpnInit(re))
    Abort("%s: %s\n",progname,rpnGetErrorString(re));

  /* Compile the script. */
  if (script_from_file) {
    if (!(rpnCompileFile(re,scriptfile)))
      Abort("%s: error compiling script file <%s>: %s!\n",argv[0],
	    scriptfile,rpnGetErrorString(re));
  }
  else {
    if (!(rpnCompile(re,script)))
      Abort("%s: error compiling expressiong <%s>: %s!\n",argv[0],script,
	    rpnGetErrorString(re));
  }

  /* Set up output parameters */
  if (outfile_flg) {
    Output = mri_copy_dataset( outfile, Input[0].ds );
    hist_add_cl( Output, argc, argv );
    safe_copy(keybuf,chunkname);
    safe_concat(keybuf,".datatype");
    if (!mri_has(Output, keybuf) 
	|| strncmp( mri_get_string( Output, keybuf ), "float", 5 )) {
      /* We will coerce the output datatype to 4-byte floats if it is
       * not one of the floating types.
       */
      mri_set_string( Output, keybuf, "float32" );
    }
    /* If the -complex command line argument is set, we must make sure
     * that the output file is complex.
     */
    if (complex_flg && !Input[0].is_complex) {
      char* dimstr;
      snprintf(keybuf,KEYBUF_SIZE-1,"%s.dimensions",chunkname);
      keybuf[KEYBUF_SIZE-1]= '\0';
      if (!mri_has(Output,keybuf))
	Abort("%s: output file inexplicably has no dimensions string!\n",
	      progname);
      dimstr= mri_get_string(Output,keybuf);
      if (strchr(dimstr,'v')) {
	if (strchr(dimstr,'v') != dimstr) /* v is not first */
	  Abort("%s: output file needs to be complex but already has dim v!\n",
		progname);
	else {
	  snprintf(keybuf,KEYBUF_SIZE-1,"%s.extent.v",chunkname);
	  if (mri_get_int(Output,keybuf)!=1) {
	    Abort("%s: output file needs to be complex but has nontrivial dim v!\n",
		  progname);
	  }
	  else {
	    mri_set_int(Output,keybuf,2);
	  }
	}
      }
      else {
	char valbuf[KEYBUF_SIZE];
	valbuf[KEYBUF_SIZE-1]= '\0';
	snprintf(valbuf,KEYBUF_SIZE-1,"v%s",dimstr);
	mri_set_string(Output,keybuf,valbuf);
	snprintf(keybuf,KEYBUF_SIZE-1,"%s.extent.v",chunkname);
	mri_set_int(Output,keybuf,2);
      }
    }
  }
  else Output= NULL;

  voxels_moved= 0;
  while (voxels_moved < Input[0].length) {
    voxels_this_chunk= ((Input[0].length - voxels_moved) > RPN_CHUNKSIZE) ?
      RPN_CHUNKSIZE : Input[0].length - voxels_moved;
    if (outfile_flg) {
      out_chunk= rpnRun(re, voxels_this_chunk, voxels_moved);
      if (!out_chunk) Abort("%s: execution error after %lld voxels: %s!\n",
			    progname,voxels_moved,rpnGetErrorString(re));
      if (complex_flg) {
	/* Real and imaginary values for an output voxel are in
	 * first and second positions in the stack; out_chunk
	 * actually points to the imaginary part.  Allocate a
	 * buffer on the heap on the first pass through.
	 */
	static double* obuf= NULL;
	long i;
	if (!obuf) {
	  if (!(obuf= (double*)malloc(2*RPN_CHUNKSIZE*sizeof(double))))
	    Abort("%s: unable to allocate %d bytes!\n",
		  2*RPN_CHUNKSIZE*sizeof(double));
	}
	for (i=0; i<voxels_this_chunk; i++) {
	  obuf[2*i]= *(out_chunk+i-RPN_CHUNKSIZE);
	  obuf[2*i+1]= *(out_chunk+i);
	}
	mri_set_chunk(Output, chunkname, 2*voxels_this_chunk, 2*voxels_moved,
		      MRI_DOUBLE, obuf);
      }
      else {
	mri_set_chunk(Output, chunkname, voxels_this_chunk, voxels_moved,
		      MRI_DOUBLE, out_chunk);
      }
    }
    else {
      if (!rpnRun(re, voxels_this_chunk, voxels_moved))
	Abort("%s: execution error after %lld voxels: %s!\n",
	      progname,voxels_moved,rpnGetErrorString(re));
    }
    voxels_moved += voxels_this_chunk;
  }

  /* Write and close data-sets */
  for (i=0; i<n_input_files; i++) ds_close(Input+i);
  if (outfile_flg) mri_close_dataset( Output );

  rpnDestroyEngine(re);
  
  if (verbose_flg) Message( "#      Math application complete.\n" );

  return 0;
}

