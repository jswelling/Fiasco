/************************************************************
 *                                                          *
 *  mri_complex_to_scalar.c                                 *
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
/*************************************************************

  DESCRIPTION OF MRI_COMPLEX_TO_SCALAR

  mri_complex_to_scalar takes a complex pgh MRI dataset (of type
  vxyzt, where images.extent.v = 2) and generates a corresponding
  scalar field.  The output dataset is of type float32.

**************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"

static char rcsid[] = "$Id: mri_complex_to_scalar.c,v 1.5 2003/08/07 19:42:30 bakalj Exp $";

#define CHUNKSIZE 4096
#define DEFAULT_CHUNK_NAME "images"
#define KEYBUF_SIZE 512

static char* progname;

static int r_flg=0;
static int i_flg=0;
static int m_flg=0;
static int p_flg=0;
static int u_flg= 0;
static int n_complex_pairs= 0;
static int run_length= 0; /* length of first non-trivial dim after v */

static int verbose_flg= 0;

static void safe_copy(char* str1, char* str2) 
{
  strncpy(str1, str2, KEYBUF_SIZE);
  str1[KEYBUF_SIZE-1]= '\0';
}

static void safe_concat(char* str1, char* str2) 
{
  strncat(str1, str2, (KEYBUF_SIZE-strlen(str1))-1);
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

static void parse_input_layout(MRI_Dataset* Input, char* chunk)
{
  char* dimstr;
  char* dim_runner;
  char keybuf[KEYBUF_SIZE];
  /* Check that program will function on data-set, and calculate
   * total dataset size.
   */
  if( !mri_has( Input, chunk ) )
    Abort( "%s has no chunk \"%s\"\n", progname, chunk );
  safe_copy(keybuf,chunk);
  safe_concat(keybuf,".dimensions");
  if ( !mri_has( Input, keybuf ) )
    Abort( "%s: input dataset is missing %s.dimensions tag!\n",
	   progname,chunk);
  dimstr= strdup(mri_get_string(Input,keybuf));
  if (*dimstr != 'v')
    Abort( "%s operates only on data with first dimension v.", progname );
  safe_copy(keybuf,chunk);
  safe_concat(keybuf,".extent.v");
  if( !mri_has( Input, keybuf ) ||
      ( mri_get_int( Input, keybuf )  != 2 ) )
    Abort( "%s operates only on complex-valued data (v = 2).", progname );
  dim_runner= dimstr+1;
  n_complex_pairs= run_length= 1;
  while (*dim_runner) {
    if (run_length==1) run_length= safe_get_extent(Input,chunk,dim_runner);
    n_complex_pairs *= safe_get_extent(Input,chunk,dim_runner);
    dim_runner++;
  }
}

int main( int argc, char** argv ) 
{

  MRI_Dataset *Input = NULL, *Output = NULL;
  int i;
  char infile[512], outfile[512];
  char chunk[512];
  char keybuf[KEYBUF_SIZE];
  float* inbuf;
  float* outbuf;
  float* in_runner;
  float* out_runner;
  int n_to_go;
  int n_this_block;
  int pairs_this_run;
  long long in_offset;
  long long out_offset;
  double last_phase= 0.0;
  int unwind_cycles= 0;

  progname= argv[0];

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Deprecate old options */
  if (cl_present( "r" ))
     Abort ("Option r has been expanded to real.  Please see help file.\n");
  if (cl_present( "i" ))
     Abort ("Option i has been expanded to imaginary|ima.  Please see help file.\n");
  if (cl_present( "m" ))
     Abort ("Option m has been expanded to magnitude|mag.  Please see help file.\n");
  if (cl_present( "p" ))
     Abort ("Option p has been expanded to phase|pha.  Please see help file.\n");
  if (cl_present( "u" ))
     Abort ("Option u has been expanded to phase_unwound|phu.  Please see help file.\n");

  /* Check switches.  Make sure value is 0 or 1 for convenience in
   * checking for redundant switches.
   */
  r_flg= (cl_present("real") ? 1 : 0);
  i_flg= (cl_present("imaginary|ima") ? 1 : 0);
  m_flg= (cl_present("magnitude|mag") ? 1 : 0);
  p_flg= (cl_present("phase|pha") ? 1 : 0);
  u_flg= (cl_present("phase_unwound|phu") ? 1 : 0);
  cl_get("chunk|chu|c", "%option %s[%]", DEFAULT_CHUNK_NAME, chunk);
  verbose_flg= cl_present("v");

  /* Get filenames */
  if (!cl_get("", "%s", infile)) {
    fprintf(stderr,"%s: Input file name not given.",progname);
    Help( "usage" );
    exit(-1);
  }
  if (!cl_get("", "%s", outfile)) {
    fprintf(stderr,"%s: Output file name not given.",progname);
    Help( "usage" );
    exit(-1);
  }
  if (cl_cleanup_check()) {
    int i;
    fprintf(stderr,"%s: invalid argument in command line:\n    ",progname);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(-1);
  }

  /*** End command-line parsing ***/

  /* Print version number */
  if (verbose_flg) Message( "# %s\n", rcsid );

  if (r_flg + i_flg + m_flg + p_flg + u_flg > 1)
    Abort("%s: two exclusive flags were used in the command line.\n",progname);
  if (r_flg + i_flg + m_flg + p_flg + u_flg == 0) m_flg= 1;

  /* Open input dataset */
  if( !strcmp( infile, outfile ) )
    Abort( "Input and output files must be distinct." );
  Input = mri_open_dataset( infile, MRI_READ );

  parse_input_layout(Input, chunk);

  Output = mri_copy_dataset( outfile, Input );
  hist_add_cl( Output, argc, argv );
  safe_copy(keybuf,chunk);
  safe_concat(keybuf,".datatype");
  mri_set_string( Output, keybuf, "float32" );
  safe_copy(keybuf,chunk);
  safe_concat(keybuf,".extent.v");
  mri_set_int( Output, keybuf, 1);

  if (!(outbuf= (float*)malloc(CHUNKSIZE*sizeof(float))))
    Abort("%s: unable to allocate %d bytes!\n",progname,sizeof(float));

  n_to_go= n_complex_pairs;
  in_offset= out_offset= 0;
  pairs_this_run= 0;
  while (n_to_go) {
    n_this_block= (CHUNKSIZE>n_to_go) ? n_to_go : CHUNKSIZE;
    inbuf= mri_get_chunk(Input,chunk,2*n_this_block,in_offset,MRI_FLOAT);
    in_runner= inbuf;
    out_runner= outbuf;
    for (i=0; i<n_this_block; i++) {

      if (pairs_this_run==run_length) {
	/* new run starting */
	pairs_this_run= 0;
	last_phase= 0.0;
	unwind_cycles= 0;
      }

      if (r_flg) {
	/* copy real part to output */
	*out_runner++= *in_runner;
	in_runner += 2;
      }
      else if (i_flg) {
	/* copy imaginary part to output */
	*out_runner++= *(in_runner+1);
	in_runner += 2;
      }
      else if (m_flg) {
	/* copy magnitude to output */
	float magsqr;
	magsqr= *in_runner * *in_runner;
	in_runner++;
	magsqr += *in_runner * *in_runner;
	in_runner++;
#ifdef __GNUC__
	*out_runner++= (float)sqrt(magsqr);
#else
	*out_runner++= sqrtf(magsqr);
#endif
      }
      else if (p_flg) {
	/* copy phase to output */
#ifdef __GNUC__
	*out_runner++= (float)atan2( *(in_runner+1), *in_runner );
#else
	*out_runner++= atan2f( *(in_runner+1), *in_runner );
#endif
	in_runner += 2;
      }
      else if (u_flg) {
	double dphase;
	double phase= atan2( *(in_runner+1), *in_runner );

	phase += unwind_cycles*2.0*M_PI;
	dphase= phase-last_phase;

	if (fabs(dphase+2*M_PI)<fabs(dphase)) {
	  unwind_cycles += 1;
	  phase += 2*M_PI;
	}
	else if (fabs(dphase-2*M_PI)<fabs(dphase)) {
	  unwind_cycles -= 1;
	  phase -= 2*M_PI;
	}
	
	last_phase= phase;
	*out_runner++= (float)phase;
	in_runner += 2;
      }
      pairs_this_run++;
    }
    mri_set_chunk(Output,chunk,n_this_block,out_offset,MRI_FLOAT,outbuf);
    in_offset += 2*n_this_block;
    out_offset += n_this_block;
    n_to_go -= n_this_block;
  }

  /* Write and close data-sets */
  mri_close_dataset( Input );
  mri_close_dataset( Output );
  
  if (verbose_flg) Message( "#      Image conversion complete.\n" );

  return 0;
}

