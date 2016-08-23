/************************************************************
 *                                                          *
 *  mri_splus_filter.c                                     *
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

  DESCRIPTION OF MRI_SPLUS_FILTER

  mri_splus_filter takes a pgh MRI dataset (of type txyz, or
  vtxyz if images.extent.v = 1) and filters it through Splus.
  Under the direction of Splus, multiple xyz datasets are
  constructed.  The output dataset is of type float32.

  mri_splus_filter [-splus splus_exe] [-dlo dlo.o] [-init init.S] 
                   -script splus_script.S infile

    -splus splus_exe   gives the full path of the Splus executable
                       (default "/usr/statlocal/bin/Splus" )
    -dlo dlo.o         gives the path of the dynamically loaded 
                       routines (default  "splus_binary_pipes.o").
    -init init.S       gives the path of the Splus routines to
                       access the dynamically loaded routines
                       (default "bio_init.S")
    -script script.S   gives the Splus script that defines
                       the mapping between input and output

    infile             is the input dataset

  The Splus script must be constructed to take the input data in
  voxel order (read from a binary pipe by the slave_splus mechanism).
  It must write only the following to Splus's stdout: an integer
  giving the number of output datasets, followed on separate lines
  by the names to be given to those datasets.  After each name is
  given, it must write the data for that file to the binary output
  pipe (as per the slave_splus mechanism).

**************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"
#include "slave_splus.h"

static char rcsid[] = "$Id: mri_splus_filter.c,v 1.7 2003/08/07 20:52:00 bakalj Exp $";

#define DEFAULT_SPLUS_EXE "/usr/statlocal/bin/Splus"
#define DEFAULT_DLO_PATH "splus_binary_pipes.o"
#define DEFAULT_INIT_PATH "bio_init.S"

#define REPLY_LENGTH 512

int main( int argc, char* argv[] ) 
{
  char *fake_argv[2];
  SlaveSplus* sp= NULL;
  MRI_Dataset *Input = NULL, *Output = NULL;
  int r_flg=0;
  int i_flg=0;
  int m_flg=0;
  int p_flg=0;
  int i;
  int j;
  int k;
  int t;
  char scriptname[512];
  char splusname[512];
  char dloname[512];
  char initname[512];
  char infile[512], outfile[512];
  char reply[REPLY_LENGTH];
  float* chunk_out= NULL;
  float* chunk_in= NULL;
  int dx;
  int dy;
  int dz;
  int dt;
  char* dimstr;
  long total_voxels;
  long voxels_moved;

  /* Print version number */
  Message( "# %s\n", rcsid );

  /* Check to see if help was requested */
  if( ( argc > 1 ) && !strcmp( argv[1], "-help" ) )
    {
      if( argc == 2 )
	Help( "selecttopic" );
      else
	Help( argv[2] );
    }

  Abort("%s: mri_splus_filter has been deprecated. \n", argv[0]);

  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Get filenames */
  if (!cl_get("script","%option %s",scriptname)) {
    fprintf(stderr,"%s: Splus script name not given.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  cl_get("splus","%option %s[%]",DEFAULT_SPLUS_EXE,splusname);
  cl_get("dlo","%option %s[%]",DEFAULT_DLO_PATH,dloname);
  cl_get("init","%option %s[%]",DEFAULT_INIT_PATH,initname);
  if (!cl_get("", "%s", infile)) {
    fprintf(stderr,"%s: Input file name not given.",argv[0]);
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
  
  /* Open input dataset */
  if( !strcmp( infile, outfile ) )
    Abort( "Input and output files must be distinct." );
  Input = mri_open_dataset( infile, MRI_READ );

  /* Check that program will function on data-set */
  if( !mri_has( Input, "images" ) )
    Abort( "%s operates only on standard images.", argv[0] );
  if ( !mri_has( Input, "images.dimensions" ) )
    Abort( "%s requires images.dimensions.", argv[0]);
  dimstr= mri_get_string( Input, "images.dimensions" );
  if ( !strcmp(dimstr,"vtxyz") ) {
    if (!mri_has(Input,"images.extent.v") 
	|| (mri_get_int(Input,"images.extent.v") != 1))
      Abort( "%s requires that vtxyz data have vector length 1.",argv[0]);
  }

  /* Set parameters in local variables */
  if( !mri_has( Input, "images.datatype" ) )
    Abort( "%s: Header missing images.datatype key.", argv[0] );
  if( !mri_has( Input, "images.extent.t" ) ||
      !mri_has( Input, "images.extent.x" ) ||
      !mri_has( Input, "images.extent.y" ) ||
      !mri_has( Input, "images.extent.z" ) )
    Abort( "%s: images.extent key(s) missing from header.", argv[0] );
  dt = mri_get_int( Input, "images.extent.t" );
  dx = mri_get_int( Input, "images.extent.x" );
  dy = mri_get_int( Input, "images.extent.y" );
  dz = mri_get_int( Input, "images.extent.z" );
  if( ( dt <= 0 ) || ( dx <= 0 ) || ( dy <= 0 ) || ( dz <= 0 ) )
    Abort( "%s: images.extent key(s) is non-positive.",argv[0] );

  /* Spawn the Splus child */
  fake_argv[0]= splusname;
  fake_argv[1]= NULL;
  if (!(sp= spawn_slave_splus(splusname,dloname,initname))) {
    Abort("%s: spawn failed!\n",argv[0]);
  }

  /* Have the child read our script */
  fprintf(sp->in,"source(\"%s\")\n",scriptname);
  fflush(sp->in);

  /* Send the data */
  total_voxels= dt*dx*dy*dz;
  voxels_moved= 0;
  for (k=0; k<dz; k++) {
    for (j=0; j<dy; j++) {
      chunk_in= mri_get_chunk(Input, "images", dt*dx, voxels_moved, 
			      MRI_FLOAT);
      if( write(sp->din, (void*)chunk_in, dt*dx*sizeof(float)) == -1 ) {
	perror("binary write error");
	Abort("%s: error writing binary data at row %d, slice %d",
	       argv[0],j,k);
      }
      voxels_moved += dt*dx;
    }
  }

  /* Get the number of output files back */
  if (!fgets(reply,REPLY_LENGTH,sp->out)) {
    Error("%s: Failed to get output dataset count from splus!");
  }
  fprintf(stderr,"reply was <%s>\n",reply);

#ifdef never
  Output = mri_copy_dataset( outfile, Input );
  hist_add_cl( Output, argc, argv );
  mri_set_string( Output, "images.datatype", "float32" );
  mri_set_string( Output, "images.dimensions", "xyzt" );

  image_out= (float*)malloc(dx*dy*sizeof(float));

  for (t=0; t<dt; t++) {
    for (k=0; k<dz; k++) {
      image_in= mri_get_image(Input, t, k, MRI_COMPLEX_FLOAT);
      in_runner= image_in;
      out_runner= image_out;
      if (r_flg) {
	/* copy real part to output */
	for (j=0; j<dy; j++)
	  for (i=0; i<dx; i++) {
	    *out_runner++= *in_runner;
	    in_runner += 2;
	  }
      }
      else if (i_flg) {
	/* copy imaginary part to output */
	in_runner++; /* shift to imaginary components */
	for (j=0; j<dy; j++)
	  for (i=0; i<dx; i++) {
	    *out_runner++= *in_runner;
	    in_runner += 2;
	  }
      }
      else if (m_flg) {
	/* copy magnitude to output */
	float magsqr;
	for (j=0; j<dy; j++)
	  for (i=0; i<dx; i++) {
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
      }
      else if (p_flg) {
	/* copy phase to output */
	for (j=0; j<dy; j++)
	  for (i=0; i<dx; i++) {
#ifdef __GNUC__
	    *out_runner++= (float)atan2( *(in_runner+1), *in_runner );
#else
	    *out_runner++= atan2f( *(in_runner+1), *in_runner );
#endif
	    in_runner += 2;
	  }
      }
      mri_set_image( Output, t, k, MRI_FLOAT, image_out );
    }
  }

  /* Write and close data-sets */
  mri_close_dataset( Input );
  mri_close_dataset( Output );
#endif

  /* Kill the child splus */
  if (sp) kill_slave_splus(sp);
  
  Message( "#      Data filtering complete.\n" );
  return 0;

}

