/************************************************************
 *                                                          *
 *  smooth_missing.c                                        *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1999 Department of Statistics,         *
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
 *  Original programming by Joel Welling 3/1999             *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"

static char rcsid[] = "$Id: smooth_missing.c,v 1.4 2003/02/07 21:28:33 welling Exp $";


static void smooth_ts( float* corr_ts, float* ts, unsigned char** missing,
		       int z, int dt )
{
  int t;
  int t_start;
  int t_end;
  int loop;
  int in_run= 0;

  for (t=0; t<dt; t++) {
    if (missing[t][z]) {
      if (!in_run) { /* start of run */
	in_run= 1;
	t_start= t;
      }
    }
    else { /* not a missing datum */
      if (in_run) { /* fill this gap */
	t_end= t;
	if (t_start==0) { /* gap starts at t=0 */
	  for (loop=0; loop<t_end; loop++) corr_ts[loop]= ts[t];
	}
	else {
	  float slope= (ts[t_end]-ts[t_start-1])/(t_end + 1 - t_start);
	  for (loop=t_start; loop<t_end; loop++)
	    corr_ts[loop]= ts[t_start-1] + ((loop-t_start)+1)*slope;
	}
	in_run= 0;
      }
      corr_ts[t]= ts[t]; /* this datum is valid */
    }
  }

  /* Finish final run */
  if (in_run) {
    t_end= dt-1;
    if (t_start==0) {
      /* all missing; fill with zero */
      for (loop=0; loop<dt; loop++) corr_ts[loop]= 0.0;
    }
    else {
      /* we have no valid last datapoint */
      for (loop= t_start; loop<dt; loop++) corr_ts[loop]= ts[t_start-1];
    }
  }
}

int main( int argc, char* argv[] ) 
{
  MRI_Dataset *Input = NULL, *Output = NULL;
  char infile[512], hdrfile[512], outfile[512];
  long dv, dx, dy, dz, dt;
  float *ts = NULL;
  float* corr_ts= NULL;
  long x, y, t, z;
  unsigned char** missing= NULL;
  int offset;

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


  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Get filenames */
  cl_get( "dataout|d", "%option %s[%]", ".dat", outfile );
  cl_get( "headerout|h", "%option %s[%]", "smooth_missing.mri", hdrfile );
  cl_get( "input|i", "%option %s[%]", "input.mri", infile );

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
  if( !strcmp( infile, hdrfile ) ) {
    Error( "%s: Input and output files must be distinct,\n",argv[0]);
    Help("usage");
  }
    
  Input = mri_open_dataset( infile, MRI_READ );

  /* Check that program will function on data-set */
  if( !mri_has( Input, "images" ) ) {
    Error( "%s operates only on standard images.", argv[0] );
    Help("usage");
  }
  if( !mri_has( Input, "images.dimensions" ) || 
      ( strcmp( mri_get_string( Input, "images.dimensions" ), "vtxyz" ) &&
	strcmp( mri_get_string( Input, "images.dimensions" ), "txyz" ) ) ) {
    Error( "%s operates only on images in order (v)txyz.\n", argv[0] );
    Help("usage");
  }
  if ( !strcmp(mri_get_string( Input, "images.dimensions" ), "vtxyz") ) {
    if( !mri_has( Input, "images.extent.v" ) ||
	( ( dv = mri_get_int( Input, "images.extent.v" ) ) != 1 ) ) {
      Error( "%s operates only on real-valued images (v = 1).\n", argv[0] );
      Help("usage");
    }
  }
  else dv= 1;

  /* Set parameters in local variables */
  if( !mri_has( Input, "images.extent.t" ) ||
      !mri_has( Input, "images.extent.x" ) ||
      !mri_has( Input, "images.extent.y" ) ||
      !mri_has( Input, "images.extent.z" ) ) 
    Abort( "images.extent key(s) missing from header.\n" );
  dt = mri_get_int( Input, "images.extent.t" );
  dx = mri_get_int( Input, "images.extent.x" );
  dy = mri_get_int( Input, "images.extent.y" );
  dz = mri_get_int( Input, "images.extent.z" );
  if( ( dt <= 0 ) || ( dx <= 0 ) || ( dy <= 0 ) || ( dz <= 0 ) )
    Abort( "images.extent key(s) is non-positive.\n" );

  /* Make sure time sequence is long enough to de-trend */
  if( dt < 3 )
    Abort( "Time series is too short (less than 3) to smooth_missing.\n" );

  /* Set output dataset */
  Output = mri_copy_dataset( hdrfile, Input );
  hist_add_cl( Output, argc, argv );
  mri_set_string( Output, "images.file", outfile );
  mri_set_string( Output, "images.datatype", "float32" );

  /* We'll need the missing info */
  missing= get_missing(Input);

  /* Allocate image and parameter storage */
  corr_ts = (float *) emalloc( dt * sizeof(float) );

  /* Loop through the time series. Slices very most slowly, which
   * certainly makes life easier. 
   */
  for( z = 0; z < dz; z++ ) {
    /* Loop through pixels in image */
    for( y = 0; y < dy; y++ )
      for( x = 0; x < dx; x++ )
	{
	  /* Read a pixel's time series */
	  offset = dt * ( dx * ( dy * z + y ) + x );
	  ts = (float*) mri_get_chunk( Input, "images", dt, 
				       offset, MRI_FLOAT );
	  
	  smooth_ts( corr_ts, ts, missing, z, dt );

	  mri_set_chunk( Output, "images", dt, offset, MRI_FLOAT, corr_ts );
	  
	}
  }
  
  /* Write and close data-sets */
  mri_close_dataset( Input );
  mri_close_dataset( Output );

  /* Clean up */
  free(corr_ts);
  
  Message( "#      Missing value smoothing complete.\n" );
  exit(0);

}

