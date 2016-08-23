/************************************************************
 *                                                          *
 *  detrend.c                                               *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1995 Department of Statistics,         *
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
 *  Original programming by Chris Genovese  5-95            *
 *     5-96: Pittsburgh Format, Mark Fitzgerald             *
 ************************************************************/
/*************************************************************

  DESCRIPTION OF DETREND.C

  detrend is used to remove linear temporal trends from a sequence of
  images


**************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"

static char rcsid[] = "$Id: detrend.c,v 1.11 2003/09/25 20:27:59 welling Exp $";


int main( int argc, char* argv[] ) 
{

  MRI_Dataset *Input = NULL, *Output = NULL, *Parput = NULL;
  char infile[512], hdrfile[512], outfile[512];
  char phfile[512];
  char* pdfile;
  long dv, dx, dy, dz, dt;
  float *ts = NULL, *times = NULL, *corr_ts = NULL, halft;
  float **beta0 = NULL, **sebeta0 = NULL, **beta1 = NULL, **sebeta1 = NULL;
  long x, y, t, z;
  long offset;
  double sumxx, sumy, sumxy, stdvres, sqrtdt;

  /* Print version number */
  Message( "# %s\n", rcsid );

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Deprecate old options */

   if (cl_present( "dataout|d" ))
     Abort ("Option dataout|d has been replaced by infile outfile format.  Please see help file.\n");
   if (cl_present( "headerout|h" ))
     Abort ("Option headerout|h has been replaced by infile outfile format.  Please see help file.\n");
   if (cl_present( "input|i" ))
     Abort ("Option input|i has been replaced by infile outfile format.  Please see help file.\n");
   if (cl_present( "pardata|pard" ))       
     Abort ("Option pardata|pard has been replaced by estimates|est|e.  Please see help file.\n");
   if (cl_present( "parhead|parh" ))       
     Abort ("Option parhead|parh has been replaced by estimates|est|e.  Please see help file.\n");


  /* Get filenames */

  cl_get( "estimates|est|e", "%option %s[%]", "detpar", phfile );

  if(!cl_get("", "%s", infile)) {
    fprintf(stderr, "%s: Input file name not given.\n", argv[0]);
    exit(-1);
  }
  if(!cl_get("", "%s", hdrfile)) {
    fprintf(stderr, "%s: Output file name not given.\n", argv[0]);
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

  pdfile = ".dat";
  strcpy(outfile,".dat");

  /* Open input dataset */
  if( !strcmp( infile, hdrfile ) || !strcmp( infile, phfile ) ||
      !strcmp( hdrfile, phfile ) )
    Error( "Input, output, and parameter header files must be distinct (%s %s %s.",
	   infile, hdrfile, phfile );
  Input = mri_open_dataset( infile, MRI_READ );

  /* Check that program will function on data-set */
  if( !mri_has( Input, "images" ) )
    Abort( "%s operates only on standard images.", argv[0] );
  if( !mri_has( Input, "images.dimensions" ) || 
      ( strcmp( mri_get_string( Input, "images.dimensions" ), "vtxyz" ) &&
	strcmp( mri_get_string( Input, "images.dimensions" ), "txyz" ) ) )
    Abort( "%s operates only on images in order (v)txyz.", argv[0] );
  if ( !strcmp(mri_get_string( Input, "images.dimensions" ), "vtxyz") ) {
    if( !mri_has( Input, "images.extent.v" ) ||
	( ( dv = mri_get_int( Input, "images.extent.v" ) ) != 1 ) ) {
      Abort( "%s operates only on real-valued images (v = 1).", argv[0] );
    }
  }
  else dv= 1;
  /* Check laboriously for output data file name conflicts */
  {
    int check_input= 0;
    char* in_fname= NULL;
    int check_output= 0;
    int check_pdfile= 0;
    if ( mri_has(Input, "images.file" ) ) {
      in_fname= mri_get_string( Input, "images.file" );
      if (*in_fname != '.') check_input= 1;
    }
    if (*outfile != '.') check_output= 1;
    if (*pdfile != '.') check_pdfile= 1;
    if ((check_input && check_output && !strcmp(in_fname,outfile))
	|| (check_input && check_pdfile && !strcmp(in_fname,pdfile))
	|| (check_output && check_pdfile && !strcmp(outfile,pdfile)))
      Abort( "Input, output, and/or parameter data files must be distinct (%s %s %s).",
	     mri_get_string( Input, "images.file" ), outfile, pdfile );
  }

  /* Set parameters in local variables */
  if( !mri_has( Input, "images.extent.t" ) ||
      !mri_has( Input, "images.extent.x" ) ||
      !mri_has( Input, "images.extent.y" ) ||
      !mri_has( Input, "images.extent.z" ) )
    Abort( "images.extent key(s) missing from header." );
  dt = mri_get_int( Input, "images.extent.t" );
  dx = mri_get_int( Input, "images.extent.x" );
  dy = mri_get_int( Input, "images.extent.y" );
  dz = mri_get_int( Input, "images.extent.z" );
  if( ( dt <= 0 ) || ( dx <= 0 ) || ( dy <= 0 ) || ( dz <= 0 ) )
    Abort( "images.extent key(s) is non-positive." );

  /* Make sure time sequence is long enough to de-trend */
  if( dt < 3 )
    Abort( "Time series is too short (less than 3) to detrend." );

  /* Set output dataset */
  Output = mri_copy_dataset( hdrfile, Input );
  hist_add_cl( Output, argc, argv );
  mri_set_string( Output, "images.file", outfile );
  mri_set_string( Output, "images.datatype", "float32" );

  /* Set parameter dataset */
  Parput = mri_copy_dataset( phfile, Input );
  mri_set_string( Parput, "images.file", ".dat" );
  mri_set_string( Parput, "images.datatype", "float32" );
  mri_set_string( Parput, "images.dimensions", "vxyzt" );
  mri_set_int( Parput, "images.extent.t", 4 );

  /* Allocate image and parameter storage */
  ts = (float *) emalloc( dt * sizeof(float) );
  times = (float *) emalloc( dt * sizeof(float) );
  corr_ts = (float *) emalloc( dt * sizeof(float) );
  beta0 = Matrix( dy, dx, float );
  sebeta0 = Matrix( dy, dx, float );
  beta1 = Matrix( dy, dx, float );
  sebeta1 = Matrix( dy, dx, float );

  /* DETREND */

  /* Establish time sequence and pre-calculate fixed sum */
  halft = (float) ( dt - 1 ) / 2.0;
  sumxx = 0.0;
  for( t = 0; t < dt; t++ )
    {
      times[t] = (float) t - halft;
      sumxx += ( times[t] * times[t] );
    }
  sqrtdt = sqrt( (double) dt );

  /* Loop through slices */
  for( z = 0; z < dz; z++ )
    {
      /* Loop through pixels in image */
      for( y = 0; y < dy; y++ )
	for( x = 0; x < dx; x++ )
	  {
	    /* Read a pixel's time series */
	    offset = dt * ( dx * ( dy * z + y ) + x );
	    ts = (float*) mri_get_chunk( Input, "images", dt, 
					 offset, MRI_FLOAT );
	    
	    /* Calculate sums for regression */
	    sumy = sumxy = 0.0;
	    for( t = 0; t < dt; t++ )
	      {
		sumy += ts[t];
		sumxy += ( times[t] * ts[t] );
	      }

	    /* Calculate regression coefficients */
	    beta0[y][x] = sumy / (float) dt;
	    beta1[y][x] = sumxy / sumxx;

	    /* Calculate detrended time series and estimate standard error */
	    stdvres = 0.0;
	    for( t = 0; t < dt; t++ )
	      {
		corr_ts[t] = ts[t] - beta1[y][x] * times[t];
		stdvres += pow( ( corr_ts[t] - beta0[y][x] ), 2.0 );
	      }
	    mri_set_chunk( Output, "images", dt, offset, MRI_FLOAT, corr_ts );
	    stdvres /= (double) ( dt - 2 );
	    
	    /* Calculate standard errors for coefficients */
	    sebeta0[y][x] = stdvres / sqrtdt;
	    sebeta1[y][x] = stdvres / sumxx;

	  }

      /* Write out a slice's worth of regression coefficients */
      /*   and standard error thereof                         */
      mri_set_image( Parput, 0, z, MRI_FLOAT, *beta0 );
      mri_set_image( Parput, 1, z, MRI_FLOAT, *sebeta0 );
      mri_set_image( Parput, 2, z, MRI_FLOAT, *beta1 );
      mri_set_image( Parput, 3, z, MRI_FLOAT, *sebeta1 );
      
    }
	    
  
  /* Write and close data-sets */
  mri_close_dataset( Input );
  mri_close_dataset( Output );
  mri_close_dataset( Parput );
  
  Message( "#      De-trending complete.\n" );
  exit(0);

}

