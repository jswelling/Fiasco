/************************************************************
 *                                                          *
 *  phadj.c                                                 *
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
 *  Original programming by Mark Fitzgerald  2-95           *
 *     5-96: Pittsburgh Format, Mark Fitzgerald             *
 ************************************************************/
/*************************************************************

  DESCRIPTION OF PHADJ.C

  phadj.c implements a global phase change across complex
    images in order to maximize phase similarity between a
    sequence of images

  phadj.m [-input Input-header-file] [-headerout Output-header-file]
          [-dataout Output-data-file] [-parameters Parameter-file]
	  [-fixed Fixed-image-number]

**************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"

static char rcsid[] = "$Id: phadj.c,v 1.7 2003/02/07 21:32:13 welling Exp $";


int main( int argc, char* argv[] ) 
{

  MRI_Dataset *Input = NULL, *Output = NULL;
  char infile[512], hdrfile[512];
  char outfile[512], parfile[512];
  long fixed_image, temp_fixed_image;
  long dx, dy, dz, dt;
  FILE *fp = NULL;
  FComplex **img = NULL, **corr_img = NULL;
  float **fixed_ph = NULL, **adj = NULL;
  unsigned char **missing = NULL;
  long xmin, xmax, ymin, ymax, phadj_area;
  float phs, cosph, sinph;
  long x, y, t, z;
  DComplex sum;
  FComplex mean;

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
  cl_get( "headerout|h", "%option %s[%]", "phadj.mri", hdrfile );
  cl_get( "input|i", "%option %s[%]", "input.mri", infile );
  cl_get( "parameters|p", "%option %s[%]", "phadj.par", parfile );
  cl_get( "fixed|f", "%option %ld[%]", 0, &fixed_image );

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
  if( !strcmp( infile, hdrfile ) )
    Abort( "Input and output files must be distinct." );
  Input = mri_open_dataset( infile, MRI_READ );

  /* Check that program will function on data-set */
  
  if( !mri_has( Input, "images" ) )
    Abort( "%s operates only on standard images.", argv[0] );
  if( !mri_has( Input, "images.dimensions" ) || 
      strcmp( mri_get_string( Input, "images.dimensions" ), "vxyzt" ) )
    Abort( "%s operates only on images in order vxyzt.", argv[0] );
  if( !mri_has( Input, "images.extent.v" ) ||
      mri_get_int( Input, "images.extent.v" ) != 2 )
    Abort( "%s operates only on complex-valued images (v = 2).", argv[0] );
  if( !mri_has( Input, "images.file" ) ||
      !strcmp( outfile, mri_get_string( Input, "images.file" ) ) )
    Abort( "Input and output image files are the same (%s).", outfile );

  /* Set output dataset */
  Output = mri_copy_dataset( hdrfile, Input );
  hist_add_cl( Output, argc, argv );
  mri_set_string( Output, "images.datatype", "float32" );
  mri_set_string( Output, "images.file", outfile );

  /* Read/Create missing image indicators */
  missing = get_missing( Output );

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


  /* Calculate area of region over which */
  /*   phase differences are calculated  */
  xmin = (long) ( ( 3 * dx ) / 8 );
  xmax = (long) ( ( 5 * dx ) / 8 );
  ymin = (long) ( ( 3 * dy ) / 8 );
  ymax = (long) ( ( 5 * dy ) / 8 );
  phadj_area = ( xmax - xmin ) * ( ymax - ymin );


  /* Allocate image and parameter storage */
  img = (FComplex **) emalloc( dy * sizeof(FComplex *) );
  corr_img = Matrix( dy, dx, FComplex );
  fixed_ph = Matrix( (long) ( ymax - ymin ), 
		     (long) ( xmax - xmin ), float );
  adj = Matrix( dt, dz, float );
  

  /* Loop through slices */
  for( z = 0; z < dz; z++ )
    {

      temp_fixed_image = fixed_image;

      /* If image is missing find a non-missing one to */
      /*   which to calibrate                          */
      if( missing[temp_fixed_image][z] )
	{
	  temp_fixed_image++;
	  if( temp_fixed_image == dt )
	    temp_fixed_image = 0;
	  if( temp_fixed_image == fixed_image )
	    Abort( "All images for slice %ld are missing.  %s failing.",
		   z, argv[0] );
	}

      /* Get image to be considered fixed and */
      /*   calculate fixed phases             */
      *img = (FComplex *)
	mri_get_image( Input, temp_fixed_image, z, MRI_COMPLEX_FLOAT );
      realign_matrix( (void **) img, dy, 
		      (long) ( dx * sizeof(FComplex) ) );
      for( y = ymin; y < ymax; y++ )
	for( x = xmin; x < xmax; x++ )
	  fixed_ph[y-ymin][x-xmin] = Phase( img[y][x] );


      /* Loop through images */
      for( t = 0; t < dt; t++ ) 
	{

	  /* Get image and calculate mean phase difference */
	  
	  *img = (FComplex *)
	    mri_get_image( Input, t, z, MRI_COMPLEX_FLOAT );
	  realign_matrix( (void **) img, dy, 
			  (long) ( dx * sizeof(FComplex) ) );
	  sum.real = sum.imag = 0.0;
	  for( y = ymin; y < ymax; y++ )
	    for( x = xmin; x < xmax; x++ )
	      {
		phs = Phase( img[y][x] ) - fixed_ph[y-ymin][x-xmin];
		sum.real += cos( phs );
		sum.imag += sin( phs );
	      }

	  mean.real = (float) ( sum.real / (double) phadj_area );
	  mean.imag = (float) ( sum.imag / (double) phadj_area );

	  /* Calculate phase adjustment parameter */
	  /*   and cosine/sine theereof           */
	  adj[t][z] = -1.0 * Phase( mean );
	  cosph = cos( adj[t][z] );
	  sinph = sin( adj[t][z] );

	  /* Rotate phase of each complex number in the image */
	  /*   by the opposite of the mean phase difference   */
	  for( y = 0; y < dy; y++ )
	    for( x = 0; x < dx; x++ )
	      {
		corr_img[y][x].real = cosph * img[y][x].real -
		  sinph * img[y][x].imag;
		corr_img[y][x].imag = cosph * img[y][x].imag +
		  sinph * img[y][x].real;
	      }

	  /* Set corrected output image */
	  mri_set_image( Output, t, z, MRI_COMPLEX_FLOAT, *corr_img );

	}

    }


  /* Write and close data-sets */
  mri_close_dataset( Input );
  mri_close_dataset( Output );
  
  /* Write parameter file */
  fp = efopen( parfile, "w" );
  for( t = 0; t < dt; t++ )
    for( z = 0; z < dz; z++ )
      fprintf( fp, "%15.6f\n", adj[t][z] );
  efclose( fp );
  
  Message( "#      Phase adjustment complete.\n" );
  exit(0);

}

