/************************************************************
 *                                                          *
 *  mnadj.c                                                 *
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

  DESCRIPTION OF MNADJ.C

  mnadj.c is used to coerce all of the images to have
    (approximately) identical means, to adjust for 
    global drifts in signal intensity

  mnadj.m [-input Input-header-file] [-headerout Output-header-file]
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
#include "misc.h"

static char rcsid[] = "$Id: meanc.c,v 1.11 2003/12/02 23:20:41 welling Exp $";


int main( int argc, char* argv[] ) 
{

  MRI_Dataset *Input = NULL, *Output = NULL;
  char infile[512], hdrfile[512];
  char parfile[512];
  long fixed_image, temp_fixed_image;
  long dv, dx, dy, dz, dt;
  FILE *fp = NULL;
  FComplex **c_img = NULL, **c_corr_img = NULL;
  float **img = NULL, **corr_img = NULL, **adj = NULL;
  unsigned char **missing = NULL;
  long xmin, xmax, ymin, ymax, mnadj_area;
  float fixed_mean, mean;
  long x, y, t, z;
  double sum;
  int all_images_missing= 0;
  int allow_neg_means= 0;


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
  cl_get( "headerout|h", "%option %s[%]", "mnadj.mri", hdrfile );
  cl_get( "input|i", "%option %s[%]", "input.mri", infile );
  cl_get( "parameters|p", "%option %s[%]", "mnadj.par", parfile );
  cl_get( "fixed|f", "%option %ld[%]", 0, &fixed_image );
  allow_neg_means= cl_present("n|allow_negative_means");

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
  if( mri_has( Input, "images.dimensions" ) )
    {
      dv = ( !strcmp( mri_get_string( Input, "images.dimensions" ), 
		      "vxyzt" ) )?
	mri_get_int( Input, "images.extent.v" ): 
	( !strcmp( mri_get_string( Input, "images.dimensions" ), "xyzt" ) )? 
	1: 0;
      if( ( dv < 1 ) || ( dv > 2 ) )
	Abort( "%s takes only reals or complex numbers of the form (v)xyzt.",
	       argv[0] );
    }
  else
    Abort( "%s does not have the images.dimensions key.", infile );

  /* Set output dataset */
  Output = mri_copy_dataset( hdrfile, Input );
  hist_add_cl( Output, argc, argv );
  mri_set_string( Output, "images.datatype", "float32" );

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

  /* Make sure fixed image is in bounds */
  if( ( fixed_image < 0 ) || ( fixed_image >= dt ) )
    fixed_image = 0;

  /* Allocate image and parameter storage */
  if( dv == 1 )
    {
      img = (float **) emalloc( dy * sizeof(float *) );
      corr_img = Matrix( dy, dx, float );
    }
  else
    {
      c_img = (FComplex **) emalloc( dy * sizeof(FComplex *) );
      c_corr_img = Matrix( dy, dx, FComplex );
    }
  adj = Matrix( dt, dz, float );
  

  /* Calculate area of region over which mean is calculated */
  xmin = (long) ( ( 3 * dx ) / 8 );
  xmax = (long) ( ( 5 * dx ) / 8 );
  ymin = (long) ( ( 3 * dy ) / 8 );
  ymax = (long) ( ( 5 * dy ) / 8 );
  mnadj_area = ( xmax - xmin ) * ( ymax - ymin );


  /* Loop through slices */
  for( z = 0; z < dz; z++ )
    {
      all_images_missing= 0;
      temp_fixed_image = fixed_image;

    set_fixed_image:
      /* If image is missing find a non-missing one to */
      /*   which to calibrate                          */
      if( missing[temp_fixed_image][z] )
	{
	  temp_fixed_image++;
	  if( temp_fixed_image == dt )
	    temp_fixed_image = 0;
	  if ( temp_fixed_image == fixed_image ) {
	    Warning(1, "%s: All images for slice %ld are missing.\n",
		    argv[0],z );
	    all_images_missing= 1;
	  }
	}

      if (all_images_missing) {
	for( t = 0; t < dt; t++ ) 
	  {
	    
	    /* Copy image unchanged to output */
	    if( dv == 1 )
	      {
		*img = (float *)
		  mri_get_image( Input, t, z, MRI_FLOAT );
		mri_set_image( Output, t, z, MRI_FLOAT, *img );
	      }
	    else
	      {
		*c_img = (FComplex *)
		  mri_get_image( Input, t, z, MRI_COMPLEX_FLOAT );
		mri_set_image( Output, t, z, MRI_COMPLEX_FLOAT, *c_img );
	      }
	}
      }
      else {
	/* Get image to be considered fixed and calculate */
	/*   its mean (or mean modulus for complex data)  */
	sum = 0;
	if( dv == 1 )
	  {
	    *img = (float *)
	      mri_get_image( Input, temp_fixed_image, z, MRI_FLOAT );
	    realign_matrix( (void **) img, dy, 
			    (long) ( dx * sizeof(float) ) );
	    for( y = ymin; y < ymax; y++ )
	      for( x = xmin; x < xmax; x++ )
		sum += img[y][x];
	  }
	else
	  {
	    *c_img = (FComplex *)
	      mri_get_image( Input, temp_fixed_image, z, MRI_COMPLEX_FLOAT );
	    realign_matrix( (void **) c_img, dy, 
			    (long) ( dx * sizeof(FComplex) ) );
	    for( y = ymin; y < ymax; y++ )
	      for( x = xmin; x < xmax; x++ )
		sum += Modulus( c_img[y][x] );
	  }
	
	/* If sum of image (or modulus of image) is non-positive,      */
	/*   data has been corrupted --- declare missing and try again */
	if( sum == 0.0 || (sum < 0.0 && !allow_neg_means) )
	  {
	    Warning( 1,
		     "Image %ld for slice %ld has non-positive mean.  Declared missing.\n",
		     temp_fixed_image, z );
	    missing[temp_fixed_image][z] = (unsigned char) 1;
	    goto set_fixed_image;
	  }
	
	/* Set mean to which to calibrate for this slice */
	fixed_mean = (float) ( sum / (double) mnadj_area );
	
	
	/* Loop through images */
	for( t = 0; t < dt; t++ ) 
	  {
	    adj[t][z]= 1.0;
	    
	    /* Get image and calculate its mean */
	    sum = 0;
	    if( dv == 1 )
	      {
		*img = (float *)
		  mri_get_image( Input, t, z, MRI_FLOAT );
		realign_matrix( (void **) img, dy, 
				(long) ( dx * sizeof(float) ) );
		for( y = ymin; y < ymax; y++ )
		  for( x = xmin; x < xmax; x++ )
		    sum += img[y][x];
	      }
	    else
	      {
		*c_img = (FComplex *)
		  mri_get_image( Input, t, z, MRI_COMPLEX_FLOAT );
		realign_matrix( (void **) c_img, dy, 
				(long) ( dx * sizeof(FComplex) ) );
		for( y = ymin; y < ymax; y++ )
		  for( x = xmin; x < xmax; x++ )
		    sum += Modulus( c_img[y][x] );
	      }
	    mean = (float) ( sum / (double) mnadj_area );
	    
	    /* If sum of image (or modulus of image) is non-positive */
	    /*   data has been corrupted --- declare missing         */
	    if( sum == 0.0 || ( sum < 0.0 && !allow_neg_means) )
	      {
		if (!missing[t][z])
		  Warning( 1, 
			 "Image %ld for slice %ld has non-positive mean.  Declared missing.\n",
			   t, z );
		missing[t][z] = (unsigned char) 1;
		
		/* Set output image to zeroes */
		if( dv == 1 )
		  for( y = 0; y < dy; y++ )
		    for( x = 0; x < dx; x++ )
		      corr_img[y][x] = 0.0;
		else
		  for( y = 0; y < dy; y++ )
		    for( x = 0; x < dx; x++ )
		      c_corr_img[y][x].real = c_corr_img[y][x].imag = 0.0;
		
	      }
	    else 
	      {
		/* Calculate multiplicative adjustment parameter */
		adj[t][z] = fixed_mean / mean;
		
		/* Multiply adjustment across entire image */
		if( dv == 1 )
		  for( y = 0; y < dy; y++ )
		    for( x = 0; x < dx; x++ )
		      corr_img[y][x] = img[y][x] * adj[t][z];
		else
		  for( y = 0; y < dy; y++ )
		    for( x = 0; x < dx; x++ )
		      {
			c_corr_img[y][x].real = c_img[y][x].real * adj[t][z];
			c_corr_img[y][x].imag = c_img[y][x].imag * adj[t][z];
		      }
		
		/* Set corrected output image */
		if( dv == 1 )
		  mri_set_image( Output, t, z, MRI_FLOAT, *corr_img );
		else
		  mri_set_image( Output, t, z, MRI_COMPLEX_FLOAT, *c_corr_img );
	      }
	}
      }

    }

  /* Write out missing chunk, since it may have changed */
  mri_set_chunk( Output, "missing", (long) ( dt * dz ), 0,
		 MRI_UNSIGNED_CHAR, *missing );

  /* Write and close data-sets */
  mri_close_dataset( Input );
  mri_close_dataset( Output );
  
  /* Write parameter file */
  fp = efopen( parfile, "w" );
  fprintf(fp,"##Format: order:z_fastest type:raw names:(meanc)\n");
  for( t = 0; t < dt; t++ )
    for( z = 0; z < dz; z++ )
      fprintf( fp, "%15.6f\n", adj[t][z] );
  efclose( fp );
  
  Message( "#      Mean adjustment complete.\n" );
  exit(0);

}

