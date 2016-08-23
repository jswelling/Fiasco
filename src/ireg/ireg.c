/************************************************************
 *                                                          *
 *  ireg.c                                                  *
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

  DESCRIPTION OF IREG.C

  ireg.c is used to perform motion correction on a set
    of images, given estimates of the motion (as those
    produced by estireg.c)

  ireg.m [-input Input-header-file] [-headerout Output-header-file]
           [-dataout Output-data-file]
	   [-parameters Registration-parameter-file]

**************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"

static char rcsid[] = "$Id: ireg.c,v 1.10 2003/02/07 21:32:12 welling Exp $";


int main( int argc, char* argv[] ) 
{

  MRI_Dataset *Input = NULL, *Output = NULL;
  char infile[512], hdrfile[512], outfile[512], parfile[512];
  long dv, dx, dy, dz, dt;
  FILE *fp = NULL;
  FComplex **c_image = NULL, **c_corr_image = NULL;
  float **image = NULL, **corr_image = NULL;
  long x, y, t, z;
  char scanline[512];
  RegPars **par = NULL;
  long linenum, numread, tt, zz;
  RegPars tmp;
  char which_space;

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
  cl_get( "headerout|h", "%option %s[%]", "ireg.mri", hdrfile );
  cl_get( "input|i", "%option %s[%]", "input.mri", infile );
  cl_get( "parameters|p", "%option %s[%]", "reg.par", parfile );

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
  mri_set_string( Output, "images.file", outfile );
  mri_set_string( Output, "images.dimensions", "vxyzt" );
  mri_set_int( Output, "images.extent.v", (int) dv );

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

  /* Allocate image and parameter storage */
  c_image = Matrix( dy, dx, FComplex );
  c_corr_image = Matrix( dy, dx, FComplex );
  if( dv == 1 )
    {
      image = (float **) emalloc( dy * sizeof(float *) );
      corr_image = Matrix( dy, dx, float );
    }
  par = Matrix( dt, dz, RegPars );
  

  /* Initialize registration parameters */
  /*   to zero and unset */
  for( t = 0; t < dt; t++ )
    for( z = 0; z < dz; z++ ) {
      par[t][z].x_shift = par[t][z].y_shift =
	par[t][z].rotation = 0.0;
      par[t][z].valid= 0;
    }

  /* Read in estimated registration parameters */
  fp = efopen( parfile, "r" );
  linenum = -1;
  while( !feof( fp ) && !ferror( fp ) )
    {

      /* Scan a line */
      do {
	if (feof(fp)) break;
	linenum++;
	fscanf( fp, "%510[^\n]%*[\n]", scanline );
      } while (scanline[0]=='#'); /* ignore comments */
      numread = sscanf( scanline, "%ld%ld%f%f%f%*f", &tt, &zz, 
			&(tmp.x_shift), &(tmp.y_shift), &(tmp.rotation) );
      if (feof(fp)&&scanline[0]=='#') continue;
      if( numread < 5 )
	{
	  Warning( 1, "Line %ld of %s is too short (%ld) -- Ignoring.\n",
		   linenum, parfile, numread );
	  Error( "Line %ld of %s is too short (%ld) -- Ignoring.\n",
		   linenum, parfile, numread );
	  continue;
	}
      
      /* Check to see if image and slice numbers are in bounds */
      if( ( tt < 0 ) || ( tt >= dt ) || ( zz < 0 ) || ( zz >= dz ) )
	{
	  Warning( 1,
		   "Image/Slice number out of bounds (%ld/%ld) for line %d of %s -- Ignoring.\n",
		   tt, zz, linenum, parfile );
	  continue;
	}
      
      /* Put parameters into appropriate storage */
      par[tt][zz].x_shift = tmp.x_shift;
      par[tt][zz].y_shift = tmp.y_shift;
      par[tt][zz].rotation = tmp.rotation;
      par[tt][zz].valid = 1;
    }
  efclose( fp );

  /* Check to see that we have values for all slices */
  for( t = 0; t < dt; t++ )
    for( z = 0; z < dz; z++ ) {
      if (!par[t][z].valid)
	Abort("No parameters for slice %d image %d!\n",z,t);
    }

  /* If real-valued input, set imaginary part */
  /*   of complex image to zero               */
  if( dv == 1 )
    {
      which_space = 'i';
      for( y = 0; y < dy; y++ )
	for( x = 0; x < dx; x++ )
	  c_image[y][x].imag = 0.0;
    }
  else
    which_space ='k';


  /* Loop through images */
  for( t = 0; t < dt; t++ )
    {
      for( z = 0; z < dz; z++ ) 
	{
	  
	  /* Get image and re-align matrix pointers */
	  if( dv == 1 )
	    {
	      *image = (float *) 
		mri_get_image( Input, t, z, MRI_FLOAT );
	      realign_matrix( (void **) image, dy, 
			      (long) ( dx * sizeof(float) ) );
	      
	      /* Convert to complex-valued for registration */
	      for( y = 0; y < dy; y++ )
		for( x = 0; x < dx; x++ )
		  c_image[y][x].real = image[y][x];
	      
	    }
	  else
	    {
	      *c_image = (FComplex *)
		mri_get_image( Input, t, z, MRI_COMPLEX_FLOAT );
	      realign_matrix( (void **) c_image, dy,
			      (long) ( dx * sizeof(FComplex) ) );
	    }
	  
	  /* Perform registration */
	  fourier_shift_rot( par[t][z], c_image, c_corr_image,
			     dx, dy, which_space );
	  
	  /* Write out registered image */
	  if( dv == 1 )
	    {
	      /* If i-space, take modulus */
	      for( y = 0; y < dy; y++ )
		for( x = 0; x < dx; x++ )
		  {
		    corr_image[y][x] = Modulus( c_corr_image[y][x] );
		  }
	      mri_set_image( Output, t, z, MRI_FLOAT, *corr_image );
	    }
	  else
	    {
	      mri_set_image( Output, t, z, MRI_COMPLEX_FLOAT, *c_corr_image );
	    }

	}

      /* Print progress report after completing a set of slices */
      if( !t )
	Message( "      " );
      if( ( t == ( dt - 1 ) ) || !( ( t + 1 ) % 60 ) )
	{
	  Message( "# %ld\n", (long) ( t + 1 ) );
	  if( t != ( dt - 1 ) )
	    Message( "      " );
	}
      else
	Message( "#" );
      
    }
      
  /* Write and close data-sets */
  mri_close_dataset( Input );
  mri_close_dataset( Output );
  
  Message( "#      Image registration complete.\n" );
  exit(0);

}

