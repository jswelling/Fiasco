/************************************************************
 *                                                          *
 *  outlier.c                                               *
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

  DESCRIPTION OF OUTLIER.C

  outlier.c is used to identify and "pull in" outliers

  outlier.m [-input Input-header-file] [-headerout Output-header-file]
             [-dataout Output-data-file] [-parameters Parameter-file]
             [-cutoff Stdvs] [-badimage Proportion]

**************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"

static char rcsid[] = "$Id: outlier.c,v 1.10 2004/01/16 19:18:18 welling Exp $";

static char* progname= NULL;

#define KEYBUF_SIZE 512

#define IQR_TO_STDV (1.0/1.35)

static void safe_copy(char* str1, char* str2) {
  strncpy(str1, str2, KEYBUF_SIZE);
  str1[KEYBUF_SIZE-1]= '\0';
}

static void safe_concat(char* str1, char* str2) {
  strncat(str1, str2, (KEYBUF_SIZE-strlen(str1))-1);
}

static int safe_get_extent(MRI_Dataset* ds, char* chunk, char* dim)
{
  char key_buf[KEYBUF_SIZE];
  char dim_buf[4];
  dim_buf[0]= *dim;
  dim_buf[1]= '\0';
  safe_copy(key_buf,chunk);
  safe_concat(key_buf,".extent.");
  safe_concat(key_buf,dim_buf);
  if (mri_has(ds,key_buf)) return mri_get_int(ds,key_buf);
  else Abort("%s: input file missing tag %s!\n",progname,key_buf);
  return 0; /* not reached */
}

static void check_sizes(const char* fname, MRI_Dataset* ds, 
			const char* field,
			const int dv, const int dx, const int dy, 
			const int dz, const int dt)
{
    if( mri_has( ds, "images.dimensions" ) )
      {
	int mydv;
	int mydt;
	char* dimstr= mri_get_string( ds, "images.dimensions" );
	if (!strcmp(dimstr,"vxyzt")) {
	  mydv= safe_get_extent(ds,"images","v");
	  mydt= safe_get_extent(ds,"images","t");
	}
	else if (!strcmp(dimstr,"xyzt")) {
	  mydv= 1;
	  mydt= safe_get_extent(ds,"images","t");
	}
	else if (!strcmp(dimstr,"vxyz")) {
	  mydv= safe_get_extent(ds,"images","v");
	  mydt= 1;
	}
	else if (!strcmp(dimstr,"xyz")) {
	  mydt= 1;
	  mydv= 1;
	}
	else {
	  Abort("%s: the %s dataset must be (v)xyz(t)\n",progname,field);
	}
	if (mydv != dv) 
	  Abort("%s: %s file v dim should be %d!\n",progname,field,dv);
	if (mydv != dt) 
	  Abort("%s: %s file t dim should be %d!\n",progname,field,dt);
	if (safe_get_extent(ds,"images","x") != dx) 
	  Abort("%s: %s file x dim should be %d!\n",progname,field,dx);
	if (safe_get_extent(ds,"images","y") != dy) 
	  Abort("%s: %s file y dim should be %d!\n",progname,field,dy);
	if (safe_get_extent(ds,"images","z") != dz) 
	  Abort("%s: %s file z dim should be %d!\n",progname,field,dz);
      }
    else
      Abort( "%s does not have the images.dimensions key.", fname );

}

static void calc_mean_stdv( MRI_Dataset* Input, unsigned char** missing,
			    double **mean, double **stdv, float** img,
			    int z, int dx, int dy, int dz, int dt,
			    long* countOut )
{
  long count, x, y, t;

  /* Initialize sums to zero */
  count = 0;
  for( y = 0; y < dy; y++ )
    for( x = 0; x < dx; x++ )
      mean[y][x] = stdv[y][x] = 0.0;
  
  /* Loop through images and sum up */
  for( t = 0; t < dt; t++ )
    {
      if( !missing[t][z] )
	{
	  count++;
	  *img = (float *)mri_get_image( Input, t, z, MRI_FLOAT );
	  realign_matrix( (void **) img, dy, (long) ( dx * sizeof(float) ) );
	  for( y = 0; y < dy; y++ )
	    for( x = 0; x < dx; x++ )
	      {
		mean[y][x] += (double) img[y][x];
		stdv[y][x] += (double) ( img[y][x] * img[y][x] );
	      }
	}
    }
  
  /* Convert sums to appropriate averages */
  if( count < 2 )
    {
      Warning( 1,
	       "Less than 2 non-missing images for slice %ld --- leaving unchanged",
	       z );
    }
  else
    {
      for( y = 0; y < dy; y++ )
	for( x = 0; x < dx; x++ )
	  {
	    mean[y][x] /= (double) count;
	    stdv[y][x] = 
	      ( stdv[y][x] - (double) count *
		mean[y][x] * mean[y][x] ) / (double) ( count - 1 );
	    stdv[y][x] = ( stdv[y][x] <= 0.0 )? 
	      0.0: sqrt( stdv[y][x] );
	  }
      
    }
  *countOut= count;
}

static void calc_mean_stdv_complex( MRI_Dataset* Input, 
				    unsigned char** missing,
				    DComplex **c_mean, DComplex **c_stdv, 
				    FComplex **c_img,
				    int z, int dx, int dy, int dz, int dt,
				    long* countOut)
{
  long count, x, y, t;

  /* Initialize sums to zero */
  count = 0;
  for( y = 0; y < dy; y++ )
    for( x = 0; x < dx; x++ )
      c_mean[y][x].real = c_mean[y][x].imag = 
	c_stdv[y][x].real = c_stdv[y][x].imag = 0.0;
  
  /* Loop through images and sum up */
  for( t = 0; t < dt; t++ )
    {
      if( !missing[t][z] )
	{
	  count++;
	  *c_img = (FComplex *)
	    mri_get_image( Input, t, z, MRI_COMPLEX_FLOAT );
	  realign_matrix( (void **) c_img, dy, 
			  (long) ( dx * sizeof(FComplex) ) );
	  for( y = 0; y < dy; y++ )
	    for( x = 0; x < dx; x++ )
	      {
		c_mean[y][x].real += (double) c_img[y][x].real;
		c_mean[y][x].imag += (double) c_img[y][x].imag;
		c_stdv[y][x].real += (double) 
		  ( c_img[y][x].real * c_img[y][x].real );
		c_stdv[y][x].imag += (double) 
		  ( c_img[y][x].imag * c_img[y][x].imag );
	      }
	}
    }
      
  /* Convert sums to appropriate averages */
  if( count < 2 )
    {
      Warning( 1,
	       "Less than 2 non-missing images for slice %ld --- leaving unchanged",
	       dz );
    }
  else
    {
      for( y = 0; y < dy; y++ )
	for( x = 0; x < dx; x++ )
	  {
	    c_mean[y][x].real /= (double) count;
	    c_mean[y][x].imag /= (double) count;
	    c_stdv[y][x].real = 
	      ( c_stdv[y][x].real - (double) count *
		c_mean[y][x].real * c_mean[y][x].real ) / 
	      (double) ( count - 1 );
	    c_stdv[y][x].real = ( c_stdv[y][x].real < 0.0 )?
	      0.0: sqrt( c_stdv[y][x].real );
	    c_stdv[y][x].imag = 
	      ( c_stdv[y][x].imag - (double) count *
		c_mean[y][x].imag * c_mean[y][x].imag ) / 
	      (double) ( count - 1 );
	    c_stdv[y][x].imag = ( c_stdv[y][x].imag < 0.0 )?
	      0.0: sqrt( c_stdv[y][x].imag );
	  }
    }
  *countOut= count;
}

int main( argc, argv ) 
     int argc;
     char **argv;
{

  MRI_Dataset *Input = NULL, *Output = NULL;
  MRI_Dataset *Mean = NULL, *Stdv = NULL;
  char infile[512], hdrfile[512], outfile[512], parfile[512];
  char meanfile[512], stdvfile[512], iqrfile[512];
  float cutoff, maxproportion, fimsize;
  long dv, dx, dy, dz, dt;
  FILE *fp = NULL;
  FComplex **c_img = NULL, **c_corr_img = NULL;
  float **img = NULL, **corr_img = NULL;
  DComplex **c_mean = NULL, **c_stdv = NULL;
  double **mean = NULL, **stdv = NULL;
  long **num_outs = NULL;
  unsigned char **missing = NULL;
  long count, x, y, t, z;
  float ndiff, slope;
  FComplex c_ndiff, shift;
  int mean_supplied_flag= 0;
  int stdv_supplied_flag= 0;
  int iqr_supplied_flag= 0;

  /* Print version number */
  Message( "# %s\n", rcsid );

  progname= argv[0];

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
  cl_get( "headerout|h", "%option %s[%]", "outlier.mri", hdrfile );
  cl_get( "input|i", "%option %s[%]", "input.mri", infile );
  cl_get( "parameters|p", "%option %s[%]", "outlier.par", parfile );
  cl_get( "cutoff|c", "%option %f[%]", 3.5, &cutoff );
  cl_get( "badimage|b", "%option %f[%]", 0.02, &maxproportion );
  mean_supplied_flag= cl_get( "median|mean", "%option %s", meanfile );
  stdv_supplied_flag= cl_get( "stdv", "%option %s", stdvfile );
  iqr_supplied_flag= cl_get( "iqr", "%option %s", iqrfile );

  if (cl_cleanup_check()) {
    int i;
    fprintf(stderr,"%s: invalid argument in command line:\n    ",argv[0]);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(-1);
  }

  /*** End command-line parsing ***/

  /* Check for valid cut-off specification */
  if( cutoff <= 0 )
    Abort( "Invalid cut-off specification (%ld).  Must be positive.",
	   cutoff );

  /* Check for conflicts */
  if (stdv_supplied_flag && iqr_supplied_flag) 
    Abort( "It is an error to supply both stdv and iqr!\b");

  /* Open input dataset */
  if( !strcmp( infile, hdrfile ) )
    Abort( "Input and output files must be distinct." );
  Input = mri_open_dataset( infile, MRI_READ );

  /* Open parameter write-file and write format info */
  fp = efopen( parfile, "w" );
  fprintf(fp,"##Format: order:z_fastest, type:raw\n");
  fprintf(fp,"##Format: names:(outlier)\n");

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

  /* Check subsidiary files for compatibility */
  if (mean_supplied_flag) {
    if (!strcmp(meanfile, infile) || !strcmp(meanfile,hdrfile))
      Abort("%s: Mean/Median file is not distinct from input and output\n",
	    argv[0]);
    Mean= mri_open_dataset( meanfile, MRI_READ );
    check_sizes(meanfile, Mean, "mean/median",dv,dx,dy,dz,1);
  }
  if (stdv_supplied_flag) {
    if (!strcmp(stdvfile, infile) || !strcmp(stdvfile,hdrfile))
      Abort("%s: stdv file is not distinct from input and output\n",
	    argv[0]);
    Stdv= mri_open_dataset( stdvfile, MRI_READ );
    check_sizes(stdvfile, Stdv, "stdv",dv,dx,dy,dz,1);
  }
  if (iqr_supplied_flag) {
    if (!strcmp(iqrfile, infile) || !strcmp(iqrfile,hdrfile))
      Abort("%s: iqr file is not distinct from input and output\n",
	    argv[0]);
    Stdv= mri_open_dataset( iqrfile, MRI_READ );
    check_sizes(stdvfile, Stdv, "iqr",dv,dx,dy,dz,1);
  }

  /* Set output dataset */
  Output = mri_copy_dataset( hdrfile, Input );
  hist_add_cl( Output, argc, argv );
  mri_set_string( Output, "images.datatype", "float32" );
  mri_set_string( Output, "images.file", outfile );

  /* Read/Create missing image indicators */
  missing = get_missing( Output );

  /* Allocate image and parameter storage */
  if( dv == 1 )
    {
      img = (float **) emalloc( dy * sizeof(float *) );
      corr_img = Matrix( dy, dx, float );
      mean = Matrix( dy, dx, double );
      stdv = Matrix( dy, dx, double );
    }
  else
    {
      c_img = (FComplex **) emalloc( dy * sizeof(FComplex *) );
      c_corr_img = Matrix( dy, dx, FComplex );
      c_mean = Matrix( dy, dx, DComplex );
      c_stdv = Matrix( dy, dx, DComplex );
    }
  num_outs = Matrix( dt, dz, long );

  /* Initialize number of outliers per image to 0 */
  for( t = 0; t < dt; t++ )
    for( z = 0; z < dz; z++ )
      num_outs[t][z] = 0;

  /* Loop through slices */
  for( z = 0; z < dz; z++ ) {
    count= 0;

    /* GET MEANS AND STDVS, BY CALCULATION AND/OR FROM INPUT */

    if (!mean_supplied_flag || !(stdv_supplied_flag||iqr_supplied_flag)) {
      /* We have to calculate them ourselves */
      if (dv==1) 
	calc_mean_stdv( Input, missing, mean, stdv, img,
			z, dx, dy, dz, dt, &count );
      else 
	calc_mean_stdv_complex( Input, missing, c_mean, c_stdv, c_img,
				z, dx, dy, dz, dt, &count );      
    }

    if (dv == 1) {
      if (mean_supplied_flag) {
	/* overwrite mean with values from file */
	mri_read_chunk( Mean, "images", dx*dy, z*dx*dy, MRI_DOUBLE,
			mean[0] );
	realign_matrix( (void **) mean, dy, (long) ( dx * sizeof(double) ) );
      }
      if (stdv_supplied_flag) {
	/* overwrite stdv with values from file */
	mri_read_chunk( Stdv, "images", dx*dy, z*dx*dy, MRI_DOUBLE,
			stdv[0] );
	realign_matrix( (void **) stdv, dy, (long) ( dx * sizeof(double) ) );
      }
      if (iqr_supplied_flag) {
	/* overwrite stdv with values from file */
	mri_read_chunk( Stdv, "images", dx*dy, z*dx*dy, MRI_DOUBLE,
			stdv[0] );
	realign_matrix( (void **) stdv, dy, (long) ( dx * sizeof(double) ) );
	/* rescale appropriately, since cutoffs are specified in stdv's */
	for (x=0; x<dx; x++)
	  for (y=0; y<dy; y++) stdv[y][x] *= IQR_TO_STDV;
      }
    }
    else {
      if (mean_supplied_flag) {
	/* overwrite mean with values from file */
	mri_read_chunk( Mean, "images", 2*dx*dy, z*2*dx*dy, MRI_DOUBLE,
			c_mean[0] );
	realign_matrix( (void **) c_mean, dy, (long) (dx*sizeof(DComplex)) );
      }
      if (stdv_supplied_flag) {
	/* overwrite stdv with values from file */
	mri_read_chunk( Stdv, "images", 2*dx*dy, z*2*dx*dy, MRI_DOUBLE,
			c_stdv[0] );
	realign_matrix( (void **) c_stdv, dy, (long) (dx*sizeof(DComplex) ) );
      }
      if (iqr_supplied_flag) {
	/* overwrite stdv with values from file */
	mri_read_chunk( Stdv, "images", 2*dx*dy, z*2*dx*dy, MRI_DOUBLE,
			c_stdv[0] );
	realign_matrix( (void **) c_stdv, dy, (long) (dx*sizeof(DComplex) ) );
	/* rescale appropriately, since cutoffs are specified in stdv's */
	for (x=0; x<dx; x++)
	  for (y=0; y<dy; y++) {
	    c_stdv[y][x].real *= IQR_TO_STDV;
	    c_stdv[y][x].imag *= IQR_TO_STDV;
	  }
      }
    }
    
    if (count==0) count= dt;

    /* MAKE OUTLIER CORRECTIONS */
    
    /* Loop through images again */
    for( t = 0; t < dt; t++ )
      {
	
	if( dv == 1 )
	  {
	    *img = (float *)
	      mri_get_image( Input, t, z, MRI_FLOAT );
	    realign_matrix( (void **) img, dy, 
			    (long) ( dx * sizeof(float) ) );
	    if( missing[t][z] || ( count < 2 ) )
	      {
		mri_set_image( Output, t, z, MRI_FLOAT, *img );
	      }
	    else
	      {
		for( y = 0; y < dy; y++ )
		  for( x = 0; x < dx; x++ )
		    {
		      /* Calculate number of standard deviations */
		      /*   that observation lies from the mean   */
		      ndiff = ( !stdv[y][x] )? 0.0:
			( img[y][x] - mean[y][x] ) / stdv[y][x];
		      
		      /* If observation is beyond the cutoff  */
		      /*   number of standard deviations from */
		      /*   the mean, then pull the obervation */
		      /*   in to the cutoff value             */
		      if( ndiff < ( -1.0 * cutoff ) )
			{
			  corr_img[y][x] = mean[y][x] - cutoff * stdv[y][x];
			  corr_img[y][x] = ( corr_img[y][x] < 0.0 )?
			    0.0: corr_img[y][x];
			  num_outs[t][z]++;
			}
		      else if( ndiff > cutoff )
			{
			  corr_img[y][x] = mean[y][x] + cutoff * stdv[y][x];
			  num_outs[t][z]++;
			}
		      else
			corr_img[y][x] = img[y][x];
		    }
		
		/* Write out corrected image */
		mri_set_image( Output, t, z, MRI_FLOAT, *corr_img );
		
	      }
	  }
	else
	  {
	    *c_img = (FComplex *)
	      mri_get_image( Input, t, z, MRI_COMPLEX_FLOAT );
	    realign_matrix( (void **) c_img, dy, 
			    (long) ( dx * sizeof(FComplex) ) );
	    if( missing[t][z] || ( count < 2 ) )
	      {
		mri_set_image( Output, t, z, MRI_COMPLEX_FLOAT, *c_img );
	      }
	    else
	      {
		for( y = 0; y < dy; y++ )
		  for( x = 0; x < dx; x++ )
		    {
		      /* Calculate number of standard deviations */
		      /*   that observation lies from the mean   */
		      c_ndiff.real = ( !c_stdv[y][x].real )? 0.0:
			( c_img[y][x].real - c_mean[y][x].real ) / 
			c_stdv[y][x].real;
		      c_ndiff.imag = ( !c_stdv[y][x].imag )? 0.0:
			( c_img[y][x].imag - c_mean[y][x].imag ) / 
			c_stdv[y][x].imag;
		      
		      /* If the normalized (centered and scaled)  */
		      /*   observation does not have modulus less */
		      /*   than cutoff, then pull the observation */
		      /*   straight toward the mean, onto the     */
		      /*   ellipse defined by the cutoff value    */
		      if( Modulus( c_ndiff ) > cutoff )
			{
			  num_outs[t][z]++;
			  if( c_ndiff.real )
			    {
			      slope = c_ndiff.imag * c_stdv[y][x].imag /
				( c_ndiff.real * c_stdv[y][x].real );
			      shift.real = cutoff *
				sqrt( pow( c_stdv[y][x].real, 2.0 ) *
				      pow( c_stdv[y][x].imag, 2.0 ) /
				      ( pow( c_stdv[y][x].imag, 2.0 ) +
					pow( c_stdv[y][x].real, 2.0 ) *
					pow( slope, 2.0 ) ) );
			      if( c_ndiff.real < 0.0 )
				shift.real *= -1.0;
			      shift.imag = slope * shift.real;
			      c_corr_img[y][x].real = c_mean[y][x].real +
				shift.real;
			      c_corr_img[y][x].imag = c_mean[y][x].imag +
				shift.imag;
			    }
			  else
			    {
			      c_corr_img[y][x].real = 0.0;
			      c_corr_img[y][x].imag = 
				( c_ndiff.imag < 0.0 )?
				( c_mean[y][x].imag -
				  cutoff * c_stdv[y][x].imag ):
				( c_mean[y][x].imag +
				  cutoff * c_stdv[y][x].imag );
			    }
			}
		      else
			{
			  c_corr_img[y][x].real = c_img[y][x].real;
			  c_corr_img[y][x].imag = c_img[y][x].imag;
			}
		    }
		
		/* Write out corrected image */
		mri_set_image( Output, t, z, MRI_COMPLEX_FLOAT, 
			       *c_corr_img );
		
	      }
	    
	  }
      }

    /* END OUTLIER CORRECTIONS */

  }

  /* Declare bad images as missing */
  fimsize = (float) ( dx * dy );
  for( t = 0; t < dt; t++ )
    for( z = 0; z < dz; z++ )
      {
	if( ( (float) num_outs[t][z] / fimsize ) > maxproportion )
	  missing[t][z] = (unsigned char) 1;
      }
  mri_set_chunk( Output, "missing", (long) ( dt * dz ), 0,
		 MRI_UNSIGNED_CHAR, *missing );
  
      
  /* Write and close data-sets */
  mri_close_dataset( Input );
  mri_close_dataset( Output );
  
  /* Write and close parameter file */
  for( t = 0; t < dt; t++ )
    for( z = 0; z < dz; z++ )
      fprintf( fp, "%10ld\n", num_outs[t][z] );
  efclose( fp );
  
  /* Clean up */
  if (dv==1) {
    free(img);
    FreeMatrix(corr_img);
    FreeMatrix(mean);
    FreeMatrix(stdv)
  }
  else {
    free(c_img);
    FreeMatrix(c_corr_img);
    FreeMatrix(c_mean);
    FreeMatrix(c_stdv)
  }

  Message( "#      Outlier correction complete.\n" );
  return 0;

}

