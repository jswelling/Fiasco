/************************************************************
 *                                                          *
 *  stats.c                                                   *
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
 *  Original programming by Audris Mockus 2-95              *
 *     5-96: Pittsburgh Format, Mark Fitzgerald             *
 ************************************************************/
/*************************************************************

  DESCRIPTION OF STAT.C

  stat.c calculates the mean and standard deviation for all
    experimental conditions, as well as all pair-wise, 
    pooled, two-sample t-statistics

  stat.m [-condition Experimental-conditions-file]
           [-input Input-header-file]
           [-meanprefix MPath] [-stdvprefix SPath]
           [-tsprefix TPath] [-maxtpairs MaxPairs]

**************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"

static char rcsid[] = "$Id: stats.c,v 1.21 2004/12/09 22:38:39 welling Exp $";

#define MY_MAX_T 999999.0

#define MAX_TTESTS 50

static void delete_missing_chunk( MRI_Dataset* ds )
{
  const char* key;
  mri_iterate_over_keys(ds);
  while ((key = mri_next_key(ds)) != NULL)
    if (!strncmp(key,"missing",strlen("missing"))) mri_remove(ds,key);
}

int main( int argc, char** argv ) 
{

  MRI_Dataset *Input = NULL, *POutput = NULL, *ROutput = NULL;
  FILE *fp = NULL;
  char condfile[512], infile[512];
  char meanfile[512], stdvfile[512], tsfile[512];
  long dx, dy, dz, dt;
  long num_conds, **conds = NULL, tmpcond;
  double **image = NULL, **fimage = NULL;
  double ****mean = NULL, ****stdv = NULL;
  double ***grandMean= NULL;
  double ***grandStdv= NULL;
  long **counts = NULL;
  long *grandCounts= NULL;
  long *dof= NULL;
  unsigned char **missing = NULL;
  long ec, ec2, linenum, m, x, y, t, z;
  char scanline[512], outfile[562];
  double pooled_stdv, mean_diff, total_var;
  int num_t_pairs;
  int max_ttests;


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
  cl_get( "condition|c", "%option %s[%]", "newsplit", condfile );
  cl_get( "input|i", "%option %s[%]", "input.mri", infile );
  cl_get( "meanprefix|m", "%option %s[%]", "", meanfile );
  cl_get( "stdvprefix|s", "%option %s[%]", "", stdvfile );
  cl_get( "tsprefix|t", "%option %s[%]", "", tsfile );
  cl_get( "maxtpairs|p", "%option %d[%]", MAX_TTESTS, &max_ttests );

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
  Input = mri_open_dataset( infile, MRI_READ );

  /* Check that program will function on data-set */
  if( !mri_has( Input, "images" ) )
    Abort( "%s operates only on standard images.", argv[0] );
  if( mri_has( Input, "images.dimensions" ) )
    {
      if( strcmp( mri_get_string( Input, "images.dimensions" ), "xyzt" ) &&
	  ( strcmp( mri_get_string( Input, "images.dimensions" ), "vxyzt" )
	    || !mri_has( Input, "images.extent.v" ) ||
	    ( mri_get_int( Input, "images.extent.v" ) != 1 ) ) )
	Abort( "%s works only on xyzt or vxyzt images with v = 1.",
	       argv[0] );
    }
  else
    Abort( "%s does not have the images.dimensions key.", infile );

  /* Set parameters in local variables; grab input missing info */
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
  missing = get_missing( Input );

  /* Open experimental conditions file */
  fp = efopen( condfile, "r" );
  
  /* Allocate image and experimental conditions storage */
  conds = Matrix( dt, dz, long );
  image = Matrix( dy, dx, double );
  fimage = Matrix( dy, dx, double );

  /* Initialize experimental condition memberships for images to missing */
  for( t = 0; t < dt; t++ )
    for( z = 0; z < dz; z++ )
      conds[t][z] = 0;

  /* Read in experimental condition memberships for images */
  linenum = 0;
  while( !feof( fp ) )
    {
      linenum++;
      fscanf( fp, "%510[^\n]%*[\n]", scanline );
      m = sscanf( scanline, "%ld %ld %ld", &t, &z, &tmpcond );
      if( ( m == 3 ) && ( t >= 0 ) && ( t < dt ) && 
	  ( z >= 0 ) && ( z < dz ) )
	{
	  conds[t][z] = ( missing[t][z] )? 0: tmpcond;
	}
      else
	{
	  Warning( 1,
		   "Invalid image/slice number: line %ld of %s.  Ignoring.\n",
		   linenum, condfile );
	}
      continue;
    }
  
  /* Done reading conditions */
  efclose( fp );

  /* Find total number of conditions */
  num_conds = 0;
  for( t = 0; t < dt; t++ )
    for( z = 0; z < dz; z++ )
      num_conds = ( conds[t][z] > num_conds )? conds[t][z]: num_conds;
  num_conds++;

  /* Allocate space to calculate means and standard deviations */
  mean = Matrix( num_conds, dz, double ** );
  stdv = Matrix( num_conds, dz, double ** );
  counts = Matrix( num_conds, dz, long );
  if (!(dof=(long*)malloc(dz*sizeof(long))))
    Abort("%s: unable to allocate %d bytes!\n",argv[0],dz*sizeof(long));
  for( ec = 0; ec < num_conds; ec++ )
    for( z = 0; z < dz; z++ )
      {
	mean[ec][z] = Matrix( dy, dx, double );
	stdv[ec][z] = Matrix( dy, dx, double );
      }
  if (!(grandMean= (double***)malloc(dz*sizeof(double**))))
    Abort("%s: unable to allocate %d bytes!\n",argv[0],dz*sizeof(double**));
  if (!(grandStdv= (double***)malloc(dz*sizeof(double**))))
    Abort("%s: unable to allocate %d bytes!\n",argv[0],dz*sizeof(double**));
  if (!(grandCounts= (long*)malloc(dz*sizeof(long))))
    Abort("%s: unable to allocate %d bytes!\n",argv[0],dz*sizeof(long));
  for (z=0; z<dz; z++) {
    grandMean[z]= Matrix(dy,dx,double);
    grandStdv[z]= Matrix(dy,dx,double);
  }

  /* CALCULATE MEANS AND STANDARD DEVIATIONS */
  /*    FOR EACH EXPERIMENTAL CONDITION      */

  /* Initialize sums and counts to zero */
  for( ec = 0; ec < num_conds; ec++ )
    for( z = 0; z < dz; z++ )
      {
	counts[ec][z] = 0;
	for( y = 0; y < dy; y++ )
	  for( x = 0; x < dx; x++ )
	    mean[ec][z][y][x] = stdv[ec][z][y][x] = 0.0;
      }
  for (z=0; z<dz; z++) {
    grandCounts[z]= 0;
    for( y = 0; y < dy; y++ )
      for( x = 0; x < dx; x++ )
	grandMean[z][y][x]= grandStdv[z][y][x]= 0.0;
    
  }

  /* Loops through images and slices */
  for( t = 0; t < dt; t++ )
    for( z = 0; z < dz; z++ )
      {

	/* Read next image */
	*image = (double *)
	  mri_get_image( Input, t, z, MRI_DOUBLE );
	realign_matrix( (void **) image, dy,
			(long) ( dx * sizeof(double) ) );
	
	/* Pixelwise, add image to mean
	 * of appropriate experimental condition            
	 */
	for( y = 0; y < dy; y++ )
	  for( x = 0; x < dx; x++ )
	    {
	      mean[conds[t][z]][z][y][x] += (double) image[y][x];
	      if (conds[t][z] != 0)
		grandMean[z][y][x] += (double)image[y][x];
	    }

	/* Keep count of images per condition per slice */
	counts[conds[t][z]][z]++;
	if (conds[t][z]!=0) grandCounts[z]++;

      }

  /* Calculate means from sums */
  /*   (ignoring missing condition henceforth)         */
  for( ec = 1; ec < num_conds; ec++ )
    for( z = 0; z < dz; z++ )
      {
	if( counts[ec][z] )
	  {

	    /* Calculate means */
	    for( y = 0; y < dy; y++ )
	      for( x = 0; x < dx; x++ )
		{
		  mean[ec][z][y][x] /= (double) counts[ec][z];
		}
	  }
	else
	  {
	    Warning( 1, "No data for condition %ld, slice %ld.\n", ec, z );
	  }
      }
  for (z=0; z<dz; z++) {
    if (grandCounts[z]) {
      for( y = 0; y < dy; y++ )
	for( x = 0; x < dx; x++ )
	  {
	    grandMean[z][y][x] /= (double) grandCounts[z];
	  }
    }
    else {
      Warning( 1, "No data for slice %ld in any condition!\n", z );
    }
  }

  /* Loops through images and slices again */
  for( t = 0; t < dt; t++ )
    for( z = 0; z < dz; z++ )
      {

	/* Read next image */
	*image = (double *)
	  mri_get_image( Input, t, z, MRI_DOUBLE );
	realign_matrix( (void **) image, dy,
			(long) ( dx * sizeof(double) ) );
	
	/* Pixelwise, add image and squared deviation to corresponding */
	/*   sum of appropriate experimental condition            */
	for( y = 0; y < dy; y++ )
	  for( x = 0; x < dx; x++ )
	    {
	      double val= (double)image[y][x] - mean[conds[t][z]][z][y][x];
	      stdv[conds[t][z]][z][y][x] += val*val;
	      if (conds[t][z] != 0) {
		val= (double)image[y][x] - grandMean[z][y][x];
		grandStdv[z][y][x] += val*val;
	      }
	    }
      }

  /* Calculate standard deviations from sums */
  /*   (ignoring missing condition henceforth)         */
  for( ec = 1; ec < num_conds; ec++ )
    for( z = 0; z < dz; z++ )
      {
	if( counts[ec][z] )
	  {

	    if( counts[ec][z] > 1 )
	      {
		/* Calculate standard deviations */
		for( y = 0; y < dy; y++ )
		  for( x = 0; x < dx; x++ )
		    {
		      stdv[ec][z][y][x] = 
			sqrt(stdv[ec][z][y][x] / ((double)counts[ec][z] - 1));
		    }
	      }
	    else
	      {
		/* Can't calculate standard deviation with a single image */
		Warning( 1,
			 "Standard deviation invalid for condition %ld, slice %ld --- only one image.\n",
			 ec, z );
		for( y = 0; y < dy; y++ )
		  for( x = 0; x < dx; x++ )
		    stdv[ec][z][y][x] = 0.0;
	      }

	  }
	else
	  {
	    /* Already emitted a warning */
	  }

      }
  for (z=0; z<dz; z++) {
    if (grandCounts[z]) {
      if (grandCounts[z]>1) {
	/* Calculate standard deviations */
	for( y = 0; y < dy; y++ )
	  for( x = 0; x < dx; x++ )
	    {
	      grandStdv[z][y][x] = 
		sqrt(grandStdv[z][y][x] / ((double)grandCounts[z] - 1));
	    }
      }
      else {
	/* Can't calculate standard deviation with a single image */
	Warning( 1,
		 "Grand standard deviation invalid for slice %ld --- only one image.\n",
		 z );
	for( y = 0; y < dy; y++ )
	  for( x = 0; x < dx; x++ )
	    grandStdv[z][y][x] = 0.0;
      }
    }
    else {
      /* Already emitted a warning */
    }
  }

  /* END MEAN/STANDARD DEVIATION CALCULATION */



  /* WRITE OUT MEANS AND STANDARD DEVIATIONS */

  /* Create a prototype output dataset (the first mean) */
  sprintf( outfile, "%sGrandMean.mri", meanfile );
  POutput = mri_copy_dataset( outfile, Input );
  hist_add_cl( POutput, argc, argv );
  mri_create_chunk( POutput, "counts" );
  delete_missing_chunk( POutput );
  mri_set_string( POutput, "counts.datatype", "int32" );
  mri_set_string( POutput, "counts.dimensions", "z" );
  mri_set_int( POutput, "counts.extent.z", (int) dz );
  mri_set_chunk( POutput, "counts", (int) dz, 0, MRI_INT, grandCounts );
  mri_create_chunk( POutput, "dof" );
  mri_set_string( POutput, "dof.datatype", "int32" );
  mri_set_string( POutput, "dof.dimensions", "z" );
  mri_set_int( POutput, "dof.extent.z", (int) dz );
  mri_create_chunk( POutput, "images" );
  mri_set_string( POutput, "images.file", ".dat" );
  mri_set_string( POutput, "images.datatype", "float64" );
  mri_set_string( POutput, "images.dimensions", "xyzt" );
  mri_set_int( POutput, "images.extent.x", dx );
  mri_set_int( POutput, "images.extent.y", dy );
  mri_set_int( POutput, "images.extent.z", dz );
  mri_set_int( POutput, "images.extent.t", 1 );

  /* Degrees of freedom for mean are 1 */
  for (z=0; z<dz; z++) dof[z]= 1;
  mri_set_chunk( POutput, "dof", (int) dz, 0, MRI_INT, dof );

  /* Write out grand means */
  for( z = 0; z < dz; z++ ) {
    for( y = 0; y < dy; y++ )
      for( x = 0; x < dx; x++ )
	{
	  fimage[y][x] = grandMean[z][y][x];
	}
    mri_set_image( POutput, 0, z, MRI_DOUBLE, *fimage );
  }

  /* Write out grand stdv if there are enough conditions that it matters */
  if (ec>1) {
    
    /* Name properly and copy prototype */
    sprintf( outfile, "%sGrandStdv.mri", stdvfile);
    ROutput = mri_copy_dataset( outfile, POutput );
    
    /* Write out standard deviation */
    for( z = 0; z < dz; z++ )
      {
	for( y = 0; y < dy; y++ )
	  for( x = 0; x < dx; x++ )
	    {
	      fimage[y][x] = grandStdv[z][y][x];
	    }
	mri_set_image( ROutput, 0, z, MRI_DOUBLE, *fimage );
      }
    mri_set_chunk( ROutput, "counts", (int)dz, 0, MRI_INT, grandCounts );
    /* dof for stdv is counts-1, assuming counts >= 1 */
    for (z=0; z<dz; z++) 
      dof[z]= ((grandCounts[z] > 0) ? (grandCounts[z]-1) : 0);
    mri_set_chunk( ROutput, "dof", (int) dz, 0, MRI_INT, dof );
    
    /* Close mean dataset */
    mri_close_dataset( ROutput );
  }

  /* Write out the means and standard deviations */
  for( ec = 1; ec < num_conds; ec++ )
    {
      /* Write out mean dataset */

      /* Name properly and copy prototype */
      sprintf( outfile, "%sMean_%ld.mri", meanfile, ec );
      ROutput = mri_copy_dataset( outfile, POutput );
      
      /* Write out means */
      for( z = 0; z < dz; z++ )
	{
	  /* Convert mean to double */
	  for( y = 0; y < dy; y++ )
	    for( x = 0; x < dx; x++ )
	      {
		fimage[y][x] = (double) mean[ec][z][y][x];
	      }
	  mri_set_image( ROutput, 0, z, MRI_DOUBLE, *fimage );
	}
      mri_set_chunk( ROutput, "counts", (int)dz, 0, MRI_INT, counts[ec] );
      /* dof is 1, the same for all means */

      /* Close mean dataset */
      mri_close_dataset( ROutput );
      

      /* Write out standard deviation dataset */
      
      /* Name properly and copy prototype */
      sprintf( outfile, "%sStdv_%ld.mri", stdvfile, ec );
      ROutput = mri_copy_dataset( outfile, POutput );
      
      /* Write out standard deviation */
      for( z = 0; z < dz; z++ )
	{
	  /* Convert mean to double */
	  for( y = 0; y < dy; y++ )
	    for( x = 0; x < dx; x++ )
	      {
		fimage[y][x] = (double) stdv[ec][z][y][x];
	      }
	  mri_set_image( ROutput, 0, z, MRI_DOUBLE, *fimage );
	}
      mri_set_chunk( ROutput, "counts", (int)dz, 0, MRI_INT, counts[ec] );
      /* dof for stdv is counts-1, assuming counts >= 1 */
      for (z=0; z<dz; z++) 
	dof[z]= ((counts[ec][z] > 0) ? (counts[ec][z]-1) : 0);
      mri_set_chunk( ROutput, "dof", (int) dz, 0, MRI_INT, dof );

      
      /* Close stdv dataset */
      mri_close_dataset( ROutput );
    }

  if (ec>1) {
    mri_close_dataset( POutput ); /* flush grand mean to disk */
  }
  else {
    /* The grand mean is of no interest with only one condition */
    mri_destroy_dataset( POutput );
  }

  /* END WRITE-OUT OF MEANS AND STANDARD DEVIATIONS */

  num_t_pairs= (num_conds*(num_conds-1))/2;
  if (num_t_pairs > max_ttests) {
    Message("%s: Skipping pairwise T-statistics; %d pairs is too many!\n",
	    argv[0],num_t_pairs);
  }
  else {
    
    /* WRITE OUT ALL PAIR-WISE T-STATISTICS */
    
    /* Create a prototype output dataset (the first pairwise t-statistic) */
    sprintf( outfile, "%sTmap_0-0.mri", tsfile );
    POutput = mri_copy_dataset( outfile, Input );
    hist_add_cl( POutput, argc, argv );
    delete_missing_chunk( POutput );
    mri_create_chunk( POutput, "counts1" );
    mri_set_string( POutput, "counts1.datatype", "int32" );
    mri_set_string( POutput, "counts1.dimensions", "z" );
    mri_set_int( POutput, "counts1.extent.z", (int) dz );
    mri_set_chunk( POutput, "counts1", (int) dz, 0, MRI_INT, counts[0] );
    mri_create_chunk( POutput, "counts2" );
    mri_set_string( POutput, "counts2.datatype", "int32" );
    mri_set_string( POutput, "counts2.dimensions", "z" );
    mri_set_int( POutput, "counts2.extent.z", (int) dz );
    mri_set_chunk( POutput, "counts2", (int) dz, 0, MRI_INT, counts[0] );
    mri_create_chunk( POutput, "dof" );
    mri_set_string( POutput, "dof.datatype", "int32" );
    mri_set_string( POutput, "dof.dimensions", "z" );
    mri_set_int( POutput, "dof.extent.z", (int) dz );
    for (z=0; z<dz; z++) dof[z]= 0;
    mri_set_chunk( POutput, "dof", (int) dz, 0, MRI_INT, dof );
    mri_create_chunk( POutput, "images" );
    mri_set_string( POutput, "images.file", "Tmap.0-0.dat" );
    mri_set_string( POutput, "images.datatype", "float64" );
    mri_set_string( POutput, "images.dimensions", "xyzt" );
    mri_set_int( POutput, "images.extent.x", dx );
    mri_set_int( POutput, "images.extent.y", dy );
    mri_set_int( POutput, "images.extent.z", dz );
    mri_set_int( POutput, "images.extent.t", 1 );
    
    /* Write out dummy t-statistics, just to complete prototype */
    for( z = 0; z < dz; z++ )
      mri_set_image( POutput, 0, z, MRI_DOUBLE, *fimage );
    
    /* Write out the t-statistics */
    for( ec = 1; ec < num_conds; ec++ )
      for( ec2 = (long) ( ec + 1 ); ec2 < num_conds; ec2++ )
	{
	  /* Write out tmap dataset */
	  
	  /* Name properly and copy prototype */
	  sprintf( outfile, "%sTmap_%ld-%ld.mri", tsfile, ec, ec2 );
	  ROutput = mri_copy_dataset( outfile, POutput );
	  
	  /* Loop through slices */
	  for( z = 0; z < dz; z++ )
	    {
	      /* Calculate t-statistics */
	      if( counts[ec][z] && counts[ec2][z] &&
		  ( ( counts[ec][z] + counts[ec2][z] ) > 2 ) )
		{
		  /* Enough data to calculate t-statistics properly */
		  for( y = 0; y < dy; y++ )
		    for( x = 0; x < dx; x++ )
		      {
			/* pooled_stdv = Normalized pooled standard deviation */
			pooled_stdv = (double)
			  sqrt( ( ( (stdv[ec][z][y][x] * stdv[ec][z][y][x] * 
				     ( counts[ec][z] - 1 )) +
				    (stdv[ec2][z][y][x] * stdv[ec2][z][y][x] *
				     ( counts[ec2][z] - 1 )) ) /
				  (double) ( counts[ec][z] + 
					     counts[ec2][z] - 2 ) ) * 
				(double) ( 1.0 / counts[ec][z] + 
					   1.0 / counts[ec2][z] ) );
			mean_diff = mean[ec][z][y][x] - mean[ec2][z][y][x];
			fimage[y][x] = ( pooled_stdv )? 
			  (double) ( mean_diff / pooled_stdv ):
			  ( mean_diff > 0 )? MY_MAX_T: ( mean_diff < 0 )?
			  ( -1.0 * MY_MAX_T ): 0.0;
		      }
		  dof[z]= counts[ec][z] + counts[ec2][z] - 2;
		}
	      else
		{
		  /* Not enough data: t-statistics is just indicators */
		  for( y = 0; y < dy; y++ )
		    for( x = 0; x < dx; x++ )
		      fimage[y][x] = ( ( mean[ec][z][y][x] - 
					 mean[ec2][z][y][x] ) > 0 )?
			MY_MAX_T: ( ( mean[ec][z][y][x] - 
				      mean[ec2][z][y][x] ) < 0 )?
			( -1.0 * MY_MAX_T ): 0.0;
		  dof[z]= 0;
		}
	      
	      /* Write out t-statistics image */
	      mri_set_image( ROutput, 0, z, MRI_DOUBLE, *fimage );
	      
	    }
	  mri_set_chunk( ROutput, "counts1", (int) dz, 0, MRI_INT, 
			 counts[ec] );
	  mri_set_chunk( ROutput, "counts2", (int) dz, 0, MRI_INT, 
			 counts[ec2] );
	  mri_set_chunk( ROutput, "dof", (int) dz, 0, MRI_INT, dof );
    
	  /* Close t-statistics dataset */
	  mri_close_dataset( ROutput );
	  

	  /* Write out percent signal change dataset */
	  
	  /* Name properly and copy prototype */
	  sprintf( outfile, "%sPctSC_%ld-%ld.mri", tsfile, ec, ec2 );
	  ROutput = mri_copy_dataset( outfile, POutput );
	  
	  /* Loop through slices */
	  for( z = 0; z < dz; z++ )
	    {
	      /* Percent signal change has 1 DOF */
	      dof[z]= 1;

	      /* Calculate percent change */
	      if( counts[ec][z] && counts[ec2][z] &&
		  ( ( counts[ec][z] + counts[ec2][z] ) > 2 ) )
		{
		  /* Enough data to calculate percent change properly */
		  for( y = 0; y < dy; y++ )
		    for( x = 0; x < dx; x++ )
		      {

			mean_diff = mean[ec][z][y][x] - mean[ec2][z][y][x];

			fimage[y][x]= (mean[ec2][z][y][x]) ?
			  (100.0*mean_diff/mean[ec2][z][y][x]) :
			  ((mean_diff>0.0) ? 100.0 : 
			   ((mean_diff<0.0) ? -100.0 : 0.0));
		      }
		}
	      else
		{
		  /* Not enough data: percent change is just indicators */
		  for( y = 0; y < dy; y++ )
		    for( x = 0; x < dx; x++ )
		      fimage[y][x] = ( ( mean[ec][z][y][x] - 
					 mean[ec2][z][y][x] ) > 0 )?
			100.0: ( ( mean[ec][z][y][x] - 
				      mean[ec2][z][y][x] ) < 0 )?
			( -100.0 ): 0.0;
		}
	      
	      /* Write out percent change image */
	      mri_set_image( ROutput, 0, z, MRI_DOUBLE, *fimage );
	      
	    }
	  mri_set_chunk( ROutput, "counts1", (int) dz, 0, MRI_INT, 
			 counts[ec] );
	  mri_set_chunk( ROutput, "counts2", (int) dz, 0, MRI_INT, 
			 counts[ec2] );
	  mri_set_chunk( ROutput, "dof", (int) dz, 0, MRI_INT, dof );
	  
	  /* Close percent change dataset */
	  mri_close_dataset( ROutput );

	  
	  /* Write out error in percent signal change dataset */
	  
	  /* Name properly and copy prototype */
	  sprintf( outfile, "%sEPctSC_%ld-%ld.mri", tsfile, ec, ec2 );
	  ROutput = mri_copy_dataset( outfile, POutput );
	  
	  /* Loop through slices */
	  for( z = 0; z < dz; z++ )
	    {
	      /* Calculate percent change */
	      if( counts[ec][z] && counts[ec2][z] &&
		  ( ( counts[ec][z] + counts[ec2][z] ) > 2 ) )
		{
		  /* Enough data to calculate percent change properly */
		  for( y = 0; y < dy; y++ )
		    for( x = 0; x < dx; x++ )
		      {

			/* Following is accurate to leading order */
			total_var= (((stdv[ec][z][y][x]*stdv[ec][z][y][x])
				    /(double)counts[ec][z])
				    +((stdv[ec2][z][y][x]*stdv[ec2][z][y][x])
				    /(double)counts[ec2][z]));
			fimage[y][x]= 
			  100.0 * sqrt(total_var) / mean[ec2][z][y][x];
		      }
		  dof[z]= counts[ec][z] + counts[ec2][z] - 2;
		}
	      else
		{
		  /* Not enough data: percent change is just indicators */
		  for( y = 0; y < dy; y++ )
		    for( x = 0; x < dx; x++ )
		      fimage[y][x] = 100.0;
		  dof[z]= 0;
		}
	      
	      /* Write out percent change image */
	      mri_set_image( ROutput, 0, z, MRI_DOUBLE, *fimage );
	      
	    }
	  mri_set_chunk( ROutput, "counts1", (int) dz, 0, MRI_INT, 
			 counts[ec] );
	  mri_set_chunk( ROutput, "counts2", (int) dz, 0, MRI_INT, 
			 counts[ec2] );
	  mri_set_chunk( ROutput, "dof", (int) dz, 0, MRI_INT, dof );

	  /* Close percent change dataset */
	  mri_close_dataset( ROutput );
	  
	}
    
    /* Finished with prototype */
    mri_destroy_dataset( POutput );

    /* END WRITE-OUT OF ALL PAIR-WISE T-STATISTICS */
  }

  /* Done with input */
  mri_close_dataset( Input );

  Message("#      Means, standard deviations, and t-maps complete.\n" );
  exit(0);

  return 0; 
}

