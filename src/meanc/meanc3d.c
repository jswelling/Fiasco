/************************************************************
 *                                                          *
 *  meanc3d.c                                                 *
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
 *     3-05: meanc3d is derived from meanc, Joel Welling    *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"

static char rcsid[] = "$Id: meanc3d.c,v 1.2 2005/07/30 02:14:01 welling Exp $";

#define LOC(i,j,k,dx,dy,dz) ((k*dy + j)*dx+i)

typedef long long MissingCode;

const char* progname;
static int debug_flag= 0;
static int verbose_flag= 0;

static MissingCode makeAllSlicesMissingCode( long dz, long zmin, long zmax )
{
  long long result= 0;
  int z;
  if (dz>8*sizeof(MissingCode))
    Abort("%s: can't handle %d slices! ( %d max )\n",progname,
	  dz,8*sizeof(MissingCode));
  for (z=zmin; z<=zmax; z++) result |= 1<<z;
  if (debug_flag) fprintf(stderr,"All slices missing is code %016llx\n",
			  result);
  return (MissingCode)result;
}

static MissingCode makeMissingCode( unsigned char** missing, long dt,
				    long dz, long t)
{
  long long result= 0;
  int z;
  if (dz>8*sizeof(MissingCode))
    Abort("%s: can't handle %d slices! ( %d max )\n",progname,
	  dz,8*sizeof(MissingCode));
  for (z=0; z<dz; z++) 
    if (missing[t][z]) result |= 1<<z;
  if (debug_flag) fprintf(stderr,"Time %d has missing code %016llx\n",
			  t,result);
  return (MissingCode)result;
}

static int isSliceMissing( MissingCode code, long z )
{
  return( (code & (1<<z)) != 0 );
}

static int checkAllSlicesMissing( MissingCode code, 
				  MissingCode allSlicesMissingCode )
{
  return ( (code & allSlicesMissingCode) == allSlicesMissingCode );
}

static double get_mean(float* img, long dx, long dy, long dz,
		       long xmin, long xmax, long ymin,
		       long ymax, long zmin, long zmax,
		       MissingCode code)
{
  long count= 0;
  double sum= 0.0;
  int i,j,k;
  for ( k=zmin;k<=zmax;k++ ) {
    if (!(isSliceMissing(code,k))) {
      for ( j=ymin;j<=ymax;j++ )
	for ( i=xmin; i<=xmax; i++) {
	  count++;
	  sum+= img[LOC(i,j,k,dx,dy,dz)];
	}
    }
  }
  return sum/(double)count;
}

static double get_mean_complex(FComplex* img, long dx, long dy, long dz,
			       long xmin, long xmax, long ymin,
			       long ymax, long zmin, long zmax,
			       MissingCode code)
{
  long count= 0;
  double sum= 0.0;
  int i,j,k;
  for ( k=zmin;k<=zmax;k++ )
    if (!(isSliceMissing(code,k))) {
      for ( j=ymin;j<=ymax;j++ )
	for ( i=xmin; i<=xmax; i++) {
	  count++;
	  sum+= Modulus(img[LOC(i,j,k,dx,dy,dz)]);
	}
    }
  return sum/(double)count;
}

int main( int argc, char* argv[] ) 
{

  MRI_Dataset *Input = NULL, *Output = NULL;
  char infile[512], hdrfile[512];
  char parfile[512];
  long fixed_t, temp_fixed_t;
  long dv=0, dx, dy, dz, dt;
  FILE *fp = NULL;
  float *fixed_img=NULL, *img=NULL, *corr_img = NULL, *adj = NULL;
  unsigned char **missing = NULL;
  long xmin, xmax, ymin, ymax, zmin, zmax, meanc3d_area;
  int fixed_t_valid= 0;
  double fixed_mean, mean;
  long x, y, t, z;
  long blkSize;
  double sum;
  int allow_neg_means= 0;
  MissingCode currentMissingCode;
  MissingCode allSlicesMissingCode;

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
  cl_get( "headerout|h", "%option %s[%]", "meanc3d.mri", hdrfile );
  cl_get( "input|i", "%option %s[%]", "input.mri", infile );
  cl_get( "parameters|p", "%option %s[%]", "meanc3d.par", parfile );
  cl_get( "fixed|f", "%option %ld[%]", 0, &fixed_t );
  debug_flag= cl_present("dbg|debug");
  verbose_flag= cl_present("v|verbose");
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

  /* Print version number */
  if (verbose_flag) Message( "# %s\n", rcsid );

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
  blkSize= dv*dx*dy*dz;

  /* Make sure fixed image is in bounds */
  if( ( fixed_t < 0 ) || ( fixed_t >= dt ) )
    fixed_t = 0;

  /* Allocate image and parameter storage */
  if (!(corr_img= (float*)malloc(blkSize*sizeof(float))))
    Abort("%s: unable to allocate %d bytes!\n", argv[0],blkSize*sizeof(float));
  if (!(adj=(float*)malloc(dt*sizeof(float))))
    Abort("%s: unable to allocate %d bytes!\n", argv[0],dt*sizeof(float));

  /* Calculate area of region over which mean is calculated */
  xmin = (long) ( ( 3 * dx ) / 8 );
  xmax = (long) ( ( 5 * dx ) / 8 );
  ymin = (long) ( ( 3 * dy ) / 8 );
  ymax = (long) ( ( 5 * dy ) / 8 );
  zmin = (long) ( ( 3 * dz ) / 8 );
  zmax = (long) ( ( 5 * dz ) / 8 );

  allSlicesMissingCode= makeAllSlicesMissingCode(dz,zmin,zmax);

  /* Find a suitable standard image */
  temp_fixed_t= fixed_t;
  fixed_t_valid= 0;
  do {
    int found_a_missing_image= 0;
    for (z=zmin;z<=zmax;z++) 
      if (missing[temp_fixed_t][z]) {
	found_a_missing_image= 1;
	break;
      }
    if (!found_a_missing_image) {
      fixed_t_valid= 1;
      break;
    }
    temp_fixed_t++;
    if (temp_fixed_t==dt) temp_fixed_t= 0;
  } while (temp_fixed_t!=fixed_t);
  if (!fixed_t_valid)
    Abort("%s: all images have missing slices in central region!\n",
	  argv[0]);
  fixed_t= temp_fixed_t;
  if (verbose_flag)
    Message("# Settled on image %d as fixed image.\n",fixed_t);

  /* Grab fixed image */
  fixed_img= (float*)mri_get_chunk(Input,"images",blkSize,
				   fixed_t*blkSize, MRI_FLOAT);
  currentMissingCode= makeMissingCode(missing,dt,dz,fixed_t);
  if (dv==1) {
    fixed_mean= get_mean(fixed_img,dx,dy,dz,xmin,xmax,ymin,ymax,zmin,zmax,
			 currentMissingCode);
  }
  else {
    fixed_mean= get_mean_complex((FComplex*)fixed_img,dx,dy,dz,
				 xmin,xmax,ymin,ymax,zmin,zmax,
				 currentMissingCode);
  }

  /* Loop through times */
  for (t=0; t<dt; t++) {
    long long offset= t*blkSize;
    int img_valid;
    MissingCode thisMissingCode= makeMissingCode(missing,dt,dz,t);

    img= (float*)mri_get_chunk(Input, "images", blkSize, offset, MRI_FLOAT);
    memcpy(corr_img, img, blkSize*sizeof(float));
    
    if (!checkAllSlicesMissing(thisMissingCode,allSlicesMissingCode)) {
      long i;
      if (dv==1) 
	mean= get_mean(img,dx,dy,dz,xmin,xmax,ymin,ymax,zmin,zmax,
		       thisMissingCode);
      else
	mean= get_mean_complex((FComplex*)img,dx,dy,dz,
			       xmin,xmax,ymin,ymax,zmin,zmax,
			       thisMissingCode);
      if (mean<0.0 && !allow_neg_means) {
	Warning(1,"Image %d mean is negative (%lf)\n",t,mean);
	for (z=0; z<dz; z++) missing[t][z]= 1; 
      }
      else {
	if (mean==0.0) Abort("Invalid mean at t=%d\n",t);
	if (thisMissingCode != currentMissingCode) {
	  if (debug_flag) 
	    fprintf(stderr,"Regenerating fixed mean for new code %016llx\n",
		    thisMissingCode);
	  currentMissingCode= thisMissingCode;
	  if (dv==1) {
	    fixed_mean= get_mean(fixed_img,dx,dy,dz,xmin,xmax,ymin,ymax,
				 zmin,zmax,currentMissingCode);
	  }
	  else {
	    fixed_mean= get_mean_complex((FComplex*)fixed_img,dx,dy,dz,
					 xmin,xmax,ymin,ymax,zmin,zmax,
					 currentMissingCode);
	  }
	}
	adj[t]= fixed_mean/mean;
	for (i=0; i<blkSize; i++) corr_img[i] *= adj[t];
      }
    }
    else {
      /* This completely missing image has already been copied 'correctly' */
      adj[t]= 1.0;
    }
    mri_set_chunk(Output,"images",blkSize,offset,MRI_FLOAT,corr_img);
  }

  /* Write out missing chunk, since it may have changed. */
  mri_set_chunk( Output, "missing", (long) ( dt * dz ), 0,
		 MRI_UNSIGNED_CHAR, *missing );

  /* Write and close data-sets */
  mri_close_dataset( Input );
  mri_close_dataset( Output );
  
  /* Write parameter file */
  fp = efopen( parfile, "w" );
  fprintf(fp,"##Format: order:t_only type:raw names:(meanc3d)\n");
  for( t = 0; t < dt; t++ )
      fprintf( fp, "%15.6f\n", adj[t] );
  efclose( fp );
  
  if (verbose_flag) Message( "#      Mean adjustment complete.\n" );
  exit(0);

}

