/************************************************************
 *                                                          *
 *  iwarp.c                                                *
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
 *     2-99: 3D version, Joel Welling                       *
 *     3-04: ireg3d becomes iwarp, Joel Welling             *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"

static char rcsid[] = "$Id: iwarp.c,v 1.6 2005/04/06 19:01:35 welling Exp $";

/* Access for 3D arrays */
#define MEM_BCK(matrix,nx,ny,nz,x,y,z) matrix[((((z)*ny)+(y))*nx)+(x)]
#define MEM(matrix,nx,ny,nz,x,y,z) matrix[((((x)*ny)+(y))*nz)+(z)]

typedef struct warppar_struct {
  Transform v;
  double mse;
} WarpPar;

static char* progname= NULL;
static int debug_flag= 0;
static int verbose_flag= 0;
static InterpolatorType interpolatorType= INTRP_LINEAR; /* default */

static void initializeWarpPar( WarpPar* p )
{
  trans_identity(p->v);
  p->mse= -1.0;
}

static void copyWarpPar( WarpPar* out, const WarpPar* in )
{
  trans_copy(out->v, in->v);
  out->mse= in->mse;
}

static int parseOneWarp(char* line, WarpPar* p, long* t)
{
  int numRead= sscanf(line, 
    "%ld%lg%lg%lg%lg%lg%lg%lg%lg%lg%lg%lg%lg%lg%lg%lg%lg%lg\n",
		      t,
		      &(p->v[0]), &(p->v[1]), &(p->v[2]), &(p->v[3]),
		      &(p->v[4]), &(p->v[5]), &(p->v[6]), &(p->v[7]),
		      &(p->v[8]), &(p->v[9]), &(p->v[10]), &(p->v[11]),
		      &(p->v[12]), &(p->v[13]), &(p->v[14]), &(p->v[15]),
		      &(p->mse));
  return (numRead==18);
}

static void load_warp_params(char* parfile, WarpPar* par, long dt)
{
  FILE *fp = NULL;
  char scanline[512];
  long linenum, tt;
  WarpPar p;

  fp = efopen( parfile, "r" );
  linenum = -1;
  while( !feof( fp ) && !ferror( fp ) )
    {
      linenum++;

      /* Scan a line, ignoring comments (which begin with '#') */
      if (fgets(scanline, sizeof(scanline), fp)
	  && strlen(scanline)>0 && scanline[0] != '#') {
	if (!parseOneWarp(scanline, &p, &tt)) {
	    Warning( 1, "Line %ld of %s contained an error -- Ignoring.\n",
		     linenum, parfile );
	}
	else if( ( tt < 0 ) || ( tt >= dt ) ) {
	  Warning( 1,
		   "Image out of bounds (%ld) for line %ld of %s -- Ignoring.\n",
		   tt, linenum, parfile );
	}
	else {
	  copyWarpPar(&(par[tt]),&p);
	  if (debug_flag)
	    fprintf(stderr,"%d: loaded (%g %g %g ...), mse %g\n",
		    tt,par[tt].v[0], par[tt].v[1], par[tt].v[2], par[tt].mse);
	}
      }
    }
  efclose( fp );
}

int main( int argc, char* argv[] ) 
{
  MRI_Dataset *Input = NULL, *Output = NULL;
  char infile[512], hdrfile[512], parfile[512];
  char modeString[512];
  long dv, dx, dy, dz, dt;
#ifdef never
  FComplex *c_image = NULL, *c_corr_image = NULL, *c_image_in= NULL;
  float *image_in = NULL, *corr_image = NULL;
#endif
  double* image_in= NULL, *corr_image= NULL;
  char* check = NULL;
  long x, y, t, z;
  WarpPar* par= NULL;
  long block_offset;
  int block_size;
  const char* dimstr= NULL;
  double length_x, length_y, length_z;
  double xvoxel, yvoxel, zvoxel;
  Interpolator* interpolator= NULL;

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

  if (cl_present("trilinear"))
    Abort("Option -trilinear has been replaced by -interp linear.  Please see help file.\n");

  if (cl_present("closest"))
    Abort("Option -closest has been replaced by -interp closest.  Please see help file.\n");

  /* Get filenames */
  cl_get( "headerout|h", "%option %s[%]", "iwarp.mri", hdrfile );
  cl_get( "input|i", "%option %s[%]", "input.mri", infile );
  cl_get( "parameters|p", "%option %s[%]", "warp.par", parfile );
  if (!cl_get( "x|xvoxel", "%option %lf", &xvoxel )) {
    fprintf(stderr,"%s: required argument xvoxel omitted.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  if (!cl_get( "y|yvoxel", "%option %lf", &yvoxel )) {
    fprintf(stderr,"%s: required argument yvoxel omitted.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  if (!cl_get( "z|zvoxel", "%option %lf", &zvoxel )) {
    fprintf(stderr,"%s: required argument zvoxel omitted.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  debug_flag= cl_present("debug");
  verbose_flag= cl_present("verbose|v");
  if (cl_get("interp","%option %s",modeString)) {
    interpolatorType= intrp_typeFromName(modeString);
    if (interpolatorType==INTRP_UNKNOWN) {
      fprintf(stderr,"%s: unknown mode <%s>.\n",argv[0],modeString);
      Help("usage");
      exit(-1);
    }
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

  /* Print version number */
  if (verbose_flag)
    Message( "# %s\n", rcsid );

  /* Turn on diagnostics if requested, and clear fshrot3d's
   * operation counters.
   */
  intrp_warpClearCounts();
#ifdef never
  if (trilinear_flag) {
    linwarp_set_debug(debug_flag);
    linwarp_clear_counts();
  }
  else if (closest_flag) {
    clwarp_set_debug(debug_flag);
    clwarp_clear_counts();
  }
#endif

  /* Open input dataset */
  if( !strcmp( infile, hdrfile ) )
    Abort( "Input and output files must be distinct." );
  Input = mri_open_dataset( infile, MRI_READ );

  /* Check that program will function on data-set */
  if( !mri_has( Input, "images" ) )
    Abort( "%s operates only on standard images.", argv[0] );
  if( mri_has( Input, "images.dimensions" ) )
    {
      dimstr= mri_get_string( Input, "images.dimensions" );
      if ( !strcmp( dimstr,"vxyzt") || !strcmp(dimstr,"vxyz") ) {
	if ( mri_has( Input, "images.extent.v" ) ) {
	  dv= mri_get_int(Input,"images.extent.v");
	  if( dv != 1 )
	    Abort( "%s takes only reals numbers of the form (v)xyz(t).",
		   argv[0] );
	}
	else
	  Abort( "%s: %s does not have the images.extent.v key.",
		 argv[0], infile );
      }
      else if ( !strcmp( dimstr,"xyzt") || !strcmp(dimstr,"xyz") ) {
	dv= 1;
      }
      else {
	Abort( "%s takes only reals of the form (v)xyzt.",
	       argv[0] );
      }
    }
  else
    Abort( "%s does not have the images.dimensions key.", infile );

  /* Set parameters in local variables */
  if( !mri_has( Input, "images.extent.x" ) ||
      !mri_has( Input, "images.extent.y" ) ||
      !mri_has( Input, "images.extent.z" ) )
    Abort( "images.extent key(s) missing from header." );
  dx = mri_get_int( Input, "images.extent.x" );
  dy = mri_get_int( Input, "images.extent.y" );
  dz = mri_get_int( Input, "images.extent.z" );
  if ( !strcmp(dimstr,"vxyzt") || !strcmp(dimstr,"xyzt") ) {
    if (mri_has(Input, "images.extent.t"))
      dt = mri_get_int( Input, "images.extent.t" );
    else Abort("%s: %s is missing the images.extent.t key.",
	       argv[0],infile);
  }
  else {
    dt= 1;
  }
  if( ( dt <= 0 ) || ( dx <= 0 ) || ( dy <= 0 ) || ( dz <= 0 ) )
    Abort( "images.extent key(s) is non-positive." );

  /* Set output dataset */
  Output = mri_copy_dataset( hdrfile, Input );
  hist_add_cl( Output, argc, argv );
  mri_set_string( Output, "images.datatype", "float32" );
  mri_set_string( Output, "images.dimensions", "vxyzt" );
  mri_set_int( Output, "images.extent.v", (int) dv );
  mri_set_int( Output, "images.extent.t", (int) dt );

  /* Allocate image and parameter storage */
#ifdef never
  if (!(c_image= (FComplex*)malloc(dx*dy*dz*sizeof(FComplex)))) 
    Abort("%s: unable to allocate %d bytes!\n",
	  argv[0],dx*dy*dz*sizeof(FComplex));
  if (!(c_corr_image= (FComplex*)malloc(dx*dy*dz*sizeof(FComplex)))) 
    Abort("%s: unable to allocate %d bytes!\n",
	  argv[0],dx*dy*dz*sizeof(FComplex));
  if ( trilinear_flag || closest_flag ) {
    if (!(check= (char*)malloc(dx*dy*dz*sizeof(char))))
      Abort("%s: unable to allocate %d bytes!\n",
	    argv[0],dx*dy*dz*sizeof(char));
  }
  if (!(corr_image= (float*)malloc(dx*dy*dz*sizeof(float))))
    Abort("%s: unable to allocate %d bytes!\n",
	  argv[0],dx*dy*dz*sizeof(float));
#endif
  if (!(corr_image= (double*)malloc(dx*dy*dz*sizeof(double))))
    Abort("%s: unable to allocate %d bytes!\n",
	  argv[0],dx*dy*dz*sizeof(double));
  if (!(check= (char*)malloc(dx*dy*dz*sizeof(char))))
    Abort("%s: unable to allocate %d bytes!\n",
	  argv[0],dx*dy*dz*sizeof(char));
  if (!(par= (WarpPar*)malloc(dt*sizeof(WarpPar))))
    Abort("%s: unable to allocate %d bytes!\n",
	  dt*sizeof(WarpPar));
  
  /* Initialize warp parameters to something reasonable in case 
   * they are missing 
   */
  for( t = 0; t < dt; t++ ) {
    initializeWarpPar(&(par[t]));
  }

  /* Read in estimated warp parameters */
  load_warp_params( parfile, par, dt );

#ifdef never
  /* Set imaginary part of complex image to zero */
  for(x=0; x<dx; x++)
    for(y=0; y<dy; y++)
      for (z = 0; z<dz; z++)
	MEM(c_image,dx,dy,dz,x,y,z).imag= 0.0;
#endif

  /* Set up voxel edge dimensions. */
  length_x= dx*xvoxel;
  length_y= dy*yvoxel;
  length_z= dz*zvoxel;

  /* Loop through images */
  block_offset= 0;
  block_size= dx*dy*dz;
  interpolator= intrp_createInterpolator3DByType(interpolatorType,
						 dx,dy,dz,1);
  if (debug_flag) {
    interpolator->setInt(interpolator,INTRP_OPT_DEBUG,1);
    interpolator->dumpSelf(interpolator,stderr);
  }
  for( t = 0; t < dt; t++ ) {
	  
    image_in= (double*)mri_get_chunk(Input, "images", block_size,
				     block_offset, MRI_DOUBLE);
    intrp_warpApply( interpolator, par[t].v, image_in, corr_image, check, 
		     dx, dy, dz, 1, length_x, length_y, length_z );
    mri_set_chunk( Output, "images", block_size, block_offset, MRI_DOUBLE,
		   corr_image );
    
    /* Print progress report after completing an image */
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
    
    /* And step to the next block. */
    block_offset += block_size;
  }
      
  /* Write and close data-sets */
  mri_close_dataset( Input );
  mri_close_dataset( Output );
  
  /* Clean up */
  interpolator->destroySelf(interpolator);
#ifdef never
  free(c_corr_image);
  free(c_image);
#endif
  free(corr_image);
  if (check) free(check);
  free(par);

  /* Write out shear counts */
  if (verbose_flag) {
    long count;
    intrp_warpGetCounts(&count);
    Message("# Total transforms by %s interpolation: %ld\n",
	    intrp_nameFromType(interpolatorType),count );
#ifdef never
    if (trilinear_flag) {
      int lin_count;
      linwarp_get_counts(&lin_count);
      Message("# Total transforms by linear interpolation: %d\n",
	      lin_count );
    }
    else if (closest_flag) {
      int count;
      clwarp_get_counts(&count);
      Message("# Total transforms by zeroth-order interpolation: %d\n",
	      count );
    }
    else {
      /* Other warp method counts output would go here */
    }
#endif
    Message( "#      Image warping complete.\n" );
  }

  return 0;

}

