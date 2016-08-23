/************************************************************
 *                                                          *
 *  displace.c                                              *
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

  DESCRIPTION OF DISPLACE.C

  displace.c takes a file of the format output by estireg.c 
    and, for each set of registration parameters, calculates
    the mean distance that a pixel is displaced by such 
    a registration

  displace.m [-input Registration-parameter-file]
             [-parameters Output-parameter-file]
             [-xdimension Nx] [-ydimension Ny]

**************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"

static char rcsid[] = "$Id: displace.c,v 1.10 2005/09/03 03:53:15 welling Exp $";

#define KEYBUF_SIZE 512

static char* progname= NULL;

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
  else Abort("%s: input missing tag %s!\n",progname,key_buf);
  return 0; /* not reached */
}

static int weight_valid(MRI_Dataset* wt, 
			const int has_xdim, const int dx,
			const int has_ydim, const int dy)
{
  char key_buf[KEYBUF_SIZE];
  char* dimstr;

  if (!mri_has(wt,"images.dimensions")) {
    Error("%s: weight file has no images.dimensions tag!\n",progname);
    return 0;
  }
  dimstr= mri_get_string(wt,"images.dimensions");
  if (!strncmp(dimstr,"xyz",3)) {
    /* This will definitely work. */
  }
  else if (!strncmp(dimstr,"vxyz",4)) {
    if (safe_get_extent(wt,"images","v") != 1) {
      Abort("%s: weight dataset must have v extent 1!\n",progname);
      return 0;
    }
  }
  else {
    Error("%s: weight dataset must have dimensions (v)xyz(...)!\n",progname);
    return 0;
  }

  if (has_xdim && safe_get_extent(wt,"images","x")!=dx) {
    Error("%s: weight file x extent doesn't match command line value!\n",
	  progname);
    return 0;
  }
  if (has_ydim && safe_get_extent(wt,"images","y")!=dy) {
    Error("%s: weight file y extent doesn't match command line value!\n",
	  progname);
    return 0;
  }

  return 1;
}

int main( int argc, char* argv[] ) 
{
  FILE *ifp = NULL, *ofp = NULL;
  MRI_Dataset* wt= NULL;
  char infile[512], parfile[512], wtfile[512];
  char scanline[512];
  long linenum, numread, dx, dy, dz, x, y, z, z_last, t, t_last;
  long dz_weight, weight_offset;
  RegPars pars;
  float *xloc = NULL, *yloc = NULL;
  float *weights= NULL;
  float radrot, cosrot, sinrot;
  float shiftx, shifty, newx, newy;
  double sum, weight_sum;
  int has_xdim;
  int has_ydim;
  int has_weight;

  progname= argv[0];

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

  /* Deprecate old options */
  if (cl_present( "input|i" ))
     Abort ("Option input|i has been replaced by infile format.  Please see help file.\n");
  if (cl_present( "x" ))
     Abort ("Option x has been expanded to xdm.  Please see help file.\n");
  if (cl_present( "y" ))
     Abort ("Option y has been expanded to ydm.  Please see help file.\n");
  if (cl_present( "parameters|par|p" ))
     Abort ("Option parameters|par|p has been replaced by estimates|est|e.  Please see help file.\n");

  /* Get filenames */
  cl_get( "estimates|est|e", "%option %s[%]", "displace.par", parfile );
  has_xdim= cl_get( "xdimension|xdm", "%option %ld", &dx );
  if (!has_xdim) dx= 128;
  has_ydim= cl_get( "ydimension|ydm", "%option %ld", &dy );
  if (!has_ydim) dy= 64;
  has_weight= cl_get( "weight|wgt|w", "%option %s", wtfile );

  if(!cl_get("", "%s", infile)) {
    fprintf(stderr, "%s: Input file name not given.\n", argv[0]);
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

  if( (has_xdim && ( dx <= 0 )) || (has_ydim && ( dy <= 0 )) )
    Abort( "Dimensions not valid (%ld,%ld).", dx, dy );

  if (!((has_xdim && has_ydim) || has_weight))
    Abort( "Dimension info not provided (explicitly or by weight file).");
  
  if (has_weight) {
    wt= mri_open_dataset(wtfile,MRI_READ);
    if (!weight_valid(wt,has_xdim,dx,has_ydim,dy))
      Abort("%s: weight information is not useable.\n",argv[0]);
    if (!has_xdim) dx= safe_get_extent(wt,"images","x");
    if (!has_ydim) dy= safe_get_extent(wt,"images","y");
    dz_weight= dz= safe_get_extent(wt,"images","z");
  }
  else {
    wt= NULL;
    dz_weight= 1;
    dz= 0;
  }
  
  /* Allocate space */
  if (!(xloc= (float*)malloc( dx*sizeof(float) )))
    Abort("%s: Unable to allocate %d bytes!", argv[0], dx*sizeof(float));
  if (!(yloc= (float*)malloc( dy*sizeof(float) )))
    Abort("%s: Unable to allocate %d bytes!", argv[0], dy*sizeof(float));
  if (!(weights= (float*)malloc( dx*dy*dz_weight*sizeof(float) )))
    Abort("%s: Unable to allocate %d bytes!", argv[0], 
	  dx*dy*dz_weight*sizeof(float));
  
  /* Set fixed location vectors.  Distances are in pixels. */
  for( x = 0; x < dx; x++ )
    xloc[x] = (float) ( x - (long) ( dx / 2 ) );    
  for( y = 0; y < dy; y++ )
    yloc[y] = (float) ( y - (long) ( dy / 2 ) );    
  
  /* Initialize the weights. */
  if (has_weight) {
    mri_read_chunk( wt, "images", dx*dy*dz_weight, 0, MRI_FLOAT, weights );
  }
  else {
    for (x=0; x<dx; x++)
      for (y=0; y<dy; y++) 
	for (z=0; z<dz_weight; z++) {
	  weights[((z*dy)+y)*dx+x]= 1.0;
	}
  }
  
  /* Open output parameter file */
  ofp = efopen( parfile, "w" );
  
  /* Read in estimated registration parameters */
  ifp = efopen( infile, "r" );
  linenum = -1;
  z_last= -1;
  t_last= -1;
  while( !feof( ifp ) && !ferror( ifp ) ) {
    linenum++;
    
    /* Scan a line */
    fscanf( ifp, "%512[^\n]%*[\n]", scanline );
    if (scanline[0]=='#') continue; /* skip comments */
    numread = sscanf( scanline, "%ld%ld%f%f%f%*f",
		      &t,&z,
		      &(pars.x_shift), &(pars.y_shift), &(pars.rotation) );
    if( ( numread < 5 ) ) {
      if( !feof( ifp ) ) {
	Warning( 1,
		 "Line %ld of %s is too short (%ld).  Displace set to -1.0.\n",
		 linenum, infile, numread );
	fprintf( ofp, "%15.9lf\n", -1.0 );
      }
    }
    
    else {
      /* Calculate displacement statistic */
	  
      /* Convert rotation to radians and calculate trig functions */
      radrot = pars.rotation / 180.0 * M_PI;
      cosrot = cos( radrot );
      sinrot = sin( radrot );
      
      /* Loop through image locations */
      sum = 0.0;
      weight_sum = 0.0;
      if (has_weight)
	weight_offset= z*dx*dy;
      else
	weight_offset= 0;
      for( y = 0; y < dy; y++ )
	for( x = 0; x < dx; x++ ) {
	  double xDelta;
	  double yDelta;
	  /* Calculate registered image position */
	  newx = cosrot * (xloc[x]+pars.x_shift) 
	    - sinrot * (yloc[y]+pars.y_shift);
	  newy = sinrot * (xloc[x]+pars.x_shift) 
	    + cosrot * (yloc[y]+pars.y_shift);
	  
	  /* Calculate distance between registered  */
	  /*   and original position and add to sum */
	  xDelta= newx - xloc[x];
	  yDelta= newy - yloc[y];

	  sum += 
	    weights[weight_offset+y*dx+x]*sqrt(xDelta*xDelta + yDelta*yDelta);
	  weight_sum += weights[weight_offset+y*dx+x];
	}
      
      /* Calculate mean displacement and write to parameter file */
      if (weight_sum != 0.0)
	sum /= weight_sum;
      fprintf( ofp, "%15.9lf\n", sum );
	  
    }
  }
  
  /* Close files */
  fclose( ifp );
  fclose( ofp );

  Message( "#      Displacement statistics calculated.\n" );
  exit(0);

}

