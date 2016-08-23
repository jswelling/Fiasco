/************************************************************
 *                                                          *
 *  matte.c                                                 * 
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
 *  Original programming by Jennifer Bakal June 2002        *
 ************************************************************/
/*************************************************************

  DESCRIPTION OF MATTE.C

  matte.c is used to apply one image file on top of another.

  matte -inmap matte-file infile outfile

**************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"

static char rcsid[] = "$Id: matte.c,v 1.5 2003/07/25 21:03:02 bakalj Exp $";

static char* progname;


#define MINGRAY_DEFAULT 50
#define MAXGRAY_DEFAULT 204
#define RED 0
#define GREEN 1
#define BLUE 2
#define ALPHA 3
#define MAX_COLORS 256
/* #define JENN */

/* Note: ifdef JENN comments out useful debugging print statements
   too long to type in every time */

int
main( argc, argv ) 
     int argc;
     char **argv;
{

  MRI_Dataset *Top_Input = NULL, *Bottom_Input = NULL, *Output = NULL;
  char top_infile[512], bottom_infile[512], outfile[512];
  float lowthresh, highthresh;
  long mingray, maxgray;
  long dv, dx, dy, dz, dt, bottom_dv;
  long top_imgSize, bottom_imgSize;
  float *top = NULL, *bottom = NULL, *matted = NULL;
  unsigned char cmingray, cmaxgray, crange, midgray;
  long x, y, t, z;
  float bottom_min, bottom_max, bottom_range, top_min, top_max, top_range;
  char *m_dimstr, *i_dimstr;
  int i;


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

  /* Get filenames */

  cl_get( "inmap|inm", "%option %s[%]", "map.mri", top_infile );

  if(!cl_get("", "%s", bottom_infile)) {
    fprintf(stderr, "%s: Input file name not given.\n", progname);
    exit(-1);
  }
  if(!cl_get("", "%s", outfile)) {
    fprintf(stderr, "%s: Output file name not given.\n", progname);
    exit(-1);
  }

  if (cl_cleanup_check()) {
    int i;
    fprintf(stderr,"%s: invalid argument in command line:\n    ",progname);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(-1);
  }

  /*** End command-line parsing ***/

  /* Open input datasets */
  if( !strcmp( top_infile, outfile ) || !strcmp( bottom_infile, outfile ) )
    Abort( "Input and output files must be distinct." );
  Top_Input = mri_open_dataset( top_infile, MRI_READ );
  Bottom_Input = mri_open_dataset( bottom_infile, MRI_READ );

  /* Check that program will function on data-set */
  if( !mri_has( Top_Input, "images" ) || !mri_has( Bottom_Input, "images" ) )
    Abort( "%s operates only on standard images.", progname );

  if (!mri_has( Top_Input, "images.dimensions" ) || !mri_has(Bottom_Input, "images.dimensions"))
    Abort("%s: input data missing images.dimensions info\n",progname);

  m_dimstr= mri_get_string(Top_Input, "images.dimensions");
  if (!(!strcmp(m_dimstr,"vxyzt") 
	&& mri_has(Top_Input,"images.extent.v") 
	&& mri_has(Top_Input,"images.extent.t"))
      && !(!strcmp(m_dimstr,"xyzt")
	   && mri_has(Top_Input,"images.extent.t"))
      && !(!strcmp(m_dimstr,"vxyz")
	   && mri_has(Top_Input,"images.extent.v"))
      && !(!strcmp(m_dimstr,"xyz"))) {
    Abort( "%s operates only on real-valued images in order vxyz(t).\n",
	   progname);
  }

  i_dimstr= mri_get_string(Bottom_Input, "images.dimensions");
  if (!(!strcmp(i_dimstr,"vxyzt") 
	&& mri_has(Bottom_Input,"images.extent.v") 
	&& mri_has(Bottom_Input,"images.extent.t"))
      && !(!strcmp(i_dimstr,"xyzt")
	   && mri_has(Bottom_Input,"images.extent.t"))
      && !(!strcmp(i_dimstr,"vxyz")
	   && mri_has(Bottom_Input,"images.extent.v"))
      && !(!strcmp(i_dimstr,"xyz"))) {
    Abort( "%s operates only on real-valued images in order (v)xyz(t).\n",
	   progname);
  }

  if( (*outfile != '.') /* just an extension */
      && (!mri_has( Top_Input, "images.file" ) ||
	  !mri_has( Bottom_Input, "images.file" ) ||
	  !strcmp( outfile, mri_get_string( Top_Input, "images.file" ) ) ||
	  !strcmp( outfile, mri_get_string( Bottom_Input, "images.file" ) ) ) )
    Abort( "%s: Input and output image files are the same (%s).\n", progname, outfile );
  

  /* Set output dataset */
  Output = mri_copy_dataset( outfile, Top_Input );
  hist_add_cl( Output, argc, argv );
  mri_set_string( Output, "images.datatype", "float32" );
  mri_set_string( Output, "images.file", ".dat");


  /* Set parameters in local variables */
  if( !mri_has( Top_Input, "images.extent.v" ) || 
      !mri_has( Top_Input, "images.extent.x" ) ||
      !mri_has( Top_Input, "images.extent.y" ) ||
      !mri_has( Top_Input, "images.extent.z" ) ||
      !mri_has( Bottom_Input, "images.extent.x" ) ||
      !mri_has( Bottom_Input, "images.extent.y" ) ||
      !mri_has( Bottom_Input, "images.extent.z" ) )
    Abort( "images.extent key(s) missing from header." );
  dt = mri_has(Top_Input, "images.extent.t") ?
    mri_get_int( Top_Input, "images.extent.t" ) : 1;
  dv =  mri_has(Top_Input, "images.extent.v") ?
    mri_get_int( Top_Input, "images.extent.v" ) : 1;
  if( dv != 4 )
    Abort( "%s: top image must have v=4 with RGB Alpha info.\n", progname);
  dx = mri_get_int( Top_Input, "images.extent.x" );
  dy = mri_get_int( Top_Input, "images.extent.y" );
  dz = mri_get_int( Top_Input, "images.extent.z" );
  if( ( dt <= 0 ) || ( dx <= 0 ) || ( dy <= 0 ) || ( dz <= 0 ) )
    Abort( "images.extent key(s) is(are) non-positive." );
  if( ((mri_has(Bottom_Input, "images.extent.t") 
	&& mri_get_int( Bottom_Input, "images.extent.t" ) != dt )) ||
      ( mri_get_int( Bottom_Input, "images.extent.x" ) != dx ) ||
      ( mri_get_int( Bottom_Input, "images.extent.y" ) != dy ) ||
      ( mri_get_int( Bottom_Input, "images.extent.z" ) != dz ) )
       Abort( "Dimension extents don't match in %s and %s.\n",
	   top_infile, bottom_infile );
  bottom_dv = mri_has(Bottom_Input, "images.extent.v") ?
    mri_get_int( Bottom_Input, "images.extent.v" ) : 1;
  if((bottom_dv!=1) && (bottom_dv!=3) && (bottom_dv!=4)){
    Abort( "%s: bottom image must have v=1(gray), v=3(RGB), or v=4(RGB Alpha) format.\n", progname);
  }

#ifdef JENN
fprintf (stderr, "map dimensions: v = %d, x = %d, y = %d, z = %d, t = %d; \n input v = %d \n", dv, dx, dy, dz, dt, bottom_dv);
#endif



  /* Allocate image and parameter storage */
  top_imgSize= dx*dy*dv;
  bottom_imgSize= dx*dy*bottom_dv;
  if (!(matted = (float*)malloc(top_imgSize*sizeof(float))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,top_imgSize*sizeof(float));

  /* Loop through images getting ranges */
  
  bottom_min = bottom_max = 0;
  top_min = top_max = 0;

  for( t = 0; t < dt; t++ ) {
    for( z = 0; z < dz; z++ ) {
      
      top = mri_get_chunk(Top_Input, "images", top_imgSize, (t*dz + z)*top_imgSize, MRI_FLOAT);
      bottom = mri_get_chunk(Bottom_Input, "images", bottom_imgSize, (t*dz + z)*bottom_imgSize, MRI_FLOAT);

      for ( y = 0; y < dy; y++) {
	for ( x = 0; x < dx; x++) {

	  /* get bottom range */
	    switch (bottom_dv){
	    case 1:
	      bottom_min = ( bottom[(y*dx)+x] < bottom_min ) ? bottom[(y*dx)+x] : bottom_min;
	      bottom_max = ( bottom[(y*dx)+x] > bottom_max ) ? bottom[(y*dx)+x] : bottom_max;
	      break;
	    case 3:
	    case 4:
	      for (i=0; i<3; i++){
		  bottom_min = ( bottom[(((y*dx)+x)*bottom_dv)+i] < bottom_min ) ? bottom[(((y*dx)+x)*bottom_dv)+i] : bottom_min;
		  bottom_max = ( bottom[(((y*dx)+x)*bottom_dv)+i] > bottom_max ) ? bottom[(((y*dx)+x)*bottom_dv)+i] : bottom_max;
	      }
	      break;
	    }

	  /* get top range */

	    for (i=0; i<3; i++){
	      top_min = ( top[(((y*dx)+x)*dv)+i] < top_min ) ? top[(((y*dx)+x)*dv)+i] : top_min;
	      top_max = ( top[(((y*dx)+x)*dv)+i] > top_max ) ? top[(((y*dx)+x)*dv)+i] : top_max;
	    }
#ifdef JENN
	    fprintf(stderr, "bottom location = %d, i = %d, bottom value = %f, bottom_min = %f, bottom_max = %f, top_min = %f, top_max = %f\n",((y*dx)+x)*bottom_dv,i,bottom[(((y*dx)+x)*bottom_dv)+i],bottom_min, bottom_max);
	    fprintf(stderr, "top location = %d, i = %d, top value = %f, top_min = %f, top_max = %f, top_min = %f, top_max = %f\n",((y*dx)+x)*dv,i,top[(((y*dx)+x)*dv)+i],top_min, top_max);
#endif
	}
      }
    }
  }

  bottom_range = bottom_max - bottom_min;
  top_range = top_max - top_min;
  if (top_range == 0) top_range = bottom_range;  /* otherwise picture becomes all zeros */
  /* this way, picture is essentially bottom_picture, translated somewhat by bottom_min */

#ifdef JENN
  fprintf(stderr, "bottom_min = %f, bottom_max = %f, bottom_range = %f\n", bottom_min, bottom_max, bottom_range);
  fprintf(stderr, "top_min = %f, top_max = %f, top_range = %f\n", top_min, top_max, top_range);
#endif
	  
  /* Loop through images doing matte*/
  for( t = 0; t < dt; t++ ) {
    for( z = 0; z < dz; z++ ) {

      top = mri_get_chunk(Top_Input, "images", top_imgSize, (t*dz + z)*top_imgSize, MRI_FLOAT);
      bottom = mri_get_chunk(Bottom_Input, "images", bottom_imgSize, (t*dz + z)*bottom_imgSize, MRI_FLOAT);

      for ( y = 0; y < dy; y++) {
	for ( x = 0; x < dx; x++) {
	  switch (bottom_dv){
	  case 1:
	    matted[((y*dx)+x)*dv + RED] = (top[((y*dx)+x)*dv + ALPHA])*(top[((y*dx)+x)*dv + RED]) + (1 - top[((y*dx)+x)*dv + ALPHA])*((bottom[((y*dx)+x)]-bottom_min)/bottom_range*top_range);
	    matted[((y*dx)+x)*dv + GREEN] = (top[((y*dx)+x)*dv + ALPHA])*(top[((y*dx)+x)*dv + GREEN]) + (1 - top[((y*dx)+x)*dv + ALPHA])*((bottom[((y*dx)+x)]-bottom_min)/bottom_range*top_range);       
	    matted[((y*dx)+x)*dv + BLUE] = (top[((y*dx)+x)*dv + ALPHA])*(top[((y*dx)+x)*dv + BLUE]) + (1 - top[((y*dx)+x)*dv + ALPHA])*((bottom[((y*dx)+x)]-bottom_min)/bottom_range*top_range);
	    matted[((y*dx)+x)*dv + ALPHA] = (top[((y*dx)+x)*dv + ALPHA])*(top[((y*dx)+x)*dv + ALPHA]) + (1 - top[((y*dx)+x)*dv + ALPHA])*((bottom[((y*dx)+x)]-bottom_min)/bottom_range*top_range);
	    break;
	  case 3:
	    matted[((y*dx)+x)*dv + RED] = (top[((y*dx)+x)*dv + ALPHA])*(top[((y*dx)+x)*dv + RED]) + (1 - top[((y*dx)+x)*dv + ALPHA])*((bottom[((y*dx)+x)*bottom_dv + RED]-bottom_min)/bottom_range*top_range);
	    matted[((y*dx)+x)*dv + GREEN] = (top[((y*dx)+x)*dv + ALPHA])*(top[((y*dx)+x)*dv + GREEN]) + (1 - top[((y*dx)+x)*dv + ALPHA])*((bottom[((y*dx)+x)*bottom_dv + GREEN]-bottom_min)/bottom_range*top_range);
	    matted[((y*dx)+x)*dv + BLUE] = (top[((y*dx)+x)*dv + ALPHA])*(top[((y*dx)+x)*dv + BLUE]) + (1 - top[((y*dx)+x)*dv + ALPHA])*((bottom[((y*dx)+x)*bottom_dv + BLUE]-bottom_min)/bottom_range*top_range);
	    matted[((y*dx)+x)*dv + ALPHA] = 1;
	    break;
	  case 4:
	    matted[((y*dx)+x)*dv + RED] = (top[((y*dx)+x)*dv + ALPHA])*(top[((y*dx)+x)*dv + RED]) + (1 - top[((y*dx)+x)*dv + ALPHA])*((bottom[((y*dx)+x)*bottom_dv + RED]-bottom_min)/bottom_range*top_range);
	    matted[((y*dx)+x)*dv + GREEN] = (top[((y*dx)+x)*dv + ALPHA])*(top[((y*dx)+x)*dv + GREEN]) + (1 - top[((y*dx)+x)*dv + ALPHA])*((bottom[((y*dx)+x)*bottom_dv + GREEN]-bottom_min)/bottom_range*top_range);
	    matted[((y*dx)+x)*dv + BLUE] = (top[((y*dx)+x)*dv + ALPHA])*(top[((y*dx)+x)*dv + BLUE]) + (1 - top[((y*dx)+x)*dv + ALPHA])*((bottom[((y*dx)+x)*bottom_dv + BLUE]-bottom_min)/bottom_range*top_range);
	    matted[((y*dx)+x)*dv + ALPHA] = (top[((y*dx)+x)*dv + ALPHA])*(top[((y*dx)+x)*dv + ALPHA]) + (1 - top[((y*dx)+x)*dv + ALPHA])*((bottom[((y*dx)+x)*bottom_dv + ALPHA]-bottom_min)/bottom_range*top_range);
	    break;
	  }
#ifdef JENN 
fprintf(stderr, "x=%d, y=%d, z=%d, top: red=%f, green=%f, blue=%f, alpha=%f \n", x, y, z, top[((y*dx)+x)*dv + RED], top[((y*dx)+x)*dv + GREEN], top[((y*dx)+x)*dv + BLUE], top[((y*dx)+x)*dv + ALPHA]);
 switch (bottom_dv){
 case 1:
   fprintf(stderr, "x=%d, y=%d, z=%d, bottom: red=%f, green=%f, blue=%f, alpha=%f \n", x, y, z, bottom[((y*dx)+x)], bottom[((y*dx)+x)], bottom[((y*dx)+x)], bottom[((y*dx)+x)]);
   break;
 case 3:
   fprintf(stderr, "x=%d, y=%d, z=%d, bottom: red=%f, green=%f, blue=%f, alpha=%f \n", x, y, z, bottom[((y*dx)+x)*bottom_dv+RED], bottom[((y*dx)+x)*bottom_dv+GREEN], bottom[((y*dx)+x)*bottom_dv+BLUE]);
   break;
 case 4:
   fprintf(stderr, "x=%d, y=%d, z=%d, bottom: red=%f, green=%f, blue=%f, alpha=%f \n", x, y, z, bottom[((y*dx)+x)*bottom_dv+RED], bottom[((y*dx)+x)*bottom_dv+GREEN], bottom[((y*dx)+x)*bottom_dv+BLUE], bottom[((y*dx)+x)*bottom_dv+ALPHA]);
   break;
 }
   fprintf(stderr, "x=%d, y=%d, z=%d, matted: red=%f, green=%f, blue=%f, alpha=%f \n", x, y, z, matted[((y*dx)+x)*dv + RED], matted[((y*dx)+x)*dv + GREEN], matted[((y*dx)+x)*dv + BLUE], matted[((y*dx)+x)*dv + ALPHA]);

#endif
	}
      }


	mri_set_chunk( Output, "images", top_imgSize, (t*dz + z)*top_imgSize, MRI_FLOAT, matted );

    }
  }

  free(matted);


  /* Write and close data-sets */
  mri_close_dataset( Top_Input );
  mri_close_dataset( Bottom_Input );
  mri_close_dataset( Output );
  
  Message( "#      Matte complete.\n" );
  exit(0);

}
