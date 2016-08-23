/************************************************************
 *                                                          *
 *  color_overlay.c                                                 * 
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
 *  Original programming by Jennifer Bakal July 2002        *
 ************************************************************/
/*************************************************************

  DESCRIPTION OF COLORIZE.C

  colorize.c is used to assign colors to a map based on a 
  set of thresholds                       

  colorize.m -color_map color_table infile outfile

**************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"
#include "math.h"

static char rcsid[] = "$Id: colorize.c,v 1.7 2004/09/10 17:30:20 welling Exp $";

static char* progname;



#define RED 0
#define GREEN 1
#define BLUE 2
#define ALPHA 3
#define BOUND 4
#define MAX_COLORS 256

typedef struct color_row_struct {
  float bound;
  float alpha;
  int red;
  int green;
  int blue;
} ColorRow;

  ColorRow Color_Table[MAX_COLORS];


static int load_color_table(float *color_data, int rows) {
  int i;

  for(i=0;i<rows;i++){
    Color_Table[i].red = (int) rint(color_data[5*i+RED]);
    Color_Table[i].green = (int) rint(color_data[5*i+GREEN]);
    Color_Table[i].blue = (int) rint(color_data[5*i+BLUE]);
    Color_Table[i].alpha = color_data[5*i+ALPHA];
    Color_Table[i].bound = color_data[5*i+BOUND];

  }

  return 1;
}

static int get_color_ref(float* t, int rows){
  int i;
  if (Color_Table[0].bound > *t) {
    *t= Color_Table[0].bound;
    return 0;
  }
  for (i=0; i<rows; i++)
    if ((Color_Table[i].bound <= *t) && (*t <= Color_Table[i+1].bound))
      return i;
  *t= Color_Table[rows-1].bound;
  return rows-1;
}

int
main( argc, argv ) 
     int argc;
     char **argv;
{

  MRI_Dataset *CInput = NULL, *Input = NULL, *Output = NULL;
  char color_infile[512], infile[512], outfile[512];
  float lowthresh, highthresh;
  long mingray, maxgray;
  long dv, dx, dy, dz, dt, c_dv, c_dc, o_dv;
  long imgSize, c_imgSize;
  float *map = NULL, *colors = NULL, *colorized = NULL;
  double *image= NULL;
  unsigned char cmingray, cmaxgray, crange, midgray;
  long x, y, t, z;
  float min, max, range;
  char* dimstr;
  float thresh;
  int ref,table;

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

   cl_get( "colormap|col", "%option %s[%]", "color.mri", color_infile );

  if(!cl_get("", "%s", infile)) {
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


  /* Open input dataset */
  if( !strcmp( infile, outfile ) )
    Abort( "Input and output files must be distinct." );


  CInput = mri_open_dataset( color_infile, MRI_READ );
  Input = mri_open_dataset( infile, MRI_READ );

  /* Check that program will function on data-set */
  if( !mri_has( Input, "images" ) )
    Abort( "%s operates only on standard images.", progname );

  if (!mri_has( Input, "images.dimensions" )) {
    Abort("%s: input data missing images.dimensions info.\n",progname);
  }

  if (!mri_has(CInput, "color_table"))
    Abort("%s: color map needs to have a color_table chunk.\n", progname);

  dimstr= mri_get_string(Input, "images.dimensions");
  if (!(!strcmp(dimstr,"vxyzt") 
	&& mri_has(Input,"images.extent.v") 
	&& mri_has(Input,"images.extent.t"))
      && !(!strcmp(dimstr,"xyzt")
	   && mri_has(Input,"images.extent.t"))
      && !(!strcmp(dimstr,"vxyz")
	   && mri_has(Input,"images.extent.v"))
      && !(!strcmp(dimstr,"xyz"))) {
    Abort( "%s operates only on real-valued images in order (v)xyz(t).\n",
	   progname);
  }

  if( (*outfile != '.') /* just an extension */
      && (
	  !mri_has( Input, "images.file" ) ||
	  
	  !strcmp( outfile, mri_get_string( Input, "images.file" ) ) ) )
    Abort( "Input and output image files are the same (%s).", outfile );
  
  /* Set parameters in local variables */
  if( !mri_has( Input, "images.extent.x" ) ||
      !mri_has( Input, "images.extent.y" ) ||
      !mri_has( Input, "images.extent.z" ) )
    Abort( "%s: images.extent key(s) missing from input data.\n", progname);
  if( !mri_has( CInput, "color_table.extent.v" ) || 
      !mri_has( CInput, "color_table.extent.c" ))
    Abort( "%s: color_table.extent key(s) missing from color map.\n",
	   progname);

  dt = mri_has(Input, "images.extent.t") ?
    mri_get_int( Input, "images.extent.t" ) : 1;
  dv = mri_has(Input, "images.extent.v") ?
    mri_get_int( Input, "images.extent.v" ) : 1;
  dx = mri_get_int( Input, "images.extent.x" );
  dy = mri_get_int( Input, "images.extent.y" );
  dz = mri_get_int( Input, "images.extent.z" );
  if( ( dt <= 0 ) || ( dx <= 0 ) || ( dy <= 0 ) || ( dz <= 0 ) )
    Abort( "images.extent key(s) is(are) non-positive." );

  c_dv = mri_get_int(CInput, "color_table.extent.v");
  c_dc = mri_get_int(CInput, "color_table.extent.c");
  if ((c_dv != 5))
    Abort( "%s: %s needs to have 5 columns: RGB alpha threshold. \n", 
	   progname, color_infile);

  colors = mri_get_chunk(CInput, "color_table", c_dv*c_dc, 0, MRI_FLOAT);

  table = load_color_table(colors, c_dc);

  /* Set output dataset */
  Output = mri_copy_dataset( outfile, Input );
  hist_add_cl( Output, argc, argv );
  mri_set_string( Output, "images.dimensions", "vxyzt");
  mri_set_string( Output, "images.datatype", "float32" );
  mri_set_string( Output, "images.file", ".dat");
  o_dv = 4;
  mri_set_int( Output, "images.extent.v", o_dv);
  mri_set_string( Output, "images.description.v", "RGBA" );
  mri_set_int( Output, "images.extent.t", dt);

  /* Create color table chunk in output dataset */
  mri_create_chunk( Output, "color_table" );
  mri_set_string( Output, "color_table.datatype", "float32");
  mri_set_string( Output, "color_table.dimensions", "vc");
  mri_set_int( Output, "color_table.extent.v", c_dv);
  mri_set_int( Output, "color_table.extent.c", c_dc);

  /* fprintf(stderr, "about to set color_table chunk in output\n"); */
  mri_set_chunk( Output, "color_table", c_dv*c_dc, 0, MRI_FLOAT, colors);


	/* Allocate image and parameter storage */
	imgSize= dx*dy*dv;
	c_imgSize = dx*dy*o_dv;

	if (!(colorized = (float*)malloc(c_imgSize*sizeof(float))))
	  Abort("%s: unable to allocate %d bytes!\n",
		progname,c_imgSize*sizeof(float));

  /* Loop through images */
  for( t = 0; t < dt; t++ ) {
    for( z = 0; z < dz; z++ ) 
      {

	image = mri_get_chunk(Input, "images", imgSize, (t*dz + z)*imgSize, MRI_DOUBLE);

	for ( y = 0; y < dy; y++) {
	
	  for ( x = 0; x < dx; x++) {
	    thresh = (float)image[(y*dx)+x];
	    ref = get_color_ref(&thresh, c_dc);
	    if (ref == -1) {
	      colorized[((y*dx)+x)*o_dv + RED] = 0;
	      colorized[((y*dx)+x)*o_dv + GREEN] = 0;
	      colorized[((y*dx)+x)*o_dv + BLUE] = 0;
	      colorized[((y*dx)+x)*o_dv + ALPHA] = 0;
	    }
	    else {
	      float denom= Color_Table[ref+1].bound - Color_Table[ref].bound;
	      if (denom==0.0) {
		/* This condition can occur if two adjacent color table
		 * values have the same bound.  This trick is used to
		 * produce sudden color changes.
		 */
		colorized[((y*dx)+x)*o_dv + RED] = Color_Table[ref].red;
		colorized[((y*dx)+x)*o_dv + GREEN] = Color_Table[ref].green;
		colorized[((y*dx)+x)*o_dv + BLUE] = Color_Table[ref].blue;
		colorized[((y*dx)+x)*o_dv + ALPHA] = Color_Table[ref].alpha;
	      }
	      else {

		colorized[((y*dx)+x)*o_dv + RED] = (((Color_Table[ref+1].bound - thresh) * Color_Table[ref].red) + ((thresh - Color_Table[ref].bound) * Color_Table[ref+1].red))/denom;

		colorized[((y*dx)+x)*o_dv + GREEN] = (((Color_Table[ref+1].bound - thresh) * Color_Table[ref].green) + ((thresh - Color_Table[ref].bound) * Color_Table[ref+1].green))/denom;
		colorized[((y*dx)+x)*o_dv + BLUE] = (((Color_Table[ref+1].bound - thresh) * Color_Table[ref].blue) + ((thresh - Color_Table[ref].bound) * Color_Table[ref+1].blue))/denom;
		/* colorized[((y*dx)+x)*o_dv+ ALPHA] = (((Color_Table[ref].alpha)));*/
		colorized[((y*dx)+x)*o_dv + ALPHA] = (((Color_Table[ref+1].bound - thresh) * Color_Table[ref].alpha) + ((thresh - Color_Table[ref].bound) * Color_Table[ref+1].alpha))/denom;
	      }
	    }

	    /* a useful debugging print statement that takes too long to
	     * type in each time.
	     */

#ifdef never

fprintf(stderr, "x=%d, y=%d, z=%d, thresh=%f, ref=%d, red=%f, green=%f, blue=%f, alpha=%f \n", x, y, z, thresh, ref, colorized[((y*dx)+x)*o_dv + RED], colorized[((y*dx)+x)*o_dv + GREEN], colorized[((y*dx)+x)*o_dv + BLUE], colorized[((y*dx)+x)*o_dv + ALPHA]);

#endif
	  }
	}


	/* fprintf(stderr, "about to set images chunk for slice %d in output.\n", z);*/
	mri_set_chunk( Output, "images", c_imgSize, (t*dz + z)*c_imgSize, MRI_FLOAT, colorized );
      }
  }

  free(colorized);


    /* Write and close data-sets */
  mri_close_dataset( Input );
  mri_close_dataset( CInput );
  mri_close_dataset( Output );
  
  Message( "#      Colorization complete.\n" );
  exit(0);

  }
