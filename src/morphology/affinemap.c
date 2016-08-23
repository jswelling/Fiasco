/************************************************************
 *                                                          *
 *  affinemap.c                                               *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1997 Department of Statistics,         *
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
 *  Original programming by Joel Welling, 1/01              *
 ************************************************************/

#include <stdlib.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "misc.h"
#include "stdcrg.h"

#define KEYBUF_SIZE 512

/* Notes-
 * -output is getting its history from proto; that's inappropriate.
 * -Need to allow for voxel aspect ratio in Fourier method; read it
 *  from the input file?
 */

static char rcsid[] = "$Id: affinemap.c,v 1.2 2003/06/28 00:39:53 welling Exp $";

static MRI_Dataset *Input = NULL, *Output = NULL;
static char selected_dim[512]= "t";
static int offset;
static int range;
static int data_changed= 0;
static float fillval= 0.0;
static char* progname;
static int verbose_flg= 0;
static int fourier_flg= 0;

static double xCtl[3][2]; /* points in output space */
static double xpCtl[3][2]; /* points in input space */
static double A,B,C,D,E,F; /* affine transformation parameters */
static float *in= NULL;
static float *out= NULL;
static FComplex *in_kspace= NULL;
static FComplex *in2_kspace= NULL;
static int dxIn, dyIn, dxOut, dyOut;

/* Linear interpolation, up to 3 steps.  Scale factors are:
 *
 * a is 0.0 on l side of l-r direction
 * b is 0.0 on d side of d-u direction
 * c is 0.0 on f side of f-b direction
 */
#define LIN( l, r, a )  ( ( (a)<=0.0 ) ? (l) : ( ((a)>=1.0) ? (r) : \
                                              ( (l)*(1.0-(a)) + (r)*(a) ) ) )
#define LIN2( ld, rd, lu, ru, a, b ) LIN( LIN(ld,rd,a), LIN(lu,ru,a), b )
#define LIN3( ldf, rdf, luf, ruf, ldb, rdb, lub, rub, a, b, c ) \
  LIN( LIN2( ldf, rdf, luf, ruf, a, b ), LIN2( ldb, rdb, lub, rub, a, b ), c )

static int checkDatasetDims( char* dsname, MRI_Dataset* ds, char* chunk, 
			       char* dims_required )
{
  char* dimstr;
  char buf[256];
  int i;
  int offset;

  if (strlen(chunk)>200) 
    Abort("%s: chunk name <%s> too long!\n",progname,chunk);

  if( !mri_has( ds, chunk ) ) {
    Message( "%s: dataset <%s> has no \"%s\" chunk!", progname, dsname, chunk);
    return 0;
  }

  sprintf(buf,"%s.dimensions",chunk);
  if ( !mri_has( ds, buf ) ) {
    Message( "%s: dataset <%s> has no dimension string for chunk \"%s\"!\n",
	     progname, dsname, chunk );
    return 0;
  }
  else dimstr= mri_get_string( ds, buf );

  offset= 0;
  for (i=0; i<strlen(dimstr); i++) {
    int thisextent;
    char* tchar;
    sprintf(buf,"%s.extent.%c",chunk,dimstr[i]);
    if (!mri_has(ds,buf)) {
      Message("%s: dataset <%s> has no info for %s!\n",buf);
      return 0;
    }
    thisextent= mri_get_int(ds,buf);
    if ((tchar=strchr(dims_required,dimstr[i])) != NULL) {
      /* Make sure this is farther along dims_required than the last one */
      if (tchar-dims_required >= offset) offset= tchar-dims_required;
      else {
	/* Order of strings doesn't match */
	Message( "%s: dataset <%s> dimension order doesn't match \"%s\"!\n",
		 progname, dsname, dims_required );
	return 0;
      }
    }
    else {
      /* not in required string- extent must be 1 */
      if (thisextent != 1) {
	Message( "%s: dataset <%s> has an inappropriate dimension %c!\n",
		 progname, dsname, dimstr[i] );
	return 0;
      }
    }
  }

  /* Were all of the required dims present? */
  if (offset!=strlen(dims_required)-1) {
      Message( "%s: dataset <%s> doesn't include dimensions \"%s\"!\n",
	       progname, dims_required );
      return 0;    
  }

  return 1;
}


static void loadImages()
{
  int i;

  in= mri_get_chunk(Input,"images",dxIn*dyIn, 0, MRI_FLOAT);
  if (!(out=(float*)malloc(dxOut*dyOut*sizeof(float))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,dxOut*dyOut*sizeof(float));
  
  for (i=0; i<dxOut*dyOut; i++) out[i]= 0.0;

  if (fourier_flg) {

    if (!(in_kspace=(FComplex*)malloc(dxIn*dyIn*sizeof(FComplex))))
      Abort("%s: unable to allocate %d bytes!\n",
	    progname,dxIn*dyIn*sizeof(FComplex));

    if (!(in2_kspace=(FComplex*)malloc(dxIn*dyIn*sizeof(FComplex))))
      Abort("%s: unable to allocate %d bytes!\n",
	    progname,dxIn*dyIn*sizeof(FComplex));
  }
}

static void map( double x, double y, double* xp, double* yp )
{
  *xp= A*x + B*y + C;
  *yp= D*x + E*y + F;
}

static void calcCoeffs()
{
  /* These results are from Maple */

  double x1= xCtl[0][0];
  double y1= xCtl[0][1];
  double xp1= xpCtl[0][0];
  double yp1= xpCtl[0][1];
  double x2= xCtl[1][0];
  double y2= xCtl[1][1];
  double xp2= xpCtl[1][0];
  double yp2= xpCtl[1][1];
  double x3= xCtl[2][0];
  double y3= xCtl[2][1];
  double xp3= xpCtl[2][0];
  double yp3= xpCtl[2][1];
  double denom= -x2*y3 + x1*y3 - x1*y2 + x3*y2 - x3*y1 + x2*y1;

  if (denom==0.0) Abort("%s: degenerate transformation!\n",progname);

  B= -(-x1*xp3 + x2*xp3 + x3*xp1 - x2*xp1 + x1*xp2 - x3*xp2)/denom;

  E= (-x2*yp3 + x3*yp2 + x1*yp3 - x1*yp2 - x3*yp1 + x2*yp1)/denom;

  A= (y2*xp3 - y1*xp3 + y1*xp2 - xp2*y3 + xp1*y3 - xp1*y2)/denom;

  D= -(-y2*yp3 + y1*yp3 - y1*yp2 - yp1*y3 + yp1*y2 + yp2*y3)/denom;

  C= -(xp1*x2*y3 - y2*x3*xp1 + y2*x1*xp3 - xp2*x1*y3 - y1*x2*xp3 + y1*x3*xp2)
    /denom;

  F= (yp2*x1*y3 - y1*x3*yp2 + y1*x2*yp3 - yp1*x2*y3 - y2*x1*yp3 + y2*x3*yp1)
    /denom;

  if (verbose_flg) {
    double xt, yt;
    Message("A= %f, B= %f, C= %f, \nD= %f, E= %f, F= %f\n",A,B,C,D,E,F);
    map(x1,y1,&xt,&yt);
    Message("test: (%f %f) vs (%f %f)\n",xp1,yp1,xt,yt);
    map(x2,y2,&xt,&yt);
    Message("      (%f %f) vs (%f %f)\n",xp2,yp2,xt,yt);
    map(x3,y3,&xt,&yt);
    Message("      (%f %f) vs (%f %f)\n",xp3,yp3,xt,yt);
  }
}

static void applyAffineTrilin()
{
  int i, j;
  double x, y;
  double xp, yp;

  calcCoeffs();

  for (j=0; j<dyOut; j++) 
    for (i=0; i<dxOut; i++) {
      x= i;
      y= j;
      map( x, y, &xp, &yp );
      if (xp<0.0 || yp<0.0 || xp>(dxIn-1) || yp>(dyIn-1)) {
	out[ (int)y*dxOut + (int)x ]= 0.0;
      }
      else {
	double floorX= floor(xp);
	double deltaX= xp - floorX;
	double floorY= floor(yp);
	double deltaY= yp - floorY;
	out[ (int)y*dxOut + (int)x ]= 
	  LIN2( in[((int)floorY)*dxIn + (int)floorX],
		in[((int)floorY)*dxIn + (int)floorX + 1],
		in[((int)floorY+1)*dxIn + (int)floorX],
		in[((int)floorY+1)*dxIn + (int)floorX + 1],
		deltaX, deltaY );
      }
    }
}

static void applyAffineFourier()
{
  int i, j;
  double x, y;
  double xp, yp;
  float xChop;
  float yChop;

  calcCoeffs();

  /* Take FFT of input image.  Note the flip in the Y direction,
   * which compensates for the difference between Fiasco 3D physical
   * coordiantes (used in fft3d() and fshrot3d_*()) and the pixel
   * index coordinates used here.
   */
  for (j=0; j<dyIn; j++)
    for (i=0; i<dxIn; i++) {
      in_kspace[(dyIn*i)+(dyIn-(j+1))].real= in[(dxIn*j)+i];
      in_kspace[(dyIn*i)+(dyIn-(j+1))].imag= 0.0;
    }

  fft3d(in_kspace, dxIn, dyIn, 1, +1, "xy");
  
  /* Chop in X and Y to compensate for Fourier center not
   * being in the same place as coordinate origin.
   */
  xChop= yChop= 1.0;
  for (j=0; j<dyIn; j++) {
    for (i=0; i<dxIn; i++) {
      in_kspace[(dxIn*j)+i].real *= xChop*yChop;
      in_kspace[(dxIn*j)+i].imag *= xChop*yChop;      
      xChop *= -1.0;
    }
    yChop *= -1.0;
  }

  for (j=0; j<dyOut; j++) 
    for (i=0; i<dxOut; i++) {
      x= i;
      y= j;
      map( x, y, &xp, &yp );

      /* For each voxel, we phase shift a copy of the input
       * image to move the point of interest to the origin,
       * then Fourier transform that point back to image space
       * to make our sample.  We try to minimize the number of
       * FFTs by not transforming the whole image all the way 
       * back.  Remember that fft3d() expects its inputs in Fortran
       * order, so the X and Y indices have been reversed in the
       * fft3d() calls.
       */
      memcpy(in2_kspace, in_kspace, dxIn*dyIn*sizeof(FComplex));
      fshrot3d_set_shift_phases(-xp, -yp, 0.0, in2_kspace, 
				dxIn, dyIn, 1, 1.0, 1.0, 1.0, 
				1, 1);
      fft3d(in2_kspace, dxIn, dyIn, 1, -1, "xy");

      /* I have no idea why the following alternative code block
       * is slower than FFTing the whole thing, but it apparently is!
       */
#ifdef never
      fft3d(in2_kspace, dxIn, dyIn, 1, -1, "x");
      fft3d(in2_kspace+((dxIn/2)*dyIn), 1, dyIn, 1, -1, "y");
#endif
      /* End of mysteriously bad code block. */

      out[ (j*dxOut) + i ]= in2_kspace[((dxIn/2)*dyIn) + (dyIn/2)].real;
    }
}

static void getCoords()
{
  int which;
  double x, y;

  printf("Give coordinates as floats separated by spaces, one pair per line.\n");
  for (which=0; which<3; which++) {
    int pt1_ok= 0;
    int pt2_ok= 0;
    while (!pt1_ok) {
      if (feof(stdin)||ferror(stdin)) Abort("%s: out of input!\n",progname);
      printf("input point %d>", which+1);
      fflush(stdout); fflush(stdin);
      if (scanf("%lf %lf",&x,&y) == 2) 
	pt1_ok= 1;
      else printf("try again...\n");
    }
    xpCtl[which][0]= x;
    xpCtl[which][1]= y;
    while (!pt2_ok) {
      if (feof(stdin)||ferror(stdin)) Abort("%s: out of input!\n",progname);
      printf("output point %d>", which+1);
      fflush(stdout); fflush(stdin);
      if (scanf("%lf %lf",&x,&y) == 2) pt2_ok= 1;
      else printf("try again...\n");
    }
    xCtl[which][0]= x;
    xCtl[which][1]= y;
  }
  printf("Input complete.\n"); fflush(stdout);
  if (verbose_flg) {
    for (which=0; which<3; which++)
      Message("point %d: (%f %f) -> (%f %f)\n",
	      which+1, xpCtl[which][0], xpCtl[which][1], 
	      xCtl[which][0], xCtl[which][1]);
  }
}

int main(int argc, char** argv ) 
{

  char infile[512], outfile[512], protofile[512];
  char* this_key;
  char this_chunk[KEYBUF_SIZE];
  MRI_Dataset* Proto;
  int dz;
  int dt;

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
  if (!cl_get("h|headerout", "%option [%s]",outfile)) {
    fprintf(stderr,"%s: output file name not given.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  if (!cl_get("i|input", "%option [%s]",infile)) {
    fprintf(stderr,"%s: input file name not given.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  if (!cl_get("p|protofile", "%option [%s]",protofile)) {
    fprintf(stderr,"%s: prototype file name not given.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  verbose_flg= cl_present("v|verbose");
  fourier_flg= cl_present("fourier");
  if (cl_cleanup_check()) {
    int i;
    fprintf(stderr,"%s: invalid argument in command line:\n    ",argv[0]);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(-1);
  }

  /*** End command-line parsing ***/
  
  /* Open input and proto datasets */
  if( !strcmp( infile, outfile ) )
    Abort( "Input and output files must be distinct." );
  Input = mri_open_dataset( infile, MRI_READ );
  Proto = mri_open_dataset( protofile, MRI_READ );

  /* Check the signature of input and output files */
  if (!checkDatasetDims("input",Input,"images","xy")
      || !checkDatasetDims("prototype",Proto,"images","xy")) {
    Abort( "An input dataset was not of the correct type.\n" );
  }

  /* Build the output dataset */
  Output = mri_copy_dataset( outfile, Proto );
  hist_add_cl( Output, argc, argv );
  mri_close_dataset( Proto );

  /* Read coords */
  getCoords();

  /* Load images */
  dxIn= mri_get_int(Input,"images.extent.x");
  dyIn= mri_get_int(Input,"images.extent.y");
  dxOut= mri_get_int(Output,"images.extent.x");
  dyOut= mri_get_int(Output,"images.extent.y");
  loadImages();

  /* Map */
  if (fourier_flg) applyAffineFourier();
  else applyAffineTrilin();

  /* write output */
  mri_set_chunk(Output, "images", dxOut*dyOut, 0, MRI_FLOAT, out);

  /* Write and close data-sets */
  mri_close_dataset( Input );
  mri_close_dataset( Output );
  
  Message( "#      All done.\n" );

  return 0;
}

