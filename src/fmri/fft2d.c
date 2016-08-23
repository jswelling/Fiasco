/************************************************************
 *                                                          *
 *  fft2d.c                                                 *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *                                                          *
 *  Copyright (c) 1999 Department of Statistics             *
 *                     Carnegie Mellon University           *
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
 * fft2d performs a Fourier transform of an image            *
 *   row-wise, column-wise, or both                          *
 *                                                           *
 *   data is the complex-valued image matrix (i/o)           *
 *   nx and ny are the fastest- and slowest-varying          *
 *     dimensions of the matrix data, respectively           *
 *   forward is positive for forward transform and           *
 *     negative for inverse transform                        *
 *   rowcol is '2' for full transform,                       *
 *     'x' for transform across fastest-varying dimension,   *
 *     or 'y' for transform across slowest-varying dimension *
 *   from and to are relevant only to 'x' or 'y' transforms, *
 *     giving the first and last row (column) to be          *
 *     transformed, respectively                             *
 *************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#ifdef FFTW3
#include <fftw3.h>
#else
#include <fftw.h>
#endif
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"
#include "par.h"

/* Notes-
   -y varies fastest within plan_temp.  Use plan_temp[x*ydim+y].
 */

#define FFTW_WISDOM_ENV "F_WISDOM_FILE"

static char rcsid[] = "$Id: fft2d.c,v 1.16 2005/07/30 02:12:55 welling Exp $";

static fftw_plan plan_x_f= NULL;
static fftw_plan plan_x_b= NULL;
static fftw_plan plan_y_f= NULL;
static fftw_plan plan_y_b= NULL;
#ifdef FFTW3
#define FFTW_SUCCESS 1
#define c_re(c) ((c)[0])
#define c_im(c) ((c)[1])
typedef int fftw_direction;
static fftw_complex* plan_temp= NULL;
static fftw_plan plan_2d_f= NULL;
static fftw_plan plan_2d_b= NULL;
#else
static FFTW_COMPLEX* plan_temp= NULL;
static fftwnd_plan plan_2d_f= NULL;
static fftwnd_plan plan_2d_b= NULL;
#endif
static int plan_xdim= 0;
static int plan_ydim= 0;
static int wisdom_loaded= 0;
static int plan_temp_xdim= 0;
static int plan_temp_ydim= 0;

static char* wisdom_fname()
{
  static char* name= NULL;
  char* name_base;
  int instance;
  const char* arch;

  if (!name) {
    if (!(name_base= getenv(FFTW_WISDOM_ENV))) return NULL; /* env not set */
    instance= par_instance();
    arch= par_arch();
    if (instance<0) {
      if (!(name= (char*)malloc(strlen(name_base)+strlen(arch)+16)))
	Abort("fft2d: wisdom_fname: unable to allocate %d bytes!",
	      strlen(name_base)+strlen(arch)+16);
      sprintf(name,"%s.%s",name_base,arch);
    }
    else {
      if (!(name= (char*)malloc(strlen(name_base)+strlen(arch)+64)))
	Abort("fft2d: wisdom_fname: unable to allocate %d bytes!",
	      strlen(name_base)+strlen(arch)+64);
      sprintf(name,"%s.%s.%d",name_base,arch,instance);
    }
  }
  return name;
}

static int check_wisdom()
{
  if (!wisdom_loaded) {
    char* fname;
    FILE* f;

    if (!(fname= wisdom_fname())) return 0; /* can't get file name */
    if (access(fname,R_OK)) return 0; /* can't read file */
    if (!(f= fopen(fname,"r"))) return 0; /* can't open for reading */
    if (fftw_import_wisdom_from_file(f) != FFTW_SUCCESS) {
      (void)fclose(f);
      return 0;
    }
    (void)fclose(f);
    wisdom_loaded= 1;
  }
  return 1;
}

static void save_wisdom()
{
  char* fname;
  FILE* f;

  if (!(fname= wisdom_fname())) return; /* can't get file name */
  if (!(f= fopen(fname,"w"))) return; /* can't open for writing */
  (void)fftw_export_wisdom_to_file(f);
  (void)fclose(f);
}

static int check_temp( int xdim, int ydim )
{
  /* return value non-zero means it was reallocated */
  if ((xdim != plan_temp_xdim) || (ydim != plan_temp_ydim)) {
#ifdef FFTW3
    if (plan_temp) fftw_free(plan_temp);
    if (!(plan_temp= 
	  (fftw_complex*)fftw_malloc(xdim*ydim*sizeof(fftw_complex))))
      Abort("fft2d: check_temp: unable to allocate %d bytes!\n",
	    xdim*ydim*sizeof(fftw_complex));
#else
    if (plan_temp) free(plan_temp);
    if (!(plan_temp= (FFTW_COMPLEX*)malloc(xdim*ydim*sizeof(FFTW_COMPLEX))))
      Abort("fft2d: check_temp: unable to allocate %d bytes!\n",
	    xdim*ydim*sizeof(FFTW_COMPLEX));
#endif
    plan_temp_xdim= xdim;
    plan_temp_ydim= ydim;
    return 1;
  }
  else return 0;
}

static void check_plan_x( int xdim, int ydim, fftw_direction xdir )
{
  if (check_temp(xdim,ydim) 
      || ((xdir==FFTW_FORWARD && !plan_x_f) 
	  || ((xdir==FFTW_BACKWARD && !plan_x_b)))
      || (plan_xdim != xdim) || (plan_ydim != ydim)) {

    (void)check_wisdom();

    if (xdir==FFTW_FORWARD) {
      if (plan_x_f) fftw_destroy_plan(plan_x_f);
#ifdef FFTW3
      plan_x_f= fftw_plan_many_dft( 1, &xdim, ydim,
				    plan_temp, NULL, ydim, 1,
				    plan_temp, NULL, ydim, 1,
				    FFTW_FORWARD, FFTW_MEASURE );
#else
      plan_x_f= fftw_create_plan_specific( xdim, xdir, 
					   FFTW_MEASURE | FFTW_IN_PLACE 
					   | FFTW_USE_WISDOM, 
					   plan_temp, ydim, NULL, 0 );
#endif
    }
    else {
#ifdef FFTW3
      plan_x_b= fftw_plan_many_dft( 1, &xdim, ydim,
				    plan_temp, NULL, ydim, 1,
				    plan_temp, NULL, ydim, 1,
				    FFTW_BACKWARD, FFTW_MEASURE );
#else
      if (plan_x_b) fftw_destroy_plan(plan_x_b);
      plan_x_b= fftw_create_plan_specific( xdim, xdir, 
					   FFTW_MEASURE | FFTW_IN_PLACE 
					   | FFTW_USE_WISDOM, 
					   plan_temp, ydim, NULL, 0 );
#endif
    }

    save_wisdom();
    plan_xdim= xdim;
    plan_ydim= ydim;
  }
}

static void check_plan_y( int xdim, int ydim, fftw_direction ydir )
{
  if (check_temp(xdim,ydim) 
      || ((ydir==FFTW_FORWARD && !plan_y_f) 
	  || ((ydir==FFTW_BACKWARD && !plan_y_b)))
      || (plan_xdim != xdim) || (plan_ydim != ydim)) {

    (void)check_wisdom();

    if (ydir==FFTW_FORWARD) {
      if (plan_y_f) fftw_destroy_plan(plan_y_f);
#ifdef FFTW3
      plan_y_f= fftw_plan_many_dft( 1, &ydim, xdim,
				    plan_temp, NULL, 1, ydim,
				    plan_temp, NULL, 1, ydim,
				    FFTW_FORWARD, FFTW_MEASURE );
#else
      plan_y_f= fftw_create_plan_specific( ydim, ydir, 
					   FFTW_MEASURE | FFTW_IN_PLACE 
					   | FFTW_USE_WISDOM, 
					   plan_temp, 1, NULL, 0 );
#endif
    }
    else {
      if (plan_y_b) fftw_destroy_plan(plan_y_b);
#ifdef FFTW3
      plan_y_b= fftw_plan_many_dft( 1, &ydim, xdim,
				    plan_temp, NULL, 1, ydim,
				    plan_temp, NULL, 1, ydim,
				    FFTW_BACKWARD, FFTW_MEASURE );
#else
      plan_y_b= fftw_create_plan_specific( ydim, ydir, 
					   FFTW_MEASURE | FFTW_IN_PLACE 
					   | FFTW_USE_WISDOM, 
					   plan_temp, 1, NULL, 0 );
#endif
    }

    save_wisdom();
    plan_xdim= xdim;
    plan_ydim= ydim;
  }
}

static void check_plan_2d( int xdim, int ydim, fftw_direction dir )
{
  if (check_temp(xdim,ydim) 
      || ((dir==FFTW_FORWARD && !plan_2d_f) 
	  || ((dir==FFTW_BACKWARD && !plan_2d_b)))
      || (plan_xdim != xdim) || (plan_ydim != ydim)) {

    (void)check_wisdom();

    if (dir==FFTW_FORWARD) {
#ifdef FFTW3
      if (plan_2d_f) fftw_destroy_plan(plan_2d_f);
      plan_2d_f= fftw_plan_dft_2d( xdim, ydim,
				   plan_temp, plan_temp,
				   FFTW_FORWARD, FFTW_MEASURE );
#else
      if (plan_2d_f) fftwnd_destroy_plan(plan_2d_f);
      plan_2d_f= fftw2d_create_plan_specific( xdim, ydim, dir, 
					      FFTW_MEASURE | FFTW_IN_PLACE 
					      | FFTW_USE_WISDOM, 
					      plan_temp, 1, NULL, 0 );
#endif
    }
    else {
#ifdef FFTW3
      if (plan_2d_b) fftw_destroy_plan(plan_2d_b);
      plan_2d_b= fftw_plan_dft_2d( xdim, ydim,
				   plan_temp, plan_temp,
				   FFTW_BACKWARD, FFTW_MEASURE );
#else
      if (plan_2d_b) fftwnd_destroy_plan(plan_2d_b);
      plan_2d_b= fftw2d_create_plan_specific( xdim, ydim, dir, 
					      FFTW_MEASURE | FFTW_IN_PLACE 
					      | FFTW_USE_WISDOM, 
					      plan_temp, 1, NULL, 0 );
#endif
    }

    save_wisdom();
    plan_xdim= xdim;
    plan_ydim= ydim;
  }
}

static void data_to_temp_2d( FComplex** data, long nx, long ny )
{
  long nx_half= nx/2;
  long ny_half= ny/2;
  long nx_bhalf= (nx+1)/2;
  long ny_bhalf= (ny+1)/2;
  long i;
  long j;

  for (i=0; i<nx_bhalf; i++) {
    for (j=0; j<ny_bhalf; j++) {
      c_re(plan_temp[(i*ny)+j])= data[j+ny_half][i+nx_half].real;
      c_im(plan_temp[(i*ny)+j])= data[j+ny_half][i+nx_half].imag;
    }
    for (j=0; j<ny_half; j++) {
      c_re(plan_temp[(i*ny)+j+ny_bhalf])= data[j][i+nx_half].real;
      c_im(plan_temp[(i*ny)+j+ny_bhalf])= data[j][i+nx_half].imag;
    }
  }
  for (i=0; i<nx_half; i++) {
    for (j=0; j<ny_bhalf; j++) {
      c_re(plan_temp[((i+nx_bhalf)*ny)+j])= data[j+ny_half][i].real;
      c_im(plan_temp[((i+nx_bhalf)*ny)+j])= data[j+ny_half][i].imag;
    }
    for (j=0; j<ny_half; j++) {
      c_re(plan_temp[((i+nx_bhalf)*ny)+j+ny_bhalf])= data[j][i].real;
      c_im(plan_temp[((i+nx_bhalf)*ny)+j+ny_bhalf])= data[j][i].imag;
    }
  }
}

static void data_from_temp_2d( FComplex** data, long nx, long ny, 
			       double scale )
{
  long nx_half= nx/2;
  long ny_half= ny/2;
  long nx_bhalf= (nx+1)/2;
  long ny_bhalf= (ny+1)/2;
  long i;
  long j;

  for (i=0; i<nx_bhalf; i++) {
    for (j=0; j<ny_bhalf; j++) {
      data[j+ny_half][i+nx_half].real= c_re(plan_temp[(i*ny)+j]) * scale;
      data[j+ny_half][i+nx_half].imag= c_im(plan_temp[(i*ny)+j]) * scale;
    }
    for (j=0; j<ny_half; j++) {
      data[j][i+nx_half].real= c_re(plan_temp[(i*ny)+j+ny_bhalf]) * scale;
      data[j][i+nx_half].imag= c_im(plan_temp[(i*ny)+j+ny_bhalf]) * scale;
    }
  }
  for (i=0; i<nx_half; i++) {
    for (j=0; j<ny_bhalf; j++) {
      data[j+ny_half][i].real= c_re(plan_temp[((i+nx_bhalf)*ny)+j]) * scale;
      data[j+ny_half][i].imag= c_im(plan_temp[((i+nx_bhalf)*ny)+j]) * scale;
    }
    for (j=0; j<ny_half; j++) {
      data[j][i].real= c_re(plan_temp[((i+nx_bhalf)*ny)+j+ny_bhalf]) * scale;
      data[j][i].imag= c_im(plan_temp[((i+nx_bhalf)*ny)+j+ny_bhalf]) * scale;
    }
  }
}

static void data_to_temp_x( FComplex** data, long nx, long ny )
{
  long nx_half= nx/2;
  long ny_half= ny/2;
  long nx_bhalf= (nx+1)/2;
  long ny_bhalf= (ny+1)/2;
  long i;
  long j;

  for (i=0; i<nx_bhalf; i++) {
    for (j=0; j<ny; j++) {
      c_re(plan_temp[(i*ny)+j])= data[j][i+nx_half].real;
      c_im(plan_temp[(i*ny)+j])= data[j][i+nx_half].imag;
    }
  }
  for (i=0; i<nx_half; i++) {
    for (j=0; j<ny; j++) {
      c_re(plan_temp[((i+nx_bhalf)*ny)+j])= data[j][i].real;
      c_im(plan_temp[((i+nx_bhalf)*ny)+j])= data[j][i].imag;
    }
  }
}

static void data_from_temp_x( FComplex** data, long nx, long ny, double scale )
{
  long nx_half= nx/2;
  long ny_half= ny/2;
  long nx_bhalf= (nx+1)/2;
  long ny_bhalf= (ny+1)/2;
  long i;
  long j;
  
  for (i=0; i<nx_bhalf; i++) {
    for (j=0; j<ny; j++) {
      data[j][i+nx_half].real= c_re(plan_temp[(i*ny)+j]) * scale;
      data[j][i+nx_half].imag= c_im(plan_temp[(i*ny)+j]) * scale;
    }
  }
  for (i=0; i<nx_half; i++) {
    for (j=0; j<ny; j++) {
      data[j][i].real= c_re(plan_temp[((i+nx_bhalf)*ny)+j]) * scale;
      data[j][i].imag= c_im(plan_temp[((i+nx_bhalf)*ny)+j]) * scale;
    }
  }
}

static void data_to_temp_y( FComplex** data, long nx, long ny )
{
  long nx_half= nx/2;
  long ny_half= ny/2;
  long nx_bhalf= (nx+1)/2;
  long ny_bhalf= (ny+1)/2;
  long i;
  long j;

  for (i=0; i<nx; i++) {
    for (j=0; j<ny_bhalf; j++) {
      c_re(plan_temp[(i*ny)+j])= data[j+ny_half][i].real;
      c_im(plan_temp[(i*ny)+j])= data[j+ny_half][i].imag;
    }
    for (j=0; j<ny_half; j++) {
      c_re(plan_temp[(i*ny)+j+ny_bhalf])= data[j][i].real;
      c_im(plan_temp[(i*ny)+j+ny_bhalf])= data[j][i].imag;
    }
  }
}

static void data_from_temp_y( FComplex** data, long nx, long ny, double scale )
{
  long nx_half= nx/2;
  long ny_half= ny/2;
  long nx_bhalf= (nx+1)/2;
  long ny_bhalf= (ny+1)/2;
  long i;
  long j;

  for (i=0; i<nx; i++) {
    for (j=0; j<ny_bhalf; j++) {
      data[j+ny_half][i].real= c_re(plan_temp[(i*ny)+j]) * scale;
      data[j+ny_half][i].imag= c_im(plan_temp[(i*ny)+j]) * scale;
    }
    for (j=0; j<ny_half; j++) {
      data[j][i].real= c_re(plan_temp[(i*ny)+j+ny_bhalf]) * scale;
      data[j][i].imag= c_im(plan_temp[(i*ny)+j+ny_bhalf]) * scale;
    }
  }
}

void fft2d( FComplex** data, long nx, long ny, long sign, 
	    char rowcol, long from, long to )
{
  double scale;
  fftw_direction fftw_dir;

  /* Make sure sign flag is valid */
  if ((sign != -1)&&(sign != 1))
    Abort("fft2d: invalid sign: %d\n",sign);

  if (sign==1) fftw_dir= FFTW_BACKWARD;
  else fftw_dir= FFTW_FORWARD;

  switch( rowcol )
    {

    /* 2D Fourier transform */
    case '2':
      /* Check for plan and temp space */
      check_plan_2d((int)nx,(int)ny,fftw_dir);

      /* Pack everything into temp space */
      data_to_temp_2d(data, nx, ny);

      /* Do the FFT */
#ifdef FFTW3
      if (fftw_dir==FFTW_FORWARD)
	fftw_execute(plan_2d_f);
      else
	fftw_execute(plan_2d_b);
#else
      if (fftw_dir==FFTW_FORWARD)
	fftwnd(plan_2d_f, 1, plan_temp, 1, nx*ny, NULL, 0, 0);
      else
	fftwnd(plan_2d_b, 1, plan_temp, 1, nx*ny, NULL, 0, 0);
#endif
      
      /* Unpack everything from temp space, correcting scale */
      scale = 1.0/sqrt( (float) ( nx * ny ) );
      data_from_temp_2d(data, nx, ny, scale);
      break;
      

    /* 1D Fourier transform across x-direction */
    case 'x':
    case 'X':
      
      /* Check bounds */
      if( ( from > to ) || ( from < 0 ) || ( to > ny ) )
        Abort( "fft2d 'x': From-to parameters out of range (%d,%d).",
               from, to );


      /* Check for plan and temp space */
      check_plan_x((int)nx,(int)ny,fftw_dir);

      /* Pack everything into temp space */
      data_to_temp_x(data, nx, ny);

      /* Do the FFT */
#ifdef FFTW3
      if (fftw_dir==FFTW_FORWARD)
	fftw_execute(plan_x_f);
      else
	fftw_execute(plan_x_b);
#else
      if (fftw_dir==FFTW_FORWARD)
	fftw(plan_x_f, to-from, plan_temp+from, ny, 1, NULL, 0, 0);
      else
	fftw(plan_x_b, to-from, plan_temp+from, ny, 1, NULL, 0, 0);
#endif
      
      /* Unpack everything from temp space, correcting scale */
      scale = 1.0/sqrt( (float)nx );
      data_from_temp_x(data, nx, ny, scale);

      break;
      
      
    /* 1D Fourier transform across y-direction */
    case 'y':
    case 'Y':
      
      /* Check bounds */
      if( ( from > to ) || ( from < 0 ) || ( to > nx ) )
      Abort( "fft2d 'y': From-to parameters out of range (%d,%d).",
             from, to );

      /* Check for plan and temp space */
      check_plan_y((int)nx,(int)ny,fftw_dir);
      /* Pack everything into temp space */
      data_to_temp_y(data, nx, ny);

      /* Do the FFT */
#ifdef FFTW3
      if (fftw_dir==FFTW_FORWARD)
	fftw_execute(plan_y_f);
      else
	fftw_execute(plan_y_b);
#else
      if (fftw_dir==FFTW_FORWARD)
	fftw(plan_y_f, to-from, plan_temp+(from*ny), 1, ny, NULL, 0, 0);
      else
	fftw(plan_y_b, to-from, plan_temp+(from*ny), 1, ny, NULL, 0, 0);
#endif
      
      /* Unpack everything from temp space, correcting scale */
      scale = 1.0/sqrt( (float)ny );
      data_from_temp_y(data, nx, ny, scale);
      
      break;
      
    default:
      Abort( "FFT2D: Unrecognized row-column indicator (%c).", rowcol );
      
    }
  
  return;
}

