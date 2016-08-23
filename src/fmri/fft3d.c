/************************************************************
 *                                                          *
 *  fft3d.c                                                 *
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
 *  Original programming by Joel Welling, 9/98              *
 ************************************************************/

/*************************************************************
 * fft3d performs a Fourier transform of an image            *
 *   in one, two, or all three dimensions.                   *
 *                                                           *
 *   data is the complex-valued image matrix (i/o)           *
 *   nx, ny, and nz are the dimensions of the matrix data,   *
 *     stored such that z varies fastest.                    *
 *   forward is positive for forward transform and           *
 *     negative for inverse transform                        *
 *   rowcol is "3" or "xyz" for full transform,              *
 *     "xy" for transform across x and y,                    *
 *     "yz" for transform across y and z,                    *
 *     "xz" for transform across x and z,                    *
 *     "x" for transform across x,                           *
 *     "y" for transform across y,                           *
 *     "z" for transform across z,                           *
 * The entire dataset is transformed, that is to say,        *
 * if "x" is specified a total of dy*dz transforms in        *
 * the x direction will be performed                         *
 *************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#ifdef FFTW3
#include <fftw3.h>
#else
#include <fftw.h>
#endif
#include <string.h>
#include <strings.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"
#include "par.h"

/* Notes-
   -z varies fastest within plan_temp.  Use plan_temp[(((x*ydim)+y)*zdim)+z].
 */

#define FFTW_WISDOM_ENV "F_WISDOM_FILE"

#define MEM(matrix,nx,ny,nz,x,y,z) matrix[((((x)*ny)+(y))*nz)+(z)]

static char rcsid[] = "$Id: fft3d.c,v 1.14 2007/03/21 23:50:20 welling Exp $";

static int plan_xdim= 0;
static int plan_ydim= 0;
static int plan_zdim= 0;

#ifdef FFTW3
#define FFTW_SUCCESS 1
#define c_re(c) ((c)[0])
#define c_im(c) ((c)[1])
typedef int fftw_direction;
static fftw_plan plan_x_f= NULL;
static fftw_plan plan_x_b= NULL;
static fftw_plan plan_y_f= NULL;
static fftw_plan plan_y_b= NULL;
static fftw_plan plan_z_f= NULL;
static fftw_plan plan_z_b= NULL;
static fftw_plan plan_xy_f= NULL;
static fftw_plan plan_xy_b= NULL;
static fftw_plan plan_yz_f= NULL;
static fftw_plan plan_yz_b= NULL;
static fftw_plan plan_xz_f= NULL;
static fftw_plan plan_xz_b= NULL;
static fftw_plan plan_xyz_f= NULL;
static fftw_plan plan_xyz_b= NULL;
static fftw_complex* plan_temp= NULL;
#else
static fftw_plan plan_x_f= NULL;
static fftw_plan plan_x_b= NULL;
static fftw_plan plan_y_f= NULL;
static fftw_plan plan_y_b= NULL;
static fftw_plan plan_z_f= NULL;
static fftw_plan plan_z_b= NULL;
static fftwnd_plan plan_xy_f= NULL;
static fftwnd_plan plan_xy_b= NULL;
static fftwnd_plan plan_yz_f= NULL;
static fftwnd_plan plan_yz_b= NULL;
static fftwnd_plan plan_xyz_f= NULL;
static fftwnd_plan plan_xyz_b= NULL;
static FFTW_COMPLEX* plan_temp= NULL;
#endif

static int wisdom_loaded= 0;
static int plan_temp_xdim= 0;
static int plan_temp_ydim= 0;
static int plan_temp_zdim= 0;

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
	Abort("fft3d: wisdom_fname: unable to allocate %d bytes!",
	      strlen(name_base)+strlen(arch)+16);
      sprintf(name,"%s.%s",name_base,arch);
    }
    else {
      if (!(name= (char*)malloc(strlen(name_base)+strlen(arch)+64)))
	Abort("fft3d: wisdom_fname: unable to allocate %d bytes!",
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

static int check_temp( int xdim, int ydim, int zdim )
{
  /* return value non-zero means it was reallocated */
  if ((xdim != plan_temp_xdim) || (ydim != plan_temp_ydim)
      || (zdim != plan_temp_zdim)) {
#ifdef FFTW3
    if (plan_temp) fftw_free(plan_temp);
    if (!(plan_temp= 
	  (fftw_complex*)fftw_malloc(xdim*ydim*zdim*sizeof(fftw_complex))))
      Abort("fft3d: check_temp: unable to allocate %d bytes!\n",
	    xdim*ydim*zdim*sizeof(fftw_complex));
#else
    if (plan_temp) free(plan_temp);
    if (!(plan_temp= 
	  (FFTW_COMPLEX*)malloc(xdim*ydim*zdim*sizeof(FFTW_COMPLEX))))
      Abort("fft3d: check_temp: unable to allocate %d bytes!\n",
	    xdim*ydim*zdim*sizeof(FFTW_COMPLEX));
#endif
    plan_temp_xdim= xdim;
    plan_temp_ydim= ydim;
    plan_temp_zdim= zdim;
    return 1;
  }
  else return 0;
}

static void check_plan_x( int xdim, int ydim, int zdim, fftw_direction xdir )
{
  if (check_temp(xdim,ydim,zdim) 
      || ((xdir==FFTW_FORWARD && !plan_x_f) 
	  || ((xdir==FFTW_BACKWARD && !plan_x_b)))
      || (plan_xdim != xdim) || (plan_ydim != ydim) || (plan_zdim != zdim)) {
#ifdef FFTW3
      fftw_iodim dims[1];
      fftw_iodim howmany_dims[2];
      dims[0].n= xdim;
      dims[0].is= dims[0].os= ydim*zdim;
      howmany_dims[0].n= ydim;
      howmany_dims[0].is= howmany_dims[0].os= zdim;
      howmany_dims[1].n= zdim;
      howmany_dims[1].is= howmany_dims[1].os= 1;
#endif

    (void)check_wisdom();

    if (xdir==FFTW_FORWARD) {
      if (plan_x_f) fftw_destroy_plan(plan_x_f);
#ifdef FFTW3
      plan_x_f= fftw_plan_guru_dft( 1, dims, 2, howmany_dims,
				    plan_temp, plan_temp,
				    FFTW_FORWARD, FFTW_MEASURE );
#else
      plan_x_f= fftw_create_plan_specific( xdim, xdir, 
					   FFTW_MEASURE | FFTW_IN_PLACE 
					   | FFTW_USE_WISDOM, 
					   plan_temp, ydim*zdim, NULL, 0 );
#endif
    }
    else {
      if (plan_x_b) fftw_destroy_plan(plan_x_b);
#ifdef FFTW3
      plan_x_b= fftw_plan_guru_dft( 1, dims, 2, howmany_dims,
				    plan_temp, plan_temp,
				    FFTW_BACKWARD, FFTW_MEASURE );
#else
      plan_x_b= fftw_create_plan_specific( xdim, xdir, 
					   FFTW_MEASURE | FFTW_IN_PLACE 
					   | FFTW_USE_WISDOM, 
					   plan_temp, ydim*zdim, NULL, 0 );
#endif
    }

    save_wisdom();
    plan_xdim= xdim;
    plan_ydim= ydim;
    plan_zdim= zdim;
  }
}

static void check_plan_y( int xdim, int ydim, int zdim, fftw_direction ydir )
{
  if (check_temp(xdim,ydim,zdim) 
      || ((ydir==FFTW_FORWARD && !plan_y_f) 
	  || ((ydir==FFTW_BACKWARD && !plan_y_b)))
      || (plan_xdim != xdim) || (plan_ydim != ydim) || (plan_zdim != zdim)) {
#ifdef FFTW3
      fftw_iodim dims[1];
      fftw_iodim howmany_dims[2];
      dims[0].n= ydim;
      dims[0].is= dims[0].os= zdim;
      howmany_dims[0].n= xdim;
      howmany_dims[0].is= howmany_dims[0].os= ydim*zdim;
      howmany_dims[1].n= zdim;
      howmany_dims[1].is= howmany_dims[1].os= 1;
#endif

    (void)check_wisdom();

    if (ydir==FFTW_FORWARD) {
      if (plan_y_f) fftw_destroy_plan(plan_y_f);
#ifdef FFTW3
      plan_y_f= fftw_plan_guru_dft( 1, dims, 2, howmany_dims,
				    plan_temp, plan_temp,
				    FFTW_FORWARD, FFTW_MEASURE );
#else
      plan_y_f= fftw_create_plan_specific( ydim, ydir, 
					   FFTW_MEASURE | FFTW_IN_PLACE 
					   | FFTW_USE_WISDOM, 
					   plan_temp, zdim, NULL, 0 );
#endif
    }
    else {
      if (plan_y_b) fftw_destroy_plan(plan_y_b);

#ifdef FFTW3
      plan_y_b= fftw_plan_guru_dft( 1, dims, 2, howmany_dims,
				    plan_temp, plan_temp,
				    FFTW_BACKWARD, FFTW_MEASURE );
#else
      plan_y_b= fftw_create_plan_specific( ydim, ydir, 
					   FFTW_MEASURE | FFTW_IN_PLACE 
					   | FFTW_USE_WISDOM, 
					   plan_temp, zdim, NULL, 0 );
#endif
    }

    save_wisdom();
    plan_xdim= xdim;
    plan_ydim= ydim;
    plan_zdim= zdim;
  }
}

static void check_plan_z( int xdim, int ydim, int zdim, fftw_direction zdir )
{
  if (check_temp(xdim,ydim,zdim) 
      || ((zdir==FFTW_FORWARD && !plan_z_f) 
	  || ((zdir==FFTW_BACKWARD && !plan_z_b)))
      || (plan_xdim != xdim) || (plan_ydim != ydim) || (plan_zdim != zdim)) {
#ifdef FFTW3
      fftw_iodim dims[1];
      fftw_iodim howmany_dims[2];
      dims[0].n= zdim;
      dims[0].is= dims[0].os= 1;
      howmany_dims[0].n= xdim;
      howmany_dims[0].is= howmany_dims[0].os= ydim*zdim;
      howmany_dims[1].n= ydim;
      howmany_dims[1].is= howmany_dims[1].os= zdim;
#endif

    (void)check_wisdom();

    if (zdir==FFTW_FORWARD) {
      if (plan_z_f) fftw_destroy_plan(plan_z_f);
#ifdef FFTW3
      plan_z_f= fftw_plan_guru_dft( 1, dims, 2, howmany_dims,
				    plan_temp, plan_temp,
				    FFTW_FORWARD, FFTW_MEASURE );
#else
      plan_z_f= fftw_create_plan_specific( zdim, zdir, 
					   FFTW_MEASURE | FFTW_IN_PLACE 
					   | FFTW_USE_WISDOM, 
					   plan_temp, 1, NULL, 0 );
#endif
    }
    else {
      if (plan_z_b) fftw_destroy_plan(plan_z_b);
#ifdef FFTW3
      plan_z_b= fftw_plan_guru_dft( 1, dims, 2, howmany_dims,
				    plan_temp, plan_temp,
				    FFTW_BACKWARD, FFTW_MEASURE );
#else      
      plan_z_b= fftw_create_plan_specific( zdim, zdir, 
					   FFTW_MEASURE | FFTW_IN_PLACE 
					   | FFTW_USE_WISDOM, 
					   plan_temp, 1, NULL, 0 );
#endif
    }

    save_wisdom();
    plan_xdim= xdim;
    plan_ydim= ydim;
    plan_zdim= zdim;
  }
}

static void check_plan_xy( int xdim, int ydim, int zdim, fftw_direction dir )
{
  if (check_temp(xdim,ydim,zdim) 
      || ((dir==FFTW_FORWARD && !plan_xy_f) 
	  || ((dir==FFTW_BACKWARD && !plan_xy_b)))
      || (plan_xdim != xdim) || (plan_ydim != ydim) || (plan_zdim != zdim)) {
#ifdef FFTW3
      fftw_iodim dims[2];
      fftw_iodim howmany_dims[1];
      dims[0].n= xdim;
      dims[0].is= dims[0].os= ydim*zdim;
      dims[1].n= ydim;
      dims[1].is= dims[1].os= zdim;
      howmany_dims[0].n= zdim;
      howmany_dims[0].is= howmany_dims[0].os= 1;
#endif

    (void)check_wisdom();

    if (dir==FFTW_FORWARD) {
#ifdef FFTW3
      if (plan_xy_f) fftw_destroy_plan(plan_xy_f);
      plan_xy_f= fftw_plan_guru_dft( 2, dims, 1, howmany_dims,
				    plan_temp, plan_temp,
				    FFTW_FORWARD, FFTW_MEASURE );
      if (!plan_xy_f) Abort("plan_xy_f came back null!\n");
#else
      if (plan_xy_f) fftwnd_destroy_plan(plan_xy_f);
      plan_xy_f= fftw2d_create_plan_specific( xdim, ydim, dir, 
					      FFTW_MEASURE | FFTW_IN_PLACE 
					      | FFTW_USE_WISDOM, 
					      plan_temp, zdim, NULL, 0 );
#endif
    }
    else {
#ifdef FFTW3
      if (plan_xy_b) fftw_destroy_plan(plan_xy_b);
      plan_xy_b= fftw_plan_guru_dft( 2, dims, 1, howmany_dims,
				    plan_temp, plan_temp,
				    FFTW_BACKWARD, FFTW_MEASURE );
#else
      if (plan_xy_b) fftwnd_destroy_plan(plan_xy_b);
      plan_xy_b= fftw2d_create_plan_specific( xdim, ydim, dir, 
					      FFTW_MEASURE | FFTW_IN_PLACE 
					      | FFTW_USE_WISDOM, 
					      plan_temp, zdim, NULL, 0 );
#endif
    }

    save_wisdom();
    plan_xdim= xdim;
    plan_ydim= ydim;
    plan_zdim= zdim;
  }
}

static void check_plan_yz( int xdim, int ydim, int zdim, fftw_direction dir )
{
  if (check_temp(xdim,ydim,zdim) 
      || ((dir==FFTW_FORWARD && !plan_yz_f) 
	  || ((dir==FFTW_BACKWARD && !plan_yz_b)))
      || (plan_xdim != xdim) || (plan_ydim != ydim) || (plan_zdim != zdim)) {
#ifdef FFTW3
      fftw_iodim dims[2];
      fftw_iodim howmany_dims[1];
      dims[0].n= ydim;
      dims[0].is= dims[0].os= zdim;
      dims[1].n= zdim;
      dims[1].is= dims[1].os= 1;
      howmany_dims[0].n= xdim;
      howmany_dims[0].is= howmany_dims[0].os= ydim*zdim;
#endif

    (void)check_wisdom();

    if (dir==FFTW_FORWARD) {
#ifdef FFTW3
      if (plan_yz_f) fftw_destroy_plan(plan_yz_f);
      plan_yz_f= fftw_plan_guru_dft( 2, dims, 1, howmany_dims,
				    plan_temp, plan_temp,
				    FFTW_FORWARD, FFTW_MEASURE );
#else
      if (plan_yz_f) fftwnd_destroy_plan(plan_yz_f);
      plan_yz_f= fftw2d_create_plan_specific( ydim, zdim, dir, 
					      FFTW_MEASURE | FFTW_IN_PLACE 
					      | FFTW_USE_WISDOM, 
					      plan_temp, 1, NULL, 0 );
#endif
    }
    else {
#ifdef FFTW3
      if (plan_yz_b) fftw_destroy_plan(plan_yz_b);
      plan_yz_b= fftw_plan_guru_dft( 2, dims, 1, howmany_dims,
				    plan_temp, plan_temp,
				    FFTW_BACKWARD, FFTW_MEASURE );
#else
      if (plan_yz_b) fftwnd_destroy_plan(plan_yz_b);
      plan_yz_b= fftw2d_create_plan_specific( ydim, zdim, dir, 
					      FFTW_MEASURE | FFTW_IN_PLACE 
					      | FFTW_USE_WISDOM, 
					      plan_temp, 1, NULL, 0 );
#endif
    }

    save_wisdom();
    plan_xdim= xdim;
    plan_ydim= ydim;
    plan_zdim= zdim;
  }
}

#ifdef FFTW3
static void check_plan_xz( int xdim, int ydim, int zdim, fftw_direction dir )
{
  if (check_temp(xdim,ydim,zdim) 
      || ((dir==FFTW_FORWARD && !plan_yz_f) 
	  || ((dir==FFTW_BACKWARD && !plan_yz_b)))
      || (plan_xdim != xdim) || (plan_ydim != ydim) || (plan_zdim != zdim)) {
      fftw_iodim dims[2];
      fftw_iodim howmany_dims[1];
      dims[0].n= xdim;
      dims[0].is= dims[0].os= ydim*zdim;
      dims[1].n= zdim;
      dims[1].is= dims[1].os= 1;
      howmany_dims[0].n= ydim;
      howmany_dims[0].is= howmany_dims[0].os= zdim;

    (void)check_wisdom();

    if (dir==FFTW_FORWARD) {
      if (plan_xz_f) fftw_destroy_plan(plan_xz_f);
      plan_xz_f= fftw_plan_guru_dft( 2, dims, 1, howmany_dims,
				    plan_temp, plan_temp,
				    FFTW_FORWARD, FFTW_MEASURE );
    }
    else {
      if (plan_xz_b) fftw_destroy_plan(plan_xz_b);
      plan_xz_b= fftw_plan_guru_dft( 2, dims, 1, howmany_dims,
				    plan_temp, plan_temp,
				    FFTW_BACKWARD, FFTW_MEASURE );
    }

    save_wisdom();
    plan_xdim= xdim;
    plan_ydim= ydim;
    plan_zdim= zdim;
  }
}
#endif

static void check_plan_xyz( int xdim, int ydim, int zdim, fftw_direction dir )
{
  if (check_temp(xdim,ydim,zdim) 
      || ((dir==FFTW_FORWARD && !plan_xyz_f) 
	  || ((dir==FFTW_BACKWARD && !plan_xyz_b)))
      || (plan_xdim != xdim) || (plan_ydim != ydim) || (plan_zdim != zdim)) {
#ifdef FFTW3
      fftw_iodim dims[3];
      dims[0].n= xdim;
      dims[0].is= dims[0].os= ydim*zdim;
      dims[1].n= ydim;
      dims[1].is= dims[1].os= zdim;
      dims[2].n= zdim;
      dims[2].is= dims[2].os= 1;
#endif

    (void)check_wisdom();

    if (dir==FFTW_FORWARD) {
#ifdef FFTW3
      if (plan_xyz_f) fftw_destroy_plan(plan_xyz_f); 
      plan_xyz_f= fftw_plan_guru_dft( 3, dims, 0, NULL,
				      plan_temp, plan_temp,
				      FFTW_FORWARD, FFTW_MEASURE );
#else
      if (plan_xyz_f) fftwnd_destroy_plan(plan_xyz_f);      
      plan_xyz_f= fftw3d_create_plan_specific( xdim, ydim, zdim, dir, 
					       FFTW_MEASURE | FFTW_IN_PLACE 
					       | FFTW_USE_WISDOM, 
					       plan_temp, 1, NULL, 0 );
#endif
    }
    else {
#ifdef FFTW3
      if (plan_xyz_b) fftw_destroy_plan(plan_xyz_b);      
      plan_xyz_b= fftw_plan_guru_dft( 3, dims, 0, NULL,
				      plan_temp, plan_temp,
				      FFTW_BACKWARD, FFTW_MEASURE );
#else
      if (plan_xyz_b) fftwnd_destroy_plan(plan_xyz_b);
      plan_xyz_b= fftw3d_create_plan_specific( xdim, ydim, zdim, dir, 
					       FFTW_MEASURE | FFTW_IN_PLACE 
					       | FFTW_USE_WISDOM, 
					       plan_temp, 1, NULL, 0 );
#endif
    }

    save_wisdom();
    plan_xdim= xdim;
    plan_ydim= ydim;
    plan_zdim= zdim;
  }
}

static void data_to_temp_x( FComplex* data, long nx, long ny, long nz )
{
  long nx_half= nx/2;
  long nx_bhalf= (nx+1)/2;
  long i;
  long j;
  long k;
  
  for (i=0; i<nx_bhalf; i++) {
    for (j=0; j<ny; j++) {
      for (k=0; k<nz; k++) {
	c_re(MEM(plan_temp,nx,ny,nz,i,j,k))= 
	  MEM(data,nx,ny,nz,i+nx_half,j,k).real;
	c_im(MEM(plan_temp,nx,ny,nz,i,j,k))= 
	  MEM(data,nx,ny,nz,i+nx_half,j,k).imag;
      }
    }
  }
  for (i=0; i<nx_half; i++) {
    for (j=0; j<ny; j++) {
      for (k=0; k<nz; k++) {
	c_re(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j,k))= 
	  MEM(data,nx,ny,nz,i,j,k).real;
	c_im(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j,k))= 
	  MEM(data,nx,ny,nz,i,j,k).imag;
      }
    }
  }
}

static void data_from_temp_x( FComplex* data, long nx, long ny, long nz, 
			      double scale )
{
  long nx_half= nx/2;
  long nx_bhalf= (nx+1)/2;
  long i;
  long j;
  long k;
  
  for (i=0; i<nx_bhalf; i++) {
    for (j=0; j<ny; j++) {
      for (k=0; k<nz; k++) {
	MEM(data,nx,ny,nz,i+nx_half,j,k).real=
	  c_re(MEM(plan_temp,nx,ny,nz,i,j,k))*scale;
	MEM(data,nx,ny,nz,i+nx_half,j,k).imag=
	  c_im(MEM(plan_temp,nx,ny,nz,i,j,k))*scale;
      }
    }
  }
  for (i=0; i<nx_half; i++) {
    for (j=0; j<ny; j++) {
      for (k=0; k<nz; k++) {
	MEM(data,nx,ny,nz,i,j,k).real=
	  c_re(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j,k))*scale;
	MEM(data,nx,ny,nz,i,j,k).imag=
	  c_im(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j,k))*scale;
      }
    }
  }
}

static void data_to_temp_y( FComplex* data, long nx, long ny, long nz )
{
  long ny_half= ny/2;
  long ny_bhalf= (ny+1)/2;
  long i;
  long j;
  long k;
  
  for (i=0; i<nx; i++) {
    for (j=0; j<ny_bhalf; j++) {
      for (k=0; k<nz; k++) {
	c_re(MEM(plan_temp,nx,ny,nz,i,j,k))= 
	  MEM(data,nx,ny,nz,i,j+ny_half,k).real;
	c_im(MEM(plan_temp,nx,ny,nz,i,j,k))= 
	  MEM(data,nx,ny,nz,i,j+ny_half,k).imag;
      }
    }
    for (j=0; j<ny_half; j++) {
      for (k=0; k<nz; k++) {
	c_re(MEM(plan_temp,nx,ny,nz,i,j+ny_bhalf,k))= 
	  MEM(data,nx,ny,nz,i,j,k).real;
	c_im(MEM(plan_temp,nx,ny,nz,i,j+ny_bhalf,k))= 
	  MEM(data,nx,ny,nz,i,j,k).imag;
      }
    }
  }
}

static void data_from_temp_y( FComplex* data, long nx, long ny, long nz,
			      double scale )
{
  long ny_half= ny/2;
  long ny_bhalf= (ny+1)/2;
  long i;
  long j;
  long k;
  
  for (i=0; i<nx; i++) {
    for (j=0; j<ny_bhalf; j++) {
      for (k=0; k<nz; k++) {
	MEM(data,nx,ny,nz,i,j+ny_half,k).real=
	  c_re(MEM(plan_temp,nx,ny,nz,i,j,k))*scale;
	MEM(data,nx,ny,nz,i,j+ny_half,k).imag=
	  c_im(MEM(plan_temp,nx,ny,nz,i,j,k))*scale;
      }
    }
    for (j=0; j<ny_half; j++) {
      for (k=0; k<nz; k++) {
	MEM(data,nx,ny,nz,i,j,k).real=
	  c_re(MEM(plan_temp,nx,ny,nz,i,j+ny_bhalf,k))*scale;
	MEM(data,nx,ny,nz,i,j,k).imag=
	  c_im(MEM(plan_temp,nx,ny,nz,i,j+ny_bhalf,k))*scale;
      }
    }
  }
}

static void data_to_temp_z( FComplex* data, long nx, long ny, long nz )
{
  long nz_half= nz/2;
  long nz_bhalf= (nz+1)/2;
  long i;
  long j;
  long k;
  
  for (i=0; i<nx; i++) {
    for (j=0; j<ny; j++) {
      for (k=0; k<nz_bhalf; k++) {
	c_re(MEM(plan_temp,nx,ny,nz,i,j,k))=
	  MEM(data,nx,ny,nz,i,j,k+nz_half).real;
	c_im(MEM(plan_temp,nx,ny,nz,i,j,k))=
	  MEM(data,nx,ny,nz,i,j,k+nz_half).imag;
      }
      for (k=0; k<nz_half; k++) {
	c_re(MEM(plan_temp,nx,ny,nz,i,j,k+nz_bhalf))=
	  MEM(data,nx,ny,nz,i,j,k).real;
	c_im(MEM(plan_temp,nx,ny,nz,i,j,k+nz_bhalf))=
	  MEM(data,nx,ny,nz,i,j,k).imag;
      }
    }
  }
}

static void data_from_temp_z( FComplex* data, long nx, long ny, long nz, 
			      double scale )
{
  long nz_half= nz/2;
  long nz_bhalf= (nz+1)/2;
  long i;
  long j;
  long k;
  
  for (i=0; i<nx; i++) {
    for (j=0; j<ny; j++) {
      for (k=0; k<nz_bhalf; k++) {
	MEM(data,nx,ny,nz,i,j,k+nz_half).real=
	  c_re(MEM(plan_temp,nx,ny,nz,i,j,k))*scale;
	MEM(data,nx,ny,nz,i,j,k+nz_half).imag=
	  c_im(MEM(plan_temp,nx,ny,nz,i,j,k))*scale;
      }
      for (k=0; k<nz_half; k++) {
	MEM(data,nx,ny,nz,i,j,k).real=
	  c_re(MEM(plan_temp,nx,ny,nz,i,j,k+nz_bhalf))*scale;
	MEM(data,nx,ny,nz,i,j,k).imag=
	  c_im(MEM(plan_temp,nx,ny,nz,i,j,k+nz_bhalf))*scale;
      }
    }
  }
}

static void data_to_temp_xy( FComplex* data, long nx, long ny, long nz )
{
  long nx_half= nx/2;
  long nx_bhalf= (nx+1)/2;
  long ny_half= ny/2;
  long ny_bhalf= (ny+1)/2;
  long i;
  long j;
  long k;
  
  for (i=0; i<nx_bhalf; i++) {
    for (j=0; j<ny_bhalf; j++)
      for (k=0; k<nz; k++) {
	c_re(MEM(plan_temp,nx,ny,nz,i,j,k))=
	  MEM(data,nx,ny,nz,i+nx_half,j+ny_half,k).real;
	c_im(MEM(plan_temp,nx,ny,nz,i,j,k))=
	  MEM(data,nx,ny,nz,i+nx_half,j+ny_half,k).imag;
      }
    for (j=0; j<ny_half; j++)
      for (k=0; k<nz; k++) {
	c_re(MEM(plan_temp,nx,ny,nz,i,j+ny_bhalf,k))=
	  MEM(data,nx,ny,nz,i+nx_half,j,k).real;
	c_im(MEM(plan_temp,nx,ny,nz,i,j+ny_bhalf,k))=
	  MEM(data,nx,ny,nz,i+nx_half,j,k).imag;
      }
  }
  for (i=0; i<nx_half; i++) {
    for (j=0; j<ny_bhalf; j++)
      for (k=0; k<nz; k++) {
	c_re(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j,k))=
	  MEM(data,nx,ny,nz,i,j+ny_half,k).real;
	c_im(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j,k))=
	  MEM(data,nx,ny,nz,i,j+ny_half,k).imag;
      }
    for (j=0; j<ny_half; j++)
      for (k=0; k<nz; k++) {
	c_re(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j+ny_bhalf,k))=
	  MEM(data,nx,ny,nz,i,j,k).real;
	c_im(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j+ny_bhalf,k))=
	  MEM(data,nx,ny,nz,i,j,k).imag;
      }
  }
}

static void data_from_temp_xy( FComplex* data, long nx, long ny, long nz,
			       double scale )
{
  long nx_half= nx/2;
  long nx_bhalf= (nx+1)/2;
  long ny_half= ny/2;
  long ny_bhalf= (ny+1)/2;
  long i;
  long j;
  long k;
  
  for (i=0; i<nx_bhalf; i++) {
    for (j=0; j<ny_bhalf; j++)
      for (k=0; k<nz; k++) {
	MEM(data,nx,ny,nz,i+nx_half,j+ny_half,k).real=
	  c_re(MEM(plan_temp,nx,ny,nz,i,j,k))*scale;
	MEM(data,nx,ny,nz,i+nx_half,j+ny_half,k).imag=
	  c_im(MEM(plan_temp,nx,ny,nz,i,j,k))*scale;
      }
    for (j=0; j<ny_half; j++)
      for (k=0; k<nz; k++) {
	MEM(data,nx,ny,nz,i+nx_half,j,k).real=
	  c_re(MEM(plan_temp,nx,ny,nz,i,j+ny_bhalf,k))*scale;
	MEM(data,nx,ny,nz,i+nx_half,j,k).imag=
	  c_im(MEM(plan_temp,nx,ny,nz,i,j+ny_bhalf,k))*scale;
      }
  }
  for (i=0; i<nx_half; i++) {
    for (j=0; j<ny_bhalf; j++)
      for (k=0; k<nz; k++) {
	MEM(data,nx,ny,nz,i,j+ny_half,k).real=
	  c_re(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j,k))*scale;
	MEM(data,nx,ny,nz,i,j+ny_half,k).imag=
	  c_im(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j,k))*scale;
      }
    for (j=0; j<ny_half; j++)
      for (k=0; k<nz; k++) {
	MEM(data,nx,ny,nz,i,j,k).real=
	  c_re(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j+ny_bhalf,k))*scale;
	MEM(data,nx,ny,nz,i,j,k).imag=
	  c_im(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j+ny_bhalf,k))*scale;
      }
  }
}

static void data_to_temp_yz( FComplex* data, long nx, long ny, long nz )
{
  long ny_half= ny/2;
  long ny_bhalf= (ny+1)/2;
  long nz_half= nz/2;
  long nz_bhalf= (nz+1)/2;
  long i;
  long j;
  long k;
  
  for (i=0; i<nx; i++) {
    for (j=0; j<ny_bhalf; j++) {
      for (k=0; k<nz_bhalf; k++) {
	c_re(MEM(plan_temp,nx,ny,nz,i,j,k))=
	  MEM(data,nx,ny,nz,i,j+ny_half,k+nz_half).real;
	c_im(MEM(plan_temp,nx,ny,nz,i,j,k))=
	  MEM(data,nx,ny,nz,i,j+ny_half,k+nz_half).imag;
      }
      for (k=0; k<nz_half; k++) {
	c_re(MEM(plan_temp,nx,ny,nz,i,j,k+nz_bhalf))=
	  MEM(data,nx,ny,nz,i,j+ny_half,k).real;
	c_im(MEM(plan_temp,nx,ny,nz,i,j,k+nz_bhalf))=
	  MEM(data,nx,ny,nz,i,j+ny_half,k).imag;
      }
    }
    for (j=0; j<ny_half; j++) {
      for (k=0; k<nz_bhalf; k++) {
	c_re(MEM(plan_temp,nx,ny,nz,i,j+ny_bhalf,k))=
	  MEM(data,nx,ny,nz,i,j,k+nz_half).real;
	c_im(MEM(plan_temp,nx,ny,nz,i,j+ny_bhalf,k))=
	  MEM(data,nx,ny,nz,i,j,k+nz_half).imag;
      }
      for (k=0; k<nz_half; k++) {
	c_re(MEM(plan_temp,nx,ny,nz,i,j+ny_bhalf,k+nz_bhalf))=
	  MEM(data,nx,ny,nz,i,j,k).real;
	c_im(MEM(plan_temp,nx,ny,nz,i,j+ny_bhalf,k+nz_bhalf))=
	  MEM(data,nx,ny,nz,i,j,k).imag;
      }
    }
  }
}

static void data_from_temp_yz( FComplex* data, long nx, long ny, long nz,
			       double scale )
{
  long ny_half= ny/2;
  long ny_bhalf= (ny+1)/2;
  long nz_half= nz/2;
  long nz_bhalf= (nz+1)/2;
  long i;
  long j;
  long k;
  
  for (i=0; i<nx; i++) {
    for (j=0; j<ny_bhalf; j++) {
      for (k=0; k<nz_bhalf; k++) {
	MEM(data,nx,ny,nz,i,j+ny_half,k+nz_half).real=
	  c_re(MEM(plan_temp,nx,ny,nz,i,j,k))*scale;
	MEM(data,nx,ny,nz,i,j+ny_half,k+nz_half).imag=
	  c_im(MEM(plan_temp,nx,ny,nz,i,j,k))*scale;
      }
      for (k=0; k<nz_half; k++) {
	MEM(data,nx,ny,nz,i,j+ny_half,k).real=
	  c_re(MEM(plan_temp,nx,ny,nz,i,j,k+nz_bhalf))*scale;
	MEM(data,nx,ny,nz,i,j+ny_half,k).imag=
	  c_im(MEM(plan_temp,nx,ny,nz,i,j,k+nz_bhalf))*scale;
      }
    }
    for (j=0; j<ny_half; j++) {
      for (k=0; k<nz_bhalf; k++) {
	MEM(data,nx,ny,nz,i,j,k+nz_half).real=
	  c_re(MEM(plan_temp,nx,ny,nz,i,j+ny_bhalf,k))*scale;
	MEM(data,nx,ny,nz,i,j,k+nz_half).imag=
	  c_im(MEM(plan_temp,nx,ny,nz,i,j+ny_bhalf,k))*scale;
      }
      for (k=0; k<nz_half; k++) {
	MEM(data,nx,ny,nz,i,j,k).real=
	  c_re(MEM(plan_temp,nx,ny,nz,i,j+ny_bhalf,k+nz_bhalf))*scale;
	MEM(data,nx,ny,nz,i,j,k).imag=
	  c_im(MEM(plan_temp,nx,ny,nz,i,j+ny_bhalf,k+nz_bhalf))*scale;
      }
    }
  }
}

static void data_to_temp_xz( FComplex* data, long nx, long ny, long nz )
{
  long nx_half= nx/2;
  long nx_bhalf= (nx+1)/2;
  long nz_half= nz/2;
  long nz_bhalf= (nz+1)/2;
  long i;
  long j;
  long k;

  for (i=0; i<nx_bhalf; i++) 
    for (j=0; j<ny; j++) {
      for (k=0; k<nz_bhalf; k++) {
	c_re(MEM(plan_temp,nx,ny,nz,i,j,k))=
	  MEM(data,nx,ny,nz,i+nx_half,j,k+nz_half).real;
	c_im(MEM(plan_temp,nx,ny,nz,i,j,k))=
	  MEM(data,nx,ny,nz,i+nx_half,j,k+nz_half).imag;
      }
      for (k=0; k<nz_half; k++) {
	c_re(MEM(plan_temp,nx,ny,nz,i,j,k+nz_bhalf))=
	  MEM(data,nx,ny,nz,i+nx_half,j,k).real;
	c_im(MEM(plan_temp,nx,ny,nz,i,j,k+nz_bhalf))=
	  MEM(data,nx,ny,nz,i+nx_half,j,k).imag;
      }
    }
  for (i=0; i<nx_half; i++) 
    for (j=0; j<ny; j++) {
      for (k=0; k<nz_bhalf; k++) {
	c_re(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j,k))=
	  MEM(data,nx,ny,nz,i,j,k+nz_half).real;
	c_im(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j,k))=
	  MEM(data,nx,ny,nz,i,j,k+nz_half).imag;
      }
      for (k=0; k<nz_half; k++) {
	c_re(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j,k+nz_bhalf))=
	  MEM(data,nx,ny,nz,i,j,k).real;
	c_im(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j,k+nz_bhalf))=
	  MEM(data,nx,ny,nz,i,j,k).imag;
      }
    }
}

static void data_from_temp_xz( FComplex* data, long nx, long ny, long nz,
			       double scale )
{
  long nx_half= nx/2;
  long nx_bhalf= (nx+1)/2;
  long nz_half= nz/2;
  long nz_bhalf= (nz+1)/2;
  long i;
  long j;
  long k;

  for (i=0; i<nx_bhalf; i++) 
    for (j=0; j<ny; j++) {
      for (k=0; k<nz_bhalf; k++) {
	MEM(data,nx,ny,nz,i+nx_half,j,k+nz_half).real=
	  c_re(MEM(plan_temp,nx,ny,nz,i,j,k))*scale;
	MEM(data,nx,ny,nz,i+nx_half,j,k+nz_half).imag=
	  c_im(MEM(plan_temp,nx,ny,nz,i,j,k))*scale;
      }
      for (k=0; k<nz_half; k++) {
	MEM(data,nx,ny,nz,i+nx_half,j,k).real=
	  c_re(MEM(plan_temp,nx,ny,nz,i,j,k+nz_bhalf))*scale;
	MEM(data,nx,ny,nz,i+nx_half,j,k).imag=
	  c_im(MEM(plan_temp,nx,ny,nz,i,j,k+nz_bhalf))*scale;
      }
    }
  for (i=0; i<nx_half; i++) 
    for (j=0; j<ny; j++) {
      for (k=0; k<nz_bhalf; k++) {
	MEM(data,nx,ny,nz,i,j,k+nz_half).real=
	  c_re(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j,k))*scale;
	MEM(data,nx,ny,nz,i,j,k+nz_half).imag=
	  c_im(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j,k))*scale;
      }
      for (k=0; k<nz_half; k++) {
	MEM(data,nx,ny,nz,i,j,k).real=
	  c_re(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j,k+nz_bhalf))*scale;
	MEM(data,nx,ny,nz,i,j,k).imag=
	  c_im(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j,k+nz_bhalf))*scale;
      }
    }
}

static void data_to_temp_xyz( FComplex* data, long nx, long ny, long nz )
{
  long nx_half= nx/2;
  long nx_bhalf= (nx+1)/2;
  long ny_half= ny/2;
  long ny_bhalf= (ny+1)/2;
  long nz_half= nz/2;
  long nz_bhalf= (nz+1)/2;
  long i;
  long j;
  long k;

  for (i=0; i<nx_bhalf; i++) {
    for (j=0; j<ny_bhalf; j++) {
      for (k=0; k<nz_bhalf; k++) {
	c_re(MEM(plan_temp,nx,ny,nz,i,j,k))=
	  MEM(data,nx,ny,nz,i+nx_half,j+ny_half,k+nz_half).real;
	c_im(MEM(plan_temp,nx,ny,nz,i,j,k))=
	  MEM(data,nx,ny,nz,i+nx_half,j+ny_half,k+nz_half).imag;
      }
      for (k=0; k<nz_half; k++) {
	c_re(MEM(plan_temp,nx,ny,nz,i,j,k+nz_bhalf))=
	  MEM(data,nx,ny,nz,i+nx_half,j+ny_half,k).real;
	c_im(MEM(plan_temp,nx,ny,nz,i,j,k+nz_bhalf))=
	  MEM(data,nx,ny,nz,i+nx_half,j+ny_half,k).imag;
      }
    }
    for (j=0; j<ny_half; j++) {
      for (k=0; k<nz_bhalf; k++) {
	c_re(MEM(plan_temp,nx,ny,nz,i,j+ny_bhalf,k))=
	  MEM(data,nx,ny,nz,i+nx_half,j,k+nz_half).real;
	c_im(MEM(plan_temp,nx,ny,nz,i,j+ny_bhalf,k))=
	  MEM(data,nx,ny,nz,i+nx_half,j,k+nz_half).imag;
      }
      for (k=0; k<nz_half; k++) {
	c_re(MEM(plan_temp,nx,ny,nz,i,j+ny_bhalf,k+nz_bhalf))=
	  MEM(data,nx,ny,nz,i+nx_half,j,k).real;
	c_im(MEM(plan_temp,nx,ny,nz,i,j+ny_bhalf,k+nz_bhalf))=
	  MEM(data,nx,ny,nz,i+nx_half,j,k).imag;
      }
    }
  }
  for (i=0; i<nx_half; i++) {
    for (j=0; j<ny_bhalf; j++) {
      for (k=0; k<nz_bhalf; k++) {
	c_re(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j,k))=
	  MEM(data,nx,ny,nz,i,j+ny_half,k+nz_half).real;
	c_im(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j,k))=
	  MEM(data,nx,ny,nz,i,j+ny_half,k+nz_half).imag;
      }
      for (k=0; k<nz_half; k++) {
	c_re(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j,k+nz_bhalf))=
	  MEM(data,nx,ny,nz,i,j+ny_half,k).real;
	c_im(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j,k+nz_bhalf))=
	  MEM(data,nx,ny,nz,i,j+ny_half,k).imag;
      }
    }
    for (j=0; j<ny_half; j++) {
      for (k=0; k<nz_bhalf; k++) {
	c_re(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j+ny_bhalf,k))=
	  MEM(data,nx,ny,nz,i,j,k+nz_half).real;
	c_im(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j+ny_bhalf,k))=
	  MEM(data,nx,ny,nz,i,j,k+nz_half).imag;
      }
      for (k=0; k<nz_half; k++) {
	c_re(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j+ny_bhalf,k+nz_bhalf))=
	  MEM(data,nx,ny,nz,i,j,k).real;
	c_im(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j+ny_bhalf,k+nz_bhalf))=
	  MEM(data,nx,ny,nz,i,j,k).imag;
      }
    }
  }
}

static void data_from_temp_xyz( FComplex* data, long nx, long ny, long nz,
				double scale )
{
  long nx_half= nx/2;
  long nx_bhalf= (nx+1)/2;
  long ny_half= ny/2;
  long ny_bhalf= (ny+1)/2;
  long nz_half= nz/2;
  long nz_bhalf= (nz+1)/2;
  long i;
  long j;
  long k;

  for (i=0; i<nx_bhalf; i++) {
    for (j=0; j<ny_bhalf; j++) {
      for (k=0; k<nz_bhalf; k++) {
	MEM(data,nx,ny,nz,i+nx_half,j+ny_half,k+nz_half).real=
	  c_re(MEM(plan_temp,nx,ny,nz,i,j,k))*scale;
	MEM(data,nx,ny,nz,i+nx_half,j+ny_half,k+nz_half).imag=
	  c_im(MEM(plan_temp,nx,ny,nz,i,j,k))*scale;
      }
      for (k=0; k<nz_half; k++) {
	MEM(data,nx,ny,nz,i+nx_half,j+ny_half,k).real=
	  c_re(MEM(plan_temp,nx,ny,nz,i,j,k+nz_bhalf))*scale;
	MEM(data,nx,ny,nz,i+nx_half,j+ny_half,k).imag=
	  c_im(MEM(plan_temp,nx,ny,nz,i,j,k+nz_bhalf))*scale;
      }
    }
    for (j=0; j<ny_half; j++) {
      for (k=0; k<nz_bhalf; k++) {
	MEM(data,nx,ny,nz,i+nx_half,j,k+nz_half).real=
	  c_re(MEM(plan_temp,nx,ny,nz,i,j+ny_bhalf,k))*scale;
	MEM(data,nx,ny,nz,i+nx_half,j,k+nz_half).imag=
	  c_im(MEM(plan_temp,nx,ny,nz,i,j+ny_bhalf,k))*scale;
      }
      for (k=0; k<nz_half; k++) {
	MEM(data,nx,ny,nz,i+nx_half,j,k).real=
	  c_re(MEM(plan_temp,nx,ny,nz,i,j+ny_bhalf,k+nz_bhalf))*scale;
	MEM(data,nx,ny,nz,i+nx_half,j,k).imag=
	  c_im(MEM(plan_temp,nx,ny,nz,i,j+ny_bhalf,k+nz_bhalf))*scale;
      }
    }
  }
  for (i=0; i<nx_half; i++) {
    for (j=0; j<ny_bhalf; j++) {
      for (k=0; k<nz_bhalf; k++) {
	MEM(data,nx,ny,nz,i,j+ny_half,k+nz_half).real=
	  c_re(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j,k))*scale;
	MEM(data,nx,ny,nz,i,j+ny_half,k+nz_half).imag=
	  c_im(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j,k))*scale;
      }
      for (k=0; k<nz_half; k++) {
	MEM(data,nx,ny,nz,i,j+ny_half,k).real=
	  c_re(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j,k+nz_bhalf))*scale;
	MEM(data,nx,ny,nz,i,j+ny_half,k).imag=
	  c_im(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j,k+nz_bhalf))*scale;
      }
    }
    for (j=0; j<ny_half; j++) {
      for (k=0; k<nz_bhalf; k++) {
	MEM(data,nx,ny,nz,i,j,k+nz_half).real=
	  c_re(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j+ny_bhalf,k))*scale;
	MEM(data,nx,ny,nz,i,j,k+nz_half).imag=
	  c_im(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j+ny_bhalf,k))*scale;
      }
      for (k=0; k<nz_half; k++) {
	MEM(data,nx,ny,nz,i,j,k).real=
	  c_re(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j+ny_bhalf,k+nz_bhalf))*scale;
	MEM(data,nx,ny,nz,i,j,k).imag=
	  c_im(MEM(plan_temp,nx,ny,nz,i+nx_bhalf,j+ny_bhalf,k+nz_bhalf))*scale;
      }
    }
  }
}

void fft3d( FComplex* data, long nx, long ny, long nz, long sign, 
	    char* rowcol)
{
  double scale;
  fftw_direction fftw_dir;

  /* Make sure sign flag is valid */
  if ((sign != -1)&&(sign != 1))
    Abort("fft3d: invalid sign: %d\n",sign);

  if (sign==1) fftw_dir= FFTW_BACKWARD;
  else fftw_dir= FFTW_FORWARD;

  if (!strcasecmp(rowcol,"3") || !strcasecmp(rowcol,"xyz")) {
    /* Check for plan and temp space */
    check_plan_xyz((int)nx,(int)ny,(int)nz,fftw_dir);
    
    /* Pack everything into temp space */
    data_to_temp_xyz(data, nx, ny, nz);
    
    /* Do the FFT */
#ifdef FFTW3
    if (fftw_dir==FFTW_FORWARD)
      fftw_execute(plan_xyz_f);
    else
      fftw_execute(plan_xyz_b);
#else
    if (fftw_dir==FFTW_FORWARD)
      fftwnd(plan_xyz_f, 1, plan_temp, 1, 1, NULL, 0, 0);
    else 
      fftwnd(plan_xyz_b, 1, plan_temp, 1, 1, NULL, 0, 0);
#endif
    
    /* Unpack everything from temp space, correcting scale */
    scale = 1.0/sqrt( (float)(nx*ny*nz) );
    data_from_temp_xyz(data, nx, ny, nz, scale);
  }
  else if (!strcasecmp(rowcol,"xy")) {
    /* Check for plan and temp space */
    check_plan_xy((int)nx,(int)ny,(int)nz,fftw_dir);
    
    /* Pack everything into temp space */
    data_to_temp_xy(data, nx, ny, nz);
    
    /* Do the FFT */
#ifdef FFTW3
    if (fftw_dir==FFTW_FORWARD)
      fftw_execute(plan_xy_f);
    else
      fftw_execute(plan_xy_b);
#else
    if (fftw_dir==FFTW_FORWARD) 
      fftwnd(plan_xy_f, nz, plan_temp, nz, 1, NULL, 0, 0);
    else
      fftwnd(plan_xy_b, nz, plan_temp, nz, 1, NULL, 0, 0);
#endif
    
    /* Unpack everything from temp space, correcting scale */
    scale = 1.0/sqrt( (float)(nx*ny) );
    data_from_temp_xy(data, nx, ny, nz, scale);
  }
  else if (!strcasecmp(rowcol,"yz")) {
    /* Check for plan and temp space */
    check_plan_yz((int)nx,(int)ny,(int)nz,fftw_dir);
    
    /* Pack everything into temp space */
    data_to_temp_yz(data, nx, ny, nz);
    
    /* Do the FFT */
#ifdef FFTW3
    if (fftw_dir==FFTW_FORWARD)
      fftw_execute(plan_yz_f);
    else
      fftw_execute(plan_yz_b);
#else
    if (fftw_dir==FFTW_FORWARD) 
      fftwnd(plan_yz_f, nx, plan_temp, 1, ny*nz, NULL, 0, 0);
    else 
      fftwnd(plan_yz_b, nx, plan_temp, 1, ny*nz, NULL, 0, 0);
#endif
    
    /* Unpack everything from temp space, correcting scale */
    scale = 1.0/sqrt( (float)(ny*nz) );
    data_from_temp_yz(data, nx, ny, nz, scale);
  }
  else if (!strcasecmp(rowcol,"xz")) {
    /* Check for plan and temp space */
#ifdef FFTW3
    check_plan_xz((int)nx,(int)ny,(int)nz,fftw_dir);
#else
    check_plan_x((int)nx,(int)ny,(int)nz,fftw_dir);
    check_plan_z((int)nx,(int)ny,(int)nz,fftw_dir);
#endif
    
    /* Pack everything into temp space */
    data_to_temp_xz(data, nx, ny, nz);
    
    /* Do the FFT */
#ifdef FFTW3
    if (fftw_dir==FFTW_FORWARD)
      fftw_execute(plan_xz_f);
    else
      fftw_execute(plan_xz_b);
#else
    if (fftw_dir==FFTW_FORWARD) {
      fftw(plan_x_f, ny*nz, plan_temp, ny*nz, 1, NULL, 0, 0);
      fftw(plan_z_f, nx*ny, plan_temp, 1, nz, NULL, 0, 0);
    }
    else {
      fftw(plan_x_b, ny*nz, plan_temp, ny*nz, 1, NULL, 0, 0);
      fftw(plan_z_b, nx*ny, plan_temp, 1, nz, NULL, 0, 0);
    }
#endif
    
    /* Unpack everything from temp space, correcting scale */
    scale = 1.0/sqrt( (float)(nx*nz) );
    data_from_temp_xz(data, nx, ny, nz, scale);
  }
  else if (!strcasecmp(rowcol,"x")) {
    /* Check for plan and temp space */
    check_plan_x((int)nx,(int)ny,(int)nz,fftw_dir);
    
    /* Pack everything into temp space */
    data_to_temp_x(data, nx, ny, nz);
    
    /* Do the FFT */
#ifdef FFTW3
    if (fftw_dir==FFTW_FORWARD)
      fftw_execute(plan_x_f);
    else
      fftw_execute(plan_x_b);
#else
    if (fftw_dir==FFTW_FORWARD)
      fftw(plan_x_f, ny*nz, plan_temp, ny*nz, 1, NULL, 0, 0);
    else
      fftw(plan_x_b, ny*nz, plan_temp, ny*nz, 1, NULL, 0, 0);
#endif
    
    /* Unpack everything from temp space, correcting scale */
    scale = 1.0/sqrt( (float)nx );
    data_from_temp_x(data, nx, ny, nz, scale);
  }
  else if (!strcasecmp(rowcol,"y")) {
    int i;

    /* Check for plan and temp space */
    check_plan_y((int)nx,(int)ny,(int)nz,fftw_dir);
    
    /* Pack everything into temp space */
    data_to_temp_y(data, nx, ny, nz);
    
    /* Do the FFT */
#ifdef FFTW3
    if (fftw_dir==FFTW_FORWARD)
      fftw_execute(plan_y_f);
    else
      fftw_execute(plan_y_b);
#else
    if (fftw_dir==FFTW_FORWARD)
      for (i=0; i<nx; i++) 
	fftw(plan_y_f, nz, plan_temp+(i*ny*nz), nz, 1, NULL, 0, 0);
    else
      for (i=0; i<nx; i++) 
	fftw(plan_y_b, nz, plan_temp+(i*ny*nz), nz, 1, NULL, 0, 0);
#endif
    
    /* Unpack everything from temp space, correcting scale */
    scale = 1.0/sqrt( (float)ny );
    data_from_temp_y(data, nx, ny, nz, scale);
  }
  else if (!strcasecmp(rowcol,"z")) {
    /* Check for plan and temp space */
    check_plan_z((int)nx,(int)ny,(int)nz,fftw_dir);
    
    /* Pack everything into temp space */
    data_to_temp_z(data, nx, ny, nz);
    
    /* Do the FFT */
#ifdef FFTW3
    if (fftw_dir==FFTW_FORWARD)
      fftw_execute(plan_z_f);
    else
      fftw_execute(plan_z_b);
#else
    if (fftw_dir==FFTW_FORWARD)
      fftw(plan_z_f, nx*ny, plan_temp, 1, nz, NULL, 0, 0);
    else
      fftw(plan_z_b, nx*ny, plan_temp, 1, nz, NULL, 0, 0);
#endif
    
    /* Unpack everything from temp space, correcting scale */
    scale = 1.0/sqrt( (float)nz );
    data_from_temp_z(data, nx, ny, nz, scale);
  }
  else {
    Abort( "FFT3D: Unrecognized row-column indicator (%s).", rowcol );
  }
  
  return;
}

