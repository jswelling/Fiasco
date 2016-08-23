/************************************************************
 *                                                          *
 *  fshrot.c                                                *
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
 * fourier_shift_rot performs a translation and rotation     *
 *   of a complex image using the Fourier shift theorem      *
 *   and three-pass shearing for rotations                   *
 * In the image domain, the translation is applied           *
 *   prior to the rotation                                   *
 * For details, see "Improved Registration using             *
 *   Fourier Interpolation" by Eddy, Fitzgerald, Noll,       *
 *   Journal of Magnetic Resonance in Medicine (1996)        *
 *                                                           *
 * par is the RegPars structure containing the amount        *
 *   of translation and rotation to apply                    *
 * orig_image is the 2D, complex array containing the        *
 *   image to be translated and rotated                      *
 * moved_image is the returned 2D, complex array containing  *
 *   translated, rotated version of orig_image               *
 * nx and ny and the slowest-varying and fastest-varying     *
 *   dimensions of orig_image (and moved_image) respectively *
 * domain indicates the domain of the image to which the     *
 *   movement is to be applied:                              *
 *   'f' or 'k' if orig_image is a Fourier (k-space) image   *
 *   't' or 'i' if orig_image is a regular (i-space) image   * 
 *************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"

static char rcsid[] = "$Id: fshrot.c,v 1.10 2007/03/21 23:50:20 welling Exp $";

void fourier_shift_rot( RegPars par, FComplex** orig_image, 
			FComplex** moved_image, long nx, long ny, 
			char domain )
{
  short fourier_flag;
  FComplex **temp_image1 = NULL, **temp_image2 = NULL;
  double rotate, shift_mult, temp_shift_mult;
  long halfx, halfy;
  long nx_mod, ny_mod;
  double angle, temp_angle, cos_ang, sin_ang;
  long x, y;

  /* Determine domain */
  switch( domain )
    {
    case 'i': case 'I': case 't': case 'T':
      fourier_flag = 0;
      break;
    case 'f': case 'F': case 'k': case 'K':
      fourier_flag = 1;
      break;
    default:
      Abort( "Domain argument to fourier_shift_rot unrecognized: (%c).",
	     domain );
    }
      

  /* Calculate midpoints */
  halfx = nx / 2;
  halfy = ny / 2;
  if (nx % 2) nx_mod= nx-1;
  else nx_mod= nx;
  if (ny % 2) ny_mod= ny-1;
  else ny_mod= ny;

  /* Allocate space for intermediate images */
  temp_image1 = Matrix( ny, nx, FComplex );
  temp_image2 = Matrix( ny, nx, FComplex );

  if( !fourier_flag && ( par.x_shift || par.y_shift ) )
    {
      /* Translation---before rotation only if image domain data */
      /*   and if translation is non-zero                        */

      /* Copy original image into local storage */
      memcpy( *temp_image1, *orig_image, 
	      (size_t) ( sizeof(FComplex) * nx * ny ) );

      /* Perform full 2D FFT */
      fft2d( temp_image1, nx, ny, 1, '2', 0, 0 );

      /* Apply linear phase change to transformed data */
      for( y = 0; y < ny; y++ )
	{
	  /* Calculate y-translation contribution to phase change */
	  temp_angle = -2.0 * M_PI * par.y_shift *
	    (float) ( y - halfy ) / (float) (ny_mod);

	  for( x = 0; x < nx; x++ )
	    {
	      /* Calculate phase shift to apply for particular location */
	      angle = temp_angle - 2.0 * M_PI * par.x_shift *
		(float) ( x - halfx ) / (float) (nx_mod);
	      cos_ang = cos( angle );
	      sin_ang = sin( angle );

	      /* Apply phase shift */
	      temp_image2[y][x].real = cos_ang * temp_image1[y][x].real -
		sin_ang * temp_image1[y][x].imag;
	      temp_image2[y][x].imag = cos_ang * temp_image1[y][x].imag +
		sin_ang * temp_image1[y][x].real;
	    }
	}

      /* Perform full 2D inverse FFT */
      fft2d( temp_image2, nx, ny, -1, '2', 0, 0 );

    }
  else
    {
      /* Don't translate, so just copy original image */
      /*   into appropriate local storage             */
      memcpy( *temp_image2, *orig_image, 
	      (size_t) ( sizeof(FComplex) * nx * ny ) );
    }
      

  if( par.rotation )
    {
      /* Apply rotation if non-zero */

      /* SHEARING ROTATION */

      /* Convert rotation from degrees to radians */
      rotate = par.rotation * M_PI / 180.0;

      /*** First Shear ***/

      /* Calculate constant (with respect to Fourier location) */
      /*   of exponent for linear shift term                   */
      shift_mult = -2.0 * M_PI * tan( ( rotate / 2.0 ) );

      /* Inverse Fourier transform across x-dimension */
      fft2d( temp_image2, nx, ny, -1, 'x', 0, ny );

      /* Apply linear phase change */
      for( y = 0; y < ny; y++ )
	{
	  /* Update shift multiplier for this column */
	  temp_shift_mult = shift_mult * (float) ( y - halfy );

	  for( x = 0; x < nx; x++ )
	    {
	      /* Calculate phase shift for specific Fourier location */
	      angle = temp_shift_mult * (float) ( x - halfx );

	      /* Scale angle appropriately, depending on which    */
	      /*   dimension is in Fourier space and which is not */
	      if( fourier_flag )
		angle /= (float) (ny_mod);
	      else
		angle /= (float) (nx_mod);

	      /* Apply phase shift */
	      cos_ang = cos( angle );
	      sin_ang = sin( angle );
	      temp_image1[y][x].real = cos_ang * temp_image2[y][x].real -
		sin_ang * temp_image2[y][x].imag;
	      temp_image1[y][x].imag = cos_ang * temp_image2[y][x].imag +
		sin_ang * temp_image2[y][x].real;
		
	    }
	}
      
      /* Fourier transform across x-dimension */
      fft2d( temp_image1, nx, ny, 1, 'x', 0, ny );

      /*** Second Shear ***/
      
      /* Calculate constant (with respect to Fourier location) */
      /*   of exponent for linear shift term                   */
      shift_mult = -2.0 * M_PI * sin( (float) (-1.0 * rotate) );

      /* Inverse Fourier transform across y-dimension */
      fft2d( temp_image1, nx, ny, -1, 'y', 0, nx );

      /* Apply linear phase change */
      for( x = 0; x < nx; x++ )
	{
	  /* Update shift multiplier for this column */
	  temp_shift_mult = shift_mult * (float) ( x - halfx );

	  for( y = 0; y < ny; y++ )
	    {
	      /* Calculate phase shift for specific Fourier location */
	      angle = temp_shift_mult * (float) ( y - halfy );

	      /* Scale angle appropriately, depending on which    */
	      /*   dimension is in Fourier space and which is not */
	      if( fourier_flag )
		angle /= (float) (nx_mod);
	      else
		angle /= (float) (ny_mod);

	      /* Apply phase shift */
	      cos_ang = cos( angle );
	      sin_ang = sin( angle );
	      temp_image2[y][x].real = cos_ang * temp_image1[y][x].real -
		sin_ang * temp_image1[y][x].imag;
	      temp_image2[y][x].imag = cos_ang * temp_image1[y][x].imag +
		sin_ang * temp_image1[y][x].real;
		
	    }
	}
      
      /* Fourier transform across y-dimension */
      fft2d( temp_image2, nx, ny, 1, 'y', 0, nx );

      /*** Third Shear ***/

      /* Calculate constant (with respect to Fourier location) */
      /*   of exponent for linear shift term                   */
      shift_mult = -2.0 * M_PI * tan( ( rotate / 2.0 ) );

      /* Inverse Fourier transform across x-dimension */
      fft2d( temp_image2, nx, ny, -1, 'x', 0, ny );

      /* Apply linear phase change */
      for( y = 0; y < ny; y++ )
	{
	  /* Update shift multiplier for this column */
	  temp_shift_mult = shift_mult * (float) ( y - halfy );

	  for( x = 0; x < nx; x++ )
	    {
	      /* Calculate phase shift for specific Fourier location */
	      angle = temp_shift_mult * (float) ( x - halfx );

	      /* Scale angle appropriately, depending on which    */
	      /*   dimension is in Fourier space and which is not */
	      if( fourier_flag )
		angle /= (float) (ny_mod);
	      else
		angle /= (float) (nx_mod);

	      /* Apply phase shift */
	      cos_ang = cos( angle );
	      sin_ang = sin( angle );
	      temp_image1[y][x].real = cos_ang * temp_image2[y][x].real -
		sin_ang * temp_image2[y][x].imag;
	      temp_image1[y][x].imag = cos_ang * temp_image2[y][x].imag +
		sin_ang * temp_image2[y][x].real;
		
	    }
	}
      
      /* Fourier transform across x-dimension */
      fft2d( temp_image1, nx, ny, 1, 'x', 0, ny );

      /* END SHEARING */

    }

  else

    {
      /* Since there is zero rotation, simply copy image */
      /*   into appropriate local storage                */
      memcpy( *temp_image1, *temp_image2,
	      (size_t) ( sizeof(FComplex) * nx * ny ) );
    }


  if( fourier_flag && ( par.x_shift || par.y_shift ) )
    {
      /* Translation---phase change applied after rotation for */
      /*   Fourier domain data (if translation is non-zero)    */
      /*   [Translation is still intrinsically prior to        */
      /*    the rotation in image space]                       */

      /* We are in k space already, so no need to do full 2D fft */

      /* Apply linear phase change to transformed data */
      for( y = 0; y < ny; y++ )
	{
	  /* Calculate y-translation contribution to phase change */
	  temp_angle = -2.0 * M_PI * par.y_shift *
	    (float) ( y - halfy ) / (float) (ny_mod);

	  for( x = 0; x < nx; x++ )
	    {
	      /* Calculate phase shift to apply for particular location */
	      angle = temp_angle + -2.0 * M_PI * par.x_shift *
		(float) ( x - halfx ) / (float) (nx_mod);
	      cos_ang = cos( angle );
	      sin_ang = sin( angle );

	      /* Apply phase shift */
	      moved_image[y][x].real = cos_ang * temp_image1[y][x].real -
		sin_ang * temp_image1[y][x].imag;
	      moved_image[y][x].imag = cos_ang * temp_image1[y][x].imag +
		sin_ang * temp_image1[y][x].real;
	    }
	}

      /* Still in k space, so no need to perform full 2D fft */

    }
  else
    {
      /* Don't translate, so just copy original image */
      /*   into appropriate local storage             */
      memcpy( *moved_image, *temp_image1,
	      (size_t) ( sizeof(FComplex) * nx * ny ) );
    }
      
  /* Free local storage */
  FreeMatrix( temp_image1 );
  FreeMatrix( temp_image2 );

  return;

}

