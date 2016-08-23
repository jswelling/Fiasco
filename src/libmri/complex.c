/*
 *	complex.c - Example MRI library program
 *
 *	This program creates a dataset consisting of complex-valued
 *	image data.   It generates bogus images one by one and places
 *	them into the "images" chunk of the new dataset.
 *
 *	Copyright (c) 1996 Pittsburgh Supercomputing Center
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
 *
 *	HISTORY:
 *		4/96 - Written by Greg Hood
 */
#include <stdio.h>
#include <math.h>
#include "mri.h"

#define PI	3.14159265358979323

typedef struct Complex {
  short real;
  short imag;
} Complex;

static char rcsid[] = "$Id: complex.c,v 1.3 1999/07/07 20:11:27 welling Exp $";

main ()
{
  MRI_Dataset *ds;
  int t, s;
  int x, y;
  double r, scale;
  Complex image[16][16];

  printf("Creating mri file complex.mri\n");

  /* create a new dataset */
  ds = mri_open_dataset("complex.mri", MRI_WRITE);

  /* create images chunk with 3 times, 4 slices and 16x16 voxels
     in each image; each voxel will be complex-valued */
  mri_create_chunk(ds, "images");
  mri_set_string(ds, "images.datatype", "int16");
  mri_set_string(ds, "images.dimensions", "vxyzt");
  mri_set_int(ds, "images.extent.v", 2);
  mri_set_int(ds, "images.extent.x", 16);
  mri_set_int(ds, "images.extent.y", 16);
  mri_set_int(ds, "images.extent.z", 4);
  mri_set_int(ds, "images.extent.t", 3);

  /* inform the mri library that we want the
     images chunk to go into its own file
     with extension .dat, which in this case
     would be test.dat */
  mri_set_string(ds, "images.file", ".dat");

  /* create bogus images and save them */
  for (t = 0; t < 3; ++t)
    for (s = 0; s < 4; ++s)
      {
	scale = cos (PI * s / 10.0) * PI / 40.0;
	for (x = 0; x < 16; ++x)
	  for (y = 0; y < 16; ++y)
	    {
	      image[x][y].real = x + y;
	      image[x][y].imag = x - y;
	    }
	mri_set_image(ds, t, s, MRI_COMPLEX_SHORT, &image[0][0]);
      }

  /* make sure everything is written out */
  mri_close_dataset(ds);
}
