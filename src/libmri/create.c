/*
 *	create.c - Example MRI library program
 *
 *	This program creates a new dataset from scratch.  It generates
 *	bogus images one by one and places them into the "images" chunk
 *	of the new dataset.
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
 *		3/96 - Written by Greg Hood
 */
#include <stdio.h>
#include <math.h>
#include "mri.h"

#define PI	3.14159265358979323

static char rcsid[] = "$Id: create.c,v 1.3 1999/07/07 20:11:29 welling Exp $";

main ()
{
  MRI_Dataset *ds;
  int t, s;
  int x, y;
  double r, scale;
  short image[64][64];

  printf("Creating mri file test.mri\n");

  /* create a new dataset */
  ds = mri_open_dataset("test.mri", MRI_WRITE);

  /* create images chunk with 3 times, 4 slices and 64x64 voxels
     in each image */
  mri_create_chunk(ds, "images");
  mri_set_string(ds, "images.datatype", "int16");
  mri_set_string(ds, "images.dimensions", "xyzt");
  mri_set_int(ds, "images.extent.x", 64);
  mri_set_int(ds, "images.extent.y", 64);
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
	for (x = 0; x < 64; ++x)
	  for (y = 0; y < 64; ++y)
	    {
	      r = sqrt(2.0 * (x - 32) * (x - 32) + (y - 32) * (y - 32));
	      image[y][x] = 4000.0 * cos(scale * r) * (t + 5) / 6.0;
	      if (image[y][x] < 0.0)
		image[y][x] = 0.0;
	      image[y][x] = x+y;
	    }
	mri_set_image(ds, t, s, MRI_SHORT, &image[0][0]);
      }

  /* make sure everything is written out */
  mri_close_dataset(ds);
}
