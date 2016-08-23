/*
 *	imean.c - Example MRI library program
 *
 *	This program computes the mean voxel value for each image
 *	in a dataset, then saves these values in a new chunk which is
 *	added to the dataset.  This illustrates how new chunks may
 *	created and added to existing datasets.
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
#include "mri.h"

static char rcsid[] = "$Id: imean.c,v 1.3 1999/07/07 20:11:34 welling Exp $";

int main ()
{
  MRI_Dataset *ds;
  int n_times, n_slices, n_pixels;
  int t, s, i;
  short *image;
  double total;
  float *means;

  ds = mri_open_dataset("test.mri", MRI_MODIFY);  /* load in the test data set;
						     allow modifications to it */

  /* determine how many elements are present */
  n_times = mri_get_int(ds, "images.extent.t");
  n_slices = mri_get_int(ds, "images.extent.z");
  n_pixels = mri_get_int(ds, "images.extent.x") *
             mri_get_int(ds, "images.extent.y");

  /* allocate an array to temporarily hold the means */
  means = (float *) malloc(n_times * n_slices * sizeof(float));

  /* compute the mean pixel value */
  for (t = 0; t < n_times; ++t)
    for (s = 0; s < n_slices; ++s)
    {
      total = 0.0;
      image = mri_get_image(ds, t, s, MRI_SHORT);
      for (i = 0; i < n_pixels; ++i)
	total += image[i];
      means[t*n_slices + s] = total / n_pixels;
    }

  /* store the mean values in the dataset */
  mri_create_chunk(ds, "means");
  mri_set_string(ds, "means.datatype", "float32");
  mri_set_string(ds, "means.dimensions", "zt");
  mri_set_int(ds, "means.extent.z", n_slices);
  mri_set_int(ds, "means.extent.t", n_times);
  mri_set_chunk(ds, "means", n_times*n_slices, 0, MRI_FLOAT, means);

  /* write the dataset back out */
  mri_close_dataset(ds);
}

