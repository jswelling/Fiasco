/*
 *	mean.c - Example MRI library program
 *
 *	This program computes the mean voxel value over all images
 *	in a dataset, then saves this value in the header.  This illustrates
 *	how additional fields may be added to existing datasets.
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

static char rcsid[] = "$Id: mean.c,v 1.3 1999/07/07 20:11:44 welling Exp $";

int main ()
{
  MRI_Dataset *ds;
  int n_times, n_slices, n_pixels;
  int t, s, i;
  short *image;
  double total, mean;

  ds = mri_open_dataset("test", MRI_MODIFY);	/* load in the test data set;
						   allow modifications to it */

  /* determine how many elements are present */
  n_times = mri_get_int(ds, "images.extent.t");
  n_slices = mri_get_int(ds, "images.extent.z");
  n_pixels = mri_get_int(ds, "images.extent.x") *
             mri_get_int(ds, "images.extent.y");

  /* compute the mean pixel value */
  total = 0.0;
  for (t = 0; t < n_times; ++t)
    for (s = 0; s < n_slices; ++s)
    {
      image = mri_get_image(ds, t, s, MRI_SHORT);
      for (i = 0; i < n_pixels; ++i)
	total += image[i];
    }
  mean = total / (n_times * n_slices * n_pixels);

  /* print out the result */
  printf("mean is %f\n",  mean);
  
  /* save the mean back in the dataset */
  mri_set_float(ds, "mean", mean);

  /* write the dataset back out */
  mri_close_dataset(ds);
}

