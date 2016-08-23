/*
 *	mpull.c - Example MRI library program
 *
 *	This program pulls the first image out of the test dataset
 *	as an array of floats, as opposed to the chunk's actual
 *	datatype.  It prints the upper left 5 x 5 corner's
 *	voxel values.
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

static char rcsid[] = "$Id: mpull.c,v 1.3 1999/07/07 20:11:46 welling Exp $";

main ()
{
  MRI_Dataset *ds;
  int x, y;
  int row_length, col_length;
  float *image;

  /* open the dataset for reading only */
  ds = mri_open_dataset("test.mri", MRI_READ);

  image = mri_get_image(ds, 0, 0, MRI_FLOAT);
  row_length = mri_get_int(ds, "images.extent.x");
  col_length = mri_get_int(ds, "images.extent.y");
  for (y = 0; y < 5 && y < col_length; ++y)
    {
      for (x = 0; x < 5 && x < row_length; ++x)
	printf("%f ", image[y*row_length + x]);
      printf("\n");
    }

  mri_close_dataset(ds);
}
