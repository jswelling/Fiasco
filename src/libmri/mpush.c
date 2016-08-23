/*
 *	mpush.c - Example MRI library program
 *
 *	This program creates a single image composed of
 *	32 bit integers.  It creates a short int array and
 *	then uses the automatic data conversion to put
 *	this into an int32 chunk.
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

static char rcsid[] = "$Id: mpush.c,v 1.3 1999/07/07 20:11:48 welling Exp $";

main ()
{
  MRI_Dataset *ds;
  int x, y;
  short image[16][16];

  printf("Creating dataset test.mri\n");

  /* open the dataset for writing */
  ds = mri_open_dataset("test.mri", MRI_WRITE);

  mri_create_chunk(ds, "images");
  mri_set_string(ds, "images.datatype", "int32");
  mri_set_string(ds, "images.dimensions", "xyz");
  mri_set_int(ds, "images.extent.x", 16);
  mri_set_int(ds, "images.extent.y", 16);
  mri_set_int(ds, "images.extent.z", 1);
  for (y = 0; y < 16; ++y)
    for (x = 0; x < 16; ++x)
      image[y][x] = x * (y - 1);
  mri_set_image(ds, 0, 0, MRI_SHORT, &image[0][0]);
  mri_close_dataset(ds);
}
