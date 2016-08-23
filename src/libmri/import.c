/*
 *	import.c - Example MRI library program
 *
 *	This program imports a chunk of data from an external
 *	data file "complex.dat".  It then copies this data into the
 *	test dataset.
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

static char rcsid[] = "$Id: import.c,v 1.3 1999/07/07 20:11:36 welling Exp $";

main ()
{
  MRI_Dataset *ds;
  int t, s;
  int x, y;
  double r, scale;
  short image[64][64];

  /* open the dataset */
  ds = mri_open_dataset("test.mri", MRI_MODIFY);

  /* import an images chunk with 3 times, 4 slices and 16x16
     complex voxels from the file complex.dat */
  mri_create_chunk(ds, "images");
  mri_set_string(ds, "images.datatype", "int16");
  mri_set_string(ds, "images.dimensions", "vxyzt");
  mri_set_int(ds, "images.extent.v", 2);
  mri_set_int(ds, "images.extent.x", 16);
  mri_set_int(ds, "images.extent.y", 16);
  mri_set_int(ds, "images.extent.z", 4);
  mri_set_int(ds, "images.extent.t", 3);
  mri_set_string(ds, "images.file", "complex.dat");
  mri_set_string(ds, "images.order", "external");
  mri_set_int(ds, "images.offset", 0);

  /* this call causes all of the above changes to become
     "official"; a read or write to the chunk would also
     do the same thing */
  mri_update_chunk(ds, "images");

  /* now change the location of the chunk; this effectively
     causes a copy to be made into the test dataset */
  mri_set_string(ds, "images.file", ".dat");
  mri_remove(ds, "images.order");

  /* make sure everything is written out */
  mri_close_dataset(ds);
}
