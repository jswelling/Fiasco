/*
 *	single.c - Example MRI library program
 *
 *	This program takes an existing dataset and makes sure that
 *	the "images" chunk is placed back into the "test.mri" header
 *	file.  If this was the only chunk in the dataset, the
 *	resulting dataset will consist of just one file.
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

static char rcsid[] = "$Id: single.c,v 1.3 1999/07/07 20:11:59 welling Exp $";

main ()
{
  MRI_Dataset *ds;

  printf("Combining test.mri and test.dat into a single file: test.mri\n");
  ds = mri_open_dataset("test.mri", MRI_MODIFY);
  mri_set_string(ds, "images.file", ".mri");
  mri_close_dataset(ds);
}
