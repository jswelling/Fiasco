/*
 *	endian.c - Example MRI library program
 *
 *	This program takes an existing dataset, asks the user what
 *	endianness the binary "images" data should be in,
 *	and makes sure that it is written out with the specified
 *	endianness.
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

static char rcsid[] = "$Id: endian.c,v 1.3 1999/07/07 20:11:30 welling Exp $";

main ()
{
  MRI_Dataset *ds;
  int answer;

  printf("Enter 1 to save images in little-endian format.\n");
  printf("Enter 0 to save images in big-endian format.\n");
  printf("Choice: ");
  scanf("%d", &answer);

  ds = mri_open_dataset("test.mri", MRI_MODIFY);
  mri_set_int(ds, "images.little_endian", answer);
  mri_close_dataset(ds);
}
