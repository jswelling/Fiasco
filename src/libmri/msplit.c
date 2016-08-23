/*
 *	msplit.c - Example MRI library program
 *
 *	This program takes an existing dataset and makes sure that
 *	the "images" chunk is split off by itself into a separate
 *	data file "test.dat" containing pure binary data.  If the
 *	"images" chunk was already there, no actual reading or
 *	writing takes place.
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

static char rcsid[] = "$Id: msplit.c,v 1.3 1999/07/07 20:11:54 welling Exp $";

main ()
{
  MRI_Dataset *ds;

  printf("Splitting test.mri into test.mri and test.dat\n");
  ds = mri_open_dataset("test.mri", MRI_MODIFY);
  mri_set_string(ds, "images.file", ".dat");
  mri_close_dataset(ds);
}
