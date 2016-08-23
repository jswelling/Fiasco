/*
 *	halve.c - Example MRI library program
 *
 *	This program halves the x dimension of a dataset while
 *	copying it to another
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
 *		5/96 - Written by Greg Hood
 */
#include <stdio.h>
#include "mri.h"

static char rcsid[] = "$Id: halve.c,v 1.3 1999/07/07 20:11:33 welling Exp $";

int main ()
{
  char name[256];
  MRI_Dataset *ds;
  char new_name[256];
  MRI_Dataset *nds;
  int xwidth;

  printf("Dataset to be copied: ");
  scanf("%s", name);
  ds = mri_open_dataset(name, MRI_READ);	/* load in the original */
  printf("Name of new dataset: ");
  scanf("%s", new_name);

  /* copy */
  nds = mri_copy_dataset(new_name, ds);
  mri_set_string(nds, "images.file", ".hdat");
  xwidth = mri_get_int(ds, "images.extent.x")/2;
  mri_set_int(nds, "images.extent.x", xwidth);

  /* close both the datasets */
  mri_close_dataset(ds);
  mri_close_dataset(nds);
}

