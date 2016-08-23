/*
 *	pgh_write.c - routines for writing out PGH files
 *
 *	Copyright (c) 1996  Pittsburgh Supercomputing Center
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
 *	HISTORY
 *		4/96	Written by Greg Hood (PSC)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bio.h"
#include "mri.h"
#include "fmri.h"
#include "gift.h"

static char rcsid[] = "$Id: pgh_write.c,v 1.4 2000/10/06 04:55:17 welling Exp $";

Filename pgh_write_file;
MRI_Dataset *pgh_write_ds;
Boolean pgh_write_split_data = TRUE;

void
PghStartWriting (Filename basename)
{
  Filename header;
  unsigned char buf[512];
  unsigned char *d;
  FILE *f;
  int i;
  int len;
  char *type;
  char dims[16];
  char key[64];
  char *dim_chars = "vxyzt";
  
  len = strlen(basename);
  if (len >= 4 &&
      strcmp(&basename[len-4], ".mri") == 0)
    strcpy(pgh_write_file, basename);
  else
    sprintf(pgh_write_file, "%s.mri", basename);

  if ((pgh_write_ds = mri_open_dataset(pgh_write_file, MRI_WRITE)) == NULL)
    Abort("Can't open dataset %s for writing\n", pgh_write_file);
  hist_add(pgh_write_ds,"*gift*");

  mri_create_chunk(pgh_write_ds, "images");
  switch (fh.data_type)
    {
    case GIFT_UINT8:
      type = "uint8";
      break;
    case GIFT_INT16:
      type = "int16";
      break;
    case GIFT_INT32:
      type = "int32";
      break;
    case GIFT_FLOAT32:
      type = "float32";
      break;
    case GIFT_FLOAT64:
      type = "float64";
      break;
    default:
      Abort("PghStartWriting: invalid datatype\n");
      break;
    }
  mri_set_string(pgh_write_ds, "images.datatype", type);
  dims[0] = '\0';
  if (fh.dim[0].n != 1)
    strcat(dims, "v");
  strcat(dims, "xyz");
  if (fh.dim[4].n != 1)
    strcat(dims, "t");
  mri_set_string(pgh_write_ds, "images.dimensions", dims);
  for (i = 0; i < fh.n_dims; ++i)
    {
      if ((i == 0 || i > 3) && fh.dim[i].n == 1)
	continue;
      sprintf(key, "images.extent.%c", dim_chars[i]);
      mri_set_int(pgh_write_ds, key, fh.dim[i].n);
      if (fh.dim[i].min != 0)
	{
	  sprintf(key, "images.min.%c", dim_chars[i]);
	  mri_set_int(pgh_write_ds, key, fh.dim[i].min);
	}
      if (fh.dim[i].stride != 1)
	{
	  sprintf(key, "images.stride.%c", dim_chars[i]);
	  mri_set_int(pgh_write_ds, key, fh.dim[i].stride);
	}
      if (fh.dim[i].size != 0.0)
	{
	  sprintf(key, "images.size.%c", dim_chars[i]);
	  mri_set_float(pgh_write_ds, key, fh.dim[i].size);
	}
    }
  if (pgh_write_split_data)
    mri_set_string(pgh_write_ds, "images.file", ".dat");

  for (i = 0; i < fh.n_images; ++i)
    if (fh.corrupt[i])
      {
	mri_create_chunk(pgh_write_ds, "images.corrupt");
	mri_set_string(pgh_write_ds, "images.corrupt.dimensions", "zt");
	mri_set_int(pgh_write_ds, "images.corrupt.extent.z", fh.dim[3].n);
	mri_set_int(pgh_write_ds, "images.corrupt.extent.t", fh.dim[3].n);
	mri_set_chunk(pgh_write_ds, "images.corrupt", fh.n_images, 0, MRI_UNSIGNED_CHAR, fh.corrupt);
	break;
      }
}

void
PghWriteImage (int time,
	       int slice)
{
  int t, s;
  MRI_ArrayType type;

  t = (time - fh.dim[4].min) / fh.dim[4].stride;
  s = (slice - fh.dim[3].min) / fh.dim[3].stride;
  switch (fh.data_type)
    {
    case GIFT_UINT8:
      type = MRI_UNSIGNED_CHAR;
      break;
    case GIFT_INT16:
      type = MRI_SHORT;
      break;
    case GIFT_INT32:
      type = MRI_INT;
      break;
    case GIFT_FLOAT32:
      type = MRI_FLOAT;
      break;
    case GIFT_FLOAT64:
      type = MRI_DOUBLE;
      break;
    default:
      Abort("PghWriteImage: invalid datatype\n");
      break;
    }
  mri_set_image(pgh_write_ds, t, s, type, image);
}

void
PghEndWriting ()
{
  mri_close_dataset(pgh_write_ds);
}
