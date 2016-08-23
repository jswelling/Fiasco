/*
 *	pgh_read.c - routines for reading in PGH files
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

static char rcsid[] = "$Id: pgh_read.c,v 1.5 2000/11/02 20:15:02 welling Exp $";

Filename pgh_read_file;
MRI_Dataset *pgh_read_ds;
int pgh_read_orig_time_min;
int pgh_read_orig_time_stride;
int pgh_read_orig_slice_min;
int pgh_read_orig_slice_stride;

Boolean
PghCheckFormat (Filename basename,
		FileList files)
{
  FileListElement *fle;
  Filename header;
  FILE *f;
  int i;
  Boolean pgh_format;
  char buf[16];
  
  /* for PGH, check for .mri extension and !format = pgh */
  header[0] = '\0';
  for (fle = files; fle != NULL; fle = fle->next)
    if (HasExtension(fle->name, ".mri"))
      {
	strcpy(header, fle->name);
	break;
      }	
  if (fle == NULL)
    /* no .mri file was found */
    return(FALSE);

  if ((f = fopen(header, "r")) == NULL)
    Abort("Can't open header file %s\n", header);

  pgh_format = TRUE;
  if (fread(buf, 1, 16, f) != 16 ||
      strncmp(buf, "!format = pgh", 13) != 0)
    pgh_format = FALSE;

  fclose(f);
  return(pgh_format);
}

void
PghStartReading (Filename basename,
		 FileList files)
{
  FileList fle;
  FILE *f;
  int i;
  int len;
  char *type;
  char *dims;
  char dim[2];
  unsigned char *ptr;

  /* check that only one header file exists */
  pgh_read_file[0] = '\0';
  for (fle = files; fle != NULL; fle = fle->next)
    if (HasExtension(fle->name, ".mri"))
      {
	if (pgh_read_file[0] != '\0')
	  Abort("More than one .mri file was found beginning with %s\n",
		basename);
	strcpy(pgh_read_file, fle->name);
      }

  if (pgh_read_file[0] == '\0')
    /* no .mri file was found */
    Abort("Cannot find .mri file!\n");
  
  if ((pgh_read_ds = mri_open_dataset(pgh_read_file, MRI_READ)) == NULL)
    Abort("Can't open header file %s\n", pgh_read_file);

  fh.class = GIFT_I_SPACE;
  if (!mri_has(pgh_read_ds, "images.datatype"))
    Abort("Dataset does not contain images.datatype key\n");
  type = mri_get_string(pgh_read_ds, "images.datatype");
  if (strcmp(type, "uint8") == 0)
    fh.data_type = GIFT_UINT8;
  else if (strcmp(type, "int16") == 0)
    fh.data_type = GIFT_INT16;
  else if (strcmp(type, "int32") == 0)
    fh.data_type = GIFT_INT32;
  else if (strcmp(type, "float32") == 0)
    fh.data_type = GIFT_FLOAT32;
  else if (strcmp(type, "float64") == 0)
    fh.data_type = GIFT_FLOAT64;
  else
    Abort("Unrecognized images.datatype\n");

  if (mri_has(pgh_read_ds, "images.dimensions"))
    dims = mri_get_string(pgh_read_ds, "images.dimensions");
  else
    dims = "xyzt";

  fh.n_dims = strlen(dims);
  for (i = 0; i < fh.n_dims; ++i)
    {
      dim[0] = dims[i];
      dim[1] = '\0';
      switch (dims[i])
	{
	case 'v':
	  fh.dim[i].type = GIFT_V_DIMENSION;
	  break;
	case 'x':
	  fh.dim[i].type = GIFT_X_DIMENSION;
	  break;
	case 'y':
	  fh.dim[i].type = GIFT_Y_DIMENSION;
	  break;
	case 'z':
	  fh.dim[i].type = GIFT_Z_DIMENSION;
	  break;
	case 't':
	  fh.dim[i].type = GIFT_T_DIMENSION;
	  break;
	default:
	  Abort("Unknown dimension in dataset: %c\n", dims[i]);
	  break;
	}
      if (mri_has(pgh_read_ds, mri_cat("images.stride.", dim)))
	fh.dim[i].stride = mri_get_int(pgh_read_ds, mri_cat("images.stride", dim));
      else
	fh.dim[i].stride = 1;
      if (mri_has(pgh_read_ds, mri_cat("images.min.", dim)))
	fh.dim[i].min = mri_get_int(pgh_read_ds, mri_cat("images.min.", dim));
      else
	fh.dim[i].min = 0;
      if (mri_has(pgh_read_ds, mri_cat("images.extent.", dim)))
	fh.dim[i].n = mri_get_int(pgh_read_ds, mri_cat("images.extent.", dim));
      else
	fh.dim[i].n = 1;
      fh.dim[i].max = (fh.dim[i].n - 1) * fh.dim[i].stride + fh.dim[i].min;
      if (mri_has(pgh_read_ds, mri_cat("images.size.", dim)))
	fh.dim[i].size = mri_get_float(pgh_read_ds, mri_cat("images.size.", dim));
      else
	fh.dim[i].size = 0.0;
    }

  if (fh.dim[0].type != GIFT_V_DIMENSION)
    {
      /* insert the V dimension with a single index */
      for (i = fh.n_dims; i > 0; --i)
	fh.dim[i] = fh.dim[i-1];
      ++fh.n_dims;
      ++fh.n_image_dims;
      fh.dim[0].type = GIFT_V_DIMENSION;
      fh.dim[0].min = 0;
      fh.dim[0].max = 0;
      fh.dim[0].stride = 1;
      fh.dim[0].n = 1;
      fh.dim[0].size = 0.0;
    }

  if (fh.n_dims < 3 ||
      fh.n_dims > 5)
    Abort("Invalid # of dimensions in PGH dataset\n");

  if (fh.n_dims == 3)
    {
      /* insert the Z dimension with a single index */
      fh.dim[3].type = GIFT_Z_DIMENSION;
      fh.dim[3].min = 0;
      fh.dim[3].max = 0;
      fh.dim[3].stride = 1;
      fh.dim[3].n = 1;
      fh.dim[3].size = 0.0;
      fh.n_dims = 4;
    }

  if (fh.n_dims == 4)
    {
      /* insert the T dimension with a single index */
      fh.dim[4].type = GIFT_T_DIMENSION;
      fh.dim[4].min = 0;
      fh.dim[4].max = 0;
      fh.dim[4].stride = 1;
      fh.dim[4].n = 1;
      fh.dim[4].size = 0.0;
      fh.n_dims = 5;
    }

  /* check that the dimensions are in proper order */
  if (fh.dim[0].type != GIFT_V_DIMENSION ||
      fh.dim[1].type != GIFT_X_DIMENSION ||
      fh.dim[2].type != GIFT_Y_DIMENSION ||
      fh.dim[3].type != GIFT_Z_DIMENSION ||
      fh.dim[4].type != GIFT_T_DIMENSION)
    Abort("gift requires dimensions to be in vxyzt order.\n");

  fh.n_items_per_image = fh.dim[0].n * fh.dim[1].n * fh.dim[2].n;
  fh.n_images = fh.dim[3].n * fh.dim[4].n;

  fh.corrupt = (char *) malloc(fh.n_images);
  if (mri_has(pgh_read_ds, "images.corrupt") &&
      strcmp(mri_get_string(pgh_read_ds, "images.corrupt"), "[chunk]") == 0)
    {
      ptr = mri_get_chunk(pgh_read_ds, "images.corrupt", fh.n_images, 0, MRI_UNSIGNED_CHAR);
      memcpy(fh.corrupt, ptr, fh.n_images);
    }
  else
    memset(fh.corrupt, 0, fh.n_images);

  pgh_read_orig_time_min = fh.dim[4].min;
  pgh_read_orig_time_stride = fh.dim[4].stride;
  pgh_read_orig_slice_min = fh.dim[3].min;
  pgh_read_orig_slice_stride = fh.dim[3].stride;

  /* adjust the header to the specified subset of images */
  if (slice_start >= 0 && slice_start > fh.dim[3].min)
    fh.dim[3].min = slice_start;
  if (slice_end >= 0 && slice_end < fh.dim[3].max)
    fh.dim[3].max = slice_end;
  if (slice_stride > 0)
    fh.dim[3].stride = slice_stride;
  fh.dim[3].n = (fh.dim[3].max - fh.dim[3].min) / fh.dim[3].stride + 1;
  if (time_start >= 0 && time_start > fh.dim[4].min)
    fh.dim[4].min = time_start;
  if (time_end >= 0 && time_end < fh.dim[4].max)
    fh.dim[4].max = time_end;
  if (time_stride > 0)
    fh.dim[4].stride = time_stride;
  fh.dim[4].n = (fh.dim[4].max - fh.dim[4].min) / fh.dim[4].stride + 1;

  image = (char *) malloc(fh.n_items_per_image * GiftBytesPerItem(fh.data_type));
}

void
PghReadImage (int time,
	      int slice)
{
  void *ptr;
  int t, s;
  MRI_ArrayType type;

  t = (time - pgh_read_orig_time_min) / pgh_read_orig_time_stride;
  s = (slice - pgh_read_orig_slice_min) / pgh_read_orig_slice_stride;
  
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
      Abort("PghReadImage: Invalid data type\n");
      break;
    }
  ptr = mri_get_image (pgh_read_ds, t, s, type);
  memcpy(image, ptr, fh.n_items_per_image * GiftBytesPerItem(fh.data_type));
}

void
PghEndReading ()
{
  mri_close_dataset(pgh_read_ds);
  if (delete_input)
    {
      pgh_read_ds =  mri_open_dataset(pgh_read_file, MRI_MODIFY);
      mri_destroy_dataset(pgh_read_ds);
    }
}
