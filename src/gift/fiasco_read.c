/*
 *	fiasco_read.c - routines for reading in FIASCO style files
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
 *		1/96	Written by Greg Hood (PSC)
 *		4/96	Updated many function and variable names (ghood)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fmri.h"
#include "gift.h"
#include "bio.h"
#include "fiasco.h"

static char rcsid[] = "$Id: fiasco_read.c,v 1.5 2000/11/02 20:15:50 welling Exp $";

Filename fiasco_read_header;
char fiasco_read_filename[512];
int fiasco_read_type;
int fiasco_read_typesize;
int fiasco_read_ndims;
int fiasco_read_vec;
int fiasco_read_dims[10];
int fiasco_read_nim;
int fiasco_read_volume_size;
int fiasco_read_slice_size;
int fiasco_read_size;
Boolean fiasco_read_compressed;
Filename fiasco_read_file;

Boolean
FiascoCheckFormat (Filename basename,
		   FileList files)
{
  FileList fl;
  int len;
  unsigned char buf[FIASCO_Z];
  Filename header;
  FILE *f;
  int i;
  Boolean compressed;
  Boolean fiasco_format;

  /* for FIASCO check that the header has reasonable dimensions and types */
  header[0] = '\0';
  for (fl = files; fl != NULL; fl = fl->next)
    if (HasExtension(fl->name, ".hdr"))
      {
	compressed = FALSE;
	strcpy(header, fl->name);
	break;
      }
    else if (HasExtension(fl->name, ".hdr.Z"))
      {
	compressed = TRUE;
	Uncompress(header, fl->name);
	break;
      }
  if (fl == NULL)
    /* no .hdr file was found */
    return(FALSE);

  if ((f = fopen(header, "r")) == NULL)
    Abort("Can't open header file %s\n", header);

  
  fiasco_format = TRUE;
  if (fread(buf, 1, FIASCO_Z, f) != FIASCO_Z ||
      (i = BRdInt32(&buf[FIASCO_TYPE])) < FIASCO_MRI_CHAR ||
       i > FIASCO_MRI_SHORT ||
      (i = BRdInt32(&buf[FIASCO_NDIMS])) < 2 || i > 4 ||
      (i = BRdInt32(&buf[FIASCO_VEC])) < 1 || i > 2 ||
      (i = BRdInt32(&buf[FIASCO_X])) < 1 || i > 10000 ||
      (i = BRdInt32(&buf[FIASCO_Y])) < 1 || i > 10000)
    fiasco_format = FALSE;

  fclose(f);
  if (compressed)
    Delete(header);
  return(fiasco_format);
}
      
void
FiascoStartReading (Filename basename,
		    FileList files)
{
  FileList fl;
  int len;
  Filename header;
  FILE *f;
  int nim;
  int i, j;
  int index;
  char *missing;
  Boolean compressed;
  char filename[512];
  char *p;

  /* check that only one header file exists */
  header[0] = '\0';
  for (fl = files; fl != NULL; fl = fl->next)
    if (HasExtension(fl->name, ".hdr"))
      {
	if (header[0] != '\0')
	  Abort("More than one .hdr file was found beginning with %s\n",
		basename);
	strcpy(fiasco_read_header, fl->name);
	strcpy(header, fl->name);
	strcpy(fiasco_read_filename, fl->name);
	compressed = FALSE;
      }
    else if (HasExtension(fl->name, ".hdr.Z"))
      {
	if (header[0] != '\0')
	  Abort("More than one .hdr file was found beginning with %s\n",
		basename);
	strcpy(fiasco_read_header, fl->name);
	Uncompress(header, fl->name);
	compressed = TRUE;
      }

  if (header[0] == '\0')
    /* no .hdr file was found */
    Abort("Cannot find header file!\n");

  if ((f = fopen(header, "r")) == NULL)
    Abort("Can't open header file %s\n", header);

  bio_error = FALSE;
  fiasco_read_type = FRdInt32(f);
  fiasco_read_ndims = FRdInt32(f);
  fiasco_read_vec = FRdInt32(f);
  if (bio_error)
    Abort("Can't read header.\n");
  switch (fiasco_read_type)
    {
    case FIASCO_MRI_CHAR:
      fiasco_read_typesize = 1;
      fh.data_type = GIFT_UINT8;
      break;
    case FIASCO_MRI_LONG:
      fiasco_read_typesize = 4;
      fh.data_type = GIFT_INT32;
      break;
    case FIASCO_MRI_FLOAT:
      fiasco_read_typesize = 4;
      fh.data_type = GIFT_FLOAT32;
      break;
    case FIASCO_MRI_DOUBLE:
      fiasco_read_typesize = 8;
      fh.data_type = GIFT_FLOAT64;
      break;
    case FIASCO_MRI_SHORT:
      fiasco_read_typesize = 2;
      fh.data_type = GIFT_INT16;
      break;
    default:
      Abort("Invalid type field in header.\n");
    }
  if (fiasco_read_ndims <= 1 || fiasco_read_ndims > 4)
    Abort("Invalid # of dimensions in header %s.\n", header);
  FRdInt32Array(f, fiasco_read_dims, fiasco_read_ndims);
  if (bio_error)
    Abort("Can't read header dims.\n");
  if (fiasco_read_ndims < 3)
    fiasco_read_dims[2] = 1;
  if (fiasco_read_ndims < 4)
    fiasco_read_dims[3] = 1;
  nim = 1;
  for (i = 2; i < fiasco_read_ndims; ++i)
    nim *= fiasco_read_dims[i];
  if (fiasco_read_vec < 1)
    Abort("vec field in header must be greater than 0.\n");

  if (fread(filename, 1, 512, f) != 512)
    Abort("Can't read file name.\n");
  if (MAX_FILENAME_LENGTH <= 512)
    filename[MAX_FILENAME_LENGTH-1] = '\0';

  if (strchr(filename, '/') != NULL)
    strcpy(fiasco_read_filename, filename);
  else
    {
      if ((p = strrchr(fiasco_read_filename, '/')) != NULL)
	*(p+1) = '\0';
      strcat(fiasco_read_filename, filename);
    }

  /* read in missing array */
  missing = (char *) malloc(nim);
  if (fread(missing, 1, nim, f) != nim)
    Abort("Can't read in missing array.\n");
  fclose(f);
  if (compressed)
    Delete(header);

  /* transfer to canonical header */
  fh.class = GIFT_I_SPACE;
  fh.n_dims = 5;
  fh.n_image_dims = 3;

  fh.dim[0].type = GIFT_V_DIMENSION;
  fh.dim[0].min = 0;
  fh.dim[0].max = fiasco_read_vec - 1;
  fh.dim[0].stride = 1;
  fh.dim[0].n = fiasco_read_vec;
  fh.dim[0].size = 0.0;
  fh.dim[1].type = GIFT_X_DIMENSION;
  fh.dim[1].min = 0;
  fh.dim[1].max = fiasco_read_dims[0] - 1;
  fh.dim[1].stride = 1;
  fh.dim[1].n = fiasco_read_dims[0];
  fh.dim[1].size = 0.0;
  fh.dim[2].type = GIFT_Y_DIMENSION;
  fh.dim[2].min = 0;
  fh.dim[2].max = fiasco_read_dims[1] - 1;
  fh.dim[2].stride = 1;
  fh.dim[2].n = fiasco_read_dims[1];
  fh.dim[2].size = 0.0;
  fh.dim[3].type = GIFT_Z_DIMENSION;
  fh.dim[3].min = 0;
  if (slice_start >= 0)
    fh.dim[3].min = slice_start;
  fh.dim[3].max = fiasco_read_dims[2] - 1;
  if (slice_end >= 0 && slice_end < fh.dim[3].max)
    fh.dim[3].max = slice_end;
  fh.dim[3].stride = 1;
  if (slice_stride > 0)
    fh.dim[3].stride = slice_stride;
  fh.dim[3].n = (fh.dim[3].max - fh.dim[3].min) / fh.dim[3].stride + 1;
  fh.dim[3].size = 0.0;
  fh.dim[4].type = GIFT_T_DIMENSION;
  fh.dim[4].min = 0;
  if (time_start >= 0)
    fh.dim[4].min = time_start;
  fh.dim[4].max = fiasco_read_dims[3] - 1;
  if (time_end >= 0 && time_end < fh.dim[4].max)
    fh.dim[4].max = time_end;
  fh.dim[4].stride = 1;
  if (time_stride > 0)
    fh.dim[4].stride = time_stride;
  fh.dim[4].n = (fh.dim[4].max - fh.dim[4].min) / fh.dim[4].stride + 1;
  fh.dim[4].size = 0.0;

  if (fh.dim[3].max < fh.dim[3].min)
    Abort("No slices found in selected range.\n");
  if (fh.dim[4].max < fh.dim[4].min)
    Abort("No times found in selected range.\n");
  fh.n_images = fh.dim[3].n * fh.dim[4].n;
  fh.n_items_per_image = fh.dim[0].n * fh.dim[1].n * fh.dim[2].n;

  fiasco_read_slice_size = fiasco_read_typesize * fiasco_read_vec *
    fiasco_read_dims[0] * fiasco_read_dims[1];
  fiasco_read_volume_size = fiasco_read_slice_size * fiasco_read_dims[2];

  fh.corrupt = (char *) malloc(fh.n_images);
  index = 0;
  for (i = fh.dim[4].min; i <= fh.dim[4].max; i += fh.dim[4].stride)
    for (j = fh.dim[3].min; j <= fh.dim[3].max; j += fh.dim[3].stride)
      fh.corrupt[index++] = missing[i * fiasco_read_dims[2] + j];
  free(missing);

  /* open the data file */
  if (HasExtension(fiasco_read_filename, ".Z"))
    {
      Uncompress(fiasco_read_file, fiasco_read_filename);
      fiasco_read_compressed = TRUE;
    }
  else
    {
      StringCopy(fiasco_read_file, fiasco_read_filename, MAX_FILENAME_LENGTH);
      fiasco_read_compressed = FALSE;
    }
		 
  if ((input = fopen(fiasco_read_file, "r")) == NULL)
    Abort("Can't open data file %s\n", fiasco_read_file);
  image = (char *) malloc(fh.n_items_per_image * GiftBytesPerItem(fh.data_type));
}

void
FiascoReadImage (int time,
		 int slice)
{
  long offset;

  /* calculate the offset and seek to the correct image */
  offset = time*fiasco_read_volume_size + slice*fiasco_read_slice_size;
  if (fseek(input, offset, SEEK_SET) != 0)
    Abort("Cannot seek to image (time %d, slice %d) in input file.\n",
	  time, slice);
  /* since the different data types are all supported, we may just 
     read directly into the image array */
  ReadImage(time, slice);
}

void
FiascoEndReading ()
{
  char com[sizeof(Filename) + 64];

  /* close the file */
  fclose(input);
  if (fiasco_read_compressed)
    Delete(fiasco_read_file);
  if (delete_input)
    {
      sprintf(com, "rm -f %s %s",
	      fiasco_read_header,
	      fiasco_read_filename);
      system(com);
    }

}

