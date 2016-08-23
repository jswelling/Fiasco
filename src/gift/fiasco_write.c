/*
 *	fiasco_write.c - routines for writing out FIASCO style files
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
 */

#include <stdio.h>
#include <string.h>
#include "fmri.h"
#include "bio.h"
#include "gift.h"
#include "fiasco.h"

static char rcsid[] = "$Id: fiasco_write.c,v 1.4 2000/11/02 20:15:44 welling Exp $";

Filename fiasco_write_filename;
int fiasco_write_typesize;
int fiasco_write_type;
int fiasco_write_vec;
int fiasco_write_dim;
int fiasco_write_dims[10];


void
FiascoStartWriting (Filename basename)
{
  Filename header;
  char buf[1024];
  FILE *f;
  int i;
  char filename[512];
  char *p;

  sprintf(header, "%s.hdr", basename);
  sprintf(fiasco_write_filename, "%s.dat", basename);

  if ((f = fopen(header, "w")) == NULL)
    Abort("Can't open header file %s for writing\n", header);

  switch (fh.data_type)
    {
    case GIFT_UINT8:
      fiasco_write_typesize = 1;
      fiasco_write_type = FIASCO_MRI_CHAR;
      break;
    case GIFT_INT16:
      fiasco_write_typesize = 2;
      fiasco_write_type = FIASCO_MRI_SHORT;
      break;
    case GIFT_INT32:
      fiasco_write_typesize = 4;
      fiasco_write_type = FIASCO_MRI_LONG;
      break;
    case GIFT_FLOAT32:
      fiasco_write_typesize = 4;
      fiasco_write_type = FIASCO_MRI_FLOAT;
      break;
    case GIFT_FLOAT64:
      fiasco_write_typesize = 8;
      fiasco_write_type = FIASCO_MRI_DOUBLE;
      break;
    default:
      Abort("Cannot convert intermediate type to fiasco_write_type.\n");
    }

  fiasco_write_vec = fh.dim[0].n;
  fiasco_write_dim = 0;
  for (i = 1; i < fh.n_dims; ++i)
    fiasco_write_dims[fiasco_write_dim++] = fh.dim[i].n;

  bio_error = 0;
  FWrInt32(f, fiasco_write_type);
  FWrInt32(f, fiasco_write_dim);
  FWrInt32(f, fiasco_write_vec);
  FWrInt32Array(f, fiasco_write_dims, fiasco_write_dim);

  memset(filename, 0, 512);
  p = strrchr(fiasco_write_filename, '/');
  strcpy(filename, p != NULL ? p+1 : fiasco_write_filename);

  if (bio_error ||
      fwrite(filename, 1, 512, f) != 512 ||
      fwrite(fh.corrupt, 1, fh.n_images, f) != fh.n_images)
    Abort("Cannot write to file %s\n", header);
  fclose(f);
  if (compress_output)
    Compress(header, 500 + fh.n_images);

  /* open the data file */
  if ((output = fopen(fiasco_write_filename, "w")) == NULL)
    Abort("Can't open data file %s for writing\n", fiasco_write_filename);
}

void
FiascoWriteImage (int time,
		     int slice)
{
  /* since the different data types are all supported, we may just
     write the image array directly out */
  WriteImage(time, slice);
}

void
FiascoEndWriting ()
{
  fclose(output);
  if (compress_output)
    Compress(fiasco_write_filename, fh.n_images*fh.n_items_per_image*GiftBytesPerItem(fh.data_type));
}
