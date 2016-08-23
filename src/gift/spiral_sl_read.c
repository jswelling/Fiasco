/*
 *	spiral_sl_read.c - routines for reading in SPIRAL_SL files
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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fmri.h"
#include "bio.h"
#include "gift.h"
#include "spiral_sl.h"

static char rcsid[] = "$Id: spiral_sl_read.c,v 1.5 2000/11/02 20:15:18 welling Exp $";

Boolean spiral_sl_read_compressed;
int spiral_sl_read_coil = -1;
FileList spiral_sl_read_files_head;
FileList spiral_sl_read_files_tail;
int spiral_sl_read_image_size;
int spiral_sl_read_images_left;
Filename spiral_sl_read_basename;
Filename spiral_sl_read_file_name;
FileList spiral_sl_read_files;

Boolean
SpiralSLCheckFormat (Filename basename,
		     FileList files)
{
  FileListElement *fle;
  int i, j;
  char *p;

  /* for SPIRAL_SL check for file names starting with "sl",
     and followed by at least two numbers */
  for (fle = files; fle != NULL; fle = fle->next)
    {
      p = strrchr(fle->name, '/');
      if (p == NULL)
	p = fle->name;
      else
	++p;
      if (strncmp(p, "sl", 2) == 0 &&
	  sscanf(&p[2], "%d.%d", &i, &j) == 2)
	return(TRUE);
    }
  /* no sl file was found */
  return(FALSE);
}

void
SpiralSLStartReading (Filename basename,
		      FileList files)
{
  FileListElement *fle;
  Boolean first;
  int s;		/* slice */
  int c;		/* coil */
  int t;		/* time */
  int n;
  int resolution;
  int min_time, max_time;
  int min_slice, max_slice;
  Filename first_image;
  FILE *f;
  long length;
  char *p;

  min_time = 999999999;
  max_time = -999999999;
  min_slice = 999999999;
  max_slice = -999999999;
  first = TRUE;
  for (fle = files; fle != NULL; fle = fle->next)
    {
      p = strrchr(fle->name, '/');
      if (p == NULL)
	p = fle->name;
      else
	++p;

      if (spiral_sl_read_coil >= 0)
	{
	  if (sscanf(p, "sl%d.%d.%d%n", &s, &c, &t, &n) != 3 ||
	      c != spiral_sl_read_coil)
	    continue;
	}
      else
	{
	  if (sscanf(p, "sl%d.%d%n", &s, &t, &n) != 2)
	    continue;
	}
      if (p[n] != '\0' && strcmp(&p[n], ".Z") != 0)
	continue;
      if (first)
	{
	  spiral_sl_read_compressed = HasExtension(p, ".Z");
	  if (spiral_sl_read_compressed)
	    Uncompress(first_image, fle->name);
	  else
	    strcpy(first_image, fle->name);
	  if ((f = fopen(first_image, "r")) == NULL)
	    Abort("Cannot open sl file %s\n", first_image);
	  fseek(f, 0L, SEEK_END);
	  length = ftell(f);
	  fclose(f);
	  resolution = (int) sqrt((double) (length / 2));
	  spiral_sl_read_image_size = 2 * resolution * resolution;
	  if (spiral_sl_read_image_size != length)
	    Abort("Cannot determine image resolution of sl files\n");
	  if (spiral_sl_read_compressed)
	    Delete(first_image);
	  strcpy(spiral_sl_read_basename, fle->name);
	  p = strrchr(spiral_sl_read_basename, '/');
	  if (p != NULL)
	    *(p+1) = '\0';
	  else
	    spiral_sl_read_basename[0] = '\0';
	  first = FALSE;
	  }
      else
	if (spiral_sl_read_compressed != HasExtension(p, ".Z"))
	  Abort("Invalid mixture of compressed and uncompressed SPIRAL_SL files.\n");
      if (s < min_slice)
	min_slice = s;
      if (s > max_slice)
	max_slice = s;
      if (t < min_time)
	min_time = t;
      if (t > max_time)
	max_time = t;
      }
  if (first)
    /* no sl file was found */
    Abort("Can't find appropriate sl files.\n");

  fh.class = GIFT_I_SPACE;
  fh.data_type = GIFT_INT16;
  fh.n_dims = 5;
  fh.n_image_dims = 3;

  fh.dim[0].type = GIFT_V_DIMENSION;
  fh.dim[0].min = 0;
  fh.dim[0].max = 0;
  fh.dim[0].stride = 1;
  fh.dim[0].n = 1;
  fh.dim[0].size = 0.0;

  fh.dim[1].type = GIFT_X_DIMENSION;
  fh.dim[1].min = 0;
  fh.dim[1].max = resolution-1;
  fh.dim[1].stride = 1;
  fh.dim[1].n = resolution;
  fh.dim[1].size = 0.0;

  fh.dim[2].type = GIFT_Y_DIMENSION;
  fh.dim[2].min = 0;
  fh.dim[2].max = resolution-1;
  fh.dim[2].stride = 1;
  fh.dim[2].n = resolution;
  fh.dim[2].size = 0.0;

  fh.dim[3].type = GIFT_Z_DIMENSION;
  fh.dim[3].min = min_slice;
  if (slice_start >= fh.dim[3].min)
    fh.dim[3].min = slice_start;
  fh.dim[3].max = max_slice;
  if (slice_end >= 0 && slice_end < fh.dim[3].max)
    fh.dim[3].max = slice_end;
  fh.dim[3].stride = 1;
  if (slice_stride > 0)
    fh.dim[3].stride = slice_stride;
  fh.dim[3].n = (fh.dim[3].max - fh.dim[3].min) / fh.dim[3].stride + 1;
  fh.dim[3].size = 0.0;

  fh.dim[4].type = GIFT_T_DIMENSION;
  fh.dim[4].min = min_time;
  if (time_start >= fh.dim[4].min)
    fh.dim[4].min = time_start;
  fh.dim[4].max = max_time;
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

  fh.corrupt = (char *) malloc(fh.n_images);
  memset(fh.corrupt, 0, fh.n_images);


  if (spiral_sl_read_compressed)
    {
      /* build a list of the files so we can uncompress several at a time */
      spiral_sl_read_files_head = NULL;
      spiral_sl_read_files_tail = NULL;
      for (t = fh.dim[4].min; t <= fh.dim[4].max; t += fh.dim[4].stride)
	for (s = fh.dim[3].min; s <= fh.dim[3].max; s += fh.dim[3].stride)
	  {
	    if (spiral_sl_read_coil >= 0)
	      sprintf(spiral_sl_read_file_name, "%ssl%1d.%1d.%.3d.Z",
		      spiral_sl_read_basename, s, spiral_sl_read_coil, t);
	    else
	      sprintf(spiral_sl_read_file_name, "%ssl%1d.%.3d.Z",
		      spiral_sl_read_basename, s, t);
	    AppendToFileList(&spiral_sl_read_files_head,
			     &spiral_sl_read_files_tail,
			     spiral_sl_read_file_name,
			     spiral_sl_read_image_size);
	  }
      input = NULL;
      spiral_sl_read_images_left = 0;
      spiral_sl_read_files = spiral_sl_read_files_head;
    }

  image = (char *) malloc(spiral_sl_read_image_size);
}

void
SpiralSLReadImage (int time,
		   int slice)
{
  if (spiral_sl_read_compressed)
    {
      if (spiral_sl_read_images_left == 0)
	{
	  if (input != NULL)
	    {
	      fclose(input);
	      Delete(spiral_sl_read_file_name);
	    }
	  spiral_sl_read_images_left = UncompressBatch(spiral_sl_read_file_name, &spiral_sl_read_files) /
	    spiral_sl_read_image_size;
	  if ((input = fopen(spiral_sl_read_file_name, "r")) == NULL)
	    Abort("Can't open input file %s\n", spiral_sl_read_file_name);
	}
      ReadImage(time, slice);
      --spiral_sl_read_images_left;
    }
  else
    {
      if (spiral_sl_read_coil >= 0)
	sprintf(spiral_sl_read_file_name, "%ssl%1d.%1d.%.3d",
		spiral_sl_read_basename, slice, spiral_sl_read_coil, time);
      else
	sprintf(spiral_sl_read_file_name, "%ssl%1d.%.3d",
		spiral_sl_read_basename, slice, time);
      if ((input = fopen(spiral_sl_read_file_name, "r")) == NULL)
	Abort("Can't open input file %s\n", spiral_sl_read_file_name);
      ReadImage(time, slice);
      fclose(input);
    }
}

void
SpiralSLEndReading ()
{
  FileListElement *fle, *nfle;
  char com[sizeof(Filename) + 64];

  if (spiral_sl_read_compressed)
    {
      if (input != NULL)
	{
	  fclose(input);
	  Delete(spiral_sl_read_file_name);
	}
      for (fle = spiral_sl_read_files_head; fle != NULL; fle = nfle)
	{
	  nfle = fle->next;
	  free(fle);
	}
    }
  if (delete_input)
    {
      sprintf(com, "rm -f %ssl*", spiral_sl_read_basename);
      system(com);
    }
}
