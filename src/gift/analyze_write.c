/*
 *	analyze_write.c - routines for writing out ANALYZE files
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
#include <sys/stat.h>
#include "fmri.h"
#include "bio.h"
#include "gift.h"
#include "analyze.h"

static char rcsid[] = "$Id: analyze_write.c,v 1.6 2000/11/02 20:15:55 welling Exp $";

Filename analyze_write_basename;
Filename analyze_write_name;
Boolean analyze_write_split_times = TRUE;
unsigned char *analyze_write_hdr;
int analyze_write_type;
int analyze_write_file_size;
int analyze_write_time;
int analyze_write_sequence;

void
AnalyzeStartWriting (Filename basename)
{
  int bp;
  int glmin, glmax;
  char *pos;

  if (fh.dim[0].n != 1)
    if (fh.data_type != GIFT_FLOAT32 || fh.dim[0].n != 2)
      Abort("ANALYZE format cannot handle more than 1 item per voxel except\nfor complex values\n");

  analyze_write_hdr = (unsigned char *) malloc(ANALYZE_HEADER_SIZE);
  bio_big_endian_output = big_endian_output;
  BWrInt32(&analyze_write_hdr[ANALYZE_SIZEOF_HDR], ANALYZE_HEADER_SIZE);
  BWrInt32(&analyze_write_hdr[ANALYZE_EXTENTS], 16384);
  analyze_write_hdr[ANALYZE_REGULAR] = 'r';
  BWrInt16(&analyze_write_hdr[ANALYZE_N_DIMS], 4);
  BWrInt16(&analyze_write_hdr[ANALYZE_X],
	   (fh.dim[1].max - fh.dim[1].min) / fh.dim[1].stride + 1);
  BWrInt16(&analyze_write_hdr[ANALYZE_Y],
	   (fh.dim[2].max - fh.dim[2].min) / fh.dim[2].stride + 1);
  BWrInt16(&analyze_write_hdr[ANALYZE_Z],
	   (fh.dim[3].max - fh.dim[3].min) / fh.dim[3].stride + 1);
  BWrInt16(&analyze_write_hdr[ANALYZE_T],
	   analyze_write_split_times ?
	   1 :
	   (fh.dim[4].max - fh.dim[4].min) / fh.dim[4].stride + 1);
  switch (fh.data_type)
    {
    case GIFT_UINT8:
      analyze_write_type = ANALYZE_DATATYPE_UINT8;
      bp = 8;
      glmin = 0;
      glmax = 255;
      break;

    case GIFT_INT16:
      analyze_write_type = ANALYZE_DATATYPE_INT16;
      bp = 16;
      glmin = 0;
      glmax = 32767;
      break;

    case GIFT_INT32:
      analyze_write_type = ANALYZE_DATATYPE_INT32;
      bp = 32;
      glmin = 0;
      glmax = 2147483647;
      break;

    case GIFT_FLOAT32:
      if (fh.dim[0].n == 2)
	{
	  analyze_write_type = ANALYZE_DATATYPE_COMPLEX64;
	  bp = 64;
	}
      else
	{
	  analyze_write_type = ANALYZE_DATATYPE_FLOAT32;
	  bp = 32;
	}
      glmin = 0;
      glmax = 2147483647;
      break;

    case GIFT_FLOAT64:
      analyze_write_type = ANALYZE_DATATYPE_FLOAT64;
      bp =  64;
      glmin = 0;
      glmax = 2147483647;
      break;

    default:
      Abort("Invalid intermediate datatype.\n");
    }

  BWrInt16(&analyze_write_hdr[ANALYZE_DATATYPE], analyze_write_type);
  BWrInt16(&analyze_write_hdr[ANALYZE_BITPIX], bp);
  BWrFloat32(&analyze_write_hdr[ANALYZE_VOXEL_X_SIZE],
	     fh.dim[1].size);
  BWrFloat32(&analyze_write_hdr[ANALYZE_VOXEL_Y_SIZE],
	     fh.dim[2].size);
  BWrFloat32(&analyze_write_hdr[ANALYZE_VOXEL_Z_SIZE],
	     fh.dim[3].size);
  BWrFloat32(&analyze_write_hdr[ANALYZE_T_STEP],
	     fh.dim[4].size);
  BWrInt32(&analyze_write_hdr[ANALYZE_GLMAX], glmax);
  BWrInt32(&analyze_write_hdr[ANALYZE_GLMIN], glmin);

  output = NULL;
  strcpy(analyze_write_basename, basename);
  if ((pos = strrchr(analyze_write_basename, '/')) != NULL)
    {
      *pos = '\0';
      mkdir(analyze_write_basename, 0755);
      *pos = '/';
      if (*(pos+1) == '\0')
	strcat(analyze_write_basename, "func");
    }
  if (!analyze_write_split_times)
    {
      sprintf(analyze_write_name, "%s.hdr", analyze_write_basename);
      if ((output = fopen(analyze_write_name, "w")) == NULL)
	Abort("Can't open %s for writing.\n", analyze_write_name);
      if (fwrite(analyze_write_hdr, 1, ANALYZE_HEADER_SIZE, output) != ANALYZE_HEADER_SIZE)
	Abort("Can't write ANALYZE header %s\n", analyze_write_name);
      fclose(output);
      if (compress_output)
	Compress(analyze_write_name, ANALYZE_HEADER_SIZE);
      sprintf(analyze_write_name, "%s.img", analyze_write_basename);
      if ((output = fopen(analyze_write_name, "w")) == NULL)
	Abort("Can't open %s for writing.\n", analyze_write_name);
      analyze_write_file_size = fh.n_images * fh.n_items_per_image * GiftBytesPerItem(fh.data_type);
    }
  else
    {
      analyze_write_file_size = fh.dim[3].n * fh.n_items_per_image * GiftBytesPerItem(fh.data_type);
      analyze_write_time = -999999999;
      analyze_write_sequence = 1;
    }
}

void
AnalyzeWriteImage (int time,
		   int slice)
{
  int i;

  if (analyze_write_split_times && time != analyze_write_time)
    {
      if (output != NULL)
	{
	  fclose(output);
	  if (compress_output)
	    Compress(analyze_write_name, analyze_write_file_size);
	}
      sprintf(analyze_write_name, "%s.%.5d.hdr", analyze_write_basename, analyze_write_sequence);
      if ((output = fopen(analyze_write_name, "w")) == NULL)
	Abort("Can't open %s for writing.\n", analyze_write_name);
      if (fwrite(analyze_write_hdr, 1, ANALYZE_HEADER_SIZE, output) != ANALYZE_HEADER_SIZE)
	Abort("Can't write ANALYZE header %s\n", analyze_write_name);
      fclose(output);
      if (compress_output)
	Compress(analyze_write_name, ANALYZE_HEADER_SIZE);
      sprintf(analyze_write_name, "%s.%.5d.img", analyze_write_basename, analyze_write_sequence);
      if ((output = fopen(analyze_write_name, "w")) == NULL)
	Abort("Can't open %s for writing.\n", analyze_write_name);
      analyze_write_time = time;
      ++analyze_write_sequence;
    }
  WriteImage(time, slice);
}

void
AnalyzeEndWriting ()
{
  if (output != NULL)
    {
      fclose(output);
      if (compress_output)
	Compress(analyze_write_name, analyze_write_file_size);
    }
  free(analyze_write_hdr);
}
