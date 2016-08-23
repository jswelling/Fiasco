/*
 *	analyze_read.c - routines for reading in ANALYZE files
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
#include "fmri.h"
#include "bio.h"
#include "gift.h"
#include "analyze.h"

static char rcsid[] = "$Id: analyze_read.c,v 1.6 2000/11/02 20:16:01 welling Exp $";

Boolean analyze_read_numbered;
Filename analyze_read_basename;
int analyze_read_dt;		/* datatype */
int analyze_read_x;
int analyze_read_y;
int analyze_read_z;
int analyze_read_t;
float analyze_read_vx;
float analyze_read_vy;
float analyze_read_vz;
float analyze_read_vt;
unsigned char *analyze_read_tmp_buf;
FileList analyze_read_hdr_files_head;
FileList analyze_read_hdr_files_tail;
FileList analyze_read_img_files_head;
FileList analyze_read_img_files_tail;
Boolean analyze_read_headers_compressed;
Boolean analyze_read_images_compressed;
Filename analyze_read_hdr_file_name;
Filename analyze_read_img_file_name;
int analyze_read_file_size;
int analyze_read_image_size;
int analyze_read_time;
int analyze_read_times_left;
FileList analyze_read_image_files;

Boolean
AnalyzeCheckFormat (Filename basename,
		    FileList files)
{
  Boolean analyze_read_format;
  FileListElement *fle;
  Boolean compressed;
  Filename header;
  FILE *f;
  unsigned char buf[ANALYZE_HEADER_SIZE];
  
  /* for ANALYZE check that the header has ANALYZE_SIZEOF_HDR containing 348,
     ANALYZE_EXTENTS containing 16384, and ANALYZE_REGULAR containing 'r' */

  for (fle = files; fle != NULL; fle = fle->next)
    if (HasExtension(fle->name, ".hdr"))
      {
	compressed = FALSE;
	strcpy(header, fle->name);
	break;
      }
    else if (HasExtension(fle->name, ".hdr.Z"))
      {
	compressed = TRUE;
	Uncompress(header, fle->name);
	break;
      }
  if (fle == NULL)
    /* no .hdr file was found */
    return(FALSE);

  if ((f = fopen(header, "r")) == NULL)
    Abort("Can't open header file %s\n", header);

  analyze_read_format = TRUE;
  bio_big_endian_input = big_endian_input;
  if (fread(buf, 1, ANALYZE_HEADER_SIZE, f) != ANALYZE_HEADER_SIZE ||
      BRdInt32(&buf[ANALYZE_SIZEOF_HDR]) != ANALYZE_HEADER_SIZE ||
      BRdInt32(&buf[ANALYZE_EXTENTS]) != 16384 ||
      buf[ANALYZE_REGULAR] != 'r')
    analyze_read_format = FALSE;

  fclose(f);
  if (compressed)
    Delete(header);
  return(analyze_read_format);
}

void
AnalyzeStartReading (Filename basename,
		     FileList files)
{
  Boolean first;
  FileListElement *fle;
  FILE *f;
  unsigned char buf[ANALYZE_HEADER_SIZE];
  int x, y, z, t, dt;
  float vx, vy, vz, vt;
  int len;
  FileList fl;
  int n;
  int i;
  Filename header;
  int min_time, max_time;
  int ts;
  Boolean compressed;
  Boolean numbered;
  int analyze_read_img_file_size;
  

  /* first go through and make sure all the header files
     are either uncompressed or compressed, but not a
     mixture of both; also collect information about
     the minimum and maximum times represented by
     the files */
  first = TRUE;
  min_time = 999999999;
  max_time = -999999999;
  for (fle = files; fle != NULL; fle = fle->next)
    if ((compressed = HasExtension(fle->name, ".hdr.Z")) || HasExtension(fle->name, ".hdr"))
      {
	numbered = FALSE;
	len = strlen(fle->name);
	i = len - (compressed ? 7 : 5);
	if (i >= 0 && fle->name[i] >= '0' && fle->name[i] <= '9')
	  {
	    numbered = TRUE;
	    while (i > 0 && fle->name[i-1] >= '0' && fle->name[i-1] <= '9')
	      --i;
	    sscanf(&(fle->name[i]), "%d", &t);
	  }
	if (first)
	  {
	    strcpy(analyze_read_basename, fle->name);
	    if (numbered)
	      analyze_read_basename[i] = '\0';
	    else
	      analyze_read_basename[i+1] = '\0';
	    analyze_read_numbered = numbered;
	    analyze_read_headers_compressed = compressed;
	    strcpy(analyze_read_img_file_name, fle->name);
	    analyze_read_img_file_name[strlen(analyze_read_img_file_name) - (compressed ? 5 : 3)] = '\0';
	    strcat(analyze_read_img_file_name, "img");
	    if ((f = fopen(analyze_read_img_file_name, "r")) != NULL)
	      {
		fclose(f);
		analyze_read_images_compressed = FALSE;
	      }
	    else
	      {
		strcat(analyze_read_img_file_name, ".Z");
		if ((f = fopen(analyze_read_img_file_name, "r")) != NULL)
		  {
		    fclose(f);
		    analyze_read_images_compressed = TRUE;
		  }
		else
		  Abort("Cannot find corresponding .img files.\n");
	      }
	    first = FALSE;
	  }
	else
	  {
	    if (compressed != analyze_read_headers_compressed)
	      Abort("Invalid mixture of compressed and uncompressed ANALYZE header files.\n");
	    if (!numbered || !analyze_read_numbered)
	      Abort("Only one unnumbered ANALYZE file can be converted at a time.\n");
	  }
	if (analyze_read_numbered)
	  {
	    if (t < min_time)
	      min_time = t;
	    if (t > max_time)
	      max_time = t;
	  }
      }
  if (first)
    Abort("Cannot find any ANALYZE header files\n");
    
  /* make a list of all the header files we will have to use */
  analyze_read_hdr_files_head = NULL;
  analyze_read_hdr_files_tail = NULL;
  if (analyze_read_numbered)
    {
      if (time_start >= 0 && time_start > min_time)
	min_time = time_start;
      if (time_end >= 0 && time_end < max_time)
	max_time = time_end;
      if (min_time > max_time)
	Abort("Cannot find ANALYZE header files in desired time range.\n");
      if (time_stride > 0)
	ts = time_stride;
      else
	ts = 1;

      for (t = min_time; t <= max_time; t += ts)
	{
	  sprintf(analyze_read_hdr_file_name, "%s%.5d.hdr%s",
		  analyze_read_basename, t,
		  analyze_read_headers_compressed ? ".Z" : "");
	  AppendToFileList(&analyze_read_hdr_files_head,
			   &analyze_read_hdr_files_tail,
			   analyze_read_hdr_file_name,
			   ANALYZE_HEADER_SIZE);
	}
    }
  else
    {
      sprintf(analyze_read_hdr_file_name, "%s.hdr%s",
	      analyze_read_basename, 
	      analyze_read_headers_compressed ? ".Z" : "");
      AppendToFileList(&analyze_read_hdr_files_head,
		       &analyze_read_hdr_files_tail,
		       analyze_read_hdr_file_name,
		       ANALYZE_HEADER_SIZE);
    }
		     

  /* go through header files and make sure they are consistent */
  first = TRUE;
  fl = analyze_read_hdr_files_head;
  while (fl != NULL)
    {
      if (analyze_read_headers_compressed)
	n = UncompressBatch(header, &fl) / ANALYZE_HEADER_SIZE;
      else
	{
	  strcpy(header, fl->name);
	  fl = fl->next;
	  n = 1;
	}
      if ((f = fopen(header, "r")) == NULL)
	Abort("Can't open header file %s\n", header);
      bio_big_endian_input = big_endian_input;

      for (i = 0; i < n; ++i)
	{
	  if (fread(buf, 1, ANALYZE_HEADER_SIZE, f) != ANALYZE_HEADER_SIZE)
	    Abort("Can't read header in file %s\n", fl->name);

	  if (BRdInt32(&buf[ANALYZE_SIZEOF_HDR]) != ANALYZE_HEADER_SIZE ||
	      BRdInt32(&buf[ANALYZE_EXTENTS]) != 16384 ||
	      buf[ANALYZE_REGULAR] != 'r')
	    Abort("Invalid header detected in file %s\n", fl->name);

	  if (BRdInt16(&buf[ANALYZE_N_DIMS]) != 4)
	    Abort("Unimplemented: ANALYZE file with n_dims != 4\n");

	  x = BRdInt16(&buf[ANALYZE_X]);
	  y = BRdInt16(&buf[ANALYZE_Y]);
	  z = BRdInt16(&buf[ANALYZE_Z]);
	  t = BRdInt16(&buf[ANALYZE_T]);
	  dt = BRdInt16(&buf[ANALYZE_DATATYPE]);
	  vx = BRdFloat32(&buf[ANALYZE_VOXEL_X_SIZE]);
	  vy = BRdFloat32(&buf[ANALYZE_VOXEL_Y_SIZE]);
	  vz = BRdFloat32(&buf[ANALYZE_VOXEL_Z_SIZE]);
	  vt = BRdFloat32(&buf[ANALYZE_T_STEP]);

	  if (first)
	    {
	      if (analyze_read_numbered && t != 1)
		Abort("Numbered ANALYZE files must have only time in them.\n");

	      analyze_read_dt = dt;
	      analyze_read_x = x;
	      analyze_read_y = y;
	      analyze_read_z = z;
	      analyze_read_t = t;
	      analyze_read_vx = vx;
	      analyze_read_vy = vy;
	      analyze_read_vz = vz;
	      analyze_read_vt = vt;
	      analyze_read_image_size = x * y;
	      switch (dt)
		{
		case ANALYZE_DATATYPE_INT16:
		  analyze_read_image_size *= 2;
		  break;
		case ANALYZE_DATATYPE_INT32:
		case ANALYZE_DATATYPE_FLOAT32:
		  analyze_read_image_size *= 4;
		  break;
		case ANALYZE_DATATYPE_COMPLEX64:
		case ANALYZE_DATATYPE_FLOAT64:
		  analyze_read_image_size *= 8;
		  break;
		default:
		  break;
		}
	      analyze_read_img_file_size = t * z * analyze_read_image_size;
	      first = FALSE;
	    }
	  else
	    {
	      if (dt != analyze_read_dt ||
		  x != analyze_read_x ||
		  y != analyze_read_y ||
		  z != analyze_read_z ||
		  t != analyze_read_t ||
		  vx != analyze_read_vx ||
		  vy != analyze_read_vy ||
		  vz != analyze_read_vz ||
		  vt != analyze_read_vt)
		Abort("Incompatible ANALYZE headers.\n");
	    }
	}

      fclose(f);

      if (analyze_read_headers_compressed)
	Delete(header);
    }

  /* make a list of all the image files we will have to use */
  analyze_read_img_files_head = NULL;
  analyze_read_img_files_tail = NULL;
  if (analyze_read_numbered)
    for (t = min_time; t <= max_time; t += ts)
      {
	sprintf(analyze_read_img_file_name, "%s%.5d.img%s",
		analyze_read_basename, t,
		analyze_read_images_compressed ? ".Z" : "");
	AppendToFileList(&analyze_read_img_files_head,
			 &analyze_read_img_files_tail,
			 analyze_read_img_file_name,
			 analyze_read_img_file_size);
      }
  else
    {
      min_time = 1;
      if (time_start >= 0 && time_start > min_time)
	min_time = time_start;
      max_time = analyze_read_t;
      if (time_end >= 0 && time_end < max_time)
	max_time = time_end;
      if (time_stride > 0)
	ts = time_stride;
      else
	ts = 1;
      if (min_time > max_time)
	Abort("Cannot find ANALYZE images in desired time range.\n");
      sprintf(analyze_read_img_file_name, "%s.img%s",
	      analyze_read_basename, 
	      analyze_read_images_compressed ? ".Z" : "");
      AppendToFileList(&analyze_read_img_files_head,
		       &analyze_read_img_files_tail,
		       analyze_read_img_file_name,
		       analyze_read_img_file_size);
    }

  fh.class = GIFT_I_SPACE;
  fh.n_dims = 5;
  fh.n_image_dims = 3;

  fh.dim[0].type = GIFT_V_DIMENSION;
  fh.dim[0].min = 0;
  fh.dim[0].max = 0;
  fh.dim[0].stride = 1;
  fh.dim[0].n = 1;
  fh.dim[0].size = 0.0;

  switch (analyze_read_dt)
    {
    case ANALYZE_DATATYPE_UINT8:
      fh.data_type = GIFT_INT16;
      break;
    case ANALYZE_DATATYPE_INT16:
      fh.data_type = GIFT_INT16;
      break;
    case ANALYZE_DATATYPE_INT32:
      fh.data_type = GIFT_INT32;
      break;
    case ANALYZE_DATATYPE_FLOAT32:
      fh.data_type = GIFT_FLOAT32;
      break;
    case ANALYZE_DATATYPE_COMPLEX64:
      fh.data_type = GIFT_FLOAT32;
      fh.dim[0].max = 1;
      fh.dim[0].n = 2;
      break;
    case ANALYZE_DATATYPE_FLOAT64:
      fh.data_type = GIFT_FLOAT64;
      break;
    default:
      Abort("Cannot work with ANALYZE datatype (%d)\n",
	    analyze_read_dt);
    }

  fh.dim[1].type = GIFT_X_DIMENSION;
  fh.dim[1].min = 0;
  fh.dim[1].max = analyze_read_x - 1;
  fh.dim[1].stride = 1;
  fh.dim[1].n = analyze_read_x;
  fh.dim[1].size = vx;

  fh.dim[2].type = GIFT_Y_DIMENSION;
  fh.dim[2].min = 0;
  fh.dim[2].max = analyze_read_y - 1;
  fh.dim[2].stride = 1;
  fh.dim[2].n = analyze_read_y;
  fh.dim[2].size = vy;

  fh.dim[3].type = GIFT_Z_DIMENSION;
  fh.dim[3].min = 0;
  if (slice_start >= 0 && slice_start > fh.dim[3].min)
    fh.dim[3].min = slice_start;
  fh.dim[3].max = analyze_read_z - 1;
  if (slice_end >= 0 && slice_end < fh.dim[3].max)
    fh.dim[3].max = slice_end;
  fh.dim[3].stride = 1;
  if (slice_stride > 0)
    fh.dim[3].stride = slice_stride;
  fh.dim[3].n = (fh.dim[3].max - fh.dim[3].min) / fh.dim[3].stride + 1;
  fh.dim[3].size = vz;

  fh.dim[4].type = GIFT_T_DIMENSION;
  fh.dim[4].min = min_time;
  fh.dim[4].max = max_time;
  fh.dim[4].stride = ts;
  fh.dim[4].n = (fh.dim[4].max - fh.dim[4].min) / fh.dim[4].stride + 1;
  fh.dim[4].size = vt;

  if (fh.dim[3].max < fh.dim[3].min)
    Abort("No slices found in selected range.\n");
  fh.n_images = fh.dim[3].n * fh.dim[4].n;
  fh.n_items_per_image = fh.dim[0].n * fh.dim[1].n * fh.dim[2].n;

  fh.corrupt = (char *) malloc(fh.n_images);
  memset(fh.corrupt, 0, fh.n_images);

  input = NULL;
  analyze_read_time = -999999999;
  analyze_read_times_left = 0;
  analyze_read_image_files = analyze_read_img_files_head;

  if (analyze_read_dt == ANALYZE_DATATYPE_UINT8)
    analyze_read_tmp_buf = (unsigned char *) malloc(fh.n_items_per_image);
  image = (char *) malloc(fh.n_items_per_image * GiftBytesPerItem(fh.data_type));
}

void
AnalyzeReadImage (int time,
		  int slice)
{
  int i;
  long offset;
  
  if (time != analyze_read_time)
    {
      if (analyze_read_times_left <= 0)
	{
	  if (input != NULL)
	    {
	      fclose(input);
	      if (analyze_read_images_compressed)
		Delete(analyze_read_img_file_name);
	    }
	  if (analyze_read_images_compressed)
	    analyze_read_file_size = UncompressBatch(analyze_read_img_file_name, &analyze_read_image_files);
	  else
	    {
	      strcpy(analyze_read_img_file_name, analyze_read_image_files->name);
	      analyze_read_file_size = analyze_read_image_files->size;
	      analyze_read_image_files = analyze_read_image_files->next;
	    }
	  analyze_read_times_left = analyze_read_file_size / analyze_read_image_size / analyze_read_z;
	  if ((input = fopen(analyze_read_img_file_name, "r")) == NULL)
	    Abort("Can't open input file %s\n", analyze_read_img_file_name);
	}
      analyze_read_time = time;
      --analyze_read_times_left;
    }

  /* calculate the image offset from the beginning of the file and
     seek to the beginning of the image */
  offset = (analyze_read_numbered ? 0 : (time - 1) * analyze_read_z * analyze_read_image_size) +
    slice * analyze_read_image_size;
  if (fseek(input, offset, SEEK_SET) != 0)
    Abort("Cannot seek to image (time %d, slice %d) in input file.\n",
	  time, slice);
    
  if (analyze_read_dt == ANALYZE_DATATYPE_UINT8)
    {
      if (fread(analyze_read_tmp_buf, 1, fh.n_items_per_image, input) != fh.n_items_per_image)
	Abort("Can't read input file %s\n", analyze_read_img_file_name);
      for (i = 0; i < fh.n_items_per_image; ++i)
	((short *) image)[i] = analyze_read_tmp_buf[i];
    }
  else
    ReadImage(time, slice);
}

void
AnalyzeEndReading ()
{
  FileListElement *fle, *nfle;
  char com[sizeof(Filename) + 64];

  if (input != NULL)
    {
      fclose(input);
      if (analyze_read_images_compressed)
	Delete(analyze_read_img_file_name);
    }
    
  for (fle = analyze_read_hdr_files_head; fle != NULL; fle = nfle)
    {
      nfle = fle->next;
      free(fle);
    }
  for (fle = analyze_read_img_files_head; fle != NULL; fle = nfle)
    {
      nfle = fle->next;
      free(fle);
    }
       
  if (analyze_read_dt == ANALYZE_DATATYPE_UINT8)
    free(analyze_read_tmp_buf);
  
  if (delete_input)
    {
      sprintf(com, "rm -f %s*", analyze_read_basename);
      system(com);
    }
}

