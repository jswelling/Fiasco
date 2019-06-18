/************************************************************
 *                                                          *
 *  png_reader.c                                         *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 2007 Department of Statistics,         *
 *                        Carnegie Mellon University        *
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
 *                                                          *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "mri.h"
#if (SGI64 || SGI5 || SGIMP)
#include <bstring.h>
#endif

#include "bio.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"
#include "smartreader.h"

static char rcsid[] = "$Id: png_reader.c,v 1.5 2007/06/25 22:25:58 welling Exp $";

/* Notes-
 * -get time and date from text if necessary?
 * -time and text can come after image; we'll miss them if they do
 */

#ifdef USE_PNG

#include <png.h>

/* A reader for PNG files */
typedef struct png_data_struct {
  png_structp png_ptr;
  png_infop info_ptr;
  png_infop end_info;
  png_time modtime;
  long long currentOffset;
  unsigned char* buf;
} PngData;

static void pngDestroySelf( FileHandler* self )
{
  PngData* data= (PngData*)(self->hook);
  if (data->buf) free(data->buf);
  png_destroy_read_struct(&(data->png_ptr), &(data->info_ptr), 
			  &(data->end_info));
  baseDestroySelf(self);
}

static void custom_read_data(png_structp png_ptr,
			     png_bytep bufptr, png_size_t length)
{
  FileHandler* self= (FileHandler*)(png_get_io_ptr(png_ptr));
  PngData* data= (PngData*)(self->hook);
  if (!(self->file)) {
    FH_REOPEN(self);
    if (fseek(self->file, data->currentOffset, SEEK_SET))
      perror("Reopened file and cannot seek!");
  }
  data->currentOffset += fread(bufptr, 1, length, self->file);
}

static void pngProcessHeader( FileHandler* self, KVHash* info, 
			      SList* chunkStack )
{
  PngData* data= (PngData*)(self->hook);
  KVHash* defs= kvGetHash(info,"definitions");
  png_uint_32 width;
  png_uint_32 height;
  int bit_depth;
  int color_type;
  int interlace_type;
  int compression_type;
  int filter_method;
  int nChannels;
  png_uint_32 resX;
  png_uint_32 resY;
  png_timep timep;
  png_text* comments= NULL;
  int nComments;
  char str[64];
  char str2[256];
  int i;

  FH_REOPEN(self);
  data->currentOffset= 0;

  if (setjmp(png_jmpbuf(data->png_ptr))) {
        png_destroy_read_struct(&(data->png_ptr), &(data->info_ptr),
           &(data->end_info));
	FH_CLOSE(self);
        Abort("%s: error within libpng on file %s!\n",progname,self->fileName);
    }

  png_set_read_fn( data->png_ptr, self, custom_read_data );
  png_read_info(data->png_ptr, data->info_ptr);
  png_get_IHDR(data->png_ptr, data->info_ptr, &width, &height,
	       &bit_depth, &color_type, &interlace_type, &compression_type,
	       &filter_method);
  nChannels= png_get_channels(data->png_ptr, data->info_ptr);
  if (png_get_tIME(data->png_ptr, data->info_ptr, &(timep))
      == PNG_INFO_tIME) {
    bcopy(timep,&(data->modtime),sizeof(png_time)); 
    snprintf(str,sizeof(str),"%02d/%02d/%04d",
	     (int)data->modtime.month,
	     (int)data->modtime.day,
	     (int)data->modtime.year);
    kvDefString(info,"date",str);
    kvDefString(defs,"date","date modified");
    snprintf(str,sizeof(str),"%02d:%02d:%02d",
	     (int)data->modtime.hour,
	     (int)data->modtime.minute,
	     (int)data->modtime.second);
    kvDefString(info,"time",str);
    kvDefString(defs,"time","time modified");
  }
  else {
    if (debug) fprintf(stderr,"time info in %s is not valid\n",
		       self->fileName);
    bzero(&(data->modtime),sizeof(png_time));
  }
  if (nChannels==1) kvDefString(info,"dimstr","xy");
  else {
    kvDefString(info,"dimstr","vxy");
    kvDefInt(info,"dv",nChannels);
  }
  kvDefInt(info,"dx",(long)width);
  kvDefInt(info,"dy",(long)height);
  kvDefInt(info,"png_color_type",color_type);
  kvDefInt(info,"png_bit_depth",bit_depth);
  kvDefInt(info,"png_interlace_type",interlace_type);
  kvDefInt(info,"png_compression_type",compression_type);
  kvDefInt(info,"png_filter_method",filter_method);

  if (color_type == PNG_COLOR_TYPE_PALETTE)
    png_set_palette_to_rgb(data->png_ptr);
  
  /* We can coerce libpng into inflating 1, 2, or 4 bit data to 8,
   * so we only need to worry about the 8 and 16 bit depth cases.
   */
  if (color_type == PNG_COLOR_TYPE_GRAY &&
      bit_depth < 8) png_set_expand_gray_1_2_4_to_8(data->png_ptr);
  if (bit_depth==16) {
    kvDefInt(info,"datatype_in",SRDR_UINT16);
    Warning(1,"%s: PNG image pixel data is unsigned shorts; this format is not completely supported!\n",progname);
  }
  else {
    kvDefInt(info,"datatype_in",SRDR_UINT8);
  }
  kvDefInt(info,"handler_datatype_out",kvGetInt(info,"datatype_in"));
 
  nComments= png_get_text(data->png_ptr, data->info_ptr, &comments, NULL);
  for (i=0; i<nComments; i++) {
    char* here;
    snprintf(str,sizeof(str),"png_txt_%s",comments[i].key);
    strncpy(str2,comments[i].text,sizeof(str2));
    here= str2+strlen(str2)-1;
    while (isspace(*here) && here>=str2) *here--= '\0';
    kvDefString(info,str,str2);
  }

  resX= png_get_x_pixels_per_meter(data->png_ptr,data->info_ptr);
  resY= png_get_y_pixels_per_meter(data->png_ptr,data->info_ptr);
  if (resX != 0) {
    kvDefDouble(info,"voxel_x",1000.0/resX);
    kvDefString(defs,"voxel_x","X voxel size including gap (mm)");
  }
  if (resY != 0) {
    kvDefDouble(info,"voxel_y",1000.0/resY);
    kvDefString(defs,"voxel_y","Y voxel size including gap (mm)");
  }

  png_read_update_info(data->png_ptr, data->info_ptr);

  FH_CLOSE(self);
}

static void pngRead( FileHandler* self, KVHash* info,
		     long long offset, long n,
		     SRDR_Datatype datatype, void* obuf )
{
  PngData* data= (PngData*)(self->hook);
  SRDR_Datatype srdrType= kvGetInt(info,"datatype_in");
  const char* dimstr= kvGetString(info,"dimstr");

  if (datatype != srdrType)
    Abort("%s: pngRead: a datatype translation is needed! (%s to %s)\n",
	  progname,srdrTypeName[srdrType],srdrTypeName[datatype]);

  if (!(data->buf)) {
    /* Image bytes have not yet been read in */
    int rowbytes= png_get_rowbytes(data->png_ptr, data->info_ptr);
    int bytesPerPixel= rowbytes/kvGetInt(info,"dx");
    int nrows= kvGetInt(info,"dy");
    int i;
    unsigned char** rowPtrs= NULL;

    if (setjmp(png_jmpbuf(data->png_ptr))) {
      png_destroy_read_struct(&(data->png_ptr), &(data->info_ptr),
			      &(data->end_info));
      FH_CLOSE(self);
      Abort("%s: error within libpng on file %s!\n",progname,self->fileName);
    }

    if (!(data->buf= (unsigned char*)malloc(rowbytes*nrows)))
      Abort("%s: pngRead: unable to allocate %d bytes!\n",
	    rowbytes*nrows);
    if (!(rowPtrs=(unsigned char**)malloc(nrows*sizeof(unsigned char*))))
      Abort("%s: pngRead: unable to allocate %d bytes!\n",
	    nrows*sizeof(unsigned char*));
    for (i=0; i<nrows; i++) rowPtrs[i]= &(data->buf[i*rowbytes]);

    png_read_image(data->png_ptr, rowPtrs);
    png_read_end(data->png_ptr, data->end_info);
    if (png_get_tIME(data->png_ptr, data->info_ptr, NULL)
	== PNG_INFO_tIME) {
      fprintf(stderr,"*now* there is time info!\n");
    }
    else fprintf(stderr,"time info still invalid\n");
    free(rowPtrs);
    FH_CLOSE(self);
  }

  switch (srdrType) {
  case SRDR_UINT8:
    {
      bcopy(data->buf+offset, obuf, n);
    }
    break;
  case SRDR_UINT16:
    {
      bcopy(data->buf+offset, obuf, 2*n);
    }
    break;
  default:
    Abort("%s: png_reader: internal error; unexpected reader datatype!\n",
	  progname);
  }
}

static int pngCompare( FileHandler* f1, FileHandler* f2 )
{
  PngData* data1= (PngData*)(f1->hook);
  PngData* data2= (PngData*)(f2->hook);

  if (data1->modtime.year<data2->modtime.year) return -1;
  else if (data1->modtime.year>data2->modtime.year) return 1;
  else if (data1->modtime.month<data2->modtime.month) return -1;
  else if (data1->modtime.month>data2->modtime.month) return 1;
  else if (data1->modtime.day<data2->modtime.day) return -1;
  else if (data1->modtime.day>data2->modtime.day) return 1;
  else if (data1->modtime.hour<data2->modtime.hour) return -1;
  else if (data1->modtime.hour>data2->modtime.hour) return 1;
  else if (data1->modtime.minute<data2->modtime.minute) return -1;
  else if (data1->modtime.minute>data2->modtime.minute) return 1;
  else if (data1->modtime.second<data2->modtime.second) return -1;
  else if (data1->modtime.second>data2->modtime.second) return 1;
  else return (baseCompare(f1,f2));
}

FileHandler* pngFactory(char* fname, KVHash* info)
{
  FileHandler* result= baseFactory(fname);
  PngData* data;

  result->typeName= strdup("PngDataHandler");

  result->destroySelf= pngDestroySelf;
  result->read= pngRead;
  result->processHeader= pngProcessHeader;
  result->compareThisType= pngCompare;

  if (!(data= (PngData*)malloc(sizeof(PngData))))
    Abort("%s: unable to allocate %d bytes!\n",progname,sizeof(PngData));
  result->hook= data;
  data->png_ptr= png_create_read_struct( PNG_LIBPNG_VER_STRING,
					 NULL, NULL, NULL );
  if (!(data->png_ptr)) {
    Abort("Could not create png_struct!\n");
  }
  data->info_ptr= png_create_info_struct(data->png_ptr);
  if (!(data->info_ptr)) {
    png_destroy_read_struct(&(data->png_ptr), NULL, NULL);
    Abort("Could not create png info_ptr!\n");
  }
  data->end_info= png_create_info_struct(data->png_ptr);
  if (!(data->end_info)) {
    png_destroy_read_struct(&(data->png_ptr),&(data->info_ptr),NULL);
    Abort("Could not create png end_info ptr!\n");
  }
  data->buf= NULL;

  if (setjmp(png_jmpbuf(data->png_ptr))) {
        png_destroy_read_struct(&(data->png_ptr), &(data->info_ptr),
           &(data->end_info));
	/* No need to close the file, since it isn't currently open */
        return (NULL);
    }

  return result;
}

int pngTester(const char* filename)
{
  FILE* f= NULL;
  unsigned char hdr[8];
  int match= 0;

  if ((f = fopen(filename,"r"))==NULL) return 0;
  if (fread(hdr,1,sizeof(hdr),f) != sizeof(hdr)) {
    perror("Error reading file");
    fclose(f);
    return 0;
  }
  match= !png_sig_cmp(hdr,0,sizeof(hdr)); 
  if (fclose(f)) perror("Error closing file");
  return match;
}

#else /* ifdef USE_PNG */

FileHandler* pngFactory(char* fname, KVHash* info)
{
  FileHandler* result= baseFactory("NotARealFile");
  return result;
}

int pngTester(const char* filename)
{
  return 0;
}

#endif

