/************************************************************
 *                                                          *
 *  analyze_reader.c                                       *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1995 Department of Statistics,         *
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
 *  Original programming by Mark Fitzgerald  5-96           *
 *  Modified to exclusively use libmri calls for output,    *
 *       Greg Hood (PSC), 9-98                              *
 *  Modified to read header files, and to use LX2           *
 *       resampling stuff, Joel Welling (PSC/Stats), 5-1999 *
 ************************************************************/
/*
 * This is based in large part on the ANALYZE header structure
 * from the Medical Image Format FAQ, which in turn came from
 * Ellis Workman elw@mayo.edu.
 */
/* Notes-
 * -I'm not handling any history info, including what might sometimes
 *  specify a different filename for the data.
 * -I'm not handling color tables.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#if (SGI64 || SGI5 || SGIMP)
#include <bstring.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <math.h>

#include "mri.h"
#include "bio.h"
#include "fmri.h"
#include "array.h"
#include "stdcrg.h"
#include "misc.h"
#include "smartreader.h"
#include "nr_sub.h"
#include "rcn.h"

static char rcsid[] = "$Id: analyze_reader.c,v 1.6 2004/10/26 23:54:08 welling Exp $";

#define FRZ_ANALYZE_HDR_KEY_SZ 40
#define FRZ_ANALYZE_IMAGE_DIM_SZ 108
#define FRZ_ANALYZE_DATA_HST_SZ 200
#define FRZ_ANALYZE_TOT_SZ (FRZ_ANALYZE_HDR_KEY_SZ + FRZ_ANALYZE_IMAGE_DIM_SZ + FRZ_ANALYZE_DATA_HST_SZ)

#define FRZ_ANALYZE_SIZEOF_HDR_OFF 0
#define FRZ_ANALYZE_DATA_TYPE_OFF 4
#define FRZ_ANALYZE_DATA_TYPE_LEN 10
#define FRZ_ANALYZE_DB_NAME_OFF 14
#define FRZ_ANALYZE_DB_NAME_LN 18
#define FRZ_ANALYZE_EXTENTS_OFF 32
#define FRZ_ANALYZE_SESSION_ERROR_OFF 36
#define FRZ_ANALYZE_REGULAR_OFF 38
#define FRZ_ANALYZE_HKEY_UN0_OFF 39

#define FRZ_ANALYZE_DIM_OFF 40
#define FRZ_ANALYZE_DIM_LN 8
#define FRZ_ANALYZE_VOX_UNITS_OFF 56
#define FRZ_ANALYZE_VOX_UNITS_LN 4
#define FRZ_ANALYZE_CAL_UNITS_OFF 60
#define FRZ_ANALYZE_CAL_UNITS_LN 8
#define FRZ_ANALYZE_DATATYPE_OFF 70
#define FRZ_ANALYZE_PIXDIM_OFF 76
#define FRZ_ANALYZE_PIXDIM_LN 8
#define FRZ_ANALYZE_VOX_OFFSET_OFF 108
#define FRZ_ANALYZE_COMPRESSED_OFF 132

#define FRZ_ANALYZE_DESCRIP_OFF 148
#define FRZ_ANALYZE_DESCRIP_LN 80
#define FRZ_ANALYZE_AUXFILE_OFF 228
#define FRZ_ANALYZE_AUXFILE_LN 24

#define FRZ_NIFTI_MAGIC_OFF 344
#define FRZ_NIFTI_MAGIC_LN 4

static int
IntBRdFloat32 (unsigned char *addr)
{
  return((int) Round(BRdFloat32(addr)));
}

static int scan_ANALYZE_header(KVHash* info, const char* readfile)
{
  KVHash* defs;
  KVHash* extNames;
  FILE *fphead;
  unsigned char hdr[FRZ_ANALYZE_TOT_SZ];
  char buf[512];      /* scratch space */
  short dim[FRZ_ANALYZE_DIM_LN];
  float pixdim[FRZ_ANALYZE_PIXDIM_LN];
  float vox_offset;
  int ierror= 0;
  int i;
  int dtype;
  
  /* This bit is from the GE code io_signa_lx.c, with mods for portability */
  ierror= 0;
  if ((fphead = fopen(readfile,"r"))!=NULL)
    {
      if (fseek(fphead, (long) 0, SEEK_SET)) {
	perror("Error seeking header");
	ierror=1;
      }
      else {
	if (fread(hdr,sizeof(char),FRZ_ANALYZE_TOT_SZ, fphead)
	    != FRZ_ANALYZE_TOT_SZ) {
	  perror("Error reading header");
	  ierror=1;
	}
	else {
	  if (fclose(fphead)) {
	    perror("Error closing header");
	    ierror=1;
	  }
	}
      }
    }
  else {
    perror("Error opening header");
    ierror= 1;
  }
  if (ierror) return 0;
  
  /* Definitions */
  defs= kvGetHash(info,"definitions");
  extNames= kvGetHash(info,"external_names");
  
  /* We test to see if the data does have the expected endian order.
   */
  if (BRdInt32(hdr+FRZ_ANALYZE_SIZEOF_HDR_OFF)!=FRZ_ANALYZE_TOT_SZ) {
    /* Oops, try it the other way! */
    bio_big_endian_input= (bio_big_endian_input ? 0 : 1);
    if (BRdInt32(hdr+FRZ_ANALYZE_SIZEOF_HDR_OFF)!=FRZ_ANALYZE_TOT_SZ) 
      Abort("%s: internal error: analyze_reader unexpectedly can't find file length!\n",
	    progname);
  }
  
  kvDefBoolean(info,"big_endian_input",bio_big_endian_input);
  if (debug) fprintf(stderr,"Header indicates %s input\n",
		     bio_big_endian_input ? "bigendian" : "littleendian");
  
  /* If the image is compressed, we're out of luck */
  if (BRdInt32(hdr+FRZ_ANALYZE_COMPRESSED_OFF) != 0)
    Warning(1,"%s: header indicates compressed image data!\n",progname);

#ifdef never
  /* The doc says this should be 16384, but it doesn't seem to be
   * reliably set.
   */
  fprintf(stderr,"extents= %d\n",BRdInt32(hdr+FRZ_ANALYZE_EXTENTS_OFF));
#endif

  strncpy(buf,(char*)hdr+FRZ_ANALYZE_VOX_UNITS_OFF,FRZ_ANALYZE_VOX_UNITS_LN);
  if (strlen(buf)>0) kvDefString(info,"ANALYZE_vox_units",buf);
  strncpy(buf,(char*)hdr+FRZ_ANALYZE_CAL_UNITS_OFF,FRZ_ANALYZE_CAL_UNITS_LN);
  if (strlen(buf)>0) kvDefString(info,"ANALYZE_cal_units",buf);
  strncpy(buf,(char*)hdr+FRZ_ANALYZE_DESCRIP_OFF,FRZ_ANALYZE_DESCRIP_LN);
  if (strlen(buf)>0) kvDefString(info,"ANALYZE_description",buf);
  strncpy(buf,(char*)hdr+FRZ_ANALYZE_DB_NAME_OFF,FRZ_ANALYZE_DB_NAME_LN);
  if (strlen(buf)>0) kvDefString(info,"ANALYZE_db_name",buf);
  strncpy(buf,(char*)hdr+FRZ_NIFTI_MAGIC_OFF,FRZ_NIFTI_MAGIC_LN);
  if (strlen(buf)>0) kvDefString(info,"ANALYZE_nifti_magic",buf);
  
  

  dtype= BRdInt16(hdr+FRZ_ANALYZE_DATATYPE_OFF);
  switch (dtype) {
  case 1: /* binary */
    Abort("%s: analyze_reader: binary datatype not supported!\n",progname);
    break;
  case 2: /* uchar */
    kvDefInt(info,"datatype_in",SRDR_UINT8);
    break;
  case 4: /* short */
    kvDefInt(info,"datatype_in",SRDR_INT16);
    break;
  case 8: /* int */
    kvDefInt(info,"datatype_in",SRDR_INT32);
    break;
  case 16: /* float */
    kvDefInt(info,"datatype_in",SRDR_FLOAT32);
    break;
  case 32: /* complex */
    kvDefInt(info,"datatype_in",SRDR_FLOAT32);
    kvDefInt(info,"dv",2);
    break;
  case 64: /* double */
    kvDefInt(info,"datatype_in",SRDR_FLOAT64);
    break;
  case 128: /* RGB */
    kvDefInt(info,"datatype_in",SRDR_UINT8);
    kvDefInt(info,"dv",3);
    break;
  default: /* invalid or none */
    Abort("%s: analyze_reader: unknown datatype not supported!\n",progname);
  }

  kvDefInt(info,"handler_datatype_out",kvGetInt(info,"datatype_in"));

  for (i=0; i<FRZ_ANALYZE_DIM_LN; i++) {
    dim[i]= BRdInt16(hdr+FRZ_ANALYZE_DIM_OFF + 2*i);
  }
  if (dim[0]>7) dim[0]= 7; /* worry about invalid file */

  for (i=0; i<FRZ_ANALYZE_PIXDIM_LN; i++) {
    pixdim[i]= BRdFloat32(hdr+FRZ_ANALYZE_PIXDIM_OFF + 4*i);
  }

  if (kvLookup(info,"dv") != NULL) {
    strncpy(buf,"vxyztabc",dim[0]+1);
    buf[dim[0]+1]= '\0';
  }
  else {
    strncpy(buf,"xyztabc",dim[0]);
    buf[dim[0]]= '\0';
  }
  kvDefString(info,"dimstr",buf);

  if (dim[0]>=1) {
    kvDefInt(info,"dx",dim[1]);
    kvDefString(info,"description.x","gridded image-space");
  }
  if (dim[0]>=2) {
    kvDefInt(info,"dy",dim[2]);
    kvDefString(info,"description.y","gridded image-space");
  }
  if (dim[0]>=3) {
    kvDefInt(info,"dz",dim[3]);
    kvDefString(info,"description.z","gridded image-space");
  }
  if (dim[0]>=4) {
    kvDefInt(info,"dt",dim[4]);
    kvDefString(info,"description.t","gridded image-space");
  }
  if (dim[0]>=5) kvDefInt(info,"da",dim[5]);
  if (dim[0]>=6) kvDefInt(info,"db",dim[6]);
  if (dim[0]>=7) kvDefInt(info,"dc",dim[7]);

  if (pixdim[1]!=0.0) {
    kvDefDouble(info,"voxel_x",pixdim[1]);
    kvDefString(defs,"voxel_x","X voxel size (mm)");
  }
  if (pixdim[2]!=0.0) {
    kvDefDouble(info,"voxel_y",pixdim[2]);
    kvDefString(defs,"voxel_y","Y voxel size (mm)");
  }
  if (pixdim[3]!=0.0) {
    kvDefDouble(info,"voxel_z",pixdim[3]);
    kvDefString(defs,"voxel_z","Z voxel size including gap (mm)");
  }
  if (pixdim[4]!=0.0) {
    kvDefDouble(info,"voxel_t",pixdim[2]);
    kvDefString(defs,"voxel_t","t voxel size");
  }

  /* I have no idea why vox_offset is a float. */
  vox_offset= BRdFloat32(hdr+FRZ_ANALYZE_VOX_OFFSET_OFF);
  if (debug) fprintf(stderr,"vox_offset: %f\n",vox_offset);
  if (vox_offset>0.0) {
    long long offset= (long long)rint(vox_offset);
    kvDefLong(info,"start_offset",offset);
  }
  else if (vox_offset<0.0) {
    long long skip= (long long)rint( -vox_offset );
    kvDefLong(info,"sliceskip",skip);
  }
  else {
    /* We set no offsets if vox_offset is strictly 0.0 */
  }

  return 1;
}

static void processHeader( FileHandler* self, KVHash* info, SList* cStack )
{
  /* Call the base class method */
  baseProcessHeader( self, info, cStack );

  if (!scan_ANALYZE_header(info, self->fileName))
    Abort("analyze_reader: unable to read or parse header from <%s>!\n",
	  self->fileName);
}

static void analyzeClose( FileHandler* self )
{
  /* The file pointer is pointed to the .img file, but the process
   * to close it is the same as for the base case.
   */
  baseClose(self);
}

static void analyzeReopen( FileHandler* self )
{
  /* We actually want to reopen the img file rather than the hdr. */
  if (!(self->file)) {
    char buf[256];
    char* here;
    strncpy(buf, self->fileName, sizeof(buf));
    buf[sizeof(buf)-1]= '\0';
    if (!(here= strrchr(buf,'.')))
      Abort("%s: Unexpectedly can't find the extension in <%s>!\n",
	    progname, buf);
    if (strcmp(here,".hdr"))
      Abort("%s: Extension of an ANALYZE header is not .hdr!\n",progname);
    strcpy(here,".img"); /* just overwrite the letters */
    if (!(self->file= fopen(buf,"r")))
      Abort("%s: unable to open file <%s> for reading!\n",
	    progname,self->fileName);
  }
}

static void analyzeRead( FileHandler* self, KVHash* info,
			 long long offset, long n,
			 SRDR_Datatype datatype_out, void* obuf )
{

  /* We just need to point ourselves at the .BRIK file */
  analyzeReopen(self);

  baseRead( self, info, offset, n, datatype_out, obuf );
}

FileHandler* analyzeFactory(char* fname, KVHash* info)
{
  FileHandler* result= baseFactory(fname);
  result->processHeader= processHeader;
  result->read= analyzeRead;
  result->close= analyzeClose;
  result->reopen= analyzeReopen;
  result->typeName= strdup( "ANALYZE Image" );

  return result;
}

static long long getFileSize( const char* fname )
{
  struct stat s;
  if (strcmp(fname,"NotARealFile")) {
    if (stat(fname,&s))
      Abort("%s: base_reader: stat failed on %s: %s!\n",
	    progname, fname, strerror(errno));
    return (long long)s.st_size;
  }
  else return 0;
}

static int analyzeHeaderTest( const char* filename )
{
  FILE* f;
  char buf[4];
  int ierror= 0;

  if (getFileSize(filename) != FRZ_ANALYZE_TOT_SZ) return 0;
  else {
    if ((f = fopen(filename,"r"))!=NULL)
      {
	if (fseek(f, (long) FRZ_ANALYZE_SIZEOF_HDR_OFF, SEEK_SET)) {
	  perror("Error seeking header");
	  ierror=1;
	}
	else {
	  if (fread(buf,sizeof(char),4, f) != 4) {
	    perror("Error reading header");
	    ierror=1;
	  }
	  else {
	    if (fclose(f)) {
	      perror("Error closing header");
	      ierror=1;
	    }
	  }
	}
      }
    else {
      perror("Error opening header");
      ierror= 1;
    }
    
    if (ierror) return 0;
    else {
      if (BRdInt32((unsigned char*)buf)==FRZ_ANALYZE_TOT_SZ) return 1;
      else {
	/* Try the other byte order */
	bio_big_endian_input= (bio_big_endian_input ? 0 : 1);
	if (BRdInt32((unsigned char*)buf)==FRZ_ANALYZE_TOT_SZ) return 1;
	else return 0;
      }
    }
  }
}

int analyzeTester(const char* filename)
{
  return analyzeHeaderTest(filename);
}

