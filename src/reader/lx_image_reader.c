/************************************************************
 *                                                          *
 *  lx_image_reader.c                                       *
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
 * Much of the header-parsing algorithm is taken from Ifile.c 
 * in the AFNI package.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#if (SGI64 || SGI5 || SGIMP)
#include <bstring.h>
#endif

#include "mri.h"
#include "bio.h"
#include "fmri.h"
#include "array.h"
#include "stdcrg.h"
#include "misc.h"
#include "smartreader.h"
#include "nr_sub.h"
#include "rcn.h"

static char rcsid[] = "$Id: lx_image_reader.c,v 1.15 2004/09/10 00:29:48 welling Exp $";

#define FRZ_LXI_VALID_LOGO "IMGF"
#define FRZ_LXI_VALID_LOGO_BACKWARDS "FGMI"
#define FRZ_LXI_VALID_LOGO_OFF 0
#define FRZ_LXI_HEADER_SIZE_BYTES 7904
#define FRZ_LXI_SKIP_OFF 4
#define FRZ_LXI_XDIM_OFF 8
#define FRZ_LXI_YDIM_OFF 12
#define FRZ_LXI_BPP_OFF 16
#define FRZ_LXI_COMPRESSED_OFF 20
#define FRZ_LXI_GEMS_HDR_OFF 148
#define FRZ_LXI_GEMS_VOXEL_X_OFF 50
#define FRZ_LXI_GEMS_VOXEL_Y_OFF 54
#define FRZ_LXI_GEMS_VOXEL_Z_OFF 26
#define FRZ_LXI_GEMS_TR_OFF 194
#define FRZ_LXI_GEMS_TE_OFF 202
#define FRZ_LXI_GEMS_DATE_OFF 1402
#define FRZ_LXI_GEMS_TIME_OFF 1412
#define FRZ_LXI_GEMS_TLC_OFF 154
#define FRZ_LXI_GEMS_TRC_OFF 166
#define FRZ_LXI_GEMS_BRC_OFF 178

static int
IntBRdFloat32 (unsigned char *addr)
{
  return((int) Round(BRdFloat32(addr)));
}

static int scan_LXI_header(KVHash* info, const char* readfile)
{
  KVHash* defs;
  KVHash* extNames;
  FILE *fphead;
  unsigned char hdr[FRZ_LXI_HEADER_SIZE_BYTES];
  char buf[256];      /* scratch space */
  int ierror= 0;
  int xdim;
  int gems_off;
  unsigned char* gemshdr;
  double voxel[3];
  double tlc[3];
  double trc[3];
  double brc[3];
  double edge_lr[3];
  double edge_bt[3];
  double norm[3];
  float temp[3];
  int i;
  int axis_sum;
  
  /* This bit is from the GE code io_signa_lx.c, with mods for portability */
  ierror= 0;
  if ((fphead = fopen(readfile,"r"))!=NULL)
    {
      if (fseek(fphead, (long) 0, SEEK_SET)) {
	perror("Error seeking header");
	ierror=1;
      }
      else {
	if (fread(hdr,sizeof(char),FRZ_LXI_HEADER_SIZE_BYTES, fphead)
	    != FRZ_LXI_HEADER_SIZE_BYTES) {
	  if (ferror(fphead)) {
	    perror("Error reading header");
	  }
	  else Abort("%s: file <%s> is shorter than expected!\n",
		  progname, readfile);
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

  /* We test to see if the data does have the expected endian order-
   * x dimension should be some reasonable.
   */
  xdim= BRdInt32(hdr+FRZ_LXI_XDIM_OFF);
  if (xdim<0 || xdim>8192) {
    /* Oops, try it the other way! */
    bio_big_endian_input= (bio_big_endian_input ? 0 : 1);
    xdim= BRdInt32(hdr+FRZ_LXI_XDIM_OFF);
  }
  kvDefBoolean(info,"big_endian_input",bio_big_endian_input);
  if (debug) fprintf(stderr,"Header indicates %s input\n",
		     bio_big_endian_input ? "bigendian" : "littleendian");

  /* If the image is compressed, we're out of luck */
  if (BRdInt32(hdr+FRZ_LXI_COMPRESSED_OFF) != 1)
    Warning(1,"%s: header indicates compressed image data!\n",progname);

  kvDefLong(info,"start_offset",BRdInt32(hdr+FRZ_LXI_SKIP_OFF));
  kvDefString(info,"dimstr","xyz");
  /* image space, so no row flipping! */
  kvDefBoolean(info,"rowflip",0); 
  kvDefString(info,"rowflip_pattern","none");
  kvDefString(defs,"rowflip","EPI row reversal needed");
  kvDefString(defs,"rowflip_pattern","EPI row reversal pattern");
  kvDefInt(info,"dx",xdim);
  kvDefString(info,"description.x","gridded image-space");
  kvDefInt(info,"dy",BRdInt32(hdr+FRZ_LXI_YDIM_OFF));
  kvDefString(info,"description.y","gridded image-space");
  kvDefInt(info,"dz",1);
  kvDefString(info,"description.z","gridded image-space");
  /* Type depends on bits per pixel */
  switch (BRdInt32(hdr+FRZ_LXI_BPP_OFF)) {
  case 8: kvDefInt(info,"datatype_in",SRDR_UINT8); break;
  case 16: kvDefInt(info,"datatype_in",SRDR_INT16); break;
  case 32: kvDefInt(info,"datatype_in",SRDR_INT32); break;
  default: Abort("%s: lx_image_reader can't handle %d bits per pixel!\n",
		 progname,BRdInt32(hdr+FRZ_LXI_BPP_OFF));
  }

  if (kvGetLong(info,"start_offset") < FRZ_LXI_HEADER_SIZE_BYTES) {
    Warning(1,"This appears to be only a partial header.\n");
  }
  else {
   
    gems_off= BRdInt32(hdr+FRZ_LXI_GEMS_HDR_OFF);
    if (debug) fprintf(stderr,"gems_off= %d\n",gems_off);
    if (gems_off<=0 || gems_off+256>=FRZ_LXI_HEADER_SIZE_BYTES)
      Abort("%s: lx_image_reader: internal header offset of %d looks invalid!\n",
	    progname,gems_off);
    gemshdr= hdr+gems_off;
    
    kvDefInt(info,"TR",BRdInt32(gemshdr + FRZ_LXI_GEMS_TR_OFF));
    kvDefString(defs,"TR","TR (us)");
    kvDefInt(info,"TE",BRdInt32(gemshdr + FRZ_LXI_GEMS_TE_OFF));
    kvDefString(defs,"TE","TE (us)");
    
    strncpy(buf,(char*)gemshdr + FRZ_LXI_GEMS_DATE_OFF,10);
    buf[10]= '\0';
    kvDefString(info,"date",buf);
    kvDefString(defs,"date","scan date");
    strncpy(buf,(char*)gemshdr + FRZ_LXI_GEMS_TIME_OFF,8);
    buf[8]= '\0';
    kvDefString(info,"time",buf);
    kvDefString(defs,"time","scan time");
    
    /* We'll grab the voxel sizes now, but they are subject to
     * shuffling based on the stated corner order.
     */
    voxel[0]= BRdFloat32(gemshdr + FRZ_LXI_GEMS_VOXEL_X_OFF);
    voxel[1]= BRdFloat32(gemshdr + FRZ_LXI_GEMS_VOXEL_Y_OFF);
    voxel[2]= BRdFloat32(gemshdr + FRZ_LXI_GEMS_VOXEL_Z_OFF);
    
    BRdFloat32Array(gemshdr + FRZ_LXI_GEMS_TLC_OFF, temp, 3);
    for (i=0; i<3; i++) tlc[i]= temp[i]; /* type conversion */
    BRdFloat32Array(gemshdr + FRZ_LXI_GEMS_TRC_OFF, temp, 3);
    for (i=0; i<3; i++) trc[i]= temp[i]; /* type conversion */
    BRdFloat32Array(gemshdr + FRZ_LXI_GEMS_BRC_OFF, temp, 3);
    for (i=0; i<3; i++) brc[i]= temp[i]; /* type conversion */

    /* The LX internal coordinate system is oriented differently from
     * Fiasco's.  LX has X running right-left in the images and Z
     * decreasing toward the feet; Fiasco has X left-right and Z
     * increasing toward the feet.
     */
    tlc[0] *= -1.0;
    trc[0] *= -1.0;
    brc[0] *= -1.0;
    tlc[2] *= -1.0;
    trc[2] *= -1.0;
    brc[2] *= -1.0;
    
    subtractVec3( edge_lr, trc, tlc );
    subtractVec3( edge_bt, trc, brc );
    crossVec3( norm, edge_lr, edge_bt );
    normalizeVec3( norm );
    
    if (debug) {
      fprintf(stderr,"corners: tlc= (%f, %f, %f)\n", tlc[0], tlc[1], tlc[2]);
      fprintf(stderr,"         trc= (%f, %f, %f)\n", trc[0], trc[1], trc[2]);
      fprintf(stderr,"         brc= (%f, %f, %f)\n", brc[0], brc[1], brc[2]);
      fprintf(stderr,"         norm=(%f, %f, %f)\n", norm[0], norm[1], norm[2]);
    }
    
    axis_sum= 0;
    /* X direction is the direction of greatest change between TRC and TLC */
    if (fabs(edge_lr[0])>fabs(edge_lr[1])) {
      if (fabs(edge_lr[0])>fabs(edge_lr[2])) {
	/* component 0 longest */
	kvDefDouble(info,"voxel_x",voxel[0]);
	axis_sum += 0;
      }
      else {
	/* component 2 longest */
	kvDefDouble(info,"voxel_x",voxel[2]);
	axis_sum += 2;
      }
    }
    else {
      if (fabs(edge_lr[1])>fabs(edge_lr[2])) {
	/* component 1 longest */
	kvDefDouble(info,"voxel_x",voxel[1]);
	axis_sum += 1;
      }
      else {
	/* component 2 longest */
	kvDefDouble(info,"voxel_x",voxel[2]);
	axis_sum += 2;
      }
    }
    
    /* Y direction is the direction of greatest change between BRC and TRC */
    if (fabs(edge_bt[0])>fabs(edge_bt[1])) {
      if (fabs(edge_bt[0])>fabs(edge_bt[2])) {
	/* component 0 longest */
	kvDefDouble(info,"voxel_y",voxel[0]);
	axis_sum += 0;
      }
      else {
	/* component 2 longest */
	kvDefDouble(info,"voxel_y",voxel[2]);
	axis_sum += 2;
      }
    }
    else {
      if (fabs(edge_bt[1])>fabs(edge_bt[2])) {
	/* component 1 longest */
	kvDefDouble(info,"voxel_y",voxel[1]);
	axis_sum += 1;
      }
      else {
	/* component 2 longest */
	kvDefDouble(info,"voxel_y",voxel[2]);
	axis_sum += 2;
      }
    }
    
    /* And the left-over direction is Z. */
    if (axis_sum > 3)
      Abort("%s: disordered axes in lx_image_reader!\n",progname);
    kvDefDouble(info,"voxel_z",voxel[3-axis_sum]);
    
    kvDefString(defs,"voxel_x","X voxel size including gap (mm)");
    kvDefString(defs,"voxel_y","Y voxel size including gap (mm)");
    kvDefString(defs,"voxel_z","Z voxel size including gap (mm)");
    
    defVec3( info, "slice_tlc", tlc );
    defVec3( info, "slice_trc", trc );
    defVec3( info, "slice_brc", brc );
    defVec3( info, "slice_norm", norm );
    
    kvDefString(defs,"slice_tlc.0", "slice top left corner");
    kvDefString(extNames,"slice_tlc.0", "tlc.0");
    kvDefString(defs,"slice_tlc.1", "slice top left corner");
    kvDefString(extNames,"slice_tlc.1", "tlc.1");
    kvDefString(defs,"slice_tlc.2", "slice top left corner");
    kvDefString(extNames,"slice_tlc.2", "tlc.2");
    
    kvDefString(defs,"slice_trc.0", "slice top right corner");
    kvDefString(extNames,"slice_trc.0", "trc.0");
    kvDefString(defs,"slice_trc.1", "slice top right corner");
    kvDefString(extNames,"slice_trc.1", "trc.1");
    kvDefString(defs,"slice_trc.2", "slice top right corner");
    kvDefString(extNames,"slice_trc.2", "trc.2");
    
    kvDefString(defs,"slice_brc.0", "slice bottom right corner");
    kvDefString(extNames,"slice_brc.0", "brc.0");
    kvDefString(defs,"slice_brc.1", "slice bottom right corner");
    kvDefString(extNames,"slice_brc.1", "brc.1");
    kvDefString(defs,"slice_brc.2", "slice bottom right corner");
    kvDefString(extNames,"slice_brc.2", "brc.2");
  }

  return 1;
}

static void processHeader( FileHandler* self, KVHash* info, SList* cStack )
{
  /* Call the base class method */
  baseProcessHeader( self, info, cStack );

  if (!scan_LXI_header(info, self->fileName))
    Abort("lx_image_reader: unable to read or parse header from <%s>!\n",
	  self->fileName);
}

FileHandler* lxImageFactory(char* fname, KVHash* info)
{
  FileHandler* result= baseFactory(fname);
  result->processHeader= processHeader;
  result->typeName= strdup( "GE LX Image" );

  return result;
}

static int lxImageHeaderTest( const char* filename )
{
  FILE* f;
  char buf[8];
  int ierror= 0;

  if ((f = fopen(filename,"r"))!=NULL)
    {
      if (fseek(f, (long) FRZ_LXI_VALID_LOGO_OFF, SEEK_SET)) {
	perror("Error seeking header");
	ierror=1;
      }
      else {
	if (fread(buf,sizeof(char),sizeof(FRZ_LXI_VALID_LOGO), f)
	    != sizeof(FRZ_LXI_VALID_LOGO)) {
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
  else 
    return (!strncmp(buf,FRZ_LXI_VALID_LOGO,strlen(FRZ_LXI_VALID_LOGO))
	    || !strncmp(buf,FRZ_LXI_VALID_LOGO_BACKWARDS,
			strlen(FRZ_LXI_VALID_LOGO)));
}

int lxImageTester(const char* filename)
{
  return lxImageHeaderTest(filename);
}

