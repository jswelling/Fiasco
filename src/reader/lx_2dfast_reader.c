/************************************************************
 *                                                          *
 *  lx_2dfast_reader.c                                     *
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
/* This routine uses functionality from the "epirecon" code
 * supplied with the GE LX2 scanner system, author 
 * Bryan J. Mock (GE Medical Systems).
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "mri.h"
#if (SGI64 || SGI5 || SGIMP)
#include <bstring.h>
#endif

#include "bio.h"
#include "fmri.h"
#include "array.h"
#include "stdcrg.h"
#include "misc.h"
#include "smartreader.h"
#include "nr_sub.h"
#include "frozen_header_info.h"

/* LX header files */
#if defined(USE_RDBM_LX2)
#define ADD_VSUFFIX(x) x ## _lx2
#elif defined(USE_RDBM_PRELX)
#define ADD_VSUFFIX(x) x ## _prelx
#elif defined(USE_RDBM_EXCITE)
#define ADD_VSUFFIX(x) x ## _excite
#else
#define ADD_VSUFFIX(x) x ## _cnv4
#endif

static char rcsid[] = "$Id: lx_2dfast_reader.c,v 1.4 2004/10/26 23:55:58 welling Exp $";

#define BLANKY_NUM_DEFAULT 2
#define BLANKY_NUM_MAX 10

#define GE_RESAMPLE_SCRIPT "ge_ramp_resample.csh"

static int
IntBRdFloat32 (unsigned char *addr)
{
  return((int) Round(BRdFloat32(addr)));
}

void ADD_VSUFFIX(scan2dfastHeader)( KVHash* info, 
				     unsigned char* rdbhead, 
				     unsigned char* acq_tab,
				     unsigned char* examhead, 
				     unsigned char* serieshead,
				     unsigned char* imagehead )
{
  /* Read from header.  A two-step translation then happens:
   * from the header fields to the language of Mock's "epirecon",
   * then to the language of this program.
   */
  KVHash* defs= kvGetHash(info,"definitions");
  
  float bandwidth;       /* reciever bandwidth from header */
  char plane[4];         /* if pl = L/R > axial or A/P > sag, */
  /* S/I > coronal */
  int xres, yres;        /* acq. resolution */
  int rcxres,rcyres;     /* reconstructed resolution (with fovar=1.0) */
  float fovar;            /* Field of View Aspect Ratio from header */
  int numreps;            /* total reps and rep number (nex in header) */
  int frame;             /* size of each image points calc. from */
  /* header info */
  int ssp_size;          /* Adding these together (+ header) should */
  int pt_size = 0;        /* For extended Dynamic Range data */
  int rhraw_size;        /* give expected file size in bytes        */
  float hnw;             /* homodyne correction transition width in */
  /* header */
  float fw, fr;           /* fermi width and radius from header */
  char pl;                /* orientation of slice */
  short swapfp;           /* swapped freq/phase encode direction 1 = */
  /* yes - determines if images are reoriented */
  short rot, tpose;
  int blanky_num;         /* # of y lines to blank @ top & bottom of image */
  int bviews = 0;         /* # of baseline views */
  int ileaves = 1;        /* # of interleaves if multi-shot EPI */
  int fast_rec = 0;       /* fast_reciever check */
  int vpsht;              /* Ky views per shot */
  int interm_xres;        /* reconstructed x resolution (with fovar != 1.0)*/
  
  /* Offset of the actual sample data */
  kvDefLong(info,"start_offset",
	    FRZ_POOL_HEADER_SIZE + 
	    (2 * kvGetInt(info,"pt_size")
	     * (kvGetInt(info,"bviews")+kvGetInt(info,"bl_save"))));

  /* Set acquired matrix size xres,yres, reconstruction size rcxres, */
  /* rcyres, number of slices, frame size of each image, and number of */
  /* reps (nex) from rawheader */
  xres = BRdInt16(rdbhead+FRZ_RDBHEAD_RDB_HDR_DA_XRES_OFF);
  yres = BRdInt16(rdbhead+FRZ_RDBHEAD_RDB_HDR_DA_YRES_OFF)-1;
  fovar = BRdFloat32(rdbhead+FRZ_RDBHEAD_RDB_HDR_PHASE_SCALE_OFF);
  
  /* Other variables in the nomenclature of Mock's "epirecon" */
  rcxres = BRdInt16(rdbhead+FRZ_RDBHEAD_RDB_HDR_RC_XRES_OFF);
  interm_xres = (int)((float)rcxres/fovar);
  rcyres = BRdInt16(rdbhead+FRZ_RDBHEAD_RDB_HDR_RC_YRES_OFF);
  numreps = BRdFloat32(imagehead+FRZ_IMAGEHEAD_NEX_OFF);
  /* allow for ext dyn. range data */
  pt_size = kvGetInt(info,"pt_size");
  frame = xres*yres*pt_size;
  ssp_size = BRdInt32(rdbhead+FRZ_RDBHEAD_RDB_HDR_SSPSAVE_OFF);
  rhraw_size = 
    BRdInt32(rdbhead+FRZ_RDBHEAD_RDB_HDR_RAW_PASS_SIZE_OFF);
  hnw = BRdFloat32(rdbhead+FRZ_RDBHEAD_RDB_HDR_NTRAN_OFF);
  fw = BRdFloat32(rdbhead+FRZ_RDBHEAD_RDB_HDR_FERMI_WIDTH_OFF);
  fr = BRdFloat32(rdbhead+FRZ_RDBHEAD_RDB_HDR_FERMI_RADIUS_OFF);
  pl = *((char*)imagehead + FRZ_IMAGEHEAD_LOC_RAS_OFF);
  swapfp = BRdInt16(imagehead+FRZ_IMAGEHEAD_SWAPPF_OFF);
  bandwidth = BRdFloat32(imagehead+FRZ_IMAGEHEAD_VBW_OFF);
  blanky_num = BRdInt16(rdbhead+FRZ_RDBHEAD_RDB_HDR_SLBLANK_OFF);
  bviews = kvGetInt(info,"bviews");
  fast_rec = BRdInt32(rdbhead+FRZ_RDBHEAD_RDB_HDR_FAST_REC_OFF);
  ileaves = BRdInt16(rdbhead+FRZ_RDBHEAD_RDB_HDR_ILEAVES_OFF);
  
  /* set number of blank lines to 2 if out of bounds */
  if( (blanky_num == 0) || (blanky_num > BLANKY_NUM_MAX) )
    blanky_num = BLANKY_NUM_DEFAULT;
  
  /* number of views per shot */
  if(ileaves != 0) vpsht= yres/ileaves;
  else vpsht= yres; /* presumably single shot data */
  
  /* BJM: safety check in case BW field in image header is NULL */
  if(bandwidth == 62.0)
    bandwidth+= 0.5;
  else if (bandwidth == 0.0) /* assume GE epibold.e */
    bandwidth = BRdFloat32(rdbhead+FRZ_RDBHEAD_RDB_HDR_USER16_OFF);
  else if (bandwidth == 0.0) { /* if STILL Zero.. */
    Warning(1,"%s: Warning: Reciever bandwidth was ZERO\n", progname);
    if (kvLookup(info,"bandpassdir")) {
      Warning(1,
	      "%s: Not performing band pass asymmetry correction\n",progname);
      kvDeleteAll(info,"bandpassdir");
    }
  }
  kvDefDouble(info,"bandwidth",bandwidth);
  kvDefString(defs,"bandwidth","receiver bandwdith in MHz");
  
  kvDefDouble(info,"fermi_width",fw);
  kvDefString(defs,"fermi_width","Fermi filter width");
  kvDefDouble(info,"fermi_radius",fr);
  kvDefString(defs,"fermi_radius","Fermi filter radius");

  kvDefDouble(info,"fov_x",BRdFloat32(imagehead+FRZ_IMAGEHEAD_DFOV_OFF));
  kvDefString(defs,"fov_x","X field of view (mm)");
  
  if(fovar == 0.5) { 
    kvDefDouble(info,"fov_y",
		BRdFloat32(imagehead+FRZ_IMAGEHEAD_DFOV_RECT_OFF)/fovar);
  }
  else {
    kvDefDouble(info,"fov_y",
		BRdFloat32(imagehead+FRZ_IMAGEHEAD_DFOV_RECT_OFF));
  }
  kvDefString(defs,"fov_y","Y field of view (mm)");

  kvDefDouble(info,"voxel_x",
	      BRdFloat32(imagehead+FRZ_IMAGEHEAD_PIXSIZE_X_OFF));
  kvDefString(defs,"voxel_x","X voxel size (mm)");
  kvDefDouble(info,"voxel_y",
	      BRdFloat32(imagehead+FRZ_IMAGEHEAD_PIXSIZE_Y_OFF));
  kvDefString(defs,"voxel_y","Y voxel size (mm)");
  kvDefDouble(info,"image_x",BRdFloat32(imagehead+FRZ_IMAGEHEAD_DIM_X_OFF));
  kvDefDouble(info,"image_y",BRdFloat32(imagehead+FRZ_IMAGEHEAD_DIM_Y_OFF));
  kvDefInt(info,"overscan",BRdInt16(rdbhead+FRZ_RDBHEAD_RDB_HDR_HNOVER_OFF));
  
  if (pt_size<4) kvDefInt(info,"datatype_in",SRDR_INT16);
  else kvDefInt(info,"datatype_in",SRDR_INT32);
  kvDefInt(info,"dv",2);
  kvDefString(info,"description.v","complex real/imaginary");

  kvDefInt(info,"dy",vpsht);
  kvDefString(info,"description.y","gridded k-space index");
  kvDefInt(info,"dy_base",rcyres);
  kvDefString(defs,"dy_base","y samples after reconstruction and clipping");
  if (rcyres == yres) {
    kvDefBoolean(info,"partialk",0);
  }
  else {
    kvDefBoolean(info,"partialk",0);
  }
  kvDefString(defs,"partialk","partial-k completion needed");

  if (ileaves != 0) kvDefInt(info,"ds",ileaves);
  else kvDefInt(info,"ds",1);
  kvDefInt(info,"dt",numreps);
  kvDefString(info,"description.t","gridded image-space index");

  kvDefInt(info,"ncoils",
	   BRdInt16(rdbhead + FRZ_RDBHEAD_RDB_HDR_DAB_0_STOP_RCV_OFF) -
	   BRdInt16(rdbhead + FRZ_RDBHEAD_RDB_HDR_DAB_0_START_RCV_OFF) + 1);
  kvDefString(defs,"ncoils","# of coils");
  kvDefInt(info,"coil_record_length",
	   BRdInt32(rdbhead + FRZ_RDBHEAD_RDB_HDR_RAW_PASS_SIZE_OFF) 
	   / kvGetInt(info,"ncoils"));
  kvDefString(defs,"coil_record_length","# of bytes in one coil record");
  kvDefInt(info,"baseline_length",4*kvGetInt(info,"ndat"));
  kvDefString(defs,"baseline_length","# of bytes in one baseline record");
  if (kvGetInt(info,"ncoils")>1)
    kvDefInt(info,"dc",kvGetInt(info,"ncoils"));

  /* For some reason, these files are broken up into groups of at most
   * 7 slices.
   */
  if (kvGetInt(info,"nslices") > 7)
    kvDefInt(info,"dz",7);
  else 
    kvDefInt(info,"dz",kvGetInt(info,"nslices"));
  kvDefString(info,"description.z","gridded image-space index");
  
  if (rcxres != xres) {
    /* ramp resampling required */
    if (!kvLookup(info,"rampfile")) 
      Warning(1,"%s: scan_2dfast_header: ramp resampling needed, but no ramp file given!\n",
	      progname);
    kvDefInt(info,"dx_resampled",rcxres);
    kvDefInt(info,"dq",xres);
    kvDefString(info,"description.q","ungridded k-space index");
    kvDefBoolean(info,"resample",1);
    kvDefString(defs,"resample","regridding for ramp sampling required");
    kvDefString(info,"resample_method",GE_RESAMPLE_SCRIPT);    
    if (kvLookup(info,"dc") && kvGetInt(info,"dc")>1) {
      if (kvLookup(info,"ds") && kvGetInt(info,"ds")>1) {
	kvDefString(info,"dimstr","vqyszct");
      }
      else {
	kvDefString(info,"dimstr","vqyzct");
      }
    }
    else {
      if (kvLookup(info,"ds") && kvGetInt(info,"ds")>1) {
	kvDefString(info,"dimstr","vqyszt");
      }
      else {
	kvDefString(info,"dimstr","vqyzt");
      }
    }
  }
  else {
    kvDefInt(info,"dx",xres);
    kvDefString(info,"description.x","gridded k-space index");
    if (kvLookup(info,"dc") && kvGetInt(info,"dc")>1) {
      if (kvLookup(info,"ds") && kvGetInt(info,"ds")>1) {
	kvDefString(info,"dimstr","vxyszct");
      }
      else {
	kvDefString(info,"dimstr","vxyzct");
      }
    }
    else {
      if (kvLookup(info,"ds") && kvGetInt(info,"ds")>1) {
	kvDefString(info,"dimstr","vxyszt");
      }
      else {
	kvDefString(info,"dimstr","vxyzt");
      }
    }
  }

  kvDefBoolean(info,"reorder",1);
  kvDefString(info,"reorder_pattern","even/odd");
  kvDefLong(info,"skip",0);
  kvDefLong(info,"sliceskip",0);  

  kvDefBoolean(info,"rowflip",0);
  kvDefString(defs,"rowflip","EPI row reversal needed");
  kvDefString(info,"rowflip_pattern","none");
  kvDefString(defs,"rowflip_pattern","EPI row reversal pattern");

  if (kvGetInt(info,"overscan")==0) {
    kvDefBoolean(info,"ychop",1);
    kvDefInt(info,"kspace_ctr.y",0);
  }
  else {
    kvDefBoolean(info,"ychop",1);
    kvDefInt(info,"kspace_ctr.y",yres-kvGetInt(info,"overscan"));
  }
  kvDefString(defs,"kspace_ctr.y", "this row crosses the k-space origin");
    
  kvDefInt(info,"kspace_ctr.x",xres/2);
  kvDefString(defs,"kspace_ctr.x", "this column crosses the k-space origin");
  
}

