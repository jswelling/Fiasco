/************************************************************
 *                                                          *
 *  lx_splx_reader.c                                     *
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
/* This routine includes functionality from the "gspnnn" series of
 * spiral reconstruction codes, by Andy Stenger and others.
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
#include "rttraj.h"
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

static char rcsid[] = "$Id: lx_splx_reader.c,v 1.16 2004/10/26 23:55:58 welling Exp $";

static int
IntBRdFloat32 (unsigned char *addr)
{
  return((int) Round(BRdFloat32(addr)));
}


void ADD_VSUFFIX(scanSplxHeader)( KVHash* info, 
				  unsigned char* rdbhead, 
				  unsigned char* acq_tab,
				  unsigned char* examhead, 
				  unsigned char* serieshead,
				  unsigned char* imagehead )
{
  int itmp;
  int nodd;
  long slice_record_length;
  KVHash* defs= kvGetHash(info,"definitions");
  KVHash* extNames= kvGetHash(info,"external_names");
  
  kvDefInt(info,"nph1",1);
  kvDefString(defs,"nph1","# of baselines per slice");
  kvDefString(extNames,"nph1",""); /* supress output */
  itmp= IntBRdFloat32(rdbhead + FRZ_RDBHEAD_RDB_HDR_USER1_OFF);
  if (itmp<1) itmp= 1;
  kvDefInt(info,"nphmult",itmp);
  kvDefString(defs,"nphmult","# of images (phases) per baseline");
  kvDefString(extNames,"nphmult",""); /* supress output */
  kvDefInt(info,"nimages", kvGetInt(info,"nph1")*kvGetInt(info,"nphmult"));
  kvDefString(defs,"nimages","# of images in 1 Pfile (per coil)");
  kvDefString(extNames,"nimages",""); /* supress output */
  kvDefInt(info,"ndatfr",BRdInt16(rdbhead + FRZ_RDBHEAD_RDB_HDR_FRAME_SIZE_OFF));
  kvDefString(defs,"ndatfr","# of samples in each acq frame");
  kvDefInt(info,"npr",IntBRdFloat32(rdbhead + FRZ_RDBHEAD_RDB_HDR_USER4_OFF));
  kvDefString(defs,"npr","# of projections (spirals)");
  kvDefInt(info,"chop",1);
  kvDefString(defs,"chop","may be 1 or -1");
  kvDefDouble(info,"nex",1);
  kvDefString(defs,"nex","# of excitations");
  kvDefInt(info,"densamp",5);
  kvDefString(defs,"densamp","???");

  if (!strncmp(kvGetString(info,"pulse_seq"),"splxr",5)) {
    /* user16 is nzpe; always interleaved */
    kvDefInt(info,"nzpe",
	     IntBRdFloat32(rdbhead + FRZ_RDBHEAD_RDB_HDR_USER16_OFF));
    kvDefInt(info,"sliceorder",0);
  }
  else if (!strncmp(kvGetString(info,"pulse_seq"),"splxcnv",7)) {
    /* user16 *should* be slice order, but isn't. nzpe is 1. */
    kvDefInt(info,"nzpe",1);
    kvDefInt(info,"sliceorder",0);
  }
  else if (!strncmp(kvGetString(info,"pulse_seq"),"splxvh3",7)) {
    /* Always interleaved; nzpe is 1. */
    kvDefInt(info,"nzpe",1);
    kvDefInt(info,"sliceorder",0);
  }
  else {
    /* user16 is slice order, nzpe is 1 */
    kvDefInt(info,"sliceorder",
	     IntBRdFloat32(rdbhead + FRZ_RDBHEAD_RDB_HDR_USER16_OFF));
    kvDefInt(info,"nzpe",1);
  }
  kvDefString(defs,"nzpe","Z projections (3d acquisition)");
  kvDefString(defs,"sliceorder","slice order: 0=interleaved, 1=sequential");

  if (kvGetInt(info,"sliceorder")==0) {
    kvDefInt(info,"reorder",1);
    kvDefString(info,"reorder_pattern","even/odd");
  }
  else {
    kvDefInt(info,"reorder",1);
    kvDefString(info,"reorder_pattern","even/odd");
  }

  kvDefInt(info,"gtype",IntBRdFloat32(rdbhead + FRZ_RDBHEAD_RDB_HDR_USER5_OFF));
  kvDefString(defs,"gtype","k-space trajectory (0=spiral, 1=reverse)");
  if (kvGetInt(info,"gtype")==1) {
    kvDefBoolean(info,"reverse",1);
    kvDefString(defs,"reverse","true for reversed sample order");
  }
  kvDefInt(info,"opxres",IntBRdFloat32(rdbhead + FRZ_RDBHEAD_RDB_HDR_USER3_OFF));
  kvDefString(defs,"opxres","nominal output resolution in image space");
  kvDefDouble(info,"slewrate",BRdFloat32(rdbhead + FRZ_RDBHEAD_RDB_HDR_USER7_OFF));
  kvDefString(defs,"slewrate","design slew rate, (mT/m/ms)");
  kvDefInt(info,"ngap",IntBRdFloat32(rdbhead + FRZ_RDBHEAD_RDB_HDR_USER10_OFF));
  kvDefString(defs,"ngap","??? for concat readout");
  kvDefInt(info,"concat",IntBRdFloat32(rdbhead + FRZ_RDBHEAD_RDB_HDR_USER11_OFF));
  kvDefString(defs,"concat","??? for concat readout");
  kvDefInt(info,"ndat",kvGetInt(info,"ndatfr")*(kvGetInt(info,"concat")+1));
  kvDefString(defs,"ndat","# of data items per projection");
  kvDefDouble(info,"fast_rec_lpf",BRdFloat32(rdbhead + FRZ_RDBHEAD_RDB_HDR_USER12_OFF));
  kvDefString(defs,"fast_rec_lpf","fast rec cutoff in kHz");
  kvDefDouble(info,"mapdel",BRdFloat32(rdbhead + FRZ_RDBHEAD_RDB_HDR_USER15_OFF));
  kvDefString(defs,"mapdel","field mapping offset (us)");
  kvDefDouble(info,"samp_time",BRdFloat32(rdbhead + FRZ_RDBHEAD_RDB_HDR_USER13_OFF));
  kvDefString(defs,"samp_time","sample time (usec)");
  kvDefDouble(info,"fsgcm",BRdFloat32(rdbhead + FRZ_RDBHEAD_RDB_HDR_USER6_OFF));
  kvDefString(defs,"fsgcm","design gradient max, G/cm");
  kvDefInt(info,"risetime",(int)rint(kvGetDouble(info,"fsgcm")*10000.0
				     / kvGetDouble(info,"slewrate")));
  kvDefString(defs,"risetime","????");
  kvDefDouble(info,"opfov",BRdFloat32(rdbhead + FRZ_RDBHEAD_RDB_HDR_USER0_OFF));
  kvDefString(defs,"opfov","field-of-view (cm)");

  kvDefDouble(info,"voxel_x",
	      10.0*kvGetDouble(info,"opfov")/kvGetInt(info,"opxres"));
  kvDefString(defs,"voxel_x","X voxel size (mm)");
  kvDefDouble(info,"voxel_y", kvGetDouble(info,"voxel_x"));
  kvDefString(defs,"voxel_y","Y voxel size (mm)");

  kvDefDouble(info,"ts",kvGetDouble(info,"samp_time")*1e-6);
  kvDefString(defs,"ts","sample spacing in time (seconds)");
  kvDefDouble(info,"gts",4.0);
  kvDefString(defs,"gts","gradient time spacing (msec)");
  kvDefString(extNames,"gts","grad_time_spacing");

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

  kvDefInt(info,"nfiles",1);
  kvDefString(defs,"nfiles","# of input files");
  kvDefInt(info,"total_images",
	   kvGetInt(info,"nimages")*kvGetInt(info,"nfiles"));
  kvDefString(defs,"total_images","# of images over all P files");
  
  if (kvGetBoolean(info,"reverse")) {
    kvDefDouble(info,"pix_shifth",
		BRdFloat32(rdbhead + FRZ_RDBHEAD_RDB_HDR_YOFF_OFF)/128.0);
    kvDefDouble(info,"pix_shiftv",
		-BRdFloat32(rdbhead + FRZ_RDBHEAD_RDB_HDR_XOFF_OFF)/128.0);
  }
  else {
    kvDefDouble(info,"pix_shifth",
		BRdFloat32(rdbhead + FRZ_RDBHEAD_RDB_HDR_YOFF_OFF)
		/ BRdInt16(rdbhead + FRZ_RDBHEAD_RDB_HDR_IM_SIZE_OFF));
    kvDefDouble(info,"pix_shiftv",
		-BRdFloat32(rdbhead + FRZ_RDBHEAD_RDB_HDR_XOFF_OFF)
		/ BRdInt16(rdbhead + FRZ_RDBHEAD_RDB_HDR_IM_SIZE_OFF));
  }

  /*
   **************
   * Finally we have the necessary info to proceed.
   **************
   */

  /* nzpe is number of Fourier samples in slice direction, a new
   * feature with the "gsp3d" series of codes.  If it is unset, that
   * indicates an older code, which will be a 2D acquisition.
   */
  kvDefString(info,"chunkname","samples");
  if (kvGetInt(info,"nzpe") == 1
      || kvGetInt(info,"nzpe") == 0) {
    /* 2D acquisition */
    kvDefString(info,"dimstr","vpstzc");
  }
  else {
    /* 3D acquisition */
    kvDefString(info,"dimstr","vpsktzc");
    kvDefInt(info,"dk",kvGetInt(info,"nzpe"));
    kvDefString(defs,"dk","# of Fourier samples in Z");
    kvDefString(info,"description.k","gridded k-space");
    
    Warning(1,
	 "%s: scan_splx_header: this looks like 3D acquisition (nxpe= %d)!\n",
	    progname,kvGetInt(info,"nzpe"));
  }

  /* Some calculations are needed regarding the data layout.
   * Believe it or not, the defs below make things clearer 
   */
#define NPR kvGetInt(info,"npr")
#define NZPE kvGetInt(info,"nzpe")
#define CONCAT kvGetInt(info,"concat")
#define NDATFR kvGetInt(info,"ndatfr")
#define NPHMULT kvGetInt(info,"nphmult")

  nodd= NPHMULT*NPR*NZPE*(CONCAT+1) % 2;
  kvDefLong(info, "slice_record_length",
	    4 * NDATFR * ( ((CONCAT+1) * NPHMULT * NPR * NZPE) + (nodd+1) ));

#undef NPR
#undef NZPE
#undef CONCAT
#undef NDATFR
#undef NPHMULT

  kvDefString(defs,"slice_record_length","# of bytes in one slice record");
  kvDefInt(info,"dv",2);
  kvDefString(info,"description.v","complex real/imaginary");
  kvDefInt(info,"dp",kvGetInt(info,"ndat"));
  kvDefString(defs,"dp","# of complex samples");
  kvDefString(info,"description.p","ungridded k-space");
  kvDefInt(info,"ds",kvGetInt(info,"npr"));
  kvDefString(defs,"ds","# of shots");
  kvDefString(info,"description.s","discrete");
  kvDefInt(info,"dc",kvGetInt(info,"ncoils"));
  kvDefString(defs,"dc","# of coils");
  kvDefString(info,"description.c","discrete");
  kvDefInt(info,"dz",kvGetInt(info,"nslices"));
  kvDefString(info,"description.z","gridded image-space");
  kvDefInt(info,"dt",kvGetInt(info,"total_images"));
  kvDefString(info,"description.t","gridded image-space");
  
  /* Offset of the actual sample data */
  kvDefLong(info,"pool_header_size", FRZ_POOL_HEADER_SIZE);
  kvDefString(defs,"pool_header_size", "size of LX file header");
  kvDefLong(info,"start_offset",
	    FRZ_POOL_HEADER_SIZE + 
	    (2 * kvGetInt(info,"pt_size") * kvGetInt(info,"ndatfr")));
  
  kvDefLong(info,"skip.t", 
	    kvGetLong(info,"slice_record_length")
	    -(srdrTypeSize[kvGetInt(info,"datatype_in")]
	      *kvGetInt(info,"dv")*kvGetInt(info,"dp")*kvGetInt(info,"ds")
	      *kvGetInt(info,"dt")));
  kvDefLong(info,"skip.z",
	    kvGetInt(info,"coil_record_length")
	    - (kvGetInt(info,"dz")*kvGetLong(info,"slice_record_length")));
    
}

void ADD_VSUFFIX(addSplxTrajChunk)(FileHandler* self, KVHash* info, 
				   SList* cStack)
{
  /* Add chunks for x and y sample locations */
  KVHash* defs= kvGetHash(info,"definitions");
  ChunkHandlerPair* chPair= NULL;
  KVHash* subInfo= kvFactory(KV_DEFAULT_SIZE);
  float*** t2k;
  float* sampx= NULL;
  float* sampy= NULL;
  float* xrunner;
  float* yrunner;
  int res;
  int ds= kvGetInt(info,"ds");
  int dp= kvGetInt(info,"dp");
  int dc= kvGetInt(info,"dc");
  int dz= kvGetInt(info,"dz");
  int s;
  int p;
  int c;
  int z;
  
  /* Calculate the trajectory.  The extra space on the end of the
   * last dimension of t2k is to cover for the fact that 
   * getrttrajghg() doesn't know exactly when to stop, and
   * may write an extra sample into the t2k arrays.
   */
  t2k = Alloc3DFloatArray(2, ds, dp+2);
  res = getrttrajghg(kvGetInt(info,"opxres"), ds, 
		     kvGetDouble(info,"ts"), 
		     kvGetDouble(info,"gts")*1e-6,
		     kvGetDouble(info,"fsgcm"), kvGetDouble(info,"opfov"),
		     kvGetDouble(info,"slewrate"), 
		     kvGetDouble(info,"gts")*21000*1e-6, 
		     kvGetInt(info,"gtype"), 
		     t2k[0], t2k[1]);

  /* Sometimes the trajectory calculation calculates one too few
   * samples.  If that happens, we'll use a linear interpolation 
   * from the last sample to fill in the gap.
   */
  if (res+ds < ds*dp) {
    if (res+ds == ds*dp - ds) {
      Warning(1,"%s: extrapolating short spiral trajectory!\n",progname);
      for (s=0; s<ds; s++) {
	t2k[0][s][dp-1]= 2*t2k[0][s][dp-2] - t2k[0][s][dp-3];
	t2k[1][s][dp-1]= 2*t2k[1][s][dp-2] - t2k[1][s][dp-3];
      }
    }
#ifdef never
    else Abort("%s: Spiral traj routine returned too few samples (%d)!\n",
	       progname,res);
#endif
    else {
      Warning(1,"%s: Spiral traj routine returned too few samples (%d)!\n",
	    progname,res);
      if (ds != 1) Abort("%s: Short trajectory hack is not implemented for multishot!\n",
			 progname);
      for (p=res+1; p<dp; p++) {
	t2k[0][0][p]= 0.0;
	t2k[1][0][p]= 0.0;
      }
    }
  }
  
  /* make copies of the trajectory to hand off.  For reversed spirals
   * we store the trajectory points in reverse order, to correspond with
   * the storage order of the samples.
   */
  if (!(sampx= (float*)malloc(ds*dp*dc*dz*sizeof(float))))
    Abort("%s: unable to allocate %d bytes!\n",progname,ds*dp*dc*dz*sizeof(float));
  if (!(sampy= (float*)malloc(ds*dp*dc*dz*sizeof(float))))
    Abort("%s: unable to allocate %d bytes!\n",progname,ds*dp*dc*dz*sizeof(float));
  xrunner= sampx;
  yrunner= sampy;
  for (z=0; z<dz; z++)
    for (c=0; c<dc; c++)
      for (s=0; s<ds; s++) {
	if (kvGetBoolean(info,"reverse")) {
	  for (p=0; p<dp; p++) {
	    *xrunner++= t2k[0][s][dp - (p+1)];
	    *yrunner++= t2k[1][s][dp - (p+1)];
	  }
	}
	else {
	  for (p=0; p<dp; p++) {
	    *xrunner++= t2k[0][s][p];
	    *yrunner++= t2k[1][s][p];
	  }
	}
      }
  
  /* and clean up */
  Free3DFloatArray(t2k);
  
  initInfoHash(subInfo);
  kvDefString(subInfo,"chunkname","sample_kxloc");
  kvDefString(subInfo,"chunkfile",".sampxloc");
  kvDefString(subInfo,"dimstr","pscz");
  kvDefInt(subInfo,"dp",kvGetInt(info,"dp"));
  kvDefString(subInfo,"description.p","ungridded k-space");
  kvDefInt(subInfo,"ds",kvGetInt(info,"ds"));
  kvDefString(subInfo,"description.s","discrete");
  kvDefInt(subInfo,"dc",kvGetInt(info,"dc"));
  kvDefString(subInfo,"description.c","discrete");
  kvDefInt(subInfo,"dz",kvGetInt(info,"dz"));
  kvDefString(subInfo,"description.z","gridded image-space");
  kvDefLong(subInfo,"start_offset",0);
  kvDefInt(subInfo,"datatype_in",SRDR_FLOAT32);
  kvDefInt(subInfo,"handler_datatype_out",SRDR_FLOAT32);
  
  if (!(chPair= (ChunkHandlerPair*)malloc(sizeof(ChunkHandlerPair))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,sizeof(ChunkHandlerPair));
  chPair->info= subInfo;
  chPair->handler= ramDataHandlerFactory(sampx, dp*ds*dc*dz, SRDR_FLOAT32);
  slist_push(cStack,chPair);
  
  chPair= NULL;
  subInfo= kvFactory(KV_DEFAULT_SIZE);
  
  initInfoHash(subInfo);
  kvDefString(subInfo,"chunkname","sample_kyloc");
  kvDefString(subInfo,"chunkfile",".sampyloc");
  kvDefString(subInfo,"dimstr","pscz");
  kvDefInt(subInfo,"dp",kvGetInt(info,"dp"));
  kvDefString(subInfo,"description.p","ungridded k-space");
  kvDefInt(subInfo,"ds",kvGetInt(info,"ds"));
  kvDefString(subInfo,"description.s","discrete");
  kvDefInt(subInfo,"dc",kvGetInt(info,"dc"));
  kvDefString(subInfo,"description.c","discrete");
  kvDefInt(subInfo,"dz",kvGetInt(info,"dz"));
  kvDefString(subInfo,"description.z","gridded image-space");
  kvDefLong(subInfo,"start_offset",0);
  kvDefInt(subInfo,"datatype_in",SRDR_FLOAT32);
  kvDefInt(subInfo,"handler_datatype_out",SRDR_FLOAT32);
  
  if (!(chPair= (ChunkHandlerPair*)malloc(sizeof(ChunkHandlerPair))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,sizeof(ChunkHandlerPair));
  chPair->info= subInfo;
  chPair->handler= ramDataHandlerFactory(sampy, dp*ds*dc*dz, SRDR_FLOAT32);
  slist_push(cStack,chPair);
}

