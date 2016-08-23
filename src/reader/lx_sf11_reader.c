/************************************************************
 *                                                          *
 *  lx_sf11_reader.c                                     *
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

#define SWAP(a,b) {tempr=(a);(a)=(b);(b)=tempr;}

static char rcsid[] = "$Id: lx_sf11_reader.c,v 1.10 2004/10/26 23:55:58 welling Exp $";

static int
IntBRdFloat32 (unsigned char *addr)
{
  return((int) Round(BRdFloat32(addr)));
}

static void
CalculateLocation (float *p1,
		   float *p2,
		   float *p3,
		   int rot,
		   int trans)
{
  /* 
    top left     - p1
    bottom left  - p2
    top right    - p3
    bottom right - p4
  */
  float p4[3], tempr;
  int i,n;

  /* calc missing corner (bottom right) */
  for (i=0; i<3; i++)
    p4[i] = p3[i] + (p2[i] - p1[i]);

  if (trans)
    for (i=0; i<3; i++)
      SWAP(p2[i],p3[i]);

  for (n=0; n<rot; n++)
    for (i=0; i<3; i++)
      {
	SWAP(p1[i],p3[i]);
	SWAP(p3[i],p4[i]);
	SWAP(p4[i],p2[i]);
      }

  /* reformat to match image header */
  for (i=0; i<3; i++)
    {
      SWAP(p2[i],p3[i]);
      SWAP(p3[i],p4[i]);
    }
}

static void calcSliceGeometry( KVHash* info,
			       unsigned char* rdbhead, 
			       unsigned char* acq_tab,
			       unsigned char* examhead, 
			       unsigned char* serieshead,
			       unsigned char* imagehead )
{
  float ctr[3];
  float p1[3];
  float p2[3];
  float p3[3];
  float n2[3];
  float n3[3];
  int i;
  int rotation;
  int transpose;
  float s_i_offset;
  float tempr2;
  float tempr3;
  float pix_shifth;
  float pix_shiftv;
  KVHash* defs= kvGetHash(info,"definitions");
  KVHash* extNames= kvGetHash(info,"external_names");

  ctr[0] = BRdFloat32(imagehead + FRZ_IMAGEHEAD_CTR_R_OFF);
  ctr[1] = BRdFloat32(imagehead + FRZ_IMAGEHEAD_CTR_A_OFF);
  ctr[2] = BRdFloat32(imagehead + FRZ_IMAGEHEAD_CTR_S_OFF);

  for (i = 0; i < 3; ++i)
    {
      p1[i]= 
	BRdFloat32(acq_tab + FRZ_DATAACQ_GW_POINT1_0_OFF + i*sizeof(float));
      p2[i]= 
	BRdFloat32(acq_tab + FRZ_DATAACQ_GW_POINT2_0_OFF + i*sizeof(float));
      p3[i]= 
	BRdFloat32(acq_tab + FRZ_DATAACQ_GW_POINT3_0_OFF + i*sizeof(float));
    }
  rotation = BRdInt16(rdbhead + FRZ_RDBHEAD_RDB_HDR_ROTATION_OFF);
  transpose = BRdInt16(rdbhead + FRZ_RDBHEAD_RDB_HDR_TRANSPOSE_OFF);
  CalculateLocation(p1, p2, p3, rotation, transpose);
  s_i_offset = BRdFloat32(imagehead + FRZ_IMAGEHEAD_TLHC_S_OFF) - p1[2];
  ctr[2] -= s_i_offset;

  /* calculate in-plane normal vectors n2,n3 */
  tempr2 = tempr3 = 0.0;
  for (i = 0; i < 3; i++)
    {
      n2[i] = p2[i] - p1[i];
      tempr2 += n2[i]*n2[i];
      n3[i] = p2[i] - p3[i];
      tempr3 += n3[i]*n3[i];
    }
  tempr2 = sqrt(tempr2);
  tempr3 = sqrt(tempr3);
  for (i = 0; i < 3; i++)
    {
      n2[i] /= tempr2;
      n3[i] /= tempr3;
    }

  pix_shifth = -(n2[0]*ctr[0] +  n2[1]*ctr[1] +  n2[2]*ctr[2])/tempr2;
  pix_shiftv = (n3[0]*ctr[0] +  n3[1]*ctr[1] +  n3[2]*ctr[2])/tempr2;
  kvDefDouble(info,"pix_shifth", pix_shifth);
  kvDefDouble(info,"pix_shiftv", pix_shiftv);
  kvDefString(defs,"pix_shifth","horizontal pixel shift (frac of FOV)");
  kvDefString(defs,"pix_shiftv","vertical pixel shift (frac of FOV)");
}


void ADD_VSUFFIX(scanSf11Header)( KVHash* info, 
				  unsigned char* rdbhead, 
				  unsigned char* acq_tab,
				  unsigned char* examhead, 
				  unsigned char* serieshead,
				  unsigned char* imagehead )
{
  int nodd;
  long slice_record_length;
  KVHash* defs= kvGetHash(info,"definitions");
  KVHash* extNames= kvGetHash(info,"external_names");
  
  kvDefInt(info,"nph1",
	   IntBRdFloat32(rdbhead + FRZ_RDBHEAD_RDB_HDR_USER4_OFF));
  kvDefString(defs,"nph1","# of baselines per slice");
  kvDefString(extNames,"nph1",""); /* supress output */
  kvDefInt(info,"nphmult",
	   IntBRdFloat32(rdbhead + FRZ_RDBHEAD_RDB_HDR_USER10_OFF));
  kvDefString(defs,"nphmult","# of images (phases) per baseline");
  kvDefString(extNames,"nphmult",""); /* supress output */

  /* slices get counted differently from standard LX */
  kvDefInt(info,"nslices",
	   BRdInt16(rdbhead+FRZ_RDBHEAD_RDB_HDR_NSLICES_OFF)
	   / (float)kvGetInt(info,"nph1"));
  

  kvDefInt(info,"nimages", kvGetInt(info,"nph1")*kvGetInt(info,"nphmult"));
  kvDefString(defs,"nimages","# of images in 1 Pfile (per coil)");
  kvDefString(extNames,"nimages",""); /* supress output; we have dt */
  kvDefInt(info,"ndatfr",
	   BRdInt16(rdbhead + FRZ_RDBHEAD_RDB_HDR_FRAME_SIZE_OFF));
  kvDefString(defs,"ndatfr","# of samples in each acq frame");
  kvDefInt(info,"ndat",kvGetInt(info,"ndatfr"));
  kvDefString(defs,"ndat","# of data items per projection");
  kvDefInt(info,"npr",IntBRdFloat32(rdbhead + FRZ_RDBHEAD_RDB_HDR_USER5_OFF));
  kvDefString(defs,"npr","# of projections (spirals)");
  kvDefInt(info,"chop",IntBRdFloat32(rdbhead+FRZ_RDBHEAD_RDB_HDR_USER7_OFF));
  kvDefString(defs,"chop","may be 1 or -1");
  kvDefDouble(info,"nex",BRdFloat32(rdbhead+FRZ_RDBHEAD_RDB_HDR_USER8_OFF));
  kvDefString(defs,"nex","# of excitations");
  kvDefInt(info,"densamp",
	   IntBRdFloat32(rdbhead+FRZ_RDBHEAD_RDB_HDR_USER0_OFF));
  kvDefString(defs,"densamp","???");
  /* This one seems to always want reordering */
  kvDefInt(info,"sliceorder",0);
  kvDefString(defs,"sliceorder","slice order: 0=interleaved, 1=sequential");
  kvDefInt(info,"reorder",1);
  kvDefString(info,"reorder_pattern","even/odd");
  kvDefInt(info,"gtype",0);
  kvDefString(defs,"gtype","k-space trajectory (0=spiral, 1=reverse)");
  kvDefDouble(info,"samp_time",
	      BRdFloat32(rdbhead + FRZ_RDBHEAD_RDB_HDR_USER2_OFF));
  kvDefString(defs,"samp_time","sample time (usec)");
  kvDefDouble(info,"fsgcm",
	      BRdFloat32(rdbhead + FRZ_RDBHEAD_RDB_HDR_USER12_OFF));
  kvDefString(defs,"fsgcm","design gradient max, G/cm");
  kvDefInt(info,"risetime",
	   IntBRdFloat32(rdbhead + FRZ_RDBHEAD_RDB_HDR_USER11_OFF));
  kvDefString(defs,"risetime","????");
  kvDefDouble(info,"opfov",
	      0.1*BRdFloat32(rdbhead + FRZ_RDBHEAD_RDB_HDR_USER6_OFF));
  kvDefString(defs,"opfov","field-of-view (cm)");

  kvDefDouble(info,"ts",kvGetDouble(info,"samp_time")*1e-6);
  kvDefString(defs,"ts","sample spacing in time (seconds)");
  kvDefDouble(info,"gts",
	      BRdFloat32(rdbhead + FRZ_RDBHEAD_RDB_HDR_USER3_OFF));
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
  
  calcSliceGeometry(info, rdbhead, acq_tab, examhead, serieshead, imagehead);

  /*
   **************
   * Finally we have the necessary info to proceed.
   **************
   */

  /* Some calculations are needed regarding the data layout.
   * Believe it or not, the defs below make things clearer 
   */
#define NPR kvGetInt(info,"npr")
#define NDAT kvGetInt(info,"ndat")
#define NPHMULT kvGetInt(info,"nphmult")
#define NPH1 kvGetInt(info,"nph1")
#define PT_SIZE kvGetInt(info,"pt_size")

  kvDefInt(info,"datatype_in",SRDR_INT16);
  kvDefString(info,"chunkname","samples");
  kvDefString(info,"dimstr","vpstbzc");
  kvDefInt(info,"dv",2);
  kvDefString(info,"description.v","compex real/imaginary");
  kvDefInt(info,"dp",NDAT);
  kvDefString(defs,"dp","# of complex samples");
  kvDefString(info,"description.p","ungridded k-space");
  kvDefInt(info,"ds",NPR);
  kvDefString(defs,"ds","# of shots");
  kvDefString(info,"description.s","discrete");
  kvDefInt(info,"dt",NPHMULT);
  kvDefString(defs,"dt","number of images per baseline");
  kvDefString(info,"description.t","gridded image-space");
  kvDefInt(info,"db",NPH1);
  kvDefString(defs,"db","# of baselines per slice");
  kvDefString(info,"description.b","discrete");
  kvDefInt(info,"dz",kvGetInt(info,"nslices"));
  kvDefString(info,"description.z","gridded image-space");
  kvDefInt(info,"dc",kvGetInt(info,"ncoils"));
  kvDefString(defs,"dc","# of coils");
  kvDefString(info,"description.c","discrete");

  nodd= NPHMULT*NPR % 2;

  kvDefLong(info, "baseline_record_length",
	    2*PT_SIZE * NDAT * ( NPHMULT*NPR + 1 + nodd ));
  kvDefString(defs,"baseline_record_length","# of bytes per baseline record");
	    
  kvDefLong(info, "slice_record_length",
	    NPH1 * kvGetLong(info,"baseline_record_length"));
  kvDefString(defs,"slice_record_length","# of bytes per slice record");

  /* Offset of the actual sample data is header plus 1 baseline */
  kvDefLong(info,"start_offset",
	    FRZ_POOL_HEADER_SIZE + (2 * PT_SIZE * NDAT));

  /* samples (p dimension) and shots (s dimension) are contiguous and
   * hence require no skips */

  /* baseline records (t dimension) have an extra empty field at the end
   * if nodd is true, plus we have to skip the baseline at the beginning
   * of the next baseline record.
   */
  kvDefLong(info,"skip.t", (nodd + 1) * 2 * PT_SIZE * NDAT);

  /* There appears to be no skip necessary between baseline records 
   * (the b dimension).  We may need a skip between coils, however.
   */

  kvDefLong(info,"skip.z",
	    kvGetInt(info,"coil_record_length")
	    - (kvGetInt(info,"dz")*kvGetLong(info,"slice_record_length")));
    
#undef NPR
#undef NDATFR
#undef NPHMULT
#undef NPH1
#undef PT_SIZE

}

void ADD_VSUFFIX(addSf11TrajChunk)(FileHandler* self, KVHash* info, 
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
  int nsamples;
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
   * some versions of getrttraj...() doesn't know exactly when to stop, 
   * and may write an extra sample into the t2k arrays.
   */
  t2k = Alloc3DFloatArray(2, ds, dp+2);
  nsamples = getrttrajvd(dp, ds, 
			 kvGetDouble(info,"ts"), 
			 kvGetDouble(info,"gts")*1e-6, 0.0, 
			 kvGetDouble(info,"fsgcm"), 
			 10.0*kvGetDouble(info,"opfov"),
			 kvGetInt(info,"risetime"), 
			 kvGetInt(info,"densamp"), 
			 1.0, 1.0, 
			 t2k[0], t2k[1], &res);

  kvDefInt(info,"opxres",res);
  kvDefDouble(info,"voxel_x",
	      kvGetDouble(info,"opfov")/kvGetInt(info,"opxres"));
  kvDefString(defs,"voxel_x","X voxel size (mm)");
  kvDefDouble(info,"voxel_y", kvGetDouble(info,"voxel_x"));
  kvDefString(defs,"voxel_y","Y voxel size (mm)");

  /* Sometimes the trajectory calculation calculates one too few
   * samples.  If that happens, we'll use a linear interpolation 
   * from the last sample to fill in the gap.
   */
  if (nsamples+ds < ds*dp) {
    if (nsamples+ds == ds*dp - ds) {
      Warning(1,"%s: extrapolating short spiral trajectory!\n",progname);
      for (s=0; s<ds; s++) {
	t2k[0][s][dp-1]= 2*t2k[0][s][dp-2] - t2k[0][s][dp-3];
	t2k[1][s][dp-1]= 2*t2k[1][s][dp-2] - t2k[1][s][dp-3];
      }
    }
#ifdef never
    else Abort("%s: Spiral traj routine returned too few samples (%d)!\n",
	       progname,nsamples);
#endif
    else {
      Warning(1,"%s: Spiral traj routine returned too few samples (%d)!\n",
	    progname,nsamples);
      if (ds != 1) 
	Abort("%s: Short trajectory hack is not implemented for multishot!\n",
	      progname);
      for (p=nsamples+1; p<dp; p++) {
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

