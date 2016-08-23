/************************************************************
 *                                                          *
 *  lx_reader.c                                             *
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
/* Notes-
   ---SPLX---
   -I'm still not reading all the damn coordinate info from the file.
    Some of the missing stuff is actually used in reconstruction.
   -splx should set skip and sliceskip, or epibold should not.
   ---EPIBOLD---
   -Is their homodyne recon really equiv to our partialk?
   -Header provides info necessary to reorient sagital and coronal images,
    but not oblique.
   -What's the blank lines stuff in epirecon?
   -Am I supporting partial FOV scanning?  Warn user either way.
   -I *think* twoshot is in the input with successive slices for
    the same shot consecutively (then moving on to the next shot).
   -What's "dc zipper" effect that baseline data is supposed to get
    rid of for fast receiver case?
   ---2DFAST---
   -Gee, I really should learn something about this pulse sequence.
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
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

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

static char rcsid[] = "$Id: lx_reader.c,v 1.28 2004/10/26 23:55:58 welling Exp $";

/* We actually seem to have some floating point precision issues in
 * distinguishing the rdbm header revision numbers.
 */
#define REVISION_EPSILON 0.00001

/* Number of times to try shifting reading frame to find a valid header */
#define MAX_FRAME_SHIFTS 5

/* This is the range of the FFT of the first image with autoscaling on */
#define AUTOSCALE_RANGE 8192.0

static int
IntBRdFloat32 (unsigned char *addr)
{
  return((int) Round(BRdFloat32(addr)));
}

static void calcVolumeBounds(KVHash* info)
{
  double slice_tlc[3];
  double slice_brc[3];
  double slice_trc[3];
  double slice_blc[3];
  double slice_ctr[3];
  double tlf[3];
  double trf[3];
  double tlb[3];
  double trb[3];
  double blf[3];
  double brf[3];
  double blb[3];
  double brb[3];
  double delta[3];
  double x_norm[3];
  double y_norm[3];
  double z_norm[3];
  double xscale;
  double yscale;

  /* This routine will probably have to be updated as we learn more 
   * about their notion of coordinate systems.
   */

  getVec3(info,"slice_tlc",slice_tlc);
  getVec3(info,"slice_trc",slice_trc);
  getVec3(info,"slice_brc",slice_brc);
  getVec3(info,"slice_ctr",slice_ctr);

  subtractVec3(delta,slice_ctr,slice_trc);
  xplusbyVec3(slice_blc,slice_trc,delta,2.0);
  defVec3(info,"slice_blc",slice_blc);

  subtractVec3(delta,slice_trc,slice_tlc);
  xscale= normVec3(delta);
  if (xscale==0.0) {
    Warning(1,"%s: lx_reader: nonsense slice corner information!\n",progname);
    return;
  }
  multVec3( x_norm, delta, 1.0/xscale );

  subtractVec3(delta,slice_trc,slice_brc);
  yscale= normVec3(delta);
  if (yscale==0.0) {
    Warning(1,"%s: lx_reader: nonsense slice corner information!\n",progname);
    return;
  }
  multVec3( y_norm, delta, 1.0/yscale );

  crossVec3( z_norm, x_norm, y_norm );

#ifdef never
  fprintf(stderr,"Norms: (%f %f %f) (%f %f %f) (%f %f %f)\n",
	  x_norm[0], x_norm[1], x_norm[2],
	  y_norm[0], y_norm[1], y_norm[2],
	  z_norm[0], z_norm[1], z_norm[2]);
#endif
  defVec3(info,"slice_norm",z_norm);

  copyVec3(blb, slice_tlc);
  copyVec3(brb, slice_trc);
  copyVec3(brf, slice_brc);
  copyVec3(blf, slice_blc);
  if ( kvLookup(info,"dx") != NULL && kvGetInt(info,"dx")!=0 ) {
    if (!kvLookup(info,"voxel_x"))
      kvDefDouble(info, "voxel_x", xscale/((double)kvGetInt(info,"dx") - 1.0));
  }
  if ( kvLookup(info,"dy") != NULL && kvGetInt(info,"dy")!=0 ) {
    if (!kvLookup(info,"voxel_y"))
      kvDefDouble(info, "voxel_y", yscale/((double)kvGetInt(info,"dy") - 1.0));
  }

  /* voxel_z was defined directly from the slice info in the header */
  /* nslices is dz */
  xplusbyVec3(tlf, blf, z_norm, 
	      (kvGetInt(info,"nslices")-1)*kvGetDouble(info,"voxel_z"));
  xplusbyVec3(trf, brf, z_norm, 
	      (kvGetInt(info,"nslices")-1)*kvGetDouble(info,"voxel_z"));
  xplusbyVec3(tlb, blb, z_norm, 
	      (kvGetInt(info,"nslices")-1)*kvGetDouble(info,"voxel_z"));
  xplusbyVec3(trb, brb, z_norm, 
	      (kvGetInt(info,"nslices")-1)*kvGetDouble(info,"voxel_z"));

  defVec3(info, "trf", trf);
  defVec3(info, "trb", trb);
  defVec3(info, "tlb", tlb);
  defVec3(info, "tlf", tlf);
  defVec3(info, "brf", brf);
  defVec3(info, "brb", brb);
  defVec3(info, "blb", blb);
  defVec3(info, "blf", blf);

}

static int scan_LX_data_header(KVHash* info, const char* readfile)
{
  KVHash* defs;
  KVHash* extNames;
  FILE *fphead;
  unsigned char header[FRZ_RDB_HEADER_SIZE_BYTES];
  unsigned char* rdbhead;
  unsigned char* acq_tab;
  unsigned char* examhead;
  unsigned char* serieshead;
  unsigned char* imagehead;
  char* here;
  int i;
  char buf[256];      /* scratch space */
  char pulse_seq[64]; /* pulse sequence ID string */
  float revnum;       /* header revision number */
  int ierror= 0;
  int exnum= 0;
  int num_frame_shifts= 0;
  double slice_tlc[3];
  double slice_trc[3];
  double slice_blc[3];
  double slice_brc[3];
  double slice_ctr[3];

  /* Definitions */
  defs= kvGetHash(info,"definitions");
  extNames= kvGetHash(info,"external_names");

  /* This bit is from the GE code io_signa_lx.c, with mods for portability */
  ierror= 0;
  if ((fphead = fopen(readfile,"r"))!=NULL)
    {
      if (fseek(fphead, (long) 0, SEEK_SET)) {
	perror("Error seeking header");
	ierror=1;
      }
      else {
	if (fread(header,sizeof(char),FRZ_RDB_HEADER_SIZE_BYTES, fphead)
	    != FRZ_RDB_HEADER_SIZE_BYTES) {
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

  rdbhead= header+FRZ_RDB_HDR_OFF;
  acq_tab= header+FRZ_RDB_DATAACQ_OFF;
  examhead= header+FRZ_RDB_EXAMDATATYPE_OFF;
  serieshead= header+FRZ_RDB_SERIESDATATYPE_OFF;
  imagehead= header+FRZ_RDB_MRIMAGEDATATYPE_OFF;

  /* We test to see if the data does have the expected endian order-
   revision number should be in the neighborhood of 7.0 */
  revnum= BRdFloat32(rdbhead+FRZ_RDBHEAD_RDB_HDR_RDBM_REV_OFF);
  if (revnum<1.0 || revnum>1000.0) {
    /* Oops, try it the other way! */
    bio_big_endian_input= (bio_big_endian_input ? 0 : 1);
    revnum= BRdFloat32(rdbhead+FRZ_RDBHEAD_RDB_HDR_RDBM_REV_OFF);
  }
  kvDefBoolean(info,"big_endian_input",bio_big_endian_input);
  if (debug) fprintf(stderr,"Header indicates %s input\n",
		     bio_big_endian_input ? "bigendian" : "littleendian");

  /* Sometimes there seems to be an extra, empty DATAACQ block. If so,
   * shift the reading frame accordingly.
   */
  ierror= 0;
  if ((fphead = fopen(readfile,"r"))!=NULL) {
    while ((exnum=BRdInt16(examhead+FRZ_EXAMHEAD_EX_NO_OFF)) == 0) {
      long dataacq_end= FRZ_RDB_DATAACQ_OFF + FRZ_RDB_DATAACQ_SIZE;
      num_frame_shifts++;
      Message("Shifting reading frame! (%d)\n",num_frame_shifts);
      if (fseek(fphead, (long)(dataacq_end
			       + num_frame_shifts*FRZ_RDB_DATAACQ_SIZE),
		SEEK_SET)) {
	perror("Error seeking header");
	ierror=1;
	break;
      }
	
      if (fread(header+dataacq_end, sizeof(char),
		FRZ_RDB_HEADER_SIZE_BYTES - dataacq_end,
		fphead)
	  != FRZ_RDB_HEADER_SIZE_BYTES-dataacq_end) {
	perror("Error reading header");
	ierror= 1;
	break;
      }

      if (num_frame_shifts >= MAX_FRAME_SHIFTS) {
	Error("Could not find a valid header after %d frame shifts!\n",
	      num_frame_shifts);
	ierror= 1;
	break;
      }
    }
    if (fclose(fphead)) {
      perror("Error closing header");
      ierror=1;
    }
  }
  else {
    perror("Error reopening header");
    ierror= 1;
  }
  if (ierror != 0) return 0;

  /*
   * Sometimes the GE header logo is INVALID rather than GE_MED_NMR.
   * Make a note of this and continue.
   *
   */
   if (!(strncmp((char*)rdbhead+FRZ_RDBHEAD_RDB_HDR_LOGO_OFF, 
		 FRZ_RDBHEAD_RDB_INVALID_LOGO, 
		 FRZ_RDBHEAD_RDB_HDR_LOGO_SIZE))) {
     kvDefBoolean(info,"GE_hdr_invalid",1);
     kvDefString(defs,"GE_hdr_invalid","GE RDB_LOGO is marked INVALID?");
   }
      
  /* 
   * these values seem true for all pulse sequences 
   */
  kvDefDouble(info,"table_delta",
	      BRdFloat32(imagehead+FRZ_IMAGEHEAD_TBLDLTA_OFF));
  kvDefString(defs,"table_delta","Table delta (mm)");

  kvDefString(extNames,"slice_thickness","slthick");
  kvDefDouble(info,"rdbm_rev", revnum);
  kvDefInt(info,"nslices",BRdInt16(rdbhead+FRZ_RDBHEAD_RDB_HDR_NSLICES_OFF));
  kvDefString(defs,"nslices","number of slices");
  kvDefString(extNames,"nslices",""); /* supress export */
  kvDefInt(info, "exnum", exnum); /* exam number */
  kvDefString(defs,"exnum","exam number");
  kvDefInt(info,"flip", 
	   BRdInt16(imagehead+FRZ_IMAGEHEAD_MR_FLIP_OFF)); /* flip angle */
  kvDefString(defs,"flip","flip angle");
  kvDefDouble(info,"slice_thickness",
	      BRdFloat32(imagehead+FRZ_IMAGEHEAD_SLTHICK_OFF));
  kvDefString(defs,"slice_thickness","slice thickness (mm)");
  kvDefString(extNames,"slice_thickness","slthick");
  kvDefDouble(info,"slice_gap",
	      BRdFloat32(imagehead+FRZ_IMAGEHEAD_SCANSPACING_OFF)); 
  kvDefString(defs,"slice_gap","slice gap (mm)");
  kvDefDouble(info,"voxel_z", 
	      kvGetDouble(info,"slice_thickness")
	      + kvGetDouble(info,"slice_gap"));
  kvDefString(defs,"voxel_z","Z voxel size including gap (mm)");
  kvDefDouble(info,"fov_z",
	      kvGetInt(info,"nslices")*kvGetDouble(info,"voxel_z")
	      - kvGetDouble(info,"slice_gap"));
  kvDefString(defs,"fov_z","Z field of view (mm)");
  kvDefInt(info,"TR",BRdInt32(imagehead+FRZ_IMAGEHEAD_TR_OFF));
  kvDefString(defs,"TR","TR (us)");
  kvDefInt(info,"TE",BRdInt32(imagehead+FRZ_IMAGEHEAD_TE_OFF));
  kvDefString(defs,"TE","TE (us)");
  strncpy(buf, (char*)(rdbhead+FRZ_RDBHEAD_RDB_HDR_SCAN_DATE_OFF), 10);
  buf[10]= '\0';
  kvDefString(info,"date",buf);
  kvDefString(defs,"date","scan date");
  strncpy(buf, (char*)(rdbhead+FRZ_RDBHEAD_RDB_HDR_SCAN_TIME_OFF), 8); 
  buf[8]= '\0';
  kvDefString(info,"time",buf);
  kvDefString(defs,"time","scan time");
#ifdef never
  /* Patient name info has generally been blanked out.  We don't want
   * to propogate it in any case.
   */
  strncpy(buf, (char*)(examhead+FRZ_EXAMHEAD_PATNAME_OFF),
	  FRZ_EXAMHEAD_PATNAME_SIZE);
  buf[FRZ_EXAMHEAD_PATNAME_SIZE]= '\0';
  fprintf(stderr,"******Patient Name*******<%s>\n",buf);
#endif

  /* Find the pulse sequence descriptor.  We've seen some which don't 
   * start at the beginning of the field, so we need to feel around a
   * bit to try to find the right region.
   */
  here= (char*)(imagehead + FRZ_IMAGEHEAD_PSDNAME_OFF);
  for (i=0; i<33; i++) {
    if (isprint(*here)) break;
    here++;
  }
  if (strrchr(here,'/')) 
    here= strrchr(here,'/')+1; /* skip leading file path */
  strncpy(pulse_seq, here, 33-i);
  pulse_seq[33-i] = '\0';
  if (debug) fprintf(stderr,"Pulse sequence name shifted by %d of %d chars\n",
		     i,33);
  kvDefString(info,"pulse_seq",pulse_seq);
  kvDefString(defs,"pulse_seq","pulse sequence");
  /* Determine Image Orientation - use imagehead.plane integer istead */
  /* of loc raster character */
  if(BRdInt16(imagehead+FRZ_IMAGEHEAD_PLANE_OFF) == 2)
    sprintf(buf,"AX");
  else if(BRdInt16(imagehead+FRZ_IMAGEHEAD_PLANE_OFF) == 4)
    sprintf(buf,"SAG");
  else if(BRdInt16(imagehead+FRZ_IMAGEHEAD_PLANE_OFF) == 8)
    sprintf(buf,"COR");
  else {
    Warning(1,"%s: Unable to detect image orientation\n",progname);
    Warning(1,"%s: Assuming AXIAL images; verify orientation!\n",progname);
    sprintf(buf,"AX");
  }
  kvDefString(info,"plane",buf);
  kvDefString(defs,"plane","scan plane");
  kvDefDouble(info,"autoscale_range",AUTOSCALE_RANGE);
  kvDefString(defs,"autoscale_range","range for autoscaling");

  kvDefInt(info,"index_rots", 
	   BRdInt16(rdbhead+FRZ_RDBHEAD_RDB_HDR_ROTATION_OFF));
  kvDefInt(info,"index_transpose",
	   BRdInt16(rdbhead+FRZ_RDBHEAD_RDB_HDR_TRANSPOSE_OFF));

  slice_tlc[0]= BRdFloat32(imagehead + FRZ_IMAGEHEAD_TLHC_R_OFF);
  slice_tlc[1]= BRdFloat32(imagehead + FRZ_IMAGEHEAD_TLHC_A_OFF);
  slice_tlc[2]= BRdFloat32(imagehead + FRZ_IMAGEHEAD_TLHC_S_OFF);
  slice_trc[0]= BRdFloat32(imagehead + FRZ_IMAGEHEAD_TRHC_R_OFF);
  slice_trc[1]= BRdFloat32(imagehead + FRZ_IMAGEHEAD_TRHC_A_OFF);
  slice_trc[2]= BRdFloat32(imagehead + FRZ_IMAGEHEAD_TRHC_S_OFF);
  slice_brc[0]= BRdFloat32(imagehead + FRZ_IMAGEHEAD_BRHC_R_OFF);
  slice_brc[1]= BRdFloat32(imagehead + FRZ_IMAGEHEAD_BRHC_A_OFF);
  slice_brc[2]= BRdFloat32(imagehead + FRZ_IMAGEHEAD_BRHC_S_OFF);
  slice_ctr[0]= BRdFloat32(imagehead + FRZ_IMAGEHEAD_CTR_R_OFF);
  slice_ctr[1]= BRdFloat32(imagehead + FRZ_IMAGEHEAD_CTR_A_OFF);
  slice_ctr[2]= BRdFloat32(imagehead + FRZ_IMAGEHEAD_CTR_S_OFF);

  /* The LX internal coordinate system is oriented differently from
   * Fiasco's.  LX has X running right-left in the images and Z
   * decreasing toward the feet; Fiasco has X left-right and Z
   * increasing toward the feet.
   */
  slice_tlc[0] *= -1.0;
  slice_trc[0] *= -1.0;
  slice_brc[0] *= -1.0;
  slice_ctr[0] *= -1.0;
  slice_tlc[2] *= -1.0;
  slice_trc[2] *= -1.0;
  slice_brc[2] *= -1.0;
  slice_ctr[2] *= -1.0;    

  defVec3(info,"slice_tlc",slice_tlc);
  defVec3(info,"slice_trc",slice_trc);
  defVec3(info,"slice_brc",slice_brc);
  defVec3(info,"slice_ctr",slice_ctr);

  kvDefString(defs,"slice_tlc.0","slice 0 top left corner R (RAS coords)");
  kvDefString(defs,"slice_tlc.1","slice 0 top left corner A (RAS coords)");
  kvDefString(defs,"slice_tlc.2","slice 0 top left corner S (RAS coords)");
  
  kvDefString(defs,"slice_trc.0","slice 0 top right corner R (RAS coords)");
  kvDefString(defs,"slice_trc.1","slice 0 top right corner A (RAS coords)");
  kvDefString(defs,"slice_trc.2","slice 0 top right corner S (RAS coords)");
  
  kvDefString(defs,"slice_brc.0","slice 0 bottom right corner R (RAS coords)");
  kvDefString(defs,"slice_brc.1","slice 0 bottom right corner A (RAS coords)");
  kvDefString(defs,"slice_brc.2","slice 0 bottom right corner S (RAS coords)");
  
  kvDefString(defs,"slice_blc.0","slice 0 bottom left corner R (RAS coords)");
  kvDefString(defs,"slice_blc.1","slice 0 bottom left corner A (RAS coords)");
  kvDefString(defs,"slice_blc.2","slice 0 bottom left corner S (RAS coords)");
  
  kvDefString(defs,"slice_ctr.0","slice 0 center R (RAS coords)");
  kvDefString(defs,"slice_ctr.1","slice 0 center A (RAS coords)");
  kvDefString(defs,"slice_ctr.2","slice 0 center S (RAS coords)");

  kvDefInt(info,"pt_size",
	   BRdInt16(rdbhead+FRZ_RDBHEAD_RDB_HDR_POINT_SIZE_OFF));
  kvDefString(defs,"pt_size","bytes per value");
  
  kvDefInt(info,"bviews",
	   BRdInt16(rdbhead+FRZ_RDBHEAD_RDB_HDR_BASELINE_VIEWS_OFF));
  kvDefString(defs,"bviews","# of baseline views");

  kvDefInt(info,"bl_save",0);
  kvDefString(defs,"bl_save","1 for Genesis, 0 for LX");

  /* The rest of the interpretation depends on the pulse sequence */
  if (!strncmp(pulse_seq,"epi",3)) {
    /* Interpret as EPIBOLD */
    ADD_VSUFFIX(scanEpiboldHeader)(info, rdbhead, acq_tab, examhead, 
				   serieshead, imagehead);
  }
  else if (!strncmp(pulse_seq,"splx",4)) {
    /* Variant of SPLX */
    ADD_VSUFFIX(scanSplxHeader)(info, rdbhead, acq_tab, examhead, 
				serieshead, imagehead);
  }
  else if (strstr(pulse_seq,"sf11") != NULL) {
    /* Obsolete multi-shot spiral */
    ADD_VSUFFIX(scanSf11Header)(info, rdbhead, acq_tab, examhead,
				serieshead, imagehead);
  }
  else if (strstr(pulse_seq,"sprl518") != NULL) {
    /* Gary Glover spiral? */
    Message("We think this is a Gary Glover-type spiral!\n");

    /* The format seems *mostly* like splx, except that some
     * fields seem invalid.  We'll scan it and then delete the
     * most-recently-added definitions of those fields.
     */
    ADD_VSUFFIX(scanSplxHeader)(info, rdbhead, acq_tab, examhead, 
				serieshead, imagehead);
    kvDelete(info,"slice_tlc.0"); kvDelete(info,"slice_tlc.1"); 
    kvDelete(info,"slice_tlc.2");
    kvDelete(info,"slice_trc.0"); kvDelete(info,"slice_trc.1"); 
    kvDelete(info,"slice_trc.2");
    kvDelete(info,"slice_brc.0"); kvDelete(info,"slice_brc.1"); 
    kvDelete(info,"slice_brc.2");
    kvDelete(info,"slice_ctr.0"); 
    kvDelete(info,"slice_ctr.1"); kvDelete(info,"slice_ctr.2");

    /* There is also a data offset correction based on odd-ness of dt */
    if (kvGetInt(info,"dt")%2) {
      /* skip both dummy frame and next baseline frame */
      kvDefLong(info,"skip.t",
		2*(2 * kvGetInt(info,"pt_size") * kvGetInt(info,"ndatfr")));
    }
    else {
      /* skip only next baseline frame */
      kvDefLong(info,"skip.t",
		(2 * kvGetInt(info,"pt_size") * kvGetInt(info,"ndatfr")));
    }
  }
  else if (!strncmp(pulse_seq,"2dfast",3)) {
    /* Interpret as 2dfast structural scan */
    ADD_VSUFFIX(scan2dfastHeader)(info, rdbhead, acq_tab, examhead, 
				  serieshead, imagehead);
  }
  else {
    Message("Unknown pulse sequence name <%s>!\n",pulse_seq);
    Message("I'll pretend it's spiral, because they are the most inconsistent!\n");
    ADD_VSUFFIX(scanSplxHeader)(info, rdbhead, acq_tab, examhead, 
				serieshead, imagehead);
  }

  /* Now that we know as much about the scan as possible, calculate
   * the bounds of the scan volume.
   */
  calcVolumeBounds(info);

  /* If we did any frame shifting to find the header info, adjust the
   * starting offset appropriately.
   */
  kvDefLong(info,"start_offset", 
	    kvGetLong(info,"start_offset") 
	    + num_frame_shifts*FRZ_RDB_DATAACQ_SIZE);
  return 1;
}

static void addExtraChunks( FileHandler* self, KVHash* info, SList* cStack )
{
  KVHash* defs= kvGetHash(info,"definitions");

  if (!strncmp(kvGetString(info,"pulse_seq"),"splx",4)
      || strstr(kvGetString(info,"pulse_seq"),"sprl518")) 
    ADD_VSUFFIX(addSplxTrajChunk)(self, info, cStack);

  if (strstr(kvGetString(info,"pulse_seq"),"sf11")) 
    ADD_VSUFFIX(addSf11TrajChunk)(self, info, cStack);  

  if (kvLookup(info,"bandpassdir")) 
    ADD_VSUFFIX(addBandpassChunk)(self, info, cStack);

  if (kvLookup(info,"phasereffile")) 
    ADD_VSUFFIX(addPhaseRefChunk)(self, info, cStack);

  if (kvLookup(info,"rampfile"))
    ADD_VSUFFIX(addRampSampleChunk)(self, info, cStack);
}

static void processHeader( FileHandler* self, KVHash* info, SList* cStack )
{
  /* Call the base class method */
  baseProcessHeader( self, info, cStack );

  if (!scan_LX_data_header(info, self->fileName))
    Abort("lx_reader: unable to read or parse header from <%s>!\n",
	  self->fileName);

  addExtraChunks( self, info, cStack );
}

#ifdef never
static void fakeRead( FileHandler* self, KVHash* info,
		      long long offset, long n,
		      SRDR_Datatype datatype_out,
		      void* obuf ) {
  static int count= 0;
  int count2= 0;
  short* s= (short*)obuf;
  int i;
  int thisRun= 0;
  int maxRun= 0;
  int maxRunStart;
  int thisRunStart;
  baseRead( self, info, offset, n, datatype_out, obuf );
  for (i=0; i<n; i++) 
    if (s[i]==0) {
      count2++;
      thisRun++;
    }
    else {
      if (thisRun>maxRun) {
	maxRun= thisRun;
	maxRunStart= thisRunStart;
      }
      thisRun= 0;
      thisRunStart= i;
    }
  if (thisRun>maxRun) maxRun= thisRun;
  fprintf(stderr,"Just read %ld from %lld; count= %d; %d zeros, max run %d from %d\n",
	  n,offset,count++,count2,maxRun,maxRunStart);
}
#endif

FileHandler* ADD_VSUFFIX(lxFactory)(char* fname, KVHash* info)
{
  FileHandler* result= baseFactory(fname);
  result->processHeader= processHeader;
  result->typeName= strdup( nameof_filetype(FILE_LX) );
#ifdef never
  result->read= fakeRead;
#endif
  return result;
}

int ADD_VSUFFIX(lxTester)(const char* filename)
{
  /* NOTE- this implicitly  makes the assumption that the location and 
   *       value of the GE LX rdb header logo don't change between 
   *       revisions of rdbm.h .  So far, that seems to be true.
   */
  if (check_filetype(filename)!=FILE_LX) return 0;
  else {

    /* At this point we know it's some type of GE LX, but we don't know
     * if this version is compiled for the correct rdbm.h revision.  
     * Check the revision number in the header against the compiled-in
     * version.
     */
    unsigned char buf[sizeof(float)];
    int ierror= 0;
    FILE* fphead= NULL;
    float revnum;

    /* This bit is from the GE code io_signa_lx.c, with mods for portability */
    if ((fphead = fopen(filename,"r"))!=NULL)
      {
	if (fseek(fphead, FRZ_RDBHEAD_RDB_HDR_RDBM_REV_OFF, SEEK_SET)) {
	  perror("Error seeking header");
	  ierror=1;
	}
	else {
	  if (fread(buf, sizeof(float), 1, fphead) != 1) {
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

    /* We test to see if the data does have the expected endian order-
       revision number should be reasonable */
    revnum= BRdFloat32(buf);
    if (revnum<1.0 || revnum>1000.0) {
      /* Oops, try it the other way! */
      bio_big_endian_input= (bio_big_endian_input ? 0 : 1);
      revnum= BRdFloat32(buf);
    }

    /* And finally, we test against the compiled-in revision number */
    if (debug) 
      fprintf(stderr,"Checking rdbm revision %f against %f\n",
	      revnum, FRZ_RDBHEAD_RDB_HDR_RDBM_REVISION);
    if (revnum < FRZ_RDBHEAD_RDB_HDR_RDBM_REVISION - REVISION_EPSILON
	|| revnum >= FRZ_RDBHEAD_RDB_HDR_RDBM_REVISION+1.0) {
      return 0;
    }
  }

  return 1;
}

