/************************************************************
 *                                                          *
 *  siemens_kspace_reader.c                                 *
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
 ************************************************************/
/* This module is based heavily on a Matlab routine provided by
 * Lawrence L. Wald (wald@nmr.mgh.harvard.edu).
 */
/* Notes-
   -Offsets for timeSinceLastRF, readoutOffCtr are tenuous at best.
   -Multiple channel input is not implemented.  This will require a
    new dimension in output and changes to Line structure and routines.
   -There is slice orientation info in the line headers, but I don't
    know how to make inferences from it.
   -I never check to make sure all the mdh's have the same line length,
    k-space center, etc.
   -lx_reader seems to think that the distance between top and bottom is
    (dz-1)*voxel_z.  
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
#include "siemens_kspace_header_info.h"

static char rcsid[] = "$Id: siemens_kspace_reader.c,v 1.32 2006/06/22 06:33:21 welling Exp $";

#define SIEMENS_RESAMPLE_SCRIPT "siemens_ramp_resample.csh"

typedef struct mdh_struct {
  long long dmaLength;
  long measUID;
  long long scanCounter;
  long long timeStamp;
  long long pmuTimeStamp;
  long evalInfoMask;
  long samplesInScan;
  long usedChannels;
  long lineNum;
  long acqNum;
  long sliceNum;
  long partitionNum;
  long echoNum;
  long phaseNum;
  long repNum;
  long setNum;
  long segNum;
  long freeNum;
  long cutoffPre;
  long cutoffPost;
  long kspaceCtrCol;
  float readoutOffCtr;
  long long timeSinceLastRF;
  long kspaceCtrLine;
  long kspaceCtrPartition;
  long freeParams[SMNSKSPC_NUM_FREE_PARAMS];
  float slicePos[3]; /* saggital, coronal, transverse; DICOM coords */
  Quat quaternion;
  long long channelId;
} MeasurementDataHeader;

typedef struct line_struct {
  long long offset; /* real offset of data */
  long n; /* number of data elements */
  long line;
  long acq;
  long slice;
  long partition;
  long rep;
  long set;
  long seg;
  long nav_id;
  long mask;
  float slicePos[3]; /* DICOM coords */
  float sortingDist; /* used only when sorting slices */
} Line;

/* This function implements the coordinate transformation from
 * Siemens internal coords to Fiasco coords.
 */
static void transformSiemensToFiascoVec3(double* v)
{
  /* Slice normal and tlc require special treatment because Fiasco
   * coordinates differ from Siemens coordinates:
   *
   * Siemens-to-Fiasco transformation matrix is:
   *  /  1  0  0 \
   *  |  0  1  0 |
   *  \  0  0  1 /
   */
#ifdef never
  v[0] *= -1.0;
#endif
  v[1] *= -1.0;
  v[2] *= -1.0;
}

static void calcVolumeBounds(KVHash* info)
{
  double slice_ctr[3];
  double blf[3];
  double brf[3];
  double blb[3];
  double brb[3];
  double tlf[3];
  double trf[3];
  double tlb[3];
  double trb[3];
  double delta[3];
  double x_norm[3];
  double y_norm[3];
  double z_norm[3];
  double tmpVec[3];
  double fov_x;
  double fov_y;
  double fov_z;

  /* This routine will probably have to be updated as we learn more 
   * about their notion of coordinate systems.
   */

  if (testVec3(info,"slice_ctr") && testVec3(info,"slice_norm")
      && testVec3(info,"freqEncodeDir")
      && kvLookup(info,"fov_x") && kvLookup(info,"fov_y")
      && kvLookup(info,"fov_z")) {
    fov_x= kvGetDouble(info,"fov_x");
    fov_y= kvGetDouble(info,"fov_y");
    fov_z= kvGetDouble(info,"fov_z");
    getVec3(info,"slice_ctr",slice_ctr);
    getVec3(info,"slice_norm",z_norm);
    getVec3(info,"freqEncodeDir",x_norm);
    crossVec3(y_norm,z_norm,x_norm);

    Report("Norms: (%f %f %f) (%f %f %f) (%f %f %f)\n",
	   x_norm[0], x_norm[1], x_norm[2],
	   y_norm[0], y_norm[1], y_norm[2],
	   z_norm[0], z_norm[1], z_norm[2]);
    
    xplusbyVec3(tmpVec,slice_ctr,x_norm,-0.5*fov_x);
    xplusbyVec3(blb,tmpVec,y_norm,0.5*fov_y);
    xplusbyVec3(blf,blb,y_norm,-fov_y);
    xplusbyVec3(brb,blb,x_norm,fov_x);
    xplusbyVec3(brf,brb,y_norm,-fov_y);
    multVec3(delta,z_norm,fov_z);
    xplusbyVec3(tlf,blf,delta,1.0);
    xplusbyVec3(trf,brf,delta,1.0);
    xplusbyVec3(tlb,blb,delta,1.0);
    xplusbyVec3(trb,brb,delta,1.0);

    defVec3(info, "tlf", tlf);
    defVec3(info, "trf", trf);
    defVec3(info, "tlb", tlb);
    defVec3(info, "trb", trb);
    defVec3(info, "brf", brf);
    defVec3(info, "brb", brb);
    defVec3(info, "blb", blb);
    defVec3(info, "blf", blf);
  }
  else {
    if (debug) fprintf(stderr,"Not enough info to calculate volume bounds!\n");
  }

}

static long fixShortSign( long i )
{
  /* This stupid standard stores everything as unsigned shorts, which
   * are not supported on most architectures.  Any negative value
   * that comes in here is actually a 16 bit unsigned positive value.
   */
  return (i>=0) ? i : ~i+1;
}

static long long fixLongSign( long long i )
{
  /* For want of a bit, the longword was lost.  For loss of the
   * longword, the offset was lost.  For loss of the offset, the
   * parse was lost.
   */
  return (i>=0) ? i : ~i+1;
}

static int parse_mdh( FILE* file, long long offset, 
		       MeasurementDataHeader* mdh )
{
  unsigned char buf[SMNSKSPC_LINE_HDR_BYTES];
  int i;
  int bytesRead;
  

  if (bigfile_fseek(file, offset, SEEK_SET))
    Abort("%s: unable to seek to offset %lld in input!\n",offset);

  if ((bytesRead=fread(buf, sizeof(char), SMNSKSPC_LINE_HDR_BYTES, file))
      != SMNSKSPC_LINE_HDR_BYTES) {
    if (bytesRead==0) return 0;
    else Abort("%s: read short line (length %d) from %lld: %s\n",
	       progname, bytesRead, offset, strerror(errno));
  }

  mdh->dmaLength= fixLongSign(BRdInt32(buf+SMNSKSPC_DMALENGTH_OFF));
  mdh->measUID= BRdInt32(buf+SMNSKSPC_MEASUID_OFF);
  mdh->scanCounter= fixLongSign(BRdInt32(buf+SMNSKSPC_SCANCOUNTER_OFF));
  mdh->timeStamp= fixLongSign(BRdInt32(buf+SMNSKSPC_TIMESTAMP_OFF));
  mdh->pmuTimeStamp= fixLongSign(BRdInt32(buf+SMNSKSPC_PMUTIMESTAMP_OFF));
  mdh->evalInfoMask= BRdInt32(buf+SMNSKSPC_EVALINFOMASK_OFF);

  mdh->samplesInScan= fixShortSign(BRdInt16(buf+SMNSKSPC_SAMPLESINSCAN_OFF));
  mdh->usedChannels= fixShortSign(BRdInt16(buf+SMNSKSPC_USEDCHANNELS_OFF));
  mdh->lineNum= fixShortSign(BRdInt16(buf+SMNSKSPC_LINENUM_OFF));
  mdh->acqNum= fixShortSign(BRdInt16(buf+SMNSKSPC_ACQNUM_OFF));
  mdh->sliceNum= fixShortSign(BRdInt16(buf+SMNSKSPC_SLICENUM_OFF));
  mdh->partitionNum= fixShortSign(BRdInt16(buf+SMNSKSPC_PARTITIONNUM_OFF));
  mdh->echoNum= fixShortSign(BRdInt16(buf+SMNSKSPC_ECHONUM_OFF));
  mdh->phaseNum= fixShortSign(BRdInt16(buf+SMNSKSPC_PHASENUM_OFF));
  mdh->repNum= fixShortSign(BRdInt16(buf+SMNSKSPC_REPNUM_OFF));
  mdh->segNum= fixShortSign(BRdInt16(buf+SMNSKSPC_SEGNUM_OFF));
  mdh->setNum= fixShortSign(BRdInt16(buf+SMNSKSPC_SETNUM_OFF));
  mdh->freeNum= fixShortSign(BRdInt16(buf+SMNSKSPC_FREENUM_OFF));

  mdh->cutoffPre= fixShortSign(BRdInt16(buf+SMNSKSPC_CUTOFF_PRE_OFF));
  mdh->cutoffPost= fixShortSign(BRdInt16(buf+SMNSKSPC_CUTOFF_POST_OFF));

  mdh->kspaceCtrCol= fixShortSign(BRdInt16(buf+SMNSKSPC_KSPACECTRCOL_OFF));

  mdh->readoutOffCtr= BRdFloat32(buf+SMNSKSPC_READOUTOFFCTR_OFF);
				
  mdh->timeSinceLastRF= 
    fixLongSign(BRdInt32(buf+SMNSKSPC_TIMESINCELASTRF_OFF));

  mdh->kspaceCtrLine= 
    fixShortSign(BRdInt16(buf+SMNSKSPC_KSPACECTRLINE_OFF));
  mdh->kspaceCtrPartition= 
    fixShortSign(BRdInt16(buf+SMNSKSPC_KSPACECTRPARTITION_OFF));

  for (i=0; i<SMNSKSPC_NUM_FREE_PARAMS; i++)
    mdh->freeParams[i]= 
      fixShortSign(BRdInt16(buf+SMNSKSPC_FREEPARAMS_OFF + 2*i));

  mdh->slicePos[0]= BRdFloat32(buf+SMNSKSPC_SLICEPOS_SAG_OFF);
  mdh->slicePos[1]= BRdFloat32(buf+SMNSKSPC_SLICEPOS_COR_OFF);
  mdh->slicePos[2]= BRdFloat32(buf+SMNSKSPC_SLICEPOS_TRA_OFF);

  mdh->quaternion.x= BRdFloat32(buf+SMNSKSPC_QUATERNION_OFF);
  mdh->quaternion.y= BRdFloat32(buf+SMNSKSPC_QUATERNION_OFF+4);
  mdh->quaternion.z= BRdFloat32(buf+SMNSKSPC_QUATERNION_OFF+8);
  mdh->quaternion.w= BRdFloat32(buf+SMNSKSPC_QUATERNION_OFF+12);

  mdh->channelId= fixLongSign(BRdInt32(buf+SMNSKSPC_CHANNELID_OFF));
  return 1;
}

static void dump_mdh( FILE* ofile, MeasurementDataHeader* mdh )
{
  int i;

  fprintf(ofile,
	  "dmaLength %lld, scanCounter %lld, timeStamp %lld, pmuTimeStamp %lld\n",
	  mdh->dmaLength, mdh->scanCounter, mdh->timeStamp, mdh->pmuTimeStamp);
  fprintf(ofile,"ACQEND %d, PHASECOR %d, FIRSTSCANINSLICE %d, LASTSCANINSLICE %d, REFLECT %d\n",
	  SMNSKSPC_EVALINFO(ACQEND,mdh->evalInfoMask),
	  SMNSKSPC_EVALINFO(PHASECOR,mdh->evalInfoMask),
	  SMNSKSPC_EVALINFO(FIRSTSCANINSLICE,mdh->evalInfoMask),
	  SMNSKSPC_EVALINFO(LASTSCANINSLICE,mdh->evalInfoMask),
	  SMNSKSPC_EVALINFO(REFLECT,mdh->evalInfoMask));
  fprintf(ofile,
	  "measUID %ld, evalInfoMask 0x%08x, samplesInScan %ld, usedChannels %ld\n",
	  mdh->measUID, mdh->evalInfoMask, mdh->samplesInScan, 
	  mdh->usedChannels);
  fprintf(ofile,
	  "lineNum %ld, acqNum %ld, sliceNum %ld, partitionNum %ld\n",
	  mdh->lineNum, mdh->acqNum, mdh->sliceNum, mdh->partitionNum);
  fprintf(ofile,
	  "echoNum %ld, phaseNum %ld, repNum %ld, setNum %ld\n",
	  mdh->echoNum,mdh->phaseNum,mdh->repNum,mdh->setNum);
  fprintf(ofile,
	  "segNum %ld, freeNum %ld, cutoffPre %ld, cutoffPost %ld\n",
	  mdh->segNum, mdh->freeNum, mdh->cutoffPre, mdh->cutoffPost);
  fprintf(ofile,
	  "kspaceCtrCol %ld, kspaceCtrLine %ld, kspaceCtrPartition %ld\n",
	  mdh->kspaceCtrCol, mdh->kspaceCtrLine, mdh->kspaceCtrPartition);
  fprintf(ofile,"timeSinceLastRF %lld\n",
	  mdh->timeSinceLastRF);
  fprintf(ofile,
	  "slicePos (%f, %f, %f), readoutOffCtr %f\n",
	  mdh->slicePos[0],mdh->slicePos[1],mdh->slicePos[2], 
	  mdh->readoutOffCtr);
  fprintf(ofile,"Quaternion (%f, %f, %f, %f)\n",
	  mdh->quaternion.x,mdh->quaternion.y,mdh->quaternion.z,
	  mdh->quaternion.w);
  fprintf(ofile,"freeParams: ");
  for (i=0; i<SMNSKSPC_NUM_FREE_PARAMS; i++) 
    fprintf(ofile,"%d ",mdh->freeParams[i]);
  fprintf(ofile,"\n");
  fprintf(ofile,"timeSinceLastRF %lld, channelId %lld\n",
	  mdh->timeSinceLastRF, mdh->channelId);
}

static int reality_check_mdh( MeasurementDataHeader* mdh )
{
  /* dmaLength should equal mdh length + samplesInScan*sizeof(sample) */
  return ( mdh->dmaLength%(SMNSKSPC_LINE_HDR_BYTES + 8*mdh->samplesInScan)
	   == 0 );
}

static int reality_check_file( FILE* f )
{
  MeasurementDataHeader mdh;
  long long offset= SMNSKSPC_HEADER_SIZE_BYTES;

  while (!bigfile_fseek(f, offset, SEEK_SET)) {
    if (!parse_mdh( f, offset, &mdh )) return 1; /* end of file */
    if (!reality_check_mdh( &mdh )) return 0; /* not expected format */
    offset += mdh.dmaLength;
  }
  return 1;
}

static void line_dump( FILE* ofile, Line* thisLine )
{
  fprintf(ofile,"%8lld: ln %d, ac %d, sl %d, prt %d, rep %d, set %d, seg %d, mask %d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d -> 0x%08x\n",
	  thisLine->offset, thisLine->line, thisLine->acq,
	  thisLine->slice, thisLine->partition, thisLine->rep,
	  thisLine->set, thisLine->seg,
	  (SMNSKSPC_EVALINFO(ACQEND,thisLine->mask)!=0),
	  (SMNSKSPC_EVALINFO(RTFEEDBACK,thisLine->mask)!=0),
	  (SMNSKSPC_EVALINFO(HPFEEDBACK,thisLine->mask)!=0),
	  (SMNSKSPC_EVALINFO(ONLINE,thisLine->mask)!=0),
	  (SMNSKSPC_EVALINFO(OFFLINE,thisLine->mask)!=0),
	  (SMNSKSPC_EVALINFO(REFPHASESTABSCAN,thisLine->mask)!=0),
	  (SMNSKSPC_EVALINFO(PHASESTABSCAN,thisLine->mask)!=0),
	  (SMNSKSPC_EVALINFO(D3FFT,thisLine->mask)!=0),
	  (SMNSKSPC_EVALINFO(SIGNREV,thisLine->mask)!=0),
	  (SMNSKSPC_EVALINFO(PHASEFFT,thisLine->mask)!=0),
	  (SMNSKSPC_EVALINFO(SWAPPED,thisLine->mask)!=0),
	  (SMNSKSPC_EVALINFO(POSTSHAREDLINE,thisLine->mask)!=0),
	  (SMNSKSPC_EVALINFO(PHASECOR,thisLine->mask)!=0),
	  (SMNSKSPC_EVALINFO(ZEROLINE,thisLine->mask)!=0),
	  (SMNSKSPC_EVALINFO(ZEROPARTITION,thisLine->mask)!=0),
	  (SMNSKSPC_EVALINFO(REFLECT,thisLine->mask)!=0),
	  (SMNSKSPC_EVALINFO(NOISEADJSCAN,thisLine->mask)!=0),
	  (SMNSKSPC_EVALINFO(SHARENOW,thisLine->mask)!=0),
	  (SMNSKSPC_EVALINFO(LASTMEASUREDLINE,thisLine->mask)!=0),
	  (SMNSKSPC_EVALINFO(FIRSTSCANINSLICE,thisLine->mask)!=0),
	  (SMNSKSPC_EVALINFO(LASTSCANINSLICE,thisLine->mask)!=0),
	  (SMNSKSPC_EVALINFO(TREFFECTIVEBEGIN,thisLine->mask)!=0),
	  (SMNSKSPC_EVALINFO(TREFFECTIVEEND,thisLine->mask)!=0),
	  thisLine->mask
	  );
}

static Line* line_create( MeasurementDataHeader* mdh, long long offset )
{
  int i;
  Line* result;
  if (!(result= (Line*)malloc(sizeof(Line))))
    Abort("%s: unable to allocate %d bytes!\n",progname,sizeof(Line));

  result->offset= offset+SMNSKSPC_LINE_HDR_BYTES;
  result->n= 2*mdh->samplesInScan;
  result->line= mdh->lineNum;
  result->acq= mdh->acqNum;
  result->slice= mdh->sliceNum;
  result->partition= mdh->partitionNum;
  result->rep= mdh->repNum;
  result->set= mdh->setNum;
  result->seg= mdh->segNum;
  result->mask= mdh->evalInfoMask;
  for (i=0; i<3; i++) result->slicePos[i]= mdh->slicePos[i];
  result->sortingDist= 0.0;

  /* This helps in keeping navigator pulses in the right order */
  if (SMNSKSPC_EVALINFO(PHASECOR,mdh->evalInfoMask))
    result->nav_id= 2*result->acq + (1-result->seg);
  else {
    result->nav_id= 0;
  }

  return result;
}

static void line_destroy( void* l )
{
  free(l);
}

static long long line_count_range( Line* l1, Line* l2 )
{
  /* If l1 and l2 are the first and last of a complete set of
   * lines, how many lines are there in all?  We count only those
   * fields expected to make a difference to the Fiasco output
   * dimension string.
   */
  long long result= 1;

  result *= (l2->line - l1->line) + 1;    /* dimension y */
  result *= (l2->slice - l1->slice) + 1;  /* dimension z */
  result *= (l2->partition - l1->partition) + 1; /* equiv to slice 
						    for structurals */
  result *= (l2->rep - l1->rep) + 1;      /* dimension t */
  result *= (l2->nav_id - l1->nav_id) + 1; /* navigator dimension */

  return result;
}

static int line_compare( const void* p1, const void* p2 )
{
  Line* l1= *(Line**)p1;
  Line* l2= *(Line**)p2;

#define COMPARE_THIS(a,b) { if (a < b) return -1; else if (a > b) return 1; }

  /* repetition corresponds to typical Fiasco t */
  COMPARE_THIS(l1->rep, l2->rep);

  /* I'm guessing that partition changes more slowly than slice.
   * It corresponds to the slice number for structural scans. 
   */
  COMPARE_THIS(l1->partition, l2->partition);

  /* Slice corresponds to z (surprisingly enough) */
  COMPARE_THIS(l1->slice, l2->slice);

  /* line number corresponds to y */
  COMPARE_THIS(l1->line, l2->line);

  /* The following keeps the navigator scans in the right order */
  COMPARE_THIS(l1->nav_id, l2->nav_id);

  /* I think there might sometimes be more than one acquisition in
   * a given rep, partition, slice, line.
   */
  COMPARE_THIS(l1->acq, l2->acq);

  return 0;

#undef COMPARE_THIS
}

static int sliceSortLineCompare( const void* p1, const void* p2 )
{
  Line** linep1= (Line**)p1;
  Line** linep2= (Line**)p2;
#ifdef never
  fprintf(stderr,"compare %d %f vs %d %f\n",
	  (*linep1)->slice, (*linep1)->sortingDist,
	  (*linep2)->slice, (*linep2)->sortingDist);
#endif
  if ((*linep1)->sortingDist<(*linep2)->sortingDist) return -1;
  else if ((*linep1)->sortingDist>(*linep2)->sortingDist) return 1;
  else return 0;
}

static void relabel_lines( KVHash* info, SList* lineList )
{
  static double zDir[]= {0.0,0.0,1.0};
  static double sliceNorm[3];
  int dz= kvGetInt(info,"dz"); /* number of slices */
  int dy= 1; /* number of lines per slice; corrected below */
  int reverseY= 0;
  double sliceCtr[3];
  SList* sliceSortList= NULL;
  int* sliceOrderTable= NULL;
  int prevSlice;
  Line* line;

  /* If we have info about the scan orientation, we can improve the
   * guesses for what needs to be reversed.
   */
  if (kvLookup(info,"dy")) dy= kvGetInt(info,"dy");
  if (testVec3(info,"sliceAcqDir")) getVec3(info,"sliceAcqDir",sliceNorm);
  else copyVec3(sliceNorm, zDir); /* Best guess! */
  flipToPositiveHemisphereVec3(sliceNorm);
  defVec3(info,"slice_norm",sliceNorm);

  /* Infer slice reorder pattern by sorting the locations of the
   * slices along the slice normal direction.
   */
  slist_totop(lineList);
  sliceSortList= slist_create();
  prevSlice= -1;
  while (!slist_atend(lineList)) {
    Line* thisLine= (Line*)slist_get(lineList);
    if (debug) line_dump(stderr,thisLine);
    if (thisLine->rep != 0 || thisLine->partition != 0) break;
    if (thisLine->slice != prevSlice && thisLine->acq==0) {
      int i;
      double pos[3];
      for (i=0; i<3; i++) pos[i]= thisLine->slicePos[i];
      transformSiemensToFiascoVec3(pos);
      thisLine->sortingDist= dotVec3(sliceNorm,pos);
      slist_append(sliceSortList,thisLine);
      prevSlice= thisLine->slice;
    }
    slist_next(lineList);
  }
  if (slist_count(sliceSortList) != 0) {
    const char* name= NULL;
    int i;

    if (slist_count(sliceSortList)!=dz)
      Abort("%s: scan line list has inconsistent structure (%d vs. %d)!\n",
	    progname,slist_count(sliceSortList),dz);

    if (!(sliceOrderTable= (int*)malloc(dz*sizeof(int))))
      Abort("%s: unable to allocate %d bytes!\n",progname, dz*sizeof(int));

    slist_sort(sliceSortList,sliceSortLineCompare);

    slist_totop(sliceSortList);
    i= 0;
    while (!slist_atend(sliceSortList)) {
      Line* thisLine= (Line*)slist_get(sliceSortList);
      sliceOrderTable[thisLine->slice]= i;
      slist_next(sliceSortList);
      i++;
    }
    /* Note the slice order pattern and the location of the bottom slice.
     */
    if ((name= slp_findSlicePatternNameFromTable(dz,sliceOrderTable)) != NULL)
      kvDefString(info,"reorder_pattern",name);
    slist_totop(sliceSortList);
    line= (Line*)slist_get(sliceSortList);
#ifdef never
    line= (Line*)slist_getlast(sliceSortList);
#endif
    for (i=0; i<3; i++) sliceCtr[i]= line->slicePos[i];
    transformSiemensToFiascoVec3(sliceCtr);
    defVec3(info,"slice_ctr",sliceCtr);
  }
  slist_destroy(sliceSortList,NULL);

  if (testVec3(info,"phaseEncodeDir")) {
    double phaseEncodeDir[3];
    getVec3(info,"phaseEncodeDir",phaseEncodeDir);
    if (fabs(dotVec3(phaseEncodeDir,zDir))>=1.0/sqrt(2.0))
      reverseY= 1;
    else reverseY= 0;
  }

  /* We have no way to reverse X here; could change rowflip order later */
  slist_totop(lineList);
  while (!slist_atend(lineList)) { 
    Line* thisLine= (Line*)slist_get(lineList);
    int z= thisLine->slice;
    int y= thisLine->line;
    slist_next(lineList);
    thisLine->slice= sliceOrderTable[z];
    if (reverseY) {
      /* This never happens if dy is unknown (defaulting to 1) */
      thisLine->line= dy-(y+1);
    }
    else {
      thisLine->line= y;
    }
  }
}

static int reality_check_linelist( SList* l )
{
  Line* firstLine;
  Line* lastLine;
  long long expectedLines;
  long long gotLines;

  if (slist_empty(l)) return 1; /* list is empty, which is OK */
  slist_totop(l);
  firstLine= slist_get(l);
  lastLine= slist_getlast(l);
  expectedLines= line_count_range(firstLine,lastLine);
  gotLines= slist_count(l);
  if (debug) 
    fprintf(stderr,"reality_check_linelist: expected %lld, got %lld\n",
	    expectedLines,gotLines);
  return (expectedLines==gotLines);
}

static void destroySelf( FileHandler* self )
{
  /* We need to free any line list we happen to have */
  SList* lineList= (SList*)self->hook;
  if (lineList) {
    slist_destroy(lineList,line_destroy);
  }
  baseDestroySelf(self);
}

static void walk_mdh_structures( KVHash* info, FileHandler* self, FILE* f,
				 SList** lineList, SList** navList )
{
  MeasurementDataHeader mdh;
  long long offset= SMNSKSPC_HEADER_SIZE_BYTES;
  int i= 0;

  slist_totop(*lineList);
  slist_totop(*navList);

  while (offset < self->totalLengthBytes) {
    Line* thisLine;
#ifdef never
    fprintf(stderr,"%%%%%%%%%%%%%%%%%% %d at offset %lld %%%%%%%%%%%\n",
	    i++, offset);
#endif
    if (!parse_mdh( f, offset, &mdh )) 
      Abort("%s: unexpected end of file on %s!\n",
	    progname, self->fileName);
#ifdef never
    dump_mdh( stderr, &mdh );
#endif

    if (!(SMNSKSPC_EVALINFO(ACQEND,mdh.evalInfoMask))) { /* ignore end tag */
      if (SMNSKSPC_EVALINFO(PHASECOR,mdh.evalInfoMask))
	slist_append( *navList, line_create(&mdh,offset) );
      else slist_append( *lineList, line_create(&mdh,offset) );
    }

    offset += SMNSKSPC_LINE_HDR_BYTES + 8*mdh.samplesInScan;
  }

  slist_sort(*navList, line_compare);
  slist_sort(*lineList, line_compare);

  if (!reality_check_linelist(*navList))
    Abort("%s: some navigator lines seem to be missing or repeated! (or not EPI)\n",
	  progname);
  if (!reality_check_linelist(*lineList))
    Abort("%s: some phase encode lines seem to be missing or repeated! (or not EPI)\n",
	  progname);
}

static char* pick_rowflip_pattern( SList* lineList, long dy, 
				   const char* chunkName )
{
  Line* thisLine;
  long y;
  int revEven= 0;
  int revOdd= 0;

  /* As nearly as I can tell, when the REFLECT flag is set in
   * the info mask the scan actually goes from low to high X
   * in k-space, which to Fiasco means that the data is *not*
   * reflected.  Thus we need to take the logical 'not' of
   * the flag.
   */

  slist_totop(lineList);
  for (y=0; y<dy; y++) {
    thisLine= (Line*)slist_get(lineList);
    slist_next(lineList);
    if (y==0) revEven= (!SMNSKSPC_EVALINFO(REFLECT,thisLine->mask));
    else if (y==1) revOdd= (!SMNSKSPC_EVALINFO(REFLECT,thisLine->mask));
    else {
      if (y%2) {
	if (revOdd != (!SMNSKSPC_EVALINFO(REFLECT,thisLine->mask)))
	  Abort("%s: unexpected line reversal pattern for chunk %s!\n",
		progname, chunkName);
      }
      else {
	if (revEven != (!SMNSKSPC_EVALINFO(REFLECT,thisLine->mask)))
	  Abort("%s: unexpected line reversal pattern for chunk %s!\n",
		progname, chunkName);
      }
    }
  }
  if (revEven) {
    if (revOdd) return( strdup("both") );
    else return( strdup("even") );
  }
  else {
    if (revOdd) return( strdup("odd") );
    else return( strdup("none") );
  }
}

static void infer_scan_type_from_list( KVHash* info,
				       Line* firstLine, Line* lastLine,
				       MeasurementDataHeader* mdh )
{
  KVHash* defs= kvGetHash(info,"definitions");
  /* We have the first and last entries of a scan, and we know
   * it is not a navigator.  We have also established (when the
   * line list was reality-checked) that it contains a full scan.
   */
  static double xDir[]= {1.0,0.0,0.0};
  static double yDir[]= {0.0,1.0,0.0};
  static double zDir[]= {0.0,0.0,1.0};
  double sliceAcqDir[3];
  int i;
  for (i=0; i<3; i++) 
    sliceAcqDir[i]= firstLine->slicePos[i] - lastLine->slicePos[i];
  transformSiemensToFiascoVec3(sliceAcqDir);
  if (normVec3(sliceAcqDir) == 0.0) {
    /* One slice; we can't do a thing. */
    return;
  }
  else { 
    double phaseEncodeDir[3];
    double freqEncodeDir[3];
    double xDot;
    double yDot;
    double zDot;
    normalizeVec3(sliceAcqDir);
    xDot= fabs(dotVec3(sliceAcqDir,xDir));
    yDot= fabs(dotVec3(sliceAcqDir,yDir));
    zDot= fabs(dotVec3(sliceAcqDir,zDir));
    
    /* I'm sure the quaternion in mdh can give much better info, but
     * until we understand it...
     */
    if (zDot>xDot) { 
      if (zDot>yDot) {
	/* Axial, possibly oblique */
	for (i=0; i<3; i++) freqEncodeDir[i]= xDir[i];
      }
      else {
	/* yDot greatest; Coronal */
	for (i=0; i<3; i++) freqEncodeDir[i]= xDir[i];
      }
    }
    else {
      if (yDot>xDot) {
	/* yDot largest; Coronal */
	for (i=0; i<3; i++) freqEncodeDir[i]= xDir[i];
      }
      else {
	/* xDot largest; Saggital */
	for (i=0; i<3; i++) freqEncodeDir[i]= -yDir[i];
      }
    } 
    crossVec3( phaseEncodeDir, sliceAcqDir, freqEncodeDir );
    normalizeVec3( phaseEncodeDir );
    defVec3(info,"freqEncodeDir",freqEncodeDir);
    kvDefString(defs,"freqEncodeDir.0","freq encoding dir (Fiasco coords)");
    defVec3(info,"phaseEncodeDir",phaseEncodeDir);
    kvDefString(defs,"phaseEncodeDir.0","phase encoding dir (Fiasco coords)");
    defVec3(info,"sliceAcqDir",sliceAcqDir);
    kvDefString(defs,"sliceAcqDir.0","slice acq progress dir (Fiasco coords)");
  }
}

static void infer_key_value_pairs_from_list( KVHash* info, SList* lineList,
					     MeasurementDataHeader* mdh )
{
  KVHash* defs= kvGetHash(info,"definitions");
  KVHash* extNames= kvGetHash(info,"external_names");
  Line* firstLine= NULL;
  Line* lastLine= NULL;
  char* rowflip= NULL;
  long nrows= 0;

  slist_totop(lineList);
  firstLine= (Line*)slist_get(lineList);
  lastLine= (Line*)slist_getlast(lineList);

  if (mdh->usedChannels != 1)
    Abort("%s: siemens_kspace_reader: multi-channel input not implemented (see notes in code)!\n",
	  progname);
  
  if (firstLine->nav_id == lastLine->nav_id) {
    kvDefString(info,"dimstr","vxyzt");
    nrows= (lastLine->line - firstLine->line) + 1;
    kvDefInt(info,"dy", nrows);
    kvDefString(info,"description.y","gridded k-space");
    kvDefInt(info,"kspace_ctr.y",mdh->kspaceCtrLine);
    kvDefString(defs,"kspace_ctr.y", "this row crosses the k-space origin");
    if (mdh->kspaceCtrLine == 0) {
      kvDefBoolean(info,"xchop",1);
      kvDefBoolean(info,"partialk",0);
      kvDefInt(info,"dy_base",nrows);
    }
    else if (mdh->kspaceCtrLine == kvGetInt(info,"dy")/2) {
      kvDefBoolean(info,"xchop",0);
      kvDefBoolean(info,"partialk",0);
      kvDefInt(info,"dy_base",nrows);
    }
    else {
      kvDefBoolean(info,"xchop",0);
      kvDefBoolean(info,"partialk",1);
      kvDefInt(info,"dy_base",2*mdh->kspaceCtrLine);
    }
    kvDefString(defs,"partialk","partial-k completion needed");
    kvDefString(defs,"dy_base","y samples after reconstruction and clipping");
    infer_scan_type_from_list( info, firstLine, lastLine, mdh );
  }
  else {
    kvDefString(info,"dimstr","vxnzt");
    nrows= (lastLine->nav_id - firstLine->nav_id) + 1;
    kvDefInt(info,"dn",nrows);
    kvDefString(info,"description.n","discrete");        
  }
  kvDefInt(info,"dv",2);
  kvDefString(info,"description.v","complex real/imaginary");
  kvDefInt(info,"dx",mdh->samplesInScan);
  kvDefString(info,"description.x","gridded k-space");
  if (lastLine->slice != firstLine->slice) 
    kvDefInt(info,"dz",(lastLine->slice - firstLine->slice) + 1);
  else
    kvDefInt(info,"dz",(lastLine->partition - firstLine->partition) + 1);
  kvDefString(info,"description.z","gridded image-space");
  kvDefInt(info,"dt",(lastLine->rep - firstLine->rep) + 1);
  kvDefString(info,"description.t","gridded image-space");
  kvDefInt(info,"datatype_in",SRDR_FLOAT32);
  kvDefInt(info,"handler_datatype_out",SRDR_FLOAT32);

  /* We will be handling slice reordering when we sort the scan lines.
   * The specific reordering pattern depends on the scan orientation.
   */
  kvDefInt(info,"sliceorder",1);
  kvDefString(defs,"sliceorder","slice order: 0=interleaved, 1=sequential");
  kvDefBoolean(info,"reorder",0);

  rowflip= pick_rowflip_pattern(lineList, nrows,
				kvGetString(info,"chunkname"));
  kvDefBoolean(info,"rowflip",( strcmp(rowflip,"none")!= 0) );
  kvDefString(defs,"rowflip","EPI row reversal needed");
  kvDefString(info,"rowflip_pattern",rowflip);
  kvDefString(defs,"rowflip_pattern","EPI row reversal pattern");
  free(rowflip);

  /* Everything works much better if the image values aren't so tiny */
  kvDefBoolean(info,"autoscale",1);
  kvDefDouble(info,"autoscale_range",5000.0);

  kvDefInt(info,"kspace_ctr.x",mdh->kspaceCtrCol);
  kvDefString(defs,"kspace_ctr.x", "this column crosses the k-space origin");
  if (mdh->kspaceCtrCol == 0)
    kvDefBoolean(info,"ychop",1);
  else if (mdh->kspaceCtrCol != mdh->samplesInScan/2)
    Abort("%s: siemens_kspace_reader.c: k-space center is at col %d!\n",
	  mdh->kspaceCtrCol);

  kvDefLong(info,"start_offset",0); /* we will hide all offsets */
  kvDefLong(info,"skip.x",0); /* trigger breaking of read after each line */
}

static void listRead( FileHandler* self, KVHash* info,
		      long long offset, long n,
		      SRDR_Datatype datatype_out,
		      void* obuf ) 
{
  SList* lineList= (SList*)self->hook;
  Line* thisLine;
  
  if (slist_atend(lineList))
    Abort("%s: listRead: ran out of input lines at offset %lld!\n",
	  progname,offset);
  
  thisLine= (Line*)slist_get(lineList);
  slist_next(lineList);
  
  if (n != thisLine->n)
    Abort("%s: listRead: internal error; asked for length %d, not %d!\n",
	  progname,n, thisLine->n);

  FH_REOPEN(self); /* just in case */
  
  /* Ignore the given offset and use our own */
  baseRead(self, info, thisLine->offset, n, datatype_out, obuf);
}

static char* strip( char* s )
{
  char* tail;

  if (!s) return NULL;
  while (isspace(*s)) s++;
  if (*s=='\0') return s;
  tail= s+(strlen(s)-1);
  while (isspace(*tail)) tail--;
  *(tail+1)='\0';
  return s;
}

static void parse_asc_entry( KVHash* info, char* s, int lineNum )
{
  char* key;
  char* value;
  char* cp;

  key= strip(strtok_r(s, "=", &cp));
  value= strip(strtok_r(NULL, "=", &cp));
  if (value) { /* ignore anything that isn't a key-value pair */
    int m_mode= (!strncmp(key,"m_",2));
    char* keyTail= strrchr(key,'.');
    if (keyTail) keyTail++; /* skip the '.' */
    else keyTail= key;
    if (*value=='[') value++; /* clip leading '[' */
    if (value[strlen(value)-1]==']') 
      value[strlen(value)-1]= '\0'; /* clip trailing ']' */
    value= strip(value);
    if (m_mode) {
      char buf[256];
      int i= 0;
      char* vcp;
      char* tok= strtok_r(value," ",&vcp);
      if (!strncmp(keyTail,"m_",2)) keyTail += 2;
      if (!strncmp(keyTail,"ax",2)) {
	/* A complex pair, not in a list */
	double real, imag;
	real= atof(tok);
	tok= strtok_r(NULL,"+",&vcp);
	if (!(tok= strchr(tok,'i')))
	  Abort("%s: siemens_kspace_reader: tried to read a bad complex value!\n",
		progname);
	imag= atof(tok);
	snprintf(buf,sizeof(buf),"%s.real",key);
	kvDefDouble(info,buf,real);
	snprintf(buf,sizeof(buf),"%s.imag",key);
	kvDefDouble(info,buf,imag);
      }
      else if (!strncmp(keyTail,"d",1) 
	       || !strncmp(keyTail,"fl",2)
	       || !strncmp(keyTail,"afl",3)
	       ) {
	/* double value */
	while (tok) {
	  snprintf(buf,sizeof(buf),"%s.%d",key,i);
	  kvDefDouble(info,buf,atof(tok));
	  tok= strtok_r(NULL," ",&vcp);
	  i++;
	}
      }
      else if (!strncmp(keyTail,"b",1)) {
	/* bool value */
	while (tok) {
	  snprintf(buf,sizeof(buf),"%s.%d",key,i);
	  kvDefBoolean(info,buf,atol(tok));
	  tok= strtok_r(NULL," ",&vcp);
	  i++;
	}
      }
      else if (!strncmp(keyTail,"t",1) /* ascii text */
	       || !strncmp(keyTail,"at",2) /* ascii time */
	       ) {
	while (tok) {
	  snprintf(buf,sizeof(buf),"%s.%d",key,i);
	  kvDefString(info,buf,tok);
	  tok= strtok_r(NULL," ",&vcp);
	  i++;
	}
      }
      else if (!strncmp(keyTail,"uc",2)
	       || !strncmp(keyTail,"ush",3)
	       || !strncmp(keyTail,"ui",2)
	       || !strncmp(keyTail,"ul",2)
	       || !strncmp(keyTail,"un",2)
	       || !strncmp(keyTail,"al",2)
	       || !strncmp(keyTail,"i",1)
	       || !strncmp(keyTail,"n",1)
	       || !strncmp(keyTail,"l",1)
	       || !strncmp(keyTail,"e",1) /* enumeration? */
	       ) {
	/* Unsigned int of some sort */
	while (tok) {
	  int base= 0;
	  /* This one particular tag seems to be in hex, but without
	   * the leading "0x".
	   */
	  if (!strcmp(keyTail,"alCoilPlugId")) base= 16;
	  snprintf(buf,sizeof(buf),"%s.%d",key,i);
	  kvDefLong(info,buf,strtol(tok,NULL,base));
	  tok= strtok_r(NULL," ",&vcp);
	  i++;
	}
      }
      else Warning(1,"%s: siemens_kspace_reader: unparsed key <%s>!\n",
		   progname,key);
    }
    else {
      if (*keyTail=='a') keyTail += 1;
      if (!strncmp(keyTail,"d",1) 
	  || !strncmp(keyTail,"fl",2)
	  ) {
	/* double value */
	kvDefDouble(info,key,atof(value));
      }
      else if (!strncmp(keyTail,"l",1)) {
	/* long value */
	kvDefInt(info,key,atol(value));
      }
      else if (!strncmp(keyTail,"b",1)) {
	/* bool value */
	kvDefBoolean(info,key,atol(value));
      }
      else if (!strncmp(keyTail,"t",1)) {
	/* ASCII text in quotes */
	value++;
	value[strlen(value)-1]= '\0'; 
	kvDefString(info,key,value);
      }
      else if (!strncmp(keyTail,"uc",2)
	       || !strncmp(keyTail,"ush",3)
	       || !strncmp(keyTail,"ui",2)
	       || !strncmp(keyTail,"ul",2)
	       || !strncmp(keyTail,"un",2)
	       || !strncmp(keyTail,"i",1)
	       || !strncmp(keyTail,"n",1)
	       || !strncmp(keyTail,"e",1) /* enumeration? */
	       ) {
	/* Unsigned int of some sort */
	kvDefLong(info,key,strtol(value,NULL,0));
      }
      else Warning(1,"%s: siemens_kspace_reader: unparsed key <%s>!\n",
		   progname,key);
    }
  }
}

static void parse_asc_file( KVHash* info, char* fname )
{
  KVHash* auxInfo= kvFactory(KV_DEFAULT_SIZE);
  KVHash* defs= kvGetHash(info,"definitions");
  KVHash* extNames= kvGetHash(info,"external_names");
  FILE* f= NULL;
  char buf[256];
  int lineNum= 1;
  double c0[3], c1[3], norm[3], sep[3], sepMag, normMag, vox_z, gap;

  if (!(f= fopen(fname,"r")))
    Abort("%s: siemens_kspace_reader: auxiliary file <%s> is missing or not readable!\n",
	  progname, fname);

  while (!feof(f)) {
    if (!fgets(buf, sizeof(buf), f)) {
      if (ferror(f)) {
	Abort("%s: siemens_kspace_reader: error on line %d of %s: %s\n",
		progname, lineNum, fname, strerror(errno));
	break;
      }
    }
    if (buf[0]!='\0' && buf[strlen(buf)-1]=='\n')
      buf[strlen(buf)-1]= '\0'; /* strip newline */
    parse_asc_entry( auxInfo, buf, lineNum );
    lineNum++;
  }

  if (fclose(f))
    Warning(1,"%s: siemens_kspace_reader: warning, cannot close <%s>!\n",
	    progname, fname);

#define COPY_DEF( name1, type1, name2, type2, def ) \
  {  \
     if (kvLookup(auxInfo,name1)) \
       kvDef##type2 (info,name2,kvGet##type1 (auxInfo,name1)); \
     kvDefString(defs,name2,def); \
  }
#define DEFAULT_TO_ZERO( name2 ) \
  if (!kvLookup(info,name2)) kvDefDouble(info,name2,0.0);

  /* Let us now grab our favorite parts */
  COPY_DEF( "alTE[0]", Long, "TE", Int, "TE (us)" );

  COPY_DEF( "alTR[0]", Long, "TR", Int, "TR (us)" );

  COPY_DEF( "dFlipAngleDegrees", Double, "flip", Double, "flip angle" );

  COPY_DEF( "m_aflRegridADCDuration.0", Double, "regridADCDuration", 
	    Double, "ADC open time for regridding (us)" );

  COPY_DEF( "m_alRegridDelaySamplesTime.0", Long, "regridSampleDelay", 
	    Double, "sample delay time for regridding (us)" );

  COPY_DEF( "m_alRegridDestSamples.0", Long, "regridDestSamples", 
	    Int, "samples after regridding" );

  COPY_DEF( "m_alRegridDestSamples.0", Long, "dx_resampled", 
	    Int, "samples after regridding" );

  COPY_DEF( "sKSpace.lBaseResolution", Long, "dx_base", 
	    Int, "x samples after clipping" );

  COPY_DEF( "m_alRegridFlattopTime.0", Long, "regridFlattopTime", 
	    Double, "ramp flat top time for regridding (us)" );

  COPY_DEF( "m_alRegridRampdownTime.0", Long, "regridRampdownTime", 
	    Double, "ramp down time for regridding (us)" );

  COPY_DEF( "m_alRegridRampupTime.0", Long, "regridRampupTime", 
	    Double, "ramp up time for regridding (us)" );

  COPY_DEF( "m_alRegridMode.0", Long, "regridMode", 
	    Int, "regridding mode" );

  COPY_DEF( "m_flKSpaceFilterWidth.0", Double, "kspaceFilterWidth", 
	    Double, "kspace filter width" );

  COPY_DEF( "m_aflMagneticFieldStrength.0", Double, "fieldStrength", 
	    Int, "magnetic field strength (Tesla)" );

  COPY_DEF( "m_lAbsTablePosition.0", Long, "table_delta", 
	    Double, "table position (mm)" );

  COPY_DEF( "m_lEchoSpacing.0", Long, "echoSpacing", 
	    Double, "echo spacing (us)" );

  COPY_DEF( "m_tICEProgramName.0", String, "ICEProgramName", 
	    String, "Siemens controlling ICE program" );

  COPY_DEF( "m_tPatientPosition.0", String, "patientPosition", 
	    String, "patient position" );

  COPY_DEF( "m_tSequenceString.0", String, "pulse_seq", 
	    String, "pulse sequence" );

  COPY_DEF( "m_tScanningSequence.0", String, "scanningSequence", 
	    String, "Siemens scanning sequence" );

  COPY_DEF( "m_tScanOptions.0", String, "scanOptions", 
	    String, "Siemens scan options" );

  COPY_DEF( "m_tSequenceVariant.0", String, "pulse_seq_variant", 
	    String, "pulse sequence variant" );

  COPY_DEF( "sAdjVolume.dPhaseFOV", Double, "fov_y", 
	    Double, "Y field of view (mm)" );
  if (!kvLookup(info,"fov_y"))
    COPY_DEF( "sSliceArray.asSlice[0].dPhaseFOV", Double, "fov_y", 
	      Double, "Y field of view (mm)" );

  COPY_DEF( "sAdjVolume.dReadoutFOV", Double, "fov_x", 
	    Double, "X field of view (mm)" );
  if (!kvLookup(info,"fov_x"))
    COPY_DEF( "sSliceArray.asSlice[0].dReadoutFOV", Double, "fov_x", 
	      Double, "X field of view (mm)" );

  COPY_DEF( "sSliceArray.asSlice[0].dThickness", Double, "slice_thickness", 
	    Double, "Slice thickness (mm)" );
  kvDefString( extNames, "slice_thickness", "slthick" );

  /* We now get slice_norm and slice_ctr from the linelist */
  /*
    getVec3_anat(auxInfo,"sSliceArray.asSlice[0].sNormal.d",slice_norm);
    transformSiemensToFiascoVec3(slice_norm);
    defVec3(info,"slice_norm",slice_norm);

    getVec3_anat(auxInfo,"sSliceArray.asSlice[0].sPosition.d",slice_ctr);
    transformSiemensToFiascoVec3(slice_ctr);
    defVec3(info,"slice_ctr",slice_ctr);
  */

#undef COPY_DEF
#undef DEFAULT_TO_ZERO

  /* Can we easily get the sample time? */
  if (kvLookup(info,"regridADCDuration") 
      && (kvLookup(info,"dx") || kvLookup(info,"dq"))) {
    long dx;
    if (kvLookup(info,"dx")) dx= kvGetInt(info,"dx");
    else dx= kvGetInt(info,"dq");
    kvDefDouble(info,"samp_time",kvGetDouble(info,"regridADCDuration")/dx);
    kvDefString(defs,"samp_time","sample time (usec)");
  }

  /* We need to do some algebra to get the slice gap and fov_z*/
  if (getVec3_anat(auxInfo,"sSliceArray.asSlice[0].sPosition.d",c0) 
      && getVec3_anat(auxInfo, "sSliceArray.asSlice[1].sPosition.d", c1)) {
    transformSiemensToFiascoVec3(c0);
    transformSiemensToFiascoVec3(c1);
    if (!getVec3(info,"slice_norm",norm))
      Abort("%s: siemens_kspace_reader internal error: slice_norm undefined!\n",
	    progname);
    subtractVec3(sep,c1,c0);
    normMag= normVec3( norm );
    sepMag= normVec3( sep );
    vox_z= sepMag/normMag;
    gap= vox_z - kvGetDouble(info,"slice_thickness");
    kvDefDouble(info,"slice_gap",gap);
    kvDefString(defs,"slice_gap","slice gap (mm)");
    kvDefDouble(info,"voxel_z",vox_z);
  }
  else {
    /* Let's leave the slice gap undefined, since we can't calculate it
     * kvDefDouble(info,"slice_gap",0.0);
     * kvDefString(defs,"slice_gap","slice gap (mm)");
     */
    kvDefDouble(info,"voxel_z",kvGetDouble(info,"slice_thickness"));
  }
  
  if (kvLookup(info,"slice_gap")) {
    kvDefDouble(info,"fov_z",
		kvGetInt(info,"dz")*kvGetDouble(info,"slice_thickness")
		+ (kvGetInt(info,"dz")-1)*kvGetDouble(info,"slice_gap"));
  }
  else {
    kvDefDouble(info,"fov_z",
		kvGetInt(info,"dz")*kvGetDouble(info,"voxel_z"));
  }
  kvDefString(defs,"fov_z","Z field of view (mm)");

  if (kvLookup(info,"fov_x") && kvLookup(info,"dx_base")) {
    kvDefDouble(info,"voxel_x",
		kvGetDouble(info,"fov_x")/kvGetInt(info,"dx_base"));
    kvDefString(defs,"voxel_x","X voxel size (mm)");
  }
  if (kvLookup(info,"fov_y") && kvLookup(info,"dy_base")) {
    kvDefDouble(info,"voxel_y",
		kvGetDouble(info,"fov_y")/kvGetInt(info,"dy_base"));
    kvDefString(defs,"voxel_y","Y voxel size (mm)");
  }

  if (kvLookup(info,"regridRampupTime") 
      && (kvGetDouble(info,"regridRampupTime") != 0.0)) {
    /* Ramp resampling is needed */
    char* c;
    kvDefBoolean(info,"resample",1);
    kvDefString(defs,"resample","regridding for ramp sampling required");
    strncpy(buf,kvGetString(info,"dimstr"),sizeof(buf));
    if (!(c=strchr(buf,'x')))
      Abort("%s: siemens_kspace_reader.c: internal error; no dimension x!\n",
	    progname);
    *c= 'q';
    kvDefString(info,"dimstr",buf);
    kvDefInt(info,"dq",kvGetInt(info,"dx"));
    kvDefString(info,"description.q","ungridded k-space");
    kvDelete(info,"dx");
    kvDelete(info,"description.x");
    if (kvLookup(info,"skip.x")) {
      kvDefLong(info,"skip.q",kvGetLong(info,"skip.x"));
      kvDelete(info,"skip.x");
    }
    kvDefString(info,"resample_method",SIEMENS_RESAMPLE_SCRIPT);
  }
}

static void processHeader( FileHandler* self, KVHash* info, SList* cStack )
{
  SList* lineList= slist_create();
  SList* navList= slist_create();
  FILE *fphead;
  ChunkHandlerPair* chPair= NULL;
  FileHandler* navHandler= NULL;
  KVHash* navInfo= kvFactory(KV_DEFAULT_SIZE);
  MeasurementDataHeader mdh;
  int i;
  int ierror= 0;

  /* Call the base class method */
  baseProcessHeader( self, info, cStack );

  /* This bit is from the GE code io_signa_lx.c, with mods for portability */
  ierror= 0;
  if ((fphead = fopen(self->fileName,"r"))!=NULL)
    {
      if (fseek(fphead, SMNSKSPC_HEADER_OFF, SEEK_SET)) {
	perror("Error seeking header");
	ierror=1;
      }
      else {
	/* Parity check */
	unsigned char buf[4];
	if (fread(buf, sizeof(long), 1, fphead) != 1) {
	  perror("Error reading header");
	  ierror=1;
	}
	else {
	  long firstWord= BRdInt32(buf);
	  if (firstWord != SMNSKSPC_HEADER_SIZE_BYTES) {
	    /* Oops, try it the other way! */
	    bio_big_endian_input= (bio_big_endian_input ? 0 : 1);
	    firstWord= BRdInt32(buf);
	    if (firstWord != SMNSKSPC_HEADER_SIZE_BYTES) 
	      Abort("%s: internal error: unexpectedly found wrong first word in file\n",progname);
	  }
	  kvDefBoolean(info,"big_endian_input",bio_big_endian_input);
	}

	/* this will find all the acquisition lines in the file */
	walk_mdh_structures( info, self, fphead, &lineList, &navList );

	/* we need one representative mdh to infer some file attributes */
	(void)parse_mdh(fphead, SMNSKSPC_HEADER_SIZE_BYTES, &mdh);
	if (debug) {
	  fprintf(stderr,"Prototype MeasurementDataHeader:\n");
	  dump_mdh(stderr,&mdh);
	}

	if (fclose(fphead)) {
	  perror("Error closing header");
	  ierror=1;
	}
      }
    }
  else {
    perror("Error opening header");
    ierror= 1;
  }
  if (ierror)
    Abort("siemens_kspace_reader: unable to read or parse header from <%s>!\n",
	  self->fileName);

  /* Now that we have all the data, we can make inferences about the
   * file structure.  We also need to implement slice reordering
   * based on those inferred values.  This is done by re-sorting on
   * corrected slice numbers.
   */
  infer_key_value_pairs_from_list( info, lineList, &mdh );
  relabel_lines( info, lineList );
  slist_sort(lineList, line_compare);

  /* If there is an auxiliary file, try to parse it */
  if (kvLookup(info,"auxfile"))
    parse_asc_file(info, kvGetString(info,"auxfile"));

  /* This routine may do nothing if we don't have enough info */
  calcVolumeBounds(info);

  /* Set up to handle the main chunk ("images") */
  slist_totop(lineList);
  self->hook= lineList;

  /* If there is a navigator chunk, modify a base-class 
     FileHandler to handle it */
  if (!slist_empty(navList)) {
    double sliceAcqDir[3];
    double phaseEncodeDir[3];
    double freqEncodeDir[3];
    initInfoHash(navInfo);
    kvDefString(navInfo,"chunkname","navigator");
    kvDefString(navInfo,"chunkfile",".nav");
    getVec3(info,"sliceAcqDir",sliceAcqDir);
    defVec3(navInfo,"sliceAcqDir",sliceAcqDir);
    getVec3(info,"phaseEncodeDir",phaseEncodeDir);
    defVec3(navInfo,"phaseEncodeDir",phaseEncodeDir);
    getVec3(info,"freqEncodeDir",freqEncodeDir);
    defVec3(navInfo,"freqEncodeDir",freqEncodeDir);
    infer_key_value_pairs_from_list( navInfo, navList, &mdh );
    relabel_lines( navInfo, navList ); /* causing slice sorting; see above */
    slist_sort(navList, line_compare);
    slist_totop(navList);
    navHandler= baseFactory(self->fileName);
    navHandler->destroySelf= destroySelf;
    navHandler->typeName= strdup( "SiemensNavigator" );
    navHandler->read= listRead;
    navHandler->hook= navList;
    if (!(chPair= (ChunkHandlerPair*)malloc(sizeof(ChunkHandlerPair))))
      Abort("%s: unable to allocate %d bytes!\n",
	    progname,sizeof(ChunkHandlerPair));
    chPair->info= navInfo;
    chPair->handler= navHandler;
    slist_push(cStack,chPair);
  }
}

FileHandler* siemensKspaceFactory(char* fname, KVHash* info)
{
  FileHandler* result= baseFactory(fname);
  result->processHeader= processHeader;
  result->destroySelf= destroySelf;
  result->typeName= strdup( "SiemensKspace" );
  result->read= listRead;
  return result;
}

int siemensKspaceTester(const char* filename)
{
  char buf[4];
  int ierror= 0;
  int match= 0;
  FILE* fphead= NULL;

  /* This bit is from the GE code io_signa_lx.c, with mods for portability */
  if ((fphead = fopen(filename,"r"))!=NULL)
    {
      if (fseek(fphead, SMNSKSPC_HEADER_OFF, SEEK_SET)) {
	perror("Error seeking header");
	ierror=1;
      }
      else {
	if (fread(buf, sizeof(long), 1, fphead) != 1) {
	  perror("Error reading header");
	  ierror=1;
	}
      }
    }
  else {
    perror("Error opening header");
    ierror= 1;
  }
  if (ierror) {
    if (fphead) fclose(fphead);
    return 0;
  }

  /* Try both endian orders */

  /* We test to see if the data does have the expected endian order-
     first word should be known value */
  if (BRdInt32((unsigned char*)buf) == SMNSKSPC_HEADER_SIZE_BYTES) {
    match= 1;
  }
  else {
    bio_big_endian_input= (bio_big_endian_input ? 0 : 1);
    if (BRdInt32((unsigned char*)buf) == SMNSKSPC_HEADER_SIZE_BYTES) {
      match= 1;
    }
  }

  if (match) match= reality_check_file(fphead);
  if (fphead) fclose(fphead);

  return match;
}

