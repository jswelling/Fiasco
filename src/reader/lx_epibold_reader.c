/************************************************************
 *                                                          *
 *  lx_epibold_reader.c                                     *
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

static char rcsid[] = "$Id: lx_epibold_reader.c,v 1.13 2004/10/26 23:55:58 welling Exp $";

#define BLANKY_NUM_DEFAULT 2
#define BLANKY_NUM_MAX 10

#define GE_RESAMPLE_SCRIPT "ge_ramp_resample.csh"

static int
IntBRdFloat32 (unsigned char *addr)
{
  return((int) Round(BRdFloat32(addr)));
}

void ADD_VSUFFIX(scanEpiboldHeader)( KVHash* info, 
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
  if (rcxres != xres) {
    /* ramp resampling required */
    if (!kvLookup(info,"rampfile")) 
      Warning(1,"%s: scan_epibold_header: ramp resampling needed, but no ramp file given!\n",
	      progname);
    kvDefString(info,"dimstr","vqyzt");
    kvDefInt(info,"dx_resampled",rcxres);
    kvDefInt(info,"dq",xres);
    kvDefString(info,"description.q","ungridded k-space");
    kvDefBoolean(info,"resample",1);
    kvDefString(defs,"resample","regridding for ramp sampling required");
    kvDefString(info,"resample_method",GE_RESAMPLE_SCRIPT);    
  }
  else {
    kvDefString(info,"dimstr","vxyzt");
    kvDefInt(info,"dx",xres);
    kvDefString(info,"description.x","gridded k-space");
  }

  kvDefInt(info,"dy",vpsht);
  kvDefString(info,"description.y","gridded k-space");
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
  kvDefString(info,"description.z","gridded image-space");
  kvDefInt(info,"dz",kvGetInt(info,"nslices"));
  kvDefString(info,"description.y","gridded image-space");
  
  kvDefBoolean(info,"reorder",1);
  kvDefString(info,"reorder_pattern","even/odd");
  kvDefLong(info,"skip",0);
  kvDefLong(info,"sliceskip",0);  

  kvDefBoolean(info,"rowflip",1);
  kvDefString(defs,"rowflip","EPI row reversal needed");
  kvDefString(info,"rowflip_pattern","odd");
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

void ADD_VSUFFIX(addBandpassChunk)(FileHandler* self, 
				   KVHash* info, SList* cStack)
{
  /* Add chunks for bandpass directory information */
  KVHash* defs= kvGetHash(info,"definitions");
  char buf[1024];
  char linebuf[256];
  ChunkHandlerPair* chPair= NULL;
  KVHash* subInfo= kvFactory(KV_DEFAULT_SIZE);
  
  double bandw= kvGetDouble(info,"bandwidth");
  FILE *f;
  int index;
  float v1, v2, v3;
  int nbpc;
  float x;
  int i;
  int j;
  float* table;
  float* tablerunner;
  
  /* Determine bp correction file to read based on bandwidth */
  strncpy(buf, kvGetString(info,"bandpassdir"),512);
  if (bandw > 62.5) /* fast receiver */
    if (bandw <= 100.0) strcat(buf,"/bcrvf1.dat");
    else if (bandw <= 200.0) strcat(buf,"/bcrvf2.dat");
    else if (bandw <= 300.0)  strcat(buf,"/bcrvf3.dat");
    else if (bandw <= 400.0)  strcat(buf,"/bcrvf4.dat");
    else strcat(buf,"/bcrvf5.dat");
  else strcat(buf,"/bcrvf1.dat"); /* digital receiver */
  kvDefString(info,"bandpassfile",buf);
  kvDefString(defs,"bandpassfile",
	      "file from which bandpass filter was read");
  
  /* Count lines in the bandpass file */
  if (!(f= fopen(buf,"r")))
    Abort("%s: Couldn't open file <%s> to read bandpass info!\n",
	  progname,buf);
  i= 0;
  while (!feof(f)) {
    if (fgets(linebuf, sizeof(linebuf), f)==NULL) break;
    if (strlen(linebuf)!=0) i++;
  }
  nbpc= i;
  (void)fclose(f);
  
  if (!(table=(float*)malloc(2*nbpc*sizeof(float)))) 
    Abort("%s: unable to allocate %d bytes!\n",2*nbpc*sizeof(float));
  
  /* Read the table for content, and translate the relevant bits to
   * real-imaginary representation.
   */
  if (!(f= fopen(buf,"r")))
    Abort("%s: Couldn't open file <%s> to read bandpass info (pass 2)!\n",
	  progname,buf);
  tablerunner= table;
  for (i=0; i<nbpc; i++) {
    float x;
    float mag;
    float phase;
    float val_r;
    float val_i;
    if (fscanf(f,"%f %f %f",&x, &mag, &phase)!=3) 
      Abort("%s: bandpass table <%s> is corrupted at line %d!\n",
	    progname, buf, i);
    
    val_r= mag*cos(phase*(M_PI/180));
    *tablerunner++= val_r;
    val_i= mag*sin(phase*(M_PI/180));
    *tablerunner++= val_i;
  }
  (void)fclose(f);
  
  initInfoHash(subInfo);
  kvDefString(subInfo,"chunkname","bandpass");
  kvDefString(subInfo,"chunkfile",".band");
  kvDefString(subInfo,"dimstr","vl");
  kvDefInt(subInfo,"dv",2);
  kvDefString(subInfo,"description.v","complex real/imaginary");
  kvDefInt(subInfo,"dl",nbpc);
  kvDefLong(subInfo,"start_offset",0);
  kvDefInt(subInfo,"datatype_in",SRDR_FLOAT32);
  kvDefInt(subInfo,"handler_datatype_out",SRDR_FLOAT32);
  
  if (!(chPair= (ChunkHandlerPair*)malloc(sizeof(ChunkHandlerPair))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,sizeof(ChunkHandlerPair));
  chPair->info= subInfo;
  chPair->handler= 
    ramDataHandlerFactory(table, 2*nbpc, SRDR_FLOAT32);
  slist_push(cStack,chPair);
}

void ADD_VSUFFIX(addPhaseRefChunk)(FileHandler* self, 
				   KVHash* info, SList* cStack)
{
  /* Add phase reference dataset.  This is actually two chunks,
   * for the constant and linear correction terms. */
  KVHash* defs= kvGetHash(info,"definitions");
  ChunkHandlerPair* chPair= NULL;
  KVHash* subInfo= kvFactory(KV_DEFAULT_SIZE);
  FileHandler* subHandler= rawFactory(kvGetString(info,"phasereffile"), info);
  
  initInfoHash(subInfo);
  kvDefString(subInfo,"chunkname","phaseref_const");
  kvDefString(subInfo,"chunkfile",".poff");
  kvDefString(subInfo,"dimstr","ys");
  kvDefInt(subInfo,"dy",kvGetInt(info,"dy"));
  kvDefString(subInfo,"description.y","gridded k-space");
  kvDefInt(subInfo,"ds",kvGetInt(info,"ds"));
  kvDefString(subInfo,"description.y","gridded k-space");
  kvDefLong(subInfo,"start_offset",0);
  kvDefInt(subInfo,"datatype_in",SRDR_FLOAT32);
  kvDefInt(subInfo,"handler_datatype_out",SRDR_FLOAT32);
  
  if (!(chPair= (ChunkHandlerPair*)malloc(sizeof(ChunkHandlerPair))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,sizeof(ChunkHandlerPair));
  chPair->info= subInfo;
  chPair->handler= subHandler;
  slist_push(cStack,chPair);
  
  chPair= NULL;
  subInfo= kvFactory(KV_DEFAULT_SIZE);
  subHandler= rawFactory(kvGetString(info,"phasereffile"),info);
  
  initInfoHash(subInfo);
  kvDefString(subInfo,"chunkname","phaseref_linear");
  kvDefString(subInfo,"chunkfile",".pl");
  kvDefString(subInfo,"dimstr","ys");
  kvDefInt(subInfo,"dy",kvGetInt(info,"dy"));
  kvDefString(subInfo,"description.y","gridded k-space");
  kvDefInt(subInfo,"ds",kvGetInt(info,"ds"));
  kvDefString(subInfo,"description.y","gridded k-space");
  kvDefLong(subInfo,"start_offset",512*sizeof(float));
  kvDefInt(subInfo,"datatype_in",SRDR_FLOAT32);
  kvDefInt(subInfo,"handler_datatype_out",SRDR_FLOAT32);
  
  if (!(chPair= (ChunkHandlerPair*)malloc(sizeof(ChunkHandlerPair))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,sizeof(ChunkHandlerPair));
  chPair->info= subInfo;
  chPair->handler= subHandler;
  slist_push(cStack,chPair);
}

void ADD_VSUFFIX(addRampSampleChunk)(FileHandler* self, 
				     KVHash* info, SList* cStack)
{
  KVHash* defs= kvGetHash(info,"definitions");
  /* Add ramp sampling dataset */
  ChunkHandlerPair* chPair= NULL;
  KVHash* subInfo= kvFactory(KV_DEFAULT_SIZE);
  FileHandler* subHandler= rawFactory(kvGetString(info,"rampfile"), info);
  
  initInfoHash(subInfo);
  kvDefString(subInfo,"chunkname","rampsample");
  kvDefString(subInfo,"chunkfile",".ramp");
  kvDefString(subInfo,"dimstr","qx");
  kvDefInt(subInfo,"dq",kvGetInt(info,"dq"));
  kvDefString(subInfo,"description.q","ungridded k-space");
  kvDefInt(subInfo,"dx",kvGetInt(info,"dx_resampled"));
  kvDefString(subInfo,"description.x","gridded k-space");
  kvDefLong(subInfo,"start_offset",0);
  kvDefInt(subInfo,"datatype_in",SRDR_FLOAT32);
  kvDefInt(subInfo,"handler_datatype_out",SRDR_FLOAT32);
  
  if (!(chPair= (ChunkHandlerPair*)malloc(sizeof(ChunkHandlerPair))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,sizeof(ChunkHandlerPair));
  chPair->info= subInfo;
  chPair->handler= subHandler;
  slist_push(cStack,chPair);
}
