/*
 *	header.c
 *
 *    Routines for reading a P file header used as input to the
 *    "spiral" program.  For more information about
 *    the "spiral program, see spiral.c.
 *
 *    Copyright (c) 1998 by Douglas C. Noll and the University of Pittsburgh and
 *	  the Pittsburgh Supercomputing Center 
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
 *    HISTORY
 *	1/98 - split off from spiral.c (Greg Hood, PSC)
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "fmri.h"
#include "bio.h"
#include "array.h"
#include "acct.h"
#include "spiral.h"
#include "stdcrg.h"
#include "frozen_header_info.h"

static char rcsid[] = "$Id: header.c,v 1.20 2003/03/29 00:15:20 welling Exp $";

static void CalculateLocation (float *p1, float *p2, float *p3, int rot, int trans);
static int IntBRdFloat32 (unsigned char *addr);

void
ReadFileHeader (const Filename input_file, Context *c, unsigned char *hdr)
{
  FILE *f;
  unsigned char local_hdr[FRZ_RDB_HEADER_SIZE_BYTES];
  unsigned char *h;
  float s_i_offset;
  float tempr2,tempr3;
  float p1[3],p2[3],p3[3];
  float n2[3],n3[3];
  int i, k, n;
  float tmp1, tmp2;
  char com[sizeof(Filename)+128];
  float x_size;		/* the physical size of a voxel in mm */
  float y_size;
  float z_size;
  float x_spacing;	/* the physical spacing between adjacent slices */
  float y_spacing;
  float z_spacing;
  float revnum;         /* read header revision number to check endian-ness */

  if (hdr != NULL)
    h = hdr;
  else
    h = &local_hdr[0];

  Acct(READOPEN);
  /* open the file */
#ifdef AFS
  if (strncmp(input_file, "/afs", 4) == 0)
    {
      sprintf(com, "fs flushvolume -path %s", input_file);
      system(com);
    }
#endif
  if ((f = fopen(input_file, "r")) == NULL)
    Abort("Can't open %s!\n", input_file);
  
  /* read the header */
  Acct(READING);
  if (fread(h, FRZ_RDB_HEADER_SIZE_BYTES, 1, f) != 1)
    Abort("Can't read header of %s\n", input_file);
  Acct(READCLOSE);
  fclose(f);
  Acct(PROCESSING);

  /* Check endian-ness */
  /* We test to see if the data does have the expected endian order-
   revision number should be in the neighborhood of 7.0 */
  bio_big_endian_input = c->big_endian_input;
  revnum= BRdFloat32(&h[FRZ_RDBHEAD_RDB_HDR_RDBM_REV_OFF]);
  if (revnum<1.0 || revnum>1000.0) {
    /* Oops, try it the other way! */
    bio_big_endian_input= (bio_big_endian_input ? 0 : 1);
    c->big_endian_input= bio_big_endian_input;
    revnum= BRdFloat32(&h[FRZ_RDBHEAD_RDB_HDR_RDBM_REV_OFF]);
    Message("Found unexpected Pfile byte order; compensating!\n");
  }

  /* extract the context info */
  c->nph1 = 1;
  c->nphmult = IntBRdFloat32(&h[FRZ_RDBHEAD_RDB_HDR_USER1_OFF]);
  if (c->nphmult < 1)
    c->nphmult = 1;
  if (c->gen_map || c->lin_map) c->nimages= 2;
  else c->nimages = c->nph1 * c->nphmult;
  c->nslices = BRdInt16(&h[FRZ_RDBHEAD_RDB_HDR_NSLICES_OFF]) / c->nph1; 
  c->ndatfr = BRdInt16(&h[FRZ_RDBHEAD_RDB_HDR_FRAME_SIZE_OFF]);
  c->npr = IntBRdFloat32(&h[FRZ_RDBHEAD_RDB_HDR_USER4_OFF]);
  c->chop = 1;
  c->nex = 1;
  c->densamp = 5;
  /* Slice order: 0=interleaved, 1=sequential */
  c->sliceorder = IntBRdFloat32(&h[FRZ_RDBHEAD_RDB_HDR_USER16_OFF]); 
  /* glover cv's */
  c->gtype = IntBRdFloat32(&h[FRZ_RDBHEAD_RDB_HDR_USER5_OFF]);
  if (c->gtype==1) c->reverse= 1;
  else c->reverse= 0;
  c->opxres = IntBRdFloat32(&h[FRZ_RDBHEAD_RDB_HDR_USER3_OFF]);
  c->slewrate = BRdFloat32(&h[FRZ_RDBHEAD_RDB_HDR_USER7_OFF]);
  c->ngap = IntBRdFloat32(&h[FRZ_RDBHEAD_RDB_HDR_USER10_OFF]); /* concat readout */
  c->concat = IntBRdFloat32(&h[FRZ_RDBHEAD_RDB_HDR_USER11_OFF]); /*concat readout*/
  c->ndat = c->ndatfr*(c->concat+1);
  c->fast_rec_lpf = 
    BRdFloat32(&h[FRZ_RDBHEAD_RDB_HDR_USER12_OFF]); /* fast rec cutoff in kHz */
  if (!(c->mapdel_set)) {
    c->mapdel_set= 1; /* this routine isn't smart enough to know it's invalid */
    c->mapdel = 
      BRdFloat32(&h[FRZ_RDBHEAD_RDB_HDR_USER15_OFF]); /* field mapping offset (us) */
  }
  c->samp_time = BRdFloat32(&h[FRZ_RDBHEAD_RDB_HDR_USER13_OFF]); /* in (us) */
  c->fsgcm = BRdFloat32(&h[FRZ_RDBHEAD_RDB_HDR_USER6_OFF]);
  c->risetime = (int)rint(c->fsgcm*10000.0/c->slewrate);
  c->opfov = BRdFloat32(&h[FRZ_RDBHEAD_RDB_HDR_USER0_OFF]);
  c->ts = c->samp_time * 1e-6; /* in sec */
  c->gts = 4 * 1e-6; /* 4 us gradient spacing in sec */
  c->ncoils = ( BRdInt16(&h[FRZ_RDBHEAD_RDB_HDR_DAB_0_STOP_RCV_OFF]) -
		BRdInt16(&h[FRZ_RDBHEAD_RDB_HDR_DAB_0_START_RCV_OFF]) + 1 );
  c->coil_record_length = 
    BRdInt32(&h[FRZ_RDBHEAD_RDB_HDR_RAW_PASS_SIZE_OFF]) / c->ncoils;
  strncpy(c->time, (char*)&h[FRZ_RDBHEAD_RDB_HDR_SCAN_TIME_OFF], 8);
  c->time[8] = '\0';
  strncpy(c->date, (char*)&h[FRZ_RDBHEAD_RDB_HDR_SCAN_DATE_OFF], 10);
  c->date[10] = '\0';
  c->slthick = 
    BRdFloat32(&h[FRZ_RDB_MRIMAGEDATATYPE_OFF + FRZ_IMAGEHEAD_SLTHICK_OFF]);
  c->tr = 
    BRdInt32(&h[FRZ_RDB_MRIMAGEDATATYPE_OFF + FRZ_IMAGEHEAD_TR_OFF]);
  c->te =
    BRdInt32(&h[FRZ_RDB_MRIMAGEDATATYPE_OFF + FRZ_IMAGEHEAD_TE_OFF]);
  c->flip = 
    BRdInt16(&h[FRZ_RDB_MRIMAGEDATATYPE_OFF + FRZ_IMAGEHEAD_MR_FLIP_OFF]);
  c->spacing = 
    BRdFloat32(&h[FRZ_RDB_MRIMAGEDATATYPE_OFF 
		 + FRZ_IMAGEHEAD_SCANSPACING_OFF]);
  c->tlc[0] = 
    BRdFloat32(&h[FRZ_RDB_MRIMAGEDATATYPE_OFF + FRZ_IMAGEHEAD_TLHC_R_OFF]);
  c->tlc[1] =
    BRdFloat32(&h[FRZ_RDB_MRIMAGEDATATYPE_OFF + FRZ_IMAGEHEAD_TLHC_A_OFF]);
  c->tlc[2] =
    BRdFloat32(&h[FRZ_RDB_MRIMAGEDATATYPE_OFF + FRZ_IMAGEHEAD_TLHC_S_OFF]);
  c->trc[0] =
    BRdFloat32(&h[FRZ_RDB_MRIMAGEDATATYPE_OFF + FRZ_IMAGEHEAD_TRHC_R_OFF]);
  c->trc[1] =
    BRdFloat32(&h[FRZ_RDB_MRIMAGEDATATYPE_OFF + FRZ_IMAGEHEAD_TRHC_A_OFF]);
  c->trc[2] =
    BRdFloat32(&h[FRZ_RDB_MRIMAGEDATATYPE_OFF + FRZ_IMAGEHEAD_TRHC_S_OFF]);
  c->brc[0] =
    BRdFloat32(&h[FRZ_RDB_MRIMAGEDATATYPE_OFF + FRZ_IMAGEHEAD_BRHC_R_OFF]);
  c->brc[1] =
    BRdFloat32(&h[FRZ_RDB_MRIMAGEDATATYPE_OFF + FRZ_IMAGEHEAD_BRHC_A_OFF]);
  c->brc[2] =
    BRdFloat32(&h[FRZ_RDB_MRIMAGEDATATYPE_OFF + FRZ_IMAGEHEAD_BRHC_S_OFF]);
  strncpy(c->psd, 
	  (char*)&h[FRZ_RDB_MRIMAGEDATATYPE_OFF + FRZ_IMAGEHEAD_PSDNAME_OFF],
	  33);
  c->psd[33] = '\0';

  c->baseline_length = c->ndat*4;
  c->slice_record_length = 
    c->nimages*c->npr*c->ndat*4 + c->nph1*c->baseline_length;
  if (((c->nphmult * c->npr) % 2) == 1)
    c->slice_record_length += c->nph1*c->ndat*4;
  c->os_res = c->over_samp * c->res;
  c->total_images = c->nimages * c->nfiles;

  c->rotation = BRdInt16(&h[FRZ_RDBHEAD_RDB_HDR_ROTATION_OFF]);
  c->transpose = BRdInt16(&h[FRZ_RDBHEAD_RDB_HDR_TRANSPOSE_OFF]);

  if (c->verbosity >= VERBOSITY_MINIMAL)
    Report("Processing raw file %s\n", input_file);
  if (c->verbosity >= VERBOSITY_FILE)
    {
      Report("\tnumber of time points = %d (%d*%d)\n", c->nimages,
	     c->nph1, c->nphmult);
      Report("\tnumber of slices = %d\n", c->nslices);
      Report("\tslice order is %s\n",c->sliceorder?"sequential":"interleaved");
      Report("\tnumber of samples in readout = %d\n", c->ndat);
      Report("\tnumber of samples in each acquired frame = %d\n", c->ndatfr);
      Report("\tfield of view (cm?) = %.1f\n", c->opfov);
      Report("\tnumber of spirals = %d\n", c->npr);
      Report("\ttime of scan = %s\n",c->time);
      Report("\tdate of scan = %s\n",c->date);
      Report("\n");
      Report("\ttotal # of P files = %d\n", c->nfiles);
      Report("\ttotal # of time points = %d\n", c->total_images);
      Report("\n");
      Report("\tpulse sequence descriptor: <%s>\n",c->psd);
      Report("\ttr = %d us, te= %d us, flip = %d degrees\n",
	     c->tr, c->te, c->flip);
      Report("\tfsgcm = %.1f, slewrate = %.1f, risetime = %d\n",
	     c->fsgcm, c->slewrate, c->risetime); 
      Report("\tsamp_time = %f us\n", c->samp_time);
      Report("\tconcat = %d, ngap = %d\n",c->concat, c->ngap); 
      Report("\topxres = %d\n",c->opxres);
      Report("\tfast_rec_lpf = %f kHz, mapdel = %f us\n",
	     c->fast_rec_lpf,c->mapdel);
      Report("\tts = %f, gts = %f\n",c->ts*1e6,c->gts*1e6);
      Report("\tnumber of coils = %d\n", c->ncoils);
      Report("\tslice thickness = %f mm, slice spacing = %f mm\n",
	     c->slthick, c->spacing);
      Report("\trot = %d\n",c->rotation);
      Report("\ttrans = %d\n",c->transpose);
      Report("\taps tg = %d\n",
	     BRdInt32(&h[FRZ_RDBHEAD_RDB_HDR_PS_APS_TG_OFF]));
      Report("\tmps tg = %d\n",
	     BRdInt32(&h[FRZ_RDBHEAD_RDB_HDR_PS_MPS_TG_OFF]));
      Report("\taps freq = %d\n",
	     BRdInt32(&h[FRZ_RDBHEAD_RDB_HDR_PS_APS_FREQ_OFF]));
      Report("\txshift = %f\n",
	     BRdFloat32(&h[FRZ_RDBHEAD_RDB_HDR_XOFF_OFF]));
      Report("\tyshift = %f\n",
	     BRdFloat32(&h[FRZ_RDBHEAD_RDB_HDR_YOFF_OFF]));
    }
  if (strcmp(c->psd, "splx1"))
    Error("This Pfile is of type <%s> rather than <splx1>!\n",c->psd);
  
  c->ctr[0] = 
    BRdFloat32(&h[FRZ_RDB_MRIMAGEDATATYPE_OFF + FRZ_IMAGEHEAD_CTR_R_OFF]);
  c->ctr[1] =
    BRdFloat32(&h[FRZ_RDB_MRIMAGEDATATYPE_OFF + FRZ_IMAGEHEAD_CTR_A_OFF]);
  c->ctr[2] =
    BRdFloat32(&h[FRZ_RDB_MRIMAGEDATATYPE_OFF + FRZ_IMAGEHEAD_CTR_S_OFF]);

  for (i = 0; i < 3; ++i)
    {
      p1[i]= BRdFloat32(&h[FRZ_RDB_DATAACQ_OFF + FRZ_DATAACQ_GW_POINT1_0_OFF
			  + i*4]);
      p2[i]= BRdFloat32(&h[FRZ_RDB_DATAACQ_OFF + FRZ_DATAACQ_GW_POINT2_0_OFF
			  + i*4]);
      p3[i]= BRdFloat32(&h[FRZ_RDB_DATAACQ_OFF + FRZ_DATAACQ_GW_POINT3_0_OFF
			  + i*4]);
    }
  CalculateLocation(p1, p2, p3, c->rotation, c->transpose);
  s_i_offset = c->tlc[2] - p1[2];
  c->ctr[2] -= s_i_offset;

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
  if (c->loc_shift)
    {
      if (c->reverse) {
	c->pix_shifth= 
	  BRdFloat32(&h[FRZ_RDBHEAD_RDB_HDR_YOFF_OFF])/128.0;
	c->pix_shiftv= 
	  BRdFloat32(&h[FRZ_RDBHEAD_RDB_HDR_XOFF_OFF])/128.0;
      }
      else {
	c->pix_shifth= 
	  BRdFloat32(&h[FRZ_RDBHEAD_RDB_HDR_YOFF_OFF])
	  / BRdInt16(&h[FRZ_RDBHEAD_RDB_HDR_IM_SIZE_OFF]);
	c->pix_shiftv= 
	  BRdFloat32(&h[FRZ_RDBHEAD_RDB_HDR_XOFF_OFF])
	  / BRdInt16(&h[FRZ_RDBHEAD_RDB_HDR_IM_SIZE_OFF]);
      }
      if (c->verbosity > VERBOSITY_FILE) 
	{
	  Report("  shifth = %f, shiftv = %f (in mm)\n", 
		 c->pix_shifth*tempr2, c->pix_shiftv*tempr2);
	  
	  for (n = 0; n < c->nslices; n++)
	    {
	      ReadSliceCoords(c, h, n, p1, p2, p3, NULL);
	      Report("slice %d    R/L       A/P       S/I\n",n+1);
	      Report(" TL:  %f, %f, %f\n",p1[0],p1[1],p1[2]);
	      Report(" TR:  %f, %f, %f\n",p2[0],p2[1],p2[2]);
	      Report(" BR:  %f, %f, %f\n",p3[0],p3[1],p3[2]);
	    }    
	}
    }
  else
    c->pix_shifth = c->pix_shiftv = 0.0;

  /* compute the rotation & scaling factors */
  c->factxx = -1.0/c->zoom; c->factxy = 0.0; c->factyx = 0.0; 
  c->factyy = -1.0/c->zoom; 
  if (c->transpose)
    {
      tmp1 = c->factyx; tmp2 = c->factyy;
      c->factyx = c->factxx; c->factyy = c->factxy;
      c->factxx = tmp1; c->factxy = tmp2;
    }
  for (k = 0; k < c->rotation; k++)
    {
      tmp1 = c->factyx; tmp2 = c->factyy;
      c->factyx = -c->factxx; c->factyy = -c->factxy;
      c->factxx = tmp1; c->factxy = tmp2;
    }

  /* allocate the t2k and kdens arrays */
  c->t2k = Alloc3DFloatArray(2, c->npr, c->ndat);
  c->kdens = (float *) malloc(c->ndat * sizeof(float));
  c->refl= Alloc2DFloatArray(c->nslices, 3);

  if (!(c->ref_missing= (unsigned char*)malloc(c->nslices))) 
    Abort("Unable to allocate %d bytes!\n",c->nslices);
  for (k=0; k<c->nslices; k++) c->ref_missing[k]= 0;

  if (c->hc_cor)
    {
      c->refim_in = Alloc3DFloatArray(c->nslices, c->res, c->res);
      c->sampim = Alloc3DFloatArray(c->nslices, c->os_res, c->os_res);
    }


  /* set the starting and ending slices */
  if (c->slice > 0)
    c->start_slice = c->end_slice = c->slice - 1;
  else
    {
      c->start_slice = 0;
      c->end_slice = c->nslices - 1;
    }

}

void ReadSliceCoords(Context* c, unsigned char* h, int islice,
		     float p1[3], float p2[3], float p3[3],
		     float* s_i_offset_out)
{
  int i;
  float s_i_offset;
  float p1_0[3], p2_0[3], p3_0[3];
  short rotation;
  short transpose;

  if (islice<0 || islice>=c->nslices)
    Abort("ReadSliceCoords: slice %d is out of range!\n",islice);

  rotation = BRdInt16(&h[FRZ_RDBHEAD_RDB_HDR_ROTATION_OFF]);
  transpose = BRdInt16(&h[FRZ_RDBHEAD_RDB_HDR_TRANSPOSE_OFF]);

  /* Need to do the calculation for slice 0 to get s_i_offset */
  for (i = 0; i < 3; i++)
    {
      p1_0[i]= BRdFloat32(&h[FRZ_RDB_DATAACQ_OFF 
			    + FRZ_DATAACQ_GW_POINT1_0_OFF
			    + i*4]);
      p2_0[i]= BRdFloat32(&h[FRZ_RDB_DATAACQ_OFF 
			    + FRZ_DATAACQ_GW_POINT2_0_OFF
			    + i*4]);
      p3_0[i]= BRdFloat32(&h[FRZ_RDB_DATAACQ_OFF 
			    + FRZ_DATAACQ_GW_POINT3_0_OFF
			    + i*4]);
    }
  CalculateLocation(p1_0, p2_0, p3_0, rotation, transpose);
  s_i_offset = c->tlc[2] - p1_0[2];

  if (s_i_offset_out != NULL) *s_i_offset_out= s_i_offset;

  /* Now do the calculation in earnest */
  for (i = 0; i < 3; i++)
    {
      p1[i]= BRdFloat32(&h[FRZ_RDB_DATAACQ_OFF 
			  + FRZ_DATAACQ_GW_POINT1_0_OFF
			  + islice*FRZ_RDB_SLICE_INFO_ENTRY_SIZE 
			  + i*4]);
      p2[i]= BRdFloat32(&h[FRZ_RDB_DATAACQ_OFF 
			  + FRZ_DATAACQ_GW_POINT2_0_OFF
			  + islice*FRZ_RDB_SLICE_INFO_ENTRY_SIZE 
			  + i*4]);
      p3[i]= BRdFloat32(&h[FRZ_RDB_DATAACQ_OFF 
			  + FRZ_DATAACQ_GW_POINT3_0_OFF
			  + islice*FRZ_RDB_SLICE_INFO_ENTRY_SIZE 
			  + i*4]);
    }
  CalculateLocation(p1, p2, p3, rotation, transpose);
  p1[2] += s_i_offset;
  p2[2] += s_i_offset;
  p3[2] += s_i_offset;

}

void
CheckHeaderInfo (const Filename input_file, Context *c)
{
  FILE *f;
  unsigned char h[FRZ_RDB_HEADER_SIZE_BYTES];
  Context cTmp;

  cTmp.verbosity= c->verbosity;
  ReadFileHeader(input_file, &cTmp, NULL);

  if (c->nph1 != cTmp.nph1 ||
      (c->nphmult != cTmp.nphmult && (c->nphmult > 1 || cTmp.nphmult> 1)) ||
      c->nimages != c->nph1 * c->nphmult ||
      c->nslices != cTmp.nslices ||
      c->sliceorder != cTmp.sliceorder ||
      c->ndatfr != cTmp.ndatfr ||
      c->gtype != cTmp.gtype ||
      c->opxres != cTmp.opxres ||
      c->ngap != cTmp.ngap ||
      c->concat != cTmp.concat ||
      c->slewrate != cTmp.slewrate ||
      c->fast_rec_lpf != cTmp.fast_rec_lpf ||
      c->mapdel != cTmp.mapdel ||
      c->samp_time != cTmp.samp_time ||
      c->ndat != cTmp.ndat ||
      c->npr != cTmp.npr ||
      c->chop != cTmp.chop ||
      c->risetime != cTmp.risetime ||
      c->densamp != cTmp.densamp ||
      c->fsgcm != cTmp.fsgcm ||
      c->opfov != cTmp.opfov ||
      c->ts != cTmp.ts ||
      c->gts != cTmp.gts ||
      c->ncoils != cTmp.ncoils ||
      c->coil_record_length != cTmp.coil_record_length)
    Abort("The headers of the input files are inconsistent.\n");

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

static int
IntBRdFloat32 (unsigned char *addr)
{
  return((int) Round(BRdFloat32(addr)));
}
