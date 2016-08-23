/*
 *	worker_utils.c
 *
 *    Worker routines for the parallel master/worker configuration of
 *    the "spiral" program.  All of the computationally-intensive
 *    routines are included in this file.  For more information about
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
#include <math.h>
#include <unistd.h>
#ifdef AFS
#include <sys/types.h>
#include <sys/wait.h>
#include <errno.h>
#endif
#include "mri.h"
#include "fmri.h"
#include "par.h"
#include "bio.h"
#include "array.h"
#include "acct.h"
#include "stdcrg.h"
#include "frozen_header_info.h"
#include "spiral.h"

#define NWEIGHTS	1024	/* # of times kaiser function is sampled */
#define SEGSIZE         3000    /* us */ 

static char rcsid[] = "$Id: worker_utils.c,v 1.11 2005/07/07 20:09:00 welling Exp $";

/* GLOBAL VARIABLES FOR WORKER_UTILS */
static float ***prd = NULL;  /* [2][npr][ndat]	      projection data (input)*/
static float ***grim = NULL; /* [2][os_res][os_res]  gridded image data */
static float **ws = NULL;    /* [os_res][os_res]     weighting array  */
static short int *in_buf = NULL; /* [2*ndat]input buffer for projection data */
static float weight[NWEIGHTS+1]; /* array which holds a sampled kaiser 
                                    function -- it is used
	 			    so we can do table-lookup 
                                    instead of computing the kaiser
				    function over and over */
static float maxk;	      /* maximum magnitude in the t2k array */
static float ***im = NULL;    /* [2][os_res][os_res]  image data */
static short int **mag;	      /* [os_res][os_res] magnitude of gridded data */
static float **outim = NULL;  /* [res][res]	output image data (float32) */
static short int **sh_outim = NULL; /* [res][res] output image data (int16) */
static float **totalim = NULL;/* [res][res] combined coil output image data */
static short int **mask = NULL;	    /* [res][res] circular mask for 
                                       cropping output */
static short int **imphase = NULL;  /* [res][res] image phase */
static float *wc = NULL;            /* [res] weighting coefficients */
static float **refim_out = NULL;    /* [res][res] reference image */
static float **refim_mag = NULL;    /* [res][res] reference image magnitude */
static float **fbuf = NULL;	    /* [2][os_res] FIX what? buffer */
static MRI_Dataset *ods = NULL;	    /* output datset */

/* FORWARD DECLARATIONS */
static void AllocWorkerMemory ();
static void DeallocWorkerMemory ();
static void WorkerPrecomputation ();
static void ComputeMulticoilReference ();
static void ComputeReference (int coil_num);
static void ComputeMulticoilImage ();
static void ComputeImage (int coil_num);
static long CalcOffset( int coil_num, int slice_sequence_number );
static void LoadProjections (int coil_num);
static void LoadMRIProjections (int coil_num);
static void Refocus ();
static void FixViews ();
static void PerformGridding (int coil_num);
static void FermiFilter1 ();
static void FermiFilter2 ();
static void FlipGrim ();
static void WriteRaw (int coil_num);
static void WriteSampleInfo (float **sampim);
static void WriteMagnitude (int coil_num);
static void WriteSamples (int coil_num);
static void MakePhase1Filter ();
static void MakePhase2Filter ();
static void WriteReference ();
static void TransformImage (int slice_num, int coil_num, int image_num);
static void BesselKaiserCorrectImage ();
static void WindowImageData (int middle_sample, int sample_window);
static void CopyImageData ();
static void ApplyFourierTransform ();
static void MakeFinalImage (float ***finalim, float ***im, int segment_num);
static void CopyFinalImage (float ***im, float ***finalim);
static void WriteSingleCoilImage (int coil_num);
static void WriteImage ();
static void WriteImagePhase (int coil_num);
static float kaiser (float l, float b, float u);
static float bessi0 (float x);
static void CheckAvailableInfo();
#ifdef AFS
static int MySystem( const char* command );
static int CheckForAFS( const char* fname );
static void FlushAFS( const char* fname );
#endif

void CheckAvailableInfo()
{
  /* Because we are trying to process multiple different Pfile
   * formats, we may get into a situation where some needed
   * piece of information is unavailable.  This routine is
   * just a collection of checks for this sort of situation.
   */
  if ((c.lin_cor || c.lin_map || c.hc_cor || c.gen_map)
      && !c.mapdel_set) {
    Abort("Attempted to generate or use a reference map without knowing map delay!\n");
  }
}

void DeallocWorkerMemory()
{
  Free3DFloatArray(prd);
  Free3DFloatArray(grim);
  Free2DFloatArray(ws);
  Free3DFloatArray(im);
  Free2DFloatArray(outim);
  Free2DShortArray(sh_outim);
  Free2DFloatArray(totalim);
  Free2DShortArray(mask);
  Free2DShortArray(imphase);
  if (wc != NULL)
    {
      free(wc);
      wc = NULL;
    }
  Free2DFloatArray(refim_out);
  Free2DFloatArray(refim_mag);
  Free2DFloatArray(fbuf);
  if (in_buf != NULL)
    {
      free(in_buf);
      in_buf = NULL;
    }
  if (r.slice_sampim != NULL)
    {
      Free2DFloatArray(r.slice_sampim);
      r.slice_sampim = NULL;
    }
}

void AllocWorkerMemory()
{
  /* reallocate the worker arrays */
  DeallocWorkerMemory();

  prd = Alloc3DFloatArray(2, c.npr, c.ndat);
  grim = Alloc3DFloatArray(2, c.os_res, c.os_res);
  ws = Alloc2DFloatArray(c.os_res, c.os_res);
  in_buf = (short int *) malloc(2*sizeof(short int)*c.ndat);
  im = Alloc3DFloatArray(2, c.os_res, c.os_res);
  outim = Alloc2DFloatArray(c.res, c.res);
  sh_outim = Alloc2DShortArray(c.res, c.res);
  totalim = Alloc2DFloatArray(c.res, c.res);
  mask = Alloc2DShortArray(c.res, c.res);
  imphase = Alloc2DShortArray(c.res, c.res);
  wc = (float *) malloc(c.res*sizeof(float));
  refim_out = Alloc2DFloatArray(c.res, c.res);
  refim_mag = Alloc2DFloatArray(c.res, c.res);
  fbuf = Alloc2DFloatArray(2, c.os_res); 

  /* allocate the sampim array in the Result */
  if (c.hc_cor)
    r.slice_sampim = Alloc2DFloatArray(c.os_res, c.os_res);

}

void WorkerShutdown()
{
  if (ods != NULL) mri_close_dataset(ods);
  DeallocWorkerMemory();
}

void WorkerInitialize()
{
  static int old_res = 0, old_os_res = 0;
  Filename ods_name;

  CheckAvailableInfo();

  if (c.res != old_res ||
      c.os_res != old_os_res)
    {
      /* reopen the output dataset */
      if (ods != NULL)
	{
	  Acct(WRITECLOSE);
	  mri_close_dataset(ods);
	}
      Acct(WRITEOPEN);
      ConstructFilename(ods_name, c.output_directory, c.output_name);
      /* open the output raw dataset */
      ods = mri_open_dataset(ods_name, MRI_MODIFY_DATA);
      Acct(PROCESSING);

      AllocWorkerMemory();
      WorkerPrecomputation();

      old_res = c.res;
      old_os_res = c.os_res;
    }
}

void WorkerPrecomputation()
{
  int i, j;
  int hr;
  float wctmp, wcmax;

  /* precompute the wc array that is used for kaiser-bessel correction */
  hr = c.res / 2;
  wcmax = sinh(sqrt(c.gridb*c.gridb))/sqrt(c.gridb*c.gridb);
  for (i = 0; i < c.res; i++)
    {
      wctmp = 
	PI*PI*c.grid_len*c.grid_len*(i-hr)*(i-hr)/(c.os_res*c.os_res/4) 
	- c.gridb*c.gridb;
      if(wctmp == 0.0)
	wc[i] = wcmax;
      else if (wctmp < 0.0)
	wc[i] = sqrt(-wctmp)/sinh(sqrt(-wctmp))*wcmax;
      else
	wc[i] = sqrt(wctmp)/sin(sqrt(wctmp))*wcmax;
    }
  /* precompute the mask array that is used for kaiser-bessel correction */
  for (i = 0; i < c.res; ++i)
    for (j = 0; j < c.res; ++j)
      mask[i][j] = (hypot((float)(i-hr),(float)(j-hr)) < .51*c.res);
  
  /* construct the convolution function table --- kaiser-bessel */
  for (i = 0; i < NWEIGHTS+1; i++)
    weight[i] = kaiser((float)NWEIGHTS, c.gridb, (float)(i-NWEIGHTS/2));
}

void
GenerateReference ()
{
  r.recon_failed= 0; /* set non-zero later if necessary */
  r.slice_num= t.slice_num;
  r.image_num= -1;
  for (t.image_num = 0; t.image_num < 2; ++t.image_num)
    if (c.ncoils == 1)
      ComputeReference(0);
    else
      ComputeMulticoilReference();
  WriteReference();
}

void
GenerateImage ()
{
  r.slice_num= t.slice_num;
  r.image_num= t.overall_image_num;
  r.recon_failed= 0; /* set non-zero later if necessary */

  /* Deal with the hopefully-rare case where the reference dataset is
   * unavailable for this particular slice.
   */
  if ((c.lin_cor || c.hc_cor) && c.ref_missing[t.slice_num]) {
    r.recon_failed= 1;
    return;
  }

  if (c.ncoils == 1)
    ComputeImage(0);
  else
    ComputeMulticoilImage();
  WriteImage();
}

static void
ComputeMulticoilReference ()
{
  int i, j, k;
  float **total_out;
  float **total_mag;

  total_out = Alloc2DFloatArray(c.res, c.res);
  total_mag = Alloc2DFloatArray(c.res, c.res);

  for (i = 0; i < c.res; i++)
    for (j = 0; j < c.res; j++)
      {
	refim_out[i][j] = 0.0;
	refim_mag[i][j] = 0.0;
      }

  for (k = 0; k < c.ncoils; k++)
    {
      ComputeReference(k);
      for (i = 0; i < c.res; ++i)
	for (j = 0; j < c.res; ++j)
	  {
	    total_out[i][j] += refim_out[i][j] * refim_mag[i][j];
	    total_mag[i][j] += refim_mag[i][j];
	  }
    }
  for (i = 0; i < c.res; i++)
    for (j = 0; j < c.res; j++)
      refim_out[i][j] = total_out[i][j] / (total_mag[i][j] + 0.001);

  Free2DFloatArray(total_out);
  Free2DFloatArray(total_mag);
}

static void
ComputeReference (int coil_num)
{
  if (c.verbosity >= VERBOSITY_TIMESTEP)
    Report("Co %2d, Sl %2d, Im %2d: Ld.",
	   coil_num, t.slice_num, t.image_num+1);
  if (c.input_is_mri) {
    LoadMRIProjections(coil_num);
  }
  else
    LoadProjections(coil_num);
  if (c.write_samples)
    WriteSamples(coil_num);
  Refocus();
  FixViews();
  if (c.verbosity >= VERBOSITY_TIMESTEP)
    Report("Rm.");
  PerformGridding(coil_num);
  if (c.verbosity >= VERBOSITY_TIMESTEP)
    Report("Fl.");
  FermiFilter1();
  FermiFilter2();
  if (c.reverse)
    FlipGrim(); /* flip the raw gridded image corner-for-corner */
  if (c.write_raw)
    WriteRaw(coil_num);
  if (c.write_mag)
    WriteMagnitude(coil_num);
  TransformImage(t.slice_num, coil_num, t.image_num);
  if (t.image_num == 0)
    MakePhase1Filter();
  else
    MakePhase2Filter();
  if (c.verbosity >= VERBOSITY_TIMESTEP)
    Report("\n");
}

static void
ComputeMulticoilImage ()
{
  float imt, imc, imm, imx;
  int i, j, k;
  int hr;

  for (i = 0; i < c.res; ++i)
    for (j = 0; j < c.res; ++j)
      totalim[i][j] = 0.0;

  hr = c.res / 2;
  for (k = 0; k < c.ncoils; ++k)
    {
      ComputeImage(k);
      if (r.recon_failed) return;

      if (c.do_recon) {
	if (c.all_coils)
	  {
	    WriteSingleCoilImage(k);
	    if (c.write_phase)
	      WriteImagePhase(k);
	  }
	
	/* combine the image in outim with that in totalim */
	imx = 0.0;
	for (i = 0; i < c.res; ++i)
	  for (j = 0; j < c.res; ++j)
	    if (mask[i][j])
	      totalim[i][j] += outim[i][j] * outim[i][j];
      }
    }

  if (c.do_recon)
    for (i = 0; i < c.res; ++i)
      for (j = 0; j < c.res; ++j)
	outim[i][j] = sqrt(totalim[i][j]);
}

static void
ComputeImage (int coil_num)
{
  if (c.verbosity >= VERBOSITY_TIMESTEP)
    Report("Co %2d, Sl %2d, Im %2d: Ld.", coil_num, t.slice_num, t.image_num+1);
  if (c.input_is_mri)    
    LoadMRIProjections(coil_num);
  else
    LoadProjections(coil_num);
  if (c.write_samples)
    WriteSamples(coil_num);
  if (c.do_recon) {
    if (!r.recon_failed) Refocus();
    if (!r.recon_failed) FixViews();
    if (c.verbosity >= VERBOSITY_TIMESTEP)
      Report("Rm.");
    if (!r.recon_failed) PerformGridding(coil_num);
    if (c.verbosity >= VERBOSITY_TIMESTEP)
      Report("Fl.");
    if (!r.recon_failed) FermiFilter1();
    if (!r.recon_failed) FermiFilter2();
    if (c.reverse)
      FlipGrim(); /* flip the raw gridded image corner-for-corner */
    if (c.write_raw)
      WriteRaw(coil_num);
    if (c.write_mag)
      WriteMagnitude(coil_num);
    if (!r.recon_failed) TransformImage (t.slice_num, coil_num, t.image_num);
    if (!r.recon_failed) BesselKaiserCorrectImage ();
    if (c.write_phase)
      WriteImagePhase(coil_num);
  }
  if (c.verbosity >= VERBOSITY_TIMESTEP)
    Report("\n");
}

static long
CalcOffset(int coil_num, int slice_sequence_number)
{
  /* This routine calculates file offsets for slices, based on 
   * context and task info.  This process is a lot more obscure than
   * common sense would suggest it should be.
   */
  long offset;
  long coil_offset;
  long nproj;
  int nodd;

#ifdef never
  offset = FRZ_RDB_HEADER_SIZE_BYTES +
    coil_num * c.coil_record_length +
    slice_sequence_number * c.slice_record_length +
    (t.image_num / c.nphmult + 1) * c.baseline_length +
    t.image_num * c.npr * c.ndat * 4;
  if (((c.nphmult * c.npr) % 2) == 1)
    offset += (t.image_num / c.nphmult) * c.ndat * 4;
#endif

  nodd= (c.nphmult*c.npr*(c.concat+1))%2;
  coil_offset = c.coil_record_length * coil_num;
  nproj = slice_sequence_number * ((c.concat+1)*c.nphmult*c.npr + 1 + nodd)
    + (c.concat+1)*t.image_num*c.npr + 1;

  offset = FRZ_RDB_HEADER_SIZE_BYTES + nproj*c.ndatfr*4 + coil_offset;

  return offset; 
}

static void
LoadProjections (int coil_num)
{
  int n;
  int j;
  int chopper = 1;
  FILE *f;
  int slice_sequence_number;
  long offset;
  char com[sizeof(Filename)+128];
  int nodd;

  Acct(READOPEN);
#ifdef AFS
  if (CheckForAFS(".")) FlushAFS(".");
#endif
  if ((f = fopen(t.filename, "r")) == NULL)
    Abort("Could not open input file %s\n", t.filename);

  Acct(READING);
  /* all odd slices are located before all even slices
     in the P file (assuming a numbering scheme that
     begins with 1); hence we need to compute a sequence
     number in order to pull the right slice out of the
     file; note that our internal slice numbers begin
     with 0) */
  if (c.sliceorder) {
    slice_sequence_number = t.slice_num;
  }
  else {
    if (t.slice_num % 2 == 0)
      slice_sequence_number = t.slice_num / 2;
    else
      slice_sequence_number = (c.nslices + t.slice_num) / 2;
  }

  offset= CalcOffset(coil_num, slice_sequence_number); 
  if (fseek(f, offset, SEEK_SET) != 0)
    Abort("Could not load projections\n");

  bio_error = FALSE;
  bio_big_endian_input = c.big_endian_input;
  for (n = 0; n < c.npr; ++n)
    {
      FRdInt16Array(f, in_buf, 2*c.ndat);
      if (bio_error)
	Abort("Could not load projections\n");
      if (c.reverse) {
	/* reverse sample order on the fly */
	for (j = 0; j < c.ndat; ++j)
	  {
	    prd[1][n][c.ndat-j-1] = in_buf[2*j] * chopper;
	    prd[0][n][c.ndat-j-1] = in_buf[2*j+1] * chopper;
	  }
      }
      else {
	for (j = 0; j < c.ndat; ++j)
	  {
	    prd[0][n][j] = in_buf[2*j] * chopper;
	    prd[1][n][j] = in_buf[2*j+1] * chopper;
	  }
      }
      chopper *= c.chop;
    }

  Acct(READCLOSE);
  fclose(f);
  Acct(PROCESSING);
}

static void
LoadMRIProjections (int coil_num)
{
  int n;
  int j;
  long overall_offset;
  char com[sizeof(Filename)+128];
  short* buf;
  short* runner;
  MRI_Dataset* ds;

  /* Unfortunately the data order in the Pgh MRI file doesn't match
   * that in the prd[][][] array, so we must read into a buffer and
   * reorder.
   */
  if (!(buf= (short*)malloc(2*c.ndat*sizeof(short))))
    Abort("Unable to allocate %d bytes!\n",2*c.ndat*sizeof(short));

  Acct(READOPEN);

  if (!(ds= mri_open_dataset(t.filename, MRI_READ)))
    Abort("Could not open input file %s\n", t.filename);

  /* Note that the dataset was checked for structure when all the context
   * values were read from it.  Thus we don't have to check the presence
   * of the expected keys.
   */

  /* Note if this slice is missing */
  if (mri_has(ds,"missing") 
      && !strcmp(mri_get_string(ds,"missing"),"[chunk]")) {
    unsigned char this_slice_missing;
    long long offset= 0;
    if (!strcmp(mri_get_string(ds,"missing.dimensions"),"zt")) {
      offset= t.image_num*c.nslices + t.slice_num;
    }
    else if (!strcmp(mri_get_string(ds,"missing.dimensions"),"tz")) {
      offset= t.slice_num*c.nimages + t.image_num;
    }
    else Abort("Missing info order %s is unsupported!\n",
	       mri_get_string(ds,"missing.dimensions"));
    mri_read_chunk(ds,"missing", 1, offset, MRI_UNSIGNED_CHAR,
		   &this_slice_missing);
    if (this_slice_missing != 0) r.recon_failed= 1;
  }

  Acct(READING);

  /* read the chunk */
  for (n=0; n<c.npr; n++) {
    /* can't use t.overall_image_num here because it's always 0 when
     * generating references.
     */
    overall_offset = 
      ((((((((t.file_index*c.nimages + t.image_num) 
	     * c.nslices)
	    +t.slice_num)
	   *c.ncoils) 
	  + coil_num)
	 *c.npr)
	+ n)*c.ndat*2);
    mri_read_chunk(ds, "samples", 2*c.ndat, overall_offset,
		  MRI_SHORT, buf);
    runner= buf;
    if (c.reverse) {
      /* reverse sample order on the fly */
      for (j=0; j<c.ndat; j++) {
	prd[1][n][c.ndat-j-1]= *runner++;
	prd[0][n][c.ndat-j-1]= *runner++;
      }      
    }
    else {
      for (j=0; j<c.ndat; j++) {
	prd[0][n][j]= *runner++;
	prd[1][n][j]= *runner++;
      }
    }
  }

  Acct(READCLOSE);
  mri_close_dataset(ds);
  Acct(PROCESSING);

  free(buf);
}

static void
Refocus ()
{
  int i, j, j2;
  double theta;
  float tr, ti, cs, sn, pt2;

  pt2 = (c.fast_rec_lpf*1e3*c.ts);
  for (j=0; j<c.ndatfr; j++)
    {
      theta = (-2.0*PI*(pt2 + c.ph_twist/c.ndat) - c.refl[t.slice_num][0])*j;
      cs = cos(theta);
      sn = sin(theta);
      for (i=0; i<c.npr; i++)
	{
	  tr = prd[0][i][j];
	  ti = prd[1][i][j];
	  prd[0][i][j] = tr*cs - ti*sn;
	  prd[1][i][j] = tr*sn + ti*cs;
	}
    }
  if (c.concat)
    {
      for (j=0; j<c.ndatfr; j++) {
	j2 = j+c.ndatfr+c.ngap;
	theta = (-2.0*PI*(pt2 + c.ph_twist/c.ndat) 
		 - c.refl[t.slice_num][0])*j2;
	cs = cos(theta);
	sn = sin(theta);
	for (i=0; i<c.npr; i++) {
	  tr = prd[0][i][j2]; 
	  ti = prd[1][i][j2];
	  prd[0][i][j2] = tr*cs - ti*sn;
	  prd[1][i][j2] = tr*sn + ti*cs;
	}
      }
      for (i=0; i<c.npr; i++) {
	for (j=0; j<c.ngap; j++) {
	  prd[0][i][c.ndatfr+j] = 
	    0.5*(prd[0][i][c.ndatfr-1] + prd[0][i][c.ndatfr+c.ngap]);
	  prd[1][i][c.ndatfr+j] = 
	    0.5*(prd[1][i][c.ndatfr-1] + prd[1][i][c.ndatfr+c.ngap]); 
	}
      }
    }
}

static void
FixViews ()
{
  int i, j;
  float p0r, p0i, prr, pri, cc, dd, tmp;
  FILE *warn, *f;
  float cs, sn, theta;
  long offset;
  int slice_sequence_number;
  float *tmp_p0, *tmp_p1, *tmp_r0, *tmp_r1;
  short *tmp_s;

  /* remove mag and phase variations (e.g. Hu's method) */
  p0r = p0i = 0.0;
  for (i=0; i<c.npr; i++)
    for (j = 0; j <= c.samp_delay; j++)
      { 
	p0r += prd[0][i][j];
	p0i += prd[1][i][j];
      }
  p0r /= c.npr;
  p0i /= c.npr;

  for (i = 0; i < c.npr; i++)
    {
      prr = pri = 0.0;
      for (j = 0; j <= c.samp_delay; j++)
	{ 
	  prr += prd[0][i][j];
	  pri += prd[1][i][j];
	}

      /* phase correction angle */
      tmp = atan2(p0i,p0r) - atan2(pri,prr); 
      if (tmp >  PI)
	tmp -= 2.0*PI;
      if (tmp < -PI)
	tmp += 2.0*PI;
      cc = cos(tmp * c.ph_factor);
      dd = sin(tmp * c.ph_factor);

      /* magnitude correction term   */
      if (prr == 0.0 && pri == 0.0) {
	Error("Signal is 0 at k-space origin (image %d slice %d)!\n",
	      t.overall_image_num, t.slice_num);
	r.recon_failed= 1;
	return;
      }
      else
	tmp = hypot(p0r, p0i) / hypot(prr, pri);
      cc *= tmp*c.mag_factor + (1.0 - c.mag_factor);
      dd *= tmp*c.mag_factor + (1.0 - c.mag_factor);
      for (j = 0; j < c.ndat; j++)
	{ 
	  prr = prd[0][i][j];
	  pri = prd[1][i][j];
	  tmp = cc*prr - dd*pri;
	  prd[1][i][j] = dd*prr + cc*pri;
	  prd[0][i][j] = tmp;
	}
    }
}

static void FlipGrim()
{
  /* This routine flips the raw gridded image corner-for-corner */
  float*** tmp_grim= Alloc3DFloatArray(2,c.os_res,c.os_res);
  int i;
  int j;
  
  for (i=0; i<c.os_res; i++)
    for (j=0; j<c.os_res; j++) {
      tmp_grim[0][i][j]= grim[0][c.os_res-i-1][c.os_res-j-1];
      tmp_grim[1][i][j]= grim[1][c.os_res-i-1][c.os_res-j-1];
    }
  
  for (i=0; i<c.os_res; i++)
    for (j=0; j<c.os_res; j++) {
      grim[0][i][j]= tmp_grim[0][i][j];
      grim[1][i][j]= tmp_grim[1][i][j];
    }
  
  Free3DFloatArray(tmp_grim);
}

static void
PerformGridding (int coil_num)
{
  int i, j;
  int lx, ly;
  float kx, ky;
  float w;
  float pr, pi;
  float mkr;
  float dkx,dky;
  float dwin;
  float w2;
  float wx;
  float tmpd;
  float rot1, rot2;
  float rotfact;
  int maxj;
  int min_lx, max_lx;
  int min_ly, max_ly;
  int count;

  for (i=0; i<c.os_res; i++)
    for (j=0; j<c.os_res; j++)
      {
	grim[0][i][j] = 0.0;
	grim[1][i][j] = 0.0;
	ws[i][j] = 0.0;
      }
  if (t.overall_image_num == 0 && coil_num == 0 && c.hc_cor)
    {
      r.sampim_valid = 1;
      for (i=0; i<c.os_res; i++)
	for (j=0; j<c.os_res; j++)
	  r.slice_sampim[i][j] = 0.0;
    }
  else
    r.sampim_valid= 0;

  maxk = 0.0;
  dwin = c.wind;
  rotfact = 2.0*PI / c.res / c.over_samp;
  maxj = c.ndat - c.samp_delay - 2;
  for (i = 0; i < c.npr; i++)
    for (j=0; j <= maxj; j++)
      { 
	w2 = c.kdens[j];
	kx = c.t2k[0][i][j];
	ky = c.t2k[1][i][j];
	if ( (mkr = hypot(kx,ky)) > maxk)
	  maxk = mkr;
	pr = prd[0][i][j+c.samp_delay];
	pi = prd[1][i][j+c.samp_delay];
	dkx = (c.factxx*kx + c.factxy*ky)*c.over_samp;
	dky = (c.factyx*kx + c.factyy*ky)*c.over_samp;
	if (c.reg_file[0] != '\0')
	  { 
	    rot1 = cos(rotfact*dkx*(c.pix_shifth*c.res - t.reg_xs + c.lr_shift));
	    rot2 = -sin(rotfact*dkx*(c.pix_shifth*c.res - t.reg_xs + c.lr_shift));
	    tmpd = pr*rot1 + pi*rot2;
	    pi = pi*rot1 - pr*rot2;
	    pr = tmpd;

	    rot1 = cos(rotfact*dky*(c.pix_shiftv*c.res - t.reg_ys + c.tb_shift));
	    rot2 = -sin(rotfact*dky*(c.pix_shiftv*c.res - t.reg_ys + c.tb_shift));
	    tmpd = pr*rot1 + pi*rot2;
	    pi = pi*rot1 - pr*rot2;
	    pr = tmpd;

	    rot1 = cos(-PI/180.0 * t.reg_rot);
	    rot2 = -sin(-PI/180.0 * t.reg_rot);
	    tmpd = dkx*rot1 + dky*rot2;
	    dky = dky*rot1 - dkx*rot2;
	    dkx = tmpd; 
	  }
	else if (c.lr_shift != 0 || c.tb_shift != 0 || c.loc_shift)
	  { 
	    rot1 = cos(rotfact*dkx*(c.pix_shifth*c.res + c.lr_shift));
	    rot2 = -sin(rotfact*dkx*(c.pix_shifth*c.res + c.lr_shift));
	    tmpd = pr*rot1 + pi*rot2;
	    pi = pi*rot1 - pr*rot2;
	    pr = tmpd;

	    rot1 = cos(rotfact*dky*(c.pix_shiftv*c.res + c.tb_shift));
	    rot2 = -sin(rotfact*dky*(c.pix_shiftv*c.res + c.tb_shift));
	    tmpd = pr*rot1 + pi*rot2;
	    pi = pi*rot1 - pr*rot2;
	    pr = tmpd;
	  }
	dkx += (c.res/2 - c.refl[t.slice_num][2]*j)*c.over_samp; 
	dky += (c.res/2 - c.refl[t.slice_num][1]*j)*c.over_samp;

	/* compute loop bounds to avoid checking within loop (borrowed
	   from Mark Hahn's newgsp.C) */
	min_lx = ceil(dkx - dwin);
	if (min_lx < 0)
	  min_lx = 0;

	max_lx = floor(dkx + dwin);
	if (max_lx >= c.os_res)
	  max_lx = c.os_res - 1;

	for (lx = min_lx; lx <= max_lx; lx++)
	  {
	    wx = weight[ Round((((dkx-lx)/dwin)*NWEIGHTS/2 )) + NWEIGHTS/2] *w2;

	    min_ly = ceil(dky - dwin);
	    if (min_ly < 0)
	      min_ly = 0;

	    max_ly = floor(dky + dwin);
	    if (max_ly >= c.os_res)
	      max_ly = c.os_res - 1;

	    for (ly = min_ly; ly <= max_ly; ly++)
	      {
		w = wx*weight[ Round((((dky-ly)/dwin)*NWEIGHTS/2 )) + NWEIGHTS/2];
		if (t.overall_image_num == 0 && coil_num == 0 && c.hc_cor)
		  r.slice_sampim[lx][ly] += w * j;
		grim[0][lx][ly] += pr * w;
		grim[1][lx][ly] += pi * w;
		ws[lx][ly] += w;
	      }
	  }
      }

  if (t.overall_image_num == 0 && coil_num == 0 && c.hc_cor)
    {
      for (i=0; i<c.os_res; i++)
	for (j=0; j<c.os_res; j++)
	  r.slice_sampim[i][j] /= (0.001+ws[i][j]);
      if (c.write_raw)
	WriteSampleInfo(r.slice_sampim);
      /* copy into the Context so that we have correct
	 values for the reconstruction of this slice */
      memcpy(&c.sampim[r.slice_num][0][0],
	     &r.slice_sampim[0][0],
	     c.os_res * c.os_res * sizeof(float));
    }
}

static void
FermiFilter1 ()
{
  float ww;
  float rad;
  int i, j;
  int offs;

  rad = (c.res / 2) * c.over_samp;
  if (maxk*c.over_samp/c.zoom < rad)
    rad = maxk*c.over_samp/c.zoom;
  offs = (c.res / 2) * c.over_samp;
  for (i = 0; i < c.os_res; ++i)
    for (j = 0; j < c.os_res; ++j)
      {
	ww = 1.0 / (1.0 + exp((hypot((float)(offs-i), (float)(offs-j)) - rad) *
			      6.4 / rad));
	ws[i][j] = ww / (0.001 + ws[i][j]);
      }
}

static void
FermiFilter2 ()
{
  int i, j;

  for (i = 0; i < c.os_res; ++i)
    for (j = 0; j < c.os_res; ++j)
      {
	grim[0][i][j] *= ws[i][j];
	grim[1][i][j] *= ws[i][j];
      }
}

static void
WriteRaw (int coil_num)
{
  int overall_offset;

  overall_offset = ((coil_num * c.total_images + t.overall_image_num) *
		    c.nslices + t.slice_num) * 2 * c.os_res * c.os_res;
  
  Acct(WRITING);
  mri_set_chunk(ods, "raw", 2*c.os_res*c.os_res, overall_offset,
		MRI_FLOAT, &grim[0][0][0]);
  Acct(PROCESSING);
}

static void
WriteSampleInfo (float **sampim)
{
  int overall_offset;

  overall_offset = t.slice_num * c.os_res * c.os_res;
  
  Acct(WRITING);
  mri_set_chunk(ods, "sampim", c.os_res*c.os_res, overall_offset,
		MRI_FLOAT, &sampim[0][0]);
  Acct(PROCESSING);
}

static void
WriteSamples (int coil_num)
{
  int n, j;
  int overall_offset;
  short* buf;
  short* runner;

  if (!(buf= (short*)malloc(2*c.ndat*sizeof(short))))
    Abort("Unable to allocate %d bytes!\n",2*c.ndat*sizeof(short));

  /* write the chunk */
  for (n=0; n<c.npr; n++) {
    runner= buf;
    if (c.reverse) {
      /* Samples were reversed on input, so we must un-
       * reverse them on output.
       */
      for (j=0; j<c.ndat; j++) {
	*runner++= (short)rint(prd[1][n][c.ndat-j-1]);
	*runner++= (short)rint(prd[0][n][c.ndat-j-1]);
      }      
    }
    else {
      for (j=0; j<c.ndat; j++) {
	*runner++= (short)rint(prd[0][n][j]);
	*runner++= (short)rint(prd[1][n][j]);
      }
    }
    /* can't use t.overall_image_num here because it's always 0 when
     * generating references.
     */
    overall_offset = 
      ((((((((t.file_index*c.nimages + t.image_num) 
	     * c.nslices)
	    +t.slice_num)
	   *c.ncoils) 
	  + coil_num)
	 *c.npr)
	+ n)*c.ndat*2);
    Acct(WRITING);
    mri_set_chunk(ods, "samples", 2*c.ndat, overall_offset,
		  MRI_SHORT, buf);
    Acct(PROCESSING);
  }

  free(buf);
}

static void
WriteMagnitude (int coil_num)
{
  int i, j;
  float imm, imx;
  int overall_offset;
  short **mag;

  /* construct the output buffer from the magnitudes of the values
     in the grim array */
  mag = Alloc2DShortArray(c.os_res, c.os_res);
  imx = 0.0;
  for (i = 0; i < c.os_res; i++)
    for (j = 0; j < c.os_res; j++)
      {
	imm = hypot(grim[0][j][i], grim[1][j][i]);
	if (imm > imx)
	  imx = imm;
	mag[i][j] = imm;
      }
  if (c.verbosity >= VERBOSITY_TIMESTEP)
    Report("max=%3d", (int) imx);

  /* write the chunk */
  overall_offset = ((coil_num * c.total_images + t.overall_image_num) *
		    c.nslices + t.slice_num) * c.os_res * c.os_res;
  Acct(WRITING);
  mri_set_chunk(ods, "magnitude", c.os_res*c.os_res, overall_offset,
		MRI_SHORT, &mag[0][0]);
  Acct(PROCESSING);

  Free2DShortArray(mag);
}

static void
MakePhase1Filter ()
{
  int i, j, k, offs;
  int hfsize;	/* half-width of filter-size */
  float sumr,sumi;
  int j_start, j_end;
  int k_start, k_end;

  hfsize = (c.filter_sz - 1)/2;
  offs = (c.over_samp-1)*c.res/2;

  for (i = 0; i < c.res; i++)
  {
    j_start = MAX(0, offs - hfsize);
    j_end = MIN(c.res + offs + hfsize, c.os_res-1);
    for (j = j_start; j <= j_end; j++)
      {
	fbuf[0][j] = 0.0;
	fbuf[1][j] = 0.0;
	k_start = MAX(0, i + offs - hfsize);
	k_end = MIN(i + offs + hfsize, c.os_res-1);
	for (k = k_start; k <= k_end; k++)
	  {
	    fbuf[0][j] += im[0][k][j];
	    fbuf[1][j] += im[1][k][j];
	  }
      }
    for (j = 0; j < c.res; j++)
      {
	sumr = sumi = 0.;
	k_start = MAX(0, j + offs - hfsize);
	k_end = MIN(j + offs + hfsize, c.os_res-1);
	for (k = k_start; k <= k_end; k++)
	  {
	    sumr += fbuf[0][k];
	    sumi += fbuf[1][k];
	  }
	refim_out[i][j] = atan2(sumi,sumr);
      }
  }
}

static void
MakePhase2Filter ()
{
  int i, j, k, offs;
  int hfsize;
  float sumr, sumi;
  float phaseoff;
  int j_start, j_end;
  int k_start, k_end;

  /* phaseoff is value to offset freq map if phase twist (ph_twist) != 0 */
  phaseoff = -2.0*PI * (c.mapdel / c.samp_time) * c.ph_twist / c.ndat;
  hfsize = (c.filter_sz - 1) / 2;
  offs = (c.over_samp - 1) * c.res/2;

  for (i = 0; i < c.res; i++)
    {
      j_start = MAX(0, offs - hfsize);
      j_end = MIN(c.res + offs + hfsize, c.os_res-1);
      for (j = j_start; j <= j_end; j++)
	{
	  fbuf[0][j] = 0.0;
	  fbuf[1][j] = 0.0;
	  k_start = MAX(0, i + offs - hfsize);
	  k_end = MIN(i + offs + hfsize, c.os_res-1);
	  for (k = k_start; k <= k_end; k++)
	    {
	      fbuf[0][j] += im[0][k][j];
	      fbuf[1][j] += im[1][k][j];
	    }
	}
      for (j = 0; j < c.res; j++)
	{
	  sumr = sumi = 0.0;
	  k_start = MAX(0, j + offs - hfsize);
	  k_end = MIN(j + offs + hfsize, c.os_res - 1);
	  for (k = k_start; k <= k_end; k++)
	    {
	      sumr += fbuf[0][k];
	      sumi += fbuf[1][k];
	    }
	  refim_mag[i][j] = hypot(sumi,sumr);
	  refim_out[i][j] -= (atan2(sumi,sumr) - phaseoff); 
	}
    }
}

static void
WriteReference ()
{
  int i, j;
  int hr;	/* half the resolution */
  int rm1;	/* resolution minus 1 */	
  float max;
  double tmpxr, tmpyr, tmpxi, tmpyi, tmpc, tmpcr, tmpci;
  float out[3];

  if (r.recon_failed) {
    unsigned char one= 1;

    /* Note the recon failure in the "missing" chunk, so that users of
     * this reference file can know to ignore those slices.
     */
    /* This sample is missing */
    Acct(WRITING);
    mri_set_chunk(ods, "missing", 1, t.slice_num,
		  MRI_UNSIGNED_CHAR, &one);
    Acct(PROCESSING);

    if (c.lin_map) {
      for (i=0; i<3; i++) out[i]= 0.0;
      
      Acct(WRITING);
      mri_set_chunk(ods, "linear", 3, 3*t.slice_num, MRI_FLOAT, out);
      Acct(PROCESSING);	
    }
    
    if (c.gen_map) {
	for (i = 0; i < c.res; i++)
	  for (j = 0; j < c.res; j++) {
	      refim_out[i][j]= 0.0;
	  }
	
	Acct(WRITING);
	mri_set_chunk(ods, "general", c.res*c.res, t.slice_num*c.res*c.res,
		      MRI_FLOAT, &refim_out[0][0]);
	Acct(PROCESSING);
    }
  }
  else {
    unsigned char zero= 0;

    /* This sample is not missing */
    Acct(WRITING);
    mri_set_chunk(ods, "missing", 1, t.slice_num,
		  MRI_UNSIGNED_CHAR, &zero);
    Acct(PROCESSING);
    

    if (c.lin_map)
      {
	hr = c.res / 2;
	rm1 = c.res - 1;
	
	max = 0.0;
	for (i = 0; i < c.res; i++)
	  for (j = 0; j < c.res; j++) 
	    if (hypot((float)(i-hr), (float)(j-hr)) < hr &&
		refim_mag[i][j] > max)
	      max = refim_mag[i][j];
	max /= 4.0;
	
	tmpxr = tmpyr = tmpxi = tmpyi = tmpcr = tmpci = 0.0;
	/* find linear terms */
	for (i = 0; i < rm1; i++) 
	  for (j = 0; j < rm1; j++) 
	    if (refim_mag[i][j] > max)
	      {
		tmpxr += refim_mag[i][j] * refim_mag[i][j+1]
		  * cos(refim_out[i][j]-refim_out[i][j+1]);
		tmpxi += refim_mag[i][j] * refim_mag[i][j+1]
		  * sin(refim_out[i][j]-refim_out[i][j+1]);
		tmpyr += refim_mag[i][j] * refim_mag[i+1][j]
		  * cos(refim_out[i][j]-refim_out[i+1][j]);
		tmpyi += refim_mag[i][j] * refim_mag[i+1][j]
		  * sin(refim_out[i][j]-refim_out[i+1][j]);
	      }
	
	tmpxr = -atan2(tmpxi,tmpxr);
	tmpyr = -atan2(tmpyi,tmpyr);
	/* uncomment to zero linear terms */
	/* tmpxr = tmpyr = 0.; */ 
	
	/* find constant term with linears removed */
	for (i = 0; i < rm1; i++) 
	  for (j = 0; j < rm1; j++) 
	    if (refim_mag[i][j] > max)
	      {
		tmpcr += refim_mag[i][j]
		  * cos(refim_out[i][j] - tmpxr*(j-hr) - tmpyr*(i-hr));
		tmpci += refim_mag[i][j]
		  * sin(refim_out[i][j] - tmpxr*(j-hr) - tmpyr*(i-hr));
	      }
	
	tmpc = atan2(tmpci,tmpcr);
	
	out[0] = tmpc * c.samp_time / c.mapdel;
	out[1] = tmpxr * c.res * c.samp_time / c.mapdel / 2.0 /PI;
	out[2] = tmpyr * c.res * c.samp_time / c.mapdel / 2.0 /PI;
	
	Acct(WRITING);
	mri_set_chunk(ods, "linear", 3, 3*t.slice_num, MRI_FLOAT, out);
	Acct(PROCESSING);
	
	for (i = 0; i < rm1; i++) 
	  for (j = 0; j < rm1; j++)
	    {
	      tmpcr = refim_mag[i][j]
		* cos(refim_out[i][j] - tmpxr*(j-hr) - tmpyr*(i-hr) - tmpc);
	      tmpci = refim_mag[i][j]
		* sin(refim_out[i][j] - tmpxr*(j-hr) - tmpyr*(i-hr) - tmpc);
	      refim_out[i][j] = atan2(tmpci,tmpcr);
	    }
      }
    
    if (c.gen_map)
      {
	for (i = 0; i < c.res; i++)
	  for (j = 0; j < c.res; j++)
	    {
	      if (refim_out[i][j] > PI)
		refim_out[i][j] -= 2*PI;
	      if (refim_out[i][j] < -PI)
		refim_out[i][j] += 2*PI;
	      if (refim_out[i][j] > 0.5*PI)
		refim_out[i][j] = 0.5*PI;
	      if (refim_out[i][j] < -0.5*PI)
		refim_out[i][j] = -0.5*PI;
	    }
	
	Acct(WRITING);
	mri_set_chunk(ods, "general", c.res*c.res, t.slice_num*c.res*c.res,
		      MRI_FLOAT, &refim_out[0][0]);
	Acct(PROCESSING);
      }
  }

  if (c.verbosity >= VERBOSITY_TIMESTEP)
    Report("Wr.");
}

static void
TransformImage (int slice_num,
		int coil_num,
		int image_num)
{
  int nsegs;
  int segment_num;
  int sample_window;
  int middle_sample;
  float ***finalim;

  if (c.verbosity >= VERBOSITY_TIMESTEP)
    Report("FT.");
  if (c.hc_cor)
    {
      /* segmented hc recon */
      nsegs = floor(1.0 * c.ndat * c.samp_time / SEGSIZE);
      sample_window = SEGSIZE / c.samp_time;
      finalim = Alloc3DFloatArray(2, c.os_res, c.os_res);
      for (segment_num = 0; segment_num <= nsegs; segment_num++)
	{
	  if (c.verbosity >= VERBOSITY_TIMESTEP)
	    Report("%d.", segment_num+1);
	  middle_sample = segment_num * sample_window;
	  WindowImageData(middle_sample, sample_window);
	  ApplyFourierTransform();
	  MakeFinalImage(finalim, im, segment_num);
	}
      CopyFinalImage(im, finalim);
      Free3DFloatArray(finalim);
    }
  else
    {
      CopyImageData();
      ApplyFourierTransform();
    }
}

static void
BesselKaiserCorrectImage ()
{
  int i, j;
  int offs;
  float imi, imr;
  int hr;

  hr = c.res / 2;
  offs = (c.over_samp-1)*hr;

  for (i = 0; i < c.res; i++)
    for (j = 0; j < c.res; j++)
      if (mask[i][j])
	{
	  imr = im[0][-j+c.res+offs-1][-i+c.res+offs-1]*wc[i]*wc[j];
	  imi = im[1][-j+c.res+offs-1][-i+c.res+offs-1]*wc[i]*wc[j];
	  outim[i][j] = sqrt(imi*imi+imr*imr);
	}
      else
	outim[i][j] = 0.0;
}

static void
WindowImageData (int middle_sample,
		 int sample_window)
{
  float start_sample, end_sample;
  float radius;
  float ww;
  int i, j;

  start_sample = middle_sample - sample_window;
  end_sample = middle_sample + sample_window;

  for (i = 0; i < c.os_res; ++i)
    for (j = 0; j < c.os_res; ++j)
      {
	radius = c.sampim[t.slice_num][i][j];
	if (radius > start_sample && radius < end_sample)
	  ww = 0.5 * (1.0 + cos(PI*(radius - middle_sample) / sample_window));
	else
	  ww = 0.0;
	im[0][i][j] = ww*grim[0][i][j];
	im[1][i][j] = ww*grim[1][i][j];
      }
}

static void
CopyImageData ()
{
  memcpy(&im[0][0][0], &grim[0][0][0], 2*c.os_res*c.os_res*sizeof(float));
}

static void
ApplyFourierTransform ()         /* formerly ft_image */
{
  int i, j;
  static FComplex** w= NULL;
  static int w_dim= 0;

  if (w_dim != c.os_res && w) FreeMatrix(w);

  if (!w) {
    w= Matrix(c.os_res, c.os_res, FComplex);
    w_dim= c.os_res;
  }

  for (i=0; i<c.os_res; i++)
    for (j=0; j<c.os_res; j++) {
      w[i][j].real= im[0][i][j];
      w[i][j].imag= im[1][i][j];
    }
  fft2d(w, c.os_res, c.os_res, 1, '2', 0, c.os_res-1);
  for (i=0; i<c.os_res; i++)
    for (j=0; j<c.os_res; j++) {
      im[0][i][j]= c.os_res*w[i][j].real;
      im[1][i][j]= c.os_res*w[i][j].imag;
    }

}

static void
MakeFinalImage (float ***finalim,
		float ***im,
		int segment_num)
{
  int i, j;
  int offs;
  float phase_factor;

  phase_factor = -segment_num;
  offs = (c.over_samp-1)*c.res/2;
  if (segment_num == 0)
    for (i=0; i<c.res; i++)
      {
	memcpy(&finalim[0][i][0],
	       &im[0][i+offs][offs],
	       c.res * sizeof(float));
	memcpy(&finalim[1][i][0],
	       &im[1][i+offs][offs],
	       c.res * sizeof(float));
      }
  else
    for (i = 0; i < c.res; i++)
      for (j = 0; j < c.res; j++)
	{
	  double theta= 
	    phase_factor*c.refim_in[t.slice_num][i][j]*SEGSIZE/c.mapdel;
	  finalim[0][i][j] += 
	    im[0][i+offs][j+offs]*cos(theta) -
	    im[1][i+offs][j+offs]*sin(theta);
	  finalim[1][i][j] += 
	    im[0][i+offs][j+offs]*sin(theta) +
	    im[1][i+offs][j+offs]*cos(theta);
	}
}

static void
CopyFinalImage (float ***im, float ***finalim)
{
  int i;
  int offs;

  offs = (c.over_samp-1)*c.res/2;
  for (i = 0; i < c.res; i++)
    {
      memcpy(&im[0][i+offs][offs],
	     &finalim[0][i][0],
	     c.res * sizeof(float));
      memcpy(&im[1][i+offs][offs],
	     &finalim[1][i][0],
	     c.res * sizeof(float));
    }
}

static void
WriteSingleCoilImage (int coil_num)
{
  char chunk_name[64];
  int overall_offset;
  float **f_output;
  short **s_output;
  int i, j;

  Acct(WRITING);
  if (c.verbosity >= VERBOSITY_TIMESTEP)
    Report("Wr.");
  sprintf(chunk_name, "coil%d", coil_num);
  overall_offset = (t.overall_image_num * c.nslices + (t.slice_num - c.start_slice)) * c.res * c.res;

  if (c.output_float)
    {
      f_output = Alloc2DFloatArray(c.res, c.res);
      for (i = 0; i < c.res; ++i)
	for (j = 0; j < c.res; ++j)
	  f_output[i][j] = c.output_scale * outim[i][j] / (c.os_res * c.os_res);
      mri_set_chunk(ods, chunk_name, c.res * c.res, overall_offset, MRI_FLOAT, &f_output[0][0]);
      Free2DFloatArray(f_output);
    }
  else
    {
      s_output = Alloc2DShortArray(c.res, c.res);
      for (i = 0; i < c.res; ++i)
	for (j = 0; j < c.res; ++j)
	  s_output[i][j] = c.output_scale * outim[i][j] / (c.os_res * c.os_res);
      mri_set_chunk(ods, chunk_name, c.res * c.res, overall_offset, MRI_SHORT, &s_output[0][0]);
      Free2DShortArray(s_output);
    }
  Acct(PROCESSING);
}

static void
WriteImage ()
{
  int i, j;
  float v, maxv;

  /* If there is no image, best not to write it */
  if (!c.do_recon) return;

  maxv = 0.0;
  if (r.recon_failed) {
    for (i = 0; i < c.res; ++i)
      for (j = 0; j < c.res; ++j)
	{
	  outim[i][j] = 0.0;
	  sh_outim[i][j] = 0;
	}
  }
  else {
    for (i = 0; i < c.res; ++i)
      for (j = 0; j < c.res; ++j)
	{
	  v = c.output_scale * outim[i][j] / (c.os_res * c.os_res);
	  outim[i][j] = v;
	  sh_outim[i][j] = v;
	  if (v > maxv)
	    maxv = v;
	}
  }
  if (c.verbosity >= VERBOSITY_TIMESTEP)
    Report("max=%3d", (int) maxv);

  Acct(WRITING);
  if (c.verbosity >= VERBOSITY_TIMESTEP)
    Report("Wr.");
  if (c.output_float)
    mri_set_image(ods, t.overall_image_num, t.slice_num - c.start_slice, MRI_FLOAT, &outim[0][0]);
  else
    mri_set_image(ods, t.overall_image_num, t.slice_num - c.start_slice, MRI_SHORT, &sh_outim[0][0]);
  Acct(PROCESSING);
  if (c.verbosity >= VERBOSITY_TIMESTEP)
    Report("\n");
}

static void
WriteImagePhase (int coil_num)
{
  int i, j;
  int hr;
  int offs, imnum;
  float imi, imr;
  int overall_offset;

  hr = c.res / 2;
  offs = (c.over_samp-1)*hr;
  for (i = 0; i < c.res; i++)
    for (j = 0; j < c.res; j++)
      if (mask[i][j])
	{
	  imr = im[0][-j+c.res+offs-1][-i+c.res+offs-1];
	  imi = im[1][-j+c.res+offs-1][-i+c.res+offs-1];
	  imphase[i][j] = 1000 * atan2(imi,imr);
	}
      else
        imphase[i][j] = 0;

  overall_offset = (((coil_num*c.total_images + t.overall_image_num)
		     * c.nslices + t.slice_num 
		     - c.start_slice)) * c.res * c.res;
  Acct(WRITING);
  mri_set_chunk(ods, "phase", c.res*c.res, overall_offset,
		MRI_SHORT, &imphase[0][0]);
  Acct(PROCESSING);
}

static float
kaiser (float l,
	float b,
	float u)
{
  float f;

  if (fabs(u) <= 0.5*fabs(l))
    {
      f = bessi0(b*sqrt(1.0-(2.0*u/l)*(2.0*u/l)));
      return f;
    }
  else
    return(0.0);
}

static float
bessi0 (float x)
{
  float ax,ans;
  double y;
  
  if ((ax=fabs(x)) < 3.75)
    {
      y=x/3.75;
      y*=y;
      ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
    }
  else
    {
      y=3.75/ax;
      ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
			+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
			+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
			+y*0.392377e-2))))))));
    }
  return(ans);
}

#ifdef AFS

static int MySystem (const char *command) {
  int pid, status;
  extern char **environ;
  extern int errno;
  
  if (command == 0)
    return 1;
  pid = fork();
  if (pid == -1)
    return -1;
  if (pid == 0) {
    char *argv[4];
    argv[0] = "sh";
    argv[1] = "-c";
    argv[2] = (char*)command;
    argv[3] = 0;
    execve("/bin/sh", argv, environ);
    exit(127);
  }
  do {
    if (waitpid(pid, &status, 0) == -1) {
      if (errno != EINTR)
        return -1;
    } else
      return status;
  } while(1);
}


static int CheckForAFS( const char* fname )
{
  char buf[256];
  int retval;
  snprintf(buf,sizeof(buf),"fs whereis %s >/dev/null 2>&1",(char*)fname);
  retval= MySystem(buf);
  return retval;
}

static void FlushAFS( const char* path )
{
  char buf[256];
  snprintf(buf,sizeof(buf),"fs flushvolume -path %s >/dev/null 2>&1",path);
  MySystem(buf);
}

#endif
