/*
 *	mriheader.c
 *
 *    Routines to read data from a Pgh MRI file which contains the
 *    same information as a P file header for the "spiral" program.
 *    For more information about the "spiral" program, see spiral.c.
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
 *      9/01 - modified to handle Pgh MRI files
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mri.h"
#include "fmri.h"
#include "bio.h"
#include "array.h"
#include "acct.h"
#include "spiral.h"

static char rcsid[] = "$Id: mriheader.c,v 1.11 2005/07/07 20:07:12 welling Exp $";

static void CalculateLocation (float *p1, float *p2, float *p3, int rot, int trans);

void 
CreateMRIFile (const Filename output_file, int argc, char** argv, Context *c)
{
  MRI_Dataset *ds;
  int ishot;
  int icoil;
  int islice;

  Acct(WRITEOPEN);
  ds = mri_open_dataset(output_file, MRI_WRITE);
  hist_add_cl(ds,argc,argv);
  
  /* create the samples chunk */
  Acct(WRITING);
  mri_create_chunk(ds, "samples");
  mri_set_string(ds, "samples.datatype", "int16");
  mri_set_string(ds, "samples.dimensions", "vpsczt");
  mri_set_int(ds, "samples.extent.v", 2);
  mri_set_string(ds, "samples.description.v", "complex real/imaginary");
  mri_set_int(ds, "samples.extent.z", c->nslices);
  mri_set_string(ds, "samples.description.z", "gridded image-space");
  mri_set_int(ds, "samples.extent.t", c->total_images);
  mri_set_string(ds, "samples.description.t", "gridded image-space");
  mri_set_int(ds, "samples.extent.c", c->ncoils);
  mri_set_string(ds, "samples.description.c", "discrete");
  mri_set_int(ds, "samples.extent.p", c->ndat);
  mri_set_string(ds, "samples.description.p", "ungridded k-space");
  mri_set_int(ds, "samples.extent.s", c->npr);
  mri_set_string(ds, "samples.description.s", "discrete");
  
  mri_set_string(ds, "samples.psd", c->psd);
  mri_set_int(ds, "samples.nph1", c->nph1);
  mri_set_int(ds, "samples.nphmult", c->nphmult);
  mri_set_int(ds, "samples.nimages", c->nimages);
  mri_set_int(ds, "samples.ndatfr", c->ndatfr);
  mri_set_int(ds, "samples.chop", c->chop);
  mri_set_int(ds, "samples.nex", c->nex);
  mri_set_int(ds, "samples.densamp", c->densamp);
  mri_set_int(ds, "samples.sliceorder", c->sliceorder);
  mri_set_int(ds, "samples.gtype", c->gtype);
  mri_set_int(ds, "samples.reverse", c->reverse);
  mri_set_int(ds, "samples.opxres", c->opxres);
  mri_set_float(ds, "samples.slewrate", c->slewrate);
  mri_set_int(ds, "samples.ngap", c->ngap);
  mri_set_int(ds, "samples.concat", c->concat);
  mri_set_float(ds, "samples.fast_rec_lpf", c->fast_rec_lpf);
  if (c->mapdel_set)
    mri_set_float(ds, "samples.mapdel", c->mapdel);
  mri_set_float(ds, "samples.samp_time", c->samp_time);
  mri_set_float(ds, "samples.fsgcm", c->fsgcm);
  mri_set_int(ds, "samples.risetime", c->risetime);
  mri_set_float(ds, "samples.opfov", c->opfov);
  mri_set_float(ds, "samples.grad_time_spacing", c->gts*1e6);
  mri_set_string(ds, "samples.time", c->time);
  mri_set_string(ds, "samples.date", c->date);
  mri_set_float(ds, "samples.slthick", c->slthick);
  mri_set_int(ds, "samples.tr", c->tr);
  mri_set_int(ds, "samples.te", c->te);
  mri_set_int(ds, "samples.flip", c->flip);
  mri_set_float(ds, "samples.spacing", c->spacing);
  mri_set_float(ds, "samples.tlc.0", c->tlc[0]);
  mri_set_float(ds, "samples.tlc.1", c->tlc[1]);
  mri_set_float(ds, "samples.tlc.2", c->tlc[2]);
  mri_set_float(ds, "samples.trc.0", c->trc[0]);
  mri_set_float(ds, "samples.trc.1", c->trc[1]);
  mri_set_float(ds, "samples.trc.2", c->trc[2]);
  mri_set_float(ds, "samples.brc.0", c->brc[0]);
  mri_set_float(ds, "samples.brc.1", c->brc[1]);
  mri_set_float(ds, "samples.brc.2", c->brc[2]);
  mri_set_float(ds, "samples.ctr.0", c->ctr[0]);
  mri_set_float(ds, "samples.ctr.1", c->ctr[1]);
  mri_set_float(ds, "samples.ctr.2", c->ctr[2]);
  mri_set_int(ds,"samples.index_rots", c->rotation);
  mri_set_int(ds,"samples.index_transpose", c->transpose);
  if (c->loc_shift) {
    mri_set_float(ds, "samples.pix_shifth", c->pix_shifth);
    mri_set_float(ds, "samples.pix_shiftv", c->pix_shiftv);
  }
  mri_set_string(ds, "samples.file", ".dat");

  mri_create_chunk(ds, "sample_kxloc");
  mri_set_string(ds, "sample_kxloc.datatype", "float32");
  mri_set_string(ds, "sample_kxloc.dimensions", "pscz");
  mri_set_int(ds, "sample_kxloc.extent.z", c->nslices);
  mri_set_string(ds, "sample_kyloc.description.z", "gridded image-space");
  mri_set_int(ds, "sample_kxloc.extent.c", c->ncoils);
  mri_set_string(ds, "sample_kxloc.description.c", "discrete");
  mri_set_int(ds, "sample_kxloc.extent.p", c->ndat);
  mri_set_string(ds, "sample_kxloc.description.p", "ungridded k-space");
  mri_set_int(ds, "sample_kxloc.extent.s", c->npr);
  mri_set_string(ds, "sample_kxloc.description.s", "discrete");
  mri_set_string(ds, "sample_kxloc.file", ".sampxloc");
  mri_create_chunk(ds, "sample_kyloc");
  mri_set_string(ds, "sample_kyloc.datatype", "float32");
  mri_set_string(ds, "sample_kyloc.dimensions", "pscz");
  mri_set_int(ds, "sample_kyloc.extent.z", c->nslices);
  mri_set_string(ds, "sample_kyloc.description.z", "gridded image-space");
  mri_set_int(ds, "sample_kyloc.extent.c", c->ncoils);
  mri_set_string(ds, "sample_kyloc.description.c", "discrete");
  mri_set_int(ds, "sample_kyloc.extent.p", c->ndat);
  mri_set_string(ds, "sample_kyloc.description.p", "ungridded k-space");
  mri_set_int(ds, "sample_kyloc.extent.s", c->npr);
  mri_set_string(ds, "sample_kyloc.description.s", "discrete");
  mri_set_string(ds, "sample_kyloc.file", ".sampyloc");
  
  for (islice=0; islice<c->nslices; islice++)
    for (icoil=0; icoil<c->ncoils; icoil++)
      for (ishot=0; ishot<c->npr; ishot++) {
	int offset = ((((islice*c->ncoils)+icoil)*c->npr)+ishot)*c->ndat;
	if (c->reverse) {
	  /* By convention, Pgh MRI files store the trajectory in the
	   * same order as the samples.  Since these samples are
	   * reversed, we must reverse the trajectory also.
	   */
	  float* xvals;
	  float* yvals;
	  int p;
	  if (!(xvals= (float*)malloc(c->ndat*sizeof(float))))
	    Abort("Unable to allocate %d bytes!\n",c->ndat*sizeof(float));
	  if (!(yvals= (float*)malloc(c->ndat*sizeof(float))))
	    Abort("Unable to allocate %d bytes!\n",c->ndat*sizeof(float));
	  for (p=0; p<c->ndat; p++) {
	    xvals[p]= c->t2k[0][ishot][c->ndat-(p+1)];
	    yvals[p]= c->t2k[1][ishot][c->ndat-(p+1)];
	  }
	  mri_set_chunk(ds, "sample_kxloc", c->ndat, offset,
			MRI_FLOAT, xvals);
	  mri_set_chunk(ds, "sample_kyloc", c->ndat, offset,
			MRI_FLOAT, yvals);
	  free(xvals);
	  free(yvals);
	}
	else {
	  mri_set_chunk(ds, "sample_kxloc", c->ndat, offset,
			MRI_FLOAT, c->t2k[0][ishot]);
	  mri_set_chunk(ds, "sample_kyloc", c->ndat, offset,
			MRI_FLOAT, c->t2k[1][ishot]);
	}
      }
  
  Acct(WRITECLOSE);
  mri_close_dataset(ds);
  Acct(PROCESSING);
}

void
ReadMRIFileHeader (const Filename input_file, Context *c)
{
  MRI_Dataset *ds;
  char com[sizeof(Filename)+128];
  int k;
  float tmp1, tmp2;
  int islice;
  int icoil;
  int ishot;

  Acct(READOPEN);
  /* open the file */
#ifdef AFS
  if (strncmp(input_file, "/afs", 4) == 0)
    {
      sprintf(com, "fs flushvolume -path %s", input_file);
      system(com);
    }
#endif
  if (!(ds= mri_open_dataset(input_file, MRI_READ)))
    Abort("Can't open %s\n", input_file);
  if (!mri_has(ds,"samples"))
    Abort("Input dataset %s has no 'samples' chunk!\n",input_file);
  if (!mri_has(ds,"samples.dimensions"))
    Abort("Input dataset %s has no 'samples.dimensions' tag!\n",input_file);
  if (strcmp(mri_get_string(ds,"samples.dimensions"),"vpsczt"))
    Abort("Input dataset %s sample dimensions are not vpsczt!\n",input_file);
  if (!mri_has(ds,"sample_kxloc"))
    Abort("Input dataset %s has no 'sample_kxloc' chunk!\n",input_file);
  if (!mri_has(ds,"sample_kxloc.dimensions"))
    Abort("Input dataset %s has no 'sample_kxloc.dimensions' tag!\n",
	  input_file);
  if (strcmp(mri_get_string(ds,"sample_kxloc.dimensions"),"pscz"))
    Abort("Input dataset %s kxloc dimensions are not pscz!\n",input_file);
  if (!mri_has(ds,"sample_kyloc"))
    Abort("Input dataset %s has no 'sample_kyloc' chunk!\n",input_file);
  if (!mri_has(ds,"sample_kyloc.dimensions"))
    Abort("Input dataset %s has no 'sample_kyloc.dimensions' tag!\n",
	  input_file);
  if (strcmp(mri_get_string(ds,"sample_kyloc.dimensions"),"pscz"))
    Abort("Input dataset %s kyloc dimensions are not pscz!\n",input_file);
  
  /* extract the context info */
  c->nslices= mri_get_int(ds,"samples.extent.z");
  c->ncoils= mri_get_int(ds,"samples.extent.c");
  c->ndat= mri_get_int(ds,"samples.extent.p");
  c->npr = mri_get_int(ds,"samples.extent.s");
  
  if (c->nphmult < 1)
    c->nphmult = 1;
  if (c->gen_map || c->lin_map) {
    c->nimages= c->total_images= 2;
  }
  else {
    c->total_images= mri_get_int(ds,"samples.extent.t");
    if (mri_has(ds,"samples.nimages"))
      c->nimages = mri_get_int(ds,"samples.nimages");
    else c->nimages= c->total_images;
  }
  if (mri_has(ds,"samples.nph1"))
    c->nph1 = mri_get_int(ds,"samples.nph1");
  else 
    c->nph1 = 1;
  if (mri_has(ds,"samples.nphmult"))
    c->nphmult = mri_get_int(ds,"samples.nphmult");
  else c->nphmult = c->total_images;
  if (mri_has(ds,"samples.ndatfr"))
    c->ndatfr = mri_get_int(ds,"samples.ndatfr");
  else
    c->ndatfr = mri_get_int(ds,"samples.extent.p");
  c->chop = 1;
  c->nex = 1;
  c->densamp = 5;
  /* Slice order: 0=interleaved, 1=sequential */
  if (mri_has(ds,"samples.sliceorder"))
    c->sliceorder = mri_get_int(ds,"samples.sliceorder");
  else if (mri_has(ds,"samples.reorder") && 
	   mri_get_int(ds,"samples.reorder")==1)
    c->sliceorder= 0;
  else c->sliceorder= 1;
  /* glover cv's */
  if (mri_has(ds,"samples.gtype"))
    c->gtype = mri_get_int(ds,"samples.gtype");
  else c->gtype= 0; /* doesn't matter since trajectory loaded from file */
  if (mri_has(ds,"samples.reverse"))
    c->reverse = mri_get_int(ds,"samples.reverse");
  else c->reverse= 0;
  c->opxres = mri_get_int(ds,"samples.opxres");
  if (mri_has(ds,"samples.slewrate"))
    c->slewrate = mri_get_float(ds,"samples.slewrate");
  else c->slewrate= 0.0; /* slewrate is unused because trajectory is known */
  if (mri_has(ds,"samples.ngap"))
    c->ngap = mri_get_int(ds,"samples.ngap"); /* concat readout */
  else c->ngap= 0;
  if (mri_has(ds,"samples.concat"))
    c->concat = mri_get_int(ds,"samples.concat"); /*concat readout*/
  else c->concat= 0;
  if (mri_has(ds,"samples.fast_rec_lpf"))
    c->fast_rec_lpf = 
      mri_get_float(ds,"samples.fast_rec_lpf"); /* fast rec cutoff in kHz */
  else c->fast_rec_lpf= 0.0;
  if (mri_has(ds,"samples.mapdel")) {
    if (c->mapdel_set) {
      Warning(1,"Map delay from file is overridden (%f vs. %f)\n",
	      mri_get_float(ds,"samples.mapdel"),c->mapdel);
    }
    else {
      c->mapdel_set= 1;
      c->mapdel = 
	mri_get_float(ds,"samples.mapdel"); /* field mapping offset (us) */
    }
  }
  else {
    if (!c->mapdel_set) c->mapdel = 0.0;
  }
  c->samp_time = mri_get_float(ds,"samples.samp_time"); /* in (us) */
  if (mri_has(ds,"samples.fsgcm"))
    c->fsgcm = mri_get_float(ds,"samples.fsgcm");
  else c->fsgcm= 0.0; /* only used for calculating trajectory */
  if (mri_has(ds,"samples.risetime"))
    c->risetime = mri_get_int(ds,"samples.risetime");
  else c->risetime= 0.0; /* only used for calculating trajectory */
  c->opfov = mri_get_float(ds,"samples.opfov");
  c->ts = c->samp_time*1e-6; /* in sec */
  if (mri_has(ds,"samples.grad_time_spacing"))
    c->gts = mri_get_float(ds,"samples.grad_time_spacing")*1e-6; /* gradient spacing in sec */
  else c->gts= 0.0; /* only used for calculating trajectory */
  c->coil_record_length = 0; /* There are no coil records for MRI files */
  if (mri_has(ds,"samples.time"))
    strncpy(c->time, mri_get_string(ds, "samples.time"), 8);
  else strcpy(c->time,"UNKNOWN");
  c->time[8] = '\0';
  if (mri_has(ds,"samples.date")) 
    strncpy(c->date, mri_get_string(ds, "samples.date"), 10);
  else strcpy(c->date,"UNKNOWN");
  c->date[10] = '\0';
  if (mri_has(ds,"samples.slthick"))
    c->slthick = mri_get_float(ds, "samples.slthick");
  else if (mri_has(ds,"samples.voxel_z"))
    c->slthick = mri_get_float(ds, "samples.voxel_z");
  else c->slthick= 0.0;
  c->tr = mri_get_int(ds, "samples.tr");
  c->te = mri_get_int(ds, "samples.te");
  if (mri_has(ds,"samples.flip"))
    c->flip = mri_get_int(ds, "samples.flip");
  else c->flip= 0.0;
  if (mri_has(ds,"samples.spacing"))
    c->spacing = mri_get_float(ds, "samples.spacing");
  else if (mri_has(ds,"samples.slice_gap"))
    c->spacing = mri_get_float(ds, "samples.slice_gap");
  else c->spacing= 0.0;
  if (mri_has(ds,"samples.tlc.0")) 
    c->tlc[0] = mri_get_float(ds, "samples.tlc.0");
  else c->tlc[0]= 0.0;
  if (mri_has(ds,"samples.tlc.1")) 
    c->tlc[1] = mri_get_float(ds, "samples.tlc.1");
  else c->tlc[1]= 0.0;
  if (mri_has(ds,"samples.tlc.2")) 
    c->tlc[2] = mri_get_float(ds, "samples.tlc.2");
  else c->tlc[2]= 0.0;
  if (mri_has(ds,"samples.trc.0")) 
    c->trc[0] = mri_get_float(ds, "samples.trc.0");
  else c->trc[0]= 0.0;
  if (mri_has(ds,"samples.trc.1")) 
    c->trc[1] = mri_get_float(ds, "samples.trc.1");
  else c->trc[1]= 0.0;
  if (mri_has(ds,"samples.trc.2")) 
    c->trc[2] = mri_get_float(ds, "samples.trc.2");
  else c->trc[2]= 0.0;
  if (mri_has(ds,"samples.brc.0")) 
    c->brc[0] = mri_get_float(ds, "samples.brc.0");
  else c->brc[0]= 0.0;
  if (mri_has(ds,"samples.brc.1")) 
    c->brc[1] = mri_get_float(ds, "samples.brc.1");
  else c->brc[1]= 0.0;
  if (mri_has(ds,"samples.brc.2")) 
    c->brc[2] = mri_get_float(ds, "samples.brc.2");
  else c->brc[2]= 0.0;
  if (mri_has(ds,"samples.ctr.0")) 
    c->ctr[0] = mri_get_float(ds, "samples.ctr.0");
  else c->ctr[0]= 0.0;
  if (mri_has(ds,"samples.ctr.1")) 
    c->ctr[1] = mri_get_float(ds, "samples.ctr.1");
  else c->ctr[1]= 0.0;
  if (mri_has(ds,"samples.ctr.2")) 
    c->ctr[2] = mri_get_float(ds, "samples.ctr.2");
  else c->ctr[2]= 0.0;

  if (mri_has(ds,"samples.psd"))
    strncpy(c->psd, mri_get_string(ds, "samples.psd"), 33);
  else if (mri_has(ds,"samples.pulse_seq"))
    strncpy(c->psd, mri_get_string(ds, "samples.pulse_seq"), 33);
  else strncpy(c->psd, "UNKNOWN", 33);
  c->psd[33] = '\0';
  if (mri_has(ds,"samples.index_rots"))
    c->rotation= mri_get_int(ds, "samples.index_rots");
  else c->rotation= 0;
  if (mri_has(ds,"samples.index_transpose"))
    c->transpose= mri_get_int(ds, "samples.index_transpose");
  else c->transpose= 0;
  if (mri_has(ds,"samples.pix_shifth"))
    c->pix_shifth= mri_get_float(ds, "samples.pix_shifth");
  else c->pix_shifth= 0.0;
  if (mri_has(ds,"samples.pix_shiftv"))
    c->pix_shiftv= mri_get_float(ds, "samples.pix_shiftv");
  else c->pix_shiftv= 0.0;

  c->baseline_length = c->ndat*4;
  c->slice_record_length = 
    c->nimages*c->npr*c->ndat*4 + c->nph1*c->baseline_length;
  if (((c->nphmult * c->npr) % 2) == 1)
    c->slice_record_length += c->nph1*c->ndat*4;
  c->os_res = c->over_samp * c->res;

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
      if (c->mapdel_set) 
	Report("\tfast_rec_lpf = %f kHz, mapdel = %f us\n",
	       c->fast_rec_lpf,c->mapdel);
      else
	Report("\tfast_rec_lpf = %f kHz, mapdel not set\n",
	       c->fast_rec_lpf);
      Report("\tts = %f, gts = %f\n",c->ts*1e6,c->gts*1e6);
      Report("\tnumber of coils = %d\n", c->ncoils);
      Report("\tslice thickness = %f mm, slice spacing = %f mm\n",
	     c->slthick, c->spacing);
      Report("\trot = %d\n",c->rotation);
      Report("\ttrans = %d\n",c->transpose);
    }
  
  if (c->loc_shift)
    {
      if (c->verbosity > VERBOSITY_FILE) 
	{
	  Report("  shifth = %f, shiftv = %f (in mm)\n", 
		 c->pix_shifth, c->pix_shiftv);
	}
    }
  
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
  if (!(c->kdens = (float *) malloc(c->ndat * sizeof(float))))
    Abort("Unable to allocate %d bytes!\n",c->ndat * sizeof(float));
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

  /* Import trajectory info.  We are actually importing all
   * slice's data into the same place, because (currently) spiral
   * recon uses the same trajectory info for every slice.
   */
  for (islice=0; islice<c->nslices; islice++) {
    for (icoil=0; icoil<c->ncoils; icoil++) {
      for (ishot=0; ishot<c->npr; ishot++) {
	long long offset = 
	  ((((islice*c->ncoils)+icoil)*c->npr)+ishot)*c->ndat;
	if (c->reverse) {
	  /* Pgh MRI files store trajectory info in the same order as
	   * the samples, so if the data is reversed, we need to reverse
	   * these too.
	   */
	  int p;
	  float* xvals= mri_get_chunk(ds, "sample_kxloc", c->ndat, 
				      offset, MRI_FLOAT);
	  float* yvals= mri_get_chunk(ds, "sample_kyloc", c->ndat, 
				      offset, MRI_FLOAT);
	  for (p=0; p<c->ndat; p++) {
	    c->t2k[0][ishot][p]= xvals[c->ndat-(p+1)];
	    c->t2k[1][ishot][p]= yvals[c->ndat-(p+1)];
	  }
	}
	else {
	  mri_read_chunk(ds, "sample_kxloc", c->ndat, offset,
			 MRI_FLOAT, c->t2k[0][ishot]);
	  mri_read_chunk(ds, "sample_kyloc", c->ndat, offset,
			 MRI_FLOAT, c->t2k[1][ishot]);
	}
      }
    }
  }

  /* Close the dataset so that workers can reopen it later */
  mri_close_dataset(ds);
}

void
CheckMRIHeaderInfo (const Filename input_file, Context *c)
{
  FILE *f;
  Context cTmp;

  cTmp.verbosity= c->verbosity;
  ReadMRIFileHeader(input_file, &cTmp);

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
      c->mapdel_set != cTmp.mapdel_set ||
      (c->mapdel_set && (c->mapdel != cTmp.mapdel)) ||
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

