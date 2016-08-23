/*
 *	commun.c
 *
 *    Communication routines for the parallel master/worker configuration of
 *    the "spiral" program.  These are based on PVM, and are entirely
 *    ignored if "PVM" is not defined.  For more information about
 *    the "spiral program, see spiral.c.
 *
 *    Copyright (c) 1998 by the Pittsburgh Supercomputing Center
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
#include <stdio.h>
#include <stdlib.h>
#include "fmri.h"
#include "array.h"
#include "spiral.h"
#include "par.h"

static char rcsid[] = "$Id: commun.c,v 1.22 2003/04/22 21:43:51 welling Exp $";

void
PackContext ()
{
  par_pkstr(c.input_directory);
  par_pkstr(c.reference_directory);
  par_pkint(c.big_endian_input);
  par_pkstr(c.tmp_directory);
  par_pkstr(c.output_directory);
  par_pkint(c.big_endian_output);
  par_pkint(c.input_is_mri);
  par_pkstr(c.output_name);
  par_pkstr(c.output_data);
  par_pkint(c.nfiles);
  par_pkint(c.ncoils);
  par_pkint(c.nslices);
  par_pkint(c.sliceorder);
  par_pkint(c.nimages);
  par_pkint(c.nph1);
  par_pkint(c.nphmult);
  par_pkint(c.npr);
  par_pkint(c.ndat);
  par_pkint(c.ndatfr);
  par_pkint(c.total_images);
  par_pkint(c.gtype);
  par_pkint(c.opxres);
  par_pkint(c.ngap);
  par_pkint(c.concat);
  par_pkint(c.chop);
  par_pkint(c.risetime);
  par_pkint(c.densamp);
  par_pkdouble(c.ts);
  par_pkdouble(c.gts);
  par_pkdouble(c.fsgcm);
  par_pkdouble(c.opfov);
  par_pkdouble(c.slewrate);
  par_pkdouble(c.fast_rec_lpf);
  par_pkdouble(c.mapdel);
  par_pkint(c.mapdel_set);
  par_pkfloat(c.pix_shifth);
  par_pkfloat(c.pix_shiftv);
  par_pkint(c.coil_record_length);
  par_pkint(c.slice_record_length);
  par_pkint(c.baseline_length);
  par_pkint(c.res);
  par_pkint(c.os_res);
  par_pkint(c.slice);
  par_pkint(c.samp_delay);
  par_pkint(c.samp_cor);
  par_pkfloat(c.ph_twist);
  par_pkint(c.lr_shift);
  par_pkint(c.tb_shift);
  par_pkint(c.loc_shift);
  par_pkfloat(c.zoom);
  par_pkfloat(c.mag_factor);
  par_pkfloat(c.ph_factor);
  par_pkint(c.start_slice);
  par_pkint(c.end_slice);
  par_pkint(c.lin_cor);
  par_pkint(c.lin_map);
  par_pkint(c.hc_cor);
  par_pkint(c.gen_map);
  par_pkstr(c.reg_file);
  par_pkint(c.reg_2x);
  par_pkint(c.write_raw);
  par_pkint(c.write_mag);
  par_pkint(c.write_samples);
  par_pkint(c.write_phase);
  par_pkint(c.do_recon);
  par_pkint(c.over_samp);
  par_pkfloat(c.grid_len);
  par_pkfloat(c.gridb);
  par_pkfloat(c.wind);
  par_pkfloat(c.factxx);
  par_pkfloat(c.factxy);
  par_pkfloat(c.factyx);
  par_pkfloat(c.factyy);
  par_pkfloat(c.samp_time);
  par_pkint(c.filter_sz);
  par_pkint(c.all_coils);
  par_pkint(c.output_float);
  par_pkfloat(c.output_scale);
  par_pkint(c.verbosity);
  par_pkint(c.reverse);
  par_pkfloatarray(&c.t2k[0][0][0], 2*c.npr*c.ndat);
  par_pkfloatarray(c.kdens, c.ndat);
  if (c.hc_cor) {
    par_pkfloatarray(&c.refim_in[0][0][0], c.nslices*c.res*c.res);
    par_pkfloatarray(&c.sampim[0][0][0], c.nslices*c.os_res*c.os_res);
  }
  par_pkfloatarray(&c.refl[0][0], 3*c.nslices);
  par_pkbytearray(c.ref_missing, c.nslices);
}

void
UnpackContext ()
{
  static int previously_allocated = FALSE;
  static int previously_allocated_hc = FALSE;
  static int previously_allocated_lc = FALSE;
  static int previously_allocated_rm = FALSE;

  par_upkstr(c.input_directory);
  par_upkstr(c.reference_directory);
  c.big_endian_input= par_upkint();
  par_upkstr(c.tmp_directory);
  par_upkstr(c.output_directory);
  c.big_endian_output= par_upkint();
  c.input_is_mri= par_upkint();
  par_upkstr(c.output_name);
  par_upkstr(c.output_data);
  c.nfiles= par_upkint();
  c.ncoils= par_upkint();
  c.nslices= par_upkint();
  c.sliceorder= par_upkint();
  c.nimages= par_upkint();
  c.nph1= par_upkint();
  c.nphmult= par_upkint();
  c.npr= par_upkint();
  c.ndat= par_upkint();
  c.ndatfr= par_upkint();
  c.total_images= par_upkint();
  c.gtype= par_upkint();
  c.opxres= par_upkint();
  c.ngap= par_upkint();
  c.concat= par_upkint();
  c.chop= par_upkint();
  c.risetime= par_upkint();
  c.densamp= par_upkint();
  c.ts= par_upkdouble();
  c.gts= par_upkdouble();
  c.fsgcm= par_upkdouble();
  c.opfov= par_upkdouble();
  c.slewrate= par_upkdouble();
  c.fast_rec_lpf= par_upkdouble();
  c.mapdel= par_upkdouble();
  c.mapdel_set= par_upkint();
  c.pix_shifth= par_upkfloat();
  c.pix_shiftv= par_upkfloat();
  c.coil_record_length= par_upkint();
  c.slice_record_length= par_upkint();
  c.baseline_length= par_upkint();
  c.res= par_upkint();
  c.os_res= par_upkint();
  c.slice= par_upkint();
  c.samp_delay= par_upkint();
  c.samp_cor= par_upkint();
  c.ph_twist= par_upkfloat();
  c.lr_shift= par_upkint();
  c.tb_shift= par_upkint();
  c.loc_shift= par_upkint();
  c.zoom= par_upkfloat();
  c.mag_factor= par_upkfloat();
  c.ph_factor= par_upkfloat();
  c.start_slice= par_upkint();
  c.end_slice= par_upkint();
  c.lin_cor= par_upkint();
  c.lin_map= par_upkint();
  c.hc_cor= par_upkint();
  c.gen_map= par_upkint();
  par_upkstr(c.reg_file);
  c.reg_2x= par_upkint();
  c.write_raw= par_upkint();
  c.write_mag= par_upkint();
  c.write_samples= par_upkint();
  c.write_phase= par_upkint();
  c.do_recon= par_upkint();
  c.over_samp= par_upkint();
  c.grid_len= par_upkfloat();
  c.gridb= par_upkfloat();
  c.wind= par_upkfloat();
  c.factxx= par_upkfloat();
  c.factxy= par_upkfloat();
  c.factyx= par_upkfloat();
  c.factyy= par_upkfloat();
  c.samp_time= par_upkfloat();
  c.filter_sz= par_upkint();
  c.all_coils= par_upkint();
  c.output_float= par_upkint();
  c.output_scale= par_upkfloat();
  c.verbosity= par_upkint();
  c.reverse= par_upkint();

  /* deallocate old t2k and kdens arrays */
  if (previously_allocated)
    {
      Free3DFloatArray(c.t2k);
      free(c.kdens);
      previously_allocated = FALSE;
    }

  if (previously_allocated_hc)
    {
      Free3DFloatArray(c.refim_in);
      Free3DFloatArray(c.sampim);
      previously_allocated_hc = FALSE;
    }

  if (previously_allocated_lc)
    {
      Free2DFloatArray(c.refl);
      previously_allocated_lc = FALSE;
    }

  if (previously_allocated_rm)
    {
      free(c.ref_missing);
      previously_allocated_rm = FALSE;
    }

  /* allocate the t2k and kdens arrays */
  c.t2k = Alloc3DFloatArray(2, c.npr, c.ndat);
  c.kdens = (float *) malloc(c.ndat * sizeof(float));
  previously_allocated = TRUE;

  par_upkfloatarray(&c.t2k[0][0][0], 2*c.npr*c.ndat);
  par_upkfloatarray(c.kdens, c.ndat);

  if (c.hc_cor)
    {
      /* allocate the arrays */
      c.refim_in = Alloc3DFloatArray(c.nslices, c.res, c.res);
      c.sampim = Alloc3DFloatArray(c.nslices, c.os_res, c.os_res);
      previously_allocated_hc = TRUE;
      par_upkfloatarray(&c.refim_in[0][0][0], c.nslices*c.res*c.res);
      par_upkfloatarray(&c.sampim[0][0][0], c.nslices*c.os_res*c.os_res);
    }

  c.refl= Alloc2DFloatArray(c.nslices, 3);
  previously_allocated_lc = TRUE;
  par_upkfloatarray(&c.refl[0][0], 3*c.nslices);

  if (!(c.ref_missing= (unsigned char*)malloc(c.nslices)))
    Abort("Could not allocate %d bytes!\n",c.nslices);
  previously_allocated_rm= TRUE;
  par_upkbytearray(c.ref_missing, c.nslices);
}

void
PackTask ()
{
  par_pkint(t.file_index);
  par_pkstr(t.filename);
  par_pkint(t.slice_num);
  par_pkint(t.image_num);
  par_pkint(t.overall_image_num);
  par_pkfloat(t.reg_xs);
  par_pkfloat(t.reg_ys);
  par_pkfloat(t.reg_rot);
}

void
UnpackTask ()
{
  t.file_index= par_upkint();
  par_upkstr(t.filename);
  t.slice_num= par_upkint();
  t.image_num= par_upkint();
  t.overall_image_num= par_upkint();
  t.reg_xs= par_upkfloat();
  t.reg_ys= par_upkfloat();
  t.reg_rot= par_upkfloat();
}

void
PackResult ()
{
  par_pkint(r.image_num);
  par_pkint(r.slice_num);
  par_pkint(r.sampim_valid);
  par_pkint(r.recon_failed);
  if (r.sampim_valid)
    par_pkfloatarray(&r.slice_sampim[0][0], c.os_res*c.os_res);
}

void
UnpackResult ()
{
  static int previously_allocated = FALSE;

  r.image_num= par_upkint();
  r.slice_num= par_upkint();
  r.sampim_valid= par_upkint();
  r.recon_failed= par_upkint();
  if (r.sampim_valid)
    {
      if (!previously_allocated)
	{
	  r.slice_sampim = Alloc2DFloatArray(c.os_res, c.os_res);
	  previously_allocated = TRUE;
	}
      par_upkfloatarray(&r.slice_sampim[0][0], c.os_res*c.os_res);
    }
}

