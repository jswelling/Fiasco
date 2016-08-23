#!/bin/csh -exf
# siemens_ramp_resample.csh
#/************************************************************
# *                                                          *
# *  Permission is hereby granted to any individual or       *
# *  institution for use, copying, or redistribution of      *
# *  this code and associated documentation, provided        *
# *  that such code and documentation are not sold for       *
# *  profit and the following copyright notice is retained   *
# *  in the code and documentation:                          *
# *     Copyright (c) 2003 Department of Statistics,         *
# *                        Carnegie Mellon University        *
# *                                                          *
# *  This program is distributed in the hope that it will    *
# *  be useful, but WITHOUT ANY WARRANTY; without even the   *
# *  implied warranty of MERCHANTABILITY or FITNESS FOR A    *
# *  PARTICULAR PURPOSE.  Neither Carnegie Mellon University *
# *  nor any of the authors assume any liability for         *
# *  damages, incidental or otherwise, caused by the         *
# *  installation or use of this software.                   *
# *                                                          *
# *  CLINICAL APPLICATIONS ARE NOT RECOMMENDED, AND THIS     *
# *  SOFTWARE HAS NOT BEEN EVALUATED BY THE UNITED STATES    *
# *  FDA FOR ANY CLINICAL USE.                               *
# *                                                          *
# *                                                          *
# *  Original programming by Joel Welling                    *
# ************************************************************/
#
echo '$Id: siemens_ramp_resample.csh,v 1.6 2003/09/10 20:17:23 bakalj Exp $'
#
echo '#'`date`$0

# This script will be called with two parameters, the un-resampled
# file and the name to be given the resampled file.  

set infile = ${1}
set outfile = ${2}

echo "Ramp resampling ${infile} using the Siemens algorithm"

# Set things up to be tolerant of use of x or q in input dim string
set dimstr = \
  `mri_printfield -field images.dimensions -nofail ${infile}`
if ( ${dimstr} =~ *qy* ) then
  set perm_dimstr = `echo $dimstr | sed 's/yz//g' | sed 's/q/yzq/g'`
  set remap_dimstr = $perm_dimstr
  set resamp_dimstr = `echo $dimstr | sed 's/q/x/g'`
  set qdim = `mri_printfield -field images.extent.q ${infile}`
else
  set perm_dimstr = `echo $dimstr | sed 's/yz//g' | sed 's/x/yzx/g'`
  set remap_dimstr = `echo $perm_dimstr | sed 's/q/x/g'`  
  set resamp_dimstr = $dimstr
  set qdim = `mri_printfield -field images.extent.x ${infile}`
endif
if ( ${resamp_dimstr} != "vxyzt" ) then
  echo "****Error: Unexpected dimension string for ${infile}"
  exit -1
endif

# Extract ramp resampling data from the input file

set gradShape = trapezoidal
set xdim = $qdim
set fov_x = \
  `mri_printfield -field images.fov.x ${infile}`
set rampUp = \
  `mri_printfield -field images.regridRampupTime ${infile}`
set rampDown = \
  `mri_printfield -field images.regridRampdownTime ${infile}`
set flatTop = \
  `mri_printfield -field images.regridFlattopTime ${infile}`
set sampDelay = \
  `mri_printfield -field images.regridSampleDelay ${infile}`
set durADC = \
  `mri_printfield -field images.regridADCDuration ${infile}`
set lag = "half"

# Intermediate file storage
if ( ! -d data ) mkdir data
set rawlocs = 'data/resample_rawlocs'
set junk = 'data/resample_tmp1'
set junk_stretch = 'data/resample_tmp2'
set junk_stretch_complex = 'data/resample_tmp3'
set fake_samples = 'data/resample_fake_samples'
set rawlocs_zero = 'data/resample_zerolocs'
set scaledlocs = 'data/resample_scaledlocs'
set samples = 'data/resample_samples'
set samples_resampled = 'data/resample_tmp4'
set samples_resampled_fftx = 'data/resample_tmp5'
set samples_resampled_p = 'data/resample_tmp6'
set resample_sum = 'data/resample_tmp7'
set resample_sum_m = 'data/resample_tmp8'
set resample_sum_m_v = 'data/resample_tmp9'
set resample_sum_m_vq = 'data/resample_tmp10'
set junk2 = 'data/resample_tmp11'
set junk3 = 'data/resample_tmp12'
set resample = 'data/resample_matrix'

# This absurd bit is necessary because we can't do a floating point
# division in csh.
mri_counter -l 1 ${junk2}
set vox_x = \
  `mri_rpn_math -out ${junk3} '0,'$fov_x','$xdim',/,1,if_print_1' ${junk2}`

#
# We need to generate the Siemens ramp resampling matrix on the fly.
#
# Generate k-space locations for samples
pulsecalc epireadout $gradShape $xdim $rampUp $flatTop $rampDown $sampDelay \
  $durADC $lag | mri_from_ascii -ord vq -l 1:${xdim} -ind q ${rawlocs}
mri_remap -order pscz -length ${xdim}:1:1:1 ${rawlocs} 
mri_rpn_math -out ${rawlocs_zero} '0' ${rawlocs}
mri_rpn_math -out ${scaledlocs} '$1,'${xdim}',1,-,*,'${xdim}',2,/,-' ${rawlocs}

# Make up a diagonal matrix of 1's (scaled appropriately for FFT)
mri_counter -d q -l ${xdim} ${junk}
mri_remap -order vpsczt -length 1:${xdim}:1:1 ${junk}
mri_interp -d t -con -len ${xdim} ${junk} ${junk_stretch}
mri_interp -d v -con -len 2 ${junk_stretch} ${junk_stretch_complex}
mri_rpn_math -out ${fake_samples} \
    '0,1,'${xdim}',/,sqrt,$v,0,==,$p,$t,==,*,if_keep' ${junk_stretch_complex}

# Assemble a dataset which corresponds to a single peak at each
# of the sampled locations
if ( -e ${samples}.mri ) mri_destroy_dataset ${samples}
mri_copy_chunk -chunk images -chunk_out samples ${fake_samples} ${samples}
mri_copy_chunk -chunk images -chunk_out sample_kxloc ${scaledlocs} ${samples}
mri_copy_chunk -chunk images -chunk_out sample_kyloc ${rawlocs_zero} ${samples}

# Find the resampled signal strengths associated with each sample location
parallel.run.csh slow_ft -phscale 0.0 -xvoxel ${vox_x} -yvoxel ${vox_x} \
    -xres ${xdim} -yres 1 ${samples} ${samples_resampled}
mri_fft -d x -cpx ${samples_resampled} ${samples_resampled_fftx}
mri_remap -order vqx -length 2:${xdim}:${xdim} \
    ${samples_resampled_fftx} 

# Reorder the matrix, and re-weight the samples so that each output
# location gets a total contribution weight of 1.
mri_permute -order vxq ${samples_resampled_fftx} ${samples_resampled_p} 
mri_remap -order vqx -length 2:${xdim}:${xdim} \
    ${samples_resampled_p} 
mri_subsample -d q -l 1 -sum ${samples_resampled_p} ${resample_sum}
mri_complex_to_scalar ${resample_sum} ${resample_sum_m}
mri_interp -d v -len 2 -con ${resample_sum_m} ${resample_sum_m_v}
mri_interp -d q -len ${qdim} -con ${resample_sum_m_v} ${resample_sum_m_vq}
mri_rpn_math -out ${resample} '$1,$2,/' \
  ${samples_resampled_p} ${resample_sum_m_vq}

# Reorder data, apply the resampling matrix, and reorder it back
mri_permute -order ${perm_dimstr} ${infile} ${infile}_p 
mri_remap -order vyzqt -length ':::'${qdim}':' ${infile}_p 
time mri_matmult -v -complex -out ${outfile}_p ${infile}_p ${resample}
mri_destroy_dataset ${infile}_p
mri_permute -order ${resamp_dimstr} ${outfile}_p ${outfile} 
mri_destroy_dataset ${outfile}_p

# Clean up
mri_destroy_dataset ${rawlocs}
mri_destroy_dataset ${scaledlocs}
mri_destroy_dataset ${rawlocs_zero}
mri_destroy_dataset ${junk}
mri_destroy_dataset ${junk2}
mri_destroy_dataset ${junk3}
mri_destroy_dataset ${junk_stretch}
mri_destroy_dataset ${junk_stretch_complex}
mri_destroy_dataset ${fake_samples}
mri_destroy_dataset ${samples}
mri_destroy_dataset ${samples_resampled}
mri_destroy_dataset ${samples_resampled_fftx}
mri_destroy_dataset ${samples_resampled_p}
mri_destroy_dataset ${resample_sum}
mri_destroy_dataset ${resample_sum_m}
mri_destroy_dataset ${resample_sum_m_v}
mri_destroy_dataset ${resample_sum_m_vq}

