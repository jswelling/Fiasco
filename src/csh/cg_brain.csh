#!/bin/csh -efx
# cg_brain.csh
#/************************************************************
# *                                                          *
# *  Permission is hereby granted to any individual or       *
# *  institution for use, copying, or redistribution of      *
# *  this code and associated documentation, provided        *
# *  that such code and documentation are not sold for       *
# *  profit and the following copyright notice is retained   *
# *  in the code and documentation:                          *
# *     Copyright (c) 1997                                   *
# *                        Department of Statistics,         *
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
echo '#'`date` $0
echo '#$Id: cg_brain.csh,v 1.19 2003/09/17 22:52:13 welling Exp $'

if(! -d stat)mkdir stat
if(! -d ps) mkdir ps
if (`mri_printfield -fld images.dimensions $1` == "xyzt") then
	mri_remap -chunk images -order vxyzt $1
endif
set i_own_permute = 0
if ( ( ! -r $3 ) && ( ! -r $3.mri ) ) then
  set i_own_permute = 1
  mri_permute -memlimit 32000000 -chunk images -order vtxyz $1 $3
          
endif
echo '#'`date`
smooth_missing -i $3.mri -h ${3}_s.mri
if ( $i_own_permute ) then
  mri_destroy_dataset $3
endif

# Generate vmpfx input parameter file
generate_vmpfx_in -proto $F_VMPFX_PROTO -nimage $F_NIMAGE -nslice \
	$F_NSLICE -droot ${3:t}_s -dpath ${3:h}/ -split $F_SPLIT_NEW \
	-cond $F_SPLIT_COND -fixed $F_VMPFX_FIXED -iai $F_IAI > $F_VMPFX_INPUT

# Run vmpfx
vmpfx $F_VMPFX_INPUT

# Clean up some storage
mri_destroy_dataset ${3}_s

# Move vmpfx results to appropriate directories
mv vmpfx_results.out vmpfx_results.log out

# Add up all the parts to figure out how many results we will generate
# Total is estimates plus standard errors for a whole bunch of parameters
set ofile = out/vmpfx_results.out
set n_estimates = `awk '/estimates/ { print $2 ; exit }' < $ofile`
set n_errors = `awk '/standard.errors/ { print $2 ; exit }' < $ofile`
set n_diag_bytes = `awk '/diagnostics/ { print $3 ; exit }' < $ofile`
set n_null_est = `awk '/null.estimates/ { print $2 ; exit }' < $ofile`
set n_null_err = `awk '/null.standard.errors/ { print $2 ; exit }' < $ofile`
set n_null_diag_bytes = `awk '/null.diagnostics/ { print $3 ; exit }' < $ofile`
set n_mp_base = `awk '/baseline/ { print $2 ; exit }' < $ofile`
set n_mp_drift_coef = `awk '/drift.coef/ { print $2 ; exit }' < $ofile`
set n_mp_noise_prec = `awk '/noise.precision/ { print $2 ; exit }' < $ofile`
set n_mp_shape = `awk '/shape/ { print $2 ; exit }' < $ofile`
set n_mp_resp = `awk '/responsiveness/ { print $2 ; exit }' < $ofile`
set n_cycle = `awk '/ mp / { print $3 ; exit }' < $ofile`
@ n_total = $n_estimates + $n_errors

# Scream and die if the output format of vmpfx has changed.  (This
# only catches a tiny fraction of possible changes, of course)
@ n_cycle_temp = ( 8 * ( $n_estimates + $n_errors + $n_null_est \
	+ $n_null_err) ) + $n_diag_bytes + $n_null_diag_bytes 
if ( $n_cycle != $n_cycle_temp ) then
	echo 'Unexpected format for vmpfx output!'
        exit -1
endif
if ( $n_estimates != $n_errors ) then
	echo 'Unexpected format for vmpfx subrecord output!'
	exit -1
endif
if ( $n_null_diag_bytes != 24 ) then
	echo 'Unexpected format for vmpfx null diagnostic output!'
        exit -1
endif
if ( $n_estimates != ( $n_mp_base + $n_mp_drift_coef + \
	+ $n_mp_noise_prec + $n_mp_shape + $n_mp_resp )) then
	echo 'Unexpected format for vmpfx estimates output!'
	exit -1
endif

# Copy the binary data representing estimates and standard errors, and
# Assemble the output file
@ n_data_bytes = $n_total * 8
@ bayes_fac_offset = $n_cycle - 12
set vdim = `mri_printfield -field images.extent.v $3`
set xdim = `mri_printfield -field images.extent.x $3`
set ydim = `mri_printfield -field images.extent.y $3`
set zdim = `mri_printfield -field images.extent.z $3`
if ( `fiasco_getarch.csh` == LINUX ) then
  set whichend = "-littleendian"
else
  set whichend = "-bigendian"
endif

@ vskip = $n_cycle - $n_data_bytes
smartreader -i vmpfx_results.bin -out $2 -ignoreheader \
  -dimorder vxyz -dims ${n_total}:${xdim}:${ydim}:${zdim} -type double \
  -offset 0 -def skip.v=$vskip $whichend

@ vskip = $n_cycle - 4
smartreader -i vmpfx_results.bin -out vmpfx_bfac_tmp -ignoreheader \
  -dimorder vxyz -dims 1:${xdim}:${ydim}:${zdim} -type float \
  -offset $bayes_fac_offset -def skip.v=$vskip $whichend \
  -def chunkfile=".bfac_dat"
mri_remap -order xyz vmpfx_bfac_tmp

mri_copy_chunk -chunk images -chunk_out lnbfac vmpfx_bfac_tmp $2
mri_destroy_dataset vmpfx_bfac_tmp

rm vmpfx_results.bin

echo '#'`date`
