#!/bin/csh -efx
# spline_detrend.csh
#/************************************************************
# *                                                          *
# *  Permission is hereby granted to any individual or       *
# *  institution for use, copying, or redistribution of      *
# *  this code and associated documentation, provided        *
# *  that such code and documentation are not sold for       *
# *  profit and the following copyright notice is retained   *
# *  in the code and documentation:                          *
# *     Copyright (c) 1999 Department of Statistics,         *
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
echo '#$Id: spline_detrend.csh,v 1.13 2003/10/04 22:00:57 bakalj Exp $'
if(! -d stat)mkdir stat
if(! -d ps) mkdir ps
if (`mri_printfield -field ${F_DETREND_CHUNK}.dimensions $1` == "xyzt") then
	mri_remap -chunk ${F_DETREND_CHUNK} -order vxyzt $1
endif
mri_permute -memlimit 32000000 -chunk ${F_DETREND_CHUNK} -order vtxyz $1 ${1}_p
echo '#'`date`
smooth_missing -i ${1}_p.mri -h ${1}_s.mri
mri_destroy_dataset ${1}_p
echo '#'`date`

# Generate vmpfx input parameter file for null model
generate_vmpfx_in -proto $F_VMPFX_NPROTO -nimage $F_NIMAGE -nslice \
	$F_NSLICE -droot ${1:t}_s -dpath ${1:h}/ \
	-iai $F_IAI -null > $F_VMPFX_NINPUT

# Run vmpfx
echo '#'`date`
vmpfx $F_VMPFX_NINPUT
echo '#'`date`

# Move vmpfx results to appropriate directories
mv vmpfx_results.out out/spline_detrend_vmpfx.out
mv vmpfx_results.log out/spline_detrend_vmpfx.log
mv vmpfx_results.bin ${3}_d.dat

# Reality check the vmpfx output
set ofile = out/spline_detrend_vmpfx.out
set tot_size = `awk '/max.posterior/ { print $2 ; exit }' < $ofile`
set res_units = `awk '/record/ { print $2 ; exit }' < $ofile`
set res_count = `awk '/residuals/ { print $2 ; exit }' < $ofile`
set res_size = `awk '/residuals/ { print $3 ; exit }' < $ofile`
@ prod = $res_count * $res_size * $res_units
if ( $res_size != 8 ) then
	echo 'Unexpected data type for vmpfx output!'
	exit -1
endif
if ( $tot_size != $prod ) then
	echo 'Unexpected format for vmpfx output!'
	exit -1
endif

# Prepare a .mri file for the output
set images_ext_x = `mri_printfield -field images.extent.x ${1}_s`
set images_ext_y = `mri_printfield -field images.extent.y ${1}_s`
set images_ext_z = `mri_printfield -field images.extent.z ${1}_s`
set images_ext_v = `mri_printfield -field images.extent.v ${1}_s`
set images_ext_t = `mri_printfield -field images.extent.t ${1}_s`
set images_dims  = `mri_printfield -field images.dimensions ${1}_s`
echo '\!format = pgh' > ${3}_d.mri
echo '\!version = 1.0' >> ${3}_d.mri
echo 'images = [chunk]' >> ${3}_d.mri
echo 'images.datatype = float64' >> ${3}_d.mri
echo 'images.dimensions=' $images_dims >> ${3}_d.mri
echo 'images.file= .dat' >> ${3}_d.mri
echo 'images.extent.x=' $images_ext_x >> ${3}_d.mri
echo 'images.extent.y=' $images_ext_y >> ${3}_d.mri
echo 'images.extent.z=' $images_ext_z >> ${3}_d.mri
echo 'images.extent.t=' $images_ext_t >> ${3}_d.mri
echo 'images.extent.v= 1' >> ${3}_d.mri
echo 'images.offset= 0' >> ${3}_d.mri
echo 'images.size=' `(ls -l ${3}_d.dat | awk '{ print $5 }')` >> ${3}_d.mri
if ( ${PVM_ARCH} == LINUX ) then
        echo 'images.little_endian= 1' >> ${3}_d.mri
endif

mri_destroy_dataset ${1}_s

# Restore missing data
mri_copy_chunk -chunk missing $1 ${3}_d

# Convert the output to type float, and throw out the double version
mri_type_convert -float ${3}_d ${3}_f
mri_destroy_dataset ${3}_d

# Come up with a mean dataset, and add it back into the residuals
# It's actually faster to do the permutation than to interpolate
# in time series order.
mri_subsample -d t -l 1 -mean $1 ${3}_m
mri_remap -order vxyzt ${3}_m
mri_interp -con -d t -len $images_ext_t ${3}_m ${3}_l
mri_destroy_dataset ${3}_m
mri_permute -memlimit 32000000 -order vtxyz ${3}_l ${3}_M 
mri_destroy_dataset ${3}_l
mri_rpn_math -out $3 '$1,$2,+' ${3}_f ${3}_M
mri_destroy_dataset ${3}_f
mri_destroy_dataset ${3}_M

# Permute the output back to image order
echo '#'`date`
mri_permute -memlimit 32000000 -chunk ${F_DETREND_CHUNK} -order vxyzt $3 $2
echo '#'`date`



