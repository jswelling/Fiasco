#!/bin/csh -efx
# smooth_detrend.csh
#/************************************************************
# *                                                          *
# *  Permission is hereby granted to any individual or       *
# *  institution for use, copying, or redistribution of      *
# *  this code and associated documentation, provided        *
# *  that such code and documentation are not sold for       *
# *  profit and the following copyright notice is retained   *
# *  in the code and documentation:                          *
# *     Copyright (c) 1995,1996 Department of Statistics,    *
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
# *  Original programming by Bill Eddy                       *
# ************************************************************/
#
echo '# '`date` $0
echo '# $Id: smooth_detrend.csh,v 1.8 2003/09/09 16:25:45 bakalj Exp $'
if(! -d stat)mkdir stat
if(! -d ps) mkdir ps
if (`mri_printfield -field ${F_DETREND_CHUNK}.dimensions $1` == "xyzt") then
	mri_remap -chunk ${F_DETREND_CHUNK} -order vxyzt $1
endif

# Permute to time order

mri_permute -memlimit 32000000 -chunk ${F_DETREND_CHUNK} -order vtxyz $1 $3

# Use a simple rule of thumb to figure out kernel width

set tdim = `mri_printfield -field ${F_DETREND_CHUNK}'.extent.t' $3`
@ band = ( $tdim / $F_SMDETREND_DNM )

# Build the smoothed version, and a mean dataset.
mri_smooth -d t -bandwidth $band $3 tmp_smooth_detrend
mri_subsample -d t -l 1 -mean $3 tmp_smooth_detrend2
mri_permute -memlimit 32000000 -order vxyzt tmp_smooth_detrend2 tmp_smooth_detrend3
mri_destroy_dataset tmp_smooth_detrend2

# Subtract smoothed from original
mri_rpn_math -out $2 '$1,$2,-' $3 tmp_smooth_detrend
mri_destroy_dataset tmp_smooth_detrend

# Reorder the result
mri_permute -memlimit 32000000 -chunk ${F_DETREND_CHUNK} -order vxyzt $2 $3 

# Add mean component back in.  We have to do this in vxyzt order to
# avoid creating yet another big temporary dataset.
mri_rpn_math -out $2 '$1,$2,+' $3 tmp_smooth_detrend3
mri_destroy_dataset tmp_smooth_detrend3

# Make a permuted copy for compatibility with detrend.csh
mri_permute -memlimit 32000000 -chunk ${F_DETREND_CHUNK} -order vtxyz $2 $3

echo '# '`date`' Smoothing detrend complete.'







