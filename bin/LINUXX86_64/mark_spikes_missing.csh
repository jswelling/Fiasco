#!/bin/tcsh -efx
# mark_spikes_missing.csh
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
# This file takes 3 inputs:
# 1) the input dataset
# 2) the output dataset
# 3) the spike count threshold above which to mark additional slices missing
#

echo '#'`date` $0
echo '#$Id: mark_spikes_missing.csh,v 1.3 2003/12/02 19:48:41 welling Exp $'

# Read some dimensions from the Pgh MRI file
set zextent = `mri_printfield -field 'images.extent.z' -nofail $1`
if (dummy$zextent == dummy) then
    set zextent = `mri_printfield -field 'samples.extent.z' $1`
    endif
set textent = `mri_printfield -field 'images.extent.t' -nofail $1`
if (dummy$textent == dummy) then
    set textent = `mri_printfield -field 'samples.extent.t' $1`
    endif
#echo textent = ${textent}, zextent = ${zextent}. 

# Get the missing chunk from the input data
if ( -e tmp_mark_missing.mri ) mri_destroy_dataset tmp_mark_missing
mri_copy_chunk -c missing $1 tmp_mark_missing
 
# Make a chunk with the spike data in it
mri_from_ascii -c missing -ord vtz -ind tz -l 1:${textent}:${zextent} \
    tmp2_mark_missing < `most_recent_parms.py despike`
mri_permute -order vzt -c missing tmp2_mark_missing tmp3_mark_missing 
mri_remap -order zt -c missing tmp3_mark_missing 

# Modify the missing data
mri_rpn_math -out tmp4_mark_missing -c missing '$1,1,$2,'$3',<,if_keep' \
  tmp_mark_missing tmp3_mark_missing
 
# Change back to characters
mri_type_convert -char -c missing tmp4_mark_missing tmp5_mark_missing
 
# Move new missing info into the output file
mri_copy_dataset $1 $2
mri_copy_chunk -c missing -replace tmp5_mark_missing $2
 
# Clean up
foreach fname (tmp*_mark_missing.mri)
  mri_destroy_dataset $fname
#end
 
# Write to summary file
set step = `depath.csh $0`
echo $step.$$ >> $F_SUMM_MISSING
count_missing.csh $2.mri >> $F_SUMM_MISSING


