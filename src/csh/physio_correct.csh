#!/bin/csh -efx
# physio_correct.csh
#/************************************************************
# *                                                          *
# *  Permission is hereby granted to any individual or       *
# *  institution for use, copying, or redistribution of      *
# *  this code and associated documentation, provided        *
# *  that such code and documentation are not sold for       *
# *  profit and the following copyright notice is retained   *
# *  in the code and documentation:                          *
# *     Copyright (c) 1995 Department of Statistics,         *
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
echo '#'`date` $0
echo '#$Id: physio_correct.csh,v 1.19 2006/03/10 01:23:04 welling Exp $'
#
# This script takes four parameters.  The first is the input file
# name, the second is the output file name, the third is the file
# containing valid "missing" data, and the fourth is the
# file of physiology data.  If the fourth parameter is not given,
# rename the input to the output.
#
if (dummy$4 == dummy) then
  mri_copy_dataset $1 $2
else

# Check for existence of directories
  set cardio_f = $F_PHYS_CARDIO
  set resp_f = $F_PHYS_RESP
  set param_f = $F_PHYS_PARAM
  if ((dummy${2:h} != dummy) && (! -d ${2:h})) mkdir ${2:h}
  if ((dummy${cardio_f:h} != dummy) && (! -d ${cardio_f:h})) \
	mkdir ${cardio_f:h}
  if ((dummy${resp_f:h} != dummy) && (! -d ${resp_f:h})) \
	mkdir ${resp_f:h}
  if ((dummy${param_f:h} != dummy) && (! -d ${param_f:h})) \
	mkdir ${param_f:h}

# Step 1: subsample the physio data file, and shuffle to slicing order
  @ totaltimes = $F_NIMAGE * $F_NSLICE
  mri_subsample -d t -l $totaltimes -w ${F_PHYS_WINDOW} -${F_PHYS_SUBOP} $4 \
	physio_tmp_subsample
  set pattern = `mri_printfield -fld images.reorder_pattern $1 -nofail`
  if ( dummy$pattern == dummy ) then
    set pat_string = ""
  else
    set pat_string = "-reorder $pattern"
  endif
  mri_scan_fold -zdm $F_NSLICE ${pat_string} \
	physio_tmp_subsample physio_tmp_scan
  mri_subset -d v -l 1 -s 1 physio_tmp_scan $F_PHYS_RESP
  mri_subset -d v -l 1 -s 2 physio_tmp_scan $F_PHYS_CARDIO
  mri_destroy_dataset physio_tmp_subsample
  mri_destroy_dataset physio_tmp_scan

# Step 2: Prepare permuted physiology data, plus a factor to act as
# a surrogate for time.
  mri_permute -memlimit 32000000 -order tz -chunk physio \
	$F_PHYS_CARDIO ${F_PHYS_CARDIO}_p
  mri_permute -memlimit 32000000 -order tz -chunk physio \
	$F_PHYS_RESP ${F_PHYS_RESP}_p
  mri_counter -c physio -d t -l ${F_NIMAGE} physio_timer_tz
  mri_remap -c physio -o tz physio_timer_tz 
  mri_interp -d z -len ${F_NSLICE} -con physio_timer_tz physio_timer
  mri_destroy_dataset physio_timer_tz

# Step 3: Run the regression
  mri_permute -memlimit 32000000 -order vtxyz $1 ${1}_p
  mri_copy_chunk -chunk missing -replace $3 ${1}_p
  mri_glm -var -ssqr -output ${2}_p -est ${F_PHYS_PARAM}:images \
	${1}_p \
	${F_PHYS_RESP}_p:physio \
	${F_PHYS_CARDIO}_p:physio \
	physio_timer:physio
  mri_destroy_dataset ${F_PHYS_RESP}_p
  mri_destroy_dataset ${F_PHYS_CARDIO}_p
  mri_destroy_dataset ${1}_p
  mri_destroy_dataset physio_timer

# Step 4: Construct the output, with constant term restored
  mri_subset -d v -s 0 -l 2 ${F_PHYS_PARAM} physio_const
  mri_remap -c physio -o vxyzt physio_const
  mri_permute -memlimit 32000000 -order vxyzt ${2}_p ${2}_noconst
  mri_destroy_dataset ${2}_p
  mri_rpn_math -out $2 '$1,$2,+' ${2}_noconst physio_const
  mri_destroy_dataset ${2}_noconst
  mri_destroy_dataset physio_const

endif

