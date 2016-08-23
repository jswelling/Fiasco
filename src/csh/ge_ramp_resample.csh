#!/bin/csh -efx
# ge_ramp_resample.csh
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
echo '$Id: ge_ramp_resample.csh,v 1.2 2003/09/10 20:02:50 bakalj Exp $'
#
echo '#'`date`$0

# This script will be called with two parameters, the un-resampled
# file and the name to be given the resampled file.  

set infile = ${1}
set outfile = ${2}

echo "Ramp resampling ${infile} using the GE algorithm"

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

# Extract the chunk containing the GE ramp resampling table
set rampfile = ${infile}_ramp
set rampdim = \
  `mri_printfield -field rampsample.dimensions -nofail ${infile}`
if ( dummy${rampdim} == dummy ) then
  echo "Ramp sampling is implied, but no ramp chunk was provided!"
  exit -1
endif
mri_copy_chunk -chunk rampsample -chunk_out images -replace ${infile} ${rampfile}

# Reorder data, apply the resampling matrix, and reorder it back
mri_permute -order ${perm_dimstr} ${infile} ${infile}_p 
mri_delete_chunk -chunk rampsample ${infile}_p 
mri_remap -order vyzqt -length ':::'${qdim}':' ${infile}_p 
time mri_matmult -v -out ${outfile}_p ${infile}_p ${rampfile}
mri_destroy_dataset ${infile}_p
mri_permute -order ${resamp_dimstr} ${outfile}_p ${outfile}
mri_destroy_dataset ${outfile}_p


