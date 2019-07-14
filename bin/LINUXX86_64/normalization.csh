#!/bin/csh -ef
# normalization.csh
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
echo '#$Id: normalization.csh,v 1.8 2007/06/22 00:57:01 welling Exp $'

set normscript = "generated_scripts/apply_normalize.gen_csh"
set coregscript = "generated_scripts/apply_coreg_inplane.gen_csh"

# Do we have the necessary info
if ( ! -e ${1}.mri ) then
  echo "Normalization prototype does not exist!"
  exit -1
endif
if ( ! -e $coregscript ) then
  echo "Coregistration information is missing!"
  exit -1
endif

# Check needed dirs
if ( ! -d coregistered_func) mkdir coregistered_func
if ( ! -d generated_scripts) mkdir generated_scripts

# Build the generated normalization script
normalize_prep.py -v --oscript $normscript $1

# Make separate counts by condition, generating a list of files to
# normalize in the process.
pushd stat > /dev/null
set meanList = ()
set stdvList = ()
set countsList = ()
set rawMeanList = ()
set done = 0
set count = 1
while ( ! $done ) 
  if ( -e Mean_${count}.mri ) then
    set rawMeanList = ( $rawMeanList Mean_${count}.mri )
    @ count = $count + 1
  else
    set done = 1
  endif
end
if ( ${#rawMeanList} < 10 ) then
  foreach ds ( $rawMeanList )
    set meanList = ( $meanList ${ds:r} )
    set counts_name = `echo ${ds:r} | sed 's/Mean/Counts/g'`
    set countsList = ( $countsList $counts_name )
    set stdv_name = `echo ${ds:r} | sed 's/Mean/Stdv/g'`
    set stdvList = ( $stdvList $stdv_name )
    mri_copy_chunk -chunk counts -chunk_out images -chunk_file .dat $ds nrm_tmp_1
    mri_remap -order xyz nrm_tmp_1
    mri_interp -constant -d x -len 64 nrm_tmp_1 nrm_tmp_2
    mri_interp -constant -d y -len 64 nrm_tmp_2 nrm_tmp_3
    mri_rpn_math -out $counts_name '$2' $ds nrm_tmp_3
    mri_destroy_dataset nrm_tmp_1
    mri_destroy_dataset nrm_tmp_2
    mri_destroy_dataset nrm_tmp_3
  end
else
  echo "### There are ${#rawMeanList} condition means and stdvs- too many! ###"
endif
popd > /dev/null

# Normalize appropriate files
set dsList = ( GrandMean GrandStdv $meanList $stdvList $countsList )
foreach ds ( $dsList )
  set stat_ds = stat/${ds:r}
  set aligned_ds = coregistered_func/aligned_${ds:r}
  set norm_ds = normalized_func/norm_${ds:r}
  if ( ! -e ${stat_ds}.mri ) then
    echo "$stat_ds does not exist; skipping!"
    continue
  else
    echo "Processing $stat_ds"
  endif
  set needs_align = 0
  if ( ! -e ${aligned_ds}.mri ) then
    set needs_align = 1
  else
    set test_newer = \
      `test_in_subshell.csh test_newer.csh ${aligned_ds}.mri ${stat_ds}.mri`
    if ( $test_newer ) set needs_align = 1
    set test_newer = \
      `test_in_subshell.csh test_newer.csh ${aligned_ds}.mri ${coregscript}`
    if ( $test_newer ) set needs_align = 1
  endif
  if ( $needs_align ) then
    echo "Coregistering $ds"
    $coregscript $stat_ds $aligned_ds tmp_normalization > /dev/null
    mri_destroy_dataset tmp_normalization
  endif
  set needs_norm = 0
  if ( ! -e ${norm_ds}.mri ) then
    set needs_norm = 1
  else
    set test_newer = \
      `test_in_subshell.csh test_newer.csh ${norm_ds}.mri ${aligned_ds}.mri`
    if ( $test_newer ) set needs_norm = 1
    set test_newer = \
      `test_in_subshell.csh test_newer.csh ${norm_ds}.mri ${normscript}`
    if ( $test_newer ) set needs_norm = 1
  endif
  if ( $needs_norm ) then
    echo "Normalizing $ds"
    $normscript $aligned_ds $norm_ds > /dev/null
  endif
end

