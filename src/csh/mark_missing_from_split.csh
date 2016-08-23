#! /bin/csh -exf
# mark_missing_from_split.csh
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
# ************************************************************/
#
echo '#'`date` $0
echo '#$Id: mark_missing_from_split.csh,v 1.6 2005/06/05 06:52:22 welling Exp $'

# Parse split file if an output is missing, or if it is 
# more recent than an output
set split_new = $F_SPLIT_NEW
if (! -d ${split_new:h} ) mkdir ${split_new:h}
set split_cond = $F_SPLIT_COND
if (! -d ${split_cond:h} ) mkdir ${split_cond:h}
set do_split = 0
if (! -r $F_SPLIT_NEW || ! -r $F_SPLIT_COND ) then
  set do_split = 1
else
  set ages = `ls -t $F_SPLIT_FILE $F_SPLIT_NEW $F_SPLIT_COND`
  if ($ages[3] != $F_SPLIT_FILE) set do_split = 1
endif
if ($do_split) then
  intsplit -splitin $F_SPLIT_FILE -code $F_SPLIT_COND \
           -splitout $F_SPLIT_NEW $1
endif

# Snag relevant info from the input DS
set zdim = `mri_printfield -field 'images.extent.z' -nofail $1`
if (dummy$zdim == dummy) then
  set zdim = `mri_printfield -field 'samples.extent.z' $1`
endif
set tdim = `mri_printfield -field 'images.extent.t' -nofail $1`
if (dummy$tdim == dummy) then
  set tdim = `mri_printfield -field 'samples.extent.t' $1`
endif
set missingchunk = \
  "`mri_printfield -field missing -nofail $1`"
if ( dummy"${missingchunk}" == 'dummy[chunk]' ) then
  set has_missing = 1
else
  set has_missing = 0
endif

# Make up a dataset with the split information, shaped like the
# appropriate missing chunk
mri_from_ascii -order tz -l ${tdim}:${zdim} -c missing -ind tz \
  mark_missing_tmp1 < $F_SPLIT_NEW
mri_permute -c missing -order zt mark_missing_tmp1 mark_missing_tmp2

# Insert or merge in the new missing information
if ( $has_missing ) then
  mri_type_convert -char -c missing mark_missing_tmp2 mark_missing_tmp3
  mri_rpn_math -c missing -out mark_missing_tmp4 '$1,0,==,1,$2,0,!=,if_keep' \
    mark_missing_tmp3 $1
  mri_type_convert -char -c missing mark_missing_tmp4 mark_missing_tmp5
  mri_copy_chunk -c missing -replace mark_missing_tmp5 $1
  foreach fname (mark_missing_tmp1 mark_missing_tmp2 mark_missing_tmp3 \
                 mark_missing_tmp4 mark_missing_tmp5)
    mri_destroy_dataset $fname
  end
else
  mri_type_convert -char -c missing mark_missing_tmp2 mark_missing_tmp3
  mri_rpn_math -c missing -out mark_missing_tmp4 '$1,0,==' mark_missing_tmp3
  mri_type_convert -char -c missing mark_missing_tmp4 mark_missing_tmp5
  mri_copy_chunk -c missing mark_missing_tmp5 $1
  foreach fname (mark_missing_tmp1 mark_missing_tmp2 mark_missing_tmp3 \
                 mark_missing_tmp4 mark_missing_tmp5)
    mri_destroy_dataset $fname
  end
endif

# Note the fact that things have been marked missing
set lclsum = $F_SUMM_MISSING
set sumdir = ${lclsum:h}
if (! -d $sumdir ) mkdir $sumdir
set step = `depath.csh $0`
echo $step.$$ >> $F_SUMM_MISSING
count_missing.csh $1.mri >> $F_SUMM_MISSING


