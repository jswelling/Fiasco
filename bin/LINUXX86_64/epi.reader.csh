#!/bin/csh -efx
# epi.reader.csh
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
echo '$Id: epi.reader.csh,v 1.29 2006/03/10 01:23:04 welling Exp $'
#
echo '#'`date`$0

set outfile = $argv[$#argv]

# Read data, based on pfile if F_READER_DIMS is empty.  Each specified
# input file is read separately and finally concatenated.

@ nextToLast = $#argv - 1
set rawfiles = ( $argv[1-$nextToLast] )
set count = 0
set infiles = ()
set rdr_read_fname = ${outfile}_${F_READER_DIMORDER}
set dim_words = `echo $F_READER_DIMS`
foreach rawfile ( $rawfiles )
  #cat $1.Z | uncompress > trash
  set thisInFile = tmp_${count}
  if ( $#dim_words ) then

    set smart_dims = ${F_READER_VECLEN}:`echo $F_READER_DIMS | sed 's/ /:/g'`
    smartreader -verbose \
        -input $rawfile -dims ${smart_dims} -out $thisInFile \
	-offset $F_READER_OFFSET -reorder $F_READER_REORDER \
	-skip $F_READER_SKIP -sliceskip $F_READER_SSKIP -type $F_READER_TYPE \
	-dataorder $F_READER_DIMORDER \
	-tag $F_DIR $F_READER_ENDIAN $F_READER_OPTS

  else

    smartreader -verbose -input $rawfile -out $thisInFile \
	-tag $F_DIR -autoscale $F_READER_OPTS

  endif
  #rm trash

  set tdim = `mri_printfield -field images.extent.t $thisInFile`
  echo $count | mri_from_ascii -order t -len 1 -chunk acq_blocks acq
  mri_interp -d t -len $tdim -constant acq acq_stretch
  mri_destroy_dataset acq
  mri_copy_chunk -chunk acq_blocks acq_stretch $thisInFile
  mri_destroy_dataset acq_stretch
  
  set infiles = ( $infiles $thisInFile )
  @ count = $count + 1
end
if ( $count == 1 ) then
  mri_copy_dataset $infiles $rdr_read_fname
else
  mri_paste -d t -outfile $rdr_read_fname $infiles
endif
foreach ds ( $infiles )
  mri_destroy_dataset $ds
end

# EPI datasets requiring ramp resampling use the dimension q in place
# of x.
set dimstr = `mri_printfield -field images.dimensions ${rdr_read_fname}`
if ( $dimstr =~ *qy* ) then
  set x = q
else
  set x = x
endif

# If this was read from multiple blocks of several images each, the
# final dimensions may be ta (because of the behavior of smartreader
# with the -multi option).  We need to remap these datasets.
set tdim = `mri_printfield -field images.extent.t ${rdr_read_fname}`
set zdim = `mri_printfield -field images.extent.z ${rdr_read_fname}`
if ( $dimstr =~ *ta ) then
  set adim = `mri_printfield -field images.extent.a ${rdr_read_fname}`
  @ tdim_new = $tdim * $adim
  set tdim = $tdim_new
  set dimstr_new = `echo $dimstr | sed 's/ta/t/g'`
  set dimstr = $dimstr_new
  set extents = \
    `echo $dimstr | sed "s/t/${tdim}/g" | sed 's/[^0123456789]/:/g'`
  mri_remap -order $dimstr -length $extents ${rdr_read_fname}
endif

# If there is no v dimension, add a trivial one.
if ( $dimstr !~ v* ) then
  set dimstr_new = v${dimstr}
  set dimstr = $dimstr_new
  mri_remap -order $dimstr ${rdr_read_fname}
endif
# Note the v dimension
set vdim = `mri_printfield -field images.extent.v ${rdr_read_fname}`

# Check that the user has specified a number of images and slices, etc,
# equal to the number actually found
if (! { dataset_matches_env.py -v ${rdr_read_fname} } ) then
  echo "***ERROR*** environment variables conflict with file info"
endif

# Chop in X, Y, or both as required
set xchop = \
  `mri_printfield -field images.xchop -nofail ${rdr_read_fname}`
set ychop = \
  `mri_printfield -field images.ychop -nofail ${rdr_read_fname} `
set rdr_chop_fname = ${outfile}_chop
if ( dummy${xchop} == dummy1 ) then
  if ( dummy${ychop} == dummy1 ) then
    mri_rpn_math -out ${rdr_chop_fname} \
      '$1,1,-1,$ydim,$y,1,+,-,2,%,if_keep,1,-1,$'$x',2,%,if_keep,*,*' \
      ${rdr_read_fname}
    mri_setfield -field images.ychop -delete ${rdr_chop_fname}
  else
    mri_rpn_math -out ${rdr_chop_fname} \
      '$1,1,-1,$'$x',2,%,if_keep,*' ${rdr_read_fname}
  endif
  mri_destroy_dataset ${rdr_read_fname}
  mri_setfield -field images.xchop -delete ${rdr_chop_fname}
else
  if ( dummy${ychop} == dummy1 ) then
    mri_rpn_math -out ${rdr_chop_fname} \
      '$1,1,-1,$ydim,$y,1,+,-,2,%,if_keep,*' ${rdr_read_fname}
    mri_destroy_dataset ${rdr_read_fname}
    mri_setfield -field images.ychop -delete ${rdr_chop_fname}
  else
    set rdr_chop_fname = ${rdr_read_fname}
  endif
endif

if ( $?F_READER_SPACE ) then
  set rdr_space = $F_READER_SPACE
else
  if ( $vdim == 1 ) then
    set rdr_space = i
  else
    set rdr_space = k
  endif
endif

if ( $dimstr != "v${x}yzt" ) then
  echo 'Permuting data from ' $F_READER_DIMORDER ' to v'$x'yzt'
  set rdr_permute_fname = ${outfile}_${rdr_space}
  mri_permute -memlimit 32000000 -order v${x}yzt ${rdr_chop_fname}.mri \
          ${rdr_permute_fname}.mri
  mri_destroy_dataset ${rdr_chop_fname}
else
  set rdr_permute_fname = ${rdr_chop_fname}   
endif

# We need to flip real images to k-space for the convenience of
# the downstream processing stream.  If we do so, however, we can
# skip deghosting.
#
# F_READER_SPACE is not always correctly set, so don't rely on it.
if ( $rdr_space == "i" || $vdim == 1) then
  echo "Converting fdata from image space to k-space"
  set rdr_kspace_fname = ${outfile}_kspace
  mri_fft -d ${x}y -fwd -cpx ${rdr_permute_fname} ${rdr_kspace_fname}
  mri_setfield -field images.needs_deghost -value 0 ${rdr_kspace_fname}
  mri_destroy_dataset ${rdr_permute_fname}
else
  set rdr_kspace_fname = ${rdr_permute_fname}
  mri_setfield -field images.needs_deghost -value 1 ${rdr_kspace_fname}
endif

# Do scan reordering if required
set reorder = \
  `mri_printfield -field images.reorder -nofail ${rdr_kspace_fname}`
if (${#reorder} != 0) then
  if (${reorder} != 0) then
    set rdr_reorder_fname = ${outfile}_reorder
    mri_scan_fold -v ${rdr_kspace_fname} ${rdr_reorder_fname}
    mri_destroy_dataset ${rdr_kspace_fname}
    mri_setfield -field images.reorder -value 0 ${rdr_reorder_fname}
  else
    set rdr_reorder_fname = ${rdr_kspace_fname}
  endif
else
  set rdr_reorder_fname = ${rdr_kspace_fname}
endif

# Do autoscale if required
set rdr_out_fname = ${outfile}
set xdim = \
  `mri_printfield -field images.extent.${x} ${rdr_reorder_fname}`
set ydim = `mri_printfield -field images.extent.y ${rdr_reorder_fname}`
set autoscale = \
  `mri_printfield -field images.autoscale -nofail ${rdr_reorder_fname}`
if (${#autoscale} != 0) then
  if (${autoscale} != 0) then
    set autoscale_range = \
      `mri_printfield -field images.autoscale_range -nofail ${rdr_reorder_fname}`
    if (${#autoscale_range} == 0) then
      set autoscale_range = 1
    endif
    mri_subset -d t -l 1 -s 0 ${rdr_reorder_fname} rdr_autoscale_tmp
    mri_fft -d ${x}y -mod rdr_autoscale_tmp rdr_autoscale_tmp2
    mri_destroy_dataset rdr_autoscale_tmp
    @ wdim = ${xdim} * ${ydim} * ${zdim}
    mri_remap -order w -length ${wdim} rdr_autoscale_tmp2
    mri_subsample -d w -l 1 -max rdr_autoscale_tmp2 rdr_autoscale_tmp3
    mri_destroy_dataset rdr_autoscale_tmp2
    set val = `mri_rpn_math -out rdr_autoscale_tmp4 '0,$1,1,if_print_1' rdr_autoscale_tmp3`
    mri_destroy_dataset rdr_autoscale_tmp3
    mri_destroy_dataset rdr_autoscale_tmp4
    mri_rpn_math -out ${rdr_out_fname} '$1,'${autoscale_range}','${val}',/,*' \
      ${rdr_reorder_fname}
    mri_destroy_dataset ${rdr_reorder_fname}
    mri_setfield -field images.autoscale -value 0 ${rdr_out_fname}
  else
    mri_copy_dataset ${rdr_reorder_fname} ${rdr_out_fname}
    mri_destroy_dataset ${rdr_reorder_fname}
  endif
else
  mri_copy_dataset ${rdr_reorder_fname} ${rdr_out_fname}
  mri_destroy_dataset ${rdr_reorder_fname}
endif
