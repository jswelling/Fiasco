#! /bin/csh -exf
# phase_lock.csh
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
echo '#$Id: phase_lock.csh,v 1.9 2003/11/24 21:55:13 welling Exp $'

################
# Notes:
# 
################

set imagechunk = "`mri_printfield -fld images -nofail $1 `"
if ( dummy"${imagechunk}" == 'dummy[chunk]' ) then
  set mainchunk = 'images'
else
  set mainchunk = 'samples'
  set sampleschunk = "`mri_printfield -fld samples -nofail $1 `"
  if ( dummy"${sampleschunk}" != 'dummy[chunk]' ) then
    echo "$0 : input file does not contain images or samples!"
    exit -1
  endif
endif
set missingchunk = "`mri_printfield -fld missing -nofail $1`"
if ( dummy"${missingchunk}" == 'dummy[chunk]' ) then
  set has_missing = 1
endif

# Find out some dimensions
set dimstr = `mri_printfield -field ${mainchunk}.dimensions $1`
set vdim = `mri_printfield -field ${mainchunk}.extent.v $1`
set zdim = `mri_printfield -field ${mainchunk}.extent.z $1`
set tdim = `mri_printfield -field ${mainchunk}.extent.t $1`

if ( ${mainchunk} == 'images' ) then
  if ( ${dimstr} != 'vxyzt' ) then
    echo "$0 can only handle images chunk in vxyzt order!"
    exit -1
  endif
  set xdim = `mri_printfield -field ${mainchunk}.extent.x $1`
  set ydim = `mri_printfield -field ${mainchunk}.extent.y $1`
  @ pdim = ( $xdim * $ydim )
  set sdim = 1
  set cdim = 1
else
  if ( ${dimstr} != 'vpsczt' ) then
    echo "$0 can only handle samples chunk in vpsczt order!"
    exit -1
  endif
  set pdim = `mri_printfield -field ${mainchunk}.extent.p $1`
  set sdim = `mri_printfield -field ${mainchunk}.extent.s $1`
  set cdim = `mri_printfield -field ${mainchunk}.extent.c $1`
endif

#
# Make some scratch space
#
if ( ${?F_TEMP} ) then
  set tmpdir = ${F_TEMP}/phase_lock_$$
else
  set tmpdir = ./phase_lock_$$
endif
if (! -e $tmpdir) mkdir $tmpdir

#
# Find the sample at the origin.  For EPI we assume it's just the
# center sample; for spiral we check the trajectory.
#
if ( ${mainchunk} == 'images' ) then
  @ pcenter = ( ( $ydim / 2 ) * $xdim + ( $xdim / 2 ) )
else
  mri_copy_chunk -chunk sample_kxloc -chunk_out images $1 $tmpdir/kxloc
  mri_copy_chunk -chunk sample_kyloc -chunk_out images $1 $tmpdir/kyloc

  set dz = 0
  while ( $dz < $zdim )
    set dc = 0
    while ( $dc < $cdim )
      set ds = 0
      while ( $ds < $sdim )
        set whichp = `mri_rpn_math -out $tmpdir/junk '0,$p,$s,'$ds',==,$c,'$dc',==,*,$z,'$dz',==,*,$1,0,==,*,$2,0,==,*,if_print_1' $tmpdir/kxloc $tmpdir/kyloc | tail -1 `
        if ( ${?origin_sample} ) then
          if ( $origin_sample != $whichp ) then
            echo "##ERROR##: k-space sample locations which vary by " \
                 "slice, coil, or shot are not supported."
            exit -1
          endif
        else
          set pcenter = $whichp
        endif
        @ ds = $ds + 1
      end
      @ dc = $dc + 1
    end
    @ dz = $dz + 1
  end
  echo "origin sample is " $pcenter
endif

# We'll need a scratch copy of the complex data in time series order.
if ( ${mainchunk} == 'images' ) then
  mri_permute -order vztxy $1 $tmpdir/in_p 
  mri_remap -order vztpsc -length 2:${zdim}:${tdim}:${pdim}:${sdim}:${cdim} \
    $tmpdir/in_p 
else
  mri_permute -chunk samples -order vztpsc $1 $tmpdir/in_p 
endif
mri_complex_to_scalar -chunk ${mainchunk} -mag $tmpdir/in_p $tmpdir/raw_mags
mri_complex_to_scalar -chunk ${mainchunk} -phu $tmpdir/in_p $tmpdir/raw_phases
foreach dset ( $tmpdir/raw_mags $tmpdir/raw_phases )
  if ( $has_missing ) mri_delete_chunk -chunk missing $dset
  if ( $mainchunk == 'samples' ) then
    mri_delete_chunk -chunk sample_kxloc $dset
    mri_delete_chunk -chunk sample_kyloc $dset
  endif
end
mri_destroy_dataset $tmpdir/in_p

if ( $F_PHLOCK_ALG == "weighted" ) then
# Construct an average phase for each slice, using the magnitude
# at each voxel as a weight.
  mri_rpn_math -out $tmpdir/product -chunk ${mainchunk} '$1,$2,*' \
      $tmpdir/raw_mags $tmpdir/raw_phases
  mri_permute -chunk ${mainchunk} -order vpsczt \
    $tmpdir/product $tmpdir/product_p 
  mri_permute -chunk ${mainchunk} -order vpsczt \
    $tmpdir/raw_mags $tmpdir/raw_mags_p 
  mri_subsample -d p -l 1 -sum \
     $tmpdir/product_p $tmpdir/prod_sum 
  mri_subsample -d p -l 1 -sum \
     $tmpdir/raw_mags_p $tmpdir/mag_sum 
  mri_rpn_math -chunk ${mainchunk} -out $tmpdir/rephase_raw '$1,$2,/' \
      $tmpdir/prod_sum $tmpdir/mag_sum
  mri_remap -chunk ${mainchunk} -order vztpsc $tmpdir/rephase_raw 
else if ( ${F_PHLOCK_ALG} == "center" ) then
# The phase of the center voxel is used.
  mri_subset -d p -l 1 -s $pcenter \
    $tmpdir/raw_phases $tmpdir/rephase_raw
else
  echo "Unknown algorithm ${F_PHLOCK_ALG} requested in $0!"
  exit -1
endif

# Smooth the estimates
mri_smooth -d t -bandwidth $F_PHLOCK_BAND \
    $tmpdir/rephase_raw $tmpdir/rephase

# Construct output with corrected phases.  We do a little song and dance 
# with an extra pair of remaps and a permute to cut down on pasting time.
mri_rpn_math -out $tmpdir/result_p_r -chunk ${mainchunk} '$1,$2,$3,-,cos,*' \
    $tmpdir/raw_mags $tmpdir/raw_phases $tmpdir/rephase
mri_rpn_math -out $tmpdir/result_p_i -chunk ${mainchunk} '$1,$2,$3,-,sin,*' \
    $tmpdir/raw_mags $tmpdir/raw_phases $tmpdir/rephase
mri_destroy_dataset $tmpdir/raw_mags
mri_destroy_dataset $tmpdir/raw_phases
mri_remap -order ztpscv -chunk ${mainchunk} $tmpdir/result_p_r 
mri_remap -order ztpscv -chunk ${mainchunk} $tmpdir/result_p_i
mri_paste -d v -out $tmpdir/result_p_p \
    $tmpdir/result_p_r $tmpdir/result_p_i
mri_destroy_dataset $tmpdir/result_p_r
mri_destroy_dataset $tmpdir/result_p_i
mri_permute -order vztpsc -chunk ${mainchunk} \
    $tmpdir/result_p_p $tmpdir/result_p 
mri_destroy_dataset $tmpdir/result_p_p
mri_permute -order vpsczt -chunk ${mainchunk} $tmpdir/result_p $tmpdir/result
mri_destroy_dataset $tmpdir/result_p
mri_copy_dataset $1 $2
mri_delete_chunk -chunk $mainchunk $2
mri_copy_chunk -chunk $mainchunk $tmpdir/result $2
mri_destroy_dataset $tmpdir/result
if ( ${mainchunk} == 'images' ) then
  mri_remap -order vxyzt -length 2:${xdim}:${ydim}:${zdim}:${tdim} $2
else
  mri_remap -chunk samples -order vpsczt \
    -length 2:${pdim}:${sdim}:${cdim}:${zdim}:${tdim} $2
endif

# Produce summary information
if(! -d par) mkdir par
set step = `depath.csh $0`
echo '##Format: order:index_tz type:raw names:(phase_lock)' \
    > par/$F_PHLOCK_RWPARMS.$$
mri_rpn_math -out $tmpdir/junk -chunk $mainchunk '0,$t,$z,$1,1,if_print_3' \
    $tmpdir/rephase_raw >> par/$F_PHLOCK_RWPARMS.$$
echo "$step par/${F_PHLOCK_RWPARMS}.$$" >> $F_SUMM_INPUT
echo '##Format: order:index_tz type:filtered names:(phase_lock)' \
    > par/$F_PHLOCK_PARMS.$$
mri_rpn_math -out $tmpdir/junk -chunk $mainchunk '0,$t,$z,$1,1,if_print_3' \
    $tmpdir/rephase >> par/$F_PHLOCK_PARMS.$$
echo "$step par/${F_PHLOCK_PARMS}.$$" >> $F_SUMM_INPUT

# Clean up
rm -r $tmpdir

