#! /bin/csh -exf
# phase_undrift.csh
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
echo '#$Id: phase_undrift.csh,v 1.7 2003/09/09 16:25:45 bakalj Exp $'

# Can we handle this dataset?
set dimstr = `mri_printfield -field samples.dimensions $1`
if ( $dimstr != "vpsczt" ) then
  echo "# Error: ${0}: input dataset $1 is not vpsczt!"
  exit -1
endif

# Find out some dimensions
set vdim = `mri_printfield -field samples.extent.v $1`
set pdim = `mri_printfield -field samples.extent.p $1`
set sdim = `mri_printfield -field samples.extent.s $1`
set cdim = `mri_printfield -field samples.extent.c $1`
set zdim = `mri_printfield -field samples.extent.z $1`
set tdim = `mri_printfield -field samples.extent.t $1`
@ samples_per_image = ( $pdim * $sdim * $cdim * $zdim )
set TR = `mri_printfield -field samples.tr $1`
set samp_time = `mri_printfield -field samples.samp_time $1`

#
# Make some scratch space
#
if ( ${?F_TEMP} ) then
  set tmpdir = ${F_TEMP}/phase_undrift_$$
else
  set tmpdir = ./phase_undrift_$$
endif
if (! -e $tmpdir) mkdir $tmpdir
set homedir = $PWD

mri_copy_chunk -chunk sample_kxloc -chunk_out images $1 $tmpdir/kxloc
mri_copy_chunk -chunk sample_kyloc -chunk_out images $1 $tmpdir/kyloc

#
# Look through the trajectory, finding the last sample which claims to
# be taken at the origin.
#
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
               "slice are not supported."
          exit -1
        endif
      else
        set origin_sample = $whichp
      endif
      @ ds = $ds + 1
    end
    @ dc = $dc + 1
  end
  @ dz = $dz + 1
end

echo "origin sample is " $origin_sample

mri_complex_to_scalar -phu -c samples $1 $tmpdir/phases
mri_complex_to_scalar -mag -c samples $1 $tmpdir/mags
mri_subset -d p -l 1 -s $origin_sample $1 $tmpdir/samples_origin
cd $tmpdir
mri_delete_chunk -chunk sample_kxloc phases
mri_delete_chunk -chunk sample_kyloc phases
mri_delete_chunk -chunk sample_kxloc mags
mri_delete_chunk -chunk sample_kyloc mags
mri_delete_chunk -chunk sample_kxloc samples_origin 
mri_delete_chunk -chunk sample_kyloc samples_origin 

mri_remap -c samples -order vsczt samples_origin 
mri_permute -c samples -order vtscz samples_origin samples_origin_p 
mri_complex_to_scalar -phu -c samples samples_origin_p phases_origin_u
mri_smooth -d t -bandwidth $F_PHUND_BAND phases_origin_u phases_smooth

# Make estimates of the change of phase per image, by doing numerical
# derivatives in an amazingly obscure way
mri_smooth -smoother_type shift -bdw 1 -d t phases_smooth phases_p1
mri_smooth -smoother_type shift -bdw 2 -d t phases_smooth phases_p2
mri_smooth -smoother_type shift -bdw -1 -d t phases_smooth phases_m1
mri_smooth -smoother_type shift -bdw -2 -d t phases_smooth phases_m2

mri_rpn_math -out phases_ctr -c samples '$2,$1,-,0.5,*' phases_m1 phases_p1
mri_rpn_math -out phases_fwd2 -c samples '$2,$1,-,1.0,*' phases_p1 phases_p2
mri_rpn_math -out phases_bck -c samples '$2,$1,-,1.0,*' phases_m1 phases_smooth

mri_rpn_math -out phase_change -c samples \
  '$1,$2,$t,0,==,if_keep,$3,$t,$tdim,1,-,==,if_keep' \
  phases_ctr phases_fwd2 phases_bck

# Set up to apply the correction.
mri_permute -c samples -order vsczt phase_change phase_change_p
mri_remap -chunk samples -order vscztp phase_change_p 
mri_permute -c samples -order vscztp phases phases_p
mri_permute -c samples -order vsczt phases_smooth phases_smooth_p
mri_remap -chunk samples -order vscztp phases_smooth_p 

mri_rpn_math -out phase_change_per_sample -c samples \
  '$1,'$samp_time',*,'$TR',/' phase_change_p

# Here we actually apply the correction
mri_rpn_math -out phases_corrected_p -c samples '$1,$2,$p,*,$3,+,-' \
  phases_p phase_change_per_sample phases_smooth_p

# And now we construct the output
mri_permute -c samples -order vpsczt phases_corrected_p phases_corrected
mri_rpn_math -out samples_r -c samples '$1,$2,cos,*' mags phases_corrected
mri_rpn_math -out samples_i -c samples '$1,$2,sin,*' mags phases_corrected

cd $homedir

# Final assembly of the output
set missingchunk = "`mri_printfield -fld missing -nofail $tmpdir/samples_r `"
if ( dummy"${missingchunk}" == 'dummy[chunk]' ) then
  mri_delete_chunk -chunk missing $tmpdir/samples_r
endif
set missingchunk = "`mri_printfield -fld missing -nofail $tmpdir/samples_i `"
if ( dummy"${missingchunk}" == 'dummy[chunk]' ) then
  mri_delete_chunk -chunk missing $tmpdir/samples_i
endif
mri_paste -d v -out $2 $tmpdir/samples_r $tmpdir/samples_i
mri_copy_chunk -chunk sample_kxloc $1 $2
mri_copy_chunk -chunk sample_kyloc $1 $2
set missingchunk = "`mri_printfield -fld missing -nofail $1`"
if ( dummy"${missingchunk}" == 'dummy[chunk]' ) then
  mri_copy_chunk -chunk missing $1 $2
endif

# Produce summary information
if(! -d par) mkdir par
set step = `depath.csh $0`
echo '##Format: order:index_tz type:raw names:(phase_roll)' \
    > par/$F_PHUND_RWPARMS.$$
mri_rpn_math -out $tmpdir/junk -c samples '0,$t,$z,$1,1,if_print_3' \
    $tmpdir/phases_origin_u >> par/$F_PHUND_RWPARMS.$$
echo "$step par/${F_PHUND_RWPARMS}.$$" >> $F_SUMM_INPUT
echo '##Format: order:index_tz type:filtered names:(phase_roll)' \
    > par/$F_PHUND_PARMS.$$
mri_rpn_math -out $tmpdir/junk -c samples '0,$t,$z,$1,1,if_print_3' \
    $tmpdir/phases_smooth >> par/$F_PHUND_PARMS.$$
echo "$step par/${F_PHUND_PARMS}.$$" >> $F_SUMM_INPUT

# Clean up
rm -r $tmpdir

