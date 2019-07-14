#! /bin/csh -exf
# phase_b0_homo_correct.csh.csh
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
echo '#$Id: phase_b0_homo_correct.csh,v 1.3 2003/09/09 16:25:45 bakalj Exp $'

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

#
# Make some scratch space
#
if ( ${?F_TEMP} ) then
  set tmpdir = ${F_TEMP}/phase_b0_homo_correct_$$
else
  set tmpdir = ./phase_b0_homo_correct_$$
endif
if (! -e $tmpdir) mkdir $tmpdir
set homedir = $PWD

mri_copy_chunk -chunk sample_kxloc -chunk_out images $1 $tmpdir/kxloc
mri_copy_chunk -chunk sample_kyloc -chunk_out images $1 $tmpdir/kyloc
mri_copy_chunk -chunk missing -chunk_out images $1 $tmpdir/missing

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
mri_delete_chunk -chunk missing phases
mri_delete_chunk -chunk missing mags
mri_delete_chunk -chunk missing samples_origin

mri_remap -c samples -order vsczt samples_origin 

mri_rpn_math -out phases_adj -c samples '$1,0.5,pi,*,$p,*,-' phases


# Regress across samples.  Since mri_glm requires that the regression
# dimension be called t, we have to rename p to t and t to something else
# (specifically q) to do the regression.

mri_permute -order vpsctz -c samples phases_adj phases_adj_p 
mri_remap -c samples -order tscqz \
    -length ${pdim}:${sdim}:${cdim}:${tdim}:${zdim} phases_adj_p 
mri_subset -d q -l 1 -s 0 phases_adj_p proto
mri_remap -order tz -c samples proto
mri_rpn_math -out lin_factor -c samples '$t,'${origin_sample}',-' proto
mri_rpn_math -out quad_factor -c samples '$1,dup,*' lin_factor

mri_glm -output phase_residuals_p -est param_p:images \
  phases_adj_p:samples lin_factor:samples quad_factor:samples

mri_remap -order vpsctz -chunk samples \
    -length 1:${pdim}:${sdim}:${cdim}:${tdim}:${zdim} phase_residuals_p 
mri_permute -order vpsczt -c samples phase_residuals_p phase_residuals 

mri_remap -order vsctz -length 3:${sdim}:${cdim}:${tdim}:${zdim} param_p 
mri_permute -order vsczt param_p param
foreach offset ( 0 1 2 )
  mri_subset -d v -l 1 -s ${offset} param param${offset}
  mri_smooth -d t -bandwidth $F_PHB0HOMO_BAND param${offset} param${offset}_sm
end

mri_rpn_math -out phase_corrected -c samples '$1,0.5,pi,*,$p,*,+' \
    phase_residuals
mri_rpn_math -out result_r -c samples '$1,$2,cos,*' mags phase_corrected
mri_rpn_math -out result_i -c samples '$1,$2,sin,*' mags phase_corrected

# And now we go home and construct the output
cd $homedir
mri_paste -d v -out $2 $tmpdir/result_r $tmpdir/result_i
mri_copy_chunk -chunk images -chunk_out sample_kxloc $tmpdir/kxloc $2
mri_copy_chunk -chunk images -chunk_out sample_kyloc $tmpdir/kyloc $2
mri_copy_chunk -chunk images -chunk_out missing $tmpdir/missing $2

# Produce summary information
if(! -d par) mkdir par
set step = `depath.csh $0`
echo '##Format: order:index_tz type:raw names:(b0_const,b0_lin,b0_quad)' \
    > par/$F_PHB0HOMO_RWPARMS.$$
mri_rpn_math -out $tmpdir/junk '0,$t,$z,$1,$2,$3,1,if_print_5' \
    $tmpdir/param0 $tmpdir/param1 $tmpdir/param2 \
    >> par/$F_PHB0HOMO_RWPARMS.$$
echo "$step par/${F_PHB0HOMO_RWPARMS}.$$" >> $F_SUMM_INPUT
echo '##Format: order:index_tz type:filtered names:(b0_const,b0_lin,b0_quad)' \
    > par/$F_PHB0HOMO_PARMS.$$
mri_rpn_math -out $tmpdir/junk '0,$t,$z,$1,$2,$3,1,if_print_5' \
    $tmpdir/param0_sm $tmpdir/param1_sm $tmpdir/param2_sm \
    >> par/$F_PHB0HOMO_PARMS.$$
echo "$step par/${F_PHB0HOMO_PARMS}.$$" >> $F_SUMM_INPUT

# Clean up
#rm -r $tmpdir

