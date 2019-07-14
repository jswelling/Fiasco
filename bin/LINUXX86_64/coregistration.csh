#!/bin/csh -efx
# coregistration.csh
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
echo '#$Id: coregistration.csh,v 1.6 2005/04/12 22:50:40 welling Exp $'

# Do we have the necessary info?
# Beware: input names often contain wildcard characters.
set nonomatch
if ( ! ${?F_CRG_INPL_INPUT} ) then
  echo "No coregistration inplane input specified; skipping coregistration"
  exit
endif
if ( "dummy"${F_CRG_INPL_INPUT} == "dummy" ) then
  echo "No coregistration inplane input specified; skipping coregistration"
  exit
endif

if ( ! ${?F_CRG_ANAT_INPUT} ) then
  echo "No coregistration anatomy input specified; skipping coregistration"
  exit
endif

if ( "dummy"${F_CRG_ANAT_INPUT} == "dummy" ) then
  echo "No coregistration anatomy input specified; skipping coregistration"
  exit
endif

if ( ! ${?F_CRG_FUNC} ) then
  echo "No coregistration functional dataset specified; skipping coregistration"
  exit
endif

if ( "dummy"${F_CRG_FUNC} == "dummy" ) then
  echo "No coregistration functional dataset specified; skipping coregistration"
  exit
endif
unset nonomatch

# Check needed dirs
if ( ! -d anat) mkdir anat
if ( ! -d coregistered_anat) mkdir coregistered_anat
if ( ! -d coregistered_func) mkdir coregistered_func
if ( ! -d generated_scripts) mkdir generated_scripts

# Generate the inplane and axial anatomical datasets.  If the output
# already exists and is newer than the input, we'll use the existing
# output.  This allows the user to manually do skull stripping (a later
# step) and then re-run without having everything overwritten.
set input_dir = ${F_CRG_INPL_INPUT:h}
set test_newer = \
      `test_in_subshell.csh test_newer.csh anat/${F_CRG_INPL}.mri $input_dir `
if ( $test_newer ) then
  smartreader $F_CRG_INPL_OPTS -input $F_CRG_INPL_INPUT \
      -out anat/$F_CRG_INPL -verbose
endif
# There is a problem with some inplane anatomicals acquired before
# about 11/02 on the 3T GE at UPMC MR Center in which the slice planes 
# are skewed.  This can cause smartreader to fail to identify the set of 
# slices as a volume, and hence to fail to calculate volume corners 
# appropriately.  We patch over this by copying in volume corners from 
# the (theoretically identical) functional bounding volume if necessary.
set has_corners = 0
foreach cnr ( "tlf" "trf" "tlb" "trb" "blf" "brf" "blb" "brb" )
  set fld = `mri_printfield -fld "images."$cnr".0" -nofail anat/$F_CRG_INPL`
  if ( junk$fld != junk ) then
    set has_corners = 1
    break
  endif
end
if ( ! $has_corners ) then
  foreach cnr ( "tlf" "trf" "tlb" "trb" "blf" "brf" "blb" "brb" )
    set fld1 = `mri_printfield -fld "images."$cnr".0" -nofail $F_CRG_FUNC`
    set fld2 = `mri_printfield -fld "images."$cnr".1" -nofail $F_CRG_FUNC`
    set fld3 = `mri_printfield -fld "images."$cnr".2" -nofail $F_CRG_FUNC`
    if ( junk$fld1 != junk && junk$fld2 != junk && junk$fld3 != junk ) then
      mri_setfield -all012 -fld "images."$cnr \
        -value '" '$fld1':'$fld2':'$fld3'"' anat/$F_CRG_INPL
    endif
  end
endif

set input_dir = ${F_CRG_ANAT_INPUT:h}
set test_newer = \
      `test_in_subshell.csh test_newer.csh anat/${F_CRG_ANAT}.mri $input_dir `
if ( ! $test_newer ) then
  set test_newer = \
       `test_in_subshell.csh test_newer.csh anat/${F_CRG_ANAT}.mri anat/${F_CRG_INPL}.mri `
endif
if ( $test_newer ) then
  smartreader $F_CRG_ANAT_OPTS -input $F_CRG_ANAT_INPUT \
      -out anat/${F_CRG_ANAT}_tmp -verbose
  align_by_permutation.py -v --out anat/${F_CRG_ANAT} \
      anat/${F_CRG_ANAT}_tmp anat/${F_CRG_INPL}
  mri_destroy_dataset anat/${F_CRG_ANAT}_tmp
endif

# Strip the skulls from the anatomical data
set test_newer = \
      `test_in_subshell.csh test_newer.csh anat/${F_CRG_ANAT_STRIPPED}.mri anat/${F_CRG_ANAT}.mri`
if ( $test_newer ) then
  strip_skull.py -v anat/$F_CRG_ANAT anat/$F_CRG_ANAT_STRIPPED
endif
set test_newer = \
      `test_in_subshell.csh test_newer.csh anat/${F_CRG_INPL_STRIPPED}.mri anat/${F_CRG_INPL}.mri`
if ( $test_newer ) then
  strip_skull.py -v --like=anat/$F_CRG_ANAT_STRIPPED \
  anat/$F_CRG_INPL anat/$F_CRG_INPL_STRIPPED
endif

# Coregister
set arglist = ()
if ( $?F_CRG_INPL_ALG ) then
  set arglist = ( $arglist --inplanealg $F_CRG_INPL_ALG )
endif
if ( $?F_CRG_STRCT_ALG ) then
  set arglist = ( $arglist --structalg $F_CRG_STRCT_ALG )
endif
if ( $?F_CRG_WARP_ALG ) then
  set arglist = ( $arglist --warpalg $F_CRG_WARP_ALG )
endif
coregister.py $arglist \
    -v $F_CRG_FUNC anat/$F_CRG_INPL_STRIPPED anat/$F_CRG_ANAT_STRIPPED
mv apply_coreg*.gen_csh generated_scripts
mv apply_cowarp*.gen_csh generated_scripts

# Make the coregistered anatomical for normalization
generated_scripts/apply_coreg_struct_to_inplane.gen_csh \
    anat/$F_CRG_ANAT coreg_tmp1 coreg_tmp2
generated_scripts/apply_cowarp_inplane.gen_csh coreg_tmp1 \
    coregistered_anat/$F_CRG_ANAT_COWARPED
generated_scripts/apply_cowarp_inplane.gen_csh coreg_tmp2 \
    coregistered_anat/$F_CRG_ANAT_INPL
mri_destroy_dataset coreg_tmp1
mri_destroy_dataset coreg_tmp2

# Make the coregistered functional prototype, and snag end-to-end
# mutual information on the way.
generated_scripts/apply_coreg_inplane.gen_csh $F_CRG_FUNC \
    coregistered_func/$F_CRG_FUNC_COREG coreg_tmp3
mutual_information -estimates mi_tmp.par  -nbins 256 \
    coregistered_anat/$F_CRG_ANAT_INPL coreg_tmp3
set mutual_info = `grep -v '#' mi_tmp.par | awk '{print $2}'`
mri_destroy_dataset coreg_tmp3
rm mi_tmp.par

# Make an image to help the user verify coregistration worked as expected
setenv F_PS_TITLE3 "Mutual Information $mutual_info"
coregister_makeps.py -v --out ps/coreg_test.ps \
    --func coregistered_func/$F_CRG_FUNC_COREG \
    --strct coregistered_anat/$F_CRG_ANAT_INPL


