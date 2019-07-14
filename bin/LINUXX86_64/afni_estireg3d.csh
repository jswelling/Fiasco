#!/bin/csh -efx
# afni_estireg3d.csh
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
# *                                                          *
# *   Original programming by S. Ali Imam                    *
# ************************************************************/
#
echo '#'`date` $0
echo '#$Id: afni_estireg3d.csh,v 1.5 2006/01/23 19:55:22 welling Exp $'
if (! -d par) mkdir par

# Check for help command
if ( $#argv >= 1 ) then
  if ( junk$argv[1] == junk-help ) then
    if ( $#argv >= 2 ) then
      scripthelp $0 $argv[2]
    else
      scripthelp $0
    endif
    exit
  endif
endif

# Default behavior is to align data to mean.
if ( ! ${?F_EST3D_ALIGN} ) then
    setenv F_EST3D_ALIGN "mean"
endif

set dimstr = `mri_printfield -field images.dimensions $1`
if ( $dimstr == xyzt ) then
  set perm_order = "txyz"
else
  set perm_order = "vtxyz"
endif

set xvox = \
  `mri_printfield -field images.voxel_spacing.x -nofail $1`
if (${#xvox} != 0)  then
  if ( $xvox != 0 ) then
    echo "x voxel already there"
  else
    mri_setfield -field images.voxel_spacing.x -value ${F_XVOXEL} $1
  endif
else
  mri_setfield -field images.voxel_spacing.x -value ${F_XVOXEL} $1
endif

set yvox = \
  `mri_printfield -field images.voxel_spacing.y -nofail $1`
if (${#yvox} != 0)  then
  if ( $yvox != 0 ) then
    echo "y voxel already there"
  else
    mri_setfield -field images.voxel_spacing.y -value ${F_YVOXEL} $1
  endif
else
  mri_setfield -field images.voxel_spacing.y -value ${F_YVOXEL} $1
endif

set zvox = \
  `mri_printfield -field images.voxel_spacing.z -nofail $1`
if (${#zvox} != 0)  then
  if ( $zvox != 0 ) then
    echo "z voxel already there"
  else
    mri_setfield -field images.voxel_spacing.z -value ${F_ZVOXEL} $1
  endif
else
  mri_setfield -field images.voxel_spacing.z -value ${F_ZVOXEL} $1
endif

#
# changing the target value into and afni dataset
# so we can use 3dvolreg on it.
#
pghtoafni -func $1 ${1}+orig

if ( $F_EST3D_ALIGN == "mean" ) then
#
# calculating a mean in order to transform into afni
# and perform needed tasks
#
  mri_permute -order $perm_order $1 ${1}_p
  mri_subsample -mean -d t -l 1 ${1}_p ${1}_mean
  mri_destroy_dataset ${1}_p
  mri_remap -order vxyzt ${1}_mean
    
  pghtoafni -anat ${1}_mean Mean1+orig
  set target = Mean1+orig
else if ( $F_EST3D_ALIGN == "median" ) then
#
#  calculating the median
#
  mri_permute -order $perm_order $1 ${1}_p
  mri_subsample -median -d t -l 1 ${1}_p ${1}_median
  mri_destroy_dataset ${1}_p
  mri_remap -order vxyzt ${1}_median
  
  pghtoafni -anat ${1}_median Median+orig
  set target = Median+orig
else
  pghtoafni -anat $F_EST3D_ALIGN Mean1+orig
  set target = Mean1+orig
endif

#
#  run afni's 3d volume registration
#
3dvolreg -base ${target} -dfile paramfile \
    -prefix aligned ${1}+orig

# Translate the output file if an output dataset was specified
if ( dummy$2 != dummy ) then
  smartreader -i data/aligned+orig.HEAD -out $2

  #
  # Getting rid of the extra files made by smartreader
  #
  set useless = ${2}.V*
  while ( $#useless != 0 )
    mri_delete_chunk -c afni.${useless[1]:e} ${useless[1]:r}
    rm $useless[1]
    if ( $#useless != 1 ) then
      set useless = ${2}.V*
    else 
      break
    endif
  end
endif

set xdim = `mri_printfield -fld images.extent.x $1`
set ydim = `mri_printfield -fld images.extent.y $1`

par_translate -i afni3d -o estireg3d -t ${F_NIMAGE} \
    -z ${F_NSLICE} -nx $xdim -ny $ydim -nz ${F_NSLICE} \
    -vx ${F_XVOXEL} -vy ${F_YVOXEL} -vz ${F_ZVOXEL} \
    paramfile par/${F_EST3D_PARMS}.$$

# We will pretend to be estireg3d, to better blend with our surroundings.
set step = estireg3d
echo "$step par/${F_EST3D_PARMS}.$$" >> $F_SUMM_INPUT

#
#  Cleaning up
#
if ( $F_EST3D_ALIGN == "median" ) then
  rm Median+orig*
else if ( $F_EST3D_ALIGN == "mean" ) then
  rm Mean1+orig*
else
  rm Mean1+orig*
endif

rm ${1}+orig*
rm paramfile
rm data/aligned+orig*
