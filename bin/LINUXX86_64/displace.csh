#!/bin/csh -efx
# displace.csh
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
echo '#$Id: displace.csh,v 1.9 2003/09/09 16:25:45 bakalj Exp $'

if ( dummy${1} != dummy ) then
  set xdim2d = `mri_printfield -field images.extent.x $1`
  set ydim2d = `mri_printfield -field images.extent.y $1`
  set xdim3d = $xdim2d
  set ydim3d = $ydim2d
  set zdim3d = `mri_printfield -field images.extent.z $1`
else
  set xdim2d = $F_DISPLACE_XDIM
  set ydim2d = $F_DISPLACE_YDIM
  if ( $?F_CLIP2_WIDTH ) then
    set xdim3d = $F_CLIP2_WIDTH
  endif
  if ( $?F_CLIP2_HEIGHT ) then
    set ydim3d = $F_CLIP2_HEIGHT
  endif
  set zdim3d = $F_NSLICE
endif

set dimstr = `mri_printfield -field images.dimensions $1`
if ( $dimstr == xyzt ) then
  set perm_order = "txyz"
else
  set perm_order = "vtxyz"
endif

set wt_generated = 0
if ( $?F_DISPLACE_WEIGHT ) then
  if ( dummy${F_DISPLACE_WEIGHT} != dummy ) then
    if ( -e ${F_DISPLACE_WEIGHT}.mri ) then
      set weightstr = "-weight $F_DISPLACE_WEIGHT"
    else if ( ${F_DISPLACE_WEIGHT} == "mean" ) then
      mri_permute -order ${perm_order} ${1} ${1}_p
      mri_subsample -d t -l 1 -mean ${1}_p ${1}_mean
      mri_destroy_dataset ${1}_p
      mri_remap -order xyz ${1}_mean 
      set weightstr = "-weight ${1}_mean"
      set wt_generated = ${1}_mean
    else if ( ${F_DISPLACE_WEIGHT} == "median" ) then
      mri_permute -order ${perm_order} ${1} ${1}_p
      mri_subsample -d t -l 1 -median ${1}_p ${1}_median
      mri_destroy_dataset ${1}_p
      mri_remap -order xyz ${1}_median 
      set weightstr = "-weight ${1}_median"
      set wt_generated = ${1}_median
    else
      echo "# Error: cannot find displacement weight file $F_DISPLACE_WEIGHT"
      exit -1
    endif
  else
    set weightstr = ""
  endif
else
  set weightstr = ""
endif

if(! -d par)mkdir par
set estpar = `most_recent_parms.py estireg`
if ( dummy$estpar != dummy ) then
	displace -estimates par/$F_DISPLACE_OUTPUT.$$ \
		-xdimension $xdim2d -ydimension $ydim2d $weightstr \
		$estpar
    set step = `depath.csh $0`
    echo "$step par/$F_DISPLACE_OUTPUT.$$" >> $F_SUMM_INPUT
endif
set estsmpar = `most_recent_parms.py parsm`
if ( dummy$estsmpar != dummy ) then
	displace -estimates par/$F_DISPLACE_SMOUT.$$ \
		-xdimension $xdim2d -ydimension $ydim2d $weightstr \
		$estsmpar
    set step = `depath.csh $0`
    echo "sm$step par/$F_DISPLACE_SMOUT.$$" >> $F_SUMM_INPUT
endif
set estpar3d = `most_recent_parms.py estireg3d`
if ( dummy$estpar3d != dummy ) then
    if ( $?xdim3d && $?ydim3d && $?zdim3d ) then
      displace3d -e par/${F_DISP3D_RAW_PARMS}.$$ \
	-xdimension $xdim3d -ydimension $ydim3d -zdimension $zdim3d \
	-xvoxel ${F_XVOXEL} -yvoxel ${F_YVOXEL} \
        -zvoxel ${F_ZVOXEL} $weightstr $estpar3d 
      set step = `depath.csh $0`3d_raw
      echo "$step par/${F_DISP3D_RAW_PARMS}.$$" >> $F_SUMM_INPUT
    else
      echo "3d displacement: Needed dimensions are unavailable!"
    endif
endif
set estsmpar3d = `most_recent_parms.py parsm3d`
if ( dummy$estsmpar3d != dummy ) then
    if ( $?xdim3d && $?ydim3d && $?zdim3d ) then
      displace3d -e par/${F_DISP3D_SMTH_PARMS}.$$ \
	-xdimension $xdim3d -ydimension $ydim3d -zdimension $zdim3d \
	-xvoxel ${F_XVOXEL} -yvoxel ${F_YVOXEL} -zvoxel ${F_ZVOXEL} \
        -filtered $weightstr $estsmpar3d
      set step = `depath.csh $0`3d
      echo "$step par/${F_DISP3D_SMTH_PARMS}.$$" >> $F_SUMM_INPUT
    else
      echo "3d displacement: Needed dimensions are unavailable!"
    endif
endif

if ( $wt_generated != 0 ) then
  mri_destroy_dataset $wt_generated
endif
