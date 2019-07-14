#!/bin/csh -efx
# estireg3d.csh
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
# *  Original programming by Bill Eddy                       *
# ************************************************************/
#
echo '#'`date` $0
echo '#$Id: estireg3d.csh,v 1.6 2003/10/01 20:49:13 bakalj Exp $'
if(! -d par) mkdir par 

# Default behavior is to align to the mean
if ( ! ${?F_EST3D_ALIGN} ) then
  setenv F_EST3D_ALIGN "mean"
endif

set dimstr = `mri_printfield -field images.dimensions $1`
if ( $dimstr == xyzt ) then
  set perm_order = "txyz"
else
  set perm_order = "vtxyz"
endif

if ( $F_EST3D_ALIGN == "mean" ) then
#
# Calculate mean and stdv images.  This involves creating a temporary
# fake newsplit file
#
  if (! -d stat) mkdir stat
  mri_rpn_math -out data/junk -c missing \
      '$t,$z,1.0,1.0,if_print_3,0.0' $1 > fake_newsplit
  stats -condition fake_newsplit -input $1 -meanprefix stat/unmoved \
      -stdvprefix stat/unmoved -tsprefix stat/unmoved -maxtpairs 0
  mri_destroy_dataset data/junk
  rm fake_newsplit
  set target = stat/unmovedMean_1
  set stdvtarget = stat/unmovedStdv_1
else if ( $F_EST3D_ALIGN == "median" ) then
#
# Calculate median and iqr images.
#
  mri_permute -order $perm_order $1 ${1}_p 
  mri_subsample -median -d t -l 1 ${1}_p ${1}_median
  mri_subsample -iqr -d t -l 1 ${1}_p ${1}_iqr
  mri_destroy_dataset ${1}_p
  mri_remap -order vxyzt ${1}_median 
  mri_remap -order vxyzt ${1}_iqr
  set target = ${1}_median
  set stdvtarget = ${1}_iqr
else
  set target = $F_EST3D_ALIGN
  set stdvtarget = $F_EST3D_STDVALN
endif
#
# Estimate 3D alignment parameters
#
parallel.run.csh estireg3d -xvoxel ${F_XVOXEL} -yvoxel ${F_YVOXEL} \
    -zvoxel ${F_ZVOXEL} -input $1 -align ${target} -stdv ${stdvtarget} \
    -parameters par/${F_EST3D_PARMS}.$$ -alg ${F_EST3D_ALG}
set step = `depath.csh $0`
echo "$step par/${F_EST3D_PARMS}.$$" >> $F_SUMM_INPUT
#
# Clean up
#
if ( $F_EST3D_ALIGN == "mean" ) then
  mri_destroy_dataset stat/unmovedMean_1
  mri_destroy_dataset stat/unmovedStdv_1
endif
