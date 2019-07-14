#!/bin/csh -exf
# ispace_filter.csh
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
echo '$Id: ispace_filter.csh,v 1.7 2007/06/07 21:58:04 welling Exp $'
#
echo '# '`date`' '$0

#
# Check input dataset shape
#
set dimstr = `mri_printfield -field 'images.dimensions' $1`
if ( $dimstr !~ vxy* && $dimstr !~ xy* ) then
  echo "Input dataset dimensions don't begin with vxy!"
  exit -1
endif

if ( $dimstr =~ v* ) then
  set vdim = `mri_printfield -field 'images.extent.v' $1`
  if ( $vdim != 1 ) then
    echo "Input dataset must be scalar!"
    exit -1
  endif
endif

# We need to be sure to use integer math for this; csh does so
set xdim = `mri_printfield -field images.extent.x $1`
set ydim = `mri_printfield -field images.extent.y $1`
set zdim = `mri_printfield -field images.extent.z $1`
@ xhalf = ( $xdim / 2 )
@ yhalf = ( $ydim / 2 )
@ zhalf = ( $zdim / 2 )

#
# Calculate filter dimensions.  We use dimensions specified in image
# space.  Once again we must follow the ridiculous practice of using
# mri_rpn_math to do floating point calculations from csh.
#
mri_counter -d t -l 1 ifilter_tmp1
set iradius = ${F_IFILTER_RADIUS}
set radius = \
  `mri_rpn_math -out ifilter_tmp2 '0,'${xhalf}','${iradius}',/,1,if_print_1' ifilter_tmp1`
mri_destroy_dataset ifilter_tmp1
mri_destroy_dataset ifilter_tmp2

if ( ! -d data ) mkdir data
set filter_name = data/filter_${F_IFILTER_TYPE}.$$
set scaled_filter_name = data/filter_${F_IFILTER_TYPE}_scaled.$$

#
# Create the filter.  Filter calculations depend on the grid
# dimension in the filtering direction, and so depend on the
# filter type.
#
switch ( ${F_IFILTER_TYPE} )
  # These filters are in-plane, or full 3D.  Data extents assumed
  # equal over filter directions.
  case 'none':
    breaksw;
  case 'cylinder':
  case 'boxcar':
  case 'fermi':
  case 'gaussian':
  case 'gaussianFWHM':
  case 'gaussian3D':
  case 'gaussian3DFWHM':
  case 'hanning':
  case 'hamming':
  case 'welch':
    set radius = \
      `python -c 'from math import *; print '$xdim'/(2*'$iradius')'`
  breaksw;
  # These filters are in the slice direction.
  case zGaussian:
  case zGaussianFWHM:
    set radius = \
      `python -c 'from math import *; print '$zdim'/(2*'$iradius')'`
  breaksw;
  default :
    echo "Error- Unknown filter type" ${type} "requested"
    exit -1
  breaksw;

endsw

#
# Kick the data and filter to kspace, apply the filter, and 
# pull the data back to image space.  Don't forget to re-scale the
# filter so that the value of a constant is preserved under filtering.
#
mri_fft -d xyz -cpx $1 ${1}_kspace
make_filter.csh ${F_IFILTER_TYPE} 1 ${radius} ${1}_kspace ${filter_name}
mri_complex_to_scalar -real ${filter_name} ${filter_name}_r
mri_subset -d x -l 1 -s $xhalf ${filter_name}_r ${filter_name}_rx
mri_subset -d y -l 1 -s $yhalf ${filter_name}_rx ${filter_name}_rxy
mri_subset -d z -l 1 -s $zhalf ${filter_name}_rxy ${filter_name}_rxyz
mri_rpn_math -out ${scaled_filter_name} '$1,$2,/' \
    ${filter_name} ${filter_name}_rxyz
mri_rpn_math -out ${2}_kspace '$1,$2,*' ${1}_kspace ${scaled_filter_name}
mri_destroy_dataset ${1}_kspace
mri_fft -inv -d xyz -mod ${2}_kspace $2
foreach dsname ( ${2}_kspace ${filter_name} ${filter_name}_r \
                 ${filter_name}_rx ${filter_name}_rxy ${filter_name}_rxyz )
  mri_destroy_dataset ${dsname}
end

