#!/bin/csh -ef
# make_filter.csh
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
# ************************************************************/
#
# usage: make_filter.csh type width radius proto outname
#
echo '$Id: make_filter.csh,v 1.7 2007/06/07 21:58:04 welling Exp $'
#
echo '# '`date`' '$0

if ( $#argv != 5 ) then
  echo "usage: $0 type width radius protoDS outDS"
  exit -1
endif

set type = $1
set width = $2
set radius = $3
set inDS = $4
set outDS = $5

#
# Get filter width and radius
#
set xdim = `mri_printfield -field images.extent.x ${inDS}`
set ydim = `mri_printfield -field images.extent.y ${inDS}`
set zdim = `mri_printfield -field images.extent.z ${inDS}`

# We need to be sure to use integer math for this; csh does so
@ xctr = ( $xdim / 2 )
@ yctr = ( $ydim / 2 )
@ zctr = ( $zdim / 2 )

mri_subset -d t -l 1 -s 0 ${inDS} tmp_filter_shape

echo "### Generating " ${type} "filter" ${outDS}
switch ( ${type} )
  case 'none' :
      mri_rpn_math -out ${outDS} '1' tmp_filter_shape
    breaksw;

  case 'cylinder' :
    mri_rpn_math -out ${outDS} \
      '1,0,$x,'$xctr',-,dup,*,$y,'$yctr',-,dup,*,+,sqrt,'$radius',<,if_keep' \
      tmp_filter_shape
    breaksw;

  case 'boxcar' :
    mri_rpn_math -out ${outDS} \
      '1,0,$x,'$xctr',-,'$radius',<,if_keep,1,0,$y,'$yctr',-,'$radius',<,if_keep,*' \
      tmp_filter_shape
    breaksw;

  case 'fermi' :
      mri_rpn_math -out ${outDS} \
        '1,1,$x,'$xctr',-,dup,*,$y,'$yctr',-,dup,*,+,sqrt,'$radius',-,'$width',/,exp,+,/' \
        tmp_filter_shape
    breaksw;

  case 'gaussian' :
      mri_rpn_math -out ${outDS} \
        '$x,'$xctr',-,dup,*,$y,'$yctr',-,dup,*,+,'$radius',dup,*,/,-1,*,exp,'$radius',dup,*,pi,*,2,sqrt,*,/' \
        tmp_filter_shape
    breaksw;

  case 'gaussianFWHM' :
      # Assumes xdim==ydim
      # Treats input 'radius' as if xdim/(2*radius) were intended to
      # be the full width half maximum in image space of a convolution
      # kernel of which this is the k-space version.  This is appropriate
      # for using ispace_filter with IFILTER_RADIUS being the intended
      # FWHM in pixels
      if ( $xdim != $ydim ) then
        echo "xdim does not equal ydim for gaussianFWHM in $0 "
        exit -1
      endif
      set b_kspace = \
        `python -c 'from math import *; print 4*sqrt(log(2))*'$radius'/pi'`
      mri_rpn_math -out ${outDS} \
        '$x,'$xctr',-,dup,*,$y,'$yctr',-,dup,*,+,'$b_kspace',dup,*,/,-1,*,exp,'$radius',dup,*,pi,*,2,sqrt,*,/' \
        tmp_filter_shape
    breaksw;

  case 'zGaussian' :
      mri_rpn_math -out ${outDS} \
        '$z,'$zctr',-,dup,*,'$radius',dup,*,/,-1,*,exp' \
        tmp_filter_shape
    breaksw;

  case 'zGaussianFWHM' :
      # Treats input 'radius' as if zdim/(2*radius) were intended to
      # be the full width half maximum in image space of a convolution
      # kernel of which this is the k-space version.  This is appropriate
      # for using ispace_filter with IFILTER_RADIUS being the intended
      # FWHM in pixels
      set b_kspace = \
        `python -c 'from math import *; print 4*sqrt(log(2))*'$radius'/pi'`
      echo "b_kspace is" $b_kspace
      mri_rpn_math -out ${outDS} \
        '$z,'$zctr',-,dup,*,'$b_kspace',dup,*,/,-1,*,exp' \
        tmp_filter_shape
    breaksw;

  case 'gaussian3D' :
      mri_rpn_math -out ${outDS} \
        '$x,'$xctr',-,dup,*,$y,'$yctr',-,dup,*,+,$z,'$zctr',-,dup,*,+,'$radius',dup,*,/,-1,*,exp' \
        tmp_filter_shape
    breaksw;

  case 'gaussian3DFWHM' :
      # Assumes xdim==ydim==zdim
      # Treats input 'radius' as if xdim/(2*radius) were intended to
      # be the full width half maximum in image space of a convolution
      # kernel of which this is the k-space version.  This is appropriate
      # for using ispace_filter with IFILTER_RADIUS being the intended
      # FWHM in pixels
      if ( $xdim != $ydim ) then
        echo "xdim does not equal ydim for gaussian3DFWHM in $0 "
        exit -1
      endif
      if ( $xdim != $zdim ) then
        echo "xdim does not equal zdim for gaussian3DFWHM in $0 "
        exit -1
      endif
      set b_kspace = \
        `python -c 'from math import *; print 4*sqrt(log(2))*'$radius'/pi'`
      mri_rpn_math -out ${outDS} \
        '$x,'$xctr',-,dup,*,$y,'$yctr',-,dup,*,+,$z,'$zctr',-,dup,*,+,'$b_kspace',dup,*,/,-1,*,exp' \
        tmp_filter_shape
    breaksw;

  case 'hanning' :
      mri_rpn_math -out ${outDS} \
	'0.5,2,pi,*,$x,*,$xdim,1,-,/,cos,0.5,*,-,0.5,2,pi,*,$y,*,$ydim,1,-,/,cos,0.5,*,-,*' \
        tmp_filter_shape
    breaksw;

  case 'hamming' :
      mri_rpn_math -out ${outDS} \
	'0.54,2,pi,*,$x,*,$xdim,1,-,/,cos,0.46,*,-,0.54,2,pi,*,$y,*,$ydim,1,-,/,cos,0.46,*,-,*' \
        tmp_filter_shape
    breaksw;

  case 'welch' :
      mri_rpn_math -out ${outDS} \
        '1,$x,0.5,$xdim,1,-,*,-,0.5,$xdim,1,+,*,/,dup,*,-,1,$y,0.5,$ydim,1,-,*,-,0.5,$ydim,1,+,*,/,dup,*,-,*' \
        tmp_filter_shape
    breaksw;

  default :
    echo "Error- Unknown filter type" ${type} "requested"
    exit -1
  breaksw;

endsw

mri_destroy_dataset tmp_filter_shape


