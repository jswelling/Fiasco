#!/bin/csh -exf
# epi.filter.csh
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
echo '$Id: epi.filter.csh,v 1.7 2003/09/10 16:44:02 bakalj Exp $'
#
echo '#'`date` $0

#
# Get filter width and radius
#
set width = `mri_printfield -field images.fermi_width -nofail $1`
if ( dummy${width} == dummy ) then
  set width = ${F_FILTER_WIDTH}
endif
set radius = `mri_printfield -field images.fermi_radius -nofail $1`
if ( dummy${radius} == dummy ) then
  if ( ${?F_FILTER_RADIUS} ) then
    set radius = ${F_FILTER_RADIUS}
  else
    set xdim = `mri_printfield -field images.extent.x $1`
    @ radius = ( $xdim / 2 )
  endif
endif

if ( ! -d data ) mkdir data
set filter_name = data/filter_${F_FILTER_TYPE}.$$

# Create and apply the filter
make_filter.csh ${F_FILTER_TYPE} ${width} $radius} $1 ${filter_name}
mri_rpn_math -out $2 '$1,$2,*' $1 ${filter_name}



