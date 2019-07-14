#!/bin/csh -efx
# ts.deghost.csh
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
echo '$Id: ts.deghost.csh,v 1.12 2003/05/14 18:17:04 welling Exp $'
#
echo '#'`date`$0
if(! -d par) mkdir par

# We need to treat the second shot as a set of additional slices to
# avoid messing up the smoothing of estimates.
set zdim = `mri_printfield -i $1 -field images.extent.z`
set tdim = `mri_printfield -i $1 -field images.extent.t`
set dimorder = `mri_printfield -i $1 -field images.dimensions`
if ( $dimorder != 'vxyzt' ) then
  echo "ts.deghost.csh: input is not in vxyzt order!"
  exit -1
endif
set vdim = `mri_printfield -i $1 -field images.extent.v`
if ( $vdim != 2) then
  echo "ts.deghost.csh: input is not complex!"
  exit -1
endif
@ tdim_half = ( $tdim / 2 )
@ zdim_dbl = ( $zdim * 2 )
mri_remap -file $1 -order vxyzt -extents 2:::${zdim_dbl}:${tdim_half}

epi.deghost.csh $1 $2

# Remap the input and output file back to the original dimension order
mri_remap -file $1 -order vxyzt -extents 2:::${zdim}:${tdim}
mri_remap -file $2 -order vxyzt -extents 2:::${zdim}:${tdim}

