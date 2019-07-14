#! /bin/csh -ef
# weighted_sum_squares.csh.csh
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
# $Id: weighted_sum_squares.csh,v 1.2 2003/09/10 20:01:27 bakalj Exp $

#
# Usage:
#  weighted_sum_squares.csh dataset_to_sum weight_dataset
#
# Results for all times are emitted as ascii, one per line

# Can we handle this dataset?
set dimstr = `mri_printfield -field images.dimensions $1`
if ( $dimstr != "vxyzt" ) then
  echo "# Error: ${0}: input dataset $1 is not vxyzt!"
  exit -1
endif
set wt_dimstr = `mri_printfield -field images.dimensions $2`
if ( $wt_dimstr != "xyz" ) then
  echo "# Error: ${0}: weight dataset $2 is not xyz!"
  exit -1
endif

# Find out some dimensions
set vdim = `mri_printfield -field images.extent.v $1`
set xdim = `mri_printfield -field images.extent.x $1`
set ydim = `mri_printfield -field images.extent.y $1`
set zdim = `mri_printfield -field images.extent.z $1`
set tdim = `mri_printfield -field images.extent.t $1`
@ qdim = ( $xdim * $ydim * $zdim )

if ( $vdim > 2 ) then
  echo "# Error: ${0}: input is not scalar or complex!"
  exit -1
endif

set wt_xdim = `mri_printfield -field images.extent.x $2`
set wt_ydim = `mri_printfield -field images.extent.y $2`
set wt_zdim = `mri_printfield -field images.extent.z $2`
if ( $wt_xdim != $xdim ) then
  echo "# Error: ${0}: incommensurate X dimensions!"
  exit -1
endif
if ( $wt_ydim != $ydim ) then
  echo "# Error: ${0}: incommensurate Y dimensions!"
  exit -1
endif
if ( $wt_zdim != $zdim ) then
  echo "# Error: ${0}: incommensurate Z dimensions!"
  exit -1
endif

#
# Make some scratch space
#
if ( ${?F_TEMP} ) then
  set tmpdir = ${F_TEMP}/weighted_sum_squares_$$
else
  set tmpdir = ./weighted_sum_squares_$$
endif
if (! -e $tmpdir) mkdir $tmpdir
set homedir = $PWD

if ( $vdim == 1 ) then
  mri_copy_dataset $1 $tmpdir/mags
else
  mri_complex_to_scalar -mag $1 $tmpdir/mags
endif
mri_copy_dataset $2 $tmpdir/weights

cd $tmpdir

# Calculate the weighted sum of squares
mri_rpn_math -out weighted_squares '$1,dup,*,$2,*' mags weights
mri_remap -order qt -length ${qdim}:${tdim} weighted_squares 
mri_remap -order q -length ${qdim} weights
mri_subsample -d q -l 1 -sum weighted_squares total
mri_subsample -d q -l 1 -sum weights wt_total

# Emit results
mri_rpn_math -out junk '0,$t,$1,$2,/,1,if_print_2' total wt_total

# Clean up
cd $homedir
rm -r $tmpdir

