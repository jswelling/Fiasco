#!/bin/csh -efx
# epi.clip2.csh
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
echo '#$Id: epi.clip2.csh,v 1.9 2003/09/09 16:25:45 bakalj Exp $'

set tmpfile = $F_TEMP/tmp_clip2_$$
set xdim = `mri_printfield -field images.extent.x $1`
set ydim = `mri_printfield -field images.extent.y $1`

if ( $F_CLIP2_WIDTH == "whole" ) then
  @ width = $xdim
else if ( $F_CLIP2_WIDTH == "half" ) then
  @ width = $xdim / 2
else
  @ width = $F_CLIP2_WIDTH
endif

if ( $F_CLIP2_HEIGHT == "whole" ) then
  @ height = $ydim
else if ( $F_CLIP2_HEIGHT == "half" ) then
  @ height = $ydim / 2
else
  @ height = $F_CLIP2_HEIGHT
endif

if ( $F_CLIP2_XCENTER == "center" ) then
  @ xoff = ( $xdim / 2 ) - ( $width / 2 )
else
  @ xoff = $F_CLIP2_XCENTER - ( $width / 2 )
endif

if ( $F_CLIP2_YCENTER == "center" ) then
  @ yoff = ( $ydim / 2 ) - ( $height / 2 )
else
  @ yoff = $F_CLIP2_YCENTER - ( $height / 2 )
endif

mri_subset -d x -l $width -s $xoff $1 $tmpfile
mri_subset -d y -l $height -s $yoff $tmpfile $2

mri_destroy_dataset $tmpfile

