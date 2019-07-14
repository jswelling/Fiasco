#!/bin/csh -ef
# ts.reconcile.csh
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
echo '# $Id: ts.reconcile.csh,v 1.3 2002/01/08 00:49:18 welling Exp $'

#Get dimensions
set vextent = \
`mri_printfield -input data/$F_READER_OUTPUT -field images.extent.v -nofail`
set textent = \
`mri_printfield -input data/$F_READER_OUTPUT -field images.extent.t -nofail`
set zextent = \
`mri_printfield -input data/$F_READER_OUTPUT -field images.extent.z -nofail`
set xextent = \
`mri_printfield -input data/$F_READER_OUTPUT -field images.extent.x -nofail`
set yextent = \
`mri_printfield -input data/$F_READER_OUTPUT -field images.extent.y -nofail`
set sextent = \
`mri_printfield -input data/$F_READER_OUTPUT -field images.extent.s -nofail`

if ( dummy$sextent == dummy ) then
  set interleaved_slices = 0
else
  set interleaved_slices = 1
endif

# Set slice and image numbers
setenv F_NSLICE $zextent
if ($interleaved_slices) then
  setenv F_NIMAGE $textent
else 
  @ t_half = $textent / 2
  setenv F_NIMAGE $t_half
endif

# Are clip dims compatible with image dims?
if ($interleaved_slices) then
  @ xctr = $xextent / 2
  @ yctr = $yextent
else
  @ xctr = $xextent
  @ yctr = $yextent / 2
endif
if ($F_CLIP1_XCENTER != $xctr) then
  echo "%%%%% Error! First clip not centered in X %%%%%"
  exit -1
endif
if ($F_CLIP1_YCENTER != $yctr) then
  echo "%%%%% Error! First clip not centered in Y %%%%%"
  exit -1
endif
if ($F_CLIP1_WIDTH > $xextent) then
  echo "%%%%% Error! First clip region too wide %%%%%"
  exit -1
endif
if ($F_CLIP1_HEIGHT > $yextent) then
  echo "%%%%% Error! First clip region too high %%%%%"
  exit -1
endif
if ($interleaved_slices) then
  echo "%%%%% Error! Interleaved twoshow EPI not yet implemented %%%%%"
  exit -1
endif

# Set some titling information
setenv F_PS_TITLE1 "$F_DESCRIPTION"
setenv F_PS_TITLE2 "$F_NSLICE Slices"
