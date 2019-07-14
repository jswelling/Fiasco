#!/bin/csh -efx
# epi.vhr.csh
# /************************************************************
# *                                                          *
# *  Permission is hereby granted to any individual or       *
# *  institution for use, copying, or redistribution of      *
# *  this code and associated documentation, provided        *
# *  that such code and documentation are not sold for       *
# *  profit and the following copyright notice is retained   *
# *  in the code and documentation:                          *
# *     Copyright (c) 1999 Department of Statistics,         *
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
# *  Original programming by Bill Eddy,                      *
# *  Additional programming by Scott Ziolko                  *
# ************************************************************/
##################
# 
# Usage: substitute this script for epi.recon2.csh and epi.clip2.csh,
# for example:
#
# epi.vhr.csh data/$F_PHYS_OUTPUT data/$F_CLIP2_OUTPUT data/$F_CLIP1_OUTPUT
#
# The default enhanced dimensions are calculated below;  they may be
# changed by adding lines like the following (but with other
# numbers) to your *.local.csh file.  The default is to double
# resolution in both directions.
#
# setenv F_PAD_RNG_X 256
# setenv F_PAD_RNG_Y 128
##################
#
echo '#'`date`$0
echo '#$Id: epi.vhr.csh,v 1.7 2003/09/10 20:17:23 bakalj Exp $'
#
#

# Define the names of some temporary files
setenv F_PADX_OUTPUT padx
setenv F_PADY_OUTPUT pady
setenv F_VHR_OUTPUT vhrecon

# Get the dimensions of the image in order to correctly
# offset the images while padding
set XDIM = `mri_printfield -field images.extent.x $1`
set YDIM = `mri_printfield -field images.extent.y $1`

# Make up default padding dims if not already set
if (! ${?F_PAD_RNG_X} ) then
  @ range_x = $XDIM * 2
  setenv F_PAD_RNG_X $range_x
endif
if (! ${?F_PAD_RNG_Y} ) then
  @ range_y = $YDIM * 2
  setenv F_PAD_RNG_Y $range_y
endif

# Determine offset in order to center data in padding
@ OFFSET_XDIM = $F_PAD_RNG_X / 2 - $XDIM / 2

# Pad the x-dimension to range supplied in epi.local.csh
# -d dimension -l length -s shift infile outfile
# Default padding is done with 0 in all locations
mri_pad -d x -l $F_PAD_RNG_X -s $OFFSET_XDIM $1 data/$F_PADX_OUTPUT
echo "x-dimension padding complete"

# Determine offset in order to center data in padding
@ OFFSET_YDIM = $F_PAD_RNG_Y / 2 - $YDIM / 2

# Pad the y-dimension to range supplied in epi.local.csh
# -d dimension -l length -s shift infile outfile
# Default padding is done with 0 in all locations
mri_pad -d y -l $F_PAD_RNG_Y -s $OFFSET_YDIM data/$F_PADX_OUTPUT \
	data/$F_PADY_OUTPUT
echo "y-dimension padding complete"

rm data/${F_PADX_OUTPUT}.mri data/${F_PADX_OUTPUT}.dat 
 
# Set clipping parameters
# 4 clipping parameters are needed,
# F_CLIP2_WIDTH which is the width of the padded data
# F_CLIP2_HEIGHT which is the height of the padded data
# F_CLIP2_XCENTER which specifies the location of the center
# of the data in the x-dimension
# F_CLIP2_YCENTER which specifies the location of the center
# of the data in the y-dimension

@ XCENTER = ${F_PAD_RNG_X} / 2
setenv F_CLIP2_XCENTER ${XCENTER}

@ YCENTER = ${F_PAD_RNG_Y} / 2
setenv F_CLIP2_YCENTER ${YCENTER}

# This is used to account for 64X64 as well as 128X64 raw data
if (${F_PAD_RNG_X} > ${F_PAD_RNG_Y}) then
  setenv F_CLIP2_WIDTH ${F_PAD_RNG_Y}
  setenv F_CLIP2_HEIGHT ${F_PAD_RNG_Y}
else
  setenv F_CLIP2_WIDTH ${F_PAD_RNG_X}
  setenv F_CLIP2_HEIGHT ${F_PAD_RNG_X}
endif

# reconstruct corrected, registered images
echo '#'`date`$0
echo '# Reconstructing Hi-Res Images $'
parallel.run.csh recon -input data/${F_PADY_OUTPUT}.mri \
	-headerout data/${F_VHR_OUTPUT}.mri \
	-dataout .dat -direction $F_RECON2_DIRECTION \
	-recon $F_RECON2_MODULUS
rm data/${F_PADY_OUTPUT}.mri data/${F_PADY_OUTPUT}.dat 

# Need to copy "missing" chunk from dataset used for motion estimation
# to new reconstructed dataset
mri_copy_chunk -chunk missing -replace $3 data/${F_VHR_OUTPUT}

#
# cut them down to size to speed up computations
echo 'Clipping Hi-Res Images $'
epi.clip2.csh data/${F_VHR_OUTPUT} $2

# clean up
mri_destroy_dataset data/${F_VHR_OUTPUT}
mri_destroy_dataset $1
mri_destroy_dataset $3



