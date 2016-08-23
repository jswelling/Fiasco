#!/bin/csh -efx
# threshold_mse.csh
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
echo '#$Id: threshold_mse.csh,v 1.5 2003/09/10 19:59:22 bakalj Exp $'
#
# Input parameters: 
#  1) input MRI file
#  2) output MRI file
#  3) registration parameter file
#  4) MSE threshold value as a float
#
set ifile = $1
set ofile = $2
set parfile = $3
set thresh = $4
set nt = `mri_printfield -field 'images.extent.t' $ifile`
set nz = `mri_printfield -field 'images.extent.z' $ifile`

# Create output file; transfer images chunk
mri_copy_chunk -rep -c images $ifile $ofile

# Copy out missing info for modification
mri_copy_chunk -rep -c missing $ifile tmp_tmse

# Make a file of mse values
mri_from_ascii -order vtz -c missing -l 4:${nt}:${nz} tz \
	tmp2_tmse < $parfile
mri_permute -c missing -o vzt -memlimit 32000000 tmp2_tmse tmp3_tmse 
mri_subset -d v -l 1 -s 3 tmp3_tmse tmp4_tmse

# Mask missing info with thresholded MSE data. Multiply by 1.1
# is to assure that conversion to bytes which follows doesn't
# have rounding errors.
mri_rpn_math -out tmp5_tmse -c missing \
	'1.0,$1,'${thresh}',$2,<,if_keep,1.1,*' tmp_tmse tmp4_tmse 
mri_type_convert -char -c missing tmp5_tmse tmp6_tmse

# Finish assembling output
mri_copy_chunk -rep -c missing tmp6_tmse $ofile

# Clean up
rm tmp_tmse.mri
rm tmp2_tmse.mri tmp2_tmse.dat
rm tmp3_tmse.mri tmp3_tmse.dat
rm tmp4_tmse.mri tmp4_tmse.dat
rm tmp5_tmse.mri
rm tmp6_tmse.mri

