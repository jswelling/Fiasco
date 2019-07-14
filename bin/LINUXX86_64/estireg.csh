#!/bin/csh -efx
# estireg.csh
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
echo '#$Id: estireg.csh,v 1.14 2003/09/22 04:02:33 bakalj Exp $'
if(! -d par) mkdir par 
set tmpname = ${F_ESTIREG_FIXED}
set fullfilename = ${tmpname:r}.mri

set dimstr = `mri_printfield -field images.dimensions $1`
if ( $dimstr == xyzt ) then
  set perm_order = "txyz"
else
  set perm_order = "vtxyz"
endif

if ( $F_ESTIREG_FIXED == "mean" ) then

# Align to mean of input images
  mri_permute -order $perm_order $1 ${1}_p
  mri_subsample -mean -d t -l 1 ${1}_p ${1}_mean
  mri_destroy_dataset ${1}_p
  mri_remap -order vxyzt ${1}_mean 
  parallel.run.csh estireg -input $1.mri -parameters par/$F_ESTIREG_PARMS.$$ \
                -fixed 0 -align ${1}_mean

else if ( $F_ESTIREG_FIXED == "median" ) then

# Align to mean of input images
  mri_permute -order $perm_order $1 ${1}_p
  mri_subsample -median -d t -l 1 ${1}_p ${1}_median
  mri_destroy_dataset ${1}_p
  mri_remap -order vxyzt ${1}_median 
  parallel.run.csh estireg -input $1.mri -parameters par/$F_ESTIREG_PARMS.$$ \
                -fixed 0 -align ${1}_median

else if ( -e ${fullfilename} ) then

  mri_remap -order vxyzt ${F_ESTIREG_FIXED} 
  parallel.run.csh estireg -input $1.mri -parameters par/$F_ESTIREG_PARMS.$$ \
                -fixed 0 -align ${F_ESTIREG_FIXED}

else

# Align to the selected image number
  parallel.run.csh estireg -input $1.mri -parameters par/$F_ESTIREG_PARMS.$$ \
                -fixed $F_ESTIREG_FIXED
endif
set step = `depath.csh $0`
echo "$step par/$F_ESTIREG_PARMS.$$" >> $F_SUMM_INPUT
