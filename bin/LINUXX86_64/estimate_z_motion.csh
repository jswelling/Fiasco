#!/bin/csh -efx
# estimate_z_motion.csh
#/************************************************************
# *                                                          *
# *  Permission is hereby granted to any individual or       *
# *  institution for use, copying, or redistribution of      *
# *  this code and associated documentation, provided        *
# *  that such code and documentation are not sold for       *
# *  profit and the following copyright notice is retained   *
# *  in the code and documentation:                          *
# *     Copyright (c) 1998 Department of Statistics,         *
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
# *  Original programming by Joel Welling                    *
# ************************************************************/
#
# This script makes a crude estimate of motion in the Z direction
# for a dataset.  The input data is assumed to be in image space,
# and to have dimensions vxyzt.  The z dimension gets mapped onto
# the x direction, and the x direction onto the y direction.  The
# logical place to use this script would be after the clip2 step
# of Fiasco.
echo '#'`date` $0
echo '# $Id: estimate_z_motion.csh,v 1.10 2003/09/10 20:17:23 bakalj Exp $'

# Pull in some parameters for smregpar.  This will smash in/split, so
# save and restore it.
cp in/split in/e_z_motion_split_sav
source ${FIASCO}/epi.defaults.csh
rm in/split
mv in/e_z_motion_split_sav in/split

# Generate some variables
set fname = ${1:r}
set ntimes = `mri_printfield -fld images.extent.t ${fname}`
set nslices = `mri_printfield -fld images.extent.z ${fname}`
set nlines = `mri_printfield -fld images.extent.y ${fname}`
set ncols = `mri_printfield -fld images.extent.x ${fname}`

# Assumed input is image-space data in vxyzt form
# Permute and relabel.  Because of the algorithm employed by permute,
# it is actually faster to do two permutes than one to accomplish
# this task.
mri_permute -memlimit 32000000 -chunk images -order vyztx ${fname} ${fname}_p
mri_permute -memlimit 32000000 -chunk images -order vyzxt ${fname}_p ${fname}_p2 
rm ${fname}_p.mri ${fname}_p.dat
mri_remap  -c images -order vxyzt \
	-length 1:${nlines}:${nslices}:${ncols}:${ntimes} ${fname}_p2

# Now we average down in the fake z direction, and up in the fake x dir
mri_subsample -d z -l 7 -closest ${fname}_p2 ${fname}_s
rm ${fname}_p2.mri ${fname}_p2.dat
mri_subset -d z -s 2 -l 3 ${fname}_s ${fname}_S
rm ${fname}_s.mri ${fname}_s.dat
mri_interp -d y -len $nlines -lin ${fname}_S ${fname}_i
rm ${fname}_S.mri ${fname}_S.dat

# Estimate motion on these fake slices
@ fixed = $ntimes / 2
parallel.run.csh estireg -input ${fname}_i -parameters est_z_motion.par \
  -fixed $fixed
smregpar -headerinput ${fname}_i -parameterin est_z_motion.par \
	-parameterout est_z_motion_smooth.par \
	-cutoff $F_PARSM_CUTOFF -fixed $fixed \
	-bandwidth $F_PARSM_BANDWIDTH -kernel $F_PARSM_KERNEL \
        -threshold $F_PARSM_TRANST
echo "estireg ./est_z_motion.par" > z_motion_summary.input
echo "parsm ./est_z_motion_smooth.par" >> z_motion_summary.input



