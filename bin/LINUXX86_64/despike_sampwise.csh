#! /bin/csh -exf
# despike_sampwise.csh
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
echo '#'`date` $0
echo '#$Id: despike_sampwise.csh,v 1.3 2003/09/09 16:25:45 bakalj Exp $'

#
# spk_t0_cutoff is the cutoff threshold for image 0, in terms of
#            (orig-smoothed)/(IQR)
# spk_t0_band is the smoother bandwidth used for image t0.
#
set spk_band = $F_DESPIKE_T0_BND
set spk_cutoff = $F_DESPIKE_T0_CUT

#
# Make some of the parameters easier to remember
#
set q1 = $2
set q3 = $3
set result = $4
set spike_count = $5

#
# We guess at where the spikes are by looking for rapid changes in 
# magnitude, taking only the positive changes in this case.  (Otherwise 
# we end up picking up samples adjacent to the spike because of the 
# smoothing operation).
#
mri_complex_to_scalar -mag ${1} ${1}_m
mri_complex_to_scalar -phu ${1} ${1}_p
mri_smooth -d p -bdw $spk_band ${1}_m ${1}_smooth_m
mri_rpn_math -out ${1}_trigger '$1,$2,-,$3,$4,-,/' ${1}_m ${1}_smooth_m $q3 $q1
mri_destroy_dataset ${1}_smooth_m
mri_rpn_math -out ${1}_is_a_spike '$1,'$spk_cutoff',<' ${1}_trigger
mri_destroy_dataset ${1}_trigger

#
# Through clever use of a triangular filter we happen to have handy,
# we'll produce a set of replacement values for the spikes which are
# the average of the two adjacent points.
#
mri_smooth -d p -bdw 2 -kernel triangular ${1}_m ${1}_triangle_m
mri_smooth -d p -bdw 2 -kernel triangular ${1}_p ${1}_triangle_p
mri_rpn_math -out ${1}_interp_m '$1,$2,0.5,*,-,2,*' ${1}_triangle_m ${1}_m
mri_destroy_dataset ${1}_triangle_m
mri_destroy_dataset ${1}_m
mri_rpn_math -out ${1}_interp_p '$1,$2,0.5,*,-,2,*' ${1}_triangle_p ${1}_p
mri_destroy_dataset ${1}_triangle_p
mri_destroy_dataset ${1}_p

#
# Substitute these replacement values where appropriate 
#
mri_rpn_math -out ${1}_r_interp '$1,$2,cos,*' \
    ${1}_interp_m ${1}_interp_p
mri_rpn_math -out ${1}_i_interp '$1,$2,sin,*' \
    ${1}_interp_m ${1}_interp_p
mri_destroy_dataset ${1}_interp_m
mri_destroy_dataset ${1}_interp_p
mri_paste -d v -out ${1}_interp ${1}_r_interp ${1}_i_interp
mri_destroy_dataset ${1}_r_interp
mri_destroy_dataset ${1}_i_interp
mri_interp -con -d v -len 2 ${1}_is_a_spike ${1}_is_a_spike_complex
mri_rpn_math -out $result '$1,$2,$3,if_keep' \
    ${1} ${1}_interp ${1}_is_a_spike_complex
mri_destroy_dataset ${1}_interp
mri_destroy_dataset ${1}_is_a_spike_complex

# Produce spike count information
mri_subsample -sum -d p -l 1 ${1}_is_a_spike $spike_count

# Clean up
mri_destroy_dataset ${1}_is_a_spike
