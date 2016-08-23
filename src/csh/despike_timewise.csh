#! /bin/csh -exf
# despike_timewise.csh
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
echo '#$Id: despike_timewise.csh,v 1.7 2004/07/01 23:59:54 welling Exp $'

#
# spk_cutoff is the cutoff threshold (in terms of (mag-median)/(IQR))
#            for the body of the data (all but image 0)
set spk_cutoff = $F_DESPIKE_CUTOFF

#
# Prepare median and quartiles for magnitude.
#
mri_complex_to_scalar -mag ${1} ${1}_m
mri_permute -order vtpscz ${1}_m ${1}_m_prm 
mri_subsample -d t -l 1 -median ${1}_m_prm ${1}_m_median
mri_subsample -d t -l 1 -1qr ${1}_m_prm ${1}_m_q1
mri_subsample -d t -l 1 -3qr ${1}_m_prm ${1}_m_q3
mri_remap -order vpsczt ${1}_m_median 
mri_remap -order vpsczt ${1}_m_q1 
mri_remap -order vpsczt ${1}_m_q3
mri_destroy_dataset ${1}_m_prm

#
# A given sample is considered to be a spike if its "trigger" value
# (computed below) is larger than $spk_cutoff (defined above).
#
mri_rpn_math -out ${1}_is_a_spike '$1,$2,-,$3,$4,-,/,abs,'$spk_cutoff',<' \
    ${1}_m ${1}_m_median ${1}_m_q3 ${1}_m_q1 

#
# Prepare a replacement set of unrolled phase values by linear
# interpoloation between adjacent values.  Since we are assuming
# that spikes are one sample wide, this will provide bridge values
# when a sample must be replaced.
#
mri_complex_to_scalar -phu ${1} ${1}_u
mri_smooth -d p -smoother_type linterp ${1}_u ${1}_u_linterp

#
# By substituting for the spike points, we make filtered values.
#
mri_rpn_math -out $2 -complex \
    '$1,cx_mag,$2,pop,$3,pop,$4,pop,$5,pop,cx_if_keep,cx_fmphase' \
    ${1} ${1}_u ${1}_m_median ${1}_u_linterp ${1}_is_a_spike
mri_destroy_dataset ${1}_u
mri_destroy_dataset ${1}_u_linterp

# Produce spike count information
mri_subsample -sum -d p -l 1 ${1}_is_a_spike ${3}

# Clean up
mri_destroy_dataset ${1}_is_a_spike

