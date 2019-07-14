#!/bin/csh -ef
# epi.steps.csh
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
echo '# $Id: epi.steps.csh,v 1.18 2003/08/07 21:16:06 welling Exp $'
if (! -d data) mkdir data
#
# Start the profiler is appropriate
#
if ( ${?F_WATCHER_INTERVAL} ) then
  if ( ${F_WATCHER_INTERVAL} > 0 ) then
    nohup watcher.py -i ${F_WATCHER_INTERVAL} $$ >> ${F_WATCHER_LOG} &
  endif
endif
#
# Send a note home indicating we started a run
mail_startnote.csh epi
#
# read in the data into our standard format
epi.reader.csh $F_READER_INPUT data/$F_READER_OUTPUT
mark_missing_from_split.csh data/$F_READER_OUTPUT
#
# perform deghosting
epi.deghost.csh data/$F_READER_OUTPUT data/$F_DEGHOST_OUTPUT
mri_destroy_dataset data/$F_READER_OUTPUT
#
# perform ramp resampling if necessary
ramp_resample.csh data/$F_DEGHOST_OUTPUT data/$F_RESAMP_OUTPUT
mri_destroy_dataset data/$F_DEGHOST_OUTPUT
#
# perform baseline correction
epi.baseline.csh data/$F_RESAMP_OUTPUT data/$F_BASELINE_OUTPUT
mri_destroy_dataset data/$F_RESAMP_OUTPUT
#
# partialk correction
partialk.csh data/$F_BASELINE_OUTPUT data/$F_PARTIALK_OUTPUT
mri_destroy_dataset data/$F_BASELINE_OUTPUT
#
# mean correction
meanc.csh data/$F_PARTIALK_OUTPUT data/$F_MEANC_OUTPUT
#
# reconstruct corrected images
epi.recon1.csh data/$F_MEANC_OUTPUT data/$F_RECON1_OUTPUT
mri_destroy_dataset data/$F_MEANC_OUTPUT
#
# cut them down to size to speed up computations
epi.clip1.csh data/$F_RECON1_OUTPUT data/$F_CLIP1_OUTPUT
mri_destroy_dataset data/$F_RECON1_OUTPUT
#
# estimate registration parameters
estireg.csh data/$F_CLIP1_OUTPUT
#
# smooth registration parameters
parsm.csh data/$F_CLIP1_OUTPUT
#
# implement registration in k-space
epi.ireg.csh data/$F_PARTIALK_OUTPUT data/$F_IREG_OUTPUT
mri_destroy_dataset data/$F_PARTIALK_OUTPUT
#
# Physiological data correction
physio_correct.csh data/$F_IREG_OUTPUT data/$F_PHYS_OUTPUT \
   data/$F_CLIP1_OUTPUT $F_PHYS_DATA
mri_destroy_dataset data/$F_IREG_OUTPUT
#
# reconstruct corrected, registered images
epi.recon2.csh data/$F_PHYS_OUTPUT data/$F_RECON2_OUTPUT data/$F_CLIP1_OUTPUT
mri_destroy_dataset data/$F_PHYS_OUTPUT
mri_destroy_dataset data/$F_CLIP1_OUTPUT
#
# cut them down to size to speed up computations
epi.clip2.csh data/$F_RECON2_OUTPUT data/$F_CLIP2_OUTPUT
mri_destroy_dataset data/$F_RECON2_OUTPUT
#
# perform outlier correction
outlier.csh data/$F_CLIP2_OUTPUT data/$F_OUTLIER_OUTPUT
mri_destroy_dataset data/$F_CLIP2_OUTPUT
#
#remove trend pixelwise
detrend.csh data/$F_OUTLIER_OUTPUT data/$F_DETREND_OUTPUT data/$F_DETREND_TEMP
mri_destroy_dataset data/$F_OUTLIER_OUTPUT
#
# calculate average pixel displacement
displace.csh data/$F_DETREND_OUTPUT
#
# print summary pages
summary.csh data/$F_DETREND_OUTPUT
#
# mail summary to Fiasco Headquarters
mail_summary.csh data/$F_DETREND_OUTPUT
#
# create plots of shifts, rotation, and displacement
epi.printpar.csh
#
# Do experiment statistics
$F_STATS_SCRIPT data/$F_DETREND_OUTPUT data/$F_DETREND_TEMP
#
# print postscript files
print.files.csh default all $F_PRINTER




#mv output.dat input.dat
#mv output.mri input.mri
#chnghead.m input.mri -string images.file input.dat

