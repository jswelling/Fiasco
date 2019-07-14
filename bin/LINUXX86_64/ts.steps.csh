#!/bin/csh -ef
# ts.steps.csh
#/************************************************************
# *                                                          *
# *  Permission is hereby granted to any individual or       *
# *  institution for use, copying, or redistribution of      *
# *  this code and associated documentation, provided        *
# *  that such code and documentation are not sold for       *
# *  profit and the following copyright notice is retained   *
# *  in the code and documentation:                          *
# *  Copyright (c) 1995,1996, 1997 Department of Statistics, *
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
echo '# $Id: ts.steps.csh,v 1.25 2002/01/08 01:44:31 welling Exp $'
if (! -d data) mkdir data
#
# Send a note home indicating we started a run
mail_startnote.csh twoshot
#
# read in the data into our standard format
ts.reader.csh $F_READER_INPUT data/$F_READER_OUTPUT
#
# baseline correct individual shots
ts.baseline.csh data/$F_READER_OUTPUT data/$F_BASELINE_OUTPUT
mri_destroy_dataset data/$F_READER_OUTPUT
#
# partialk correction
partialk.csh data/$F_BASELINE_OUTPUT data/$F_PARTIALK_OUTPUT
mri_destroy_dataset data/$F_BASELINE_OUTPUT
#
# perform deghosting if desired (don't forget to change the
# first parameter of baseline2.csh to data/$F_DEGHOST_OUTPUT
# ts.deghost.csh data/$F_PARTIALK_OUTPUT data/$F_DEGHOST_OUTPUT
# mri_destroy_dataset data/$F_PARTIALK_OUTPUT
#
# merge pairs and repeat baseline correction (harmless)
ts.baseline2.csh data/$F_PARTIALK_OUTPUT data/$F_BASELINE2_OUTPUT
mri_destroy_dataset data/$F_PARTIALK_OUTPUT
#
# Repeat deghosting step with merged images if desired (don't forget
# to change the first parameters of meanc.csh and ts.ireg.csh
# to data/$F_DEGHOST2_OUTPUT
# ts.deghost2.csh data/$F_BASELINE2_OUTPUT data/$F_DEGHOST2_OUTPUT
# mri_destroy_dataset data/$F_BASELINE2_OUTPUT
#
# mean correction
meanc.csh data/$F_BASELINE2_OUTPUT data/$F_MEANC_OUTPUT
#
# reconstruct corrected images
ts.recon1.csh data/$F_MEANC_OUTPUT data/$F_RECON1_OUTPUT
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
ts.ireg.csh data/$F_BASELINE2_OUTPUT data/$F_IREG_OUTPUT
mri_destroy_dataset data/$F_BASELINE2_OUTPUT
#
# Physiological data correction
physio_correct.csh data/$F_IREG_OUTPUT data/$F_PHYS_OUTPUT \
    data/$F_CLIP1_OUTPUT $F_PHYS_DATA
mri_destroy_dataset data/$F_IREG_OUTPUT
#
# reconstruct corrected, registered images
ts.recon2.csh data/$F_PHYS_OUTPUT data/$F_RECON2_OUTPUT data/$F_CLIP1_OUTPUT
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
ts.printpar.csh
#
# Do experiment statistics
$F_STATS_SCRIPT data/$F_DETREND_OUTPUT data/$F_DETREND_TEMP
#
# print postscript files
print.files.csh default all $F_PRINTER




#mv output.dat input.dat
#mv output.mri input.mri
#chnghead.m input.mri -string images.file input.dat

