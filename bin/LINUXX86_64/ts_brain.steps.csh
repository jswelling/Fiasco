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
echo '# $Id: ts_brain.steps.csh,v 1.15 2002/01/08 01:45:58 welling Exp $'
if (! -d data) mkdir data
#
# Send a note home indicating we started a run
mail_startnote.csh twoshot_brain
#
# read in the data into our standard format
ts.reader.csh $F_READER_INPUT data/$F_READER_OUTPUT
#
# baseline correct individual shots
ts.baseline.csh data/$F_READER_OUTPUT data/$F_BASELINE_OUTPUT
#
# partialk correction
partialk.csh data/$F_BASELINE_OUTPUT data/$F_PARTIALK_OUTPUT
#
# perform deghosting if desired (don't forget to change the
# first parameter of baseline2.csh to data/$F_DEGHOST_OUTPUT
# ts.deghost.csh data/$F_PARTIALK_OUTPUT data/$F_DEGHOST_OUTPUT
#
# merge pairs and repeat baseline correction (harmless)
ts.baseline2.csh data/$F_PARTIALK_OUTPUT data/$F_BASELINE2_OUTPUT
#
# Repeat deghosting step with merged images if desired (don't forget
# to change the first parameter of meanc.csh and ts.ireg.csh 
# to data/$F_DEGHOST2_OUTPUT
# ts.deghost2.csh data/$F_BASELINE2_OUTPUT data/$F_DEGHOST2_OUTPUT
#
# mean correction
meanc.csh data/$F_BASELINE2_OUTPUT data/$F_MEANC_OUTPUT
#
# reconstruct corrected images
ts.recon1.csh data/$F_MEANC_OUTPUT data/$F_RECON1_OUTPUT
#
# cut them down to size to speed up computations
epi.clip1.csh data/$F_RECON1_OUTPUT data/$F_CLIP1_OUTPUT
#
# estimate registration parameters
estireg.csh data/$F_CLIP1_OUTPUT
#
# smooth registration parameters
parsm.csh data/$F_CLIP1_OUTPUT
#
# implement registration in k-space
ts.ireg.csh data/$F_BASELINE2_OUTPUT data/$F_IREG_OUTPUT
#
# Physiological data correction
physio_correct.csh data/$F_IREG_OUTPUT data/$F_PHYS_OUTPUT \
    data/$F_CLIP1_OUTPUT $F_PHYS_DATA
#
# reconstruct corrected, registered images
ts.recon2.csh data/$F_PHYS_OUTPUT data/$F_RECON2_OUTPUT data/$F_CLIP1_OUTPUT
#
# cut them down to size to speed up computations
epi.clip2.csh data/$F_RECON2_OUTPUT data/$F_CLIP2_OUTPUT
#
# perform outlier correction
outlier.csh data/$F_CLIP2_OUTPUT data/$F_OUTLIER_OUTPUT
#
# separate into experimental conditions
split.csh data/$F_OUTLIER_OUTPUT
#
# calculate Bayesian analysis
cg_brain.csh data/$F_OUTLIER_OUTPUT data/$F_VMPFX_OUTPUT data/$F_VMPFX_TEMP
#
# extract the interesting parts of the results
cg_brain_stats.csh data/$F_VMPFX_OUTPUT
#
# calculate average pixel displacement
displace.csh
#
# print summary pages
summary.csh data/$F_VMPFX_TEMP
#
# mail summary to Fiasco Headquarters
mail_summary.csh data/$F_VMPFX_TEMP
#
# create plots of shifts, rotation, and displacement
ts.printpar.csh
#
# create postscript files
cg_brain.makeps.csh
#
# print postscript files
print.files.csh default all $F_PRINTER




#mv output.dat input.dat
#mv output.mri input.mri
#chnghead.m input.mri -string images.file input.dat
