#!/bin/csh -efx
# spiral.steps.csh
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
echo '#$Id: spiral.steps.csh,v 1.17 2003/09/25 20:16:43 welling Exp $'
if(! -d data) mkdir data
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
mail_startnote.csh spiral
#
# Translate Pfile(s) to Pgh MRI format
#
spiral.reader.csh data/$F_READER_OUTPUT
mark_missing_from_split.csh data/$F_READER_OUTPUT
#
# Noise spike removal
#
spiral.despike.csh data/$F_READER_OUTPUT data/$F_DESPIKE_OUTPUT
mri_destroy_dataset data/$F_READER_OUTPUT
#
# Phase undrift
#
phase_undrift.csh data/$F_DESPIKE_OUTPUT data/$F_PHUND_OUTPUT
mri_destroy_dataset data/$F_DESPIKE_OUTPUT
#
# generate ref files
spiral.refs.csh data/$F_PHUND_OUTPUT
#
# do first reconstruction
spiral.recon1.csh data/$F_PHUND_OUTPUT data/$F_RECON1_OUTPUT
#
# mean and global phase correction
meanc.csh data/$F_RECON1_OUTPUT data/$F_MEANC_OUTPUT
mri_destroy_dataset data/$F_RECON1_OUTPUT
#
# estimate registration parameters
estireg.csh data/$F_MEANC_OUTPUT
#
# smooth registration parameters
parsm.csh data/$F_MEANC_OUTPUT
#
# implement registration in k-space
spiral.recon2.csh \
    data/$F_PHUND_OUTPUT data/$F_RECON2_OUTPUT data/$F_MEANC_OUTPUT
mri_destroy_dataset data/$F_PHUND_OUTPUT
mri_destroy_dataset data/$F_MEANC_OUTPUT
#
# perform outlier correction
outlier.csh data/$F_RECON2_OUTPUT data/$F_OUTLIER_OUTPUT
mri_destroy_dataset data/$F_RECON2_OUTPUT
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
spiral.printpar.csh
#
# Do experiment statistics
$F_STATS_SCRIPT data/$F_DETREND_OUTPUT data/$F_DETREND_TEMP
#
# print postscript files
print.files.csh default all $F_PRINTER

