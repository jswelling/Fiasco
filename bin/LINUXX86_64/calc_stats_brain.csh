#!/bin/csh -exf
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
echo '# $Id: calc_stats_brain.csh,v 1.4 2003/09/09 16:16:25 bakalj Exp $'

#
# This script wants two inputs: the dataset to be analyzed (in image
# order), and the same dataset in permuted form (in time series order).
#

set data = $1
set tsdata = $2

#
# Do we have what we need?
#
if ( $#argv != 2 ) then
  echo '# ' $0 'Error: I need two parameters!'
  exit -1
endif
if ( ( ! -r $1 ) && ( ! -e ${1}.mri ) ) then
  echo '# ' $0 'Error: input data is missing!'
  exit -1
endif
if ( ( ! -r $2 ) && ( ! -e ${2}.mri ) ) then
  echo '# ' $0 'Error: permuted input data is missing!'
  exit -1
endif
if ( ! -r $F_SPLIT_FILE ) then
  echo '# ' $0 'Error: split file is missing!'
  exit -1
endif

# Parse split file if an output is missing, or if it is 
# more recent than an output
set split_new = $F_SPLIT_NEW
if (! -d ${split_new:h} ) mkdir ${split_new:h}
set split_cond = $F_SPLIT_COND
if (! -d ${split_cond:h} ) mkdir ${split_cond:h}
set do_split = 0
if (! -r $F_SPLIT_NEW || ! -r $F_SPLIT_COND ) then
  set do_split = 1
else
  set ages = `ls -t $F_SPLIT_FILE $F_SPLIT_NEW $F_SPLIT_COND`
  if ($ages[3] != $F_SPLIT_FILE) set do_split = 1
endif
if ($do_split) then
  intsplit -splitin $F_SPLIT_FILE -code $F_SPLIT_COND \
           -splitout $F_SPLIT_NEW $1
endif

#
# calculate Bayesian analysis
cg_brain.csh $data data/$F_VMPFX_OUTPUT data/$F_VMPFX_TEMP

#
# extract the interesting parts of the results
cg_brain_stats.csh data/$F_VMPFX_OUTPUT

#
# create postscript files
cg_brain.makeps.csh

