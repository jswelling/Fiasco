#!/bin/csh -efx
# mail_startnote.csh
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
echo '#$Id: save_summary.csh,v 1.3 2003/09/10 20:17:23 bakalj Exp $'

# Allow the user to avoid saving a summary
if ( ${?F_NO_MAIL_SUMMARY} ) exit 0

if ( ${?F_SUMMARY_DIR} ) then
  set sumdir = ${F_SUMMARY_DIR}
else  
  set sumdir = ~fiasco/summaries
endif

# Maintain a counter, to insure filename uniqueness
if (-f ${sumdir}/counter.t) then
  set count = `cat ${sumdir}/counter.t`
else
  set count = 0
endif
@ newcount = $count + 1
echo $newcount >! ${sumdir}/counter.t

# construct unique filename root
set mydate = `date +%Y_%m_%d_%H_%M_%S`
set fname = ${mydate}_${count}.html

echo # `date` Writing summary to $sumdir/$fname

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

summary  -headerinput $1.mri -list $F_SUMM_INPUT \
         -split $F_SPLIT_NEW -cond $F_SPLIT_COND \
         -out $sumdir/$fname -fixed $F_ESTIREG_FIXED \
         -missing $F_SUMM_MISSING -env -html




