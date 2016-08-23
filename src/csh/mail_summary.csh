#!/bin/csh -ef
# mail_summary.csh
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
echo '#$Id: mail_summary.csh,v 1.8 2007/03/27 01:27:35 welling Exp $'

# Allow the user to avoid sending a summary
if ( ${?F_NO_MAIL_SUMMARY} ) exit 0

# Check for the presence of 'metasend'
set failure = 0
if ( `test_in_subshell.csh which metasend` ) then
  set failure = 1
else
  set words = ( `which metasend` )
  if ( $#words != 1 ) then
    set failure = 1
  endif
endif
if ( $failure ) then
  echo '##### Cannot send summary because metasend is unavailable. #####'
  echo '##### Consider using save_summary.csh instead? #####'
  exit 0
endif

set addr = "fiascosummaries@stat.cmu.edu"

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
  intsplit -splitin $F_SPLIT_FILE  -code $F_SPLIT_COND \
           -splitout $F_SPLIT_NEW $1
endif

echo '##Mailing summary information to ' $addr
summary  -headerinput $1.mri -list $F_SUMM_INPUT \
         -split $F_SPLIT_NEW -cond $F_SPLIT_COND \
         -out tmp_summary.html -fixed $F_ESTIREG_FIXED \
         -missing $F_SUMM_MISSING -env -html
metasend -z -b -f tmp_summary.html -m text/html -s 'Fiasco summary' -t $addr
rm tmp_summary.html

