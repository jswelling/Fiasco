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
echo '#$Id: save_startnote.csh,v 1.1 2002/07/25 23:16:17 welling Exp $'

# Allow the user to avoid sending a summary
if ( ${?F_NO_MAIL_SUMMARY} ) exit 0

if ( ${?F_SUMMARY_DIR} ) then
  set sumdir = ${F_SUMMARY_DIR}
else  
  set sumdir = ~fiasco/summaries
endif
cd $sumdir

# Maintain a counter, to insure filename uniqueness
if (-f counter.t) then
  set count = `cat counter.t`
else
  set count = 0
endif
@ newcount = $count + 1
echo $newcount >! counter.t

# construct unique filename root
set mydate = `date +%Y_%m_%d_%H_%M_%S`
set fname = ${mydate}_${count}.txt

echo # `date` Writing startnote to `pwd`/$fname

cat > $fname <<EOF
args: "$*"
date: "`date`"
user: $LOGNAME
dir: $F_DIR
description: "$F_DESCRIPTION"
header: "$F_HEADER"
EOF



