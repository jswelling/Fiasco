#!/bin/csh -efx
# ts.proc.csh
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
# $Id: ts.proc.csh,v 1.7 2007/06/21 23:13:20 welling Exp $
echo '#'`date` $0
if(! -e fiasco.local.csh) cp $FIASCO/fiasco.local.csh .
source fiasco.local.csh
if (-e ts.defaults.csh) then
  source ts.defaults.csh
else 
  source $FIASCO/ts.defaults.csh
endif
if(! -e ts.local.csh) cp $FIASCO/ts.local.csh .
source ts.local.csh 
if (! -d $F_TEMP) mkdir $F_TEMP
if(! -d out) mkdir out
if(! -d ps) mkdir ps
if($F_PARALLEL) source $FIASCO/parallel.start.csh
if ($status) exit ($status)

# Make sure to use the local script if it is present, watching out
# for some common user errors.
set rawScriptName = ts.steps.csh
if ( -f $rawScriptName ) then
  if ( -x $rawScriptName ) then
    ./$rawScriptName
  else
    echo "A local $rawScriptName exists, but it is not executable!"
    exit -1
  endif
else
  # Use whatever version is first in path
  $rawScriptName
endif

if($F_PARALLEL) source $FIASCO/parallel.finish.csh
