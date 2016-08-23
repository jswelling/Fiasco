#!/bin/csh -efx
# stats_only.proc.csh
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
# $Id: stats_only.proc.csh,v 1.2 2007/06/21 23:13:20 welling Exp $
echo '#'`date` $0

# Infer the type of run from the .local.csh files
set runtype_candidates = 0
if (-e epi.local.csh) then
  set runtype = epi
  @ runtype_candidates = $runtype_candidates + 1
else if (-e spiral.local.csh) then
  set runtype = spiral
  @ runtype_candidates = $runtype_candidates + 1
else if (-e ts.local.csh) then
  set runtype = ts
  @ runtype_candidates = $runtype_candidates + 1
endif
if ( $runtype_candidates != 1 ) then
  echo '%%%%% Error- cannot determine run type! %%%%%'
  exit -1
endif

# Source appropriate stuff
if(! -e fiasco.local.csh) cp $FIASCO/fiasco.local.csh .
source fiasco.local.csh
if (-e ${runtype}.defaults.csh) then
  source ${runtype}.defaults.csh
else
  source ${FIASCO}/${runtype}.defaults.csh
endif
source ${runtype}.local.csh

if (! -d $F_TEMP) mkdir $F_TEMP
if(! -d out) mkdir out
if(! -d ps) mkdir ps

# Make sure to use the local script if it is present, watching out
# for some common user errors.
set rawScriptName = stats_only.steps.csh
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


