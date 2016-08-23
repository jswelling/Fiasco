#!/bin/csh -efx
# ge_ramp_resample.csh
#/************************************************************
# *                                                          *
# *  Permission is hereby granted to any individual or       *
# *  institution for use, copying, or redistribution of      *
# *  this code and associated documentation, provided        *
# *  that such code and documentation are not sold for       *
# *  profit and the following copyright notice is retained   *
# *  in the code and documentation:                          *
# *     Copyright (c) 2003 Department of Statistics,         *
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
# *  Original programming by Joel Welling                    *
# ************************************************************/
#
echo '$Id: ramp_resample.csh,v 1.4 2004/09/10 00:27:49 welling Exp $'
#
echo '#'`date`$0

# This script just sets up and dispatches ramp resampling if
# necessary.  Specialized scripts handle the specific cases.

set infile = ${1}
set outfile = ${2}

# ramp resampling if necessary.  
set resample = \
  `mri_printfield -field images.resample -nofail ${infile}`
if (${#resample} != 0) then
  if (${resample} != 0) then
    set resamp_scr = \
      `mri_printfield -field images.resample_method -nofail ${infile}`
    if ( ${#resamp_scr} == 0 ) then
      echo "****ERROR**** ramp resampling requested, but no method found"
      exit -1  
    endif
    if ( -e ./${resamp_scr} ) then
      echo "Ramp resampling data using ./${resamp_scr}"
      ./${resamp_scr} ${infile} ${outfile}
      echo "Ramp resampling complete"
    else if ( -e ${FIASCO}/${resamp_scr} ) then
      echo "Ramp resampling data using ${FIASCO}/${resamp_scr}"
      ${FIASCO}/${resamp_scr} ${infile} ${outfile}
      echo "Ramp resampling complete"
    else
      echo "****ERROR**** could not find resampling script ${resamp_scr}"
      exit -1  
    endif
    mri_setfield -field images.resample -value 0 ${outfile} 
    mri_setfield -field images.description.x -value 'gridded k-space' \
	${outfile}
  else
    mri_copy_dataset ${infile} ${outfile}
  endif
else
  mri_copy_dataset ${infile} ${outfile}
endif

