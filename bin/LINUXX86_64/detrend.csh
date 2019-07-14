#!/bin/csh -efx
# detrend.csh
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
echo '#$Id: detrend.csh,v 1.13 2005/06/21 21:15:40 welling Exp $'
if (! -d par) mkdir par
if (`mri_printfield -fld ${F_DETREND_CHUNK}.dimensions $1` == "xyzt") then
	mri_remap -chunk ${F_DETREND_CHUNK} -order vxyzt $1
endif
mri_permute -memlimit 32000000 -chunk ${F_DETREND_CHUNK} -order vtxyz $1.mri $2.mri
echo '#'`date`
set F_T = $3:t
detrend -estimates par/$F_DETREND_PARHEAD $2 $3
echo '#'`date`
set F_T = $2:t
mri_permute -memlimit 32000000 -chunk ${F_DETREND_CHUNK} -order vxyzt $3.mri $2.mri
echo '#'`date`







