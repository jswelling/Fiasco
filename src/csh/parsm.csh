#!/bin/csh -efx
# parsm.csh
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
echo '#$Id: parsm.csh,v 1.7 2002/08/07 00:01:48 welling Exp $'

set inparms = `most_recent_parms.py -f estireg`

smregpar -headerinput $1.mri \
           -parameterin $inparms \
           -parameterout par/$F_PARSM_PARMS.$$ -cutoff $F_PARSM_CUTOFF \
           -fixed $F_ESTIREG_FIXED -bandwidth $F_PARSM_BANDWIDTH \
           -kernel $F_PARSM_KERNEL \
           -threshold $F_PARSM_TRANST
if(! -d par) mkdir par
set step = `depath.csh $0`
echo "$step par/$F_PARSM_PARMS.$$" >> $F_SUMM_INPUT
echo $step.$$ >> $F_SUMM_MISSING
count_missing.csh $1.mri >> $F_SUMM_MISSING

