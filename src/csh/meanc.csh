#!/bin/csh -efx
# epi.meanc.csh
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
echo '# $Id: meanc.csh,v 1.8 2003/12/04 21:20:25 welling Exp $'
if(! -d par) mkdir par
meanc          -input $1.mri -headerout $2.mri \
                 -parameters par/$F_MEANC_PARMS.$$ \
                 -fixed $F_MEANC_FIXED
set step = `depath.csh $0`
echo  "$step par/$F_MEANC_PARMS.$$" >> $F_SUMM_INPUT
echo $step.$$ >> $F_SUMM_MISSING
count_missing.csh $2.mri >> $F_SUMM_MISSING
