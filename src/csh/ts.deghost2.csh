#!/bin/csh -efx
# ts.deghost2.csh
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
echo '$Id: ts.deghost2.csh,v 1.10 2003/07/18 21:27:59 welling Exp $'
#
echo '#'`date`$0
if(! -d par) mkdir par
parallel.run.csh deghost -input $1.mri -headerout $2.mri -dataout .dat \
  -parameters par/$F_DEGHOST2_RWPARMS.$$ \
  -smoothedparameters par/$F_DEGHOST2_PARMS.$$ \
  -smoother_bandwidth $F_DEGHOST2_BAND \
  -phase $F_DEGHOST2_PHASE -reverse none
set step = `depath.csh $0`
echo "${step}_raw par/$F_DEGHOST2_RWPARMS.$$" >> $F_SUMM_INPUT
echo "$step par/$F_DEGHOST2_PARMS.$$" >> $F_SUMM_INPUT



