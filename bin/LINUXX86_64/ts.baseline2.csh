#!/bin/csh -efx
# ts.baseline2.csh
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
echo '$Id: ts.baseline2.csh,v 1.10 2003/05/14 18:18:32 welling Exp $'
#
echo '#'`date`$0
if(! -d par) mkdir par
set do_flip = \
  `mri_printfield -input $1 -field images.rowflip -nofail`
if ( dummy${do_flip} != dummy ) then
  # Input dataset specifies row flipping (implemented by the executable)
  baseline2      -f $1.mri -s $2.mri -w par/$F_BASELINE2_PARMS.$$ \
          -o $F_BASELN2_OVER
else
  # Implement row flipping based on environment variables
  baseline2      -f $1.mri -s $2.mri -w par/$F_BASELINE2_PARMS.$$ \
          -r $F_BASELN2_REVERSE -o $F_BASELN2_OVER
endif
set step = `depath.csh $0`
echo "$step par/$F_BASELINE2_PARMS.$$" >> $F_SUMM_INPUT
echo $step.$$ >> $F_SUMM_MISSING
count_missing.csh $2.mri >> $F_SUMM_MISSING

