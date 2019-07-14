#!/bin/csh -efx
# epi.baseline.csh
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
echo '$Id: epi.baseline.csh,v 1.15 2003/09/09 16:25:45 bakalj Exp $'
#
echo '#'`date`$0
if(! -d par) mkdir par

set do_flip = \
  `mri_printfield -field images.rowflip -nofail $1`
if ( dummy${do_flip} != dummy ) then
  # Input dataset specifies row flipping (implemented by the executable)
  parallel.run.csh baseline  -estimates par/$F_BASELINE_RWPARMS.$$ \
	                     -smoothedestimates par/$F_BASELINE_PARMS.$$ \
                             $F_BASELINE_SMOPTS \
	                     $1.mri $2.mri
else
  # Implement row flipping based on environment variables
  if ( $F_READER_CORONAL ) then
    if ( $F_BASELINE_REVERSE == "even" ) then
      set use_reverse = "odd"
    else if ( $F_BASELINE_REVERSE == "odd" ) then
      set use_reverse = "even"
    else
      set use_reverse = $F_BASELINE_REVERSE
    endif
  else
    set use_reverse = $F_BASELINE_REVERSE
  endif
  parallel.run.csh baseline -estimates par/$F_BASELINE_RWPARMS.$$ \
	                    -smoothedestimates par/$F_BASELINE_PARMS.$$ \
                            $F_BASELINE_SMOPTS \
	                    -reverse $use_reverse \
	                    $1.mri $2.mri
endif

set step = `depath.csh $0`
echo "${step}_raw par/$F_BASELINE_RWPARMS.$$" >> $F_SUMM_INPUT
echo "$step par/$F_BASELINE_PARMS.$$" >> $F_SUMM_INPUT
echo $step.$$ >> $F_SUMM_MISSING
count_missing.csh $2.mri >> $F_SUMM_MISSING

