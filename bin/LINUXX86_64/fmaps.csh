#!/bin/csh -fx
# fmaps.csh
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

echo '#'`date`$0
echo '#$Id: fmaps.csh,v 1.10 2003/09/09 16:25:45 bakalj Exp $'
if (! -d stat) mkdir stat
cd stat
setenv F_INPUTFILE ../$1
setenv F_SPLIT_CFL ../$F_SPLIT_COND
setenv F_SPLIT_FTZ ../$F_SPLIT_NEW
mri_anova -cnd $F_SPLIT_CFL -split $F_SPLIT_FTZ $F_INPUTFILE
# If anovas failed, clean up the incorrect Fmaps files
if ($status) then
  echo "Anovas failed;  cleaning up Fmaps files."
  rm -f Fmaps*.mri
  rm -f Fmaps*.dat
endif


