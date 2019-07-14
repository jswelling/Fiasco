#!/bin/csh -efx
# partialk.csh
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
echo '$Id: partialk.csh,v 1.9 2003/09/10 20:17:23 bakalj Exp $'
#
echo '#'`date`$0
if(! -d par) mkdir par
set do_partialk = `mri_printfield -field images.partialk -nofail $1`
set do_flip = `mri_printfield -field images.rowflip -nofail $1`

if ( dummy$do_partialk != dummy ) then

  if ( dummy$do_partialk == dummy1 ) then

    if ( dummy${do_flip} != dummy ) then
      # Input dataset specifies row flipping (implemented by the executable)
      parallel.run.csh partialk -input $1 -output $2
    else
      # Implement row flipping based on environment variables
      parallel.run.csh partialk -reverse $F_PARTIALK_REVERSE -input $1 \
        -output $2
    endif

  else

    # Partial-k correction not needed
    mri_copy_dataset $1 $2

  endif

else

  if ( $F_CLIP1_YCENTER != "center" ) then
    @ newdim = 2 * $F_CLIP1_YCENTER
  else
    @ newdim = `mri_printfield -field "images.extent.y" -nofail $1`
  endif

  if ( dummy${do_flip} != dummy ) then
    # Input dataset specifies row flipping (implemented by the executable)
    parallel.run.csh partialk -dim $newdim -input $1 -output $2
  else
    # Implement row flipping based on environment variables
    parallel.run.csh partialk -reverse $F_PARTIALK_REVERSE -dim $newdim \
      -input $1 -output $2
  endif

endif


