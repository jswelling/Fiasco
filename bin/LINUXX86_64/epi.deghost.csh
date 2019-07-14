#!/bin/csh -efx
# epi.deghost.csh
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
echo '$Id: epi.deghost.csh,v 1.20 2005/06/01 19:59:36 welling Exp $'
#
echo '#'`date`$0
if(! -d par) mkdir par

# Assume that the input dataset needs deghosting unless it explicitly
# says otherwise.
set needs_deghost = `mri_printfield -fld images.needs_deghost -nofail $1`
if ( "dummy"$needs_deghost == "dummy" ) then
  set needs_deghost = 1
endif

if ( $?F_READER_SPACE ) then
  set reader_space = $F_READER_SPACE
else
  set reader_space = k
endif

if ( $reader_space == "i" ) then
  echo "Deghosting skipped because input data was in image space"
  mri_copy_dataset $1 $2
else if ( ! $needs_deghost ) then
  echo "Deghosting skipped because input data has needs_deghost=0 attribute"
  mri_copy_dataset $1 $2  
else
  set do_flip = \
    `mri_printfield -field images.rowflip -nofail $1`
  set navchunk = "`mri_printfield -fld navigator -nofail $1`"

  if ( dummy"${navchunk}" == 'dummy[chunk]' ) then

    set nav_opts = "-estimates par/$F_NAVDGH_RWPARMS.$$ "
    set nav_opts = "$nav_opts -smoothedestimates par/$F_NAVDGH_PARMS.$$ "
    set nav_opts = "$nav_opts -smoother_bandwidth $F_NAVDGH_BAND "

    if ( dummy${do_flip} != dummy ) then
      # Input dataset specifies row flipping (implemented by the executable)
    else
      # Implement row flipping based on environment variables
      set nav_opts = "$nav_opts \
         -reverse $F_NAVDGH_REVERSE -navreverse $F_NAVDGH_NAVREVERSE"
    endif

    parallel.run.csh nav_deghost -verbose $nav_opts $1 $2
    echo "****WARNING*** deleting navigator chunk"
    mri_delete_chunk -chunk navigator $2
    set step = `depath.csh $0`
    echo "${step}_nav_raw par/$F_NAVDGH_RWPARMS.$$" >> $F_SUMM_INPUT
    echo "${step}_nav par/$F_NAVDGH_PARMS.$$" >> $F_SUMM_INPUT

  else

    set deghost_opts = "-estimates par/$F_DEGHOST_RWPARMS.$$ "
    set deghost_opts = \
      " $deghost_opts -smoothedestimates par/$F_DEGHOST_PARMS.$$ "
    set deghost_opts = \
      " $deghost_opts -smoother_bandwidth $F_DEGHOST_BAND "
    set deghost_opts = \
      " $deghost_opts -phase $F_DEGHOST_PHASE "

    set deghost_est_opts = "$deghost_opts"
    set deghost_apply_opts = "$deghost_opts"

    if ( dummy${do_flip} != dummy ) then
      # Input dataset specifies row flipping (implemented by the executable)
    else
      # Implement row flipping based on environment variables.
      # We need to keep careful track of row flipping- the deghost 'apply'
      # stage is performed on the original input to baseline, so
      # it must use baseline's row flipping information.
      set deghost_apply_opts = \
        " $deghost_apply_opts -reverse $F_BASELINE_REVERSE "
    endif

    epi.baseline.csh $1 data/$F_DEGHOST_TEMP
    parallel.run.csh deghost -mode estimate $deghost_est_opts \
      data/$F_DEGHOST_TEMP dummy
    parallel.run.csh deghost -mode apply $deghost_apply_opts \
      $1 $2
    mri_destroy_dataset data/$F_DEGHOST_TEMP

    set step = `depath.csh $0`
    echo "${step}_raw par/$F_DEGHOST_RWPARMS.$$" >> $F_SUMM_INPUT
    echo "$step par/$F_DEGHOST_PARMS.$$" >> $F_SUMM_INPUT

  endif
  mri_setfield -fld 'images.needs_deghost' -value 0 $2
  echo $step.$$ >> $F_SUMM_MISSING
  count_missing.csh $2.mri >> $F_SUMM_MISSING
endif

