#!/bin/csh -efx
# spiral.recon2.csh
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
echo '#'`date` $0
echo '#$Id: spiral.recon2.csh,v 1.15 2003/09/10 20:17:23 bakalj Exp $'

# This ridiculous bit of gymnastics is necessary because the command
# parser in Spiral wants a comma-separated list, but the csh scripts
# that got us here want a space-separated list.
if ($F_PARALLEL) then
  set par_string = "-hosts "
  foreach host ( $F_PARALLEL_HOSTS )
    set par_string = "$par_string$host,"
  end
else
  set par_string = ""
endif

set inparms = `most_recent_parms.py -f parsm`

parallel.run.csh spiral $F_SPIRAL_OPT \
       -res $F_SPIRAL_RESOLUTION -samp_delay $F_SPIRAL_SAMP_DEL \
       -phase_factor $F_SPIRAL_PHFAC -phase_twist $F_SPIRAL_PHTW \
       -filter_sz $F_SPIRAL_FILT_SZ -scale $F_SPIRAL_SCALE \
       -input_directory ./ \
       -reference_directory data -tmp_directory $F_TEMP \
       -out_name $2 -output_dir ./ -print_mode \
       -reg_file ${inparms} $par_string \
       -tag $F_DIR/$F_REFP $1

# Need to copy "missing" chunk from dataset used for motion estimation
# to new reconstructed dataset
mri_copy_chunk -chunk missing -replace $3 $2

# Reference files can get marked missing at this step.
if(! -d par) mkdir par
set step = `depath.csh $0`
echo $step.$$ >> $F_SUMM_MISSING
count_missing.csh $2 >> $F_SUMM_MISSING


