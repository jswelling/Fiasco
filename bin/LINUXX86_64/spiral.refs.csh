#!/bin/csh -efx
# spiral.refs.csh
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
echo '#$Id: spiral.refs.csh,v 1.12 2003/07/18 21:27:59 welling Exp $'

# The following ridiculous bit of gymnastics is necessary because the 
# command parser in Spiral wants a comma-separated list, but the csh 
# scripts that got us here want a space-separated list.
if ($F_PARALLEL) then
  set par_string = "-hosts "
  foreach host ( $F_PARALLEL_HOSTS )
    set par_string = "$par_string$host,"
  end
else
  set par_string = ""
endif

parallel.run.csh spiral $F_SPIRAL_REFOPT \
       -res $F_SPIRAL_RESOLUTION -samp_delay $F_SPIRAL_SAMP_DEL \
       -filter_sz $F_SPIRAL_FILT_SZ -scale $F_SPIRAL_SCALE \
       -phase_factor $F_SPIRAL_PHFAC -phase_twist $F_SPIRAL_PHTW \
       -input_dir ./ -tmp_dir $F_TEMP \
       -output_dir data -print_mode $par_string \
       -tag $F_DIR/$F_REFP \
       $1



