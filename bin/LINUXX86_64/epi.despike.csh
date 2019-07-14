#! /bin/csh -exf
# epi.despike.csh
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
echo '#$Id: epi.despike.csh,v 1.6 2004/02/03 20:54:29 welling Exp $'

# Find out some dimensions
set vdim = `mri_printfield -field images.extent.v $1`
set xdim = `mri_printfield -field images.extent.x $1`
set ydim = `mri_printfield -field images.extent.y $1`
set zdim = `mri_printfield -field images.extent.z $1`
set tdim = `mri_printfield -field images.extent.t $1`

@ pdim = ( $xdim * $ydim )

# Make a working copy, and give it the needed shape.
if ( ! -d data ) mkdir data
mri_type_convert -float $1 data/raw_samples
mri_remap -order vpsczt -length ${vdim}:${pdim}:1:1:${zdim}:${tdim} \
    data/raw_samples 

# Despike the data by comparing across times.
despike_timewise.csh data/raw_samples data/samples_filtered data/spike_count
mri_destroy_dataset data/raw_samples

# Produce an output dataset
mri_copy_dataset $1 $2
mri_remap -order vxyzt -length ${vdim}:${xdim}:${ydim}:${zdim}:${tdim} \
    data/samples_filtered
mri_copy_chunk -chunk images -chunk_out images -replace data/samples_filtered $2
mri_destroy_dataset data/samples_filtered

# Produce summary information
if(! -d par) mkdir par
set step = `depath.csh $0`
echo '##Format: order:index_tz type:raw names:(spike_count)' \
    > par/$F_DESPIKE_PARMS.$$
mri_rpn_math -out data/junk '0,$t,$z,$1,1,if_print_3' \
    data/spike_count >> par/$F_DESPIKE_PARMS.$$
echo "$step par/${F_DESPIKE_PARMS}.$$" >> $F_SUMM_INPUT
mri_destroy_dataset data/junk
mri_destroy_dataset data/spike_count

