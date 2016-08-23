#! /bin/csh -exf
# spiral.despike.csh
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
echo '#$Id: spiral.despike.csh,v 1.15 2003/09/09 16:25:45 bakalj Exp $'

set tdim = `mri_printfield -field samples.extent.t $1`
@ normimages = $tdim - 1

# Make a working copy, separating t=0 from the rest of the times
mri_copy_chunk -chunk samples -chunk_out images -replace \
    $1 data/raw_samples_all
mri_type_convert -float data/raw_samples_all data/raw_samples_all_float
mri_destroy_dataset data/raw_samples_all
mri_subset -d t -l 1 -s 0 data/raw_samples_all_float data/raw_samples_t0
mri_subset -d t -l $normimages -s 1 data/raw_samples_all_float data/raw_samples
mri_destroy_dataset data/raw_samples_all_float

# Despike the bulk of the data by comparing across times.
# This produces some files that we want as side effects.
despike_timewise.csh data/raw_samples data/samples_filtered data/spike_count

if (! -f data/raw_samples_m_median.mri) then
  echo "despike_timewise.csh failed to produce magnitude median!"
  exit -1
endif
if (! -f data/raw_samples_m_q1.mri) then
  echo "despike_timewise.csh failed to produce magnitude q1!"
  exit -1
endif
if (! -f data/raw_samples_m_q3.mri) then
  echo "despike_timewise.csh failed to produce magnitude q3!"
  exit -1
endif
mri_destroy_dataset data/raw_samples

# Despike image 0 by comparing across samples.
despike_sampwise.csh data/raw_samples_t0 \
    data/raw_samples_m_q1 data/raw_samples_m_q3 \
    data/samples_filtered_t0 data/spike_count_t0
mri_destroy_dataset data/raw_samples_t0

# Produce an output dataset
mri_paste -out data/samples_filtered_all -d t \
    data/samples_filtered_t0 data/samples_filtered
mri_destroy_dataset data/samples_filtered
mri_destroy_dataset data/samples_filtered_t0
mri_copy_dataset $1 $2
mri_copy_chunk -chunk images -chunk_out samples -replace \
    data/samples_filtered_all $2
mri_destroy_dataset data/samples_filtered_all

# Produce summary information
if(! -d par) mkdir par
mri_paste -out data/spike_count_all data/spike_count_t0 data/spike_count
set step = `depath.csh $0`
mri_destroy_dataset data/spike_count
mri_destroy_dataset data/spike_count_t0
echo '##Format: order:index_tz type:raw names:(spike_count)' \
    > par/$F_DESPIKE_PARMS.$$
mri_rpn_math -out data/junk '0,$t,$z,$1,1,if_print_3' \
    data/spike_count_all >> par/$F_DESPIKE_PARMS.$$
echo "$step par/${F_DESPIKE_PARMS}.$$" >> $F_SUMM_INPUT
mri_destroy_dataset data/junk
mri_destroy_dataset data/spike_count_all

