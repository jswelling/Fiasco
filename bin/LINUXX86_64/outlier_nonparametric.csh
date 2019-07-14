#!/bin/csh -efx
# outlier.csh
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
echo '# $Id: outlier_nonparametric.csh,v 1.9 2005/12/08 18:42:30 fiasco Exp $'
if(! -d par)mkdir par
mri_rpn_math -out ${1}_inf '$1,inf,missing,if_keep' $1
mri_remap -order vxyzt ${1}_inf 
mri_permute -memlimit 50000000 -order vtxyz ${1}_inf ${1}_p
mri_destroy_dataset ${1}_inf
mri_sort -asc ${1}_p ${1}_s
mri_subsample -d t -l 1 -count ${1}_p ${1}_c
mri_destroy_dataset ${1}_p
mri_subset -d x -l 1 -s 0 ${1}_c ${1}_cx
mri_destroy_dataset ${1}_c
mri_subset -d y -l 1 -s 0 ${1}_cx ${1}_counts
mri_destroy_dataset ${1}_cx

echo 'Counts per slice used for computing median and IQR:'
mri_rpn_math -out outlier_junk '0,$z,$1,1,if_print_2' ${1}_counts

set zdim = `mri_printfield -field images.extent.z $1`
@ z = 0
while ( $z < $zdim )
    set c = `mri_rpn_math -out outlier_junk '0,$1,$z,'$z',==,if_print_1' ${1}_counts` 
    @ median_offset = $c / 2
    @ q3_offset = $c / 4
    @ q1_offset = ( 3 * $c ) / 4
    mri_subset -d z -l 1 -s $z ${1}_s ${1}_s_z
    mri_subset -d t -l 1 -s $median_offset ${1}_s_z ${1}_median_z
    mri_subset -d t -l 1 -s $q1_offset ${1}_s_z ${1}_q1_z
    mri_subset -d t -l 1 -s $q3_offset ${1}_s_z ${1}_q3_z
    mri_rpn_math -out ${1}_iqr_z '$1,$2,-' ${1}_q1_z ${1}_q3_z
    mri_destroy_dataset ${1}_q1_z
    mri_destroy_dataset ${1}_q3_z

    if ( $z == 0 ) then
      mri_copy_dataset ${1}_median_z ${1}_median
      mri_copy_dataset ${1}_iqr_z ${1}_iqr
    else
      mri_paste -d z -out ${1}_median_tmp ${1}_median ${1}_median_z
      mri_copy_dataset ${1}_median_tmp ${1}_median
      mri_destroy_dataset ${1}_median_tmp
      mri_paste -d z -out ${1}_iqr_tmp ${1}_iqr ${1}_iqr_z
      mri_copy_dataset ${1}_iqr_tmp ${1}_iqr
      mri_destroy_dataset ${1}_iqr_tmp
    endif

    @ z = $z + 1
end

mri_remap -order vxyzt ${1}_median 
mri_remap -order vxyzt ${1}_iqr 

mri_destroy_dataset ${1}_s
mri_destroy_dataset ${1}_s_z
mri_destroy_dataset ${1}_counts
mri_destroy_dataset ${1}_median_z
mri_destroy_dataset ${1}_iqr_z
mri_destroy_dataset outlier_junk

outlier -input $1 -headerout $2 \
          -parameters par/$F_OUTLIER_PARMS.$$ -cutoff $F_OUTLIER_STDVS \
	  -median ${1}_median -iqr ${1}_iqr

set step = `depath.csh $0 | sed 's/_/./g'`
echo "$step par/$F_OUTLIER_PARMS.$$" >> $F_SUMM_INPUT
echo $step.$$ >> $F_SUMM_MISSING
count_missing.csh $2.mri >> $F_SUMM_MISSING

