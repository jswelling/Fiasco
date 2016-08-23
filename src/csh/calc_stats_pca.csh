#!/bin/csh -ef
# calc_stats_pca.csh
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
echo '#'`date` $0
echo '# $Id: calc_stats_pca.csh,v 1.7 2005/02/09 21:24:23 welling Exp $'

#
# This script wants two inputs: the dataset to be analyzed (in image
# order), and the same dataset in permuted form (in time series order).
#

set data = $1
set tsdata = $2

#
# Do we have what we need?
#
if ( $#argv != 2 ) then
  echo '# ' $0 'Error: I need two parameters!'
  exit -1
endif
if ( ( ! -r $1 ) && ( ! -e ${1}.mri ) ) then
  echo '# ' $0 'Error: input data is missing!'
  exit -1
endif
if ( ( ! -r $2 ) && ( ! -e ${2}.mri ) ) then
  echo '# ' $0 'Error: permuted input data is missing!'
  exit -1
endif

#
# Create output dir if necessary
#
if ( ! -d stat ) mkdir stat

set xdim = `mri_printfield -field 'images.extent.x' $tsdata`
set ydim = `mri_printfield -field 'images.extent.y' $tsdata`
set zdim = `mri_printfield -field 'images.extent.z' $tsdata`
set tdim = `mri_printfield -field 'images.extent.t' $tsdata`

#
# Make a Grand Mean, for overlays
#
echo "##### Building GrandMean"
mri_subsample -d t -len 1 -mean $tsdata stat/GrandMean
mri_remap -order xyzt stat/GrandMean

set impute_count = 0
while ( $impute_count <= $F_PCA_NIMPUTATIONS )

  if ( $impute_count ) then

#
# This is an imputation:
# Make a copy of the input timeseries data with missing elements 
# substituted from values calculated from PCA results. tmp_pca
# was constructed in the first iteration.
#
    echo "##### Imputation $impute_count of $F_PCA_NIMPUTATIONS"
    mri_permute -order at tmp_pca_left tmp_pca_left_p
    mri_rpn_math -out tmp_pca_left_p_scaled '$1,$2,*' \
      tmp_pca_left_p tmp_pca_evals
    mri_destroy_dataset tmp_pca_left_p
    mri_matmult -out reconstructed tmp_pca_evecs tmp_pca_left_p_scaled
    mri_destroy_dataset tmp_pca_left_p_scaled
    mri_permute -order xyzt tmp_pca tmp_pca_old_p
    mri_rpn_math -out tmp_pca_p '$2,$3,missing,if_keep' \
      $data tmp_pca_old_p reconstructed
    mri_destroy_dataset reconstructed
    mri_remap -order xyzt tmp_pca_p
    mri_destroy_dataset tmp_pca_old_p
    mri_delete_chunk -chunk missing tmp_pca_p
    mri_permute -order txyz tmp_pca_p tmp_pca
    mri_destroy_dataset tmp_pca_p

  else

#
# First pass:
# Make a copy of the input timeseries data with missing elements 
# interpolated.
#
    echo "##### Smoothing missing data"
    mri_rpn_math -out tmp_pca_demean '$1,$2,-' $data stat/GrandMean
    mri_remap -order xyzt tmp_pca_demean
    mri_permute -order txyz tmp_pca_demean tmp_pca_demean_p
    smooth_missing -input tmp_pca_demean_p -headerout tmp_pca
    mri_delete_chunk -chunk missing tmp_pca
    mri_destroy_dataset tmp_pca_demean_p

  endif
#
# Remap the data and do pca
#

  @ qdim = $xdim * $ydim * $zdim
  mri_remap -order tq -len :$qdim tmp_pca
  echo "##### Performing PCA"
  pca.py -v --methodthresh $F_PCA_MTHDTHRESH \
    --leftvectors tmp_pca_left --eigenvalues tmp_pca_evals \
    --eigenvectors tmp_pca_evecs -n $F_PCA_NCOMPONENTS tmp_pca
  mri_remap -order txyz -len :${xdim}:${ydim}:${zdim} tmp_pca
  mri_remap -order xyza -len ${xdim}:${ydim}:${zdim}:${F_PCA_NCOMPONENTS} \
    tmp_pca_evecs
  mri_remap -order ta -len :${F_PCA_NCOMPONENTS} tmp_pca_left
  mri_remap -order a -len ${F_PCA_NCOMPONENTS} tmp_pca_evals

  @ impute_count = $impute_count + 1
end

#
# Disassemble the results.
#
echo "##### Assembling results"
set count = 0
while ( $count < $F_PCA_NCOMPONENTS )
  @ ind = $count + 1
  mri_subset -d a -len 1 -s $count tmp_pca_evecs stat/pca_$ind
  mri_remap -order xyz stat/pca_$ind
  mri_subset -d a -len 1 -s $count tmp_pca_left stat/pca_left_$ind
  mri_remap -order t stat/pca_left_$ind
  @ count = $count + 1
end

#
# Move eigenvalues into place
#
mri_copy_dataset tmp_pca_evals stat/pca_eigenvals
mri_remap -order a stat/pca_eigenvals

#
# I am horrified to say that I am about to do an ANOVA in csh.
#

echo "##### Performing termwise ANOVA"
@ dof = $qdim * $tdim
set dof_list = ( $dof )

mri_remap -order tq -len :$qdim tmp_pca
mri_permute -order qt tmp_pca tmp_pca_residual
mri_rpn_math -out tmp_pca_res_sqr '$1,dup,*' tmp_pca_residual
mri_subsample -d q -len 1 -sum tmp_pca_res_sqr tmp_pca_q
mri_subsample -d t -len 1 -sum tmp_pca_q tmp_pca_qt
set msto = `mri_rpn_math '$1,'$dof',/,1,if_print_1' tmp_pca_qt`
set mse_list = ( $msto )

set count = 1
while ( $count <= $F_PCA_NCOMPONENTS )
  mri_remap -order xyzt -len ${xdim}:${ydim}:${zdim}: tmp_pca_residual
  set eval = \
    `mri_rpn_math '$1,$a,'$count',1,-,==,if_print_1' stat/pca_eigenvals`
  mri_remap -order xyza stat/pca_$count
  mri_rpn_math -out tmp_pca_scaled '$1,'$eval',*' stat/pca_left_$count
  mri_remap -order at tmp_pca_scaled
  mri_matmult -out tmp_pca_term stat/pca_$count tmp_pca_scaled
  mri_remap -order xyz stat/pca_$count
  mri_remap -order t stat/pca_left_$count
  mri_rpn_math -out tmp_pca_new_residual '$1,$2,-' \
    tmp_pca_residual tmp_pca_term
  mri_copy_dataset tmp_pca_new_residual tmp_pca_residual
  mri_remap -order qt -len ${qdim}: tmp_pca_residual
  mri_rpn_math -out tmp_pca_res_sqr '$1,dup,*' tmp_pca_residual
  mri_subsample -d q -len 1 -sum tmp_pca_res_sqr tmp_pca_q
  mri_subsample -d t -len 1 -sum tmp_pca_q tmp_pca_qt
  @ dof = $dof - ( $tdim + $qdim )
  set mse = `mri_rpn_math '$1,'$dof',/,1,if_print_1' tmp_pca_qt`
  set dof_list = ( $dof_list $dof )
  set mse_list = ( $mse_list $mse )

  @ count = $count + 1
end

set count = 1
set f_list = ( )
set p_list = ( )
while ( $count <= $F_PCA_NCOMPONENTS )
  @ countplus1 = $count + 1
  set denom_mse = ${mse_list[$countplus1]}
  set num_mse = ${mse_list[$count]}
  set denom_dof = ${dof_list[$countplus1]}
  set num_dof = ${dof_list[$count]}
  set fval = \
    `mri_rpn_math $num_mse','$denom_mse',/,1,if_print_1' tmp_pca_qt`
  set f_list = ( $f_list $fval )
  set pval = \
    `mri_rpn_math '1,'$fval','$num_mse','$denom_mse',cf,-,1,if_print_1' tmp_pca_qt`
  set p_list = ( $p_list $pval )
  @ count = $count + 1
end

@ ncomp_plus_1 = $F_PCA_NCOMPONENTS + 1
echo $mse_list | \
  mri_from_ascii -ord a -len $ncomp_plus_1 stat/pca_mse_termwise
echo $f_list | \
  mri_from_ascii -ord a -len $F_PCA_NCOMPONENTS stat/pca_F_termwise
echo $p_list| \
  mri_from_ascii -ord a -len $F_PCA_NCOMPONENTS stat/pca_P_termwise

#
# Clean up
#
echo "##### Cleaning up"
foreach ds ( tmp_pca tmp_pca_demean tmp_pca_evecs tmp_pca_left \
             tmp_pca_evals tmp_pca_scaled tmp_pca_q tmp_pca_qt \
             tmp_pca_residual tmp_pca_new_residual tmp_pca_res_sqr \
             tmp_pca_term )
  if ( -e $ds.mri ) mri_destroy_dataset $ds
end

#
# create postscript files
echo "##### Making postscript"
pca_makeps.csh

