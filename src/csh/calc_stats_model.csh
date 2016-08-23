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
echo '# $Id: calc_stats_model.csh,v 1.2 2005/06/10 20:13:26 welling Exp $'

#
# This script wants two inputs: the dataset to be analyzed (in image
# order), and the same dataset in permuted form (in time series order).
#

set data = $1
set tsdata = $2

#
# Do we have what we need?
#
if ( ! $?F_MODEL_FORMULA ) then
  echo '# ' $0 'Error: F_MODEL_FORMULA is not set!'
  exit -1
endif
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
# Build temp directory
#
set tmpdir = ./model_tmp_$$
if ( ! -e $tmpdir ) mkdir $tmpdir
echo '# Temporary directory is' $tmpdir

#
# Build the raw design matrix
#
echo "# Building model matrix..."
build_model_matrix.py -v --nimages $tdim --nslices $zdim \
    --out $tmpdir/raw_design --model "$F_MODEL_FORMULA" \
    $F_SPLIT_FILE

#
# Extract design column information
#
set ncols = `mri_printfield -fld images.extent.v $tmpdir/raw_design`
set col = 0
set colNames = ()
while ( $col < $ncols )
  set name = `mri_printfield -fld images.label.v.$col $tmpdir/raw_design`
  set colNames = ( $colNames $name )
  @ col = $col + 1
end

#
# Build updated missing information.  We copy the input dataset because
# we have to mess with its missing information and we don't want to
# corrupt the original.
#
echo "# Building updated missing info"
mri_subsample -d v -len 1 -sum $tmpdir/raw_design $tmpdir/raw_summed
mri_rpn_math -out $tmpdir/design_missing_p '1,0,$1,is_finite,if_keep' \
  $tmpdir/raw_summed
mri_remap -order tz $tmpdir/design_missing_p
mri_permute -order zt $tmpdir/design_missing_p $tmpdir/design_missing
mri_copy_dataset $2 $tmpdir/permuted
mri_copy_chunk -chunk missing -chunk_out images \
  $tmpdir/permuted $tmpdir/old_missing
mri_rpn_math -out $tmpdir/all_missing '0,1,$1,$2,+,if_keep' \
  $tmpdir/design_missing $tmpdir/old_missing
mri_type_convert -char $tmpdir/all_missing $tmpdir/all_missing_char
mri_copy_chunk -chunk images -chunk_out missing -replace \
  $tmpdir/all_missing_char $tmpdir/permuted

#
# Do the regression.
#
echo "# Regressing..."
pushd $tmpdir >& /dev/null
mri_rpn_math -out filtered_design '0,$1,dup,is_finite,if_keep' raw_design
mri_glm -v -sum_of_squares -estimates glm_out:images permuted filtered_design
mri_destroy_dataset permuted

#
# Construct the results
#
echo "# Assembling results..."
mri_subset -d v -len 1 -s 0 glm_out GrandMean
set col = 0
@ offset_base = 1
while ( $col < $ncols )
  @ offset = $col + $offset_base
  mri_subset -d v -len 1 -s $offset glm_out BetaMap_$col
  @ col = $col + 1
end
@ offset_base = $offset_base + $ncols
mri_subset -d v -len 1 -s $offset_base glm_out counts
@ offset_base = $offset_base + 1
mri_subset -d v -len 1 -s $offset_base glm_out ssto
@ offset_base = $offset_base + 1
mri_subset -d v -len 1 -s $offset_base glm_out sse
@ offset_base = $offset_base + 1
set col = 0
while ( $col < $ncols )
  @ offset = $col + $offset_base
  mri_subset -d v -len 1 -s $offset glm_out ssr_$col
  @ col = $col + 1
end
mri_rpn_math -out mse '$1,$2,'$ncols',1,+,-,/' sse counts
mri_rpn_math -out msto '$1,$2,/' ssto counts
mri_rpn_math -out msr '$1,$2,-,'$ncols',1,+,/' ssto sse
mri_rpn_math -out Fmap_model '$1,$2,/' msr mse
mri_rpn_math -out Pmap_model \
  '$1,'$ncols',1,+,$2,'$ncols',1,+,-,fcf,-1,*,inv_foldp' \
  Fmap_model counts
set col = 0
while ( $col < $ncols )
  mri_rpn_math -out Fmap_$col '$1,$2,/' ssr_$col mse
  mri_rpn_math -out Pmap_$col '$1,1,$2,'$ncols',1,+,-,fcf,-1,*,inv_foldp' \
    Fmap_$col counts
  @ col = $col + 1
end

#
# Copy data out and Clean up
#
echo "# Moving results into place..."
popd >& /dev/null
foreach ds ( GrandMean Fmap_model Pmap_model counts msr mse msto )
  mri_copy_dataset $tmpdir/$ds stat/$ds
end
set col = 0
while ( $col < $ncols )
  @ colp1 = $col + 1
  mri_copy_dataset $tmpdir/BetaMap_$col stat/BetaMap_"${colNames[$colp1]}"
  mri_copy_dataset $tmpdir/Fmap_$col stat/Fmap_"${colNames[$colp1]}"
  mri_copy_dataset $tmpdir/Pmap_$col stat/Pmap_"${colNames[$colp1]}"
  @ col = $col + 1
end
rm -r $tmpdir

#
# create postscript files
echo "##### Making postscript"
makeps.csh

