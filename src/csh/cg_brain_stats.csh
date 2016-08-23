#!/bin/csh -efx
# cg_brain_stats.csh
#/************************************************************
# *                                                          *
# *  Permission is hereby granted to any individual or       *
# *  institution for use, copying, or redistribution of      *
# *  this code and associated documentation, provided        *
# *  that such code and documentation are not sold for       *
# *  profit and the following copyright notice is retained   *
# *  in the code and documentation:                          *
# *     Copyright (c) 1997                                   *
# *                        Department of Statistics,         *
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
# *  Original programming by Joel Welling                    *
# ************************************************************/
#
# This code takes the data file produced by vmpfx and extracts
# the 'interesting' parts for display.  
#
echo '#'`date` $0
echo '#$Id: cg_brain_stats.csh,v 1.16 2003/09/09 16:18:10 bakalj Exp $'

if ( ! -d $F_STAT_MPATH) mkdir $F_STAT_MPATH
if ( ! -d $F_STAT_SPATH) mkdir $F_STAT_SPATH
if ( ! -d $F_STAT_TPATH) mkdir $F_STAT_TPATH

set ofile = out/vmpfx_results.out
if (! -e $ofile) then
	echo 'Expected vmpfx output does not exist.'
	exit -1
endif

# Add up all the parts to figure out how many results we will generate
# Total is estimates plus standard errors for a whole bunch of parameters
set n_estimates = `awk '/estimates/ { print $2 ; exit }' < $ofile`
set n_errors = `awk '/standard.errors/ { print $2 ; exit }' < $ofile`
set n_baselevel = `awk '/baseline/ { print $2 ; exit }' < $ofile`
set n_drift = `awk '/drift\.coef/ { print $2 ; exit }' < $ofile`
set n_noise = `awk '/noise\.precision/ { print $2 ; exit }' < $ofile`
set n_shape = `awk '/shape/ { print $2 ; exit }' < $ofile`
set n_resp = `awk '/responsiveness/ { print $2 ; exit }' < $ofile`

# Consistency checks
if ( $n_estimates != $n_errors ) then
	echo 'Unexpected format for vmpfx output!'
	exit -1
endif
if ( $n_errors != $n_baselevel + $n_drift + $n_noise + $n_shape + $n_resp ) then
	echo 'Parsing error on vmpfx output!'
	exit -1
endif
if ( $n_estimates != ( $n_baselevel + $n_drift + \
	+ $n_noise + $n_shape + $n_resp )) then
	echo 'Unexpected format for vmpfx estimates output!'
	exit -1
endif

# Make a temporary file with just fit data
mri_copy_chunk -c images -cho images $1 brain_stat_tmp

# Build a map of probability from log of Bayes factor.  We convert the
# input to doubles so that the final output will be in doubles.
mri_copy_chunk -c lnbfac -cho images $1 brain_stat_lnb_tmp
mri_type_convert -double brain_stat_lnb_tmp brain_stat_lnb_dbl_tmp
mri_rpn_math '$1,exp,'${F_BAYES_P0}',*,dup,1.0,+,/' \
    -out $F_STAT_MPATH/Bayes.P brain_stat_lnb_dbl_tmp
mri_destroy_dataset brain_stat_lnb_tmp
mri_destroy_dataset brain_stat_lnb_dbl_tmp

# extract baselevel statistics
@ base_stde_offset = $n_estimates
mri_subset -d v -l 1 -s 0 brain_stat_tmp $F_STAT_MPATH/Baselevel.Mean
mri_subset -d v -l 1 -s $base_stde_offset brain_stat_tmp \
    $F_STAT_SPATH/Baselevel.Stdv

# extract drift1, drift2 statistics
#@ drift1_stde_offset = 1 + $n_estimates
#mri_subset -d v -l 1 -s 1 brain_stat_tmp $F_STAT_MPATH/Drift1.Mean
#mri_subset -d v -l 1 -s $drift1_stde_offset brain_stat_tmp \
#    $F_STAT_SPATH/Drift1.Stdv 
#@ drift2_stde_offset = 2 + $n_estimates
#mri_subset -d v -l 1 -s 2 brain_stat_tmp $F_STAT_MPATH/Drift2.Mean
#mri_subset -d v -l 1 -s $drift2_stde_offset brain_stat_tmp \
#    $F_STAT_SPATH/Drift2.Stdv 

# extract noise and shape statistics.  mri_subset counts from 0, not 1
@ noise_offset = $n_baselevel + $n_drift
@ shape_offset = $n_baselevel + $n_drift + $n_noise 
@ shape_stde_offset = $shape_offset + $n_estimates
mri_subset -d v -s $noise_offset -l $n_noise \
	brain_stat_tmp $F_VMPFX_NOISE
mri_subset -d v -s $shape_offset -l $n_shape \
	brain_stat_tmp $F_VMPFX_SHAPE
mri_subset -d v -s $shape_stde_offset -l $n_shape \
	brain_stat_tmp $F_VMPFX_SHPSE

# Now we want to extract stats for responses, but only for those
# conditions that are not fixed.  This is done using special 
# translation info saved earlier by generate_vmpfx_in.

# Set up a list of unfixed Fiasco indices
set unfixed_list = foo
shift unfixed_list

# need the map between brain and Fiasco conds, generated 
# by generate_vmpfx_in.
set cond_map = `grep FiascoCondMap $F_VMPFX_INPUT`
shift cond_map

# Walk the condition map, generating stats for non-fixed
# conditions.
foreach mapentry ($cond_map)
  if (${mapentry:h} == v) then
    set pair = ${mapentry:t}
    set bindex = ${pair:r}
    set findex = ${pair:e}
    set unfixed_list = ( $unfixed_list $findex )
    @ cond_offset = $n_baselevel + $n_drift + $n_shape + $n_noise + $bindex
    @ cond_stde_offset = $cond_offset + $n_estimates
    mri_subset -d v -l 1 -s $cond_offset brain_stat_tmp \
	$F_STAT_MPATH/Resp.$findex.Mean
    mri_subset -d v -l 1 -s $cond_stde_offset brain_stat_tmp \
	$F_STAT_SPATH/Resp.$findex.Stdv 
    mri_rpn_math '$1,$2,*' -out $F_STAT_TPATH/Bayes.PResp.$findex \
	    $F_STAT_MPATH/Resp.$findex.Mean $F_STAT_SPATH/Bayes.P 
  endif
end

# Generate fake Ttest data, if there will not be too much.  The cutoff
# at 10 unfixed conditions is arbitrary, but the number of Ttests
# produced is half the square of this number.
if ( $F_FAKE_TMAPS ) then
  if ($#unfixed_list <= 10) then
    while ( $#unfixed_list > 1 )
      set base_cond = $unfixed_list[1]
      shift unfixed_list
      foreach high_cond ($unfixed_list)
        mri_rpn_math '$2,$3,-,$1,*,$4,dup,*,$5,dup,*,+,sqrt,/' \
          -out $F_STAT_MPATH/Resp.${base_cond}-${high_cond}.Tmap \
          $F_STAT_MPATH/Bayes.P $F_STAT_MPATH/Resp.$base_cond.Mean \
          $F_STAT_MPATH/Resp.$high_cond.Mean \
          $F_STAT_MPATH/Resp.$base_cond.Stdv \
          $F_STAT_MPATH/Resp.$high_cond.Stdv
      end
    end
  else
    echo '**** Too many conditions to generate tmaps! ****'
  endif
endif

# Clean up the temp file
mri_destroy_dataset brain_stat_tmp




