#!/bin/csh -efx
# s.local.csh
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
#$Id: spiral.local.csh,v 1.19 2004/03/27 01:41:25 welling Exp $
#

#
# Subject information
#
setenv F_SUBJ_AGE
setenv F_SUBJ_SEX
setenv F_SUBJ_1DIAG none
setenv F_SUBJ_2DIAG none
setenv F_SUBJ_3DIAG none

#
# Input file location and names
#
# source directory
setenv F_DIR $FIASCO/../../spiraldata/
# ref file
setenv F_REFP P09216.7
# dataset name
setenv F_HEADER	sample_splx5

#
# Description
#
setenv F_DESCRIPTION "INTERESTING SPIRAL TEST DATA"

#
# Information describing the data format.
#
# number of slices
setenv F_NSLICE 2
# number of images
setenv F_NIMAGE 20
# lengths of dimensions for reconstruction
setenv F_SPIRAL_RESOLUTION 64

# Other reader options.  See the documentation for "smartreader";
# some common flags include '-multi' and '-auxfile'
#
setenv F_READER_OPTS ""

# If an axial anatomical dataset is available, complete
# the following lines.  See the documentation for "smartreader" for
# possible entries in the OPTS line.
#
setenv F_CRG_ANAT_INPUT ""
setenv F_CRG_ANAT_OPTS ""

# If an inplane anatomical dataset is available, complete
# the following lines.  See the documentation for "smartreader" for
# possible entries in the OPTS line.
#
setenv F_CRG_INPL_INPUT ""
setenv F_CRG_INPL_OPTS ""

# Voxel dimensions for the final reconstructed functional images.
# F_ZVOXEL should include the slice separation!
setenv F_XVOXEL 3.125
setenv F_YVOXEL 3.125
setenv F_ZVOXEL 4.2

# image acquisition interval in seconds
setenv F_IAI 1.0

# physiological data (defaults to nothing)
setenv F_PHYS_DATA ""
# fraction of scan time from slice selection to when the ADCs turn
# off.  For EPI, this is typically about 32 ms * F_NSLICE / TR in ms.
setenv F_PHYS_WINDOW 0.1587

# Fixed image numbers for mean correction and motion estimation.
# F_ESTIREG_FIXED can also be set to "mean" or "median", in which 
# case alignment will be done to the mean or median image rather 
# than any specific image. Don't set F_MEANC_FIXED to "mean", though!
@ F_N2 = $F_NIMAGE / 2
setenv F_MEANC_FIXED	$F_N2
setenv F_ESTIREG_FIXED	$F_N2

# Ttest and Ftest thresholds.  For data from a 1.5T magnet F_UPPERT
# and F_LOWERT should typically be 4.5 and -4.5 respectively; for
# data from a 3T magnet they should typically be 6.0 and -6.0
# respectively.  F_UPPERF should typically be the square of F_UPPERT.
setenv F_UPPERT		6
setenv F_LOWERT		-6
setenv F_UPPERF		36

# Threshold limits for P score overlays.  Uncommenting these lines will
# force Pmaps to be generated and will produce overlays of those maps
# with the given thresholds.
#setenv F_UPPERP       0.9999
#setenv F_LOWERP       0.0001

# False discovery threshold.  Uncommenting this line will force Pmaps
# to be generated, and will calculate and produce overlays with the
# given false discovery rate.
setenv F_LOWERQ 0.05

# Active voxel count threshold.  Uncommenting the first line will cause 
# Tmap overlays to be generated showing the N voxels with strongest
# activation.  Uncommenting the second line will pick these voxels 
# from both tails of the distribution, rather than just the left tail
#setenv F_PICKN_COUNT 200
#setenv F_PICKN_TWOTAILED 1

# Uncomment the following line to prevent summary information from
# this run being sent to Fiasco Headquarters.  This summary info is
# just the information in the output summary file, plus the state of
# Fiasco environment variables.
#
# setenv F_NO_MAIL_SUMMARY

#
#
# split file data
#
cat > $F_SPLIT_FILE << %%%END_SPLIT_INPUT%%%
1 bynumimages
Condition
Rest 4
Lights 4
Rest 4
Lights 4
Rest 4
%%%END_SPLIT_INPUT%%%

# Set titling information
setenv F_PS_TITLE1 "$F_DESCRIPTION"
setenv F_PS_TITLE2 "$F_NSLICE Slices"

