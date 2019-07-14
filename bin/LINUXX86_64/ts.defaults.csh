#!/bin/csh -efx
# ts.defaults.csh
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
echo '#$Id: ts.defaults.csh,v 1.46 2005/04/06 19:00:52 welling Exp $'
if(! -d in)mkdir in
#
#
# if the following variables have not already been specified, set them here
if (! $?F_TEMP)		setenv F_TEMP		/tmp
if (! $?F_DIR)		setenv F_DIR		$FIASCO/../../tsdata/
if (! $?F_PARALLEL)	setenv F_PARALLEL	0
if (! $?F_PARALLEL_HOSTS) setenv F_PARALLEL_HOSTS ""
#
# disable parallelism unless later turned on by F_PARALLEL flag in
#	parallel.start.csh
unsetenv PAR_ENABLE
#
# ts.reader.csh
#
# input file header
setenv F_HEADER		tssample
#
# input file name
setenv F_READER_INPUT   $F_DIR/$F_HEADER.ra
#
# number of slices
setenv F_NSLICE		4
#
# number of image pairs
setenv F_NIMAGE		11
#
# dimensions
# reader sees twice this many images
@ READER_IMAGES = 2 * $F_NIMAGE
setenv F_READER_DIMS	"256 128 $F_NSLICE $READER_IMAGES"
#
# header file
setenv F_READER_OUTPUT	reader
#
# beginning offset
setenv F_READER_OFFSET	0
#
# TS reording of slices?
setenv F_READER_REORDER 1
#
# interimage skip
setenv F_READER_SKIP	0
#
# interslice skip
setenv F_READER_SSKIP   0
#
# input file byte order
setenv F_READER_ENDIAN -bigendian
# 
# data order
setenv F_READER_DIMORDER vxyzt
#
# input data space (k for k-space, i for image space)
setenv F_READER_SPACE   k
#
# data type
#setenv F_READER_TYPE	MRI_INT16
setenv F_READER_TYPE	short
#
# vector length (2 for ts k-space)
setenv F_READER_VECLEN	2
#
# flag to indicate coronal imaging
setenv F_READER_CORONAL 0
#
# Other reader options
setenv F_READER_OPTS ""
#
#
# epi.despike.csh
#
setenv F_DESPIKE_OUTPUT despike
setenv F_DESPIKE_PARMS  spike_count.par
setenv F_DESPIKE_CUTOFF 4.0
#
#
# ts.baseline.csh
#
# baseline header file
setenv F_BASELINE_OUTPUT baseline
#
# parameter estimates
setenv F_BASELINE_RWPARMS basadj_raw.par
setenv F_BASELINE_PARMS basadj.par
#
# baseline reverse
setenv F_BASELINE_REVERSE odd
#
# smoother options (usually none)
setenv F_BASELINE_SMOPTS 
#
#
# partialk.csh
#
# baseline header file
setenv F_PARTIALK_OUTPUT partialk
#
# scan line reverse 
setenv F_PARTIALK_REVERSE none
#
#
# ts.deghost.csh
#
# deghost header file
setenv F_DEGHOST_OUTPUT deghost
#
# temporary file used by deghost script
setenv F_DEGHOST_TEMP deghost_tmp
#
# deghost unsmoothed estimates
setenv F_DEGHOST_RWPARMS deghost_raw.par
#
# deghost smoothed estimates
setenv F_DEGHOST_PARMS	deghost.par
#
# deghost phase
setenv F_DEGHOST_PHASE	byslice
# 
# deghost filter bandwidth
setenv F_DEGHOST_BAND 8.0
#
# reverse what lines?
setenv F_DEGHOST_REVERSE none
#
# unsmoothed estimates if navigator is present
setenv F_NAVDGH_RWPARMS  navdeghost_raw.par
#
# smoothed estimates if navigator is present
setenv F_NAVDGH_PARMS  navdeghost.par
#
# filter bandwidth if navigator is present
setenv F_NAVDGH_BAND 12.0
#
# reverse what lines in images if navigator is present?
setenv F_NAVDGH_REVERSE odd
#
# reverse what lines in navigator if navigator is present?
setenv F_NAVDGH_NAVREVERSE even
#
#
# ts.baseline2.csh
#
# merge, with second baseline correction
setenv F_BASELINE2_OUTPUT baseline2
#
# parameter estimates
setenv F_BASELINE2_PARMS basadj2.par
#
# baseline2 reverse
setenv F_BASELN2_REVERSE none
#
# baseline2 overlap (actually shift, equal to half overlap)
setenv F_BASELN2_OVER 10
#
#
# ts.deghost2.csh (deghost of merged images)
#
# deghost header file
setenv F_DEGHOST2_OUTPUT deghost2
#
# deghost unsmoothed estimates
setenv F_DEGHOST2_RWPARMS deghost2_raw.par
#
# deghost smoothed estimates
setenv F_DEGHOST2_PARMS	deghost2.par
#
# deghost phase
setenv F_DEGHOST2_PHASE	byslice
# 
# deghost filter bandwidth
setenv F_DEGHOST2_BAND 8.0
#
# phase_lock.csh
#
# phase_lock header file
setenv F_PHLOCK_OUTPUT phase_lock
# unsmoothed estimates
setenv F_PHLOCK_RWPARMS phase_lock_raw.par
# smoothed estimates
setenv F_PHLOCK_PARMS phase_lock.par
# smoother bandwidth
setenv F_PHLOCK_BAND 3.0
# algorithm: "weighted" or "center"
setenv F_PHLOCK_ALG "weighted"
#
# phase_undrift.csh
#
# phase_undrift header file
setenv F_PHUND_OUTPUT phase_undrift
# unsmoothed estimates
setenv F_PHUND_RWPARMS phase_undrift_raw.par
# smoothed estimates
setenv F_PHUND_PARMS phase_undrift.par
# smoother bandwidth
setenv F_PHUND_BAND 4.0
#
# ts.meanc.csh
#
# meanc header file
setenv F_MEANC_OUTPUT	meanc
#
# mean correction estimates
setenv F_MEANC_PARMS	meanc.par
#
# fixed image number
@ F_N2 = $F_NIMAGE / 2
setenv F_MEANC_FIXED	$F_N2
#
# epi.meanc3d.csh
#
# meanc header file
setenv F_MEANC3D_OUTPUT	${F_MEANC_OUTPUT}
#
# mean correction estimates
setenv F_MEANC3D_PARMS	meanc3d.par
#
#
#
# recon1.csh
#
#recon1 header file
setenv F_RECON1_OUTPUT	recon1
#
#inverse transform
setenv F_RECON1_DIRECTION inverse
#
#output type
setenv F_RECON1_MODULUS modulus
#
#
# ts.clip1.csh
#
# clip1 header file
setenv F_CLIP1_OUTPUT	clip1
#
# width of output image
setenv F_CLIP1_WIDTH	256
#
# height of output image
setenv F_CLIP1_HEIGHT	128
#
# location of center in x direction
setenv F_CLIP1_XCENTER	center
#
# location of center in y direction
setenv F_CLIP1_YCENTER	64
#
#
# estireg.csh
#
# registration estimates
setenv F_ESTIREG_PARMS	estireg.par
#
# base image
setenv F_ESTIREG_FIXED	$F_N2
#
# hosts to use if parallel
setenv F_ESTIREG_HOSTS	"$F_PARALLEL_HOSTS"
#
#
# estireg3d.csh
#
# registration estimates
setenv F_EST3D_PARMS  reg3d.par
#
# algorithm specification
setenv F_EST3D_ALG 'opt=praxis'
#
# alignment target
setenv F_EST3D_ALIGN "mean"
setenv F_EST3D_STDVALN "mean"
#
#
# parsm.csh
#
# -parameterout
setenv F_PARSM_PARMS	smreg.par
#
# -cutoff
setenv F_PARSM_CUTOFF	"10.0"
#
# -bandwidth
setenv F_PARSM_BANDWIDTH "3.0"
#
# -kernel
setenv F_PARSM_KERNEL	"g"
#
# -threshold
setenv F_PARSM_TRANST	"0.1"
#
#
# smoothing parameters for smregpar3d
#
# -parameterout
setenv F_PARSM3D_PARMS	smreg3d.par
# -cutoff
setenv F_PARSM3D_CUTOFF ${F_PARSM_CUTOFF}
# -bandwidth
setenv F_PARSM3D_BANDWIDTH ${F_PARSM_BANDWIDTH}
# -kernel
setenv F_PARSM3D_KERNEL ${F_PARSM_KERNEL}
# -threshold
setenv F_PARSM3D_TRANST ${F_PARSM_TRANST}
#
#
# ts.ireg.csh
# -headerout
setenv F_IREG_OUTPUT	ireg
#
#
# physio_correct.csh
#
# Input physiology dataset (defaults to empty)
setenv F_PHYS_DATA ""
#
# Corrected output data
setenv F_PHYS_OUTPUT phys_corrected
#
# subsampling method
setenv F_PHYS_SUBOP max
#
# Fraction of scan time the ADCs are on
setenv F_PHYS_WINDOW 0.1587
#
# respiration per-slice data
setenv F_PHYS_RESP in/respiration
#
# cardiac per-slice data
setenv F_PHYS_CARDIO in/cardiac
#
# regression parameter data
setenv F_PHYS_PARAM stat/physio
#
#
# ts.recon2.csh
#
#recon2 header file
setenv F_RECON2_OUTPUT	recon2
#
#inverse transform
setenv F_RECON2_DIRECTION inverse
#
#output type
setenv F_RECON2_MODULUS modulus
#
#
# ts.clip2.csh
#
# clip2 header file
setenv F_CLIP2_OUTPUT	clip2
#
# width of output image
setenv F_CLIP2_WIDTH	256
#
# height of output image
setenv F_CLIP2_HEIGHT	128
#
# location of center in x direction
setenv F_CLIP2_XCENTER	center
#
# location of center in y direction
setenv F_CLIP2_YCENTER	64
#
#
# outlier.csh
#
# -headerout
setenv F_OUTLIER_OUTPUT outlier
#
# standard deviation threshold
setenv F_OUTLIER_STDVS	3.0
#
# adjustment parameters 
setenv F_OUTLIER_PARMS	outlier.par
#
#
#detrend.csh
#
# -headerout
setenv F_DETREND_OUTPUT detrend
#
# -chunk
setenv F_DETREND_CHUNK	images
#
# -headerout
setenv F_DETREND_TEMP	permute
#
# -parhead
setenv F_DETREND_PARHEAD detpar
#
#
#spline_detrend.csh
#
setenv F_VMPFX_NPROTO    ${FIASCO}/../../src/brain/vmpfx_null_proto.t
#
setenv F_VMPFX_NINPUT    in/vmpfx_null_in.t
# image acquisition interval in seconds
setenv F_IAI 3.0
#
#
#smooth_detrend.csh
#
#  denominator used to compute bandwidth (must be an integer!)
setenv F_SMDETREND_DNM 8
#
#
#split.csh
#
# split file
setenv F_SPLIT_FILE	in/split
#
# split file data
#
cat > $F_SPLIT_FILE << %%%END_SPLIT_INPUT%%%
1 bynumimages
Condition
Rest 4
Test 7
%%%END_SPLIT_INPUT%%%

#
# condition file
setenv F_SPLIT_COND	in/conditions
#
#
# maximum number of conditions per factor
setenv F_SPLIT_MAXCOND	16
#
# output file
setenv F_SPLIT_NEW	in/newsplit
#
#
# stat.csh
#
# mean prefix
setenv F_STAT_MPATH	stat/
#
# mean prefix
setenv F_STAT_SPATH	stat/
#
# mean prefix
setenv F_STAT_TPATH	stat/
#
# max ttests
setenv F_STAT_MAXTPR    50
#
# fmaps.csh (requires no parameters)
#
#
# displace.csh
#
# displace output
setenv F_DISPLACE_OUTPUT displace.par
# x dimension
setenv F_DISPLACE_XDIM	256
# y dimension
setenv F_DISPLACE_YDIM	128
# displace smoutput
setenv F_DISPLACE_SMOUT smdisplace.par
# weight for displacement calc.  Value "" means no weighting.  Value
# can be a filename (without the .mri) or "mean" or "median".
setenv F_DISPLACE_WEIGHT "mean"
#
#
# displace3d.csh
#
# displacement, etc. for raw motion estimates
setenv F_DISP3D_RAW_PARMS displace3d_raw.par
# displacement, etc. smoothed motion estimates
setenv F_DISP3D_SMTH_PARMS displace3d_smooth.par
# See defaults for displace.csh for weighting info
#
#
# summary.csh
#
# input list
setenv F_SUMM_INPUT	par/summary.input
#
# history of missing slices
setenv F_SUMM_MISSING   par/summary.missing
#
# output file
setenv F_SUMM_OUTPUT	summary.ps
#
# creation date
setenv F_CREDATE	"`date`"
#
# short description
setenv F_DESCRIPTION	"TS TEST DATA SET"


setenv F_PS_TITLE1	"$F_DESCRIPTION"
setenv F_PS_TITLE2	"$F_NSLICE Slices"
# upper t threshold
setenv F_UPPERT		6
#lower t threshold
setenv F_LOWERT		-6
# upper F threshold
setenv F_UPPERF		36
# image offset
setenv F_OFFSET		0
#
#printpar.csh
#
# postscript output destination
setenv F_PRINTPAR_SUBDIR ps
# maximum number of plots vertically
setenv F_PRINTPAR_MAXPLOT 4
#
#
# print.files.csh
setenv F_PRINTER	none
#
#
#
# calc_stats_brain.csh
#
# -fixed
setenv F_VMPFX_FIXED NA
#
# Base probability
setenv F_BAYES_P0 0.0025
#
# vmpfx input prototype file
setenv F_VMPFX_PROTO    ${FIASCO}/../../src/brain/vmpfx_proto.t
# 
# vmpfx input file name
setenv F_VMPFX_INPUT    in/vmpfx_in.t
#
# vmpfx output file name
setenv F_VMPFX_OUTPUT vmpfx_out
#
# vmpfx output for BOLD shape, BOLD shape stde, and noise precision
setenv F_VMPFX_SHAPE stat/bold_shape
setenv F_VMPFX_SHPSE stat/bold_shape_stde
setenv F_VMPFX_NOISE stat/vmpfx_noise_prec
#
# permuted vmpfx input data
setenv F_VMPFX_TEMP permute
#
# switch to generate fake tmaps
setenv F_FAKE_TMAPS 0


#
# Coregistration
#
# Pgh MRI dataset to be used as the prototype functional data for coreg
setenv F_CRG_FUNC stat/GrandMean
#
# name of the coregistered version of F_CRG_FUNC dataset
setenv F_CRG_FUNC_COREG aligned_GrandMean
#
# Input for smartreader to generate inplane anatomical data (unset by default)
setenv F_CRG_INPL_INPUT ""
#
# smartreader options for the inplane anatomical data
setenv F_CRG_INPL_OPTS ""
#
# name of the inplane anatomical Pgh MRI dataset
setenv F_CRG_INPL inplane
#
# name of skull-stripped inplane anatomical Pgh MRI dataset
setenv F_CRG_INPL_STRIPPED stripped_inplane
#
# Input for smartreader to generate axial anatomical data (unset by default)
setenv F_CRG_ANAT_INPUT ""
#
# smartreader options for the axial anatomical data
setenv F_CRG_ANAT_OPTS ""
#
# name of the axial anatomical Pgh MRI dataset
setenv F_CRG_ANAT axial
#
# name of skull-stripped axial anatomical Pgh MRI dataset
setenv F_CRG_ANAT_STRIPPED stripped_axial
#
# name of coregistered and warped axial anatomical Pgh MRI dataset
setenv F_CRG_ANAT_COWARPED warped_axial
#
# name of coregistered, warped axial resampled to match inplanes
setenv F_CRG_ANAT_INPL warped_axial_resampled
#
# Set the following to override the default optimization algorithm
# in coregister_inplane.py, coregister_struct_to_inplane.py, and
# cowarp_inplane.py.  The defaults actually take into account
# the acquisition types of the functional and structural scans.
# setenv F_CRG_INPL_ALG whatever
# setenv F_CRG_STRCT_ALG whatever
# setenv F_CRG_WARP_ALG whatever
#
#
# calc_stats_pca.csh:
# Stats via Principal component analysis
#
# Number of components to extract
setenv F_PCA_NCOMPONENTS 10
#
# Number of iterations of PCA to impute missing values
setenv F_PCA_NIMPUTATIONS 1
#
# First two dimensions of input matrix must be smaller than this
# to use SVD rather than an eigensystems method
setenv F_PCA_MTHDTHRESH 256
#

