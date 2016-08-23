#!/bin/csh -efx
# fiasco.local.csh
#/************************************************************
# *                                                          *
# *  Permission is hereby granted to any individual or       *
# *  institution for use, copying, or redistribution of      *
# *  this code and associated documentation, provided        *
# *  that such code and documentation are not sold for       *
# *  profit and the following copyright notice is retained   *
# *  in the code and documentation:                          *
# *     Copyright (c) 1997 Pittsburgh Supercomputing Center  *
# *                                                          *
# *  Original programming by Greg Hood (PSC)                 *
# ************************************************************/
#
# $Id: fiasco.local_orig.csh,v 1.4 2003/04/10 19:33:44 welling Exp $
#
# Environment variables specific to the local machine
#	Note that these may also be specified by the user by typing
#	"setenv" commands at the shell prompt prior to running FIASCO
# 
# temporary disk storage
setenv F_TEMP		/tmp/`whoami`
setenv MRI_TMP_DIR      $F_TEMP
#

# to enable parallelism, uncomment the next line
#setenv F_PARALLEL	1

# hosts to use if parallelism is enabled (empty string spawns one worker
#	on the local host)
#setenv F_PARALLEL_HOSTS	""

# where to find Splus
#setenv SPLUS            /usr/statlocal/bin/Splus
#setenv SPLUS /usr/statlocal/bin/Splus
setenv SPLUS "/usr/bin/R --no-save"

# what to use to compress and uncompress files
setenv F_COMPRESS       gzip
setenv F_UNCOMPRESS     gunzip

# where to store FFTW "wisdom".  Comment this out to avoid saving
# wisdom.
setenv F_WISDOM_FILE in/fftw_wisdom

# Some Fiasco programs will use the following hint to choose the
# size of internal memory buffers.  
setenv F_MEMSIZE_HINT 200000000

# Script to use in calculating experimental statistics (regions of
# activation, etc) from the reconstructed images.
#
setenv F_STATS_SCRIPT calc_stats_classic.csh

# Script to use in converting Pgh MRI files to Postscript.
# mri_to_ps.csh produces separate pages for each slice;
# allslice.csh puts up to 25 slices on one page.
#
#setenv F_MRI_TO_PS mri_to_ps.csh
setenv F_MRI_TO_PS allslice.csh

