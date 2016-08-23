#! /bin/csh -exf
# spiral.reader.csh
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
echo '#$Id: spiral.reader.csh,v 1.20 2006/03/10 01:23:04 welling Exp $'

if (-e $F_DIR/$F_REFP.Z) then
  set ref_pfile = $F_REFP.Z
else if (-e $F_DIR/$F_REFP.z) then
  set ref_pfile = $F_REFP.z
else if (-e $F_DIR/$F_REFP.gz) then
  set ref_pfile = $F_REFP.gz
else
  set ref_pfile = $F_REFP
endif

# Find any other Pfiles in the data directory.  If there aren't any, we
# need to duplicate the name of the reference file to satisfy pfilorder.
set nonomatch
set PFILES = `cd $F_DIR;echo P*`
if ( "$PFILES" == 'P*') set PFILES = $ref_pfile
unset nonomatch

set ALLPFILES = `pfileorder $ref_pfile $PFILES`

set count = 0
foreach thispfile ( $ALLPFILES )

  if ( $thispfile =~ *.Z || $thispfile =~ *.z || $thispfile =~ *.gz ) then
    set thispfile = ${thispfile:r}
    set was_compressed = 1
    $F_UNCOMPRESS $thispfile
  else
    set was_compressed = 0
  endif

  smartreader -verbose -input $F_DIR/$thispfile \
    -out ${1}_tmp0 -tag $thispfile $F_READER_OPTS

  if ( $was_compressed ) $F_COMPRESS $thispfile
  set dimstr = `mri_printfield -field samples.dimensions ${1}_tmp0`
  if ( $dimstr == "vpstbzc" ) then
    set tdim = `mri_printfield -field samples.extent.t ${1}_tmp0`
    set bdim = `mri_printfield -field samples.extent.b ${1}_tmp0`
    @ newtdim = ${tdim} * ${bdim}
    mri_remap -chunk samples -order vpstzc \
	-length :::${newtdim}:: ${1}_tmp0
  endif

  set tdim = `mri_printfield -field samples.extent.t ${1}_tmp0`
  echo $count | mri_from_ascii -order t -length 1 -chunk acq_blocks acq
  mri_interp -d t -len $tdim -constant acq acq_stretch
  mri_destroy_dataset acq

  if ( -e ${1}_tmp.mri ) then
    mri_copy_chunk -chunk samples ${1}_tmp0 ${1}_tmp1
    mri_copy_chunk -chunk acq_blocks acq_stretch ${1}_tmp1
    mri_destroy_dataset ${1}_tmp0
    mri_destroy_dataset acq_stretch
    mri_paste -d t -out ${1}_tmp2 ${1}_samples ${1}_tmp1
    mri_destroy_dataset ${1}_tmp1
    foreach fname ( ${1}_tmp2.* )
      mv $fname ${1}_samples.${fname:e}
    end
  else
    mri_copy_dataset ${1}_tmp0 ${1}_tmp
    mri_delete_chunk -chunk samples ${1}_tmp 
    mri_copy_chunk -chunk samples ${1}_tmp0 ${1}_samples
    mri_copy_chunk -chunk acq_blocks acq_stretch ${1}_samples
    mri_destroy_dataset ${1}_tmp0
    mri_destroy_dataset acq_stretch
  endif

  @ count = $count + 1
end
mri_copy_chunk -chunk samples ${1}_samples ${1}_tmp
mri_copy_chunk -chunk acq_blocks -replace ${1}_samples ${1}_tmp
mri_destroy_dataset ${1}_samples
set tdim = `mri_printfield -field samples.extent.t ${1}_tmp`
set zdim = `mri_printfield -field samples.extent.z ${1}_tmp`

# Check that the user has specified a number of images and slices, etc,
# equal to the number actually found
if (! { dataset_matches_env.py -c samples -v ${1}_tmp } ) then
  echo "***ERROR*** environment variables conflict with file info"
endif

# Do relabeling and permutation if required
set dimstr = `mri_printfield -field samples.dimensions ${1}_tmp`
if ( ${dimstr} == "vpsztc" ) then
  set cdim = `mri_printfield -field samples.extent.c ${1}_tmp`
  if ( $cdim != 1 ) then
    echo "****ERROR**** multi-coil reordering is not implemented!"
    exit -1
  endif
  mri_remap -chunk samples -order vpsczt ${1}_tmp 
  set dimstr = vpsczt
endif
if ( ${dimstr} == "vpstzc" ) then
  set cdim = `mri_printfield -field samples.extent.c ${1}_tmp`
  if ( $cdim != 1 ) then
    echo "****ERROR**** multi-coil reordering is not implemented!"
    exit -1
  endif
  mri_remap -chunk samples -order vpsctz ${1}_tmp 
  set dimstr = vpsctz
endif
if ( ${dimstr} != "vpsczt" ) then
  mri_permute -chunk samples -order vpsczt ${1}_tmp ${1}_p
else
  mri_copy_dataset ${1}_tmp ${1}_p
endif
mri_destroy_dataset ${1}_tmp

# Do scan reordering if required.  Remapping acq_blocks prevents it
# from being folded.
set reorder = `mri_printfield -field samples.reorder -nofail ${1}_p`
if (${#reorder} != 0) then
  if (${reorder} != 0) then
    mri_scan_fold ${1}_p ${1}
    mri_destroy_dataset ${1}_p
  endif
else
  mri_copy_dataset ${1}_p ${1}
  mri_destroy_dataset ${1}_p
endif
