#!/bin/csh -ef
# overlay.csh
#/************************************************************
# *                                                          *
# *  Permission is hereby granted to any individual or       *
# *  institution for use, copying, or redistribution of      *
# *  this code and associated documentation, provided        *
# *  that such code and documentation are not sold for       *
# *  profit and the following copyright notice is retained   *
# *  in the code and documentation:                          *
# *     Copyright (c) 2003 Department of Statistics,         *
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
echo '#$Id: overlay.csh,v 1.2 2003/09/09 16:25:45 bakalj Exp $'
#

# Check for help command
if ( $#argv >= 1 ) then 
  if ( foo$argv[1] == "foo-help" ) then
    if ( $#argv >= 2 ) then
      scripthelp $0 $argv[2]
    else
      scripthelp $0
    endif
    exit
  endif
endif

set lowthresh = -6.0
set highthresh = 6.0
set outfile = "overlay.mri"
set inmap = "map.mri"
set inimage = "input.mri"
set mingray = 50
set maxgray = 204

while ( "${1}" =~ -* )
  if ( "${1}" == '-inmap' ) then
    set inmap = $2
    shift; shift
  else if ( "${1}" == '-lowthresh' ) then
    set lowthresh = $2
    shift; shift
  else if ( "${1}" == '-highthresh' ) then
    set highthresh = $2
    shift; shift
  else if ( "${1}" == '-inimage' ) then
    set inimage = $2
    shift; shift
  else if ( "${1}" == '-headerout' ) then
    set outfile = $2
    shift; shift
  else if ( "${1}" == '-mingray' ) then
    set mingray = $2
    shift; shift
  else if ( "${1}" == '-maxgray' ) then
    set maxgray = $2
    shift; shift
  else
    echo "${0} : Unrecognized argument ${1}"
    exit -1
  endif
end

# Make some scratch space, and go there.
if ( ${?F_TEMP} ) then
  set tmpdir = ${F_TEMP}/overlay_tmp_$$
else
  set tmpdir = ./overlay_tmp_$$
endif
if (! -e $tmpdir) mkdir $tmpdir
set homedir = `pwd`
if ( $inimage !~ /* ) then
  set inimage = $homedir/$inimage
endif
if ( $inmap !~ /* ) then
  set inmap = $homedir/$inmap
endif
if ( $outfile !~ /* ) then
  set outfile = $homedir/$outfile
endif
cd $tmpdir

mri_subsample -d x -l 1 -max $inimage max_x
mri_subsample -d x -l 1 -min $inimage min_x
mri_subsample -d y -l 1 -max max_x max_xy
mri_subsample -d y -l 1 -min min_x min_xy
mri_subsample -d z -l 1 -max max_xy max_xyz
mri_subsample -d z -l 1 -min min_xy min_xyz

mri_rpn_math -out scaled '$1,$3,-,$2,$3,-,/,192,*,32,+' $inimage max_xyz min_xyz

mri_from_ascii -c color_table -ord vc -length 5:4 cmap << CTBL_EOF
0 0 0 1 $lowthresh
0 0 0 0 $lowthresh
0 0 0 0 $highthresh
255 255 255 1 $highthresh
CTBL_EOF

colorize -col cmap $inmap over
matte -inmap over scaled $outfile

cd $homedir
rm -r $tmpdir
