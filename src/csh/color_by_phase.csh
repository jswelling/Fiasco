#! /bin/csh -ef
# color_by_phase.csh
#/************************************************************
# *                                                          *
# *  Permission is hereby granted to any individual or       *
# *  institution for use, copying, or redistribution of      *
# *  this code and associated documentation, provided        *
# *  that such code and documentation are not sold for       *
# *  profit and the following copyright notice is retained   *
# *  in the code and documentation:                          *
# *     Copyright (c) 1998,1999 Department of Statistics,    *
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
#echo '#'`date` $0
#echo '# $Id: color_by_phase.csh,v 1.11 2005/04/06 18:59:50 welling Exp $'
#
# This script generates a color image from a complex dataset, by
# using magnitude (as a fraction of min-max range) and phase as
# the hue and value of an HSV color triple, with saturation=1.
# The one argument is the name of a single complex image file.
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

# Parse command line options
set args = `getopt bdrgp $*`
if ($#args < 2 ) then
  echo "Unrecognized command line option in " $*
  scripthelp $0
  exit -1
endif
set ropt = ""
set bopt = ""
set dopt = ""
set pngopt = ""
while ($#args > 1) 
  switch ( $args[1] )
    case '-b' : 
      set bopt = "-brighten" ; shift args; breaksw;
    case '-d' : 
      set dopt = "-darken" ; shift args; breaksw;
    case '-r' : 
      set ropt = "-reverse" ; shift args; breaksw;
    case '-p' : 
      set pngopt = "-png" ; shift args; breaksw;
    case '--' : 
      shift args; breaksw;
  endsw
  if (junk${args[1]} == 'junk--') then
    shift args;
    break;
  endif
end
if ($#args != 1) then
  echo ${0:t} ' only one filename permitted'
  scripthelp $0
  exit -1
endif
set argv = $args

# Get the root name, to construct output filename
if ($1:t == ".mri") then
  set tailname = $1:t
  set fname = $tailname:r
else
  set fname = $1:t
endif

# Grab some ancillary info from the .mri file
set type    = \
        `mri_printfield -field images.datatype -nofail $1`
set idtag = \
        `mri_printfield -field images.scan.id -nofail $1`

# Use the file name for first title
setenv F_PS_TITLE1 $1

# Use the ID tag as second title if available.
if (! ${?F_PS_TITLE2} && dummy$idtag != dummy ) then
  setenv F_PS_TITLE2 "ID $idtag"
endif

set vdim = `mri_printfield -field images.extent.v -nofail $1`
set xdim = `mri_printfield -field images.extent.x $1`
set ydim = `mri_printfield -field images.extent.y $1`
set zdim = `mri_printfield -field images.extent.z -nofail $1`
if ( foo$zdim == foo ) set zdim = 1
set tdim = `mri_printfield -field images.extent.t -nofail $1`
if ( foo$tdim == foo ) set tdim = 1

if ( foo$vdim != foo2 ) then
  echo ${0:t} "Error: Input file $1 is not complex!"
  exit -1
endif

@ qdim = $xdim * $ydim * $zdim * $tdim

#
# Make some scratch space, and go there.
#
if ( ${?F_TEMP} ) then
  set tmpdir = ${F_TEMP}/color_by_phase_tmp_$$
else
  set tmpdir = ./color_by_phase_tmp_$$
endif
if (! -e $tmpdir) mkdir $tmpdir
set homedir = `pwd`
if ( $1 =~ /* ) then
  set input = $1
else
  set input = $homedir/$1
endif
cd $tmpdir

#Collect the undetermined options.  
set allopts = "-black_min 0.0 -black_max 255.0 -mosaic -title -date_created"
set allopts = "$allopts  -date_printed $pngopt $ropt $bopt $dopt"

mri_complex_to_scalar -mag $input mags
mri_remap -order vxyzt -length 1:${xdim}:${ydim}:${zdim}:${tdim} mags
mri_complex_to_scalar -pha $input phase
mri_remap -order vxyzt -length 1:${xdim}:${ydim}:${zdim}:${tdim} phase 

# If this is byte data, assume a range 0 to 255.  Otherwise 
# we rescale the range.
if ($type == "uint8") then
  mri_rpn_math -out scaled_mags '$1,255,/' mags 
else
  mri_copy_dataset mags tmp1
  mri_remap -order q -length $qdim tmp1
  mri_subsample -d q -l 1 -max tmp1 mag_max
#  Uncomment these if you really want to calculate a minimum magnitude
#  mri_subsample -d q -l 1 -min tmp1 mag_min
#  mri_rpn_math -out scaled_mags '$1,$2,-,$3,$2,-,/' mags mag_min mag_max
  mri_rpn_math -out scaled_mags '$1,$2,/' mags mag_max
endif

#
# The following implements the HSV->RGB conversion from "Computer
# Graphics Principles and Practice" by Foley, van Dam, Feiner,
# and Hughes, chapter 13.  H is the complex phase, S is 1.0, and
# V is the rescaled magnitude.
#
mri_rpn_math -out phase_zone '$1,0,2,pi,*,$1,0,>,if_keep,+,pi,/,3,*' phase
mri_rpn_math -out phase_int '$1,floor,5,min,0,max' phase_zone
mri_rpn_math -out phase_off '$1,$2,-' phase_zone phase_int

mri_rpn_math -out p '0' scaled_mags phase_off
mri_rpn_math -out q '$1,1,$2,-,*' scaled_mags phase_off
mri_rpn_math -out t '$1,$2,*' scaled_mags phase_off

mri_rpn_math -out r '$1,$4,$2,$2,$3,$1,$5,switch_6,255,*' \
    scaled_mags p q t phase_int
mri_rpn_math -out g '$2,$2,$3,$1,$1,$4,$5,switch_6,255,*' \
    scaled_mags p q t phase_int
mri_rpn_math -out b '$3,$1,$1,$4,$2,$2,$5,switch_6,255,*' \
    scaled_mags p q t phase_int
mri_paste -d v -out result r g b

if ( foo$pngopt != foo ) then
 mri_to_img $allopts result $homedir/$fname.png
else
 mri_to_img $allopts result $homedir/$fname.ps
endif

#clean up
cd $homedir
rm -r $tmpdir

