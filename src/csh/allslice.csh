#! /bin/csh -ef
# allslice.csh
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
# *  Original programming by Bill Eddy                       *
# ************************************************************/
#
# This script exists for backward compatability;  it has become
# a pretty simple call to mri_to_img.
#
# The one argument is the name of a single image file
#

# Parse command line options
set args = `getopt bdrgn:x:v $*`
if ($#args < 2 ) then
  echo "usage: " $0 '[-b -d -r -g] filename'
  echo "  -b brightens image"
  echo "  -d dims image"
  echo "  -r reverses black and white"
  echo "  -g turns off gray scale bar"
  echo "  -n min sets minimum value"
  echo "  -v gives verbose output"
  echo "  -x max sets maximum value"
  exit -1
endif
set ropt = "-reverse"
set gopt = "-gray_scale"
set bopt = ""
set dopt = ""
set minopt = ""
set maxopt = ""
set verbose_flg = 0
while ($#args > 1) 
  switch ( $args[1] )
    case '-b' : 
      set bopt = "-brighten" ; shift args; breaksw;
    case '-d' : 
      set dopt = "-darken" ; shift args; breaksw;
    case '-r' : 
      set ropt = "" ; shift args; breaksw;
    case '-g' : 
      set gopt = "" ; shift args; breaksw;
    case '-n' :
      set minopt = "-black_min $args[2]" ; shift args; shift args; breaksw;
    case '-x' :
      set maxopt = "-black_max $args[2]" ; shift args; shift args; breaksw;
    case '-v' :
      set verbose_flg = 1 ; shift args; breaksw;
    case '--' : 
      shift args; breaksw;
  endsw
  if (junk${args[1]} == 'junk--') then
    shift args;
    break;
  endif
end
if ($#args != 1) then
  echo $0 ': only one filename permitted'
  exit -1
endif
set argv = $args

if ( $verbose_flg ) then
  echo '#'`date` $0
  echo '# $Id: allslice.csh,v 1.25 2004/12/15 01:49:20 welling Exp $'
endif

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

set ctblchunk = "`mri_printfield -field color_table -nofail $1`"
if ( dummy"${ctblchunk}" == 'dummy[chunk]' ) then
  set has_ctbl = 1
else
  set has_ctbl = 0
endif

# Use the ID tag as second title if available.
if (! ${?F_PS_TITLE2} && dummy$idtag != dummy ) then
  setenv F_PS_TITLE2 "ID $idtag"
endif

#Collect the undetermined options.  No gray scale is ever drawn with
#byte data, since that's what we use for overlays; range is also set 
#to byte range.
if ( $has_ctbl || ( $type == "uint8" ) ) then
  set allopts = "-mosaic -black_min 0 -black_max 255 -title -date_created" 
  set allopts = "$allopts -date_printed $ropt $bopt $dopt"
else
  set allopts = "-mosaic -title -date_created -date_printed"
  set allopts = "$allopts $maxopt $minopt $gopt $ropt $bopt $dopt"
endif

mri_to_img $allopts $1 $fname.ps

