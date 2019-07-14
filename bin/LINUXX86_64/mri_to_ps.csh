#!/bin/csh -ef
# mri_to_ps.csh
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
# *  Original programming by Bill Eddy                       *
# *  Mods by Joel Welling                                    *
# ************************************************************/
#
echo '#$Id: mri_to_ps.csh,v 1.14 2003/11/04 22:26:52 welling Exp $'
echo '#'`date`$0

# Parse command line options
set args = `getopt bdrg $*`
if ($#args < 2) then
  echo "usage: " $0 '[-b -d -r -g] filename'
  echo "  -b brightens image"
  echo "  -d dims image"
  echo "  -r reverses black and white"
  echo "  -g turns off gray scale bar"
  exit -1
endif
set ropt = "-reverse"
set gopt = "-gray_scale"
set bopt = ""
set dopt = ""
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
    case '--' : 
      shift args; breaksw;
  endsw
  if ($args[1] == '--') then
    shift args;
    break;
  endif
end
if ($#args != 1) then
  echo $0 ': only one filename permitted'
  exit -1
endif
set argv = $args

# Get some dimensions (assumed constant)
set zdim = `mri_printfield -field images.extent.z -nofail $1`
set tdim = `mri_printfield -field images.extent.t -nofail $1`
set type = `mri_printfield -field images.datatype -nofail $1`

# Does the input picture have a colormap?
set ctblchunk = `mri_printfield -field color_table -nofail $1`
if ( ${#ctblchunk} != 0 ) then
  set has_ctbl = 1
else
  set has_ctbl = 0
endif

# make a list of ints < zdim
set nslice_list = 0
set count = 1
while ($count < $zdim)
  set nslice_list = ( $nslice_list $count )
  @ count = $count + 1
end

# make a list of ints < t
set nimage_list = 0
set count = 1
while ($count < $tdim)
  set nimage_list = ( $nimage_list $count )
  @ count = $count + 1
end

# Convert the given .mri file to postscript, slice by slice.  We do all
# the extra subsetting to reproduce the old behavior of one .ps file
# for each input slice.
if( $has_ctbl || ($type == "uint8")) then
  set argstr = "-date_created -date_printed $ropt -title -black_min 0"
  set argstr = "$argstr -black_max 255"
else
  set argstr = "-date_created -date_printed $ropt -title $gopt"
endif
if (${?F_PS_TITLE2}) then
  set title2_set = 1
else
  set title2_set = 0
endif
foreach i ($nimage_list )
  echo Converting image $i
  if (! $title2_set) then
    setenv F_PS_TITLE2 "Image $i"
  endif
  mri_subset -d t -l 1 -s $i $1 tmp_mri_to_ps_1
  foreach j ($nslice_list)
    setenv F_PS_TITLE3 "Slice $j"
    mri_subset -d z -l 1 -s $j tmp_mri_to_ps_1 tmp_mri_to_ps_2
    if ($tdim > 1) then
      mri_to_img $argstr tmp_mri_to_ps_2 $1.$i.$j.ps
    else
      mri_to_img $argstr tmp_mri_to_ps_2 $1.$i.$j.ps
    endif
  end
end

# Clean up
rm tmp_mri_to_ps_1.* tmp_mri_to_ps_2.*

