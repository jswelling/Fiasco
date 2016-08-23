#! /bin/csh -f
# make_phase_cine.csh
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
echo '#'`date` $0
echo '# $Id: make_phase_cine.csh,v 1.8 2004/08/05 20:26:23 welling Exp $'
#
# This script generates a cine loop from a (v)xyzt dataset.
# The one argument is the name of a single scalar dataset file.
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
set args = `getopt bdr $*`
if ($#args < 2 ) then
  exit -1
endif
set ropt = ""
set bopt = ""
set dopt = ""
set nopt = ""
set xopt = ""
while ($#args > 1) 
  switch ( $args[1] )
    case '-b' : 
      set bopt = "-b" ; shift args; breaksw;
    case '-d' : 
      set dopt = "-d" ; shift args; breaksw;
    case '-r' : 
      set ropt = "-r" ; shift args; breaksw;
    case '--' : 
      shift args; breaksw;
  endsw
  if (junk${args[1]} == 'junk--') then
    shift args;
    break;
  endif
end
if ($#args != 1) then
  echo ${0:t} ': only one filename permitted'
  scripthelp $0
  exit -1
endif
set argv = $args

# Check that the needed "convert" and "identifies" utilities are present
foreach util ( convert identify )
  if ( `test_in_subshell.csh which ${util} ` ) then
    echo ${0:t} 'Error: The '${util}' utility was needed but was not found!'
    exit -1
  endif
end

# Get the root name, to construct output filename
if ($1:t == ".mri") then
  set tailname = $1:t
  set fname = $tailname:r
else
  set fname = $1:t
endif

set vdim = `mri_printfield -field images.extent.v -nofail $1`
if ( foo$vdim == foo ) set vdim = 1
set xdim = `mri_printfield -field images.extent.x $1`
set ydim = `mri_printfield -field images.extent.y $1`
set zdim = `mri_printfield -field images.extent.z -nofail $1`
if ( foo$zdim == foo ) set zdim = 1
set tdim = `mri_printfield -field images.extent.t -nofail $1`
if ( foo$tdim == foo ) set tdim = 1

@ qdim = $xdim * $ydim * $zdim

if ( $vdim != 2 ) then
  echo "Input file $1 is not complex!"
  exit -1
endif

#
# Make some scratch space, and go there.
#
if ( ${?F_TEMP} ) then
  set tmpdir = ${F_TEMP}/make_phase_cine_tmp_$$
else
  set tmpdir = ./make_phase_cine_tmp_$$
endif
if (! -e $tmpdir) mkdir $tmpdir
set homedir = `pwd`
if ( $1 =~ /* ) then
  set input = $1
else
  set input = $homedir/$1
endif
cd $tmpdir
mkdir frames

#Collect the undetermined options.  
set allopts = "-p $ropt $bopt $dopt"

set t = 0

@ xlim = $xdim - 1
@ ylim = $ydim - 1

# We need the maximum input magnitude

mri_complex_to_scalar -mag $input ${fname}_m
mri_permute -order vtxyz ${fname}_m ${fname}_p 
mri_destroy_dataset ${fname}_m
mri_remap -order tq -length ${tdim}:${qdim} ${fname}_p 
mri_subsample -d t -l 1 -max ${fname}_p ${fname}_t
set missing_chunk = `mri_printfield -field missing -nofail ${fname}_t`
if ( foo${#missing_chunk} != foo ) then
  if ( ${#missing_chunk} != 0 ) then
    mri_delete_chunk -chunk missing ${fname}_t 
  endif
endif
mri_subsample -d q -l 1 -max ${fname}_t max
set maxval = `mri_rpn_math -out junk '0,$1,1,if_print_1' max`

echo "Generating $tdim frames..."
@ xdim_pad = $xdim + 2
@ ydim_pad = $ydim + 2
while ( $t < $tdim )
  echo -n "${t}... "
  if ( ! ( $t % 10 ) && $t != 0 ) echo ""
  mri_subset -d t -l 1 -s $t $input thistime
  mri_pad -d x -l $xdim_pad -s 1 -fillvalue $maxval thistime thistime_padx
  mri_pad -d y -l $ydim_pad -s 1 -fillvalue $maxval thistime_padx tweaked_image

  # We have to jigger with the output file names to make sure they
  # sort to the right order.  This version will break if there are
  # more than 10000 frames, but it should hold us for a while!
  if ( $t < 10 ) then
    set oname = frames/frame_000$t.png
  else if ( $t < 100 ) then
    set oname = frames/frame_00$t.png
  else if ( $t < 1000 ) then
    set oname = frames/frame_0$t.png
  else
    set oname = frames/frame_$t.png
  endif


  set missing_chunk = `mri_printfield -field missing -nofail tweaked_image`
  if ( foo${#missing_chunk} != foo ) then
    if ( ${#missing_chunk} != 0 ) then
      mri_delete_chunk -chunk missing tweaked_image
    endif
  endif

  color_by_phase.csh ${allopts} tweaked_image
  if ( $t == 0 ) then
    set img_x = `identify -format '%w' tweaked_image.png`
    set img_y = `identify -format '%h' tweaked_image.png`
    if ( $img_x > 640 || $img_y > 480 ) then
      echo  ${0:t} 'Error: resulting image resolution '${img_x}'x'${img_y}' \
        ' is too large!'
      cd $homedir
      rm -r $tmpdir
      exit -1
    endif
  endif
  mv tweaked_image.png ${oname}
  @ t = $t + 1
end
echo "done"

echo "Assembling mpeg..."

cd frames
convert -scale 200%! frame_*.png ${homedir}/${fname}.mpg

echo "Done."

#clean up
cd $homedir
rm -r $tmpdir

