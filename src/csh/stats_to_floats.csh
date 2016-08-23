#! /bin/csh -ef
#
# This script is called from a Fiasco directory and expects to find
# a directory named "stat" filled with Pgh MRI files.  It converts
# appropriate files from double precision to float and puts the
# results in a new directory stat_float.
#
# More precisely, it goes through the environment looking for 
# variables of the form F_STAT_?PATH and converts any Pgh MRI
# file in those directories.  All the results end up in 'stat_float'.
#

# Check for help command
if ( $#argv >= 1 ) then 
  if ( $argv[1] == "-help" ) then
    if ( $#argv >= 2 ) then
      scripthelp $0 $argv[2]
    else
      scripthelp $0
    endif
    exit
  endif
endif

# Environment variables needed: F_STAT_SPATH, F_STAT_TPATH, F_STAT_MPATH
source fiasco.local.csh
source $FIASCO/epi.defaults.csh
if ( -f epi.local.csh) source epi.local.csh
if ( -f spiral.local.csh) source spiral.local.csh
if ( -f ts.local.csh) source ts.local.csh

set outdir = stat_float

if (! -d $outdir) mkdir $outdir

set stat_paths = \
    `setenv | grep 'F_STAT_.PATH' | awk -F '=' '{ print $2 }' | sort -u`
set datasets = ( )
foreach dir ( $stat_paths )
  echo "Processing $dir"
  foreach fname ( $dir/*.mri )
    set datasets = ( $datasets ${fname:r} )
  end
end
foreach ds ( $datasets )
    if ( ${ds:t} =~ *_[ftq] ) then
      echo "Ignoring ${ds:t}"
    else
      set newname = ${outdir}/${ds:t}
      if ( ${ds:t} =~ Pmap_* ) then
	echo "Copying ${ds:t}"
	mri_copy_dataset $ds $newname
      else
        set type = `mri_printfield -field 'images.datatype' -nofail $ds`
        if ( dummy$type == dummyfloat64 ) then
	  echo "Converting ${ds:t}"
	  mri_type_convert -float $ds $newname
        else
	  echo "Copying ${ds:t}"
	  mri_copy_dataset $ds $newname
        endif
      endif
    endif
end

