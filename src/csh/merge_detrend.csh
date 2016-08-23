#! /bin/csh -exf
#
# This script is called with a bunch of directory names as parameters.
# Everything gets aligned to the first dir given.  The assumption is
# made that the total number of degrees of freedom is large enough to
# use a standard normal distribution rather than a T distribution.
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

# Environment variables needed: F_XVOXEL, F_YVOXEL, F_ZVOXEL, and whatever
# makeps.csh needs.
source fiasco.local.csh
source $FIASCO/epi.defaults.csh
if ( -f epi.local.csh) source epi.local.csh
if ( -f spiral.local.csh) source spiral.local.csh
if ( -f ts.local.csh) source ts.local.csh

if (! -d merged) mkdir merged
if (! -d merged/in) mkdir merged/in
if (! -d merged/par) mkdir merged/par
if (! -d merged/stat) mkdir merged/stat
if (! -d merged/data) mkdir merged/data
if (! -d merged/ps) mkdir merged/ps
if (! -d merged/tmp) mkdir merged/tmp

# Build aligned directories.  If necessary we use "stats" with a temporary 
# fake newsplit file to create a grand mean and stdv used for alignment.
set count = 0
foreach dirname ($*)
    echo Processing $dirname
    if (! -d $dirname/stat_aligned) mkdir $dirname/stat_aligned
    if (! -d $dirname/data_aligned) mkdir $dirname/data_aligned 
    if ( -e $dirname/stat/GrandMean.mri ) then
      set align_this = $dirname/stat/GrandMean
    else
      set align_this = $dirname/stat/Mean_1
    endif
    if ( -e $dirname/stat/GrandStdv.mri ) then
      set align_this_stdv = $dirname/stat/GrandStdv
    else
      set align_this_stdv = $dirname/stat/Stdv_1
    endif
    if ( $count == 0 ) then
      cp $dirname/stat/* $dirname/stat_aligned
      cp $dirname/data/detrend.* $dirname/data_aligned
      set first = $dirname 
      set alignto = $align_this
      set alignto_stdv = $align_this_stdv
    else
      set ct = `test_in_subshell.csh diff -w -i -q ${first}/$F_SPLIT_COND ${dirname}/$F_SPLIT_COND`
      if ( $ct != 0 ) then
        echo "Error: Order of experimental conditions in " ${dirname} " does not match " ${first}
        exit -1
      endif 
      set parfile = merged/par/reg3d_${dirname}.par
      estireg3d -xvoxel ${F_XVOXEL} -yvoxel ${F_YVOXEL} -zvoxel ${F_ZVOXEL} \
        -input $align_this \
        -align $alignto -stdv $alignto_stdv \
	-parameters $parfile
      foreach statfile ( ${dirname}/stat/*.mri ) 
        if (( $statfile !~ *_?.mri ) || ( $statfile =~ *Mean*.mri ) \
	    || ( $statfile =~ *Stdv*.mri )) then
          echo Aligning $statfile
	  mri_remap -order xyzt $statfile 
          ireg3d -xvoxel ${F_XVOXEL} -yvoxel ${F_YVOXEL} -zvoxel ${F_ZVOXEL} \
	    -input $statfile -headerout ${dirname}/stat_aligned/${statfile:t} \
	    -parameters $parfile
          endif
      end
# A stretched par file is needed to align the many images of the detrend file.
      set stretchedpar = merged/par/reg3d_stretch_${dirname}.par
      set tdim = `mri_printfield -field images.extent.t ${dirname}/data/detrend`
      if (-f merged/tmp.mri) rm merged/tmp.*
      if (-f merged/tmp2.mri) rm merged/tmp2.*
      if (-f merged/tmp3.mri) rm merged/tmp3.*
      if (-f merged/tmp4.mri) rm merged/tmp4.*
      if (-f merged/tmp5.mri) rm merged/tmp5.*
      if (-f merged/tmp6.mri) rm merged/tmp6.*
      if (-f merged/tmp7.mri) rm merged/tmp7.*
      if (-f merged/tmp8.mri) rm merged/tmp8.*
      if (-f merged/tmp9.mri) rm merged/tmp9.*
      if (-f merged/tmp10.mri) rm merged/tmp10.*
      mri_from_ascii -ord vt -l 8:1 -ind t merged/tmp < $parfile
      mri_interp -con -d t -len $tdim merged/tmp merged/tmp2
      mri_subset -d v -l 1 -s 0 merged/tmp2 merged/tmp3
      mri_subset -d v -l 1 -s 1 merged/tmp2 merged/tmp4
      mri_subset -d v -l 1 -s 2 merged/tmp2 merged/tmp5
      mri_subset -d v -l 1 -s 3 merged/tmp2 merged/tmp6
      mri_subset -d v -l 1 -s 4 merged/tmp2 merged/tmp7
      mri_subset -d v -l 1 -s 5 merged/tmp2 merged/tmp8
      mri_subset -d v -l 1 -s 6 merged/tmp2 merged/tmp9
      mri_subset -d v -l 1 -s 7 merged/tmp2 merged/tmp10
      grep '#' $parfile > $stretchedpar
      mri_rpn_math '$t,$1,$2,$3,$4,$5,$6,$7,$8,1,if_print_9' \
	 merged/tmp3 merged/tmp4 merged/tmp5 merged/tmp6 \
	 merged/tmp7 merged/tmp8 merged/tmp9 merged/tmp10 >> $stretchedpar
      ireg3d -xvoxel ${F_XVOXEL} -yvoxel ${F_YVOXEL} -zvoxel ${F_ZVOXEL} \
	  -input ${dirname}/data/detrend \
          -headerout ${dirname}/data_aligned/detrend \
	  -parameters $stretchedpar
#      ireg3d -xvoxel ${F_XVOXEL} -yvoxel ${F_YVOXEL} -zvoxel ${F_ZVOXEL} \
#	  -input ${dirname}/data/detrend \
#          -headerout ${dirname}/data_aligned/detrend \
#	  -parameters $stretchedpar -trilinear
#      mri_copy_dataset ${dirname}/data/detrend \
#         ${dirname}/data_aligned/detrend
    endif
    @ count = $count + 1
end
rm merged/tmp.* merged/tmp2.* merged/tmp3.* merged/tmp4.* merged/tmp5.* 
rm merged/tmp6.* merged/tmp7.* merged/tmp8.* merged/tmp9.* 
rm merged/tmp10.* 

# Produce scaling factors relating the different sessions, since the
# overall gain may have changed between sessions
set count = 0
foreach dirname ($*)
  set basefile = ${dirname}/stat_aligned/GrandMean
  set tfile = merged/tmp/dirmean${count}
  set xdim = `mri_printfield -field images.extent.x $basefile`
  @ xhalf = $xdim / 2
  @ xqtr = $xdim / 4
  set ydim = `mri_printfield -field images.extent.y $basefile`
  @ yhalf = $ydim / 2
  @ yqtr = $ydim / 4
  set zdim = `mri_printfield -field images.extent.z $basefile`
  mri_subset -d x -l $xhalf -s $xqtr $basefile merged/tmp/tmp
  mri_delete_chunk -chunk counts merged/tmp/tmp 
  mri_subsample -mean -d x -l 1 merged/tmp/tmp merged/tmp/tmp2
  mri_subset -d y -l $yhalf -s $yqtr merged/tmp/tmp2 merged/tmp/tmp3
  mri_subsample -mean -d y -l 1 merged/tmp/tmp3 merged/tmp/tmp4
  mri_subsample -mean -d z -l 1 merged/tmp/tmp4 $tfile
  foreach fname ( tmp tmp2 tmp3 tmp4 )
    mri_destroy_dataset merged/tmp/$fname
  end
  @ count = $count + 1
end

# Produce the merged aligned detrend file by concatenation.
set count = 0
foreach dirname ($*)
    if ( $count == 0 ) then
	mri_copy_dataset ${dirname}/data_aligned/detrend merged/data/detrend
	cp ${dirname}/in/split merged/in/split
    else
	mri_rpn_math -out merged/tmp/singledetrend '$1,$2,*,$3,/' \
	    ${dirname}/data_aligned/detrend \
	    merged/tmp/dirmean0 merged/tmp/dirmean${count}
#	mri_rpn_math -out merged/tmp/singledetrend '$1' \
#	    ${dirname}/data_aligned/detrend \
#	    merged/tmp/dirmean0 merged/tmp/dirmean${count}
	mri_paste -d t -out merged/tmp/detrend \
	    merged/data/detrend merged/tmp/singledetrend
	mri_destroy_dataset merged/tmp/singledetrend
	mv merged/tmp/detrend.* merged/data
	tail +3 ${dirname}/in/split >> merged/in/split
    endif
    @ count = $count + 1
end

# Things are set up to expect a permute file as well.
mri_permute -order vtxyz merged/data/detrend merged/data/permute

# Produce a merged split file by concatenation.
set splitlist = ()
set totalimages = 0
set nslices = `mri_printfield -fld 'images.extent.z' $1/data/detrend`
foreach dirname ($*)
  set splitlist = ( $splitlist ${dirname}/in/split )
  @ totalimages = $totalimages + \
    `mri_printfield -fld images.extent.t ${dirname}/data/detrend`
end
paste_split.py -v --out merged/in/split --nimages $totalimages \
   --nslices $nslices $splitlist

#
# The following processing steps just do stats and make plots as per 
# the normal Fiasco processing stream.
#

cd merged
setenv F_DETREND_OUTPUT detrend
setenv F_DETREND_TMP permute

# Do experiment statistics
$F_STATS_SCRIPT data/$F_DETREND_OUTPUT data/$F_DETREND_TEMP

# No displacements, because we don't have overall transformations.
#
# No unified parameters, so no summary pages
#
# Let's not mail summary to Fiasco Headquarters
#
# No plots because we don't have unified parameters.
#
# print postscript files
print.files.csh default all $F_PRINTER


