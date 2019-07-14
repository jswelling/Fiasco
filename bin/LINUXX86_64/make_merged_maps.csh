#! /bin/csh -efx
#
# This script is called with a bunch of directory names as parameters.
# Everything gets aligned to the first dir given.  Using the Fisher
# method, the assumption is made that the total number of degrees of 
# freedom is large enough to use a standard normal distribution rather 
# than a T distribution.

echo '#'`date` $0
echo '# $Id: make_merged_maps.csh,v 1.1 2004/06/09 17:13:53 bakalj Exp $'

# Check for help command
if ( $#argv >= 1 ) then 
  if ( dummy$argv[1] == "dummy-help" ) then
    if ( $#argv >= 2 ) then
      scripthelp $0 $argv[2]
    else
      scripthelp $0
    endif
    exit
  endif
endif

#Parse command line options
set args = `getopt m: $*`
set help = ""
set method = ""
set dirs = ""
while ($#args > 1)
  switch ( $args[1] )
    case '-m' :
      set method = $args[2]; shift args; shift args; breaksw;
  endsw
  if (junk${args[1]} == 'junk--') then
    shift args;
    break;
  endif
end

if ($#args < 2 || dummy$help != dummy ) then
  scripthelp $0 usage
  exit -1
endif

while ($#args >= 1)
  set dirs= ($dirs $args[1])
  shift args
end

if ($method != "fisher" && $method != "stats" && $method != "stouffer") then
  echo "Invalid method ($method) specified.  Using Fisher method."
  set method = "fisher"
else
  echo "Combining using $method method."
endif

# Environment variables needed: F_XVOXEL, F_YVOXEL, F_ZVOXEL, and whatever
# makeps.csh needs.
source fiasco.local.csh
source $FIASCO/epi.defaults.csh
if ( -f epi.local.csh) source epi.local.csh
if ( -f spiral.local.csh) source spiral.local.csh
if ( -f ts.local.csh) source ts.local.csh

if (! -d merged) mkdir merged
if (! -d merged/par) mkdir merged/par
if (! -d merged/stat) mkdir merged/stat
if (! -d merged/ps) mkdir merged/ps

set count = 0


# Build aligned directories
foreach dirname ($dirs)
    echo Processing $dirname
    if (! -d $dirname/stat_aligned) mkdir $dirname/stat_aligned
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
        -input $align_this -align $alignto -stdv $alignto_stdv \
        -parameters $parfile -alg opt=none
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
    endif
    @ count = $count + 1
end

# Produce sums and SSEs
set cond = 1;
set xdim = `mri_printfield -field images.extent.x $dirs[1]/stat_aligned/Mean_1`
set ydim = `mri_printfield -field images.extent.y $dirs[1]/stat_aligned/Mean_1`
while (1)
    set dircount = 0;
    foreach dirname ($dirs)
      # Beware! not all subdirectories may have means for all conditions!
      if ( ! -f ${dirname}/stat_aligned/Mean_${cond}.mri ) continue
      if ( $dircount == 0 ) then
        if (-f merged/tmp.mri) rm merged/tmp.*
        mri_copy_chunk -chunk counts -chunk_out images -chunk_file .dat \
	  -replace $dirname/stat_aligned/Mean_${cond} merged/tmp
	mri_type_convert -float merged/tmp merged/tmp2
	mri_remap -order xyz merged/tmp2 
        mri_interp -d x -len $xdim -con merged/tmp2 merged/tmp3
        mri_interp -d y -len $ydim -con merged/tmp3 merged/this_count
        if (-f merged/count_sum.mri) rm merged/count_sum.*
	mri_copy_chunk -chunk images merged/this_count merged/count_sum
        if (-f merged/tmp.mri) rm merged/tmp.*
        mri_copy_chunk -chunk images -replace \
          $dirname/stat_aligned/Mean_${cond} merged/tmp 
        mri_rpn_math -out merged/sum_sum '$1,$2,*' \
          merged/tmp merged/this_count
        if (-f merged/tmp.mri) rm merged/tmp.*
        mri_copy_chunk -chunk images -replace \
          $dirname/stat_aligned/Stdv_${cond} merged/tmp
	mri_rpn_math -out merged/this_sse '$1,dup,*,$2,1.0,-,*' \
          merged/tmp merged/this_count
        if (-f merged/sse_sum.mri) rm merged/sse_sum.*
	mri_copy_chunk -chunk images merged/this_sse merged/sse_sum
      else
        if (-f merged/this_mean.mri) rm merged/this_mean.*
        mri_copy_chunk  -chunk images -replace \
          $dirname/stat_aligned/Mean_${cond} merged/this_mean 
        if (-f merged/tmp.mri) rm merged/tmp.*
        mri_copy_chunk -chunk counts -chunk_out images -chunk_file .dat \
          -replace $dirname/stat_aligned/Mean_${cond} merged/tmp
	mri_type_convert -float merged/tmp merged/tmp2
	mri_remap -order xyz merged/tmp2 
        mri_interp -d x -len $xdim -con merged/tmp2 merged/tmp3
        mri_interp -d y -len $ydim -con merged/tmp3 merged/this_count
        if (-f merged/tmp.mri) rm merged/tmp.*
        mri_copy_chunk -chunk images -replace \
          $dirname/stat_aligned/Stdv_${cond} merged/tmp
	mri_rpn_math -out merged/this_sse '$1,dup,*,$2,1.0,-,*' \
          merged/tmp merged/this_count
	mri_rpn_math -out merged/tmp '$1,$2,+' merged/count_sum merged/this_count
        if (-f merged/count_sum.mri) rm merged/count_sum.*
	mri_copy_chunk -chunk images merged/tmp merged/count_sum
        mri_rpn_math -out merged/tmp '$1,$2,$3,*,+' \
          merged/sum_sum merged/this_mean merged/this_count
        if (-f merged/sum_sum.mri) rm merged/sum_sum.*
	mri_copy_chunk -chunk images merged/tmp merged/sum_sum
        mri_rpn_math -out merged/tmp '$1,$2,+' merged/sse_sum merged/this_sse
        if (-f merged/sse_sum.mri) rm merged/sse_sum.*
	mri_copy_chunk -chunk images merged/tmp merged/sse_sum
      endif
      @ dircount = $dircount + 1
        rm merged/tmp.* merged/tmp2.* merged/tmp3.*
    end
    if ( $dircount == 0 ) break  # we are out of conditions

# Build the collective counts and mean.  We have to do some trickery
# to cause the counts to be saved in the .mri rather than .dat file.
    mri_rpn_math -out merged/stat/Mean_${cond} '$1,$2,/' \
      merged/sum_sum merged/count_sum
    mri_type_convert -int merged/count_sum merged/tmp
    mri_subset -d x -l 1 -s 0 merged/tmp merged/tmp2
    mri_subset -d y -l 1 -s 0 merged/tmp2 merged/tmp3
    mri_remap -order z merged/tmp3 
    if (-f merged/counts_int.mri) rm merged/counts_int.*
    mri_copy_chunk -chunk counts \
      ${dirname}/stat_aligned/Mean_${cond} merged/counts_int
    mri_copy_chunk -chunk images -chunk_out counts -replace \
      merged/tmp3 merged/counts_int
    mri_copy_chunk -chunk counts -chunk_out counts -replace \
      merged/counts_int merged/stat/Mean_${cond}
    rm merged/tmp.* merged/tmp2.* merged/tmp3.*


# Adjust the joint SSE for the variance of the means
    set dircount = 0
    foreach dirname ($dirs) 
      if (-f merged/tmp.mri) mri_destroy_dataset merged/tmp
      # Beware! not all subdirectories may have means for all conditions!
      if ( ! -f ${dirname}/stat_aligned/Mean_${cond}.mri ) continue
      mri_copy_chunk -chunk counts -chunk_out images -chunk_file .dat \
        $dirname/stat_aligned/Mean_${cond} merged/tmp
      mri_type_convert -float merged/tmp merged/tmp2
      mri_remap -order xyz merged/tmp2
      mri_interp -d x -len $xdim -con merged/tmp2 merged/tmp3
      mri_interp -d y -len $ydim -con merged/tmp3 merged/this_count
      mri_rpn_math -out merged/tmp '$2,$3,-,dup,*,$4,*,$1,+' \
        merged/sse_sum $dirname/stat_aligned/Mean_${cond} \
        merged/stat/Mean_${cond} merged/this_count
      if (-f merged/sse_sum.mri) rm merged/sse_sum.*
      mri_copy_chunk -chunk images merged/tmp merged/sse_sum
      @ dircount = $dircount + 1
    end
    rm merged/tmp.* merged/tmp2.* merged/tmp3.*

# Build the final stdv
    mri_rpn_math -out merged/stat/Stdv_${cond} '$1,$2,1.0,-,/,sqrt' \
      merged/sse_sum merged/count_sum
    mri_copy_chunk -chunk counts -chunk_out counts -replace \
      merged/counts_int merged/stat/Stdv_${cond}

# Continue with condition loop
    @ cond = $cond + 1
end

# Build Tmaps

@ maxcond = ( $cond - 1 )
set cond1 = 1
  switch ( $method )
    case 'stats' :
      while ( $cond1 < $maxcond )
	@ cond2 = ( $cond1 + 1 )

	set mean1 = merged/stat/Mean_${cond1}
	set stdv1 = merged/stat/Stdv_${cond1}
	if (-f merged/tmp.mri) rm merged/tmp.*
	mri_copy_chunk -chunk counts -chunk_out images -chunk_file .dat \
	    $mean1 merged/tmp
	mri_type_convert -float merged/tmp merged/tmp2
	mri_remap -order xyz merged/tmp2 
	mri_interp -d x -len $xdim -con merged/tmp2 merged/tmp3
	mri_interp -d y -len $ydim -con merged/tmp3 merged/this_count1
	@ cond2 = ( $cond1 + 1 )

	while ( $cond2 <= $maxcond )

	  echo Building Tmap_${cond1}-${cond2}
          set mean2 = merged/stat/Mean_${cond2}
          set stdv2 = merged/stat/Stdv_${cond2}
	  set tmap = merged/stat/Tmap_${cond1}-${cond2}
          if (-f merged/tmp.mri) rm merged/tmp.*
          mri_copy_chunk -chunk counts -chunk_out images -chunk_file .dat \
            $mean1 merged/tmp
          mri_type_convert -float merged/tmp merged/tmp2
          mri_remap -order xyz merged/tmp2 
          mri_interp -d x -len $xdim -con merged/tmp2 merged/tmp3
          mri_interp -d y -len $xdim -con merged/tmp3 merged/this_count2

# We need the normalized pooled standard deviation.  Read a long 
# mri_rpn_math script from a file.
          cat > merged/tmp_script.t <<EOF
            \$1,dup,*,\$2,1.0,-,*,\$3,dup,*,\$4,1.0,-,*,+
            \$2,\$4,+,2,-,/
            1.0,\$2,/,1.0,\$4,/,+,*
            sqrt
EOF
	  mri_rpn_math -out merged/this_pooled_stdv -exp merged/tmp_script.t \
            $stdv1 merged/this_count1 $stdv2 merged/this_count2
          rm merged/tmp_script.t
          mri_rpn_math -out $tmap '$1,$2,-,$3,/' \
            $mean1 $mean2 merged/this_pooled_stdv
	  mri_copy_chunk -chunk counts -chunk_out counts1 \
	    $mean1 $tmap
	  mri_copy_chunk -chunk counts -chunk_out counts2 \
            $mean2 $tmap

          @ cond2 = ( $cond2 + 1 )
	  
          rm merged/tmp.* merged/tmp2.* merged/tmp3.*
        end
	@ cond1 = $cond1 + 1
      end
    breaksw;

    case 'fisher' :
      while ( $cond1 < $maxcond )
	@ cond2 = ( $cond1 + 1 )
	while ( $cond2 <= $maxcond )
	  echo Building Tmap_${cond1}-${cond2}
	  set pmap = Pmap_${cond1}-${cond2}
	  set tmap = Tmap_${cond1}-${cond2}
	
	  set filelist = ''
	  foreach dirname ($dirs) 
	    set filelist = "$filelist ${dirname}/stat_aligned/${tmap}"
	  end

	  combine.csh -o merged/stat/${pmap} -m $method $filelist

	  @ cond2 = $cond2 + 1
        end
	@ cond1 = $cond1 + 1
      end
    breaksw

    case 'stouffer' :

      while ( $cond1 < $maxcond )
	@ cond2 = ( $cond1 + 1 )
	while ( $cond2 <= $maxcond )
	  echo Building Tmap_${cond1}-${cond2}
	  set pmap = Pmap_${cond1}-${cond2}
	  set tmap = Tmap_${cond1}-${cond2}
	
	  set filelist = ''
	  foreach dirname ($dirs) 
	    set filelist = "$filelist ${dirname}/stat_aligned/${tmap}"
	  end
	  combine.csh -o merged/stat/${pmap} -m $method $filelist

	  @ cond2 = $cond2 + 1
        end
	@ cond1 = $cond1 + 1
      end
      breaksw

    endsw

foreach ds ( merged/*.mri )
  mri_destroy_dataset $ds
end

# Make the Grand Mean and Stdv
pooled_mean.csh -out merged/stat/GrandMean merged/stat/Mean_*.mri
pooled_stdv.csh -out merged/stat/GrandStdv merged/stat/Stdv_*.mri

# Make some pictures

cd merged
makeps.csh

# print postscript files
print.files.csh default all $F_PRINTER


