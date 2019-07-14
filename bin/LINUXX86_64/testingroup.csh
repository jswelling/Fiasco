#! /bin/csh -ef
#
# This script does between-group comparisons.
#
echo '#'`date` $0
echo '# $Id: testingroup.csh,v 1.10 2003/09/16 17:31:58 bakalj Exp $'

# Check for help command
if ( $#argv >= 1 ) then 
  if ( junk$argv[1] == junk-help ) then
    if ( $#argv >= 2 ) then
      scripthelp $0 $argv[2]
    else
      scripthelp $0
    endif
    exit
  endif
endif

# Parse command line options
set args = `getopt m:d:h $*`
set help = ""
set mapstuff = ""
while ($#args > 1) 
  switch ( $args[1] )
    case '-m' : 
      if ( ${?F_MAP_PATH} ) then
        setenv F_MAP_PATH $args[2]':'${F_MAP_PATH}
      else
        setenv F_MAP_PATH $args[2]
      endif
      shift args; shift args; breaksw;
    case '-d' :
      set mapstuff = $mapstuff" "$args[2]
      shift args; shift args; breaksw;
    case '-h' :
      set help = 1 ; shift args; breaksw;
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

#
# Make some scratch space, and go there.
#
if ( ${?F_TEMP} ) then
  set tmpdir = ${F_TEMP}/testingroup_tmp_$$
else
  set tmpdir = ./testingroup_tmp_$$
endif

if (! -e $tmpdir) mkdir $tmpdir
echo "Scratch directory will be " $tmpdir
set homedir = `pwd`
echo "Current directory is " $homedir
cd $tmpdir

# Crawl around, counting conditions and groups

echo '######## Counting conditions and classifying subjects by group ########'

set nsubjs = ${#args}
set nconds = 0
set ngroups = 0
set groupnames = ( )
set sbj_grp_ids = ( )

foreach subj ($args)
  set cond = 1
  while (1)
    set mean_fname = `map_name.py -d cond=$cond -d file=Mean $subj`
    # echo 'testing' $mean_fname
    if (! -e ${mean_fname}.mri) break
    set stdv_fname = `map_name.py -d cond=$cond -d file=Stdv $subj`
    if (! -e ${stdv_fname}.mri) then
      echo "ERROR: the following file seems to be missing: " ${stdv_fname}.mri
      exit -1
    endif
    set count_fname = `map_name.py -d cond=$cond -d file=Count $subj`
    if (! -e ${count_fname}.mri) then
      echo "ERROR: the following file seems to be missing: " ${count_fname}.mri
      exit -1
    endif
    @ cond = $cond + 1
  end
  if ( $cond == 1 ) then
    echo "ERROR: no mean data for subject " $subj
    exit -1
  endif
  if ( $cond > $nconds ) then
    @ nconds = $cond
  endif
  set groupid = 0
  set subjgroup = `map_name.py -t group $subj`
  # echo 'subjgroup is ' $subjgroup
  foreach group ( $groupnames )
    # echo 'testing group ' $group ' vs ' $subjgroup
    if ( $group == $subjgroup ) break;
    @ groupid= $groupid + 1
  end
  if ( $groupid + 1 > $ngroups ) then
    @ ngroups = $ngroups + 1
    set groupnames = ( $groupnames $subjgroup )
  endif
  set sbj_grp_ids = ( $sbj_grp_ids $groupid )
end
@ nconds = $nconds - 1
set lastsubj = $subj
echo 'nsubjs is ' $nsubjs
echo 'nconds is ' $nconds
echo 'ngroups is ' $ngroups
echo 'group names: ' $groupnames
echo 'group ids: ' $sbj_grp_ids

if ( $nsubjs < 2 ) then
  echo "You need at least 2 subjects!"
  exit -1
endif

if ( $ngroups != 1 ) then
  echo "All of the subjects must be in the same group!"
  exit -1
endif

if ( $nconds != 2 ) then
  echo "There must be exactly 2 experimental conditions!"
  exit -1
endif

@ tdim = $nsubjs * $nconds

# Construct the X matrix (factors).  conds_x_subjects has entries which
# correspond to each subject's group ID.
#
# Factors:
# - factor1 will be -1 for condition 0, 1 for condition 1
# - subjfactorn will be 1 for subject n, 0 otherwise, de-meaned, skipping
#   the last subject (to stay orthogonal to the mean)
#

echo '######## Generating Factors ########'

echo $sbj_grp_ids | mri_from_ascii -ord cs -l 1:${nsubjs} junk3
mri_interp -con -d c -len $nconds junk3 conds_x_subjects

mri_rpn_math -out rawfactor1 '$c,2,*,1,-' conds_x_subjects

set which_subj = 0
set rawfactorlist = ( rawfactor1 )
set effectfactor = 1
set subjfactorlist = ( )
foreach subj ($args)
  if ( $subj == $lastsubj ) break
  mri_rpn_math -out normsubjfactor$which_subj \
    '$s,'${which_subj}',==,1.0,'$nsubjs',/,-,$2,*' \
    conds_x_subjects rawfactor1
  set subjfactorlist = ( $subjfactorlist normsubjfactor${which_subj} )
  @ which_subj = $which_subj + 1  
end
@ nrawfactors = ${#rawfactorlist}
@ nsubjfactors = ${#subjfactorlist}
@ nfactors = $nrawfactors + $nsubjfactors
echo "number of fixed effect factors is " $nrawfactors
echo "number of factors associated with individual subjects is " $nsubjfactors
echo "total nfactors is " $nfactors
echo "unscaled factors are " $rawfactorlist $subjfactorlist
echo "The factor of interest is number " $effectfactor

# Remap the factors into the order needed for GLM
foreach fname ( $rawfactorlist $subjfactorlist )
  mri_remap -order t -length ${tdim} $fname
end

# Make the factor part of the operand string for mri_glm
set rawfactoropts = ( )
foreach fname ( $rawfactorlist )
  set rawfactoropts = ( $rawfactoropts ${fname}:images )
end
set subjfactoropts = ( )
foreach fname ( $subjfactorlist )
  set subjfactoropts = ( $subjfactoropts ${fname}:images )
end

# Build the scaling matrix and the vector of mean values.  There are 
# nconditions*ngroups terms in each.  The trick here is to remember that
# stdv can be 0.0;  if that happens we want to substitute 0.0 for the Y
# value.  We sum up the counts as we go, for later use in computing the
# total degrees of freedom.  
#
# The missing chunks are deleted on the fly.

echo '######## Building raw Y and scaling vectors ########'

set scale_args = ( )
set rawy_args = ( )
foreach subj ($args)
  echo "Subject " $subj
  foreach cond ( 1 2 )
    set mean_fname = `map_name.py -d cond=$cond -d file=Mean $subj`
    set stdv_fname = `map_name.py -d cond=$cond -d file=Stdv $subj`
    set counts_fname = `map_name.py -d cond=$cond -d file=Count $subj`
    set tfile1 = junk2_${subj}_${cond}
#    mri_rpn_math -out $tfile1 \
#      '0.0,$1,dup,1.0,<=,if_keep,sqrt,1.0,$2,dup,0.0,!=,if_keep,/,0.0,$2,0.0,==,if_keep' \
#      $counts_fname $stdv_fname
    mri_rpn_math -out $tfile1 '0.0,$1,dup,1.0,<=,if_keep,sqrt' \
      $counts_fname $stdv_fname
    mri_remap -order vtxyz $tfile1
    mri_delete_chunk -chunk missing $tfile1 
    set scale_args = ( $scale_args $tfile1 )
    set tfile2 = junk3_${subj}_${cond}
    mri_copy_dataset ${mean_fname} ${tfile2}
    mri_remap -order vtxyz $tfile2
    mri_delete_chunk -chunk missing $tfile2 
    set rawy_args = ( $rawy_args $tfile2 )
    if ( ( $subj == $args[1] ) && ( $cond == 1 ) ) then
      mri_copy_dataset ${counts_fname} totalcounts
      mri_remap -order xyz totalcounts
    else
      mri_rpn_math -out junk4 '$1,$2,+' totalcounts ${counts_fname}
      mri_copy_dataset junk4 totalcounts
      mri_destroy_dataset junk4
    endif
  end
end

# We can't do all the pastes in one go because there may be too many
# for the limit compiled into mri_paste, and we don't want to do them
# all separately because it involves too many redundant copies.

echo "Assembling scale vector"
mri_copy_dataset ${scale_args[1]} scale
shift scale_args
while ( $#scale_args )
  if ( $#scale_args > 20 ) then
    set these_args = ( $scale_args[1-20] )
    set scale_args = ( $scale_args[21-] )
  else
    set these_args = ( $scale_args )
    set scale_args = ( )
  endif
  mri_paste -d t -out tmp_paste scale $these_args
  mri_copy_dataset tmp_paste scale
  mri_destroy_dataset tmp_paste
end

echo "Assembling raw Y vector"

mri_copy_dataset ${rawy_args[1]} rawy
shift rawy_args
while ( $#rawy_args )
  if ( $#rawy_args > 20 ) then
    set these_args = ( $rawy_args[1-20] )
    set rawy_args = ( $rawy_args[21-] )
  else
    set these_args = ( $rawy_args )
    set rawy_args = ( )
  endif
  mri_paste -d t -out tmp_paste rawy $these_args
  mri_copy_dataset tmp_paste rawy
  mri_destroy_dataset tmp_paste
end

# Save some space
rm junk1_* junk2_* junk3_* 

# Remap these two massive datasets into the order needed by GLM.
echo '######## Remapping Y, stdv and scaling vectors ########'
set xdim = `mri_printfield -field images.extent.x rawy `
set ydim = `mri_printfield -field images.extent.y rawy `
set zdim = `mri_printfield -field images.extent.z rawy `
@ fakezdim = $xdim * $ydim * $zdim
mri_remap -order tz -length ${tdim}:${fakezdim} rawy
mri_remap -order tz -length ${tdim}:${fakezdim} scale

# Calculate the GLM over the raw (effect) factors
echo "######## Calculating the first GLM- this may take a while ########"
echo "The command being used is:"
echo mri_glm -v -ssqr -e raw_glm_results:images -out raw_residuals \
    -scale scale rawy $rawfactoropts
mri_glm -v -ssqr -e raw_glm_results:images -out raw_residuals \
    -scale scale rawy $rawfactoropts
echo '######## First GLM complete ########'

# Map the GLM results back to appropriate dimensions
set vdim = `mri_printfield -field images.extent.v raw_glm_results`
mri_remap -order vxyz -length ${vdim}:${xdim}:${ydim}:${zdim} \
    raw_glm_results 

# Calculate the GLM over the subject factors
echo "######## Calculating the second GLM- this may take a while ########"
echo "The command being used is:"
echo mri_glm -v -ssqr -e subj_glm_results:images \
    -scale scale raw_residuals $subjfactoropts
mri_glm -v -ssqr -e subj_glm_results:images \
    -scale scale raw_residuals $subjfactoropts
echo '######## Second GLM complete; constructing results ########'

# Map the GLM results back to appropriate dimensions
set vdim = `mri_printfield -field images.extent.v subj_glm_results`
mri_remap -order vxyz -length ${vdim}:${xdim}:${ydim}:${zdim} \
    subj_glm_results 
# Grab interesting parts of the GLM.  The output of the model consists
# of (nfactors+1) estimated b's, followed by the associated summed squares.
mri_subset -d v -l 1 -s $effectfactor raw_glm_results effect_b
@ ssto_offset = $nrawfactors + 2
mri_subset -d v -l 1 -s $ssto_offset raw_glm_results ssto
@ sse_offset = $nrawfactors + 3
mri_subset -d v -l 1 -s $sse_offset raw_glm_results sse
@ essqr_offset = $effectfactor + $nrawfactors + 3
mri_subset -d v -l 1 -s $essqr_offset raw_glm_results ssb
@ subj_ssto_offset = $nsubjfactors + 2
mri_subset -d v -l 1 -s $subj_ssto_offset subj_glm_results subj_ssto
@ subj_sse_offset = $nsubjfactors + 3
mri_subset -d v -l 1 -s $subj_sse_offset subj_glm_results subj_sse
mri_subset -d v -l 1 -s 0 raw_glm_results $homedir/mean

echo "######## Generating final statistics ########"

# The number of degrees of freedom in the final sse is equal to the number of 
# input cases minus the number of fit parameters, or:
# (nsubjs*nconds) - (nfactors + 1).  
# We use this to calculate the final MSE.
@ subj_sse_dof  = ( $nsubjs * $nconds) - ( $nfactors + 1 )
mri_rpn_math -out subj_mse '$1,'$subj_sse_dof',/' subj_sse

# sssg is summed squared error for subject within group; we calculate it 
# from SSTO = SSR + SSE (with a check for negative values due to floating
# point error). Degrees of freedom is the number of such columns in the 
# factor matrix, or:
# $nsubjs - 1 .
@ sssg_dof = $nsubjs - 1
mri_rpn_math -out sssg '$1,$2,-,dup,0.0,swap,0.0,>,if_keep' \
  subj_ssto subj_sse
mri_rpn_math -out mssg '$1,'$sssg_dof',/' sssg

# Hop back home
cd $homedir

# Construct the results.  
#
# Following NKNW ch. 24, the relevant test statistic is MSB/MSAB where
# factor B is the experimental condition (fixed) and factor A is the
# subject identifier (random).  This gives a pretty simple test statistic.
# We do an additional small song-and-dance to cause the Pmaps to be saved
# in double precision.
mri_rpn_math -out $homedir/effect_f '$1,$2,/' $tmpdir/ssb $tmpdir/mssg
mri_type_convert -double $homedir/effect_f $tmpdir/tmp_doubles
mri_rpn_math -out $homedir/effect_p '$1,1,'$sssg_dof',cf' $tmpdir/tmp_doubles

# Let's keep totalcounts as well; it's useful for seeing where the
# samples were actually taken.
mri_copy_dataset $tmpdir/totalcounts $homedir/totalcounts

# Clean up
echo '######## Cleaning up #######'
#rm -r $tmpdir

# And that's it.
echo '######## Done! #######'
