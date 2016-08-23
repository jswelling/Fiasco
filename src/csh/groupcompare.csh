#! /bin/csh -ef
#
# This script does between-group comparisons.
#
echo '#'`date` $0
echo '# $Id: groupcompare.csh,v 1.24 2003/09/16 17:31:58 bakalj Exp $'

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
  set tmpdir = ${F_TEMP}/groupcompare_tmp_$$
else
  set tmpdir = ./groupcompare_tmp_$$
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
set group0_subjects = ( )
set group1_subjects = ( )

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
    if (! -e ${count_fname}.mri ) then
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
  if ( $groupid == 0 ) then
    set group0_subjects = ( $group0_subjects $subj )
    set lastgrp0 = $subj
  else if ( $groupid == 1 ) then
    set group1_subjects = ( $group1_subjects $subj )
    set lastgrp1 = $subj
  endif
end
@ nconds = $nconds - 1
set lastsubj = $subj
echo 'nsubjs is ' $nsubjs
echo 'nconds is ' $nconds
echo 'ngroups is ' $ngroups
echo 'group names: ' $groupnames
echo 'group ids: ' $sbj_grp_ids
echo 'subjects in group 0: ' $group0_subjects
echo 'subjects in group 1: ' $group1_subjects

if ( $nsubjs < 2 ) then
  echo "You need at least 2 subjects!"
  exit -1
endif

if ( $ngroups != 2 ) then
  echo "The subjects must be in exactly 2 groups!"
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
# - factor2 will be -1 for group 0, 1 for group 1
# - factor3 will be (factor1 * factor2)
# - normsubjfactorn will be 1 for subject n, 0 otherwise, de-meaned within 
#   group and condition, multiplied by factor1, with the last subject of 
#   each group dropped for orthogonality
#

echo '######## Generating Factors ########'

echo $sbj_grp_ids | mri_from_ascii -ord cs -l 1:${nsubjs} junk3
mri_interp -con -d c -len $nconds junk3 conds_x_subjects
mri_rpn_math -out group0_mask '$1,0,==' conds_x_subjects
mri_rpn_math -out group1_mask '$1,1,==' conds_x_subjects
mri_subsample -d s -l 1 -sum group0_mask group0_count
mri_subsample -d s -l 1 -sum group1_mask group1_count

mri_rpn_math -out rawfactor1 '$c,2,*,1,-' conds_x_subjects

mri_rpn_math -out rawfactor2 '$1,2,*,1,-' conds_x_subjects

mri_rpn_math -out rawfactor3 '$1,$2,*' rawfactor1 rawfactor2

set rawfactorlist = ( rawfactor1 rawfactor2 rawfactor3 )
set effectfactor = 1
set groupfactor = 2
set crossfactor = 3

set subjfactorlist = ( )
set which_subj = 0
foreach subj ( $group0_subjects )
  if ( $subj == $lastgrp0 ) then
    @ which_subj = $which_subj + 1  
    break
  endif
  mri_rpn_math -out normsubjfactor$which_subj \
    '$s,'${which_subj}',==,$2,1.0,$3,/,*,-,$4,*' \
    conds_x_subjects group0_mask group0_count rawfactor1
  set subjfactorlist = ( $subjfactorlist normsubjfactor${which_subj} )
  @ which_subj = $which_subj + 1  
end
foreach subj ( $group1_subjects )
  if ( $subj == $lastgrp1 ) then
    @ which_subj = $which_subj + 1  
    break
  endif
  mri_rpn_math -out normsubjfactor$which_subj \
    '$s,'${which_subj}',==,$2,1.0,$3,/,*,-,$4,*' \
    conds_x_subjects group1_mask group1_count rawfactor1
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
echo "The factor of interest is number " $crossfactor

# Remap the factors into the order needed for GLM
foreach fname ( $rawfactorlist $subjfactorlist )
  mri_remap -order t -length ${tdim} $fname
end

# Make the factor part of the operand strings for mri_glm
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
# stdv can be 0.0;  if that happens we want to substitute 0.0 for the scale
# value.  We sum up the counts as we go, for later use in computing the
# total degrees of freedom.  Also, smoothing may have produced numbers
# of counts less than 1; set those to zero.

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
    set tfile2 = junk3_${subj}_${cond}
    mri_rpn_math -out $tfile1 '0.0,$1,dup,1.0,<=,if_keep,sqrt' \
	$counts_fname $stdv_fname
    mri_remap -order vtxyz $tfile1 
    mri_delete_chunk -chunk missing $tfile1 
    set scale_args = ( $scale_args $tfile1 )
    mri_copy_dataset ${mean_fname} $tfile2
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
set xdim = `mri_printfield -fld images.extent.x rawy`
set ydim = `mri_printfield -fld images.extent.y rawy`
set zdim = `mri_printfield -fld images.extent.z rawy`
@ fakezdim = $xdim * $ydim * $zdim
mri_remap -order tz -length ${tdim}:${fakezdim} rawy
mri_remap -order tz -length ${tdim}:${fakezdim} scale

# Calculate the GLM over the raw (effect, group, and cross) factors
echo "######## Calculating the first GLM- this may take a while ########"
echo "The command being used is:"
echo mri_glm -v -ssqr -e raw_glm_results:images -out raw_residuals \
    -scale scale rawy $rawfactoropts
mri_glm -v -ssqr -e raw_glm_results:images -out raw_residuals \
    -scale scale rawy $rawfactoropts
echo '######## First GLM complete ########'

# Map this GLM results back to appropriate dimensions
set vdim = `mri_printfield -fld images.extent.v raw_glm_results`
mri_remap -order vxyz -length ${vdim}:${xdim}:${ydim}:${zdim} \
    raw_glm_results 

# Calculate the GLM over the subject factors
echo "######## Calculating the second GLM- this may take a while ########"
echo "The command being used is:"
echo mri_glm -v -ssqr -e subj_glm_results:images \
    -scale scale raw_residuals $subjfactoropts
mri_glm -v -ssqr -e subj_glm_results:images \
    -scale scale raw_residuals $subjfactoropts
echo '######## Second GLM complete ########'

# Map the GLM results back to appropriate dimensions
set vdim = `mri_printfield -fld images.extent.v subj_glm_results`
mri_remap -order vxyz -length ${vdim}:${xdim}:${ydim}:${zdim} \
    subj_glm_results 

# Grab interesting parts of the GLM.  The output of the model consists
# of (nfactors+1) estimated b's, followed by the associated summed squares.
# The cross term is factor3 above, so its estimate is the 3rd result
# value (at offset 2), since factor2 was not used.  Starting at offset
# (nfactors+2) are the sums of squares- SSTO at offset (nfactors+2),
# SSE at (nfactors+3), factor 0 SSR at (nfactors+4), and factor3 SSR
# at (nfactors+5).
mri_subset -d v -l 1 -s $effectfactor raw_glm_results effect_b
mri_subset -d v -l 1 -s $groupfactor raw_glm_results group_b
mri_subset -d v -l 1 -s $crossfactor raw_glm_results cross_b
@ ssto_offset = $nrawfactors + 2
mri_subset -d v -l 1 -s $ssto_offset raw_glm_results ssto
@ sse_offset = $nrawfactors + 3
mri_subset -d v -l 1 -s $sse_offset raw_glm_results sse
@ essqr_offset = $effectfactor + $nrawfactors + 3
mri_subset -d v -l 1 -s $essqr_offset raw_glm_results ssb
@ gssqr_offset = $groupfactor + $nrawfactors + 3
mri_subset -d v -l 1 -s $gssqr_offset raw_glm_results ssg
@ cssqr_offset = $crossfactor + $nrawfactors + 3
mri_subset -d v -l 1 -s $cssqr_offset raw_glm_results ssc
@ subj_ssto_offset = $nsubjfactors + 2
mri_subset -d v -l 1 -s $subj_ssto_offset subj_glm_results subj_ssto
@ subj_sse_offset = $nsubjfactors + 3
mri_subset -d v -l 1 -s $subj_sse_offset subj_glm_results subj_sse
mri_subset -d v -l 1 -s 0 raw_glm_results $homedir/mean

echo "######## Generating final statistics ########"

# The number of degrees of freedom in the sse is equal to the number of 
# input cases minus the number of fit parameters, or:
# (nsubjs*nconds) - (nrawfactors + 1).  
# We use this to calculate MSE (fixed effect model).
@ sse_dof  = ( $nsubjs * $nconds) - ( $nrawfactors + 1 )
mri_rpn_math -out mse '$1,'$sse_dof',/' sse

# We'll also need the degrees of freedom in the final sse (random
# effects model) after the second regression.
@ subj_sse_dof  = ( $nsubjs * $nconds) - ( $nfactors + 1 )

# sssg is summed squared error for subject within group; we calculate it 
# from SSTO = SSR + SSE (with a check for negative values due to floating
# point error). Degrees of freedom is the number of such columns in the 
# factor matrix, or:
# $nsubjs - 2 .
@ sssg_dof = $nsubjs - 2
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

mri_rpn_math -out $homedir/cross_f '$1,$2,/' $tmpdir/ssc $tmpdir/mssg
mri_type_convert -double $homedir/cross_f $tmpdir/tmp_doubles
mri_rpn_math -out $homedir/cross_p '$1,1,'$sssg_dof',cf' $tmpdir/tmp_doubles

# Let's keep totalcounts as well; it's useful for seeing where the
# samples were actually taken.
mri_copy_dataset $tmpdir/totalcounts $homedir/totalcounts

# Clean up
echo '######## Cleaning up #######'
#rm -r $tmpdir

# And that's it.
echo '######## Done! #######'
