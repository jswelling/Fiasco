#! /bin/csh -ef
#
# This script takes a series of Stdv files as input (each of which
# also contains a counts chunk) and produces a single pooled Stdv
# file as output.
#
echo '#'`date` $0
echo '# $Id: pooled_stdv.csh,v 1.7 2004/09/09 22:35:23 welling Exp $'

# Check for help command
if ( $#argv >= 1 ) then 
  if ( junk$argv[1] == "junk-help" ) then
    if ( $#argv >= 2 ) then
      scripthelp $0 $argv[2]
    else
      scripthelp $0
    endif
    exit
  endif
endif
if ( $#argv == 0 ) then 
  scripthelp $0
  exit
endif

# Parse arguments
set outfile = "PooledStdv"
set badarg = 0
while ( junk$argv[1] =~ junk-* ) 
  switch ( $argv[1] )
    case '-out' : 
      set outfile = $argv[2]; shift argv; shift argv; breaksw;
    default :
      set badarg = 1; shift argv; breaksw;
  endsw
  if ( $#argv == 0 ) break;
end
if ( $badarg != 0 ) then
  echo $0 ': invalid command line argument'
  scripthelp $0
  exit -1
endif
if ( $#argv == 0 ) then
  echo $0 ': no stdv files to pool!'
  scripthelp $0
  exit -1
endif
#
# Make some scratch space
#
if ( ${?F_TEMP} ) then
  set tmpdir = ${F_TEMP}/pooled_stdv_$$
else
  set tmpdir = ./pooled_stdv_$$
endif
if (! -e $tmpdir) mkdir $tmpdir
echo "Scratch directory will be " $tmpdir

# Do we have the DOF chunk?  We'll check the first and trust the rest.
set dofchunk = "`mri_printfield -fld dof -nofail $1`"
if ( dummy"${dofchunk}" == 'dummy[chunk]' ) then
  set has_dof = 1
else
  set has_dof = 0
endif

# Collect all stdvs into a single file, indexed over t dimension
# We take some steps to restore the original dimensionality of 
# the counts chunk.
set original_dims = `mri_printfield -field images.dimensions $1`
set original_cnt_dims = `mri_printfield -field counts.dimensions $1`
if ( $has_dof ) then
  set original_dof_dims = `mri_printfield -field dof.dimensions $1 -nofail`
else
  set original_dof_dims = $original_cnt_dims
endif
foreach fname ( $argv )
  mri_remap -chunk images -order xyzt $fname
  mri_remap -chunk counts -order xyzt $fname
  if ( $has_dof ) mri_remap -chunk dof -order xyzt $fname
end
set fnamelist = ( $argv )
set nTerms = $#fnamelist
mri_copy_dataset ${fnamelist[1]} $tmpdir/allstdvs
shift fnamelist
while ( $#fnamelist )
  if ( $#fnamelist > 20 ) then
    set these_args = ( $fnamelist[1-20] )
    set fnamelist = ( $fnamelist[21-] )
  else
    set these_args = ( $fnamelist )
    set fnamelist = ( )
  endif
  mri_paste -d t -out $tmpdir/tmp_paste $tmpdir/allstdvs $these_args
  mri_copy_dataset $tmpdir/tmp_paste $tmpdir/allstdvs
  mri_destroy_dataset $tmpdir/tmp_paste
end
foreach fname ( $argv )
  mri_remap -chunk images -order $original_dims $fname
  mri_remap -chunk counts -order $original_cnt_dims $fname
  if ( $has_dof ) mri_remap -chunk dof -order $original_dof_dims $fname
end

# Build the pooled stdv dataset
#
# numerator_terms is (n-1)s^2 for each of the component stdv's
#
# numerator is the sum of those terms
# denominator is (sum of n's) - N where N is the number of terms
# result is the square root of the ratio of these two.  Alternately
# we'll use the DOF information if it is given.
#
mri_copy_chunk -chunk counts -chunk_out images -replace \
    $tmpdir/allstdvs $tmpdir/allcounts
mri_permute -order ztxy $tmpdir/allstdvs $tmpdir/allstdvs_p
mri_permute -order ztxy $tmpdir/allcounts $tmpdir/allcounts_p 
mri_rpn_math -out $tmpdir/numerator_terms '$1,dup,*,$2,1,-,*' \
    $tmpdir/allstdvs_p $tmpdir/allcounts_p
mri_subsample -d t -l 1 -sum $tmpdir/numerator_terms $tmpdir/numerator_p
mri_subsample -d t -l 1 -sum $tmpdir/allcounts_p $tmpdir/totalcounts_p
if ( $has_dof ) then
  mri_copy_chunk -chunk dof -chunk_out images -replace \
    $tmpdir/allstdvs $tmpdir/alldofs 
  mri_permute -order ztxy $tmpdir/alldofs $tmpdir/alldofs_p 
  mri_subsample -d t -l 1 -sum $tmpdir/alldofs_p $tmpdir/totaldof_p
else
  mri_rpn_math -out $tmpdir/totaldof_p '$1,'${nTerms}',-' \
    $tmpdir/totalcounts_p
endif
mri_rpn_math -out $tmpdir/stdvs_p '$1,$2,/,sqrt' \
    $tmpdir/numerator_p $tmpdir/totaldof_p
mri_permute -order xyzt $tmpdir/stdvs_p $outfile 
mri_permute -order xyzt $tmpdir/totalcounts_p $tmpdir/totalcounts
mri_permute -order xyzt $tmpdir/totaldof_p $tmpdir/totaldof
mri_copy_chunk -chunk images -chunk_out counts -replace \
    $tmpdir/totalcounts $outfile
mri_copy_chunk -chunk images -chunk_out dof -replace \
    $tmpdir/totaldof $outfile 
mri_remap -chunk images -order $original_dims $outfile
mri_remap -chunk counts -order $original_cnt_dims $outfile
mri_remap -chunk dof -order $original_dof_dims $outfile

# Clean up
rm -r $tmpdir



