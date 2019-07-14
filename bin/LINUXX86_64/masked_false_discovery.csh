#! /bin/csh -ef
#
# Usage: masked_false_discovery.csh qval Mask Pmap
#
#   where qval is a false discovery rate, Mask is a mask file in Pgh MRI
#   format, and Pmap is a pmap in Pgh MRI format.
#
if ( $#argv != 3 ) then
  echo "usage: $0 qval Mask Pmap"
  echo "      qval is the Q threshold (e.g. 0.01)"
  echo "      Mask is a Pgh MRI file containing a mask of the brain"
  echo "      Pmap is a Pgh MRI file containing the P scores"
  exit -1
endif
set qval = $1
set mask = $2
set pmap = $3
mri_rpn_math -out mfd_tmp 'inf,$1,$2,if_keep' $pmap $mask
set thresh = `false_discovery.py $qval mfd_tmp`
echo "Q value of" $qval "implies P threshold of" $thresh
mri_rpn_math -out mfd_tmp2 '0,$x,$y,$z,$1,dup,'$thresh',>=,if_print_4' \
    mfd_tmp > mfd_tmp.t
set count = `wc -l mfd_tmp.t`
echo $count[1] "voxels satisfy the test:"
cat mfd_tmp.t | awk '{ print $4,"(",$1,",",$2,",",$3,")" } '
rm mfd_tmp.t
mri_destroy_dataset mfd_tmp
mri_destroy_dataset mfd_tmp2
