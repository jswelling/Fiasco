#! /bin/csh -exf
#
#  This script converts an AFNI text dump of a masked dataset
#  to a Pgh MRI format mask with value 1 in the corresponding
#  region and 0 elsewhere.
#
#  usage:
#   afni_mask_to_fiasco.csh afnimask.txt proto fiascomask
#
#   where: afnimask.txt is the text file from AFNI
#          proto.mri is a Pgh MRI file of the right "shape"
#          fiascomask.mri is the output Fiasco mask file
#
set xdim = `mri_printfield -field "images.extent.x" $2`
set ydim = `mri_printfield -field "images.extent.y" $2`
set zdim = `mri_printfield -field "images.extent.z" $2`
grep -v +++ $1 | awk -F '[ ,()]' '{ print $4,$5,$6, 1 }' | \
  mri_from_ascii -ord xyz -l ${xdim}:${ydim}:${zdim} -ind xyz $3
