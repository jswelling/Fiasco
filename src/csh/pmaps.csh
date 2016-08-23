#!/bin/csh -fx
# pmaps.csh
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
# *                                                          *
# *  Original programming by Bill Eddy                       *
# ************************************************************/
#

echo '#'`date`$0
echo '#$Id: pmaps.csh,v 1.8 2003/09/09 16:25:45 bakalj Exp $'

# What shall we use for input?
pushd $F_STAT_TPATH
echo "This is a placeholder" > Tmapdummy.mri
set tmaplist = ( Tmap*.mri )
rm Tmapdummy.mri
popd
pushd stat
echo "This is a placeholder" > Fmapdummy.mri
set fmaplist = ( Fmap*.mri )
rm Fmapdummy.mri
popd

if ( $#tmaplist > 1 ) then

# Compute P scores from Tmaps.  Remember to protect against the case
# where there are no samples in a given condition, so counts == 0.
# We'll also mask out places where the pooled Stdv is 0.0, since there 
# is presumably no data there.  This masking is accomplished by inserting
# Infs in the Pmap.
  cd $F_STAT_TPATH
  pooled_stdv.csh -out pooledStdv Stdv_*.mri
  foreach i ( $tmaplist )
    if ( $i != Tmapdummy.mri && $i !~ *_?.mri ) then
      set b = $i:t
      set j = $b:r
      set p = `echo $j | sed 's/Tmap/Pmap/g'`
      set xdim = `mri_printfield -field images.extent.x $j`
      set ydim = `mri_printfield -field images.extent.y $j`
      mri_copy_chunk -chunk counts1 -chunk_out images $j tmp_counts1
      mri_remap -order xyz tmp_counts1 
      mri_interp -d x -len $xdim -con tmp_counts1 tmp2_counts1
      mri_interp -d y -len $ydim -con tmp2_counts1 tmp3_counts1
      mri_copy_chunk -chunk counts2 -chunk_out images $j tmp_counts2
      mri_remap -order xyz tmp_counts2 
      mri_interp -d x -len $xdim -con tmp_counts2 tmp2_counts2
      mri_interp -d y -len $ydim -con tmp2_counts2 tmp3_counts2
      mri_type_convert -double $j tmp_dbl_tmap
      mri_rpn_math -out $p \
        '$1,$2,$3,+,dup,1,swap,0,==,if_keep,ct,nan,$4,0.0,==,if_keep' \
        tmp_dbl_tmap tmp3_counts1 tmp3_counts2 pooledStdv
      mri_rpn_math -out ${p}_bar \
        '$1,-1,*,$2,$3,+,dup,1,swap,0,==,if_keep,ct,nan,$4,0.0,==,if_keep' \
        tmp_dbl_tmap tmp3_counts1 tmp3_counts2 pooledStdv
      mri_destroy_dataset tmp_counts1
      mri_destroy_dataset tmp2_counts1
      mri_destroy_dataset tmp3_counts1
      mri_destroy_dataset tmp_counts2
      mri_destroy_dataset tmp2_counts2
      mri_destroy_dataset tmp3_counts2
      mri_destroy_dataset tmp_dbl_tmap
    endif
  end

#else if ( $#pmaplist > 1 ) then
## At this point we should compute P scores from Fmaps, but at the moment
## (Fiasco 5.1) no DOF information is being saved with the Fmaps so we
## can't convert them to P scores.
#
#  cd stat  
#  foreach i ( $fmaplist )
#    if ( $i != Fmapdummy.mri ) then
##     Do something wonderfully useful
#    endif
#  end
#
else

  echo "No Tmaps found; could not compute Pmaps. "

endif
