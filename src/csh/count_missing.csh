#!/bin/csh -ef
# count_missing.csh
#/************************************************************
# *                                                          *
# *  Permission is hereby granted to any individual or       *
# *  institution for use, copying, or redistribution of      *
# *  this code and associated documentation, provided        *
# *  that such code and documentation are not sold for       *
# *  profit and the following copyright notice is retained   *
# *  in the code and documentation:                          *
# *     Copyright (c) 1999 Department of Statistics,         *
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
# *  Original programming by Joel Welling                    *
# ************************************************************/
#
# $Id: count_missing.csh,v 1.6 2003/09/09 16:25:45 bakalj Exp $
#
set missing_test = `mri_printfield -field missing.size -nofail $1`
if ( $?F_TEMP ) then
  set tmpdir = $F_TEMP
else
  set tmpdir = ./
endif
if ( dummy$missing_test != dummy ) then
  set tmp1 = tmp1_count_missing_$$
  set tmp2 = tmp2_count_missing_$$
  set tmp3 = tmp3_count_missing_$$
  set tmp4 = tmp4_count_missing_$$
  mri_copy_chunk -chunk missing -chunk_out images $1 $tmpdir/$tmp1
  cd $tmpdir
  mri_type_convert -int $tmp1 $tmp2
  mri_destroy_dataset $tmp1
  mri_subsample -sum -d t -l 1 $tmp2 $tmp3
  mri_destroy_dataset $tmp2
  mri_rpn_math -out $tmp4 '$z,$1,1.0,if_print_2,0.0' $tmp3
  mri_destroy_dataset $tmp3
  mri_destroy_dataset $tmp4
else
  set zdim = `mri_printfield -field images.extent.z $1`
  set z = 0
  while ( $z < $zdim ) 
    echo $z 0
    @ z = $z + 1
  end
endif
