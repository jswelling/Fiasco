#!/bin/csh -ef
# pca_makeps.csh
#/************************************************************
# *                                                          *
# *  Permission is hereby granted to any individual or       *
# *  institution for use, copying, or redistribution of      *
# *  this code and associated documentation, provided        *
# *  that such code and documentation are not sold for       *
# *  profit and the following copyright notice is retained   *
# *  in the code and documentation:                          *
# *     Copyright (c) 2004 Department of Statistics,         *
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
# *  Mods by Joel Welling                                    *
# ************************************************************/
#
echo '#$Id: pca_makeps.csh,v 1.5 2005/07/13 23:23:34 welling Exp $'
echo '#'`date`$0

if ( ! -d ps ) mkdir ps
if ( ! -d out ) mkdir out
cd stat

# Generate GrandWhatever plots
foreach fname ( Grand*.mri )
  echo Converting $fname to PS slices
  setenv F_PS_TITLE2 ${fname:r}
  ${F_MRI_TO_PS} ${fname:r}
  mv ${fname:r}*.ps ../ps
end

# Generate PCA component plots
set i = 1
while ( $i <= $F_PCA_NCOMPONENTS )
  set fname = "pca_${i}"
  echo Converting $fname to PS slices
  setenv F_PS_TITLE2 "PCA Component $i"
  ${F_MRI_TO_PS} $fname
  mv ${fname}*.ps ../ps
  @ i = $i + 1
end

#
# Threshold PCA maps to select a given number of voxels if requested.
# This is always done two-tailed.
#
if ( ${?F_PICKN_COUNT} ) then
  set twotailed = 0
  set titlestr = "$F_PICKN_COUNT vox either tail"
  set twotailed = 1
  set backImg = "GrandMean"
  set count = 1

  while ( -e pca_${count}.mri )
    set i = pca_${count}.mri
    set b = $i:t
    set j = $b:r
    set val = `pick_n_voxels.py -T --twotails $F_PICKN_COUNT $j`
    set hival = `python -c "print repr(-1*${val})"`
    echo "# Picking $F_PICKN_COUNT implies thresholding $j at val<$val or val>$hival"
    overlay.csh -inmap $j.mri -inimage ${backImg} \
         -headerout ${j}_n.mri \
         -lowthresh $val -highthresh $hival \
         -mingray 64 -maxgray 191
    setenv F_PS_TITLE2 "$j $titlestr"
    ${F_MRI_TO_PS} ${j}_n
    mv ${j}_n*.ps ../ps
    @ count = $count + 1
  end

endif

#Quit if SPLUS does not point to a valid executable
if (! $?SPLUS) then
	echo 'SPLUS not defined so plots not generated\!'
	exit 0
endif
set splus_words = `echo $SPLUS`
set splus_exe = $splus_words[1]
if (! -x $splus_exe) then
	echo 'Cannot find Splus in ' $splus_exe ' so plots not generated\!'
	exit 0
endif

# We will need some scratch files
if ( ! -d $F_TEMP ) mkdir $F_TEMP
set tmp1 = $F_TEMP/tmp1_printpca_$$.t
set tmp2 = $F_TEMP/tmp2_printpca_$$.t
set tmp3 = $F_TEMP/tmp3_printpca_$$.t
set tmp4 = $F_TEMP/tmp4_printpca_$$.t
set tmp5 = $F_TEMP/tmp5_printpca_$$.t

# Extract the eigenvalues and left eigenvectors

mri_rpn_math '$a,$1,1,if_print_2' pca_eigenvals > $tmp1
mri_rpn_math '$a,$1,1,if_print_2' pca_F_termwise > $tmp2
mri_rpn_math '$a,$1,ln,10,ln,/,1,if_print_2' pca_P_termwise > $tmp3
mri_rpn_math '$a,$1,1,if_print_2' pca_mse_termwise > $tmp4
mri_rpn_math '$t,$1,1,if_print_2' pca_left_1 > $tmp5
set i = 2
while ( $i <= $F_PCA_NCOMPONENTS )
  mri_rpn_math '$t,$1,1,if_print_2' pca_left_$i >> $tmp5
  @ i = $i + 1
end

cat $FIASCO/../../splus/script/printpca.S | \
    sed "s/ncomponents/$F_PCA_NCOMPONENTS/g" | \
    sed "s*evalsFile*$tmp1*g" | sed "s*fsFile*$tmp2*g" | \
    sed "s*psFile*$tmp3*g" | sed "s*msesFile*$tmp4*g" | \
    sed "s*evecsFile*$tmp5*g" | sed "s*headerString*$F_DESCRIPTION*g" | \
    $SPLUS>>& ../out/S.output
if( -e test.ps) mv test.ps ../ps/pcaplots.ps


foreach fname ( $tmp1 $tmp2 $tmp3 $tmp4 $tmp5 )
  if ( -e $fname ) rm $fname
end

