#!/bin/csh -ef
# makeps.csh
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
# *  Mods by Joel Welling                                    *
# ************************************************************/
#
echo '#$Id: cg_brain.makeps.csh,v 1.14 2007/06/15 18:11:14 welling Exp $'
echo '#'`date`$0
cd stat

# Generate Baselevel plots
foreach fname ( Baselevel*.mri )
  echo Converting $fname to PS slices
  setenv F_PS_TITLE2 ${fname:r}
  ${F_MRI_TO_PS} ${fname:r}
  mv ${fname:r}*.ps ../ps
end

# Generate Bayes.P plots
if ( -r Bayes.P.mri ) then
  set fname = Bayes.P.mri
  echo Converting $fname to PS slices
  setenv F_PS_TITLE2 "Prob. Activity in Any Condition"
  ${F_MRI_TO_PS} Bayes.P
  mv Bayes.P*.ps ../ps
endif

# Generate Bayes.PResp plots
foreach fname ( Bayes.PResp.*.mri )
  echo Converting $fname to PS slices
  set cond_h = ${fname:r}
  set cond = ${cond_h:e}
  setenv F_PS_TITLE2 "Prop. Signal Change, Condition $cond"
  ${F_MRI_TO_PS} ${fname:r}
  mv ${fname:r}*.ps ../ps
end


# Threshold Bayes.P to select a given number of voxels if requested.
if ( ${?F_PICKN_COUNT} ) then
  set twotailed = 0
  set titlestr = "$F_PICKN_COUNT vox left tail"
  if ( ${?F_PICKN_TWOTAILED} ) then
    if ( ${F_PICKN_TWOTAILED} ) then
      set twotailed = 1
      set titlestr = "$F_PICKN_COUNT vox either tail"
    endif
  endif

  foreach i ( Bayes.P.mri )
    set b = $i:t
    set j = $b:r
    if ( $twotailed ) then
      set val = `pick_n_voxels.py -P --twotails $F_PICKN_COUNT $j`
      set hival = `python -c "print repr(1.0-${val})"`
      echo "# Picking $F_PICKN_COUNT implies thresholding $j at P<$val or P>$hival"
      overlay.csh -inmap $j.mri -inimage Baselevel.Mean \
           -headerout ${j}_n.mri \
           -lowthresh $val -highthresh $hival 
           -mingray 64 -maxgray 191
    else
      set val = `pick_n_voxels.py -T $F_PICKN_COUNT $j`
      echo "# Picking $F_PICKN_COUNT implies thresholding $j at P<$val "
      overlay.csh -inmap $j.mri -inimage Baselevel.Mean \
           -headerout ${j}_n.mri \
           -lowthresh $val -highthresh 2.0 \
           -mingray 64 -maxgray 191
    endif
    setenv F_PS_TITLE2 "$j $titlestr"
    ${F_MRI_TO_PS} ${j}_n
    mv ${j}_n*.ps ../ps

  end

endif

# Make overlayed versions of Tmap files
if ( ${?F_LOWERT} || ${?F_UPPERT} ) then

  if ( ${?F_LOWERT} ) then
    if ( ${?F_UPPERT} ) then
      set lowert = $F_LOWERT
      set uppert = $F_UPPERT
      set titlestr = "$F_LOWERT < fakeT < $F_UPPERT"
    else
      set lowert = $F_LOWERT
      set uppert = 1.0e10
      set titlestr = "fakeT < $F_LOWERT"
    endif
  else
      set lowert = -1.0e10
      set uppert = $F_UPPERT
      set titlestr = "fakeT > $F_UPPERT"
    else
  endif

  echo "This is a placeholder" > dummy.Tmap.mri
  foreach fname (*.Tmap.mri)
    if ( $fname != dummy.Tmap.mri && $fname !~ *_?.mri ) then
      echo Generating overlays of $fname
      set b = $fname:t
      set j = $b:r
      overlay.csh -inmap $j.mri -inimage Baselevel.Mean \
         -headerout ${j}_t.mri \
         -lowthresh $lowert -highthresh $uppert \
         -mingray 64 -maxgray 191
      setenv F_PS_TITLE2 "$j $titlestr"
      ${F_MRI_TO_PS} ${j}_t
      mv ${j}_t*.ps ../ps
    endif
  end
  rm dummy.Tmap.mri

endif

# If a false discovery rate (Q) limit has been specified,
# convert P maps to postscript using that limit.
if ( ${?F_LOWERQ} ) then
  set tmpfile = makeps_tmp_$$
  mri_rpn_math -out $tmpfile '1,$1,-' Bayes.P
  set pthresh = `false_discovery.py $F_LOWERQ $tmpfile`
  echo "# FDR of $F_LOWERQ implies thresholding Bayes.P at 1.0-$pthresh"
  overlay.csh -inmap $tmpfile.mri -inimage Baselevel.Mean \
          -headerout Bayes.P_q.mri \
          -lowthresh $pthresh -highthresh 2.0 \
          -mingray 64 -maxgray 191
  setenv F_PS_TITLE2 "Bayes.P at FDR = $F_LOWERQ"
  ${F_MRI_TO_PS} Bayes.P_q
  mv Bayes.P_q*.ps ../ps
  mri_destroy_dataset $tmpfile
endif


