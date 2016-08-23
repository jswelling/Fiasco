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
# ************************************************************/
#
echo '#$Id: makeps.csh,v 1.24 2005/07/13 23:23:02 welling Exp $'
echo '#'`date`$0
cd stat
# convert means to postscript
if ( -e GrandMean.mri ) then
  setenv F_PS_TITLE2 "GrandMean"
  ${F_MRI_TO_PS} GrandMean
  mv GrandMean*.ps ../ps
  set backImg = "GrandMean"
else
  setenv F_PS_TITLE2 "Mean_1"
  ${F_MRI_TO_PS} Mean_1
  mv Mean_1*.ps ../ps
  set backImg = "Mean_1"
endif
echo '#Mean'
#
# convert std deviations to postscript
if ( -e GrandStdv.mri ) then
  setenv F_PS_TITLE2 "GrandStdv"
  ${F_MRI_TO_PS} GrandStdv
  mv GrandStdv*.ps ../ps
else if ( -e Stdv_1.mri ) then
  setenv F_PS_TITLE2 "Stdv_1"
  ${F_MRI_TO_PS} Stdv_1
  mv Stdv_1*.ps ../ps
endif
echo '#Stdv'

# If a false discovery rate (Q) limit has been specified,
# convert P maps to postscript using that limit.
if ( ${?F_LOWERQ} ) then
  echo "This is a placeholder" > Pmapdummy.mri
  foreach i (Pmap*.mri)
    if ( $i != Pmapdummy.mri && $i !~ *_?.mri ) then
      set b = $i:t
      set j = $b:r
      set pthresh = `false_discovery.py $F_LOWERQ $j`
      echo "# FDR of $F_LOWERQ implies thresholding $j at $pthresh"
      overlay.csh -inmap $j.mri -inimage ${backImg} \
             -headerout ${j}_q.mri \
             -lowthresh $pthresh -highthresh 2.0 \
             -mingray 64 -maxgray 191
      setenv F_PS_TITLE2 "$j at FDR = $F_LOWERQ"
      ${F_MRI_TO_PS} ${j}_q
      mv ${j}_q*.ps ../ps
    endif
  end
rm Pmapdummy.mri
endif

#
# Threshold Tmaps to select a given number of voxels if requested.
#
if ( ${?F_PICKN_COUNT} ) then
  set twotailed = 0
  set titlestr = "$F_PICKN_COUNT vox left tail"
  if ( ${?F_PICKN_TWOTAILED} ) then
    if ( ${F_PICKN_TWOTAILED} ) then
      set twotailed = 1
      set titlestr = "$F_PICKN_COUNT vox either tail"
    endif
  endif

  echo "This is a placeholder" > Tmapdummy.mri
  foreach i (Tmap*.mri)
    if ( $i != Tmapdummy.mri && $i !~ *_?.mri ) then
      set b = $i:t
      set j = $b:r
      if ( $twotailed ) then
        set val = `pick_n_voxels.py -T --twotails $F_PICKN_COUNT $j`
	set hival = `python -c "print repr(-1*${val})"`
        echo "# Picking $F_PICKN_COUNT implies thresholding $j at T<$val or T>$hival"
        overlay.csh -inmap $j.mri -inimage ${backImg} \
             -headerout ${j}_n.mri \
             -lowthresh $val -highthresh $hival \
             -mingray 64 -maxgray 191
      else
        set val = `pick_n_voxels.py -T $F_PICKN_COUNT $j`
        echo "# Picking $F_PICKN_COUNT implies thresholding $j at T<$val "
        overlay.csh -inmap $j.mri -inimage ${backImg} \
             -headerout ${j}_n.mri \
             -lowthresh $val -highthresh 10000 \
             -mingray 64 -maxgray 191
      endif
      setenv F_PS_TITLE2 "$j $titlestr"
      ${F_MRI_TO_PS} ${j}_n
      mv ${j}_n*.ps ../ps
    endif
  end
  rm Tmapdummy.mri

endif

#
# convert p maps to postscript
if ( ${?F_LOWERP} || ${?F_UPPERP} ) then

  if ( ${?F_LOWERP} ) then
    if ( ${?F_UPPERP} ) then
      set lowerp = $F_LOWERP
      set upperp = $F_UPPERP
      set titlestr = "$F_LOWERP,$F_UPPERP"
    else
      set lowerp = $F_LOWERP
      set upperp = 2.0
      set titlestr = "P < $F_LOWERP"
    endif
  else
      set lowerp = -1.0
      set upperp = $F_UPPERP
      set titlestr = "P > $F_UPPERP"
    else
  endif

  echo "This is a placeholder" > Pmapdummy.mri
  foreach i (Pmap*.mri)
    if ( $i != Pmapdummy.mri && $i !~ *_?.mri ) then
      set b = $i:t
      set j = $b:r
      overlay.csh -inmap $j.mri -inimage ${backImg} \
             -headerout ${j}_p.mri \
             -lowthresh $lowerp -highthresh $upperp \
             -mingray 64 -maxgray 191
      setenv F_PS_TITLE2 "$j $titlestr"
      ${F_MRI_TO_PS} ${j}_p
      mv ${j}_p*.ps ../ps
    endif
  end
  rm Pmapdummy.mri
echo '#Pmap'

endif

#
# convert t maps to postscript
if ( ${?F_LOWERT} || ${?F_UPPERT} ) then

  if ( ${?F_LOWERT} ) then
    if ( ${?F_UPPERT} ) then
      set lowert = $F_LOWERT
      set uppert = $F_UPPERT
      set titlestr = "$F_LOWERT < T < $F_UPPERT"
    else
      set lowert = $F_LOWERT
      set uppert = 1.0e10
      set titlestr = "T < $F_LOWERT"
    endif
  else
      set lowert = -1.0e10
      set uppert = $F_UPPERT
      set titlestr = "T > $F_UPPERT"
    else
  endif

  echo "This is a placeholder" > Tmapdummy.mri
  foreach i (Tmap*.mri)
    if ( $i != Tmapdummy.mri && $i !~ *_?.mri ) then
      set b = $i:t
      set j = $b:r
      overlay.csh -inmap $j.mri -inimage ${backImg} \
             -headerout ${j}_t.mri \
             -lowthresh $lowert -highthresh $uppert \
             -mingray 64 -maxgray 191
      setenv F_PS_TITLE2 "$j $titlestr"
      ${F_MRI_TO_PS} ${j}_t
      mv ${j}_t*.ps ../ps
    endif
  end
  rm Tmapdummy.mri
  echo '#Tmap'

endif

#
# convert fmaps to postscript
#
# Fmaps are positive definite, so lower bound never actually happens

if ( ${?F_UPPERF} ) then

  echo "This is a placeholder" > Fmapdummy.mri
  foreach i (Fmap*.mri)
    if ( $i != Fmapdummy.mri && $i !~ *_?.mri ) then
      if ( $i != Fmaps.mean.mri ) then
          set b = $i:t
          set j = $b:r
          overlay.csh -inmap $j.mri -inimage ${backImg} \
             -headerout ${j}_f.mri \
             -lowthresh "-10.0" -highthresh $F_UPPERF \
             -mingray 64 -maxgray 191
          setenv F_PS_TITLE2 "$j  F > $F_UPPERF"
          ${F_MRI_TO_PS} ${j}_f
          mv ${j}_f*.ps ../ps
      endif
    endif
  end
  rm Fmapdummy.mri
  echo '#Fmap'

endif

#
# FIX: THE FOLLOWING CODE NEEDS SOME WORK;
#	JUST EXIT FOR NOW
exit

cd ../par
set count = 0
foreach coef (Intercept Slope)
  foreach value (Mean Stdv)
    mri_subset -d t -l 1 -s $count detpar ${coef}_${value}
    setenv F_PS_TITLE2 "Coeficient $coef $value"
    ${F_MRI_TO_PS} ${coef}_${value}
    mv ${coef}_${value}*.ps ../ps
    rm ${coef}_${value}.mri ${coef}_${value}.dat
    @ count = $count + 1
  end
end
echo '#Detrend'



