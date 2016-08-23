#!/bin/csh -efx
# epi.printpar.csh
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
# *  Original programming by Audris Mockus and Bill Eddy     *
# ************************************************************/
#

echo '#'`date` $0
echo '#$Id: epi.printpar.csh,v 1.24 2006/04/07 17:54:20 welling Exp $'

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

#Truncate the run name if it is too long
set shortname = `echo $F_HEADER | cut -c 1-20`
if ( "$shortname" != "$F_HEADER" ) then
  set usename = "$shortname..."
else
  set usename = "$F_HEADER"
endif

# We will need some scratch files
set tmp1 = $F_TEMP/tmp1_printpar3d_$$.t
set tmp2 = $F_TEMP/tmp2_printpar3d_$$.t
set tmp3 = $F_TEMP/tmp3_printpar3d_$$.t
set tmp4 = $F_TEMP/tmp4_printpar3d_$$.t
set tmp5 = $F_TEMP/tmp5_printpar3d_$$.t
set tmp6 = $F_TEMP/tmp6_printpar3d_$$.t
set tmp7 = $F_TEMP/tmp7_printpar3d_$$.t
set tmp8 = $F_TEMP/tmp8_printpar3d_$$.t
set tmp9 = $F_TEMP/tmp9_printpar3d_$$.t
set tmp10 = $F_TEMP/tmp10_printpar3d_$$.t

foreach fname ( $tmp1 $tmp2 $tmp3 $tmp4 $tmp5 $tmp6 $tmp7 $tmp8 $tmp9 $tmp10 )
  if ( -e $fname ) rm $fname
end

# this part plots the baseline, mean and phase adjustments

set meancparms = `most_recent_parms.py meanc`
if ( dummy${meancparms} != dummy ) then
  set domeanc = T
  egrep -v '#' $meancparms > $tmp5
else
  set domeanc = F
endif

egrep -v '#' `most_recent_parms.py -f baseline` > $tmp6
egrep -v '#' `most_recent_parms.py -f baseline_raw` > $tmp7

set deghostparms = `most_recent_parms.py deghost`
if ( dummy${deghostparms} != dummy ) then
  set dodeghost = T
  egrep -v '#' $deghostparms > $tmp8
  egrep -v '#' `most_recent_parms.py -f deghost_raw` > $tmp9
else
  set dodeghost = F
endif

cat $FIASCO/../../splus/script/epiprintpar.S | sed "s/nslice/$F_NSLICE/g" | \
    sed "s/name/$usename/g" | sed "s/inc/$F_PRINTPAR_MAXPLOT/g" | \
    sed "s/nimage/$F_NIMAGE/g" | sed "s/dodeghost/$dodeghost/g" | \
    sed "s*tmpdir*$F_TEMP*g" | \
    sed "s*meancprm*$tmp5*g" | sed "s/domeanc/$domeanc/g" | \
    sed "s*baseadjrw*$tmp7*g" | \
    sed "s*baseadjsm*$tmp6*g" | \
    sed "s*deghostrw*$tmp9*g" | \
    sed "s*deghostsm*$tmp8*g" | \
    $SPLUS>>& out/S.output
if (! -d $F_PRINTPAR_SUBDIR) mkdir $F_PRINTPAR_SUBDIR
if( -e test.ps) mv test.ps $F_PRINTPAR_SUBDIR/parms.ps

# This part plots the 2D registration parameters
set regpar = `most_recent_parms.py estireg`
if ( dummy${regpar} != dummy ) then

  egrep -v '#' $regpar > $tmp1
  egrep -v '#' `most_recent_parms.py -f parsm` > $tmp2
  egrep -v '#' `most_recent_parms.py -f displace` > $tmp3
  egrep -v '#' `most_recent_parms.py -f smdisplace` > $tmp4

  cat $FIASCO/../../splus/script/printregist.S | sed "s/nslice/$F_NSLICE/g" | \
    sed "s/name/$usename/g" | sed "s/inc/$F_PRINTPAR_MAXPLOT/g" | \
    sed "s*smdisplace*${tmp4}*g" | \
    sed "s*dixplace*${tmp3}*g" | \
    sed "s*sxreg*${tmp2}*g" | \
    sed "s*rreg*${tmp1}*g" | \
    sed "s*nimage*$F_NIMAGE*g" | \
    sed "s*tmpdir*$F_TEMP*g" | \
    $SPLUS>>& out/S.output
  if( -e test.ps) mv test.ps $F_PRINTPAR_SUBDIR/regist.ps

endif

# This part plots the 3D registration parameters
set meanc3dpar = `most_recent_parms.py meanc3d`
if ( dummy${meanc3dpar} != dummy ) then
  set domeanc3d = T
  egrep -v '#' $meanc3dpar > $tmp8
else
  set domeanc3d = F
  set meanc3dpar = "none"
endif
set smreg3dpar = `most_recent_parms.py parsm3d`
if ( dummy${smreg3dpar} != dummy ) then
  set dosmreg3d = T
  egrep -v '#' $smreg3dpar > $tmp3
else
  set dosmreg3d = F
  set smreg3dpar = "none"
endif
set rdisp3dpar = `most_recent_parms.py displace3d_raw`
if ( dummy${rdisp3dpar} != dummy ) then
  set dordisp3d = T
  egrep -v '#' $rdisp3dpar > $tmp2
else
  set dordisp3d = F
  set rdisp3dpar = "none"
endif
set sdisp3dpar = `most_recent_parms.py displace3d`
if ( dummy${sdisp3dpar} != dummy ) then
  set dosdisp3d = T
  egrep -v '#' $sdisp3dpar > $tmp4
else
  set dosdisp3d = F
  set sdisp3dpar = "none"
endif
set reg3dpar = `most_recent_parms.py estireg3d`
if ( dummy${reg3dpar} != dummy ) then

  egrep -v '#' $reg3dpar > $tmp1

  cat $FIASCO/../../splus/script/printregist3d.S | \
    sed "s/nslice/$F_NSLICE/g" | \
    sed "s*name*$usename*g" | sed "s/inc/$F_PRINTPAR_MAXPLOT/g" | \
    sed "s/nimage/$F_NIMAGE/g" | sed "s*tmpdir*$F_TEMP*g" | \
    sed "s*rreg*$tmp1*g" | \
    sed "s*meanc3d*$tmp8*g" | sed "s/domnc3d/${domeanc3d}/g" | \
    sed "s*sxreg*$tmp3*g" | sed "s/dosmreg3d/${dosmreg3d}/g" | \
    sed "s*dixplace*$tmp2*g" | sed "s/dordisp3d/${dordisp3d}/g" | \
    sed "s*smdisplace*$tmp4*g" | sed "s/dosdisp3d/${dosdisp3d}/g" | \
    $SPLUS>>& out/S.output
  if( -e test.ps) mv test.ps $F_PRINTPAR_SUBDIR/regist3d.ps

endif

foreach fname ( $tmp1 $tmp2 $tmp3 $tmp4 $tmp5 $tmp6 $tmp7 $tmp8 $tmp9 $tmp10 )
  if ( -e $fname ) rm $fname
end

