#!/bin/csh -efx
# ts.printpar.csh
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
echo '#$Id: ts.printpar.csh,v 1.18 2003/07/01 22:27:50 welling Exp $'

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

#this part plots the baseline, mean and phase adjustments
set meancparms = `most_recent_parms.py -f -t meanc`
set baselineparms = `most_recent_parms.py -f -t baseline`
set baselinerawparms = `most_recent_parms.py -f -t baseline_raw`
set deghostparms = `most_recent_parms.py -t deghost`
if ( dummy${deghostparms} != dummy ) then
  set dodeghost = T
  set deghostrawparms = `most_recent_parms.py -f -t deghost_raw`
else
  set dodeghost = F
  set deghostparms = "dummy"
  set deghostrawparms = "dummy"
endif

cat $FIASCO/../../splus/script/tsprintpar.S | sed "s/nslice/$F_NSLICE/g" | \
    sed "s/name/$usename/g" | sed "s/inc/$F_PRINTPAR_MAXPLOT/g" | \
    sed "s/nimage/$F_NIMAGE/g" | sed "s/dodeghost/$dodeghost/g" | \
    sed "s*tmpdir*$F_TEMP*g" | sed "s/meancprm/$meancparms/g" | \
    sed "s/baseadjrw/$baselinerawparms/g" | \
    sed "s/baseadjsm/$baselineparms/g" | \
    sed "s/deghostrw/$deghostrawparms/g" | \
    sed "s/deghostsm/$deghostparms/g" | \
    $SPLUS>>& out/S.output
if (! -d $F_PRINTPAR_SUBDIR) mkdir $F_PRINTPAR_SUBDIR
  if( -e test.ps) mv test.ps $F_PRINTPAR_SUBDIR/parms.ps

# This part plots the 2D registration parameters
set regpar = `most_recent_parms.py -t estireg`
if ( dummy${regpar} != dummy ) then

  set smregpar = `most_recent_parms.py -f -t parsm`
  set disppar = `most_recent_parms.py -f -t smdisplace`
  set disprawpar = `most_recent_parms.py -f -t displace`
  cat $FIASCO/../../splus/script/printregist.S | sed "s/nslice/$F_NSLICE/g" | \
    sed "s/name/$usename/g" | sed "s/inc/$F_PRINTPAR_MAXPLOT/g" | \
    sed "s*smdisplace*${disppar}*g" | \
    sed "s/dixplace/${disprawpar}/g" | \
    sed "s/sxreg/$smregpar/g" | sed "s/rreg/$regpar/g" | \
    sed "s/nimage/$F_NIMAGE/g" | sed "s*tmpdir*$F_TEMP*g" | \
    $SPLUS>>& out/S.output
  if( -e test.ps) mv test.ps $F_PRINTPAR_SUBDIR/regist.ps

endif

# This part plots the 3D registration parameters
if ( -r par/${F_EST3D_PARMS} ) then
  set tfile1 = $F_TEMP/tmp1_printpar3d_$$.t
  set tfile2 = $F_TEMP/tmp2_printpar3d_$$.t
  set tfile3 = $F_TEMP/tmp3_printpar3d_$$.t
  set tfile4 = $F_TEMP/tmp4_printpar3d_$$.t

  egrep -v '#' `most_recent_parms.py -f estireg3d` > $tfile1
  egrep -v '#' `most_recent_parms.py -f displace3d` > $tfile2
  egrep -v '#' `most_recent_parms.py -f parsm3d` > $tfile3
  egrep -v '#' `most_recent_parms.py -f smdisplace3d` > $tfile4

  cat $FIASCO/../../splus/script/printregist3d.S | \
    sed "s/nslice/$F_NSLICE/g" | \
    sed "s*name*$usename*g" | sed "s/inc/$F_PRINTPAR_MAXPLOT/g" | \
    sed "s*rreg*$tfile1*g" | sed "s*dixplace*$tfile2*g" | \
    sed "s*sxreg*$tfile3*g" | sed "s*smdisplace*$tfile4*g" | \
    sed "s/nimage/$F_NIMAGE/g" | sed "s*tmpdir*$F_TEMP*g" | \
    $SPLUS>>& out/S.output
  if( -e test.ps) mv test.ps $F_PRINTPAR_SUBDIR/regist3d.ps

  rm $tfile1 $tfile2 $tfile3 $tfile4
endif

