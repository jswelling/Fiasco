#!/bin/csh -ef
# print.files.csh
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
# *  Original programming by Bill Eddy and Audris Mockus     *
# ************************************************************/
#
# first argument chooses first set of files
# summary parms reg mean stdv stat fmaps detrend
# second argument chooses last set of files
# third argument is output device
#
echo '#'`date` $0
echo '$Id: print.files.csh,v 1.21 2005/06/10 20:12:29 welling Exp $'
set output = $3
set done = $2
if (! -d print) mkdir print
cd ps

switch ($1)
  case default:
  case summary:
  echo summary
  if(-e summary.ps) then 
    set create_time = `ls -l summary.ps | awk '{print $6,$7,$8}'`
    cat summary.ps | sed "s/%%PRINTMSG%%/(Printed:`date`)/" \
      | sed "s/%%CREATIONMSG%%/(Created: $create_time)/" \
      > ../print/summary.ps
    rm summary.ps
    if ($output != "none" ) $output ../print/summary.ps
    if($done == summary)then
      breaksw
    endif
  endif

  case parms:
  echo parms
  if (-e parms.ps)then
    set create_time = `ls -l parms.ps | awk '{print $6,$7,$8}'`
    sed "s'Stoplight'Created $create_time, Printed `date`'" \
	$FIASCO/../../src/misc/sed.script > curr.script
    cat parms.ps | sed -f curr.script \
      > ../print/parms.ps
    rm parms.ps
    if ($output != "none" ) $output ../print/parms.ps
    rm curr.script
    if($done == parms)then
      breaksw
    endif
  endif

  case reg:
  echo reg
  if(-e regist.ps) then
    set create_time = `ls -l regist.ps | awk '{print $6,$7,$8}'`
    sed "s'Stoplight'Created $create_time, Printed `date`'" $FIASCO/../../src/misc/sed.script > curr.script
    cat regist.ps | sed -f curr.script \
      > ../print/regist.ps
    rm regist.ps
    if ($output != "none" ) $output ../print/regist.ps
    rm curr.script
  endif
  if(-e regist3d.ps) then
    set create_time = `ls -l regist3d.ps | awk '{print $6,$7,$8}'`
    sed "s'Stoplight'Created $create_time, Printed `date`'" $FIASCO/../../src/misc/sed.script > curr.script
    cat regist3d.ps | sed -f curr.script \
      > ../print/regist3d.ps
    rm regist3d.ps
    if ($output != "none" ) $output ../print/regist3d.ps
    rm curr.script
  endif
  if($done == reg)then
    breaksw
  endif


  case pcaplots:
  echo pcaplots
  if(-e pcaplots.ps) then
    set create_time = `ls -l pcaplots.ps | awk '{print $6,$7,$8}'`
    sed "s'Stoplight'Created $create_time, Printed `date`'" $FIASCO/../../src/misc/sed.script > curr.script
    cat pcaplots.ps | sed -f curr.script \
      > ../print/pcaplots.ps
    rm pcaplots.ps
    if ($output != "none" ) $output ../print/pcaplots.ps
    rm curr.script
  endif
  if($done == pcaplots)then
    breaksw
  endif


  case Mean:
  echo Mean
	echo "This is a placeholder" > dummyMean.ps
	foreach i ( *Mean*.ps )
	  if ( $i != dummyMean.ps ) then

                    cat $i | sed "s/%%PRINTMSG%%/(Printed: `date`) show/" \
                      > ../print/$i

                    rm $i
		    if ($output != "none" ) $output ../print/$i
          endif
	end
	rm dummyMean.ps
        if($done == mean)then
 	  breaksw
	endif

  case Stdv:
  echo Stdv
	echo "This is a placeholder" > dummyStdv.ps
	foreach i ( *Stdv*.ps )
	  if ( $i != dummyStdv.ps ) then
                    cat $i | sed "s/%%PRINTMSG%%/(Printed: `date`) show/" \
                      > ../print/$i
                    rm $i
		    if ($output != "none" ) $output ../print/$i
          endif
	end
	rm dummyStdv.ps
        if($done == stdv)then
 	  breaksw
	endif

  case Pmap:
  echo Pmap
	echo "This is a placeholder" > dummyPmap.ps
	foreach i ( *Pmap*.ps )
	  if ( $i != dummyPmap.ps ) then
                    cat $i | sed "s/%%PRINTMSG%%/(Printed: `date`) show/" \
                      > ../print/$i
                    rm $i
		    if ($output != "none" ) $output ../print/$i
          endif
	end
	rm dummyPmap.ps
        if($done == pmap)then
 	  breaksw
	endif

  case Tmap:
  echo Tmap
	echo "This is a placeholder" > dummyTmap.ps
	foreach i ( *Tmap*.ps )
	  if ( $i != dummyTmap.ps ) then
                    cat $i | sed "s/%%PRINTMSG%%/(Printed: `date`) show/" \
                      > ../print/$i
                    rm $i
		    if ($output != "none" ) $output ../print/$i
          endif
	end
	rm dummyTmap.ps
        if($done == tmap)then
 	  breaksw
	endif

  case detrend:
  echo Detrend
	echo "This is a placeholder" > Interceptdummy.ps
	echo "This is a placeholder" > Slopedummy.ps
	foreach i ( Intercept*.ps Slope*.ps )
	  if ( ( $i != Interceptdummy.ps ) && ( $i != Slopedummy.ps ) ) then
                    cat $i | sed "s/%%PRINTMSG%%/(Printed: `date`) show/" \
                      > ../print/$i
                    rm $i
		    if ($output != "none" ) $output ../print/$i
          endif
	end
	rm Interceptdummy.ps Slopedummy.ps
        if($done == detrend)then
 	  breaksw
	endif

  case Fmaps:
  echo Fmaps
	echo "This is a placeholder" > Fmaps.dummy.ps
	foreach i ( *Fmap*.ps )
	  if ( $i != Fmaps.dummy.ps ) then
                    cat $i | sed "s/%%PRINTMSG%%/(Printed: `date`) show/" \
                      > ../print/$i
                    rm $i
		    if ($output != "none" ) $output ../print/$i
          endif
	end
	rm Fmaps.dummy.ps
        if($done == fmaps)then
 	  breaksw
	endif

  case Bayes:
  echo Bayes
	echo "This is a placeholder" > Bayes.dummy.ps
	foreach i ( Bayes.*.ps )
	  if ( $i != Bayes.dummy.ps ) then
                    cat $i | sed "s/%%PRINTMSG%%/(Printed: `date`) show/" \
                      > ../print/$i
                    rm $i
		    if ($output != "none" ) $output ../print/$i
          endif
	end
	rm Bayes.dummy.ps
        if($done == Bayes)then
 	  breaksw
	endif
  
  case pca:
  echo pca
	echo "This is a placeholder" > pca_dummy.ps
	foreach i ( pca_*.ps )
	  if ( $i != pca_dummy.ps ) then
                    cat $i | sed "s/%%PRINTMSG%%/(Printed: `date`) show/" \
                      > ../print/$i
                    rm $i
		    if ($output != "none" ) $output ../print/$i
          endif
	end
	rm pca_dummy.ps
        if($done == Bayes)then
 	  breaksw
	endif
  
  case Coreg:
  echo Coreg
	echo "This is a placeholder" > dummycoreg.ps
	foreach i ( *coreg*.ps )
	  if ( $i != dummycoreg.ps ) then

                    cat $i | sed "s/%%PRINTMSG%%/(Printed: `date`) show/" \
                      > ../print/$i

                    rm $i
		    if ($output != "none" ) $output ../print/$i
          endif
	end
	rm dummycoreg.ps
        if($done == coreg)then
 	  breaksw
	endif

endsw
