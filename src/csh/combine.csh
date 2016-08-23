#! /bin/csh -ef
# combine.csh
#/************************************************************
# *                                                          *
# *  Permission is hereby granted to any individual or       *
# *  institution for use, copying, or redistribution of      *
# *  this code and associated documentation, provided        *
# *  that such code and documentation are not sold for       *
# *  profit and the following copyright notice is retained   *
# *  in the code and documentation:                          *
# *     Copyright (c) 1998,1999 Department of Statistics,    *
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
# *  Original programming by Nicole Lazar, Bill Eddy,        *
# *  and Joel Welling                                        *
# *  Adapted for use with Stouffer method by Jennifer Bakal  *
# ************************************************************/
#
# This script takes a collection of Tmap files as input and 
# combines them using Fisher's method. Note that the files are
# assumed to be spatially aligned!
#
echo '#'`date` $0
echo '# $Id: combine.csh,v 1.5 2004/12/09 00:22:09 welling Exp $'

# Check for help command
if ( $#argv >= 1 ) then
  if ( dummy$argv[1] == "dummy-help" ) then
    if ( $#argv >= 2 ) then
      scripthelp $0 $argv[2]
    else
      scripthelp $0
    endif
    exit
  endif
endif

# Parse command line options
set args = `getopt o:m:hf $*`
set ofile = ""
set help = ""
set fold = 0
set method = ""
while ($#args > 1) 
  switch ( $args[1] )
    case '-o' : 
      set ofile = $args[2] ; shift args; shift args; breaksw;
    case '-m' :
      set method = $args[2]; shift args; shift args; breaksw;
    case '-h' :
      set help = 1 ; shift args; breaksw;
    case '-f' :
      set fold = 1 ; shift args; breaksw;
  endsw
  if (junk${args[1]} == 'junk--') then
    shift args;
    break;
  endif
end

if ($#args < 1 || dummy$help != dummy ) then
  scripthelp $0
  exit -1
endif

if ( dummy$ofile == dummy ) then
  echo $0:t ': required output file omitted; use -o outfile'
  exit -1
endif

if (dummy$method == dummy ) then
  set method = "fisher";
  echo "No method specified. Defaulting to Fisher method."
else
  if ( $fold ) then
    echo "Combining using folded $method method."
  else
    echo "Combining using $method method."
  endif
endif

# Initialize the output with zeros, making sure to use double precision math.
# From here on in, we need to be careful not to let the type of tmp_${method}_sum
# accidentally get changed to single precision.
mri_copy_chunk -chunk images -replace $args[1] tmp_${method}_tmp
mri_type_convert -double tmp_${method}_tmp tmp_${method}_1
mri_destroy_dataset tmp_${method}_tmp
mri_rpn_math -out tmp_${method}_sum '0.0' tmp_${method}_1

# Slice size info
set xdim = `mri_printfield -field images.extent.x $args[1]`
set ydim = `mri_printfield -field images.extent.y $args[1]`

switch ($method)

  case 'fisher' :

    @ number = 0
    foreach i ($args)
	@ number = $number + 1

        echo "Preparing " $i
	mri_copy_chunk -chunk counts1 -chunk_out images -replace $i tmp_${method}_counts1
	mri_copy_chunk -chunk counts2 -chunk_out images -replace $i tmp_${method}_counts2 
	mri_rpn_math -out tmp_${method}_count '$1,$2,+,2,-' \
	    tmp_${method}_counts1 tmp_${method}_counts2
	mri_destroy_dataset tmp_${method}_counts1
	mri_destroy_dataset tmp_${method}_counts2
	mri_remap -order xyzt tmp_${method}_count
	set cnt_xdim = \
          `mri_printfield -field images.extent.x tmp_${method}_count`
	set cnt_ydim = \
          `mri_printfield -field images.extent.y tmp_${method}_count`
        if ( $cnt_xdim != $xdim ) then
          mri_interp -con -d x -len $xdim \
            tmp_${method}_count tmp_${method}_count_x
        else
          mri_copy_dataset tmp_${method}_count tmp_${method}_count_x
        endif
        if ( $cnt_ydim != $ydim ) then
          mri_interp -con -d y -len $ydim \
            tmp_${method}_count_x tmp_${method}_count_xy
        else
          mri_copy_dataset tmp_${method}_count_x tmp_${method}_count_xy
        endif

# Transform the t values to -2 log p values being careful about 0 and 1.
# Note that we make sure there is always at least 1 dof, to avoid confusing
# the cumulative T computation.
	mri_rpn_math '1,$1,$2,1,max,ct,-,dup,1.0,swap,1.0,<,if_keep,ln,-2,*' \
		-out tmp_${method}_1 $i tmp_${method}_count_xy
# add them up
	mri_rpn_math -out tmp_${method}_newsum '$1,$2,+' tmp_${method}_sum tmp_${method}_1
	mri_copy_dataset tmp_${method}_newsum tmp_${method}_sum
	mri_destroy_dataset tmp_${method}_newsum
    end

# compare against the chi-square distribution with 2*k df
    echo "Assembling " $ofile
    @ number = 2 * $number
    if ( $fold ) then
      mri_rpn_math -out $ofile '$1,'$number',fcchisqr' tmp_${method}_sum
    else
      mri_rpn_math -out $ofile '$1,'$number',cchisqr' tmp_${method}_sum
    endif
    breaksw

  case 'stouffer' :
    @ number = 0
    foreach i ($args)
	@ number = $number + 1

        echo "Preparing " $i
	mri_copy_chunk -chunk counts1 -chunk_out images -replace $i tmp_${method}_counts1
	mri_copy_chunk -chunk counts2 -chunk_out images -replace $i tmp_${method}_counts2 
	mri_rpn_math -out tmp_${method}_count '$1,$2,+,2,-,0,max' \
	    tmp_${method}_counts1 tmp_${method}_counts2
	mri_destroy_dataset tmp_${method}_counts1
	mri_destroy_dataset tmp_${method}_counts2
	mri_remap -order xyzt tmp_${method}_count
	set cnt_xdim = \
          `mri_printfield -field images.extent.x tmp_${method}_count`
	set cnt_ydim = \
          `mri_printfield -field images.extent.y tmp_${method}_count`
        if ( $cnt_xdim != $xdim ) then
          mri_interp -con -d x -len $xdim \
            tmp_${method}_count tmp_${method}_count_x
        else
          mri_copy_dataset tmp_${method}_count tmp_${method}_count_x
        endif
        if ( $cnt_ydim != $ydim ) then
          mri_interp -con -d y -len $ydim \
            tmp_${method}_count_x tmp_${method}_count_xy
        else
          mri_copy_dataset tmp_${method}_count_x tmp_${method}_count_xy
        endif

# Transform the t values to phiinv(1-p) values being careful about 0 and 1
# We can always do this using 'folded' P's to preserve accuracy.  Note
# that we set the minimum number of counts to 1 to avoid a problem with
# calculating the cumulative T on 0 degrees of freedom.
	mri_rpn_math '$1,$2,1,max,fct,-1,*,0.0,1.0,inv_fcnormal,-1,*' \
		-out tmp_${method}_1 $i tmp_${method}_count_xy
# add them up
	mri_rpn_math -out tmp_${method}_newsum '$1,$2,+' tmp_${method}_sum tmp_${method}_1
	mri_copy_dataset tmp_${method}_newsum tmp_${method}_sum
	mri_destroy_dataset tmp_${method}_newsum
    end

# divide by square root of k
    echo "Assembling " $ofile
    mri_rpn_math -out tmp_${method}_sum_div_sqrtk '$1,'$number',sqrt,/' tmp_${method}_sum

# compare against the normal distribution 
    if ( $fold ) then
      mri_rpn_math -out $ofile '$1,0,1,fcnormal' tmp_${method}_sum_div_sqrtk
    else
      mri_rpn_math -out $ofile '$1,0,1,cnormal' tmp_${method}_sum_div_sqrtk
    endif

    mri_destroy_dataset tmp_${method}_sum_div_sqrtk

    breaksw
    
  case 'binomial' :

    mri_copy_dataset tmp_${method}_sum tmp_${method}_totcount
    foreach i ($args)

        echo "Preparing " $i
	mri_copy_chunk -chunk counts1 -chunk_out images -replace \
          $i tmp_${method}_counts1
	mri_copy_chunk -chunk counts2 -chunk_out images -replace \
          $i tmp_${method}_counts2 
	mri_rpn_math -out tmp_${method}_count '$1,$2,+' \
	  tmp_${method}_counts1 tmp_${method}_counts2
	mri_destroy_dataset tmp_${method}_counts1
	mri_destroy_dataset tmp_${method}_counts2
	mri_remap -order xyzt tmp_${method}_count
	set cnt_xdim = \
          `mri_printfield -field images.extent.x tmp_${method}_count`
	set cnt_ydim = \
          `mri_printfield -field images.extent.y tmp_${method}_count`
        if ( $cnt_xdim != $xdim ) then
          mri_interp -con -d x -len $xdim \
            tmp_${method}_count tmp_${method}_count_x
        else
          mri_copy_dataset tmp_${method}_count tmp_${method}_count_x
        endif
        if ( $cnt_ydim != $ydim ) then
          mri_interp -con -d y -len $ydim \
            tmp_${method}_count_x tmp_${method}_count_xy
        else
          mri_copy_dataset tmp_${method}_count_x tmp_${method}_count_xy
        endif
	mri_rpn_math '$1,0.0,<,$2,0,!=,*' -out tmp_${method}_1 \
          $i tmp_${method}_count_xy

	mri_rpn_math -out tmp_${method}_newsum '$1,$2,+' \
          tmp_${method}_sum tmp_${method}_1 
	mri_copy_dataset tmp_${method}_newsum tmp_${method}_sum
	mri_destroy_dataset tmp_${method}_newsum

	mri_rpn_math -out tmp_${method}_newtotcount \
          '$1,$2,0,!=,$3,0.0,!=,*,$3,is_finite,*,+' \
	  tmp_${method}_totcount tmp_${method}_count_xy $i
	mri_copy_dataset tmp_${method}_newtotcount tmp_${method}_totcount
	mri_destroy_dataset tmp_${method}_newtotcount

    end

#   Calculate a P score based on a binomial distribution
    echo "Assembling " $ofile
    if ( $fold ) then
      mri_rpn_math -out $ofile \
        '$1,$2,dup,1.0,swap,0.0,==,if_keep,0.5,fcbinom' \
        tmp_${method}_sum tmp_${method}_totcount
    else
      mri_rpn_math -out $ofile \
        '$1,$2,dup,1.0,swap,0.0,==,if_keep,0.5,cbinom' \
	tmp_${method}_sum tmp_${method}_totcount
    endif

    mri_destroy_dataset tmp_${method}_totcount

    breaksw
    
  endsw


mri_destroy_dataset tmp_${method}_1
mri_destroy_dataset tmp_${method}_sum
mri_destroy_dataset tmp_${method}_count
mri_destroy_dataset tmp_${method}_count_x
mri_destroy_dataset tmp_${method}_count_xy
