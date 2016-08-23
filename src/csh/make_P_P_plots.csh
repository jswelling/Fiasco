#! /bin/csh -ef
# make_P_P_plots.csh
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
# *  Original programming by Bill Eddy                       *
# ************************************************************/
#
echo '#'`date` $0
echo '# $Id: make_P_P_plots.csh,v 1.4 2004/02/05 20:02:07 welling Exp $'
#
# This script generates arrays of P-P or Q-P plots
# usage: make_P_P_plots.csh Pfile1 Pfile2 Pfile3 ...
#

# Default range
set range = 1000

# Parse command line options
set plotout = "plot_P_P.ps"
set fracfit = 2

set errflag = 0
while ( $#argv >= 1 && dummy$argv[1] =~ dummy-* )
  switch ( $1 )
    case '-out':
      if ( $#argv < 2 ) then
        set errflag = 1
      else
        set plotout = $2
        shift
      endif
      breaksw
    case '-fitdenom':
      if ( $#argv < 2 ) then
        set errflag = 1
      else
        set fracfit = $2
        shift
      endif
      breaksw
    case '-help':
      if ( $#argv > 1 ) then
        scripthelp $0 $argv[2]
      else
        scripthelp $0
      endif
      exit 0
      breaksw
    default:
      set errflag = 1
      breaksw
  endsw
  shift
end
if ( $#argv < 1 ) then
  set errflag = 1
endif
if ($errflag != 0 || $#argv < 1 ) then
  scripthelp $0 usage
  exit -1
endif

#
# Make some scratch space, and go there.
#
if ( ${?F_TEMP} ) then
  set tmpdir = ${F_TEMP}/tmp_make_P_P_plots_$$
else
  set tmpdir = ./tmp_make_P_P_plots_$$
endif
if (! -e $tmpdir) mkdir $tmpdir
set homedir = `pwd`
cd $tmpdir
mkdir tfiles

#
# Extract all the needed text files
#
echo "# Generating P score text files"
set quoted_names = ()
foreach pfile ( $* )
  if ( ${pfile} =~ /* ) then
    set input = ${pfile}
  else
    set input = $homedir/${pfile}
  endif
  set dsname = ${pfile:r}
  set safename = `echo $dsname | sed 's%/%_%g'`
  set ofile = tfiles/${safename}_pscores.t
  set quoted_names = ( $quoted_names '"'$safename'",' )
  echo "# Processing $dsname"
  mri_copy_dataset $input in
  set xdim = `mri_printfield -fld images.extent.x in`
  set ydim = `mri_printfield -fld images.extent.y in`
  set zdim = `mri_printfield -fld images.extent.z in`
  @ qdim = $xdim * $ydim * $zdim
  mri_remap -order vq -len 1:$qdim in 
  mri_sort -ascending in sorted
  set n_valid = \
    `mri_rpn_math -out junk '0,$1,dup,is_finite,if_print_1' sorted | wc -l `
  if ( $n_valid != $qdim ) then
    mri_subset -d q -len $n_valid -shift 0 sorted valid
    set qdim = $n_valid
  else
    mri_copy_dataset sorted valid
  endif
  mri_rpn_math -out x '$q,1,+' valid
  mri_rpn_math -out junk '0,$q,$1,1,if_print_2' valid > $ofile
  @ q_fitfrac = $n_valid / $fracfit
  @ q_low = ( $n_valid - $q_fitfrac ) / 2
  mri_subset -d q -len $q_fitfrac -shift $q_low sorted sorted_middle
  mri_subset -d q -len $q_fitfrac -shift $q_low x x_middle
  mri_remap -order tqz -len ${q_fitfrac}:1:1 sorted_middle
  mri_remap -order tz -len ${q_fitfrac}:1 x_middle
  mri_glm -est glm_est:images sorted_middle:images x_middle:images
  set intercept = `mri_rpn_math -out junk '0,$1,$v,0,==,if_print_1' glm_est`
  set slope = `mri_rpn_math -out junk '0,$1,$v,1,==,if_print_1' glm_est`
  set fitfile = tfiles/${safename}_fit.t
  echo $intercept $slope > $fitfile
end

#
# Here is the R script we'll use
#
cat <<SCRIPT_EOF > plot_stuff.R
namelist <- c(${quoted_names})

postscript(horizontal=T,file="test.ps")
op <- par(mfrow=c(2,3),pty="s");

for (name in namelist) {
    data <- read.table(paste(name,"_pscores.t",sep=""))
    fit <- read.table(paste(name,"_fit.t",sep=""))
    plot(data,pch='.',main=name,ylim=c(0.0,1.0),
	xlab="",ylab="")
    abline(fit[1,1],fit[1,2])
}

print("All done!");
par(op);

SCRIPT_EOF

#
# Generate the plot .ps files
#
cd tfiles
if ( ${?SPLUS} ) then
  set rcommand = $SPLUS
else 
  set rcommand = "R --no-save"
endif
echo "# Generating plot"
${rcommand} < ../plot_stuff.R >& ../R.log
cp test.ps $homedir/$plotout
echo "# Done."

#clean up
cd $homedir
#rm -r $tmpdir

