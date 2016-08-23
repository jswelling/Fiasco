#! /bin/csh -ef
# plot_pixels.csh
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
echo '# $Id: plot_pixels.csh,v 1.4 2005/09/01 00:41:10 welling Exp $'
#
# This script generates huge arrays of plots of voxel time series.
# usage: plot_pixels.csh mrifile ixmin ixmax iymin iymax izmin izma
#

# Default range
set range = 1000

# Default label vals to none
set lblfile = ""

goto skipprintusage
printusage:
  echo "usage: $0 [-r range] [-l lbl] ixmin ixmax iymin iymax izmin izmax "
  echo "          mrifile1 mrifile2 ..."
  echo "  range sets the Y scale of each plot"
  echo "  lbl is a file of scalar values used to label the plots"
  echo "    (for example an F statistic)"
  echo "  ixmin to ixmax inclusive is the X range"
  echo "  iymin to iymax inclusive is the Y range"
  echo "  izmin to izmax inclusive is the Z range (usually izmin=izmax)"
  echo "  mrifile1, mrifile2, etc. are scalar files to plot"
  echo "     (all must have the same dimensions; only one file is required) "
  exit -1
skipprintusage:

# Parse command line options
set args = `getopt r:l: $*`
if ($#args < 8 ) then
  goto printusage
endif
while ($#args >= 8) 
  switch ( $args[1] )
    case '-r' : 
      set range = $args[2] ; shift args; shift args; breaksw;
    case '-l' :
      set lblfile = $args[2] ; shift args; shift args; breaksw;
    case '--' : 
      shift args; breaksw;
  endsw
  if (junk${args[1]} == 'junk--') then
    shift args;
    break;
  endif
end
if ($#args < 7) then
  goto printusage
endif

set ixmin = $args[1]
set ixmax = $args[2]
set iymin = $args[3]
set iymax = $args[4]
set izmin = $args[5]
set izmax = $args[6]
set mrifiles = ( $args[7-] )

if ( ( 1 + $ixmax - $ixmin ) % 8 ) then
  echo "ixmax - ixmin must be a multiple of 8!"
  exit -1
endif
if ( ( 1 + $iymax - $iymin ) % 8 ) then
  echo "iymax - iymin must be a multiple of 8!"
  exit -1
endif

set vdim = `mri_printfield -field images.extent.v -nofail ${mrifiles[1]}`
if ( foo$vdim == foo ) set vdim = 1
set xdim = `mri_printfield -field images.extent.x ${mrifiles[1]}`
set ydim = `mri_printfield -field images.extent.y ${mrifiles[1]}`
set zdim = `mri_printfield -field images.extent.z ${mrifiles[1]}`
set tdim = `mri_printfield -field images.extent.t -nofail ${mrifiles[1]}`
if ( foo$tdim == foo ) set tdim = 1

if ($vdim != 1) then
  echo "Dataset ${mrifiles} is not scalar!"
  exit -1
endif

# Get the root name, to construct output filename
set basename = $mrifiles[1]
if ( ${basename:e} == "mri") then
  set tailname = ${basename:t}
  set basename = $basename:r
else
  set basename = ${basename:t}
endif

#
# Make some scratch space, and go there.
#
if ( ${?F_TEMP} ) then
  set tmpdir = ${F_TEMP}/plot_pixels_tmp_$$
else
  set tmpdir = ./plot_pixels_tmp_$$
endif
if (! -e $tmpdir) mkdir $tmpdir
set homedir = `pwd`
set input = ()
foreach fname ( $mrifiles )
  if ( ${fname} =~ /* ) then
    set input = ( $input $fname )
  else
    set input = ( $input $homedir/$fname )
  endif
end
if ( dummy${lblfile} != dummy ) then
  if ( ${lblfile} =~ /* ) then
    set fulllblfile = $lblfile
  else
    set fulllblfile = $homedir/${lblfile}
  endif
endif
cd $tmpdir
mkdir tfiles

#
# Extract all the needed text files
#
@ cnt = 0
foreach fname ( $input )
  echo "Generating time series text files for ${fname:t}"
  @ iz = $izmin
  while ( $iz <= $izmax )
    mri_subset -d z -l 1 -s ${iz} $fname t1
    @ iy = $iymin
    while ( $iy <= $iymax )
      mri_subset -d y -l 1 -s ${iy} t1 t2
      @ ix = $ixmin
      while ( $ix <= $ixmax )
        echo $ix $iy $iz
        mri_subset -d x -l 1 -s ${ix} t2 t3
        set ofile = tfiles/data_${cnt}_${ix}_${iy}_${iz}.t
        echo '##Format: order:index_t type:raw names:(value)'  > $ofile
        mri_rpn_math '$t,$1,1,if_print_2' t3 >> $ofile
        echo '# end' >> $ofile
        @ ix = $ix + 1
      end
      @ iy = $iy + 1
    end
    @ iz = $iz + 1
  end
  @ cnt = $cnt + 1
end

if ( dummy${lblfile} != dummy ) then
  echo "Generating time series label data files"
  @ iz = $izmin
  while ( $iz <= $izmax )
    mri_subset -d z -l 1 -s ${iz} $fulllblfile t1
    @ iy = $iymin
    while ( $iy <= $iymax )
      mri_subset -d y -l 1 -s ${iy} t1 t2
      @ ix = $ixmin
      while ( $ix <= $ixmax )
        echo $ix $iy $iz
        mri_subset -d x -l 1 -s ${ix} t2 t3
        set ofile = tfiles/label_${ix}_${iy}_${iz}.t
        mri_rpn_math '$1,1,if_print_1' t3 >> $ofile
        @ ix = $ix + 1
      end
      @ iy = $iy + 1
    end
    @ iz = $iz + 1
  end  
endif

#
# Here is the R script we'll use
#
cat > plot_stuff_fixed_range.R << SCRIPT_EOF

xmin <- fake_xmin
xmax <- fake_xmax
ymin <- fake_ymin
ymax <- fake_ymax
z <- fake_z
range <- fake_range

ntot <- ( (xmax-xmin+1) * (ymax-ymin+3) );

postscript(horizontal=T,
        file=paste("data_",xmin,"-",xmax,"_",ymin,"-",ymax,"_",z,".ps",
        sep=""));
nf <- layout( matrix(1:ntot,ymax-ymin+1,xmax-xmin+1,byrow= TRUE) );
op <- par(xaxt="n",yaxt="n",mgp=c(1,0,0),mar=c(1.1,1.1,1.7,1.1));

for (y in ymin:ymax) {
  for (x in xmin:xmax) {
    print(paste(x,y,z));
    lblFName <- paste("label_",x,"_",y,"_",z,".t",sep="");
    if (!file.access(lblFName,4)) {
      label <- paste(x,y,z,":",sprintf("%.2f",read.table(lblFName)[1,1]));
    }
    else {
      label <- paste(x,y,z);
    }
    cnt <- 0
    while (!file.access(paste("data_",cnt,"_",xmin,"_",
                               ymin,"_",z,".t",sep=""))) {
      data <- read.table(paste("data_",cnt,"_",x,"_",y,"_",z,".t",sep=""));
      valmin <- mean(data[2])-range/2;
      valmax <- mean(data[2])+range/2;
      if ( cnt==0 ) {
        plot(data,type="l",xlab="",ylab="",main=label, ylim=c(valmin,valmax));
      }
      else {
        lines(data,lty=cnt+1);
      }
      cnt <- cnt+1;
    }
  }
}

print("All done!");
par(op);
SCRIPT_EOF

#
# Generate the plot .ps files
#
set plotdir = ${homedir}/plots
if (! -d $plotdir) mkdir $plotdir
cd tfiles
set scriptname = ../plot_stuff_fixed_range.R
@ iz = $izmin
while ( $iz <= $izmax )
  @ iy = $iymin
  while ( $iy <= $iymax )
    @ iytop = $iy + 7
    @ ix = $ixmin
    while ( $ix <= $ixmax )
      @ ixtop = $ix + 7
      cat $scriptname | sed "s/fake_xmin/$ix/g" | \
        sed "s/fake_xmax/$ixtop/g" | sed "s/fake_ymin/$iy/g" | \
        sed "s/fake_ymax/$iytop/g" | sed "s/fake_z/$iz/g" | \
        sed "s/fake_range/$range/g" | \
        R --no-save 
      mv  data_${ix}-${ixtop}_${iy}-${iytop}_${iz}.ps \
          $plotdir/${basename}_${ix}-${ixtop}_${iy}-${iytop}_${iz}.ps
      echo $ix $iy $iz
      @ ix = $ix + 8
    end
    @ iy = $iy + 8
  end
  @ iz = $iz + 1
end

echo "Done."

#clean up
cd $homedir
rm -r $tmpdir

