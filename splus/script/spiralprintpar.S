#!/usr/statlocal/bin/Splus
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
# *  Original programming by Bill Eddy                       *
# ************************************************************/
#
# $Id: spiralprintpar.S,v 1.9 2005/11/01 20:00:44 welling Exp $
#
postscript("test.ps",horizontal=T);
shift <- c("Spikes","Mean Adjustment","Outliers");
par(mfrow=c(inc,3));

if (exists("file.access")) {
  # This is R, not S
  access <- file.access;
}

if (!access("meancprm")) {
  mc <- matrix(scan("meancprm"));
} else {
  mc <- matrix(0,nrow=nimage*nslice,ncol=1);
}

if (!access("spikeprm")) {
  spk <- matrix(scan("spikeprm"),ncol=3,byrow=T)[,3];
} else {
  spk <- matrix(0,nrow=nimage*nslice,ncol=1);
}

if (!access("outlierprm")) {
  ol <- matrix(scan("outlierprm"));
} else {
  ol <- matrix(0,nrow=nimage*nslice,ncol=1);
}

lclimgcnt <- length(mc)/nslice
for (pagenum in 0:(((nslice-1)/inc))){
  for (rownum in 0:(min(c(inc-1,nslice-pagenum*inc-1)))){
      h <- rownum+pagenum*inc
      sub <- (0:(lclimgcnt-1))*nslice + h+1
      rowsub <- (0:(lclimgcnt-1)) + (h*nimage) + 1
      plot(1:lclimgcnt,spk[rowsub],xlab="Image number",ylab=shift[1],main=paste("name Slice",h),type='h')
      plot(1:lclimgcnt,mc[sub],xlab="Image number",ylab=shift[2],main=paste("name Slice",h),pch='.')
#      lines(supsmu((1:lclimgcnt),mc[sub]))
      plot(1:lclimgcnt,ol[sub],xlab="Image number",ylab=shift[3],main=paste("name Slice",h),type='h')
  }
}
