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
# $Id: epiprintpar.S,v 1.12 2005/06/10 20:15:00 welling Exp $
#
postscript("test.ps",horizontal=T)
shift <- c("Baseline.1","Baseline.2","Mean Adjustment","Deghost")

baseraw <- matrix(scan("baseadjrw"),ncol=2,byrow=T)
basesmt <- matrix(scan("baseadjsm"),ncol=2,byrow=T)
nvars <- 2
if (domeanc) {
  means <- matrix(scan("meancprm"))
  nvars <- nvars + 1
} 
if (dodeghost) {
  degraw <- matrix(scan("deghostrw"))
  degsmt <- matrix(scan("deghostsm"))
  nvars <- nvars + 1
} 
par(mfrow=c(3,nvars))
images <- length(basesmt[,1])/nslice

for (pagenum in 0:(((nslice-1)/inc))){
  for (rownum in 0:(min(c(inc-1,nslice-pagenum*inc-1)))){
    h <- rownum+pagenum*inc
    sub <- (0:(images-1))*nslice+h+1
    subt <- (1:images) + h*images

    plot(1:images,baseraw[sub,1],xlab="Image number",ylab=shift[1],
	 main=paste("name Slice",h),pch='.')
    lines(1:images,basesmt[sub,1])

    plot(1:images,baseraw[sub,2],xlab="Image number",ylab=shift[2],
	 main=paste("name Slice",h),pch='.')
    lines(1:images,basesmt[sub,2])


    if (domeanc) {
	 plot(1:images,means[sub,1],xlab="Image number",ylab=shift[3],
	      main=paste("name Slice",h),pch='.')
#        lines(supsmu((1:images),means[sub,1]))
    }

    if (dodeghost) {
      plot(1:images,degraw[subt,1],xlab="Image number",ylab=shift[4],
	   main=paste("name Slice",h),pch='.')
      lines(1:images,degsmt[subt,1])
    }
  }
}
