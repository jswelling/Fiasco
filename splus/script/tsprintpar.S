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
# $Id: tsprintpar.S,v 1.7 2004/02/24 23:29:10 welling Exp $
#
postscript("test.ps",horizontal=T)
shift <- c("Shot 1 Base Real","Shot 1 Base Imag","Shot 2 Base Real",
	"Shot 2 Base Imag","Mean Adjustment")
par(mfrow=c(inc,5))

baseraw <- matrix(scan("par/baseadjrw"),ncol=2,byrow=T)
basesmt <- matrix(scan("par/baseadjsm"),ncol=2,byrow=T)
means <- matrix(scan("par/meancprm"))

images <- length(basesmt[,1])/(2*nslice)
for (pagenum in 0:(((nslice-1)/inc))){
  for (rownum in 0:(min(c(inc-1,nslice-pagenum*inc-1)))){
    h <- rownum+pagenum*inc

    for (shot in 1:2) {
      for (plotnum in 1:2){
        sub <- (0:(images-1))*(2*nslice) + (shot-1)*nslice + h + 1
	plot(1:images,baseraw[sub,plotnum],xlab="Image number",
	     ylab=shift[2*(shot-1)+plotnum],main=paste("name Slice",h),
	     pch='.')
        lines(1:images,basesmt[sub,plotnum])
      }
    }

    sub <- (0:(images-1))*nslice + h + 1
    plot(1:images,means[sub,1],xlab="Image number",ylab=shift[5],
	 main=paste("name Slice",h),pch='.')
   # lines(supsmu((1:images),means[sub,1]))

  }
}
