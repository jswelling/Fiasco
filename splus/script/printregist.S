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
# $Id: printregist.S,v 1.8 2005/07/30 02:09:28 welling Exp $
#
#
postscript("test.ps",horizontal=T)
shift<-c("X Shift in Pixels","Y Shift in Pixels",
	"Rotation in Degrees","Average Displacement in Pixels")
p <- 4
if ( nimage > 400 ) p <- 2
par(mfrow=c(inc,p))
x <- matrix(scan("rreg"),ncol=6,byrow=T)[,3:5]
x<-cbind(x,scan("dixplace"))
y<-matrix(scan("sxreg"),ncol=6,byrow=T)[,3:5]
y<-cbind(y,scan("smdisplace"))
images<-length(x[,1])/nslice
for (pagenum in 0:(((nslice-1)/inc))){
  for (rownum in 0:(min(c(inc-1,nslice-pagenum*inc-1)))){
    for (plotnum in 1:4){
      h<-rownum+pagenum*inc
      sub<-(0:(images-1))*nslice+h+1
#      plot(1:images,x[h*images+(1:images),plotnum],xlab="Image number",ylab=shift[plotnum],main=paste("name Slice",h))
      plot(1:images,x[sub,plotnum],xlab="Image number",ylab=shift[plotnum],main=paste("name Slice",h),pch='.')
#    lines(1:images,y[h*images+(1:images),plotnum])
    lines(1:images,y[sub,plotnum])
    }
  }
}
