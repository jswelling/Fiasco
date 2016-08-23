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
# $Id: printregist3d.S,v 1.7 2006/10/23 21:04:59 welling Exp $
#
#
postscript("test.ps",horizontal=T)

if (exists("file.access")) {
  # This is R, not S
  access <- file.access;
}
	
shift<-c("3D Quat.x", "3D Quat.y", "3D Quat.z", "3D Quat.w",
"3D X shift (voxels)", "3D Y shift (voxels)", "3D Z shift (voxels)", "3D MSE",
"3D Net Rotation", "3D Net Translation (mm)", "3D Displacement (mm)",
"3D Meanc")

par(mfrow=c(3,4))
x <- matrix(scan("rreg"),ncol=9,byrow=T)[,2:9]

if (dosmreg3d) {
  y <- matrix(scan("sxreg"),ncol=9,byrow=T)[,2:9]
} else {
  y <- matrix(0,nrow=nimage,ncol=8)
}

if (dordisp3d) {
  x <- cbind(x,matrix(scan("dixplace"),ncol=4,byrow=T)[,2:4])
} else {
  x <- cbind(x,matrix(0,ncol=3,nrow=nimage))
}

if (dosdisp3d) {
  y <- cbind(y,matrix(scan("smdisplace"),ncol=4,byrow=T)[,2:4])
} else {
  y <- cbind(y,matrix(0,ncol=3,nrow=nimage))
}

if (domnc3d) {
  mc <- matrix(scan("meanc3d"));
} else {
  mc <- matrix(0,nrow=nimage,ncol=1);
}
x <- cbind(x,mc)
y <- cbind(y,mc)

images <- length(x[,1])
for (plotnum in 12:12) {
  plot(1:images,x[,plotnum],xlab="Image number",ylab=shift[plotnum],
        main="name",pch='.')
}
for (plotnum in 9:11){
  plot(1:images,x[,plotnum],xlab="Image number",ylab=shift[plotnum],
	main="name",pch='.')
  lines(1:images,y[,plotnum])
}
for (plotnum in 1:7){
  plot(1:images,x[,plotnum],xlab="Image number",ylab=shift[plotnum],
	main="name",pch='.')
  lines(1:images,y[,plotnum])
}
for (plotnum in 8:8){
  plot(1:images,x[,plotnum],xlab="Image number",ylab=shift[plotnum],
	main="name",pch='.')
  lines(1:images,y[,plotnum])
}

