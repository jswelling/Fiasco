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
# $Id: printpca.S,v 1.2 2004/12/18 03:21:12 welling Exp $
#
postscript("test.ps",horizontal=T)
shift <- c("Baseline.1","Baseline.2","Mean Adjustment","Deghost")

evals <- matrix(scan("evalsFile"),ncol=2,byrow=T)
Fs <- matrix(scan("fsFile"),ncol=2,byrow=T)
Ps <- matrix(scan("psFile"),ncol=2,byrow=T)
MSEs <- matrix(scan("msesFile"),ncol=2,byrow=T)
evecs <- matrix(scan("evecsFile"),ncol=2,byrow=T)
hdrstring <- "headerString"

nComp <- ncomponents
nTimes <- length(evecs)/(nComp*2)
nPlots <- nComp + 1
rowsPerPage <- 4
colsPerRow <- 4

par(mfrow=c(rowsPerPage, colsPerRow))

plot(1:nComp, evals[1:nComp,2], xlab="Component number", ylab="Magnitude",
     main="Eigenvalues")

title(main=hdrstring,outer=TRUE,line=-1)

plot(1:(nComp+1), MSEs[1:(nComp+1),2], xlab="Component number", ylab="MSE",
     main="Termwise MSEs incl. GrandMean")

plot(1:nComp, Fs[1:nComp,2], xlab="Component number", ylab="F Score",
     main="Termwise F Scores")

plot(1:nComp, Ps[1:nComp,2], xlab="Component number", ylab="log(P)",
     main="Termwise P Scores")

par(mfg=c(2,1))
for ( i in 1:nComp ) {
  plot(1:nTimes, evecs[(((i-1)*nTimes)+1):(i*nTimes),2], 
       type="l", main=paste("Component ",i), xlab="Image", ylab="Val")
}


