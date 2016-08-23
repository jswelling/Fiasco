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
# $Id: anova.S,v 1.3 2004/02/24 23:29:10 welling Exp $
#

#interpretsplit <- function(spltfile)
#{
#  splt <- read.table(spltfile)
#  dm <- dim(splt)
#  nfac <- dm[2] - 1
#  facs <- as.matrix(splt[,1:nfac])
#  dimnames(facs) <- list(NULL,NULL)
#  noim <- as.numeric(splt[,dm[2]])
#  names(noim) <- NULL
#  nim <- sum(noim)
#  output <- list(grps=matrix(nrow=dm[1],ncol=2),
#		 desgn=data.frame())
#  output$grps[,2] <- noim
#  glst <- facs
#  ngr <- 0
#  for(i in 1:dm[1])
#    {
#		idd <- F
#		if(prod(facs[i,]))
#		{
#			if(ngr>0)
#			{
#				for(j in 1:ngr)
#				{
#					if(all(facs[i,]==glst[j,]))
#					{
#						idd <- T
#						output$grps[i,1] <- j
#						break
#					}
#				}
#			}
#		} else {
#			idd <- T
#			output$grps[i,1] <- 0
#		}	
#		if(!idd)
#		{
#			ngr <- ngr+1
#			glst[ngr,] <- facs[i,]
#			output$grps[i,1] <- ngr
#		}
#	}
#	dmat <- matrix(nrow=nim,ncol=nfac)
#	imno <- 0
#	for(i in 1:dm[1])
#	{
#		if(!(prod(facs[i,])))
#		{
#			facs[i,] <- rep(NA,nfac)
#		}
#	}
#	for(i in 1:dm[1])
#	{
#		dmat[(imno+1):(imno+noim[i]),] <-
#		  matrix(rep(facs[i,],rep(noim[i],nfac)),ncol=nfac)
#		imno <- imno + noim[i]
#	}
#	output$desgn <- as.data.frame(dmat)
#	output
#}

interpretsplit <- function(spltfile,nslices)
{
  splt <- read.table(spltfile)
  cvals <- splt[,3]
  cvalsBySlice <- t(matrix(cvals,ncol=nslices))
  cvalsBySlice
}

runAnovas <- function (x, tdim, xdim, ydim, zdim, formula)
{
  x<-x[,1:length(attr(terms(formula), "variables"))-1]
  m <- length(terms(formula)) # or: dim(attr(terms(formula),"factors"))[2]

  FValues <- array(0, dim=c(xdim, ydim, zdim ,m))
  inBuf <- array(0, dim=c(tdim, xdim))
  outBuf <- array(0, dim=c(xdim, m))

  for (sl in 0:(zdim-1))
   for (y in 0:(ydim-1))
   {
     inBuf <- BioReadFloats(inBuf)
     RES <- t(matrix(inBuf, nrow=xdim))
     x.new <- data.frame(t(x[sl,]),I(RES))
     names(x.new) <- c(names(x), "RES")

     wx <- aov(formula, x.new)
     swx <- summary(wx)
     for (i in 1:xdim)
       FValues [i, y+1, sl+1, ] <- as.matrix(swx[[i]])[1:m,4]
   }

  return (FValues)
}

x <- interpretsplit(getenv("F_SPLIT_FTZ"),6)
formula <- eval(parse(getenv("F_FORMULA_FTZ")))
ans <- runAnovas(x, 6, 64, 64, 4, formula)
ans[is.na(ans)]<- -1   # This finds all invalid results and substitutes -1
dims <- dim (ans)
dims[4] # emits number of following datasets
for (i in 1:dims[4]){
  paste(getenv("F_HEADER"),".F.",
	       labels(attr(terms(formula),"factors"))[[2]][i], sep="")
  BioWriteFloats(ans[,,,i])
} 
