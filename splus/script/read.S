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
# $Id: read.S,v 1.2 2004/02/24 23:29:10 welling Exp $
#
attach(paste(getenv("FIASCO"),"/../../splus/S/PVM_ARCH",sep=""),pos=1)
readHead <- function (name)
{
	 type <- ndim <- vec <-  0
	 dims <- rep(0, 100)
	 storage.mode (type) <- "integer"
	 storage.mode (ndim) <- "integer"
	 storage.mode (vec) <- "integer"
	 storage.mode (dims) <- "integer"
         tmp <- .C("readFieldHead", 
         "name"=name, 
         "type"=type,
	 "ndim"=ndim, 
         "veclen" =vec,
         "dims"=dims)
	tmp
}

writeField <- function (name, field,type,ndim,veclen,dims)
{
	 storage.mode (type) <- "integer"
	 storage.mode (ndim) <- "integer"
	 storage.mode (veclen) <- "integer"
	 storage.mode (dims) <- "integer"
	 storage.mode (field) <- "single"
         .C("writeField", 
         "name"=name, 
         "type"=type,
	 "ndim"=ndim, 
         "veclen" =veclen,
         "dims"=dims,
	 "val"=field)
	return (0)
}

readField <- function (iL=rep(0,4), iU=rep(1,4))
{
	 res <- array (0, prod((iU-iL) + 1)) 
	 storage.mode (iL) <- "integer"
	 storage.mode (iU) <- "integer"
	 storage.mode (res) <- "single"
         tmp <- .C("readField", 
         "iL"=iL,
         "iU"=iU,
         "value"=res)
	tmp
}

runAnovas <- function (name, x, formula)
{
  head <- readHead (name)
  adim <- head$dims [1]
  cdim <- head$dims [4]
  m <- length(terms(formula)) # or: dim(attr(terms(formula),"factors"))[2]

  FValues <- array(0, dim=c(head$dims [1], head$dims [2],head$dims [3],m))
  x<-x[,1:length(attr(terms(formula), "variables"))-1]

  for (sl in 0:(head$dims [3]-1))
   for (y in 0:(head$dims [2]-1))
   {
     RES <- t(matrix(readField ( iL=c(0, y, sl, 0),
               iU=c(adim-1, y, sl, cdim-1))$value, nrow=adim)) 
     x.new <- data.frame(x,I(RES))
     names(x.new) <- c(names(x), "RES")

     wx <- aov(formula, x.new)
     swx <- summary(wx)

     for (i in 1:adim)
       FValues [i, y+1, sl+1, ] <- as.matrix(swx[[i]])[1:m,4]
   }

  return (FValues)
}

interpretsplit <- function(spltfile)
{
	splt <- read.table(spltfile)
	dm <- dim(splt)
	nfac <- dm[2] - 1
	facs <- as.matrix(splt[,1:nfac])
	dimnames(facs) <- list(NULL,NULL)
	noim <- as.numeric(splt[,dm[2]])
	names(noim) <- NULL
	nim <- sum(noim)
	output <- list(grps=matrix(nrow=dm[1],ncol=2),
	  desgn=data.frame())
	output$grps[,2] <- noim
	glst <- facs
	ngr <- 0
	for(i in 1:dm[1])
	{
		idd <- F
		if(prod(facs[i,]))
		{
			if(ngr>0)
			{
				for(j in 1:ngr)
				{
					if(all(facs[i,]==glst[j,]))
					{
						idd <- T
						output$grps[i,1] <- j
						break
					}
				}
			}
		} else {
			idd <- T
			output$grps[i,1] <- 0
		}	
		if(!idd)
		{
			ngr <- ngr+1
			glst[ngr,] <- facs[i,]
			output$grps[i,1] <- ngr
		}
	}
	dmat <- matrix(nrow=nim,ncol=nfac)
	imno <- 0
	for(i in 1:dm[1])
	{
		if(!(prod(facs[i,])))
		{
			facs[i,] <- rep(NA,nfac)
		}
	}
	for(i in 1:dm[1])
	{
		dmat[(imno+1):(imno+noim[i]),] <-
		  matrix(rep(facs[i,],rep(noim[i],nfac)),ncol=nfac)
		imno <- imno + noim[i]
	}
	output$desgn <- as.data.frame(dmat)
	output
}


