#
#	Makefile for splus scripts
#
#	Copyright (c) 1997 Department of Statistics, Carnegie Mellon University
#
#	HISTORY
#		12/96	Written by William F. Eddy (CMU)
#

PKG          = script
PKG_MAKETARGETS = S_files

SSOURCE= anova.S epiprintpar.S printregist.S printregist3d.S \
	spiralprintpar.S tsprintpar.S read.S
MISCFILES= 
MAKEFILES= Makefile

include ../Makefile.splus

S_files: read.S
	@echo '%%%% read.S is no longer compiled by default %%%%'
#	cat read.S | sed "s/PVM_ARCH/$(ARCH)/g" | $(SPLUSEXE)





