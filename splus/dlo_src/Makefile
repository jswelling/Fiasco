#
#	Makefile for splus
#
#	Copyright (c) 1997 Department of Statistics, Carnegie Mellon University
#
#	HISTORY
#		12/96	Written by William F. Eddy (CMU)
#

PKG          = dlo_src

SRC= read_mri.c
CSOURCE= read.c read_mri.c splus_libbio.c splus_libmri.c 
HFILES= mr.h splus_libbio.h splus_libmri.h
ALL_MAKEFILES= Makefile
PKG_MAKETARGETS = splus_dlos
PKG_MAKEOBJS = $(DLO)/read_mri.o $(DLO)/splus_libmri.o $(DLO)/splus_libbio.o

include ../Makefile_splus

splus_dlos: $(PKG_MAKEOBJS)

$(DLO)/read_mri.o: read_mri.c
	@echo "%%%% Compiling $(@F:.o=.c) %%%%" 
	@$(CC) -c $(CFLAGS) -o $@ read_mri.c

$(DLO)/splus_libmri.o: splus_libmri.c
	@echo "%%%% Compiling $(@F:.o=.c) %%%%" 
	@$(CC) -c $(CFLAGS) -o $@ splus_libmri.c

$(DLO)/splus_libbio.o: splus_libbio.c
	@echo "%%%% Compiling $(@F:.o=.c) %%%%" 
	@$(CC) -c $(CFLAGS) -o $@ splus_libbio.c
