#
#	Makefile for the C version of cdflib
#
#	Copyright (c) 1999  Pittsburgh Supercomputing Center
#
#	HISTORY
#		4/99	Written by Joel Welling
#

PKG          = libcdf
PKG_EXPORTS  = dcdflib.h
PKG_MAKELIBS = $L/libdcdf.a
PKG_MAKEBINS = cdftester

PKG_LIBS     = $(LAPACK_LIBS)

MAKEFILES= Makefile
CSOURCE = ipmpar.c dcdflib.c cdftester.c
HFILES= dcdflib.h dcdflib_private.h
DOCFILES= README libdcdf_help.help

include ../Makefile.pkg

LIB_OBJ = $O/dcdflib.o $O/ipmpar.o

$L/libdcdf.a: $(LIB_OBJ)
	@echo "%%%% Building libdcdf.a %%%%"
	@$(AR) $(ARFLAGS) $L/libdcdf.a $(LIB_OBJ)
	@$(RANLIB) $L/libdcdf.a

$O/dcdflib.o: dcdflib.c
	$(CC_RULE)

$O/ipmpar.o: ipmpar.c
	$(CC_RULE)

cdftester: $O/cdftester.o $(LIB_OBJ)
	@echo "%%%% Linking $(@F) %%%%"
	@$(LD) -o cdftester $(LFLAGS) $O/cdftester.o $(LIB_OBJ) $(PKG_LIBS) -lm

$O/cdftester.o: cdftester.c
	$(CC_RULE)

releaseprep:
	echo "no release prep from " `pwd`

