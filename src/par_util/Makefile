#
#	Makefile for mri utilities
#
#	Copyright (c) 2000 Department of Statistics, Carnegie Mellon University
#
#	HISTORY
#		6/00	Written by Joel S. Welling (CMU)
#

PKG          = par_util
PKG_MAKEBINS = $(CB)/par_translate

PKG_LIBS     = -ldcdf -lfmri -lmri -lpar -lbio -lacct -lmisc \
	-lcrg $(LAPACK_LIBS) -lm

MAKEFILES= Makefile
CSOURCE= par_translate.c
HFILES= parameters.h
DOCFILES= par_translate_help.help

include ../Makefile.pkg

LIBFILES= $L/libmri.a $L/libpar.a $L/libbio.a $L/libarray.a $L/libmisc.a \
	$L/libacct.a

$O/par_translate.o: par_translate.c
	$(CC_RULE)

$O/par_translate_help.o: par_translate_help.help
	$(HELP_RULE)

$(CB)/par_translate: $O/par_translate.o $O/par_translate_help.o 
	@echo "%%%% Linking $(@F) %%%%" 
	@$(LD) $(LFLAGS) -o $B/$(@F) $O/par_translate.o \
		$O/par_translate_help.o $(LIBS)

releaseprep:
	echo "no release prep from " `pwd`

