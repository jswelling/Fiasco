#
#	Makefile for mri utilities
#
#	Copyright (c) 1997 Department of Statistics, Carnegie Mellon University
#
#	HISTORY
#		9/97	Written by Joel S. Welling (CMU)
#

PKG          = physio
PKG_MAKEBINS = $(CB)/cardiopeaks $(CB)/resppeaks $(CB)/syncpeaks $(CB)/eventpeaks

PKG_LIBS     = -lfmri -lmri -lpar -lbio -lacct -lcrg -lmisc $(LAPACK_LIBS) -lm

ALL_MAKEFILES= Makefile
CSOURCE= cardiopeaks.c resppeaks.c polyfit.c syncpeaks.c eventpeaks.c
HFILES= polyfit.h
DOCFILES= cardiopeaks_help.help resppeaks_help.help syncpeaks_help.help eventpeaks_help.help

include ../Makefile_pkg

LIBFILES= $L/libmri.a $L/libpar.a $L/libbio.a $L/libarray.a $L/libmisc.a \
	$L/libacct.a
HDRS= fmri.h mri.h par.h bio.h misc.h acct.h stdcrg.h

$O/cardiopeaks.o: cardiopeaks.c
	$(CC_RULE)

$O/resppeaks.o: resppeaks.c
	$(CC_RULE)

$O/polyfit.o: polyfit.c
	$(CC_RULE)

$O/syncpeaks.o: syncpeaks.c
	$(CC_RULE)

$O/eventpeaks.o: eventpeaks.c
	$(CC_RULE)

$O/cardiopeaks_help.o: cardiopeaks_help.help
	$(HELP_RULE)

$O/resppeaks_help.o: resppeaks_help.help
	$(HELP_RULE)

$O/syncpeaks_help.o: syncpeaks_help.help
	$(HELP_RULE)

$O/eventpeaks_help.o: eventpeaks_help.help
	$(HELP_RULE)

$(CB)/cardiopeaks: $O/cardiopeaks.o $O/polyfit.o $O/cardiopeaks_help.o
	@echo "%%%% Linking $(@F) %%%%" 
	@$(LD) $(LFLAGS) -o $B/cardiopeaks $O/cardiopeaks.o $O/polyfit.o \
		$O/cardiopeaks_help.o $(LIBS)

$(CB)/resppeaks: $O/resppeaks.o $O/polyfit.o $O/resppeaks_help.o
	@echo "%%%% Linking $(@F) %%%%" 
	@$(LD) $(LFLAGS) -o $B/resppeaks $O/resppeaks.o $O/polyfit.o \
		$O/resppeaks_help.o $(LIBS)

$(CB)/syncpeaks: $O/syncpeaks.o $O/polyfit.o $O/syncpeaks_help.o
	@echo "%%%% Linking $(@F) %%%%" 
	@$(LD) $(LFLAGS) -o $B/syncpeaks $O/syncpeaks.o $O/polyfit.o \
		$O/syncpeaks_help.o $(LIBS)

$(CB)/eventpeaks: $O/eventpeaks.o $O/polyfit.o $O/eventpeaks_help.o
	@echo "%%%% Linking $(@F) %%%%" 
	@$(LD) $(LFLAGS) -o $B/eventpeaks $O/eventpeaks.o $O/polyfit.o \
		$O/eventpeaks_help.o $(LIBS)

releaseprep:
	echo "no release prep from " `pwd`

