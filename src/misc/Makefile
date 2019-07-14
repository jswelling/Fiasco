#
#	Makefile for tools that don't comfortably fit elsewhere.
#
#	Copyright (c) 1996  Pittsburgh Supercomputing Center
#
#	HISTORY
#		1/96	Written by Greg Hood (PSC)
#		8/96	Made architecture-independent (Greg Hood, PSC)
#

PKG          = misc
PKG_EXPORTS  = 
PKG_MAKELIBS = 
PKG_MAKEBINS = $(CB)/prefix $(CB)/pulsecalc $(CB)/blind_minimize \
	     $(CB)/parallel_mode $(CB)/image_entropy \
	     $(CB)/mutual_information $(CB)/fiat_version
PKG_LIBS     = -ldcdf -lfmri -lmri -lpar -lbio -lacct \
	     -lcrg -lmisc -lpulse $(LAPACK_LIBS) -lm

ALL_MAKEFILES= Makefile
CSOURCE= pulsecalc.c blind_minimize.c parallel_mode.c image_entropy.c \
	mutual_information.c fiat_version.c
HFILES= 
DOCFILES= blind_minimize_help.help image_entropy_help.help \
	mutual_information_help.help
SCRIPTFILES= prefix.csh
MISCFILES= sed.script

LIBFILES= $L/libfmri.a $L/libmri.a $L/libpar.a $L/libbio.a \
	  $L/libarray.a $L/libmisc.a $L/libacct.a $L/libcrg.a

include ../Makefile_pkg

$O/fiat_version.o: fiat_version.c
	$(CC_RULE)

$(CB)/fiat_version: $O/fiat_version.o
	$(SINGLE_LD)

$O/parallel_mode.o: parallel_mode.c 
	$(CC_RULE)

$(CB)/parallel_mode: $O/parallel_mode.o
	$(SINGLE_LD)

$O/blind_minimize.o: blind_minimize.c 
	$(CC_RULE)

$O/blind_minimize_help.o: blind_minimize_help.help \
			  ../fmri/praxis_help.help \
			  ../fmri/nelmin_help.help
	$(HELP_RULE)

$(CB)/blind_minimize: $O/blind_minimize.o $O/blind_minimize_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/image_entropy.o: image_entropy.c 
	$(CC_RULE)

$O/image_entropy_help.o: image_entropy_help.help
	$(HELP_RULE)

$(CB)/image_entropy: $O/image_entropy.o $O/image_entropy_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mutual_information.o: mutual_information.c 
	$(CC_RULE)

$O/mutual_information_help.o: mutual_information_help.help
	$(HELP_RULE)

$(CB)/mutual_information: $O/mutual_information.o $O/mutual_information_help.o \
		$(LIBFILES)
	$(SINGLE_HELP_LD)

$O/pulsecalc.o: pulsecalc.c
	$(CC_RULE)

$(CB)/pulsecalc: $O/pulsecalc.o $(LIBFILES)
	$(SINGLE_LD)

$(CB)/prefix: prefix.csh
	@echo "%%%% Copying in $(@F) %%%%"
	@cp -f prefix.csh $(CB)/prefix
	@chmod a+x $B/prefix

releaseprep:
	echo "no release prep from " `pwd`

