#
#	Makefile for morphology
#
#	Copyright (c) 1996 Department of Statistics, Carnegie Mellon University
#
#	HISTORY
#		1/01	Written by Joel S. Welling (CMU)
#

PKG          = morphology
PKG_MAKEBINS = $(CB)/affinemap $(CB)/follow $(CB)/pullback $(CB)/morph \
	$(CB)/watershed

PKG_LIBS     = -lfmri -lmri -lpar -lbio -lacct -lmisc -lcrg $(LAPACK_LIBS) -lm

ALL_MAKEFILES= Makefile
CSOURCE= affinemap.c follow.c pullback.c morph.c watershed.c
DOCFILES= affinemap_help.help follow_help.help pullback_help.help \
	  morph_help.help watershed_help.help

include ../Makefile_pkg

PKG_OBJS = $O/affinemap.o $O/follow.o $O/pullback.o $O/morph.o $O/watershed.o

LIBFILES= $L/libmri.a $L/libpar.a $L/libbio.a $L/libarray.a $L/libmisc.a \
	$L/libacct.a
HDRS= fmri.h mri.h par.h bio.h misc.h acct.h stdcrg.h

$O/affinemap.o: affinemap.c
	$(CC_RULE)

$O/affinemap_help.o: affinemap_help.help
	$(HELP_RULE)

$(CB)/affinemap: $O/affinemap.o $O/affinemap_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/follow.o: follow.c
	$(CC_RULE)

$O/follow_help.o: follow_help.help
	$(HELP_RULE)

$(CB)/follow: $O/follow.o $O/follow_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/watershed.o: watershed.c
	$(CC_RULE)

$O/watershed_help.o: watershed_help.help
	$(HELP_RULE)

$(CB)/watershed: $O/watershed.o $O/watershed_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/pullback.o: pullback.c
	$(CC_RULE)

$O/pullback_help.o: pullback_help.help
	$(HELP_RULE)

$(CB)/pullback: $O/pullback.o $O/pullback_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/morph.o: morph.c
	$(CC_RULE)

$O/morph_help.o: morph_help.help
	$(HELP_RULE)

$(CB)/morph: $O/morph.o $O/morph_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/flow.o: flow.c
	$(CC_RULE)

$O/flow_help.o: flow_help.help
	$(HELP_RULE)

$(CB)/flow: $O/flow.o $O/flow_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

releaseprep:
	echo "no release prep from " `pwd`

