#
#	Makefile for partialk
#
#	Copyright (c) 1997  Pittsburgh Supercomputing Center
#
#	HISTORY
#		12/97	Adapted from baseline Makefile by Joel Welling, CMU
#

PKG          = partialk
PKG_MAKEBINS = $(CB)/partialk

PKG_LIBS     = -lfmri -lmri -lpar -lbio -lacct -lmisc -lcrg $(LAPACK_LIBS) -lm

ALL_MAKEFILES= Makefile
CSOURCE= partialk.c
DOCFILES= partialk_help.help

include ../Makefile_pkg

LIBFILES= $L/libmri.a $L/libpar.a $L/libbio.a $L/libarray.a $L/libmisc.a \
	$L/libacct.a
HDRS= fmri.h mri.h par.h bio.h misc.h acct.h stdcrg.h

$O/partialk.o: partialk.c
	$(CC_RULE)

$O/partialk_help.o: partialk_help.help
	$(HELP_RULE)

$(CB)/partialk: $O/partialk.o $O/partialk_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

releaseprep:
	echo "no release prep from " `pwd`

