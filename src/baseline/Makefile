#
#	Makefile for baseline
#
#	Copyright (c) 1996  Pittsburgh Supercomputing Center
#
#	HISTORY
#		9/96	Written by Greg Hood (PSC)
#

PKG          = baseline
PKG_MAKEBINS = $(CB)/baseline

PKG_LIBS     = -lfmri -lmri -lpar -lbio -lacct -lmisc -lcrg $(LAPACK_LIBS) -lm

LIBFILES= $L/libmri.a $L/libpar.a $L/libbio.a $L/libarray.a $L/libmisc.a \
	$L/libacct.a
HDRS= fmri.h mri.h par.h bio.h misc.h acct.h stdcrg.h

ALL_MAKEFILES= Makefile
CSOURCE= baseline.c
DOCFILES= baseline_help.help

include ../Makefile_pkg

$O/baseline.o: baseline.c
	$(CC_RULE)

$O/baseline_help.o: baseline_help.help
	$(HELP_RULE)

$(CB)/baseline: $O/baseline.o $O/baseline_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

releaseprep:
	echo "no release prep from " `pwd`

