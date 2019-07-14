#
#	Makefile for summary
#
#	Copyright (c) 1996 Department of Statistics, Carnegie Mellon University
#
#	HISTORY
#		12/96	Written by William F. Eddy (CMU)
#

PKG          = summary
PKG_MAKEBINS = $(CB)/summary

PKG_LIBS     = -lfmri -lmri -lpar -lbio -lacct -lmisc -lcrg $(LAPACK_LIBS) -lm

ALL_MAKEFILES= Makefile
CSOURCE= summary.c
DOCFILES= summary_help.help

include ../Makefile_pkg

LIBFILES= $L/libmri.a $L/libpar.a $L/libbio.a $L/libarray.a $L/libmisc.a \
	$L/libacct.a
HDRS= fmri.h mri.h par.h bio.h misc.h acct.h stdcrg.h

$O/summary.o: summary.c
	$(CC_RULE)

$O/summary_help.o: summary_help.help
	$(HELP_RULE)

$(CB)/summary: $O/summary.o $O/summary_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

releaseprep:
	echo "no release prep from " `pwd`

