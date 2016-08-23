#
#	Makefile for stats
#
#	Copyright (c) 1996 Department of Statistics, Carnegie Mellon University
#
#	HISTORY
#		12/96	Written by William F. Eddy (CMU)
#

PKG          = stats
PKG_MAKEBINS = $(CB)/stats

PKG_LIBS     = -lfmri -lmri -lpar -lbio -lacct -lmisc -lcrg $(LAPACK_LIBS) -lm

MAKEFILES= Makefile
CSOURCE= stats.c
DOCFILES= stats_help.help

include ../Makefile.pkg

LIBFILES= $L/libmri.a $L/libpar.a $L/libbio.a $L/libarray.a $L/libmisc.a \
	$L/libacct.a
HDRS= fmri.h mri.h par.h bio.h misc.h acct.h stdcrg.h

$O/stats.o: stats.c
	$(CC_RULE)

$O/stats_help.o: stats_help.help
	$(HELP_RULE)

$(CB)/stats: $O/stats.o $O/stats_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

releaseprep:
	echo "no release prep from " `pwd`

