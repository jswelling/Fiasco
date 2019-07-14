#
#	Makefile for outlier
#
#	Copyright (c) 1996 Department of Statistics, Carnegie Mellon University
#
#	HISTORY
#		12/96	Written by William F. Eddy (CMU)
#

PKG          = outlier
PKG_MAKEBINS = $(CB)/outlier

PKG_LIBS     = -lfmri -lmri -lpar -lbio -lacct -lmisc -lcrg $(LAPACK_LIBS) -lm

ALL_MAKEFILES= Makefile
CSOURCE= outlier.c
DOCFILES= outlier_help.help

include ../Makefile_pkg

LIBFILES= $L/libmri.a $L/libpar.a $L/libbio.a $L/libarray.a $L/libmisc.a \
	$L/libacct.a
HDRS= fmri.h mri.h par.h bio.h misc.h acct.h stdcrg.h

$O/outlier.o: outlier.c
	$(CC_RULE)

$O/outlier_help.o: outlier_help.help
	$(HELP_RULE)

$(CB)/outlier: $O/outlier.o $O/outlier_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

releaseprep:
	echo "no release prep from " `pwd`

