#
#	Makefile for detrend
#
#	Copyright (c) 1996 Department of Statistics, Carnegie Mellon University
#
#	HISTORY
#		12/96	Written by William F. Eddy (CMU)
#

PKG          = detrend
PKG_MAKEBINS = $(CB)/detrend

PKG_LIBS     = -lfmri -lmri -lpar -lbio -lacct -lmisc -lcrg $(LAPACK_LIBS) -lm

LIBFILES= $L/libmri.a $L/libpar.a $L/libbio.a $L/libarray.a $L/libmisc.a \
	$L/libacct.a
HDRS= fmri.h mri.h par.h bio.h misc.h acct.h stdcrg.h

ALL_MAKEFILES= Makefile
CSOURCE= detrend.c
DOCFILES= detrend_help.help

include ../Makefile_pkg

$O/detrend.o: detrend.c
	$(CC_RULE)

$O/detrend_help.o: detrend_help.help
	$(HELP_RULE)

$(CB)/detrend: $O/detrend.o $O/detrend_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

releaseprep:
	echo "no release prep from " `pwd`

