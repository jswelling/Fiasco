#
#	Makefile for intsplit
#
#	Copyright (c) 1996 Department of Statistics, Carnegie Mellon University
#
#	HISTORY
#		12/96	Written by William F. Eddy (CMU)
#

PKG          = intsplit
PKG_MAKEBINS = $(CB)/intsplit

PKG_LIBS     = -lfmri -lmri -lpar -lbio -lacct -lmisc -lcrg $(LAPACK_LIBS) -lm

MAKEFILES= Makefile
CSOURCE= intsplit.c
DOCFILES= intsplit_help.help

include ../Makefile.pkg

LIBFILES= $L/libmri.a $L/libpar.a $L/libbio.a $L/libarray.a $L/libmisc.a \
	$L/libacct.a
HDRS= fmri.h mri.h par.h bio.h misc.h acct.h stdcrg.h

$O/intsplit.o: intsplit.c
	$(CC_RULE)

$O/intsplit_help.o: intsplit_help.help
	$(HELP_RULE)

$(CB)/intsplit: $O/intsplit.o $O/intsplit_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

releaseprep:
	echo "no release prep from " `pwd`

