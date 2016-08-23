#
#	Makefile for baseline2
#
#	Copyright (c) 1996  Pittsburgh Supercomputing Center
#
#	HISTORY
#		9/96	Written by Greg Hood (PSC)
#

PKG          = baseline2
PKG_MAKEBINS = $(CB)/baseline2

PKG_LIBS     = -lfmri -lmri -lpar -lbio -lacct -lmisc -lcrg $(LAPACK_LIBS) -lm

LIBFILES= $L/libmri.a $L/libpar.a $L/libbio.a $L/libarray.a $L/libmisc.a \
	$L/libacct.a
HDRS= fmri.h mri.h par.h bio.h misc.h acct.h stdcrg.h

MAKEFILES= Makefile
CSOURCE= baseline2.c
DOCFILES= baseline2_help.help

include ../Makefile.pkg

$O/baseline2.o: baseline2.c
	$(CC_RULE)

$O/baseline2_help.o: baseline2_help.help
	$(HELP_RULE)

$(CB)/baseline2: $O/baseline2.o $O/baseline2_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

releaseprep:
	echo "no release prep from " `pwd`

