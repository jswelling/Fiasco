#
#	Makefile for deghost
#
#	Copyright (c) 1996 Department of Statistics, Carnegie Mellon University
#
#	HISTORY
#		12/96	Written by William F. Eddy (CMU)
#

PKG          = deghost
PKG_MAKEBINS = $(CB)/deghost $(CB)/nav_deghost

PKG_LIBS     = -lfmri -lmri -lpar -lbio -lacct -lmisc -lcrg \
	     $(LAPACK_LIBS) -lm

LIBFILES= $L/libmri.a $L/libpar.a $L/libbio.a $L/libarray.a $L/libmisc.a \
	$L/libacct.a $L/libfmri.a
HDRS= fmri.h mri.h par.h bio.h misc.h acct.h stdcrg.h

ALL_MAKEFILES= Makefile
CSOURCE= deghost.c nav_deghost.c
DOCFILES= deghost_help.help nav_deghost_help.help

include ../Makefile_pkg

$O/deghost.o: deghost.c
	$(CC_RULE)

$O/deghost_help.o: deghost_help.help
	$(HELP_RULE)

$(CB)/deghost: $O/deghost.o $O/deghost_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/nav_deghost.o: nav_deghost.c
	$(CC_RULE)

$O/nav_deghost_help.o: nav_deghost_help.help
	$(HELP_RULE)

$(CB)/nav_deghost: $O/nav_deghost.o $O/nav_deghost_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

releaseprep:
	echo "no release prep from " `pwd`

