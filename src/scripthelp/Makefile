#
#	Makefile for scripthelp
#
#	Copyright (c) 2001 Department of Statistics, Carnegie Mellon University
#
#	HISTORY
#		6/01  Modified from existing code by Joel Welling
#

PKG          = scripthelp
PKG_MAKEBINS = $(CB)/scripthelp

PKG_LIBS     = -lfmri -lmri -lpar -lbio -lacct -lcrg -lmisc $(LAPACK_LIBS) -lm

MAKEFILES= Makefile
CSOURCE= scripthelp.c
DOCFILES= scripthelp_help.help

include ../Makefile.pkg

LIBFILES= $L/libmri.a $L/libpar.a $L/libbio.a $L/libarray.a $L/libmisc.a \
	$L/libacct.a
HDRS= fmri.h mri.h par.h bio.h misc.h acct.h stdcrg.h

$O/scripthelp.o: scripthelp.c
	$(CC_RULE)

$O/scripthelp_help.o: scripthelp_help.help
	$(HELP_RULE)

$(CB)/scripthelp: $O/scripthelp.o $O/scripthelp_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

releaseprep:
	echo "no release prep from " `pwd`

