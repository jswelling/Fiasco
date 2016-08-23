#
#	Makefile for phadj
#
#	Copyright (c) 1996 Department of Statistics, Carnegie Mellon University
#
#	HISTORY
#		12/96	Written by William F. Eddy (CMU)
#

PKG          = phadj
PKG_MAKEBINS = $(CB)/phadj

PKG_LIBS     = -lfmri -lmri -lpar -lbio -lacct -lmisc -lcrg $(LAPACK_LIBS) -lm

MAKEFILES= Makefile
CSOURCE= phadj.c
DOCFILES= phadj_help.help

include ../Makefile.pkg

LIBFILES= $L/libmri.a $L/libpar.a $L/libbio.a $L/libarray.a $L/libmisc.a \
	$L/libacct.a
HDRS= fmri.h mri.h par.h bio.h misc.h acct.h stdcrg.h

$O/phadj.o: phadj.c
	$(CC_RULE)

$O/phadj_help.o: phadj_help.help
	$(HELP_RULE)

$(CB)/phadj: $O/phadj.o $O/phadj_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

releaseprep:
	echo "no release prep from " `pwd`

