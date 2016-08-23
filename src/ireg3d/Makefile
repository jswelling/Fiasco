#
#	Makefile for ireg3d
#
#	Copyright (c) 1996 Department of Statistics, Carnegie Mellon University
#
#	HISTORY
#		12/96	Written by William F. Eddy (CMU)
#

PKG          = ireg3d
PKG_MAKEBINS = $(CB)/ireg3d $(CB)/iwarp

PKG_LIBS     = -lfmri -lmri -lpar -lbio -lacct -lmisc -lcrg $(LAPACK_LIBS) -lm

MAKEFILES= Makefile
CSOURCE= ireg3d.c iwarp.c
DOCFILES= ireg3d_help.help iwarp_help.help

include ../Makefile.pkg

LIBFILES= $L/libmri.a $L/libpar.a $L/libbio.a $L/libarray.a $L/libmisc.a \
	$L/libacct.a $L/libfmri.a
HDRS= fmri.h mri.h par.h bio.h misc.h acct.h stdcrg.h

$O/ireg3d.o: ireg3d.c
	$(CC_RULE)

$O/ireg3d_help.o: ireg3d_help.help
	$(HELP_RULE)

$(CB)/ireg3d: $O/ireg3d.o $O/ireg3d_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/iwarp.o: iwarp.c
	$(CC_RULE)

$O/iwarp_help.o: iwarp_help.help
	$(HELP_RULE)

$(CB)/iwarp: $O/iwarp.o $O/iwarp_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

releaseprep:
	echo "no release prep from " `pwd`

