#
#	Makefile for meanc
#
#	Copyright (c) 1996 Department of Statistics, Carnegie Mellon University
#
#	HISTORY
#		12/96	Written by William F. Eddy (CMU)
#

PKG          = meanc
PKG_MAKEBINS = $(CB)/meanc $(CB)/meanc3d

PKG_LIBS     = -lfmri -lmri -lpar -lbio -lacct -lmisc -lcrg $(LAPACK_LIBS) -lm

ALL_MAKEFILES= Makefile
CSOURCE= meanc.c meanc3d.c
DOCFILES= meanc_help.help meanc3d_help.help

include ../Makefile_pkg

LIBFILES= $L/libmri.a $L/libpar.a $L/libbio.a $L/libarray.a $L/libmisc.a \
	$L/libacct.a
HDRS= fmri.h mri.h par.h bio.h misc.h acct.h stdcrg.h

$O/meanc.o: meanc.c
	$(CC_RULE)

$O/meanc_help.o: meanc_help.help
	$(HELP_RULE)

$(CB)/meanc: $O/meanc.o $O/meanc_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/meanc3d.o: meanc3d.c
	$(CC_RULE)

$O/meanc3d_help.o: meanc3d_help.help
	$(HELP_RULE)

$(CB)/meanc3d: $O/meanc3d.o $O/meanc3d_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

releaseprep:
	echo "no release prep from " `pwd`

