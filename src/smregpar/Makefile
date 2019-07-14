#
#	Makefile for smregpar
#
#	Copyright (c) 1996 Department of Statistics, Carnegie Mellon University
#
#	HISTORY
#		12/96	Written by William F. Eddy (CMU)
#

PKG          = smregpar
PKG_MAKEBINS = $(CB)/smregpar

PKG_LIBS     = -lfmri -lmri -lpar -lbio -lacct -lmisc -lcrg $(LAPACK_LIBS) -lm

ALL_MAKEFILES= Makefile
CSOURCE= smregpar.c
DOCFILES= smregpar_help.help

include ../Makefile_pkg

LIBFILES= $L/libmri.a $L/libpar.a $L/libbio.a $L/libarray.a $L/libmisc.a \
	$L/libacct.a
HDRS= fmri.h mri.h par.h bio.h misc.h acct.h stdcrg.h

$O/smregpar.o: smregpar.c
	$(CC_RULE)

$O/smregpar_help.o: smregpar_help.help
	$(HELP_RULE)

$(CB)/smregpar: $O/smregpar.o $O/smregpar_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

releaseprep:
	echo "no release prep from " `pwd`

