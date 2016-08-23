#
#	Makefile for smregpar3d
#
#	Copyright (c) 1996 Department of Statistics, Carnegie Mellon University
#
#	HISTORY
#		12/96	Written by William F. Eddy (CMU)
#                4/00   modified smregpar -> smregpar3d Joel Welling (CMU)
#

PKG          = smregpar3d
PKG_MAKEBINS = $(CB)/smregpar3d

PKG_LIBS     = -lfmri -lmri -lpar -lbio -lacct -lmisc -lcrg $(LAPACK_LIBS) -lm

MAKEFILES= Makefile
CSOURCE= smregpar3d.c
DOCFILES= smregpar3d_help.help

include ../Makefile.pkg

LIBFILES= $L/libmri.a $L/libpar.a $L/libbio.a $L/libarray.a $L/libmisc.a \
	$L/libacct.a $L/libfmri.a
HDRS= fmri.h mri.h par.h bio.h misc.h acct.h stdcrg.h

$O/smregpar3d.o: smregpar3d.c
	$(CC_RULE)

$O/smregpar3d_help.o: smregpar3d_help.help
	$(HELP_RULE)

$(CB)/smregpar3d: $O/smregpar3d.o $O/smregpar3d_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

releaseprep:
	echo "no release prep from " `pwd`

