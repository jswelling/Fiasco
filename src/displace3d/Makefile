#
#	Makefile for displace3d
#
#	Copyright (c) 1996 Department of Statistics, Carnegie Mellon University
#
#	HISTORY
#		12/96	Written by William F. Eddy (CMU)
#                5/00   Modified for 3D by Joel Welling (CMU)
#

PKG          = displace3d
PKG_MAKEBINS = $(CB)/displace3d

PKG_LIBS     = -lfmri -lmri -lpar -lbio -lacct -lmisc -lcrg $(LAPACK_LIBS) -lm

LIBFILES= $L/libmri.a $L/libpar.a $L/libbio.a $L/libarray.a $L/libmisc.a \
	$L/libacct.a
HDRS= fmri.h mri.h par.h bio.h misc.h acct.h stdcrg.h

ALL_MAKEFILES= Makefile
CSOURCE= displace3d.c
DOCFILES= displace3d_help.help

include ../Makefile_pkg

$O/displace3d.o: displace3d.c
	$(CC_RULE)

$O/displace3d_help.o: displace3d_help.help
	$(HELP_RULE)

$(CB)/displace3d: $O/displace3d.o $O/displace3d_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

releaseprep:
	echo "no release prep from " `pwd`

