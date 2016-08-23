#
#	Makefile for displace
#
#	Copyright (c) 1996 Department of Statistics, Carnegie Mellon University
#
#	HISTORY
#		12/96	Written by William F. Eddy (CMU)
#

PKG          = displace
PKG_MAKEBINS = $(CB)/displace

PKG_LIBS     = -lfmri -lmri -lpar -lbio -lacct -lmisc -lcrg $(LAPACK_LIBS) -lm

LIBFILES= $L/libmri.a $L/libpar.a $L/libbio.a $L/libarray.a $L/libmisc.a \
	$L/libacct.a
HDRS= fmri.h mri.h par.h bio.h misc.h acct.h stdcrg.h

MAKEFILES= Makefile
CSOURCE= displace.c
DOCFILES= displace_help.help

include ../Makefile.pkg

$O/displace.o: displace.c
	$(CC_RULE)

$O/displace_help.o: displace_help.help
	$(HELP_RULE)

$(CB)/displace: $O/displace.o $O/displace_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

releaseprep:
	echo "no release prep from " `pwd`

