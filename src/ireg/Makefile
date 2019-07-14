#
#	Makefile for ireg
#
#	Copyright (c) 1996 Department of Statistics, Carnegie Mellon University
#
#	HISTORY
#		12/96	Written by William F. Eddy (CMU)
#

PKG          = ireg
PKG_MAKEBINS = $(CB)/ireg 

PKG_LIBS     = -lfmri -lmri -lpar -lbio -lacct -lmisc -lcrg $(LAPACK_LIBS) -lm

ALL_MAKEFILES= Makefile
CSOURCE= ireg.c
DOCFILES= ireg_help.help

include ../Makefile_pkg

LIBFILES= $L/libmri.a $L/libpar.a $L/libbio.a $L/libarray.a $L/libmisc.a \
	$L/libacct.a $L/libfmri.a
HDRS= fmri.h mri.h par.h bio.h misc.h acct.h stdcrg.h

$O/ireg.o: ireg.c
	$(CC_RULE)

$O/ireg_help.o: ireg_help.help
	$(HELP_RULE)

$(CB)/ireg: $O/ireg.o $O/ireg_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$(CB)/ireg_chirp: $O/ireg.o $O/ireg_help.o $O/../fmri/chirprot.o $(LIBFILES)
	$(LD) $(LFLAGS) -o $B/$(@F) $O/ireg.o $O/ireg_help.o \
		$O/../fmri/chirprot.o $(LIBS)

releaseprep:
	echo "no release prep from " `pwd`

