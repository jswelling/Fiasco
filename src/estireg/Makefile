#
#	Makefile for estireg
#
#	Copyright (c) 1996  Pittsburgh Supercomputing Center
#
#	HISTORY
#		9/96	Written by Greg Hood (PSC)
#

PKG          = estireg
PKG_LIBS     = -lfmri -lmri -lpar -lbio -lacct -larray -lcrg -lmisc \
	$(LAPACK_LIBS) -lm
PKG_OBJS     = $O/estireg.o $O/estireg_help.o
PKG_MAKEBINS = $(CB)/estireg

LIBFILES= $L/libmri.a $L/libpar.a $L/libbio.a $L/libarray.a $L/libmisc.a \
	$L/libacct.a
HDRS= mri.h par.h bio.h array.h misc.h acct.h

MAKEFILES= Makefile
CSOURCE= estireg.c 
DOCFILES= estireg_help.help

include ../Makefile.pkg

$O/estireg.o: estireg.c
	$(CC_RULE)

$O/estireg_help.o: estireg_help.help
	$(HELP_RULE)

$(CB)/estireg: $(PKG_OBJS) $(LIBFILES)
	$(MULTI_LD)

$O/cross.o: cross.c
	$(CC_RULE)

$O/cross_help.o: cross_help.help
	$(HELP_RULE)

$(CB)/cross: $O/cross.o $O/cross_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/extract.o: extract.c 
	$(CC_RULE)

$O/extract_help.o: extract_help.help
	$(HELP_RULE)

$(CB)/extract: $O/extract.o $O/extract_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

releaseprep:
	echo "no release prep from " `pwd`

