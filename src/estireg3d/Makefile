#
#	Makefile for estireg3d
#
#	Copyright (c) 1996 Department of Statistics, Carnegie Mellon University
#
#	HISTORY
#		12/96	Written by William F. Eddy (CMU)
#

PKG          = estireg3d
PKG_MAKEBINS = $(CB)/estireg3d $(CB)/estiwarp

PKG_LIBS     = -lfmri -lmri -lpar -lbio -lacct -lmisc -lcrg $(LAPACK_LIBS) -lm

ALL_MAKEFILES= Makefile
CSOURCE= estireg3d.c estiwarp.c algorithm.c estireg_utils.c
HFILES= algorithm.h estireg_utils.h
DOCFILES= estireg3d_help.help estiwarp_help.help

include ../Makefile_pkg

LIBFILES= $L/libmri.a $L/libpar.a $L/libbio.a $L/libarray.a $L/libmisc.a \
	$L/libacct.a $L/libfmri.a
HDRS= fmri.h mri.h par.h bio.h misc.h acct.h stdcrg.h

COMMONOBJS= $O/algorithm.o $O/estireg_utils.o

$O/algorithm.o: algorithm.c
	$(CC_RULE)

$O/estireg_utils.o: estireg_utils.c
	$(CC_RULE)

$O/estireg3d.o: estireg3d.c
	$(CC_RULE)

$O/estireg3d_help.o: estireg3d_help.help
	$(HELP_RULE)

$(CB)/estireg3d: $O/estireg3d.o $O/estireg3d_help.o $(COMMONOBJS) $(LIBFILES)
	@echo "%%%% Linking $(@F) %%%%"
	@$(LD) $(LFLAGS) -o $(CB)/$(@F) $O/estireg3d.o $(COMMONOBJS) \
		$O/estireg3d_help.o $(LIBS)

$O/estiwarp.o: estiwarp.c
	$(CC_RULE)

$O/estiwarp_help.o: estiwarp_help.help
	$(HELP_RULE)

$(CB)/estiwarp: $O/estiwarp.o $O/estiwarp_help.o $(COMMONOBJS) $(LIBFILES)
	@echo "%%%% Linking $(@F) %%%%"
	@$(LD) $(LFLAGS) -o $(CB)/$(@F) $O/estiwarp.o $(COMMONOBJS) \
		$O/estiwarp_help.o $(LIBS)

releaseprep:
	echo "no release prep from " `pwd`

