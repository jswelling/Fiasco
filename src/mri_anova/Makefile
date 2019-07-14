#
#	Makefile for mri_anova
#
#	Copyright (c) 1998 Department of Statistics, Carnegie Mellon University
#
#	HISTORY
#		9/98	Written by Joel S. Welling (CMU)
#

PKG          = mri_anova
PKG_MAKEBINS = $(CB)/mri_anova 

PKG_LIBS     = -lfmri -lmri -lpar -lbio -lacct -lmisc -lcrg $(LAPACK_LIBS) -lm

ALL_MAKEFILES= Makefile
CSOURCE= mri_anova.c 
HFILES= 
DOCFILES= mri_anova_help.help
MISCFILES= 

include ../Makefile_pkg

PKG_OBJS = 

LIBFILES= $L/libfmri.a $L/libmri.a $L/libpar.a $L/libbio.a $L/libarray.a \
	$L/libmisc.a $L/libacct.a
HDRS= fmri.h mri.h par.h bio.h misc.h acct.h stdcrg.h slave_splus.h

$O/mri_anova.o: mri_anova.c
	$(CC_RULE)

$O/mri_anova_help.o: mri_anova_help.help
	$(HELP_RULE)

$(CB)/mri_anova: $O/mri_anova.o $O/mri_anova_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

releaseprep:
	echo "no release prep from " `pwd`

