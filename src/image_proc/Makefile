#
#	Makefile for image proc
#
#	Copyright (c) 2002 Department of Statistics, Carnegie Mellon University
#
#	HISTORY
#		May 2002 Written by Jennifer Bakal (CMU)
#

PKG          = image_proc
PKG_MAKEBINS = $(CB)/mri_to_img $(CB)/matte $(CB)/colorize

PKG_LIBS     = -lfmri -lmri -lpar -lbio -lacct -lcrg -lmisc -lm \
	-ldcdf $(LAPACK_LIBS)

MAKEFILES= Makefile
CSOURCE= mri_to_img.c matte.c colorize.c
DOCFILES= mri_to_img_help.help matte_help.help colorize_help.help

include ../Makefile.pkg

LIBFILES= $L/libmri.a $L/libpar.a $L/libbio.a $L/libarray.a $L/libmisc.a \
	$L/libacct.a
HDRS= fmri.h mri.h par.h bio.h misc.h acct.h stdcrg.h

$O/mri_to_img.o: mri_to_img.c
	$(CC_RULE)

$O/mri_to_img_help.o: mri_to_img_help.help
	$(HELP_RULE)

$(CB)/mri_to_img: $O/mri_to_img.o $O/mri_to_img_help.o
	$(SINGLE_HELP_LD)

$O/matte.o: matte.c
	$(CC_RULE)

$O/matte_help.o: matte_help.help
	$(HELP_RULE)

$(CB)/matte: $O/matte.o $O/matte_help.o
	$(SINGLE_HELP_LD)

$O/colorize.o: colorize.c
	$(CC_RULE)

$O/colorize_help.o: colorize_help.help
	$(HELP_RULE)

$(CB)/colorize: $O/colorize.o $O/colorize_help.o
	$(SINGLE_HELP_LD)

releaseprep:
	echo "no release prep from " `pwd`

