#
#	Makefile for gift format converter
#
#	Copyright (c) 1996  Pittsburgh Supercomputing Center
#
#	HISTORY
#		9/96	Written by Greg Hood (PSC)
#

PKG          = gift
PKG_LIBS     = -lmri -lfmri -lbio -lcrg -lmisc $(LAPACK_LIBS) -lm
PKG_MAKEBINS = $(CB)/gift
PKG_OBJS     = $O/gift.o $O/fiasco_read.o $O/fiasco_write.o $O/pgh_read.o \
	$O/pgh_write.o $O/analyze_read.o $O/analyze_write.o \
	$O/spiral_sl_read.o $O/spiral_sl_write.o 

MAKEFILES= Makefile
CSOURCE= analyze_read.c analyze_write.c fiasco_read.c fiasco_write.c \
	gift.c pgh_read.c pgh_write.c spiral_sl_read.c \
	spiral_sl_write.c
HFILES= analyze.h fiasco.h gift.h pgh.h spiral.h spiral_sl.h
DOCFILES= README gift_help.help specs.txt

include ../Makefile.pkg

HFILES = gift.h fiasco.h analyze.h spiral_sl.h pgh.h spiral.h

$(CB)/gift: $(PKG_OBJS) $O/gift_help.o
	$(MULTI_HELP_LD)

$O/gift_help.o: gift_help.help
	$(HELP_RULE)

$O/gift.o: gift.c
	$(CC_RULE)

$O/fiasco_read.o: fiasco_read.c
	$(CC_RULE)

$O/fiasco_write.o: fiasco_write.c
	$(CC_RULE)

$O/pgh_read.o: pgh_read.c
	$(CC_RULE)

$O/pgh_write.o: pgh_write.c
	$(CC_RULE)

$O/analyze_read.o: analyze_read.c
	$(CC_RULE)

$O/analyze_write.o: analyze_write.c
	$(CC_RULE)

$O/spiral_sl_read.o: spiral_sl_read.c
	$(CC_RULE)

$O/spiral_sl_write.o: spiral_sl_write.c
	$(CC_RULE)

releaseprep:
	echo "no release prep from " `pwd`

