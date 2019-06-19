#
# Makefile for libcrg library
#	Adapted for PSC by Greg Hood, Sept. 1996
# Author:              Christopher R. Genovese
# Last Modified:       January 1996
#

PKG          = libcrg
PKG_EXPORTS  = stdcrg.h slist.h helphelp.h
PKG_MAKELIBS = $L/libcrg.a
PKG_MAKEBINS = $(CB)/buildhelp

PKG_CFLAGS   = -DDIAGNOSTIC=1 -DDEBUG=0
PKG_LIBS     = -lmisc -lcrg
PKG_OBJS     = $O/buildhelp.o $O/slist.o $O/helphelp.o $O/buildhelp_help.o

MAKEFILES= Makefile macros.m4
CSOURCE= buildhelp.c clproc.c efunc.c help.c helpblank.c inpfproc.c \
	sorting.c slist.c helphelp.c
HFILES= stdcrg.h slist.h helphelp.h
DOCFILES= buildhelp_help.help

include ../Makefile.pkg

$L/libcrg.a: $O/sorting.o $O/inpfproc.o $O/clproc.o $O/help.o $O/efunc.o \
		$O/slist.o $O/helphelp.o
	@echo "%%%% Building $(@F) %%%%"
	@$(AR) $(ARFLAGS) $L/libcrg.a $O/sorting.o $O/inpfproc.o $O/clproc.o \
              $O/help.o $O/efunc.o $O/slist.o $O/helphelp.o
	@$(RANLIB) $L/libcrg.a

$O/sorting.o: sorting.c
	$(CC_RULE)

$O/inpfproc.o: inpfproc.c
	$(CC_RULE)

$O/clproc.o: clproc.c
	$(CC_RULE)

$O/help.o: help.c
	$(CC_RULE)

$O/efunc.o: efunc.c
	$(CC_RULE)

#
# buildhelp must be built in two steps- stage1 makes a version 
# with stub help routines, which is then used to *create* the help
# routines for the final (stage 2) buildhelp.  The stage1 version
# gets deleted after its job is done.
#

$O/buildhelp_help.o: buildhelp_help.help $(CB)/buildhelp_stage1
	@echo "%%%% Generating $(@F) using buildhelp_stage1 %%%%" ;
	@$(CB)/buildhelp_stage1 buildhelp_help.help $O/tmp
	@$(CC) -c $(CFLAGS) -o $O/buildhelp_help.o $O/tmpHelp.c
#	rm $O/tmpHelp.c

$(CB)/buildhelp_stage1: $O/buildhelp.o $O/helpblank.o $O/slist.o \
		$(PKG_MAKELIBS)
	@echo "%%%% Linking buildhelp_stage1 %%%%"
	@$(LD) $(LFLAGS) -o $(CB)/buildhelp_stage1 $O/buildhelp.o \
		$O/helpblank.o $O/slist.o $O/helphelp.o $(PKG_LIBS) -lm

$O/buildhelp.o: buildhelp.c
	$(CC_RULE)

$O/slist.o: slist.c
	$(CC_RULE)

$O/helphelp.o: helphelp.c
	$(CC_RULE)

$O/helpblank.o: helpblank.c
	$(CC_RULE)

$(CB)/buildhelp: $(CB)/buildhelp_stage1 $(PKG_OBJS) $(PKG_MAKELIBS)
	@-rm -f $(CB)/buildhelp_stage1
	@$(LD) $(LFLAGS) -o $B/$(@F) $(PKG_OBJS) $(PKG_LIBS) -lm

releaseprep:
	echo "no release prep from " `pwd`

