#
#	Makefile for FMRI utility routines
#
#	Copyright (c) 1996  Pittsburgh Supercomputing Center
#
#	HISTORY
#		1/96	Written by Greg Hood (PSC)
#		8/96	Made architecture-independent (Greg Hood, PSC)
#

PKG          = util
PKG_EXPORTS  = acct.h array.h bio.h mdbg.h misc.h par.h errors.h \
	       pulse.h rttraj.h
PKG_MAKELIBS = $L/libacct.a $L/libarray.a $L/libbio.a $L/libmdbg.a \
               $L/libmisc.a $L/libpar.a $L/libpulse.a $L/librttraj.a
PKG_MAKEBINS = 
PKG_LIBS     = -lcdf -lfmri -lmri -lpar -lbio -lacct \
	     -lcrg -lmisc -lrttraj $(LAPACK_LIBS) -lm

MAKEFILES= Makefile
CSOURCE= libacct.c libarray.c libbio.c libmdbg.c libmisc.c libpar.c \
	ptest.c libpulse.c librttraj.c
HFILES= acct.h array.h bio.h mdbg.h misc.h par.h errors.h pulse.h rttraj.h
DOCFILES= 
SCRIPTFILES= 

LIBFILES= $L/libmri.a $L/libpar.a $L/libbio.a $L/libarray.a $L/libmisc.a \
	$L/libacct.a $L/libpulse.a $L/librttraj.a

include ../Makefile.pkg

$L/libacct.a: $O/libacct.o
	@echo "%%%% Building $(@F) %%%%"
	@$(AR) $(ARFLAGS) $L/libacct.a $O/libacct.o
	@$(RANLIB) $L/libacct.a

$O/libacct.o: libacct.c
	$(CC_RULE)

$L/libarray.a: $O/libarray.o	
	@echo "%%%% Building $(@F) %%%%"
	@$(AR) $(ARFLAGS) $L/libarray.a $O/libarray.o	
	@$(RANLIB) $L/libarray.a

$O/libarray.o: libarray.c
	$(CC_RULE)

$L/libbio.a: $O/libbio.o
	@echo "%%%% Building $(@F) %%%%"
	@$(AR) $(ARFLAGS) $L/libbio.a $O/libbio.o
	@$(RANLIB) $L/libbio.a

$O/libbio.o: libbio.c
	$(CC_RULE)

$L/libmdbg.a: $O/libmdbg.o
	@echo "%%%% Building $(@F) %%%%"
	@$(AR) $(ARFLAGS) $L/libmdbg.a $O/libmdbg.o
	@$(RANLIB) $L/libmdbg.a

$O/libmdbg.o: libmdbg.c
	$(CC_RULE)

$L/libmisc.a: $O/libmisc.o
	@echo "%%%% Building $(@F) %%%%"
	@$(AR) $(ARFLAGS) $L/libmisc.a $O/libmisc.o
	@$(RANLIB) $L/libmisc.a

$O/libmisc.o: libmisc.c
	$(CC_RULE)

$L/libpar.a: $O/libpar.o
	@echo "%%%% Building $(@F) %%%%"
	@$(AR) $(ARFLAGS) $L/libpar.a $O/libpar.o
	@$(RANLIB) $L/libpar.a

$O/libpar.o: libpar.c
	$(CC_RULE)

$L/libpulse.a: $O/libpulse.o
	@echo "%%%% Building $(@F) %%%%"
	@$(AR) $(ARFLAGS) $L/libpulse.a $O/libpulse.o
	@$(RANLIB) $L/libpulse.a

$O/libpulse.o: libpulse.c
	$(CC_RULE)

$L/librttraj.a: $O/librttraj.o
	@echo "%%%% Building $(@F) %%%%"
	@$(AR) $(ARFLAGS) $L/librttraj.a $O/librttraj.o
	@$(RANLIB) $L/librttraj.a

$O/librttraj.o: librttraj.c
	$(CC_RULE)

releaseprep:
	echo "no release prep from " `pwd`

