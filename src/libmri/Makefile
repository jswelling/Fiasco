#
#	Makefile for the mri interface library
#
#	Copyright (c) 1996  Pittsburgh Supercomputing Center
#
#	HISTORY
#		3/96	Written by Greg Hood (PSC)
#		9/96	Integrated Makefile with other FMRI software
#

PKG          = libmri
PKG_EXPORTS  = mri.h
PKG_MAKELIBS = $L/libmri.a
PKG_MAKEBINS = 

# other currently inactive targets for PKG_MAKEBINS:
# $(CB)/create $(CB)/mean $(CB)/imean $(CB)/endian $(CB)/single 
# $(CB)/msplit $(CB)/mcopy $(CB)/complex $(CB)/import $(CB)/mpull 
# $(CB)/mpush $(CB)/halve

PKG_LIBS     = -lmri -lbio -lcrg -lmisc -lm

MAKEFILES= Makefile
CSOURCE= complex.c create.c endian.c halve.c imean.c import.c libmri.c \
	mcopy.c mean.c mpull.c mpush.c msplit.c single.c 
HFILES= mri.h
DOCFILES= README mri-c.doc mri-pgh.doc ref.doc 

include ../Makefile.pkg

$L/libmri.a: $O/libmri.o
	@echo "%%%% Building $(@F) %%%%"
	@$(AR) $(ARFLAGS) $L/libmri.a $O/libmri.o
	@$(RANLIB) $L/libmri.a

$O/libmri.o: libmri.c
	$(CC_RULE)

$(CB)/create: $O/create.o $L/libmri.a
	$(SINGLE_LD)

$O/create.o: create.c
	$(CC_RULE)

$(CB)/mean: $O/mean.o $L/libmri.a
	$(SINGLE_LD)

$O/mean.o: mean.c
	$(CC_RULE)

$(CB)/imean: $O/imean.o $L/libmri.a
	$(SINGLE_LD)

$O/imean.o: imean.c
	$(CC_RULE)

$(CB)/endian: $O/endian.o $L/libmri.a
	$(SINGLE_LD)

$O/endian.o: endian.c
	$(CC_RULE)

$(CB)/single: $O/single.o $L/libmri.a
	$(SINGLE_LD)

$O/single.o: single.c
	$(CC_RULE)

$(CB)/msplit: $O/msplit.o $L/libmri.a
	$(SINGLE_LD)

$O/msplit.o: msplit.c
	$(CC_RULE)

$(CB)/mcopy: $O/mcopy.o $L/libmri.a
	$(SINGLE_LD)

$O/mcopy.o: mcopy.c
	$(CC_RULE)

$(CB)/complex: $O/complex.o $L/libmri.a
	$(SINGLE_LD)

$O/complex.o: complex.c
	$(CC_RULE)

$(CB)/import: $O/import.o $L/libmri.a
	$(SINGLE_LD)

$O/import.o: import.c
	$(CC_RULE)

$(CB)/mpull: $O/mpull.o $L/libmri.a
	$(SINGLE_LD)

$O/mpull.o: mpull.c
	$(CC_RULE)

$(CB)/mpush: $O/mpush.o $L/libmri.a
	$(SINGLE_LD)

$O/mpush.o: mpush.c
	$(CC_RULE)

$(CB)/halve: $O/halve.o $L/libmri.a
	$(SINGLE_LD)

$O/halve.o: halve.c
	$(CC_RULE)

releaseprep:
	echo "no release prep from " `pwd`

