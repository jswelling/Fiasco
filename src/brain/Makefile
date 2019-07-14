#
#	Makefile for Chris Genovese "brain" reconstruction programs
#
#	Copyright (c) 1996  Pittsburgh Supercomputing Center
#
#	HISTORY
#		1/96	Written by Greg Hood (PSC)
#		8/96	Made architecture-independent (Greg Hood, PSC)
#

PKG          = brain
PKG_LIBS     = -lfmri -lcrg -lmri -lpar -lbio -lacct -lmisc \
	-larray $(LAPACK_LIBS) -lm
PKG_MAKEBINS = $(CB)/vmpfx $(CB)/generate_vmpfx_in $(CB)/byte_extract \
	$(CB)/smooth_missing

ALL_MAKEFILES= Makefile
CSOURCE= addendum.c generate_vmpfx_in.c vmpfx.c amoeba.c info.c \
	drift.c lbfgs.c qnewton-lbfgs.c tmperror.c byte_extract.c \
	smooth_missing.c
HFILES= lapack.h minqnb.h tmperror.h
DOCFILES= vmpfx_help.help generate_vmpfx_in_help.help byte_extract_help.help \
	smooth_missing_help.help README
MISCFILES= vmpfx_proto.t vmpfx_null_proto.t ramps.txt

include ../Makefile_pkg

# Available conditional compiler options:
#  Really present:
#
#    Optimization: no flags (compiles using AMOEBA and QUASINEWTON routines)
#    Model:        NON_ADDITIVE_PASTE - does not look fully implemented
#    Computation:  NO_INLINE_DRIFTPROF - probably should never be on
#
#  Mentioned by Chris Genovese:
#
#    Optimization:   AMOEBA, POWELL, QUASINEWTON
#    Model:          SMOOTH_BELLS, ADDITIVE_PASTE
#    Computation:    COMPUTE_SE, SPLUS_COMPATIBLE
#
# OPTOBJ should be set according to the Optimization flag used
# as shown here.
#

# Voxelwise Maximum Posterior, Fixed Effects Model

OPTOBJ= $O/amoeba.o $O/qnewton-lbfgs.o $O/lbfgs.o
PKG_OBJS     = $O/vmpfx.o $O/drift.o $O/info.o $(OPTOBJ) \
		$O/tmperror.o $(LIBFILES)

CFLAGS_T1 = $(CFLAGS:-O3=-g)
CFLAGS_T2 = $(CFLAGS_T1:-O2=-g)
CFLAGS_T3 = $(CFLAGS_T2:-O1=-g)
CFLAGS_NOOPT = $(CFLAGS_T3:-O=-g)

#Controls how often to output diagnostic data during processing
MONITOR_CNT= 1024

$O/vmpfx.o: vmpfx.c
	@echo %%%% Compiling $(@F:.o=.c)
	@$(CC) -c $(CFLAGS) -DMONITOR=$(MONITOR_CNT) -o $@ $(@F:.o=.c)

$O/vmpfx_help.o: vmpfx_help.help
	$(HELP_RULE)

$O/drift.o: drift.c
	$(CC_RULE)

$O/info.o: info.c
	$(CC_RULE)

$O/amoeba.o: amoeba.c
	$(CC_RULE)

$O/qnewton-lbfgs.o: qnewton-lbfgs.c
	$(CC_RULE)

$O/tmperror.o: tmperror.c
	$(CC_RULE)

$O/lbfgs.o: lbfgs.c
	@echo %%%% Compiling $(@F:.o=.c)
	@$(CC) -c $(CFLAGS_NOOPT) -o $@ $(@F:.o=.c)

$(CB)/vmpfx: $O/vmpfx.o $O/vmpfx_help.o $O/drift.o $O/info.o $(OPTOBJ) \
		$O/tmperror.o $(LIBFILES)
	$(MULTI_HELP_LD)

$O/generate_vmpfx_in.o: generate_vmpfx_in.c
	$(CC_RULE)

$O/generate_vmpfx_in_help.o: generate_vmpfx_in_help.help
	$(HELP_RULE)

$(CB)/generate_vmpfx_in: $O/generate_vmpfx_in.o $O/generate_vmpfx_in_help.o \
		$(LIBFILES)
	$(SINGLE_HELP_LD)

$O/smooth_missing.o: smooth_missing.c
	$(CC_RULE)

$O/smooth_missing_help.o: smooth_missing_help.help
	$(HELP_RULE)

$(CB)/smooth_missing: $O/smooth_missing.o $O/smooth_missing_help.o ${LIBFILES}
	$(SINGLE_HELP_LD)

$O/byte_extract.o: byte_extract.c
	$(CC_RULE)

$O/byte_extract_help.o: byte_extract_help.help
	$(HELP_RULE)

$(CB)/byte_extract: $O/byte_extract.o $O/byte_extract_help.o \
		$(LIBFILES)
	$(SINGLE_HELP_LD)

releaseprep:
	echo "no release prep from " `pwd`

# Other objects not supported in this release:
#
# Voxelwise Full Posterior, Fixed Effects Model
#
#vfpfx: vfpfx.o drift.o $(LIBFILES)
#	$(MULTI_HELP_LD)
#
# Query Tool
#
#query: query.o $(LIBFILES)
#	$(SINGLE_HELP_LD)


