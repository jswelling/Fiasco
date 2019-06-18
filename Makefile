#
#	Top-level Makefile for FMRI software
#
#	Copyright (c) 1996 Pittsburgh Supercomputing Center
#	Written by Greg Hood, PSC, Sept. 1996
#

SHELL        = /bin/sh
ARCH         := $(shell src/fiat_scripts/fiasco_getarch.csh)
COMPILE_ARCH = $(ARCH)

include Makefile.common

FMRI         = $(TOPDIR)
DEFINES      = ARCH=$(ARCH) COMPILE_ARCH=$(COMPILE_ARCH) FMRI="$(FMRI)"

all depend .DEFAULT:
	$(MAKE) $(DEFINES) -f Makefile.common $@ 2>&1 | \
		grep -v 'recipe for target'

