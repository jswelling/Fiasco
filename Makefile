#
#	Top-level Makefile for FMRI software
#
#	Copyright (c) 1996 Pittsburgh Supercomputing Center
#	Written by Greg Hood, PSC, Sept. 1996
#

SHELL        = /bin/sh
FMRI         = `pwd`
ARCH         = `src/fiat_scripts/fiasco_getarch.csh`
COMPILE_ARCH = `src/fiat_scripts/fiasco_getarch.csh`
DEFINES      = ARCH=$(ARCH) COMPILE_ARCH=$(COMPILE_ARCH) FMRI="$(FMRI)"

all depend .DEFAULT:
	@$(MAKE) $(DEFINES) -f Makefile.$(ARCH) $@ 2>&1 | \
		grep -v "commands for target"

