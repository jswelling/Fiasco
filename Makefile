#
#	Top-level Makefile for FMRI software
#
#	Copyright (c) 1996 Pittsburgh Supercomputing Center
#	Written by Greg Hood, PSC, Sept. 1996
#

SHELL        = /bin/sh
ARCH         := $(shell src/fiat_scripts/fiasco_getarch.csh)
COMPILE_ARCH = $(ARCH)

mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
current_dir := $(notdir $(patsubst %/,%,$(dir $(mkfile_path))))
FMRI         = $(dir $(mkfile_path))
DEFINES      = ARCH=$(ARCH) COMPILE_ARCH=$(COMPILE_ARCH) FMRI="$(FMRI)"

include Makefile_common
