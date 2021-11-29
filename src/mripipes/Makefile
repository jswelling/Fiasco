#
#	Makefile for the mri interface library
#
#	Copyright (c) 2004  Pittsburgh Supercomputing Center
#
#	HISTORY
#		3/96	Written by Greg Hood (PSC)
#		9/96	Integrated Makefile with other FMRI software
#

PKG          = mripipes
PKG_EXPORTS  = mripipes.h
PKG_MAKELIBS = $L/libmripipestools.a
PKG_MAKEBINS = $(CB)/_mripipes.$(SHR_EXT) $(CB)/mripipes.py \
	build_envs.bash

PKG_CFLAGS   = -I$(PYTHON_INCLUDE)
PKG_LIBS     = -lfmri -lmri -lbio -lcrg -lmisc -ldcdf $(LAPACK_LIBS) -lm

ALL_MAKEFILES= Makefile
CSOURCE= mripipes.c mripipes_wrap.c test_tool.c mri_file_input_tool.c \
	devnull_tool.c passthru_tool.c mri_file_output_tool.c \
	matmult_tool.c rpn_math_tool.c zero_src_tool.c \
	subset_tool.c pad_tool.c block_map_tool.c \
	func_2rows_unblocked_tool.c special_tool.c
HFILES= mripipes.h
DOCFILES=  
#
# mripipes.i is a SWIG input file to generate mripipes_wrap.c .  It
# won't normally be needed, but is included for completeness.
#
MISCFILES= mripipes.i mripipes.py mripipes_setup.py install_py_modules.bash

TOOLS_OBJ= $O/test_tool.o $O/mri_file_input_tool.o \
	$O/devnull_tool.o $O/passthru_tool.o $O/mri_file_output_tool.o \
	$O/matmult_tool.o $O/rpn_math_tool.o $O/zero_src_tool.o \
	$O/subset_tool.o $O/pad_tool.o $O/block_map_tool.o \
	$O/func_2rows_unblocked_tool.o $O/special_tool.o

.PHONY: build_envs.bash

include ../Makefile_pkg

$L/libmripipestools.a: $(TOOLS_OBJ)
	@echo "%%%% Building $(@F) %%%%"
	@$(AR) $(ARFLAGS) $L/libmripipestools.a $(TOOLS_OBJ)
	@$(RANLIB) $L/libmripipestools.a

$O/test_tool.o: test_tool.c
	$(CC_RULE)

$O/mri_file_input_tool.o: mri_file_input_tool.c
	$(CC_RULE)

$O/devnull_tool.o: devnull_tool.c
	$(CC_RULE)

$O/passthru_tool.o: passthru_tool.c
	$(CC_RULE)

$O/mri_file_output_tool.o: mri_file_output_tool.c
	$(CC_RULE)

$O/matmult_tool.o: matmult_tool.c
	$(CC_RULE)

$O/rpn_math_tool.o: rpn_math_tool.c
	$(CC_RULE)

$O/zero_src_tool.o: zero_src_tool.c
	$(CC_RULE)

$O/subset_tool.o: subset_tool.c
	$(CC_RULE)

$O/pad_tool.o: pad_tool.c
	$(CC_RULE)

$O/block_map_tool.o: block_map_tool.c
	$(CC_RULE)

$O/func_2rows_unblocked_tool.o: func_2rows_unblocked_tool.c
	$(CC_RULE)

$O/special_tool.o: special_tool.c
	$(CC_RULE)

$(CB)/_mripipes.$(SHR_EXT) $(CB)mripipes.py: \
$(CSOURCE) $(HFILES) $L/libmisc.a
	@echo "%%%% Building python module ${@F} %%%%"
	@env CFLAGS='$(CFLAGS)' LFLAGS='-L$(L) $(PKG_LIBS)' \
		SRCFILES='$(CSOURCE)' \
		python mripipes_setup.py build --build-lib $(CB)

mripipes_wrap.c mripipes.py: mripipes.i
	@echo "%%%% using SWIG to regenerate python mripipes wrapper! %%%%"
	$(SWIG) -python mripipes.i
	cp mripipes.py $(CB)

build_envs.bash:
	# This target captures the build environment in a sourceable file
	echo CFLAGS='"'$(FMRI_CFLAGS) -I$(PYTHON_INCLUDE)'"' \
		> build_envs.bash
	echo LFLAGS='"'-L$(L) $(FMRI_LFLAGS) $(FMRI_LIBS) $(PKG_LIBS)'"' \
		>> build_envs.bash
	echo SRCFILES='"'$(CSOURCE)'"' >> build_envs.bash

releaseprep:
	echo "no release prep from " `pwd`

