#
#	Makefile for Mark Fitzgerald's fmri library
#
#	Copyright (c) 1996  Pittsburgh Supercomputing Center
#
#	HISTORY
#		9/96	Written by Greg Hood (PSC)
#

PKG          = fmri
PKG_EXPORTS  = fmri.h glm.h smoother.h parsesplit.h quaternion.h fshrot3d.h \
	linrot3d.h history.h filetypes.h frozen_header_info.h \
	frozen_header_info_cnv4.h frozen_header_info_lx2.h \
	frozen_header_info_prelx.h frozen_header_info_excite.h \
	windaq_header_info.h kvhash.h \
	siemens_kspace_header_info.h optimizer.h linwarp.h rpn_engine.h \
	entropy.h fexceptions.h closest_warp.h spline.h interpolator.h \
	fiat.h slicepattern.h mriu.h kalmanfilter.h
PKG_MAKELIBS = $L/libfmri.a
PKG_MAKEBINS = $(CB)/smoother_tester $(CB)/quat_tester \
	$(CB)/quaternion.py $(CB)/_quaternion.$(SHR_EXT) \
	$(CB)/optimizer_tester $(CB)/exception_tester $(CB)/fft3d_tester \
	$(CB)/slicepattern_tester $(CB)/glm_tester \
	$(CB)/fiasco_numpy.py $(CB)/_fiasco_numpy.$(SHR_EXT) \

PKG_LIBS     = -lfmri -lmri -lpar -lbio -lacct -lmisc -lcrg $(LAPACK_LIBS) -lm
PKG_CFLAGS   = -I$(PYTHON_INCLUDE)

MAKEFILES= Makefile
CSOURCE= fft2d.c fmri.c fshrot.c polyt.c glm.c glm_irls.c smoother.c \
	smoother_tester.c parsesplit.c fft3d.c quaternion.c quat_tester.c \
	fshrot3d.c linrot3d.c praxis.c history.c filetypes.c kvhash.c \
	bvls.c fmin.c quaternion_wrap.c optimizer.c optimizer_tester.c \
	linwarp.c rpn_engine.c entropy.c fexceptions.c exception_tester.c \
	closest_warp.c spline.c interpolator.c fft3d_tester.c slicepattern.c \
	slicepattern_tester.c mriu.c fiasco_numpy_wrap.c  glm_tester.c \
	kalmanfilter.c
HFILES= fmri.h lapack.h glm.h smoother.h parsesplit.h quaternion.h \
	fshrot3d.h linrot3d.h history.h frozen_header_info.h \
	frozen_header_info_cnv4.h frozen_header_info_lx2.h \
	frozen_header_info_prelx.h frozen_header_info_excite.h \
	siemens_kspace_header_info.h \
	windaq_header_info.h filetypes.h kvhash.h optimizer.h linwarp.h \
	rpn_engine.h entropy.h fexceptions.h closest_warp.h mriu.h \
	spline.h interpolator.h fiat.h slicepattern.h kalmanfilter.h
DOCFILES= smoother_help.help fft2d_help.help fft3d_help.help \
	fshrot3d_help.help linrot3d_help.help praxis_help.help \
	nelmin_help.help coordsys_help.help fmin_help.help \
	linwarp_help.help entropy_help.help rpn_engine_help.help \
	closest_warp_help.help interpolator_help.help
#
# quaternion.i is a SWIG input file to generate quaternion.py and 
# quaternion_wrap.c .  It won't normally be needed, but is included 
# for completeness.  Same goes for fiasco_numpy.i.  That set of 
# utilities requires Python numpy, which isn't yet in the default
# Python build, so actual construction is commented out.
#
MISCFILES= quaternion.i quaternion_setup.py quaternion.py \
	numpy.i fiasco_numpy.i fiasco_numpy_setup.py fiasco_numpy.py \
	test_kalman.py test_interpolator.py test_logistic_regression.py \
	test_optimizer.py

FMRI_OBJ= $O/fmri.o $O/fft2d.o $O/fshrot.o $O/polyt.o $O/glm.o \
	$O/glm_irls.o $O/mriu.o \
	$O/smoother.o $O/parsesplit.o $O/fft3d.o $O/quaternion.o \
	$O/fshrot3d.o $O/linrot3d.o $O/praxis.o $O/history.o \
	$O/filetypes.o $O/kvhash.o $O/bvls.o $O/fmin.o $O/optimizer.o \
	$O/linwarp.o $O/rpn_engine.o $O/entropy.o $O/fexceptions.o \
	$O/closest_warp.o $O/spline.o $O/interpolator.o $O/slicepattern.o \
	$O/kalmanfilter.o

include ../Makefile.pkg

$L/libfmri.a: $(FMRI_OBJ)
	@echo "%%%% Building $(@F) %%%%"
	@$(AR) $(ARFLAGS) $L/libfmri.a $(FMRI_OBJ)
	@$(RANLIB) $L/libfmri.a

$O/kalmanfilter.o: kalmanfilter.c
	$(CC_RULE)

$O/mriu.o: mriu.c
	$(CC_RULE)

$O/slicepattern.o: slicepattern.c
	$(CC_RULE)

$O/slicepattern_tester.o: slicepattern_tester.c
	$(CC_RULE)

$(CB)/slicepattern_tester: $O/slicepattern_tester.o $L/libfmri.a
	$(SINGLE_LD)

$O/spline.o: spline.c
	$(CC_RULE)

$O/interpolator.o: interpolator.c
	$(CC_RULE)

$O/entropy.o: entropy.c
	$(CC_RULE)

$O/fexceptions.o: fexceptions.c
	$(CC_RULE)

$O/rpn_engine.o: rpn_engine.c
	$(CC_RULE)

$O/kvhash.o: kvhash.c
	$(CC_RULE)

$O/filetypes.o: filetypes.c
	$(CC_RULE)

$O/fmri.o: fmri.c
	$(CC_RULE)

$O/fft2d.o: fft2d.c
	$(CC_RULE)

$O/fft3d.o: fft3d.c
	$(CC_RULE)

$O/fshrot.o: fshrot.c
	$(CC_RULE)

$O/fshrot3d.o: fshrot3d.c
	$(CC_RULE)

$O/linrot3d.o: linrot3d.c
	$(CC_RULE)

$O/linwarp.o: linwarp.c
	$(CC_RULE)

$O/closest_warp.o: closest_warp.c
	$(CC_RULE)

$O/polyt.o: polyt.c
	$(CC_RULE)

$O/praxis.o: praxis.c
	$(CC_RULE)

$O/bvls.o: bvls.c
	$(CC_RULE)

$O/fmin.o: fmin.c
	$(CC_RULE)

$O/history.o: history.c
	$(CC_RULE)

$O/glm.o: glm.c
	$(CC_RULE)

$O/glm_irls.o: glm_irls.c
	$(CC_RULE)

$O/glm_tester.o: glm_tester.c 
	$(CC_RULE)

$(CB)/glm_tester: $O/glm_tester.o $L/libfmri.a $(LIBFILES)
	@echo %%%% Linking glm_tester %%%%
	@$(LD) $(LFLAGS) -o $B/$(@F) $O/glm_tester.o $(LIBS)

$O/smoother.o: smoother.c
	$(CC_RULE)

$O/smoother_tester.o: smoother_tester.c
	$(CC_RULE)

$(CB)/smoother_tester: $O/smoother_tester.o $L/libfmri.a
	$(SINGLE_LD)

$O/exception_tester.o: exception_tester.c
	$(CC_RULE)

$(CB)/exception_tester: $O/exception_tester.o $L/libfmri.a
	$(SINGLE_LD)

$O/optimizer.o: optimizer.c
	$(CC_RULE)

$O/optimizer_tester.o: optimizer_tester.c
	$(CC_RULE)

$(CB)/optimizer_tester: $O/optimizer_tester.o $L/libfmri.a
	$(SINGLE_LD)

$O/parsesplit.o: parsesplit.c
	$(CC_RULE)

$O/quaternion.o: quaternion.c
	$(CC_RULE)

$O/quat_tester.o: quat_tester.c
	$(CC_RULE)

$(CB)/quat_tester: $O/quat_tester.o $L/libfmri.a
	$(SINGLE_LD)

$O/fft3d_tester.o: fft3d_tester.c
	$(CC_RULE)

$(CB)/fft3d_tester: $O/fft3d_tester.o $L/libfmri.a
	$(SINGLE_LD)

$O/chirprot.o: chirprot.c
	$(CC_RULE)

quaternion_wrap.c quaternion.py: quaternion.i
	@echo '%%%% Using SWIG to regenerate python quaternion wrapper! ####'
	$(SWIG) -python quaternion.i

$(CB)/_quaternion.$(SHR_EXT) $(CB)/quaternion.py: \
		quaternion.py quaternion_setup.py quaternion_wrap.c \
		quaternion.c quaternion.h 
	@echo "%%%% Building python module ${@F} %%%%"
	@env CFLAGS='$(CFLAGS)' \
		LFLAGS='-L$(L) -lmisc $(LAPACK_LIBS)' \
		python quaternion_setup.py build --build-lib $(CB) 

fiasco_numpy_wrap.c fiasco_numpy.py: fiasco_numpy.i
	@echo '%%%% Using SWIG to regenerate python fiasco_numpy wrapper! ####'
	$(SWIG) -python fiasco_numpy.i

$(CB)/_fiasco_numpy.$(SHR_EXT) $(CB)/fiasco_numpy.py: \
		fiasco_numpy.py fiasco_numpy_setup.py fiasco_numpy_wrap.c \
		glm.c glm.h \
		optimizer.c optimizer.h $L/libfmri.a
	@echo "%%%% Building python module ${@F} %%%%"
	@env CFLAGS='$(CFLAGS)' \
		LFLAGS='-L$(L) -lmisc -lfmri $(LAPACK_LIBS)' \
		python fiasco_numpy_setup.py build --build-lib $(CB) 

releaseprep:
	echo "no release prep from " `pwd`

