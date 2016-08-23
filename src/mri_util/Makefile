#
#	Makefile for mri utilities
#
#	Copyright (c) 1997 Department of Statistics, Carnegie Mellon University
#
#	HISTORY
#		9/97	Written by Joel S. Welling (CMU)
#

PKG          = mri_util
PKG_MAKEBINS = $(CB)/mri_complex_to_scalar \
	$(CB)/mri_rpn_math $(CB)/mri_splus_filter $(CB)/mri_subset \
	$(CB)/mri_glm $(CB)/mri_interp \
	$(CB)/mri_from_pdigm $(CB)/mri_subsample $(CB)/mri_scan_fold \
	$(CB)/mri_change_dims $(CB)/mri_counter $(CB)/mri_type_convert \
	$(CB)/mri_pad $(CB)/mri_smooth $(CB)/mri_fft \
	$(CB)/mri_from_ascii $(CB)/mri_history $(CB)/mri_copy_chunk \
	$(CB)/mri_delete_chunk $(CB)/mri_paste $(CB)/mri_sort \
	$(CB)/mri_copy_dataset $(CB)/mri_destroy_dataset \
	$(CB)/mri_remap $(CB)/mri_printfield $(CB)/mri_permute \
	$(CB)/mri_setfield $(CB)/mri_matmult $(CB)/mri_esa \
	$(CB)/mri_resample $(CB)/mri_describe $(CB)/mri_svd \
	$(CB)/mri_kalman

PKG_LIBS     = -lfmri -ldcdf -lmri -lpar -lbio -lacct \
	-lcrg -lmisc $(LAPACK_LIBS) -lm

MAKEFILES= Makefile
CSOURCE= mri_complex_to_scalar.c mri_rpn_math.c mri_splus_filter.c \
	slave_splus.c slave_splus_tester.c \
	splus_binary_pipes.c mri_subset.c mri_glm.c \
	mri_interp.c mri_from_pdigm.c mri_subsample.c mri_scan_fold.c \
	mri_change_dims.c mri_counter.c mri_type_convert.c \
	mri_pad.c mri_paste.c mri_smooth.c mri_fft.c mri_from_ascii.c \
	mri_history.c mri_copy_chunk.c mri_delete_chunk.c mri_sort.c \
	mri_copy_dataset.c mri_destroy_dataset.c mri_remap.c \
	mri_printfield.c mri_permute.c tumble.c mri_setfield.c \
	mri_matmult.c mri_esa.c mri_resample.c mri_describe.c \
	mri_svd.c mri_kalman.c
HFILES= slave_splus.h tumble.h
DOCFILES= mri_complex_to_scalar_help.help mri_splus_filter_help.help \
	mri_rpn_math_help.help \
	mri_subset_help.help mri_glm_help.help mri_interp_help.help \
	mri_from_pdigm_help.help mri_subsample_help.help \
	mri_scan_fold_help.help mri_change_dims_help.help \
	mri_counter_help.help mri_type_convert_help.help \
	mri_pad_help.help mri_paste_help.help mri_smooth_help.help \
	mri_fft_help.help mri_from_ascii_help.help mri_history_help.help \
	mri_copy_chunk_help.help mri_delete_chunk_help.help \
	mri_sort_help.help mri_copy_dataset_help.help \
	mri_destroy_dataset_help.help mri_remap_help.help \
	mri_printfield_help.help mri_permute_help.help \
	mri_setfield_help.help mri_matmult_help.help \
	mri_esa_help.help mri_resample_help.help mri_describe_help.help \
	mri_svd_help.help mri_kalman_help.help

MISCFILES= bio_init.S

include ../Makefile.pkg

PKG_OBJS = $O/mri_splus_filter.o $O/slave_splus.o

LIBFILES= $L/libfmri.a $L/libmri.a $L/libpar.a $L/libbio.a \
	  $L/libarray.a $L/libmisc.a $L/libacct.a $L/libcrg.a

HDRS= fmri.h mri.h par.h bio.h misc.h acct.h stdcrg.h slave_splus.h

$(CB)/mri_kalman: $O/mri_kalman.o $O/mri_kalman_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mri_kalman.o: mri_kalman.c
	$(CC_RULE)

$O/mri_kalman_help.o: mri_kalman_help.help
	$(HELP_RULE)

$(CB)/mri_describe: $O/mri_describe.o $O/mri_describe_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mri_describe.o: mri_describe.c
	$(CC_RULE)

$O/mri_describe_help.o: mri_describe_help.help
	$(HELP_RULE)

$(CB)/mri_svd: $O/mri_svd.o $O/mri_svd_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mri_svd.o: mri_svd.c
	$(CC_RULE)

$O/mri_svd_help.o: mri_svd_help.help
	$(HELP_RULE)

$(CB)/mri_resample: $O/mri_resample.o $O/mri_resample_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mri_resample.o: mri_resample.c
	$(CC_RULE)

$O/mri_resample_help.o: mri_resample_help.help
	$(HELP_RULE)

$(CB)/mri_esa: $O/mri_esa.o $O/mri_esa_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mri_esa.o: mri_esa.c
	$(CC_RULE)

$O/mri_esa_help.o: mri_esa_help.help
	$(HELP_RULE)

$(CB)/mri_matmult: $O/mri_matmult.o $O/mri_matmult_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mri_matmult.o: mri_matmult.c
	$(CC_RULE)

$O/mri_matmult_help.o: mri_matmult_help.help
	$(HELP_RULE)

$(CB)/mri_printfield: $O/mri_printfield.o $O/mri_printfield_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mri_printfield.o: mri_printfield.c
	$(CC_RULE)

$O/mri_printfield_help.o: mri_printfield_help.help
	$(HELP_RULE)

$(CB)/mri_setfield: $O/mri_setfield.o $O/mri_setfield_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mri_setfield.o: mri_setfield.c
	$(CC_RULE)

$O/mri_setfield_help.o: mri_setfield_help.help
	$(HELP_RULE)

$(CB)/mri_remap: $O/mri_remap.o $O/mri_remap_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mri_remap.o: mri_remap.c
	$(CC_RULE)

$O/mri_remap_help.o: mri_remap_help.help
	$(HELP_RULE)

$(CB)/mri_copy_chunk: $O/mri_copy_chunk.o $O/mri_copy_chunk_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mri_copy_chunk.o: mri_copy_chunk.c
	$(CC_RULE)

$O/mri_copy_chunk_help.o: mri_copy_chunk_help.help
	$(HELP_RULE)

$(CB)/mri_delete_chunk: $O/mri_delete_chunk.o $O/mri_delete_chunk_help.o \
		$(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mri_delete_chunk.o: mri_delete_chunk.c
	$(CC_RULE)

$O/mri_delete_chunk_help.o: mri_delete_chunk_help.help
	$(HELP_RULE)

$O/mri_permute.o: mri_permute.c
	$(CC_RULE)

$O/tumble.o: tumble.c
	$(CC_RULE)

$O/mri_permute_help.o: mri_permute_help.help
	$(HELP_RULE)

$(CB)/mri_permute: $O/mri_permute.o $O/tumble.o $O/mri_permute_help.o \
		$(LIBFILES)
	@echo %%%% Linking mri_permute %%%%
	@$(LD) $(LFLAGS) -o $B/$(@F) $O/mri_permute.o $O/tumble.o \
	       $O/mri_permute_help.o $(LIBS)

$O/mri_glm.o: mri_glm.c
	$(CC_RULE)

$O/mri_glm_help.o: mri_glm_help.help
	$(HELP_RULE)

$(CB)/mri_glm: $O/mri_glm.o $O/mri_glm_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mri_subset.o: mri_subset.c
	$(CC_RULE)

$O/mri_subset_help.o: mri_subset_help.help
	$(HELP_RULE)

$(CB)/mri_subset: $O/mri_subset.o \
		$O/mri_subset_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mri_smooth.o: mri_smooth.c
	$(CC_RULE)

$O/mri_smooth_help.o: mri_smooth_help.help
	$(HELP_RULE)

$(CB)/mri_smooth: $O/mri_smooth.o \
		$O/mri_smooth_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mri_fft.o: mri_fft.c
	$(CC_RULE)

$O/mri_fft_help.o: mri_fft_help.help
	$(HELP_RULE)

$(CB)/mri_fft: $O/mri_fft.o \
		$O/mri_fft_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mri_from_ascii.o: mri_from_ascii.c
	$(CC_RULE)

$O/mri_from_ascii_help.o: mri_from_ascii_help.help
	$(HELP_RULE)

$(CB)/mri_from_ascii: $O/mri_from_ascii.o \
		$O/mri_from_ascii_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mri_history.o: mri_history.c
	$(CC_RULE)

$O/mri_history_help.o: mri_history_help.help
	$(HELP_RULE)

$(CB)/mri_history: $O/mri_history.o \
		$O/mri_history_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mri_interp.o: mri_interp.c
	$(CC_RULE)

$O/mri_interp_help.o: mri_interp_help.help
	$(HELP_RULE)

$(CB)/mri_interp: $O/mri_interp.o \
		$O/mri_interp_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mri_change_dims.o: mri_change_dims.c
	$(CC_RULE)

$O/mri_change_dims_help.o: mri_change_dims_help.help
	$(HELP_RULE)

$(CB)/mri_change_dims: $O/mri_change_dims.o \
		$O/mri_change_dims_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mri_counter.o: mri_counter.c
	$(CC_RULE)

$O/mri_counter_help.o: mri_counter_help.help
	$(HELP_RULE)

$(CB)/mri_counter: $O/mri_counter.o $O/mri_counter_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mri_type_convert.o: mri_type_convert.c
	$(CC_RULE)

$O/mri_type_convert_help.o: mri_type_convert_help.help
	$(HELP_RULE)

$(CB)/mri_type_convert: $O/mri_type_convert.o $O/mri_type_convert_help.o \
		$(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mri_subsample.o: mri_subsample.c
	$(CC_RULE)

$O/mri_subsample_help.o: mri_subsample_help.help
	$(HELP_RULE)

$(CB)/mri_subsample: $O/mri_subsample.o \
		$O/mri_subsample_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mri_pad.o: mri_pad.c
	$(CC_RULE)

$O/mri_pad_help.o: mri_pad_help.help
	$(HELP_RULE)

$(CB)/mri_pad: $O/mri_pad.o $O/mri_pad_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mri_sort.o: mri_sort.c
	$(CC_RULE)

$O/mri_sort_help.o: mri_sort_help.help
	$(HELP_RULE)

$(CB)/mri_sort: $O/mri_sort.o $O/mri_sort_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mri_paste.o: mri_paste.c
	$(CC_RULE)

$O/mri_paste_help.o: mri_paste_help.help
	$(HELP_RULE)

$(CB)/mri_paste: $O/mri_paste.o $O/mri_paste_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mri_scan_fold.o: mri_scan_fold.c
	$(CC_RULE)

$O/mri_scan_fold_help.o: mri_scan_fold_help.help
	$(HELP_RULE)

$(CB)/mri_scan_fold: $O/mri_scan_fold.o \
		$O/mri_scan_fold_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mri_from_pdigm.o: mri_from_pdigm.c
	$(CC_RULE)

$O/mri_from_pdigm_help.o: mri_from_pdigm_help.help
	$(HELP_RULE)

$(CB)/mri_from_pdigm: $O/mri_from_pdigm.o \
		$O/mri_from_pdigm_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mri_complex_to_scalar.o: mri_complex_to_scalar.c
	$(CC_RULE)

$O/mri_rpn_math.o: mri_rpn_math.c
	$(CC_RULE)

$O/mri_complex_to_scalar_help.o: mri_complex_to_scalar_help.help
	$(HELP_RULE)

$(CB)/mri_complex_to_scalar: $O/mri_complex_to_scalar.o \
		$O/mri_complex_to_scalar_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mri_rpn_math_help.o: mri_rpn_math_help.help
	$(HELP_RULE)

$(CB)/mri_rpn_math: $O/mri_rpn_math.o $O/mri_rpn_math_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$(CB)/mri_destroy_dataset: $O/mri_destroy_dataset.o \
		$O/mri_destroy_dataset_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mri_destroy_dataset.o: mri_destroy_dataset.c
	$(CC_RULE)

$O/mri_destroy_dataset_help.o: mri_destroy_dataset_help.help
	$(HELP_RULE)

$(CB)/mri_copy_dataset: $O/mri_copy_dataset.o \
		$O/mri_copy_dataset_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

$O/mri_copy_dataset.o: mri_copy_dataset.c
	$(CC_RULE)

$O/mri_copy_dataset_help.o: mri_copy_dataset_help.help
	$(HELP_RULE)

$O/mri_splus_filter_help.o: mri_splus_filter_help.help
	$(HELP_RULE)

$O/mri_splus_filter.o: mri_splus_filter.c
	$(CC_RULE)

$O/slave_splus.o: slave_splus.c
	$(CC_RULE)

$(CB)/mri_splus_filter: $O/mri_splus_filter.o $O/slave_splus.o \
		$O/mri_splus_filter_help.o $(LIBFILES)
	$(MULTI_HELP_LD)

releaseprep:
	echo "no release prep from " `pwd`

