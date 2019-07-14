#
#	Makefile for reader
#
#	Copyright (c) 1996 Department of Statistics, Carnegie Mellon University
#
#	HISTORY
#		12/96	Written by William F. Eddy (CMU)
#

PKG          = reader
PKG_MAKEBINS = $(CB)/reader $(CB)/smartreader $(CB)/pghtoafni \
	$(CB)/pghtonifti

PKG_LIBS     = -lfmri -lmri -lpar -lbio -lacct -lmisc -larray -lcrg \
               -lrttraj $(LAPACK_LIBS) -lm
PKG_OBJS     = $O/reader.o $O/epi_correction.o $O/rcn.o \
		$O/nr_sub.o 
SMARTREADER_OBJS = $O/smartreader.o $O/smartreader_help.o \
	$O/base_reader.o $O/windaq_reader.o $O/raw_reader.o \
	$O/lx_reader_excite.o $O/lx_epibold_reader_excite.o \
	$O/lx_splx_reader_excite.o $O/lx_sf11_reader_excite.o \
	$O/lx_reader_cnv4.o $O/lx_epibold_reader_cnv4.o \
	$O/lx_splx_reader_cnv4.o $O/lx_sf11_reader_cnv4.o \
	$O/lx_reader_lx2.o $O/lx_epibold_reader_lx2.o \
	$O/lx_splx_reader_lx2.o $O/lx_sf11_reader_lx2.o \
	$O/lx_reader_prelx.o $O/lx_epibold_reader_prelx.o \
	$O/lx_splx_reader_prelx.o $O/lx_sf11_reader_prelx.o \
	$O/lx_2dfast_reader_prelx.o $O/lx_2dfast_reader_cnv4.o \
	$O/lx_2dfast_reader_lx2.o $O/lx_2dfast_reader_excite.o \
	$O/ram_reader.o $O/convert_reader.o \
	$O/lx_image_reader.o $O/multi_reader.o $O/afni_reader.o \
	$O/dicom_reader.o $O/dicom_transfer_syntax.o $O/dicom_uid_dict.o \
	$O/dicom_parser.o $O/desmith_reader.o \
	$O/siemens_kspace_reader.o $O/analyze_reader.o \
	$O/pghmri_reader.o $O/ushort_reader.o $O/fiff_reader.o \
	$O/son_reader.o $O/png_reader.o $O/tiff_reader.o $O/fits_reader.o \
	$O/nifti_reader.o \
	$O/wildcard.o $O/tables.o \
	$O/smart_utils.o $O/vec3.o 

PGHTOAFNI_OBJS = $O/pghtoafni.o $O/pghmri_reader.o $O/base_reader.o \
	       $O/convert_reader.o $O/pghtoafni_help.o $O/tables.o \
	       $O/smart_utils.o $O/vec3.o

PGHTONIFTI_OBJS = $O/pghtonifti.o $O/pghmri_reader.o $O/base_reader.o \
	       $O/convert_reader.o $O/pghtonifti_help.o $O/tables.o \
	       $O/smart_utils.o $O/vec3.o

ALL_MAKEFILES= Makefile
CSOURCE= reader.c epi_correction.c rcn.c nr_sub.c \
	smartreader.c lx_reader.c lx_epibold_reader.c \
	lx_splx_reader.c lx_sf11_reader.c windaq_reader.c \
	raw_reader.c base_reader.c ram_reader.c \
	convert_reader.c lx_image_reader.c desmith_reader.c \
	multi_reader.c afni_reader.c wildcard.c \
	dicom_reader.c dicom_transfer_syntax.c dicom_uid_dict.c \
	dicom_parser.c siemens_kspace_reader.c \
	analyze_reader.c pghmri_reader.c pghtoafni.c \
	tables.c smart_utils.c vec3.c lx_2dfast_reader.c \
	ushort_reader.c fiff_reader.c son_reader.c png_reader.c \
	tiff_reader.c nifti_reader.c fits_reader.c pghtonifti.c
HFILES= rcn.h control_data.h nr_sub.h smartreader.h dicom_dict.h \
	dicom_uid_dict.h dicom_transfer_syntax.h dicom_parser.h \
	afni_defs.h nifti1.h
DOCFILES= reader_help.help smartreader_help.help pghtoafni_help.help \
	pghtonifti_help.help

include ../Makefile_pkg

LIBFILES= $L/libmri.a $L/libfmri.a $L/libpar.a $L/libbio.a $L/libarray.a \
	  $L/libmisc.a $L/libacct.a $L/librttraj.a
HDRS= fmri.h mri.h par.h bio.h misc.h acct.h stdcrg.h

CC_FOR_EXCITE = @$(CC) -c $(CFLAGS) -DUSE_RDBM_EXCITE -o $@
CC_FOR_CNV4 = @$(CC) -c $(CFLAGS) -o $@
CC_FOR_LX2 = @$(CC) -c $(CFLAGS) -DUSE_RDBM_LX2 -o $@
CC_FOR_PRELX = @$(CC) -c $(CFLAGS) -DUSE_RDBM_PRELX -o $@

$O/reader.o: reader.c
	$(CC_RULE)

$O/reader_help.o: reader_help.help
	$(HELP_RULE)

$O/epi_correction.o: epi_correction.c
	$(CC_RULE)

$O/rcn.o: rcn.c
	$(CC_RULE)

$O/nr_sub.o: nr_sub.c
	$(CC_RULE)

$(CB)/reader: $(PKG_OBJS) $O/reader_help.o $(LIBFILES)
	$(MULTI_HELP_LD)

$O/pghtoafni.o: pghtoafni.c
	$(CC_RULE)

$O/pghtoafni_help.o: pghtoafni_help.help
	$(HELP_RULE)

$(CB)/pghtoafni: $(PGHTOAFNI_OBJS) $(LIBFILES)
	@echo "%%%% Linking $(@F) %%%%"
	@$(LD) $(LFLAGS) -o $B/$(@F) $(PGHTOAFNI_OBJS) $(LIBS)

$O/pghtonifti.o: pghtonifti.c
	$(CC_RULE)

$O/pghtonifti_help.o: pghtonifti_help.help
	$(HELP_RULE)

$(CB)/pghtonifti: $(PGHTONIFTI_OBJS) $(LIBFILES)
	@echo "%%%% Linking $(@F) %%%%"
	@$(LD) $(LFLAGS) -o $B/$(@F) $(PGHTONIFTI_OBJS) $(LIBS)

$O/smartreader.o: smartreader.c
	$(CC_RULE)

$O/smartreader_help.o: smartreader_help.help
	$(HELP_RULE)

$O/wildcard.o: wildcard.c
	$(CC_RULE)

$O/tables.o: tables.c
	$(CC_RULE)

$O/smart_utils.o: smart_utils.c
	$(CC_RULE)

$O/vec3.o: vec3.c
	$(CC_RULE)

$O/png_reader.o: png_reader.c
	$(CC_RULE)

$O/fits_reader.o: fits_reader.c
	$(CC_RULE)

$O/tiff_reader.o: tiff_reader.c
	$(CC_RULE)

$O/fiff_reader.o: fiff_reader.c
	$(CC_RULE)

$O/son_reader.o: son_reader.c
	$(CC_RULE)

$O/ushort_reader.o: ushort_reader.c
	$(CC_RULE)

$O/multi_reader.o: multi_reader.c
	$(CC_RULE)

$O/afni_reader.o: afni_reader.c
	$(CC_RULE)

$O/lx_reader_excite.o: lx_reader.c \
		../fmri/frozen_header_info.h \
		../fmri/frozen_header_info_excite.h
	@echo "%%%% Compiling lx_reader.c to $(@F) %%%%"
	$(CC_FOR_EXCITE) lx_reader.c

$O/lx_reader_cnv4.o: lx_reader.c \
		../fmri/frozen_header_info.h \
		../fmri/frozen_header_info_cnv4.h
	@echo "%%%% Compiling lx_reader.c to $(@F) %%%%"
	$(CC_FOR_CNV4) lx_reader.c

$O/lx_reader_lx2.o: lx_reader.c  \
		../fmri/frozen_header_info.h \
		../fmri/frozen_header_info_lx2.h
	@echo "%%%% Compiling lx_reader.c to $(@F) %%%%"
	$(CC_FOR_LX2) lx_reader.c

$O/lx_reader_prelx.o: lx_reader.c  \
		../fmri/frozen_header_info.h \
		../fmri/frozen_header_info_prelx.h
	@echo "%%%% Compiling lx_reader.c to $(@F) %%%%"
	$(CC_FOR_PRELX) lx_reader.c

$O/lx_2dfast_reader_excite.o: lx_2dfast_reader.c  \
		../fmri/frozen_header_info.h \
		../fmri/frozen_header_info_excite.h
	@echo "%%%% Compiling lx_2dfast_reader.c to $(@F) %%%%"
	$(CC_FOR_EXCITE) lx_2dfast_reader.c

$O/lx_2dfast_reader_cnv4.o: lx_2dfast_reader.c  \
		../fmri/frozen_header_info.h \
		../fmri/frozen_header_info_cnv4.h
	@echo "%%%% Compiling lx_2dfast_reader.c to $(@F) %%%%"
	$(CC_FOR_CNV4) lx_2dfast_reader.c

$O/lx_2dfast_reader_prelx.o: lx_2dfast_reader.c  \
		../fmri/frozen_header_info.h \
		../fmri/frozen_header_info_prelx.h
	@echo "%%%% Compiling lx_2dfast_reader.c to $(@F) %%%%"
	$(CC_FOR_PRELX) lx_2dfast_reader.c

$O/lx_2dfast_reader_lx2.o: lx_2dfast_reader.c  \
		../fmri/frozen_header_info.h \
		../fmri/frozen_header_info_lx2.h
	@echo "%%%% Compiling lx_2dfast_reader.c to $(@F) %%%%"
	$(CC_FOR_LX2) lx_2dfast_reader.c

$O/lx_epibold_reader_excite.o: lx_epibold_reader.c  \
		../fmri/frozen_header_info.h \
		../fmri/frozen_header_info_excite.h
	@echo "%%%% Compiling lx_epibold_reader.c to $(@F) %%%%"
	$(CC_FOR_EXCITE) lx_epibold_reader.c

$O/lx_epibold_reader_cnv4.o: lx_epibold_reader.c  \
		../fmri/frozen_header_info.h \
		../fmri/frozen_header_info_cnv4.h
	@echo "%%%% Compiling lx_epibold_reader.c to $(@F) %%%%"
	$(CC_FOR_CNV4) lx_epibold_reader.c

$O/lx_epibold_reader_lx2.o: lx_epibold_reader.c \
		../fmri/frozen_header_info.h \
		../fmri/frozen_header_info_lx2.h
	@echo "%%%% Compiling lx_epibold_reader.c to $(@F) %%%%"
	$(CC_FOR_LX2) lx_epibold_reader.c

$O/lx_epibold_reader_prelx.o: lx_epibold_reader.c \
		../fmri/frozen_header_info.h \
		../fmri/frozen_header_info_prelx.h
	@echo "%%%% Compiling lx_epibold_reader.c to $(@F) %%%%"
	$(CC_FOR_PRELX) lx_epibold_reader.c

$O/lx_splx_reader_excite.o: lx_splx_reader.c \
		../fmri/frozen_header_info.h \
		../fmri/frozen_header_info_excite.h
	@echo "%%%% Compiling lx_splx_reader.c to $(@F) %%%%"
	$(CC_FOR_EXCITE) lx_splx_reader.c

$O/lx_splx_reader_cnv4.o: lx_splx_reader.c \
		../fmri/frozen_header_info.h \
		../fmri/frozen_header_info_cnv4.h
	@echo "%%%% Compiling lx_splx_reader.c to $(@F) %%%%"
	$(CC_FOR_CNV4) lx_splx_reader.c

$O/lx_splx_reader_lx2.o: lx_splx_reader.c \
		../fmri/frozen_header_info.h \
		../fmri/frozen_header_info_lx2.h
	@echo "%%%% Compiling lx_splx_reader.c to $(@F) %%%%"
	$(CC_FOR_LX2) lx_splx_reader.c

$O/lx_splx_reader_prelx.o: lx_splx_reader.c \
		../fmri/frozen_header_info.h \
		../fmri/frozen_header_info_prelx.h
	@echo "%%%% Compiling lx_splx_reader.c to $(@F) %%%%"
	$(CC_FOR_PRELX) lx_splx_reader.c

$O/lx_sf11_reader_excite.o: lx_sf11_reader.c \
		../fmri/frozen_header_info.h \
		../fmri/frozen_header_info_excite.h
	@echo "%%%% Compiling lx_sf11_reader.c to $(@F) %%%%"
	$(CC_FOR_EXCITE) lx_sf11_reader.c

$O/lx_sf11_reader_cnv4.o: lx_sf11_reader.c \
		../fmri/frozen_header_info.h \
		../fmri/frozen_header_info_cnv4.h
	@echo "%%%% Compiling lx_sf11_reader.c to $(@F) %%%%"
	$(CC_FOR_CNV4) lx_sf11_reader.c

$O/lx_sf11_reader_lx2.o: lx_sf11_reader.c \
		../fmri/frozen_header_info.h \
		../fmri/frozen_header_info_lx2.h
	@echo "%%%% Compiling lx_sf11_reader.c to $(@F) %%%%"
	$(CC_FOR_LX2) lx_sf11_reader.c

$O/lx_sf11_reader_prelx.o: lx_sf11_reader.c \
		../fmri/frozen_header_info.h \
		../fmri/frozen_header_info_prelx.h
	@echo "%%%% Compiling lx_sf11_reader.c to $(@F) %%%%"
	$(CC_FOR_PRELX) lx_sf11_reader.c

$O/windaq_reader.o: windaq_reader.c
	$(CC_RULE)

$O/desmith_reader.o: desmith_reader.c
	$(CC_RULE)

$O/raw_reader.o: raw_reader.c
	$(CC_RULE)

$O/base_reader.o: base_reader.c
	$(CC_RULE)

$O/ram_reader.o: ram_reader.c
	$(CC_RULE)

$O/convert_reader.o: convert_reader.c
	$(CC_RULE)

$O/lx_image_reader.o: lx_image_reader.c
	$(CC_RULE)

$O/analyze_reader.o: analyze_reader.c
	$(CC_RULE)

$O/nifti_reader.o: nifti_reader.c
	$(CC_RULE)

$O/dicom_reader.o: dicom_reader.c
	$(CC_RULE)

$O/dicom_transfer_syntax.o: dicom_transfer_syntax.c
	$(CC_RULE)

$O/dicom_uid_dict.o: dicom_uid_dict.c
	$(CC_RULE)

$O/dicom_parser.o: dicom_parser.c
	$(CC_RULE)

$O/siemens_kspace_reader.o: siemens_kspace_reader.c
	$(CC_RULE)

$O/pghmri_reader.o: pghmri_reader.c
	$(CC_RULE)

$(CB)/smartreader: $(SMARTREADER_OBJS) $(LIBFILES)
	@echo "%%%% Linking $(@F) %%%%"
	@$(LD) $(LFLAGS) -o $(CB)/$(@F) $(SMARTREADER_OBJS) $(LIBS)

releaseprep:
	echo "no release prep from " `pwd`






