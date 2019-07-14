#
#       Makefile for spiral reconstruction program
#
#       Copyright (c) 1998  Pittsburgh Supercomputing Center
#
#       HISTORY
#               1/96    Written by Greg Hood (PSC)
#               8/96    Made architecture-independent (Greg Hood, PSC)
#		8/98	Combined sgrid & srecon into spiral
#

PKG          = spiral
PKG_LIBS     = -lfmri -lcrg -lmri -lpar -lbio -lacct -lmisc -larray \
	     -lrttraj $(LAPACK_LIBS) -lm
PKG_MAKEBINS = $(CB)/checkpfile $(CB)/pfileorder $(CB)/spiral \
	$(CB)/srecon $(CB)/sgrid $(CB)/spiral_reader $(CB)/slow_ft

ALL_MAKEFILES= Makefile
CSOURCE= checkpfile.c commun.c header.c \
	pfileorder.c spiral.c worker.c worker_utils.c sgrid.c srecon.c \
	spiral_reader.c mriheader.c slow_ft.c dirichlet.c vpolygon.c
HFILES= rdb.h checkpfile.h spiral.h dirichlet.h vpolygon.h
DOCFILES= pfileorder_help.help spiral_help.help sgrid_help.help \
	srecon_help.help spiral_reader_help.help slow_ft_help.help

include ../Makefile_pkg

PKG_OBJS     = $O/commun.o $O/header.o $O/spiral.o \
	     $O/worker.o $O/worker_utils.o $O/mriheader.o

LIBFILES= $L/libfmri.a $L/libpar.a $L/libbio.a $L/libarray.a \
	$L/libmisc.a $L/libacct.a $L/librttraj.a
HDRS= mri.h par.h bio.h array.h misc.h acct.h

$O/worker.o: worker.c
	$(CC_RULE)

$O/worker_utils.o: worker_utils.c
	$(CC_RULE)

$O/header.o: header.c
	$(CC_RULE)

$O/mriheader.o: mriheader.c
	$(CC_RULE)

$O/commun.o: commun.c
	$(CC_RULE)

$O/spiral.o: spiral.c
	$(CC_RULE)

$O/spiral_reader.o: spiral_reader.c
	$(CC_RULE)

$O/spiral_reader_help.o: spiral_reader_help.help
	$(HELP_RULE)

$O/spiral_help.o: spiral_help.help
	$(HELP_RULE)

$(CB)/spiral: $(PKG_OBJS) $O/spiral_help.o $(LIBFILES)
	$(MULTI_HELP_LD)

$(CB)/spiral_reader: $O/spiral_reader.o $O/header.o $O/worker_utils.o \
		     $O/mriheader.o $O/spiral_reader_help.o
	@echo %%%% Linking spiral_reader %%%%
	@$(LD) $(LFLAGS) -o $(CB)/spiral_reader \
	      $O/spiral_reader.o $O/header.o $O/worker_utils.o \
	      $O/mriheader.o $O/spiral_reader_help.o $(LIBS)

$O/dirichlet.o: dirichlet.c
	$(CC_RULE)

$O/vpolygon.o: vpolygon.c
	$(CC_RULE)

$O/slow_ft.o: slow_ft.c
	$(CC_RULE)

$O/slow_ft_help.o: slow_ft_help.help
	$(HELP_RULE)

$(CB)/slow_ft: $O/slow_ft.o $O/dirichlet.o $O/vpolygon.o $O/slow_ft_help.o
	@echo %%%% Linking slow_ft %%%%
	@$(LD) $(LFLAGS) -o $(CB)/slow_ft \
		$O/slow_ft.o $O/dirichlet.o $O/slow_ft_help.o \
		$O/vpolygon.o $(LIBS)

$O/pfileorder.o: pfileorder.c
	$(CC_RULE)

$O/pfileorder_help.o: pfileorder_help.help
	$(HELP_RULE)

$(CB)/pfileorder: $O/pfileorder.o $O/pfileorder_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

#$O/tinypfile.o: tinypfile.c
#	$(CC_RULE)

#$(CB)/tinypfile: $O/tinypfile.o $O/header.o $(LIBFILES)
#	$(LD) $(LFLAGS) -o $(CB)/tinypfile $O/tinypfile.o $O/header.o $(LIBS)

$O/checkpfile.o: checkpfile.c
	$(CC_RULE)

$(CB)/checkpfile: $O/checkpfile.o $O/header.o $(LIBFILES)
	@echo %%%% Linking checkpfile %%%%
	@$(LD) $(LFLAGS) -o $(CB)/checkpfile $O/checkpfile.o $O/header.o \
		$(LIBS)

$O/sgrid.o: sgrid.c
	$(CC_RULE)

$O/srecon.o: srecon.c
	$(CC_RULE)

$O/sgrid_help.o: sgrid_help.help
	$(HELP_RULE)

$O/srecon_help.o: srecon_help.help
	$(HELP_RULE)

$(CB)/sgrid: $O/sgrid.o $O/sgrid_help.o $(LIBFILES)
	@echo %%%% Linking sgrid %%%%
	@$(LD) $(LFLAGS) -o $(CB)/sgrid $O/sgrid.o \
		$O/sgrid_help.o $(LIBS)

$(CB)/srecon: $O/srecon.o $O/srecon_help.o $(LIBFILES)
	$(SINGLE_HELP_LD)

releaseprep:
	echo "no release prep from " `pwd`

