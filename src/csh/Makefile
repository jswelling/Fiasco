# csh/Makefile
#/************************************************************
# *                                                          *
# *  Permission is hereby granted to any individual or       *
# *  institution for use, copying, or redistribution of      *
# *  this code and associated documentation, provided        *
# *  that such code and documentation are not sold for       *
# *  profit and the following copyright notice is retained   *
# *  in the code and documentation:                          *
# *     Copyright (c) 1996 Department of Statistics,         *
# *                        Carnegie Mellon University        *
# *                                                          *
# *  Original programming by Bill Eddy                       *
# ************************************************************/
#

PKG	= csh
PKG_MAKETARGETS = copy_csh make_docs
ALL_MAKEFILES= Makefile
SCRIPTFILES = FIASCO *.csh *.py
SCRIPTSWITHDOCS = FIASCO groupcompare.csh make_merged_maps.csh \
	merge_detrend.csh false_discovery.py testingroup.csh \
	pooled_mean.csh pooled_stdv.csh stats_to_floats.csh \
	make_cine.csh make_phase_cine.csh color_by_phase.csh \
	overlay.csh physio_correct_triggered.py \
        dataset_matches_env.py pick_n_voxels.py customize_local.py \
	coregister_inplane.py make_P_P_plots.csh \
	coregister_struct_to_inplane.py coregister.py \
	multi_runner.py renumber_ifiles.py convert_ifiles.py \
        make_group_anat.csh afni_estireg3d.csh combine.csh \
	calc_stats_over_roi.py icbm_to_tlrc.py tlrc_to_icbm.py \
	pca.py align_by_permutation.py coregister_makeps.py \
	cowarp_inplane.py paste_split.py build_model_matrix.py \
	pullback_roi_mask.py plot_roi_time_series.py strip_skull.py \
	epoch_subset.py afni_tshift.py plot_roi_event_response.py

include ../Makefile_pkg

copy_csh:
	-rm -f fiasco.local.csh
	awk '{ if (index($$0,"setenv SPLUS")==1) { \
			print "#"$$0; print "setenv SPLUS \"$(SPLUSEXE)\""} \
		else { print $$0}} ' < fiasco.local_orig.csh > fiasco.local.csh
	chmod go+x fiasco.local.csh $(SCRIPTFILES)
	-chmod u+x fiasco.local.csh $(SCRIPTFILES)
	-cp -f $(SCRIPTFILES) $B

make_docs:
	@( for i in $(SCRIPTSWITHDOCS) ; do \
	  echo "%%%% Generating $$i.html %%%%" ; \
          env PATH="$B:$(PATH)" $B/$$i -help htmldoc > $D//$$i.html 2> /dev/null ; \
	  : ; \
	done )

releaseprep:
	if test ! \( -f fiasco.local_orig.csh \); \
		then cp fiasco.local.csh fiasco.local_orig.csh; \
		else echo 'fiasco.local_orig.csh already exists' ; fi
	chmod +xr $(SCRIPTFILES)

