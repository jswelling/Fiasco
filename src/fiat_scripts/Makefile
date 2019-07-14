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

PKG	= fiat_scripts
PKG_MAKETARGETS = cover_matching_scripts copy_csh make_docs
ALL_MAKEFILES= Makefile
SCRIPTFILES = *.csh 
SCRIPTSWITHDOCS = 

include ../Makefile_pkg

cover_matching_scripts:
	@( for i in \( $(SCRIPTFILES) \); do \
	  if test \( -f ../csh/$$i \); \
		then ( echo "%%%% Hiding ../csh/$$i %%%%" ; \
			mv ../csh/$$i ../csh/$$i.old ); fi; \
	done )

copy_csh:
	@echo "%%%% Exporting scripts from $(PKG) %%%%"
	@-cp -f $(SCRIPTFILES) $B

make_docs:
	@( for i in $(SCRIPTSWITHDOCS) dummy ; do \
	  if [ $$i != dummy ] ; \
	     then \
	       echo "%%%% Generating $$i.html %%%%" ; \
               PATH="$B:$(PATH)" $B/$$i -help htmldoc \
	         > $D//$$i.html 2> /dev/null ; fi;\
	done )

releaseprep:
	chmod +xr $(SCRIPTFILES)

