#!/bin/csh -f
# parallel.start.csh
#/************************************************************
# *                                                          *
# *  Permission is hereby granted to any individual or       *
# *  institution for use, copying, or redistribution of      *
# *  this code and associated documentation, provided        *
# *  that such code and documentation are not sold for       *
# *  profit and the following copyright notice is retained   *
# *  in the code and documentation:                          *
# *     Copyright (c) 1997 Pittsburgh Supercomputing Center  *
# *                                                          *
# *  This program is distributed in the hope that it will    *
# *  be useful, but WITHOUT ANY WARRANTY; without even the   *
# *  implied warranty of MERCHANTABILITY or FITNESS FOR A    *
# *  PARTICULAR PURPOSE.  Neither Carnegie Mellon University *
# *  nor any of the authors assume any liability for         *
# *  damages, incidental or otherwise, caused by the         *
# *  installation or use of this software.                   *
# *                                                          *
# *  CLINICAL APPLICATIONS ARE NOT RECOMMENDED, AND THIS     *
# *  SOFTWARE HAS NOT BEEN EVALUATED BY THE UNITED STATES    *
# *  FDA FOR ANY CLINICAL USE.                               *
# *                                                          *
# *                                                          *
# *  Original programming by Nigel Goddard & Greg Hood (PSC) *
# ************************************************************/
#
echo '# $Id: parallel.start.csh,v 1.1 2005/03/09 01:47:18 welling Exp $'
echo '# '`date` $0

set oldecho = ${?echo}
unset echo

echo "performing parallel execution checks..."     
#---------------------------------------------------------------------------*
# Phase 1: check .login and .cshrc don't produce output
#---------------------------------------------------------------------------*
  set tmp = "$F_TEMP/par.$$.tmp"
  echo "checking .cshrc and .login..."
  if ( -e $HOME/.cshrc ) then
    csh -f $HOME/.cshrc > $tmp
    if ( ! -z $tmp) then
	unset echo
        echo 'Your ~/.cshrc produces output for non-interactive shells\!'
        echo 'Edit your .cshrc to remedy this, e.g. add'
        echo '     if (! $?prompt ) then'
        echo '         ...'
        echo '     endif'
        echo 'around the output statements'
        echo 'To check, verify that "(csh -f ~/.cshrc)" prints nothing'
        exit 1
    endif
    rm -f $tmp
  endif

  if ( -e $HOME/.login ) then
    csh -f $HOME/.login > $tmp
    if ( ! -z $tmp) then
	unset echo
        echo 'Your ~/.login produces output for non-interactive shells\!'
        echo 'Edit your .login to remedy this, e.g. add'
        echo '     if (! $?prompt ) then'
        echo '         ...'
        echo '     endif'
        echo 'around the output statements'
        echo 'To check, verify that "(csh -f ~/.login)" prints nothing'
        rm -f $tmp
        exit 1
    endif
    rm -f $tmp
  endif

#---------------------------------------------------------------------------*
# Phase 2: do specifics for the current parallel mode
#---------------------------------------------------------------------------*

set mode = `parallel_mode`

if ( -x ./parallel.${mode}_start.csh ) then
  source ./parallel.${mode}_start.csh
  if ( $status ) exit $status
else if ( -x ${FIASCO}/parallel.${mode}_start.csh ) then
  source ${FIASCO}/parallel.${mode}_start.csh
  if ( $status ) exit $status
endif

#---------------------------------------------------------------------------*
# Phase 3: enable parallelism
#---------------------------------------------------------------------------*

# enable parallelism
  setenv PAR_ENABLE 1
  setenv PAR_HOSTS "$F_PARALLEL_HOSTS"


#---------------------------------------------------------------------------*
# All done; clean up
#---------------------------------------------------------------------------*

if ( $oldecho ) set echo

