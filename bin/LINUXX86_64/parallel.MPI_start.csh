#!/bin/csh -f
# parallel.MPI_start.csh
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
echo '# $Id: parallel.MPI_start.csh,v 1.2 2005/10/12 18:05:38 welling Exp $'
echo '# '`date` $0

echo "performing MPI-specific checks..."     

#
# Not much to do; mpirun takes care of most logistics for us.
#

#---------------------------------------------------------------------------*
# Phase 1: check .rhosts or .ssh set appropriately
#---------------------------------------------------------------------------*

echo "checking .ssh files set appropriately..."
if ( ! -f ~/.ssh/known_hosts && ! -f ~/.ssh/known_hosts2) then
      unset echo
      echo
      echo "It looks like your account is not set up to use ssh."
      echo "The easiest way to do this is just to use the ssh"
      echo "command to connect once to each of the hosts in your"
      echo "parallel computing cluster."
      exit 1
endif
set rherr=0
foreach m ($parallel_hosts)
      echo "Checking that $m is a known host..."
      set status = `test_in_subshell.csh grep $m ~/.ssh/known_hosts`
      set status2 = `test_in_subshell.csh grep $m ~/.ssh/known_hosts2`
      if ( $status && $status2 ) then
	unset echo
	echo
	echo "You need to use ssh to connect to $m, just once"
	echo "so its identity will become known."
	exit 1
       endif
end
echo "checking that $HOST is allowed to connect to the remote hosts..."
set ssh_host_auth = 0
if ( -e ~/.ssh/authorized_keys ) then
  if ( ! `test_in_subshell.csh grep $HOST ~/.ssh/authorized_keys` ) then
    set ssh_host_auth = 1
  endif
else if ( -e ~/.ssh/authorized_keys2 ) then
  if ( ! `test_in_subshell.csh grep $user@$HOST ~/.ssh/authorized_keys2` ) \
    then
    set ssh_host_auth = 1
  endif
endif
if ( ! $ssh_host_auth ) then
	unset echo
	echo
	echo "You need to use ssh-keygen to prove your identity"
        echo "to the remote hosts.  To do this, do the commands:"
	echo "   ssh-keygen -t dsa"
	echo "   cd ${HOME}/.ssh"
	echo "   cat id_dsa.pub >> authorized_keys2"
	exit 1
endif
