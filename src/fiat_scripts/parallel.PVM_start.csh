#!/bin/csh -f
# parallel.PVM_start.csh
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
echo '# $Id: parallel.PVM_start.csh,v 1.1 2005/03/09 01:47:18 welling Exp $'
echo '# '`date` $0

echo "performing PVM-specific checks..."     

#---------------------------------------------------------------------------*
# Phase 1: setup                                                            *
#---------------------------------------------------------------------------*

  #construct list of hosts from F_PARALLEL_HOSTS
  setenv parallel_hosts	""
  foreach m ($F_PARALLEL_HOSTS)
	set hname = `echo $m | awk -F: '{print $1}'`
	set parallel_hosts = "$parallel_hosts $hname"
  end

  # How are we doing communication?
  if (${?PVM_RSH}) then
    set pvm_rsh = ${PVM_RSH}
  else
    set pvm_rsh = rsh
  endif

#---------------------------------------------------------------------------*
# Phase 2: check pvm environment variables                                  *
#---------------------------------------------------------------------------*
  if ( ${?PVM_ROOT} ) then
        setenv PVM ${PVM_ROOT}/lib/pvm
        if (! -f $PVM) then
                echo "can't find pvm executable: $PVM"
                exit 1
        endif
  else
	echo "You must set PVM_ROOT in your .cshrc!\n"
	exit 1
  endif

  if (! ${?PVM_ARCH} ) then
	if (! -e ${PVM_ROOT}/lib/pvmgetarch) then
		echo "File ${PVM_ROOT}/lib/pvmgetarch is missing."
		exit 1
	endif
	if (! -x ${PVM_ROOT}/lib/pvmgetarch) then
		echo "Cannot execute ${PVM_ROOT}/lib/pvmgetarch"
		exit 1
	endif
	setenv PVM_ARCH `${PVM_ROOT}/lib/pvmgetarch`
  endif
  echo "PVM_ARCH is: ${PVM_ARCH}"

  if (! $?PVM_EXPORT) then
        setenv PVM_EXPORT DISPLAY
  endif
  if ( $?PVM_DEBUGGER) then
	setenv PVM_EXPORT "PVM_DEBUGGER:$PVM_EXPORT"
  endif
  echo "PVM_EXPORT is: $PVM_EXPORT"

#---------------------------------------------------------------------------*
# Phase 3: check PVM access to executables                                  *
#---------------------------------------------------------------------------*
  echo "checking executables are accessible to PVM..."
  #make PVM_ROOT directories as needed
  if ( ! -d ~/pvm3) then
	echo "Making pvm3 subdirectory in your home directory"
        mkdir -p ~/pvm3
  endif
  setenv MY_PVM_ROOT ~/pvm3

# check and set afs permissions
  echo "checking filesystem type"

  set tmp = "$F_TEMP/par.$$.tmp"
  set lead=`(cd $HOME; pwd) | sed s+/+' '+g | awk '{print $1}'`
  if ($lead == "afs") then
        set status = `test_in_subshell.csh fs la $HOME >& $tmp`
        if ($status) then
                echo "unable to read acls on $HOME"
        	rm -f $tmp
                exit 1
        endif
	set tmp2 = "$F_TEMP/par.$$.tmp2"
	set status = `test_in_subshell.csh grep system:anyuser $tmp >& $tmp2`
        if ($status) then
		unset echo
		echo
		echo "You must make your home directory world readable in order"
		echo "to use parallel FIASCO on an AFS filesystem."
		echo "To do this use the command: fs sa $HOME system:anyuser rl"
        	rm -f $tmp
        	rm -f $tmp2
		exit 1
        endif
	set status = `test_in_subshell.csh grep rl $tmp2`
        if ($status) then
		unset echo
		echo
		echo "You must make your home directory world readable in order"
		echo "to use parallel FIASCO on an AFS filesystem."
		echo "To do this use the command: fs sa $HOME system:anyuser rl"
        	rm -f $tmp
        	rm -f $tmp2
		exit 1
        endif
        rm -f $tmp2

	echo "Making your ~/pvm3/bin/${PVM_ARCH} directory world-readable"
	foreach d (~/pvm3 ~/pvm3/bin ~/pvm3/bin/${PVM_ARCH})
	        fs sa $d "system:anyuser" rl
	end
	rm -f $tmp

  else
        echo "Don't yet know how to check non-afs access permissions..."
  endif


#---------------------------------------------------------------------------*
# Phase 4: check .rhosts or .ssh set appropriately
#---------------------------------------------------------------------------*
if ( $pvm_rsh =~ *rsh ) then
  echo "checking .rhosts file set appropriately..."
  if ( ! -f ~/.rhosts ) then
        unset echo
	echo
	echo "You should create an .rhosts file in your home directory"
	echo "containing the list of machines you plan to use in your"
	echo "parallel configuration."
	exit 1
  endif
  set rherr=0
  foreach m ($parallel_hosts)
	echo "Checking for $m in .rhosts file"
        set status = `test_in_subshell.csh grep $m ~/.rhosts`
        if ( $status) then
		unset echo
		echo
		echo "You need to add the line:"
		echo "$m $USER"
		echo "to your ~/.rhosts	file."
		set rherr=1
        endif
  end
  if ( $rherr) then
    exit 1
  endif
else if ( $pvm_rsh =~ *ssh ) then
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
endif

#---------------------------------------------------------------------------*
# Phase 5: kill off any pvm daemons that are already running  (strictly, this
#		should not be necessary, but it is because old pvm daemons
#		occasionally get into a wedged state and really foul things up)
#---------------------------------------------------------------------------*
  foreach m ($parallel_hosts)
    set d = `${pvm_rsh} $m "ps -u $USER" </dev/null | awk '$4~/pvmd/ {print $1}'`
    foreach p ($d)
	${pvm_rsh} $m "kill $p" </dev/null
    end
  end

#---------------------------------------------------------------------------*
# Phase 6: start pvm daemon
#---------------------------------------------------------------------------*

# start the local pvm daemon, obtain configuration
  echo conf | $PVM >& $tmp
  echo version | $PVM
  grep ' running' $tmp >& /dev/null
  if ($status) then
        echo "$PVM failed to start pvm daemon on" `hostname`
        rm -f $tmp
        exit 1
  endif
  rm -f $tmp



