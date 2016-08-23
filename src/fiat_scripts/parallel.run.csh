#!/bin/csh -f
# parallel.run.csh
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
echo '# $Id: parallel.run.csh,v 1.1 2005/03/09 01:47:18 welling Exp $'
echo '# '`date` $0

# Bail if parallelism isn't enabled; it saves time and trouble
if ( ! ${?PAR_ENABLE} ) then
  $*
  exit
endif

set mode = `parallel_mode`
switch ( $mode )
# For the non-parallel case, just run it.
  case "NONE":
    $*
    breaksw

# For PVM, parallelism happens auto-magically.
  case "PVM":
    set progname = $1
    set full_progname = `which $progname`
    set tailname = ${progname:t}
    shift

    if ( ${?MY_PVM_ROOT} ) then
      set my_pvm_root = ${MY_PVM_ROOT}
    else 
      set my_pvm_root = ${HOME}/pvm3
    endif
    if ( ! -d ${my_pvm_root}/bin ) then
      echo "Making ${my_pvm_root}/bin"
      mkdir ${my_pvm_root}/bin
    endif
    set linkdir = ${my_pvm_root}/bin/${PVM_ARCH}
    if ( ! -d $linkdir ) then
      echo "Making $linkdir"
      mkdir $linkdir
    endif
    set linkname = $linkdir/${tailname}.link.$$
    if ( -x ${my_pvm_root}/bin/${progname} ) then
      echo "WARNING: there is an old link to ${progname} in ${linkdir}!"
      echo "         (this is harmless but suspicious)"
    endif
    # link executable to where PVM can find it
    echo ln -s ${full_progname} ${linkname}
    ln -s ${full_progname} ${linkname}
    # run the program via the new link
    ${linkname} $*
    # clean up the link
    echo rm -f ${linkname}
    rm -f ${linkname}
    breaksw

# Under MPI, an appropriate mpirun command must be built.
  case "MPI":
    set progname = $1
    set full_progname = `which $progname`
    shift
    set par_args = ""
    if (${?PAR_ENABLE}) then
      set par_args = "$par_args -PAR_ENABLE=$PAR_ENABLE"
    endif
    if (${?PAR_DEBUG}) then
      set par_args = "$par_args -PAR_DEBUG=$PAR_DEBUG"
    endif
    if (${?PAR_GROUP}) then
      set par_args = "$par_args -PAR_GROUP=$PAR_GROUP"
    endif
    if (${?PAR_NOSPAWN}) then
      set par_args = "$par_args -PAR_NOSPAWN=$PAR_NOSPAWN"
    endif
    if (${?PAR_CWD}) then
      set par_args = "$par_args -PAR_CWD=$PAR_CWD"
    endif
    if (${?PAR_VERBOSE}) then
      set par_args = "$par_args -PAR_VERBOSE=$PAR_VERBOSE"
    endif
    set tmp = /tmp/par_run_tmp_$$
    touch $tmp
    #
    # Parse PAR_HOSTS according to the following rules:
    #
    #  Entries of the form "hostname" request one additional process 
    #  and add the given host to the pool.  
    #
    #  Entries of the form "hostname:n" where n is some integer add 
    #  hostname to the pool and request n additional procs; mpich will 
    #  supposedly set things up so that exactly n procs end up on that 
    #  host.
    #
    #  A single entry of the form MPI:n will request n procs but not
    #  specify any specific hosts, so the MPI implementation can decide
    #  where to put the processes.
    #
    @ nprocs = 0
    @ nhosts = 0
    @ maxcluster = 1
    if ( ${?PAR_HOSTS} ) then
      foreach worker ( $PAR_HOSTS )
        echo $worker >> $tmp
	if ( $worker =~ *:* ) then
	  set theseprocs = `echo $worker | awk -F: '{print $2}'`
	  set thismach = `echo $worker | awk -F: '{print $1}'`
	  @ nprocs = $nprocs + $theseprocs
	  if ( $thismach != MPI ) then
            @ nhosts = $nhosts + 1
	  endif
	  if ( $theseprocs > $maxcluster ) set maxcluster = $theseprocs
	else
	  @ nprocs = $nprocs + 1
	  @ nhosts = $nhosts + 1
	endif
      end
      unsetenv PAR_HOSTS
      setenv MPI_MAX_CLUSTER_SIZE $maxcluster
    endif
    if ( $nhosts > 0 ) then
      echo mpirun -np $nprocs -machinefile $tmp $full_progname $* $par_args
      mpirun -np $nprocs -machinefile $tmp $full_progname $* $par_args
    else
      echo mpirun -np $nprocs $full_progname $* $par_args
      mpirun -np $nprocs $full_progname $* $par_args
    endif
    rm $tmp
    breaksw

  default:
    echo "Unknown parallel mode ${mode}!"
    exit -1
    breaksw
endsw
