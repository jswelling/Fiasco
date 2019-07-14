#!/bin/csh -f
# parallel.PVM_finish.csh
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
# $Id: parallel.PVM_finish.csh,v 1.1 2005/03/09 01:47:18 welling Exp $
echo '#'`date` $0

#---------------------------------------------------------------------------*
# Do PVM-specific cleanup
#---------------------------------------------------------------------------*

  echo "halting PVM daemon and cleaning up ..."
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
  set retcode = `test_in_subshell.csh "echo halt | $PVM >& /dev/null"`
  echo "PVM daemon shut down; return code $retcode"

# clean up old PVM log files
  echo "cleaning log files from the hosts"
  if ( ${?PVM_RSH} ) then
    set my_rsh = $PVM_RSH
  else
    set my_rsh = rsh
  endif
  foreach m ($F_PARALLEL_HOSTS)
    set hname = `echo $m | awk -F: '{print $1}'`
    echo "  cleaning logs from $hname"
    set d = `($my_rsh $hname "ls -l /tmp" </dev/null |grep pvm |grep $USER |awk '{print $9}'; exit 0)`
    foreach f ($d)
	($my_rsh $hname "rm -f /tmp/$f" </dev/null; exit 0)
    end
  end
