#!/bin/csh -f
# parallel.finish.csh
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
# $Id: parallel.finish.csh,v 1.1 2005/03/09 01:47:18 welling Exp $

set oldecho = ${?echo}
unset echo

echo '#'`date` $0

#---------------------------------------------------------------------------*
# Do method-specific cleanup
#---------------------------------------------------------------------------*

set mode = `parallel_mode`

echo "Doing parallel shutdown for mode ${mode}, if any"

if ( -x ./parallel.${mode}_finish.csh ) then
  source ./parallel.${mode}_finish.csh
else if ( -x ${FIASCO}/parallel.${mode}_finish.csh ) then
  source ${FIASCO}/parallel.${mode}_finish.csh
endif

#---------------------------------------------------------------------------*
# general clean up
#---------------------------------------------------------------------------*
unsetenv PAR_ENABLE

echo "Parallel shutdown complete."
if ( $oldecho ) set echo

