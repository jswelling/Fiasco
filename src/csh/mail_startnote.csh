#!/bin/csh -efx
# mail_startnote.csh
#/************************************************************
# *                                                          *
# *  Permission is hereby granted to any individual or       *
# *  institution for use, copying, or redistribution of      *
# *  this code and associated documentation, provided        *
# *  that such code and documentation are not sold for       *
# *  profit and the following copyright notice is retained   *
# *  in the code and documentation:                          *
# *     Copyright (c) 1995 Department of Statistics,         *
# *                        Carnegie Mellon University        *
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
# *  Original programming by Bill Eddy                       *
# ************************************************************/
#
echo '#'`date` $0
echo '#$Id: mail_startnote.csh,v 1.3 2000/08/21 22:32:54 welling Exp $'

# Allow the user to avoid sending a summary
if ( ${?F_NO_MAIL_SUMMARY} ) exit 0

set addr = "fiascosummaries@stat.cmu.edu"
echo '##Mailing start note information to ' $addr
Mail -s 'Fiasco start note' $addr <<EOF 
args: "$*"
date: "`date`"
user: $LOGNAME
dir: $F_DIR
description: "$F_DESCRIPTION"
header: "$F_HEADER"
EOF


