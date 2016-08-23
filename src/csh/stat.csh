#!/bin/csh -efx
# stat.csh
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
echo '#$Id: stat.csh,v 1.5 2001/11/22 00:04:42 welling Exp $'
if ( ! -d $F_STAT_MPATH) mkdir $F_STAT_MPATH
if ( ! -d $F_STAT_SPATH) mkdir $F_STAT_SPATH
if ( ! -d $F_STAT_TPATH) mkdir $F_STAT_TPATH
stats -condition $F_SPLIT_NEW -input $1.mri \
      -meanprefix $F_STAT_MPATH -stdvprefix $F_STAT_SPATH \
      -tsprefix $F_STAT_TPATH -maxtpairs $F_STAT_MAXTPR
