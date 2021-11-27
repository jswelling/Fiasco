#!/usr/bin/env python
# convert_ifiles.py
#/************************************************************
# *                                                          *
# *  Permission is hereby granted to any individual or       *
# *  institution for use, copying, or redistribution of      *
# *  this code and associated documentation, provided        *
# *  that such code and documentation are not sold for       *
# *  profit and the following copyright notice is retained   *
# *  in the code and documentation:                          *
# *     Copyright (c) 1995,1996 Department of Statistics,    *
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
# ************************************************************/
#
import sys
import os
import os.path
import string
import getopt
import re
if "FIASCO" in os.environ:
    sys.path.append(os.environ["FIASCO"])
from fiasco_utils import *

idString= '# $Id: convert_ifiles.py,v 1.2 2004/04/08 22:21:36 welling Exp $'

ifileRegex= re.compile('I\.\d\d\d')
gzIfileRegex= re.compile('I\.\d\d\d\.gz')
zIfileRegex= re.compile('I\.\d\d\d\.Z')

##############################
#
# Main
#
##############################

# Check for "-help"
if len(sys.argv)>1:
    if sys.argv[1] == "-help":
        if len(sys.argv)>2:
            os.system( "scripthelp %s %s"%(sys.argv[0],sys.argv[2]) )
        else:
            os.system( "scripthelp %s"%sys.argv[0] )
        sys.exit()

try:
    (opts,pargs) = getopt.getopt(sys.argv[1:],"vd",[])
except:
    errorMessage("%s: Invalid command line parameter" % sys.argv[0])
    describeSelf()
    sys.exit(1)

Message(idString)

# Parse args
for a,b in opts:
    if a=="-v":
        setVerbose(1)
    if a=="-d":
        setDebug(1)

smartreaderFlags= ""
if getDebug() or getVerbose():
    smartreaderFlags= smartreaderFlags+ "-verbose "
flags= " "
if getDebug():
    flags= flags+"-d "
if getVerbose():
    flags= flags+"-v "

# Make the subdirectories of links, doing unpacking, etc. as need.
os.system("renumber_ifiles.py %s"%flags)

# Convert the resulting files to Pgh MRI format
for series in xrange(1,10):
    linkdir= "series%d_links"%series
    tfile= "series%d_with_guide"%series
    ofile= "series%d"%series
    if os.access(linkdir,os.R_OK):
        Message("#####################################################")
        Message("Converting series %d"%series)
        safeRun("smartreader %s -i %s/I.#### -multi -out %s"%\
                (smartreaderFlags,linkdir,tfile))
        tmpDS= MRIDataset(tfile)
        imageChunk= tmpDS.getChunk("images")
        if imageChunk.hasValue("extent.t"):
            tdim= imageChunk.getDim("t")
        else:
            tdim= 1
        if tdim > 1:
            keept = tdim - 1
            Message("Dropping guide volume from %d; %d times remain"%\
                    (series,keept))
            safeRun("mri_subset -d t -len %d -s 1 %s %s"%(keept,tfile,ofile))
        else:
            Message("Series %d is a single volume"%series)
            safeRun("mri_copy_dataset %s %s"%(tfile,ofile))
        safeRun("mri_destroy_dataset %s"%tfile)

