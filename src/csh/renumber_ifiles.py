#!/usr/bin/env python
# renumber_ifiles.py
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

idString= '# $Id: renumber_ifiles.py,v 1.2 2004/04/08 22:21:36 welling Exp $'

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

for series in xrange(1,10):
    if os.access("00%d"%series,os.R_OK):
        newdir = "series%d_links"%series
        os.mkdir(newdir)
        count = 0
        for block in [ 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 1, 3, 5, 7, 9 ]:
            subdir = "%02d%d"%(block,series)
            uncompressMsgFlag= None
            if os.access(subdir,os.R_OK):
                fList= os.listdir(subdir);
                fList.sort();
                for f in fList:
                    linkMe= None
                    if ifileRegex.match(f):
                        linkMe= f;
                    if gzIfileRegex.match(f):
                        if not uncompressMsgFlag:
                            Message("Uncompressing some files in %s"%subdir)
                            uncompressMsgFlag= 1
                        safeRun("gunzip %s"%os.path.join(subdir,f))
                        linkMe= f[:5]
                    if zIfileRegex.match(f):
                        if not uncompressMsgFlag:
                            Message("Uncompressing some files in %s"%subdir)
                            uncompressMsgFlag= 1
                        safeRun("uncompress %s"%os.path.join(subdir,f))
                        linkMe= f[:5]
                    if linkMe:
                        os.symlink(os.path.join("..",subdir,linkMe),
                                   os.path.join(newdir,"I.%04d"%count))
                        count= count+1
        Message("Series %s contains %d images total"%(series,count))
