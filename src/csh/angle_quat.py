#! /usr/bin/env python
#
#
# ************************************************************
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
# ************************************************************
#

import sys
import os
import os.path
import string
import getopt
if "FIASCO" in os.environ:
    sys.path.append(os.environ["FIASCO"])
from fiasco_utils import *


# Check for "-help"
if len(sys.argv)>1:
    if sys.argv[1] == "-help":
	if len(sys.argv)>2:
	    os.system( "scripthelp %s %s"%(sys.argv[0],sys.argv[2]) );
	else:
	    os.system( "scripthelp %s"%sys.argv[0] );
        sys.exit();

try:
    (opts,pargs) = getopt.getopt(sys.argv[1:],"x:y:z:",["xangle=","yangle=","zangle="])
except:
    errorMessage("%s: Invalid command line parameter" % sys.argv[0])
    describeSelf();
    sys.exit()

xangle=0.0
yangle=0.0
zangle=0.0

for a,b in opts:
    if a in ("-x","--xangle"):
        xangle=float(b)
    if a in ("-y","--yangle"):
        yangle=float(b)
    if a in ("-z","--zangle"):
        zangle=float(b)

Q=quat_from_euler_RzRyRx(None,zangle,yangle,xangle)
print(Q.x, Q.y, Q.z, Q.w)
