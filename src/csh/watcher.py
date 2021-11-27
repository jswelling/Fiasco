#! /usr/bin/env python
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
# Notes:
#
import sys
import os
import os.path
import string
import getopt
import time
if "FIASCO" in os.environ:
    sys.path.append(os.environ["FIASCO"])
import fiasco_utils

def describeSelf():
    print("usage [-i intervalInSeconds] pidToWatch")

idString= "$Id: watcher.py,v 1.1 2003/08/07 17:39:44 welling Exp $"

sleepInterval= 30

try:
    (opts,pargs) = getopt.getopt(sys.argv[1:],"i:",[])
except:
    fiasco_utils.errorMessage("%s: Invalid command line parameter" % sys.argv[0])
    describeSelf()
    sys.exit(1)

if len(pargs) != 1:
    describeSelf()
    sys.exit(1)

for a,b in opts:
    if a=="-i":
        sleepInterval= string.atoi(b)

topPID = string.atoi(pargs[0])
myPID = os.getpid()

while 1:
    psout= fiasco_utils.readCmdOutputToList("ps -elf")

    linesByPid= {}
    linesByParent= {}

    psout.pop(0) # get rid of the headers
    for ln in psout:
        toks= string.split(ln,None,14)
        thisPID= string.atoi(toks[3])
        if thisPID != myPID:
            linesByPid[thisPID]= toks
            linesByParent[string.atoi(toks[4])]= toks

    if not linesByPid.has_key(topPID):
        break
    
    thisPID= topPID
    thisLine= []
    while 1:
        if linesByParent.has_key(thisPID):
            thisLine= linesByParent[thisPID]
            thisPID= string.atoi(thisLine[3])
        else:
            break
    if os.access("/proc/%d/stat"%thisPID,os.R_OK):
        statfile= open("/proc/%d/stat"%thisPID)
        statinfo= string.split(statfile.readline())
        statfile.close()
        if len(thisLine)>0:
            print("cmd=<%s> status=%s"%(string.strip(thisLine[len(thisLine)-1]),statinfo[2]))
        else:
            print("cmd=%s status=%s"%(statinfo[1],statinfo[2]))
    else:
        if len(thisLine)>0:
            print(thisLine[len(thisLine)-1])
    sys.stdout.flush()
    time.sleep(sleepInterval)

