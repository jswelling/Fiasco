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
# *     Copyright (c) 2006 Department of Statistics,         *
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
#################
# Notes-
#################

import sys
import os
import os.path
import string
import getopt
from math import *
if "FIASCO" in os.environ:
    sys.path.append(os.environ["FIASCO"])
from fiasco_utils import *
import mripipes

def initCallback(dim,dimstr,fastBlk,upstream_extent,extent1,extent2,
                 slowBlk,hookDict):
##     print("initCallback Args: %s <%s> %d %d %d %d %d %s"%\
##           (dim,dimstr,fastBlk,upstream_extent,extent1,extent2,slowBlk,hookDict))
    if (len(hookDict['blockOffsetList'])%extent2 != 0):
        print("Length of block offset list does not match extent!")
        return 0
    hookDict['fastBlk']=fastBlk*extent1
    hookDict['slowBlk']=slowBlk
    hookDict['extent']= extent2
    hookDict['upstreamBlk']= fastBlk*upstream_extent
    return 1

def myCallback(size,offset,hookDict):
    blockOffsetList= hookDict['blockOffsetList']
    extent= hookDict['extent']
    fastBlk= hookDict['fastBlk']
    upstreamBlk= hookDict['upstreamBlk']
    n_fast_blks= offset/fastBlk
    fast_blk_offset= offset%fastBlk
    n_full_extents= n_fast_blks/extent
    extent_offset= n_fast_blks%extent
##     print("callback: %d %d %d %d"%\
##           (n_fast_blks,fast_blk_offset,n_full_extents,extent_offset))
    return (fastBlk-fast_blk_offset,
            n_full_extents*upstreamBlk
            + blockOffsetList[extent_offset]+fast_blk_offset)

#################
#
# Main
#
#################

infoDict= {}

# Check for "-help"
if len(sys.argv)>1:
    if sys.argv[1] == "-help":
        if len(sys.argv)>2:
            os.system( "scripthelp %s %s"%(sys.argv[0],sys.argv[2]) );
        else:
            os.system( "scripthelp %s"%sys.argv[0] );
        sys.exit();

extent= None
dim= None
newdim= None
offsetFname= None
inFname= None
outFname= None

try:
    (opts,pargs) = getopt.getopt(sys.argv[1:],"v",
                                 ["len=","dim=","newdim=",
                                  "debug","offsetfile="])
except:
    print("%s: Invalid command line parameter" % sys.argv[0])
    describeSelf();
    sys.exit()

#Check calling syntax; parse args
if len(pargs) != 2 :
    describeSelf()
    sys.exit(1)

for a,b in opts:
    if a=="-v":
        setVerbose(1)
    if a=="--debug":
        setDebug(1)
    if a=="--len":
        extent= int(b)
    if a=="--dim":
        dim= b
    if a=="--newdim":
        newdim= b
    if a=="--offsetfile":
        offsetFname=b

if not extent:
    sys.exit("The required --len option is not present")
if not dim:
    sys.exit("The required --dim option is not present")
if not newdim:
    sys.exit("The required --newdim option is not present")
if not offsetFname:
    sys.exit("The required --offsetfile option is not present")
if dim[0]==newdim[0]:
    sys.exit("dim and newdim must differ");

inFname= os.path.abspath(pargs[0])
outFname= os.path.abspath(pargs[1])

verboseMessage("Filtering %s into %s based on block list %s"%\
               (pargs[0],pargs[1],offsetFname))

offsetFile= open(offsetFname,"r")
lines= offsetFile.readlines();
offsetFile.close();
offsetList= [ int(x) for x in lines ]

hookDict= {'blockOffsetList':offsetList}

a= mripipes.createArena();
t1= mripipes.createMRIFileInputTool(a,inFname)
t2= mripipes.createBlockMapTool(a,dim,newdim,extent,len(offsetList),
                                (initCallback,myCallback,hookDict))
if getDebug(): t2.setDebug()
if getVerbose(): t2.setVerbose()
t3= mripipes.createMRIFileOutputTool(a,outFname)
t3.getSink(0).connect(t2.getSource(0))
t2.getSink(0).connect(t1.getSourceByName('images'))
a.init()
a.execute()
    
