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
if os.environ.has_key("FIASCO"):
    sys.path.append(os.environ["FIASCO"])
from fiasco_utils import *

idString= "$Id: paste_split.py,v 1.4 2006/05/04 23:02:58 welling Exp $"

#################
#
# Main
#
#################

# Check for "-help"
if len(sys.argv)>1:
    if sys.argv[1] == "-help":
	if len(sys.argv)>2:
	    os.system( "scripthelp %s %s"%(sys.argv[0],sys.argv[2]) );
	else:
	    os.system( "scripthelp %s"%sys.argv[0] );
	sys.exit();

try:
    (opts,pargs) = getopt.getopt(sys.argv[1:],"vd",\
                                 ["out=","nslices=","nimages="])
except:
    print "%s: Invalid command line parameter" % sys.argv[0]
    describeSelf();
    sys.exit()

#Check calling syntax; parse args
if len(pargs) < 1 :
    describeSelf()
    sys.exit(1)

outFname= None
nSlices= 0
nImages= 0
for a,b in opts:
    if a=="-v":
        setVerbose(1)
    if a=="-d":
        setDebug(1)
    if a=="--out":
        outFname= b
    if a=="--nslices":
        nSlices= int(b)
    if a=="--nimages":
        nImages= int(b)

if outFname==None:
    sys.exit("Required output file name not given")
if nSlices==0:
    sys.exit("Requied number of slices not given or invalid")
if nImages==0:
    sys.exit("Requied number of images not given or invalid")

fileLines= []
condDict= {}
condStringTable= []
factors= []
factorLevels= {}
factorIndexDict= {}
factorLevelDicts= {}
condMatrix= CondMatrix(nImages,nSlices)

# Slurp the files, and scan for overall set of factors
for inFname in pargs:
    preparseSplitFile( inFname, fileLines, factorIndexDict, factors )

parseAllSplitFiles( fileLines, factors, factorLevels, factorIndexDict,
                    factorLevelDicts, condMatrix, condDict,
                    condStringTable, (False) )
            
verboseMessage("Factors: %s"%factors)
verboseMessage("condDict: %s"%condDict)

ofile= open(outFname,"w")
ofile.write("%d\n"%len(factors))
for fac in factors:
    ofile.write("%s "%fac)
ofile.write("\n")
for thisCondString in condStringTable:
    ofile.write("%s 0 0\n"%thisCondString)
for t in xrange(0,nImages):
    zStart= 0
    for z in xrange(0,nSlices):
        if condMatrix.getCond(t,z)==condMatrix.getCond(t,zStart):
            pass
        else:
            if zStart==z-1:
                ofile.write("%s %d %d\n"%\
                            (condStringTable[condMatrix.getCond(t,zStart)],
                             t,zStart))
            else:
                ofile.write("%s %d %d-%d\n"%\
                            (condStringTable[condMatrix.getCond(t,zStart)],
                             t,zStart,z-1))
            zStart= z
    if zStart==nSlices-1:
        ofile.write("%s %d %d\n"%\
                    (condStringTable[condMatrix.getCond(t,zStart)],
                     t,zStart))
    else:
        ofile.write("%s %d %d-%d\n"%\
                    (condStringTable[condMatrix.getCond(t,zStart)],
                     t,zStart,nSlices-1))

ofile.close()

# Clean up
# Nothing to do here.

