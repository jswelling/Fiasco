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
if os.environ.has_key("FIASCO"):
    sys.path.append(os.environ["FIASCO"])
from fiasco_utils import *

idString= "$Id: dataset_matches_env.py,v 1.3 2003/10/02 01:14:48 welling Exp $"

chunk= "images" # by default

intPairs= [ ("F_NIMAGE", "extent.t"), ("F_NSLICE", "extent.z") ]
floatPairs= [ ("F_XVOXEL", "voxel_spacing.x"), \
              ("F_YVOXEL", "voxel_spacing.y"), \
              ("F_ZVOXEL", "voxel_spacing.z") ]

def compareFloat( envName, file, chunk, key ):
    valStr= getFieldNofail(file,chunk,key);
    if valStr == None:
        debugMessage( "%s has no key %s.%s"%(file,chunk,key) )
        return 1
    else:
        if os.environ.has_key(envName):
            val1= string.atof(os.environ[envName])
            val2= string.atof(valStr)
            return ( val1 == val2 )
        else:
            debugMessage( "Environment variable %s is not set!"%envName);
            return 0

def compareInt( envName, file, chunk, key ):
    valStr= getFieldNofail(file,chunk,key);
    if valStr == None:
        debugMessage( "%s has no key %s.%s"%(file,chunk,key) )
        return 1
    else:
        if os.environ.has_key(envName):
            val1= string.atoi(os.environ[envName])
            val2= string.atoi(valStr)
            return ( val1 == val2 )
        else:
            debugMessage( "Environment variable %s is not set!"%envName);
            return 0

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
    (opts,pargs) = getopt.getopt(sys.argv[1:],"vdc:",[])
except:
    errorMessage("%s: Invalid command line parameter" % sys.argv[0])
    describeSelf()
    sys.exit(1)

Message(idString)

# Check calling syntax; parse args
if len(pargs) != 1 and len(pargs) != 3:
    describeSelf()
    sys.exit(1)

for a,b in opts:
    if a=="-v":
        setVerbose(1)
    if a=="-d":
        setDebug(1)
    if a=="-c":
        chunk= b

inDS= pargs[0]

for envStr, key in intPairs:
    debugMessage("Checking %s against %s.%s"%(envStr,chunk,key))
    if not compareInt( envStr, inDS, chunk, key ):
        verboseMessage("%s does not match %s.%s!"%(envStr,chunk,key))
        sys.exit(1)

for envStr, key in floatPairs:
    debugMessage("Checking %s against %s.%s"%(envStr,chunk,key))
    if not compareFloat( envStr, inDS, chunk, key ):
        verboseMessage("%s does not match %s.%s!"%(envStr,chunk,key))
        sys.exit(1)

debugMessage("Checking F_READER_SPACE against %s.extent.v"%chunk)
vDimStr= getFieldNofail(inDS, chunk, "extent.v")
if vDimStr==None:
    vDim= 1
else:
    vDim= int(vDimStr)
if os.environ.has_key("F_READER_SPACE"):
    space= os.environ["F_READER_SPACE"]
    if space == "i":
        if vDim != 1:
            verboseMessage("F_READER_SPACE is i but data is not scalar!")
            sys.exit(1)
    elif space == "k":
        if vDim != 2:
            verboseMessage("F_READER_SPACE is k but data is not complex!")
            sys.exit(1)
    else:
        verboseMessage("F_READER_SPACE is neither i nor k!")
        sys.exit(1)

                       
                 
