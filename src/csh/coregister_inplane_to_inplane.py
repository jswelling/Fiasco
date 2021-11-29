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

import sys
import os
import os.path
import string
import getopt
from math import *
if "FIASCO" in os.environ:
    sys.path.append(os.environ["FIASCO"])
from fiasco_utils import *

idString= "$Id: coregister_inplane_to_inplane.py,v 1.4 2006/05/04 23:02:58 welling Exp $"

tol= 0.001

generatedScriptHdr="""\
#!/bin/csh -ef
#
# This is an automatically generated script!  It's only good
# for this particular special case.  Don't copy it around or
# edit it unless you *really* know what you're doing!
#
# Inputs:
#   $1: input inplane MRI file name (functional or matching anatomical)
#   $2: output inplane MRI file, identical to $1 but with corners
#       expressed in inplane coordinates
#
# NOTE: the input file must have consistent voxel size and
#       bounds information.
#
mri_copy_dataset $1 $2
"""

##############################
#
# Main
#
##############################

# Check for "-help"
if len(sys.argv)>1:
    if sys.argv[1] == "-help":
	if len(sys.argv)>2:
	    os.system( "scripthelp %s %s"%(sys.argv[0],sys.argv[2]) );
	else:
	    os.system( "scripthelp %s"%sys.argv[0] );
	sys.exit();

try:
    (opts,pargs) = getopt.getopt(sys.argv[1:],"vd",["oscript="])
except:
    print("%s: Invalid command line parameter" % sys.argv[0])
    describeSelf();
    sys.exit()

#Check calling syntax; parse args
if len(pargs) != 1 :
    describeSelf()
    sys.exit(1)

oscript= "apply_coreg_inplane_to_inplane.gen_csh"
for a,b in opts:
    if a=="-v":
        setVerbose(1)
    if a=="-d":
        setDebug(1)
    if a=="--oscript":
        oscript= b

inplaneBBox= bboxFromFilename(pargs[0])
inplaneChunk= inplaneBBox.chunk

#Check reasonableness of input
dimstr= inplaneBBox.chunk.getValue('dimensions');
if dimstr != "xyz":
    if dimstr == "xyzt":
	if inplaneChunk.getDim("t") != 1:
	    sys.exit("Input file must have t extent 1!")
    elif dimstr == "vxyzt":
        if inplaneChunk.getDim("t") != 1:
	    sys.exit("Input file must have t extent 1!")
        if inplaneChunk.getDim("v") != 1:
	    sys.exit("Input file must have v extent 1!")
    elif dimstr == "vxyz":
	if inplaneChunk.getDim("v") != 1:
	    sys.exit("Input file must have v extent 1!")
    else:
        sys.exit("Input file must have dimensions (v)xyz(t)!")

# Where are we?
homedir= os.getcwd()

# All we need to do is to reset corner coordinates.
inplaneBBox.printBounds("Inplane bounding box, scanner coordinates")
inplaneBBox.setCorners(Vec4(1.0,0.0,0.0), Vec4(0.0,1.0,0.0), \
                       Vec4(0.0,0.0,1.0))
inplaneBBox.printBounds("Inplane bounding box, inplane coordinates")

# Write a shell script that implements this change of coordinates
scriptName= os.path.join(homedir,oscript);
ofile= open(scriptName,"w")
ofile.write(generatedScriptHdr)
for corner in ('tlf', 'trf', 'tlb', 'trb', 'blf', 'brf', 'blb', 'brb', 'ctr'):
    vtx= inplaneBBox.__dict__[corner]
    ofile.write("mri_setfield -fld images.%s -all012 -value ' %f,%f,%f' $2\n"%\
                (corner,vtx[0], vtx[1], vtx[2]));
ofile.close()
safeRun("chmod ug+x %s"%scriptName)

# Clean up
# None to do for this script

