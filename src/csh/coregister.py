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

idString= "$Id: coregister.py,v 1.13 2006/05/04 23:02:58 welling Exp $"

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
                                 ["warpalg=","structalg=",\
                                  "inplanealg="])
except:
    print("%s: Invalid command line parameter" % sys.argv[0])
    describeSelf();
    sys.exit()

#Check calling syntax; parse args
if len(pargs) != 3 :
    describeSelf()
    sys.exit(1)

warpAlgString= None
structAlgString= None
inplAlgString= None
for a,b in opts:
    if a=="-v":
        setVerbose(1)
    if a=="-d":
        setDebug(1)
    if a=="--warpalg":
        warpAlgString= b;
    if a=="--structalg":
        structAlgString= b;
    if a=="--inplanealg":
        inplAlgString= b;

funcBBox= bboxFromFilename(os.path.abspath(pargs[0]))
funcChunk= funcBBox.chunk
funcVox= funcBBox.vox
inplaneBBox= bboxFromFilename(os.path.abspath(pargs[1]))
inplaneChunk= inplaneBBox.chunk
inplaneVox= inplaneBBox.vox
strctBBox= bboxFromFilename(os.path.abspath(pargs[2]))
strctChunk= strctBBox.chunk
strctVox= strctBBox.vox

#Check reasonableness of input
for thisChunk in (funcBBox.chunk, inplaneBBox.chunk, strctBBox.chunk):
    dimstr= thisChunk.getValue('dimensions');
    if dimstr != "xyz":
        if dimstr == "xyzt":
            if thisChunk.getDim("t") != 1:
                sys.exit("Input file must have t extent 1!")
        elif dimstr == "vxyzt":
            if thisChunk.getDim("t") != 1:
                sys.exit("Input file must have t extent 1!")
            if thisChunk.getDim("v") != 1:
                sys.exit("Input file must have v extent 1!")
        elif dimstr == "vxyz":
            if thisChunk.getDim("v") != 1:
                sys.exit("Input file must have v extent 1!")
        else:
            sys.exit("Input file must have dimensions (v)xyz(t)!")

# Move to a temporary directory 
tmpdir= makeTempDir('tmp_coregister')
homedir= os.getcwd()
os.chdir(tmpdir)

flags= ""
if (getVerbose()):
    flags= flags+"-v "

if (getDebug()):
    flags= flags+"-d "

verboseMessage("#### Coregistering functionals to inplanes ####")
if inplAlgString==None:
    safeRun("coregister_inplane.py %s %s %s"%\
            (flags,funcBBox.chunk.ds.fname,inplaneBBox.chunk.ds.fname))
else:
    safeRun("coregister_inplane.py %s --alg %s %s %s"%\
            (flags,inplAlgString,funcBBox.chunk.ds.fname,\
             inplaneBBox.chunk.ds.fname))
verboseMessage("#### Coregistering axial structurals to inplanes ####")
if structAlgString==None:
    safeRun("coregister_struct_to_inplane.py %s %s %s"%\
            (flags,inplaneBBox.chunk.ds.fname,strctBBox.chunk.ds.fname))
else:
    safeRun("coregister_struct_to_inplane.py %s --alg %s %s %s"%\
            (flags,structAlgString,\
             inplaneBBox.chunk.ds.fname,strctBBox.chunk.ds.fname))
verboseMessage("#### Generating transformation to inplane coordinates ####")
safeRun("coregister_inplane_to_inplane.py %s %s"%\
        (flags,inplaneBBox.chunk.ds.fname))

# Use the newly generated scripts to produce datasets for cowarping
verboseMessage("#### Generating datasets for cowarping ####")
safeRun("apply_coreg_inplane.gen_csh %s alignedFunc resampFunc"%\
        funcBBox.chunk.ds.fname);
scaleX= (float(funcChunk.getDim("x")-1)*funcBBox.vox[0])/\
        (float(inplaneChunk.getDim("x")-1)*inplaneBBox.vox[0])
scaleY= (float(funcChunk.getDim("y")-1)*funcBBox.vox[1])/\
        (float(inplaneChunk.getDim("y")-1)*inplaneBBox.vox[1])
xLow= 0.5*(1.0-scaleX)*float(inplaneChunk.getDim('x')-1)
xHigh= float(inplaneChunk.getDim('x')-1)-xLow
yLow= 0.5*(1.0-scaleY)*float(inplaneChunk.getDim('y')-1)
yHigh= float(inplaneChunk.getDim('y')-1)-yLow
safeRun("mri_resample -d x -len %d -start %f -end %f -interp bspline %s tmp"%\
        (funcBBox.chunk.getDim('x'),xLow,xHigh,
         inplaneBBox.chunk.ds.fname))
safeRun("mri_resample -d y -len %d -start %f -end %f -interp bspline tmp lowresInplane"%\
        (funcBBox.chunk.getDim('y'),yLow,yHigh))
lowresInplaneBBox= bboxFromFilename("lowresInplane")
lowresInplaneBBox.setVox(funcBBox.vox)
lowresInplaneBBox.setCorners(Vec4(1.0,0.0,0.0), Vec4(0.0,1.0,0.0), \
                             Vec4(0.0,0.0,1.0))
lowresInplaneBBox.exportBounds()

# Generate the cowarp script
verboseMessage("#### Cowarping inplanes to functionals ####")
if warpAlgString==None:
    safeRun("cowarp_inplane.py %s alignedFunc lowresInplane"%flags)
else:
    safeRun("cowarp_inplane.py %s --alg %s alignedFunc lowresInplane"%\
            (flags,warpAlgString))

verboseMessage("#### Cleaning up ####")
# Move all the generated scripts back to the home directory
safeRun("cp *.gen_csh %s"%homedir)

# Clean up
os.chdir(homedir)
if not getDebug():
    removeTmpDir(tmpdir)


