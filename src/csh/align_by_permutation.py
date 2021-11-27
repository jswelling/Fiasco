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

idString= "$Id: align_by_permutation.py,v 1.6 2007/02/06 21:45:42 welling Exp $"

def checkInputStructure( chunk, unsafeFlag ):
    dimstr= chunk.getValue('dimensions');
    if dimstr != "xyz":
        if dimstr == "xyzt":
            if chunk.getDim("t") != 1 and not unsafeFlag:
                sys.exit("Input file %s must have t extent 1!"%\
                         os.path.basename(chunk.ds.fname))
        elif dimstr == "vxyzt":
            if chunk.getDim("t") != 1 and not unsafeFlag:
                sys.exit("Input file %s must have t extent 1!"%\
                         os.path.basename(chunk.ds.fname))
            if chunk.getDim("v") != 1:
                sys.exit("Input file %s must have v extent 1!"%\
                         os.path.basename(chunk.ds.fname))
        elif dimstr == "vxyz":
            if chunk.getDim("v") != 1:
                sys.exit("Input file %s must have v extent 1!"%\
                         os.path.basename(chunk.ds.fname))
        else:
            sys.exit("Input file %s must have dimensions (v)xyz(t)!"%\
                     os.path.basename(chunk.ds.fname))
    

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
    (opts,pargs) = getopt.getopt(sys.argv[1:],"vd",["out=","unsafe"])
except:
    print("%s: Invalid command line parameter" % sys.argv[0])
    describeSelf();
    sys.exit()

#Check calling syntax; parse args
if len(pargs) != 2 :
    describeSelf()
    sys.exit(1)

outDSName= None
unsafeFlag= 0
for a,b in opts:
    if a=="-v":
        setVerbose(1)
    if a=="-d":
        setDebug(1)
    if a=="--out":
        outDSName= b
    if a=="--unsafe":
        unsafeFlag= 1

if outDSName==None:
    sys.exit("Required output dataset name not given.")

inDS= MRIDataset(os.path.abspath(pargs[0]))
inChunk= inDS.getChunk('images')
protoDS= MRIDataset(os.path.abspath(pargs[1]))
protoChunk= protoDS.getChunk('images')

#Check reasonableness of input
checkInputStructure(inChunk,unsafeFlag)
checkInputStructure(protoChunk,unsafeFlag)

# Create a temporary directory 
tmpdir= makeTempDir('tmp_align_by_permutation')
homedir= os.getcwd()

# Get relevant dimensions
xdim= inChunk.getDim("x");
ydim= inChunk.getDim("y");
zdim= inChunk.getDim("z");
dimstr= inChunk.getValue('dimensions');

inBBox= BBox(inChunk)
protoBBox= BBox(protoChunk)
if getVerbose():
    inBBox.printBounds("Input bounding box:")
    protoBBox.printBounds("Prototype bounding box:")

inRHCoordTest= inBBox.zedge.dot(inBBox.xedge.cross(inBBox.yedge))
protoRHCoordTest= protoBBox.zedge.dot(protoBBox.xedge.cross(protoBBox.yedge))
if inRHCoordTest*protoRHCoordTest < 0.0 and not unsafeFlag:
    sys.exit("Input and prototype coord systems don't have same handedness!")

inAxes= { "x":inBBox.xedge.clone(), \
          "y":inBBox.yedge.clone(), \
          "z":inBBox.zedge.clone() }
for v1 in ["x","y","z"]: inAxes[v1].normalize()
protoAxes= { "x":protoBBox.xedge.clone(), \
             "y":protoBBox.yedge.clone(), \
             "z":protoBBox.zedge.clone() }
for v1 in ["x","y","z"]: protoAxes[v1].normalize()

becomesMap= {}
usedToBeMap= {}
needsReversed= {}
for v1 in ["x","y","z"]: 
    largestDot= 0.0;
    comp= None
    for v2 in ["x","y","z"]:
        val= inAxes[v1].dot(protoAxes[v2])
        if math.fabs(val)>math.fabs(largestDot):
            largestDot= val
            comp= v2
    debugMessage("%s matches %s, dot %f"%(v1,comp,largestDot))
    becomesMap[v1]= comp
    needsReversed[v1]= ( largestDot < 0 )
debugMessage("becomesMap: %s"%repr(becomesMap))
for v1 in becomesMap.keys():
    usedToBeMap[becomesMap[v1]]= v1
debugMessage("usedToBeMap: %s"%repr(usedToBeMap))
debugMessage("needsReversed: %s"%repr(needsReversed))
debugMessage("inAxes: %s"%repr(inAxes))
debugMessage("protoAxes: %s"%repr(protoAxes))
newDimstr= usedToBeMap['x']+usedToBeMap['y']+usedToBeMap['z']
newExtents= "%d:%d:%d"%(inChunk.getDim(usedToBeMap['x']),\
                        inChunk.getDim(usedToBeMap['y']),\
                        inChunk.getDim(usedToBeMap['z']))
if dimstr.startswith('v'):
    newDimstr= 'v'+newDimstr
    newExtents= ":"+newExtents
if dimstr.endswith('t'):
    newDimstr= newDimstr+'t'
    newExtents= newExtents+":"
debugMessage("dimstr <%s> becomes <%s>, extents <%s>"%\
             (dimstr,newDimstr,newExtents))

# Flip the axis vectors as appropriate
for v1 in ['x','y','z']:
    if needsReversed[v1]: inAxes[v1]= -1.0*inAxes[v1]

# We will now use the needsReversed info to determine which data
# dimensions need to be flipped.  There is a trick here, since the
# Y data dimension is opposite the Y coordinate dimension in Fiasco
# coordinates.  Thus we first dink with the needsReversed info
# to correct for this.
if becomesMap['y'] != 'y':
    needsReversed[usedToBeMap['y']]= ( not needsReversed[usedToBeMap['y']] )
    needsReversed['y']= ( not needsReversed['y'] )
debugMessage("needsReversed after correction for data order: %s"%\
             repr(needsReversed))

# Handle axis reversals via the double-fft trick
currentDSName= inDS.fname
if needsReversed['x']:
    if needsReversed['y']:
        # use xy fft
        safeRun("mri_fft -d xy -fwd -cpx %s %s"%\
                (currentDSName,os.path.join(tmpdir,"tmp1")))
        safeRun("mri_fft -d xy -fwd -mod %s %s"%\
                (os.path.join(tmpdir,"tmp1"),os.path.join(tmpdir,"tmp2")))
        needsReversed['y']= 0
    else:
        # use x fft
        safeRun("mri_fft -d x -fwd -cpx %s %s"%\
                (currentDSName,os.path.join(tmpdir,"tmp1")))
        safeRun("mri_fft -d x -fwd -mod %s %s"%\
                (os.path.join(tmpdir,"tmp1"),os.path.join(tmpdir,"tmp2")))
    currentDSName= os.path.join(tmpdir,"tmp2")
    needsReversed['x']= 0
    if not dimstr.startswith('v'):
        safeRun("mri_remap -order %s %s"%(dimstr,currentDSName))
if needsReversed['y']:
    if needsReversed['z']:
        # use yz fft
        safeRun("mri_fft -d yz -fwd -cpx %s %s"%\
                (currentDSName,os.path.join(tmpdir,"tmp3")))
        safeRun("mri_fft -d yz -fwd -mod %s %s"%\
                (os.path.join(tmpdir,"tmp3"),os.path.join(tmpdir,"tmp4")))
        needsReversed['z']= 0
    else:
        # use y fft
        safeRun("mri_fft -d y -fwd -cpx %s %s"%\
                (currentDSName,os.path.join(tmpdir,"tmp3")))
        safeRun("mri_fft -d y -fwd -mod %s %s"%\
                (os.path.join(tmpdir,"tmp3"),os.path.join(tmpdir,"tmp4")))
    currentDSName= os.path.join(tmpdir,"tmp4")
    needsReversed['y']= 0
    if not dimstr.startswith('v'):
        safeRun("mri_remap -order %s %s"%(dimstr,currentDSName))
if needsReversed['z']:
    # use z fft
    safeRun("mri_fft -d z -fwd -cpx %s %s"%\
            (currentDSName,os.path.join(tmpdir,"tmp5")))
    safeRun("mri_fft -d z -fwd -mod %s %s"%\
            (os.path.join(tmpdir,"tmp5"),os.path.join(tmpdir,"tmp6")))
    currentDSName= os.path.join(tmpdir,"tmp6")
    needsReversed['z']= 0
    if not dimstr.startswith('v'):
        safeRun("mri_remap -order %s %s"%(dimstr,currentDSName))
debugMessage("inAxes now %s"%repr(inAxes))

if dimstr != newDimstr:
    safeRun("mri_permute -order %s %s %s"%(newDimstr,currentDSName,outDSName))
    safeRun("mri_remap -order %s -len %s %s"%(dimstr,newExtents,outDSName))
else:
    safeRun("mri_copy_dataset %s %s"%(currentDSName,outDSName))
outDS= MRIDataset(outDSName)
outChunk= outDS.getChunk('images')
outBBox= BBox(outChunk)
outBBox.setCtr(inBBox.ctr)
outBBox.setVox([inChunk.getFloat("voxel_spacing.%s"%usedToBeMap['x']),\
                inChunk.getFloat("voxel_spacing.%s"%usedToBeMap['y']),\
                inChunk.getFloat("voxel_spacing.%s"%usedToBeMap['z'])])
outBBox.setCorners(inAxes[usedToBeMap['x']],\
                   inAxes[usedToBeMap['y']],\
                   inAxes[usedToBeMap['z']])
if getVerbose():
    outBBox.printBounds("Output bounding box:")
outBBox.exportBounds()

# Clean up
os.chdir(homedir)
if not getDebug():
    removeTmpDir(tmpdir)


