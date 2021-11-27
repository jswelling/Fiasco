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

idString= "$Id: cowarp_inplane.py,v 1.5 2006/05/04 23:02:58 welling Exp $"

tol= 0.001

generatedScriptHdr="""\
#!/bin/csh -ef
#
# This is an automatically generated script!  It's only good
# for this particular special case.  Don't copy it around or
# edit it unless you *really* know what you're doing!
#
# Inputs:
#   $1: input structural MRI file name
#   $2: output structural MRI, warped to match shape and
#       location of functionals
#
# NOTE: the input file must have consistent voxel size and
#       bounds information.
#
cat > tmp_cowarp_inpl.par << EOF
"""
generatedScriptBody="""\
EOF
set xvox = `mri_printfield -fld images.voxel_spacing.x $1`
set yvox = `mri_printfield -fld images.voxel_spacing.y $1`
set zvox = `mri_printfield -fld images.voxel_spacing.z $1`
iwarp -interp bspline -x $xvox -y $yvox -z $zvox -i $1 -h $2 \
  -p tmp_cowarp_inpl.par
rm tmp_cowarp_inpl.par
"""

def pickMagicPadSize( sz ):
    for magic in magicSizes:
        if magic>=sz:
            return magic-sz;

def getMean(chunk):
    totDim= chunk.getDim('x')*chunk.getDim('y')*chunk.getDim('z')
    safeRun("mri_copy_dataset %s tmp_mean_1"%chunk.ds.fname)
    safeRun("mri_remap -order q -length %d tmp_mean_1"%totDim)
    safeRun("mri_subsample -d q -length 1 -mean tmp_mean_1 tmp_mean_2")
    lines= readCmdOutputToList("mri_rpn_math '$1,1,if_print_1' tmp_mean_2")
    return float(lines[0])

def makeGradMag3D(ch,vox,dsname):
    safeRun("mri_smooth -d x -smoother_type ddx %s tmp_x"%ch.ds.fname)
    safeRun("mri_smooth -d y -smoother_type ddx %s tmp_y"%ch.ds.fname)
    safeRun("mri_smooth -d z -smoother_type ddx %s tmp_z"%ch.ds.fname)
    safeRun("mri_rpn_math -out %s '$1,dup,*,%f,dup,*,/,$2,dup,*,%f,dup,*,/,+,$3,dup,*,%f,dup,*,/,+,sqrt' tmp_x tmp_y tmp_z"%(dsname,vox[0],vox[1],vox[2]))
    newDS= MRIDataset(dsname)
    newChunk= newDS.getChunk('images')
    return BBox(newChunk)

def checkBBoxCornersMatch( bbox1, bbox2 ):
    for corner in ["tlf","blf","trf","brf","tlb","blb","trb","brb"]:
        if bbox1.__dict__.has_key(corner) and bbox2.__dict__.has_key(corner):
            loc1= bbox1.__dict__[corner]
            loc2= bbox1.__dict__[corner]
            if (loc1-loc2).mag() > tol:
                sys.exit("Bounding boxes don't align at corner %s!"%corner)

myAlgWords= ["gradmag","meanc"]

def simplifyAlgString( str ):
    words= string.split(str,",")
    goodWords= []
    for word in words:
        if not word in myAlgWords: goodWords.append(word)
    return string.join(goodWords,",")

def algRequiresGradMag( str ):
    return str.find("gradmag")>=0

def algRequiresMeanCorrect( str ):
    return str.find("meanc")>=0

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
    (opts,pargs) = getopt.getopt(sys.argv[1:],"vd",["oscript=","alg="])
except:
    print("%s: Invalid command line parameter" % sys.argv[0])
    describeSelf();
    sys.exit()

#Check calling syntax; parse args
if len(pargs) != 2 :
    describeSelf()
    sys.exit(1)

oscript= "apply_cowarp_inplane.gen_csh"
clAlgString= None
for a,b in opts:
    if a=="-v":
        setVerbose(1)
    if a=="-d":
        setDebug(1)
    if a=="--oscript":
        oscript= b
    if a=="--alg":
        clAlgString= b

funcDS= MRIDataset(os.path.abspath(pargs[0]))
funcChunk= funcDS.getChunk('images')
strctDS= MRIDataset(os.path.abspath(pargs[1]))
strctChunk= strctDS.getChunk('images')

# Move to a temporary directory 
tmpdir= makeTempDir('tmp_cowarp_inplane')
homedir= os.getcwd()
os.chdir(tmpdir)

# Get relevant dimensions
xdim= funcChunk.getDim("x");
ydim= funcChunk.getDim("y");
zdim= funcChunk.getDim("z");

#Check reasonableness of input
dimstr= funcChunk.getValue('dimensions');
if dimstr != "xyz":
    if dimstr == "xyzt":
        if funcChunk.getDim("t") != 1:
            sys.exit("Input file must have t extent 1!")
    elif dimstr == "vxyzt":
        if funcChunk.getDim("t") != 1:
            sys.exit("Input file must have t extent 1!")
        if funcChunk.getDim("v") != 1:
            sys.exit("Input file must have v extent 1!")
    elif dimstr == "vxyz":
        if funcChunk.getDim("v") != 1:
            sys.exit("Input file must have v extent 1!")
    else:
        sys.exit("Input file must have dimensions (v)xyz(t)!")

if clAlgString != None:
    algString= clAlgString
elif funcChunk.isT1Weighted() == strctChunk.isT1Weighted():
    algString= "opt=praxis,obj=mse,meanc,inplane"
else:
    algString= "opt=praxis,obj=jointentropy,inplane"

funcBBox= BBox(funcChunk)
strctBBox= BBox(strctChunk)
if getVerbose():
    funcBBox.printBounds("Aligned functional bounding box:")
    strctBBox.printBounds("Inplane structural bounding box:")
checkBBoxCornersMatch(funcBBox,strctBBox)

funcVox= [ funcChunk.getFloat('voxel_spacing.x'), \
           funcChunk.getFloat('voxel_spacing.y'), \
           funcChunk.getFloat('voxel_spacing.z') ];
verboseMessage("functional voxel: "+str(funcVox))
strctVox= [ strctChunk.getFloat('voxel_spacing.x'), \
           strctChunk.getFloat('voxel_spacing.y'), \
           strctChunk.getFloat('voxel_spacing.z') ];
verboseMessage("structural voxel: "+str(strctVox))
for i in range(0,3):
    if abs(funcVox[i]-strctVox[i]) > tol:
        sys.exit("Structural and functional voxel sizes do not match!")

moveThisBBox= strctBBox
alignToBBox= funcBBox

if algRequiresGradMag(algString):
    verboseMessage("Producing grad magnitudes...")
    moveThisBBox= moveThisBBox.chunk.makeGradMag3D(strctVox,'strctGrad')
    alignToBBox= alignToBBox.chunk.makeGradMag3D(funcVox,'funcGrad')
    
if algRequiresMeanCorrect(algString):
    verboseMessage("Mean correcting...")
    moveThisMean= moveThisBBox.chunk.getMean()
    alignToMean= alignToBBox.chunk.getMean()
    safeRun("mri_rpn_math -out alignto '$1,%f,*' %s"%\
            ((moveThisMean/alignToMean),alignToBBox.chunk.ds.fname))
    alignToBBox= bboxFromFilename("alignto")

flags= "-x %f -y %f -z %f"%tuple(funcVox)
if getVerbose():
    flags= flags+' -v'
if getDebug():
    flags= flags+' -debug'
if algString != None:
    flags= flags+" -alg %s"%simplifyAlgString(algString)
safeRun("mri_remap -order xyzt %s"%moveThisBBox.chunk.ds.fname)
verboseMessage("Beginning co-warp")
safeRun("parallel.run.csh estiwarp %s -input %s -align %s -p warp.par"%\
        (flags, moveThisBBox.chunk.ds.fname, alignToBBox.chunk.ds.fname))
verboseMessage("Finished co-warp")

# All done- now we write a shell script that implements this transformation.
scriptName= os.path.join(homedir,oscript)
ofile= open(scriptName,"w")
ofile.write(generatedScriptHdr)
ifile= open("warp.par","r")
parLines= ifile.readlines()
ifile.close()
for line in parLines:
    ofile.write(line)
ofile.write(generatedScriptBody)
ofile.close()
safeRun("chmod ug+x %s"%scriptName)

# Clean up
os.chdir(homedir)
if not getDebug():
    removeTmpDir(tmpdir)


