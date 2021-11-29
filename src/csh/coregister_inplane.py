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

idString= "$Id: coregister_inplane.py,v 1.19 2006/05/04 23:02:58 welling Exp $"

magicSizes= [ 2, 4, 8, 16, 32, 48, 64, 96, 128, 192, 256 ]

maxPermittedCtrShift= 10.0   # Max allowed sep btwn func and inplane ctrs, mm

generatedScriptHdr="""\
#!/bin/csh -ef
#
# This is an automatically generated script!  It's only good
# for this particular special case.  Don't copy it around or
# edit it unless you *really* know what you're doing!
#
# Switches: (must appear first)
#   -closest : use 'closest' interpolation rather than Fourier interpolation
#
# Inputs:
#   $1: input functional MRI file name
#   $2: output functional MRI, translated and rotated to
#       orientation of the inplane
#   $3: $2 resampled up to the in-plane structural resolution
#
set modeflag = '-interp fourier'
if ( "$1" == '-closest' ) then
  shift argv
  set modeflag = '-interp closest'
endif
cat > tmp_coreg_inpl.par << EOF
"""
generatedScriptBody="""\
EOF
ireg3d $modeflag -x %f -y %f -z %f -i $1 -h $2 -p tmp_coreg_inpl.par
mri_resample -interp closest -d x -len %d -start %f -end %f $2 tmp_coreg_inpl_1
mri_resample -interp closest -d y -len %d -start %f -end %f tmp_coreg_inpl_1 $3

mri_setfield -field images.tlf -all012 -value ' %f,%f,%f' $2
mri_setfield -field images.trf -all012 -value ' %f,%f,%f' $2
mri_setfield -field images.tlb -all012 -value ' %f,%f,%f' $2
mri_setfield -field images.trb -all012 -value ' %f,%f,%f' $2
mri_setfield -field images.blf -all012 -value ' %f,%f,%f' $2
mri_setfield -field images.brf -all012 -value ' %f,%f,%f' $2
mri_setfield -field images.blb -all012 -value ' %f,%f,%f' $2
mri_setfield -field images.brb -all012 -value ' %f,%f,%f' $2
mri_setfield -field images.ctr -all012 -value ' %f,%f,%f' $2
mri_setfield -field images.voxel_spacing -allxyz -value ' %f,%f,%f' $3
mri_setfield -field images.voxel_size -allxyz -value ' %f,%f,%f' $3
mri_setfield -field images.tlf -all012 -value ' %f,%f,%f' $3
mri_setfield -field images.trf -all012 -value ' %f,%f,%f' $3
mri_setfield -field images.tlb -all012 -value ' %f,%f,%f' $3
mri_setfield -field images.trb -all012 -value ' %f,%f,%f' $3
mri_setfield -field images.blf -all012 -value ' %f,%f,%f' $3
mri_setfield -field images.brf -all012 -value ' %f,%f,%f' $3
mri_setfield -field images.blb -all012 -value ' %f,%f,%f' $3
mri_setfield -field images.brb -all012 -value ' %f,%f,%f' $3
mri_setfield -field images.ctr -all012 -value ' %f,%f,%f' $3

mri_destroy_dataset tmp_coreg_inpl_1
rm tmp_coreg_inpl.par
"""

def pickMagicPadSize( sz ):
    for magic in magicSizes:
        if magic>=sz:
            return magic-sz;

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
    (opts,pargs) = getopt.getopt(sys.argv[1:],"vd",["oscript=","maxctrsep=",\
                                                    "alg="])
except:
    print("%s: Invalid command line parameter" % sys.argv[0])
    describeSelf();
    sys.exit()

#Check calling syntax; parse args
if len(pargs) != 2 :
    describeSelf()
    sys.exit(1)

oscript= "apply_coreg_inplane.gen_csh"
clAlgString= None
for a,b in opts:
    if a=="-v":
        setVerbose(1)
    if a=="-d":
        setDebug(1)
    if a=="--oscript":
        oscript= b
    if a=="--maxctrsep":
        maxPermittedCtrShift= float(b)
    if a=="--alg":
        clAlgString= b


funcDS= MRIDataset(os.path.abspath(pargs[0]))
funcChunk= funcDS.getChunk('images')
strctDS= MRIDataset(os.path.abspath(pargs[1]))
strctChunk= strctDS.getChunk('images')

# Move to a temporary directory 
tmpdir= makeTempDir('tmp_coregister_inplane')
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
elif funcChunk.isT1Weighted()==strctChunk.isT1Weighted():
    algString= "obj=mse,opt=praxis,meanc"
else:
    algString= "opt=nelmin,obj=jointentropy"

funcBBox= BBox(funcChunk)
strctBBox= BBox(strctChunk)
if getVerbose():
    funcBBox.printBounds("Functional bounding box:")
    strctBBox.printBounds("Inplane structural bounding box:")

# By definition, the bounding volume of the aligned functional now
# matches that of the inplane anatomical dataset.  We also want to
# think of these datasets as being in inplane coordinates, such that
# the Fiasco Z direction corresponds to the slice normal direction.
# We implement this by setting the BBoxes of the structural and
# functional appropriately, so that they will be written out
# that way when we generate the scripts. Note that we do not export
# these bounds, since we want the input files to keep their original
# corners!
funcBBox.setCorners(Vec4(1.0,0.0,0.0), Vec4(0.0,1.0,0.0), Vec4(0.0,0.0,1.0))
strctBBox.setCorners(Vec4(1.0,0.0,0.0), Vec4(0.0,1.0,0.0), Vec4(0.0,0.0,1.0))

funcTlf= funcBBox.tlf
funcBlf= funcBBox.blf
funcTrf= funcBBox.trf
funcBrf= funcBBox.brf
funcTlb= funcBBox.tlb
funcBlb= funcBBox.blb
funcTrb= funcBBox.trb
funcBrb= funcBBox.brb
funcCtr= funcBBox.ctr
strctCtr= strctBBox.ctr
verboseMessage("Separation between centers is %f mm"%(funcCtr-strctCtr).mag())
ctrShift= funcCtr - strctCtr
verboseMessage("Center shift: %s"%ctrShift)
if ctrShift.mag() > maxPermittedCtrShift:
    sys.exit("Structural and functional centers don't line up!")

funcVox= [ funcChunk.getFloat('voxel_spacing.x'), \
           funcChunk.getFloat('voxel_spacing.y'), \
           funcChunk.getFloat('voxel_spacing.z') ];
verboseMessage("functional voxel: "+str(funcVox))
strctVox= [ strctChunk.getFloat('voxel_spacing.x'), \
           strctChunk.getFloat('voxel_spacing.y'), \
           strctChunk.getFloat('voxel_spacing.z') ];
verboseMessage("structural voxel: "+str(strctVox))
xStart= (funcTlf[0]-funcCtr[0])/strctVox[0] + (strctChunk.getDim('x')/2);
xEnd= (funcTrf[0]-funcCtr[0])/strctVox[0] + (strctChunk.getDim('x')/2);
safeRun("mri_resample -d x -len %d -start %f -end %f %s resampx"%\
        (funcChunk.getDim('x'),xStart,xEnd,strctChunk.ds.fname))
yStart= (funcTlf[1]-funcCtr[1])/strctVox[1] + (strctChunk.getDim('y')/2);
yEnd= (funcTlb[1]-funcCtr[1])/strctVox[1] + (strctChunk.getDim('y')/2);
safeRun("mri_resample -d y -len %d -start %f -end %f resampx strctResamp"%\
        (funcChunk.getDim('y'),yStart,yEnd))
strctResampBBox= bboxFromFilename('strctResamp')

pad= pickMagicPadSize(funcChunk.getDim('z'))

safeRun("mri_pad -d z -len %d -shift %d %s funcPad"%\
        (zdim+pad, pad/2, funcChunk.ds.fname))
safeRun("mri_pad -d z -len %d -shift %d %s strctResampPad"%\
        (zdim+pad, pad/2, strctResampBBox.chunk.ds.fname))

funcPadBBox= bboxFromFilename("funcPad")
strctResampPadBBox= bboxFromFilename("strctResampPad")
moveThisBBox= funcPadBBox
alignToBBox= strctResampPadBBox

if algRequiresGradMag(algString):
    verboseMessage("Producing grad magnitudes...")
    moveThisBBox= moveThisBBox.chunk.makeGradMag3D(funcVox,'funcGrad')
    alignToBBox= alignToBBox.chunk.makeGradMag3D(funcVox,'strctGrad')
    
if algRequiresMeanCorrect(algString):
    verboseMessage("Mean correcting...")
    moveThisMean= moveThisBBox.chunk.getMean()
    alignToMean= alignToBBox.chunk.getMean()
    safeRun("mri_rpn_math -out alignto '$1,%f,*' %s"%\
            ((moveThisMean/alignToMean),alignToBBox.chunk.ds.fname))
    alignToBBox= bboxFromFilename("alignto")

flags= "-x %f -y %f -z %f"%tuple(funcVox)
if getDebug():
    flags= flags+" -debug"
if algString != None:
    flags= flags+" -alg %s"%simplifyAlgString(algString)
verboseMessage("Beginning 3D alignment")
safeRun("parallel.run.csh estireg3d %s -input %s -align %s -p rot.par"%\
        (flags,moveThisBBox.chunk.ds.fname,alignToBBox.chunk.ds.fname))
verboseMessage("Finished 3D alignment")

vals= None
try:
    ifile= open("rot.par","r")
    parLines= ifile.readlines()
    ifile.close()
    for line in parLines:
        if line[0] != '#':
            vals= map(float,string.split(line))
            break
except:
    pass
if vals == None:
    sys.exit("Estireg3d exited without producing an alignment estimate!")

safeRun("displace3d -estimates disp.par -xvx %f -yvx %f -zvx %f "%\
        (funcVox[0], funcVox[1], funcVox[2]) \
        + "-weight %s rot.par"%alignToBBox.chunk.ds.fname);
dispVals= None
try:
    ifile= open("disp.par","r")
    dispLines= ifile.readlines()
    ifile.close()
    for line in dispLines:
        if line[0] != '#':
            dispVals= map(float,string.split(line))
            break
except:
    pass
if dispVals == None:
    sys.exit("displace3d exited without producing a displacement estimate!")

# Report the transformation we've estimated
Q_inv= Quat(vals[1],vals[2],vals[3],vals[4])
invAligningTrans= Transform()
quat_to_trans(invAligningTrans,Q_inv,vals[5],vals[6],vals[7])
aligningTrans= Transform()
if not trans_inverse(aligningTrans,invAligningTrans):
    sys.exit("Failed to invert the aligning transformation!")
Q= Quat(0.0,0.0,0.0,1.0)
trans_to_quat(Q, aligningTrans)
rawShift= Vec4( aligningTrans[3], aligningTrans[7], \
                aligningTrans[11], aligningTrans[15] )
unscaledShift= Vec4(rawShift[0]*funcVox[0],\
                     rawShift[1]*funcVox[1],\
                     rawShift[2]*funcVox[2])

Message("Geometrical transformation for functional/inplane relative motion:")
Message("Quaternion : %s"%Q)
Message("Translation  (mm): %s"%unscaledShift)
Message("Rot, trans, and mean disp from this transform: %f %f %f\n"%\
        (dispVals[1],dispVals[2],dispVals[3]));

#
# Uncomment this if you want to copy intermediate files back to
# the home directory
#
##safeRun("ireg3d -x %f -y %f -z %f -i %s -h rot1 -p rot.par"%\
##        (funcVox[0],funcVox[1],funcVox[2],funcChunk.ds.fname))
        
##safeRun("mri_resample -d x -len %d -start 0.0 -end %f rot1 fresampx"%\
##        (strctChunk.getDim('x'),funcChunk.getDim('x')-1.0))
##safeRun("mri_resample -d y -len %d -start 0.0 -end %f fresampx funcResamp"%\
##        (strctChunk.getDim('y'),funcChunk.getDim('y')-1.0))
##safeRun("mri_setfield -field %s.voxel_spacing -allxyz -value ' %f,%f,%f' funcResamp"%\
##        (funcChunk.name,strctVox[0],strctVox[1],strctVox[2]))
##safeRun("mri_setfield -field %s.voxel_size -allxyz -value ' %f,%f,%f' funcResamp"%\
##        (funcChunk.name,strctVox[0],strctVox[1],strctVox[2]))
##safeRun("cp *.par %s"%homedir)
##safeRun("mri_copy_dataset rot1 %s/funcAlignedInplane"%homedir)
##safeRun("mri_copy_dataset strctResamp %s/strctDownsamp"%homedir)
##safeRun("mri_copy_dataset funcResamp %s/funcResamp"%homedir)
##safeRun("mri_copy_dataset funcGradPad %s/funcGradPad"%homedir)
##safeRun("mri_copy_dataset strctResampGradPad %s/strctResampGradPad"%homedir)

# All done- now we write a shell script that implements this transformation.
scriptName= os.path.join(homedir,oscript);
ofile= open(scriptName,"w")
ofile.write(generatedScriptHdr)
ifile.close()
for line in parLines:
    ofile.write(line)

# Args for ireg3d
argList= funcVox[0:3]
# Args for mri_resample -d x
start= ( strctBBox.tlf[0]-funcBBox.tlf[0] )/funcVox[0];
end= (funcChunk.getDim('x')-1.0) \
     + ( strctBBox.trf[0]-funcBBox.trf[0] )/funcVox[0];
argList= argList + [ strctChunk.getDim('x'), start, end ]

# Args for mri_resample -d y; remember sign change because Y is
# sampled in reverse order
start= -( strctBBox.tlb[1]-funcBBox.tlb[1] )/funcVox[0];
end= (funcChunk.getDim('y')-1.0) \
     -( strctBBox.tlf[1]-funcBBox.tlf[1] )/funcVox[0];
argList= argList + [ strctChunk.getDim('y'), start, end ]

# Remember that these are the corners in inplane coordinates, not the
# scanner coordinate corners which actually appear in the inplane struct
# input file.
for corner in ('tlf', 'trf', 'tlb', 'trb', 'blf', 'brf', 'blb', 'brb', 'ctr'):
    vtx= funcBBox.__dict__[corner]
    argList= argList+[ vtx[0], vtx[1], vtx[2] ]
argList = argList+strctVox[0:3] 
argList = argList+strctVox[0:3] 
# Remember that these are the corners in inplane coordinates, not the
# scanner coordinate corners which actually appear in the inplane struct
# input file.
for corner in ('tlf', 'trf', 'tlb', 'trb', 'blf', 'brf', 'blb', 'brb', 'ctr'):
    vtx= strctBBox.__dict__[corner]
    argList= argList+[ vtx[0], vtx[1], vtx[2] ]

ofile.write(generatedScriptBody%tuple(argList))
ofile.close()
safeRun("chmod ug+x %s"%scriptName)

# Clean up
os.chdir(homedir)
if not getDebug():
    removeTmpDir(tmpdir)


