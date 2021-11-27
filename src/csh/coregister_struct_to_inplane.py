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

idString= "$Id: coregister_struct_to_inplane.py,v 1.18 2006/05/04 23:02:58 welling Exp $"

magicSizes= [ 2, 4, 8, 16, 32, 48, 64, 96, 128, 192, 256, 320, 384 ]

parFileFormat= """\
##Format: order:index_t, type:raw
##Format: names:(3d_qbar_x,3d_qbar_y,3d_qbar_z,3d_qbar_w,3d_deltabarx,3d_deltabary,3d_deltabarz,mse)
# Input file: none
# Alignment file: none
# Stdv file: none
# voxel size x= %f, y= %f, z= %f
0 %f %f %f %f %f %f %f %f
"""

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
#   $1: input structural MRI file name
#   $2: output structural MRI, translated and rotated to
#       orientation of the inplane dataset, inplane coordinates
#   $3: $2 resampled to match in-plane structural, inplane coordinates
#
set modeflag = '-interp fourier'
if ( "$1" == '-closest' ) then
  shift argv
  set modeflag = '-interp closest'
endif
cat > tmp_coreg_strct_to_inplane.par << EOF
"""
generatedScriptBody="""\
EOF
ireg3d $modeflag -x %f -y %f -z %f -i $1 -h $2 -p tmp_coreg_strct_to_inplane.par
mri_setfield -field images.tlf -all012 -value ' %f,%f,%f' $2
mri_setfield -field images.trf -all012 -value ' %f,%f,%f' $2
mri_setfield -field images.tlb -all012 -value ' %f,%f,%f' $2
mri_setfield -field images.trb -all012 -value ' %f,%f,%f' $2
mri_setfield -field images.blf -all012 -value ' %f,%f,%f' $2
mri_setfield -field images.brf -all012 -value ' %f,%f,%f' $2
mri_setfield -field images.blb -all012 -value ' %f,%f,%f' $2
mri_setfield -field images.brb -all012 -value ' %f,%f,%f' $2
mri_setfield -field images.ctr -all012 -value ' %f,%f,%f' $2
mri_resample -d x -len %d -start %f -end %f $2 tmp_coreg_stoi_1
mri_resample -d y -len %d -start %f -end %f tmp_coreg_stoi_1 tmp_coreg_stoi_2
mri_resample -d z -len %d -start %f -end %f tmp_coreg_stoi_2 $3
mri_setfield -field images.tlf -all012 -value ' %f,%f,%f' $3
mri_setfield -field images.trf -all012 -value ' %f,%f,%f' $3
mri_setfield -field images.tlb -all012 -value ' %f,%f,%f' $3
mri_setfield -field images.trb -all012 -value ' %f,%f,%f' $3
mri_setfield -field images.blf -all012 -value ' %f,%f,%f' $3
mri_setfield -field images.brf -all012 -value ' %f,%f,%f' $3
mri_setfield -field images.blb -all012 -value ' %f,%f,%f' $3
mri_setfield -field images.brb -all012 -value ' %f,%f,%f' $3
mri_setfield -field images.ctr -all012 -value ' %f,%f,%f' $3
mri_setfield -field images.voxel_spacing -allxyz -value ' %f,%f,%f' $3
mri_setfield -field images.voxel_size -allxyz -value ' %f,%f,%f' $3
mri_destroy_dataset tmp_coreg_stoi_1
mri_destroy_dataset tmp_coreg_stoi_2
rm tmp_coreg_strct_to_inplane.par
"""

rot1ParHeader="""\
##Format: order:index_t, type:raw
##Format: names:(3d_qbar_x,3d_qbar_y,3d_qbar_z,3d_qbar_w,3d_deltabarx,3d_deltabary,3d_deltabarz,mse)
# Estimate for rot to bring structurals to inplane alignment
"""

def pickMagicPadSize( sz ):
    for magic in magicSizes:
        if magic>=sz:
            return magic-sz;

def calcBBoxAxes( bbox ):
    xvec= bbox.brf - bbox.blf
    xvec.normalize()
    yvec= bbox.blb - bbox.blf
    yvec.normalize()
    zvec= bbox.tlf - bbox.blf
    zvec.normalize()
    return (xvec, yvec, zvec)

def checkDots(lbl,qTrans,vec,basisX,basisY,basisZ):
    dots= []
    for b in (basisX, basisY, basisZ):
        vDup= vec.clone()
        trans_vec_mult(qTrans,vDup)
        dots.append(vDup.dot(b))
    debugMessage("%s: %f %f %f "%(lbl,dots[0],dots[1],dots[2]))
    
def makeAligningRotation( fromBBox, toBBox ):
    fromX, fromY, fromZ = calcBBoxAxes(fromBBox)
    toX, toY, toZ = calcBBoxAxes(toBBox)
    debugMessage( "Generating aligning transformation:")
    debugMessage( "  fromX: %s"%fromX )
    debugMessage( "  fromY: %s"%fromY )
    debugMessage( "  fromZ: %s"%fromZ )
    debugMessage( "  toX: %s"%toX )
    debugMessage( "  toY: %s"%toY )
    debugMessage( "  toZ: %s"%toZ )
    axis1= fromZ.cross( toZ )
    cosTheta= fromZ.dot( toZ )
    sinTheta= axis1.mag()
    if axis1.mag() != 0.0:
        axis1.normalize()
    theta= atan2(sinTheta,cosTheta)
    debugMessage( "  theta= %f radians = %f degrees"%(theta,(180.0/pi)*theta) )
    debugMessage( "  axis is %s"%axis1 )
    Q1= quat_from_axis_angle(None,axis1[0],axis1[1],axis1[2],theta)
    Q1Trans= Transform() 
    quat_to_trans(Q1Trans,Q1,0.0,0.0,0.0);
    rotFromX= fromX.clone()
    trans_vec_mult(Q1Trans,rotFromX)
    debugMessage("  rotFromX: %s"%rotFromX)
    axis2= rotFromX.cross(toX) # direction same as toZ
    debugMessage("  axis2: %s"%axis2)
    if axis2.mag() == 0.0:
        phi= 0.0
        Q2= quat_identity(None)
    else:
        cosPhi= rotFromX.dot(toX)
        sinPhi= axis2.mag()
        axis2.normalize()
        phi= atan2(sinPhi,cosPhi)
        Q2= quat_from_axis_angle(None,axis2[0],axis2[1],axis2[2],phi)
    Q= quat_copy(None,Q1)
    quat_mult_right(Q,Q2)
    QTrans= Transform()
    quat_to_trans(QTrans,Q,0.0,0.0,0.0)
    debugMessage( "  phi= %f radians = %f degrees"%(phi,(180.0/pi)*phi) )
    if getDebug():
        identityTrans= Transform()
        checkDots("      raw X",identityTrans,fromX,toX,toY,toZ)
        checkDots("  aligned X",QTrans,fromX,toX,toY,toZ)
        checkDots("      raw Y",identityTrans,fromY,toX,toY,toZ)
        checkDots("  aligned Y",QTrans,fromY,toX,toY,toZ)
        checkDots("      raw Z",identityTrans,fromZ,toX,toY,toZ)
        checkDots("  aligned Z",QTrans,fromZ,toX,toY,toZ)
    return Q

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

oscript= "apply_coreg_struct_to_inplane.gen_csh"
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

inplaneBBox= bboxFromFilename(os.path.abspath(pargs[0]))
inplaneChunk= inplaneBBox.chunk
inplaneVox= inplaneBBox.vox
strctBBox= bboxFromFilename(os.path.abspath(pargs[1]))
strctChunk= strctBBox.chunk
strctVox= strctBBox.vox

#Check reasonableness of input
for thisChunk in (inplaneBBox.chunk, strctBBox.chunk):
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

if clAlgString != None:
    algString= clAlgString
elif inplaneChunk.isT1Weighted() == strctChunk.isT1Weighted():
    algString= "opt=praxis,optol=0.01,obj=mse,meanc"
else:
    algString= "opt=praxis,optol=0.01,obj=jointentropy"

# Move to a temporary directory 
tmpdir= makeTempDir('tmp_coreg_stoi')
homedir= os.getcwd()
os.chdir(tmpdir)

(strctX,strctY,strctZ)= calcBBoxAxes(strctBBox)
scannerToStrctTrans= Transform([strctX[0],strctX[1],strctX[2],0.0,\
                                strctY[0],strctY[1],strctY[2],0.0,\
                                strctZ[0],strctZ[1],strctZ[2],0.0,\
                                0.0,0.0,0.0,1.0])

inplaneBBox.printBounds("Inplane bounding box")
strctBBox.printBounds("Structural bounding box")
inplaneCtr= inplaneBBox.ctr
strctCtr= strctBBox.ctr
ctrShiftScannerCoords= inplaneCtr - strctCtr
Message( "Center shift in scanner coordinates: %s"%ctrShiftScannerCoords )

ctrShift= ctrShiftScannerCoords.clone()
trans_vec_mult(scannerToStrctTrans,ctrShift)
Message( "Center shift in structural coordinates: %s"%ctrShift )
Message( "Separation between centers is %f mm"%ctrShift.mag() )

Q1= makeAligningRotation(strctBBox, inplaneBBox)
Q1_bar= quat_copy(None,Q1)
quat_conjugate(Q1_bar)
scaledCtrShift= Vec4( ctrShift[0]/strctVox[0], \
                      ctrShift[1]/strctVox[1], \
                      ctrShift[2]/strctVox[2] )

Message("Geometrical transformation to bring structurals into inplane space:")
Message("Quaternion 1: %s"%Q1)
Message("Translation 1 (unscaled): %s"%ctrShift)
Message("Translation 1 (scaled):   %s"%scaledCtrShift)

safeRun("mri_copy_dataset %s strctDup"%strctBBox.chunk.ds.fname)
safeRun("mri_copy_dataset %s inplaneDup"%inplaneBBox.chunk.ds.fname)
safeRun("mri_remap -order xyzt strctDup")
f= open("rot1.par","w")
f.write(rot1ParHeader)
f.write("0 %f %f %f %f %f %f %f 1.0\n"%\
        (Q1.x,Q1.y,Q1.z,Q1.w,\
         scaledCtrShift[0],scaledCtrShift[1],scaledCtrShift[2]))
f.close()
safeRun("ireg3d -x %f -y %f -z %f -i strctDup -h rot1 -p rot1.par -real"%\
        (strctVox[0],strctVox[1],strctVox[2]))

rot1BBox= bboxFromFilename("rot1")
rot1BBox.setCtr( inplaneCtr )
rot1BBox.setVox(strctVox)
rot1BBox.setCorners( Vec4(1.0,0.0,0.0), Vec4(0.0,1.0,0.0), Vec4(0.0,0.0,1.0) )
rot1BBox.printBounds("Bounds for rotated structurals,inplane coordinates")
rot1BBox.exportBounds()
rot1Chunk= rot1BBox.chunk

inplaneDupBBox= bboxFromFilename("inplaneDup")
inplaneDupBBox.setCtr(inplaneCtr)
inplaneDupBBox.setVox(inplaneVox)
inplaneDupBBox.setCorners( Vec4(1.0,0.0,0.0), Vec4(0.0,1.0,0.0), \
                           Vec4(0.0,0.0,1.0))
inplaneDupBBox.printBounds("Bounds for inplanes, inplane coordinates")
inplaneDupBBox.exportBounds()

rot1Tlf= rot1BBox.tlf
rot1Blf= rot1BBox.blf
rot1Trf= rot1BBox.trf
rot1Brf= rot1BBox.brf
rot1Tlb= rot1BBox.tlb
rot1Blb= rot1BBox.blb
rot1Trb= rot1BBox.trb
rot1Brb= rot1BBox.brb
rot1Ctr= copy.copy(inplaneCtr)
rot1Vox= rot1BBox.vox

#
# rot1 gives us the transformation to bring the structurals to the
# inplane coordinates, as calculated directly from the scanner data.
# We now estimate an additional small rotation and shift to bring
# the rot1 data to the structural data, to account for head motion
# between the two scans.
#
# We must pick a target size in Z which is good for FFTs and larger 
# than the number of actual slices in the inplane dataset.
#

zPad= pickMagicPadSize(inplaneBBox.zdim+1)+1

#
# Start by downsampling the inplanes in x and y, and downsampling the
# rotated structurals in z, so that everyone has the same (maximum)
# voxel size and the corresponding bounding box.  Remember that rot1BBox
# has its corner coords in the in-plane coordinate system.
#

xStart= (rot1Tlf[0]-rot1Ctr[0])/inplaneVox[0] + (inplaneChunk.getDim('x')/2);
xEnd= (rot1Trf[0]-rot1Ctr[0])/inplaneVox[0] + (inplaneChunk.getDim('x')/2);
safeRun("mri_resample -d x -len %d -start %f -end %f %s inplResampX"%\
        (rot1BBox.chunk.getDim('x'), xStart, xEnd, inplaneChunk.ds.fname))

yStart= (rot1Tlf[1]-rot1Ctr[1])/inplaneVox[1] + (inplaneChunk.getDim('y')/2);
yEnd= (rot1Tlb[1]-rot1Ctr[1])/inplaneVox[1] + (inplaneChunk.getDim('y')/2);
safeRun("mri_resample -d y -len %d -start %f -end %f inplResampX inplResampXY"%\
        (rot1BBox.chunk.getDim('y'), yStart, yEnd))

safeRun("mri_pad -d z -len %d -shift %d inplResampXY inplResampXYZ"%\
        (inplaneBBox.zdim+zPad, (zPad+1)/2))

inplResampVox= [ strctVox[0], strctVox[1], inplaneVox[2] ]

inplResampXYZBBox= bboxFromFilename('inplResampXYZ')
inplResampXYZBBox.setCtr( inplaneCtr )
inplResampXYZBBox.setVox( inplResampVox )
inplResampXYZBBox.setCorners( Vec4(1.0,0.0,0.0), Vec4(0.0,1.0,0.0), \
                             Vec4(0.0,0.0,1.0) )
inplResampXYZBBox.printBounds("Bounds for resampled inplanes,inplane coordinates")
inplResampXYZBBox.exportBounds()

inplaneTlf= inplResampXYZBBox.tlf
inplaneBlf= inplResampXYZBBox.blf
inplaneTrf= inplResampXYZBBox.trf
inplaneBrf= inplResampXYZBBox.brf
inplaneTlb= inplResampXYZBBox.tlb
inplaneBlb= inplResampXYZBBox.blb
inplaneTrb= inplResampXYZBBox.trb
inplaneBrb= inplResampXYZBBox.brb

zStart= ((inplaneBlf[2]-inplaneCtr[2])/rot1Vox[2]) + (rot1Chunk.getDim('z')/2);
zEnd= ((inplaneTlf[2]-inplaneCtr[2])/rot1Vox[2]) + (rot1Chunk.getDim('z')/2);
safeRun("mri_resample -d z -len %d -start %f -end %f %s rot1ResampZ"%\
        (inplResampXYZBBox.chunk.getDim('z'), zStart, zEnd, \
         rot1Chunk.ds.fname))

rot1ResampZBBox= bboxFromFilename('rot1ResampZ')
rot1ResampZBBox.setCtr( inplaneCtr )
rot1ResampZBBox.setVox( inplResampVox )
rot1ResampZBBox.setCorners( Vec4(1.0,0.0,0.0), Vec4(0.0,1.0,0.0), \
                            Vec4(0.0,0.0,1.0) )
rot1ResampZBBox.printBounds("Bounds for resampled structurals,inplane coordinates")
rot1ResampZBBox.exportBounds()

moveThisBBox= rot1ResampZBBox
alignToBBox= inplResampXYZBBox

if algRequiresGradMag(algString):
    verboseMessage("Producing grad magnitudes...")
    moveThisBBox= moveThisBBox.chunk.makeGradMag3D(moveThisBBox.vox,\
                                                   'strctGrad')
    alignToBBox= alignToBBox.chunk.makeGradMag3D(alignToBBox.vox,'inplGrad')
    
if algRequiresMeanCorrect(algString):
    verboseMessage("Mean correcting...")
    moveThisMean= moveThisBBox.chunk.getMaskedMean(inplResampXYZBBox.chunk,40)
    alignToMean= alignToBBox.chunk.getMaskedMean(inplResampXYZBBox.chunk,40)
    safeRun("mri_rpn_math -out alignto '$1,%f,*' %s"%\
            ((moveThisMean/alignToMean),alignToBBox.chunk.ds.fname))
    alignToBBox= bboxFromFilename("alignto")

flags= "-x %f -y %f -z %f"%tuple(inplResampVox)
if getDebug():
    flags= flags+" -debug"
if algString != None:
    flags= flags+" -alg %s"%simplifyAlgString(algString)
verboseMessage("Beginning 3D alignment")
safeRun("parallel.run.csh estireg3d %s -input %s -align %s -p rot2.par"%\
        (flags,moveThisBBox.chunk.ds.fname,alignToBBox.chunk.ds.fname))
verboseMessage("Finished 3D alignment")

vals= None
try:
    ifile= open("rot2.par","r")
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
        (inplResampVox[0], inplResampVox[1], inplResampVox[2]) \
        + "-weight %s rot2.par"%alignToBBox.chunk.ds.fname);
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

Q2= Quat(vals[1],vals[2],vals[3],vals[4])
Q2_bar= quat_copy(None,Q2)
quat_conjugate(Q2_bar)
rawShift2= Vec4(vals[5],vals[6],vals[7])
unscaledShift2= Vec4(rawShift2[0]*inplResampVox[0],\
                     rawShift2[1]*inplResampVox[1],\
                     rawShift2[2]*inplResampVox[2]);
rescaledShift2= Vec4(unscaledShift2[0]/strctVox[0],
                     unscaledShift2[1]/strctVox[1],
                     unscaledShift2[2]/strctVox[2]);

Message("Geometrical transformation for structural/inplane relative motion:")
Message("Quaternion 2: %s"%Q2)
Message("Translation 2 (unscaled): %s"%unscaledShift2)
Message("Translation 2 (scaled):   %s"%rescaledShift2)
Message("Rot, trans, and mean disp from this transform: %f %f %f\n"%\
        (dispVals[1],dispVals[2],dispVals[3]));

# Merge first and second transformations. Because of the Fiasco
# convention that the parameter file for aligning A to B gives
# the transformation from B to A, we are calculating:
#   Tnet Tnet = ( (T2 R2)' (T1 R1)' )'
#             = T1 R1 T2 R2
#             = T1 ( R1 T2 R1' ) R1 R2
# where ' means inverse.  We'll use quaternion algebra.
#
t2Quat= quat_identity(None)
t2Quat.x= unscaledShift2[0]
t2Quat.y= unscaledShift2[1]
t2Quat.z= unscaledShift2[2]
t2Quat.w= 0.0
quat_mult_right(t2Quat,Q1_bar)
quat_mult_left(Q1,t2Quat)
rotatedT2= Vec4(t2Quat.x,t2Quat.y,t2Quat.z)
netTrans= ctrShift + rotatedT2
Message("unscaled net translation: %s"%netTrans)
rescaledNetTrans= Vec4(netTrans[0]/strctVox[0], \
                     netTrans[1]/strctVox[1], \
                     netTrans[2]/strctVox[2])
netQ= quat_copy(None,Q1)
quat_mult_right(netQ,Q2)
netQ_bar= quat_copy(None,netQ)
quat_conjugate(netQ_bar)
Message("Net transformation to align structural to inplane:")
Message("Net Quaternion: %s"%netQ)
Message("Net Translation (unscaled): %s"%netTrans)
Message("Net Translation (scaled):   %s"%rescaledNetTrans)

# Generate a parameter file for the net rotation
f= open("netRot.par","w")
f.write(parFileFormat%( strctVox[0], strctVox[1], strctVox[2], \
                        netQ.x, netQ.y, netQ.z, netQ.w, \
                        rescaledNetTrans[0], rescaledNetTrans[1], \
                        rescaledNetTrans[2], 1.0 ))
f.close()

safeRun("ireg3d -x %f -y %f -z %f -i strctDup -h alignedStruct -p netRot.par -real"%\
        (strctVox[0],strctVox[1],strctVox[2]))
alignedStrctBBox= bboxFromFilename('alignedStruct')
alignedStrctBBox.setCtr( inplaneCtr )
alignedStrctBBox.setCorners( Vec4(1.0,0.0,0.0), Vec4(0.0,1.0,0.0), \
                            Vec4(0.0,0.0,1.0) )
alignedStrctBBox.printBounds("Bounds for aligned structurals,inplane coordinates")
alignedStrctBBox.exportBounds()

# All done- now we write a shell script that implements this transformation.
# This block consists mainly of filling out the blanks in the generated
# script.
scriptName= os.path.join(homedir,oscript);
ofile= open(scriptName,"w")
ofile.write(generatedScriptHdr)
ifile= open("netRot.par","r")
parLines= ifile.readlines()
ifile.close()
for line in parLines:
    ofile.write(line)
argList= strctVox[0:3]
for corner in ('tlf', 'trf', 'tlb', 'trb', 'blf', 'brf', 'blb', 'brb', 'ctr'):
    vtx= alignedStrctBBox.__dict__[corner]
    argList= argList+[ vtx[0], vtx[1], vtx[2] ]
    
xStart= (inplaneDupBBox.tlf[0]-inplaneDupBBox.ctr[0])/strctVox[0] \
        + (alignedStrctBBox.chunk.getDim('x')/2);
xEnd= (inplaneDupBBox.trf[0]-inplaneDupBBox.ctr[0])/strctVox[0] \
        + (alignedStrctBBox.chunk.getDim('x')/2);
argList= argList+ [ inplaneDupBBox.chunk.getDim('x'), xStart, xEnd ]

yStart= (inplaneDupBBox.tlf[1]-inplaneDupBBox.ctr[1])/strctVox[1] \
        + (alignedStrctBBox.chunk.getDim('y')/2);
yEnd= (inplaneDupBBox.tlb[1]-inplaneDupBBox.ctr[1])/strctVox[1] \
        + (alignedStrctBBox.chunk.getDim('y')/2);
argList= argList+ [ inplaneDupBBox.chunk.getDim('y'), yStart, yEnd ]

zStart= (inplaneDupBBox.blf[2]-inplaneDupBBox.ctr[2])/strctVox[2] \
        + (alignedStrctBBox.chunk.getDim('z')/2);
zEnd= (inplaneDupBBox.tlf[2]-inplaneDupBBox.ctr[2])/strctVox[2] \
        + (alignedStrctBBox.chunk.getDim('z')/2);
argList= argList+ [ inplaneDupBBox.chunk.getDim('z'), zStart, zEnd ]

for corner in ('tlf', 'trf', 'tlb', 'trb', 'blf', 'brf', 'blb', 'brb', 'ctr'):
    vtx= inplaneDupBBox.__dict__[corner]
    argList= argList+[ vtx[0], vtx[1], vtx[2] ]

argList= argList+inplaneVox[0:3]    
argList= argList+inplaneVox[0:3]
ofile.write(generatedScriptBody%tuple(argList))
ofile.close()
safeRun("chmod ug+x %s"%scriptName)

#safeRun("cp *.par %s"%homedir)
#safeRun("mri_copy_dataset rot1 %s/rot1"%homedir)
#safeRun("mri_copy_dataset alignto %s/alignto"%homedir)

# Clean up
os.chdir(homedir)
if not getDebug():
    removeTmpDir(tmpdir)


