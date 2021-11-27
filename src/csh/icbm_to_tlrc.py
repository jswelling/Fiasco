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
if "FIASCO" in os.environ:
    sys.path.append(os.environ["FIASCO"])
from fiasco_utils import *

idString= "$Id: icbm_to_tlrc.py,v 1.3 2006/05/04 23:02:58 welling Exp $"

###########
# This program uses a 2-affine-transformation warp to reshape ICBM data
# in Talairach space to match the Talairach atlas.  The particular method
# is from http://www.mrc-cbu.cam.ac.uk/Imaging/Common/mnispace.shtml
# (Matthew Brett 5/8/99 and 14/2/02)
###########
# Transform to warp ICBM data in Talairach space to match Talairach atlas,
# for points at and above the anterior commisure
aboveACTransform= Transform([0.9900, 0.0, 0.0, 0.0, \
                             0.0, 0.9688, 0.0460, 0.0,\
                             0.0, -0.0485, 0.9189, 0.0,\
                             0.0, 0.0, 0.0, 1.0])
# Transform to warp ICBM data in Talairach space to match Talairach atlas,
# for points below the anterior commisure
belowACTransform= Transform([0.9900, 0.0, 0.0, 0.0,\
                             0.0, 0.9688, 0.0420, 0.0,\
                             0.0, -0.0485, 0.8390, 0.0,\
                             0.0, 0.0, 0.0, 1.0])

# The transform P maps Pgh MRI 3D coordinates (right handed with origin
# at the Fourier center of rotation) to Talairach coordinates.  P happens
# to be its own inverse, not for any deep reason.
P= Transform([1.0, 0.0, 0.0,  0.0,\
              0.0,-1.0, 0.0, 15.0,\
              0.0, 0.0,-1.0, 10.0,\
              0.0, 0.0, 0.0,  1.0])
PInverse= P

# Standard dimensions of ICBM and Talairach grids
icbmGridSize= [181,217,181]
tlrcGridSize= [161,191,151]
dimNames= ["x","y","z"]

def createIwarpParFile( fname, trans, tdim ):
    f= file(fname,"w")
    for t in xrange(tdim):
        f.write("%d "%t)
        for i in xrange(16):
            f.write("%f "%trans[i])
        f.write("1.0\n")
    f.close()
    

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

outDSName= None
interpMethod= "linear"

try:
    (opts,pargs) = getopt.getopt(sys.argv[1:],"vd",["closest","trilinear"])
except:
    print("%s: Invalid command line parameter" % sys.argv[0])
    describeSelf();
    sys.exit()

for a,b in opts:
    if a=="-v":
        setVerbose(1)
    if a=="-d":
        setDebug(1)
    if a=="--trilinear":
        interpMethod= "linear"
    if a=="--closest":
        interpMethod= "closest"
    
    

if len(pargs) < 1:
    sys.exit("Required input dataset name not specified")
if len(pargs) < 2:
    sys.exit("Required output dataset name not specified")
elif len(pargs) > 2:
    sys.exit("Too many dataset names specified")
inDSName= pargs[0]
outDSName= pargs[1]

# Is the input dataset an ICBM volume?
ds= MRIDataset(inDSName)
if not ds.hasChunk("images"):
    sys.exit("%s has no images chunk!"%inDSName)
images= ds.getChunk("images")
if not images.hasValue("dimensions"):
    sys.exit("%s images chunk has no dimensions!"%inDSName)
if images.getValue("dimensions")[0:3] != "xyz":
    sys.exit("%s images chunk has unexpected dimensions %s!"%\
             (inDSName,images.getValue("dimensions")))
for i in xrange(3):
    if images.getDim(dimNames[i]) != icbmGridSize[i]:
        sys.exit("%s images chunk dimensions are not those of an ICBM dataset!"%inDSName)
vox= getVox(images)
if vox[0]!=1.0 or vox[1]!=1.0 or vox[2]!=1.0:
    sys.exit("%s images chunk voxel size is not that of an ICBM dataset!"%inDSName)

# Make a temp directory
tmpdir= makeTempDir('tmp_icbm_to_tlrc')
homedir= os.getcwd()

#
# The actual affine transformations to be applied above and below the AC
# are P*A*Pinv, where A is the appropriate ICBM-to-Talairach affine trans
# in Talairach coordinates.  This produces transforms we can use on the
# data in Pgh MRI 3D coordinates, which is what 'iwarp' expects.
#
fullAboveACTrans= Transform();
trans_mult_right(fullAboveACTrans,P)
trans_mult_right(fullAboveACTrans,aboveACTransform)
trans_mult_right(fullAboveACTrans,PInverse)
debugMessage("fullAboveACTrans follows:")
debugMessage("%s"%fullAboveACTrans)

fullBelowACTrans= Transform();
trans_mult_right(fullBelowACTrans,P)
trans_mult_right(fullBelowACTrans,belowACTransform)
trans_mult_right(fullBelowACTrans,PInverse)
debugMessage("fullBelowACTrans follows:")
debugMessage("%s"%fullBelowACTrans)

# Generate iwarp parameter files
if images.hasValue("extent.t"):
    tdim= images.getDim("t")
else:
    tdim= 1

# Run iwarp, once for each half
upperParName= os.path.join(tmpdir,"upper.par")
lowerParName= os.path.join(tmpdir,"lower.par")
upperMRIName= os.path.join(tmpdir,"upper")
lowerMRIName= os.path.join(tmpdir,"lower")
createIwarpParFile(upperParName, fullAboveACTrans, tdim)
createIwarpParFile(lowerParName, fullBelowACTrans, tdim)
safeRun("iwarp -interp %s -xvoxel %f -yvoxel %f -zvoxel %f -i %s -p %s -h %s"%\
        (interpMethod,vox[0],vox[1],vox[2],inDSName,upperParName,upperMRIName))
safeRun("iwarp -interp %s -xvoxel %f -yvoxel %f -zvoxel %f -i %s -p %s -h %s"%\
        (interpMethod,vox[0],vox[1],vox[2],inDSName,lowerParName,lowerMRIName))

# We now need the location of the AC in Pgh MRI 3D coordinates.  Rather
# than just writing it down, let's calculate it, so as to avoid an
# error if we change the coordinate transformations in the future.
# We also need to remember that the Pgh MRI 3D coordinate system
# has its origin at the Fourier center, but Pgh MRI grid coordinates
# have their Z origin at the superior edge.
acCoords= Vec4(0.0,0.0,0.0)
trans_vec_mult(PInverse,acCoords)
debugMessage("AC coordinates in Pgh MRI 3D space: %s"%acCoords)
cutoff= acCoords[2]+(images.getDim("z")/2)
debugMessage("Breakpoint between the two affine transforms is at Z=%f"%cutoff)

# The 'upper' data goes on the superior side of the AC, but in Pgh MRI
# coordinates Z values increase in the inferior direction.  Thus we really
# want to use the upper values when z is *less* than or equal to cutoff.
mergedMRIName= os.path.join(tmpdir,"merged")
safeRun("mri_rpn_math -out %s '$1,$2,$z,%f,<=,if_keep' %s %s"%\
        (mergedMRIName,cutoff,upperMRIName,lowerMRIName))

# Clip down to standard Talairach dimensions
clipXName= os.path.join(tmpdir,"clip_x")
safeRun("mri_subset -d x -len %d -s %d %s %s"%\
        (tlrcGridSize[0],int((icbmGridSize[0]-tlrcGridSize[0])/2),
         mergedMRIName,clipXName))
clipXYName= os.path.join(tmpdir,"clip_xy")
safeRun("mri_subset -d y -len %d -s %d %s %s"%\
        (tlrcGridSize[1],int((icbmGridSize[1]-tlrcGridSize[1])/2),
         clipXName,clipXYName))
safeRun("mri_subset -d z -len %d -s %d %s %s"%\
        (tlrcGridSize[2],int((icbmGridSize[2]-tlrcGridSize[2])/2),
        clipXYName,outDSName))

# Clean up
os.chdir(homedir)
if not getDebug():
    removeTmpDir(tmpdir)
