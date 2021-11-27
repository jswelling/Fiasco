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

idString= "$Id: coregister_makeps.py,v 1.5 2006/05/04 23:02:58 welling Exp $"

cmapString= """
0 0 0 0 %s
0 255 255 1 %s
0 128 255 1 %s
0 0 255 1 %s
0 0 128 1 %s
0 0 0 1 %s
"""

def getMax(chunk):
    totDim= chunk.getDim('x')*chunk.getDim('y')*chunk.getDim('z')
    safeRun("mri_copy_dataset %s tmp_mean_1"%chunk.ds.fname)
    safeRun("mri_remap -order q -length %d tmp_mean_1"%totDim)
    safeRun("mri_subsample -d q -length 1 -max tmp_mean_1 tmp_mean_2")
    lines= readCmdOutputToList("mri_rpn_math '$1,1,if_print_1' tmp_mean_2")
    return float(lines[0])

def checkInputStructure( chunk ):
    dimstr= chunk.getValue('dimensions');
    if dimstr != "xyz":
        if dimstr == "xyzt":
            if chunk.getDim("t") != 1:
                sys.exit("Input file %s must have t extent 1!"%\
                         os.path.basename(chunk.ds.fname))
        elif dimstr == "vxyzt":
            if chunk.getDim("t") != 1:
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

def pickSliceNums( nSlices ):
    return ( int(nSlices/3), int(nSlices/2), int((2*nSlices)/3) )

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
    (opts,pargs) = getopt.getopt(sys.argv[1:],"vd",["out=","func=","strct="])
except:
    sys.stderr.write("%s: Invalid command line parameter\n" % sys.argv[0])
    describeSelf();
    sys.exit()

#Check calling syntax; parse args
if len(pargs) != 0 :
    describeSelf()
    sys.exit(1)

outFile= None
funcFile= None
strctFile= None
for a,b in opts:
    if a=="-v":
        setVerbose(1)
    if a=="-d":
        setDebug(1)
    if a=="--out":
        outFile= os.path.abspath(b)
    if a=="--func":
        funcFile= b
    if a=="--strct":
        strctFile= b

if outFile==None: sys.exit("Required output postscript file name not given")
if funcFile==None: sys.exit("Required input functional file name not given")
if strctFile==None: sys.exit("Required input structural file name not given")

funcDS= MRIDataset(os.path.abspath(funcFile))
funcChunk= funcDS.getChunk('images')
strctDS= MRIDataset(os.path.abspath(strctFile))
strctChunk= strctDS.getChunk('images')
checkInputStructure(funcChunk)
checkInputStructure(strctChunk)
funcBBox= BBox(funcChunk)
strctBBox= BBox(strctChunk)

# Move to a temporary directory 
tmpdir= makeTempDir('tmp_coregister_makeps')
homedir= os.getcwd()
os.chdir(tmpdir)

# Get relevant dimensions
funcXdim= funcChunk.getDim("x")
funcYdim= funcChunk.getDim("y")
funcZdim= funcChunk.getDim("z")
strctXdim= strctChunk.getDim("x")
strctYdim= strctChunk.getDim("y")
strctZdim= strctChunk.getDim("z")

scaleX= (float(strctXdim-1)*strctBBox.vox[0])/\
        (float(funcXdim-1)*funcBBox.vox[0])
scaleY= (float(strctYdim-1)*strctBBox.vox[1])/\
        (float(funcYdim-1)*funcBBox.vox[1])
scaleZ= (float(strctZdim-1)*strctBBox.vox[2])/\
        (float(funcZdim-1)*funcBBox.vox[2])
xLow= 0.5*(1.0-scaleX)*float(funcXdim-1)
xHigh= float(funcXdim-1)-xLow
yLow= 0.5*(1.0-scaleY)*float(funcYdim-1)
yHigh= float(funcYdim-1)-yLow
zLow= 0.5*(1.0-scaleZ)*float(funcZdim-1)
zHigh= float(funcZdim-1)-zLow

safeRun("mri_resample -d x -len %d -start %f -end %f %s stretch_x"%\
        (strctXdim,xLow,xHigh,funcDS.fname))
safeRun("mri_resample -d y -len %d -start %f -end %f stretch_x stretch_xy"%\
        (strctYdim,yLow,yHigh))
if funcDS.hasChunk("missing"):
    safeRun("mri_delete_chunk -chunk missing stretch_xy")
safeRun("mri_resample -d z -len %d -start %f -end %f stretch_xy stretch_xyz"%\
        (strctZdim,zLow,zHigh))

funcMax= getMax(funcChunk)

writeStringToCmdInput(cmapString%(0.15*funcMax, 0.15*funcMax,
                                0.5*funcMax,0.5*funcMax,
                                funcMax,funcMax),
                    "mri_from_ascii -c color_table -ord vc -length 5:6 cmap")
safeRun("colorize -col cmap stretch_xyz func_over")
safeRun("mri_remap -order vxyz func_over")

# Mask  the structural before overlaying to improve contrast
safeRun("mri_permute -order xyzv func_over func_over_p")
safeRun("mri_subset -d v -len 1 -s 3 func_over_p mask")
safeRun("mri_remap -order vxyz mask")
safeRun("mri_rpn_math -out strct_masked '$1,0.0,$2,if_keep' %s mask"%\
        strctDS.fname)
safeRun("matte -inmap func_over strct_masked overlaid")
safeRun("mri_subsample -d z -len 9 -closest overlaid overlaid_z")
safeRun("mri_to_img -ps -reverse -mosaic overlaid_z %s"%outFile)

xFOV= strctXdim*strctChunk.getFloat("voxel_spacing.x")
yFOV= strctYdim*strctChunk.getFloat("voxel_spacing.y")
zFOV= strctZdim*strctChunk.getFloat("voxel_spacing.z")
xCoverage= int(strctXdim*(zFOV/xFOV))
yCoverage= int(strctYdim*(zFOV/yFOV))

sliceImageList= []
(s1,s2,s3)= pickSliceNums( strctZdim )
count= 1
for slice in [s1,s2,s3]:
    safeRun("mri_subset -d z -len 1 -s %d overlaid z%d"%(slice,count))
    safeRun("mri_remap -order vxyz -len 4:%d:%d:1 z%d"%\
            (strctXdim,strctYdim,count))
    sliceImageList.append("z%d"%count)
    count= count+1

(s1,s2,s3)= pickSliceNums( strctXdim )
count= 1
for slice in [s1,s2,s3]:
    safeRun("mri_subset -d x -len 1 -s %d overlaid x%d"%(slice,count))
    safeRun("mri_remap -order vxyz -len 4:%d:%d:1 x%d"%\
            (strctYdim,strctZdim,count))
    safeRun("mri_interp -d x -len %d -linear x%d x%d_fix"%\
            (strctXdim,count,count))
    safeRun("mri_interp -d y -len %d -linear x%d_fix x%d_fix2"%\
            (yCoverage,count,count))
    safeRun("mri_pad -d y -len %d -shift %d x%d_fix2 x%d_fix3"%
            (strctYdim,int((strctYdim-yCoverage)/2),count,count))
    sliceImageList.append("x%d_fix3"%count)
    count= count+1

(s1,s2,s3)= pickSliceNums( strctYdim )
count= 1
for slice in [s1,s2,s3]:
    safeRun("mri_subset -d y -len 1 -s %d overlaid y%d"%(slice,count))
    safeRun("mri_remap -order vxyz -len 4:%d:%d:1 y%d"%\
            (strctXdim,strctZdim,count))
    safeRun("mri_interp -d y -len %d -linear y%d y%d_fix"%\
            (yCoverage,count,count))
    safeRun("mri_pad -d y -len %d -shift %d y%d_fix y%d_fix2"%
            (strctYdim,int((strctYdim-yCoverage)/2),count,count))
    sliceImageList.append("y%d_fix2"%count)
    count= count+1

argList= "mri_paste -d z -out images"
for img in sliceImageList:
    argList= argList + " " + img
safeRun(argList)

argList= "-ps -reverse -mosaic -title -date_created -date_printed"
argList= argList + " -black_min 0 -black_max 255"
if not os.environ.has_key('F_PS_TITLE1'):
    if funcChunk.hasValue("scan.id"):
        os.environ['F_PS_TITLE1']= "ID %s"%funcChunk.getValue("scan.id")
    else:
        os.environ['F_PS_TITLE1']= os.path.basename(funcChunk.ds.fname)
os.environ['F_PS_TITLE2']= "Coregistration Test"
safeRun("mri_to_img %s images %s"%(argList,outFile))

# Clean up
os.chdir(homedir)
if not getDebug():
    removeTmpDir(tmpdir)


