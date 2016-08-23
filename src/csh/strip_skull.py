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

idString= "$Id: strip_skull.py,v 1.15 2006/05/04 23:02:58 welling Exp $"

afniSkullStripCmd= "3dIntracranial"
afni3dSkullStripCmd= "3dSkullStrip"

def outputAlreadyExists( inChunk, outFname ):
    fullInFname= inChunk.ds.fname+'.mri'
    fullOutFname= outFname+'.mri'
    if not os.access(fullOutFname,os.F_OK):
        return 0;
    else:
        return isNewer( fullOutFname, fullInFname )
    
def afniSkullStripApplicable( infoDict, inChunk, outFname ):
    return 1 # always usable

def afni3dSkullStripStrip( infoDict, inChunk, outFname ):
    checkExeFound( afni3dSkullStripCmd, 'AFNI' )
    safeRun("mri_type_convert -chunk %s -short %s struct"%\
            (inChunk.name,inChunk.ds.fname))
    if inChunk.hasValue('afni.SCENE_DATA.0'):
        afniSpace= inChunk.getValue( 'afni.SCENE_DATA.0' )
    else:
        afniSpace= '+orig'
    safeRun("pghtoafni -anat struct struct%s"%afniSpace)
    argString= "-input struct%s -prefix strippedStruct"%afniSpace
    outLines= readCmdOutputToList("%s %s"%(afni3dSkullStripCmd,argString))
    for line in outLines:
        line= string.strip(line)
    safeRun("smartreader -i strippedStruct%s.HEAD -out tmp_stripped"%\
            afniSpace)
    safeRun("mri_rpn_math -out %s '$2' struct tmp_stripped"%outFname)
    safeRun("mri_setfield -field images.skull_strip_info -value \"%s\" %s"%\
            ("afni3dSkullStrip()",outFname))
    safeRun("mri_history -add \"%s\" %s"%(string.join(sys.argv),outFname))

def afniSkullStrip( infoDict, inChunk, outFname ):
    checkExeFound( afniSkullStripCmd, 'AFNI' )
    safeRun("mri_type_convert -chunk %s -short %s struct"%\
            (inChunk.name,inChunk.ds.fname))
    if inChunk.hasValue('afni.SCENE_DATA.0'):
        afniSpace= inChunk.getValue( 'afni.SCENE_DATA.0' )
    else:
        afniSpace= '+orig'
    safeRun("pghtoafni -anat struct struct%s"%afniSpace)
    argString= "-anat struct%s -prefix strippedStruct"%afniSpace
    if infoDict.has_key('min'):
        argString= "%s -min_val %f"%(argString,infoDict['min'])
    if infoDict.has_key('max'):
        argString= "%s -max_val %f"%(argString,infoDict['max'])
    if infoDict.has_key('minConn'):
        argString= "%s -min_conn %d"%(argString,infoDict['minConn'])
    if infoDict.has_key('maxConn'):
        argString= "%s -max_conn %d"%(argString,infoDict['maxConn'])
    outLines= readCmdOutputToList("%s %s"%(afniSkullStripCmd,argString))
    min= None
    max= None
    minConn= None
    maxConn= None
    for line in outLines:
        line= string.strip(line)
        if line.startswith("min value"):
            min= float(string.split(line)[-1])
        if line.startswith("max value"):
            max= float(string.split(line)[-1])
        if line.startswith("min conn"):
            minConn= int(string.split(line)[-1])
        if line.startswith("max conn"):
            maxConn= int(string.split(line)[-1])
    safeRun("smartreader -i strippedStruct%s.HEAD -out tmp_stripped"%\
            afniSpace)
    safeRun("mri_rpn_math -out %s '$2' struct tmp_stripped"%outFname)
    safeRun("mri_setfield -field images.skull_strip_info -value \"%s\" %s"%\
            ("afniSkullStrip(%s,%s,%s,%s)"%(min,max,minConn,maxConn),\
            outFname))
    safeRun("mri_history -add \"%s\" %s"%(string.join(sys.argv),outFname))

def afniSkullStripLikeThisOneApplicable( infoDict, inChunk, outFname ):
    if infoDict.has_key('like'):
        return 1
    else:
        return 0

def afniSkullStripLikeThisOne( infoDict, inChunk, outFname ):
    if not infoDict.has_key('like'):
        raise Exception,"no 'like' info provided"
    os.chdir(infoDict['homedir'])
    likeFName= os.path.abspath(infoDict['like'])
    os.chdir(infoDict['tmpdir'])
    likeDS= MRIDataset(likeFName)
    likeChunk= likeDS.getChunk('images')
    likeString= likeChunk.getValue('skull_strip_info')
    verboseMessage("Found strip info <%s>"%likeString)
    if likeString.startswith('afniSkullStrip(') and \
       likeString.endswith(')'):
        cloneDict= infoDict.copy()
        argList= string.split(likeString[len('afniSkullStrip('):-1],',')
        cloneDict['min']= float(argList[0])
        cloneDict['max']= float(argList[1])
        cloneDict['minConn']= int(argList[2])
        cloneDict['maxConn']= int(argList[3])
        afniSkullStrip(cloneDict, inChunk, outFname)
    elif likeString.startswith('afni3dSkullStrip(') and \
       likeString.endswith(')'):
        cloneDict= infoDict.copy()
        argList= string.split(likeString[len('afniSkullStrip('):-1],',')
        afni3dSkullStrip(cloneDict, inChunk, outFname)

mthdList= [ (afniSkullStripApplicable, afni3dSkullStripStrip),
            (afniSkullStripApplicable, afniSkullStrip),
            (afniSkullStripLikeThisOneApplicable, afniSkullStripLikeThisOne) ]

def realityCheck(inChunk,outChunk):
    global tmpdir
    xdim= inChunk.getDim('x')
    ydim= inChunk.getDim('y')
    zdim= inChunk.getDim('z')
    totdim= xdim*ydim*zdim
    safeRun("mri_copy_dataset %s copy"%inChunk.ds.fname)
    safeRun("mri_remap -order q -length %d copy"%totdim)
    safeRun("mri_subsample -d q -len 1 -max copy max")
    lines= readCmdOutputToList("mri_rpn_math '$1,1,if_print_1' max")
    max= float(lines[0])
    thresh= 0.03*max # this scaling constant is very iffy
    debugMessage("realityCheck: Input max: %f -> threshold %f"%(max,thresh))
    safeRun("mri_rpn_math -out copy_thresh '$1,%f,<=' copy"%thresh)
    safeRun("mri_subsample -d q -len 1 -sum copy_thresh copy_sum")
    safeRun("mri_rpn_math -out stripped_thresh '$1,%f,<=' %s"%(thresh,outChunk.ds.fname))
    safeRun("mri_remap -order q -length %d stripped_thresh"%totdim)
    safeRun("mri_subsample -d q -length 1 -sum stripped_thresh stripped_sum");
    lines= readCmdOutputToList("mri_rpn_math '$1,$2,1,if_print_2' copy_sum stripped_sum")
    words= string.split(lines[0])
    liveVoxelsInInput= int(words[0])
    liveVoxelsInOutput= int(words[1])
    debugMessage("realityCheck: live voxel counts %d and %d"%\
                 (liveVoxelsInInput,liveVoxelsInOutput))
    if liveVoxelsInOutput < 0.2*liveVoxelsInInput:
        raise Exception("Too many voxels stripped away! (%d/%d remain)"%\
                        (liveVoxelsInOutput,liveVoxelsInInput))

#################
#
# Main
#
#################

infoDict= {}

# Check for "-help"
if len(sys.argv)>1:
    if sys.argv[1] == "-help":
	if len(sys.argv)>2:
	    os.system( "scripthelp %s %s"%(sys.argv[0],sys.argv[2]) );
	else:
	    os.system( "scripthelp %s"%sys.argv[0] );
	sys.exit();

try:
    (opts,pargs) = getopt.getopt(sys.argv[1:],"vd",["like="])
except:
    print "%s: Invalid command line parameter" % sys.argv[0]
    describeSelf();
    sys.exit()

#Check calling syntax; parse args
if len(pargs) != 2 :
    describeSelf()
    sys.exit(1)

for a,b in opts:
    if a=="-v":
        setVerbose(1)
    if a=="-d":
        setDebug(1)
    if a=="--like":
        infoDict["like"]= b

inFname= os.path.abspath(pargs[0])
outFname= os.path.abspath(pargs[1])
inDS= MRIDataset(inFname)
inChunk= inDS.getChunk('images')

#Check reasonableness of input
dimstr= inChunk.getValue('dimensions');
if dimstr != "xyz":
    if dimstr == "xyzt":
        if inChunk.getDim("t") != 1:
            sys.exit("Input file must have t extent 1!")
    elif dimstr == "vxyzt":
        if inChunk.getDim("t") != 1:
            sys.exit("Input file must have t extent 1!")
        if inChunk.getDim("v") != 1:
            sys.exit("Input file must have v extent 1!")
    elif dimstr == "vxyz":
        if inChunk.getDim("v") != 1:
            sys.exit("Input file must have v extent 1!")
        else:
            sys.exit("Input file must have dimensions (v)xyz(t)!")

# Move to a temporary directory 
tmpdir= os.path.abspath(makeTempDir('tmp_strip_skull'))
homedir= os.getcwd()
infoDict['tmpdir']= tmpdir
infoDict['homedir']= homedir
os.chdir(tmpdir)

flags= ""
if (getVerbose()):
    flags= flags+"-v "
if (getDebug()):
    flags= flags+"-d "
infoDict['flags']= flags

if not outputAlreadyExists(inChunk,outFname):
    for (mthdApplicable, mthd) in mthdList:
        if mthdApplicable(infoDict, inChunk, outFname):
            try:
                # start in clean directory
                for file in os.listdir(tmpdir):
                    kidPath= os.path.join(tmpdir,file)
                    if os.path.isdir(kidPath):
                        removeTmpDir(kidPath)
                    else:
                        os.remove(kidPath)
                outDS= None
                verboseMessage("Trying stripping method %s"%mthd.__name__)
                mthd(infoDict,inChunk,outFname)
                outDS= MRIDataset(outFname)
                outChunk= outDS.getChunk("images")
                realityCheck(inChunk,outChunk)
                break
            except Exception, inst:
                verboseMessage("Method %s failed (%s); cleaning up"%\
                               (mthd.__name__,inst))
                if outDS != None:
                    safeRun("mri_destroy_dataset %s"%outDS.fname)
                    outDS= None
    if outDS==None:
        sys.exit("All known methods for skull stripping have failed.")

else:
    verboseMessage("Output file %s exists and is new."%outFname)

# Clean up
os.chdir(homedir)
if not getDebug():
    removeTmpDir(tmpdir)


