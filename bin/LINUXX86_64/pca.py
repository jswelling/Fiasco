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
if os.environ.has_key("FIASCO"):
    sys.path.append(os.environ["FIASCO"])
from fiasco_utils import *

idString= "$Id: pca.py,v 1.6 2006/05/04 23:02:58 welling Exp $"

def pickUniqueIndex( dimstr ):
    for c in string.ascii_letters:
        if dimstr.find(c)<0:
            return c
    sys.exit("Ran out of possible dimension indices!")

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

inDSName= None
evalDSName= None
evecDSName= None
leftvecDSName= None
nComponents= None

# default matrix size to switch to a covariance-based method
methodThreshold= 256 

try:
    (opts,pargs) = getopt.getopt(sys.argv[1:],"vdn:",\
                                 ["eigenvectors=","eigenvalues=",\
                                  "leftvectors=","methodthresh="])
except:
    print "%s: Invalid command line parameter" % sys.argv[0]
    describeSelf();
    sys.exit()

for a,b in opts:
    if a=="-v":
        setVerbose(1)
    if a=="-d":
        setDebug(1)
    if a=="-n":
        nComponents= int(b)
    if a=="--eigenvectors":
        evecDSName= b
    if a=="--eigenvalues":
        evalDSName= b
    if a=="--leftvectors":
        leftvecDSName= b
    if a=="--methodthresh":
        methodThreshold= int(b)

if len(pargs) < 1:
    sys.exit("Required input dataset name not specified")
elif len(pargs) > 1:
    sys.exit("Too many dataset names specified")
inDSName= pargs[0]

if evecDSName==None:
    sys.exit("Required eigenvector dataset name not specified")

fullInDSName= os.path.abspath(inDSName)
fullEvecDSName= os.path.abspath(evecDSName)
if evalDSName != None:
    fullEvalDSName= os.path.abspath(evalDSName)
else:
    fullEvalDSName= None
if leftvecDSName != None:
    fullLeftvecDSName= os.path.abspath(leftvecDSName)
else:
    fullLeftvecDSName= None

# Do the input and output paths exist?
if not dsExists(fullInDSName):
    sys.exit("Input dataset %s does not exist!"%inDSName)
checkPathExists( fullEvecDSName )
if evalDSName != None:
    checkPathExists( fullEvalDSName )
if leftvecDSName != None:
    checkPathExists( fullLeftvecDSName )

# Is the input dataset acceptable?
inDS= MRIDataset(inDSName)
if not inDS.hasChunk("images"):
    sys.exit("%s has no images chunk!"%inDSName)
images= inDS.getChunk("images")
if not images.hasValue("dimensions"):
    sys.exit("%s images chunk has no dimensions!"%inDSName)
dimstr= images.getValue("dimensions")
if len(dimstr)<2:
    sys.exit("Input dataset must be at least 2D")

tChar= dimstr[0]
tDim= images.getDim(tChar)
sChar= dimstr[1]
sDim= images.getDim(sChar)
loopString= dimstr[2:]
loopStringSpacer= len(loopString)*':'

if nComponents == None:
    nComponents= sDim

verboseMessage("Input dataset %s chunk images: dims %d %d, loop dims <%s>"%\
               (inDSName,tDim,sDim,loopString))

# Make a temp directory
tmpdir= makeTempDir('tmp_pca')
homedir= os.getcwd()

#
# We use one of three methods, depending on the situation.  If
# both dimensions are smaller than some threshold, we simply do
# an SVD and declare victory.
#
# If one or more dimensions are above the threshold, we use a
# method based on the covariance matrix.
#
# If the first dimension is larger than the second, calculate xTx
# (x-transpose x) and solve the eigensystem.  If the second
# dimension is larger, calculate xxT, solve the eigensystem, and
# then substitute back to get the eigenvectors of x.
#
os.chdir(tmpdir)
uChar= pickUniqueIndex(dimstr)
if tDim <= methodThreshold and sDim <= methodThreshold:
    verboseMessage("Using singular value decomposition")
    smallDim = min(tDim,sDim)
    svdArgList = ""
    if fullLeftvecDSName != None:
        svdArgList= svdArgList + "-umatrix g "
    if fullEvalDSName != None:
        svdArgList= svdArgList + "-wvector l "
    if fullEvecDSName != None:
        svdArgList= svdArgList + "-vmatrix dt "
    safeRun("mri_svd %s %s"%(svdArgList,fullInDSName))
    if fullEvalDSName != None:
        safeRun("mri_remap -order %c%s -len %d%s l"%\
                (uChar,loopString,smallDim,loopStringSpacer));
        if smallDim > nComponents:
            safeRun("mri_subset -d %c -len %d -s 0 l %s"%\
                    (tChar, nComponents, fullEvalDSName))
        else:
            safeRun("mri_copy_dataset l %s"%fullEvalDSName)
    if fullLeftvecDSName != None:
        safeRun("mri_remap -order %c%c%s -len :%d%s g"%\
                (tChar,uChar,loopString,smallDim,loopStringSpacer));
        if smallDim > nComponents:
            safeRun("mri_subset -d %c -len %d -s 0 g %s"%\
                    (uChar, nComponents, fullLeftvecDSName))
        else:
            safeRun("mri_copy_dataset g %s"%fullLeftvecDSName)
    if fullEvecDSName != None:
        safeRun("mri_remap -order %c%c%s -len %d:%s dt"%\
                (uChar,sChar,loopString,sDim,loopStringSpacer))
        if smallDim > nComponents:
            safeRun("mri_subset -d %c -len %d -s 0 dt dt_subset"%\
                    (uChar,nComponents))
        else:
            safeRun("mri_copy_dataset dt dt_subset")
        safeRun("mri_permute -order %c%c%s dt_subset %s"%\
                (sChar,uChar,loopString,fullEvecDSName))

elif tDim < sDim:
    verboseMessage("Using x-x-transpose method")
    if fullLeftvecDSName == None:
        fullLeftvecDSName= os.path.join(tmpdir,"g")
    if fullEvalDSName == None:
        fullEvalDSName= os.path.join(tmpdir,"l")
    safeRun("mri_permute -o %c%c%s %s xt"%\
            (sChar, tChar, loopString, fullInDSName))
    safeRun("mri_remap -o %c%c%s -len :%d%s xt"%\
            (sChar,uChar,loopString,tDim,loopStringSpacer))
    safeRun("mri_matmult -out xxt %s xt"%fullInDSName)
    safeRun("mri_esa -descend -num %d -eigenvectors %s xxt ll"%\
            (nComponents,fullLeftvecDSName))
    safeRun("mri_remap -o %c%s -len %d: ll"%(uChar,loopString,nComponents))
    safeRun("mri_rpn_math -out %s '$1,sqrt' ll"%fullEvalDSName)
    safeRun("mri_permute -o %c%c%s %s gt"%\
            (uChar,tChar,loopString,fullLeftvecDSName))
    safeRun("mri_matmult -out dt gt %s"%fullInDSName)
    safeRun("mri_remap -order %c%c%s -len %d:%s dt"%\
            (uChar,sChar,loopString,nComponents,loopStringSpacer))
    safeRun("mri_rpn_math -out dt_scaled '$1,$2,/' dt %s"%fullEvalDSName);
    safeRun("mri_permute -o %c%c%s dt_scaled %s"%\
            (sChar,uChar,loopString,fullEvecDSName))

else:
    verboseMessage("Using x-transpose-x method")
    if fullEvalDSName == None and fullLeftvecDSName != None:
        fullEvalDSName= os.path.join(tmpdir,"l")
    safeRun("mri_permute -o %c%c%s %s xt"%\
            (sChar, tChar, loopString, fullInDSName))
    safeRun("mri_remap -o %c%c%s -len %d:%s xt"%\
            (uChar,tChar,loopString,sDim,loopStringSpacer))
    safeRun("mri_matmult -out xtx xt %s"%fullInDSName)
    safeRun("mri_esa -descend -num %d -eigenvectors %s xtx ll"%\
            (nComponents,fullEvecDSName))
    safeRun("mri_remap -order %c%c%s -len %d:%d%s %s"%\
            (sChar,uChar,loopString,sDim,nComponents,loopStringSpacer,\
             fullEvecDSName))
    if fullEvalDSName != None: 
        safeRun("mri_remap -o %c%s -len %d%s ll"%\
                (uChar,loopString,nComponents,loopStringSpacer))
        safeRun("mri_rpn_math -out %s '$1,sqrt' ll"%fullEvalDSName)
    if fullLeftvecDSName != None:
        safeRun("mri_matmult -out g %s %s"%\
                (fullInDSName, fullEvecDSName))
        safeRun("mri_permute -order %c%c%s g gt"%\
                (uChar,tChar,loopString))
        safeRun("mri_rpn_math -out gt_scaled '$1,$2,/' gt %s"%fullEvalDSName)
        safeRun("mri_permute -order %c%c%s gt_scaled %s"%
                (tChar,uChar,loopString,fullLeftvecDSName))

# Clean up
os.chdir(homedir)
if not getDebug():
    removeTmpDir(tmpdir)
