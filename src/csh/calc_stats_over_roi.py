#! /usr/bin/env python
#
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

import sys
import os
import os.path
import string
import getopt
if "FIASCO" in os.environ:
    sys.path.append(os.environ["FIASCO"])
from fiasco_utils import *

idString= "$Id: calc_stats_over_roi.py,v 1.4 2007/02/14 00:04:29 welling Exp $"

basePath= os.environ['PWD']

def maybeBuildNormalizedAnat(dsName):
    global basePath
    normPath= os.path.join(basePath, "normalized_anat","norm_warped_%s"%dsName)
    coregPath= os.path.join(basePath,"coregistered_anat", \
                            "warped_%s"%dsName)
    origPath= os.path.join(basePath,"anat",dsName)
    coregScript= os.path.join(basePath,"generated_scripts",\
                              "apply_coreg_struct_to_inplane.gen_csh")
    warpScript= os.path.join(basePath,"generated_scripts",\
                             "apply_cowarp_inplane.gen_csh")
    normScript= os.path.join(basePath,"generated_scripts",\
                              "apply_normalize.gen_csh")
    if not dsExists(origPath):
        sys.exit("Required input dataset <%s> not found!"%origPath)
    if not dsExists(coregPath) or isNewer("%s.mri"%origPath,
                                          "%s.mri"%coregPath):
        verboseMessage("Coregistering %s"%dsName);
        safeRun("%s %s tmp_coreg tmp2_coreg > /dev/null 2>&1"%\
                (coregScript,origPath))
        safeRun("%s tmp_coreg %s > /dev/null 2>&1"%(warpScript,coregPath))
        safeRun("mri_destroy_dataset tmp_coreg")
        safeRun("mri_destroy_dataset tmp2_coreg")
    if not dsExists(normPath) or isNewer("%s.mri"%coregPath,"%s.mri"%normPath):
        verboseMessage("Normalizing %s"%dsName);
        safeRun("%s -anat %s %s > /dev/null 2>&1"%\
                (normScript,coregPath,normPath))
    return normPath

def maybeBuildNormalizedFunc(dsName):
    global basePath
    normPath= os.path.join(basePath, "normalized_func","norm_%s"%dsName)
    coregPath= os.path.join(basePath,"coregistered_func", \
                            "aligned_%s"%dsName)
    origPath= os.path.join(basePath,"stat",dsName)
    coregScript= os.path.join(basePath,"generated_scripts",\
                              "apply_coreg_inplane.gen_csh")
    normScript= os.path.join(basePath,"generated_scripts",\
                              "apply_normalize.gen_csh")
    if not dsExists(origPath):
        sys.exit("Required input dataset <%s> not found!"%origPath)
    if not dsExists(coregPath) or isNewer("%s.mri"%origPath,
                                          "%s.mri"%coregPath):
        verboseMessage("Coregistering %s"%dsName);
        safeRun("%s %s %s tmp_coreg > /dev/null 2>&1"%\
                (coregScript,origPath,coregPath))
        safeRun("mri_destroy_dataset tmp_coreg")
    if not dsExists(normPath) or isNewer("%s.mri"%coregPath,"%s.mri"%normPath):
        verboseMessage("Normalizing %s"%dsName);
        safeRun("%s %s %s > /dev/null 2>&1"%(normScript,coregPath,normPath))
    return normPath

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
    (opts,pargs) = getopt.getopt(sys.argv[1:],"vdl:h:r:",["basepath="])
except:
    print("%s: Invalid command line parameter" % sys.argv[0])
    describeSelf();
    sys.exit()

# These default thresholds are against a typical axial anatomical scan,
# which has values in the range 0 < v < 1000 or so.
lowThresh= 0.0
highThresh= 1000000
subj= None
roiMaskName= None
for a,b in opts:
    if a=="-v":
        setVerbose(1)
    if a=="-d":
        setDebug(1)
    if a=="-l":
        lowThresh= float(b)
    if a=="-h":
        highThresh= float(b)
    if a=="-r":
        roiMaskName= b
    if a=="--basepath":
        basePath= b

if roiMaskName==None:
    sys.exit("Required raw ROI mask dataset not specified")

if len(pargs) != 1:
    sys.exit("Failed to specify which statistic to sample")
sampleThis= pargs[0]

# Make a temp directory and go there
tmpdir= makeTempDir('tmp_calc_stats_over_roi')
homedir= os.getcwd()
os.chdir(tmpdir)

# Make sure we have the needed input fields
fullStrippedAxialDSName= maybeBuildNormalizedAnat("stripped_axial")
fullStatDSName= maybeBuildNormalizedFunc(sampleThis)

# Calculate median and IQR
valsInROI= \
   readCmdOutputToList("mri_rpn_math '$1,dup,is_finite,$2,*,$3,%f,<=,*,$3,%f,>=,*,if_print_1' %s %s %s"%\
                       (lowThresh,highThresh,fullStatDSName,\
                        os.path.join(homedir,roiMaskName),\
                        fullStrippedAxialDSName))
debugMessage("Found %d valid samples"%len(valsInROI))
if len(valsInROI)==0:
    print("0.0 0.0 0.0 0.0 0.0 0")
else:
    floatsInROI= [float(x) for x in valsInROI]
    floatsInROI.sort()
    n= len(valsInROI)
    q1Loc= int(n/4)
    medianLoc= int(n/2)
    q3Loc= int((3*n)/4)
    print("%r %r %r %r %r %d"%\
          (floatsInROI[0], floatsInROI[q1Loc], floatsInROI[medianLoc],
           floatsInROI[q3Loc], floatsInROI[-1], n))
    
# Clean up
os.chdir(homedir)
if not getDebug():
    removeTmpDir(tmpdir)

