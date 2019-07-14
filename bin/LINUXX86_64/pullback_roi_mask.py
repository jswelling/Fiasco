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

idString= "$Id: pullback_roi_mask.py,v 1.4 2006/05/04 23:02:58 welling Exp $"

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

srcdir= None
thresh= 1.0
strippedThresh= 1.0
try:
    (opts,pargs) = getopt.getopt(sys.argv[1:],"vd",
                                 ["srcdir=","thresh=","strippedthresh="])
except:
    print "%s: Invalid command line parameter" % sys.argv[0]
    describeSelf();
    sys.exit()

for a,b in opts:
    if a=="-v":
        setVerbose(1)
    if a=="-d":
        setDebug(1)
    if a=="--srcdir":
        srcdir= os.path.abspath(b)
    if a=="--thresh":
        thresh= float(b)
    if a=="--strippedthresh":
        strippedThresh= float(b)

if srcdir==None:
    sys.exit("Required srcdir value not specified")
if len(pargs) != 2:
    sys.exit("a single input mask dataset argument is required")

inDS= os.path.abspath(pargs[0])
if not dsExists(inDS):
    sys.exit("Input dataset %s not found!"%inDS)

# Check for needed scripts and input datasets in the source directory
protoDS= os.path.join(srcdir,"stat","GrandMean")
strippedAxialDS= os.path.join(srcdir,"anat","stripped_axial")
coregStrippedAxialDS= os.path.join(srcdir,"coregistered_anat","warped_stripped_axial")
normStrippedAxialDS= os.path.join(srcdir,"normalized_anat",
                                  "norm_warped_stripped_axial")
for dsNamePair in [ (protoDS,"GrandMean"),
                    (strippedAxialDS,"stripped_axial") ]:
    dsName,strName= dsNamePair
    if not dsExists(dsName):
        sys.exit("Expected %s dataset not found below source dir!"%strName)
coregScript= os.path.join(srcdir,"generated_scripts",
                          "apply_coreg_inplane.gen_csh")
normScript= os.path.join(srcdir,"generated_scripts",
                         "apply_normalize.gen_csh")
coregStructScript= os.path.join(srcdir,"generated_scripts",
                                "apply_coreg_struct_to_inplane.gen_csh")
cowarpScript= os.path.join(srcdir,"generated_scripts",
                           "apply_cowarp_inplane.gen_csh")
for scriptNamePair in [ (coregScript,"apply_coreg_inplane.gen_csh"),
                        (normScript,"apply_normalize.gen_csh"),
                        (coregStructScript,
                         "apply_coreg_struct_to_inplane.gen_csh"),
                        (cowarpScript,"apply_cowarp_inplane.gen_csh")]:
    scriptName, strName= scriptNamePair
    if not os.access(scriptName, os.X_OK):
        sys.exit("Expected %s script not found "%strName
                 +"below source dir, or not executable!")

outDS= os.path.abspath(pargs[1])
checkPathExists(outDS)

# Make a temp directory
tmpdir= makeTempDir('tmp_pullback_roi_mask')
homedir= os.getcwd()
os.chdir(tmpdir)

# Make an address dataset in Talairach space
safeRun("mri_copy_dataset %s proto"%protoDS)
ch= MRIDataset("proto").getChunk("images")
xdim= ch.getDim('x')
ydim= ch.getDim('y')
zdim= ch.getDim('z')
qdim= xdim*ydim*zdim
safeRun("mri_remap -preserve -order q -len %d proto"%qdim)
safeRun("mri_rpn_math -out addr '$q' proto")
safeRun("mri_remap -order xyz -len %d:%d:%d addr"%(xdim,ydim,zdim))
safeRun("%s -closest addr addr_coreg junk"%coregScript)
safeRun("mri_copy_dataset addr addr_coreg")
safeRun("%s addr_coreg addr_tlrc"%normScript)

# Update the normalized stripped axial dataset
if not dsExists(normStrippedAxialDS) \
   or isNewer("%s.mri"%strippedAxialDS,"%s.mri"%normStrippedAxialDS) \
   or isNewer(coregStructScript,"%s.mri"%normStrippedAxialDS) \
   or isNewer(cowarpScript,"%s.mri"%normStrippedAxialDS) \
   or isNewer(normScript,"%s.mri"%normStrippedAxialDS):
   
    if not dsExists(coregStrippedAxialDS) \
       or isNewer("%s.mri"%strippedAxialDS,"%s.mri"%coregStrippedAxialDS) \
       or isNewer(coregStructScript,"%s.mri"%coregStrippedAxialDS) \
       or isNewer(cowarpScript,"%s.mri"%coregStrippedAxialDS):
        safeRun("%s %s stripped_coreg junk"%(coregStructScript,strippedAxialDS))
        safeRun("%s stripped_coreg %s"%(cowarpScript,coregStrippedAxialDS))
        
    safeRun("%s %s %s"%(normScript,coregStrippedAxialDS,normStrippedAxialDS))

# Assemble the mask in Talairach space
safeRun("mri_rpn_math -out mask_masked '0,$1,$2,%f,<=,if_keep' %s %s"%\
        (strippedThresh,inDS,normStrippedAxialDS))

# Pull the mask back to orig space
safeRun("pullback -qdim %d -addr addr_tlrc mask_masked mask_raw"%qdim)
safeRun("mri_remap -order xyz -len %d:%d:%d mask_raw"%(xdim,ydim,zdim))
safeRun("mri_rpn_math -out %s '$2,%f,<=' %s mask_raw"%\
        (outDS,thresh,protoDS))
        
# Clean up
os.chdir(homedir)
if not getDebug():
    removeTmpDir(tmpdir)
