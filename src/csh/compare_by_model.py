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
import array
import threading
import popen2
from math import *
if os.environ.has_key("FIASCO"):
    sys.path.append(os.environ["FIASCO"])
from fiasco_utils import *

idString= "$Id: compare_by_model.py,v 1.3 2006/05/04 23:02:58 welling Exp $"

################################
# Notes-
################################

#################
#
# Main
#
#################

# Check for "-help"
if len(sys.argv)>1:
    if sys.argv[1] == "-help":
	if len(sys.argv)>2:
	    os.system( "scripthelp %s %s"%(sys.argv[0],sys.argv[2]) );
	else:
	    os.system( "scripthelp %s"%sys.argv[0] );
	sys.exit();

try:
    (opts,pargs) = getopt.getopt(sys.argv[1:],"vd",\
                                 ["model=", "ordered","unordered",
                                  "contrasts=", "split=","scale="])
except:
    print "%s: Invalid command line parameter" % sys.argv[0]
    describeSelf();
    sys.exit()

#Check calling syntax; parse args
if len(pargs) != 1 :
    errorMessage("Required split file name not given.")
    describeSelf()
    sys.exit(1)
inFname= os.path.abspath( pargs[0] )

scaleFname= None
formula= None
splitFnameList= []
orderedFlag= True
orderedSet= False
contrastType= None
for a,b in opts:
    if a=="-v":
        setVerbose(1)
    elif a=="-d":
        setDebug(1)
    elif a=="--split":
        splitFnameList.append( os.path.abspath(b) )
    elif a=="--scale":
        scaleFname= os.path.abspath(b)
    elif a=="--model":
        formula= b
    elif a=="--ordered":
        orderedFlag= True
        orderedSet= True
    elif a=="--unordered":
        orderedFlag= False
        orderedSet= True
    elif a=="--contrasts":
        contrastType= b
    
verboseMessage(idString)

if formula==None:
    sys.exit("Required model formula not given")
if len(splitFnameList)==0:
    sys.exit("Required split file name not given")

# Make a temporary directory, and go there
tmpdir= makeTempDir('tmp_compare_by_model')
homedir= os.getcwd()
os.chdir(tmpdir)

dz= getDim(inFname,"images","z")
dt= getDim(inFname,"images","t")

verboseMessage("Assembling model matrix")
argString= ""
if getVerbose(): argString += "-v "
if getDebug(): argString += "-v " # -d makes too much output
argString += "--nslices %d "%dz
argString += "--nimages %d "%dt
argString += "--model '%s' "%formula
if contrastType != None:
    argString += "--contrasts %s "%contrastType
if orderedSet:
    if orderedFlag: argString += "--ordered "
    else: argString += "--unordered "
for fname in splitFnameList:
    argString += "%s "%fname
safeRun("build_model_matrix.py --out xmatrix %s"%argString)
xmatrixDS= MRIDataset("xmatrix")
xmatrixChunk= xmatrixDS.getChunk("images")
nCols= xmatrixChunk.getDim("v")

factorNameList= ["intercept"]
for col in xrange(nCols):
    factorName= xmatrixChunk.getValue("label.v.%d"%col)
    factorNameList.append(factorName)

verboseMessage("Collecting missing data")
safeRun("mri_copy_dataset %s inData"%inFname)
safeRun("mri_subsample -d v -len 1 -sum xmatrix xm_sum")
safeRun("mri_remap -order tz xm_sum")
safeRun("mri_rpn_math -out design_missing_p '1,0,$1,is_finite,if_keep' xm_sum")
safeRun("mri_permute -order zt design_missing_p design_missing")
if xmatrixDS.hasChunk('missing'):
    safeRun("mri_copy_chunk -chunk missing -chunk_out images inData old_mis")
    safeRun("mri_rpn_math -out all_mis '0,1,$1,$2,+,if_keep' design_missing old_mis")
    safeRun("mri_type_convert -char all_mis all_mis_char")
    safeRun("mri_copy_chunk -chunk images -chunk_out missing -replace all_mis_char inData")
else:
    safeRun("mri_type_convert -char design_missing all_mis_char")
    safeRun("mri_copy_chunk -chunk images -chunk_out missing -replace all_mis_char inData")
safeRun("mri_rpn_math -out xmatrix_filt '0,$1,dup,is_finite,if_keep' xmatrix")

verboseMessage("Performing regression")
argString= "-variance "
if getVerbose(): argString += "-verbose "
#if getDebug(): argString += "-debug " # produces too much output!
if scaleFname != None: argString += "-scale %s:images"%scaleFname
safeRun("mri_glm %s -estimates est:images inData:images xmatrix_filt:images"%\
        argString)

verboseMessage("Assembling result datasets:")
nFactors= len(factorNameList)
for facOffset in xrange(nFactors):
    verboseMessage("    Beta and variance for %s"%factorNameList[facOffset])
    betaName= os.path.join(homedir,"beta_%s"%factorNameList[facOffset])
    varName=  os.path.join(homedir,"var_%s"%factorNameList[facOffset])
    safeRun("mri_subset -d v -len 1 -s %d est %s"%(facOffset,betaName))
    safeRun("mri_subset -d v -len 1 -s %d est %s"%(facOffset+nFactors,
                                                   varName))

# Clean up
os.chdir(homedir)
if not getDebug():
    removeTmpDir(tmpdir)


