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
#  -afniTest needs to check for presence of adwarp in path!
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

idString= "$Id: normalize_prep.py,v 1.8 2006/05/04 23:02:58 welling Exp $"

afniNormalizeCmd= "adwarp"

afniGeneratedScript="""\
#!/bin/csh -f
#
# This is an automatically generated script!  It's only good
# for this particular special case.  Don't copy it around or
# edit it unless you *really* know what you're doing!
#
# Switches: (must appear first)
#   -anat : treat input as anatomical rather than functional
#
# Inputs: (after any switches)
#   $1: input functional MRI dataset
#   $2: output normalized functional MRI dataset
#
set typeflag = '-func'
if ( "$1" == '-anat' ) then
  shift argv
  set typeflag = '-anat'
endif
set prototype = %s
if ( ! -e ${prototype}.HEAD ) then
  echo "Expected prototype file ${prototype}.HEAD not found!"
  exit -1
endif
if ( ! -e ${prototype}.BRIK ) then
  echo "Expected prototype file ${prototype}.BRIK not found!"
  exit -1
endif
pghtoafni $typeflag $1 tmp_normalize+orig
%s -func NN -apar ${prototype} -dpar tmp_normalize+orig
smartreader -input tmp_normalize+tlrc.HEAD -out $2
rm tmp_normalize+orig.HEAD tmp_normalize+orig.BRIK
rm tmp_normalize+tlrc.HEAD tmp_normalize+tlrc.BRIK

"""

def generalPrep():
    global homedir
    path1= os.path.join(homedir,'normalized_anat')
    path2= os.path.join(homedir,'normalized_func')
    for path in [ path1, path2 ]:
        if not os.access(path,os.W_OK):
            os.makedirs(path)

def afniPrep(inFname):
    global homedir
    global tmpdir
    global oscript

    checkExeFound(afniNormalizeCmd,'AFNI')

    baseName= os.path.basename(inFname)
    headName= inFname+'+tlrc.HEAD'
    origName= "%s+orig"%baseName
    warpName= "norm_%s"%baseName
    finalName= os.path.join(homedir,"normalized_anat","norm_%s"%baseName)
    safeRun("pghtoafni -anat %s %s"%(inFname,os.path.join(tmpdir,origName)))
    os.chdir(tmpdir)
    safeRun("%s -apar %s -dpar %s -prefix %s"%\
            (afniNormalizeCmd,os.path.join(homedir,headName),\
             origName,os.path.basename(warpName)))
    afniTlrcName= "%s+tlrc"%os.path.basename(warpName)
    safeRun("cp %s.* %s"%(afniTlrcName,\
                          os.path.join(homedir,"normalized_anat")))
    safeRun("smartreader -i %s+tlrc.HEAD -out %s"%(warpName,finalName))
    # Now we write a shell script that implements this transformation.
    scriptName= os.path.join(homedir,oscript);
    if (os.access(scriptName,os.F_OK)):
        os.remove(scriptName);
    ofile= open(scriptName,"w")
    ofile.write(afniGeneratedScript%(os.path.join(homedir,"normalized_anat",\
                                                  afniTlrcName),
                                     afniNormalizeCmd))
    ofile.close()
    safeRun("chmod ug+x %s"%scriptName)

def afniTest( inFname ):
    # If there is a file with that name + '+tlrc.HEAD', they're
    # using AFNI methods.
    headFname= inFname+'+tlrc.HEAD'
    if os.access(headFname,os.R_OK):
        return 1;
    else:
        return 0;

methodList= [(afniTest, afniPrep)] # we'll have more someday

#################
#
# Main
#
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

idString= "$Id: normalize_prep.py,v 1.8 2006/05/04 23:02:58 welling Exp $"

# Check for "-help"
if len(sys.argv)>1:
    if sys.argv[1] == "-help":
	if len(sys.argv)>2:
	    os.system( "scripthelp %s %s"%(sys.argv[0],sys.argv[2]) );
	else:
	    os.system( "scripthelp %s"%sys.argv[0] );
	sys.exit();

try:
    (opts,pargs) = getopt.getopt(sys.argv[1:],"vd",["oscript="])
except:
    print("%s: Invalid command line parameter" % sys.argv[0])
    describeSelf();
    sys.exit()

#Check calling syntax; parse args
if len(pargs) != 1 :
    describeSelf()
    sys.exit(1)

oscript= "apply_normalize.gen_csh"
for a,b in opts:
    if a=="-v":
        setVerbose(1)
    if a=="-d":
        setDebug(1)
    if a=="--oscript":
        oscript= b

inFname= pargs[0]
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

# Create a temporary directory 
tmpdir= makeTempDir('tmp_normalize_prep')
homedir= os.getcwd()

# General setup
generalPrep()

# Figure out what game we're playing
done= 0
for test, mthd in methodList:
    if test(inFname):
        try:
            mthd(inFname)
            done= 1
        except Exception, inst:
            sys.exit("Cannot normalize from %s: %s"%(inFname,inst))
if not done:
    sys.exit("Cannot normalize from %s: can't find needed inputs"%inFname)

# Clean up
os.chdir(homedir)
if not getDebug():
    removeTmpDir(tmpdir)


