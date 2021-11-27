#! /usr/bin/env python
import sys
import os
import string
import getopt
import copy
import math
if "FIASCO" in os.environ:
    sys.path.append(os.environ["FIASCO"])
from fiasco_utils import *

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
    (opts,pargs) = getopt.getopt(sys.argv[1:],"vd",[])
except:
    print ("%s: Invalid command line parameter" % sys.argv[0])
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

inFname= os.path.abspath(pargs[0])
outFname= os.path.abspath(pargs[1])
inDS= MRIDataset(inFname)
inChunk= inDS.getChunk('images')

#Check that the needed executable is in the path
checkExeFound("3dTshift","AFNI")

#Check reasonableness of input
if not inChunk.getValue('reorder_pattern'):
    sys.exit("Input file is missing the required reorder_pattern tag.")
verboseMessage("Reorder pattern <%s>"%inChunk.getValue('reorder_pattern'))
dimstr= inChunk.getValue('dimensions');
if dimstr == "xyzt":
    if inChunk.getDim("t") == 1:
        sys.exit("Input file must have t extent greater than 1!")
elif dimstr == "vxyzt":
    if inChunk.getDim("t") == 1:
        sys.exit("Input file must have t extent greater than 1!")
    if inChunk.getDim("v") != 1:
        sys.exit("Input file must have v extent 1!")
else:
    sys.exit("Input file must have dimensions (v)xyz(t)!")

# Move to a temporary directory
tmpdir= os.path.abspath(makeTempDir('tmp_afni_tshift'))
homedir= os.getcwd()
os.chdir(tmpdir)

flags= ""
afniFlags= ""
if (getVerbose()):
    flags= flags+"-verbose "
    afniFlags= afniFlags+"-verbose "

# We have to smooth over missing data before time shifting 
if dimstr=="xyzt":
    safeRun("mri_permute -order txyz %s permuted"%inFname)
else:
    safeRun("mri_permute -order vtxyz %s permuted"%inFname)
    safeRun("mri_remap -order txyz permuted")
safeRun("smooth_missing -i permuted -h nomissing_p")
safeRun("mri_destroy_dataset permuted")
safeRun("mri_permute -order xyzt nomissing_p nomissing")
safeRun("mri_destroy_dataset nomissing_p")
safeRun("pghtoafni %s -func nomissing raw+orig"%flags)
safeRun("3dTshift -Fourier %s raw+orig"%afniFlags)
safeRun("smartreader %s -i tshift+orig.HEAD -out tmp"%flags)
tmpDS= MRIDataset('tmp')
tmpChunk= tmpDS.getChunk('images')
safeRun("mri_rpn_math %s -out %s '$2' %s tmp"%(flags,outFname,inFname))
safeRun("mri_setfield -fld images.reorder_pattern -delete %s"%\
        outFname)
safeRun("mri_history -add \"%s\" %s"%(string.join(sys.argv),outFname))    

# Clean up
os.chdir(homedir)
if not getDebug():
    removeTmpDir(tmpdir)
