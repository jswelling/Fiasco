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

idString= "$Id: pick_n_voxels.py,v 1.3 2006/05/04 23:02:58 welling Exp $"

def buildExtentStr(fname,dimstr,totpix):
    extstr= ""
    for i in range(len(dimstr)):
        thischar= dimstr[i]
        if thischar=="q":
            extstr= extstr + "%d:"%totpix
        else:
            extstr= extstr + "%d:"%getDim(fname,thischar) 
    extstr= extstr[0:-1] # remove trailing ':'
    return extstr

##############################
#
# Main
#
##############################

maskDS= None
maskChunk= None
foldFlag= 0
tFlag= 0 # input is Tmap
pFlag= 0 # input is Pmap
fFlag= 0 # input is Fmap

# Check for "-help"
if len(sys.argv)>1:
    if sys.argv[1] == "-help":
        if len(sys.argv)>2:
            os.system( "scripthelp %s %s"%(sys.argv[0],sys.argv[2]) );
        else:
            os.system( "scripthelp %s"%sys.argv[0] );
        sys.exit();

try:
    (opts,pargs) = getopt.getopt(sys.argv[1:],"vdTPF",["mask=","twotails"])
except:
    print("%s: Invalid command line parameter" % sys.argv[0])
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
    if a=="-T":
        tFlag= 1
    if a=="-P":
        pFlag= 1
    if a=="-F":
        fFlag= 1
    if a=="--mask":
        maskDS= MRIDataset(b)
        maskChunk= maskDS.getChunk('images')
    if a=="--twotails":
        foldFlag= 1

nTypeFlags= tFlag + pFlag + fFlag
if nTypeFlags == 0:
    sys.exit("Either -T, -P, or -F must be specified!")

if nTypeFlags>1:
    sys.exit("Only specify one of -T, -P, or -F!")

if fFlag and foldFlag:
    sys.exit("The --twotails flag is incompatible with -F!");

nVox= string.atoi(pargs[0])
inDS= MRIDataset(pargs[1])
inChunk= inDS.getChunk('images')

#Make up a temporary directory
tmpdir= makeTempDir('tmp_pickn')
homedir= os.getcwd()

#Get relevant dimensions
xdim= inChunk.getDim("x");
ydim= inChunk.getDim("y");
zdim= inChunk.getDim("z");
totpix= xdim*ydim*zdim
dimstr= inChunk.getValue('dimensions');

#Check reasonableness of input
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
if maskChunk != None:
    if maskChunk.getValue('dimensions') != dimstr:
        sys.exit("Dimension order of mask doesn't match that of input!")
    for dim in ['v', 'x', 'y', 'z', 't']:
        if maskChunk.hasValue(dim):
            if maskChunk.getDim(dim) != inChunk.getDim(dim):
                sys.exit("%s extent of mask doesn't match that of input!"%dim)

#Generate sorted file
os.chdir(tmpdir)
if maskDS == None:
    if foldFlag:
        if pFlag: # Pmap case
            safeRun("mri_rpn_math -out unsorted '1.0,$1,-,$1,dup,0.5,>,if_keep' %s"%\
                    os.path.join(homedir,inDS.fname))
        elif tFlag: # Tmap case
            safeRun("mri_rpn_math -out unsorted '$1,-1,1,$1,0.0,>,if_keep,*' %s"%\
                    os.path.join(homedir,inDS.fname))
        else: # Fmap case
            sys.exit("Unexpectedly tried to fold an Fmap!")
    else:
        safeRun("mri_copy_chunk -chunk images -replace %s unsorted"%os.path.join(homedir,inDS.fname))
else:
    if foldFlag:
        if pFlag: # Pmap case
            safeRun("mri_rpn_math -out unsorted '1.0,$1,-,$1,dup,0.5,>,if_keep,1.0,$2,0,==,if_keep' %s %s"%\
                    (os.path.join(homedir,inDS.fname),\
                     os.path.join(homedir,maskDS.fname)))
        elif tFlag: # Tmap case
            safeRun("mri_rpn_math -out unsorted '$1,-1,1,$1,0.0,>,if_keep,*,0.0,$2,0,==,if_keep' %s %s"%\
                    (os.path.join(homedir,inDS.fname),\
                     os.path.join(homedir,maskDS.fname)))
        else: # Fmap case
            sys.exit("Unexpectedly tried to fold an Fmap!")
    else:
        if pFlag: # Pmap case
            safeRun("mri_rpn_math -out unsorted '$1,1.0,$2,0,==,if_keep' %s %s"%\
                    (os.path.join(homedir,inDS.fname),\
                     os.path.join(homedir,maskDS.fname)))
        elif tFlag: # Tmap case
            safeRun("mri_rpn_math -out unsorted '$1,0.0,$2,0,==,if_keep' %s %s"%\
                    (os.path.join(homedir,inDS.fname),\
                     os.path.join(homedir,maskDS.fname)))
        else: # Fmap case
            safeRun("mri_rpn_math -out unsorted '$1,0.0,$2,0,==,if_keep' %s %s"%\
                    (os.path.join(homedir,inDS.fname),\
                     os.path.join(homedir,maskDS.fname)))
            

safeRun("mri_remap -chunk images -order vq -length 1:%d unsorted"%totpix)
if fFlag:
    safeRun("mri_sort unsorted sorted") # get right tail
else:
    safeRun("mri_sort -asc unsorted sorted") # get left tail
safeRun("mri_subset -d q -len 1 -s %d sorted selected"%(nVox-1))
vals= readCmdOutputToList("mri_rpn_math -out junk '0,$1,1,if_print_1' selected")
thresh= string.atof(vals[0])
print(repr(thresh))

# Clean up
os.chdir(homedir)
removeTmpDir(tmpdir)


