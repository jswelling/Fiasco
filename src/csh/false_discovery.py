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
import string
import getopt

idString= "$Id: false_discovery.py,v 1.14 2006/05/04 23:02:58 welling Exp $"

def safeRun(cmd):
    cmdout= os.popen(cmd)
    if cmdout.close() != None :
	print("Command failed: <%s>"%cmd)
	sys.exit(1)

def removeTmpDir(path):
    for file in os.listdir(path):
	os.remove("%s/%s"%(path,file));
    os.rmdir(path);

def describeSelf():
    os.system("scripthelp %s usage"%sys.argv[0]);

def getDim(mrifile,index):
    cmdout= os.popen("mri_printfield -field images.extent.%s %s" % (index,mrifile))
    xstr= cmdout.read()
    if cmdout.close() != None :
	print("mri_printfield failed for images.extent.%s on %s!"%(index,mrifile))
	sys.exit(1)
    return string.atoi(string.strip(xstr));

def getDimStr(mrifile):
    cmdout= os.popen("mri_printfield -field images.dimensions %s" % mrifile)
    dimstr= cmdout.read()
    if cmdout.close() != None :
	print("mri_printfield failed for images.dimensions on %s!"%mrifile)
	sys.exit(1)
    return string.strip(dimstr);

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

def cFunc(totpix, samplesIndependent):
    if samplesIndependent:
	return 1;
    else:
	result= 0;
	for i in xrange(1,totpix):
	    result= result + (1.0/i)
	return result;

##############################
#
# Main
#
##############################

samplesIndependent= 0;

# Check for "-help"
if len(sys.argv)>1:
    if sys.argv[1] == "-help":
	if len(sys.argv)>2:
	    os.system( "scripthelp %s %s"%(sys.argv[0],sys.argv[2]) );
	else:
	    os.system( "scripthelp %s"%sys.argv[0] );
	sys.exit();

try:
    (opts,pargs) = getopt.getopt(sys.argv[1:],"",["independent"])
except:
    print("%s: Invalid command line parameter" % sys.argv[0])
    describeSelf();
    sys.exit()

#Check calling syntax; parse args
if len(pargs) != 2 :
    describeSelf()
    sys.exit(1)

for a,b in opts:
    if a=="--independent":
	samplesIndependent= 1;

qrate= string.atof(pargs[0])
infile= pargs[1]

#Make up a temporary directory
if "F_TEMP" in os.environ:
    tmpdir= "%s/tmp_fd_%d"%(os.environ["F_TEMP"],os.getpid());
else:
    tmpdir= "./tmp_fd_%d"%os.getpid();
#print("temp directory is %s"%tmpdir;)
os.makedirs(tmpdir);

#Get relevant dimensions
xdim= getDim(infile,"x");
ydim= getDim(infile,"y");
zdim= getDim(infile,"z");
totpix= xdim*ydim*zdim
dimstr= getDimStr(infile);

#Check reasonableness of input
if dimstr != "xyz":
    if dimstr == "xyzt":
	if getDim(infile,"t") != 1:
	    print("Input file must have t extent 1!")
	    sys.exit(1)
    elif dimstr == "vxyzt":
	if getDim(infile,"t") != 1:
	    print("Input file must have t extent 1!")
	    sys.exit(1)
	if getDim(infile,"v") != 1:
	    print("Input file must have v extent 1!")
	    sys.exit(1)
    elif dimstr == "vxyz":
	if getDim(infile,"v") != 1:
	    print("Input file must have v extent 1!")
	    sys.exit(1)
    else:
        print("Input file must have dimensions (v)xyz(t)!")
        sys.exit(1)

#Generate sorted file
newdimstr= "vq"
newextstr= "1:%d"%totpix
safeRun("mri_copy_chunk -chunk images -replace %s %s/tmp_fd"%\
	(infile,tmpdir))
safeRun("mri_remap -chunk images -order %s -length %s %s/tmp_fd "%\
	(newdimstr,newextstr,tmpdir))
safeRun("mri_sort -asc %s/tmp_fd %s/tmp_fd_2"%(tmpdir,tmpdir))

#Walk through the file of sorted p scores, finding the first entry to
#exceed the cutoff.  We drop all NaN and Inf values at this point
cmdout= os.popen("mri_rpn_math -out %s/tmp_fd_3 '$1,dup,is_finite,if_print_1,0.0' %s/tmp_fd_2"%(tmpdir,tmpdir))
allLines= cmdout.readlines()
if cmdout.close() != None :
    print("mri_rpn_math failed!")
    sys.exit(1)
validpix= len(allLines);
#print("totpix %d; validpix %d"%(totpix,validpix))
i= 1
scale= qrate/cFunc(validpix,samplesIndependent);
#print("%d: scale is %f "%(samplesIndependent,scale);)
rejectNullThresh= 0.0
for line in allLines:
    try:
	val= string.atof(string.strip(line))
    except ValueError:
	continue
    cutoff= (scale*i)/validpix

#    print("%d: %g vs %g"%(i,cutoff,val))
    if cutoff >= val:
	rejectNullThresh= val

    i= i+1

# We use repr() to get more significant digits in output.
print(repr(rejectNullThresh))

# Clean up
removeTmpDir(tmpdir);


