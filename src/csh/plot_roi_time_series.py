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
import threading
import subprocess
from subprocess import PIPE
if "FIASCO" in os.environ:
    sys.path.append(os.environ["FIASCO"])
from fiasco_utils import *

idString= "$Id: plot_roi_time_series.py,v 1.4 2007/07/20 15:46:22 welling Exp $"
runningsumWindow= 6
autoBarLocList= [6.5, 23, 29.5, 46, 52.5, 69, 75.5, 92, 98.5, 115, 121.5,
                 138, 143, 149.5, 166, 172.5, 189, 195.5, 212, 218.5, 235,
                 241.5, 258, 264.5, 281, 286, 292.5, 309, 315.5, 332, 338.5,
                 355, 361.5, 378, 384.5, 401, 407.5, 424]
barLocList= None

def parseBarLocations(txt):
    print("Parsing <%s>"%txt)
    substrings= txt.split(",")
    return [ float(x) for x in substrings ]

def slurpThisPipe(p,lines):
    while 1:
        line= p.readline()
        if len(line)==0: break
        lines.append(line)

def runChildProcess(rCommand, voxValListList, rsList, bars, plotName,
                    plotTitle):
    stdoutLines= []
    stderrLines= []
    debugMessage("Starting <%s>"%rCommand)
    childHook= subprocess.Popen(rCommand, shell=True,
                                stdin=PIPE, stdout=PIPE, stderr=PIPE,
                                close_fds=True)
    stdoutThread= threading.Thread(group=None,target=slurpThisPipe,\
                                   args=(childHook.fromchild,stdoutLines))
    stdoutThread.start()
    stderrThread= threading.Thread(group=None,target=slurpThisPipe,\
                                   args=(childHook.childerr,stderrLines))
    stderrThread.start()
    toChild= childHook.tochild
    debugMessage("Child process started.")

    toChild.write('postscript(horizontal=T,file="%s");\n'%plotName)
    toChild.write('par(mfrow=c(2,1));\n')

    # We need max and min
    min= max= rsList[0]
    for val in rsList:
        if val<min: min=val
        if val>max: max=val
    for valList in voxValListList:
        for val in valList:
            sval= val
            if sval<min: min=sval
            if sval>max: max=sval
    
    first= 1
    toChild.write("rsdata <- c(\n")
    for val in rsList:
        if first:
            toChild.write("%g\n"%val)
            first= 0
        else:
            toChild.write(",%g\n"%val)
    toChild.write(");\n")

    # Top plot
    subTitle="%d step running sum of mean"%runningsumWindow
    toChild.write('plot(rsdata,type="l",xlab="Image number",ylab="Intensity",main="%s",sub="%s",ylim=c(%g,%g));\n'%(plotTitle,subTitle,min,max))
    if bars != None:
        for val in bars:
            toChild.write("lines(c(%g,%g),c(%g,%g),lty=2);\n"%\
                          (val,val,min,max))

    # Bottom plot
    for i in xrange(len(voxValListList)):
        toChild.write("curve <- c(\n");
        first= 1
        for val in voxValListList[i]:
            if first:
                toChild.write("%g\n"%val)
                first= 0
            else:
                toChild.write(",%g\n"%val)
        toChild.write(");\n")
        if i==0:
            subTitle= "%d contributing voxels"%len(voxValListList)
            toChild.write('plot(curve,type="l",lty=%d,xlab="Image number",ylab="Intensity",sub="%s",ylim=c(%g,%g));\n'%\
                          ((i+1),subTitle,min,max))
        else:
            toChild.write("lines(curve,lty=%d);\n"%(i+1))
    if bars != None:
        for val in bars:
            toChild.write("lines(c(%g,%g),c(%g,%g),lty=2);\n"%\
                          (val,val,min,max))
    toChild.close()
    status= childHook.wait()
    debugMessage("Child exited with status %d"%status)
    stdoutThread.join()
    stderrThread.join()
    if len(stdoutLines)!=0 and len(stderrLines)==0 and getDebug():
        errorMessage("These lines were written to stdout:")
        for line in stdoutLines:
            errorMessage("<%s>"%line)
    if len(stderrLines)!=0:
        errorMessage("R wrote the following lines to stderr!")
        for line in stderrLines:
            errorMessage("<%s>"%line)
        if len(stdoutLines)!=0:
            errorMessage("These lines were written to stdout:")
            for line in stdoutLines:
                errorMessage("<%s>"%line)
        sys.exit("R shouldn't write anything to stderr.")
    if status != 0:
        sys.exit("R exited abnormally, with status %d!"%status)

##############################
#
# Main
#
##############################

# Figure out what slave to use
rCommand= "R --slave --no-save --no-restore-history --no-restore-data"
#rCommand= "R --no-save --no-restore-history --no-restore-data"
#rCommand= "cat"
if "SPLUS" in os.environ:
    rCommand= os.environ['SPLUS']
if rCommand.find('R') >= 0:
    if rCommand.find('--no-save')<0: rCommand += ' --no-save'
    if rCommand.find('--slave')<0 and not getDebug(): rCommand += ' --slave'
elif rCommand.find('Splus')>=0:
    sys.exit("%s currently doesn't work with Splus- try R!"%\
             os.path.basename(sys.argv[0]))

# Check for "-help"
if len(sys.argv)>1:
    if sys.argv[1] == "-help":
        if len(sys.argv)>2:
            os.system( "scripthelp %s %s"%(sys.argv[0],sys.argv[2]) );
        else:
            os.system( "scripthelp %s"%sys.argv[0] );
        sys.exit();

try:
    (opts,pargs) = getopt.getopt(sys.argv[1:],"vdlr:",["out=","title=",
                                                       "autobars","bars=",
                                                       "window="])
except:
    print("%s: Invalid command line parameter" % sys.argv[0])
    describeSelf();
    sys.exit()

# These default thresholds are against a typical axial anatomical scan,
# which has values in the range 0 < v < 1000 or so.
roiMaskName= None
plotName= "plot.ps"
plotTitle= ""
for a,b in opts:
    if a=="-v":
        setVerbose(1)
    if a=="-d":
        setDebug(1)
    if a=="-r":
        roiMaskName= os.path.abspath(b)
    if a=="--out":
        plotName= b
    if a=="--title":
        plotTitle= b
    if a=="--autobars":
        barLocList= autoBarLocList
    if a=="--bars":
        barLocList= parseBarLocations(b)
    if a=="--window":
        runningsumWindow= int(b)

if roiMaskName==None:
    sys.exit("Required ROI mask dataset not specified")

if len(pargs) != 1:
    sys.exit("Failed to specify which dataset to sample")
sampleThis= os.path.abspath(pargs[0])

# What have we got?
roiDS= MRIDataset(roiMaskName)
roiChunk= roiDS.getChunk('images')
sampleDS= MRIDataset(sampleThis)
sampleChunk= sampleDS.getChunk('images')
roiDims= roiChunk.getValue('dimensions')
if roiDims.find('t')>=0:
    roiTDim= roiChunk.getDim('t')
else:
    roiTDim= 1
if roiDims[-1]=='t':
    roiDims=roiDims[:-1]
tDim= sampleChunk.getDim('t')
qlen= 1
origExtents= ""
for dim in roiDims:
    ext= roiChunk.getDim(dim)
    qlen *= ext
    origExtents= "%s:%s"%(origExtents,ext)
origExtents= origExtents[1:]
    
# Make a temp directory and go there
tmpdir= makeTempDir('tmp_plot_roi_time_series')
homedir= os.getcwd()
os.chdir(tmpdir)

# Create a missing-free version of the data
safeRun("smooth_missing -i %s -headerout data_p"%sampleThis)

# Remap it to match the mask
safeRun("mri_remap -order tq -len %d:%d data_p"%(tDim,qlen))

# Get a list of the q's in the mask
safeRun("mri_rpn_math -out boolmask '$1,0,!=' %s"%roiMaskName)
safeRun("mri_remap -order qt -len %d:%d boolmask"%(qlen,roiTDim))
lines= readCmdOutputToList("mri_rpn_math '$q,$1,if_print_1' boolmask")
voxList= []
for line in lines:
    voxList.append(int(line.strip()))
nVoxels= len(voxList)

voxValListList= []
for vox in voxList:
    voxValList= []
    safeRun("mri_subset -d q -len 1 -s %d data_p onevox"%vox)
    lines= readCmdOutputToList("mri_rpn_math '$1,1,if_print_1' onevox")
    for line in lines:
        voxValList.append(float(line.strip()))
    voxValListList.append(voxValList)

# Create the running sum list
meanList= []
for i in xrange(tDim):
    val= 0.0
    for list in voxValListList:
        val += list[i]
    meanList.append(val/float(nVoxels))
rsList= []
val= 0.0
for i in xrange(tDim):
    val += meanList[i]
    if i>=runningsumWindow:
        val -= meanList[i-runningsumWindow]
    rsList.append(val)
for i in xrange(runningsumWindow):
    rsList[i] /= float(i+1)
for i in xrange(runningsumWindow,tDim):
    rsList[i] /= float(runningsumWindow)
## for i in xrange(tDim):
##     print("%d: vals %s -> mean %g -> running sum %g"%\
##           (i,voxValListList[0][i],meanList[i],rsList[i]))

# Create the plot
runChildProcess(rCommand, voxValListList, rsList, barLocList, "test.ps",
                plotTitle)
safeRun("mv test.ps %s"%os.path.join(homedir,plotName))
        
# Clean up
os.chdir(homedir)
if not getDebug():
    removeTmpDir(tmpdir)

