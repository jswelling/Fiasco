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
from __future__ import print_function
import sys
import os
import os.path
import string
import getopt
import threading
import subprocess
import types
if 'FIASCO' in os.environ:
    sys.path.append(os.environ["FIASCO"])
from fiasco_utils import *
import customize_local

idString= "$Id: plot_roi_event_response.py,v 1.6 2007/08/16 20:44:45 welling Exp $"

# These are the defaults for some command line options
eventIdentifierPhrase= "Stim"
eventTimeToLive= 10.0   # seconds
timeBinsPerSec= 2
plotName= "plot.ps"

barLocList= None

def parseBarLocations(txt):
    verboseMessage("Parsing bar location string <%s>"%txt)
    substrings= txt.split(",")
    return [ float(x) for x in substrings ]

def slurpThisPipe(p,lines):
    while 1:
        line= p.readline()
        if len(line)==0: break
        lines.append(line)

def runChildProcess(rCommand, allEventTupleDict, trendlineTupleDict,
                    bars, plotName, plotTitle, plotSubTitle):
    stdoutLines= []
    stderrLines= []
    clrList= ["black","red","blue","green"]
    debugMessage("Starting <%s>"%rCommand)
    p= subprocess.Popen(rCommand, shell=True)
    childHook= (p.stdout, p.stdin, p.stderr)
    stdoutThread= threading.Thread(group=None,target=slurpThisPipe,\
                                   args=(childHook.fromchild,stdoutLines))
    stdoutThread.start()
    stderrThread= threading.Thread(group=None,target=slurpThisPipe,\
                                   args=(childHook.childerr,stderrLines))
    stderrThread.start()
    toChild= childHook.tochild

    if (allEventTupleDict != None) \
           and (allEventTupleDict.keys() != trendlineTupleDict.keys()):
        exit("The two coordinate lists given have different keys!")

    debugMessage("Child process started.")

    toChild.write('nowhere= file("/dev/null","w");\n')
    toChild.write('sink(nowhere,type="message");\n')
    toChild.write('require(gplots,warn.conflicts=FALSE,quietly=TRUE);\n')
    toChild.write('sink(NULL,type="message");\n')
    toChild.write('postscript(horizontal=T,file="%s");\n'%plotName)

    # If the Y coordinates of the trend lines are 2-component tuples,
    # it's a sign that we are in IQR mode and are supposed to draw
    # a band rather than a single trend line.
    (x,y)= trendlineTupleDict[trendlineTupleDict.keys()[0]][0]
    iqrMode= ( type(y)==types.TupleType )

    # We need max and min
    min= max= None
    for name in trendlineTupleDict.keys():
        tupleList= trendlineTupleDict[name]
        if iqrMode:
            if not max or not min:
                (x,ytuple)= tupleList[0]
                (y1,y2)= ytuple
                min= max= y1
            for (x,ytuple) in tupleList:
                (y1,y2)= ytuple
                if min>y1: min= y1
                if max<y1: max= y1
                if min>y2: min= y2
                if max<y2: max= y2
        else:
            if not max or not min:
                (x,y)= tupleList[0]
                min= max= y
            for (x,y) in tupleList:
                if min>y: min= y
                if max<y: max= y
    if allEventTupleDict != None:
        for name in allEventTupleDict.keys():
            tupleList= allEventTupleDict[name]
            for (x,y) in tupleList:
                if min>y: min= y
                if max<y: max= y

    names= trendlineTupleDict.keys()
    if iqrMode:
        toChild.write('nameList <- c("%s IQR"\n'%names[0])
        for n in names[1:]:
            toChild.write(',"%s IQR"\n'%n)
        toChild.write(');\n')
    else:
        toChild.write('nameList <- c("%s"\n'%names[0])
        for n in names[1:]:
            toChild.write(',"%s"\n'%n)
        toChild.write(');\n')
    
    for i in xrange(len(trendlineTupleDict.keys())):
        if i==0: toChild.write('clrList <- c("%s"\n'%clrList[i])
        else: toChild.write(',"%s"\n'%clrList[i%len(clrList)])
    toChild.write(');\n')
                            
    counter= 1
    for name in trendlineTupleDict.keys():
        tupleList= trendlineTupleDict[name]
        first= 1
        toChild.write("xvals <- c(\n")
        for (x,y) in tupleList:
            if first:
                toChild.write("%g\n"%x)
                first= 0
            else:
                toChild.write(",%g\n"%x)
        toChild.write(");\n")
        if iqrMode:
            first= 1
            toChild.write("y1vals <- c(\n")
            for (x,ytuple) in tupleList:
                (y1,y2)= ytuple
                if first:
                    toChild.write("%g\n"%y1)
                    first= 0
                else:
                    toChild.write(",%g\n"%y1)
            toChild.write(");\n")
            first= 1
            toChild.write("y2vals <- c(\n")
            for (x,ytuple) in tupleList:
                (y1,y2)= ytuple
                if first:
                    toChild.write("%g\n"%y2)
                    first= 0
                else:
                    toChild.write(",%g\n"%y2)
            toChild.write(");\n")
            if counter==1:
                toChild.write('plot(xvals,y1vals,col="%s",type="l",xlab="time(sec)",ylab="scaled intensity",main="%s",sub="%s",ylim=c(%g,%g));\n'%\
                              (clrList[(counter-1)%len(clrList)],
                               plotTitle,plotSubTitle,min,max))
                toChild.write('lines(xvals,y2vals,col="%s",lty=%d);\n'%\
                              (clrList[(counter-1)%len(clrList)],counter))
            else:
                toChild.write('lines(xvals,y1vals,col="%s",lty=%d);\n'%\
                              (clrList[(counter-1)%len(clrList)],counter))
                toChild.write('lines(xvals,y2vals,col="%s",lty=%d);\n'%\
                              (clrList[(counter-1)%len(clrList)],counter))
        else:
            first= 1
            toChild.write("yvals <- c(\n")
            for (x,y) in tupleList:
                if first:
                    toChild.write("%g\n"%y)
                    first= 0
                else:
                    toChild.write(",%g\n"%y)
            toChild.write(");\n")
            if counter==1:
                toChild.write('plot(xvals,yvals,col="%s",type="l",xlab="time(sec)",ylab="scaled intensity",main="%s",sub="%s",ylim=c(%g,%g));\n'%\
                              (clrList[(counter-1)%len(clrList)],
                               plotTitle,plotSubTitle,min,max))
            else:
                toChild.write('lines(xvals,yvals,col="%s",lty=%d);\n'%\
                              (clrList[(counter-1)%len(clrList)],counter))

        if allEventTupleDict!=None:
            tupleList= allEventTupleDict[name]
            first= 1
            toChild.write("xvals <- c(\n")
            for (x,y) in tupleList:
                if first:
                    toChild.write("%g\n"%x)
                    first= 0
                else:
                    toChild.write(",%g\n"%x)
            toChild.write(");\n")
            first= 1
            toChild.write("yvals <- c(\n")
            for (x,y) in tupleList:
                if first:
                    toChild.write("%g\n"%y)
                    first= 0
                else:
                    toChild.write(",%g\n"%y)
            toChild.write(");\n")
            toChild.write('points(xvals,yvals,col="%s",pch=%d);\n'%\
                          (clrList[(counter-1)%len(clrList)],counter-1))

        counter += 1

    if allEventTupleDict != None:
        toChild.write('smartlegend(x="right",y="top",nameList,col=clrList,lty=1:%d,pch=0:%d);\n'%\
                      ((counter-1),(counter-2)))
    else:
        toChild.write('smartlegend(x="right",y="top",nameList,col=clrList,lty=1:%d);\n'%\
                      (counter-1))

    
    if bars != None:
        for val in bars:
            toChild.write("lines(c(%g,%g),c(%g,%g),lty=2);\n"%\
                          (val,val,min,max))

    toChild.close()
    status= 0
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

def median( vlist ):
    vlist.sort()
    return vlist[len(vlist)/2]

def q1q3( vlist ):
    # Returns a tuple
    vlist.sort()
    return (vlist[len(vlist)/4],vlist[(3*len(vlist))/4],)

def mean( vlist ):
    sum= 0.0
    for v in vlist: sum += v
    sum /= float(len(vlist))
    return sum

def unpackState( tuple ):
    # This routine unpackes a tuple entry as found in a state list
    if len(tuple)==3:
        (name,dur,start)= tuple
    else:
        (name,dur,start,comment)= tuple
    return (name,start,start+dur)

def rescale(vlist):
    vals= [y for (x,y) in vlist]
    m= mean(vals)
    if m != 0.0:
        return [(x,y/m) for (x,y) in vlist]
    else: return vlist

def updateEventTupleLists(liveEventTupleList,allEventTupleDict,state):
    (name,start,end)= unpackState(state)
    name= name.strip()
    if name.find(eventIdentifierPhrase)>=0:
        debugMessage("Hit trigger event %s starting at time %g"%(name,start))
        if not allEventTupleDict.has_key(name):
            allEventTupleDict[name]= []
        liveEventTupleList.append((name,start,[]))

def cullLiveEventTupleList(liveEventTupleList,timeToLive,time):
    for t in liveEventTupleList:
        (name,startTime,sampList)= t
        if startTime+timeToLive <= time:
            debugMessage("Retiring %s from live events list"%name)
            # The collection of samples is rescaled and appended to
            # the full list of samples associated with the event
            allEventTupleDict[name] += rescale( sampList )
            liveEventTupleList.remove(t)

def recordSampleInEventTupleLists(liveEventTupleList, val, time):
    for i in xrange(len(liveEventTupleList)):
        t= liveEventTupleList[i]
        (name,start,sampList)= t
##         print "State %s sees value %g at time %g - %g = %g"%\
##               (name,val,time,start,time-start)
        sampList.append( ((time-start), val) )

def binTime( time ):
    return round(timeBinsPerSec*time)/timeBinsPerSec
    
def writeForGnuplot( format, name, tList ):
    fname= os.path.join(homedir,format%name)
    file= open(fname,"w")
    for t in tList:
        file.write("%g %g\n"%t)
    file.close()
    verboseMessage("Wrote %s (%d entries)"%(fname,len(tList)))
    

##############################
#
# Main
#
##############################

# Figure out what slave to use
rCommand= "R --slave --no-save --no-restore-history --no-restore-data"
#rCommand= "R --no-save --no-restore-history --no-restore-data"
#rCommand= "cat"
if 'SPLUS' in os.environ:
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
                                                       "bars=","bps=",
                                                       "ttl=","phrase=",
                                                       "iqr","nosamp"])
except:
    print("%s: Invalid command line parameter" % sys.argv[0])
    describeSelf();
    sys.exit()

roiMaskName= None
iqrMode= 0
nosampMode= 0
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
    if a=="--bars":
        barLocList= parseBarLocations(b)
    if a=="--bps":
        timeBinsPerSec= float(b)
    if a=="--ttl":
        eventTimeToLive= float(b)
    if a=="--phrase":
        eventIdentifierPhrase= b
    if a=="--iqr":
        iqrMode= 1
    if a=="--nosamp":
        nosampMode= 1

if roiMaskName==None:
    sys.exit("Required ROI mask dataset not specified")

sampleThis= None
transDict= {}
for arg in pargs:
    if string.find(arg,"=") >= 0:
        (key,val) = map(string.strip,string.split(arg,"="))
        transDict[key]= val
    else:
        if sampleThis == None:
            sampleThis= os.path.abspath(arg)
        else:
            sys.exit("%s: command line contains repeated input filenames"
                     %sys.argv[0])

if sampleThis==None:
    sys.exit("Failed to specify which dataset to sample")

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

# Populate the dictionaries full of methods that are used to
# communicate with customize_local. transDict provides info
# to get the process started.
verboseMessage("Seeking stimulus sequence data...")
(dbDict, envDict, catDict, hookDict)= \
         customize_local.populateMethodTables(transDict)

#
# The buildStateList method takes one argument, the integer acq number.
# A list of states is returned.  Each entry is a tuple of length 3 or
# 4.  If length is 3, the tuple is (stateName,stateDuration,stateStart).
# If length is 4, a comment string is also present.
#
buildStateList= hookDict['BuildStateList']
if buildStateList==None:
    sys.exit("customize_local knows no stimulus timing method for this subject, run and task")
verboseMessage("...done")
debugMessage("buildStateList method is %s"%buildStateList)

# Make a temp directory and go there
tmpdir= makeTempDir('tmp_plot_roi_event_response')
homedir= os.getcwd()
os.chdir(tmpdir)

# Create a missing-free version of the data
verboseMessage("Smoothing missing data...")
safeRun("smooth_missing -i %s -headerout data_p"%sampleThis)
verboseMessage("...done")

# Remap it to match the mask
safeRun("mri_remap -order tq -len %d:%d data_p"%(tDim,qlen))

# Get a list of the q's in the mask
verboseMessage("Identifying mask voxels...")
safeRun("mri_rpn_math -out boolmask '$1,0,!=' %s"%roiMaskName)
safeRun("mri_remap -order qt -len %d:%d boolmask"%(qlen,roiTDim))
lines= readCmdOutputToList("mri_rpn_math '$q,$1,if_print_1' boolmask")
voxList= []
for line in lines:
    voxList.append(int(line.strip()))
nVoxels= len(voxList)
verboseMessage("...done (%d voxels)"%nVoxels)

verboseMessage("Collecting sample values...")
voxValListList= []
for vox in voxList:
    voxValList= []
    safeRun("mri_subset -d q -len 1 -s %d data_p onevox"%vox)
    lines= readCmdOutputToList("mri_rpn_math '$1,1,if_print_1' onevox")
    for line in lines:
        voxValList.append(float(line.strip()))
    voxValListList.append(voxValList)
verboseMessage("...done")

# Invert the list to produce a dict of lists, indexed by time
samples= {}
for list in voxValListList:
    for i in xrange(tDim):
        if not samples.has_key(i):
            samples[i]= []
        samples[i].append( list[i] )

# Produce one useful summary value per time
verboseMessage("Producing summary values...")
summaryValList= []
for i in xrange(tDim):
    summaryValList.append( median(samples[i]) )
verboseMessage("...done")
    
# Get a list of acquisition numbers
if sampleDS.hasChunk("acq_blocks"):
    acqList= readCmdOutputToList("mri_rpn_math -c acq_blocks '$1,1,if_print_1' %s"%\
                        sampleThis)
    acqList= [int(s) for s in acqList]
else:
    acqList= tDim*[0]
acqOffsetList= [0]
base= 0
for i in xrange(1,tDim):
    if acqList[i]!=acqList[i-1]:
        base= i
    acqOffsetList.append(i-base)

# For each timepoint, find the time since the last interesting timepoint.
# Interesting timepoints are those that correspond to the beginning of
# a trial.
verboseMessage("Associating samples with events...")
TR= customize_local.getTR()
stateList= None
currentState= None
currentStateStart= None
currentStateEnd= None
currentStateName= None
liveEventTupleList= [] # format is [(name,startTime,sampleList), ...]
allEventTupleDict= {}  # format is { name:sampleList, ...}
currentAcq= -1 # impossible value
for i in xrange(1,tDim):
    acq= acqList[i]
    if acq != currentAcq:
        stateList= buildStateList(acq)
        (currentStateName, currentStateStart, currentStateEnd)= \
                          unpackState(stateList[0])
        updateEventTupleLists(liveEventTupleList,allEventTupleDict,
                              stateList[0])
        currentAcq= acq
    acqOffset= acqOffsetList[i]
    # assume times shifted to half-way through the image acquisition
    timeWithinAcq= TR*(acqOffset+0.5)
    cullLiveEventTupleList(liveEventTupleList,eventTimeToLive,timeWithinAcq)
    # Skip forward to the current state
    while timeWithinAcq>currentStateEnd:
        if len(stateList)>1:
            stateList= stateList[1:]
            (currentStateName, currentStateStart, currentStateEnd)= \
                               unpackState(stateList[0])
        else:
            # Overrunning end of stateList- bad experimental design!
            currentStateName= 'NA'
            currentStateStart= currentStateEnd
            currentStateEnd += 100.0
        updateEventTupleLists(liveEventTupleList,allEventTupleDict,
                              stateList[0])
    recordSampleInEventTupleLists(liveEventTupleList,summaryValList[i],
                                  timeWithinAcq)
    debugMessage("t= %g: state %s (%g to %g)"%\
                   (timeWithinAcq,currentStateName,
                    currentStateStart,currentStateEnd))
verboseMessage("...done")

## for name in allEventTupleDict.keys():
##     sampTupleList= allEventTupleDict[name]
##     writeForGnuplot( "samples_%s.ascii", name, sampTupleList )

# Group samples into bins, and make a trend line from summary statistics
# of each bin.  In IQR mode, the summary statistics points include tuples
# of Q1 and Q3 rather than single values.
verboseMessage("Binning samples...")
trendlineTupleDict= {}
for name in allEventTupleDict.keys():
    sampTupleList= allEventTupleDict[name]
    binDict= {}
    for (x,y) in sampTupleList:
        binnedTime= binTime(x)
        if not binDict.has_key(binnedTime):
            binDict[binnedTime]= []
        binDict[binnedTime].append(y)
    trendTupleList= []
    timeList= binDict.keys()
    timeList.sort()
    for time in timeList:
        debugMessage("%s bin %g: %d samples"%(name,time,len(binDict[time])))
        if iqrMode:
            trendTupleList.append((time,q1q3(binDict[time])))
        else:
            trendTupleList.append((time,median(binDict[time])))
    trendlineTupleDict[name]= trendTupleList
verboseMessage("...done")

## for name in trendlineTupleDict.keys():
##     pointTupleList= trendlineTupleDict[name]
##     writeForGnuplot( "trend_%s.ascii", name, pointTupleList )

# Create the plot
subTitle= "%d voxels, %d samples, %d bins, event lifetime %3.2g sec"%\
          (nVoxels,len(summaryValList),len(binDict.keys()),eventTimeToLive)
if nosampMode:
    runChildProcess(rCommand, None, trendlineTupleDict,
                    barLocList, "test.ps", plotTitle,subTitle)
else:
    runChildProcess(rCommand, allEventTupleDict, trendlineTupleDict,
                    barLocList, "test.ps", plotTitle,subTitle)
safeRun("mv test.ps %s"%os.path.join(homedir,plotName))
        
# Clean up
os.chdir(homedir)
if not getDebug():
    removeTmpDir(tmpdir)

