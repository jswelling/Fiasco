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
# *     Copyright (c) 1995,1996 Department of Statistics,    *
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
# Notes:
#
import sys
import os
import os.path
import string
import getopt
if "FIASCO" in os.environ:
    sys.path.append(os.environ["FIASCO"])
from fiasco_utils import *

idString= "$Id: physio_correct_triggered.py,v 1.12 2006/03/10 01:23:04 welling Exp $"

def calcBlockDuration(start, end, pulses):
    return pulses[end][0]+2-pulses[start][0]

def printBlockList(title, blocks):
    verboseMessage("%s:"%title)
    blocks.reverse()
    for i in xrange(len(blocks)):
        verboseMessage("  block %d: %d to %d; %d pulses, %d samples"%(i,blocks[i][0],blocks[i][1],blocks[i][2],blocks[i][3]))
    blocks.reverse()

def checkEnvInt(val, name):
    if name in os.environ:
        result= string.atoi(os.environ[name])
    else:
        result= val
    return result

def checkEnvFloat(val, name):
    if name in os.environ:
        result= string.atof(os.environ[name])
    else:
        result= val
    return result

def checkEnvString(val, name):
    if name in os.environ:
        result= os.environ[name]
    else:
        result= val
    return result

##############################
#
# Main
#
##############################

# Set defaults, testing environment for overriding values
syncDS= checkEnvString("data/sync","F_PHYS_SYNC")
syncThreshDS= checkEnvString("data/sync_thresh","F_PHYS_TRIGGER")
cleanSyncThreshDS= checkEnvString("data/clean_sync_thresh","F_PHYS_CLEANTRIGGER")
cardioDS= checkEnvString("data/cardio","F_PHYS_CARDIO")
subsampCardioDS= checkEnvString("data/cardio_subsampled","F_PHYS_SSCARDIO")
respDS= checkEnvString("data/resp","F_PHYS_RESP")
subsampRespDS= checkEnvString("data/resp_subsampled","F_PHYS_SSRESP")
timeDS= checkEnvString("data/time","F_PHYS_TIME")
subsampTimeDS= checkEnvString("data/time_subsampled","F_PHYS_SSTIME")
paramDS= checkEnvString("stat/physio","F_PHYS_PARAM")
sampleChunk= checkEnvString("samples","F_PHYS_SAMPCHUNK")
sampTime= checkEnvFloat(10,"F_PHYS_SAMPTIME")          # sampling interval in milliseconds
#sampTime= checkEnvFloat(2,"F_PHYS_SAMPTIME")           # sampling interval in milliseconds
nDisDAqs= checkEnvInt(7,"F_PHYS_NDISDAQS")             # how many discarded data acquisitions expected?
#nDisDAqs= checkEnvInt(4,"F_PHYS_NDISDAQS")             # how many discarded data acquisitions expected?
syncOffset= checkEnvInt(0,"F_PHYS_CHTRIG")             # which channel is sync?
cardioOffset= checkEnvInt(1, "F_PHYS_CHCARDIO")        # which channel is cardio?
respOffset= checkEnvInt(2, "F_PHYS_CHRESP")            # which channel is respiration?
#syncOffset= checkEnvInt(1,"F_PHYS_CHSYNC")             # which channel is sync?
#cardioOffset= checkEnvInt(6, "F_PHYS_CHCARDIO")        # which channel is cardio?
#respOffset= checkEnvInt(7, "F_PHYS_CHRESP")            # which channel is respiration?
syncThresh= checkEnvFloat(-500.0, "F_PHYS_SYNCTHRESH") # threshold height for finding sync signals

# Check for "-help"
if len(sys.argv)>1:
    if sys.argv[1] == "-help":
	if len(sys.argv)>2:
	    os.system( "scripthelp %s %s"%(sys.argv[0],sys.argv[2]) )
	else:
	    os.system( "scripthelp %s"%sys.argv[0] )
	sys.exit()

try:
    (opts,pargs) = getopt.getopt(sys.argv[1:],"vd",["samptime=","ndisdaq=","syncchannel=","cardiochannel=","respchannel=","syncthresh="])
except:
    errorMessage("%s: Invalid command line parameter" % sys.argv[0])
    describeSelf()
    sys.exit(1)

Message(idString)

# Check calling syntax; parse args
if len(pargs) != 4 and len(pargs) != 3:
    describeSelf()
    sys.exit(1)

for a,b in opts:
    if a=="-v":
        setVerbose(1)
    if a=="-d":
	setDebug(1)
    if a=="--samptime":
        sampTime= string.atoi(b)
    if a=="--ndisdaq":
        nDisDAqs= string.atoi(b)
    if a=="--syncchannel":
        syncOffset= string.atoi(b)
    if a=="--cardiochannel":
        cardioOffset= string.atoi(b)
    if a=="--respchannel":
        respOffset= string.atoi(b)
    if a=="--syncthresh":
        syncThresh= string.atof(b)

inDS= pargs[0]
outDS= pargs[1]
missingDS= pargs[2]
if len(pargs) < 4:
    physioDS= ""
else:
    physioDS= pargs[3]

# If no physio dataset is given, just copy input to output and exit.
if len(physioDS) == 0:
    Message("No physio data given; copying %s to %s"%(inDS,outDS))
    safeRun("mri_copy_dataset %s %s"%(inDS,outDS))
    sys.exit(0)

# Get relevant dimensions
tdim= getDim(inDS,sampleChunk,"t")
zdim= getDim(inDS,sampleChunk,"z")
TR= string.atof(getField(inDS,sampleChunk,"tr"))/1000.0   # in milliseconds
physioTdim= getDim(physioDS,"images","t")
dimstr= getDimStr(inDS,sampleChunk)

#Check reasonableness of input
if dimstr[len(dimstr)-2:] != "zt":
    sys.exit("Input file %s must have rightmost dimensions zt!"%inDS)

# Make up a temporary directory
tmpdir= makeTempDir("tmp_physio")

# Create directories to output datasets if necessary
for thisDS in ( outDS, syncDS, syncThreshDS, cleanSyncThreshDS, cardioDS, subsampCardioDS, respDS, subsampRespDS, timeDS, subsampTimeDS, paramDS ):
    checkPathExists( thisDS )

# Estimate upper limit of sync pulse separation within blocks
blockThresh= int(0.5 + (1.3*TR/(zdim*sampTime)))

verboseMessage("tdim= %d, zdim= %d; product %d; TR= %f"%(tdim,zdim,tdim*zdim,TR))
verboseMessage("total physio samples: %d"%physioTdim)
verboseMessage("Estimated block threshold %d"%blockThresh)

# Strip out the separate parts of the physio signal
safeRun("mri_subset -d v -l 1 -s %d %s %s"%(syncOffset,physioDS,syncDS))
safeRun("mri_subset -d v -l 1 -s %d %s %s"%(cardioOffset,physioDS,cardioDS))
safeRun("mri_subset -d v -l 1 -s %d %s %s"%(respOffset,physioDS,respDS))

# Trigger on the negative-going side of the sync component
safeRun("mri_smooth -d t -smoother_type ddx %s %s/junk "%(syncDS,tmpdir))
safeRun("mri_rpn_math -out %s '$1,%f,>' %s/junk"%(syncThreshDS,syncThresh,tmpdir))

# Read the trigger sample values into a list
syncThreshSamples= map(int,readCmdOutputToList("mri_rpn_math -out %s/junk '0,$1,1,if_print_1' %s"%(tmpdir,syncThreshDS)))

#
# For every non-zero sample, find the number of samples
# between it and the previous non-zero sample.  The entries
# in syncPulses are [indexOfPulse, samplesSinceLastPulse]
#
syncPulses= []
samps= 0
for i in xrange(len(syncThreshSamples)):
    if syncThreshSamples[i] == 0:
        samps = samps + 1
    else:
        syncPulses.append( (i, samps) )
        samps= 0

#
# blocks will be a list of lists, each of which
# is (firstPulse,lastPulse, nPulses, nSamples).  The last
# block in time is first in the list.
#
blocks = []
state= "outside_block"
for i in xrange(len(syncPulses)-1,0,-1):
    if syncPulses[i][1] >= blockThresh:
        if state == "in_block":
            blocks.append( (i, blockBack, blockBack+1-i, calcBlockDuration(i,blockBack,syncPulses)) )
            state= "outside_block"
        else:
            continue
    else:
        if state == "in_block":
            continue
        else:
            blockBack= i;
            state= "in_block"

if state == "in_block":
    blocks.append( (0, blockBack, blockBack+1, calcBlockDuration(0,blockBack,syncPulses)) )
    state= "outside_block"
    
if getVerbose():
    printBlockList("raw blocks", blocks);

#
# Find the relevant (non-DisDAq) part of each block, counting backwards
#
timesFound= 0
trimmedBlocks= []
for thisBlock in blocks:
    timesThisBlock= int( (thisBlock[2] + zdim/2)/zdim ) - nDisDAqs
    if timesThisBlock < 0:
        sys.exit("Tripped over an unexpectedly short data block at physio sample %d!"%thisBlock[0])
    if thisBlock[2] != (timesThisBlock + nDisDAqs)*zdim:
        deficit= (timesThisBlock + nDisDAqs)*zdim - thisBlock[2]
        errorMessage("Warning: block beginning at sample %d is missing %d trigger pulses!"%(thisBlock[0],deficit))
    lastPulse= thisBlock[1]
    firstPulse= lastPulse + 1 - timesThisBlock*zdim
    timesFound = timesFound + timesThisBlock
    trimmedBlocks.append( (firstPulse, lastPulse, timesThisBlock*zdim, calcBlockDuration(firstPulse,lastPulse,syncPulses)) )
    if timesFound >= tdim:
        break

if getVerbose():
    printBlockList("trimmed blocks",trimmedBlocks)

#
# Construct a list of cleaned-up sync markers
#
cleanSamples= []
for x in syncThreshSamples:
    cleanSamples.append(0)
for block in trimmedBlocks:
    for x in xrange( block[0],block[1]+1 ):
        cleanSamples[syncPulses[x][0]]= 1
writeListToCmdInput(cleanSamples,"mri_from_ascii -order t -len %d %s"%(len(cleanSamples),cleanSyncThreshDS))

#
# Construct time dataset
#
time= []
for x in xrange(len(cleanSamples)):
    time.append(x)
writeListToCmdInput(time,"mri_from_ascii -order t -len %d %s"%(len(cleanSamples),timeDS))

#
# Subsample and fold the time series datasets
#
for thisDS in [ cardioDS, respDS, timeDS ]:
    safeRun("mri_rpn_math -out %s/junk '0,$1,$2,if_print_1' %s %s | mri_from_ascii -order t -len %d %s/%s"%(tmpdir,thisDS,cleanSyncThreshDS,zdim*tdim,tmpdir,os.path.basename(thisDS)))

reorderPattern= getField(inDS,sampleChunk,"reorder_pattern")
for thisDS in [ cardioDS, respDS, timeDS ]:
    safeRun("mri_scan_fold -zdm %d -reorder %s %s/%s %s"%\
            (zdim,reorderPattern,tmpdir,os.path.basename(thisDS),thisDS))

#
# Regress the physio signal out of the input dataset, adding
# the constant term back in.
#
if dimstr[0] != 'v':
    dimstr= "v"+dimstr
otherDims= string.replace(dimstr,"t","")
otherDims= string.replace(otherDims,"v","")
permDims= "vt"+otherDims
debugMessage("Dimensions <%s>, permuted dimensions <%s>"%(dimstr,permDims))
safeRun("mri_permute -c %s -order %s %s %s_p"%(sampleChunk,permDims,inDS,inDS))
safeRun("mri_copy_chunk -chunk missing -replace %s %s_p"%(missingDS,inDS))
for thisDS in [ cardioDS, respDS, timeDS ]:
    safeRun("mri_permute -order tz %s %s_p"%(thisDS,thisDS))
safeRun("mri_glm -v -var -ssqr -output %s_p -est %s:%s %s_p:%s %s_p:images %s_p:images %s_p:images"%(outDS,paramDS,sampleChunk,inDS,sampleChunk,respDS,cardioDS,timeDS))
for thisDS in [ inDS, cardioDS, respDS, timeDS ]:
    safeRun("mri_destroy_dataset %s_p"%thisDS)
safeRun("mri_subset -d v -s 0 -l 2 %s %s/phys_const"%(paramDS,tmpdir))
safeRun("mri_remap -order %s -c %s %s/phys_const"%(dimstr,sampleChunk,tmpdir))
safeRun("mri_permute -order %s -c %s %s_p %s_noconst"%(dimstr,sampleChunk,outDS,outDS))
safeRun("mri_destroy_dataset %s_p"%outDS)
safeRun("mri_rpn_math -out %s -c %s '$1,$2,+' %s_noconst %s/phys_const"%(outDS,sampleChunk,outDS,tmpdir))
safeRun("mri_destroy_dataset %s_noconst"%outDS)

# Clean up
removeTmpDir(tmpdir)


