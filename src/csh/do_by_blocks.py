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
#  - Needs documentation in scripthelp!
#
#
import sys
import os
import os.path
import string
import getopt
import re
if os.environ.has_key("FIASCO"):
    sys.path.append(os.environ["FIASCO"])
from fiasco_utils import *

idString= "$Id: do_by_blocks.py,v 1.9 2007/07/10 16:55:29 welling Exp $"

zdim= None
tdim= None
paramRegex= re.compile('.*\.par\.(\d+)')
datRegex= re.compile('.*\.dat')
mriRegex= re.compile('.*\.mri')
summaryInputRegex= re.compile('.*summary\.input')
summaryMissingRegex= re.compile('.*summary\.missing')
anyRegex= re.compile('.*')
summary_missing_file= None
summary_missing_step= None
first_created_mri= None
parDict= {}
parFormatDict= {}

def parZFastestMethod( fname, parList ):
    # This one is easy- just dump in order
    ofile= file(fname,'a')
    for parLine in parList:
        for tok in parLine:
            ofile.write(tok)
        ofile.write("\n")
    ofile.close()

def parIndexTZMethod( fname, parList ):
    global tdim, zdim
    slices= []
    for i in xrange(zdim):
        slices.append([])
    for parLine in parList:
        slices[int(parLine[1])].append(parLine[2:])
    ofile= file(fname,'a')
    for z in xrange(zdim):
        t= 0
        for parLine in slices[z]:
            ofile.write("%d %d "%(t,z))
            for tok in parLine:
                ofile.write(tok+' ')
            ofile.write("\n")
            t= t+1
    ofile.close()

parFormatExportMethods= {'z_fastest':parZFastestMethod, \
                         'index_tz':parIndexTZMethod }

def exportParam(topName,parts,tmpdir):
    global parDict, parFormatDict
    verboseMessage("Exporting %s as a param file!"%topName);
    firstPart= parts[0]
    f= open(os.path.join(tmpdir,firstPart),"r")
    lines= f.readlines()
    f.close()
    ofile= file(topName,"w")
    debugMessage("Creating %s"%topName)
    parFormatDict[topName]= "unknown"
    keepTheseLines= []
    for line in lines:
        if line[0:9]=='##Format:':
            if ofile != None:
                ofile.write(line) #Transfer format info
                if string.find(line[9:],'order:')>=0:
                    pairs= string.split(line[9:])
                    for p in pairs:
                        toks= string.split(p,':')
                        if string.strip(toks[0])=='order':
                            format= toks[1].split(',')[0].strip()
                            parFormatDict[topName]= format
                            debugMessage("Param file %s has format %s"%\
                                         (topName,format))
        else:
            if line[0] != '#':
                keepTheseLines.append( string.split(line) )
    ofile.close()
    for thisPart in parts[1:]:
        f= open(os.path.join(tmpdir,thisPart),"r")
        lines= f.readlines()
        f.close()
        for line in lines:
            if line[0] != '#':
                keepTheseLines.append( string.split(line) )
    parDict[topName]= keepTheseLines

def exportDat(topName,parts,tmpdir):
    verboseMessage("Ignoring %s as a .dat file!"%topName);

def exportMri(topName,parts,tmpdir):
    global first_created_mri
    verboseMessage("Exporting %s as a .mri file!"%topName);
    if len(parts)>1:
        cmd= "mri_paste -d t -out %s "%topName
        for p in parts:
            cmd= "%s %s"%(cmd,os.path.join(tmpdir,p))
    else: 
        cmd= "mri_copy_dataset %s %s "%\
             (os.path.join(tmpdir,parts[0]),topName)
    safeRun(cmd)
    if first_created_mri==None:
        first_created_mri= topName

def exportSummaryInput(topName,parts,tmpdir):
    verboseMessage("Handling %s as a summary.input file!"%topName);
    firstPart= parts[0]
    f= open(os.path.join(tmpdir,firstPart),"r")
    lines= f.readlines()
    f.close()
    if len(lines)>0:
        ofile= file(key,"a")
        for l in lines:
            ofile.write("%s.%d\n"%(l[0:string.rfind(l,".")],os.getpid()))
        ofile.close()

def exportSummaryMissing(topName,parts,tmpdir):
    global first_created_mri, summary_missing_file, summary_missing_step
    verboseMessage("Handling %s as a summary.input file!"%topName);
    firstPart= parts[0]
    f= open(os.path.join(tmpdir,firstPart),"r")
    lines= f.readlines()
    f.close()
    if len(lines)>0:
        l= lines[0]
        summary_missing_file= topName
        summary_missing_step= "%s.%d"%(l[0:string.rfind(l,".")],os.getpid())

def exportAny(topName,parts,tmpdir):
    verboseMessage("Ignoring %s as a any file!"%topName);

exportMethodPairs= [ (paramRegex,exportParam), (mriRegex,exportMri), \
                     (datRegex,exportDat), \
                     (summaryInputRegex,exportSummaryInput), \
                     (summaryMissingRegex,exportSummaryMissing), \
                     (anyRegex,exportAny) ] # 'any' must be last!

def checkEnvInt(val, name):
    if os.environ.has_key(name):
        result= string.atoi(os.environ[name])
    else:
        result= val
    return result

def checkEnvFloat(val, name):
    if os.environ.has_key(name):
        result= string.atof(os.environ[name])
    else:
        result= val
    return result

def checkEnvString(val, name):
    if os.environ.has_key(name):
        result= os.environ[name]
    else:
        result= val
    return result

def takeCensus(path, censusDict):
    for file in os.listdir(path):
        kidPath= os.path.join(path,file)
        if os.path.isdir(kidPath):
            censusDict= takeCensus(kidPath,censusDict)
        else:
            censusDict[kidPath]= os.stat(kidPath)[8]
            debugMessage("%s has time value %d"%(kidPath,censusDict[kidPath]))
    return censusDict

##############################
#
# Main
#
##############################

mainChunk= "samples"

# Check for "-help"
if len(sys.argv)>1:
    if sys.argv[1] == "-help":
	if len(sys.argv)>2:
	    os.system( "scripthelp %s %s"%(sys.argv[0],sys.argv[2]) )
	else:
	    os.system( "scripthelp %s"%sys.argv[0] )
	sys.exit()

try:
    (opts,pargs) = getopt.getopt(sys.argv[1:],"vdc:",[])
except:
    errorMessage("%s: Invalid command line parameter" % sys.argv[0])
    describeSelf()
    sys.exit(1)

Message(idString)

# Parse flags
for a,b in opts:
    if a=="-v":
        setVerbose(1)
    if a=="-d":
	setDebug(1)
    if a=="-c":
        mainChunk= b

# Check calling syntax; parse args
if len(pargs) < 2:
    describeSelf()
    sys.exit(1)

inDSList= []
outDSList= []

cmdScript= pargs[0]
verboseMessage("Command script is <%s>"%cmdScript)

fullCmd= "%s"%cmdScript
pargs= pargs[1:]
for thisArg in pargs:
    if os.path.isabs(thisArg):
        sys.exit("This script only works with inputs in relative paths.  %s is an absolute path!"%thisArg)
    fullCmd= "%s %s"%(fullCmd,thisArg)
    if dsExists(thisArg):
        debugMessage("Treating %s as input"%thisArg)
        inDSList.append(thisArg)
    else:
        debugMessage("Treating %s as output"%thisArg)
        outDSList.append(thisArg)

#Check reasonableness of input
tdim= 0
zdim= 0
prevDS= None
if len(inDSList) == 0:
    sys.exit("None of the given datasets exist!");
for thisDS in inDSList:
    ds= MRIDataset(thisDS)
    if not ds.hasChunk("acq_blocks"):
        sys.exit("Dataset %s has no acq_blocks chunk!"%thisDS)
    if not ds.hasChunk(mainChunk):
        sys.exit("Dataset %s has no %s chunk!"%(thisDS,mainChunk))
    acq= ds.getChunk("acq_blocks")
    main= ds.getChunk(mainChunk)
    if acq.getValue("dimensions") != 't':
        sys.exit("Dataset %s acq_blocks dimstring is not t!"%thisDS)
    if tdim==0:
        tdim= acq.getDim("t")
    elif acq.getDim("t") != tdim:
        sys.exit("Dataset %s acq_blocks tdim doesn't match %s!"%\
                 (thisDS,prevDS))
    if tdim != main.getDim("t"):
        sys.exit("acq_blocks tdim does not match that of %s chunk!"%\
                 (mainChunk))
    if zdim == 0:
        zdim= main.getDim("z")
    elif main.getDim("z") != zdim:
        sys.exit("Dataset %s has a different z dimension from %s!"%\
                 (thisDS,prevDS))
    prevDS= thisDS
                 
debugMessage("tdim= %d, zdim= %d"%(tdim,zdim))

# Make up a temporary directory
tmpdir= makeTempDir("tmp_do_by_blocks")
homedir= os.getcwd()

# Figure out the block indices and boundaries
blockNumList= map(int,
                  readCmdOutputToList("mri_rpn_math -chunk acq_blocks -out %s/junk '0,$1,1,if_print_1' %s"%(tmpdir,inDSList[0])))
prevBlockNum= -1
prevBlockStart= -1
counter= 0
blockList= []
for blockNum in blockNumList:
    if blockNum != prevBlockNum:
        if prevBlockNum >= 0:
            subDir= "block_%d"%prevBlockNum
            blockList.append((prevBlockNum,prevBlockStart,counter-1,subDir))
        prevBlockNum= blockNum
        prevBlockStart= counter
    counter= counter+1
subDir= "block_%d"%prevBlockNum
blockList.append((prevBlockNum,prevBlockStart,counter-1,subDir))

if getDebug():
    debugMessage("###### Block list follows:")
    for block in blockList:
        debugMessage("block %d: %d to %d, subdir %s"%block)

os.chdir(tmpdir)

beforeInfo= []
for block in blockList:
    ( num, start, end, subDir ) = block
    verboseMessage("###### Setting up times %d through %d in %s"%\
                   (start,end,subDir))
    os.mkdir( subDir )
    os.chdir( subDir )
    for thisDS in inDSList + outDSList:
        checkPathExists(thisDS)
    for thisDS in inDSList:
        wholeDS= os.path.join(homedir,thisDS)
        partDS= os.path.join(os.curdir,thisDS)
        safeRun("mri_subset -d t -len %d -shift %d %s %s"%((end-start+1),\
                                                           start, \
                                                           wholeDS,partDS))
    beforeInfo.append( takeCensus(".", {}) )
    os.chdir( '..' )

afterInfo= []
for block in blockList:
    ( num, start, end, subDir ) = block
    verboseMessage("###### Running times %d through %d in %s"%\
                   (start,end,subDir))
    os.chdir( subDir )
    safeRun( fullCmd )
    afterInfo.append( takeCensus(".", {}) )
    os.chdir( '..' )
verboseMessage("###### Finished applying script to individual blocks")
                                   
changedDict= {}
createdDict= {}
for block in blockList:
    ( num, start, end, subDir ) = block
    debugMessage("###### Classifying block %d files."%num)
    thisBefore= beforeInfo[num]
    thisAfter= afterInfo[num]
    for key in thisAfter:
        if paramRegex.match(key):
            debugMessage("<%s> is a parameter file!"%key)
            topKey= "%s.%d"%(key[0: string.rfind(key,'.')], os.getpid())
        else:
            topKey= key
        if not thisBefore.has_key(key):
            debugMessage("<%s> is newly created!"%key)
            if num == 0:
                createdDict[topKey]= [ os.path.join(subDir,key) ]
            elif createdDict.has_key(topKey):
                createdDict[topKey].append( os.path.join(subDir,key) )
            else:
                sys.exit("Unexpectedly found new file %s late in the game!"%\
                         key)
        elif thisBefore[key]<thisAfter[key]:
            debugMessage("<%s> was modified!"%key)
            if num == 0:
                changedDict[topKey]= [ os.path.join(subDir,key) ]
            elif changedDict.has_key(topKey):
                changedDict[topKey].append( os.path.join(subDir,key) )
            else:
                sys.exit("Unexpectedly found new file %s late in the game!"%\
                         key)
        else:
            debugMessage("<%s> is unchanged."%key)
    
if getDebug():
    debugMessage("###### Created files:")
    for key in createdDict:
        debugMessage("%s -> %s"%(key,createdDict[key]))
    debugMessage("###### Changed files:")
    for key in changedDict:
        debugMessage("%s -> %s"%(key,createdDict[key]))

# Commence constructing the output files!
os.chdir(homedir)

# Build the output files
for key in createdDict:
    checkPathExists(key)
    for regex, mthd in exportMethodPairs:
        if regex.match(key):
            mthd(key,createdDict[key],tmpdir)
            break

for key in changedDict:
    for regex, mthd in exportMethodPairs:
        if regex.match(key):
            mthd(key,changedDict[key],tmpdir)
            break

# Generate entries for summary.missing
if summary_missing_file != None and summary_missing_step != None \
       and first_created_mri != None:
    verboseMessage("Generating %s info from %s"%\
                   (summary_missing_file,first_created_mri))
    ofile= file(summary_missing_file,"a")
    ofile.write("%s\n"%summary_missing_step)
    cmd= "count_missing.csh %s"%first_created_mri
    lines= readCmdOutputToList(cmd)
    for l in lines:
        ofile.write(l)
    ofile.close()

# Fill the exported parameter files
for fname in parDict:
    parList= parDict[fname]
    format= parFormatDict[fname]
    if parFormatExportMethods.has_key(format):
        parFormatExportMethods[format](fname,parList)
    else:
        Message("Unable to export %s: format %s is unknown."%(fname,format))

# Clean up
if not getDebug():
    removeTmpDir(tmpdir)


