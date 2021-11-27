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
if "FIASCO" in os.environ:
    sys.path.append(os.environ["FIASCO"])
from fiasco_utils import *

idString= "$Id: build_model_matrix.py,v 1.10 2006/05/04 23:02:58 welling Exp $"

################################
# Notes-
################################

def slurpThisPipe(p,lines):
    while 1:
        line= p.readline()
        if len(line)==0: break
        lines.append(line)

def runChildProcess(rCommand, factors, levelTable, factorLevels, rOptions):
    (orderedFlag, orderedContrastType, unorderedContrastType)= rOptions
    stdoutLines= []
    stderrLines= []
    debugMessage("Starting <%s>"%rCommand)
    childHook= popen2.Popen3(rCommand,True)
    stdoutThread= threading.Thread(group=None,target=slurpThisPipe,\
                                   args=(childHook.fromchild,stdoutLines))
    stdoutThread.start()
    stderrThread= threading.Thread(group=None,target=slurpThisPipe,\
                                   args=(childHook.childerr,stderrLines))
    stderrThread.start()
    toChild= childHook.tochild
    debugMessage("Child process started.")

    toChild.write("options('na.action'=na.pass)\n")
    if not orderedFlag  or orderedContrastType != "none":
        toChild.write("options(contrasts=c('contr.%s','contr.%s'))\n"%\
                      (unorderedContrastType,orderedContrastType))
    for i in xrange(len(factors)):
        if orderedFlag and (contrastType == "none"):
            toChild.write("fac%d <- c(\n"%i)
            first= True
            for lvl in levelTable[i]:
                if first:
                    toChild.write("%s\n"%lvl)
                    first= False
                else:
                    toChild.write(",%s\n"%lvl)
            toChild.write(")\n")
        else:
            if orderedFlag:
                toChild.write("fac%d <- ordered(c(\n"%i)
            else:
                toChild.write("fac%d <- factor(c(\n"%i)
            first= True
            for lvl in levelTable[i]:
                if first:
                    toChild.write("'%s'\n"%lvl)
                    first= False
                else:
                    toChild.write(",'%s'\n"%lvl)
            toChild.write("),levels=c(\n")
            first= True
            for lvl in factorLevels[factors[i]]:
                if first:
                    toChild.write("'%s'\n"%lvl)
                    first= False
                else:
                    toChild.write(",'%s'\n"%lvl)
            toChild.write("))\n")

    for i in xrange(len(factors)):
        if i==0:
            terms= "'%s'=fac%d"%(factors[i],i)
            first= False
        else:
            terms= "%s, '%s'=fac%d"%(terms,factors[i],i)
    toChild.write("df <- data.frame( %s )\n"%terms)
    toChild.write("mm <- model.matrix(%s,df)\n"%formula)
    toChild.write("print('######## dims follow')\n")
    toChild.write("print(attr(mm,'dim'))\n")
    toChild.write("print('######## row labels follow')\n")
    toChild.write("attr(mm,'dimnames')[[1]]\n")
    toChild.write("print('######## column labels follow')\n")
    toChild.write("attr(mm,'dimnames')[[2]]\n")
    toChild.write("print('######## vals follow')\n")
    toChild.write("for (i in 1:length(mm)) { print(mm[[i]]) }\n")
    toChild.write("print('######## end')\n")
    toChild.write("q()\n")
    toChild.flush()
    status= childHook.wait()
    debugMessage("Child exited with status %d"%status)
    stdoutThread.join()
    stderrThread.join()
    if len(stderrLines)!=0:
        errorMessage("R wrote the following lines to stderr!")
        for line in stderrLines:
            errorMessage("<%s>"%line)
        errorMessage("These lines were written to stdout:")
        for line in stdoutLines:
            errorMessage("<%s>"%line)
        sys.exit("R shouldn't write anything to stderr.")
    if status != 0:
        sys.exit("R exited abnormally, with status %d!"%status)

    state= 0
    valList= []
    columnLabels= []
    rowLabels= []
    for line in stdoutLines:
        line= line.strip()
        debugMessage("in state %d: <%s>"%(state,string.strip(line)))
        if line[0] != "[" : continue # gets rid of echos
        if state==0:
            if line.find("dims follow")>=0: state= 1
        elif state==1:
            words= string.split(line)
            arrayDims= (int(words[1]),int(words[2]))
            debugMessage("Got array dims <%s> at slice %d"%\
                         (str(arrayDims),z))
            state=2
        elif state==2:
            if line.find("row labels follow")>=0: state= 3
        elif state==3:
            if line.find("column labels follow")>=0: state= 4
            else: rowLabels += string.split(line)[1:]
        elif state==4:
            if line.find("vals follow")>=0: state= 5
            else: columnLabels += string.split(line)[1:]
        elif state==5:
            if line.find("end")>=0:
                state= 6
                break
            else:
                word= string.split(line)[1]
                if word=='NA': valList.append('NaN')
                else: valList.append(word)
        else: break

    if state!=6:
        sys.exit("The conversation with the slave process failed!")

    return ( arrayDims, valList, rowLabels, columnLabels )

#################
#
# Main
#
#################

# Figure out what slave to use
rCommand= "R --slave --no-save --no-restore-history --no-restore-data"
if os.environ.has_key('SPLUS'):
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
    (opts,pargs) = getopt.getopt(sys.argv[1:],"vd",\
                                 ["nslices=","nimages=","out=","model=",
                                  "ordered","unordered","contrasts="])
except:
    print("%s: Invalid command line parameter" % sys.argv[0])
    describeSelf();
    sys.exit()

#Check calling syntax; parse args
if len(pargs) < 1 :
    errorMessage("Required split file name not given.")
    describeSelf()
    sys.exit(1)

nSlices= 0
nImages= 0
outFname= None
formula= None
orderedFlag= True
contrastType= None
orderedContrastType= "none"
unorderedContrastType= "treatment"
for a,b in opts:
    if a=="-v":
        setVerbose(1)
    elif a=="-d":
        setDebug(1)
    elif a=="--nslices":
        nSlices= int(b)
    elif a=="--nimages":
        nImages= int(b)
    elif a=="--out":
        outFname= os.path.abspath(b)
    elif a=="--model":
        formula= b
    elif a=="--ordered":
        orderedFlag= True
    elif a=="--unordered":
        orderedFlag= False
    elif a=="--contrasts":
        contrastType= b

verboseMessage(idString)

if formula==None:
    sys.exit("Required model formula not given")
if outFname==None:
    sys.exit("Required output file name not given")
if nSlices==0:
    sys.exit("Required number of slices not given or invalid")
if nImages==0:
    sys.exit("Required number of images not given or invalid")
if contrastType != None:
    if orderedFlag:
        if not contrastType in ["none","poly","sum"]:
            sys.exit("Unknown or unsupported ordered contrast type %s!"%
                     contrastType)
        orderedContrastType= contrastType
    else:
        if not contrastType in ["treatment","helmert"]:
            sys.exit("Unknown or unsupported unordered contrast type %s!"%
                     contrastType)
        unorderedContrastType= contrastType

fileLines= []
condDict= {}
condStringTable= []
factors= []
factorIndexDict= {}
factorLevelDicts= {}
factorLevels= {}
condMatrix= CondMatrix(nImages,nSlices)

# Slurp the files, and scan for overall set of factors
for inFname in pargs:
    preparseSplitFile( inFname, fileLines, factorIndexDict, factors )

parseAllSplitFiles( fileLines, factors, factorLevels,
                    factorIndexDict, factorLevelDicts,
                    condMatrix, condDict, condStringTable,
                    (False) )
            
verboseMessage("Factors: %s"%factors)
for fac in factors:
    verboseMessage("Levels for %s: %s"%(fac,factorLevels[fac]))
verboseMessage("condDict: %s"%condDict)

tmpdir= makeTempDir('tmp_model_matrix')
homedir= os.getcwd()
os.chdir(tmpdir)

filesToPaste= []
for z in xrange(0,nSlices):
    levelTable= []
    for fac in factors: levelTable.append([])
    for t in xrange(0,nImages):
        if condMatrix.isUnset(t,z):
            for i in xrange(len(factors)):
                levelTable[i].append('NA')
        else:
            levels= string.split(condStringTable[condMatrix.getCond(t,z)])
            for i in xrange(len(factors)):
                levelTable[i].append(levels[i])

    verboseMessage("Slice %d..."%z)
    ( arrayDims, valList, rowLabels, columnLabels )= \
      runChildProcess(rCommand, factors, levelTable, factorLevels,
                      (orderedFlag,orderedContrastType,unorderedContrastType)) 
    
    writeListToCmdInput(valList,
                        "mri_from_ascii -order tv -length %d:%d tmp"%\
                        arrayDims)
    safeRun("mri_permute -order vt tmp slice_%d"%z)
    safeRun("mri_remap -order vtz slice_%d"%z)
    filesToPaste.append("slice_%d"%z)

verboseMessage("Doing final assembly.")
if len(filesToPaste)>1:
    loopCnt= 0
    while len(filesToPaste)>0:
        if loopCnt==0:
            args= ""
            for fname in filesToPaste[0:16]:
                args= "%s %s"%(args,fname)
            safeRun("mri_paste -d z -out tmp_pasted_%d %s"%(loopCnt,args))
            filesToPaste= filesToPaste[16:]
        else:
            args= ""
            for fname in filesToPaste[0:15]:
                args= "%s %s"%(args,fname)
            safeRun("mri_paste -d z -out tmp_pasted_%d tmp_pasted_%d %s"%\
                    (loopCnt,loopCnt-1,args))
            safeRun("mri_destroy_dataset tmp_pasted_%d"%(loopCnt-1))
            filesToPaste= filesToPaste[15:]
        loopCnt= loopCnt+1
    safeRun("mri_copy_dataset tmp_pasted_%d pasted"%(loopCnt-1))
    safeRun("mri_destroy_dataset tmp_pasted_%d"%(loopCnt-1))
else:
    safeRun("mri_copy_dataset %s pasted"%filesToPaste[0])

# Strip off the intercept column if present, for compatibility with
# mri_glm (which adds an intercept column)
if columnLabels[0]== '"(Intercept)"':
    columnLabels= columnLabels[1:]
    safeRun("mri_subset -d v -len %d -shift 1 pasted %s"%\
            (len(columnLabels),outFname))
else:
    safeRun("mri_copy_dataset pasted %s"%outFname)

# Add label strings to output, so user can tell which term is which.
# Factor level arrays are the same from slice to slice, so labels
# should be also.
verboseMessage("Column labels: %s"%columnLabels)
for i in xrange(len(columnLabels)):
    label= columnLabels[i]
    debugMessage( "label %d: <%s>"%(i,label) )
    safeRun("mri_setfield -fld images.label.v.%d -val %s %s"%\
            (i,label,outFname))

# Clean up
os.chdir(homedir)
if not getDebug():
    removeTmpDir(tmpdir)


