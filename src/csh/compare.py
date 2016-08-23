#! /usr/bin/env python
#
#/************************************************************
# *                                                          *
# *  Permission is hereby granted to any individual or       *
# *  institution for use, copying, or redistribution of      *
# *  this code and associated documentation, provided        *
# *  that such code and documentation are not sold for       *
# *  profit and the following copyright notice is retained   *
# *  in the code and documentation:                          *
# *     Copyright (c) 1998,1999 Department of Statistics,    *
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
# *                                                          *
# *  Original programming by Mihir Arjunwadkar               *
# *  adapted by Jenn Bakal                                   *
# ************************************************************/
#
# This script compares two Pgh mri files or directories of Pgh mri
# files.
#

import sys
import os
import os.path
import string
import getopt
if os.environ.has_key("FIASCO"):
    sys.path.append(os.environ["FIASCO"])
from fiasco_utils import *

idString= "$Id: compare.py,v 1.2 2005/01/28 14:16:10 bakalj Exp $"

def errHandle(verbosity,status,progname,message):
    if verbosity:
        errorMessage(message % progname)
    sys.exit(status)

#def stripEndReturn(list):
#    for i in range(len(list)):
#        list[i].strip()
#        if (list[i][-1]=="\n"):
#            list[i]=list[i][:-1]


def listText(verbosity,file):
    print "made it to listText function"
    linecount=readCmdOutputToList("grep -an '^L' %s.mri | head -1 | awk -F: '{print $1-1}'" % file)
    print linecount
    if (linecount==[]):
        print "made it in if linecount"
        linecount=readCmdOutputToList("wc -l %s.mri | awk '{print $1}'" %file)
    print linecount
    text=readCmdOutputToList("head -%s %s.mri"%(string.atoi(linecount.pop()),file))
    for i in range(len(text)):
        text[i].strip()
        if (text[i][-1]=="\n"):
            text[i]=text[i][:-1]
    return text

#safeRun("head -%s %s.mri"%(string.atoi(linecount.pop()),file))

#    safeRun("head -%s %s.mri | sed -e 's/^[        ]*//' -e 's/[   ]*$//' | grep -v '^$' | awk '{for(i=1;i<NF;i++)printf "%s ", $i; printf "%s\n", $NF}'" % (string.atoi(linecount.pop()),file))
#    return linecount


def listChunks(mri_text):
    chunk_list=[]
    for line in mri_text:
#        print line
        if line.endswith(' = [chunk]'):
            chunk_list.append(line.replace(' = [chunk]',''))
#            print chunk_list
    return chunk_list            
                

def setDatafile(file,chunk):
    datpath=[]
    if dsExists(file):
        if chunkExists(file,chunk):
            dir=os.path.dirname(file)
#            print dir
            base=os.path.basename(file)
#            print base
            datpath=(dir,base)
            chunkfile=getFieldNofail(file,chunk,"file")
#            print chunkfile

            if chunkfile=='':
                physfile=file+'.mri'
            elif chunkfile[0]=='/':
                physfile=chunkfile
            elif chunkfile[0]=='.':
                if chunkfile[0:1]=='./' or chunkfile[0:2]=='../':
                    physfile=datpath+chunkfile
                else:
                    physfile=file+chunkfile
            else:
                physfile=datpath+chunkfile

            datpath = (dir,physfile)
            
#        print physfile

            
    return datpath

##############################
#
# Main
#
##############################

progname= os.path.basename(sys.argv[0]);

# Check for "-help"
if len(sys.argv)>1:
    if sys.argv[1] == "-help":
	if len(sys.argv)>2:
	    os.system( "scripthelp %s %s"%(sys.argv[0],sys.argv[2]) );
	else:
	    os.system( "scripthelp %s"%sys.argv[0] );
	sys.exit();

try:
    (opts,pargs) = getopt.getopt(sys.argv[1:],"vdt:Tqran",["verbose","debug","threshold=","Terminate","quiet","relative","absolute","notext"])
except:
    errorMessage("%s: Invalid command line parameter" % sys.argv[0])
    describeSelf();
    sys.exit()

Message(idString)

#Check calling syntax; parse args
#if len(pargs) != 2 :
#    describeSelf()
#    sys.exit(1)

comparison=""
threshold=0
terminate=0
notext=0

print opts

for a,b in opts:
    if a in ("-v","--verbose"):
        setVerbose(1)
    if a in ("-d","--debug"):
        setDebug(1)
    if a in ("-T","--terminate"):
        terminate= 1
    if a in ("-q","--quiet"):
        setVerbose(0)
    if a in ("-t","--threshold"):
        threshold = b
    if a in ("-r","--relative"):
        comparison= "relative"
    if a in ("-a","--absolute"):
        comparison= "absolute"
    if a in ("-n","--notext"):
        notext=1

files=pargs
print files

if len(files)<2:
    errHandle(getVerbose(),2,progname,"%s: Need to specify two input files!\n")

##for i in range(len(files)):
##     files[i]=files[i].replace('.mri','')
##     files[i]=files[i].replace('\n','')
##     print files


if not(dsExists(files[0])) or not(dsExists(files[1])):
    errorMessage("%s: One or more input files is nonexistent" % sys.argv[0])

ds1=MRIDataset(files[0])
ds2=MRIDataset(files[1])

print "about to compare text"

#text comparison
if not notext:
    for k in ds1.orphans.keys():
        if not ds2.orphans.has_key(k):
            errorMessage("%s: .mri files don't match at %s" %(sys.argv[0], k))
        elif not ds2.orphans[k]==ds1.orphans[k]:
            errorMessage("%s: .mri files don't match at %s" % (sys.argv[0], k))
        else:
            errorMessage("%s: .mri files match at %s" % (sys.argv[0], k))

#compare existence of chunks
for k in ds1.chunks.keys():
    if not ds2.hasChunk(k):
        errorMessage("%s: both datasets must have identical chunks" % sys.argv[0])
        sys.exit()


#chunkwise comparison

if (comparison == ""):
    comparison = "absolute"
elif (comparison == "absolute" ):
    rpnscript = "$1,$2,-,abs"
    if threshold==0:
        threshold=1
elif (comparison == "relative" ):
    rpnscript = "$1,$2,-,abs,$1abs,$2,abs,max,1,swap,dup,if_keep,/"
    if threshold==0:
        threshold=.1
else:
    errorMessage("%s: unknown comparison type" % sys.argv[0])
    sys.exit()

#Make up a temporary directory
tmpdir= makeTempDir('tmp_compare')
homedir= os.getcwd()
os.chdir(tmpdir)

for k in ds1.chunks.keys():
    diff=readCmdOutputToList("mri_rpn_math -out blah -chunk %s '0,$1,$2,-,abs,dup,%s,<,if_print_1' %s %s | wc -l"%\
                             (k,threshold,os.path.join(homedir,ds1.fname),os.path.join(homedir,ds2.fname)))
    errorMessage("number of diffs for chunk %s = %s"%(k,diff))




sys.exit()




# Clean up
os.chdir(homedir)
removeTmpDir(tmpdir)


