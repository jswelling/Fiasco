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
import stat
import string
import getopt
import urllib
import re
import grp
if "FIASCO" in os.environ:
    sys.path.append(os.environ["FIASCO"])
from fiasco_utils import *

idString= "$Id: multi_runner.py,v 1.6 2006/05/04 23:02:58 welling Exp $"

db_uname= None
db_passwd= None
tagString= None
infileName= None
srcDir= None
destSubdir= None

pfileRegex= re.compile('P(\d+)\.(\d)')

runList= []

def setPathProtection(path):
    grpInfo= grp.getgrnam("concussion")
    os.chown(path,os.geteuid(),grpInfo[2])
    currentMode= os.stat(path).st_mode
    os.chmod(path,currentMode | stat.S_IRWXG)

def getpass(prompt = "Password: "):
    import termios, sys
    fd = sys.stdin.fileno()
    old = termios.tcgetattr(fd)
    new = termios.tcgetattr(fd)
    new[3] = new[3] & ~termios.ECHO          # lflags
    try:
        termios.tcsetattr(fd, termios.TCSADRAIN, new)
        passwd = raw_input(prompt)
        termios.tcsetattr(fd, termios.TCSADRAIN, old)
        print("")
    finally:
        termios.tcsetattr(fd, termios.TCSADRAIN, old)
    return passwd

def makeAnalysisPath(spec):
    global destSubdir
    id, task, runnum= spec
    return "/home/YODA/roushre/ANALYSIS/%03d/%s_%d/%s"%\
           (id,task,runnum,destSubdir)

def makeDataPath(spec):
    id, task, runnum= spec
    return "/home/YODA/roushre/DATA/Pfiles/%03d/%s_%d"%(id,task,runnum)

def buildAndCheckAnalysis( path, spec ):
    global db_passwd, srcDir, tagString
    id, task, runnum= spec
    statusFlag= 1
    Message("Checking %d %s %d"%spec)
    dataDir= makeDataPath(spec)
    if not os.access( dataDir, os.R_OK ):
        Message("Input data directory %s does not exist!"%dataDir)
        return 0
    pfilesFound= 0
    subdirsFound= 0
    for f in os.listdir(dataDir):
        if pfileRegex.match(f):
            pfilesFound= 1
            break
        if os.path.isdir(os.path.join(dataDir,f)):
            subdirsFound= 1
            break
    if not pfilesFound and not subdirsFound:
        Message("No Pfiles or subdirectories found in %s!"%dataDir)
        return 0
    if os.access( path, os.F_OK ):
        Message("%s already exists!"%path)
        return 0
    try:
        os.makedirs(path)
    except OSError as s:
        Message("Error creating %s: %s!"%(path,s))
    setPathProtection(path)
    homeDir= os.getcwd()
    try:
        os.chdir(path)
        os.environ['F_DB_PASSWD']= db_passwd
        argStr= "--src=%s db_uname=%s subj=%d runnum=%d task=%s tag=%s"%\
                (srcDir,db_uname,id,runnum,task,tagString)
        for fname in [ "spiral.local.csh", "spiral.steps.csh", \
                       "fiasco.local.csh", "epi.local.csh", \
                       "epi.steps.csh", "run_me.csh" ]:
            if os.access(os.path.join(srcDir,fname),os.X_OK):
                safeRun("customize_local.py %s %s"%(argStr,fname))
    except:
        Message("An error occurred while customizing files in %s!"%path)
        statusFlag= 0
    os.chdir(homeDir)
    return statusFlag

def doFiascoRun( spec ):
    id, task, runnum= spec
    path= makeAnalysisPath(spec)
    Message("Executing %d %s %d in %s"%(id,task,runnum,path))
    homedir= os.getcwd()
    try:
        os.chdir(path)
        if os.access( "./run_me.csh", os.X_OK ):
            safeRun("./run_me.csh")
        else:
            safeRun("FIASCO spiral")
        Message("This run succeeded.")
    except:
        Message("This run failed!")
    os.chdir(homedir)

##############################
#
# Main
#
##############################

# Check for "-help"
if len(sys.argv)>1:
    if sys.argv[1] == "-help":
        if len(sys.argv)>2:
            os.system( "scripthelp %s %s"%(sys.argv[0],sys.argv[2]) )
        else:
            os.system( "scripthelp %s"%sys.argv[0] )
        sys.exit()

try:
    (opts,pargs) = getopt.getopt(sys.argv[1:],"vd",\
                                 ["src=","db_uname=","db_passwd=","tag=",\
                                  "infile=","dest="])
                                                    
except:
    errorMessage("%s: Invalid command line parameter" % sys.argv[0])
    describeSelf()
    sys.exit(1)

Message(idString)

# Parse args
for a,b in opts:
    if a=="-v":
        setVerbose(1)
    if a=="-d":
        setDebug(1)
    if a=="--src":
        srcDir= b
    if a=="--db_uname":
        db_uname= b
    if a=="--db_passwd":
        db_passwd= b
    if a=="--tag":
        tagString= b
    if a=="--infile":
        infileName= b
    if a=="--dest":
        destSubdir= b

while db_uname==None or len(db_uname)==0:
    db_uname= raw_input("Database Username: ")
while db_passwd==None or len(db_passwd)==0:
    db_passwd= getpass("Database Password: ")
while tagString==None or len(tagString)==0:
    tagString= raw_input("Tag for this run (e.g. Analysis#12): ")
if srcDir==None:
    srcDir= raw_input("Source directory (defaults to FIASCO): ")
if srcDir[0] != '/':
    srcDir= os.path.join(os.getcwd(),srcDir)
while destSubdir==None or len(destSubdir)==0:
    destSubdir= raw_input("Destination subdirectory (e.g. auto_7): ")
if destSubdir[0]=='/':
    sys.exit("Destination subdirectory must be a 1-word relative path!")
while infileName==None or len(infileName)==0 \
          or not os.access(infileName,os.R_OK):
    infileName= raw_input("Table file name: ")

infile= file(infileName,"r")
lines= infile.readlines()
infile.close()

for line in lines:
    line= string.strip(line)
    if len(line) == 0:
        continue
    if line[0]=='#':
        continue  # this is a comment
    idStr,task,runnumStr = string.split(line)
    id= int(idStr)
    runnum= int(runnumStr)
    spec= (id, task, runnum)
    workPath= makeAnalysisPath(spec)
    if buildAndCheckAnalysis( workPath, spec ):
        runList.append( spec )

for spec in runList:
    doFiascoRun(spec)

print('######## All Finished! ########')
