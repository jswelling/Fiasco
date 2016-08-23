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
# This script should not be run directly; it's library routines
# for other Python scripts.
#
_idString_= "$Id: fiasco_utils.py,v 1.33 2007/04/12 23:44:35 welling Exp $"

import sys
import os
import string
import getopt
import math
import copy
import array
from quaternion import *

global _debug_
global _verbose_
_verbose_= 0
_debug_= 0

def setVerbose( val ):
    global _verbose_
    _verbose_= val

def getVerbose():
    return _verbose_

def setDebug( val ):
    global _debug_
    _debug_= val

def getDebug():
    return _debug_

def Message( str ):
    sys.stderr.write("%s\n"%str)

def verboseMessage( str ):
    if _verbose_:
        Message(str)

def debugMessage( str ):
    if _debug_:
        Message(str)

def errorMessage( str ):
    Message(str)

def makeTempDir( baseString ):
    if os.environ.has_key("F_TEMP"):
        tmpdir= "%s/%s_%d"%(os.environ["F_TEMP"],baseString,os.getpid())
    else:
        tmpdir= "./%s_%d"%(baseString,os.getpid())
    verboseMessage("temp directory is %s"%tmpdir)
    os.makedirs(tmpdir)
    return tmpdir

def removeTmpDir(path):
    debugMessage("removing temporary directory tree %s"%path)
    for file in os.listdir(path):
        kidPath= os.path.join(path,file)
        if os.path.isdir(kidPath):
            removeTmpDir(kidPath)
        else:
            os.remove(kidPath)
    os.rmdir(path)

def safeRun(cmd):
    debugMessage("running <%s>"%cmd)
    if _debug_:
        (cmdin,cmdout)= os.popen4(cmd,'t',128)
        try:
            for line in cmdout:
                debugMessage(string.strip(line))
        except StopIteration:
            pass
    else:
        cmdout= os.popen(cmd)
        cmdout.readlines() # Throw away the output
    if cmdout.close() != None :
	sys.exit("Command failed: <%s>"%cmd)

def readCmdOutputToList(cmd):
    debugMessage("running <%s>"%cmd)
    cmdout= os.popen(cmd)
    result= cmdout.readlines()
    if cmdout.close() != None :
	sys.exit("Command failed: <%s>"%cmd)
    debugMessage("read %d entries"%len(result))
    return result

def writeListToCmdInput(myList,cmd):
    debugMessage("running <%s>"%cmd)
    cmdin= os.popen( cmd, 'w')
    for x in myList:
        cmdin.write("%s "%str(x))
    if cmdin.close() != None :
        sys.exit("Command failed: <%s>"%cmd)

def writeStringToCmdInput(myString,cmd):
    debugMessage("running <%s>"%cmd)
    cmdin= os.popen( cmd, 'w')
    cmdin.write(myString)
    if cmdin.close() != None :
        sys.exit("Command failed: <%s>"%cmd)

def getField(mrifile,chunk,field):
    cmd= "mri_printfield -field %s.%s %s" % (chunk,field,mrifile)
    debugMessage("running <%s>"%cmd)
    cmdout= os.popen(cmd)
    xstr= cmdout.read()
    if cmdout.close() != None :
	sys.exit("mri_printfield failed for %s.%s on %s!"%(chunk,field,mrifile))
    return string.strip(xstr)

def getFieldNofail(mrifile,chunk,field):
    cmd= "mri_printfield -field %s.%s -nofail %s" % (chunk,field,mrifile)
    debugMessage("running <%s>"%cmd)
    cmdout= os.popen(cmd)
    xstr= cmdout.read()
    if cmdout.close() != None :
	sys.exit("mri_printfield failed for %s.%s on %s!"%(chunk,field,mrifile))
    xstr= string.strip(xstr)
    if len(xstr)>0:
        return xstr
    else:
        return None

def getDim(mrifile,chunk,index):
    return string.atoi(getField(mrifile,chunk,"extent.%s"%index))

def getDimStr(mrifile,chunk):
    return getField(mrifile,chunk,"dimensions")

def describeSelf():
    os.system("scripthelp %s usage"%sys.argv[0])

def checkPathExists( thisDS ):
    dir= os.path.dirname(thisDS)
    if dir != "":
        if not os.access(dir,os.F_OK):
            verboseMessage("Making directories to <%s>"%dir)
            os.makedirs(dir)
        if not os.access(dir,os.W_OK):
            sys.exit("Cannot write to location of %s!"%thisDS)

def dsExists( thisDS ):
    if os.path.splitext(thisDS)[1]=='.mri':
        return os.access(thisDS,os.F_OK)
    else:
        return os.access("%s.mri"%thisDS,os.F_OK)

def chunkExists( thisDS, chunk ):
    cmd= "mri_printfield -field %s -nofail %s" % (chunk,thisDS)
    debugMessage("running <%s>"%cmd)
    cmdout= os.popen(cmd)
    xstr= cmdout.read()
    if cmdout.close() != None :
	sys.exit("mri_printfield failed for %s on %s!"%(chunk,thisDS))
    xstr= string.strip(xstr)
    return (xstr=='[chunk]')

def checkExeFound( exeName, pkgName ):
    "Throw an exception if exeName (from package pkgName) is not in PATH"
    try:
        whichStr= readCmdOutputToList("which %s"%exeName)[0]
    except:
        raise Exception("%s executable %s is not in PATH"%\
                        (pkgName,exeName))
    splitWhichStr= string.split(string.strip(whichStr))
    if len(splitWhichStr) != 1:
        raise Exception("%s executable %s is not in PATH"%\
                        (pkgName,exeName))

def isNewer( path1, path2 ):
    "Return true if file path1 is newer than file path2, based on modification time"
    t1= os.stat(path1)[8]
    t2= os.stat(path2)[8]
    return ( t1 > t2 )

def parseCSV( ifile ):
    "returns a tuple containing a list of keys and a list of dicts"
    lineList= []
    possibleDelimiters= [",","\t",None] # empty string means whitespace-delimited
    lines= ifile.readlines()
    delimFound= 0
    delimForThisFile= None
    for delim in possibleDelimiters:
        wordCount= len(lines[0].split(delim))
        if wordCount<2: continue
        for line in lines[1:]:
            nwords= len(line.split(delim))
            if nwords>0 and nwords != wordCount:
                debugMessage("Delim is not <%s>\n"%delim)
                break
        else:
            delimFound= 1
            delimForThisFile= delim
            break
    if not delimFound:
        sys.exit("Cannot find the right delimiter for this CSV input!")
    debugMessage("delimForThisFile= <%s>"%delimForThisFile)
    keys= lines[0].split(delimForThisFile)
    keys= [ x.strip() for x in keys ]
    stringsAreQuoted= 1
    for key in keys:
        if len(key)>0 and (not key.startswith('"') or not key.endswith('"')):
            stringsAreQuoted= 0
    debugMessage("stringsAreQuoted= %d"%stringsAreQuoted)
    if stringsAreQuoted:
        keys = [ x[1:-1] for x in keys ]
    lines= lines[1:]
    lineNum= 1
    for line in lines:
        words= line.split(delimForThisFile)
        words= [x.strip() for x in words]
        if len(words)>0:
            dict= {}
            if len(words)!=len(keys):
                errorMessage("Line length error: %d vs %d"%\
                             (len(words),len(keys)))
                for i in xrange(len(keys)):
                    errorMessage("%d: <%s> <%s>"%(i,keys[i],words[i]))
                sys.exit("Line length error parsing CSV at line %d"%(lineNum))
            for i in xrange(len(keys)):
                if stringsAreQuoted and words[i].startswith('"') \
                       and words[i].endswith('"'):
                    dict[keys[i]]= words[i][1:-1]
                else:
                    try:
                        dict[keys[i]]= int(words[i])
                    except ValueError:
                        try:
                            dict[keys[i]]= float(words[i])
                        except ValueError:
                            dict[keys[i]]= words[i]

            lineList.append(dict)
        lineNum += 1
    return (keys, lineList)

def writeCSV( ofile, keyList, recDictList, delim=",", quoteStrings=0 ):
    if quoteStrings:
        ofile.write('"%s"'%keyList[0])
        for key in keyList[1:]:
            ofile.write('%s"%s"'%(delim,key))
    else:
        ofile.write("%s"%keyList[0])
        for key in keyList[1:]:
            ofile.write("%s%s"%(delim,key))
    ofile.write("\n")
    for dict in recDictList:
        if dict.has_key(keyList[0]):
            val= dict[keyList[0]]
        else: val= 'NA'
        if isinstance(val,int): ofile.write("%d"%val)
        elif isinstance(val,float): ofile.write("%r"%val)
        elif quoteStrings:
            ofile.write('"%s"'%val)
        else:
            ofile.write("%s"%val)
        for key in keyList[1:]:
            if dict.has_key(key):
                val= dict[key]
            else: val= 'NA'
            if isinstance(val,int): ofile.write("%s%d"%(delim,val))
            elif isinstance(val,float): ofile.write("%s%r"%(delim,val))
            elif quoteStrings:
                ofile.write('%s"%s"'%(delim,val))
            else:
                ofile.write("%s%s"%(delim,val))
        ofile.write("\n")
    debugMessage("Wrote %d recs, delim=<%s>, quoteStrings= %d"%\
                 (len(recDictList),delim,quoteStrings))

def modifyCSVTagList(oldTagList,deleteTagList,addTagList):
    "Deletes the first list, then adds the second.  Enforces policy like moving 'Comments' to the end of the list."
    if oldTagList:
        newKeys= [x for x in oldTagList]
    else: newKeys= []
    if deleteTagList:
        for x in deleteTagList:
            if x in newKeys: newKeys.remove(x)
    if addTagList:
        for x in addTagList:
            if not x in newKeys: newKeys.append(x)
    # Guarantee presence of 'FLAGS' field
    if not 'FLAGS' in newKeys:
        newKeys.append('FLAGS')
    # Move some fields to end
    for x in ['FLAGS','Comments','comments','COMMENTS']:
        if x in newKeys:
            newKeys.remove(x)
            newKeys.append(x)
    return newKeys

def checkCSVFlagTag(dict,flagStr):
    "Check for presence of flagStr (any case) in dict['FLAGS']"
    if dict.has_key('FLAGS'):
        upperFlagStr= flagStr.upper()
        words= dict['FLAGS'].split(',')
        for word in words:
            if word.upper()==upperFlagStr:
                return 1
        return 0
    else:
        return 0
    
def addCSVFlagTag(dict,flagStr):
    "Add this string to dict['FLAGS'] if it is not already present"
    if not checkCSVFlagTag(dict,flagStr):
        if dict.has_key('FLAGS'):
            words= dict['FLAGS'].split(',')
        else:
            words= []
        for junkWord in ['NONE','NA','None','NULL']:
            if junkWord in words: words.remove(junkWord)
        words.append(flagStr)
        result= None
        for word in words:
            if result: result= "%s,%s"%(result,word)
            else: result= word
        dict['FLAGS']= result

def formRoiName(roiNum, subRoiNum=0):
    "Maps number pairs like 7,0 to '7' and 7,3 to '7c' and 7,27 to '7aa'"
    nLetters= len(string.ascii_lowercase)
    if not subRoiNum or subRoiNum==0: return "%d"%roiNum
    elif subRoiNum<=nLetters:
        return "%d%s"%(roiNum,string.ascii_lowercase[subRoiNum-1])
    elif subRoiNum<nLetters*(nLetters+1):
        return "%d%s%s"%(roiNum,
                         string.ascii_lowercase[((subRoiNum-1)/nLetters)-1],
                         string.ascii_lowercase[(subRoiNum%nLetters)-1])
    else:
        sys.exit("Sub-roi number %d is too large!"%subRoiNum)
        
def splitRoiName(roiName):
    "Returns a tuple (roiNum, subRoiNum, subRoiField) from a name like '7b'"
    digits= 1
    while digits<=len(roiName):
        if roiName[0:digits].isdigit(): digits += 1
        else: break
    digits -= 1
    if not roiName[0:digits].isdigit():
        sys.exit("Can't find the leading digit substring in roi <%s>!"%roiName)
    if digits<len(roiName) and not roiName[digits:].isalpha():
        sys.exit("Can't find trailing alpha substring in roi <%s>!"%roiName)
    roiNum= int(roiName[0:digits])
    subRoiField= roiName[digits:]
    if digits==len(roiName):
        subRoiNum= 0
    elif len(subRoiField)==1:
        subRoiNum= string.ascii_lowercase.find(subRoiField.lower())
        if subRoiNum<0:
            sys.exit("Mapping failure for subroi string <%s>!"%subRoiField)
        subRoiNum += 1
    elif len(subRoiField)==2:
        letters= subRoiField.lower()
        val= string.ascii_lowercase.find(letters[0])+1
        if val<=0:
            sys.exit("Mapping failure for subroi string <%s>!"%subRoiField)
        subRoiNum = val
        val= string.ascii_lowercase.find(letters[1])+1
        if val<=0:
            sys.exit("Mapping failure for subroi string <%s>!"%subRoiField)
        subRoiNum= subRoiNum*len(string.ascii_lowercase) + val
    else:
        sys.exit("I cannot map the subroi string <%s>!"%subRoiField)
    debugMessage("<%s> -> %d <%s> -> %d %d"%\
                 (roiName,roiNum,subRoiField,roiNum,subRoiNum))
    return (roiNum,subRoiNum,subRoiField)    

class MRIChunk:
    def __init__(self,ds,name):
        self.ds= ds
        self.name= name
        self.dict= {}
    def addPair(self,key,value):
        self.dict[key]= value
    def __str__(self):
        return "chunk %s: %s"%(self.name,str(self.dict))
    def hasValue(self,key):
        return self.dict.has_key(key)
    def getValue(self,key):
        return self.dict[key]
    def getFloat(self,key):
        return float(self.dict[key])
    def getDim(self,d):
        return int(self.getValue("extent.%s"%d))
    def isT1Weighted(self):
        te= self.getFloat('te')
        debugMessage("isT1Weighted: ds %s chunk %s has TE %f msec"%\
                     (self.ds.fname,self.name,te))
        return (te<=15000)
    def getMean(self):
        totDim= self.getDim('x')*self.getDim('y')*self.getDim('z')
        safeRun("mri_copy_dataset %s tmp_mean_1"%self.ds.fname)
        safeRun("mri_remap -order q -length %d tmp_mean_1"%totDim)
        safeRun("mri_subsample -d q -length 1 -mean tmp_mean_1 tmp_mean_2")
        lines= readCmdOutputToList("mri_rpn_math '$1,1,if_print_1' tmp_mean_2")
        return float(lines[0])
    def getMaskedMean(self, maskChunk,thresh):
        safeRun("mri_rpn_math -out tmp_mask_1 '$1,%f,<' %s"%\
                (thresh,maskChunk.ds.fname))
        safeRun("mri_rpn_math -out tmp_mask_2 '0,$1,$2,if_keep' " +\
                "%s tmp_mask_1"%self.ds.fname)
        totDim= self.getDim('x')*self.getDim('y')*self.getDim('z')
        safeRun("mri_remap -order q -length %d tmp_mask_1"%totDim)
        safeRun("mri_remap -order q -length %d tmp_mask_2"%totDim)
        safeRun("mri_subsample -d q -length 1 -sum tmp_mask_1 tmp_mask_3")
        safeRun("mri_subsample -d q -length 1 -sum tmp_mask_2 tmp_mask_4")
        lines= readCmdOutputToList("mri_rpn_math '$1,$2,/,1,if_print_1' "+\
                                   "tmp_mask_4 tmp_mask_3")
        return float(lines[0])
    def makeGradMag3D(self,vox,dsFilename):
        safeRun("mri_smooth -d x -smoother_type ddx %s tmp_x"%self.ds.fname)
        safeRun("mri_smooth -d y -smoother_type ddx %s tmp_y"%self.ds.fname)
        safeRun("mri_smooth -d z -smoother_type ddx %s tmp_z"%self.ds.fname)
        rpnProg= "$1,dup,*,%f,dup,*,/,"%vox[0] \
                 + "$2,dup,*,%f,dup,*,/,+,"%vox[1] \
                 + "$3,dup,*,%f,dup,*,/,+,sqrt"%vox[2]
        safeRun("mri_rpn_math -out %s '%s' tmp_x tmp_y tmp_z"%\
                (dsFilename,rpnProg))
        return BBox(MRIDataset(dsFilename).getChunk('images'))

class MRIDataset:
    def __init__(self,fname):
        self.fname= fname
        self.chunks= {}
        self.orphans= {}
        if os.path.splitext(fname)[1]=='.mri':
            f= file(fname,"r")
        else:
            f= file("%s.mri"%fname,"r")
        lines= f.readlines()
        f.close()
        currentChunkName= ""
        prefix= currentChunkName+'.'
        for line in lines:
            if string.find(line,'\f')>=0:
                break
            (key,val)= map(string.strip,string.split(line,'=',1))
            if val == '[chunk]':
                currentChunkName= key
                prefix= currentChunkName+'.'
                self.chunks[currentChunkName]= MRIChunk(self,key)
            else:
                if prefix!="." and string.find(key,prefix)==0:
                    # still working on this chunk
                    self.chunks[currentChunkName].addPair(key[len(prefix):],val)
                else:
                    currentChunkName= ""
                    prefix= "."
                    self.orphans[key]= val

    def __str__(self):
        chStr= ""
        for ch in self.chunks:
            print str(ch)
            chStr= chStr + "<%s> "%str(ch)
        return "dataset %s: orphans %s, chunks <%s>"%\
               (self.fname,str(self.orphans),chStr)

    def hasChunk(self,chname):
        return self.chunks.has_key(chname)

    def getChunk(self,chname):
        return self.chunks[chname]

    def getOrphan(self,key):
        return self.orphans[key]

class Vec4:
    "A homogeneous coordinate 4-vector"
    def __init__(self,v1,v2,v3,v4=1):
        self.data= [v1,v2,v3,v4]
    def __str__(self): return str(self.data)
    def __repr__(self): return "Vec4(%lf,%lf,%lf,%lf)"%tuple(self.data)
    def __add__(self,other):
        return Vec4( (self.data[0]*other.data[3])+(other.data[0]*self.data[3]),
                     (self.data[1]*other.data[3])+(other.data[1]*self.data[3]),
                     (self.data[2]*other.data[3])+(other.data[2]*self.data[3]),
                     self.data[3]*other.data[3])
    def __sub__(self,other):
        return Vec4( (self.data[0]*other.data[3])-(other.data[0]*self.data[3]),
                     (self.data[1]*other.data[3])-(other.data[1]*self.data[3]),
                     (self.data[2]*other.data[3])-(other.data[2]*self.data[3]),
                     self.data[3]*other.data[3])
    def clone(self):
        return Vec4(self.data[0],self.data[1],self.data[2],self.data[3])
    def dot(self,other):
        tot= 0.0
        for i in xrange(3):
            tot= tot + self.data[i]*other.data[i]
        return tot/(self.data[3]*other.data[3])
    def mag(self): return math.sqrt(self.dot(self))
    def normalize(self):
        val= self.mag()
        for i in xrange(3):
            self.data[i]= self.data[i]/(val*self.data[3])
        self.data[3]= 1.0
    def cross(self,other):
        v1= self.data
        v2= other.data
        return Vec4(v1[1]*v2[2]-v1[2]*v2[1],
                    v1[2]*v2[0]-v1[0]*v2[2],
                    v1[0]*v2[1]-v1[1]*v2[0],
                    v1[3]*v2[3])
    def __rmul__(self,other):
        return Vec4(other*self.data[0], other*self.data[1],
                    other*self.data[2], self.data[3])
    def __getitem__(self,range): return self.data[range]
    def __setitem__(self,range,val):
        self.data[range]= val

class Transform:
    "A homogeneous 4x4 transformation matrix"
    def __init__(self,dList=[1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,\
                             0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0]):
        if len(dList) != 16:
            raise TypeError, 'bad Transform initialization list'
        self.data= copy.copy(dList)
    def __str__(self): return str(self.data)
    def __repr__(self): return "Transform(%s)"%repr(self.data)
    def __len__(self): return len(self.data)
    def __getitem__(self,range): return self.data[range]
    def __setitem__(self,range,val):
        self.data[range]= val

def getVec4( chunk, which ):
    if chunk.hasValue(which+'.0') and chunk.hasValue(which+'.1') and chunk.hasValue(which+'.2'):
        return Vec4( chunk.getFloat(which+'.0'), chunk.getFloat(which+'.1'),
                 chunk.getFloat(which+'.2') )
    else:
        return None

def getVox( chunk ):
    which= 'voxel_spacing'
    if chunk.hasValue(which+'.x') and chunk.hasValue(which+'.y') \
           and chunk.hasValue(which+'.z'):
        return [ chunk.getFloat(which+'.x'), chunk.getFloat(which+'.y'), \
                 chunk.getFloat(which+'.z') ]
    else:
        return None

class BBox:
    def __init__(self,chunk):
        self.chunk= chunk
        self.tlf= getVec4(chunk,'tlf')
        self.trf= getVec4(chunk,'trf')
        self.tlb= getVec4(chunk,'tlb')
        self.trb= getVec4(chunk,'trb')
        self.blf= getVec4(chunk,'blf')
        self.brf= getVec4(chunk,'brf')
        self.blb= getVec4(chunk,'blb')
        self.brb= getVec4(chunk,'brb')
        self.ctr= getVec4(chunk,'ctr')
        self.xdim= chunk.getDim('x')
        self.ydim= chunk.getDim('y')
        self.zdim= chunk.getDim('z')
        self.vox= getVox(chunk)
        # Try to get any missing corners
        changes= 1
        missing=1
        self.xedge= None
        self.yedge= None
        self.zedge= None
        while changes>0 and missing>0:
            missing= 0
            changes= 0
            if self.xedge==None:
                missing= missing+1
                if self.tlf != None and self.trf != None:
                    self.xedge= self.trf - self.tlf
                elif self.blf != None and self.brf != None:
                    self.xedge= self.brf - self.blf
                elif self.tlb != None and self.trb != None:
                    self.xedge= self.trb - self.tlb
                elif self.blb != None and self.brb != None:
                    self.xedge= self.brb - self.blb
                if self.xedge!=None: changes= changes+1
            if self.yedge==None:
                missing= missing+1
                if self.trb != None and self.trf != None:
                    self.yedge= self.trb - self.trf
                elif self.brb != None and self.brf != None:
                    self.yedge= self.brb - self.brf
                if self.tlb != None and self.tlf != None:
                    self.yedge= self.tlb - self.tlf
                elif self.blb != None and self.blf != None:
                    self.yedge= self.blb - self.blf
                if self.yedge!=None: changes= changes+1
            if self.zedge==None:
                missing= missing+1
                if self.tlf != None and self.blf != None:
                    self.zedge= self.tlf - self.blf
                if self.trf != None and self.brf != None:
                    self.zedge= self.trf - self.brf
                if self.tlb != None and self.blb != None:
                    self.zedge= self.tlb - self.blb
                if self.trf != None and self.brf != None:
                    self.zedge= self.trb - self.brb
                if self.zedge!=None: changes= changes+1
            if self.tlf == None:
                missing= missing+1
                if self.xedge != None and self.trf != None:
                    self.tlf= self.trf - self.xedge
                    changes= changes+1
                elif self.yedge != None and self.tlb != None:
                    self.tlf= self.tlb - self.yedge
                    changes= changes+1
                elif self.zedge != None and self.blf != None:
                    self.tlf= self.blf + self.zedge; 
                    changes= changes+1
            if self.trf == None:
                missing= missing+1
                if self.xedge != None and self.tlf != None:
                    self.trf= self.tlf + self.xedge
                    changes= changes+1
                elif self.yedge != None and self.trb != None:
                    self.trf= self.trb - self.yedge
                    changes= changes+1
                elif self.zedge != None and self.brf != None:
                    self.trf= self.brf + self.zedge; 
                    changes= changes+1
            if self.tlb == None:
                missing= missing+1
                if self.xedge != None and self.trb != None:
                    self.tlb= self.trb - self.xedge
                    changes= changes+1
                elif self.yedge != None and self.tlf != None:
                    self.tlb= self.tlf + self.yedge
                    changes= changes+1
                elif self.zedge != None and self.blb != None:
                    self.tlb= self.blb + self.zedge; 
                    changes= changes+1
            if self.trb == None:
                missing= missing+1
                if self.xedge != None and self.tlb != None:
                    self.trb= self.tlb + self.xedge
                    changes= changes+1
                elif self.yedge != None and self.trf != None:
                    self.trb= self.trf + self.yedge
                    changes= changes+1
                elif self.zedge != None and self.brb != None:
                    self.trb= self.brb + self.zedge; 
                    changes= changes+1
            if self.blf == None:
                missing= missing+1
                if self.xedge != None and self.brf != None:
                    self.blf= self.brf - self.xedge
                    changes= changes+1
                elif self.yedge != None and self.blb != None:
                    self.blf= self.blb - self.yedge
                    changes= changes+1
                elif self.zedge != None and self.tlf != None:
                    self.blf= self.tlf - self.zedge; 
                    changes= changes+1
            if self.brf == None:
                missing= missing+1
                if self.xedge != None and self.blf != None:
                    self.brf= self.blf + self.xedge
                    changes= changes+1
                elif self.yedge != None and self.brb != None:
                    self.brf= self.brb - self.yedge
                    changes= changes+1
                elif self.zedge != None and self.trf != None:
                    self.brf= self.trf - self.zedge; 
                    changes= changes+1
            if self.blb == None:
                missing= missing+1
                if self.xedge != None and self.brb != None:
                    self.blb= self.brb - self.xedge
                    changes= changes+1
                elif self.yedge != None and self.blf != None:
                    self.blb= self.blf + self.yedge
                    changes= changes+1
                elif self.zedge != None and self.tlb != None:
                    self.blb= self.tlb - self.zedge; 
                    changes= changes+1
            if self.brb == None:
                missing= missing+1
                if self.xedge != None and self.blb != None:
                    self.brb= self.blb + self.xedge
                    changes= changes+1
                elif self.yedge != None and self.brb != None:
                    self.brb= self.brf + self.yedge
                    changes= changes+1
                elif self.zedge != None and self.trf != None:
                    self.brb= self.trb - self.zedge; 
                    changes= changes+1
            missing= missing-changes
        # Here we assume that the center is the geometrical center
        # of the corners, since it would be recorded if it were
        # otherwise.
        if self.ctr == None:
            if self.tlf != None and self.brb != None:
                self.ctr= 0.5*(self.tlf+self.brb)
            elif self.trf != None and self.blb != None:
                self.ctr= 0.5*( self.trf+self.blb )
            elif self.tlb != None and self.brf != None:
                self.ctr= 0.5*( self.tlb+self.brf )
            elif self.trb != None and self.blf != None:
                self.ctr= 0.5*(self.trb+self.blf)
            elif self.tlf != None and self.blf != None \
                     and self.blb != None and self.brf != None:
                dx= self.brf-self.blf
                dy= self.blf-self.blb
                dz= self.tlf-self.blf
                self.ctr= self.blf+(0.5*(dx+dy+dz))
            else:
                sys.exit("Cannot find enough corners in dataset %s chunk %s!"%\
                         (chunk.ds.fname,chunk.name))
                    
                    

    def setCorners(self,xAxis,yAxis,zAxis):
        "Set the corner locs based on current center, size, and given unit vectorss"
        dx= (self.xdim/2)*self.vox[0] # high side
        dy= (self.ydim/2)*self.vox[1] # high side
        dz= (self.zdim/2)*self.vox[2] # high side
        self.trb= self.ctr + (dx*xAxis) + (dy*yAxis) + (dz*zAxis)
        self.tlb= self.ctr - (dx*xAxis) + (dy*yAxis) + (dz*zAxis)
        self.trf= self.ctr + (dx*xAxis) - (dy*yAxis) + (dz*zAxis)
        self.tlf= self.ctr - (dx*xAxis) - (dy*yAxis) + (dz*zAxis)
        self.brb= self.ctr + (dx*xAxis) + (dy*yAxis) - (dz*zAxis)
        self.blb= self.ctr - (dx*xAxis) + (dy*yAxis) - (dz*zAxis)
        self.brf= self.ctr + (dx*xAxis) - (dy*yAxis) - (dz*zAxis)
        self.blf= self.ctr - (dx*xAxis) - (dy*yAxis) - (dz*zAxis)

    def setCtr(self,ctr):
        self.ctr= ctr

    def setVox(self,vox):
        self.vox= copy.copy(vox)

    def printBounds(self, label):
        Message('-----------------------')
        Message(label)
        Message('  tlf: '+str(self.tlf))
        Message('  trf: '+str(self.trf))
        Message('  tlb: '+str(self.tlb))
        Message('  trb: '+str(self.trb))
        Message('  blf: '+str(self.blf))
        Message('  brf: '+str(self.brf))
        Message('  blb: '+str(self.blb))
        Message('  brb: '+str(self.brb))
        Message('  ctr: '+str(self.ctr))
        if self.chunk.hasValue('table_delta'):
            Message("   tableDelta: %f"%\
                    self.chunk.getFloat('table_delta'))
        Message('  voxel size: %s'%self.vox)
        Message('-----------------------')

    def copyBoundsFrom(self,other):
        self.trb= other.trb
        self.tlb= other.tlb
        self.trf= other.trf
        self.tlf= other.tlf
        self.brb= other.brb
        self.blb= other.blb
        self.brf= other.brf
        self.blf= other.blf
        self.ctr= other.ctr

    def exportBound(self, which, val):
        safeRun("mri_setfield -field %s.%s -all012 -value ' %f,%f,%f' %s"%\
                (self.chunk.name,which,val[0],val[1],val[2],
                 self.chunk.ds.fname))

    def exportBounds(self):
        self.exportBound("tlf",self.tlf)
        self.exportBound("trf",self.trf)
        self.exportBound("tlb",self.tlb)
        self.exportBound("trb",self.trb)
        self.exportBound("blf",self.blf)
        self.exportBound("brf",self.brf)
        self.exportBound("blb",self.blb)
        self.exportBound("brb",self.brb)
        self.exportBound("ctr",self.ctr)
        safeRun("mri_setfield -field %s.voxel_spacing -allxyz -value ' %f,%f,%f' %s"%\
                (self.chunk.name, self.vox[0], self.vox[1], self.vox[2], \
                 self.chunk.ds.fname))

def bboxFromFilename(fname):
    newDS= MRIDataset(fname)
    ch= newDS.getChunk('images')
    return BBox(ch)

class CondMatrix:
    def __init__(self,nImages,nSlices):
        self.nImages= nImages
        self.nSlices= nSlices
        self.data= array.array('i',nImages*nSlices*[-1])
    def getCond(self,t,z):
        if t>=self.nImages:
            raise Exception("Image number %d is out of range!"%t)
        if z>=self.nSlices:
            raise Exception("Slice number %d is out of range!"%t)
        return self.data[t*self.nSlices+z]
    def setCond(self,t,z,val):
        if t>=self.nImages:
            raise Exception("Image number %d is out of range!"%t)
        if z>=self.nSlices:
            raise Exception("Slice number %d is out of range!"%t)
        self.data[t*self.nSlices+z]= val
    def isUnset(self,t,z):
        if t>=self.nImages:
            raise Exception("Image number %d is out of range!"%t)
        if z>=self.nSlices:
            raise Exception("Slice number %d is out of range!"%t)
        return (self.data[t*self.nSlices+z]==-1)

def _parsePair_(word,lim,fname,lineNum):
    # debugMessage("pair: <%s>"%word)
    parts= string.split(word,'-')
    if parts[0]=='all':
        low= 0
        high= lim
    else:
        low= int(parts[0])
    if len(parts)==1:        
        high= low+1
    else:
        high= int(parts[1])+1 # add 1 because end numbers are inclusive
    if low>=high:
        Message("Disordered range pair, file %s line %d: ignored"%\
                (fname,lineNum))
        low= high= 0
    if high>lim:
        Message("Range upper bound too large, file %s line %d: set to limit"%\
                (fname,lineNum))
        high= lim
        if low>high:
            low= high
    return (low,high)

def _facListToString_( facs ):
    result= "%s"%facs[0]
    for fac in facs[1:]:
        result= "%s %s"%(result,fac)
    return result

def preparseSplitFile( inFname, fileLines, factorIndexDict, factors ):
    debugMessage("preparse of %s"%inFname)
    inFile= open(inFname,"r")
    lines= inFile.readlines()
    inFile.close()

    # Parse the set-up line of this file
    words= string.split(lines[0])
    nfac= int(words[0])
    if nfac<=0:
        sys.exit("Input file %s has invalid first line!"%inFname)
    if len(words)>1:
        bynumimagesFlag= (words[1]=="bynumimages")
    else: bynumimagesFlag= 0

    # Snag factor names for this file
    facsThisFile= string.split(lines[1])
    debugMessage("File %s factors: %s"%(inFname,facsThisFile))
    for fac in facsThisFile:
        if fac not in factors:
            factorIndexDict[fac]= len(factors)
            factors.append(fac)
    
    fileLines.append((inFname,facsThisFile,bynumimagesFlag,lines[2:]))

def parseAllSplitFiles( fileLines, factors, factorLevels,
                        factorIndexDict, factorLevelDicts,
                        condMatrix, condDict, condStringTable,
                        splitParseOpts ):

    ( spreadNAFlag )= splitParseOpts
    nSlices= condMatrix.nSlices
    nImages= condMatrix.nImages

    for fac in factors:
        factorLevelDicts[fac]= {'NA':0}
        factorLevels[fac]= []
    naMask= len(factors)*['NA']
    naString= _facListToString_(naMask)
    condDict[naString]= 0
    condStringTable.append(naString)

    # Scan for content
    baseImage= 0
    for (inFname, facsThisFile, bynumimagesFlag, lines) in fileLines:
        maxImage= 0
        currentImage= 0
        for lineNum in xrange(len(lines)):
            line= lines[lineNum]
            words= string.split(line)

            # strip comments
            for i in xrange(len(words)):
                if words[i].startswith('#'):
                    words= words[:i]
                    break

            # skip empty lines
            if len(words)==0: continue

            if len(words) < len(facsThisFile)+1:
                sys.exit("Line %d of %s is too short!"%(lineNum+2,inFname))

            # debugMessage("words for line %d: <%s>"%(lineNum,words))

            facMask= len(factors)*['NA']
            for fac in facsThisFile:
                lvl= words.pop(0)
                if not factorLevelDicts[fac].has_key(lvl):
                    debugMessage("Factor %s gets level %s at line %d"%\
                                 (fac,lvl,lineNum))
                    factorLevelDicts[fac][lvl]= len(factorLevelDicts[fac])
                    factorLevels[fac].append(lvl)
                facMask[factorIndexDict[fac]]= lvl

            if spreadNAFlag:
                # Implement the rule that one missing factor makes the
                # slice missing
                if 'NA' in facMask:
                    facMask= len(factors)*['NA']

            facString= _facListToString_(facMask)
            if not condDict.has_key(facString):
                debugMessage("Adding condition <%s> at line %d"%\
                             (facString,lineNum))
                condDict[facString]= len(condDict)
                condStringTable.append(facString)
            currentCond= condDict[facString]

            if bynumimagesFlag:
                if len(words) != 1:
                    sys.exit("Bad bynumimages record at file %s line %d"%\
                             (inFname,lineNum+2))
                startImage= currentImage
                endImage= currentImage+int(words[0])
                currentImage= endImage
                startSlice=0
                endSlice= nSlices
            else:
                (startImage,endImage)= _parsePair_(words[0],nImages,\
                                                   inFname,lineNum+2)
                if len(words)==1:
                    startSlice= 0
                    endSlice= nSlices
                elif len(words)==2:
                    (startSlice,endSlice)= _parsePair_(words[1],nSlices,\
                                                       inFname,lineNum+2)
                else:
                    sys.exit("Bad record at file %s line %d"%\
                             (inFname,lineNum+2))

            debugMessage("bounds: %d %d %d %d; baseImage %d"%\
                         (startSlice,endSlice,startImage,endImage,baseImage))

            for t in xrange(startImage,endImage):
                for z in xrange(startSlice,endSlice):
                    condMatrix.setCond(t+baseImage,z,currentCond)
            if endImage>maxImage:
                maxImage= endImage

        baseImage= maxImage


