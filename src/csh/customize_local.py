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

#
# MySQL client library is in:
# /afs/andrew.cmu.edu/system/gamma/i386_rh80/local/depot/mysql.1056149249/lib/mysql/
# This directory wants to be in LD_LIBRARY_PATH to access DB by that
# mechanism.
#

import sys
import os
import os.path
import string
import getopt
import urllib
import re
import math
import cPickle
if os.environ.has_key("FIASCO"):
    sys.path.append(os.environ["FIASCO"])
from fiasco_utils import *

idString= "$Id: customize_local.py,v 1.18 2007/07/26 19:33:11 welling Exp $"

urlKeyTransDict= { "ID":"subj_id", "RUN":"run" }
transDict= None # initialized to a dictionary in populateMethodTables
dbDict= {}
pfileRegex= re.compile('P(\d+)\.(\d)')
geLxImageRegex= re.compile('I\.\d\d\d')
siemensDicomImageRegex= re.compile('\d+-\d+-\d*\.dcm')
siemensDicomDirRegex= re.compile('D\d\d')
srcDir= None

dbCursor= None

##############################
#
# Database Access
#
##############################

def termsFromUrl( dict, name, url ):
    global transDict
    debugMessage("Reading URL <%s>"%url)
    f= urllib.urlopen(url);
    urlLines= f.readlines();
    f.close()
    in_row_flag= 0
    recList= []
    thisRec= []
    labelList= []
    for line in urlLines:
        if string.find(line,'ERROR:') >= 0:
            Message("SQL database query failed: %s"%\
                    line[string.find(line,':')+1:])
            raise Exception("Query via URL failed: <%s>"%\
                            line[string.find(line,':')+1:])
            break
        if in_row_flag:
            if string.find(line,'%%%END') >= 0:
                in_row_flag= 0
                recList.append(thisRec)
            else:
                (key,val)= string.split(string.strip(line),':',1)
#                debugMessage("Adding <%s>=><%s>"%(key,val))
                if len(recList)==0:
                    labelList.append(key)
                thisRec.append(val)
        else:
            if string.find(line,'%%%START') >= 0:
                in_row_flag= 1
                thisRec= []
    if len(recList)==1:
        # Make direct dictionary entries for the record elements
        for i in xrange(len(labelList)):
            dict[labelList[i]]= recList[0][i]
    else:
        # Make records and labels entries.
        dict[name]= recList
        dict["%s_columns"%name]= labelList

def dbAccessPrep():
    global dbCursor
    if dbCursor != None:
        dbCursor.close()
    dbCursor= None
    try:
        sys.path.append(
            '/home/IMALIK/sai/MySQL-python-0.9.2/build/lib.linux-i686-2.2')
        import MySQLdb
    except ImportError, e:
        debugMessage("Failed to set up MySQLdb: %s"%e)
    try:
        dbCursor = MySQLdb.connect(host='temper',
                                   user=transDict['db_uname'],
                                   passwd=transDict['db_passwd'],
                                   db='Concussion_Study').cursor()
    except Exception, e:
        debugMessage("Failed to connect to MySQL database: %s"%e)

def termsFromQuery( dict, dbTable, qualifierList, orderList=[] ):
    global dbCursor
    if dbCursor != None:
        dbCursor.execute("describe %s;"%dbTable)
        columns = dbCursor.fetchall()
        query= "select * from %s"%dbTable
        if len(qualifierList)>0:
            first= 1
            query= query+" where "
            for qual in qualifierList:
                key, val= qual
                if first:
                    query= query+"%s = %s"%(key,val)
                    first= 0
                else:
                    query= query+" and %s = %s"%(key,val)
        if len(orderList)>0:
            first= 1
            query= query+" order by "
            for col in orderList:
                if first:
                    query= query+"%s"%col
                    first= 0
                else:
                    query= query+", %s"%col
        query= query+";"
        debugMessage("Executing query <%s>"%query)
        dbCursor.execute(query)
        list = dbCursor.fetchall()
        debugMessage("Resulting list has %d entries."%len(list))
        if len(list)==1:
            for i in xrange(len(columns)):
                dict[ columns[i][0] ] = str(list[0][i])
        else:
            dict[dbTable]= list
            dict["%s_columns"%dbTable]= [x[0] for x in columns]
    else:
        # Fall back to the URL-based method
#        url= "http://www.stat.cmu.edu/~welling/do_query_textonly.php?"
        url= "http://lib.stat.cmu.edu/concussion/?"
        url= url+"db_uname=%s&db_passwd=%s&"%\
             (transDict["db_uname"],transDict["db_passwd"])
        url= url+"db_table=%s"%dbTable
        for qual in qualifierList:
            key, val= qual
            if urlKeyTransDict.has_key(key):
                key= urlKeyTransDict[key]
            url= url+"&%s=%s"%(key,val)
        if len(orderList)>0:
            for i in xrange(len(orderList)):
                url= url+"&orderBy%d=%s"%(i,orderList[i])
        termsFromUrl(dict, dbTable, url)

def dbAccessFinalize():
    global dbCursor
    if dbCursor != None:
        dbCursor.close()
    dbCursor= None

def fillDbDict( dict ):
    global transDict
    termsFromQuery( dict, 'CS_ID', [("ID",transDict["subj"])] )
    termsFromQuery( dict, 'CS_SCANDATA', [("ID",transDict["subj"]),
                                          ("RUN",transDict["runnum"])] )
    if string.upper(transDict['task']) == "STROOP":
        termsFromQuery( dict, 'CS_STROOP', [("ID",transDict["subj"]),
                                            ("RUN",transDict["runnum"])],
                        ['BLOCK','TRIAL'])
    if string.upper(transDict['task']) == 'NBACK' and dbDict.has_key('NBACK'):
        pattern= dict['NBACK']
        if len(pattern)==1:
            tmpDict= {}
            termsFromQuery( tmpDict, 'LU_NBKORDER', [("letter","%s"%pattern)] )
            pattern= tmpDict['ORDER']
            dbDict['NBACK']= string.replace(pattern, ' ', '')
        termsFromQuery( dict, 'CS_NBACK', [("ID",transDict["subj"]),
                                          ("RUN",transDict["runnum"])],
                        ['BLOCK','TRIAL'])
    if string.upper(transDict['task']) == 'ARROW':
        tmpDict= {}
        termsFromQuery( tmpDict, 'CS_VOLUMES', [("ID",transDict["subj"]),
                                          ("RUN",transDict["runnum"])] )
        if tmpDict.has_key('CS_VOLUMES_columns'):
            rows= tmpDict['CS_VOLUMES']
            cols= tmpDict['CS_VOLUMES_columns']
            fnameOffset= None
            for fnameOffset in xrange(len(cols)):
                if cols[fnameOffset]=='PFILE': break
            nImagesOffset= None
            for nImagesOffset in xrange(len(cols)):
                if cols[nImagesOffset]=='NUM_VOL': break
            if fnameOffset is None or nImagesOffset is None or \
                   cols[fnameOffset]!='PFILE' \
                   or cols[nImagesOffset]!='NUM_VOL':
                raise Exception("CS_VOLUMES table lacks PFILE or NUM_VOL column!")
            pfileDict= {}
            for row in rows:
                pfileDict[row[fnameOffset]]= (row[fnameOffset],\
                                              row[nImagesOffset])
        else:
            # Only one pfile is listed in the db, so result format is different
            pfileDict= {tmpDict['PFILE']:(tmpDict['PFILE'],tmpDict['NUM_VOL'])}
        dict['PFILES']= pfileDict
        termsFromQuery( dict, 'CS_ARROW', [("ID",transDict["subj"]),
                                          ("RUN",transDict["runnum"])],
                        ['PFILE','BLOCK','TRIAL'])
        # Convert the lookup table for conditions into more useful form
        tmpDict={}
        termsFromQuery(tmpDict, 'LU_COND_ARW', [])
        rows= tmpDict['LU_COND_ARW']
        cols= tmpDict['LU_COND_ARW_columns']
        newDict= {}
        for row in rows:
            newDict[row[0]]= row[1:]
        dict['LU_COND_ARW']= newDict
        dict['LU_COND_ARW_columns']= cols

##############################
#
# Convenience Functions
#
##############################

def getTR():
    global envDict
    return float(envDict['F_IAI'](""))

def getNSlices():
    global envDict
    return int(envDict['F_NSLICE'](""))

def getScannerType():
    global dbDict
    return dbDict['SCANNER_TYPE']

def getBehaviorDir():
    global envDict, catDict, transDict, dbDict
    id= int(transDict['subj']);
    task= transDict['task']
    run= transDict['runnum']
    return "/home/YODA/roushre/DATA/behavioral_data/%03d/"%\
           (int(transDict['subj']))
    
def countImageFiles( dirName ):
    count= 0
    for fname in os.listdir(dirName):
        if geLxImageRegex.match(fname):
            count= count+1
    debugMessage("count is %d for <%s>"%(count,dirName))
    return count

def countDicomFiles( dirName ):
    count= 0
    for fname in os.listdir(dirName):
        if siemensDicomImageRegex.match(fname):
            count= count+1
    debugMessage("count is %d for <%s>"%(count,dirName))
    return count

def pulseSeqHasReferenceImage():
    global envDict
    if envDict['F_RUN_TYPE']("")=="spiral":
        return 1
    else:
        return 0

def sortByPfileOrder( fnameList ):
    vList= []
    for f in fnameList:
        m= pfileRegex.match(f)
        if m:
            suffix= m.group(2) # they had better all be the same!
            vList.append(int(m.group(1)))
    if len(vList)>0: # They are indeed proper Pfiles
        vList.sort()
        oldV= vList[0]
        firstV= vList[0]
        deltaList= [ firstV ]
        for v in vList[1:]:
            deltaList.append(v-oldV)
            oldV= v
        iBiggest= 0
        biggestDelta= deltaList[0]
        for i in xrange(1,len(deltaList)):
            if deltaList[i]>biggestDelta:
                iBiggest= i
                biggestDelta= deltaList[i]
        orderedFnameList= []
        for i in xrange(len(vList)):
            index= (i+iBiggest)%len(vList)
            orderedFnameList.append("P%05d.%s"%(vList[index],suffix))
    else: # Something else, like generic directory names
        orderedFnameList= []
        for name in fnameList:
            orderedFnameList.append(name)
        orderedFnameList.sort()
    return orderedFnameList

def appendStateTuple( stateList, name, duration, start, comment=None ):
    stateList.append( ( name, duration, start, comment ) )
    return start+duration

def walkStateList( ofile, states, slicePattern, sliceTime, baseImage, \
                   lastImage, fillState, markFirstImageNA=0 ):
    debugMessage("walkStateList %d -> %d"%(baseImage,lastImage))
    now= 0.0
    stateIndex= 0
    iState= 0
    iSlice= 0
    image= baseImage
    nSlices= len(slicePattern)
    if len(states[0])==3:
        (sampleStateName, junkDuration, junkStart)= states[0]
    else:
        (sampleStateName, junkDuration, junkStart, junkComment)= states[0]
    nConditions= len(string.split(sampleStateName))
    naState= nConditions*'NA '
    while image<lastImage:
        if len(states[iState])==3:
            (stateName, stateDuration, stateStart)= states[iState]
            comment= None
        elif len(states[iState])==4:
            (stateName, stateDuration, stateStart, comment)= states[iState]
        stateEnd= stateStart+stateDuration
        if  now > stateEnd :
            iState= iState+1
            if iState == len(states):
                break
        else:
            if markFirstImageNA and image-baseImage==0:
                state= naState
            else:
                state= stateName
##            debugMessage("State is %s at time %f, slice %d image %d"%\
##                         (state, now, slicePattern[iSlice], \
##                          image))
            if comment is None:
                ofile.write("%s %d %d\n"%(state,image,slicePattern[iSlice]))
            else:
                ofile.write("%s %d %d # %s\n"%(state,image,
                                               slicePattern[iSlice],comment))
            now= now+sliceTime
            iSlice += 1
            if iSlice >= nSlices:
                iSlice= 0
            if iSlice==0:
                image= image+1
    while image<lastImage:
        if markFirstImageNA and image-baseImage==0:
            state= naState
            comment= "fill first image with NA"
        else:
            state= fillState
            comment= "fill end of image range"
##        debugMessage("Filling with state %s at time %f, slice %d image %d"%\
##                     (state, now, slicePattern[iSlice], image))
        if comment!=None:
            ofile.write("%s %d %d # %s\n"%(state,image,\
                                          slicePattern[iSlice], comment))
        else:
            ofile.write("%s %d %d\n"%(state,image,\
                                      slicePattern[iSlice]))
        now= now+sliceTime
        iSlice += 1
        if iSlice >= nSlices:
            iSlice= 0
        if iSlice==0:
            image= image+1
    
def funcSliceOrder( nSlices, scannerType ):
    pattern= []
    if scannerType=='GE':
        # reorder pattern is even/odd
        i= 0
        while i<=nSlices-1:
            pattern.append(i)
            i += 2
        i= 1
        while i<=nSlices-1:
            pattern.append(i)
            i += 2
    elif scannerType=='SIEMENS':
        # reorder pattern is reversed_odd/even
        i= 1
        while i<=nSlices-1:
            pattern.append((nSlices-1)-i)
            i += 2
        i= 0
        while i<=nSlices-1:
            pattern.append((nSlices-1)-i)
            i += 2
    else:
        raise Exception("funcSliceOrder: unknown scanner type %s!"%scannerType)
    debugMessage("Generated slice pattern: <%s>"%pattern)
    return pattern
    
##############################
#
# Base Translation Methods
#
##############################

def dirMethod( line ):
    global envDict, catDict, transDict, dbDict
    id= int(transDict['subj']);
    task= transDict['task']
    run= transDict['runnum']
    return "/home/YODA/roushre/DATA/Pfiles/%03d/%s_%s"%\
           (int(transDict['subj']),task,run)
#    return readCmdOutputToList("map_name.py -d file=pfiles -d task=%s -d runnum=%s %s"%\
#                               (transDict["task"],transDict["runnum"],\
#                                transDict["subj"]))[0]

def inplOptsMethod( line ):
    global envDict, catDict, transDict, dbDict
    return "'-multi'"

def anatOptsMethod( line ):
    global envDict, catDict, transDict, dbDict
    return "'-multi'"

def refpMethod( line ):
    global envDict, catDict, transDict, dbDict
    pfileDir= dirMethod(line)
    if not os.access(pfileDir,os.R_OK):
        raise Exception("Pfile directory %s not found or not readable!"\
                        %pfileDir)
    fileList= []
    for f in os.listdir(pfileDir):
        if pfileRegex.match(f):
            fileList.append(f)
    if len(fileList)==0:
        raise Exception("Pfile directory %s does not contain any Pfiles!"\
                        %pfileDir)
    fileList= sortByPfileOrder(fileList)
    # Good opportunity to make sure the directory and the DB match
    if dbDict.has_key('PFILES'): # some don't have this entry
        dbVersionOfFileList= sortByPfileOrder(dbDict['PFILES'].keys())
        if fileList != dbVersionOfFileList:
            raise Exception("Pfile directory %s is not consistent with the database!"\
                            %pfileDir)        
    return fileList[0]

def headerMethod( line ):
    global envDict, catDict, transDict, dbDict
    result= "'subj %s %s %s"%(transDict["subj"],transDict["task"],\
                              transDict["runnum"])
    if transDict.has_key('tag'):
        result= "%s %s'"%(result,transDict['tag'])
    else:
        result= "%s'"%(result)
    return result

def genderMethod( line ):
    global envDict, catDict, transDict, dbDict
    if int(dbDict['GENDER'])==1:
        return "F"
    else:
        return "M"

def diag1Method( line ):
    global envDict, catDict, transDict, dbDict
    if int(dbDict['PT'])==2:
        return "none"
    else:
        return "concussion"

def ageMethod( line ):
    global envDict, catDict, transDict, dbDict
    if dbDict.has_key('SDATE') and dbDict.has_key('DOB'):
        birthDate= string.strip(dbDict['DOB'])
        if len(birthDate)==0:
            return ''
        scanDate= string.strip(dbDict['SDATE'])
        if len(scanDate)==0:
            return ''
        yBirth,mBirth,dBirth= map(int,string.split(dbDict['DOB'],'-'))
        words= string.split(dbDict['SDATE'],'-')
        cleanWords= []
        for word in words:
            cleanWords.append(string.split(word)[0])
        yScan= int(cleanWords[0])
        mScan= int(cleanWords[1])
        dScan= int(cleanWords[2])
        age= yScan-yBirth
        if mBirth < mScan or ( mBirth == mScan and dBirth <= dScan ):
            age= age+1
        return str(age)
    else:
        return ''

def descriptionMethod( line ):
    global envDict, catDict, transDict, dbDict
    tmpAge= ageMethod( line )
    if len(tmpAge)<1:
        tmpAge= '??'
    return "'%s age %s, %s"%(genderMethod(line),tmpAge,headerMethod(line)[1:])

def nsliceMethod( line ):
    global envDict, catDict, transDict, dbDict
    try:
        return int( dbDict['SLICE'])
    except:
        raise Exception("Could not parse SLICE database entry!")

def xvoxelMethod( line ):
    global envDict, catDict, transDict, dbDict
    return '3.125'

def yvoxelMethod( line ):
    global envDict, catDict, transDict, dbDict
    return '3.125'

def zvoxelMethod( line ):
    global envDict, catDict, transDict, dbDict
    return '3.2'

def estiregfixedMethod( line ):
    global envDict, catDict, transDict, dbDict
    return 'median'

def est3dalignMethod( line ):
    global envDict, catDict, transDict, dbDict
    return 'median'

def lowertMethod( line ):
    global envDict, catDict, transDict, dbDict
#    return '-6'
    return '-4.5'

def uppertMethod( line ):
    global envDict, catDict, transDict, dbDict
#    return '6'
    return '4.5'

def upperfMethod( line ):
    global envDict, catDict, transDict, dbDict
    return '36'

def lowerqMethod( line ):
    global envDict, catDict, transDict, dbDict
#    return '0.00625'
    return '0.05'

def tmpdirMethod( line ):
    global envDict, catDict, transDict, dbDict
    return "/home/YODA/%s/tmp"%os.environ['LOGNAME']

def readerDimsMethod( line ):
    global envDict, catDict, transDict, dbDict
    return "''"

def baseReaderOptsMethod( line ): return "''"
        
##############################
#
# NBACK_V1 Task Translation Methods
#
##############################

def nbackV1IaiMethod( line ): return '1.5'

def nbackV1NimageMethod( line ):
    global envDict, catDict, transDict, dbDict
    try:
        return int(dbDict['NIMAGES'])
    except:
        pass # NIMAGES field in DB is not valid
    count= 0
    for c in dbDict['NBACK']:
        if string.find(string.whitespace,c)<0:
            count = count+1
    nImages= 35*count;
    if getScannerType()=='SIEMENS': nImages += count*7 # disdaqs
    return nImages

def nbackV1SplitMethod( ofile, line ):
    global envDict, catDict, transDict, dbDict
    try:
        # BIRC files start with 2 TRs of disdaq images; those images
        # are not retained in datasets from the GE scanner at the MR Ctr
        if getScannerType()=='GE': # GE
            nDisdaqs= 0
        elif getScannerType()=='SIEMENS': # Siemens
            nDisdaqs= 7
        else:
            raise Exception("nbackV1SplitMethod: unknown scanner type <%s>"%
                            getScannerType())
        if pulseSeqHasReferenceImage():
            nRefs= 1
        else:
            nRefs= 0
        debugMessage("NBACK entry is <%s>; nDisdaqs %d, nRefs %d"%\
                     (dbDict['NBACK'],nDisdaqs,nRefs))
        ofile.write("1 bynumimages\n")
        ofile.write("Condition\n")
        ofile.write("0Back 0\n")
        ofile.write("1Back 0\n")
        ofile.write("2Back 0\n")
        firstAcq= 1 # for this task, one block per acq
        for c in dbDict['NBACK']:
            if string.find(string.whitespace,c) < 0:
                if firstAcq:
                    firstAcq= 0
                    if nDisdaqs!=0:
                        ofile.write("NA %d\n"%nDisdaqs)
                    ofile.write("%dBack %d\n"%(int(c),35))
                else:
                    if nDisdaqs != 0 or nRefs != 0:
                        ofile.write("NA %d\n"%(nRefs+nDisdaqs))
                    ofile.write("%dBack %d\n"%(int(c),(35-nRefs)))
    except Exception,e:
        print "Caught exception: %s"%e
        raise Exception("Generation of NBACK V1 split failed!")

##############################
#
# NBACK_V2 Task Translation Methods
#
##############################

def nbackV2IaiMethod( line ): return '2.0'

def nbackV2NimageMethod( line ):
    global envDict, catDict, transDict, dbDict
    try:
        return int(dbDict['NIMAGES'])
    except:
        pass # NIMAGES field in DB is not valid
    # 3 acquisitions of 141, plus 2 disdaqs if BIRC
    nImages= 3*141
    if getScannerType()=='SIEMENS': nImages += 3*2 # disdaqs
    return nImages

def nbackV2BuildStateList(trialOrder,nDisdaqs):
    states= []
    start= 0.0
    # Cover disdaqs, if any
    if nDisdaqs>0:
        start= appendStateTuple( states, 'NA', nDisdaqs*getTR(), start )
    # Lag to compensate for BOLD
    start= appendStateTuple( states, 'NA', 2.0, start )    
    # pre-acq fix, plus disdaqs
    start= appendStateTuple( states, 'fixation', 2.0, start )
    for blockTrials in trialOrder:
        debugMessage("nbackV2BuildStateList %s: building block:<%s>"%\
                     ("NBACK_V2",blockTrials))
        start= appendStateTuple( states, 'fixation', 7.0, start )
        start= appendStateTuple( states, 'directions', 4.0, start )
        for (cond,reps) in blockTrials:
            for trial in xrange(reps):
                start= appendStateTuple( states, '%s'%cond, \
                                         0.5, start ) #stim
                start= appendStateTuple( states, '%s'%cond, \
                                         1.7, start ) #fix
        start= appendStateTuple( states, "fixation", 2.0, start )
    start= appendStateTuple( states, "fixation", 4.0, start ) # post-acq fix
    return states

def nbackV2SplitMethod( ofile, line ):
    global envDict, catDict, transDict, dbDict
    if dbDict['FULL_TASK_NAME'] != 'NBACK_V2':
        raise Exception("Internal error: %s != NBACK_V2"%\
                        dbDict['FULL_TASK_NAME'])
    try:
        ofile.write("1\n")
        ofile.write("Condition\n")
        ofile.write("0back 0 0\n")
        ofile.write("1back 0 0\n")
        ofile.write("2back 0 0\n")
        ofile.write("fixation 0 0\n")
        ofile.write("directions 0 0\n")
        #ofile.write("startBlock 0 0\n")
        #ofile.write("error 0 0\n")
        debugMessage("NBACK entry is <%s>"%dbDict['NBACK'])
        
        dbTable= dbDict['CS_NBACK']
        if dbTable==None:
            raise Exception("No CS_NBACK entry in database!")
        dbCols= dbDict['CS_NBACK_columns']
        if dbCols==None:
            raise Exception("Only one CS_NBACK entry in database!")
        condOffset= accOffset= blockOffset= trialOffset= None
        for i in xrange(len(dbCols)):
            if dbCols[i]=='COND': condOffset= i
            if dbCols[i]=='ACC': accOffset= i
            if dbCols[i]=='BLOCK': blockOffset= i
            if dbCols[i]=='TRIAL': trialOffset= i
        if condOffset==None or accOffset==None or blockOffset==None \
               or trialOffset==None:
            raise Exception("Can't find expected columns in CS_NBACK table!")

        trialOrder= []
        blockTrials= []
        trialsPerBlock= 15
##         for row in dbTable:
##             # The order of trials depends on getting the row list back
##             # from the database in properly sorted order.
##             block= int(row[blockOffset])-1
##             trial= int(row[trialOffset])-1
##             cond= int(row[condOffset])-1 # number 0-2 rather than 1-3
##             if trial<2:
##                 blockTrials.append(("startBlock",1))
##             elif int(row[accOffset])!=0:
##                 blockTrials.append(("%dback"%cond,1)) 
##             else:
##                 blockTrials.append(("error",1))
##             if len(blockTrials)==trialsPerBlock:
##                 trialOrder.append(blockTrials)
##                 blockTrials= []
##         debugMessage("Table entries: <%s>"%trialOrder)

        # One block consists of 15 reps (trials) in the given condition
        for c in dbDict['NBACK']:
            trialOrder.append([("%cback"%c,15)]) 

##         for i in xrange(len(trialOrder)):
##             print "%d: %s"%(i,trialOrder[i])

        nSlices= getNSlices()
        slicePattern= funcSliceOrder(nSlices, getScannerType())
        sliceTime= getTR()/nSlices
        image= 0
        block= 0
        # BIRC files start with 2 TRs of disdaq images; those images
        # are not retained in datasets from the GE scanner at the MR Ctr
        if getScannerType()=='GE': # GE
            nDisdaqs= 0
        elif getScannerType()=='SIEMENS': # Siemens
            nDisdaqs= 2
        else:
            raise Exception("nbackV2SplitMethod: unknown scanner type <%s>"%
                            getScannerType())
        imagesThisAcq= 141+nDisdaqs
        blocksThisAcq= 6
        for acq in xrange(3):
            markRefMissing= ( pulseSeqHasReferenceImage() and (block != 0) )
            states= nbackV2BuildStateList(trialOrder[block:(block+blocksThisAcq)],\
                                          nDisdaqs)
            walkStateList( ofile, states, slicePattern, sliceTime, \
                           image, image+imagesThisAcq, 'NA', markRefMissing)
            image += imagesThisAcq
            block += blocksThisAcq

    except Exception,e:
        print "Caught exception: %s"%e
        raise Exception("Generation of NBACK V2 split failed!")
    
##############################
#
# ARROW_V1 Task Translation Methods
#
##############################

def arrowsV1IaiMethod( line ): return '2.0'

def arrowsV1NimageMethod( line ):
    global envDict, catDict, transDict, dbDict
    try:
        return int(dbDict['NIMAGES'])
    except:
        pass # NIMAGES field in DB is not valid
    pfileDict= dbDict['PFILES']
    nImages= 0
    for pfileName in pfileDict.keys():
        (name,imagesThisPfile)= pfileDict[pfileName]
        nImages += int(imagesThisPfile)
        # Disdaqs are included in the database of Pfile image counts
##      if getScannerType()=='SIEMENS':
##          nImages += len(pfileDict.keys())*2 # disdaqs
    return nImages

# For easy renaming of the various arrows states
## arrowV1FactorNames= ['ConIncon','Side','Stage']

## arrowV1StateNmDict= { 'NA':['NA','NA','NA'],
##                       'ConLFix1':['Con','L','Fix1'],
##                       'InconLFix1':['Incon','L','Fix1'],
##                       'ConRFix1':['Con','R','Fix1'],
##                       'InconRFix1':['Incon','R','Fix1'],
##                       'ConLStim':['Con','L','Stim'],
##                       'InconLStim':['Incon','L','Stim'],
##                       'ConRStim':['Con','R','Stim'],
##                       'InconRStim':['Incon','R','Stim'],
##                       'ConLFix2':['Con','L','Fix2'],
##                       'InconLFix2':['Incon','L','Fix2'],
##                       'ConRFix2':['Con','R','Fix2'],
##                       'InconRFix2':['Incon','R','Fix2'],
##                       'ConLCue1':['Con','L','Cue1'],
##                       'InconLCue1':['Incon','L','Cue1'],
##                       'ConRCue1':['Con','R','Cue1'],
##                       'InconRCue1':['Incon','R','Cue1'],
##                       'PreBlockFixation':['NA','NA','NA'],
##                       'PostBlockFixation':['NA','NA','NA'],
##                       'Instructions':['NA','NA','NA'],
### Entries to cope with an old error in which 'con' was replaced with
### 'cor' in some of the pickled state list files
##                       'CorLFix1':['Con','L','Fix1'],
##                       'IncorLFix1':['Incon','L','Fix1'],
##                       'CorRFix1':['Con','R','Fix1'],
##                       'IncorRFix1':['Incon','R','Fix1'],
##                       'CorLStim':['Con','L','Stim'],
##                       'IncorLStim':['Incon','L','Stim'],
##                       'CorRStim':['Con','R','Stim'],
##                       'IncorRStim':['Incon','R','Stim'],
##                       'CorLFix2':['Con','L','Fix2'],
##                       'IncorLFix2':['Incon','L','Fix2'],
##                       'CorRFix2':['Con','R','Fix2'],
##                       'IncorRFix2':['Incon','R','Fix2'],
##                       'CorLCue1':['Con','L','Cue1'],
##                       'IncorLCue1':['Incon','L','Cue1'],
##                       'CorRCue1':['Con','R','Cue1'],
##                       'IncorRCue1':['Incon','R','Cue1'] }


####
# This block is appropriate for auto_birc_10
# Don't forget to comment out the time lag for BOLD in arrowsV1BuildStateList
####
## arrowV1FactorNames= ['Factor']

## arrowV1StateNmDict= { 'NA':['NA'],
##                       'ConLFix1':['Valid'],
##                       'InconLFix1':['Valid'],
##                       'ConRFix1':['Valid'],
##                       'InconRFix1':['Valid'],
##                       'ConLStim':['Valid'],
##                       'InconLStim':['Valid'],
##                       'ConRStim':['Valid'],
##                       'InconRStim':['Valid'],
##                       'ConLFix2':['Valid'],
##                       'InconLFix2':['Valid'],
##                       'ConRFix2':['Valid'],
##                       'InconRFix2':['Valid'],
##                       'ConLCue1':['Valid'],
##                       'InconLCue1':['Valid'],
##                       'ConRCue1':['Valid'],
##                       'InconRCue1':['Valid'],
##                       'PreBlockFixation':['Valid'],
##                       'PostBlockFixation':['Valid'],
##                       'Instructions':['Valid'],
### Entries to cope with an old error in which 'con' was replaced with
### 'cor' in some of the pickled state list files
##                       'CorLFix1':['Valid'],
##                       'IncorLFix1':['Valid'],
##                       'CorRFix1':['Valid'],
##                       'IncorRFix1':['Valid'],
##                       'CorLStim':['Valid'],
##                       'IncorLStim':['Valid'],
##                       'CorRStim':['Valid'],
##                       'IncorRStim':['Valid'],
##                       'CorLFix2':['Valid'],
##                       'IncorLFix2':['Valid'],
##                       'CorRFix2':['Valid'],
##                       'IncorRFix2':['Valid'],
##                       'CorLCue1':['Valid'],
##                       'IncorLCue1':['Valid'],
##                       'CorRCue1':['Valid'],
##                       'IncorRCue1':['Valid'] }
####
# End block for auto_birc_11
####

## ####
## # This block is appropriate for auto_birc_11
## # Don't forget to include the 4.0 sec BOLD time lag in arrowsV1BuildStateList
## ####
## arrowV1FactorNames= ['ConIncon']

## arrowV1StateNmDict= { 'NA':['NA'],
##                       'ConLFix1':['NA'],
##                       'InconLFix1':['NA'],
##                       'ConRFix1':['NA'],
##                       'InconRFix1':['NA'],
##                       'ConLStim':['Con'],
##                       'InconLStim':['Incon'],
##                       'ConRStim':['Con'],
##                       'InconRStim':['Incon'],
##                       'ConLFix2':['Con'],
##                       'InconLFix2':['Incon'],
##                       'ConRFix2':['Con'],
##                       'InconRFix2':['Incon'],
##                       'ConLCue1':['NA'],
##                       'InconLCue1':['NA'],
##                       'ConRCue1':['NA'],
##                       'InconRCue1':['NA'],
##                       'PreBlockFixation':['NA'],
##                       'PostBlockFixation':['NA'],
##                       'Instructions':['NA'],
## # Entries to cope with an old error in which 'con' was replaced with
## # 'cor' in some of the pickled state list files
##                       'CorLFix1':['NA'],
##                       'IncorLFix1':['NA'],
##                       'CorRFix1':['NA'],
##                       'IncorRFix1':['NA'],
##                       'CorLStim':['Con'],
##                       'IncorLStim':['Incon'],
##                       'CorRStim':['Con'],
##                       'IncorRStim':['Incon'],
##                       'CorLFix2':['Con'],
##                       'IncorLFix2':['Incon'],
##                       'CorRFix2':['Con'],
##                       'IncorRFix2':['Incon'],
##                       'CorLCue1':['Con'],
##                       'IncorLCue1':['Incon'],
##                       'CorRCue1':['Con'],
##                       'IncorRCue1':['Incon'] }
## ####
## # End block for auto_birc_11
## ####

####
# This block is appropriate for event-related analysis separating Con fm Incon
# Don't forget to remove the 4.0 sec BOLD time lag in arrowsV1BuildStateList
####
arrowV1FactorNames= ['ConIncon']

arrowV1StateNmDict= { 'NA':['NA'],
                      'ConLFix1':['NA'],
                      'InconLFix1':['NA'],
                      'ConRFix1':['NA'],
                      'InconRFix1':['NA'],
                      'ConLStim':['ConStim'],
                      'InconLStim':['InconStim'],
                      'ConRStim':['ConStim'],
                      'InconRStim':['InconStim'],
                      'ConLFix2':['Con'],
                      'InconLFix2':['Incon'],
                      'ConRFix2':['Con'],
                      'InconRFix2':['Incon'],
                      'ConLCue1':['NA'],
                      'InconLCue1':['NA'],
                      'ConRCue1':['NA'],
                      'InconRCue1':['NA'],
                      'PreBlockFixation':['NA'],
                      'PostBlockFixation':['NA'],
                      'Instructions':['NA'],
# Entries to cope with an old error in which 'con' was replaced with
# 'cor' in some of the pickled state list files
                      'CorLFix1':['NA'],
                      'IncorLFix1':['NA'],
                      'CorRFix1':['NA'],
                      'IncorRFix1':['NA'],
                      'CorLStim':['ConStim'],
                      'IncorLStim':['InconStim'],
                      'CorRStim':['ConStim'],
                      'IncorRStim':['InconStim'],
                      'CorLFix2':['Con'],
                      'IncorLFix2':['Incon'],
                      'CorRFix2':['Con'],
                      'IncorRFix2':['Incon'],
                      'CorLCue1':['Con'],
                      'IncorLCue1':['Incon'],
                      'CorRCue1':['Con'],
                      'IncorRCue1':['Incon'] }
####
# End block for event-related analysis separating Con from Incon
####

# This order determines the Fiasco condition number order for the arrows task
arrowV1StateNmOrder= ['ConLFix1',
                      'InconLFix1',
                      'ConRFix1',
                      'InconRFix1',
                      'ConLStim',
                      'InconLStim',
                      'ConRStim',
                      'InconRStim',
                      'ConLFix2',
                      'InconLFix2',
                      'ConRFix2',
                      'InconRFix2',
                      'ConLCue1',
                      'InconLCue1',
                      'ConRCue1',
                      'InconRCue1',
                      'PreBlockFixation',
                      'PostBlockFixation',
                      'Instructions']

def arrowsV1FilterOneRawState(rawName,durationMSec,startMSec,acc=1):
    duration= 0.001*durationMSec # convert from milliseconds to seconds
    start= 0.001*startMSec       # convert from milliseconds to seconds
    if (acc):
        lvls= arrowV1StateNmDict[rawName]
        nm= lvls[0]
        for lvl in lvls[1:]:
            nm= "%s %s"%(nm, lvl)
        return (nm, duration, None)
    else:
        naString= len(arrowV1FactorNames)*'NA '
        return (naString, duration, "subject task error")
    

def arrowsV1BuildStateList(rawStateList,nDisdaqs):
    states= []
    start= 0.0
    naString= len(arrowV1FactorNames)*'NA '
    # lag to compensate for BOLD
    #start= appendStateTuple( states, naString, 4.0, start )
    # Cover disdaqs with NA
    if nDisdaqs != 0:
        start= appendStateTuple( states, naString, nDisdaqs*getTR(),
                                 start, "disdaq" )

    # Beyond this point, everything is in the supplied raw state list
    firstState= 1
    for s in rawStateList:
        name, duration, comment = apply(arrowsV1FilterOneRawState, s)
        if firstState:
            if comment:
                comment += " (first state of acq)"
            else:
                comment= "(first state of acq)"
            firstState= 0
        start= appendStateTuple( states, name, duration, start, comment )
    return states

def arrowsV1WriteSplitHeader(ofile):
        nfactors= len(arrowV1FactorNames)
        ofile.write("%d\n"%nfactors)
        for n in arrowV1FactorNames:
            ofile.write("%s "%n)
        ofile.write("\n")
        for stateName in arrowV1StateNmOrder:
            stateParts= arrowV1StateNmDict[stateName]
            s= "%s"%stateParts[0]
            for w in stateParts[1:]:
                s= "%s %s"%(s,w)
            ofile.write("%s 0 0\n"%s)

def arrowsV1BuildStateListFromDB(acq):
    global envDict, catDict, transDict, dbDict
    # BIRC files start with 2 TRs of disdaq images; those images
    # are not retained in datasets from the GE scanner at the MR Ctr
    if getScannerType()=='GE': # GE
        nDisdaqs= 0
    elif getScannerType()=='SIEMENS': # Siemens
        nDisdaqs= 2
    else:
        raise Exception("arrowsV1SplitMethod: unknown scanner type <%s>"%
                        getScannerType())

    # Fetch the raw EPrime state information, already pickled
    # for our convenience by a preprocessing step.d
    stateFileName= "statelist_arrows_subj%03d_run%1d_acq%1d.pkl"%\
                   (int(transDict["subj"]),
                    int(transDict["runnum"]),
                    acq+1)
    fullStateFileName= os.path.join(getBehaviorDir(),
                                    stateFileName)
    if not os.access(fullStateFileName,os.R_OK):
        raise Exception("Cannot find pickled behavioral record <%s>!"%\
                        fullStateFileName)
    f= open(fullStateFileName,"rb")
    rawStatesThisAcq= cPickle.load(f)
    f.close()

    # Adjust the raw state list for the scanner
    states= arrowsV1BuildStateList(rawStatesThisAcq, nDisdaqs)
    return states

def arrowsV1SplitMethod( ofile, line ):
    global envDict, catDict, transDict, dbDict
    try:
        arrowsV1WriteSplitHeader(ofile)
        nSlices= getNSlices()
        slicePattern= funcSliceOrder(nSlices, getScannerType())
        sliceTime= getTR()/nSlices
        image= 0
        block= 0
        pfileDict= dbDict['PFILES']
        pfiles= sortByPfileOrder(pfileDict.keys())
        nAcqs= len(pfiles)
        run= transDict['runnum']
        for acq in xrange(nAcqs):
            states= arrowsV1BuildStateListFromDB(acq)

            (pfileName,imagesString)= pfileDict[pfiles[acq]]
            imagesThisAcq= int(imagesString) # includes nDisdaqs
            markRefMissing= (pulseSeqHasReferenceImage() and (acq != 0))
            walkStateList( ofile, states, slicePattern, sliceTime, \
                           image, image+imagesThisAcq,
                           len(arrowV1FactorNames)*'NA ',
                           markRefMissing )
            image += imagesThisAcq
    except Exception,e:
        print "Caught exception: %s"%e
        raise Exception("Generation of ARROWS_V1 split failed!")


##############################
#
# ARROW_V2 Task Translation Methods
#
##############################

def arrowsV2IaiMethod( line ): return '2.0'

def arrowsV2NimageMethod( line ):
    global envDict, catDict, transDict, dbDict
    try:
        return int(dbDict['NIMAGES'])
    except:
        pass # NIMAGES field in DB is not valid
    pfileDict= dbDict['PFILES']
    nImages= 0
    for pfileName in pfileDict.keys():
        (name,imagesThisPfile)= pfileDict[pfileName]
        nImages += int(imagesThisPfile)
        # Disdaqs are included in the database of Pfile image counts
##      if getScannerType()=='SIEMENS':
##          nImages += len(pfileDict.keys())*2 # disdaqs
    return nImages

# For easy renaming of the various arrows_v2 states
arrowV2StateNmDict= { 'NA':'NA',
                    'cue_Red':'red',
                    'cue_Green':'green',
                    'fix1_Red':'red',
                    'fix1_Green':'green',
                    'stim_Red_Left':'red',
                    'stim_Red_Right':'red',
                    'stim_Green_Left':'green',
                    'stim_Green_Right':'green',
                    'fix2_Red_Left':'red',
                    'fix2_Red_Right':'red',
                    'fix2_Green_Left':'green',
                    'fix2_Green_Right':'green',
                    'goodbye':'NA'}

# This order determines the Fiasco condition number order for the arrows task
arrowV2StateNmOrder= [ "cue_Green","cue_Red",
                       "fix1_Green", "fix1_Red",
                       "stim_Green_Right","stim_Green_Left",
                       "stim_Red_Right","stim_Red_Left",
                       "fix2_Green_Right", "fix2_Green_Left",
                       "fix2_Red_Right","fix2_Red_Left" ]

def arrowsV2BuildStateList(trialOrder,nDisdaqs):
    states= []
    start= 0.0
    # lag to compensate for BOLD
    start= appendStateTuple( states, 'NA', 0.0, start )
    # Cover disdaqs with NA
    if nDisdaqs != 0:
        start= appendStateTuple( states, 'NA',
                                 nDisdaqs*getTR(), start )

    for blockTrials in trialOrder:
        debugMessage("arrowsV2BuildStateList %s: building block:<%s>"%\
                     ("ARROWS_V2",blockTrials))
        for (cond,reps) in blockTrials:
            if cond=='NA':
                for trial in range(reps):
                    start= appendStateTuple( states,
                                             arrowV2StateNmDict['NA'], \
                                             20.0, start )
            else:
                if cond[0:-1]=='COR': clr= "Green"
                elif cond[0:-1]=='INCOR': clr= "Red"
                else: raise Exception(("arrowsV2BuildStateList: unknown cond %s"\
                                       +" in ARROWS_V2 task")%cond)
                if cond[-1]=="L": dir= "Left"
                elif cond[-1]=="R": dir= "Right"
                else: raise Exception(("arrowsV2BuildStateList: unknown cond %s"\
                                       +" in ARROWS_V2 task")%cond)
                for trial in xrange(reps):
                    # Cue
                    start= appendStateTuple( states,
                                             arrowV2StateNmDict['cue_'+clr],
                                             0.5, start )
                    # fixation
                    start= appendStateTuple( states,
                                             arrowV2StateNmDict['fix1_'+clr],
                                             7.5, start )
                    # Stim
                    start= appendStateTuple( states,
                                             arrowV2StateNmDict["stim_%s_%s"%(clr,dir)],
                                             0.5, start )
                    # fixation
                    start= appendStateTuple( states,
                                             arrowV2StateNmDict['fix2_%s_%s'%(clr,dir)],
                                             11.5, start )
        # post-block fixation
        start= appendStateTuple( states, arrowV2StateNmDict['goodbye'],
                                 1.0, start )

def arrowsV2SplitMethod( ofile, line ):
    global envDict, catDict, transDict, dbDict
    try:
        ofile.write("1\n")
        ofile.write("Condition\n")
        for stateName in arrowV2StateNmOrder:
            ofile.write("%s 0 0\n"%arrowV2StateNmDict[stateName])
        dbTable= dbDict['CS_ARROW']
        if dbTable==None:
            raise Exception("No CS_ARROW entry in database!")
        dbCols= dbDict['CS_ARROW_columns']
        if dbCols==None:
            raise Exception("Only one CS_ARROW entry in database!")
        condOffset= accOffset= acqOffset= blockOffset= None
        for i in xrange(len(dbCols)):
            if dbCols[i]=='COND': condOffset= i
            if dbCols[i]=='ACC': accOffset= i
            if dbCols[i]=='PFILE': acqOffset= i
            if dbCols[i]=='BLOCK': blockOffset= i
        if condOffset==None or accOffset==None or acqOffset==None or \
           blockOffset==None:
            raise Exception("Can't find expected columns in CS_ARROW table!")
        lutDict= dbDict['LU_COND_ARW']
        if lutDict==None:
            raise Exception("No LU_COND_ARW entry in database!")
        lutCols= dbDict['LU_COND_ARW_columns']
        if lutCols==None:
            raise Exception("Only one LU_COND_ARW entry in database!")
        for lutOffset in xrange(len(lutCols)):
            if lutCols[lutOffset]=='LABEL':
                break
        trialOrderDict= {}
        for row in dbTable:
            block= int(row[blockOffset])-1
            acq= int(row[acqOffset])-1
            if int(row[accOffset])!=0:
                condTuple= ("%s"%lutDict[row[condOffset]][lutOffset-1],
                            1) # one repetition of this
            else:
                condTuple= ("NA", # subject got it wrong; mark missing
                            1)    # one repetition of this
            if trialOrderDict.has_key(acq):
                blockDict= trialOrderDict[acq]
            else:
                blockDict= {}
                trialOrderDict[acq]= blockDict
            if blockDict.has_key(block):
                trialList= blockDict[block]
            else:
                trialList= []
                blockDict[block]= trialList
            # The order of trials depends on getting the row list back
            # from the database in properly sorted order.
            trialList.append(condTuple)
        debugMessage("Table entries: <%s>"%trialOrderDict)

        nSlices= getNSlices()
        slicePattern= funcSliceOrder(nSlices, getScannerType())
        sliceTime= getTR()/nSlices
        image= 0
        block= 0
        # BIRC files start with 3 TRs of disdaq images; those images
        # are not retained in datasets from the GE scanner at the MR Ctr
        if getScannerType()=='GE': # GE
            nDisdaqs= 0
        elif getScannerType()=='SIEMENS': # Siemens
            nDisdaqs= 4
        else:
            raise Exception("arrowsV2SplitMethod: unknown scanner type <%s>"%
                            getScannerType())
        pfileDict= dbDict['PFILES']
        pfiles= sortByPfileOrder(pfileDict.keys())
        for acq in xrange(len(pfiles)):
            (pfileName,imagesString)= pfileDict[pfiles[acq]]
            blockDict= trialOrderDict[acq]
            blocksThisAcq= len(blockDict)
            # pfileDict image counts come from the database, and
            # include any disdaq's
            imagesThisAcq= int(imagesString) # includes nDisdaqs
            markRefMissing= (pulseSeqHasReferenceImage() and (acq != 0))
            trialList= []
            for block in xrange(blocksThisAcq):
                trialList.append(blockDict[block])
            states= arrowsV2BuildStateList(trialList, nDisdaqs)
            walkStateList( ofile, states, slicePattern, sliceTime, \
                           image, image+imagesThisAcq, 'NA',
                           markRefMissing )
            image += imagesThisAcq
    except Exception,e:
        print "Caught exception: %s"%e
        raise Exception("Generation of ARROWS_V2 split failed!")
    
##############################
#
# STROOP Task Translation Methods
#
##############################

def stroopIaiMethod( line ): return '1.5'

def stroopNimageMethod( line ):
    global envDict, catDict, transDict, dbDict
    try:
        return int(dbDict['NIMAGES'])
    except:
        pass # NIMAGES field in DB is not valid
    nImages= 301*len(dbDict['CS_STROOP'])/150
    if getScannerType()=='SIEMENS':
        raise Exception("Number of disdaqs for STROOP task not implemented!")
    return nImages

def stroopSplitMethod( ofile, line ):
    try:
        ofile.write("1 bynumimages\n")
        ofile.write("Condition\n")
        ofile.write("C 0\n")
        ofile.write("I 0\n")
        bob = dbDict['CS_STROOP']
        if bob==None:
            raise Exception("No CS_STROOP entry in database!")
        cols = dbDict['CS_STROOP_columns']
        for offset in xrange(len(cols)):
            if cols[offset]=='COND':
                break
        debugMessage("found offset %d"%offset)
        #ofile.write("C 1\n")
        for i in xrange(len(bob)):
            if i == 150 or i == 300:
                ofile.write("NA 2\n")
            else:
                if int(bob[i][offset]) == 1:
                    ofile.write("C 2\n")
                else:
                    ofile.write("I 2\n")
    except:
        raise Exception("Generation of STROOP split failed!")

##############################
#
# GE Scanner Translation Methods
#
##############################

def geRunTypeMethod( line ): return "spiral"

def geInplInputMethod( line ):
    global envDict, catDict, transDict, dbDict
    id= int(transDict['subj']);
    dirBase= "/home/YODA/roushre/DATA/Pfiles/%03d/struct_%s/"%\
             (id,transDict['runnum'])
    for dirnum in range(1,20):
        dirPath= os.path.join(dirBase,"%03d"%dirnum)
        if os.access(dirPath,os.R_OK):
            count= countImageFiles(dirPath)
            if count==int(dbDict['SLICE']):
                return "'%s'"%os.path.join(dirPath,'I.###')
    raise Exception("No suitable inplane image set in %s"%dirBase)

def geAnatInputMethod( line ):
    global envDict, catDict, transDict, dbDict
    id= int(transDict['subj']);
    dirBase= "/home/YODA/roushre/DATA/Pfiles/%03d/struct_%s/"%\
             (id,transDict['runnum'])
    for dirnum in range(1,20):
        dirPath= os.path.join(dirBase,"%03d"%dirnum)
        if os.access(dirPath,os.R_OK):
            count= countImageFiles(dirPath)
            if count>100 and count<150:
                return "'%s'"%os.path.join(dirPath,'I.###')
    raise Exception("No suitable anat image set in %s"%dirBase)

def geReaderInputMethod( line ):
    global envDict, catDict, transDict, dbDict
    # GE scanner produces Pfiles; this clause should
    # never happen because spiral runs have no F_READER_INPUT field.
    raise Exception("readerInputMethod: EPI input from GE is not implemented")

##############################
#
# Siemens Scanner Translation Methods
#
##############################

def siemensReaderOptsMethod( line ): return "'-autoscale_range 8000.0 -multi'"

def siemensRunTypeMethod( line ): return "epi"

def siemensInplInputMethod( line ):
    global envDict, catDict, transDict, dbDict
    id= int(transDict['subj']);
    dirBase= "/home/YODA/roushre/DATA/Pfiles/%03d/struct_%s/"%\
             (id,transDict['runnum'])
    for dirnum in range(1,25):
        dirPath= os.path.join(dirBase,"%02d"%dirnum)
        if os.access(dirPath,os.R_OK):
            count= countDicomFiles(dirPath)
            if count==int(dbDict['SLICE']):
                return "'%s'"%os.path.join(dirPath,'^.dcm')
    for dirnum in range(1,25):
        dirPath= os.path.join(dirBase,"D%02d"%dirnum)
        if os.access(dirPath,os.R_OK):
            count= countDicomFiles(dirPath)
            if count==int(dbDict['SLICE']):
                return "'%s'"%os.path.join(dirPath,'^.dcm')
    raise Exception("No suitable inplane image set in %s"%dirBase)

def siemensAnatInputMethod( line ):
    global envDict, catDict, transDict, dbDict
    id= int(transDict['subj']);
    dirBase= "/home/YODA/roushre/DATA/Pfiles/%03d/struct_%s/"%\
             (id,transDict['runnum'])
    for dirnum in range(1,25):
        dirPath= os.path.join(dirBase,"%02d"%dirnum)
        if os.access(dirPath,os.R_OK):
            count= countDicomFiles(dirPath)
            if count>150 and count<200:
                return "'%s'"%os.path.join(dirPath,'^.dcm')
    for dirnum in range(1,25):
        dirPath= os.path.join(dirBase,"D%02d"%dirnum)
        if os.access(dirPath,os.R_OK):
            count= countDicomFiles(dirPath)
            if count>150 and count<200:
                return "'%s'"%os.path.join(dirPath,'^.dcm')
    raise Exception("No suitable anat image set in %s"%dirBase)

def siemensReaderInputMethod( line ):
    global envDict, catDict, transDict, dbDict
    parentDir= envDict['F_DIR'](line)
    # Siemens scanner produces directories full of Dicoms
    subdirList= os.listdir(parentDir)
    usefulSubdirList= []
    for dir in subdirList:
        if siemensDicomDirRegex.match(dir):
            usefulSubdirList.append(dir)
    if len(usefulSubdirList)==0:
        raise Exception("readerInputMethod: no dirs of dcm files below %s"%\
                        parentDir)
    usefulSubdirList.sort()
    spec= '" '
    for dir in usefulSubdirList:
        spec += "$F_DIR/%s/^.dcm "%dir
    spec += '"'
    return spec


##############################
#
# Translation Tables
#
##############################

baseEnvDict= {"F_HEADER":headerMethod, "F_DIR":dirMethod, "F_REFP":refpMethod,
              "F_DESCRIPTION":descriptionMethod, 'F_SUBJ_SEX':genderMethod,
              'F_SUBJ_1DIAG':diag1Method, 'F_SUBJ_AGE':ageMethod,
              'F_XVOXEL':xvoxelMethod, 'F_YVOXEL':yvoxelMethod,
              'F_ZVOXEL':zvoxelMethod, 'F_NSLICE':nsliceMethod,
              'F_IAI':None, 'F_ESTIREG_FIXED':estiregfixedMethod,
              'F_LOWERT':lowertMethod, 'F_UPPERT':uppertMethod,
              'F_UPPERF':upperfMethod, 'F_LOWERQ':lowerqMethod,
              'F_NIMAGE':None,
#              'F_TEMP':tmpdirMethod,
              'F_EST3D_ALIGN':est3dalignMethod,
              'F_CRG_INPL_INPUT':None,
              'F_CRG_INPL_OPTS':inplOptsMethod,
              'F_CRG_ANAT_INPUT':None,
              'F_CRG_ANAT_OPTS':anatOptsMethod,
              'F_RUN_TYPE':None,
              'F_READER_INPUT':None,
              'F_READER_DIMS':readerDimsMethod,
              'F_READER_OPTS':baseReaderOptsMethod
              }

baseCatDict= {"$F_SPLIT_FILE":None}

baseHookDict= {"BuildStateList":None}

taskEnvDict= { "NBACK_V1":{'F_IAI':nbackV1IaiMethod,
                           'F_NIMAGE':nbackV1NimageMethod,
                           },
               "NBACK_V2":{'F_IAI':nbackV2IaiMethod,
                           'F_NIMAGE':nbackV2NimageMethod,
                           },
               "ARROWS_V1":{'F_IAI':arrowsV1IaiMethod,
                            'F_NIMAGE':arrowsV1NimageMethod,
                            },
               "ARROWS_V2":{'F_IAI':arrowsV2IaiMethod,
                            'F_NIMAGE':arrowsV2NimageMethod,
                            },
               "STROOP":{'F_IAI':stroopIaiMethod,
                         'F_NIMAGE':stroopNimageMethod,
                         },
               }
taskCatDict= { "NBACK_V1":{ "$F_SPLIT_FILE":nbackV1SplitMethod },
               "NBACK_V2":{ "$F_SPLIT_FILE":nbackV2SplitMethod },
               "ARROWS_V1":{ "$F_SPLIT_FILE":arrowsV1SplitMethod },
               "ARROWS_V2":{ "$F_SPLIT_FILE":arrowsV2SplitMethod },
               "STROOP":{ "$F_SPLIT_FILE":stroopSplitMethod },
               }

taskHookDict= { "NBACK_V1":{},
                "NBACK_V2":{},
                "ARROWS_V1":{ "BuildStateList":arrowsV1BuildStateListFromDB },
                "ARROWS_V2":{},
                "STROOP":{},
               }

scannerEnvDict= { "SIEMENS":{'F_RUN_TYPE':siemensRunTypeMethod,
                             'F_READER_OPTS':siemensReaderOptsMethod,
                             'F_CRG_INPL_INPUT':siemensInplInputMethod,
                             'F_CRG_ANAT_INPUT':siemensAnatInputMethod,
                             'F_READER_INPUT':siemensReaderInputMethod,
                             },
                  "GE":{'F_RUN_TYPE':geRunTypeMethod,
                        'F_CRG_INPL_INPUT':geInplInputMethod,
                        'F_CRG_ANAT_INPUT':geAnatInputMethod,
                        'F_READER_INPUT':geReaderInputMethod,
                        }
                  }

scannerCatDict= { "SIEMENS":{},
                  "GE":{}
                  }

scannerHookDict= { "SIEMENS":{},
                   "GE":{}
                   }

##############################
#
# End translation Tables
#
##############################
    
def checkInitialDefs( dict ):
    errStr= ""
    errFlag= 0
    if not dict.has_key('subj'):
        errFlag= 1
        errStr= errStr + " subj=subjectID (int)"
    if not dict.has_key('task'):
        errFlag= 1
        errStr= errStr + " task=taskname (string)"
    if not dict.has_key('runnum'):
        errFlag= 1
        errStr= errStr + " runnum=N (int)"
    if not dict.has_key('tag'):
        errFlag= 1
        errStr= errStr + " tag=tagString (string)"
    if not dict.has_key('db_uname'):
        errFlag= 1
        errStr= errStr + " db_uname=databaseUsername (string)"
    if not dict.has_key('db_passwd'):
        errFlag= 1
        errStr= errStr + " db_passwd=databasePassword (string)"
    if errFlag:
        Message("The following required key-value pairs are missing:")
        Message(errStr)
        sys.exit("Not enough initial information given!")

def overlayTaskSpecificMethods( baseEnvDict, baseCatDict, baseHookDict, task ):
    global taskEnvDict, taskCatDict;
    envDict= baseEnvDict.copy()
    if taskEnvDict.has_key(task):
        dict= taskEnvDict[task]
        for k in dict.keys():
            envDict[k]= dict[k]
    else:
        raise Exception("No task-specific environment variable methods for task %s!"%task)
    catDict= baseCatDict.copy()
    if taskCatDict.has_key(task):
        dict= taskCatDict[task]
        for k in dict.keys():
            catDict[k]= dict[k]
    else:
        raise Exception("No task-specific cat methods for task %s!"%task)
    hookDict= baseHookDict.copy()
    if taskHookDict.has_key(task):
        dict= taskHookDict[task]
        for k in dict.keys():
            hookDict[k]= dict[k]
    else:
        raise Exception("No task-specific hook methods for task %s!"%task)
    return (envDict,catDict,hookDict)

def overlayScannerSpecificMethods( baseEnvDict, baseCatDict, baseHookDict,
                                   scanner ):
    global scannerEnvDict, scannerCatDict;
    envDict= baseEnvDict.copy()
    if scannerEnvDict.has_key(scanner):
        dict= scannerEnvDict[scanner]
        for k in dict.keys():
            envDict[k]= dict[k]
    else:
        raise Exception("No scanner-specific environment variable methods for scanner %s!"%scanner)
    catDict= baseCatDict.copy()
    if scannerCatDict.has_key(scanner):
        dict= scannerCatDict[scanner]
        for k in dict.keys():
            catDict[k]= dict[k]
    else:
        raise Exception("No scanner-specific cat methods for scanner %s!"%task)
    hookDict= baseHookDict.copy()
    if scannerHookDict.has_key(scanner):
        dict= scannerHookDict[scanner]
        for k in dict.keys():
            hookDict[k]= dict[k]
    else:
        raise Exception("No scanner-specific hook methods for scanner %s!"%task)
    return (envDict,catDict,hookDict)

def deriveFullTaskName( dbDict, transDict ):
    if string.upper(transDict['task'])=='NBACK':
        if int(dbDict['DID_NBK_V1']): return 'NBACK_V1'
        else: return 'NBACK_V2'
    elif string.upper(transDict['task'])=='ARROW':
        if int(dbDict['DID_ARROW_V2']): return 'ARROWS_V2'
        else: return 'ARROWS_V1'
    else: return string.upper(transDict['task'])
    
def deriveScannerType( dbDict, transDict ):
    if int(dbDict['RF'])==1: return 'GE'
    elif int(dbDict['RF'])==2: return 'SIEMENS'
    else:
        raise Exception("Unknown RF value %s; can't guess scanner type!"%\
                    dbDict['RF'])

def populateMethodTables( localTransDict ):
    global transDict, dbDict, envDict, catDict, hookDict
    transDict= localTransDict
    # Check environment for database password if it's not already set
    if not transDict.has_key('db_passwd'):
        if os.environ.has_key('F_DB_PASSWD'):
            transDict['db_passwd']= os.environ['F_DB_PASSWD']

    # Do I have the initial definitions I need?
    checkInitialDefs(transDict)

    dbAccessPrep()
    fillDbDict( dbDict )
    dbAccessFinalize()

    # Check DB for inconsistencies
    if not dbDict.has_key("GENDER"):
        sys.exit("No subject id %s found in database!"%transDict["subj"])
    if not dbDict.has_key("SDATE"):
        sys.exit("No scan %s for subject id %s found in database!"%\
                 (transDict["runnum"],transDict["subj"]))
    if dbDict.has_key('DID_NBK_V1') and dbDict.has_key('DID_NBK_V2'):
        if int(dbDict['DID_NBK_V1']) and int(dbDict['DID_NBK_V2']):
            sys.exit("Inconsistent version of NBACK task for subject id %s!"%\
                     transDict["subj"])
    if string.upper(transDict['task'])=='NBACK' and dbDict.has_key('DID_NBK_V1') \
           and dbDict.has_key('DID_NBK_V2'):
        if not int(dbDict['DID_NBK_V1']) \
               and not int(dbDict['DID_NBK_V2']):
            sys.exit("Inconsistent version of NBACK task for subject id %s!"%\
                     transDict["subj"])
    # patch: db dosn't contain DID_ARROW_V2, and no subject has it set
    dbDict['DID_ARROW_V2']= 0
    if dbDict.has_key('DID_ARROW_V1') and dbDict.has_key('DID_ARROW_V2'):
        if int(dbDict['DID_ARROW_V1']) and int(dbDict['DID_ARROW_V2']):
            sys.exit("Inconsistent version of ARROW task for subject id %s!"%\
                     transDict["subj"])
    if string.upper(transDict['task'])=='ARROW' \
           and dbDict.has_key('DID_ARROW_V1') \
           and dbDict.has_key('DID_ARROW_V2'):
        if not int(dbDict['DID_ARROW_V1']) \
               and not int(dbDict['DID_ARROW_V2']):
            sys.exit("Inconsistent version of ARROW task for subject id %s!"%\
                     transDict["subj"])

    # Make some convenient derived names
    dbDict['FULL_TASK_NAME']= deriveFullTaskName(dbDict,transDict)
    dbDict['SCANNER_TYPE']= deriveScannerType(dbDict,transDict)

    # Modify actions for the specific scanner and task
    (envDict, catDict, hookDict)= \
              overlayScannerSpecificMethods( baseEnvDict, baseCatDict,
                                             baseHookDict, getScannerType() )
    (envDict, catDict, hookDict)= \
              overlayTaskSpecificMethods(envDict, catDict,
                                         baseHookDict,
                                         dbDict['FULL_TASK_NAME'])

    return (dbDict,envDict,catDict,hookDict)

##############################
#
# Main
#
##############################

def main(argv=None):
    # These dictionaries are shared by everyone
    global envDict, catDict, dbDict, hookDict
    # transDict is global to everyone else, spread via populateMethodTables()
    transDict= {}
    
    if argv is None:
        argv= sys.argv

    # Check for "-help"
    if len(argv)>1:
        if argv[1] == "-help":
            if len(argv)>2:
                os.system( "scripthelp %s %s"%(argv[0],argv[2]) )
            else:
                os.system( "scripthelp %s"%argv[0] )
            sys.exit()

    try:
        (opts,pargs) = getopt.getopt(argv[1:],"vd",["src="])
    except:
        errorMessage("%s: Invalid command line parameter" % argv[0])
        describeSelf()
        sys.exit(1)

    verboseMessage(idString)

    srcDir= None
    
    # Parse args
    for a,b in opts:
        if a=="-v":
            setVerbose(1)
        if a=="-d":
            setDebug(1)
        if a=="--src":
            srcDir= b

    outName= None
    for arg in pargs:
        if string.find(arg,"=") >= 0:
            (key,val) = map(string.strip,string.split(arg,"="))
            transDict[key]= val
        else:
            if outName == None:
                outName= arg
            else:
                sys.exit("%s: command line contains repeated output filenames"
                         %argv[0])

    if outName==None:
        sys.exit("%s: name of file to customize was not given!"%argv[0])

    # Fill out the working dictionaries with appropriate methods
    (dbDict,envDict,catDict,hookDict)= populateMethodTables(transDict)

    if srcDir==None:
        if not os.environ.has_key('FIASCO'):
            sys.exit("No source directory given and FIASCO environment variable is not set!")
        srcDir= os.environ['FIASCO']

    inName= os.path.join(srcDir,outName);

    debugMessage("Translating %s to %s"%(inName,outName))
    if not os.path.isfile(inName):
        sys.exit("%s does not exist!"%inName)

    f= file(inName)
    canExecute= os.access(inName,os.X_OK)
    inLines= f.readlines()
    f.close()

    ofile= file(outName,"w")
    skipMode= 0
    skipEnd= None
    try:
        for line in inLines:
            trimmedLine= string.strip(line)
            if skipMode:
                if trimmedLine[0:len(skipEnd)]==skipEnd:
                    ofile.write(line)
                    skipMode= 0
                else:
                    continue
            elif len(trimmedLine)==0:
                ofile.write(line)
            elif trimmedLine[0]=='#':
                ofile.write(line)
            else:
                words= string.split(trimmedLine)
                if words[0]=='setenv':
                    debugMessage("found setenv of <%s>"%words[1])
                    if envDict.has_key(words[1]):
                        debugMessage("found translation method for setenv <%s>"%\
                                     words[1])
                        mthd= envDict[words[1]]
                        if mthd==None:
                            raise Exception("Missing needed method to customize <%s>!"%\
                                            trimmedLine)
                        ofile.write("setenv %s %s\n"%\
                                    (words[1], mthd(trimmedLine)))
                    else:
                        ofile.write(line)
                elif words[0]=='cat' and words[1]=='>' and words[3]=='<<':
                    if catDict.has_key(words[2]):
                        mthd= catDict[words[2]]
                        if mthd==None:
                            raise Exception("Missing needed method to customize %s!"%\
                                            words[2])
                        ofile.write(line)
                        mthd(ofile,line)
                        skipMode= 1
                        skipEnd= words[4]
                    else:
                        ofile.write(line)
                else:
                    ofile.write(line)
    except Exception, e:
        sys.exit("Customization failed: %s"%e)

    ofile.close()

    if canExecute:
        verboseMessage("Making %s executable by anyone!"%outName)
        safeRun("chmod +x %s"%outName)

############
# Main hook
############

if __name__=="__main__":
    main()



