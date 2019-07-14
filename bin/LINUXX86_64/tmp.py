import sys
import os
import string


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
    print("<%s> -> %d <%s> -> %d %d"%\
                 (roiName,roiNum,subRoiField,roiNum,subRoiNum))
##     debugMessage("<%s> -> %d <%s> -> %d %d"%\
##                  (roiName,roiNum,subRoiField,roiNum,subRoiNum))
    return (roiNum,subRoiNum,subRoiField)    

