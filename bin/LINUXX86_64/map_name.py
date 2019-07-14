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
# The purpose of this script is to provide a level of indirection
# between session, group, and study names and the directory structure in
# which they are stored.
#

import sys
import os
import string
import getopt

# Maximum mapping depth
maxDepth= 20;

# What characters delimit strings?
delimiters = ['$','@']

# Dictionary of group directories; starts empty
defaultDict= { };
mapDict= {"-":defaultDict };

# Mode flags
verbose= 0;

# Default target
target= "path";

def trim( str ):
    where= len(str)-1;
    while str[where]=='\012':
        str = str[0:where]
        where = where-1
        if where<0:
            break
    string.strip(str)
    return str

def addToDict( dict, str ):
    if len(str)<1:
	return;
    brk= string.find(str,":");
    key= string.strip(str[0:brk]);
    val= string.strip(str[brk+1:]);
    dict[key]= val;

def addTriple( dictKey, key, value ):
    global mapDict
    if mapDict.has_key(dictKey):
	mapDict[dictKey][key]= value;
    else:
	newDict= { };
	newDict[key]= value;
	mapDict[dictKey]= newDict;
        if verbose:
            print "Added new dictionary <%s> to mapDict"%dictKey

def loadByColumn( lines, istart, n ):
    global mapDict
    dicts= string.split(lines[istart]);
    i= istart+1;
    while i < istart+n:
	if len(lines[i])>0:
	    vals= string.split(lines[i]);
	    key= vals[0];
	    keylist= [];
	    for junk in vals: keylist.append(key);
	    map(addTriple,keylist,dicts,vals);
	i= i+1;

def loadMap( fname ):
    global mapDict
    if verbose :
	print "Loading dictionary file <%s>"%fname
    f= open(fname)
    lines= map(trim,f.readlines())
    f.close();
    i= 0;
    while i < len(lines):
	if len(string.strip(lines[i])) > 0:
	    [mapAttr,nString]= string.split(lines[i]);
	    i= i+1;
	    if nString=="*":
		n= len(lines)-i;
	    else:
		n= string.atoi(nString);
	    if string.lower(mapAttr)=="bycolumn":
		loadByColumn(lines,i,n);
		i= i+n;
	    else:
		thisDict= { }
		if mapDict.has_key(mapAttr):
		    thisDict= mapDict[mapAttr];
		for j in xrange(0,n):
		    addToDict( thisDict, lines[i] );
		    i= i+1;
		mapDict[mapAttr]= thisDict;
                if verbose:
                    print "Added new dictionary <%s> to mapDict"%mapAttr
	else:
	    i= i+1;

def loadEnv( str ):
    # We need to walk the string backwards, so that things closer
    # to the front override those which come later.
    str= string.strip(str);
    if verbose :
	print "Loading environment %s"%str;
    while len(str)>0 :
	brk= string.rfind( str, ":" );
	if brk > 0:
	    thismap= string.strip( str[brk+1:] );
	    str= string.strip( str[0:brk] );
	    loadMap( thismap );
	else:
	    loadMap( str );
	    str= "";
    
def safeMap( dict, term ):
    if verbose:
	print "safely mapping <%s>"%term;
    if len(term)==0:
	return term;
    if dict.has_key(term):
	return dict[term];
    else:
	return term;

def expandDictionary(dict,nameOfAnotherDict):
    global mapDict
    if mapDict.has_key(nameOfAnotherDict):
        if verbose:
            print "Adding dictionary for <%s>"%nameOfAnotherDict
        for (key,val) in mapDict[nameOfAnotherDict].items():
            dict[key]= val

def maybeExpandDictionary(dict,termKey,term):
    global mapDict
    if mapDict.has_key(term):
        subDict= mapDict[term]
        if subDict.has_key(termKey):
            if subDict[termKey]==term:
                if verbose:
                    print "Adding dictionary for <%s>:<%s>"%(termKey,term)
                for (key,val) in subDict.items():
                    if key != termKey:
                        dict[key]= val

def hasDelimPairs( term ):
    for c in delimiters:
        brk= string.find(term,c)
        if brk>=0:
            brk2= string.find(term[brk+1:],c)
            if brk2>=0:
                return 1
    return 0

def defineTerm( dict, term ):
    if verbose:
        print "defining term <%s>"%term
    brk= string.find(term,"=");
    if brk > 0:
        key= string.strip(term[0:brk])
        val= string.strip(term[brk+1:])
	dict[key]= val
        maybeExpandDictionary(dict,key,val)
    else:
	dict["-"]= mapExpression(term);
        maybeExpandDictionary(dict,"-",val)

def recursiveMap( dict, term, lvl ):
    mappingPass= 1
    if lvl>maxDepth:
        raise RuntimeError,"Recursion too deep mapping %s"%term;
    if verbose:
        print "recursively mapping <%s> at level %d"%(term,lvl);
    result= term
    while hasDelimPairs(result):
        for c in delimiters:
            brk1= string.find( result, c );
            if brk1>=0:
                brk2= string.find(result[brk1+1:], c);
                if brk2<=0:
                    if lvl == 0:
                        raise RuntimeError,"Unmatched %s in %s"%(c,result);
                else:
                    brk2= brk1 + brk2 + 1; # relative to start of result
                    mapResult1= recursiveMap(dict,result[brk1+1:brk2],lvl+1)
                    mapResult2= recursiveMap(dict,result[brk2+1:],lvl+1)
                    maybeExpandDictionary(dict,result[brk1+1:brk2],mapResult1)
                    maybeExpandDictionary(dict,result[brk2+1:],mapResult2)
                    result= result[0:brk1]+ mapResult1 + mapResult2
        if verbose:
            print "level %d pass %d: <%s> has become <%s>"%\
                  (lvl,mappingPass,term,result)
        mappingPass= mappingPass+1
    if dict.has_key(result):
        result= recursiveMap(dict,dict[term],lvl+1);
        maybeExpandDictionary(dict,term,result)
    if verbose:
        print "mapped <%s> to <%s> at level %d"%(term,result,lvl)
    return result

def describeSelf():
    print """
Usage: %s [-v] [-m map] [-m map] [-m map] [-t target] [-d key1=val1] [-d key2=val2] ... mainkey

-v specifies verbose mode.
-m specifies the name of a map file (see below)
-t target sets the string to be translated; the default is "path"
-d key=value sets the value associated with the given key.  This
   key-value pair may in turn be overridden by a pair produced by mainkey.
mainkey is the key used to select a collection of key-value pairs 
   used in translating the target.  mainkey will usually be a subject ID.

An output string is generated based on map files in the map path
specified by the environment variable F_MAP_PATH and any maps 
specified on the command line.  F_MAP_PATH is a colon-separated
string of map file names.  Map files appearing in F_MAP_PATH 
supercede other map files appearing later in that string.  Maps
specified on the command line supercede those in F_MAP_PATH, with
maps appearing later on the command line superceding those which
appear earlier.

Each map file consists of a series of blocks in either of two formats.  
The first format is:

keyname n
term1 : mapofterm1
term2 : mapofterm2
...
termn : mapoftermn

The second format is:
bycolumn n
key1 key2 key3 ...
val11 val12 val13 ...
val21 val22 val23 ...
val31 val32 val33 ...
val41 val42 val43 ...
...
valn1 valn2 valn3 ...

In either case n can be an integer representing the number of additional
lines in the block, or the character "*" to indicate that the block extends
to the end of the file.  Each such block enters term:mapofterm pairs
into a dictionary with the given key name.  The first form enters all
of the given term:mapofterm pairs into a dictionary specified by keyname.
The second form treats each line separately, adding the pairs 
keyij:valij into a dictionary with the name keyi1, for all i and j.
Thus the first value in each row serves as both a dictionary name and as
the value associated with the first key, which implies that the values 
in the first column must be unique.

A dictionary is constructed by copying the dictionary with the name
"default", adding any key:value pairs from the command line, and
finally adding any definitions from the dictionary the name of which
matches mainkey.  The target is then repeatedly mapped through this
dictionary, replacing the whole target or any part delimited by '$' or '@'
characters with the corresponding value.  The string which ultimately
results is the output.

""" % sys.argv[0]

##############################
#
# Main
#
##############################

if len(sys.argv) == 1 :
    describeSelf()
    sys.exit(-1)

try:
    (opts,pargs) = getopt.getopt(sys.argv[1:],'-m:-v-d:-t:')
except:
    print "%s: Invalid command line parameter" % sys.argv[0]
    describeSelf();
    sys.exit(-1);

if len(pargs)!=1:
    print "%s: too few or too many main keys" % sys.argv[0]
    describeSelf();
    sys.exit(-1);

# Handle flags not effecting mappings
for a,b in opts:
    if a=="-v":
	# verbose mode on
	verbose= 1;
    if a=="-t":
	# set the target string
	target= b;

# Load the global translation maps
useAltScript= False
if os.environ.has_key("F_MAP_SCRIPT"):
    altScript= os.environ["F_MAP_SCRIPT"]
    if os.path.isfile(altScript):
        useAltScript= True

if os.environ.has_key("F_MAP_PATH"):
    loadEnv(os.environ["F_MAP_PATH"]);

# Handle flags which effect mappings
for a,b in opts:
    if a=="-m":
	# Overlay any local maps
	loadMap(b)

valDict= { };
if mapDict.has_key("default"):
    addDict= mapDict["default"];
    for k in addDict.keys():
	valDict[k]= addDict[k]

# Handle flags which add mappings to the default dictionary
for a,b in opts:
    if a=="-d":
	# Add this key to the default dictionary
	defineTerm(valDict,b)

mainkey= pargs[0]
if useAltScript:
    execfile(altScript)
    print altMapMethod( valDict, target, mainkey )
else:    
    # Add terms from sub-dictionary associated with mainkey to value dictionary
    if mapDict.has_key(mainkey):
        expandDictionary(valDict,mainkey)
    else:
        raise RuntimeError,"Main key %s has no associated values!"%mainkey;

    # And output the result string
    print recursiveMap(valDict,target,0);
