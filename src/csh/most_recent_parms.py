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
# This script finds the most recent entry in the F_SUMM_INPUT file
# of a given parameter type and returns the filename.
#

import sys
import os
import string
import getopt

# Mode flags
verbose= 0;
tailname= 0;
failOnMissing= 0;

def trim( str ):
    where= len(str)-1;
    while str[where]=='\012':
        str = str[0:where]
        where = where-1
        if where<0:
            break
    string.strip(str)
    return str

def match( str, target ):
    val= string.find( str, target );
    if val < 0:
	return val;
    else:
	if val > 0:
	    if str[val-1] != '.':
		return -1;
	if val + len(target) < len(str):
	    if str[val+len(target)] != '.':
		return -1;
    return val;

def describeSelf():
    print """
Usage: %s [-v] [-t] [-f] typestring

-v specifies verbose mode.
-t returns the tail of the filename only, dropping leading directories.
-f causes an error exit if the file does not exist.
typestring is the parameter type string, for example 'deghost'.

The file specified by the environment variable F_SUMM_INPUT is
scanned, and the last entry for which the first word contains the
given typestring is found.  The second word of that entry is the
filename where the corresponding set of parameters is located;
that filename is returned.  The given typestring must occur in the
first word bounded by spaces or '.' .  If no such entry is found an
empty string is returned.

""" % progname

##############################
#
# Main
#
##############################

progname= os.path.basename(sys.argv[0]);

if len(sys.argv) == 1 :
    describeSelf()
    sys.exit(-1)

try:
    (opts,pargs) = getopt.getopt(sys.argv[1:],'-v-t-f')
except:
    sys.stderr.write("%s: Invalid command line parameter\n" % progname);
    describeSelf();
    sys.exit(-1);

if len(pargs)!=1:
    sys.stderr.write("%s: too few or too many typestrings\n" % progname);
    describeSelf();
    sys.exit(-1);

typestring = pargs[0];

# Handle flags not effecting mappings
for a,b in opts:
    if a=="-v":
	# verbose mode on
	verbose= 1;
    if a=="-t":
	# tailname mode on
	tailname= 1;
    if a=="-f":
	# tailname mode on
	failOnMissing= 1;

# Find the parameter file list
if os.environ.has_key("F_SUMM_INPUT"):
    fname= string.strip(os.environ["F_SUMM_INPUT"]);
    if verbose :
	print "Loading environment %s"%fname;
    f= open(fname);
    lines= map(trim,f.readlines());
    f.close();
else:
    sys.stderr.write( "%s: F_SUMM_INPUT is not defined!\n" % progname);
    sys.exit(-1);

# Scan for the answer
result= "";

for l in lines:
    if len(l)>0:
	[a,b]= string.split(l);
	if verbose : 
	    print "<%s> <%s>; target <%s>"%(a,b,typestring);
	if match(a,typestring) >= 0:
	    result= b;

if failOnMissing:
    if len(result)==0:
	sys.stderr.write( "%s: no parameter data of type %s!\n"%(progname,typestring));
	sys.exit(-1);
    if not os.access(result,os.R_OK):
	sys.stderr.write( "%s: cannot find parameter file %s!\n"%(progname,result));
	sys.exit(-1);    

if tailname:
    loc= string.rfind(result,"/");
    if loc < 0:
	print result;
    else:
	print result[loc+1:];
else:
    print result;
