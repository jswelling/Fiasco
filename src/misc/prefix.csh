#!/bin/csh
#$Id: prefix.csh,v 1.1 2003/03/19 00:13:34 welling Exp $
if ( $1 == "-help" ) then
	echo "$0 is used only by make, in building the Fiasco executables"
	exit -1
endif
set p=$1
shift
foreach i ($*)
    echo -n $p${i} " "
end
echo
