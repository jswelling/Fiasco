#! /bin/csh -ef
if ( ! -e $1 ) exit 1
if ( ! -e $2 ) exit 0
set foo = `ls -t $1 $2`
if ( ${foo[1]} == $1 ) then
  exit 0
else
  exit 1
endif
