#! /bin/csh -f
# fiasco_update_from_cvs.csh
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
# *  Original programming by Bill Eddy                       *
# ************************************************************/
#
echo '#' `date` $0
set idstring = '$Id: fiasco_update_from_cvs.csh,v 1.12 2005/03/10 02:05:10 welling Exp $'
echo '#' $idstring

# This script updates the Fiasco directory tree from which it is called 
# from CVS. 

# Are we at the top of a Fiasco tree?
if ( ! -e Makefile.common || ! -d src || ! -d bin ) then
  echo "This script needs to be run from the top of a Fiasco tree."
  exit -1
endif

# Set the architecture from Fiasco's copy of pvmgetarch

if ( ! -x src/csh/fiasco_getarch.csh \
     && ! -x src/fiat_scripts/fiasco_getarch.csh ) then
  set arch = `pvmgetarch`
  if ( $status != 0 ) then
    echo "pvmgetarch is not in the path and I cannot find fiasco_getarch.csh!"
    exit -1
  endif
endif
if ( -x src/fiat_scripts/fiasco_getarch.csh ) then
  set arch = `src/fiat_scripts/fiasco_getarch.csh`
else 
  if ( -x src/csh/fiasco_getarch.csh ) then
    set arch = `src/csh/fiasco_getarch.csh`
  endif
endif

# Are we at the top of a Fiasco tree? Check a bit more deeply
if ( ! -d bin/${arch} || ! -f bin/${arch}/FIASCO ) then
  echo "This script needs to be run from the top of a Fiasco tree."
  exit -1
endif

# Is CVS available?
if ( ! ${?CVSROOT} && ! -e CVS/Root ) then
  echo 'No CVS root available; set CVSROOT!'
  exit -1
endif

if ( ${arch} =~ SGI* ) then
  set logfile = update_`date '+%Y-%m-%d'`.log
else
  set logfile = update_`date -I`.log
endif

echo "Update script ID string: $idstring" > $logfile
echo "Update on " `date` >> $logfile
echo '###########################################################' >> $logfile
echo '#####      cvs update -d ' $* >> $logfile
echo '###########################################################' >> $logfile
echo "Running CVS update..."
cvs update -d $* |& grep -v doc/generated | grep -v cvs.locks | \
    grep -v depend.mk. | grep -v $arch >> $logfile
echo '###########################################################' >> $logfile
echo '#####      ./configure' >> $logfile
echo '###########################################################' >> $logfile
echo "Running configure script..."
./configure >>& $logfile
echo '###########################################################' >> $logfile
echo '#####      make depend' >> $logfile
echo '###########################################################' >> $logfile
echo "Running make depend..."
make depend >>& $logfile
echo '###########################################################' >> $logfile
echo '#####      make' >> $logfile
echo '###########################################################' >> $logfile
echo "Running make..."
make >>& $logfile
echo '###########################################################' >> $logfile
echo '#####      Done' >> $logfile
echo '###########################################################' >> $logfile

echo 'Mailing results...'
Mail -s "Fiasco CVS update on ${HOST}, ${PWD} at `date`" \
    fiasco_admin@stat.cmu.edu < $logfile
echo 'Done'

