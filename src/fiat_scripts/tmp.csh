#!/bin/csh -exf
#
# This routine provides a way to find a machine type when PVM is
# not installed.  For the most part it is just a copy of "pvmgetarch"
# from PVM 3.3.4, translated from sh to csh.
#
# $Id: fiasco_getarch.csh,v 1.1 2003/07/11 22:52:15 welling Exp $

# Let the installed PVM override if available.
if ( ${?PVM_ARCH} ) then
  echo $PVM_ARCH
  exit 0
else if ( ${?PVM_ROOT} ) then
  echo `${PVM_ROOT}/lib/pvmgetarch`
  exit 0
endif

# Translated PVM code follows.  See the top-level Fiasco directory for
# copyright information.
#
# pvmgetarch.sh
#
# Generate PVM architecture string.
#
# This is a heuristic thing that may need to be tuned from time
# to time.  I don't know of a real solution to determining the
# machine type.
#
# Notes:
#   1. Local people mess with things.
#   2. It's good to try a few things for robustness.
#   3. Don't use test -x
#
# 08 Apr 1993  Robert Manchek  manchek@CS.UTK.EDU.
# 24 Aug 1994  last revision
# 28 Jul 1995  release 3.3.8
#

#
# begin section that may need to be tuned.
#
set ARCH = UNKNOWN

#
# determine the machine type from scratch
#
if ( -f /bin/uname || -f /usr/bin/uname ) then
	if ( -f /bin/uname ) then
		set os = "`/bin/uname -s`"
		set ht = "`/bin/uname -m`"
	else
		set os = "`/usr/bin/uname -s`"
		set ht = "`/usr/bin/uname -m`"
	endif

	switch ( "${os}","${ht}" )
	case "SunOS,sun3*" :
          set ARCH = SUN3 
          breaksw
	case "SunOS,sun4*" :
          set ARCH = SUN4 
          breaksw
	case "SunOS,i86pc" :
          set ARCH = X86SOL2 
          breaksw
	case "ULTRIX,RISC" :
          set ARCH = PMAX 
          breaksw
	case "ULTRIX,VAX" :
          set ARCH = UVAX 
          breaksw
	case "AIX*,*" :
          set ARCH = RS6K 
          breaksw
	case "*HP*,9000/[2345]*" :
          set ARCH = HP300 
          breaksw
	case "*HP*,9000/[78]*" :
          set ARCH = HPPA 
          breaksw
	case "IRIX,*" :
          set ARCH = SGI 
          breaksw
	case "*,alpha" :
          set ARCH = ALPHA 
          breaksw
	case "CRSOS,smp" :
          set ARCH = CRAYSMP 
          breaksw
	case "*,paragon" :
          set ARCH = PGON 
          breaksw
	case "dgux,AViiON" :
          set ARCH = DGAV 
          breaksw
	case "*,88k" :
          set ARCH = E88K 
          breaksw
	case "*,mips" :
          set ARCH = MIPS 
          breaksw
	case "*,CRAY-2" :
          set ARCH = CRAY2 
          breaksw
	case "Linux,i[3456]86" :
          set ARCH = LINUX 
          breaksw
	case "BSD/OS,i[3456]86" :
          set ARCH = BSD386 
          breaksw
	case "FreeBSD,i386" :
          set ARCH = FREEBSD 
          breaksw
	case "SUPER-UX,SX-3" :
          set ARCH = SX3 
          breaksw
	case "uts,*" :
          set ARCH = UTS2 
          breaksw
	case "realix,M88*" :
          set ARCH = M88K 
          breaksw
	case "DomainOS,DN*" :
          set ARCH = APOLLO 
          breaksw
	case "Darwin,*" :
	  set ARCH = DARWIN
	  breaksw
	endsw
endif

if ( "$ARCH"  ==  UNKNOWN ) then
	if ( -f /bin/arch ) then
		switch ("`/bin/arch`")
		case ksr1 :
                  set ARCH = KSR1
                  breaksw
		case sun2 :
                  set ARCH = SUN2
                  breaksw
		case sun3 :
                  set ARCH = SUN3
                  breaksw
		case sun4 :
                  set ARCH = SUN4
                  breaksw
		endsw
	endif
endif

if ( "$ARCH"  ==  UNKNOWN ) then

	if ( -f /usr/etc/RELDEF ) set ARCH = ATT

	if ( -f /ultrixboot ) then
		if ( -f /pcs750.bin ) then
			set ARCH = UVAX
		else
			set ARCH = PMAX
		endif
	else
		if ( -f /pcs750.bin ) set ARCH = VAX
	endif

	if ( -d /usr/alliant ) set ARCH = AFX8
	if ( -f /usr/bin/cluster ) set ARCH = BFLY
	if ( -d /usr/convex ) set ARCH = CNVX
	if ( -f /unicos ) set ARCH = CRAY
	if ( -f /hp-ux ) set ARCH = HP300
	if ( -f /usr/bin/getcube ) set ARCH = I860
	if ( -f /usr/bin/asm56000 ) set ARCH = NEXT
	if ( -f /etc/vg ) set ARCH = RS6K
	if ( -d /usr/include/caif ) set ARCH = RT
	if ( -f /bin/4d ) set ARCH = SGI
	if ( -f /dynix ) set ARCH = SYMM
	if ( -f /bin/titan ) set ARCH = TITN

	if ( -f /netbsd ) then
		switch ( "`/usr/bin/machine`" )
		case i386:   
                  set ARCH = NETBSDI386
                  breaksw
		case amiga:  
                  set ARCH = NETBSDAMIGA
                  breaksw
		case hp300:  
                  set ARCH = NETBSDHP300
                  breaksw
		case mac68k: 
                  set ARCH = NETBSDMAC68K
                  breaksw
		case pmax:   
                  set ARCH = NETBSDPMAX
                  breaksw
		case sparc:  
                  set ARCH = NETBSDSPARC
                  breaksw
		case sun3:   
                  set ARCH = NETBSDSUN3
                  breaksw
		endsw
	else if ( -f /usr/bin/machine ) then
		switch ( "`/usr/bin/machine`" )
		case i386: 
                  set ARCH = BSD386
                  breaksw
		endsw
	endif
	if ( -f /usr/bin/uxpm  ) then
          if ( { /usr/bin/uxpm } ) then
	    set ARCH = UXPM
          endif
	endif
endif

if ( "$ARCH"  ==  UNKNOWN ) then
	if ( -f /bin/uname || -f /usr/bin/uname ) then
		if ( -f /bin/uname ) then
			set os = "`/bin/uname -s`"
			set ht = "`/bin/uname -m`"
		else
			set os = "`/usr/bin/uname -s`"
			set ht = "`/usr/bin/uname -m`"
		endif

		switch ( "$os,$ht" )
		case *,i[345]86 :
                  set ARCH = SCO
                  breaksw
		endsw
	endif
endif

#
# update the machine type to derive subclasses
#
if ( "$ARCH"  ==  SUN4 ) then
	set rel = "`/bin/uname -r`"
	switch ( "$rel" )
	case 5.* :   
          set ARCH = SUN4SOL2
          breaksw
	endsw
endif
if ( "$ARCH"  ==  SUN4SOL2 ) then
	set nproc = "`/bin/mpstat | wc -l`"
	if ( $nproc > 2 ) set ARCH = SUNMP
endif
if ( "$ARCH"  ==  ALPHA ) then
	set rel = "`/usr/bin/uname -r`"
	switch ( "$rel" )
	case *3.*:
	  set nproc = "`/usr/sbin/sizer -p`"
	  if ( $nproc > 1 ) set ARCH = ALPHAMP
          breaksw
	endsw
endif
if ( "$ARCH"  ==  SGI ) then
	set rel = "`/bin/uname -r`"
	switch ( "$rel" )
	case 5.* :   
          set ARCH = SGI5
          breaksw
	case 6.* :
          set ARCH = SGI64
          breaksw
	endsw
endif
if ( "$ARCH" ==  SGI64 ) then
	set nproc = "`/usr/sbin/mpadmin -n | wc -w`"
	if ( $nproc > 1 && ${?SGIMP} ) then
          if ( $SGIMP  ==  ON ) set ARCH = SGIMP64
        endif
endif
if ( "$ARCH" ==  SGI5 ) then
	set nproc = "`/usr/sbin/mpadmin -n | wc -w`"
	if ( $nproc > 1 && "$SGIMP"  ==  ON ) set ARCH = SGIMP
endif
if ( "$ARCH"  ==  SUN4 && -f /dev/cm ) set ARCH = CM2
if ( "$ARCH"  == SUN4 && -f /dev/cmni ) set ARCH = CM5
if ( "$ARCH"  ==  CNVX ) then
	if ( { /usr/convex/getsysinfo -f native_default } ) then
		set ARCH = CNVXN
	endif
endif
if ( "$ARCH"  ==  PMAX && -d /usr/maspar ) set ARCH = MASPAR
if ( "$ARCH"  ==  RS6K ) then 
	set nproc = "`/usr/sbin/lsdev -C -c processor | wc -l`"
	if ( $nproc > 1 ) set ARCH = RS6KMP
endif
if ( "$ARCH"  ==  HPPA && -f /bin/sysinfo ) set ARCH = CSPP
if ( "$ARCH"  ==  HPPA ) then
	set nproc = "`/usr/bin/vmstat -n | wc -l`"
	if ( $nproc > 8 ) set ARCH = HPPAMP
endif
if ( "$ARCH"  ==  HPPA ) then
  if ( `/bin/uname -m` ==  9000/780 ) set ARCH = HPPA20
endif

#
# ugh, done.
#

echo $ARCH
exit

