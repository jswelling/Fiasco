/************************************************************
 *                                                          *
 *  byte_extract.c                                     *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1998 Department of Statistics,         *
 *                        Carnegie Mellon University        *
 *                                                          *
 *  This program is distributed in the hope that it will    *
 *  be useful, but WITHOUT ANY WARRANTY; without even the   *
 *  implied warranty of MERCHANTABILITY or FITNESS FOR A    *
 *  PARTICULAR PURPOSE.  Neither Carnegie Mellon University *
 *  nor any of the authors assume any liability for         *
 *  damages, incidental or otherwise, caused by the         *
 *  installation or use of this software.                   *
 *                                                          *
 *  CLINICAL APPLICATIONS ARE NOT RECOMMENDED, AND THIS     *
 *  SOFTWARE HAS NOT BEEN EVALUATED BY THE UNITED STATES    *
 *  FDA FOR ANY CLINICAL USE.                               *
 *                                                          *
 *                                                          *
 *  Original programming by Joel Welling 5-98               *
 ************************************************************/
/*************************************************************

  DESCRIPTION OF byte_extract.c
  This program takes a file of raw bytes, assumed to represent a
  series of subgroups of identical length.  From each subgroup it
  emits a specified contiguous subset of the bytes.  Input is
  from stdin; output is to stdout.

  byte_extract -cycle cycle_len -from start -n nbytes

    -length cycle_len specifies the number of bytes in each cycle.
    -from start specifies the offset to start from within each cycle
          (counting from 0)
    -n nbytes specifies the number of bytes from the given cycle to
          be output.

**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "misc.h"
#include "stdcrg.h"

static char rcsid[] = "$Id: byte_extract.c,v 1.5 2003/02/07 21:28:33 welling Exp $";

static char* progname;

int main( int argc, char *argv[] )
{
  int cycle;
  int offset;
  int nbytes;
  char* inbuf= NULL;
  int cycles_copied= 0;

  progname= argv[0];

  /* Check to see if help was requested */
  if( ( argc > 1 ) && !strcmp( argv[1], "-help" ) )
    {
      if( argc == 2 )
        Help( "selecttopic" );
      else
        Help( (char*)(argv[2]) );
    }

  /*** Parse command line ***/

  cl_scan( argc, (char**)argv );

  if (!cl_get( "cycle", "%option %d", &cycle )) {
    fprintf(stderr,"%s: -cycle option required.\n",argv[0]);
    Help( "usage" );
    exit(-1);
  }
  
  if (!cl_get( "from", "%option %d", &offset )) {
    fprintf(stderr,"%s: -from option required.\n",argv[0]);
    Help( "usage" );
    exit(-1);
  }
  
  if (!cl_get( "n", "%option %d", &nbytes )) {
    fprintf(stderr,"%s: -n option required.\n",argv[0]);
    Help( "usage" );
    exit(-1);
  }
  
  if (cl_cleanup_check()) {
    int i;
    fprintf(stderr,"%s: invalid argument in command line:\n    ",argv[0]);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(-1);
  }
  /*** End command-line parsing ***/

  if (nbytes+offset > cycle) {
    fprintf(stderr,"%s: cycle too short for specified output bytes!\n", 
	    progname);
    exit(-1);
  }
  
  if ( !(inbuf= (char*)malloc(cycle)) ) {
    fprintf(stderr,"%s: unable to allocate %d byte input buffer!\n",
	    progname,cycle);
    exit(-1);
  }
    
  cycles_copied= 0;
  while (!feof(stdin)) {
    if (fread(inbuf,1,cycle,stdin) != cycle) {
      if (feof(stdin)) break;
      else {
	fprintf(stderr,
		"%s: error reading from stdin; input not an integral number of cycles?\n",
		progname);
	exit(-1);
      }
    }
    if (fwrite(inbuf+offset, 1, nbytes, stdout) != nbytes) {
      fprintf(stderr,"%s: error writing to stdout!\n", progname);
      exit(-1);
    }
    cycles_copied++;
  }
  
  /* Just on the off chance we got an input file that was too short... */
  if (cycles_copied<1) {
    fprintf(stderr,"%s: input dataset does not contain a full cycle!\n",
	    progname);
    exit(-1);
  }

  return 0;
}
