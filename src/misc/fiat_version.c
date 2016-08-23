/************************************************************
 *                                                          *
 *  pulsecalc.c                                             *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 2003 Department of Statistics,         *
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
 *  Original programming by Joel Welling January 2003       *
 ************************************************************/
/* This supplies some simple tests of libpulse. */

#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include "fmri.h"
#include "errors.h"
#include "misc.h"

static char rcsid[] = "$Id: fiat_version.c,v 1.2 2007/03/21 23:54:26 welling Exp $";

static void describeSelf(const char* progname)
{
  fprintf(stderr,"%s: usage:\n",progname);
  fprintf(stderr,"   %s [-major | -minor | -revision | -help]\n",progname);
}

int main( int argc, char* argv[] )
{
  if (argc != 1 && argc != 2) {
    describeSelf(argv[0]);
    exit(-1);
  }
  if (argc==1) printf("%s\n",FIAT_VERSION_STRING);
  else if (!strcasecmp(argv[1],"-major"))
    printf("%d\n",FIAT_MAJOR_VERSION);
  else if (!strcasecmp(argv[1],"-minor"))
    printf("%d\n",FIAT_MINOR_VERSION);
  else if (!strcasecmp(argv[1],"-revision"))
    printf("%d\n",FIAT_REVISION);
  else if (!strcasecmp(argv[1],"-h") || !strcasecmp(argv[1],"-help"))
    describeSelf(argv[0]);
  else {
    describeSelf(argv[0]);
    exit(-1);
  }
  return 0;
}

