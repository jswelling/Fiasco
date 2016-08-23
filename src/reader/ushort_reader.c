/************************************************************
 *                                                          *
 *  ushort_reader.c                                         *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1995 Department of Statistics,         *
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
 *  Original programming by Mark Fitzgerald  5-96           *
 *  Modified to exclusively use libmri calls for output,    *
 *       Greg Hood (PSC), 9-98                              *
 *  Modified to read header files, and to use LX2           *
 *       resampling stuff, Joel Welling (PSC/Stats), 5-1999 *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "mri.h"
#if (SGI64 || SGI5 || SGIMP)
#include <bstring.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

#include "bio.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"
#include "smartreader.h"

static char rcsid[] = "$Id: ushort_reader.c,v 1.2 2007/06/18 20:38:45 welling Exp $";

static void ushortProcessHeader( FileHandler* self, KVHash* info, 
			       SList* chunkStack )
{
  /* Definitions */
  if (!kvLookup(info,"definitions"))
    kvDefHash(info, "definitions", kvFactory(KV_DEFAULT_SIZE));

  kvDefInt(info, "handler_datatype_out", SRDR_INT32);
}

static void ushortRead( FileHandler* self, KVHash* info,
			long long offset, long n,
			SRDR_Datatype datatype_out,
			void* obuf )
{
  int* intBuf= (int*)obuf;
  short* shortBuf= (short*)obuf;
  long i;

  if (datatype_out != SRDR_INT32)
    Abort("%s: ushortRead was asked to output type %s rather than %s!\n",
          progname,
          srdrTypeName[datatype_out],
          srdrTypeName[SRDR_INT32]);
 
  baseReopen(self);
 
  if (debug) fprintf(stderr,"ushortRead: reading %d ushorts at %lld\n",
                     n, offset);
  if (bigfile_fseek(self->file, offset, SEEK_SET))
    Abort("%s: unable to seek to offset %lld in %s: %s\n",
          progname, offset, self->fileName, strerror(errno));

  /* We will read the shorts into the first half of the output array,
   * which is sized for integers.  We'll then do the 1's complement
   * math, working from the back of out output buffer backwards.
   */
  FRdInt16Array(self->file, (short*)obuf, n);
  if (bio_error)
    Abort("%s: read failed on %s trying to read %d bytes starting at %lld\n",
          progname, self->fileName, n*srdrTypeSize[datatype_out], offset);
  for (i=n-1; i>=0; i--) {
    if (shortBuf[i]>=0)
      intBuf[i]= shortBuf[i];
    else
      intBuf[i]= 65536 + shortBuf[i];
  }
}

FileHandler* ushortFactory(char* fname, KVHash* info)
{
  FileHandler* result= baseFactory(fname);
  result->typeName= strdup( "raw unsigned shorts" );
  result->read= ushortRead;
  result->processHeader= ushortProcessHeader;
  return result;
}

int ushortTester(const char* filename)
{
  return 0; /* This reader type can't parse any file header */
}

