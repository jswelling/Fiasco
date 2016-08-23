/************************************************************
 *                                                          *
 *  base_reader.c                                         *
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
/*
 * These methods implement a sort of "base class" which supplies
 * functionality common to other reader types.
 */

#if defined(SUN4SOL2) || defined(LINUX) || defined(HPPA) || defined(HPPA20) || defined(ALPHA) || defined(DARWIN)
#define NO_FSEEK64
#if !defined(HPPA)
#define HAVE_FSEEKO
#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE
#endif
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#if (SGI64 || SGI5 || SGIMP)
#include <bstring.h>
#endif

#include "mri.h"
#include "bio.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"
#include "smartreader.h"

static char rcsid[] = "$Id: base_reader.c,v 1.20 2007/06/12 15:51:43 welling Exp $";

int
bigfile_fseek (FILE *f, long long offset, int whence)
{
#ifdef NO_FSEEK64
#ifdef HAVE_FSEEKO
  return(fseeko(f, (off_t) offset, whence));
#else
  if (offset < 0 || offset >= LONG_MAX)
    {
      mri_report_error(NULL, "mri_fseek: offset out-of-range\n");
      return(-1);
    }
  return(fseek(f, (long) offset, whence));
#endif
#else
  return(fseek64(f,offset,whence));
#endif
}

void baseDestroySelf( FileHandler* self )
{
  FH_CLOSE(self);
  if (self->fileName) free(self->fileName);
  if (self->hook) free(self->hook);
  if (self->typeName) free((char*)self->typeName);
  free(self);
}

void baseProcessHeader( FileHandler* self, KVHash* info, SList* chunkStack )
{
  /* Definitions */
  if (!kvLookup(info,"definitions")) 
    kvDefHash(info, "definitions", kvFactory(KV_DEFAULT_SIZE));

  /* I have nothing to contribute. */
}

static int read_binary( FILE* fpi, SRDR_Datatype type_in, int n, void* buf, 
			int big_endian_input)
{
  switch (type_in) {
  case SRDR_UINT8:
    FRdUInt8Array(fpi, (unsigned char*)buf, n);
    break;
  case SRDR_INT16:
    FRdInt16Array(fpi, (short*)buf, n);
    break;
  case SRDR_INT32:
    FRdInt32Array(fpi, (int*)buf, n);
    break;
  case SRDR_FLOAT32:
    FRdFloat32Array(fpi, (float*)buf, n);
    break;
  case SRDR_FLOAT64:
    FRdFloat64Array(fpi, (double*)buf, n);
    break;
  case SRDR_INT64:
    FRdInt64Array(fpi, (long long*)buf, n);
    break;
  default:
    Abort("%s: Internal error: unknown datatype %d\n",progname,
	  type_in);
  }
  return( bio_error==0 );
}

void baseRead( FileHandler* self, KVHash* info,
	       long long offset, long n,
	       SRDR_Datatype datatype_out,
	       void* obuf )
{
  if (datatype_out != kvGetInt(info,"datatype_in"))
    Abort("%s: baseRead was asked to do a translation (%s to %s)!\n",
	  progname,
	  srdrTypeName[kvGetInt(info,"datatype_in")],
	  srdrTypeName[datatype_out]);

  baseReopen(self);

  if (debug) fprintf(stderr,"baseRead: reading %d of type %s at %lld\n",
		     n, srdrTypeName[datatype_out], offset);
  if (bigfile_fseek(self->file, offset, SEEK_SET))
    Abort("%s: unable to seek to offset %lld in %s: %s\n",
	  progname, offset, self->fileName, strerror(errno));

  if (!read_binary(self->file, datatype_out, n, obuf, 
		   kvGetBoolean(info,"big_endian_input")))
    Abort("%s: read failed on %s trying to read %d bytes starting at %lld\n",
	  progname, self->fileName, n*srdrTypeSize[datatype_out], offset);
}

void baseClose( FileHandler* self )
{
  if (self->file) {
    if (fclose(self->file)) perror("Error closing file");
    self->file= NULL;
  }
}

void baseReopen( FileHandler* self )
{
  if (!(self->file)) {
    if (!(self->file= fopen(self->fileName,"r")))
      Abort("%s: unable to open file <%s> for reading!\n",
	    progname,self->fileName);
  }
}

int baseCompare( FileHandler* f1, FileHandler* f2 )
{
  /* Just alphabetically compare on filenames */
  return strcoll(f1->fileName, f2->fileName);
}

static long long getFileSize( const char* fname )
{
  struct stat s;
  if (strcmp(fname,"NotARealFile")) {
    if (stat(fname,&s))
      Abort("%s: base_reader: stat failed on %s: %s!\n",
	    progname, fname, strerror(errno));
    return (long long)s.st_size;
  }
  else return 0;
}

FileHandler* baseFactory(char* fname)
{
  FileHandler* result;
  if (!(result= (FileHandler*)malloc(sizeof(FileHandler)))) 
    Abort("Unable to allocate %d bytes!\n",sizeof(FileHandler));

  result->fileName= strdup(fname);
  result->totalLengthBytes= getFileSize(fname);
  result->file= NULL;
  result->processHeader= baseProcessHeader;
  result->read= baseRead;
  result->destroySelf= baseDestroySelf;
  result->close= baseClose;
  result->reopen= baseReopen;
  result->compareThisType= baseCompare;
  result->typeName= NULL;
  result->hook= NULL;
  return result;
}

