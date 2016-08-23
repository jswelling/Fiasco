/************************************************************
 *                                                          *
 *  ram_reader.c                                         *
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

#include "bio.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"
#include "smartreader.h"

static char rcsid[] = "$Id: ram_reader.c,v 1.8 2004/07/27 20:06:16 welling Exp $";

/* A bogus FileHandler that reads data from memory */
typedef struct ram_data_struct {
  void* buf;
  long length; /* in things of this type */
  SRDR_Datatype type;
} RamData;

static void ramDestroySelf( FileHandler* self )
{
  RamData* data= (RamData*)(self->hook);

  free(data->buf);
  baseDestroySelf(self);
}

static void ramProcessHeader( FileHandler* self, KVHash* info, 
			      SList* chunkStack )
{
  RamData* data= (RamData*)(self->hook);
  kvDefInt(info,"datatype_in",data->type);
  kvDefInt(info,"handler_datatype_out",data->type);
}

static void ramRead( FileHandler* self, KVHash* info,
		     long long offset, long n,
		     SRDR_Datatype datatype, void* obuf )
{
  RamData* data= (RamData*)(self->hook);

  if (datatype != data->type) 
    Abort("%s: ramRead was asked for type %s but contains type %s!\n",
	  progname, srdrTypeName[datatype], srdrTypeName[data->type]);

  if (offset<0)
    Abort("%s: ramRead was called with a negative offset!\n",progname);
  if (offset > (data->length - n)*srdrTypeSize[data->type])
    Abort("%s: read of %d from %lld would read past the end of ramRead buffer!\n",
	  progname,n,offset);

  memcpy(obuf, (char*)data->buf + offset, n*srdrTypeSize[data->type]);
}

FileHandler* ramDataHandlerFactory(void* buf, long length, SRDR_Datatype type)
{
  FileHandler* result= baseFactory("NotARealFile");
  RamData* data;

  result->typeName= strdup("RamDataHandler");

  result->destroySelf= ramDestroySelf;
  result->read= ramRead;
  result->processHeader= ramProcessHeader;
  result->totalLengthBytes= length*srdrTypeSize[type];

  if (!(data= (RamData*)malloc(sizeof(RamData))))
    Abort("%s: unable to allocate %d bytes!\n",sizeof(RamData));
  
  data->buf= buf;
  data->length= length;
  data->type= type;
  result->hook= data;

  return result;
}


