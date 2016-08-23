/************************************************************
 *                                                          *
 *  convert_reader.c                                         *
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
#include <strings.h>
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

static char rcsid[] = "$Id: convert_reader.c,v 1.5 2007/06/18 20:38:44 welling Exp $";

/* A bogus FileHandler that converts data types */

typedef struct convert_data_struct {
  FileHandler* child;
  SRDR_Datatype type_in;
  SRDR_Datatype type_out;
  void* tbuf;
  int tbufSize; /* in "things", not bytes */
} ConvertData;

static void destroySelf( FileHandler* self )
{
  ConvertData* data= (ConvertData*)(self->hook);

  FH_DESTROYSELF(data->child);
  if (data->tbufSize>0) free(data->tbuf);
  baseDestroySelf(self);
}

static void processHeader( FileHandler* self, KVHash* info, SList* chunkStack )
{
  ConvertData* data= (ConvertData*)(self->hook);

  FH_PROCESSHEADER( data->child, info, chunkStack );

  /* I have nothing more to contribute. */
}

static void uint16_to_int32_convert(void* temp_out, void* temp_in, 
				    long length, SRDR_Datatype datatype_in)
{
  short* short_in= (short*)temp_in;
  int* int_out= (int*)temp_out;
  int i;

  /* We must do 1's complement arithmetic on the incoming unsigned shorts. */
  for (i=0; i<length; i++) {
    if (short_in[i]>=0)
      int_out[i]= short_in[i];
    else
      int_out[i]= 65536 + short_in[i];
  }

}

static void float_convert(void* temp_out, void* temp_in, long length, 
			  SRDR_Datatype datatype_in)
{
  int i;
  float* fout= (float*)temp_out;

  switch (datatype_in) {
    case SRDR_UINT8:
      for (i=0; i<length; i++) fout[i]= (float)(*((unsigned char*)temp_in+i));
      break;
    case SRDR_INT16:
      for (i=0; i<length; i++) fout[i]= (float)(*((short*)temp_in+i));
      break;
    case SRDR_INT32:
      for (i=0; i<length; i++) fout[i]= (float)(*((int*)temp_in+i));
      break;
    case SRDR_FLOAT32:
      bcopy(temp_in, temp_out, length*srdrTypeSize[datatype_in]);
      break;
    case SRDR_FLOAT64:
      for (i=0; i<length; i++) fout[i]= (float)(*((double*)temp_in+i));
      break;
    default:
      Abort("%s: Internal error: unknown datatype %d\n",progname,
	    datatype_in);
  }
}

static void convertRead( FileHandler* self, KVHash* info,
			 long long offset, long n,
			 SRDR_Datatype datatype,
			 void* obuf )
{
  ConvertData* data= (ConvertData*)(self->hook);

  if (datatype != data->type_out) 
    Abort("%s: convertRead was asked for type %s but produces type %s!\n",
	  progname, srdrTypeName[datatype], srdrTypeName[data->type_out]);
  
  if (data->tbufSize < n) {
    if (data->tbuf) free(data->tbuf);
    if (!(data->tbuf=malloc(n*srdrTypeSize[data->type_out])))
      Abort("%s: unable to allocate %d bytes!\n", progname, data->tbufSize);
    data->tbufSize= n;
  }
  
  if (data->type_out == data->type_in) {
    FH_READ(data->child, info, offset, n, data->type_in, obuf);
  }
  else {
    FH_READ(data->child, info, offset, n, data->type_in, data->tbuf);
    if (data->type_out==SRDR_FLOAT32) {
      float_convert(obuf, data->tbuf, n, data->type_in);
    }
    else if (data->type_in==SRDR_UINT16 && data->type_out==SRDR_INT32) {
      uint16_to_int32_convert(obuf, data->tbuf, n, data->type_in);
    }
    else
      Abort("%s: internal error: conversion from type %s to type %s on input is not supported!\n",
	    progname, srdrTypeName[data->type_in],
	    srdrTypeName[data->type_out]);
  }
}

FileHandler* convertFactory(FileHandler* child, 
			    SRDR_Datatype type_in, SRDR_Datatype type_out)
{
  FileHandler* result= baseFactory("NotARealFile");
  ConvertData* data;
  char* typeName;
  int typeNameLength;
  
  if (type_out != type_in 
      && type_out != SRDR_FLOAT32
      && (!(type_in==SRDR_UINT16 && type_out==SRDR_INT32)))
    Abort("%s: conversion to type %s on input is not supported!\n",
	  progname, srdrTypeName[type_out]);
  
  typeNameLength= 
    strlen("Converter[,]()") 
    + strlen(srdrTypeName[type_in]) + strlen(srdrTypeName[type_out])
    + strlen(child->typeName) + 8;
  if (!(typeName=(char*)malloc(typeNameLength)))
    Abort("%s: unable to allocate %d bytes!\n",progname,typeNameLength);
  sprintf(typeName,"Converter[%s,%s](%s)",
	  srdrTypeName[type_out], srdrTypeName[type_in], child->typeName);
  
  result->typeName= typeName;
  result->totalLengthBytes= child->totalLengthBytes;
  
  result->destroySelf= destroySelf;
  result->read= convertRead;
  result->processHeader= processHeader;
  
  if (!(data= (ConvertData*)malloc(sizeof(ConvertData))))
    Abort("%s: unable to allocate %d bytes!\n",sizeof(ConvertData));
  
  data->child= child;
  data->type_in= type_in;
  data->type_out= type_out;
  data->tbuf= NULL;
  data->tbufSize= 0; /* in "things", not bytes */

  result->hook= data;

  return result;
}


