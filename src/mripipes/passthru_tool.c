/************************************************************
 *      passthru_tool.c                               *
 *                                                          *
 *	Copyright (c) 2004 Pittsburgh Supercomputing Center *
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
 *	History:                                            *
 *		5/04: Written by Joel Welling               *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mri.h>
#include <fmri.h>
#include <stdcrg.h>
#include <kvhash.h>
#include <mripipes.h>

/***********************
 * Notes-
 **********************/

static char rcsid[] = "$Id: passthru_tool.c,v 1.4 2005/08/09 23:11:06 welling Exp $"; 

#define BLOCKSIZE (1024*1024)

static DataSink* createPassthruSink( Tool* owner )
{
  DataSink* result= createBaseSink(owner);
  return result;
}

static long getUInt8Chunk(DataSource* self,
				long size, long long offset, char* buf)
{
  return self->owner->sinkArray[0]->
    source->pGetUInt8Chunk(self->owner->sinkArray[0]->source,
			   size, offset, buf);
}

static long getInt16Chunk(DataSource* self,
				long size, long long offset, short* buf)
{
  return self->owner->sinkArray[0]->
    source->pGetInt16Chunk(self->owner->sinkArray[0]->source,
			   size, offset, buf);
}

static long getInt32Chunk(DataSource* self,
				long size, long long offset, int* buf)
{
  return self->owner->sinkArray[0]->
    source->pGetInt32Chunk(self->owner->sinkArray[0]->source,
			   size, offset, buf);
}

static long getInt64Chunk(DataSource* self,
				long size, long long offset, long long* buf)
{
  return self->owner->sinkArray[0]->
    source->pGetInt64Chunk(self->owner->sinkArray[0]->source,
			   size, offset, buf);
}

static long getFloat32Chunk(DataSource* self,
				long size, long long offset, float* buf)
{
  return self->owner->sinkArray[0]->
    source->pGetFloat32Chunk(self->owner->sinkArray[0]->source,
			     size, offset, buf);
}

static long getFloat64Chunk(DataSource* self,
				long size, long long offset, double* buf)
{
  return self->owner->sinkArray[0]->
    source->pGetFloat64Chunk(self->owner->sinkArray[0]->source,
			     size, offset, buf);
}

static DataSource* createPassthruSource( Tool* owner )
{
  DataSource* result= createBaseSource(owner);
  result->pGetUInt8Chunk= getUInt8Chunk;
  result->pGetInt16Chunk= getInt16Chunk;
  result->pGetInt32Chunk= getInt32Chunk;
  result->pGetInt64Chunk= getInt64Chunk;
  result->pGetFloat32Chunk= getFloat32Chunk;
  result->pGetFloat64Chunk= getFloat64Chunk;
  return result;
}

static int init(Tool* self)
{
  DataSource* mySource= self->sourceArray[0];
  DataSource* upstreamSource= self->sinkArray[0]->source;
  if (!baseToolInit(self)) return 0;
  kvCopyUniqueExceptHashes( mySource->attr, upstreamSource->attr );
  mySource->pSetName(mySource, upstreamSource->pGetName(upstreamSource));
  return 1;
}

Tool* createPassthruTool(Arena* arena) {
  Tool* result= createBaseTool(arena);
  result->pInit= init;
  result->pAddSink( result, createPassthruSink(result) );
  result->pAddSource( result, createPassthruSource(result) );
  result->typeName= "passthru";
  return result;
}


