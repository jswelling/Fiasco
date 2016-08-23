/************************************************************
 *      block_map_tool.c                                     *
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
#include <fexceptions.h>
#include <mripipes.h>

/***********************
 * Notes-
 **********************/

static char rcsid[] = "$Id: block_map_tool.c,v 1.3 2006/07/31 23:24:15 welling Exp $"; 

#define BLOCKSIZE (1024*1024)

typedef struct block_map_data_struct {
  int dim;
  int newdim;
  long extent1;  /* output extent of dim */
  long extent2; /* output extent of newdim */
  long upstream_extent; /* input extent of dim */
  long fast_blksize;
  long long slow_blksize;
  BLOCKMAPCBFUNC cbFunc;
  BLOCKMAPINITFUNC initFunc;
  void* cbData;
} BlockMapData;

static void recalcBounds( DataSource* self, long* size, long long* offset )
{
  int errCode=0;
  BlockMapData* data= (BlockMapData*)(self->owner->hook);
  long cb_size= *size;
  long long cb_offset= *offset;

  if (self->owner->debug) fprintf(stderr,"Calling back...");
  if (!(data->cbFunc(&cb_size, &cb_offset, data->cbData)))
    pipeAbort("Callback function has signaled an error!\n");
  if (self->owner->debug) 
    fprintf(stderr,"Made it! Values returned are %ld %lld\n",
	    cb_size,cb_offset);

  *offset= cb_offset;
  if (cb_size<*size) *size=cb_size;
}


static long getUInt8Chunk(DataSource* self,
				long size, long long offset, char* buf)
{
  recalcBounds(self, &size, &offset);
  return self->owner->sinkArray[0]->
    source->pGetUInt8Chunk(self->owner->sinkArray[0]->source,
			   size, offset, buf);
}

static long getInt16Chunk(DataSource* self,
				long size, long long offset, short* buf)
{
  recalcBounds(self, &size, &offset);
  return self->owner->sinkArray[0]->
    source->pGetInt16Chunk(self->owner->sinkArray[0]->source,
			   size, offset, buf);
}

static long getInt32Chunk(DataSource* self,
				long size, long long offset, int* buf)
{
  recalcBounds(self, &size, &offset);
  return self->owner->sinkArray[0]->
    source->pGetInt32Chunk(self->owner->sinkArray[0]->source,
			   size, offset, buf);
}

static long getInt64Chunk(DataSource* self,
				long size, long long offset, long long* buf)
{
  recalcBounds(self, &size, &offset);
  return self->owner->sinkArray[0]->
    source->pGetInt64Chunk(self->owner->sinkArray[0]->source,
			   size, offset, buf);
}

static long getFloat32Chunk(DataSource* self,
				long size, long long offset, float* buf)
{
  recalcBounds(self, &size, &offset);
  return self->owner->sinkArray[0]->
    source->pGetFloat32Chunk(self->owner->sinkArray[0]->source,
			     size, offset, buf);
}

static long getFloat64Chunk(DataSource* self,
				long size, long long offset, double* buf)
{
  if (self->owner->debug) 
    fprintf(stderr,"BlockMapTool: Request %ld from %lld\n",size, offset);
  recalcBounds(self, &size, &offset);
  if (self->owner->debug) 
    fprintf(stderr,"BlockMapTool: Request changed to %ld from %lld\n",
				  size, offset);
  return self->owner->sinkArray[0]->
    source->pGetFloat64Chunk(self->owner->sinkArray[0]->source,
			     size, offset, buf);
}

static DataSource* createBlockMapSource( Tool* owner )
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

static int structureCheck( Tool* self, DataSource* upstreamSource )
{
  BlockMapData* data= (BlockMapData*)self->hook;
  const char* dimstr= getSourceDims(upstreamSource);

  /* We require that the dimension be in the input stream, and that
   * its extent in the input stream is at least as large as the fast
   * extent in the output stream.
   */
  if (!strchr(dimstr,data->dim)) 
    pipeAbort("Input data stream does not include dimension %c!\n",data->dim);
  if (data->extent1>data->upstream_extent)
    pipeAbort("Input data stream dim %c extent is too small!\n",data->dim);

  return 1;
}

static int init(Tool* self)
{
  DataSource* mySource= self->sourceArray[0];
  DataSource* upstreamSource= self->sinkArray[0]->source;
  BlockMapData* data= (BlockMapData*)self->hook;

  if (!baseToolInit(self)) return 0;
  kvCopyUniqueExceptHashes( mySource->attr, upstreamSource->attr );
  mySource->pSetName(mySource, upstreamSource->pGetName(upstreamSource));

  data->upstream_extent= getSourceDimExtent(mySource, (char)data->dim);

  /* Here we edit the attributes of the output to comply with the 
   * block_map rules.
   */
  if (structureCheck(self,upstreamSource)) {
    const char* dimstr= getSourceDims(upstreamSource);
    char buf[256];
    long i;
    if (strlen(dimstr)+3 > sizeof(buf))
      pipeAbort("Input data stream dimension string is too long!\n");
    i=0;
    while (dimstr[i]!=data->dim) { buf[i]= dimstr[i]; i++; }
    buf[i++]= data->dim;
    buf[i++]= data->newdim;
    while (dimstr[i-1]) { buf[i]= dimstr[i-1]; i++; }
    buf[i]= '\0';
    setSourceDims(mySource, buf);
    setSourceDimExtent(mySource, (char)data->dim, data->extent1);
    setSourceDimExtent(mySource, (char)data->newdim, data->extent2);
    calcSourceBlockSizes(upstreamSource, dimstr, (char)data->dim,
			 &(data->fast_blksize),&(data->slow_blksize));
    if (self->verbose) {
      fprintf(stderr,
	      "BlockMapTool: dims %c%c, extents %ld,%ld, upstream_extent %ld\n",
	      data->dim, data->newdim, data->extent1, data->extent2,
	      data->upstream_extent);
      fprintf(stderr,
	      "              fast_blksize %ld, slow_blksize %lld\n",
	      data->fast_blksize, data->slow_blksize);
    }
    if (self->debug) fprintf(stderr,"Calling back for init...");
    if (!(data->initFunc(data->dim, dimstr,
			 data->fast_blksize, data->upstream_extent,
			 data->extent1, data->extent2,
			 data->slow_blksize, data->cbData)))
      pipeAbort("Init callback function has signaled an error!\n");
    if (self->debug) 
      fprintf(stderr,"Made it!\n");
  }
  else return 0;

  return 1;
}

Tool* createBlockMapTool(Arena* arena, const char* dim, const char* newdim,
			 long extent1, long extent2,
			 BLOCKMAPINITFUNC initFunc,
			 BLOCKMAPCBFUNC cbFunc, void* cbData ) {
  Tool* result= createBaseTool(arena);
  BlockMapData* data= NULL;
  if (!(data=(BlockMapData*)malloc(sizeof(BlockMapData))))
    pipeAbort("Unable to allocate %d bytes!\n",sizeof(BlockMapData));
  if (extent1<=0) 
    pipeAbort("createBlockMapTool: invalid extent1 %d!\n",extent1);
  if (extent2<=0) 
    pipeAbort("createBlockMapTool: invalid extent2 %d!\n",extent2);
  result->hook= data;
  data->dim= *dim;
  data->newdim= *newdim;
  data->extent1= extent1;
  data->extent2= extent2;
  data->fast_blksize= 0;
  data->slow_blksize= 0;
  data->initFunc= initFunc;
  data->cbFunc= cbFunc;
  data->cbData= cbData;
  result->pInit= init;
  result->pAddSink( result, createBaseSink(result) );
  result->pAddSource( result, createBlockMapSource(result) );
  result->typeName= "block_map";
  return result;
}


