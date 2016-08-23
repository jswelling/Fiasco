/************************************************************
 *      subset_tool.c                                     *
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

static char rcsid[] = "$Id: subset_tool.c,v 1.2 2005/08/09 23:11:06 welling Exp $"; 

#define BLOCKSIZE (1024*1024)

typedef struct subset_data_struct {
  int dim;
  long shift;
  long extent;
  long upstream_extent;
  long fast_blksize;
  long long slow_blksize;
} SubsetData;

static void recalcBounds( DataSource* self, long* size, long long* offset )
{
  SubsetData* data= (SubsetData*)(self->owner->hook);

  long long n_fast_blks= *offset / data->fast_blksize;
  long fast_blk_offset= *offset - n_fast_blks*data->fast_blksize;
  long long n_full_extents = n_fast_blks / data->extent;
  long extent_offset= n_fast_blks - n_full_extents*data->extent;

  long upstream_extent_offset= extent_offset+data->shift;

  long long upstream_offset= 
    (((n_full_extents*data->upstream_extent)+upstream_extent_offset)
     * data->fast_blksize) + fast_blk_offset;

  /* The maximum size we can return is the rest of the extent (times
   * fast_blksize), because that is the location of the next break.
   */
  long max_size= 
    ((data->extent-extent_offset)*data->fast_blksize)
    - fast_blk_offset;

  *offset= upstream_offset;
  if (*size > max_size) *size= max_size;
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
    fprintf(stderr,"SubsetTool: Request %ld from %lld\n",size, offset);
  recalcBounds(self, &size, &offset);
  if (self->owner->debug) 
    fprintf(stderr,"SubsetTool: Request changed to %ld from %lld\n",
				  size, offset);
  return self->owner->sinkArray[0]->
    source->pGetFloat64Chunk(self->owner->sinkArray[0]->source,
			     size, offset, buf);
}

static DataSource* createSubsetSource( Tool* owner )
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
  SubsetData* data= (SubsetData*)self->hook;
  const char* dimstr= getSourceDims(upstreamSource);

  /* We require that the dimension be in the input stream, and that
   * the input extent be long enough.
   */
  if (!strchr(dimstr,data->dim)) return 0;

  if (data->upstream_extent < data->shift + data->extent) return 0;

  return 1;
}

static int init(Tool* self)
{
  DataSource* mySource= self->sourceArray[0];
  DataSource* upstreamSource= self->sinkArray[0]->source;
  SubsetData* data= (SubsetData*)self->hook;

  if (!baseToolInit(self)) return 0;
  kvCopyUniqueExceptHashes( mySource->attr, upstreamSource->attr );
  mySource->pSetName(mySource, upstreamSource->pGetName(upstreamSource));

  data->upstream_extent= getSourceDimExtent(mySource, (char)data->dim);

  /* Here we edit the attributes of the output to comply with the 
   * subset rules.
   */
  if (structureCheck(self,upstreamSource)) {
    setSourceDimExtent(mySource, (char)data->dim, data->extent);
    calcSourceBlockSizes(upstreamSource, getSourceDims(upstreamSource), 
			 (char)data->dim,
			 &(data->fast_blksize),&(data->slow_blksize));
    if (self->debug) {
      fprintf(stderr,"SubsetTool: dim %c, shift %ld, extent %ld\n",
	      data->dim, data->shift, data->extent);
      fprintf(stderr,
	      "SubsetTool: upstream_extent %ld, fast_blksize %ld, slow_blksize %lld\n",
	      data->upstream_extent, data->fast_blksize, data->slow_blksize);
    }
  }
  else return 0;
  return 1;
}

Tool* createSubsetTool(Arena* arena, const char* dim, const int extent,
		       const int shift ) {
  Tool* result= createBaseTool(arena);
  SubsetData* data= NULL;
  if (!(data=(SubsetData*)malloc(sizeof(SubsetData))))
    Abort("Unable to allocate %d bytes!\n",sizeof(SubsetData));
  if (extent<=0) Abort("createSubsetTool: invalid extent %d!\n",extent);
  if (shift <0) Abort("createSubsetTool: invalid shift %d!\n",shift);
  result->hook= data;
  data->dim= *dim;
  data->shift= shift;
  data->extent= extent;
  data->fast_blksize= 0;
  data->slow_blksize= 0;
  result->pInit= init;
  result->pAddSink( result, createBaseSink(result) );
  result->pAddSource( result, createSubsetSource(result) );
  result->typeName= "subset";
  return result;
}


