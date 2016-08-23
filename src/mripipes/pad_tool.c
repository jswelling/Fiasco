/************************************************************
 *      pad_tool.c                                     *
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

static char rcsid[] = "$Id: pad_tool.c,v 1.3 2005/08/09 23:11:06 welling Exp $"; 

#define BLOCKSIZE (1024*1024)

typedef struct pad_data_struct {
  int dim;
  long shift;
  long extent;
  long upstream_extent;
  long fast_blksize;
  long long slow_blksize;
  double fillValue;
} PadData;

#define TYPEFILL( type, val, len, buf, bufOffset ) \
{ \
  type *here= (type*)buf + bufOffset; \
  type *there= here + len; \
  while (here<there) *here++= (type)val; \
}

static void fill(long long bufOffset, long n, double val, int type, void* buf)
{
  switch (type) {
  case MRI_UNSIGNED_CHAR: TYPEFILL( unsigned char, val, n, buf, bufOffset ); break;
  case MRI_SHORT: TYPEFILL( short, val, n, buf, bufOffset ); break;
  case MRI_INT: TYPEFILL( int, val, n, buf, bufOffset ); break;
  case MRI_LONGLONG: TYPEFILL( long long, val, n, buf, bufOffset ); break;
  case MRI_FLOAT: TYPEFILL( float, val, n, buf, bufOffset ); break;
  case MRI_DOUBLE: TYPEFILL( double, val, n, buf, bufOffset ); break;
  }
}

#undef TYPEFILL

static long getTypeChunk(DataSource* self, int type,
			 long size, long long offset, void* buf)
{
  PadData* data= (PadData*)(self->owner->hook);
  long long n_fast_blks= offset / data->fast_blksize;
  long fast_blk_offset= offset - n_fast_blks*data->fast_blksize;
  long long n_full_extents = n_fast_blks / data->extent;
  long extent_offset= n_fast_blks - n_full_extents*data->extent;
  long n;
  long nGot= 0;
  long long baseOffset= offset;
  DataSource* src= NULL;

  /* The maximum size we can return is the rest of the extent (times
   * fast_blksize), because that is the location of the next break.
   */
  long max_size= 
    ((data->extent-extent_offset)*data->fast_blksize)
    - fast_blk_offset;

  if (self->owner->debug) {
    fprintf(stderr,
	    "PadTool: n_fast_blks= %lld, fast_blk_offset= %d, n_full_extents= %lld\n",
	    n_fast_blks, fast_blk_offset, n_full_extents);
    fprintf(stderr,
	    "PadTool: extent_offset= %d, baseOffset= %lld, size= %d, max_size= %d\n",
	    extent_offset, baseOffset, size, max_size);
  }

  if (size>max_size) size= max_size;

  if (extent_offset < data->shift) {
    /* We are in an initial block of fill */
    n= 0;
    if (fast_blk_offset>0) {
      n += data->fast_blksize - fast_blk_offset;
      fast_blk_offset= 0;
      extent_offset += 1;
    }
    if (extent_offset < data->shift) { 
      n += (data->shift - extent_offset)*data->fast_blksize;
      extent_offset= data->shift;
    }
    if (n>size) n= size;
    fill(offset-baseOffset, n, data->fillValue, type, buf);
    offset += n;
    size -= n;
  }

  if (size>0 && extent_offset - data->shift < data->upstream_extent) {
    /* We are in a block of upstream data */
    long long upstreamOffset= 
      ((n_full_extents*data->upstream_extent+(extent_offset-data->shift))
       *data->fast_blksize) + fast_blk_offset;
    n= 0;
    if (fast_blk_offset>0) {
      n += data->fast_blksize - fast_blk_offset;
      fast_blk_offset= 0;
      extent_offset += 1;
    }
    if (extent_offset - data->shift < data->upstream_extent) { 
      n += ((data->upstream_extent + data->shift - extent_offset)
	    *data->fast_blksize);
      extent_offset= data->shift + data->upstream_extent;
    }    
    if (n>size) n= size;
    src= self->owner->sinkArray[0]->source;
    switch (type) {
    case MRI_UNSIGNED_CHAR:
      nGot= src->pGetUInt8Chunk(self->owner->sinkArray[0]->source,
				n, upstreamOffset, (char*)buf);
      break;
    case MRI_SHORT:
      nGot= src->pGetInt16Chunk(self->owner->sinkArray[0]->source,
				n, upstreamOffset, (short*)buf + (offset-baseOffset));
      break;
    case MRI_INT:
      nGot= src->pGetInt32Chunk(self->owner->sinkArray[0]->source,
				n, upstreamOffset, (int*)buf + (offset-baseOffset));
      break;
    case MRI_LONGLONG:
      nGot= src->pGetInt64Chunk(self->owner->sinkArray[0]->source,
				n, upstreamOffset, 
				(long long*)buf + (offset-baseOffset));
      break;
    case MRI_FLOAT:
      nGot= src->pGetFloat32Chunk(self->owner->sinkArray[0]->source,
				  n, upstreamOffset, (float*)buf + (offset-baseOffset));
      break;
    case MRI_DOUBLE:
      nGot= src->pGetFloat64Chunk(self->owner->sinkArray[0]->source,
				  n, upstreamOffset, 
				  (double*)buf + (offset-baseOffset));
      break;
    }
    offset += nGot;
    if (nGot != n) {
      /* Upstream operations produced less than we'd like */
      if (self->owner->debug) 
	fprintf(stderr,"PadTool: end of upstream data; returning %d\n",
		(long)(offset - baseOffset));
      return (long)(offset - baseOffset);
    }
    size -= n;
  }

  /* Carry out the final block of fill */
  if (size>0 && extent_offset < data->extent) {
    n= 0;
    if (fast_blk_offset>0) {
      n += data->fast_blksize - fast_blk_offset;
      fast_blk_offset= 0;
      extent_offset += 1;
    }
    if (extent_offset < data->extent) {
      n += (data->extent - extent_offset)*data->fast_blksize;
      extent_offset= data->extent;
    }
    if (n>size) n= size;
    fill(offset-baseOffset, n, data->fillValue, type, buf);
    offset += n;
  }
  if (self->owner->debug) 
    fprintf(stderr,"PadTool: returning %d\n",(long)(offset - baseOffset));
  return (long)(offset-baseOffset);
}

static long getUInt8Chunk(DataSource* self,
				long size, long long offset, char* buf)
{
  return getTypeChunk(self, MRI_UNSIGNED_CHAR, size, offset, buf);
}

static long getInt16Chunk(DataSource* self,
				long size, long long offset, short* buf)
{
  return getTypeChunk(self, MRI_SHORT, size, offset, buf);
}

static long getInt32Chunk(DataSource* self,
				long size, long long offset, int* buf)
{
  return getTypeChunk(self, MRI_INT, size, offset, buf);
}

static long getInt64Chunk(DataSource* self,
				long size, long long offset, long long* buf)
{
  return getTypeChunk(self, MRI_LONGLONG, size, offset, buf);
}

static long getFloat32Chunk(DataSource* self,
				long size, long long offset, float* buf)
{
  return getTypeChunk(self, MRI_FLOAT, size, offset, buf);
}

static long getFloat64Chunk(DataSource* self,
				long size, long long offset, double* buf)
{
  return getTypeChunk(self, MRI_DOUBLE, size, offset, buf);
}

static DataSource* createPadSource( Tool* owner )
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
  PadData* data= (PadData*)self->hook;
  const char* dimstr= getSourceDims(upstreamSource);

  /* We require that the dimension be in the input stream, and that
   * the input extent be short enough.
   */
  if (!strchr(dimstr,data->dim)) return 0;

  if (data->upstream_extent + data->shift > data->extent) return 0;

  return 1;
}

static int init(Tool* self)
{
  DataSource* mySource= self->sourceArray[0];
  DataSource* upstreamSource= self->sinkArray[0]->source;
  PadData* data= (PadData*)self->hook;

  if (!baseToolInit(self)) return 0;
  kvCopyUniqueExceptHashes( mySource->attr, upstreamSource->attr );
  mySource->pSetName(mySource, upstreamSource->pGetName(upstreamSource));

  data->upstream_extent= getSourceDimExtent(mySource, (char)data->dim);

  /* Here we edit the attributes of the output to comply with the 
   * pad rules.
   */
  if (structureCheck(self,upstreamSource)) {
    setSourceDimExtent(mySource, (char)data->dim, data->extent);
    calcSourceBlockSizes(upstreamSource, getSourceDims(upstreamSource), 
			 (char)data->dim,
			 &(data->fast_blksize),&(data->slow_blksize));
    if (self->debug) {
      fprintf(stderr,"PadTool: dim %c, shift %ld, extent %ld\n",
	      data->dim, data->shift, data->extent);
      fprintf(stderr,
	      "PadTool: upstream_extent %ld, fast_blksize %ld, slow_blksize %lld\n",
	      data->upstream_extent, data->fast_blksize, data->slow_blksize);
    }
  }
  else return 0;
  return 0; /* NOTREACHED */
}

Tool* createPadTool(Arena* arena, const char* dim, const int extent,
		    const int shift, const double fillValue ) {
  Tool* result= createBaseTool(arena);
  PadData* data= NULL;
  if (!(data=(PadData*)malloc(sizeof(PadData))))
    Abort("Unable to allocate %d bytes!\n",sizeof(PadData));
  if (extent<=0) Abort("createPadTool: invalid extent %d!\n",extent);
  if (shift <0) Abort("createPadTool: invalid shift %d!\n",shift);
  result->hook= data;
  data->dim= *dim;
  data->shift= shift;
  data->extent= extent;
  data->fast_blksize= 0;
  data->slow_blksize= 0;
  data->fillValue= fillValue;
  result->pInit= init;
  result->pAddSink( result, createBaseSink(result) );
  result->pAddSource( result, createPadSource(result) );
  result->typeName= "pad";
  return result;
}


