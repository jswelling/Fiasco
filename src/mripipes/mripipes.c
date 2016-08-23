/************************************************************
 *	mripipes.c                                          *
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
 *		2/04: Written by Joel Welling               *
 ************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <mri.h>
#include <fmri.h>
#include <misc.h>
#include <fexceptions.h>
#include <mripipes.h>

static char rcsid[] = "$Id: mripipes.c,v 1.12 2007/06/21 23:26:51 welling Exp $";

/*************
 * Notes-
 *************/

#define SOURCE_ARRAY_LENGTH_INCR 10
#define SINK_ARRAY_LENGTH_INCR 10

static int baseSourceInit(DataSource* self)
{
  self->initialized= 1;
  return 1;
}

static int baseSourceConnect(DataSource* self, DataSink* sink)
{
  self->sink= sink;
  return 1;
}

static int baseSinkConnect(DataSink* self, DataSource* source_in)
{
  self->source= source_in;
  self->source->pConnect(self->source,self);
  return 1;
}

static int baseSinkInit(DataSink* self)
{
  if (!self->source) return 0;
  self->initialized= 1;
  return 1;
}

static void baseSinkSetName( DataSink* self, const char* name )
{
  if (self->name) free(self->name);
  self->name= strdup(name);
}

static const char* baseSinkGetName( DataSink* self )
{
  return (const char*)(self->name);
}

static void setName( DataSource* self, const char* name )
{
  if (self->name) free(self->name);
  self->name= strdup(name);
}

static const char* getName( DataSource* self )
{
  return (const char*)(self->name);
}

static void baseSourceDestroySelf( DataSource* self )
{
  if (self->hook) free(self->hook);
  if (self->name) free(self->name);
  free(self);
}

static long baseGetUInt8Chunk( DataSource* self, long size, long long offset,
			       char* buf )
{
  Abort("Tool type <%s> does not implement getUInt8Chunk!\n",
	self->owner->typeName);
  return 0;
}

static long baseGetInt16Chunk( DataSource* self, long size, long long offset,
			short* buf )
{
  Abort("Tool type <%s> does not implement getInt16Chunk!\n",
	self->owner->typeName);
  return 0;
}

static long baseGetInt32Chunk( DataSource* self, long size, long long offset,
			int* buf )
{
  Abort("Tool type <%s> does not implement getInt32Chunk!\n",
	self->owner->typeName);
  return 0;
}

static long baseGetInt64Chunk( DataSource* self, long size, long long offset,
			long long* buf )
{
  Abort("Tool type <%s> does not implement getInt64Chunk!\n",
	self->owner->typeName);
  return 0;
}

static long baseGetFloat32Chunk( DataSource* self, long size, long long offset,
			float* buf )
{
  Abort("Tool type <%s> does not implement getFloat32Chunk!\n",
	self->owner->typeName);
  return 0;
}

static long baseGetFloat64Chunk( DataSource* self, long size, long long offset,
			  double* buf )
{
  Abort("Tool type <%s> does not implement getFloat64Chunk!\n",
	self->owner->typeName);
  return 0;
}

DataSource* createBaseSource(Tool* owner)
{
  DataSource* result;
  if (!(result=(DataSource*)malloc(sizeof(DataSource)))) {
    fprintf(stderr,"malloc failed!\n");
    exit(-1);
  }
 
  result->pInit= baseSourceInit;
  result->pConnect= baseSourceConnect;
  result->pDestroySelf= baseSourceDestroySelf;
  result->pGetUInt8Chunk= baseGetUInt8Chunk;
  result->pGetInt16Chunk= baseGetInt16Chunk;
  result->pGetInt32Chunk= baseGetInt32Chunk;
  result->pGetInt64Chunk= baseGetInt64Chunk;
  result->pGetFloat32Chunk= baseGetFloat32Chunk;
  result->pGetFloat64Chunk= baseGetFloat64Chunk;
  result->pSetName= setName;
  result->pGetName= getName;
  result->sink= NULL;
  result->owner= owner;
  result->attr= kvFactory(KV_DEFAULT_SIZE);
  result->name= NULL;
  result->hook= NULL;
  result->initialized= 0;
  return result;
}

static void baseSinkDestroySelf( DataSink* self )
{
  if (self->hook) free(self->hook);
  if (self->name) free(self->name);
  free(self);
}

DataSink* createBaseSink(Tool* owner)
{
  DataSink* result;
  if (!(result=(DataSink*)malloc(sizeof(DataSink)))) {
    fprintf(stderr,"malloc failed!\n");
    exit(-1);
  }
 
  result->pInit= baseSinkInit;
  result->pConnect= baseSinkConnect;
  result->pDestroySelf= baseSinkDestroySelf;
  result->pSetName= baseSinkSetName;
  result->pGetName= baseSinkGetName;
  result->source= NULL;
  result->owner= owner;
  result->name= NULL;
  result->hook= NULL;
  result->initialized= 0;
  return result;
}

static void addSource( Tool* t, DataSource* src )
{
  if (t->sourceArrayLength <= t->nSources) {
    int newLen= t->sourceArrayLength+SOURCE_ARRAY_LENGTH_INCR;
    DataSource** newArray;
    if (!(newArray=(DataSource**)malloc(newLen*sizeof(DataSource*))))
      Abort("Unable to allocate %d bytes!\n",newLen*sizeof(DataSource*));
    if (t->sourceArray) {
      int i;
      for (i=0; i<t->nSources; i++) newArray[i]= t->sourceArray[i];
      free(t->sourceArray);
    }
    t->sourceArray= newArray;
    t->sourceArrayLength= newLen;
  }
  t->sourceArray[t->nSources++]= src;
}

static void addSink( Tool* t, DataSink* sink )
{
  if (t->sinkArrayLength <= t->nSinks) {
    int newLen= t->sinkArrayLength+SINK_ARRAY_LENGTH_INCR;
    DataSink** newArray;
    if (!(newArray=(DataSink**)malloc(newLen*sizeof(DataSink*))))
      Abort("Unable to allocate %d bytes!\n",newLen*sizeof(DataSink*));
    if (t->sinkArray) {
      int i;
      for (i=0; i<t->nSinks; i++) newArray[i]= t->sinkArray[i];
      free(t->sinkArray);
    }
    t->sinkArray= newArray;
    t->sinkArrayLength= newLen;
  }
  t->sinkArray[t->nSinks++]= sink;
}

static DataSource* getSourceByName( Tool* t, const char* nm )
{
  int i;
  for (i=0; i<t->nSources; i++)
    if (!strcmp(t->sourceArray[i]->name,nm)) return t->sourceArray[i];
  return NULL;
}

static DataSink* getSinkByName( Tool* t, const char* nm )
{
  int i;
  for (i=0; i<t->nSinks; i++)
    if (!strcmp(t->sinkArray[i]->name,nm)) return t->sinkArray[i];
  return NULL;
}

void baseToolDestroySelf( Tool* self )
{
  int i;
  for (i=0; i<self->nSources; i++) {
    DataSource* src= self->sourceArray[i];
    src->pDestroySelf(src);
  }
  for (i=0; i<self->nSinks; i++) {
    DataSink* sink= self->sinkArray[i];
    sink->pDestroySelf(sink);
  }
  if (self->hook) free(self->hook);
  free(self);
}

int baseToolInit( Tool* self ) {
  int i;
  /* Initialize all sources and sinks */
  for (i=0; i<self->nSources; i++) {
    DataSource* src= self->sourceArray[i];
    if (!src->pInit(src)) return 0;
  }
  for (i=0; i<self->nSinks; i++) {
    DataSink* sink= self->sinkArray[i];
    if (!sink->pInit(sink)) return 0;
  }
  self->initialized= 1;
  return 1;
}

int baseToolExecute( Tool* self ) {
  return 1;
}

Tool* createBaseTool(Arena* arena) {
  Tool* result;
  if (!(result=(Tool*)malloc(sizeof(Tool)))) {
    Abort("Unable to allocate %d bytes!\n",sizeof(Tool));
    fprintf(stderr,"malloc failed!\n");
    exit(-1);
  }
  result->pInit= baseToolInit;
  result->pExecute= baseToolExecute;
  result->sourceArray= NULL;
  result->nSources= 0;
  result->sourceArrayLength= 0;
  result->sinkArray= NULL;
  result->nSinks= 0;
  result->sinkArrayLength= 0;
  result->pAddSource= addSource;
  result->pAddSink= addSink;
  result->pGetSourceByName= getSourceByName;
  result->pGetSinkByName= getSinkByName;
  result->pDestroySelf= baseToolDestroySelf;
  result->hook= NULL;
  result->typeName= "baseTool";
  result->initialized= 0;
  result->verbose= 0;
  result->debug= 0;
  result->owner= arena;
  arena->pAddTool(arena,result);
  return result;
}

static void arenaDestroyToolsInList(void* p)
{
  Tool* t= (Tool*)p;
  t->pDestroySelf(t);
}

static void arenaDestroySelf(Arena* a)
{
  if (a->tools) slist_destroy(a->tools, NULL);
}

static void arenaAddTool(Arena* a, Tool* t)
{
  slist_append(a->tools,t);
}

static int recursiveInitialize(Arena* a, Tool* t)
{
  if (t->initialized) return 1;
  if (t->nSinks != 0) {
    int i;
    for (i=0; i<t->nSinks; i++) {
      if (t->sinkArray[i]->source != NULL) {
	Tool* upstreamTool= t->sinkArray[i]->source->owner;
	if (upstreamTool != NULL) {
	  if (a->debug)
	    Message("recursing up to <%s>'s sink %d's source\n",t->typeName,i);
	  if (!recursiveInitialize(a,upstreamTool)) {
	    return 0;
	  }
	}
      }
    }
  }
  if (a->verbose) Message("Initializing a %s\n",t->typeName);
  if (!t->pInit(t)) {
    if (a->verbose) Message("<%s> init failed!\n",t->typeName);
    return 0;
  }
  if (t->nSinks != 0) {
    int i;
    for (i=0; i<t->nSinks; i++) {
      if (t->sinkArray[i]->source == NULL) 
	Abort("Tool <%s> has an unconnected sink!\n",t->typeName);
    }
  }
  return 1;
}

static int arenaInit(Arena* self)
{
  slist_totop(self->tools);
  while (!slist_atend(self->tools)) {
    Tool* thisTool= (Tool*)slist_next(self->tools);
    if (thisTool->nSources==0) {
      if (!(self->drain)) self->drain= thisTool;
      else Abort("This network has more than one drain!\n");
    }
  }
  if (!recursiveInitialize(self,self->drain)) return 0;
  return 1;
}

static int arenaExecute(Arena* self)
{
  if (!self->drain)
    Abort("This network has not been initialized!\n");
  return self->drain->pExecute(self->drain);
}

Arena* createArena() {
  Arena* result;
  if (!(result=(Arena*)malloc(sizeof(Arena))))
    Abort("Unable to allocate %d bytes!\n",sizeof(Arena));
  result->tools= slist_create();
  result->drain= NULL;
  result->pDestroySelf= arenaDestroySelf;
  result->pAddTool= arenaAddTool;
  result->pInit= arenaInit;
  result->pExecute= arenaExecute;
  result->verbose= 0;
  result->debug= 0;
  return result;
}

int getSourceDimExtent( DataSource* source, const char dim )
{
  char buf[64];
  sprintf(buf,"extent.%c",dim);
  if (!kvLookup(source->attr,buf))
    Abort("Source for %s has no %s info!\n",source->owner->typeName,buf);
  return atoi(kvGetString(source->attr,buf));
}

void setSourceDimExtent( DataSource* source, const char dim, const int extent )
{
  char buf[64];
  char buf2[64];
  sprintf(buf,"extent.%c",dim);
  sprintf(buf2,"%d",extent);
  kvDefString(source->attr,buf,buf2);
}

const char* getSourceDims( DataSource* source )
{
  if (!kvLookup(source->attr,"dimensions"))
    Abort("Source for %s has no dimensions info!\n",
	  source->owner->typeName);
  return kvGetString(source->attr,"dimensions");
}

void setSourceDims( DataSource* source, const char* dimstr )
{
  const char* here= NULL;
  /* Verify that there are no duplicated dimensions */
  for (here=dimstr; *(here+1); here++)
    if (strchr(here+1,*here))
      pipeAbort("Tried to give a pipe the invalid dimension string %s!\n",
		dimstr);
  
  kvDefString(source->attr, "dimensions", dimstr);
}

int getSourceDataType(DataSource* source)
{
  const char* typestr= NULL;
  if (!(kvLookup(source->attr,"datatype"))) 
    Abort("Source for %s has no datatype info!\n",source->owner->typeName);
  typestr= kvGetString(source->attr,"datatype");
  if (!strcmp(typestr,"uint8")) return MRI_UNSIGNED_CHAR;
  else if (!strcmp(typestr,"int16")) return MRI_SHORT;
  else if (!strcmp(typestr,"int32")) return MRI_INT;
  else if (!strcmp(typestr,"int64")) return MRI_LONGLONG;
  else if (!strcmp(typestr,"float32")) return MRI_FLOAT;
  else if (!strcmp(typestr,"float64")) return MRI_DOUBLE;
  return 0;
}

void forceGetAllFloat64(DataSource* source, long size, long long offset,
			double* buf)
{
  long nToGo= size;
  while (nToGo) {
    long long nThisBlock= source->pGetFloat64Chunk(source,nToGo,offset,buf);
    if (nThisBlock==0)
      Abort("Output %s of %s provided no data\n",
	    source->name,source->owner->typeName);
    nToGo -= nThisBlock;
    buf += nThisBlock;
    offset += nThisBlock;
  }
}

void calcSourceBlockSizes(DataSource* source, const char* dimstr,
			  const char selected_dim,
			  long* fast_blocksize_out, 
			  long long* slow_blocksize_out ) {
  /* This routine will fail if selected_dim is not in dimstr! */
  long fast_blocksize;
  long long slow_blocksize;
  const char* this_dim;

  fast_blocksize= 1;
  this_dim= dimstr;
  while (*this_dim != selected_dim) 
    fast_blocksize *= getSourceDimExtent(source,*this_dim++);

  this_dim++; /* step over selected dim */

  slow_blocksize= 1;
  while (*this_dim)
    slow_blocksize *= getSourceDimExtent(source, *this_dim++);

  *fast_blocksize_out= fast_blocksize;
  *slow_blocksize_out= slow_blocksize;
}

void pipeAbort(char* fmt, ...)
{
  va_list args;
  char buf[256];

  va_start(args, fmt);
  vsnprintf(buf, sizeof(buf), fmt, args);
  va_end(args);
  fex_raiseException(EXCEPTION_BASE,buf);
}

