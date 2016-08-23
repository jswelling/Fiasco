/************************************************************
 *      mri_file_output_tool.c                               *
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
#include <mripipes.h>

/***********************
 * Notes-
 **********************/

static char rcsid[] = "$Id: mri_file_output_tool.c,v 1.4 2005/08/09 23:11:06 welling Exp $"; 

#define BLOCKSIZE (1024*1024)

typedef struct FOTdata_struct {
  char* fname;
  MRI_Dataset* ds;
  int* typeArray;
  long long* totalSizeArray;
  long long* offsetArray;
  int* liveArray;
  void** obufArray;
} FOTdata;

static long typesize[8] = { 0, 1, 2, 4, 4, 8, 4, 8 };

static char* typename[8] = { "raw", "uint8", "int16", "int32", "float32", 
		      "float64", "long", "longlong" };

static DataSink* createMRIFileSink(Tool* owner);

static int connect(DataSink* self, DataSource* source_in)
{
  DataSink* parent= (DataSink*)(self->hook);
  parent->pConnect(self, source_in);
  self->owner->pAddSink(self->owner,createMRIFileSink(self->owner));
  return 1;
}

static void destroySelf(DataSink* self )
{
  if (self->hook) {
    DataSink* parent= (DataSink*)(self->hook);
    parent->pDestroySelf(parent);
    self->hook= NULL;
  }
  free(self);
}

static DataSink* createMRIFileSink( Tool* owner )
{
  DataSink* result= createBaseSink(owner);
  result->hook= createBaseSink(owner); /* we'll need some base methods later */
  result->pConnect= connect;
  result->pDestroySelf= destroySelf;
  return result;
}

static int customInit( Tool* self )
{
  FOTdata* data= (FOTdata*)self->hook;
  DataSink* lastSink= self->sinkArray[self->nSinks - 1];
  if (!lastSink->source) {
    self->nSinks -= 1;
    lastSink->pDestroySelf(lastSink);
  }
  if (!baseToolInit(self)) return 0;
  if (!(data->ds= mri_open_dataset(data->fname,MRI_WRITE))) {
    Message("Unable to open <%s> for reading!\n", data->fname);
    return 0;
  }
  return 1;
}

static long long calcTotalSize(DataSink* sink)
{
  long long totalSize= 1;
  const char* dimstr= NULL;
  const char* here= NULL;
  long long offset= 0;  
  if (!(kvLookup(sink->source->attr,"dimensions"))) return 0;
  here= dimstr=  kvGetString(sink->source->attr,"dimensions");
  while (*here) {
    char buf[64];
    int thisExtent;
    sprintf(buf,"extent.%c",*here);
    if (!kvLookup(sink->source->attr,buf))
      Abort("%s: sink from %s has no <%s> tag!\n",
	    sink->owner->typeName,sink->source->name,buf);
    thisExtent= atoi(kvGetString(sink->source->attr,buf));
    totalSize*=thisExtent;
    here++;
  }
  return totalSize;
}

static int customExecute( Tool* self )
{
  FOTdata* data= (FOTdata*)self->hook;
  int i;
  int someAreLive= 0;
  if (!baseToolExecute(self)) return 0;

  /* Initialize tables describing the streams */
  if (data->typeArray) free(data->typeArray);
  if (data->totalSizeArray) free(data->totalSizeArray);
  if (data->offsetArray) free(data->offsetArray);
  if (data->liveArray) free(data->liveArray);
  if (!(data->typeArray= (int*)malloc(self->nSinks*sizeof(int))))
    Abort("Unable to allocate %d bytes!\n",self->nSinks*sizeof(int));
  if (!(data->totalSizeArray=
	(long long*)malloc(self->nSinks*sizeof(long long))))
    Abort("Unable to allocate %d bytes!\n",self->nSinks*sizeof(long long));
  if (!(data->offsetArray=
	(long long*)malloc(self->nSinks*sizeof(long long))))
    Abort("Unable to allocate %d bytes!\n",self->nSinks*sizeof(long long));
  if (!(data->liveArray=(int*)malloc(self->nSinks*sizeof(int))))
    Abort("Unable to allocate %d bytes!\n",self->nSinks*sizeof(int));
  if (!(data->obufArray=(void**)malloc(self->nSinks*sizeof(void*))))
    Abort("Unable to allocate %d bytes!\n",self->nSinks*sizeof(void*));

  /* Copy the stream attributes to the file, associating with relevant
   * chunks.
   */
  for (i=0; i<self->nSinks; i++) {
    DataSink* sink= self->sinkArray[i];
    if (!strcmp(sink->source->name,"orphans")) {
      KVIterator* kvi= NULL;
      data->typeArray[i]= 0;
      data->totalSizeArray[i]= 0;
      data->offsetArray[i]= 0;
      data->liveArray[i]= 0;
      kvi= kvUniqueIteratorFactory(sink->source->attr);
      while (kvIteratorHasMorePairs(kvi)) {
	KVPair* p= kvIteratorNextPair(kvi);
	mri_set_string(data->ds,kvKey(p),
		       kvGetString(sink->source->attr,kvKey(p)));
      }
      kvDestroyIterator(kvi);
    }
    else {
      char keyBuf[256];
      KVIterator* kvi= NULL;
      data->typeArray[i]= getSourceDataType(self->sinkArray[i]->source);
      data->totalSizeArray[i]= calcTotalSize(self->sinkArray[i]);
      data->offsetArray[i]= 0;
      data->liveArray[i]= (data->totalSizeArray[i] != 0) ? 1:0;
      if (!(data->obufArray[i]=malloc(BLOCKSIZE*typesize[data->typeArray[i]])))
	Abort("Unable to allocate %d bytes!\n",
	      BLOCKSIZE*typesize[data->typeArray[i]]);
      mri_create_chunk(data->ds,sink->source->name);
      kvi= kvUniqueIteratorFactory(sink->source->attr);
      while (kvIteratorHasMorePairs(kvi)) {
	KVPair* p= kvIteratorNextPair(kvi);
	snprintf(keyBuf,sizeof(keyBuf)-1,"%s.%s",
		 sink->source->name,kvKey(p));
	mri_set_string(data->ds,keyBuf,
		       kvGetString(sink->source->attr,kvKey(p)));
      }      
      kvDestroyIterator(kvi);
    }
  }

  /* Move the data, alternating between streams. */
  someAreLive= 1;
  while (someAreLive) {
    someAreLive= 0;
    for (i=0; i<self->nSinks; i++) {
      if (data->liveArray[i]) {
	DataSink* sink= self->sinkArray[i];
	long long nToGo= data->totalSizeArray[i]-data->offsetArray[i];
	long nThisBlock= (nToGo>BLOCKSIZE)? BLOCKSIZE:nToGo;
	long nGot= 0;
	switch (data->typeArray[i]) {
	case MRI_UNSIGNED_CHAR:
	  nGot= sink->source->pGetUInt8Chunk(sink->source,
					     nThisBlock, data->offsetArray[i],
					     (char*)data->obufArray[i]);
	  break;
	case MRI_SHORT:
	  nGot= sink->source->pGetInt16Chunk(sink->source,
					     nThisBlock, data->offsetArray[i],
					     (short*)data->obufArray[i]);
	  break;
	case MRI_INT:
	  nGot= sink->source->pGetInt32Chunk(sink->source,
					     nThisBlock, data->offsetArray[i],
					     (int*)data->obufArray[i]);
	  break;
	case MRI_LONGLONG:
	  nGot= sink->source->pGetInt64Chunk(sink->source,
					     nThisBlock, data->offsetArray[i],
					     (long long*)data->obufArray[i]);
	  break;
	case MRI_FLOAT:
	  nGot= sink->source->pGetFloat32Chunk(sink->source,
					       nThisBlock, 
					       data->offsetArray[i],
					       (float*)data->obufArray[i]);
	  break;
	case MRI_DOUBLE:
	  nGot= sink->source->pGetFloat64Chunk(sink->source,
					       nThisBlock, 
					       data->offsetArray[i],
					       (double*)data->obufArray[i]);
	  break;
	}
	mri_write_chunk(data->ds,sink->source->name,
			nGot, data->offsetArray[i], 
			data->typeArray[i], data->obufArray[i]);
	if (self->debug)
	  fprintf(stderr,
		  "MRIFileOutputTool: wrote %d at %lld, chunk %s type %s\n",
		  nGot, data->offsetArray[i], 
		  sink->source->name, typename[data->typeArray[i]]);
	data->offsetArray[i] += nGot;
	if (data->offsetArray[i]<data->totalSizeArray[i]) someAreLive= 1;
      }
    }
  }
  mri_close_dataset(data->ds);
  return 1;
}

static void customDestroySelf(Tool* self)
{
  FOTdata* data= (FOTdata*)self->hook;
  if (data->typeArray) free(data->typeArray);
  if (data->totalSizeArray) free(data->totalSizeArray);
  if (data->offsetArray) free(data->offsetArray);
  if (data->liveArray) free(data->liveArray);
  if (data->obufArray) {
    int i;
    for (i=0; i<self->nSinks; i++) free(data->obufArray[i]);
    free(data->obufArray);
  }
  baseToolDestroySelf(self);
}

Tool* createMRIFileOutputTool(Arena* arena, const char* fname) {
  Tool* result= createBaseTool(arena);
  FOTdata* data;
  DataSink* firstSink= createMRIFileSink(result);
  if (!(data=(FOTdata*)malloc(sizeof(FOTdata))))
    Abort("Unable to allocate %d bytes!\n",sizeof(FOTdata));
  result->hook= data;
  result->pInit= customInit;
  result->pExecute= customExecute;
  result->pDestroySelf= customDestroySelf;
  result->typeName= "MRIFileOutput";
  data->fname= strdup(fname);
  data->ds= NULL;
  data->liveArray= NULL;
  data->typeArray= NULL;
  data->totalSizeArray= NULL;
  data->offsetArray= NULL;
  data->obufArray= NULL;
  if (!(data->ds= mri_open_dataset(data->fname,MRI_WRITE))) {
    Message("Unable to open <%s> for reading!\n", data->fname);
    return 0;
  }
  result->pAddSink(result,firstSink);
  return result;
}


