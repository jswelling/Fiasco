/************************************************************
 *      mri_file_input_tool.c                               *
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

static char rcsid[] = "$Id: mri_file_input_tool.c,v 1.6 2007/06/21 23:23:28 welling Exp $"; 

typedef struct FITdata_struct {
  char* fname;
  MRI_Dataset* ds;
} FITdata;

static long orphanGetUInt8Chunk(DataSource* self,
				long size, long long offset, char* buf)
{
  Abort("Attempted to read from the orphans!\n");
  return 0;
}

static long orphanGetInt16Chunk(DataSource* self,
				long size, long long offset, short* buf)
{
  Abort("Attempted to read from the orphans!\n");
  return 0;
}

static long orphanGetInt32Chunk(DataSource* self,
				long size, long long offset, int* buf)
{
  Abort("Attempted to read from the orphans!\n");
  return 0;
}

static long orphanGetInt64Chunk(DataSource* self,
				long size, long long offset, long long* buf)
{
  Abort("Attempted to read from the orphans!\n");
  return 0;
}

static long orphanGetFloat32Chunk(DataSource* self,
				long size, long long offset, float* buf)
{
  Abort("Attempted to read from the orphans!\n");
  return 0;
}

static long orphanGetFloat64Chunk(DataSource* self,
				long size, long long offset, double* buf)
{
  Abort("Attempted to read from the orphans!\n");
  return 0;
}

static long getUInt8Chunk(DataSource* self, 
			  long size, long long offset, char* buf )
{
  (void)mri_read_chunk(((FITdata*)(self->owner->hook))->ds,self->name,
		       size, offset, MRI_UNSIGNED_CHAR, buf);
  return size;
}

static long getInt16Chunk(DataSource* self, 
			  long size, long long offset, short* buf )
{
  (void)mri_read_chunk(((FITdata*)(self->owner->hook))->ds,self->name,
		       size, offset, MRI_SHORT, buf);
  return size;
}

static long getInt32Chunk(DataSource* self, 
			  long size, long long offset, int* buf )
{
  (void)mri_read_chunk(((FITdata*)(self->owner->hook))->ds,self->name,
		       size, offset, MRI_INT, buf);
  return size;
}

static long getInt64Chunk(DataSource* self, 
			  long size, long long offset, long long* buf )
{
  (void)mri_read_chunk(((FITdata*)(self->owner->hook))->ds,self->name,
		       size, offset, MRI_LONGLONG, buf);
  return size;
}

static long getFloat32Chunk(DataSource* self, 
			    long size, long long offset, float* buf )
{
  (void)mri_read_chunk(((FITdata*)(self->owner->hook))->ds,self->name,
		       size, offset, MRI_FLOAT, buf);
  return size;
}

static long getFloat64Chunk(DataSource* self, 
			    long size, long long offset, double* buf )
{
  if (self->owner->debug) 
    fprintf(stderr,"read %ld MRI_DOUBLES from %lld, %s\n",
	    size, offset, 
	    ((FITdata*)(self->owner->hook))->fname);
  (void)mri_read_chunk(((FITdata*)(self->owner->hook))->ds,self->name,
		       size, offset, MRI_DOUBLE, buf);
  return size;
}

static DataSource* createOrphanSource( Tool* owner, 
					    const char* chunkname )
{
  DataSource* result= createBaseSource(owner);
  result->pSetName(result, chunkname);
  result->pGetUInt8Chunk= orphanGetUInt8Chunk;
  result->pGetInt16Chunk= orphanGetInt16Chunk;
  result->pGetInt32Chunk= orphanGetInt32Chunk;
  result->pGetInt64Chunk= orphanGetInt64Chunk;
  result->pGetFloat32Chunk= orphanGetFloat32Chunk;
  result->pGetFloat64Chunk= orphanGetFloat64Chunk;
  return result;
}

static DataSource* createMRIFileSource( Tool* owner, 
					    const char* chunkname )
{
  DataSource* result= createBaseSource(owner);
  result->pSetName(result, chunkname);
  result->pGetUInt8Chunk= getUInt8Chunk;
  result->pGetInt16Chunk= getInt16Chunk;
  result->pGetInt32Chunk= getInt32Chunk;
  result->pGetInt64Chunk= getInt64Chunk;
  result->pGetFloat32Chunk= getFloat32Chunk;
  result->pGetFloat64Chunk= getFloat64Chunk;
  return result;
}

Tool* createMRIFileInputTool(Arena* arena, const char* fname) {
  Tool* result= createBaseTool(arena);
  FITdata* data;
  DataSource* orphanSource= NULL;
  DataSource* currentSource= NULL;
  char currentChunkKey[256];
  const char* key;
  if (!(data=(FITdata*)malloc(sizeof(FITdata))))
    Abort("Unable to allocate %d bytes!\n",sizeof(FITdata));
  result->hook= data;
  result->typeName= "MRIFileInput";
  data->fname= strdup(fname);
  if (!(data->ds= mri_open_dataset(data->fname,MRI_READ))) {
    Message("Unable to open <%s> for reading!\n", data->fname);
    return 0;
  }
  currentChunkKey[sizeof(currentChunkKey)-1]= '\0';
  orphanSource= createMRIFileSource(result,"orphans");
  result->pAddSource(result,orphanSource);
  mri_iterate_over_keys(data->ds);
  while (key= mri_next_key(data->ds)) {
    const char* val= mri_get_string(data->ds,key);
    if (key[0]=='!') continue; /* ! denotes comment in Pgh MRI */
    if (!strcmp(val,"[chunk]")) {
      snprintf(currentChunkKey,sizeof(currentChunkKey)-1,"%s.",key);
      currentSource= createMRIFileSource(result,key);
      result->pAddSource(result,currentSource);
    }
    else {
      if (currentSource != NULL 
	  && !strncmp(currentChunkKey,key,strlen(currentChunkKey))) {
	const char* relativeKey= key+strlen(currentChunkKey);
	if (strcmp(relativeKey,"size")
	    && strcmp(relativeKey,"little_endian")
	    && strcmp(relativeKey,"offset")) /* keys maintained by libmri */
	  kvDefString(currentSource->attr,relativeKey,val);
      }
      else {
	if (strcmp(key,"size")
	    && strcmp(key,"little_endian")
	    && strcmp(key,"offset")) /* keys maintained by libmri */
	  kvDefString(orphanSource->attr,key,val);
      }
    }
  }
  return result;
}


