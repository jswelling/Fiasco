/************************************************************
 *      devnull_tool.c                                      *
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

static char rcsid[] = "$Id: devnull_tool.c,v 1.5 2004/06/05 01:21:08 welling Exp $"; 

#define BLOCKSIZE (1024*1024)

static DataSink* createDevnullSink( Tool* owner )
{
  DataSink* result= createBaseSink(owner);
  return result;
}

static int execute( Tool* self ) {
  DataSink* mySink= self->sinkArray[0];
  float *inbuf= NULL;
  long long totalSize= 1;
  const char* dimstr= kvGetString(mySink->source->attr,"dimensions");
  const char* here= dimstr;
  long long offset= 0;
  fprintf(stderr,"execute called!\n");
  if (!baseToolExecute(self)) return 0; 
  while (*here) {
    char buf[64];
    int thisExtent;
    sprintf(buf,"extent.%c",*here);
    thisExtent= atoi(kvGetString(mySink->source->attr,buf));
    totalSize*=thisExtent;
    here++;
  }
  fprintf(stderr,"totalSize is %lld\n",totalSize);
  if (!(inbuf=(float*)malloc(BLOCKSIZE*sizeof(float))))
    Abort("Unable to allocate %d bytes!\n",BLOCKSIZE*sizeof(float));
  while (totalSize>0) {
    long nThisChunk= (totalSize>BLOCKSIZE)?BLOCKSIZE:totalSize;
    nThisChunk= 
      mySink->source->pGetFloat32Chunk(mySink->source, nThisChunk, offset,
				       inbuf);
    fprintf(stderr,"Got %d: %f %f ... %f\n",nThisChunk,inbuf[0],inbuf[1],
	    inbuf[nThisChunk-1]);
    totalSize -= nThisChunk;
    offset += nThisChunk;
  }
  free(inbuf);
  return 1;
}

Tool* createDevnullTool(Arena* arena) {
  Tool* result= createBaseTool(arena);
  result->pExecute= execute;
  result->pAddSink( result, createDevnullSink(result) );
  result->typeName= "devnull";
  return result;
}


