/************************************************************
 *      func_2rows_unblocked_tool.c                         *
 *                                                          *
 *	Copyright (c) 2007 Pittsburgh Supercomputing Center *
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

static char rcsid[] = "$Id: func_2rows_unblocked_tool.c,v 1.2 2007/06/21 23:31:55 welling Exp $"; 

typedef struct func2unblk_data_struct {
  long long obufOffset;
  long obufSize;
  long obufValidLength;
  long nInputs;
  DataSource* leftSource;
  double* leftBuf;
  DataSource* rightSource;
  double* rightBuf;
  long nOutputs;
  double* obuf;
  FUNC2UNBLKCBFUNC func;
  void* cbData;
} Func2UnblkData;

static DataSink* createLeftSink( Tool* owner )
{
  DataSink* result= createBaseSink(owner);
  return result;
}

static DataSink* createRightSink( Tool* owner )
{
  DataSink* result= createBaseSink(owner);
  return result;
}

static void fillBuffer(Tool* self, long size, long long offset)
{
  /* offset is the location in the output stream of the start point
   * of the downstream tool's request.
   *
   * baseOffset is the location in the output stream of the block
   * containing offset.
   *
   * upstreamBaseOffset is the location in the input stream of the
   * offset corresponding to baseOffset.
   *
   * Fortunately for the case of fast_blksize==1 this calculation
   * is really trivial.
   */
  Func2UnblkData* data= (Func2UnblkData*)self->hook;

  /* If data->obuf is invalid or contains non-applicable data, 
   * calculate the block base offset which is closest to the 
   * requested offset and fill the obuf from that offset.
   */
  if (!data->obufValidLength || offset<data->obufOffset 
      || offset+size >=data->obufOffset+data->obufValidLength) {
    long long baseOffset= offset - (offset%data->nOutputs);
    long long upstreamBaseOffset= (baseOffset*data->nInputs)/data->nOutputs;

    data->obufValidLength= 0;

    if (data->obufSize != data->nOutputs)
      pipeAbort("func_2rows_unblocked_tool internal error: bufsize doesn't match chunksize!\n");
    forceGetAllFloat64(data->leftSource, data->nInputs, upstreamBaseOffset,
		       data->leftBuf);
    forceGetAllFloat64(data->rightSource, data->nInputs, upstreamBaseOffset,
		       data->rightBuf);
    if (!(data->func(data->leftBuf, data->rightBuf, data->nInputs, 
		     data->obuf, data->nOutputs, data->cbData)))
      pipeAbort("Func2UnblkTool fillBuffer: callback signaled error!");
      
    data->obufValidLength= data->nOutputs;
    data->obufOffset= baseOffset;
  }
  else if (self->debug) fprintf(stderr,"No new calculation needed\n");
} 

static long getFloat32Chunk(DataSource* self,
			    long size, long long offset, float* buf)
{
  Func2UnblkData* data= (Func2UnblkData*)(self->owner->hook);
  long shift;
  long n;
  long i;

  fillBuffer(self->owner, size, offset);

  /* The obuf is now guaranteed to have an applicable segment;
   * copy it to the output.
   */
  shift= (long)(offset - data->obufOffset);
  if (size>data->obufSize-shift) n= data->obufSize-shift;
  else n= size;
  for (i=0; i<n; i++) buf[i]= (float)data->obuf[i+shift];
  return n;
}

static long getFloat64Chunk(DataSource* self,
				long size, long long offset, double* buf)
{
  Func2UnblkData* data= (Func2UnblkData*)(self->owner->hook);
  long shift;
  long n;

  fillBuffer(self->owner, size, offset);

  /* The obuf is now guaranteed to have an applicable segment;
   * copy it to the output.
   */
  shift= (long)(offset - data->obufOffset);
  if (size>data->obufSize-shift) n= data->obufSize-shift;
  else n= size;
  (void)memcpy( buf, &(data->obuf[shift]), n*sizeof(double) );
  return n;
}

static DataSource* createSource( Tool* owner )
{
  DataSource* result= createBaseSource(owner);
  result->pGetFloat32Chunk= getFloat32Chunk;
  result->pGetFloat64Chunk= getFloat64Chunk;
  return result;
}

static void destroySelf(Tool* self)
{
  Func2UnblkData* data= (Func2UnblkData*)self->hook;
  if (data->obuf) free(data->obuf);
  if (data->leftBuf) free(data->leftBuf);
  if (data->rightBuf) free(data->rightBuf);
  baseToolDestroySelf(self);
}

static void restructureDims( Tool* self, DataSource* Output,
			     DataSource* left, DataSource* right )
{
  /* This routine sets the first extent of the output to match
   * data->nOutputs.
   */
  Func2UnblkData* data= (Func2UnblkData*)self->hook;
  char buf[256];
  char key_buf[32];
  const char* dimstr= getSourceDims(left);

  snprintf(key_buf,sizeof(key_buf),"extent.%c",dimstr[0]);
  snprintf(buf,sizeof(buf),"%ld",data->nOutputs);
  kvDefString(Output->attr,key_buf,buf);
}

static int structureCheck( Tool* self, DataSource* left, DataSource* right)
{
  Func2UnblkData* data= (Func2UnblkData*)self->hook;
  const char* leftDims= getSourceDims(left);
  const char* rightDims= getSourceDims(right);
  int len= strlen(leftDims);
  int i;

  if (strlen(leftDims) != strlen(rightDims))
    pipeAbort("mismatched dim strings <%s> and <%s> in func_2rows_unblocked_tool!",
	      leftDims, rightDims);
  for (i=0; i<len; i++) {
    if (leftDims[i]!=rightDims[i])
      pipeAbort("mismatched dimensions %c and %c in func_2rows_unblocked_tool!",
		leftDims[i],rightDims[i]);
    if (getSourceDimExtent(left,leftDims[i]) 
	!= getSourceDimExtent(right,rightDims[i]))
      pipeAbort("mismatched dim lengths %d and %d for dim %c in func_2rows_unblocked_tool!",
		getSourceDimExtent(left,leftDims[i]),
		getSourceDimExtent(right,rightDims[i]),
		leftDims[i]);
  }

  return 1;
}

static int init( Tool* self )
{
  Func2UnblkData* data= (Func2UnblkData*)self->hook;
  DataSource* myOutput= self->sourceArray[0];

  data->leftSource= self->pGetSinkByName(self,"left")->source;
  data->rightSource= self->pGetSinkByName(self,"right")->source;

  if (!baseToolInit(self)) return 0;
  
  /* Verify that left and right input have correct dims and types */
  if (!structureCheck(self, data->leftSource, data->rightSource)) return 0;
  
  /* Create input buffers.  structureCheck has guaranteed that
   * both inputs require the same size buffer. 
   */
  data->nInputs= getSourceDimExtent(data->leftSource, 
				    getSourceDims(data->leftSource)[0]);
  if (data->leftBuf) {
    free(data->leftBuf);
    data->leftBuf= NULL;
  }
  if (!(data->leftBuf=(double*)malloc(data->nInputs*sizeof(double))))
    pipeAbort("Unable to allocate %d bytes!\n",data->nInputs*sizeof(double));
  if (data->rightBuf) {
    free(data->rightBuf);
    data->rightBuf= NULL;
  }
  if (!(data->rightBuf=(double*)malloc(data->nInputs*sizeof(double))))
    pipeAbort("Unable to allocate %d bytes!\n",data->nInputs*sizeof(double));

  /* Create an output buffer */
  if (data->obufSize) {
    free(data->obuf);
    data->obuf= NULL;
  }
  if (!(data->obuf=(double*)malloc(data->nOutputs*sizeof(double))))
    pipeAbort("Unable to allocate %d bytes!\n",data->nOutputs*sizeof(double));
  data->obufSize= data->nOutputs;

  /* Reshape attributes for our output */
  kvCopyUniqueExceptHashes( myOutput->attr, data->leftSource->attr );
  restructureDims( self, myOutput, data->leftSource, data->rightSource );
  myOutput->pSetName(myOutput,data->leftSource->pGetName(data->leftSource));
  
  return 1;
}

Tool* createFunc2UnblkTool(Arena* arena, FUNC2UNBLKCBFUNC func, long nOut,
			   void* cbData) {
  Tool* result= createBaseTool(arena);
  Func2UnblkData* data= NULL;
  DataSink* leftSink= NULL;
  DataSink* rightSink= NULL;

  if (!(data=(Func2UnblkData*)malloc(sizeof(Func2UnblkData))))
    pipeAbort("Unable to allocate %d bytes!\n",sizeof(Func2UnblkData));
  result->hook= data;
  data->obufSize= 0;
  data->obufValidLength= 0;
  data->nInputs= 0;
  data->leftSource= NULL;
  data->leftBuf= NULL;
  data->rightSource= NULL;
  data->rightBuf= NULL;
  data->obuf= NULL;
  data->nOutputs= nOut;
  data->func= func;
  data->cbData= cbData;
  result->pDestroySelf= destroySelf;
  result->pInit= init;
  leftSink= createLeftSink(result);
  leftSink->pSetName(leftSink,"left");
  result->pAddSink( result, leftSink );
  rightSink= createRightSink(result);
  rightSink->pSetName(rightSink,"right");
  result->pAddSink( result, rightSink );
  result->pAddSource( result, createSource(result) );
  result->typeName= "func2unblk";
  return result;
}



