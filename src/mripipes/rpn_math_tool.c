/************************************************************
 *      rpn_math_tool.c                               *
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

static char rcsid[] = "$Id: rpn_math_tool.c,v 1.4 2005/08/09 23:11:06 welling Exp $"; 

typedef struct rpn_data_struct {
  RpnEngine* re;
  char* script;
  int complexFlag;
  long obufSize;
  long long obufOffset;
  long obufValidLength;
  double* obuf;
} RPNData;

static DataSink* createRPNSink(Tool* owner);

static int connect(DataSink* self, DataSource* source_in)
{
  DataSink* parent= (DataSink*)(self->hook);
  parent->pConnect(self, source_in);
  self->owner->pAddSink(self->owner,createRPNSink(self->owner));
  return 1;
}

static void sinkDestroySelf(DataSink* self )
{
  if (self->hook) {
    DataSink* parent= (DataSink*)(self->hook);
    parent->pDestroySelf(parent);
    self->hook= NULL;
  }
  free(self);
}

static DataSink* createRPNSink( Tool* owner )
{
  DataSink* result= createBaseSink(owner);
  result->hook= createBaseSink(owner); /* we'll need some base methods later */
  result->pConnect= connect;
  result->pDestroySelf= sinkDestroySelf;
  return result;
}

static void fillBuffer(Tool* self, long size, long long offset)
{
  RPNData* data= (RPNData*)self->hook;
  RpnEngine* re= data->re;
  double* here;

  /* If data->obuf is invalid or contains non-applicable data, 
   * calculate the block base offset which is closest to the 
   * requested offset and fill the obuf from that offset.
   */
  if (!data->obufValidLength || offset<data->obufOffset 
      || offset+size >=data->obufOffset+data->obufValidLength) {
    data->obufValidLength= 0;
    if (re->complexFlag) {
      long i;
      double* vals;
      long n= (size>2*RPN_CHUNKSIZE) ? RPN_CHUNKSIZE:(long)((size+1)/2);
      if (data->obufSize != 2*RPN_CHUNKSIZE)
	Abort("rpn_math_tool internal error: bufsize doesn't match 2*chunksize!\n");
      vals= rpnRun( re, n, offset );
      here= data->obuf;
      for (i=0; i<n; i++) {
	*here++= vals[i];
	*here++= vals[RPN_CHUNKSIZE+i];
      }
      data->obufValidLength= n;
    }
    else {
      double* vals;
      long n= (size>RPN_CHUNKSIZE) ? RPN_CHUNKSIZE:size;
      if (data->obufSize != RPN_CHUNKSIZE)
	Abort("rpn_math_tool internal error: bufsize doesn't match chunksize!\n");
      vals= rpnRun( re, n, offset );
      memcpy(data->obuf,vals,n*sizeof(double));
      data->obufValidLength= n;
    }
    data->obufOffset= offset;
  }
  else if (self->debug) fprintf(stderr,"No new calculation needed\n");
} 

static long getFloat32Chunk(DataSource* self,
				long size, long long offset, float* buf)
{
  RPNData* data= (RPNData*)(self->owner->hook);
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
  RPNData* data= (RPNData*)(self->owner->hook);
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
  for (i=0; i<n; i++) buf[i]= data->obuf[i+shift];
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
  RPNData* data= (RPNData*)self->hook;
  if (data->re) rpnDestroyEngine(data->re);
  if (data->script) free(data->script);
  baseToolDestroySelf(self);
}

static const char* getDimensionsCB( const int which, void* hook )
{
  Tool* owner= (Tool*)hook;
  return getSourceDims( owner->sinkArray[which]->source );
}

static const long getDimExtentCB( const int which, const char dim, 
			   void* hook )
{
  Tool* owner= (Tool*)hook;
  return getSourceDimExtent( owner->sinkArray[which]->source, dim );
}

static void inputCB( const int which, const long n, 
		     const long long offset, double* buf, 
		     void* hook )
{
  Tool* owner= (Tool*)hook;
  DataSource* src= owner->sinkArray[which]->source;
  if (owner->debug) 
    fprintf(stderr,"rpn_math_tool: inputCB: n= %ld, offset %lld\n",n,offset);
  forceGetAllFloat64(src, n, offset, buf);
}

static void inputComplexCB( const int which, const long n, 
			    const long long offset, double* buf1, 
			    double* buf2, void* hook )
{
  static double* tbuf= NULL;
  static long tbufLength= 0;
  Tool* owner= (Tool*)hook;
  DataSource* src= owner->sinkArray[which]->source;
  int i;

  if (tbufLength<2*n) {
    if (tbufLength) free(tbuf);
    if (!(tbuf= (double*)malloc(2*n*sizeof(double))))
      Abort("rpn_math_tool: inputComplexCB: unable to allocate %d bytes!\n",
	    2*n*sizeof(double));
  }
  if (owner->debug) 
    fprintf(stderr,"rpn_math_tool: inputComplexCB: n= %ld, offset %lld\n",
	    n,offset);
  forceGetAllFloat64(src, n, offset, tbuf);
  for (i=0; i<n; i++) {
    buf1[i]= tbuf[2*i];
    buf2[i]= tbuf[2*i+1];
  }
}

static int  missingCB( const long z, const long t, void* hook )
{
  return 0;
}

static int init( Tool* self )
{
  RPNData* data= (RPNData*)self->hook;
  DataSink* lastSink= self->sinkArray[self->nSinks - 1];
  DataSource* firstInput= self->sinkArray[0]->source;
  DataSource* myOutput= self->sourceArray[0];
  if (!lastSink->source) {
    self->nSinks -= 1;
    lastSink->pDestroySelf(lastSink);
  }
  if (!baseToolInit(self)) return 0;
  data->re= createRpnEngine(self->nSinks, self,
			    getDimensionsCB, getDimExtentCB,
			    inputCB, inputComplexCB, missingCB);
  rpnSetOutputFlag(data->re,1);
  rpnSetComplex(data->re,data->complexFlag);
  if (self->debug) rpnSetDebug(data->re,1);
  if (self->verbose) rpnSetVerbose(data->re,1);
  if (!rpnInit(data->re)) {
    fprintf(stderr,"%s: %s\n",self->typeName,rpnGetErrorString(data->re));
    return 0;
  }
  if (!rpnCompile(data->re,data->script)) {
    fprintf(stderr,"%s: %s\n",self->typeName,rpnGetErrorString(data->re));
    return 0;
  }
  if (data->obufSize) {
    free(data->obuf);
    data->obufSize= 0;
  }
  if (data->re->complexFlag) {
    if (!(data->obuf=(double*)malloc(2*RPN_CHUNKSIZE*sizeof(double))))
      Abort("Unable to allocate %d bytes!\n",2*RPN_CHUNKSIZE*sizeof(double));
    data->obufSize= 2*RPN_CHUNKSIZE;
  }
  else {
    if (!(data->obuf=(double*)malloc(RPN_CHUNKSIZE*sizeof(double))))
      Abort("Unable to allocate %d bytes!\n",RPN_CHUNKSIZE*sizeof(double));
    data->obufSize= RPN_CHUNKSIZE;
  }
  kvCopyUniqueExceptHashes( myOutput->attr, firstInput->attr );
  myOutput->pSetName(myOutput,firstInput->pGetName(firstInput));

  return 1;
}

Tool* createRpnMathTool(Arena* arena, const char* script_in) {
  Tool* result= createBaseTool(arena);
  RPNData* data= NULL;
  if (!(data=(RPNData*)malloc(sizeof(RPNData))))
    Abort("Unable to allocate %d bytes!\n",sizeof(RPNData));
  result->hook= data;
  data->re= NULL;
  data->obufSize= 0;
  data->obuf= NULL;
  data->obufValidLength= 0;
  data->script= strdup(script_in);
  data->complexFlag= 0;
  result->pDestroySelf= destroySelf;
  result->pInit= init;
  result->pAddSink( result, createRPNSink(result) );
  result->pAddSource( result, createSource(result) );
  result->typeName= "rpn_math";
  return result;
}

Tool* createComplexRpnMathTool(Arena* arena, const char* script_in) {
  Tool* result= createRpnMathTool(arena, script_in);
  RPNData* data= (RPNData*)(result->hook);
  data->complexFlag= 1;
  return result;
}


