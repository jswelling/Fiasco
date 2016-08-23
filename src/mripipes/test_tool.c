/************************************************************
 *	test_tool.c                                          *
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

#include <stdlib.h>
#include <string.h>
#include <mri.h>
#include <fmri.h>
#include <misc.h>
#include <mripipes.h>

static char rcsid[] = "$Id: test_tool.c,v 1.3 2004/06/02 00:52:49 welling Exp $";

static long getInt32Chunk( DataSource* self, long size, long long offset,
			   int* buf)
{
  int* result;
  int i;
  fprintf(stderr,"getInt32 called; size %d, offset %lld\n",size,offset);
  for (i=0; i<size; i++) buf[i]= i%30;
  return size;
}

DataSource* createTestSource(Tool* owner)
{
  DataSource* result= createBaseSource(owner);
  result->pGetInt32Chunk= getInt32Chunk;
  return result;
}

DataSink* createTestSink(Tool* owner)
{
  DataSink* result= createBaseSink(owner);
  return result;
}

static int initUpstreamTool( Tool* t ) {
  fprintf(stderr,"initUpstreamTool called!\n");
  return 1;
}

static int executeUpstreamTool( Tool* t ) {
  fprintf(stderr,"executeUpstreamTool called!\n");
  return 1;
}

Tool* createUpstreamTool(Arena* arena) {
  Tool* result= createBaseTool(arena);
  result->pInit= initUpstreamTool;
  result->pExecute= executeUpstreamTool;
  result->pAddSource(result, createTestSource(result));
  result->typeName= "upstream";
  return result;
}

static int initDownstreamTool( Tool* t ) {
  fprintf(stderr,"initDownstreamTool called!\n");
  return 1;
}

static int executeDownstreamTool( Tool* t ) {
  int buf[100];
  DataSink* mySink= t->sinkArray[0];
  DataSource* partnerSource= mySink->source;
  fprintf(stderr,"executeDownstreamTool called!\n");
  partnerSource->pGetInt32Chunk( partnerSource, 100, 12356,buf );
  return 1;
}

Tool* createDownstreamTool(Arena* arena) {
  Tool* result= createBaseTool(arena);
  result->pInit= initDownstreamTool;
  result->pExecute= executeDownstreamTool;
  result->pAddSink( result, createTestSink(result) );
  result->typeName= "downstream";
  return result;
}

Tool* createStreamTool(Arena* arena) {
  Tool* result= createBaseTool(arena);
  result->typeName= "stream";
  return result;
}

