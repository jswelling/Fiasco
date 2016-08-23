/************************************************************
 *	zero_src_tool.c                                          *
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

static char rcsid[] = "$Id: zero_src_tool.c,v 1.2 2005/01/03 06:14:23 welling Exp $";

static long getUInt8Chunk( DataSource* self, long size, long long offset,
			   char* buf)
{
  int* result;
  int i;
  for (i=0; i<size; i++) buf[i]= 0;
  return size;
}

static long getInt16Chunk( DataSource* self, long size, long long offset,
			   short* buf)
{
  int* result;
  int i;
  for (i=0; i<size; i++) buf[i]= 0;
  return size;
}

static long getInt32Chunk( DataSource* self, long size, long long offset,
			   int* buf)
{
  int* result;
  int i;
  for (i=0; i<size; i++) buf[i]= 0;
  return size;
}

static long getInt64Chunk( DataSource* self, long size, long long offset,
			   long long* buf)
{
  int* result;
  int i;
  for (i=0; i<size; i++) buf[i]= 0;
  return size;
}

static long getFloat32Chunk( DataSource* self, long size, long long offset,
			     float* buf)
{
  int* result;
  int i;
  for (i=0; i<size; i++) buf[i]= 0.0;
  return size;
}

static long getFloat64Chunk( DataSource* self, long size, long long offset,
			     double* buf)
{
  int* result;
  int i;
  for (i=0; i<size; i++) buf[i]= 0.0;
  return size;
}

static DataSource* createZeroSource(Tool* owner)
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

static int initZeroSrcTool( Tool* t ) {
  fprintf(stderr,"initZeroSrcTool called!\n");
  return 1;
}

static int executeZeroSrcTool( Tool* t ) {
  fprintf(stderr,"executeZeroSrcTool called!\n");
  return 1;
}

Tool* createZeroSrcTool(Arena* arena, const char* dimstr, const char* extstr) {
  Tool* result= createBaseTool(arena);
  DataSource* mySrc= createZeroSource(result);
  const char* dimHere;
  const char* extHere;

  result->pInit= initZeroSrcTool;
  result->pExecute= executeZeroSrcTool;
  result->pAddSource(result, mySrc);
  result->typeName= "ZeroSrc";

  mySrc->pSetName(mySrc,"images");
  kvDefString(mySrc->attr,"datatype","float64");
  kvDefString(mySrc->attr,"dimensions",dimstr);
  dimHere= dimstr;
  extHere= extstr;
  while (*dimHere) {
    int ext;
    char keybuf[64];
    char buf[64];
    if (!sscanf(extHere,"%d",&ext)) {
      result->pDestroySelf(result);
      return NULL;
    }
#ifdef never
    fprintf(stderr,"Ext is %d for dim %c\n",ext,*dimHere);
#endif
    snprintf(keybuf,sizeof(keybuf),"extent.%c",*dimHere);
    snprintf(buf,sizeof(buf),"%d",ext);
    kvDefString(mySrc->attr,keybuf,buf);
    dimHere++;
    extHere= strchr(extHere,':');
    if ((*dimHere && !extHere) || (! *dimHere && extHere) ) {
      result->pDestroySelf(result);
      fprintf(stderr,
	      "%s: dimension and extent string don't match or extent syntax error\n",
	      result->typeName);
      return NULL;
    }
    extHere++; /* skip the ':' */
  }
  

  return result;
}


