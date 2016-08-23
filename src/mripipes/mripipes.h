/************************************************************
 *      mripipes.h                                          *
 *                                                          *
 *      Copyright (c) 2004 Pittsburgh Supercomputing Center *
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

#include <kvhash.h>
#include <slist.h>

typedef int (*BLOCKMAPINITFUNC)(char dim, const char* dimstr,
				long long fastBlkSize, 
				long upstream_extent,
				long extent, long extent2,
				long long slowBlkSize, 
				void* hook);
typedef int (*BLOCKMAPCBFUNC)(long* size, long long* offset, void* hook);

typedef int (*FUNC2UNBLKCBFUNC)(const double* left, const double* right,
				long nIn,
				double* out, long nOut,
				void* hook);

struct arena_struct;
struct tool_struct;
struct data_sink_struct;

typedef struct data_source_struct {
  int (*pInit)(struct data_source_struct* self);
  long (*pGetUInt8Chunk)(struct data_source_struct* self, 
			 long size, long long offset, char* buf );
  long (*pGetInt16Chunk)(struct data_source_struct* self, 
			  long size, long long offset, short* buf );
  long (*pGetInt32Chunk)(struct data_source_struct* self, 
			long size, long long offset, int* buf );
  long (*pGetInt64Chunk)(struct data_source_struct* self, 
			 long size, long long offset, long long* buf );
  long (*pGetFloat32Chunk)(struct data_source_struct* self, 
			   long size, long long offset, float* buf );
  long (*pGetFloat64Chunk)(struct data_source_struct* self, 
			   long size, long long offset, double* buf );
  int (*pConnect)(struct data_source_struct* self, 
		 struct data_sink_struct* sink);
  void (*pSetName)(struct data_source_struct* self,
		   const char* name);
  void (*pDestroySelf)(struct data_source_struct* self);
  const char* (*pGetName)(struct data_source_struct* self);
  int initialized;
  struct data_sink_struct* sink;
  struct tool_struct* owner;
  KVHash* attr;
  char* name;
  void* hook;
} DataSource;

typedef struct data_sink_struct {
  int (*pInit)(struct data_sink_struct* self);
  int (*pConnect)(struct data_sink_struct* self, DataSource* source_in);
  void (*pSetName)(struct data_sink_struct* self,
		   const char* name);
  void (*pDestroySelf)(struct data_sink_struct* self);
  const char* (*pGetName)(struct data_sink_struct* self);
  int initialized;
  DataSource* source;
  struct tool_struct* owner;
  char* name;
  void* hook;
} DataSink;

typedef struct tool_struct {
  int (*pInit)(struct tool_struct* self);
  int (*pExecute)(struct tool_struct* self);
  void (*pAddSource)(struct tool_struct* self, DataSource* src);
  void (*pAddSink)(struct tool_struct* self, DataSink* sink);
  void (*pDestroySelf)(struct tool_struct* self);
  DataSource* (*pGetSourceByName)(struct tool_struct* self, const char* nm);
  DataSink* (*pGetSinkByName)(struct tool_struct* self, const char* nm);
  int initialized;
  int verbose;
  int debug;
  int nSources;
  int sourceArrayLength;
  DataSource** sourceArray;
  int nSinks;
  int sinkArrayLength;
  DataSink** sinkArray;
  const char* typeName;
  struct arena_struct* owner;
  void* hook;
} Tool;

typedef struct arena_struct {
  void (*pDestroySelf)(struct arena_struct* self);
  void (*pAddTool)(struct arena_struct* self, Tool* t);
  int (*pInit)(struct arena_struct* self);
  int (*pExecute)(struct arena_struct* self);
  SList* tools;
  Tool* drain;
  int verbose;
  int debug;
} Arena;

extern void pipeAbort(char* fmt, ...);

const char* getSourceDims( DataSource* source );
void setSourceDims( DataSource* source, const char* dimstr );
int getSourceDimExtent( DataSource* source, const char dim );
void setSourceDimExtent( DataSource* source, const char dim, const int extent);
int getSourceDataType(DataSource* source);
void forceGetAllFloat64(DataSource* source, long size, long long offset,
			double* buf);
void calcSourceBlockSizes(DataSource* source, const char* dimstr,
			  const char selected_dim,
			  long* fast_blocksize_out, 
			  long long* slow_blocksize_out );

int baseToolInit(Tool* self);
int baseToolExecute(Tool* self);
void baseToolDestroySelf(Tool* self);

Arena* createArena(void);

DataSource* createBaseSource(Tool* owner);
DataSink* createBaseSink(Tool* owner);

DataSource* createTestSource(Tool* owner);
DataSink* createTestSink(Tool* owner);

Tool* createBaseTool(Arena* arena);
Tool* createUpstreamTool(Arena* arena);
Tool* createStreamTool(Arena* arena);
Tool* createDownstreamTool(Arena* arena);
Tool* createMRIFileInputTool(Arena* arena, const char* fname);
Tool* createMRIFileOutputTool(Arena* arena, const char* fname);
Tool* createDevnullTool(Arena* arena);
Tool* createPassthruTool(Arena* arena);
Tool* createMatmultTool(Arena* arena);
Tool* createComplexMatmultTool(Arena* arena);
Tool* createRpnMathTool(Arena* arena, const char* script_in);
Tool* createComplexRpnMathTool(Arena* arena, const char* script_in);
Tool* createZeroSrcTool(Arena* arena, const char* dimstr, 
			const char* extentstr);
Tool* createSubsetTool(Arena* arena, const char* dim, const int extent,
		       const int shift );
Tool* createPadTool(Arena* arena, const char* dim, const int extent,
		    const int shift, const double fillValue );
Tool* createBlockMapTool(Arena* arena, const char* dim, const char* newdim,
			 long extent1, long extent2,
			 BLOCKMAPINITFUNC initFunc,
			 BLOCKMAPCBFUNC cbFunc, void* cbData );
Tool* createFunc2UnblkTool(Arena* arena, FUNC2UNBLKCBFUNC func, long nOutputs,
			   void* cbData);
Tool* createSpecialTool(Arena* arena);
