%module mripipes
%{
#include "fmri.h"
#include "fexceptions.h"
#include "mripipes.h"
%}

/* $Id: mripipes.i,v 1.15 2007/06/21 23:25:15 welling Exp $ */

/*******************************************************
 * Notes-
 * -BlockMapTool's python wrapper does a Py_INCREF on a
 *  tuple which is passed in in the constructor arg list.
 *  This should get DECREF'd when the tool is deleted, but
 *  I can't figure out how to cause SWIG to make this happen.
 *******************************************************/

%include "typemaps.i"

%nodefault;

typedef struct data_source_struct {
  int (*pInit)(DataSource* self);
  long (*pGetUInt8Chunk)(DataSource* self,
                         long size, long long offset, char* buf );
  long (*pGetInt16Chunk)(DataSource* self,
			 long size, long long offset, short* buf );
  long (*pGetInt32Chunk)(DataSource* self,
			 long size, long long offset, int* buf );
  long (*pGetInt64Chunk)(DataSource* self,
			 long size, long long offset, long long* buf );
  long (*pGetFloat32Chunk)(DataSource* self,
			   long size, long long offset, float* buf );
  long (*pGetFloat64Chunk)(DataSource* self,
			   long size, long long offset, double* buf );
  int (*pConnect)(DataSource* self, DataSink* sink);
  void (*pSetName)(DataSource* self, const char* name);
  void (*pDestroySelf)(DataSource* self);
  int initialized;
  DataSink* sink;
  KVHash* attr;
  char* name;
  void* hook;
} DataSource;

typedef struct data_sink_struct {
  int (*pInit)(DataSink* self);
  int (*pConnect)(DataSink* self, DataSource* source_in);
  void (*pDestroySelf)(DataSink* self);
  int initialized;
  DataSource* source;
  char* name;
  void* hook;
} DataSink;

typedef struct tool_struct {
  int (*pInit)(Tool* self);
  int (*pExecute)(Tool* self);
  void (*pAddSource)(struct tool_struct* self, DataSource* src);
  void (*pAddSink)(struct tool_struct* self, DataSink* sink);
  void (*pDestroySelf)(Tool* self);
  DataSource* (*pGetSourceByName)(struct tool_struct* self, const char* nm);
  DataSink* (*pGetSinkByName)(struct tool_struct* self, const char* nm);
  int initialized;
  int verbose;
  int debug;
  DataSource** sourceArray;
  int nSources;
  int sourceArrayLength;
  DataSink** sinkArray;
  int nSinks;
  int sinkArrayLength;
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

%makedefault;

%{
/* This function matches the prototype of the normal C callback
   function for BlockMapTool. However, we use the clientdata pointer
   for holding a reference to a Python callable object. */
static int PythonBLOCKMAPCBFUNC(long* size, long long* offset,
				void *clientdata)
{
    PyObject *initFunc, *cbFunc, *data, *arglist;
    PyObject *result;
    long     retCode = 1;
    
    // Break data into components
    initFunc= PyTuple_GetItem(clientdata,0);
    cbFunc= PyTuple_GetItem(clientdata,1);
    data= PyTuple_GetItem(clientdata,2);
    
    arglist = Py_BuildValue("lLO",*size,*offset,data);  // Build argument list
    result = PyEval_CallObject(cbFunc,arglist);         // Call Python
    Py_DECREF(arglist);                                 // Trash arglist
    if (result) {
      if (!PyArg_ParseTuple(result,"lL:PythonBLOCKMAPCBFUNC",size,offset)) {
	retCode= 0;
      }
    }
    else retCode= 0;
    Py_XDECREF(result);
    return retCode;
}

/* This function matches the prototype of the normal C Init callback
   function for BlockMapTool. However, we use the clientdata pointer
   for holding a reference to a Python callable object. */
 static int PythonBLOCKMAPINITFUNC(char dim, const char* dimstr,
				   long long fastBlkSize,
				   long upstream_extent,
				   long extent1, long extent2,
				   long long slowBlkSize,
				   void *clientdata)
{
    PyObject *initFunc, *cbFunc, *data, *arglist;
    PyObject *result;
    int retCode= 0;
    
    // Break data into components
    initFunc= PyTuple_GetItem(clientdata,0);
    cbFunc= PyTuple_GetItem(clientdata,1);
    data= PyTuple_GetItem(clientdata,2);
    
    arglist = Py_BuildValue("csLlllLO",
			    dim, dimstr, fastBlkSize, 
			    upstream_extent, extent1, extent2,
			    slowBlkSize, data);   // Build argument list
    result = PyEval_CallObject(initFunc,arglist); // Call Python
    Py_DECREF(arglist);                           // Trash arglist
    if (result) {
      retCode= PyObject_IsTrue(result);
    }
    else retCode= 0;
    Py_XDECREF(result);
    return retCode;
}
%}

/* The expected input here is a Python Tuple containing two functions
 * and a pointer for data.
 */
%typemap(in) (BLOCKMAPINITFUNC initFunc, BLOCKMAPCBFUNC cbFunc, void* cbData) {
  if (!PyTuple_Check($input)) {
    PyErr_SetString(PyExc_TypeError, "Expected a tuple!");
    return NULL;
  }
  if (PyTuple_Size($input)!=3) {
    PyErr_SetString(PyExc_TypeError, "Tuple size is not 3!");
    return NULL;
  }
  if (!PyCallable_Check(PyTuple_GetItem($input,0))) {
    PyErr_SetString(PyExc_TypeError, "First tuple element is not a function!");
    return NULL;
  }
  if (!PyCallable_Check(PyTuple_GetItem($input,1))) {
    PyErr_SetString(PyExc_TypeError, 
		    "Second tuple element is not a function!");
    return NULL;
  }

  Py_INCREF($input);
  $1= PythonBLOCKMAPINITFUNC; $2= PythonBLOCKMAPCBFUNC; $3= $input;
}

%extend DataSource {
  int init() { return self->pInit(self); }
  void setName( const char* name ) { self->pSetName(self,name); }
  const char* getName() { return self->name; }
  ~DataSource() { self->pDestroySelf(self); }
  /*
  char* getUInt8Chunk( long size, long long offset ) {
    self->pGetUInt8Chunk( self, size, offset );
  }
  char* getInt16Chunk( long size, long long offset ) {
    self->pGetInt16Chunk( self, size, offset );
  }
  char* getInt32Chunk( long size, long long offset ) {
    self->pGetInt32Chunk( self, size, offset );
  }
  char* getInt64Chunk( long size, long long offset ) {
    self->pGetInt64Chunk( self, size, offset );
  }
  char* getFloat32Chunk( long size, long long offset ) {
    self->pGetFloat32Chunk( self, size, offset );
  }
  char* getFloat64Chunk( long size, long long offset ) {
    self->pGetFloat64Chunk( self, size, offset );
  }
  */
};

%extend DataSink {
  int init() { return self->pInit(self); }
  void setName( const char* name ) { self->pSetName(self,name); }
  const char* getName() { return self->name; }
  int connect( DataSource* source ) { return self->pConnect(self,source); }
  ~DataSink() { self->pDestroySelf(self); }
}

%extend Tool {
  int init() { return self->pInit(self); }
  int execute() { return self->pExecute(self); }
  void addSource(DataSource* src) { self->pAddSource(self,src); }
  int getNSources() { return self->nSources; }
  DataSource* getSource( int i ) {
    if (i>=0 && i<self->nSources) return self->sourceArray[i];
    else return NULL;
  }
  void addSink(DataSink* sink) { self->pAddSink(self,sink); }
  int getNSinks() { return self->nSinks; }
  DataSink* getSink( int i ) {
    if (i>=0 && i<self->nSinks) return self->sinkArray[i];
    else return NULL;
  }
  DataSource* getSourceByName( const char* name ) { 
    return self->pGetSourceByName(self,name); 
  }
  DataSink* getSinkByName( const char* name ) { 
    return self->pGetSinkByName(self,name); 
  }
  void setDebug(int i=1) {
    self->debug= i;
  }
  void setVerbose(int i=1) {
    self->verbose= i;
  }
  ~Tool() { self->pDestroySelf(self); }
}

/* Below this line, anything that raises a FEX exception
 * will cause a Python exception! */
%exception {
  ExceptionContext eCtx;
  Exception* evar= NULL;
  __fex_windExceptionContext(&eCtx);
  if (!sigsetjmp(eCtx.env,1)) {
    $action 
    __fex_unwindExceptionContext(); 
  }
  else {
    if (!(evar=__fex_getCurrentException())) {
      PyErr_SetString(PyExc_RuntimeError,
		      "Error translating FEX exception from C code");
    }
    else {
      PyObject* exception= NULL;
      __fex_unwindExceptionContext();
      switch (fex_getExceptionType(evar)) {
      case EXCEPTION_IO: 
	exception= PyExc_IOError;
	break;
      case EXCEPTION_MEM:
	exception= PyExc_MemoryError;
	break;
      default:
	exception= PyExc_RuntimeError;
      }
      PyErr_SetString(exception,fex_getExceptionString(evar));
      __fex_destroyException(evar);
    }
    return NULL;
  }
}

%extend Arena {
  int init() { return self->pInit(self); }
  int execute() { return self->pExecute(self); }
  void addTool(Tool* t) { self->pAddTool(self,t); }
  ~Arena() { self->pDestroySelf(self); }
  void setDebug(int i=1) {
    self->debug= i;
  }
  void setVerbose(int i=1) {
    self->verbose= i;
  }
}

const char* getSourceDims( DataSource* source );
int getSourceDimExtent( DataSource* source, char dim );

void baseToolDestroySelf(Tool* self);
int baseToolInit(Tool* self);
int baseToolExecute(Tool* self);

Arena* createArena();

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
Tool* createRpnMathTool(Arena* arena, const char* script);
Tool* createComplexRpnMathTool(Arena* arena, const char* script);
Tool* createZeroSrcTool(Arena* arena, const char* dimstr, 
			const char* extentstr);
Tool* createSubsetTool(Arena* arena, const char* dim, const int extent,
		       const int shift );
Tool* createPadTool(Arena* arena, const char* dim, const int extent,
                    const int shift, const double fillValue=0.0 );
Tool* createBlockMapTool(Arena* arena, const char* dim, const char* newdim,
			 long extent1, long extent2,
			 BLOCKMAPINITFUNC initFunc,
			 BLOCKMAPCBFUNC cbFunc, void* cbData );
Tool* createFunc2UnblkTool(Arena* arena, FUNC2UNBLKCBFUNC func,
			   long nOutput, void* cbData);
Tool* createSpecialTool(Arena* arena);
