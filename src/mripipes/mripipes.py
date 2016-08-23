# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _mripipes

def _swig_setattr(self,class_type,name,value):
    if (name == "this"):
        if isinstance(value, class_type):
            self.__dict__[name] = value.this
            if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
            del value.thisown
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    self.__dict__[name] = value

def _swig_getattr(self,class_type,name):
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0
del types


class DataSource(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, DataSource, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, DataSource, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<C DataSource instance at %s>" % (self.this,)
    __swig_setmethods__["pInit"] = _mripipes.DataSource_pInit_set
    __swig_getmethods__["pInit"] = _mripipes.DataSource_pInit_get
    if _newclass:pInit = property(_mripipes.DataSource_pInit_get, _mripipes.DataSource_pInit_set)
    __swig_setmethods__["pGetUInt8Chunk"] = _mripipes.DataSource_pGetUInt8Chunk_set
    __swig_getmethods__["pGetUInt8Chunk"] = _mripipes.DataSource_pGetUInt8Chunk_get
    if _newclass:pGetUInt8Chunk = property(_mripipes.DataSource_pGetUInt8Chunk_get, _mripipes.DataSource_pGetUInt8Chunk_set)
    __swig_setmethods__["pGetInt16Chunk"] = _mripipes.DataSource_pGetInt16Chunk_set
    __swig_getmethods__["pGetInt16Chunk"] = _mripipes.DataSource_pGetInt16Chunk_get
    if _newclass:pGetInt16Chunk = property(_mripipes.DataSource_pGetInt16Chunk_get, _mripipes.DataSource_pGetInt16Chunk_set)
    __swig_setmethods__["pGetInt32Chunk"] = _mripipes.DataSource_pGetInt32Chunk_set
    __swig_getmethods__["pGetInt32Chunk"] = _mripipes.DataSource_pGetInt32Chunk_get
    if _newclass:pGetInt32Chunk = property(_mripipes.DataSource_pGetInt32Chunk_get, _mripipes.DataSource_pGetInt32Chunk_set)
    __swig_setmethods__["pGetInt64Chunk"] = _mripipes.DataSource_pGetInt64Chunk_set
    __swig_getmethods__["pGetInt64Chunk"] = _mripipes.DataSource_pGetInt64Chunk_get
    if _newclass:pGetInt64Chunk = property(_mripipes.DataSource_pGetInt64Chunk_get, _mripipes.DataSource_pGetInt64Chunk_set)
    __swig_setmethods__["pGetFloat32Chunk"] = _mripipes.DataSource_pGetFloat32Chunk_set
    __swig_getmethods__["pGetFloat32Chunk"] = _mripipes.DataSource_pGetFloat32Chunk_get
    if _newclass:pGetFloat32Chunk = property(_mripipes.DataSource_pGetFloat32Chunk_get, _mripipes.DataSource_pGetFloat32Chunk_set)
    __swig_setmethods__["pGetFloat64Chunk"] = _mripipes.DataSource_pGetFloat64Chunk_set
    __swig_getmethods__["pGetFloat64Chunk"] = _mripipes.DataSource_pGetFloat64Chunk_get
    if _newclass:pGetFloat64Chunk = property(_mripipes.DataSource_pGetFloat64Chunk_get, _mripipes.DataSource_pGetFloat64Chunk_set)
    __swig_setmethods__["pConnect"] = _mripipes.DataSource_pConnect_set
    __swig_getmethods__["pConnect"] = _mripipes.DataSource_pConnect_get
    if _newclass:pConnect = property(_mripipes.DataSource_pConnect_get, _mripipes.DataSource_pConnect_set)
    __swig_setmethods__["pSetName"] = _mripipes.DataSource_pSetName_set
    __swig_getmethods__["pSetName"] = _mripipes.DataSource_pSetName_get
    if _newclass:pSetName = property(_mripipes.DataSource_pSetName_get, _mripipes.DataSource_pSetName_set)
    __swig_setmethods__["pDestroySelf"] = _mripipes.DataSource_pDestroySelf_set
    __swig_getmethods__["pDestroySelf"] = _mripipes.DataSource_pDestroySelf_get
    if _newclass:pDestroySelf = property(_mripipes.DataSource_pDestroySelf_get, _mripipes.DataSource_pDestroySelf_set)
    __swig_setmethods__["initialized"] = _mripipes.DataSource_initialized_set
    __swig_getmethods__["initialized"] = _mripipes.DataSource_initialized_get
    if _newclass:initialized = property(_mripipes.DataSource_initialized_get, _mripipes.DataSource_initialized_set)
    __swig_setmethods__["sink"] = _mripipes.DataSource_sink_set
    __swig_getmethods__["sink"] = _mripipes.DataSource_sink_get
    if _newclass:sink = property(_mripipes.DataSource_sink_get, _mripipes.DataSource_sink_set)
    __swig_setmethods__["attr"] = _mripipes.DataSource_attr_set
    __swig_getmethods__["attr"] = _mripipes.DataSource_attr_get
    if _newclass:attr = property(_mripipes.DataSource_attr_get, _mripipes.DataSource_attr_set)
    __swig_setmethods__["name"] = _mripipes.DataSource_name_set
    __swig_getmethods__["name"] = _mripipes.DataSource_name_get
    if _newclass:name = property(_mripipes.DataSource_name_get, _mripipes.DataSource_name_set)
    __swig_setmethods__["hook"] = _mripipes.DataSource_hook_set
    __swig_getmethods__["hook"] = _mripipes.DataSource_hook_get
    if _newclass:hook = property(_mripipes.DataSource_hook_get, _mripipes.DataSource_hook_set)
    def init(*args): return _mripipes.DataSource_init(*args)
    def setName(*args): return _mripipes.DataSource_setName(*args)
    def getName(*args): return _mripipes.DataSource_getName(*args)
    def __del__(self, destroy=_mripipes.delete_DataSource):
        try:
            if self.thisown: destroy(self)
        except: pass

class DataSourcePtr(DataSource):
    def __init__(self, this):
        _swig_setattr(self, DataSource, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, DataSource, 'thisown', 0)
        _swig_setattr(self, DataSource,self.__class__,DataSource)
_mripipes.DataSource_swigregister(DataSourcePtr)

class DataSink(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, DataSink, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, DataSink, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<C DataSink instance at %s>" % (self.this,)
    __swig_setmethods__["pInit"] = _mripipes.DataSink_pInit_set
    __swig_getmethods__["pInit"] = _mripipes.DataSink_pInit_get
    if _newclass:pInit = property(_mripipes.DataSink_pInit_get, _mripipes.DataSink_pInit_set)
    __swig_setmethods__["pConnect"] = _mripipes.DataSink_pConnect_set
    __swig_getmethods__["pConnect"] = _mripipes.DataSink_pConnect_get
    if _newclass:pConnect = property(_mripipes.DataSink_pConnect_get, _mripipes.DataSink_pConnect_set)
    __swig_setmethods__["pDestroySelf"] = _mripipes.DataSink_pDestroySelf_set
    __swig_getmethods__["pDestroySelf"] = _mripipes.DataSink_pDestroySelf_get
    if _newclass:pDestroySelf = property(_mripipes.DataSink_pDestroySelf_get, _mripipes.DataSink_pDestroySelf_set)
    __swig_setmethods__["initialized"] = _mripipes.DataSink_initialized_set
    __swig_getmethods__["initialized"] = _mripipes.DataSink_initialized_get
    if _newclass:initialized = property(_mripipes.DataSink_initialized_get, _mripipes.DataSink_initialized_set)
    __swig_setmethods__["source"] = _mripipes.DataSink_source_set
    __swig_getmethods__["source"] = _mripipes.DataSink_source_get
    if _newclass:source = property(_mripipes.DataSink_source_get, _mripipes.DataSink_source_set)
    __swig_setmethods__["name"] = _mripipes.DataSink_name_set
    __swig_getmethods__["name"] = _mripipes.DataSink_name_get
    if _newclass:name = property(_mripipes.DataSink_name_get, _mripipes.DataSink_name_set)
    __swig_setmethods__["hook"] = _mripipes.DataSink_hook_set
    __swig_getmethods__["hook"] = _mripipes.DataSink_hook_get
    if _newclass:hook = property(_mripipes.DataSink_hook_get, _mripipes.DataSink_hook_set)
    def init(*args): return _mripipes.DataSink_init(*args)
    def setName(*args): return _mripipes.DataSink_setName(*args)
    def getName(*args): return _mripipes.DataSink_getName(*args)
    def connect(*args): return _mripipes.DataSink_connect(*args)
    def __del__(self, destroy=_mripipes.delete_DataSink):
        try:
            if self.thisown: destroy(self)
        except: pass

class DataSinkPtr(DataSink):
    def __init__(self, this):
        _swig_setattr(self, DataSink, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, DataSink, 'thisown', 0)
        _swig_setattr(self, DataSink,self.__class__,DataSink)
_mripipes.DataSink_swigregister(DataSinkPtr)

class Tool(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Tool, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Tool, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<C Tool instance at %s>" % (self.this,)
    __swig_setmethods__["pInit"] = _mripipes.Tool_pInit_set
    __swig_getmethods__["pInit"] = _mripipes.Tool_pInit_get
    if _newclass:pInit = property(_mripipes.Tool_pInit_get, _mripipes.Tool_pInit_set)
    __swig_setmethods__["pExecute"] = _mripipes.Tool_pExecute_set
    __swig_getmethods__["pExecute"] = _mripipes.Tool_pExecute_get
    if _newclass:pExecute = property(_mripipes.Tool_pExecute_get, _mripipes.Tool_pExecute_set)
    __swig_setmethods__["pAddSource"] = _mripipes.Tool_pAddSource_set
    __swig_getmethods__["pAddSource"] = _mripipes.Tool_pAddSource_get
    if _newclass:pAddSource = property(_mripipes.Tool_pAddSource_get, _mripipes.Tool_pAddSource_set)
    __swig_setmethods__["pAddSink"] = _mripipes.Tool_pAddSink_set
    __swig_getmethods__["pAddSink"] = _mripipes.Tool_pAddSink_get
    if _newclass:pAddSink = property(_mripipes.Tool_pAddSink_get, _mripipes.Tool_pAddSink_set)
    __swig_setmethods__["pDestroySelf"] = _mripipes.Tool_pDestroySelf_set
    __swig_getmethods__["pDestroySelf"] = _mripipes.Tool_pDestroySelf_get
    if _newclass:pDestroySelf = property(_mripipes.Tool_pDestroySelf_get, _mripipes.Tool_pDestroySelf_set)
    __swig_setmethods__["pGetSourceByName"] = _mripipes.Tool_pGetSourceByName_set
    __swig_getmethods__["pGetSourceByName"] = _mripipes.Tool_pGetSourceByName_get
    if _newclass:pGetSourceByName = property(_mripipes.Tool_pGetSourceByName_get, _mripipes.Tool_pGetSourceByName_set)
    __swig_setmethods__["pGetSinkByName"] = _mripipes.Tool_pGetSinkByName_set
    __swig_getmethods__["pGetSinkByName"] = _mripipes.Tool_pGetSinkByName_get
    if _newclass:pGetSinkByName = property(_mripipes.Tool_pGetSinkByName_get, _mripipes.Tool_pGetSinkByName_set)
    __swig_setmethods__["initialized"] = _mripipes.Tool_initialized_set
    __swig_getmethods__["initialized"] = _mripipes.Tool_initialized_get
    if _newclass:initialized = property(_mripipes.Tool_initialized_get, _mripipes.Tool_initialized_set)
    __swig_setmethods__["verbose"] = _mripipes.Tool_verbose_set
    __swig_getmethods__["verbose"] = _mripipes.Tool_verbose_get
    if _newclass:verbose = property(_mripipes.Tool_verbose_get, _mripipes.Tool_verbose_set)
    __swig_setmethods__["debug"] = _mripipes.Tool_debug_set
    __swig_getmethods__["debug"] = _mripipes.Tool_debug_get
    if _newclass:debug = property(_mripipes.Tool_debug_get, _mripipes.Tool_debug_set)
    __swig_setmethods__["sourceArray"] = _mripipes.Tool_sourceArray_set
    __swig_getmethods__["sourceArray"] = _mripipes.Tool_sourceArray_get
    if _newclass:sourceArray = property(_mripipes.Tool_sourceArray_get, _mripipes.Tool_sourceArray_set)
    __swig_setmethods__["nSources"] = _mripipes.Tool_nSources_set
    __swig_getmethods__["nSources"] = _mripipes.Tool_nSources_get
    if _newclass:nSources = property(_mripipes.Tool_nSources_get, _mripipes.Tool_nSources_set)
    __swig_setmethods__["sourceArrayLength"] = _mripipes.Tool_sourceArrayLength_set
    __swig_getmethods__["sourceArrayLength"] = _mripipes.Tool_sourceArrayLength_get
    if _newclass:sourceArrayLength = property(_mripipes.Tool_sourceArrayLength_get, _mripipes.Tool_sourceArrayLength_set)
    __swig_setmethods__["sinkArray"] = _mripipes.Tool_sinkArray_set
    __swig_getmethods__["sinkArray"] = _mripipes.Tool_sinkArray_get
    if _newclass:sinkArray = property(_mripipes.Tool_sinkArray_get, _mripipes.Tool_sinkArray_set)
    __swig_setmethods__["nSinks"] = _mripipes.Tool_nSinks_set
    __swig_getmethods__["nSinks"] = _mripipes.Tool_nSinks_get
    if _newclass:nSinks = property(_mripipes.Tool_nSinks_get, _mripipes.Tool_nSinks_set)
    __swig_setmethods__["sinkArrayLength"] = _mripipes.Tool_sinkArrayLength_set
    __swig_getmethods__["sinkArrayLength"] = _mripipes.Tool_sinkArrayLength_get
    if _newclass:sinkArrayLength = property(_mripipes.Tool_sinkArrayLength_get, _mripipes.Tool_sinkArrayLength_set)
    __swig_setmethods__["typeName"] = _mripipes.Tool_typeName_set
    __swig_getmethods__["typeName"] = _mripipes.Tool_typeName_get
    if _newclass:typeName = property(_mripipes.Tool_typeName_get, _mripipes.Tool_typeName_set)
    __swig_setmethods__["owner"] = _mripipes.Tool_owner_set
    __swig_getmethods__["owner"] = _mripipes.Tool_owner_get
    if _newclass:owner = property(_mripipes.Tool_owner_get, _mripipes.Tool_owner_set)
    __swig_setmethods__["hook"] = _mripipes.Tool_hook_set
    __swig_getmethods__["hook"] = _mripipes.Tool_hook_get
    if _newclass:hook = property(_mripipes.Tool_hook_get, _mripipes.Tool_hook_set)
    def init(*args): return _mripipes.Tool_init(*args)
    def execute(*args): return _mripipes.Tool_execute(*args)
    def addSource(*args): return _mripipes.Tool_addSource(*args)
    def getNSources(*args): return _mripipes.Tool_getNSources(*args)
    def getSource(*args): return _mripipes.Tool_getSource(*args)
    def addSink(*args): return _mripipes.Tool_addSink(*args)
    def getNSinks(*args): return _mripipes.Tool_getNSinks(*args)
    def getSink(*args): return _mripipes.Tool_getSink(*args)
    def getSourceByName(*args): return _mripipes.Tool_getSourceByName(*args)
    def getSinkByName(*args): return _mripipes.Tool_getSinkByName(*args)
    def setDebug(*args): return _mripipes.Tool_setDebug(*args)
    def setVerbose(*args): return _mripipes.Tool_setVerbose(*args)
    def __del__(self, destroy=_mripipes.delete_Tool):
        try:
            if self.thisown: destroy(self)
        except: pass

class ToolPtr(Tool):
    def __init__(self, this):
        _swig_setattr(self, Tool, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Tool, 'thisown', 0)
        _swig_setattr(self, Tool,self.__class__,Tool)
_mripipes.Tool_swigregister(ToolPtr)

class Arena(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Arena, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Arena, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<C Arena instance at %s>" % (self.this,)
    __swig_setmethods__["pDestroySelf"] = _mripipes.Arena_pDestroySelf_set
    __swig_getmethods__["pDestroySelf"] = _mripipes.Arena_pDestroySelf_get
    if _newclass:pDestroySelf = property(_mripipes.Arena_pDestroySelf_get, _mripipes.Arena_pDestroySelf_set)
    __swig_setmethods__["pAddTool"] = _mripipes.Arena_pAddTool_set
    __swig_getmethods__["pAddTool"] = _mripipes.Arena_pAddTool_get
    if _newclass:pAddTool = property(_mripipes.Arena_pAddTool_get, _mripipes.Arena_pAddTool_set)
    __swig_setmethods__["pInit"] = _mripipes.Arena_pInit_set
    __swig_getmethods__["pInit"] = _mripipes.Arena_pInit_get
    if _newclass:pInit = property(_mripipes.Arena_pInit_get, _mripipes.Arena_pInit_set)
    __swig_setmethods__["pExecute"] = _mripipes.Arena_pExecute_set
    __swig_getmethods__["pExecute"] = _mripipes.Arena_pExecute_get
    if _newclass:pExecute = property(_mripipes.Arena_pExecute_get, _mripipes.Arena_pExecute_set)
    __swig_setmethods__["tools"] = _mripipes.Arena_tools_set
    __swig_getmethods__["tools"] = _mripipes.Arena_tools_get
    if _newclass:tools = property(_mripipes.Arena_tools_get, _mripipes.Arena_tools_set)
    __swig_setmethods__["drain"] = _mripipes.Arena_drain_set
    __swig_getmethods__["drain"] = _mripipes.Arena_drain_get
    if _newclass:drain = property(_mripipes.Arena_drain_get, _mripipes.Arena_drain_set)
    __swig_setmethods__["verbose"] = _mripipes.Arena_verbose_set
    __swig_getmethods__["verbose"] = _mripipes.Arena_verbose_get
    if _newclass:verbose = property(_mripipes.Arena_verbose_get, _mripipes.Arena_verbose_set)
    __swig_setmethods__["debug"] = _mripipes.Arena_debug_set
    __swig_getmethods__["debug"] = _mripipes.Arena_debug_get
    if _newclass:debug = property(_mripipes.Arena_debug_get, _mripipes.Arena_debug_set)
    def init(*args): return _mripipes.Arena_init(*args)
    def execute(*args): return _mripipes.Arena_execute(*args)
    def addTool(*args): return _mripipes.Arena_addTool(*args)
    def __del__(self, destroy=_mripipes.delete_Arena):
        try:
            if self.thisown: destroy(self)
        except: pass
    def setDebug(*args): return _mripipes.Arena_setDebug(*args)
    def setVerbose(*args): return _mripipes.Arena_setVerbose(*args)

class ArenaPtr(Arena):
    def __init__(self, this):
        _swig_setattr(self, Arena, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Arena, 'thisown', 0)
        _swig_setattr(self, Arena,self.__class__,Arena)
_mripipes.Arena_swigregister(ArenaPtr)


getSourceDims = _mripipes.getSourceDims

getSourceDimExtent = _mripipes.getSourceDimExtent

baseToolDestroySelf = _mripipes.baseToolDestroySelf

baseToolInit = _mripipes.baseToolInit

baseToolExecute = _mripipes.baseToolExecute

createArena = _mripipes.createArena

createBaseSource = _mripipes.createBaseSource

createBaseSink = _mripipes.createBaseSink

createTestSource = _mripipes.createTestSource

createTestSink = _mripipes.createTestSink

createBaseTool = _mripipes.createBaseTool

createUpstreamTool = _mripipes.createUpstreamTool

createStreamTool = _mripipes.createStreamTool

createDownstreamTool = _mripipes.createDownstreamTool

createMRIFileInputTool = _mripipes.createMRIFileInputTool

createMRIFileOutputTool = _mripipes.createMRIFileOutputTool

createDevnullTool = _mripipes.createDevnullTool

createPassthruTool = _mripipes.createPassthruTool

createMatmultTool = _mripipes.createMatmultTool

createComplexMatmultTool = _mripipes.createComplexMatmultTool

createRpnMathTool = _mripipes.createRpnMathTool

createComplexRpnMathTool = _mripipes.createComplexRpnMathTool

createZeroSrcTool = _mripipes.createZeroSrcTool

createSubsetTool = _mripipes.createSubsetTool

createPadTool = _mripipes.createPadTool

createBlockMapTool = _mripipes.createBlockMapTool

createFunc2UnblkTool = _mripipes.createFunc2UnblkTool

createSpecialTool = _mripipes.createSpecialTool

