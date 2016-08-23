# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _fiasco_numpy

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


GLM_COMPLEX = _fiasco_numpy.GLM_COMPLEX
GLM_RESIDUALS = _fiasco_numpy.GLM_RESIDUALS
GLM_VARIANCES = _fiasco_numpy.GLM_VARIANCES
GLM_COVARIANCES = _fiasco_numpy.GLM_COVARIANCES
GLM_SSQR = _fiasco_numpy.GLM_SSQR
GLM_ORTHO = _fiasco_numpy.GLM_ORTHO
GLM_DEBUG = _fiasco_numpy.GLM_DEBUG
GLM_DEVIANCE = _fiasco_numpy.GLM_DEVIANCE
GLM_TYPE = _fiasco_numpy.GLM_TYPE
GLM_TYPE_BASE = _fiasco_numpy.GLM_TYPE_BASE
GLM_TYPE_LLSQ = _fiasco_numpy.GLM_TYPE_LLSQ
GLM_TYPE_LOGISTIC = _fiasco_numpy.GLM_TYPE_LOGISTIC
GLM_TYPE_POISSON = _fiasco_numpy.GLM_TYPE_POISSON
class Regressor(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Regressor, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Regressor, name)
    def __repr__(self):
        return "<C Regressor instance at %s>" % (self.this,)
    def __init__(self, *args):
        _swig_setattr(self, Regressor, 'this', _fiasco_numpy.new_Regressor(*args))
        _swig_setattr(self, Regressor, 'thisown', 1)
    def __del__(self, destroy=_fiasco_numpy.delete_Regressor):
        try:
            if self.thisown: destroy(self)
        except: pass
    def get(*args): return _fiasco_numpy.Regressor_get(*args)
    def is_settable(*args): return _fiasco_numpy.Regressor_is_settable(*args)
    def set(*args): return _fiasco_numpy.Regressor_set(*args)
    def n_params(*args): return _fiasco_numpy.Regressor_n_params(*args)
    def fit(*args): return _fiasco_numpy.Regressor_fit(*args)
    def context_valid(*args): return _fiasco_numpy.Regressor_context_valid(*args)

class RegressorPtr(Regressor):
    def __init__(self, this):
        _swig_setattr(self, Regressor, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Regressor, 'thisown', 0)
        _swig_setattr(self, Regressor,self.__class__,Regressor)
_fiasco_numpy.Regressor_swigregister(RegressorPtr)


create_base_regressor = _fiasco_numpy.create_base_regressor

create_llsq_regressor = _fiasco_numpy.create_llsq_regressor

create_logistic_regressor = _fiasco_numpy.create_logistic_regressor

create_poisson_regressor = _fiasco_numpy.create_poisson_regressor
class ScalarFunction(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, ScalarFunction, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, ScalarFunction, name)
    def __repr__(self):
        return "<C ScalarFunction instance at %s>" % (self.this,)
    def __del__(self, destroy=_fiasco_numpy.delete_ScalarFunction):
        try:
            if self.thisown: destroy(self)
        except: pass
    def value(*args): return _fiasco_numpy.ScalarFunction_value(*args)
    def reset(*args): return _fiasco_numpy.ScalarFunction_reset(*args)
    def __init__(self, *args):
        _swig_setattr(self, ScalarFunction, 'this', _fiasco_numpy.new_ScalarFunction(*args))
        _swig_setattr(self, ScalarFunction, 'thisown', 1)

class ScalarFunctionPtr(ScalarFunction):
    def __init__(self, this):
        _swig_setattr(self, ScalarFunction, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ScalarFunction, 'thisown', 0)
        _swig_setattr(self, ScalarFunction,self.__class__,ScalarFunction)
_fiasco_numpy.ScalarFunction_swigregister(ScalarFunctionPtr)


buildSimpleScalarFunction = _fiasco_numpy.buildSimpleScalarFunction
class Optimizer(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Optimizer, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Optimizer, name)
    def __del__(self, destroy=_fiasco_numpy.delete_Optimizer):
        try:
            if self.thisown: destroy(self)
        except: pass
    def __str__(*args): return _fiasco_numpy.Optimizer___str__(*args)
    def __repr__(*args): return _fiasco_numpy.Optimizer___repr__(*args)
    def getDebugLevel(*args): return _fiasco_numpy.Optimizer_getDebugLevel(*args)
    def setDebugLevel(*args): return _fiasco_numpy.Optimizer_setDebugLevel(*args)
    def go(*args): return _fiasco_numpy.Optimizer_go(*args)
    def __init__(self, *args):
        _swig_setattr(self, Optimizer, 'this', _fiasco_numpy.new_Optimizer(*args))
        _swig_setattr(self, Optimizer, 'thisown', 1)

class OptimizerPtr(Optimizer):
    def __init__(self, this):
        _swig_setattr(self, Optimizer, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Optimizer, 'thisown', 0)
        _swig_setattr(self, Optimizer,self.__class__,Optimizer)
_fiasco_numpy.Optimizer_swigregister(OptimizerPtr)


optimizerFromStringRep = _fiasco_numpy.optimizerFromStringRep

createBaseOptimizer = _fiasco_numpy.createBaseOptimizer

createNoneOptimizer = _fiasco_numpy.createNoneOptimizer

createPraxisOptimizer = _fiasco_numpy.createPraxisOptimizer

createNelminOptimizer = _fiasco_numpy.createNelminOptimizer

createNelminTOptimizer = _fiasco_numpy.createNelminTOptimizer
INTRP_CLOSEST = _fiasco_numpy.INTRP_CLOSEST
INTRP_LINEAR = _fiasco_numpy.INTRP_LINEAR
INTRP_CATMULLROM = _fiasco_numpy.INTRP_CATMULLROM
INTRP_BEZIER = _fiasco_numpy.INTRP_BEZIER
INTRP_BSPLINE = _fiasco_numpy.INTRP_BSPLINE
INTRP_UNKNOWN = _fiasco_numpy.INTRP_UNKNOWN
INTRP_OPT_DEBUG = _fiasco_numpy.INTRP_OPT_DEBUG
INTRP_OPT_FASTBLK = _fiasco_numpy.INTRP_OPT_FASTBLK
INTRP_OPT_EXTENT = _fiasco_numpy.INTRP_OPT_EXTENT
INTRP_OPT_TENSION = _fiasco_numpy.INTRP_OPT_TENSION
class Interpolator(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Interpolator, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Interpolator, name)
    def __repr__(self):
        return "<C Interpolator instance at %s>" % (self.this,)
    def __del__(self, destroy=_fiasco_numpy.delete_Interpolator):
        try:
            if self.thisown: destroy(self)
        except: pass
    def getTypeName(*args): return _fiasco_numpy.Interpolator_getTypeName(*args)
    def prep(*args): return _fiasco_numpy.Interpolator_prep(*args)
    def calc(*args): return _fiasco_numpy.Interpolator_calc(*args)
    def dumpSelf(*args): return _fiasco_numpy.Interpolator_dumpSelf(*args)
    def setInt(*args): return _fiasco_numpy.Interpolator_setInt(*args)
    def setDouble(*args): return _fiasco_numpy.Interpolator_setDouble(*args)
    def getInt(*args): return _fiasco_numpy.Interpolator_getInt(*args)
    def getDouble(*args): return _fiasco_numpy.Interpolator_getDouble(*args)
    def __init__(self, *args):
        _swig_setattr(self, Interpolator, 'this', _fiasco_numpy.new_Interpolator(*args))
        _swig_setattr(self, Interpolator, 'thisown', 1)

class InterpolatorPtr(Interpolator):
    def __init__(self, this):
        _swig_setattr(self, Interpolator, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Interpolator, 'thisown', 0)
        _swig_setattr(self, Interpolator,self.__class__,Interpolator)
_fiasco_numpy.Interpolator_swigregister(InterpolatorPtr)


intrp_typeFromName = _fiasco_numpy.intrp_typeFromName

intrp_nameFromType = _fiasco_numpy.intrp_nameFromType

intrp_warpClearCounts = _fiasco_numpy.intrp_warpClearCounts

intrp_warpGetCounts = _fiasco_numpy.intrp_warpGetCounts

warpApply = _fiasco_numpy.warpApply

createClosestInterpolator1D = _fiasco_numpy.createClosestInterpolator1D

createClosestInterpolator2D = _fiasco_numpy.createClosestInterpolator2D

createClosestInterpolator3D = _fiasco_numpy.createClosestInterpolator3D

createLinearInterpolator1D = _fiasco_numpy.createLinearInterpolator1D

createLinearInterpolator2D = _fiasco_numpy.createLinearInterpolator2D

createLinearInterpolator3D = _fiasco_numpy.createLinearInterpolator3D

createSplineInterpolator1D = _fiasco_numpy.createSplineInterpolator1D

createSplineInterpolator2D = _fiasco_numpy.createSplineInterpolator2D

createSplineInterpolator3D = _fiasco_numpy.createSplineInterpolator3D

createInterpolator1DByType = _fiasco_numpy.createInterpolator1DByType

createInterpolator2DByType = _fiasco_numpy.createInterpolator2DByType

createInterpolator3DByType = _fiasco_numpy.createInterpolator3DByType
class KalmanProcess(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, KalmanProcess, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, KalmanProcess, name)
    def __repr__(self):
        return "<C KalmanProcess instance at %s>" % (self.this,)
    def __del__(self, destroy=_fiasco_numpy.delete_KalmanProcess):
        try:
            if self.thisown: destroy(self)
        except: pass
    def dumpSelf(*args): return _fiasco_numpy.KalmanProcess_dumpSelf(*args)
    def setDebug(*args): return _fiasco_numpy.KalmanProcess_setDebug(*args)
    def getDebug(*args): return _fiasco_numpy.KalmanProcess_getDebug(*args)
    def setA(*args): return _fiasco_numpy.KalmanProcess_setA(*args)
    def getA(*args): return _fiasco_numpy.KalmanProcess_getA(*args)
    def setH(*args): return _fiasco_numpy.KalmanProcess_setH(*args)
    def getH(*args): return _fiasco_numpy.KalmanProcess_getH(*args)
    def setQ(*args): return _fiasco_numpy.KalmanProcess_setQ(*args)
    def getQ(*args): return _fiasco_numpy.KalmanProcess_getQ(*args)
    def setR(*args): return _fiasco_numpy.KalmanProcess_setR(*args)
    def getR(*args): return _fiasco_numpy.KalmanProcess_getR(*args)
    def getL(*args): return _fiasco_numpy.KalmanProcess_getL(*args)
    def getM(*args): return _fiasco_numpy.KalmanProcess_getM(*args)
    def apply(*args): return _fiasco_numpy.KalmanProcess_apply(*args)
    def __init__(self, *args):
        _swig_setattr(self, KalmanProcess, 'this', _fiasco_numpy.new_KalmanProcess(*args))
        _swig_setattr(self, KalmanProcess, 'thisown', 1)

class KalmanProcessPtr(KalmanProcess):
    def __init__(self, this):
        _swig_setattr(self, KalmanProcess, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, KalmanProcess, 'thisown', 0)
        _swig_setattr(self, KalmanProcess,self.__class__,KalmanProcess)
_fiasco_numpy.KalmanProcess_swigregister(KalmanProcessPtr)

class KalmanState(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, KalmanState, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, KalmanState, name)
    def __repr__(self):
        return "<C KalmanState instance at %s>" % (self.this,)
    def __del__(self, destroy=_fiasco_numpy.delete_KalmanState):
        try:
            if self.thisown: destroy(self)
        except: pass
    def dumpSelf(*args): return _fiasco_numpy.KalmanState_dumpSelf(*args)
    def setX(*args): return _fiasco_numpy.KalmanState_setX(*args)
    def setP(*args): return _fiasco_numpy.KalmanState_setP(*args)
    def getM(*args): return _fiasco_numpy.KalmanState_getM(*args)
    def __init__(self, *args):
        _swig_setattr(self, KalmanState, 'this', _fiasco_numpy.new_KalmanState(*args))
        _swig_setattr(self, KalmanState, 'thisown', 1)

class KalmanStatePtr(KalmanState):
    def __init__(self, this):
        _swig_setattr(self, KalmanState, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, KalmanState, 'thisown', 0)
        _swig_setattr(self, KalmanState,self.__class__,KalmanState)
_fiasco_numpy.KalmanState_swigregister(KalmanStatePtr)

class KalmanFilter(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, KalmanFilter, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, KalmanFilter, name)
    def __repr__(self):
        return "<C KalmanFilter instance at %s>" % (self.this,)
    def __del__(self, destroy=_fiasco_numpy.delete_KalmanFilter):
        try:
            if self.thisown: destroy(self)
        except: pass
    def dumpSelf(*args): return _fiasco_numpy.KalmanFilter_dumpSelf(*args)
    def getL(*args): return _fiasco_numpy.KalmanFilter_getL(*args)
    def getM(*args): return _fiasco_numpy.KalmanFilter_getM(*args)
    def getState(*args): return _fiasco_numpy.KalmanFilter_getState(*args)
    def run(*args): return _fiasco_numpy.KalmanFilter_run(*args)
    def __init__(self, *args):
        _swig_setattr(self, KalmanFilter, 'this', _fiasco_numpy.new_KalmanFilter(*args))
        _swig_setattr(self, KalmanFilter, 'thisown', 1)

class KalmanFilterPtr(KalmanFilter):
    def __init__(self, this):
        _swig_setattr(self, KalmanFilter, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, KalmanFilter, 'thisown', 0)
        _swig_setattr(self, KalmanFilter,self.__class__,KalmanFilter)
_fiasco_numpy.KalmanFilter_swigregister(KalmanFilterPtr)


createKalmanProcess = _fiasco_numpy.createKalmanProcess

createKalmanState = _fiasco_numpy.createKalmanState

createKalmanFilter = _fiasco_numpy.createKalmanFilter

