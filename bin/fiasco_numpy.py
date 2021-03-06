# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.10
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info as _swig_python_version_info
if _swig_python_version_info >= (2, 7, 0):
    def swig_import_helper():
        import importlib
        pkg = __name__.rpartition('.')[0]
        mname = '.'.join((pkg, '_fiasco_numpy')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_fiasco_numpy')
    _fiasco_numpy = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_fiasco_numpy', [dirname(__file__)])
        except ImportError:
            import _fiasco_numpy
            return _fiasco_numpy
        if fp is not None:
            try:
                _mod = imp.load_module('_fiasco_numpy', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _fiasco_numpy = swig_import_helper()
    del swig_import_helper
else:
    import _fiasco_numpy
del _swig_python_version_info
try:
    _swig_property = property
except NameError:
    pass  # Python < 2.2 doesn't have 'property'.

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_setattr_nondynamic(self, class_type, name, value, static=1):
    if (name == "thisown"):
        return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name, None)
    if method:
        return method(self, value)
    if (not static):
        if _newclass:
            object.__setattr__(self, name, value)
        else:
            self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)


def _swig_setattr(self, class_type, name, value):
    return _swig_setattr_nondynamic(self, class_type, name, value, 0)


def _swig_getattr(self, class_type, name):
    if (name == "thisown"):
        return self.this.own()
    method = class_type.__swig_getmethods__.get(name, None)
    if method:
        return method(self)
    raise AttributeError("'%s' object has no attribute '%s'" % (class_type.__name__, name))


def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except __builtin__.Exception:
    class _object:
        pass
    _newclass = 0

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
    """Proxy of C regressor_struct struct."""

    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Regressor, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Regressor, name)
    __repr__ = _swig_repr

    def __init__(self):
        """__init__(regressor_struct self) -> Regressor"""
        this = _fiasco_numpy.new_Regressor()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _fiasco_numpy.delete_Regressor
    __del__ = lambda self: None

    def get(self, f):
        """get(Regressor self, glm_feature f) -> int"""
        return _fiasco_numpy.Regressor_get(self, f)


    def is_settable(self, f):
        """is_settable(Regressor self, glm_feature f) -> int"""
        return _fiasco_numpy.Regressor_is_settable(self, f)


    def set(self, f, val):
        """set(Regressor self, glm_feature f, int val)"""
        return _fiasco_numpy.Regressor_set(self, f, val)


    def n_params(self, nfactors):
        """n_params(Regressor self, int nfactors) -> int"""
        return _fiasco_numpy.Regressor_n_params(self, nfactors)


    def fit(self, obs, factors, counts, param_out):
        """fit(Regressor self, double * obs, double const * factors, double const * counts, double * param_out)"""
        return _fiasco_numpy.Regressor_fit(self, obs, factors, counts, param_out)


    def context_valid(self, nobs, nfactors):
        """context_valid(Regressor self, int nobs, int nfactors) -> int"""
        return _fiasco_numpy.Regressor_context_valid(self, nobs, nfactors)

Regressor_swigregister = _fiasco_numpy.Regressor_swigregister
Regressor_swigregister(Regressor)


def create_base_regressor():
    """create_base_regressor() -> Regressor"""
    return _fiasco_numpy.create_base_regressor()

def create_llsq_regressor():
    """create_llsq_regressor() -> Regressor"""
    return _fiasco_numpy.create_llsq_regressor()

def create_logistic_regressor():
    """create_logistic_regressor() -> Regressor"""
    return _fiasco_numpy.create_logistic_regressor()

def create_poisson_regressor():
    """create_poisson_regressor() -> Regressor"""
    return _fiasco_numpy.create_poisson_regressor()
class ScalarFunction(_object):
    """Proxy of C ScalarFunction struct."""

    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, ScalarFunction, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, ScalarFunction, name)
    __repr__ = _swig_repr
    __swig_destroy__ = _fiasco_numpy.delete_ScalarFunction
    __del__ = lambda self: None

    def value(self, par):
        """value(ScalarFunction self, double const * par) -> double"""
        return _fiasco_numpy.ScalarFunction_value(self, par)


    def reset(self, par):
        """reset(ScalarFunction self, double const * par)"""
        return _fiasco_numpy.ScalarFunction_reset(self, par)


    def __init__(self):
        """__init__(ScalarFunction self) -> ScalarFunction"""
        this = _fiasco_numpy.new_ScalarFunction()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
ScalarFunction_swigregister = _fiasco_numpy.ScalarFunction_swigregister
ScalarFunction_swigregister(ScalarFunction)


def buildSimpleScalarFunction(valFunc):
    """buildSimpleScalarFunction(SFUNVALUEFUNC valFunc) -> ScalarFunction"""
    return _fiasco_numpy.buildSimpleScalarFunction(valFunc)
class Optimizer(_object):
    """Proxy of C Optimizer struct."""

    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Optimizer, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Optimizer, name)
    __swig_destroy__ = _fiasco_numpy.delete_Optimizer
    __del__ = lambda self: None

    def __str__(self):
        """__str__(Optimizer self) -> char const *"""
        return _fiasco_numpy.Optimizer___str__(self)


    def __repr__(self):
        """__repr__(Optimizer self) -> char *"""
        return _fiasco_numpy.Optimizer___repr__(self)


    def getDebugLevel(self):
        """getDebugLevel(Optimizer self) -> int"""
        return _fiasco_numpy.Optimizer_getDebugLevel(self)


    def setDebugLevel(self, lvl):
        """setDebugLevel(Optimizer self, int const lvl)"""
        return _fiasco_numpy.Optimizer_setDebugLevel(self, lvl)


    def go(self, f, par):
        """go(Optimizer self, ScalarFunction f, double * par) -> int"""
        return _fiasco_numpy.Optimizer_go(self, f, par)


    def __init__(self):
        """__init__(Optimizer self) -> Optimizer"""
        this = _fiasco_numpy.new_Optimizer()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
Optimizer_swigregister = _fiasco_numpy.Optimizer_swigregister
Optimizer_swigregister(Optimizer)


def optimizerFromStringRep(rep):
    """optimizerFromStringRep(char const * rep) -> Optimizer"""
    return _fiasco_numpy.optimizerFromStringRep(rep)

def createBaseOptimizer():
    """createBaseOptimizer() -> Optimizer"""
    return _fiasco_numpy.createBaseOptimizer()

def createNoneOptimizer():
    """createNoneOptimizer() -> Optimizer"""
    return _fiasco_numpy.createNoneOptimizer()

def createPraxisOptimizer(t0, h0):
    """createPraxisOptimizer(double t0, double h0) -> Optimizer"""
    return _fiasco_numpy.createPraxisOptimizer(t0, h0)

def createNelminOptimizer(stopping_val, scale, maxRestarts):
    """createNelminOptimizer(double stopping_val, double scale, int maxRestarts) -> Optimizer"""
    return _fiasco_numpy.createNelminOptimizer(stopping_val, scale, maxRestarts)

def createNelminTOptimizer(stopping_val, scale, T, maxRestarts):
    """createNelminTOptimizer(double stopping_val, double scale, double T, int maxRestarts) -> Optimizer"""
    return _fiasco_numpy.createNelminTOptimizer(stopping_val, scale, T, maxRestarts)
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
    """Proxy of C Interpolator_struct struct."""

    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Interpolator, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Interpolator, name)
    __repr__ = _swig_repr
    __swig_destroy__ = _fiasco_numpy.delete_Interpolator
    __del__ = lambda self: None

    def getTypeName(self):
        """getTypeName(Interpolator self) -> char const *"""
        return _fiasco_numpy.Interpolator_getTypeName(self)


    def prep(self, data):
        """prep(Interpolator self, double * data)"""
        return _fiasco_numpy.Interpolator_prep(self, data)


    def calc(self, calcResult, calcLoc, runLength, offset):
        """calc(Interpolator self, double * calcResult, double * calcLoc, long runLength, long offset)"""
        return _fiasco_numpy.Interpolator_calc(self, calcResult, calcLoc, runLength, offset)


    def dumpSelf(self, ofile):
        """dumpSelf(Interpolator self, FILE * ofile)"""
        return _fiasco_numpy.Interpolator_dumpSelf(self, ofile)


    def setInt(self, which, val):
        """setInt(Interpolator self, int which, long val)"""
        return _fiasco_numpy.Interpolator_setInt(self, which, val)


    def setDouble(self, which, val):
        """setDouble(Interpolator self, int which, double val)"""
        return _fiasco_numpy.Interpolator_setDouble(self, which, val)


    def getInt(self, which):
        """getInt(Interpolator self, int which) -> long"""
        return _fiasco_numpy.Interpolator_getInt(self, which)


    def getDouble(self, which):
        """getDouble(Interpolator self, int which) -> double"""
        return _fiasco_numpy.Interpolator_getDouble(self, which)


    def __init__(self):
        """__init__(Interpolator_struct self) -> Interpolator"""
        this = _fiasco_numpy.new_Interpolator()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
Interpolator_swigregister = _fiasco_numpy.Interpolator_swigregister
Interpolator_swigregister(Interpolator)


def intrp_typeFromName(name):
    """intrp_typeFromName(char const * name) -> InterpolatorType"""
    return _fiasco_numpy.intrp_typeFromName(name)

def intrp_nameFromType(type):
    """intrp_nameFromType(InterpolatorType type) -> char const *"""
    return _fiasco_numpy.intrp_nameFromType(type)

def intrp_warpClearCounts():
    """intrp_warpClearCounts()"""
    return _fiasco_numpy.intrp_warpClearCounts()

def intrp_warpGetCounts(count):
    """intrp_warpGetCounts(long * count)"""
    return _fiasco_numpy.intrp_warpGetCounts(count)

def warpApply(interp, t, orig_image, moved_image, check, nx, ny, nz, fast_blk, length_x, length_y, length_z):
    """warpApply(Interpolator interp, Transform t, double * orig_image, double * moved_image, char * check, long nx, long ny, long nz, long fast_blk, double length_x, double length_y, double length_z)"""
    return _fiasco_numpy.warpApply(interp, t, orig_image, moved_image, check, nx, ny, nz, fast_blk, length_x, length_y, length_z)

def createClosestInterpolator1D(nx, fast_blk):
    """createClosestInterpolator1D(long nx, long fast_blk) -> Interpolator"""
    return _fiasco_numpy.createClosestInterpolator1D(nx, fast_blk)

def createClosestInterpolator2D(nx, ny, fast_blk):
    """createClosestInterpolator2D(long nx, long ny, long fast_blk) -> Interpolator"""
    return _fiasco_numpy.createClosestInterpolator2D(nx, ny, fast_blk)

def createClosestInterpolator3D(nx, ny, nz, fast_blk):
    """createClosestInterpolator3D(long nx, long ny, long nz, long fast_blk) -> Interpolator"""
    return _fiasco_numpy.createClosestInterpolator3D(nx, ny, nz, fast_blk)

def createLinearInterpolator1D(nx, fast_blk):
    """createLinearInterpolator1D(long nx, long fast_blk) -> Interpolator"""
    return _fiasco_numpy.createLinearInterpolator1D(nx, fast_blk)

def createLinearInterpolator2D(nx, ny, fast_blk):
    """createLinearInterpolator2D(long nx, long ny, long fast_blk) -> Interpolator"""
    return _fiasco_numpy.createLinearInterpolator2D(nx, ny, fast_blk)

def createLinearInterpolator3D(nx, ny, nz, fast_blk):
    """createLinearInterpolator3D(long nx, long ny, long nz, long fast_blk) -> Interpolator"""
    return _fiasco_numpy.createLinearInterpolator3D(nx, ny, nz, fast_blk)

def createSplineInterpolator1D(type, nx, fast_blk):
    """createSplineInterpolator1D(SplineType type, long nx, long fast_blk) -> Interpolator"""
    return _fiasco_numpy.createSplineInterpolator1D(type, nx, fast_blk)

def createSplineInterpolator2D(type, nx, ny, fast_blk):
    """createSplineInterpolator2D(SplineType type, long nx, long ny, long fast_blk) -> Interpolator"""
    return _fiasco_numpy.createSplineInterpolator2D(type, nx, ny, fast_blk)

def createSplineInterpolator3D(type, nX, ny, nz, fast_blk):
    """createSplineInterpolator3D(SplineType type, long nX, long ny, long nz, long fast_blk) -> Interpolator"""
    return _fiasco_numpy.createSplineInterpolator3D(type, nX, ny, nz, fast_blk)

def createInterpolator1DByType(type, nx, fast_blk):
    """createInterpolator1DByType(InterpolatorType type, long nx, long fast_blk) -> Interpolator"""
    return _fiasco_numpy.createInterpolator1DByType(type, nx, fast_blk)

def createInterpolator2DByType(type, nx, ny, fast_blk):
    """createInterpolator2DByType(InterpolatorType type, long nx, long ny, long fast_blk) -> Interpolator"""
    return _fiasco_numpy.createInterpolator2DByType(type, nx, ny, fast_blk)

def createInterpolator3DByType(type, nx, ny, nz, fast_blk):
    """createInterpolator3DByType(InterpolatorType type, long nx, long ny, long nz, long fast_blk) -> Interpolator"""
    return _fiasco_numpy.createInterpolator3DByType(type, nx, ny, nz, fast_blk)
class KalmanProcess(_object):
    """Proxy of C KalmanProcess_struct struct."""

    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, KalmanProcess, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, KalmanProcess, name)
    __repr__ = _swig_repr
    __swig_destroy__ = _fiasco_numpy.delete_KalmanProcess
    __del__ = lambda self: None

    def dumpSelf(self, ofile):
        """dumpSelf(KalmanProcess self, FILE * ofile)"""
        return _fiasco_numpy.KalmanProcess_dumpSelf(self, ofile)


    def setDebug(self, val):
        """setDebug(KalmanProcess self, int val)"""
        return _fiasco_numpy.KalmanProcess_setDebug(self, val)


    def getDebug(self):
        """getDebug(KalmanProcess self) -> int"""
        return _fiasco_numpy.KalmanProcess_getDebug(self)


    def setA(self, IN_ARRAY2):
        """setA(KalmanProcess self, double * IN_ARRAY2)"""
        return _fiasco_numpy.KalmanProcess_setA(self, IN_ARRAY2)


    def getA(self, INPLACE_ARRAY2):
        """getA(KalmanProcess self, double * INPLACE_ARRAY2)"""
        return _fiasco_numpy.KalmanProcess_getA(self, INPLACE_ARRAY2)


    def setH(self, IN_ARRAY2):
        """setH(KalmanProcess self, double * IN_ARRAY2)"""
        return _fiasco_numpy.KalmanProcess_setH(self, IN_ARRAY2)


    def getH(self, INPLACE_ARRAY2):
        """getH(KalmanProcess self, double * INPLACE_ARRAY2)"""
        return _fiasco_numpy.KalmanProcess_getH(self, INPLACE_ARRAY2)


    def setQ(self, IN_ARRAY2):
        """setQ(KalmanProcess self, double * IN_ARRAY2)"""
        return _fiasco_numpy.KalmanProcess_setQ(self, IN_ARRAY2)


    def getQ(self, INPLACE_ARRAY2):
        """getQ(KalmanProcess self, double * INPLACE_ARRAY2)"""
        return _fiasco_numpy.KalmanProcess_getQ(self, INPLACE_ARRAY2)


    def setR(self, IN_ARRAY2):
        """setR(KalmanProcess self, double * IN_ARRAY2)"""
        return _fiasco_numpy.KalmanProcess_setR(self, IN_ARRAY2)


    def getR(self, INPLACE_ARRAY2):
        """getR(KalmanProcess self, double * INPLACE_ARRAY2)"""
        return _fiasco_numpy.KalmanProcess_getR(self, INPLACE_ARRAY2)


    def getL(self):
        """getL(KalmanProcess self) -> int"""
        return _fiasco_numpy.KalmanProcess_getL(self)


    def getM(self):
        """getM(KalmanProcess self) -> int"""
        return _fiasco_numpy.KalmanProcess_getM(self)


    def apply(self, state, IN_ARRAY1, missing, calcLogLikelihoodDelta, updatePandK):
        """apply(KalmanProcess self, KalmanState state, double * IN_ARRAY1, int missing, int calcLogLikelihoodDelta, int updatePandK) -> double"""
        return _fiasco_numpy.KalmanProcess_apply(self, state, IN_ARRAY1, missing, calcLogLikelihoodDelta, updatePandK)


    def __init__(self):
        """__init__(KalmanProcess_struct self) -> KalmanProcess"""
        this = _fiasco_numpy.new_KalmanProcess()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
KalmanProcess_swigregister = _fiasco_numpy.KalmanProcess_swigregister
KalmanProcess_swigregister(KalmanProcess)

class KalmanState(_object):
    """Proxy of C KalmanState_struct struct."""

    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, KalmanState, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, KalmanState, name)
    __repr__ = _swig_repr

    def dumpSelf(self, ofile):
        """dumpSelf(KalmanState self, FILE * ofile)"""
        return _fiasco_numpy.KalmanState_dumpSelf(self, ofile)


    def setX(self, INPLACE_ARRAY1):
        """setX(KalmanState self, double * INPLACE_ARRAY1)"""
        return _fiasco_numpy.KalmanState_setX(self, INPLACE_ARRAY1)


    def setP(self, INPLACE_ARRAY2):
        """setP(KalmanState self, double * INPLACE_ARRAY2)"""
        return _fiasco_numpy.KalmanState_setP(self, INPLACE_ARRAY2)


    def getM(self):
        """getM(KalmanState self) -> int"""
        return _fiasco_numpy.KalmanState_getM(self)


    def __init__(self):
        """__init__(KalmanState_struct self) -> KalmanState"""
        this = _fiasco_numpy.new_KalmanState()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
KalmanState_swigregister = _fiasco_numpy.KalmanState_swigregister
KalmanState_swigregister(KalmanState)

class KalmanFilter(_object):
    """Proxy of C KalmanFilter_struct struct."""

    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, KalmanFilter, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, KalmanFilter, name)
    __repr__ = _swig_repr
    __swig_destroy__ = _fiasco_numpy.delete_KalmanFilter
    __del__ = lambda self: None

    def dumpSelf(self, ofile):
        """dumpSelf(KalmanFilter self, FILE * ofile)"""
        return _fiasco_numpy.KalmanFilter_dumpSelf(self, ofile)


    def getL(self):
        """getL(KalmanFilter self) -> int"""
        return _fiasco_numpy.KalmanFilter_getL(self)


    def getM(self):
        """getM(KalmanFilter self) -> int"""
        return _fiasco_numpy.KalmanFilter_getM(self)


    def getState(self):
        """getState(KalmanFilter self) -> KalmanState"""
        return _fiasco_numpy.KalmanFilter_getState(self)


    def run(self, tdim, zIn, xOut, POut, calcLogLikelihood):
        """run(KalmanFilter self, long tdim, double const * zIn, double * xOut, double * POut, int calcLogLikelihood) -> double"""
        return _fiasco_numpy.KalmanFilter_run(self, tdim, zIn, xOut, POut, calcLogLikelihood)


    def __init__(self):
        """__init__(KalmanFilter_struct self) -> KalmanFilter"""
        this = _fiasco_numpy.new_KalmanFilter()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
KalmanFilter_swigregister = _fiasco_numpy.KalmanFilter_swigregister
KalmanFilter_swigregister(KalmanFilter)


def createKalmanProcess(L, M):
    """createKalmanProcess(int L, int M) -> KalmanProcess"""
    return _fiasco_numpy.createKalmanProcess(L, M)

def createKalmanState(M):
    """createKalmanState(int M) -> KalmanState"""
    return _fiasco_numpy.createKalmanState(M)

def createKalmanFilter(process):
    """createKalmanFilter(KalmanProcess process) -> KalmanFilter"""
    return _fiasco_numpy.createKalmanFilter(process)
# This file is compatible with both classic and new-style classes.


