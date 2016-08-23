# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _quaternion

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


class Quat(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Quat, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Quat, name)
    __swig_setmethods__["x"] = _quaternion.Quat_x_set
    __swig_getmethods__["x"] = _quaternion.Quat_x_get
    if _newclass:x = property(_quaternion.Quat_x_get, _quaternion.Quat_x_set)
    __swig_setmethods__["y"] = _quaternion.Quat_y_set
    __swig_getmethods__["y"] = _quaternion.Quat_y_get
    if _newclass:y = property(_quaternion.Quat_y_get, _quaternion.Quat_y_set)
    __swig_setmethods__["z"] = _quaternion.Quat_z_set
    __swig_getmethods__["z"] = _quaternion.Quat_z_get
    if _newclass:z = property(_quaternion.Quat_z_get, _quaternion.Quat_z_set)
    __swig_setmethods__["w"] = _quaternion.Quat_w_set
    __swig_getmethods__["w"] = _quaternion.Quat_w_get
    if _newclass:w = property(_quaternion.Quat_w_get, _quaternion.Quat_w_set)
    def __str__(*args): return _quaternion.Quat___str__(*args)
    def __repr__(*args): return _quaternion.Quat___repr__(*args)
    def __init__(self, *args):
        _swig_setattr(self, Quat, 'this', _quaternion.new_Quat(*args))
        _swig_setattr(self, Quat, 'thisown', 1)
    def __del__(self, destroy=_quaternion.delete_Quat):
        try:
            if self.thisown: destroy(self)
        except: pass

class QuatPtr(Quat):
    def __init__(self, this):
        _swig_setattr(self, Quat, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Quat, 'thisown', 0)
        _swig_setattr(self, Quat,self.__class__,Quat)
_quaternion.Quat_swigregister(QuatPtr)


trans_identity = _quaternion.trans_identity

trans_copy = _quaternion.trans_copy

trans_mult_right = _quaternion.trans_mult_right

trans_mult_left = _quaternion.trans_mult_left

trans_transpose = _quaternion.trans_transpose

trans_inverse = _quaternion.trans_inverse

trans_vec_mult = _quaternion.trans_vec_mult

trans_dump = _quaternion.trans_dump

trans_to_quat = _quaternion.trans_to_quat

quat_to_trans = _quaternion.quat_to_trans

quat_copy = _quaternion.quat_copy

quat_nrm_sqrt = _quaternion.quat_nrm_sqrt

quat_normalize = _quaternion.quat_normalize

quat_mult_right = _quaternion.quat_mult_right

quat_mult_left = _quaternion.quat_mult_left

quat_identity = _quaternion.quat_identity

quat_conjugate = _quaternion.quat_conjugate

quat_create = _quaternion.quat_create

quat_from_axis_angle = _quaternion.quat_from_axis_angle

quat_to_axis_angle = _quaternion.quat_to_axis_angle

quat_from_euler_RzRyRx = _quaternion.quat_from_euler_RzRyRx

quat_to_euler_RzRyRx = _quaternion.quat_to_euler_RzRyRx

