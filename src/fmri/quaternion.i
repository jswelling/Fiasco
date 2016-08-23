%module quaternion
%{
#include "quaternion.h"
%}
%include "typemaps.i"
%typemap(in) double [4] (double temp[4]) {
  int i;
  if (!PySequence_Check($input)) {
    PyErr_SetString(PyExc_ValueError,"Expected a sequence");
    return NULL;
  }
  if (PySequence_Length($input) != 4) {
    PyErr_SetString(PyExc_ValueError,"Size mismatch. Expected 4 elements");
    return NULL;
  }
  for (i = 0; i < 4; i++) {
    PyObject *o = PySequence_GetItem($input,i);
    if (PyNumber_Check(o)) {
      temp[i] = PyFloat_AsDouble(o);
    } else {
      PyErr_SetString(PyExc_ValueError,"Sequence elements must be numbers");      
      return NULL;
    }
  }
  $1 = temp;
}
%typemap(in) Vec4 (double temp[4]) {
  int i;
  for (i = 0; i < 4; i++) {
    PyObject *o = PySequence_GetItem($input,i);
    if (PyNumber_Check(o)) {
      temp[i] = PyFloat_AsDouble(o);
    } else {
      PyErr_SetString(PyExc_ValueError,"Sequence elements must be numbers");      
      return NULL;
    }
  }
  $1 = temp;
}
%typemap(in) double [16] (double temp[16]) {
  int i;
  if (!PySequence_Check($input)) {
    PyErr_SetString(PyExc_ValueError,"Expected a sequence");
    return NULL;
  }
  if (PySequence_Length($input) != 16) {
    PyErr_SetString(PyExc_ValueError,"Size mismatch. Expected 16 elements");
    return NULL;
  }
  for (i = 0; i < 16; i++) {
    PyObject *o = PySequence_GetItem($input,i);
    if (PyNumber_Check(o)) {
      temp[i] = PyFloat_AsDouble(o);
    } else {
      PyErr_SetString(PyExc_ValueError,"Sequence elements must be numbers");      
      return NULL;
    }
  }
  $1 = temp;
}
%typemap(in) Transform (double temp[16]) {
  int i;
  for (i = 0; i < 16; i++) {
    PyObject *o = PySequence_GetItem($input,i);
    if (PyNumber_Check(o)) {
      temp[i] = PyFloat_AsDouble(o);
    } else {
      PyErr_SetString(PyExc_ValueError,"Sequence elements must be numbers");      
      return NULL;
    }
  }
  $1 = temp;
}
%typemap(argout) Transform {
  int i;
  for (i = 0; i < 16; i++) {
    PyObject *o = PyFloat_FromDouble((double) $1[i]);
    PySequence_SetItem($input,i,o);
  }
}

%typemap(argout) Vec4 {
  int i;
  for (i = 0; i < 4; i++) {
    PyObject *o = PyFloat_FromDouble((double) $1[i]);
    PySequence_SetItem($input,i,o);
  }
}

%typemap(argout) double [ANY] {
  int i;
  for (i = 0; i < $1_dim0; i++) {
    PyObject *o = PyFloat_FromDouble((double) $1[i]);
    PyList_SetItem($input,i,o);
  }
}
%apply double *INOUT{ double* x, double* y, double* z, double* theta,
			 double* x_angle, double* y_angle, double* z_angle }
%include "quaternion.h"


%extend Quat {
  char *__str__() {
    static char tmp[1024];
    sprintf(tmp,"Quat(%g,%g,%g,%g)", self->x,self->y,self->z,self->w);
    return tmp;
  }
  char *__repr__() {
    static char tmp[1024];
    sprintf(tmp,"Quat(%g,%g,%g,%g)", self->x,self->y,self->z,self->w);
    return tmp;
  }
  Quat(double x, double y, double z, double w) {
    return quat_create(NULL,x,y,z,w);;
  }
}
