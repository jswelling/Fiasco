/* -*- C -*-  (not really, but good for syntax highlighting) */

/*******************************************************
 * Notes-
 * -ScalarFunction's python wrapper does a Py_INCREF on a
 *  tuple which is passed in in the constructor arg list.
 *  This should get DECREF'd when the tool is deleted, but
 *  I can't figure out how to cause SWIG to make this happen.
 *******************************************************/

%module fiasco_numpy
%{
#define SWIG_FILE_WITH_INIT
#include "numpy/old_defines.h"
#include "glm.h"
#include "optimizer.h"
#include "quaternion.h"
#include "interpolator.h"
#include "kalmanfilter.h"
%}

%include "typemaps.i"
%include "numpy.i"

%init %{ 
import_array(); 
%}

typedef enum {
  GLM_COMPLEX,
  GLM_RESIDUALS,
  GLM_VARIANCES,
  GLM_COVARIANCES,
  GLM_SSQR,
  GLM_ORTHO,
  GLM_DEBUG,
  GLM_DEVIANCE,
  GLM_TYPE
} glm_feature;

typedef enum {
  GLM_TYPE_BASE,      /* does nothing */
  GLM_TYPE_LLSQ,      /* linear least squares */
  GLM_TYPE_LOGISTIC,  /* logistic regression */
  GLM_TYPE_POISSON    /* Poisson regresson (link function is ln()) */
} glm_type;

%apply (double* IN_ARRAY2, int DIM1, int DIM2) {
  (const double* factors, int nfactors_1, int nfactors_2)
    };
%apply (double* IN_ARRAY1, int DIM1) {
    (const double* counts, int ncounts)
      }
%apply (double* INPLACE_ARRAY1, int DIM1) {
  (double* obs, int nobs),
    (double* param_out, int nparam)
    };

%exception {
$action
  if (PyErr_Occurred()) SWIG_fail;
}

%feature("autodoc","1");
typedef struct regressor_struct{
  %extend {
    Regressor() { return glm_create_llsq_regressor(); }
    ~Regressor() { self->destroy_self(self); }
    int get(glm_feature f) { return(self->get(self,f)); }
    int is_settable(glm_feature f) { return(self->is_settable(self,f)); }
    void set( glm_feature f, int val) {
      self->set(self,f,val);
      if (glm_error_msg()) {
	PyErr_SetString(PyExc_RuntimeError,glm_error_msg());
      }
    }
    int n_params(int nfactors) { return(self->n_params(self,nfactors)); }
    void fit(double* obs, int nobs,
	     const double* factors, int nfactors_1, int nfactors_2, 
	     const double* counts, int ncounts,
	     double* param_out, int nparam);
    int context_valid(int nobs, int nfactors)
    {
      return (self->context_valid(self,nobs,nfactors));
    }
  }
} Regressor;

%{
void regressor_struct_fit(struct regressor_struct* r, double* obs, int nobs,
			  const double* factors, int nfactors_1, int nfactors_2,
			  const double* counts, int ncounts,
			  double* param_out, int nparam) {
  if (nfactors_1 != nobs) {
    PyErr_SetString(PyExc_RuntimeError,
		    "factor vector length does not match obs vector");
  }
  if ((counts != NULL) && (ncounts != nobs)) {
    PyErr_SetString(PyExc_RuntimeError,
		    "counts length does not match obs vector");
  }
  if (nparam < r->n_params(r,nfactors_2)) {
    PyErr_SetString(PyExc_RuntimeError,"param_out array is too small");
    return;
  }
  if (r->fit(r, obs, factors, counts, param_out, nobs, nfactors_2)) {
    PyErr_SetString(PyExc_RuntimeError,glm_error_msg());
  }
}
%}

%rename(create_base_regressor) glm_create_base_regressor;
Regressor* glm_create_base_regressor(void);
%rename(create_llsq_regressor) glm_create_llsq_regressor;
Regressor* glm_create_llsq_regressor(void);
%rename(create_logistic_regressor) glm_create_logistic_regressor;
Regressor* glm_create_logistic_regressor(void);
%rename(create_poisson_regressor) glm_create_poisson_regressor;
Regressor* glm_create_poisson_regressor(void);

%{
  typedef double (*SFUNVALUEFUNC)(const double* vals, const int n, 
				  void* clientdata);
  typedef void (*SFUNRESETFUNC)(const double* vals, const int n, 
				void* clientdata);

/* This function matches the prototype of the normal C value callback
   function for ScalarFunction. However, we use the clientdata pointer
   for holding a reference to a Python callable object. */
  static double PythonSFunValueFunc(const double* vals, const int n,
				    void *clientdata)
{
    PyObject *valFunc, *resetFunc, *data, *valList, *arglist;
    PyObject *valArray;
    PyObject *resultObj;
    double result;
    int dimensions[1];
    
    // Break data into components
    valFunc= PyTuple_GetItem(clientdata,0);
    resetFunc= PyTuple_GetItem(clientdata,1);
    data= PyTuple_GetItem(clientdata,3);
    
    dimensions[0]= n;
    valArray= PyArray_FromDimsAndData(1,dimensions,PyArray_DOUBLE,(char*)vals);

    arglist = Py_BuildValue("OO",valArray,data);        // Build argument list
    resultObj = PyEval_CallObject(valFunc,arglist);     // Call Python
    Py_DECREF(arglist);                                 // Trash arglist
    if (resultObj) result= PyFloat_AsDouble(resultObj);
    else PyErr_SetString(PyExc_RuntimeError,"Optimizer value function failed");
    Py_XDECREF(resultObj);
    return result;
}

/* This function matches the prototype of the normal C reset callback
   function for ScalarFunction. However, we use the clientdata pointer
   for holding a reference to a Python callable object. */
  static void PythonSFunResetFunc(const double* vals, const int n,
				    void *clientdata)
{
    PyObject *valFunc, *resetFunc, *data, *valList, *arglist;
    PyObject *valArray;
    PyObject *resultObj;
    double result;
    int dimensions[1];
    
    // Break data into components
    valFunc= PyTuple_GetItem(clientdata,0);
    resetFunc= PyTuple_GetItem(clientdata,1);
    data= PyTuple_GetItem(clientdata,3);
    
    dimensions[0]= n;
    valArray= PyArray_FromDimsAndData(1,dimensions,PyArray_DOUBLE,(char*)vals);

    arglist = Py_BuildValue("OO",valArray,data);        // Build argument list
    resultObj = PyEval_CallObject(resetFunc,arglist);     // Call Python
    Py_DECREF(arglist);                                 // Trash arglist
    if (resultObj) Py_XDECREF(resultObj);
    else PyErr_SetString(PyExc_RuntimeError,"Optimizer reset function failed");
}
%}

%apply (double* INPLACE_ARRAY1, int DIM1) {
  (const double* par, const int nPar)
}

%feature("autodoc","1");
typedef struct ScalarFunction {
  %extend {
    ~ScalarFunction() { self->destroySelf(self); }
    double value( const double* par, const int nPar ) {
      return self->value(self, par, nPar);
    }
    void reset(const double* par, const int nPar) {
      self->reset(self, par, nPar);
    }
  }
} ScalarFunction;

/* The expected input here is a Python Tuple containing two functions,
 * an integer and a pointer for data.
 */
%typemap(in) (SFUNVALUEFUNC valFunc, SFUNRESETFUNC resetFunc, 
	      const int nDim, void* cbData) {
  if (!PyTuple_Check($input)) {
    PyErr_SetString(PyExc_TypeError, "Expected a tuple!");
    return NULL;
  }
  if (PyTuple_Size($input)!=4) {
    PyErr_SetString(PyExc_TypeError, "Tuple size is not 4!");
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
  if (!PyInt_Check(PyTuple_GetItem($input,2))) {
    PyErr_SetString(PyExc_TypeError, 
		    "Third tuple element is not an int!");
    return NULL;
  }
  $3= (int)PyInt_AsLong(PyTuple_GetItem($input,2));

  Py_INCREF($input);
  $1= PythonSFunValueFunc; $2= PythonSFunResetFunc; $4= $input;
}

%feature("autodoc","1");
ScalarFunction* buildSimpleScalarFunction( SFUNVALUEFUNC valFunc,
					   SFUNRESETFUNC resetFunc,
					   const int nDim, void* cbData);


%feature("autodoc","1");
typedef struct Optimizer {
/*   void (*setTol)(struct Optimizer* self, const double tol); */
/*   double (*getTol)(struct Optimizer* self); */
/*   void (*setScale)(struct Optimizer* self, const double scale); */
/*   double (*getScale)(struct Optimizer* self); */
	  
/*   int (*go)( struct Optimizer* self, ScalarFunction* f,  */
/* 	     double* par, const int nPar, double* best); */
  %extend {
    ~Optimizer() { self->destroySelf(self); }
    const char* __str__() { return self->getMethodName(self); }
    char* __repr__() { return self->getStringRep(self); }
    int getDebugLevel() { return self->getDebugLevel(self); }
    void setDebugLevel(const int lvl) { self->setDebugLevel(self,lvl); }
    %apply(double* INPLACE_ARRAY1, int DIM1) {
      (double* par, const int nPar)
    };
    %apply double *OUTPUT { double* best };
    
    int go( ScalarFunction* f, double* par, const int nPar, double* best )
    { return self->go(self,f,par,nPar,best); }
  }
} Optimizer;

%feature("autodoc","1");
Optimizer* optimizerFromStringRep( const char* rep );
Optimizer* createBaseOptimizer(void);
Optimizer* createNoneOptimizer(void);
Optimizer* createPraxisOptimizer(double t0, double h0);
Optimizer* createNelminOptimizer(double stopping_val, double scale, 
				 int maxRestarts);
Optimizer* createNelminTOptimizer(double stopping_val, double scale, double T,
				  int maxRestarts);

typedef enum { INTRP_CLOSEST, INTRP_LINEAR, INTRP_CATMULLROM, 
	       INTRP_BEZIER, INTRP_BSPLINE, INTRP_UNKNOWN } InterpolatorType;

enum { INTRP_OPT_DEBUG=0, INTRP_OPT_FASTBLK, INTRP_OPT_EXTENT, 
       INTRP_OPT_TENSION };

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
  bcopy($1,temp,sizeof(Transform));
}
%typemap(in) (FILE* ofile) {
%#if PY_MAJOR_VERSION < 3  
  if (PyFile_Check($input)) {
    $1= PyFile_AsFile($input);
  } else {
    PyErr_SetString(PyExc_ValueError,"Input must be a file");
    return NULL;
  }
%#else
  PyErr_SetString(PyExc_ValueError,"Python 3 does not support use of files in this way");
  return NULL;
%#endif
}

%typemap(in) (double* data, long sz) (PyArrayObject* array=NULL) {
  array = obj_to_array_no_conversion($input, NPY_DOUBLE);
  if (!array || !require_contiguous(array)
      || !require_native(array)) SWIG_fail;
  $1 = (double*) array_data(array);
  $2 = 1;
  {
    int i;
    for (i=0; i < array_numdims(array); ++i) $2 *= array_size(array,i);
  }
}

%typemap(in) (double* calcResult) (PyArrayObject* array=NULL) {
  array = obj_to_array_no_conversion($input, NPY_DOUBLE);
  if (!array || !require_contiguous(array)
      || !require_native(array)) SWIG_fail;
  $1 = (double*) array_data(array);
}

%typemap(in) (double* calcLoc) (double temp[16]) {
  int i;
  if (!PySequence_Check($input)) SWIG_fail;
  if (PySequence_Size($input)>16) SWIG_fail; /* arbitrary limit on dims */
  for (i=0; i<PySequence_Size($input); i++) {
    PyObject *o= PySequence_GetItem($input,i);
    if (PyNumber_Check(o)) temp[i]= PyFloat_AsDouble(o);
    else SWIG_fail;
  }
  $1= temp;
}

%feature("autodoc","1");
typedef struct Interpolator_struct {
  %extend {
    ~Interpolator() { self->destroySelf(self); }
    const char* getTypeName() { return self->typeName; }
    void prep( double* data, long sz ) 
    { 
      self->prep(self, data, sz); 
    }
    void calc( double* calcResult, double* calcLoc, long runLength, 
	       long offset )
    { 
      self->calc(self, calcResult, calcLoc, runLength, offset); 
    }
    void dumpSelf( FILE* ofile ) { self->dumpSelf(self, ofile); }
    void setInt( int which, long val ) { self->setInt(self, which, val); }
    void setDouble( int which, double val )
    { self->setDouble(self,which,val); }
    long getInt( int which ) { return(self->getInt(self,which)); }
    double getDouble( int which ) { return(self->getDouble(self,which)); }
  }
} Interpolator;

%feature("autodoc","1");
extern InterpolatorType intrp_typeFromName( const char* name );
extern const char* intrp_nameFromType( InterpolatorType type );
extern void intrp_warpClearCounts(void);
extern void intrp_warpGetCounts(long* count);

%rename(warpApply) intrp_warpApply;
void intrp_warpApply( Interpolator* interp, Transform t,
		      double* orig_image,
		      double* moved_image,
		      char* check,
		      long nx, long ny, long nz, long fast_blk,
		      double length_x, double length_y, double length_z );

%rename(createClosestInterpolator1D) intrp_createClosestInterpolator1D;
extern Interpolator* 
intrp_createClosestInterpolator1D( long nx, long fast_blk );

%rename(createClosestInterpolator2D) intrp_createClosestInterpolator2D;
extern Interpolator* 
intrp_createClosestInterpolator2D( long nx, long ny, long fast_blk );

%rename(createClosestInterpolator3D) intrp_createClosestInterpolator3D;
extern Interpolator* 
intrp_createClosestInterpolator3D( long nx, long ny, long nz, long fast_blk );

%rename(createLinearInterpolator1D) intrp_createLinearInterpolator1D;
extern Interpolator* 
intrp_createLinearInterpolator1D( long nx, long fast_blk );

%rename(createLinearInterpolator2D) intrp_createLinearInterpolator2D;
extern Interpolator* 
intrp_createLinearInterpolator2D( long nx, long ny, long fast_blk );

%rename(createLinearInterpolator3D) intrp_createLinearInterpolator3D;
extern Interpolator* 
intrp_createLinearInterpolator3D( long nx, long ny, long nz, long fast_blk );

%rename(createSplineInterpolator1D) intrp_createSplineInterpolator1D;
extern Interpolator* 
intrp_createSplineInterpolator1D( SplineType type, long nx, long fast_blk );

%rename(createSplineInterpolator2D) intrp_createSplineInterpolator2D;
extern Interpolator* 
intrp_createSplineInterpolator2D( SplineType type, long nx, long ny, 
				  long fast_blk );

%rename(createSplineInterpolator3D) intrp_createSplineInterpolator3D;
extern Interpolator* 
intrp_createSplineInterpolator3D( SplineType type, long nX, long ny, long nz, 
				  long fast_blk );

%rename(createInterpolator1DByType) intrp_createInterpolator1DByType;
extern Interpolator* 
intrp_createInterpolator1DByType( InterpolatorType type, 
				  long nx, long fast_blk );
%rename(createInterpolator2DByType) intrp_createInterpolator2DByType;
extern Interpolator* 
intrp_createInterpolator2DByType( InterpolatorType type, 
				   long nx, long ny, long fast_blk );
%rename(createInterpolator3DByType) intrp_createInterpolator3DByType;
extern Interpolator* 
intrp_createInterpolator3DByType( InterpolatorType type, 
				  long nx, long ny, long nz, long fast_blk );
						       
%feature("autodoc","1");

%apply(double* IN_ARRAY2, int DIM1, int DIM2) {
  (const double* AIn, int dim1, int dim2)
}

%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {
  (double* AOut, int dim1, int dim2)
}

typedef struct KalmanProcess_struct {
  %extend {
    ~KalmanProcess() { self->destroySelf(self); }
    void dumpSelf( FILE* ofile ) { self->dumpSelf(self,ofile); }
    void setDebug( int val ) { self->setDebug(self,val); }
    int getDebug() { return self->getDebug(self); }
    void setA( double* IN_ARRAY2, int DIM1, int DIM2 ) { 
      if (DIM1 != self->M || DIM2!= self->M) {
	PyErr_SetString(PyExc_RuntimeError,"array size mismatch");
	return;
      }
      self->setA(self,IN_ARRAY2,DIM1*DIM2); 
    }
    void getA( double* INPLACE_ARRAY2, int DIM1, int DIM2 ) { 
      if (DIM1 != self->M || DIM2 != self->M) {
	PyErr_SetString(PyExc_RuntimeError,"array size mismatch");
	return;
      }
      self->getA(self,INPLACE_ARRAY2,DIM1*DIM2); 
    }
    void setH( double* IN_ARRAY2, int DIM1, int DIM2 ) { 
      if (DIM1!=self->L || DIM2!=self->M) {
	PyErr_SetString(PyExc_RuntimeError,"array size mismatch");
	return;
      }
      self->setH(self,IN_ARRAY2,DIM1*DIM2); 
    }
    void getH( double* INPLACE_ARRAY2, int DIM1, int DIM2 ) { 
      if (DIM1 != self->L || DIM2 != self->M) {
	PyErr_SetString(PyExc_RuntimeError,"array size mismatch");
	return;
      }
      self->getH(self,INPLACE_ARRAY2,DIM1*DIM2); 
    }
    void setQ( double* IN_ARRAY2, int DIM1, int DIM2 ) { 
      if (DIM1!=self->M || DIM2!=self->M) {
	PyErr_SetString(PyExc_RuntimeError,"array size mismatch");
	return;
      }
      self->setQ(self,IN_ARRAY2,DIM1*DIM2); 
    }
    void getQ( double* INPLACE_ARRAY2, int DIM1, int DIM2 ) { 
      if (DIM1 != self->M || DIM2 != self->M) {
	PyErr_SetString(PyExc_RuntimeError,"array size mismatch");
	return;
      }
      self->getQ(self,INPLACE_ARRAY2,DIM1*DIM2); 
    }
    void setR( double* IN_ARRAY2, int DIM1, int DIM2 ) { 
      if (DIM1!=self->L || DIM2!=self->L) {
	PyErr_SetString(PyExc_RuntimeError,"array size mismatch");
	return;
      }
      self->setR(self,IN_ARRAY2,DIM1*DIM2); 
    }
    void getR( double* INPLACE_ARRAY2, int DIM1, int DIM2 ) { 
      if (DIM1 != self->L || DIM2 != self->L) {
	PyErr_SetString(PyExc_RuntimeError,"array size mismatch");
	return;
      }
      self->getR(self,INPLACE_ARRAY2,DIM1*DIM2); 
    }
    int getL() { return self->L; }
    int getM() { return self->M; }
    double apply( KalmanState* state, double* IN_ARRAY1, int DIM1,
		  int missing, int calcLogLikelihoodDelta, int updatePandK )
    { 
      if (DIM1 != self->L) {
	PyErr_SetString(PyExc_RuntimeError,"wrong z vector length");
	return 0.0;
      }
      return( self->apply(self,state,IN_ARRAY1,DIM1,missing,
			  calcLogLikelihoodDelta,
			  updatePandK)); 
    }
  }
} KalmanProcess;

typedef struct KalmanState_struct {
  %extend {
    ~KalmanProcess() { self->destroySelf(self); }
    void dumpSelf( FILE* ofile ) { self->dumpSelf(self,ofile); }
    void setX( double* INPLACE_ARRAY1, int DIM1 ) { 
      if (DIM1!=self->M) {
	PyErr_SetString(PyExc_RuntimeError,"wrong x vector length");
	return;
      }
      self->setX(self,INPLACE_ARRAY1,DIM1); 
    }
    void setP( double* INPLACE_ARRAY2, int DIM1, int DIM2 ) { 
      if (DIM1!=self->M || DIM2!=self->M) {
	PyErr_SetString(PyExc_RuntimeError,"array size mismatch");
	return;
      }
      self->setP(self,INPLACE_ARRAY2,DIM1*DIM2); 
    }
    int getM() { return self->M; }
  }
} KalmanState;

%apply ( double* IN_ARRAY2, int DIM1, int DIM2) {
  ( const double* zIn, int zDim1, int zDim2 )
};
%apply ( double* INPLACE_ARRAY2, int DIM1, int DIM2) {
  ( double* xOut, int xDim1, int xDim2 )
};
%apply ( double* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {
  ( double* POut, int pDim1, int pDim2, int pDim3 )
};
typedef struct KalmanFilter_struct {
  %extend {
    ~KalmanFilter() { self->destroySelf(self); }
    void dumpSelf( FILE* ofile ) { self->dumpSelf(self,ofile); }
    int getL() { return self->L; }
    int getM() { return self->M; }
    KalmanState* getState() { return self->state; }
    double run( long tdim,
		const double* zIn, int zDim1, int zDim2,
		double* xOut, int xDim1, int xDim2,
		double* POut, int pDim1, int pDim2, int pDim3,
		int calcLogLikelihood ) {
      if (zDim2!=self->L || zDim1!=tdim) {
	PyErr_SetString(PyExc_RuntimeError,"array size mismatch on zIn");
	return 0.0;
      }
      if (xDim2!=self->M || xDim1!=tdim) {
	PyErr_SetString(PyExc_RuntimeError,"array size mismatch on xOut");
	return 0.0;
      }
      if (pDim3!=self->M || pDim2!=self->M || pDim1!=tdim) {
	PyErr_SetString(PyExc_RuntimeError,"array size mismatch on POut");
	return 0.0;
      }
      return self->run(self, tdim, 
		       zIn, zDim1*zDim2,
		       xOut, xDim1*xDim2,
		       POut, pDim1*pDim2*pDim3,
		       calcLogLikelihood);
    }
  }
} KalmanFilter;

%rename(createKalmanProcess) klmn_createRowMajorKalmanProcess;
extern KalmanProcess* klmn_createRowMajorKalmanProcess(int L, int M);
%rename(createKalmanState) klmn_createRowMajorKalmanState;
extern KalmanState* klmn_createRowMajorKalmanState(int M);
%rename(createKalmanFilter) klmn_createKalmanFilter;
extern KalmanFilter* klmn_createKalmanFilter(KalmanProcess* process);
