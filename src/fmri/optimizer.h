/************************************************************
 *                                                          *
 *  optimizer.h                                              *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1998 Department of Statistics,         *
 *                        Carnegie Mellon University        *
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
 *  Original programming by Joel Welling, 2/2004            *
 ************************************************************/
/* This package implements smoothing methods. */

#ifndef INCL_OPTIMIZER_H
#define INCL_OPTIMIZER_H 1

typedef struct ScalarFunction {
  double (*value)( struct ScalarFunction* self, 
		   const double* par, const int nPar );
  void (*reset)( struct ScalarFunction* self, const double* par, 
		 const int nPar );
  void (*destroySelf)( struct ScalarFunction* self );
  int nPar;
  void* data;
} ScalarFunction;

typedef struct Optimizer {
  const char* (*getMethodName)(struct Optimizer* self);
  char* (*getStringRep)(struct Optimizer* self);
  void (*setDebugLevel)(struct Optimizer* self, const int lvl);
  int (*getDebugLevel)(struct Optimizer* self);
  void (*setTol)(struct Optimizer* self, const double tol);
  double (*getTol)(struct Optimizer* self);
  void (*setScale)(struct Optimizer* self, const double scale);
  double (*getScale)(struct Optimizer* self);
	  
  int (*go)( struct Optimizer* self, ScalarFunction* f, 
	     double* par, const int nPar, double* best);
  void (*destroySelf)(struct Optimizer* self);
  int debugLevel;
  void* data;
} Optimizer;

ScalarFunction* buildSimpleScalarFunction( 
	     double (*value)(const double*, const int, void*),
	     void (*reset)(const double*, const int, void*),
	     const int nDim, void* userHook);

Optimizer* optimizerFromStringRep( const char* rep );
Optimizer* createBaseOptimizer(void);
Optimizer* createNoneOptimizer(void);
Optimizer* createPraxisOptimizer(double t0, double h0);
Optimizer* createNelminOptimizer(double stopping_val, double scale, 
				 int maxRestarts);
Optimizer* createNelminTOptimizer(double stopping_val, double scale, double T,
				  int maxRestarts);

#endif  /* ifndef INCL_OPTIMIZER_H */
