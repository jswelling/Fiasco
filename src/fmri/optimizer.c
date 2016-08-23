/************************************************************
 *                                                          *
 *  optimizer.c                                              *
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "mri.h"
#include "fmri.h" /* Includes optimizer.h */
#include "lapack.h"
#include "misc.h"
#include "stdcrg.h"

static char rcsid[] = "$Id: optimizer.c,v 1.6 2005/06/01 19:54:25 welling Exp $";

static double praxis_machep= 0.0;

typedef struct SSFData {
  double (*value)(const double*, const int, void*);
  void (*reset)(const double*, const int, void*);
  void* userHook;
  int n;
} SSFData;

static void ssfDestroySelf( ScalarFunction* self )
{
  if (self->data) free(self->data);
  free(self);
}

static double ssfValue( ScalarFunction* sf, const double* par, const int nPar )
{
  SSFData* d= (SSFData*)(sf->data);
  if (nPar != d->n) 
    Abort("ssfValue: dimensionalities %d and %d do not match!\n",nPar,d->n);
  return d->value(par, d->n, d->userHook );
}

static void ssfReset( ScalarFunction* sf, const double* par, const int nPar )
{
  SSFData* d= (SSFData*)(sf->data);
  if (nPar != d->n) 
    Abort("ssfReset: dimensionalities %d and %d do not match!\n",nPar,d->n);
  d->reset(par, d->n, d->userHook );
}

ScalarFunction* buildSimpleScalarFunction( 
	     double (*value)(const double*, const int, void*),
	     void (*reset)(const double*, const int, void*),
	     const int nDim, void* userHook)
{
  ScalarFunction* result= NULL;
  SSFData* d= NULL;
  if (!(result=(ScalarFunction*)malloc(sizeof(ScalarFunction))))
    Abort("buildSimpleScalarFunction: unable to allocate %d bytes!\n",
	  sizeof(ScalarFunction));
  if (!(d=(SSFData*)malloc(sizeof(SSFData))))
    Abort("buildSimpleScalarFunction: unable to allocate %d bytes!\n",
	  sizeof(SSFData));
  d->value= value;
  d->reset= reset;
  d->n= nDim;
  d->userHook= userHook;
  result->data= d;
  result->destroySelf= ssfDestroySelf;
  result->value= ssfValue;
  result->reset= ssfReset;

  return result;
}

static const char* baseGetMethodName(Optimizer* self)
{
  static const char* name= "base";
  return name;
}

static char* baseGetStringRep(Optimizer* self)
{
  return strdup("BaseOptimizer()");
}

static void baseDestroySelf(Optimizer* self)
{
  if (self->data) free(self->data);
  free(self);
}

static void baseSetDebugLevel(Optimizer* self, const int lvl)
{
  self->debugLevel= lvl;
}

static int baseGetDebugLevel(Optimizer* self)
{
  return self->debugLevel;
}

static int baseGo(Optimizer* self, ScalarFunction* f, double* par, 
		    const int nPar, double* best)
{
  Abort("baseGo: tried to run an instance of base Optimizer class!\n");
  return 0;
}

static void baseSetTol(Optimizer* self, const double tol)
{
  /* no-op */
  if (self->debugLevel>0) 
    Message("Optimizer tol set to %lf\n",tol);
}

static double baseGetTol(Optimizer* self)
{
  return 0.0;
}

static void baseSetScale(Optimizer* self, const double scale)
{
  /* no-op */
  if (self->debugLevel>0) 
    Message("Optimizer scale set to %lf\n",scale);
}

static double baseGetScale(Optimizer* self)
{
  return 1.0;
}

Optimizer* createBaseOptimizer()
{
  Optimizer* result;

  if (!(result=(Optimizer*)malloc(sizeof(Optimizer))))
    Abort("createBaseOptimizer: unable to allocate %d bytes!\n",
	  sizeof(Optimizer));

  result->data= NULL;
  result->debugLevel= 0;
  result->getMethodName= baseGetMethodName;
  result->getStringRep= baseGetStringRep;
  result->destroySelf= baseDestroySelf;
  result->go= baseGo;
  result->setDebugLevel= baseSetDebugLevel;
  result->getDebugLevel= baseGetDebugLevel;
  result->setTol= baseSetTol;
  result->getTol= baseGetTol;
  result->setScale= baseSetScale;
  result->getScale= baseGetScale;
  
  return result;
}

static const char* noneGetMethodName(Optimizer* self)
{
  static const char* name= "none";
  return name;
}

static char* noneGetStringRep(Optimizer* self)
{
  return strdup("NoneOptimizer()");
}

static int noneGo(Optimizer* self, ScalarFunction* f, double* par, 
		    const int nPar, double* best)
{
  int i;
  *best= f->value(f,par,nPar);
  return 1;
}

Optimizer* createNoneOptimizer()
{
  Optimizer* result= createBaseOptimizer();

  result->getMethodName= noneGetMethodName;
  result->getStringRep= noneGetStringRep;
  result->go= noneGo;
  return result;
}

typedef struct PraxisData {
  double t0;
  double h0;
} PraxisData;

static const char* praxisGetMethodName(Optimizer* self)
{
  static const char* name= "praxis";
  return name;
}

static char* praxisGetStringRep(Optimizer* self)
{
  PraxisData* pd= (PraxisData*)(self->data);
  char buf[256];
  snprintf(buf,sizeof(buf),
	   "PraxisOptimizer(%lg,%lg)",pd->t0,pd->h0);
  return strdup(buf);
}

static double praxisValue( double* par, int nPar, 
			   void* p)
{
  ScalarFunction* sf= (ScalarFunction*)p;
  return sf->value(sf,par,nPar);
}

static void praxisReset( double* par, int nPar, 
			   void* p)
{
  ScalarFunction* sf= (ScalarFunction*)p;
  sf->reset(sf,par,nPar);
}

static void praxisSetTol(Optimizer* self, const double tol)
{
  PraxisData* pd= (PraxisData*)(self->data);
  baseSetTol(self,tol);
  pd->t0= tol;
}

static double praxisGetTol(Optimizer* self)
{
  PraxisData* pd= (PraxisData*)(self->data);
  return pd->t0;
}

static void praxisSetScale(Optimizer* self, const double scale)
{
  PraxisData* pd= (PraxisData*)(self->data);
  baseSetScale(self,scale);
  pd->h0= scale;
}

static double praxisGetScale(Optimizer* self)
{
  PraxisData* pd= (PraxisData*)(self->data);
  return pd->h0;
}

static int praxisGo(Optimizer* self, ScalarFunction* f, double* par, 
		    const int nPar, double* best)
{
  PraxisData* pd= (PraxisData*)(self->data);
  int prx_npar= nPar;

  *best= praxis(pd->t0, praxis_machep, pd->h0, prx_npar, self->debugLevel, 
		   par, praxisValue, praxisReset, 0.0, f);
  return 1;
}

Optimizer* createPraxisOptimizer(double t0, double h0)
{
  Optimizer* result= createBaseOptimizer();
  PraxisData* pd= NULL;

  if (result->data) free(result->data);
  if (!(pd=(PraxisData*)malloc(sizeof(PraxisData))))
    Abort("createPraxisOptimizer: unable to allocate %d bytes!\n");
  result->data= pd;
  pd->t0= t0;
  pd->h0= h0;
  if (praxis_machep==0.0) {
    praxis_machep= (2.0*DLAMCH("e"));
#ifdef never
    fprintf(stderr,"machine eps = %lg\n",praxis_machep);
#endif
  }
  if ((1.0+praxis_machep) == 1.0) 
    Abort("createPraxisOptimizer: precision test failed!\n");

  result->getMethodName= praxisGetMethodName;
  result->getStringRep= praxisGetStringRep;
  result->setTol= praxisSetTol;
  result->getTol= praxisGetTol;
  result->setScale= praxisSetScale;
  result->getScale= praxisGetScale;
  result->go= praxisGo;

  return result;
}

typedef struct NelminData {
  float* startBuf;
  float* valBuf;
  float* scaleBuf;
  double* doubleBuf;
  int bufLength;
  int nPar;
  float scale;
  float stopping_val;
  int maxRestarts;
  ScalarFunction* sf;
} NelminData;

static const char* nelminGetMethodName(Optimizer* self)
{
  static const char* name= "nelmin";
  return name;
}

static char* nelminGetStringRep(Optimizer* self)
{
  NelminData* pd= (NelminData*)(self->data);
  char buf[256];
  snprintf(buf,sizeof(buf),
	   "NelminOptimizer(%g,%g,%d)",pd->stopping_val,pd->scale,
	   pd->maxRestarts);
  return strdup(buf);
}

static float nelminValue( float* par, void* p)
{
  NelminData* d= (NelminData*)p;
  int i;
  for (i=0; i<d->nPar; i++) d->doubleBuf[i]= par[i];
  return (float)(d->sf->value(d->sf,d->doubleBuf,d->nPar));
}

static void nelminReset( float* par, void* p)
{
  NelminData* d= (NelminData*)p;
  int i;
  for (i=0; i<d->nPar; i++) d->doubleBuf[i]= par[i];
  d->sf->reset(d->sf,d->doubleBuf,d->nPar);
}

static void nelminSetTol(Optimizer* self, const double tol)
{
  NelminData* pd= (NelminData*)(self->data);
  baseSetTol(self,tol);
  pd->stopping_val= tol;
}

static double nelminGetTol(Optimizer* self)
{
  NelminData* pd= (NelminData*)(self->data);
  return pd->stopping_val;
}

static void nelminSetScale(Optimizer* self, const double scale)
{
  NelminData* pd= (NelminData*)(self->data);
  baseSetScale(self,scale);
  pd->scale= scale;
}

static double nelminGetScale(Optimizer* self)
{
  NelminData* pd= (NelminData*)(self->data);
  return pd->scale;
}

static int nelminGo(Optimizer* self, ScalarFunction* f, double* par, 
		    const int nPar, double* best)
{
  int numpar;
  int steps_per_conv_check = 3;
  int max_iter = 800;
  int num_iter=0, num_restart=0, return_cond=0;
  float mseval= -1.0;
  NelminData* pd= (NelminData*)(self->data);
  int lcl_npar= nPar;
  int i;
  int retval= 0;
  int max_restart= pd->maxRestarts;

  if (pd->bufLength<nPar) {
    if (pd->startBuf) free(pd->startBuf);
    if (pd->valBuf) free(pd->valBuf);
    if (!(pd->startBuf=(float*)malloc(nPar*sizeof(float))))
      Abort("nelminGo: cannot allocate %d bytes!\n",nPar*sizeof(float));
    if (!(pd->valBuf=(float*)malloc(nPar*sizeof(float))))
      Abort("nelminGo: cannot allocate %d bytes!\n",nPar*sizeof(float));
    if (!(pd->scaleBuf=(float*)malloc(nPar*sizeof(float))))
      Abort("nelminGo: cannot allocate %d bytes!\n",nPar*sizeof(float));
    if (!(pd->doubleBuf=(double*)malloc(nPar*sizeof(double))))
      Abort("nelminGo: cannot allocate %d bytes!\n",nPar*sizeof(double));
    pd->bufLength= nPar;
  }

  for (i=0; i<nPar; i++) {
    pd->startBuf[i]= pd->valBuf[i]= (float)par[i];
    pd->scaleBuf[i]= pd->scale;
  }
  pd->nPar= nPar;
  pd->sf= f;

  nelmin( nelminValue, nelminReset,
	  &lcl_npar, pd->startBuf, pd->valBuf, &mseval,
	  &(pd->stopping_val), pd->scaleBuf, &steps_per_conv_check,
	  &max_iter, &num_iter, &num_restart, &return_cond, &max_restart,
	  pd );
  if (return_cond != 0) {
    if (self->debugLevel) 
      Message("Nelmin optimization failed on return code %d\n",return_cond);
    if (return_cond==2)
      Warning(1,"nelminGo: Nelder-Mead convergence failed!\n");
    else 
      Warning(1,"nelminGo: Nelder-Mead call illegal value!\n");
    Warning(1,"         failed at iter = %d of %d, restart= %d of %d\n",
	    num_iter, max_iter, num_restart, max_restart);
    retval= 0;
  }
  else retval= 1;

  *best= (double)mseval;
  for (i=0; i<nPar; i++) par[i]= (double)pd->valBuf[i];
  return retval;
}

Optimizer* createNelminOptimizer(double stopping_val, double scale,
				 int maxRestarts)
{
  Optimizer* result= createBaseOptimizer();
  NelminData* pd= NULL;

  if (result->data) free(result->data);
  if (!(pd=(NelminData*)malloc(sizeof(NelminData))))
    Abort("createNelminOptimizer: unable to allocate %d bytes!\n");
  result->data= pd;
  pd->startBuf= NULL;
  pd->valBuf= NULL;
  pd->scaleBuf= NULL;
  pd->bufLength= 0;
  pd->nPar= 0;
  pd->stopping_val= (float)stopping_val;
  pd->scale= (float)scale;
  pd->maxRestarts= maxRestarts;

  result->getMethodName= nelminGetMethodName;
  result->getStringRep= nelminGetStringRep;
  result->setTol= nelminSetTol;
  result->getTol= nelminGetTol;
  result->setScale= nelminSetScale;
  result->getScale= nelminGetScale;
  result->go= nelminGo;
  
  return result;
}

typedef struct NelminTData {
  float* startBuf;
  float* valBuf;
  float* scaleBuf;
  double* doubleBuf;
  int bufLength;
  int nPar;
  float scale;
  float stopping_val;
  float Tval;
  int maxRestarts;
  ScalarFunction* sf;
} NelminTData;

static const char* nelminTGetMethodName(Optimizer* self)
{
  static const char* name= "nelmin_t";
  return name;
}

static char* nelminTGetStringRep(Optimizer* self)
{
  NelminTData* pd= (NelminTData*)(self->data);
  char buf[256];
  snprintf(buf,sizeof(buf),
	   "NelminTOptimizer(%g,%g,%g,%d)",pd->stopping_val,pd->scale,
	   pd->Tval,pd->maxRestarts);
  return strdup(buf);
}

static float nelminTValue( float* par, void* p)
{
  NelminTData* d= (NelminTData*)p;
  int i;
  for (i=0; i<d->nPar; i++) d->doubleBuf[i]= par[i];
  return (float)(d->sf->value(d->sf,d->doubleBuf,d->nPar));
}

static void nelminTReset( float* par, void* p)
{
  NelminTData* d= (NelminTData*)p;
  int i;
  for (i=0; i<d->nPar; i++) d->doubleBuf[i]= par[i];
  d->sf->reset(d->sf,d->doubleBuf,d->nPar);
}

static void nelminTSetTol(Optimizer* self, const double tol)
{
  NelminTData* pd= (NelminTData*)(self->data);
  baseSetTol(self,tol);
  pd->stopping_val= tol;
}

static double nelminTGetTol(Optimizer* self)
{
  NelminTData* pd= (NelminTData*)(self->data);
  return pd->stopping_val;
}

static void nelminTSetScale(Optimizer* self, const double scale)
{
  NelminTData* pd= (NelminTData*)(self->data);
  baseSetScale(self,scale);
  pd->scale= scale;
}

static double nelminTGetScale(Optimizer* self)
{
  NelminTData* pd= (NelminTData*)(self->data);
  return pd->scale;
}

static int nelminTGo(Optimizer* self, ScalarFunction* f, double* par, 
		    const int nPar, double* best)
{
  int numpar;
  int steps_per_conv_check = 3;
  int max_iter = 800;
  int num_iter=0, num_restart=0, return_cond=0;
  float mseval= -1.0;
  NelminTData* pd= (NelminTData*)(self->data);
  int lcl_npar= nPar;
  int i;
  int retval= 0;
  float tCritSqr= pd->Tval*pd->Tval;
  int max_restart= pd->maxRestarts;

  if (pd->bufLength<nPar) {
    if (pd->startBuf) free(pd->startBuf);
    if (pd->valBuf) free(pd->valBuf);
    if (!(pd->startBuf=(float*)malloc(nPar*sizeof(float))))
      Abort("nelminTGo: cannot allocate %d bytes!\n",nPar*sizeof(float));
    if (!(pd->valBuf=(float*)malloc(nPar*sizeof(float))))
      Abort("nelminTGo: cannot allocate %d bytes!\n",nPar*sizeof(float));
    if (!(pd->scaleBuf=(float*)malloc(nPar*sizeof(float))))
      Abort("nelminTGo: cannot allocate %d bytes!\n",nPar*sizeof(float));
    if (!(pd->doubleBuf=(double*)malloc(nPar*sizeof(double))))
      Abort("nelminTGo: cannot allocate %d bytes!\n",nPar*sizeof(double));
    pd->bufLength= nPar;
  }

  for (i=0; i<nPar; i++) {
    pd->startBuf[i]= pd->valBuf[i]= (float)par[i];
    pd->scaleBuf[i]= pd->scale;
  }
  pd->nPar= nPar;
  pd->sf= f;

  nelmin_t( nelminTValue, nelminTReset,
	  &lcl_npar, pd->startBuf, pd->valBuf, &mseval,
	  &(pd->stopping_val), pd->scaleBuf, &steps_per_conv_check,
	  &max_iter, &num_iter, &num_restart, &return_cond, &max_restart,
	  &tCritSqr, pd );
  if (return_cond != 0) {
    if (self->debugLevel) 
      Message("NelminT optimization failed on return code %d\n",return_cond);
    if (return_cond==2)
      Warning(1,"nelminTGo: Nelder-Mead convergence failed!\n");
    else 
      Warning(1,"nelminTGo: Nelder-Mead call illegal value!\n");
    Warning(1,"         failed at iter = %d of %d, restart= %d of %d\n",
	    num_iter, max_iter, num_restart, max_restart);
    retval= 0;
  }
  else retval= 1;

  *best= (double)mseval;
  for (i=0; i<nPar; i++) par[i]= (double)pd->valBuf[i];
  return retval;
}

Optimizer* createNelminTOptimizer(double stopping_val, double scale, 
				  double Tval, int maxRestarts)
{
  Optimizer* result= createBaseOptimizer();
  NelminTData* pd= NULL;

  if (result->data) free(result->data);
  if (!(pd=(NelminTData*)malloc(sizeof(NelminTData))))
    Abort("createNelminTOptimizer: unable to allocate %d bytes!\n");
  result->data= pd;
  pd->startBuf= NULL;
  pd->valBuf= NULL;
  pd->scaleBuf= NULL;
  pd->bufLength= 0;
  pd->nPar= 0;
  pd->stopping_val= (float)stopping_val;
  pd->scale= (float)scale;
  pd->Tval= Tval;
  pd->maxRestarts= maxRestarts;

  result->getMethodName= nelminTGetMethodName;
  result->getStringRep= nelminTGetStringRep;
  result->setTol= nelminTSetTol;
  result->getTol= nelminTGetTol;
  result->setScale= nelminTSetScale;
  result->getScale= nelminTGetScale;
  result->go= nelminTGo;
  
  return result;
}

Optimizer* optimizerFromStringRep( const char* rep )
{
  char* args= strchr(rep,'(');
  if ( !args || (args-rep)<1 ) return NULL;
  if (!strncmp(rep,"BaseOptimizer",strlen("BaseOptimizer"))) {
    return createBaseOptimizer();
  }
  else if (!strncmp(rep,"NoneOptimizer",strlen("NoneOptimizer"))) {
    return createNoneOptimizer();
  }
  else if (!strncmp(rep,"PraxisOptimizer",strlen("PraxisOptimizer"))) {
    double t0;
    double h0;
    int n= sscanf(args,"(%lg,%lg)",&t0,&h0);
    if (n==2) return createPraxisOptimizer(t0,h0);
    else return NULL;
  }
  else if (!strncmp(rep,"NelminOptimizer",strlen("NelminOptimizer"))) {
    float stopping_val;
    float scale;
    int maxRestarts;
    int n= sscanf(args,"(%g,%g,%d)",&stopping_val,&scale,&maxRestarts);
    if (n==3) return createNelminOptimizer(stopping_val,scale,maxRestarts);
    else return NULL;
  }
  else if (!strncmp(rep,"NelminTOptimizer",strlen("NelminTOptimizer"))) {
    float stopping_val;
    float scale;
    float Tval;
    int maxRestarts;
    int n= sscanf(args,"(%g,%g,%g,%d)",&stopping_val,&scale,&Tval,
		  &maxRestarts);
    if (n==3) return createNelminTOptimizer(stopping_val,scale,Tval,
					    maxRestarts);
    else return NULL;
  }
  else return NULL;
}


