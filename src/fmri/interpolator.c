/************************************************************
 *                                                          *
 *  interpolator.c                                          *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *                                                          *
 *  Copyright (c) 2005 Department of Statistics             *
 *                     Carnegie Mellon University           *
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
 *  Original programming by Joel Welling 2/2005             *
 ************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <assert.h>
#include "mri.h"
#include "fmri.h"
#include "interpolator.h"

static char rcsid[] = "$Id: interpolator.c,v 1.17 2008/02/12 01:04:58 welling Exp $";

/********************
 * Notes-
 * -Check that data range is respected?  A few asserts would do it.
 * -Distingush calcLinear2D() from calcClosest2D() (and similarly
 *  for 3D), and do the Closest case more efficiently?
 * -Direct calls to prep2D and prep1D would be faster than indirection
 *  through class methods.
 * -Spline interpolator will fail on creation if dim<4.  Work around that
 *  by substituting appropriate linear interpolators?
 ********************/

/* These names must correspond to opts, in order. */
static char* optNameTable[]= {"DEBUG","FASTBLK","EXTENT","TENSION"};
#define OPT_NAME_TABLE_ENTRIES (sizeof(optNameTable)/sizeof(char*))

/* These names must correspond to types, in order. */
static char* typeNameTable[]= 
  {"CLOSEST","LINEAR","CATMULLROM","BEZIER","BSPLINE", "UNKNOWN"};
#define TYPE_NAME_TABLE_ENTRIES (sizeof(typeNameTable)/sizeof(char*))

#define OFFSET(nx,ny,nz,x,y,z) ((((z)*ny)+(y))*nx + (x))

/* Tolerance for coord calculations (allow for numerical noise) */
#define EPSILON 0.000001

typedef struct Info2D_struct {
  Interpolator* interpX;
  Interpolator* interpY;
  long nx;
  long ny;
  double* buf;
  long bufSize;
} Info2D;

typedef struct Info3D_struct {
  Interpolator* interpX;
  Interpolator* interpYZ;
  long nx;
  long ny;
  long nz;
  double* buf;
  long bufSize;
} Info3D;

static long warpCount= 0;

void intrp_warpClearCounts()
{
  warpCount= 0;
}

void intrp_warpGetCounts( long* count )
{
  *count= warpCount;
}

void intrp_warpApply( Interpolator* interp, Transform t_in,
		      double* orig_image,
		      double* moved_image,
		      char* check,
		      long nx, long ny, long nz, long fast_blk,
		      double length_x, double length_y, double length_z )
{
  Vec4 pout; /* point in output space (grid-aligned) */
  Vec4 p; /* point in input space */
  long iout, jout, kout;
  long halfx= nx/2;
  long halfy= ny/2;
  long halfz= nz/2;
  Transform t;
  Transform tTemp;

#ifdef never
  double xmin= 1000.0;
  double xmax= -1000.0;
  double ymin= 1000.0;
  double ymax= -1000.0;
  double zmin= 1000.0;
  double zmax= -1000.0;
#endif

  /* Step counter */
  warpCount++;

  /* Make a version of the transform in which the voxel scale
   * factors have been sucked into the matrix.  Signs of terms 
   * are determined by relationship between grid and 3D (radiological) 
   * coords.
   */
  trans_identity(tTemp);
  tTemp[0]= length_x/(double)nx;
  tTemp[5]= -length_y/(double)ny;
  tTemp[10]= length_z/(double)nz;
  trans_copy(t, t_in);
  trans_mult_right(t,tTemp);
  tTemp[0]= 1.0/tTemp[0];
  tTemp[5]= 1.0/tTemp[5];
  tTemp[10]= 1.0/tTemp[10];
  trans_mult_left(tTemp,t);

  /* Just do it */
  interp->prep(interp, orig_image, nx*ny*nz*fast_blk);
  pout[3]= 1.0;
  for (kout=0; kout<nz; kout++) {
    pout[2]= (double)(kout-halfz); 
    for (jout=0; jout<ny; jout++) {
      pout[1]= (double)(jout-halfy);
      for (iout=0; iout<nx; iout++) {
	pout[0]= (double)(iout-halfx);
	bcopy(pout,p,sizeof(pout));
	trans_vec_mult(t, p);
	if (p[3]!=1.0) {
	  if (p[3]==0.0) {
	    Abort("interpolator.c: encountered singular projective transform!\n");
	  }
	  else {
	    p[0] /= p[3];
	    p[1] /= p[3];
	    p[2] /= p[3];
	    p[3] == 1.0;
	  }
	}
	p[0] += halfx;
	p[1] += halfy;
	p[2] += halfz;
#ifdef never
	if (p[0]<xmin) xmin= p[0];
	if (p[0]>xmax) xmax= p[0];
	if (p[1]<ymin) ymin= p[1];
	if (p[1]>ymax) ymax= p[1];
	if (p[2]<zmin) zmin= p[2];
	if (p[2]>zmax) zmax= p[2];
#endif
	if ((p[0]<(-EPSILON)) || (p[0]>(double)(nx-1)+EPSILON)
	    || (p[1]<(-EPSILON)) || (p[1]>(double)(ny-1)+EPSILON)
	    || (p[2]<(-EPSILON)) || (p[2]>(double)(nz-1)+EPSILON)) {
	  
	  /* This point maps outside the input array */
	  moved_image[OFFSET(nx,ny,nz,iout,jout,kout)]= 0.0;
	  check[OFFSET(nx,ny,nz,iout,jout,kout)]= 0; /* mark out */
	}
	else {
	  interp->calc(interp,
		       moved_image+OFFSET(nx,ny,nz,iout,jout,kout),
		       p,fast_blk,0);
	  
	  check[OFFSET(nx,ny,nz,iout,jout,kout)]= 1; /* mark in */
	} 
      }
    }
  }
#ifdef never
  fprintf(stderr,"mapping %g<x<%g, %g<y<%g, %g<z<%g\n",
	  xmin,xmax,ymin,ymax,zmin,zmax);
#endif
}


InterpolatorType intrp_typeFromName( const char* name )
{
  int i;
  for (i=0; i<TYPE_NAME_TABLE_ENTRIES; i++)
    if (!strcasecmp(name,typeNameTable[i])) return (InterpolatorType)i;
  return INTRP_UNKNOWN;
}

const char* intrp_nameFromType( InterpolatorType type )
{
  return typeNameTable[(int)type];
}

Interpolator* 
intrp_createInterpolator1DByType( InterpolatorType type, 
				  long nx, long fast_blk )
{
  switch (type) {
  case INTRP_CLOSEST:
    return intrp_createClosestInterpolator1D( nx, fast_blk );
  case INTRP_LINEAR:
    return intrp_createLinearInterpolator1D( nx, fast_blk );
  case INTRP_CATMULLROM:
    return intrp_createSplineInterpolator1D( SPL_CATMULLROM, nx, fast_blk );
  case INTRP_BEZIER:
    return intrp_createSplineInterpolator1D( SPL_BEZIER, nx, fast_blk );
  case INTRP_BSPLINE:
    return intrp_createSplineInterpolator1D( SPL_BSPLINE, nx, fast_blk );
  default:
    Abort("intrp_createInterpolator1DByType: can't create an interpolator of unknown type!\n");
  }
  return NULL; /* not reached */
}

Interpolator* 
intrp_createInterpolator2DByType( InterpolatorType type, 
				   long nx, long ny, long fast_blk )
{
  switch (type) {
  case INTRP_CLOSEST:
    return intrp_createClosestInterpolator2D( nx, ny, fast_blk );
  case INTRP_LINEAR:
    return intrp_createLinearInterpolator2D( nx, ny, fast_blk );
  case INTRP_CATMULLROM:
    return intrp_createSplineInterpolator2D( SPL_CATMULLROM, nx, ny, 
					     fast_blk );
  case INTRP_BEZIER:
    return intrp_createSplineInterpolator2D( SPL_BEZIER, nx, ny, fast_blk );
  case INTRP_BSPLINE:
    return intrp_createSplineInterpolator2D( SPL_BSPLINE, nx, ny, fast_blk );
  default:
    Abort("intrp_createInterpolator1DByType: can't create an interpolator of unknown type!\n");
  }
  return NULL; /* not reached */
}

Interpolator* 
intrp_createInterpolator3DByType( InterpolatorType type, 
				  long nx, long ny, long nz, long fast_blk )
{
  switch (type) {
  case INTRP_CLOSEST:
    return intrp_createClosestInterpolator3D( nx, ny, nz, fast_blk );
  case INTRP_LINEAR:
    return intrp_createLinearInterpolator3D( nx, ny, nz, fast_blk );
  case INTRP_CATMULLROM:
    return intrp_createSplineInterpolator3D( SPL_CATMULLROM, nx, ny, nz,
					     fast_blk );
  case INTRP_BEZIER:
    return intrp_createSplineInterpolator3D( SPL_BEZIER, nx, ny, nz, 
					     fast_blk );
  case INTRP_BSPLINE:
    return intrp_createSplineInterpolator3D( SPL_BSPLINE, nx, ny, nz, 
					     fast_blk );
  default:
    Abort("intrp_createInterpolator1DByType: can't create an interpolator of unknown type!\n");
  }
  return NULL; /* not reached */
}					       

static void baseDestroySelf( Interpolator* self )
{
  if (self->debug) 
    fprintf(stderr,"Destroying <%s> interpolator at 0x%lx\n",
	    self->typeName,(long)self);
  if (self->hook) free(self->hook);
  free(self);
}

static void baseDumpSelf( const Interpolator* self, FILE* ofile )
{
  fprintf(ofile,
	  "<%s> interpolator at 0x%lx, fast blk %ld, extent %d, data 0x%lx\n",
	  self->typeName,(long)self,self->fast_blksize,self->extent,
	  (long)self->dataField);
}

static void baseSetInt( Interpolator* self, int which, long val )
{
  if (which<0 || which>OPT_NAME_TABLE_ENTRIES) 
    Abort("interpolator: baseSetInt: unknown option %d!\n",which);
  else {
    switch (which) {
    case INTRP_OPT_DEBUG: self->debug= val; break;
    case INTRP_OPT_FASTBLK: self->fast_blksize= val; break;
    case INTRP_OPT_EXTENT: self->extent= val; break;
    default: 
      Abort("interpolator:baseSetInt: option %s is not an int!\n",
	    optNameTable[which]);
    }
  }
}

static void baseSetDouble( Interpolator* self, int which, double val )
{
  if (which<0 || which>OPT_NAME_TABLE_ENTRIES) 
    Abort("interpolator: baseSetInt: unknown option %d!\n",which);
  else {
    switch (which) {
    case INTRP_OPT_TENSION: break; /* derived class should handle it */
    default: 
      Abort("interpolator:baseSetDouble: option %s is not a double!\n",
	    optNameTable[which]);
    }
  }
}

static long baseGetInt( const Interpolator* self, int which )
{
  if (which<0 || which>OPT_NAME_TABLE_ENTRIES) 
    Abort("interpolator: baseSetInt: unknown option %d!\n",which);
  else {
    switch (which) {
    case INTRP_OPT_DEBUG: return self->debug;
    case INTRP_OPT_FASTBLK: return self->fast_blksize;
    case INTRP_OPT_EXTENT: return self->extent;
    default: 
      Abort("interpolator:baseGetInt: option %s is not an int!\n",
	    optNameTable[which]);
    }
  }
  return 0; /* not reached */
}

static double baseGetDouble( const Interpolator* self, int which )
{
  if (which<0 || which>OPT_NAME_TABLE_ENTRIES) 
    Abort("interpolator: baseSetInt: unknown option %d!\n",which);
  else {
    switch (which) {
    case INTRP_OPT_TENSION: return 0.0; /* derived class should have handled */
    default: 
      Abort("interpolator:baseGetDouble: option %s is not a double!\n",
	    optNameTable[which]);
    }
  }
  return 0.0; /* not reached */
}

static void basePrep( Interpolator* self, double* data, long sz )
{
  if (self->debug)
    fprintf(stderr,
	    "interpolator:basePrep: type %s at 0x%lx, %ld doubles starting at 0x%lx\n",
	    self->typeName,(long)self,sz,(long)data);
#ifdef never
  {
    long i;
    fprintf(stderr,"Scanning...");
    for (i=0; i<sz; i++) 
      if (data[i]<-10000.0 || data[i]>10000.0) 
	Abort("Outside range at %g!\n",data[i]);
    fprintf(stderr,"done\n");
  }
#endif
  if (sz < self->dataFieldLength)
    Abort("interpolator:basePrep: size mismatch: %ld vs. %ld!\n",sz,
	  self->dataFieldLength);
  self->dataField= data;
}

static void baseCalc( Interpolator* self, double* result, double* loc,
		      long runLength, long offset )
{
  /* Assume *loc is 0.0 and just copy input to output.  This is
   * actually useful as a fallback case when extent==1.
   */
  long i;
  for (i=0; i<runLength; i++) result[i]= self->dataField[i+offset];
}

static Interpolator* createBaseInterpolator( int extent, int fast_blksize )
{
  Interpolator* result= NULL;
  if (!(result=(Interpolator*)malloc(sizeof(Interpolator))))
    Abort("createBaseInterpolator: unable to allocate %d bytes!\n",
	  sizeof(Interpolator));

  result->prep= basePrep;
  result->calc= baseCalc;
  result->destroySelf= baseDestroySelf;
  result->dumpSelf= baseDumpSelf;
  result->setInt= baseSetInt;
  result->setDouble= baseSetDouble;
  result->getInt= baseGetInt;
  result->getDouble= baseGetDouble;
  result->typeName= "base";

  result->dataField= NULL;
  result->dataFieldLength= fast_blksize*extent;
  result->fast_blksize= fast_blksize;
  result->extent= extent;
  result->debug= 0;
  result->hook= NULL;

  return result;
}

static void calcClosest1D( Interpolator* self, double* result, double* loc,
			   long runLength, long offset )
{
  int i;
  int iloc= (int)rint(*loc);
  if (iloc<0) iloc= 0;
  if (iloc>self->extent-1) iloc= self->extent-1;
  if (self->debug)
    fprintf(stderr,"interpolator:calcClosest1D: Doing %d at %d; loc= %f\n",
	    runLength, iloc*self->fast_blksize+offset,*loc);
  for (i=0; i<runLength; i++)
    result[i]= self->dataField[iloc*self->fast_blksize + offset + i];
}

Interpolator* intrp_createClosestInterpolator1D( long nx, long fast_blksize )
{
  Interpolator* result= createBaseInterpolator( nx, fast_blksize );
  result->typeName= "closest1D";
  result->calc= calcClosest1D;
  return result;
}

static void setInt2D( Interpolator* self, int which, long val )
{
  Info2D* info= (Info2D*)self->hook;
  switch (which) {
  case INTRP_OPT_EXTENT:
    Abort("interpolator.c:setInt2D: cannot set EXTENT on 2D interpolators.\n");
  case INTRP_OPT_DEBUG:
    /* Don't pass debug to children; it makes too much output */
    baseSetInt(self,which,val);
    break; 
  default: 
    {
      info->interpX->setInt(info->interpX, which, val);
      info->interpY->setInt(info->interpY, which, val);
      baseSetInt( self, which, val );
    }
    break;
  }
}

static long getInt2D( const Interpolator* self, int which )
{
  Info2D* info= (Info2D*)self->hook;
  switch (which) {
  case INTRP_OPT_DEBUG: return self->debug;
  default: return info->interpX->getInt(info->interpX, which);
  }
}

static void setDouble2D( Interpolator* self, int which, double val )
{
  Info2D* info= (Info2D*)self->hook;
  info->interpX->setDouble(info->interpX, which, val);
  info->interpY->setDouble(info->interpY, which, val);
  baseSetDouble( self, which, val );
}

static double getDouble2D( const Interpolator* self, int which )
{
  Info2D* info= (Info2D*)self->hook;
  return info->interpX->getDouble(info->interpX, which);
}

static void destroySelf2D( Interpolator* self )
{
  Info2D* info= (Info2D*)self->hook;
  info->interpX->destroySelf(info->interpX);
  info->interpY->destroySelf(info->interpY);
  if (info->buf) free(info->buf);
  baseDestroySelf(self);
}

static void dumpSelf2D( const Interpolator* self, FILE* ofile )
{
  Info2D* info= (Info2D*)self->hook;
  baseDumpSelf(self,ofile);
  fprintf(ofile,"Sub-interpolators follow:\n");
  info->interpX->dumpSelf(info->interpX, ofile);
  info->interpY->dumpSelf(info->interpY, ofile);
}

static void prep2D( Interpolator* self, double* data, long sz )
{
  Info2D* info= (Info2D*)self->hook;
  long neededBufSize= info->nx*self->fast_blksize;
  basePrep(self, data, sz);
  if (info->bufSize != neededBufSize) {
    if (info->buf) free(info->buf);
    if (!(info->buf=(double*)malloc(neededBufSize*sizeof(double))))
      Abort("%s:%d: unable to allocate %d bytes!\n",__FILE__,__LINE__,
	    neededBufSize*sizeof(double));
    info->bufSize= neededBufSize;
    info->interpX->prep( info->interpX, info->buf, info->bufSize );
  }
  info->interpY->prep( info->interpY, data, sz );
}

static void calc2D( Interpolator* self, double* result, double* loc,
		    long runLength, long offset )
{
  Info2D* info= (Info2D*)self->hook;
  long leftXLoc= (long)floor(loc[0]);
  long rightXLoc= (long)ceil(loc[0]);
  if (leftXLoc<0) {
    leftXLoc= 0;
    if (rightXLoc<0) rightXLoc= 0;
  }
  if (rightXLoc>info->nx-1) {
    rightXLoc= info->nx-1;
    if (leftXLoc>info->nx-1) leftXLoc= info->nx-1;
  }
  if (self->debug)
    fprintf(stderr,
	    "interpolator:calc2D: %ld values, offset %ld, loc= %f %f\n",
	    runLength, offset, loc[0], loc[1]);
  if ( leftXLoc==rightXLoc ) {
    info->interpY->calc( info->interpY, result, loc+1, runLength, 
			 offset+self->fast_blksize*leftXLoc );
  }
  else {
    info->interpY->calc( info->interpY,
			 info->buf+leftXLoc*self->fast_blksize, loc+1, 
			 2*self->fast_blksize,
			 leftXLoc*self->fast_blksize );
    info->interpX->calc( info->interpX, result, loc, runLength, offset );
  }
}

Interpolator* intrp_createClosestInterpolator2D( long nx, long ny, 
						 long fast_blksize )
{
  Interpolator* result= NULL;

  if (ny==1) {
    result= intrp_createClosestInterpolator1D( nx, fast_blksize );
  }
  else {
    Info2D* info= NULL;
    result= createBaseInterpolator( nx*ny, fast_blksize );
    
    if (!(info=(Info2D*)malloc(sizeof(Info2D))))
      Abort("%s:%d: unable to allocate %d bytes!\n",
	    __FILE__,__LINE__,sizeof(Info2D));
    result->hook= info;
    
    info->nx= nx;
    info->ny= ny;
    if (nx==0)
      info->interpX= createBaseInterpolator( nx, fast_blksize );
    else
      info->interpX= intrp_createClosestInterpolator1D( nx, fast_blksize );
    info->interpY= intrp_createClosestInterpolator1D( ny, nx*fast_blksize );
    info->buf= NULL;
    info->bufSize= 0;
    result->prep= prep2D;
    result->calc= calc2D;
    result->destroySelf= destroySelf2D;
    result->dumpSelf= dumpSelf2D;
    result->setInt= setInt2D;
    result->setDouble= setDouble2D;
    result->getInt= getInt2D;
    result->getDouble= getDouble2D;
    result->typeName= "closest2D";
  }

  return result;
}

static void setInt3D( Interpolator* self, int which, long val )
{
  Info3D* info= (Info3D*)self->hook;
  switch (which) {
  case INTRP_OPT_EXTENT:
    Abort("interpolator.c:setInt3D: cannot set EXTENT on 3D interpolators.\n");
  case INTRP_OPT_DEBUG:
    /* Don't pass debug to children; it makes too much output */
    baseSetInt(self,which,val);
    break; 
  default: 
    {
      info->interpYZ->setInt(info->interpYZ, which, val);
      info->interpX->setInt(info->interpX, which, val);
      baseSetInt( self, which, val );
    }
    break;
  }
}

static long getInt3D( const Interpolator* self, int which )
{
  Info3D* info= (Info3D*)self->hook;
  switch (which) {
  case INTRP_OPT_DEBUG: return self->debug;
  default: return info->interpYZ->getInt(info->interpYZ, which);
  }
}

static void setDouble3D( Interpolator* self, int which, double val )
{
  Info3D* info= (Info3D*)self->hook;
  info->interpYZ->setDouble(info->interpYZ, which, val);
  info->interpX->setDouble(info->interpX, which, val);
  baseSetDouble( self, which, val );
}

static double getDouble3D( const Interpolator* self, int which )
{
  Info3D* info= (Info3D*)self->hook;
  return info->interpYZ->getDouble(info->interpYZ, which);
}

static void destroySelf3D( Interpolator* self )
{
  Info3D* info= (Info3D*)self->hook;
  info->interpYZ->destroySelf(info->interpYZ);
  info->interpX->destroySelf(info->interpX);
  if (info->buf) free(info->buf);
  baseDestroySelf(self);
}

static void dumpSelf3D( const Interpolator* self, FILE* ofile )
{
  Info3D* info= (Info3D*)self->hook;
  baseDumpSelf(self,ofile);
  fprintf(ofile,"Sub-interpolators follow:\n");
  info->interpYZ->dumpSelf(info->interpYZ, ofile);
  info->interpX->dumpSelf(info->interpX, ofile);
}

static void prep3D( Interpolator* self, double* data, long sz )
{
  Info3D* info= (Info3D*)self->hook;
  long neededBufSize= info->nx*self->fast_blksize;
  basePrep(self, data, sz);
  if (info->bufSize != neededBufSize) {
    if (info->buf) free(info->buf);
    if (!(info->buf=(double*)malloc(neededBufSize*sizeof(double))))
      Abort("%s:%d: unable to allocate %d bytes!\n",__FILE__,__LINE__,
	    neededBufSize*sizeof(double));
    info->bufSize= neededBufSize;
    info->interpX->prep( info->interpX, info->buf, info->bufSize );
  }
  info->interpYZ->prep( info->interpYZ, data, sz );
}

static void calc3D( Interpolator* self, double* result, double* loc,
		    long runLength, long offset )
{
  Info3D* info= (Info3D*)self->hook;
  long leftXLoc= (long)floor(loc[0]);
  long rightXLoc= (long)ceil(loc[0]);
  if (leftXLoc<0) {
    leftXLoc= 0;
    if (rightXLoc<0) rightXLoc= 0;
  }
  if (rightXLoc>info->nx-1) {
    rightXLoc= info->nx-1;
    if (leftXLoc>info->nx-1) leftXLoc= info->nx-1;
  }

  if (self->debug)
    fprintf(stderr,
	    "interpolator:calc3D: %ld values, offset %ld, loc= %f %f %f\n",
	    runLength, offset, loc[0], loc[1], loc[2]);
  if (leftXLoc==rightXLoc) {
    info->interpYZ->calc( info->interpYZ, result, loc+1, runLength,
			  offset+self->fast_blksize*leftXLoc );
  }
  else {
    info->interpYZ->calc( info->interpYZ, 
			  info->buf+leftXLoc*self->fast_blksize, loc+1, 
			  2*self->fast_blksize,
			  leftXLoc*self->fast_blksize );
    info->interpX->calc( info->interpX, result, loc, runLength, offset );
  }
}

Interpolator* intrp_createClosestInterpolator3D( long nx, long ny, long nz, 
						 long fast_blksize )
{
  Interpolator* result= NULL;

  if (nz==1) {
    result= intrp_createClosestInterpolator2D(nx,ny,fast_blksize);
  }
  else {
    Info3D* info= NULL;
    result= createBaseInterpolator( nx*ny*nz, fast_blksize );
    
    if (!(info=(Info3D*)malloc(sizeof(Info3D))))
      Abort("%s:%d: unable to allocate %d bytes!\n",
	    __FILE__,__LINE__,sizeof(Info3D));
    result->hook= info;
    
    info->nx= nx;
    info->ny= ny;
    info->nz= nz;

    if (nx==1)
      info->interpX= createBaseInterpolator( nx, fast_blksize );
    else
      info->interpX= intrp_createClosestInterpolator1D( nx, fast_blksize );

    info->interpYZ= intrp_createClosestInterpolator2D( ny, nz, 
						       nx*fast_blksize );
    info->buf= NULL;
    info->bufSize= 0;
    result->prep= prep3D;
    result->calc= calc3D;
    result->destroySelf= destroySelf3D;
    result->dumpSelf= dumpSelf3D;
    result->setInt= setInt3D;
    result->setDouble= setDouble3D;
    result->getInt= getInt3D;
    result->getDouble= getDouble3D;
    result->typeName= "closest3D";
  }

  return result;
}

static void calcLinear1D( Interpolator* self, double* result, double* loc,
			  long runLength, long offset )
{
  long i;
  long leftLoc= (long)floor(*loc);
  double delta= *loc-leftLoc;
  long rightLoc;
  if (leftLoc<0) {
    leftLoc= 0;
    delta= 0.0;
  }
  rightLoc= leftLoc+1;
  if (rightLoc>self->extent-1)
    rightLoc= self->extent-1;
  if (self->debug)
    fprintf(stderr,
	    "interpolator:calcLinear1D: %ld at %ld; loc= %f->%ld, delta= %f\n",
	    runLength, leftLoc*self->fast_blksize+offset,*loc,leftLoc,delta);
  for (i=0; i<runLength; i++) {
    double v1= self->dataField[leftLoc*self->fast_blksize + offset + i];
    double v2= self->dataField[rightLoc*self->fast_blksize + offset + i];
    result[i]= (1.0-delta)*v1 + delta*v2;
  }
}

Interpolator* intrp_createLinearInterpolator1D( long nx, long fast_blksize )
{
  Interpolator* result= createBaseInterpolator( nx, fast_blksize );
  result->typeName= "linear1D";
  result->calc= calcLinear1D;
  return result;
}

Interpolator* intrp_createLinearInterpolator2D( long nx, long ny, long fast_blksize )
{
  Interpolator* result= NULL;

  if (ny==1) {
    result= intrp_createLinearInterpolator1D( nx, fast_blksize );
  }
  else {
    Info2D* info= NULL;
    result= createBaseInterpolator( nx*ny, fast_blksize );
    
    if (!(info=(Info2D*)malloc(sizeof(Info2D))))
      Abort("%s:%d: unable to allocate %d bytes!\n",
	    __FILE__,__LINE__,sizeof(Info2D));
    result->hook= info;
    
    info->nx= nx;
    info->ny= ny;
    if (nx==1) 
      info->interpX= createBaseInterpolator( nx, fast_blksize );
    else
      info->interpX= intrp_createLinearInterpolator1D( nx, fast_blksize );
    info->interpY= intrp_createLinearInterpolator1D( ny, nx*fast_blksize );
    info->buf= NULL;
    info->bufSize= 0;
    result->prep= prep2D;
    result->calc= calc2D;
    result->destroySelf= destroySelf2D;
    result->dumpSelf= dumpSelf2D;
    result->setInt= setInt2D;
    result->setDouble= setDouble2D;
    result->getInt= getInt2D;
    result->getDouble= getDouble2D;
    result->typeName= "linear2D";
  }

  return result;
}

Interpolator* intrp_createLinearInterpolator3D( long nx, long ny, long nz, 
						long fast_blksize )
{
  Interpolator* result= NULL;

  if (nz==1) {
    result= intrp_createLinearInterpolator2D( nx, ny, fast_blksize );
  }
  else {
    Info3D* info= NULL;
    result= createBaseInterpolator( nx*ny*nz, fast_blksize );
    
    if (!(info=(Info3D*)malloc(sizeof(Info3D))))
      Abort("%s:%d: unable to allocate %d bytes!\n",
	    __FILE__,__LINE__,sizeof(Info3D));
    result->hook= info;
    
    info->nx= nx;
    info->ny= ny;
    info->nz= nz;
    if (nx==1)
      info->interpX= createBaseInterpolator( nx, fast_blksize );
    else
      info->interpX= intrp_createLinearInterpolator1D( nx, fast_blksize );
    info->interpYZ= intrp_createLinearInterpolator2D( ny, nz, 
						      nx*fast_blksize );
    info->buf= NULL;
    info->bufSize= 0;
    result->prep= prep3D;
    result->calc= calc3D;
    result->destroySelf= destroySelf3D;
    result->dumpSelf= dumpSelf3D;
    result->setInt= setInt3D;
    result->setDouble= setDouble3D;
    result->getInt= getInt3D;
    result->getDouble= getDouble3D;
    result->typeName= "linear3D";
  }

  return result;
}

static void splineSetDouble( Interpolator* self, int which, double val )
{
  switch (which) {
  case INTRP_OPT_TENSION: 
    {
      Spline* spl= (Spline*)(self->hook);
      spl_set_tension(spl,val);
    }
    break;
  default: baseSetDouble( self, which, val ); break;
  }
}

static double splineGetDouble( const Interpolator* self, int which )
{
  switch (which) {
  case INTRP_OPT_TENSION: 
    {
      Spline* spl= (Spline*)(self->hook);
      return spl_get_tension(spl);
    }
  default: return baseGetDouble( self, which );
  }
  return 0.0; /* not reached */
}

static void destroySpline1D( Interpolator* self )
{
  if (self->hook) {
    spl_destroy((Spline*)(self->hook));
    self->hook= NULL;
  }
  baseDestroySelf(self);
}

static void dumpSpline1D( const Interpolator* self, FILE* ofile )
{
  baseDumpSelf(self,ofile);
  spl_dump(ofile,(Spline*)(self->hook));
}

static void prepSpline1D( Interpolator* self, double* data, long sz )
{
  Spline* spline= (Spline*)self->hook;
  basePrep(self, data, sz);
  spl_reset(spline, self->fast_blksize, self->extent, data);
}

static void calcSpline1D( Interpolator* self, double* result, double* loc,
			  long runLength, long offset )
{
  Spline* spl= (Spline*)(self->hook);
  if (self->debug)
    fprintf(stderr,
	    "interpolator:calcSpline1D: %ld values, offset %ld, loc= %f\n",
	    runLength, offset, *loc);
  spl_calc( result, spl, runLength, offset, 
	    *loc/((double)(self->extent - 1)) );
}

Interpolator* intrp_createSplineInterpolator1D( SplineType type,
						long nx, long fast_blksize )
{
  Interpolator* result= createBaseInterpolator( nx, fast_blksize );
  result->typeName= "spline1D";
  result->setDouble= splineSetDouble;
  result->getDouble= splineGetDouble;
  result->destroySelf= destroySpline1D;
  result->dumpSelf= dumpSpline1D;
  result->prep= prepSpline1D;
  result->calc= calcSpline1D;
  result->hook= spl_create( fast_blksize, nx, type, NULL );
  return result;
}

static void calcSpline2D( Interpolator* self, double* result, double* loc,
			  long runLength, long offset )
{
  Info2D* info= (Info2D*)self->hook;
  long lvlB= (long)floor(loc[0]);
  long lvlC= (long)ceil(loc[0]);
  long lvlA, lvlD;
  if (lvlB<0) {
    lvlB= 0;
    if (lvlC<0) lvlC= 0;
  }
  if (lvlC>info->nx-1) {
    lvlC= info->nx-1;
    if (lvlB>info->nx-1) lvlB= info->nx-1;
  }
  if (lvlB>0) lvlA= lvlB-1;
  else lvlA= 0;
  if (lvlC<info->nx-1) lvlD= lvlC+1;
  else lvlD= info->nx-1;

  if (self->debug)
    fprintf(stderr,
	    "interpolator:calcSpline2D: %ld values, offset %ld, loc= %f %f\n",
	    runLength, offset, loc[0], loc[1]);

  if (lvlD-lvlA>0) {
    info->interpY->calc( info->interpY, 
			 info->buf+lvlA*self->fast_blksize, loc+1,
			 (lvlD-lvlA)*self->fast_blksize, 
			 lvlA*self->fast_blksize );
    info->interpX->calc( info->interpX, result, loc, runLength, offset );
  }
  else {
    /* All stacked up on one point */
    info->interpY->calc( info->interpY, result, loc+1, runLength, 
			 offset+self->fast_blksize*lvlB );
  }
}

Interpolator* intrp_createSplineInterpolator2D( SplineType type,
						long nx, long ny, 
						long fast_blksize )
{
  Interpolator* result= NULL;

  if (ny==1) {
    result= intrp_createSplineInterpolator1D( type, nx, fast_blksize );
  }
  else {
    Info2D* info= NULL;
    result= createBaseInterpolator( nx*ny, fast_blksize );
    
    if (!(info=(Info2D*)malloc(sizeof(Info2D))))
      Abort("%s:%d: unable to allocate %d bytes!\n",
	    __FILE__,__LINE__,sizeof(Info2D));
    result->hook= info;
    
    info->nx= nx;
    info->ny= ny;
    if (nx==1)
      info->interpX= createBaseInterpolator( nx, fast_blksize );
    else
      info->interpX= intrp_createSplineInterpolator1D( type, nx, 
						       fast_blksize );
    info->interpY= intrp_createSplineInterpolator1D( type, ny, 
						     nx*fast_blksize );
    info->buf= NULL;
    info->bufSize= 0;
    result->prep= prep2D;
    result->calc= calcSpline2D;
    result->destroySelf= destroySelf2D;
    result->dumpSelf= dumpSelf2D;
    result->setInt= setInt2D;
    result->setDouble= setDouble2D;
    result->getInt= getInt2D;
    result->getDouble= getDouble2D;
    result->typeName= "spline2D";
  }

  return result;
}

static void calcSpline3D( Interpolator* self, double* result, double* loc,
		    long runLength, long offset )
{
  Info3D* info= (Info3D*)self->hook;
  long lvlB= (long)floor(loc[0]);
  long lvlC= (long)ceil(loc[0]);
  long lvlA, lvlD;
  if (lvlB<0) {
    lvlB= 0;
    if (lvlC<0) lvlC= 0;
  }
  if (lvlC>info->nx-1) {
    lvlC= info->nx-1;
    if (lvlB>info->nx-1) lvlB= info->nx-1;
  } 
  if (lvlB>0) lvlA= lvlB-1;
  else lvlA= 0;
  if (lvlC<info->nx-1) lvlD= lvlC+1;
  else lvlD= info->nx-1; 

  if (self->debug)
    fprintf(stderr,
	    "interpolator:calcSpline3D: %ld values, offset %ld, loc= %f %f %f\n",
	    runLength, offset, loc[0], loc[1], loc[2]);
  if (lvlD-lvlA>0) {
    info->interpYZ->calc( info->interpYZ, 
			  info->buf+lvlA*self->fast_blksize, loc+1,
			  (lvlD-lvlA)*self->fast_blksize, 
			  lvlA*self->fast_blksize );
    info->interpX->calc( info->interpX, result, loc, runLength, offset );
  }
  else {
    /* All stacked up on one point */
    info->interpYZ->calc( info->interpYZ, result, loc+1, runLength, 
			  offset+self->fast_blksize*lvlB );
  }
}

Interpolator* intrp_createSplineInterpolator3D( SplineType type,
						long nx, long ny, long nz, 
						long fast_blksize )
{
  Interpolator* result= NULL;

  if (nz==1) {
    result= intrp_createSplineInterpolator2D( type, nx, ny, fast_blksize );
  }
  else {
    Info3D* info= NULL;
    result= createBaseInterpolator( nx*ny*nz, fast_blksize );
    
    if (!(info=(Info3D*)malloc(sizeof(Info3D))))
      Abort("%s:%d: unable to allocate %d bytes!\n",
	    __FILE__,__LINE__,sizeof(Info3D));
    result->hook= info;
    
    info->nx= nx;
    info->ny= ny;
    info->nz= nz;
    if (nx==1)
      info->interpX= createBaseInterpolator( nx, fast_blksize );
    else
      info->interpX= intrp_createSplineInterpolator1D( type, nx, 
						       fast_blksize );
    info->interpYZ= intrp_createSplineInterpolator2D( type,
						      ny, nz, 
						      nx*fast_blksize );
    info->buf= NULL;
    info->bufSize= 0;
    result->prep= prep3D;
    result->calc= calcSpline3D;
    result->destroySelf= destroySelf3D;
    result->dumpSelf= dumpSelf3D;
    result->setInt= setInt3D;
    result->setDouble= setDouble3D;
    result->getInt= getInt3D;
    result->getDouble= getDouble3D;
    result->typeName= "spline3D";
  }

  return result;
}

