/************************************************************
 *                                                          *
 *  morph.c                                                *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1995 Department of Statistics,         *
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
 *  Original programming by Mark Fitzgerald  2-95           *
 *     5-96: Pittsburgh Format, Mark Fitzgerald             *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <ctype.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"
#include "slist.h"
#include "misc.h"

static char rcsid[] = "$Id: morph.c,v 1.6 2007/03/21 23:58:25 welling Exp $";

#define DEFAULT_NBR_THRESH 20
#define DEFAULT_MAX_REPS 1

/* Access for 3D arrays */
#define LOC(matrix,x,y,z,dx,dy,dz) matrix[((((z)*dy)+(y))*dx)+(x)]

/* Struct to hold operations */
typedef enum { OP_ERODE, OP_DILATE, OP_LABEL, OP_FLOOD } OpType;
typedef struct step_struct {
  OpType op;
  int maxReps;
  int nbrThresh;
} StepType;

typedef struct kvol_struct {
  FComplex* data;
  long dx;
  long dy;
  long dz;
} KVol;

typedef struct long_stack_struct {
  long* base;
  long* top;
  long size;
} LongStack;

/* Volumes are accessed z-fastest in k-space for compatibility with fft3d */
#define KVOL_OFFSET( v, i, j, k ) ((((i)*v->dy)+(j))*v->dz + (k))
#define KVOL_LOC( v, i, j, k ) (v->data + KVOL_OFFSET(v,i,j,k))
#define KVOL_OFFLOC( v, offset ) (v->data + offset)
#define SET( v, offset ) {KVOL_OFFLOC(v,offset)->real= 1.0;}
#define CLEAR( v, offset ) {KVOL_OFFLOC(v,offset)->real= 0.0;}
#define ISSET( v, offset ) (KVOL_OFFLOC(v,offset)->real>=0.5)
#define ISCLEAR( v, offset ) (KVOL_OFFLOC(v,offset)->real<=0.5)

static int debug_flag= 0;
static int verbose_flag= 0;
static char* progname;
static SList* stepList= NULL;

/* Sizes which are handled well by FFTW */
static long magic[]= { 2,4,8,12,16,24,32,48,64,80,96,128,144,160,192,216,256,
384,512,768,1024,1280,1536,2048};

static const char* getOpName( OpType o )
{
  switch (o) {
  case OP_ERODE: return "erode";
  case OP_DILATE: return "dilate";
  case OP_LABEL: return "label";
  case OP_FLOOD: return "flood";
  }
  return NULL;
}

/* This is like strtok_r, but respects nested parens */
static char* mystrtokr( char* str, const char* delim, char** t )
{
  char* start;
  char* here;
  int parenLvl= 0;

  if (str) start=str;
  else start= *t;
  if (!start) return NULL;
  here= start;
  while (*here) {
    if (*here=='(') parenLvl++;
    else if (*here==')') parenLvl--;
    if (parenLvl==0 && strchr(delim,*here)) {
      /* Break the string here! */
      *here= '\0';
      *t= here+1;
      return start;
    }
    else if (parenLvl<0)
      Abort("%s: bad token <%s>\n",progname,start);
    
    ++here;
  }
  /* If we got here, we're out of string */
  *t= NULL;
  return start;
}

static void kvol_clear( KVol* v )
{
  FComplex* runner= v->data;
  FComplex* limit= runner+v->dx*v->dy*v->dz;
  while (runner<limit) {
    runner->imag= runner->real= 0.0;
    runner++;
  }
}

static KVol* kvol_create( long dx, long dy, long dz )
{
  KVol* result= NULL;
  if (debug_flag) fprintf(stderr,"Allocating KVol( %d, %d, %d )\n",dx,dy,dz);
  if (!(result=(KVol*)malloc(sizeof(KVol))))
    Abort("%s: unable to allocate %d bytes!\n",progname,sizeof(KVol));
  if (!(result->data=(FComplex*)malloc(dx*dy*dz*sizeof(FComplex))))
    Abort("%s: unable to allocate %d bytes!\n",progname,sizeof(KVol));
  result->dx= dx;
  result->dy= dy;
  result->dz= dz;
  return result;
}

static KVol* kvol_create_from_prototype( KVol* orig )
{
  KVol* result= kvol_create(orig->dx,orig->dy,orig->dz);
  return result;
}

static void kvol_destroy( KVol* v )
{
  free(v->data);
  v->data= NULL;
  free(v);
}

static void checkFormat( MRI_Dataset* ds, char* fname )
{
  if( !mri_has( ds, "images" ) )
    Abort( "%s operates only on standard images; %s won't work", 
	   progname, fname );
  if( mri_has( ds, "images.dimensions" ) )
    {
      char* dimstr= mri_get_string(ds,"images.dimensions");
      if (!strcmp( dimstr,"xyz")) {
	/* All is well */
      }
      else if (!strcmp( dimstr,"vxyz")) {
	if (!mri_has(ds,"images.extent.v"))
	  Abort("%s: %s is missing images.extent.t tag!\n",progname,fname);
	if (mri_get_int(ds,"images.extent.v") != 1)
	  Abort("%s: %s has v extent greater than one!\n",progname,fname);
      }
      else if (!strcmp( dimstr,"vxyzt")) {
	if (!mri_has(ds,"images.extent.v"))
	  Abort("%s: %s is missing images.extent.v tag!\n",progname,fname);
	if (mri_get_int(ds,"images.extent.v") != 1)
	  Abort("%s: %s has v extent greater than one!\n",progname,fname);
	if (!mri_has(ds,"images.extent.t"))
	  Abort("%s: %s is missing images.extent.t tag!\n",progname,fname);
      }
      else if (!strcmp( dimstr,"xyzt")) {
	if (!mri_has(ds,"images.extent.t"))
	  Abort("%s: %s is missing images.extent.t tag!\n",progname,fname);
      }
      else
	Abort("%s: input dataset %s must have dimensions xyz!\n",
	      progname,fname);
    }
  else
    Abort( "%s: %s does not have the images.dimensions key.", 
	   progname,fname);
  if (!mri_has(ds,"images.extent.x"))
    Abort("%s: %s is missing images.extent.x tag!\n",progname,fname);
  if (!mri_has(ds,"images.extent.y"))
    Abort("%s: %s is missing images.extent.y tag!\n",progname,fname);
  if (!mri_has(ds,"images.extent.z"))
    Abort("%s: %s is missing images.extent.z tag!\n",progname,fname);
}

static long pad_up( n )
{
  int i;
  for (i=0; i<sizeof(magic)/sizeof(long); i++) {
    if (magic[i]>n) return magic[i];
  }
  Abort("%s: dimension extent of %d is too big!\n",progname,n);
}

static void kvol_mask_copy_in_nonzero( KVol* v, float* buf, 
				       long dx, long dy, long dz)
{
  long kxoff= (v->dx - dx)/2;
  long kyoff= (v->dy - dy)/2;
  long kzoff= (v->dz - dz)/2;
  long i;
  long j;
  long k;

  for (k=0; k<dz; k++)
    for (j=0; j<dy; j++)
      for (i=0; i<dx; i++) {
	FComplex* here= KVOL_LOC(v,i+kxoff,j+kyoff,k+kzoff);
	here->imag= 0.0;
	here->real= (LOC(buf,i,j,k,dx,dy,dz) != 0.0) ? 1.0 : 0.0;
      }

}

static void kvol_mask_copy_out_nonzero( float* buf, KVol* v,
					long dx, long dy, long dz)
{
  long kxoff= (v->dx - dx)/2;
  long kyoff= (v->dy - dy)/2;
  long kzoff= (v->dz - dz)/2;
  long i;
  long j;
  long k;

  for (k=0; k<dz; k++)
    for (j=0; j<dy; j++)
      for (i=0; i<dx; i++) {
	FComplex* here= KVOL_LOC(v,i+kxoff,j+kyoff,k+kzoff);
	LOC(buf,i,j,k,dx,dy,dz)= (here->real != 0) ? 1.0 : 0.0;
      }

}

static void kvol_copy_out_nonzero( float* buf, KVol* v,
				   long dx, long dy, long dz)
{
  long kxoff= (v->dx - dx)/2;
  long kyoff= (v->dy - dy)/2;
  long kzoff= (v->dz - dz)/2;
  long i;
  long j;
  long k;

  for (k=0; k<dz; k++)
    for (j=0; j<dy; j++)
      for (i=0; i<dx; i++) {
	FComplex* here= KVOL_LOC(v,i+kxoff,j+kyoff,k+kzoff);
	LOC(buf,i,j,k,dx,dy,dz)= here->real;
      }

}

static void mark_neighborhood( KVol* v )
{
  long xCtr= v->dx/2;
  long yCtr= v->dy/2;
  long zCtr= v->dz/2;
  long i,j,k;
  FComplex save;

  save= *KVOL_LOC(v,xCtr,yCtr,zCtr);
  for (k=0; k<3; k++)
    for (j=0; j<3; j++)
      for (i=0; i<3; i++) {
	FComplex* here= KVOL_LOC(v,i+xCtr-1,j+yCtr-1,k+zCtr-1);
	here->real= 1.0;
      }
  *KVOL_LOC(v,xCtr,yCtr,zCtr)= save; /* hollow out the middle */
}

static void kvol_fft3d( KVol* v, long sign )
{
  fft3d(v->data, v->dx, v->dy, v->dz, sign, "xyz");
}

static void kvol_mult( KVol* v, KVol* multBy )
{
  FComplex* r= v->data;
  FComplex* r2= multBy->data;
  FComplex* limit= v->data + v->dx*v->dy*v->dz;

  if (v->dx != multBy->dx || v->dy != multBy->dy || v->dz != multBy->dz)
    Abort("%s: attempted to multiply inconsistent volumes!\n",progname);\

  while (r<limit) {
    float rval= r->real*r2->real - r->imag*r2->imag;
    float ival= r->real*r2->imag + r->imag*r2->real;
    r->real= rval;
    r->imag= ival;
    r++;
    r2++;
  }
}

static void kvol_logical_or( KVol* v, KVol* andBy )
{
  FComplex* r= v->data;
  FComplex* r2= andBy->data;
  FComplex* limit= v->data + v->dx*v->dy*v->dz;

  if (v->dx != andBy->dx || v->dy != andBy->dy || v->dz != andBy->dz)
    Abort("%s: attempted to multiply inconsistent volumes!\n",progname);\

  while (r<limit) {
    if (r->real<0.5 && r2->real>=0.5) r->real= 1.0;
    r++;
    r2++;
  }
}

static void kvol_copy( KVol* out, const KVol* in )
{
  FComplex* r= out->data;
  FComplex* r2= in->data;
  FComplex* limit= out->data + out->dx*out->dy*out->dz;

  if (out->dx != in->dx || out->dy != in->dy || out->dz != in->dz)
    Abort("%s: attempted to copy inconsistent volumes!\n",progname);\

  while (r<limit) *r++ = *r2++;
}

static void kvol_masked_fill( KVol* v, const KVol* mask, float val )
{
  FComplex* r= v->data;
  FComplex* r2= mask->data;
  FComplex* limit= v->data + v->dx*v->dy*v->dz;

  if (v->dx != mask->dx || v->dy != mask->dy || v->dz != mask->dz)
    Abort("%s: attempted to copy inconsistent volumes!\n",progname);\

  while (r<limit) {
    if (r2->real>=0.5) { r->real= val; r->imag= 0.0; }
    r++; r2++;
  }
}

static long kvol_findNextNonzeroOffset( KVol* v, long baseOffset )
{
  FComplex* r= v->data + baseOffset;
  FComplex* end= v->data + v->dx*v->dy*v->dz;
  while ( r<end && r->real<0.5 ) r++;
  if (r<end) return ((long)(r-v->data));
  else return -1;
}

static void kvol_set( KVol* v, long offset, float val )
{
  if (offset>=v->dx*v->dy*v->dz)
    Abort("%s: attempted to set a voxel at an out-of-range offset!\n",
	  progname);
  (v->data+offset)->real= val;
}

static void kvol_thresh( KVol* v, const float thresh )
{
  FComplex* r= v->data;
  FComplex one= {1.0, 0.0};
  FComplex zero= {0.0, 0.0};

  /* Note that we are using rounding, since we're trying to
   * implement integer arithmetic.
   */
  for (r=v->data; r<v->data + v->dx*v->dy*v->dz; r++) {
    if (r->real>thresh) *r = one;
    else *r= zero;
  };
}

static void kvol_scale( KVol* v, const float scale )
{
  FComplex* r= v->data;

  /* Note that we are using rounding, since we're trying to
   * implement integer arithmetic.
   */
  for (r=v->data; r<v->data + v->dx*v->dy*v->dz; r++) {
    r->real *= scale;
    r->imag *= scale;
  };
}

static int kvol_count_changes( const KVol* v1, const KVol* v2 )
{
  FComplex* r1= v1->data;
  FComplex* r2= v2->data;
  FComplex* limit= v1->data + v1->dx*v1->dy*v1->dz;
  long result= 0;

  if (v1->dx != v2->dx || v1->dy != v2->dy || v1->dz != v2->dz)
    Abort("%s: attempted to compare inconsistent volumes!\n",progname);\

  /* Note that we are using rounding, since we're trying to
   * implement integer arithmetic.
   */
  while (r1<limit) {
    if (fabs(r1->real - r2->real) > 0.5) result++;
    r1++;
    r2++;
  }
  return result;
}

SList* compileAlg(const char* algString, int defaultMaxReps, 
		  int defaultNbrThresh)
{
  SList* result= slist_create();
  char* str= strdup(algString);
  char* t= NULL;
  char* tok;

  tok= mystrtokr(str,",",&t);
  while (tok) {
    char* p1= strchr(tok,'(');
    char* p2= NULL;
    StepType* step= (StepType*)malloc(sizeof(StepType));
    if (!step)
      Abort("%s: unable to allocate %d bytes!\n",progname,sizeof(StepType));
    if (p1) p2= strchr(p1,',');
    if (!strncasecmp(tok,"erode",strlen("erode"))) {
      step->op= OP_ERODE;
    }
    else if (!strncasecmp(tok,"dilate",strlen("dilate"))) {
      step->op= OP_DILATE;
    }
    else if (!strncasecmp(tok,"label",strlen("label"))) {
      step->op= OP_LABEL;
    }
    else if (!strncasecmp(tok,"flood",strlen("flood"))) {
      step->op= OP_FLOOD;
    }
    else Abort("%s: unrecognized algorithm step %s!\n",progname,tok);

    if ( step->op==OP_LABEL || step->op==OP_FLOOD ) {
      step->maxReps= -1;
      step->nbrThresh= 1;
      if (p1)
	fprintf(stderr,"%s: arguments to %s operation ignored\n",
		progname,getOpName(step->op));
    }
    else {
      if (p1) step->nbrThresh= atoi((p1+1));
      else step->nbrThresh= defaultNbrThresh;
      if (p2) step->maxReps= atoi(p2+1);
      else step->maxReps= defaultMaxReps;
      if (step->nbrThresh<=0)
	Abort("%s: invalid neighbor count for algorithm step %s!\n",
	      progname,tok);      
      if (step->maxReps<=0) 
	Abort("%s: invalid number of reps for algorithm step %s!\n",
	      progname,tok);
    }
    slist_append(result,step);
    tok= mystrtokr(NULL,",",&t);
  }

  free(str);
  return result;
}

static long applyDilation( KVol* kIn, KVol* kScratch, KVol* kFootprint, 
			   KVol* kMask, int thresh )
{
  kvol_copy(kScratch,kIn);
  kvol_fft3d(kScratch,1);
  kvol_mult(kScratch,kFootprint);
  kvol_fft3d(kScratch,-1);
  kvol_thresh( kScratch, 
	       (float)thresh-0.5); /* 0.5 for rounding */
  kvol_logical_or(kScratch,kIn);  /* don't lose any voxels */
  kvol_mult(kScratch,kMask);   /* bound with the fixed mask */
  return kvol_count_changes( kIn, kScratch );
}

static long applyErosion( KVol* kIn, KVol* kScratch, KVol* kFootprint, 
			   KVol* kMask, int thresh )
{
  kvol_copy(kScratch,kIn);
  kvol_fft3d(kScratch,1);
  kvol_mult(kScratch,kFootprint);
  kvol_fft3d(kScratch,-1);
  kvol_thresh( kScratch, 
	       (float)thresh-0.5); /* 0.5 for rounding */
  kvol_mult(kScratch,kIn); /* nothing is allowed to turn on */
  return kvol_count_changes( kIn, kScratch );
}

static LongStack* stack_create( KVol* proto )
{
  long initialStackSz= (long)pow((double)(proto->dx*proto->dy*proto->dz),0.667);
  long* stack= NULL;
  LongStack* result= NULL;

  if (!(result=(LongStack*)malloc(sizeof(LongStack))))
    Abort("%s: unable to allocate %d bytes!\n",sizeof(LongStack));

  if (!(stack=(long*)malloc(initialStackSz*sizeof(long))))
    Abort("%s: unable to allocate %d bytes!\n",initialStackSz*sizeof(long));

  result->top= result->base= stack;
  result->size= initialStackSz;
  if (debug_flag)
    fprintf(stderr,"Allocated a stack of total size %ld\n",result->size);
  return(result);
}

static void stack_destroy( LongStack* stk )
{
  if (stk) {
    if (stk->base) free(stk->base);
    free(stk);
  }
}

static inline void stack_push( LongStack* stk, long val )
{
  if (stk->top>=stk->base+stk->size) {
    /* Reallocate stack */
    long newSz= 2*stk->size;
    long* newBase= NULL;
    if (!(newBase= (long*)realloc(stk->base,newSz*sizeof(long))))
      Abort("%s: unable to realloc %d bytes!\n",newSz*sizeof(long));
    stk->size= newSz;
    stk->top= newBase+(stk->top-stk->base);
    stk->base= newBase;
    if (debug_flag)
      fprintf(stderr,"Stack reallocated; size is now %ld\n",stk->size);
  }
  *stk->top++= val;
}

static inline long stack_pop( LongStack* stk )
{
  if (stk->top>stk->base) {
    return(*(--stk->top));
  }
  else return 0;
}

static inline int stack_empty( LongStack* stk )
{
  return (stk->top<=stk->base);
}

#define CHECKANDSTACK( offset ) \
{ if (offset>=0 && offset<=maxOffset \
  && ISCLEAR(kIn,offset) && ISSET(kMask,offset)) stack_push(stk,offset); }

static inline void checkAndStackNeighbors( KVol* kIn, KVol* kMask,
					   LongStack* stk, long offset)
{
  long maxOffset= KVOL_OFFSET(kIn, kIn->dx-1, kIn->dy-1, kIn->dz-1);
  long nbrOffset= offset;
  nbrOffset-=kIn->dz; CHECKANDSTACK(nbrOffset); 
  nbrOffset+=2*(kIn->dz); CHECKANDSTACK(nbrOffset);
  
  nbrOffset= offset-(kIn->dy*kIn->dz); CHECKANDSTACK(nbrOffset);
  nbrOffset-=kIn->dz; CHECKANDSTACK(nbrOffset); 
  nbrOffset+=2*(kIn->dz); CHECKANDSTACK(nbrOffset);
  
  nbrOffset= offset+(kIn->dy*kIn->dz); CHECKANDSTACK(nbrOffset);
  nbrOffset-=kIn->dz; CHECKANDSTACK(nbrOffset); 
  nbrOffset+=2*(kIn->dz); CHECKANDSTACK(nbrOffset);
  
  nbrOffset= offset+1; CHECKANDSTACK(nbrOffset);
  nbrOffset-=kIn->dz; CHECKANDSTACK(nbrOffset); 
  nbrOffset+=2*(kIn->dz); CHECKANDSTACK(nbrOffset);
  
  nbrOffset= offset+1-(kIn->dy*kIn->dz); CHECKANDSTACK(nbrOffset);
  nbrOffset-=kIn->dz; CHECKANDSTACK(nbrOffset); 
  nbrOffset+=2*(kIn->dz); CHECKANDSTACK(nbrOffset);
  
  nbrOffset= offset+1+(kIn->dy*kIn->dz); CHECKANDSTACK(nbrOffset);
  nbrOffset-=kIn->dz; CHECKANDSTACK(nbrOffset); 
  nbrOffset+=2*(kIn->dz); CHECKANDSTACK(nbrOffset);
  
  nbrOffset= offset-1; CHECKANDSTACK(nbrOffset);
  nbrOffset-=kIn->dz; CHECKANDSTACK(nbrOffset); 
  nbrOffset+=2*(kIn->dz); CHECKANDSTACK(nbrOffset);
  
  nbrOffset= offset-1-(kIn->dy*kIn->dz); CHECKANDSTACK(nbrOffset);
  nbrOffset-=kIn->dz; CHECKANDSTACK(nbrOffset); 
  nbrOffset+=2*(kIn->dz); CHECKANDSTACK(nbrOffset);
  
  nbrOffset= offset-1+(kIn->dy*kIn->dz); CHECKANDSTACK(nbrOffset);
  nbrOffset-=kIn->dz; CHECKANDSTACK(nbrOffset); 
  nbrOffset+=2*(kIn->dz); CHECKANDSTACK(nbrOffset);
}

static long applyFlood( KVol* kIn, KVol* kMask )
{
  LongStack* stk= stack_create(kIn);
  long count= 0;
  long offset;

  /* For each live voxel within the mask, push all neighbors which
   * are not live but which are within the mask onto the stack.
   */
  offset= 0;
  while (1) {
    offset= kvol_findNextNonzeroOffset(kIn,offset);
    if (offset<0) break;
    if (ISSET(kMask,offset)) checkAndStackNeighbors(kIn,kMask,stk,offset);
    ++offset;
  }

  /* Process the stack.  The fact that we are protected from wrap-around
   * by a layer of buffer cells allows us to take some liberties with
   * address calculation.
   *
   * The voxels on the stack were off when they were added.  If they
   * are on, it means they've already been hit (they can be redundantly
   * stacked).  If they are off, we turn them on, and stack their valid
   * neighbors.
   */
  while (!stack_empty(stk)) {
    offset= stack_pop(stk);
    if ISCLEAR(kIn,offset) {
      SET(kIn,offset);
      count++;
      checkAndStackNeighbors(kIn,kMask,stk,offset);
      if (debug_flag) fprintf(stderr,"Count= %ld; stack depth now %ld\n",
			      count,(long)(stk->top-stk->base));
    }
  }

  stack_destroy(stk);
  return count;
}

int main( int argc, char* argv[] ) 
{
  MRI_Dataset *Input = NULL, *Output = NULL, *Mask= NULL;
  char infile[512], hdrfile[512], maskfile[512];
  char algString[512];
  long dx, dy, dz, dt, dtMask;
  KVol* kVol0= NULL;
  KVol* kVol1= NULL;
  KVol* kFootprint= NULL;
  KVol* kMask= NULL;
  float *ioVol= NULL; /* not padded */
  float* buf;
  KVol* kBefore;
  KVol* kAfter;
  long iLoop;
  long t;
  long changes= 0;
  long opCount= 0;
  int nbrThresh= DEFAULT_NBR_THRESH;
  int maxReps= DEFAULT_MAX_REPS;
  long long offset= 0;
  long long maskOffset= 0;
  long long blocksize;
  int outputIsAMask= 1;

  progname= argv[0];

  /* Print version number */
  Message( "# %s\n", rcsid );

  /* Check to see if help was requested */
  if( ( argc > 1 ) && !strcmp( argv[1], "-help" ) )
    {
      if( argc == 2 )
	Help( "selecttopic" );
      else
	Help( argv[2] );
    }

  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Get filenames */
  if (!cl_get( "output|out", "%option %s", hdrfile )) {
    fprintf(stderr,"%s: output dataset not specified!\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  if (!cl_get( "input|i", "%option %s", infile )) {
    fprintf(stderr,"%s: input dataset not specified!\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  if (!cl_get( "mask|m", "%option %s", maskfile )) {
    fprintf(stderr,"%s: mask dataset not specified!\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  if (!cl_get( "alg", "%option %s", algString )) {
    fprintf(stderr,"%s: algorithm not specified!\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  cl_get( "thresh", "%option %d", &nbrThresh );
  cl_get( "reps", "%option %d", &maxReps );
  debug_flag= cl_present("debug");
  verbose_flag= cl_present("v|verbose");

  if (cl_cleanup_check()) {
    int i;
    fprintf(stderr,"%s: invalid argument in command line:\n    ",argv[0]);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(-1);
  }

  /*** End command-line parsing ***/

  if (nbrThresh<=0) 
    Abort("%s: Invalid neighbor count threshold %d!\n",
	  progname,nbrThresh);
  if (maxReps<=0) 
    Abort("%s: Invalid repetition limit %d!\n",
	  progname,maxReps);

  /* What are we doing? */
  stepList= compileAlg(algString,maxReps,nbrThresh);

  /* Open input dataset */
  if ( !strcmp( infile, hdrfile ) || !strcmp( hdrfile, maskfile ) )
    Abort( "Input, and mask files must be distinct from output." );
  Input = mri_open_dataset( infile, MRI_READ );
  Mask  = mri_open_dataset( maskfile, MRI_READ );

  /* Check that program will function on datasets */
  checkFormat(Input,infile);
  checkFormat(Mask,maskfile);

  dx= mri_get_int(Input,"images.extent.x");
  dy= mri_get_int(Input,"images.extent.y");
  dz= mri_get_int(Input,"images.extent.z");
  if (mri_has(Input,"images.extent.t"))
    dt= mri_get_int(Input,"images.extent.t");
  else dt= 1;
  if (mri_has(Mask,"images.extent.t"))
    dtMask= mri_get_int(Input,"images.extent.t");
  else dtMask= 1;

  if (dx<1 || dy<1 || dz<1)
    Abort("%s: invalid dimensions %d x %d x %d!\n",progname,dx,dy,dz);

  if ( mri_get_int(Mask,"images.extent.x") != dx 
       || mri_get_int(Mask,"images.extent.y") != dy 
       || mri_get_int(Mask,"images.extent.z") != dz )
    Abort("%s: input and mask are not commensurate!\n",progname);

  /* Set output dataset */
  Output = mri_copy_dataset( hdrfile, Input );
  hist_add_cl( Output, argc, argv );
  mri_set_string( Output, "images.file", ".dat" );

  /* Allocate output mask storage */
  if (!(ioVol= (float*)malloc(dx*dy*dz*sizeof(float)))) 
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,dx*dy*dz*sizeof(float));

  /* Allocate and clear k-space storage */
  kVol0= kvol_create( pad_up(dx), pad_up(dy), pad_up(dz));
  kVol1= kvol_create( pad_up(dx), pad_up(dy), pad_up(dz));
  kMask= kvol_create( pad_up(dx), pad_up(dy), pad_up(dz));
  kFootprint= kvol_create( pad_up(dx), pad_up(dy), pad_up(dz));
  kBefore= kVol0;
  kAfter= kVol1;

  /* Build the neighborhood footprint */
  mark_neighborhood(kFootprint);

  /* Flip the footprint to Fourier Transform space and adjust its scale */
  kvol_fft3d(kFootprint,1);
  kvol_scale(kFootprint,
	     sqrt(kFootprint->dx*kFootprint->dy*kFootprint->dz));

  blocksize= dx*dy*dz;
  for (t=0; t<dt; t++) {
    /* Load 'em up */
    if (verbose_flag) fprintf(stderr,"######### Starting t= %d #########\n",t);
    kvol_mask_copy_in_nonzero(kBefore,
			      mri_get_chunk(Input, "images", dx*dy*dz, 
					    offset, MRI_FLOAT),
			      dx,dy,dz);
    kvol_mask_copy_in_nonzero(kMask,
			      mri_get_chunk(Mask, "images", dx*dy*dz, 
					    maskOffset, MRI_FLOAT),
			      dx,dy,dz);
    slist_totop(stepList);
    while (!slist_atend(stepList)) {
      StepType* step= (StepType*)slist_next(stepList);
      if (debug_flag) 
	fprintf(stderr,"%s on %d neighbors, max %d reps!\n",
		getOpName(step->op),step->nbrThresh,step->maxReps);
      switch (step->op) {
      case OP_DILATE:
	{
	  for (iLoop=0; iLoop<step->maxReps; iLoop++) {
	    KVol* kTmp;
	    changes= applyDilation(kBefore,kAfter,
				   kFootprint,kMask,step->nbrThresh);
	    opCount++;
	    if (verbose_flag)
	      fprintf(stderr,"Repetition %d: %d changes\n",iLoop,changes);
	    kTmp= kBefore;
	    kBefore= kAfter;
	    kAfter= kTmp;
	    if (changes==0) break;
	  }
	  outputIsAMask= 1;
	}
	break;
      case OP_ERODE:
	{
	  for (iLoop=0; iLoop<step->maxReps; iLoop++) {
	    KVol* kTmp;
	    changes= applyErosion(kBefore,kAfter,
				  kFootprint,kMask,step->nbrThresh);
	    opCount++;
	    if (verbose_flag)
	      fprintf(stderr,"Repetition %d: %d changes\n",iLoop,changes);
	    kTmp= kBefore;
	    kBefore= kAfter;
	    kAfter= kTmp;
	    if (changes==0) break;
	  }
	  outputIsAMask= 1;
	}
	break;
      case OP_FLOOD:
	{
	  changes= applyFlood(kBefore,kMask);
	  opCount++;
	  outputIsAMask= 1;
	}
	break;
      case OP_LABEL:
	{
	  long seedOffset= 0;
	  long i;
	  int done= 0;
	  float labelVal= 1.0;
	  KVol* kSeeded= kvol_create_from_prototype(kBefore);
	  KVol* kScratch= kvol_create_from_prototype(kBefore);
	  KVol* kTmp;
	  kvol_clear(kAfter);

	  while (1) {
	    /* Find next seed loc */
	    seedOffset= kvol_findNextNonzeroOffset(kBefore,seedOffset);
	    if (seedOffset<0) break; /* this means no further seeds */

	    if (verbose_flag)
	      fprintf(stderr,"Expanding seed at offset %ld\n",
		      seedOffset,step->nbrThresh);
	    kvol_clear(kSeeded);
	    kvol_set(kSeeded,seedOffset,1.0);

	    /* Flood that loc */
	    applyFlood(kSeeded,kBefore);
	    opCount++;

	    /* Copy that loc's region into a labeled volume */
	    kvol_masked_fill( kAfter, kSeeded, labelVal );

	    /* Zero that loc's region */
	    kvol_masked_fill( kBefore, kSeeded, 0.0 );

	    /* Move on to next label value */
	    labelVal += 1.0;
	    ++seedOffset;
	  }

	  /* Swap volumes, because output is done from kBefore */
	  kTmp= kBefore;
	  kBefore= kAfter;
	  kAfter= kTmp;

	  /* Allow scalar rather than boolean output */
	  outputIsAMask= 0;

	  kvol_destroy(kSeeded);
	  kvol_destroy(kScratch);
	}
	break;
      }
      free(step);
    }

    /* Write output.  Note that last result is back
     * in kBefore already.
     */
    if (outputIsAMask)
      kvol_mask_copy_out_nonzero( ioVol, kBefore, dx, dy, dz );
    else
      kvol_copy_out_nonzero( ioVol, kBefore, dx, dy, dz );
    mri_set_chunk( Output, "images", dx*dy*dz, offset, MRI_FLOAT, ioVol );

    /* Make a mark on the screen and set up for next step */
    if (!verbose_flag) {
      if (((t+1) % 60) != 0 && t != dt-1 ) {
	Message("#");
      }
      else {
	Message("# %ld\n",(long)(t+1));
      }
    }

    offset += blocksize;
    maskOffset += blocksize;
    if (t>=dtMask) maskOffset= 0;
    kvol_clear(kBefore);
    kvol_clear(kMask);
  }

  /* Clean up */
  mri_close_dataset( Mask );
  mri_close_dataset( Input );
  mri_close_dataset( Output );
  kvol_destroy(kVol0);
  kvol_destroy(kVol1);
  kvol_destroy(kFootprint);
  free(ioVol);

  Message( "#      Morphology operations complete (%d ops).\n",
	   opCount);
  return 0;
}

