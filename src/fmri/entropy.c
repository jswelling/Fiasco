/************************************************************
 *                                                          *
 *  entropy.c                                            *
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
 *  Original programming by Chris Hefferan 6/04             *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "misc.h"
#include "stdcrg.h"

static char rcsid[] = "$Id: entropy.c,v 1.9 2005/03/08 22:24:29 welling Exp $";

EntropyContext* ent_createContext()
{
  EntropyContext* result;

  if (!(result=(EntropyContext*)malloc(sizeof(EntropyContext))))
    Abort("entropy: unable to allocate %d bytes!\n",
	  sizeof(EntropyContext));

  result->histogram= NULL;
  result->histogramArraySize= 0;
  result->debugFlag= 0;
  result->verboseFlag= 0;
  result->num_of_boxes= 0;
  result->min= 0.0;
  result->max= 0.0;
  result->min_set= 0.0;
  result->max_set= 0.0;

  return result;
}

void ent_destroyContext(EntropyContext* ctx)
{
  if (ctx->histogram) free(ctx->histogram);
  free(ctx);
}

void ent_setDebug(EntropyContext* ctx, int val)
{
  ctx->debugFlag= val;
}

void ent_setVerbose(EntropyContext* ctx, int val)
{
  ctx->verboseFlag= val;
}

int ent_getDebug(const EntropyContext* ctx)
{
  return ctx->debugFlag;
}

int ent_getVerbose(const EntropyContext* ctx)
{
  return ctx->verboseFlag;
}

void ent_setMin(EntropyContext* ctx, double val)
{
  ctx->min= val;
  ctx->min_set= 1;
}

void ent_setMax(EntropyContext* ctx, double val)
{
  ctx->max= val;
  ctx->max_set= 1;
}

double ent_getMin(const EntropyContext* ctx)
{
  return ctx->min;
}

double ent_getMax(const EntropyContext* ctx)
{
  return ctx->max;
}

int ent_getMinSet(const EntropyContext* ctx)
{
  return ctx->min_set;
}

int ent_getMaxSet(const EntropyContext* ctx)
{
  return ctx->max_set;
}

void ent_unsetMin(EntropyContext* ctx)
{
  ctx->min_set= 0;
}

void ent_unsetMax(EntropyContext* ctx)
{
  ctx->max_set= 0;
}

void ent_setNBins(EntropyContext* ctx, const long nbins)
{
  ctx->num_of_boxes= nbins;
}

long ent_getNBins(const EntropyContext* ctx)
{
  return ctx->num_of_boxes;
}

int ent_getNBinsSet(const EntropyContext* ctx)
{
  return (ctx->num_of_boxes != 0);
}

MutualInfoContext* ent_createMIContext()
{
  MutualInfoContext* result;

  if (!(result=(MutualInfoContext*)malloc(sizeof(MutualInfoContext))))
    Abort("entropy: unable to allocate %d bytes!\n",
	  sizeof(MutualInfoContext));

  result->histogram= NULL;
  result->histogramArraySize= 0;
  result->debugFlag= 0;
  result->verboseFlag= 0;
  result->num_of_boxes= 0;
  result->ec1= ent_createContext();
  result->ec2= ent_createContext();

  return result;
}

void ent_destroyMIContext(MutualInfoContext* ctx)
{
  if (ctx->histogram) free(ctx->histogram);
  if (ctx->ec1) ent_destroyContext(ctx->ec1);
  if (ctx->ec2) ent_destroyContext(ctx->ec2);
}

void ent_setMIDebug(MutualInfoContext* ctx, int val)
{
  ctx->debugFlag= val;
  ent_setDebug(ctx->ec1,val);
  ent_setDebug(ctx->ec2,val);
}

void ent_setMIVerbose(MutualInfoContext* ctx, int val)
{
  ctx->verboseFlag= val;
  ent_setVerbose(ctx->ec1,val);
  ent_setVerbose(ctx->ec2,val);
}

int ent_getMIDebug(const MutualInfoContext* ctx)
{
  return ctx->debugFlag;
}

int ent_getMIVerbose(const MutualInfoContext* ctx)
{
  return ctx->verboseFlag;
}

void ent_setMIMin1(MutualInfoContext* ctx, double val)
{
  ent_setMin(ctx->ec1,val);
}

void ent_setMIMax1(MutualInfoContext* ctx, double val)
{
  ent_setMax(ctx->ec1,val);
}

double ent_getMIMin1(const MutualInfoContext* ctx)
{
  return ent_getMin(ctx->ec1);
}

double ent_getMIMax1(const MutualInfoContext* ctx)
{
  return ent_getMax(ctx->ec1);
}

int ent_getMIMin1Set(const MutualInfoContext* ctx)
{
  return ent_getMinSet(ctx->ec1);
}

int ent_getMIMax1Set(const MutualInfoContext* ctx)
{
  return ent_getMaxSet(ctx->ec1);
}

void ent_unsetMIMin1(MutualInfoContext* ctx)
{
  ent_unsetMin(ctx->ec1);
}

void ent_unsetMIMax1(MutualInfoContext* ctx)
{
  ent_unsetMax(ctx->ec1);
}

void ent_setMIMin2(MutualInfoContext* ctx, double val)
{
  ent_setMin(ctx->ec2,val);
}

void ent_setMIMax2(MutualInfoContext* ctx, double val)
{
  ent_setMax(ctx->ec2,val);
}

double ent_getMIMin2(const MutualInfoContext* ctx)
{
  return ent_getMin(ctx->ec2);
}

double ent_getMIMax2(const MutualInfoContext* ctx)
{
  return ent_getMax(ctx->ec2);
}

int ent_getMIMin2Set(const MutualInfoContext* ctx)
{
  return ent_getMinSet(ctx->ec2);
}

int ent_getMIMax2Set(const MutualInfoContext* ctx)
{
  return ent_getMaxSet(ctx->ec2);
}

void ent_unsetMIMin2(MutualInfoContext* ctx)
{
  ent_unsetMin(ctx->ec2);
}

void ent_unsetMIMax2(MutualInfoContext* ctx)
{
  ent_unsetMax(ctx->ec2);
}

/* NOTE: unsetNBins(ctx) is accomplished by doing setNBins(0) */
void ent_setMINBins(MutualInfoContext* ctx, const long nbins)
{
  ctx->num_of_boxes= nbins;
  ent_setNBins(ctx->ec1,nbins);
  ent_setNBins(ctx->ec2,nbins);
}

long ent_getMINBins(const MutualInfoContext* ctx)
{
  return ctx->num_of_boxes;
}

int ent_getMINBinsSet(const MutualInfoContext* ctx)
{
  return (ctx->num_of_boxes != 0);
}

static int compareDouble(const void* v1, const void* v2)
{
  double d1= *(double*)v1;
  double d2= *(double*)v2;
  if (d1<d2) return -1;
  else if (d1>d2) return 1;
  else return 0;
}

static long pickNBins( long nSamples, long stride, double* samples )
{
  long num_of_boxes= 0;
  /*
  double iqr;
  double mean;
  double stdv;
  */

  /* Note- many of these rules require valid 'samples' input, and
   * most of the calling routines don't provide it!
   */

  /* Joel's stupid rule */
  /* num_of_boxes = sqrt(nSamples); */
  /* Sturges' rule */
  num_of_boxes = ceil(log((double)nSamples)/log(2.0) + 1);
  
  /* Freedman and Diaconis */
  /* 
     qsort(samples, nSamples, sizeof(double), compareDouble);
     iqr= samples[(3*nSamples)/4]-samples[(nSamples/4)];
     num_of_boxes =(long) (2.0*iqr*pow(nSamples,-1.0/3.0)); 
  */ 
  
  /* Scott (1979) */
  /*
    mean= 0.0;
    for (i=0; i<nSamples; i++) mean+=samples[i];
    mean /= (double)nSamples;
    stdv= 0.0;
    for (i=0; i<nSamples; i++) stdv += (samples[i]-mean)*(samples[i]-mean);
    stdv= sqrt(stdv/(double)(nSamples-1));
    num_of_boxes = (long) (3.5*stdev*pow(nSamples,-1.0/3.0)); 
  */
  return num_of_boxes;
}

static long scanMaskedRegionDouble( EntropyContext* ctx, const double* img,
				    const int* mask, long pixels, 
				    long stride, long maskstride )
{
  long ofInterest= 0;
  long i;
  long j;
  long first= -1;

  /** Count the number of voxels that will be used **/
  j= 0;
  for( i=0; i<pixels; i++) { 
    if(mask[j] != 0) {
      if (!ofInterest) first= i;
      ofInterest++;
    }
    j += maskstride;
  }
  if (ofInterest>0) {

    if (!ctx->max_set) {
      i= first*stride;
      ctx->max = img[i];
      i += stride;
      j = (first+1)*maskstride;
      while (i<pixels*stride) {
	if( img[i] > ctx->max && mask[j] != 0) ctx->max = img[i]; 
	i += stride;
	j += maskstride;
      }
      ctx->max_set= 1;
    }
    
    if (!ctx->min_set) {
      i= first*stride;
      ctx->min = img[i];
      i += stride;
      j = (first+1)*maskstride;
      while (i<pixels*stride) {
	if( img[i] < ctx->min && mask[j] != 0) ctx->min = img[i]; 
	i += stride;
	j += maskstride;
      }
      ctx->min_set= 1;
    }
  }
  return ofInterest;
}

static void scanRegionDouble( EntropyContext* ctx, const double* img,
			      long pixels, long stride )
{
  long i;

  if (!ctx->max_set) {
    ctx->max= img[0];
    for (i=stride; i<pixels*stride; i+=stride)
      if( img[i] > ctx->max) ctx->max = img[i]; 
    ctx->max_set= 1;
  }
  
  if (!ctx->min_set) {
    ctx->min= img[0];
    for (i=stride; i<pixels*stride; i+=stride)
      if( img[i] < ctx->min) ctx->min = img[i]; 
    ctx->min_set= 1;
  }
  
}

static long scanMaskedRegionFloat( EntropyContext* ctx, const float* img,
				    const int* mask, long pixels, 
				    long stride, long maskstride )
{
  long ofInterest= 0;
  long i;
  long j;
  long first= -1;

  /** Count the number of voxels that will be used **/
  j= 0;
  for( i=0; i<pixels; i++) { 
    if(mask[j] != 0) {
      if (!ofInterest) first= i;
      ofInterest++;
    }
    j += maskstride;
  }
  if (ofInterest>0) {

    if (!ctx->max_set) {
      i= first*stride;
      ctx->max = img[i];
      i += stride;
      j = (first+1)*maskstride;
      while (i<pixels*stride) {
	if( img[i] > ctx->max && mask[j] != 0) ctx->max = img[i]; 
	i += stride;
	j += maskstride;
      }
      ctx->max_set= 1;
    }
    
    if (!ctx->min_set) {
      i= first*stride;
      ctx->min = img[i];
      i += stride;
      j = (first+1)*maskstride;
      while (i<pixels*stride) {
	if( img[i] < ctx->min && mask[j] != 0) ctx->min = img[i]; 
	i += stride;
	j += maskstride;
      }
      ctx->min_set= 1;
    }
  }
  return ofInterest;
}

static void scanRegionFloat( EntropyContext* ctx, const float* img,
			     long pixels, long stride )
{
  long i;

  if (!ctx->max_set) {
    ctx->max= img[0];
    for (i=stride; i<pixels*stride; i+=stride)
      if( img[i] > ctx->max) ctx->max = img[i]; 
    ctx->max_set= 1;
  }
  
  if (!ctx->min_set) {
    ctx->min= img[0];
    for (i=stride; i<pixels*stride; i+=stride)
      if( img[i] < ctx->min) ctx->min = img[i]; 
    ctx->min_set= 1;
  }
  
}

double ent_calcMaskedImageEntropyDouble( EntropyContext* ctx, 
					 const double* img, const int* mask, 
					 long dx, long dy, long dz,
					 long stride, long maskstride )
{
  long a,b,c,d,e,f,g,i,h;
  
  long pixels;
  long outOfRangeCounts= 0;
  double range;
  double boxRange;
  double shannonEntropy;
  long ofInterest;

  if (stride<1) Abort("calcMaskedImageEntropyDouble: nonsense stride!\n");
  if (maskstride<1) Abort("calcMaskedImageEntropyDouble: nonsense maskstride!\n");
  pixels = dx*dy*dz;
  if ((ofInterest=scanMaskedRegionDouble(ctx, img, mask, pixels, stride, maskstride)) == 0)
    return 0.0;
  
  range = ctx->max - ctx->min;
  if (range==0.0) return 0.0;
  
  /* set up bins */ 
  if (ctx->num_of_boxes==0) /* not set from command line */
    ctx->num_of_boxes= pickNBins(ofInterest,stride,NULL);

  boxRange = range / (ctx->num_of_boxes-1);
  if (ctx->debugFlag)
    fprintf(stderr,"binning into %ld bins, max= %g, min= %g, box range %g\n",
	    ctx->num_of_boxes,ctx->max,ctx->min,boxRange);
  
  if (ctx->num_of_boxes > ctx->histogramArraySize) {
    if (ctx->debugFlag) fprintf(stderr,"Reallocating histogram array from %ld to %ld\n",
		       ctx->histogramArraySize,ctx->num_of_boxes);
    if (ctx->histogramArraySize) free(ctx->histogram);
    if (!(ctx->histogram= (long*)malloc(ctx->num_of_boxes*sizeof(long))))
      Abort("calcMaskedImageEntropyDouble: unable to allocate %d bytes!\n",
	    ctx->num_of_boxes*sizeof(long));
    ctx->histogramArraySize= ctx->num_of_boxes;
  }
  
  /* initialize histogram array */  
  for(b = 0; b < ctx->num_of_boxes; b++)
    ctx->histogram[b] = 0;
  
  /* fill histogram array */
  c= 0;
  d= 0;
  while (c<stride*pixels) {
    if(mask[d] != 0){
      long bin = (long)((img[c]-ctx->min)*(1.0/boxRange));
      if(bin>=0 && bin<ctx->num_of_boxes)
	ctx->histogram[bin]++;
      else {
	outOfRangeCounts++;
	if (ctx->debugFlag) fprintf(stderr,"miss: bin would be %ld\n",bin);
      }
    }
    c += stride;
    d += maskstride;
  }
  if (outOfRangeCounts && ctx->verboseFlag)
    Warning(1,"calcMaskedImageEntropyDouble: %ld of %ld samples out of binned range!\n",
	    outOfRangeCounts,ofInterest);
  
  /* sum histogram entries */
  /* Calculate the entropy using eqn (2) from the Mutual Information 
   * paper
   */ 
  shannonEntropy = 0.0;
  for(g = 0; g < ctx->num_of_boxes; g++)
    {
      if(ctx->histogram[g] != 0) /* only accounting for significant probabilities */ 
	{
	  double prob= (1.0/ofInterest) * (double)ctx->histogram[g];   
	  shannonEntropy += -prob*log(prob); 
	}
    }
  shannonEntropy /= log(2.0); /* rescale log_10 to log_2 */
  
  if (ctx->debugFlag) fprintf(stderr,"Entropy: %lg\n",shannonEntropy);
  return shannonEntropy;
}

double ent_calcImageEntropyDouble( EntropyContext* ctx, 
				   const double* img, 
				   long dx, long dy, long dz, long stride )
{
  long a,b,c,e,g,i;
  
  long pixels;
  long outOfRangeCounts= 0;
  double range;
  double boxRange;
  double shannonEntropy;

  if (stride<1) Abort("calcImageEntropyDouble: nonsense stride!\n");

  pixels = dx*dy*dz;
  scanRegionDouble(ctx, img, pixels, stride);
  range = ctx->max - ctx->min;
  if (range==0.0) return 0.0;
  
  /* set up bins */ 
  if (ctx->num_of_boxes==0) /* not set from command line */
    ctx->num_of_boxes= pickNBins(pixels,stride,NULL);
  boxRange = range / (ctx->num_of_boxes-1);
  if (ctx->debugFlag)
    fprintf(stderr,"binning into %ld bins, max= %g, min= %g, box range %g\n",
	    ctx->num_of_boxes,ctx->max,ctx->min,boxRange);
  
  if (ctx->num_of_boxes > ctx->histogramArraySize) {
    if (ctx->debugFlag) fprintf(stderr,
			    "Reallocating histogram array from %ld to %ld\n",
			    ctx->histogramArraySize,
				ctx->num_of_boxes*ctx->num_of_boxes);
    if (ctx->histogramArraySize) free(ctx->histogram);
    if (!(ctx->histogram= (long*)malloc(ctx->num_of_boxes*sizeof(long))))
      Abort("calcImageEntropyDouble: unable to allocate %d bytes!\n",
	    ctx->num_of_boxes*sizeof(long));
    ctx->histogramArraySize= ctx->num_of_boxes;
  }
  
  /* initialize histogram array */  
  for(b = 0; b < ctx->num_of_boxes; b++)
    ctx->histogram[b] = 0;
  
  /* fill histogram array */
  for(c=0; c < stride*pixels; c+=stride)
    {
      long bin = (long)((img[c]-ctx->min)*(1.0/boxRange));
      if (bin>=0 && bin<ctx->num_of_boxes)
	ctx->histogram[bin]++;
      else {
	outOfRangeCounts++;
	if (ctx->debugFlag) fprintf(stderr,"miss: bin would be %ld\n",bin);
      }
    }
  if (outOfRangeCounts && ctx->verboseFlag)
    Warning(1,"calcImageEntropyDouble: %ld of %ld samples out of binned range!\n",
	    outOfRangeCounts,pixels);

  /* sum histogram entries */
  /* Calculate the entropy using eqn (2) from the Mutual Information 
   * paper
   */ 
  shannonEntropy = 0.0;
  for(g = 0; g < ctx->num_of_boxes; g++)
    {
      if(ctx->histogram[g] != 0) /* only accounting for significant probabilities */ 
	{
	  double prob= (1.0/pixels) * (double)ctx->histogram[g];   
	  shannonEntropy += -prob*log(prob); 
	}
    }
  shannonEntropy /= log(2.0); /* rescale log_10 to log_2 */
  
  if (ctx->debugFlag) fprintf(stderr,"Entropy: %lg\n",shannonEntropy);
  return shannonEntropy;
}

double ent_calcMaskedJointEntropyDouble( MutualInfoContext* mc,
					 const double* img1,
					 const double* img2,
					 const int* mask,
					 long dx, long dy, long dz, 
					 long stride1, long stride2,
					 long maskstride )
{
  long a,b,c,d,e,f,g;
  
  long pixels= dx*dy*dz;
  double min1;
  double min2;
  double range1;
  double range2;
  double boxRange1;
  double boxRange2;
  double jointShannonEntropy= 0.0;
  long ofInterest= 0;
  long outOfRangeCounts= 0;

  if (stride1<1) Abort("calcMaskedJointEntropyDouble: nonsense stride1!\n");
  if (stride2<1) Abort("calcMaskedJointEntropyDouble: nonsense stride2!\n");
  if (maskstride<1) Abort("calcMaskedJointEntropyDouble: nonsense maskstride!\n");

  /** Count the number of voxels that will be used **/
  ofInterest= scanMaskedRegionDouble(mc->ec1,img1,mask,pixels,stride1,maskstride);
  if (ofInterest==0) return 0.0;
  (void)scanMaskedRegionDouble(mc->ec2,img2,mask,pixels,stride2,maskstride);

  min1= ent_getMin(mc->ec1);
  min2= ent_getMin(mc->ec2);
  range1 = ent_getMax(mc->ec1) - ent_getMin(mc->ec1);
  range2 = ent_getMax(mc->ec2) - ent_getMin(mc->ec2);
  if (mc->debugFlag) 
    fprintf(stderr,"min1= %lg, min2= %lg, range1= %lg, range2= %lg\n",
	    min1,min2,range1,range2);
      
  if (range1 == 0.0 && range2 == 0.0) return 0.0;

  /* set up bins */ 
  if (mc->num_of_boxes==0) 
    mc->num_of_boxes= pickNBins(ofInterest,stride1,NULL);

  boxRange1 = range1 / (mc->num_of_boxes - 1);
  boxRange2 = range2 / (mc->num_of_boxes - 1);
  
  if (mc->num_of_boxes*mc->num_of_boxes > mc->histogramArraySize) {
    if (mc->debugFlag) 
      fprintf(stderr,"Reallocating histogram array from %ld to %ld\n",
	      mc->histogramArraySize,mc->num_of_boxes*mc->num_of_boxes);
    if (mc->histogramArraySize) free(mc->histogram);
    if (!(mc->histogram= 
	  (long*)malloc(mc->num_of_boxes*mc->num_of_boxes*sizeof(long))))
      Abort("ent_calcMaskedJointEntropyDouble: unable to allocate %d bytes!\n",
	    mc->num_of_boxes*mc->num_of_boxes*sizeof(long));
    mc->histogramArraySize= mc->num_of_boxes* mc->num_of_boxes;
  }
 
  /* initialize histogram array */  
  for(b = 0; b < mc->num_of_boxes * mc->num_of_boxes; b++)    
     mc->histogram[b] = 0;
     
  /* fill histogram array */
  e= 0;
  f= 0;
  g= 0;
  while ( e<stride1*pixels ) {
    if(mask[g] != 0) { 
      long bin = (long)((img1[e]-min1)*(1.0/boxRange1));
      long bin2 = (long)((img2[f]-min2)*(1.0/boxRange2));
      if (bin>=0 && bin<mc->num_of_boxes
	  && bin2>=0 && bin2<mc->num_of_boxes) {
	long location = (long)(bin*mc->num_of_boxes) + bin2; 
	mc->histogram[location]++;
      }
      else {
	outOfRangeCounts++;
	if (mc->debugFlag) fprintf(stderr,"miss: bins would be %ld %ld\n",bin,bin2);
      }
    }
    e += stride1;
    f += stride2;
    g += maskstride;
  }
  if (outOfRangeCounts && mc->verboseFlag)
    Warning(1,"calcMaskedJointEntropyDouble: %ld of %ld samples out of binned range!\n",
	    outOfRangeCounts,ofInterest);

  if (mc->debugFlag) {
    FILE* f;
    fprintf(stderr,"Writing %ld-by-%ld histogram of raw longs to debug.raw\n",
	    mc->num_of_boxes, mc->num_of_boxes);
    f= fopen("debug.raw","w");
    (void)fwrite(mc->histogram,sizeof(long),mc->num_of_boxes*mc->num_of_boxes,
		 f);
    (void)fclose(f);
  }
  
  /* sum histogram entries */
  jointShannonEntropy = 0.0;
  
  /* Calculate the entropy using eqn (2) from the Mutual Information 
   * paper
   */ 
  
  for(d = 0; d < mc->num_of_boxes*mc->num_of_boxes; d++)
    {
      if(mc->histogram[d] != 0) 
	{
	  double prob= (1.0/ofInterest) * (double)mc->histogram[d];   
	  jointShannonEntropy += -prob*log(prob); 
	}
      
    }
  jointShannonEntropy /= log(2.0); /* rescale log_10 to log_2 */
  if (mc->debugFlag) fprintf(stderr,"Joint entropy: %lg\n",jointShannonEntropy);
  return jointShannonEntropy;
}

double ent_calcJointEntropyDouble( MutualInfoContext* mc,
				   const double* img1,
				   const double* img2,
				   long dx, long dy, long dz, 
				   long stride1, long stride2 )
{
  long a,b,c,d,e,f;
  
  long pixels= dx*dy*dz;
  double min1;
  double min2;
  double range1;
  double range2;
  double boxRange1;
  double boxRange2;
  double jointShannonEntropy= 0.0;
  long outOfRangeCounts= 0;

  if (stride1<1) Abort("calcJointEntropyDouble: nonsense stride1!\n");
  if (stride2<1) Abort("calcJointEntropyDouble: nonsense stride2!\n");

  scanRegionDouble(mc->ec1, img1, pixels, stride1);
  scanRegionDouble(mc->ec2, img1, pixels, stride2);

  min1= ent_getMin(mc->ec1);
  min2= ent_getMin(mc->ec2);
  range1 = ent_getMax(mc->ec1) - ent_getMin(mc->ec1);
  range2 = ent_getMax(mc->ec2) - ent_getMin(mc->ec2);
  if (mc->debugFlag) 
    fprintf(stderr,"min1= %lg, min2= %lg, range1= %lg, range2= %lg\n",
	    min1,min2,range1,range2);
      
  if (range1 == 0.0 && range2 == 0.0) return 0.0;

  /* set up bins */ 
  if (mc->num_of_boxes==0) 
    mc->num_of_boxes= pickNBins(pixels,stride1,NULL);

  boxRange1 = range1 / (mc->num_of_boxes - 1);
  boxRange2 = range2 / (mc->num_of_boxes - 1);
  
  if (mc->num_of_boxes*mc->num_of_boxes > mc->histogramArraySize) {
    if (mc->debugFlag) 
      fprintf(stderr,"Reallocating histogram array from %ld to %ld\n",
	      mc->histogramArraySize,mc->num_of_boxes*mc->num_of_boxes);
    if (mc->histogramArraySize) free(mc->histogram);
    if (!(mc->histogram= 
	  (long*)malloc(mc->num_of_boxes*mc->num_of_boxes*sizeof(long))))
      Abort("ent_calcJointEntropyDouble: unable to allocate %d bytes!\n",
	    mc->num_of_boxes*mc->num_of_boxes*sizeof(long));
    mc->histogramArraySize= mc->num_of_boxes* mc->num_of_boxes;
  }
 
  /* initialize histogram array */  
  for(b = 0; b < mc->num_of_boxes * mc->num_of_boxes; b++)    
     mc->histogram[b] = 0;
     
  /* fill histogram array */
  f= 0;
  c= 0;
  while (c<stride1*pixels) {
    long bin = (long)((img1[c]-min1)*(1.0/boxRange1));
    long bin2 = (long)((img2[f]-min2)*(1.0/boxRange2));
    if (bin>=0 && bin<mc->num_of_boxes
	&& bin2>=0 && bin2<mc->num_of_boxes) {
      long location = (long)(bin*mc->num_of_boxes) + bin2; 
      mc->histogram[location]++;
    }
    else {
      outOfRangeCounts++;
      if (mc->debugFlag) fprintf(stderr,"miss: bins would be %ld %ld\n",bin,bin2);
    }
    c += stride1;
    f += stride2;
  }
  if (outOfRangeCounts && mc->verboseFlag)
    Warning(1,"calcJointEntropyDouble: %ld of %ld samples out of binned range!\n",
	    outOfRangeCounts,pixels);

  if (mc->debugFlag) {
    FILE* f;
    fprintf(stderr,"Writing %ld-by-%ld histogram of raw longs to debug.raw\n",
	    mc->num_of_boxes, mc->num_of_boxes);
    f= fopen("debug.raw","w");
    (void)fwrite(mc->histogram,sizeof(long),mc->num_of_boxes*mc->num_of_boxes,
		 f);
    (void)fclose(f);
  }
  
  /* sum histogram entries */
  jointShannonEntropy = 0.0;
  
  /* Calculate the entropy using eqn (2) from the Mutual Information 
   * paper
   */ 
  for(d = 0; d < mc->num_of_boxes*mc->num_of_boxes; d++)
    {
      if(mc->histogram[d] != 0) 
	{
	  double prob= (1.0/pixels) * (double)mc->histogram[d];   
	  jointShannonEntropy += -prob*log(prob); 
	}
    }
  jointShannonEntropy /= log(2.0); /* rescale log_10 to log_2 */
  if (mc->debugFlag) fprintf(stderr,"Joint entropy: %lg\n",jointShannonEntropy);
  return jointShannonEntropy;
}

double ent_calcMaskedMutualInformationDouble( MutualInfoContext* mc,
					      const double* img1, 
					      const double* img2, 
					      const int* mask, 
					      long dx, long dy, long dz,
					      long stride1, long stride2,
					      long maskstride )
{
  double mask1Entropy = ent_calcMaskedImageEntropyDouble(mc->ec1, 
							 img1, mask, dx, dy, 
							 dz, stride1, maskstride);
  double mask2Entropy = ent_calcMaskedImageEntropyDouble(mc->ec2, 
							 img2, mask, dx, dy, 
							 dz, stride2, maskstride);
  double mutualInformation = 
    mask1Entropy + mask2Entropy 
    - ent_calcMaskedJointEntropyDouble(mc, img1, img2, mask, dx, dy, dz, 
				       stride1, stride2, maskstride);

  return mutualInformation;
}

double ent_calcMutualInformationDouble( MutualInfoContext* mc,
					const double* img1, 
					const double* img2, 
					long dx, long dy, long dz,
					long stride1, long stride2 )
{
  double img1Entropy = ent_calcImageEntropyDouble(mc->ec1, 
						  img1, dx, dy, dz, stride1);
  double img2Entropy = ent_calcImageEntropyDouble(mc->ec2, 
						  img2, dx, dy, dz, stride2);
  double mutualInformation = 
    img1Entropy + img2Entropy - ent_calcJointEntropyDouble(mc, img1, img2,
							   dx, dy, dz, 
							   stride1, stride2);

  return mutualInformation;
}

double ent_calcMaskedImageEntropyFloat( EntropyContext* ctx, 
					 const float* img, const int* mask, 
					 long dx, long dy, long dz,
					 long stride, long maskstride )
{
  long a,b,c,d,e,f,g,i,h;
  
  long pixels;
  long outOfRangeCounts= 0;
  double range;
  double boxRange;
  double shannonEntropy;
  long ofInterest;

  if (stride<1) Abort("calcMaskedImageEntropyFloat: nonsense stride!\n");
  if (maskstride<1) Abort("calcMaskedImageEntropyFloat: nonsense maskstride!\n");
  pixels = dx*dy*dz;
  if ((ofInterest=scanMaskedRegionFloat(ctx, img, mask, pixels, stride, maskstride)) == 0)
    return 0.0;
  
  range = ctx->max - ctx->min;
  if (range==0.0) return 0.0;
  
  /* set up bins */ 
  if (ctx->num_of_boxes==0) /* not set from command line */
    ctx->num_of_boxes= pickNBins(ofInterest,stride,NULL);

  boxRange = range / (ctx->num_of_boxes-1);
  if (ctx->debugFlag)
    fprintf(stderr,"binning into %ld bins, max= %g, min= %g, box range %g\n",
	    ctx->num_of_boxes,ctx->max,ctx->min,boxRange);
  
  if (ctx->num_of_boxes > ctx->histogramArraySize) {
    if (ctx->debugFlag) fprintf(stderr,"Reallocating histogram array from %ld to %ld\n",
		       ctx->histogramArraySize,ctx->num_of_boxes);
    if (ctx->histogramArraySize) free(ctx->histogram);
    if (!(ctx->histogram= (long*)malloc(ctx->num_of_boxes*sizeof(long))))
      Abort("calcMaskedImageEntropyFloat: unable to allocate %d bytes!\n",
	    ctx->num_of_boxes*sizeof(long));
    ctx->histogramArraySize= ctx->num_of_boxes;
  }
  
  /* initialize histogram array */  
  for(b = 0; b < ctx->num_of_boxes; b++)
    ctx->histogram[b] = 0;
  
  /* fill histogram array */
  c= 0;
  d= 0;
  while (c<stride*pixels) {
    if(mask[d] != 0){
      long bin = (long)((img[c]-ctx->min)*(1.0/boxRange));
      if(bin>=0 && bin<ctx->num_of_boxes)
	ctx->histogram[bin]++;
      else {
	outOfRangeCounts++;
	if (ctx->debugFlag) fprintf(stderr,"miss: bin would be %ld\n",bin);
      }	
    }
    c += stride;
    d += maskstride;
  }
  if (outOfRangeCounts && ctx->verboseFlag)
    Warning(1,"calcMaskedImageEntropyFloat: %ld of %ld samples out of binned range!\n",
	    outOfRangeCounts,ofInterest);
  
  /* sum histogram entries */
  /* Calculate the entropy using eqn (2) from the Mutual Information 
   * paper
   */ 
  shannonEntropy = 0.0;
  for(g = 0; g < ctx->num_of_boxes; g++)
    {
      if(ctx->histogram[g] != 0) /* only accounting for significant probabilities */ 
	{
	  double prob= (1.0/ofInterest) * (double)ctx->histogram[g];   
	  shannonEntropy += -prob*log(prob); 
	}
    }
  shannonEntropy /= log(2.0); /* rescale log_10 to log_2 */
  
  if (ctx->debugFlag) fprintf(stderr,"Entropy: %lg\n",shannonEntropy);
  return shannonEntropy;
}

double ent_calcImageEntropyFloat( EntropyContext* ctx, 
				   const float* img, 
				   long dx, long dy, long dz, long stride )
{
  long a,b,c,e,g,i;
  
  long pixels;
  long outOfRangeCounts= 0;
  double range;
  double boxRange;
  double shannonEntropy;

  if (stride<1) Abort("calcImageEntropyFloat: nonsense stride!\n");

  pixels = dx*dy*dz;
  scanRegionFloat(ctx, img, pixels, stride);
  range = ctx->max - ctx->min;
  if (range==0.0) return 0.0;
  
  /* set up bins */ 
  if (ctx->num_of_boxes==0) /* not set from command line */
    ctx->num_of_boxes= pickNBins(pixels,stride,NULL);
  boxRange = range / (ctx->num_of_boxes-1);
  if (ctx->debugFlag)
    fprintf(stderr,"binning into %ld bins, max= %g, min= %g, box range %g\n",
	    ctx->num_of_boxes,ctx->max,ctx->min,boxRange);
  
  if (ctx->num_of_boxes > ctx->histogramArraySize) {
    if (ctx->debugFlag) fprintf(stderr,
			    "Reallocating histogram array from %ld to %ld\n",
			    ctx->histogramArraySize,
				ctx->num_of_boxes*ctx->num_of_boxes);
    if (ctx->histogramArraySize) free(ctx->histogram);
    if (!(ctx->histogram= (long*)malloc(ctx->num_of_boxes*sizeof(long))))
      Abort("calcImageEntropyFloat: unable to allocate %d bytes!\n",
	    ctx->num_of_boxes*sizeof(long));
    ctx->histogramArraySize= ctx->num_of_boxes;
  }
  
  /* initialize histogram array */  
  for(b = 0; b < ctx->num_of_boxes; b++)
    ctx->histogram[b] = 0;
  
  /* fill histogram array */
  for(c=0; c < stride*pixels; c+=stride)
    {
      long bin = (long)((img[c]-ctx->min)*(1.0/boxRange));
      if (bin>=0 && bin<ctx->num_of_boxes)
	ctx->histogram[bin]++;
      else {
	outOfRangeCounts++;
	if (ctx->debugFlag) fprintf(stderr,"miss: bin would be %ld\n",bin);
      }
    }
  if (outOfRangeCounts && ctx->verboseFlag)
    Warning(1,"calcImageEntropyFloat: %ld of %ld samples out of binned range!\n",
	    outOfRangeCounts,pixels);

  /* sum histogram entries */
  /* Calculate the entropy using eqn (2) from the Mutual Information 
   * paper
   */ 
  shannonEntropy = 0.0;
  for(g = 0; g < ctx->num_of_boxes; g++)
    {
      if(ctx->histogram[g] != 0) /* only accounting for significant probabilities */ 
	{
	  double prob= (1.0/pixels) * (double)ctx->histogram[g];   
	  shannonEntropy += -prob*log(prob); 
	}
    }
  shannonEntropy /= log(2.0); /* rescale log_10 to log_2 */
  
  if (ctx->debugFlag) fprintf(stderr,"Entropy: %lg\n",shannonEntropy);
  return shannonEntropy;
}

double ent_calcMaskedJointEntropyFloat( MutualInfoContext* mc,
					 const float* img1,
					 const float* img2,
					 const int* mask,
					 long dx, long dy, long dz, 
					 long stride1, long stride2,
					 long maskstride )
{
  long a,b,c,d,e,f,g;
  
  long pixels= dx*dy*dz;
  double min1;
  double min2;
  double range1;
  double range2;
  double boxRange1;
  double boxRange2;
  double jointShannonEntropy= 0.0;
  long ofInterest= 0;
  long outOfRangeCounts= 0;

  if (stride1<1) Abort("calcMaskedJointEntropyFloat: nonsense stride1!\n");
  if (stride2<1) Abort("calcMaskedJointEntropyFloat: nonsense stride2!\n");
  if (maskstride<1) Abort("calcMaskedJointEntropyFloat: nonsense maskstride!\n");

  /** Count the number of voxels that will be used **/
  ofInterest= scanMaskedRegionFloat(mc->ec1,img1,mask,pixels,stride1,maskstride);
  if (ofInterest==0) return 0.0;
  (void)scanMaskedRegionFloat(mc->ec2,img2,mask,pixels,stride2,maskstride);

  min1= ent_getMin(mc->ec1);
  min2= ent_getMin(mc->ec2);
  range1 = ent_getMax(mc->ec1) - ent_getMin(mc->ec1);
  range2 = ent_getMax(mc->ec2) - ent_getMin(mc->ec2);
  if (mc->debugFlag) 
    fprintf(stderr,"min1= %lg, min2= %lg, range1= %lg, range2= %lg\n",
	    min1,min2,range1,range2);
      
  if (range1 == 0.0 && range2 == 0.0) return 0.0;

  /* set up bins */ 
  if (mc->num_of_boxes==0) 
    mc->num_of_boxes= pickNBins(ofInterest,stride1,NULL);

  boxRange1 = range1 / (mc->num_of_boxes - 1);
  boxRange2 = range2 / (mc->num_of_boxes - 1);
  
  if (mc->num_of_boxes*mc->num_of_boxes > mc->histogramArraySize) {
    if (mc->debugFlag) 
      fprintf(stderr,"Reallocating histogram array from %ld to %ld\n",
	      mc->histogramArraySize,mc->num_of_boxes*mc->num_of_boxes);
    if (mc->histogramArraySize) free(mc->histogram);
    if (!(mc->histogram= 
	  (long*)malloc(mc->num_of_boxes*mc->num_of_boxes*sizeof(long))))
      Abort("ent_calcMaskedJointEntropyFloat: unable to allocate %d bytes!\n",
	    mc->num_of_boxes*mc->num_of_boxes*sizeof(long));
    mc->histogramArraySize= mc->num_of_boxes* mc->num_of_boxes;
  }
 
  /* initialize histogram array */  
  for(b = 0; b < mc->num_of_boxes * mc->num_of_boxes; b++)    
     mc->histogram[b] = 0;
     
  /* fill histogram array */
  e= 0;
  f= 0;
  g= 0;
  while ( e<stride1*pixels ) {
    if(mask[g] != 0) { 
      long bin = (long)((img1[e]-min1)*(1.0/boxRange1));
      long bin2 = (long)((img2[f]-min2)*(1.0/boxRange2));
      if (bin>=0 && bin<mc->num_of_boxes
	  && bin2>=0 && bin2<mc->num_of_boxes) {
	long location = (long)(bin*mc->num_of_boxes) + bin2; 
	mc->histogram[location]++;
      }
      else {
	outOfRangeCounts++;
	if (mc->debugFlag) fprintf(stderr,"miss: bins would be %ld %ld\n",bin,bin2);
      }
    }
    e += stride1;
    f += stride2;
    g += maskstride;
  }
  if (outOfRangeCounts && mc->verboseFlag)
    Warning(1,"calcMaskedJointEntropyFloat: %ld of %ld samples out of binned range!\n",
	    outOfRangeCounts,ofInterest);

  /* sum histogram entries */
  jointShannonEntropy = 0.0;
  
  /* Calculate the entropy using eqn (2) from the Mutual Information 
   * paper
   */ 
  
  for(d = 0; d < mc->num_of_boxes*mc->num_of_boxes; d++)
    {
      if(mc->histogram[d] != 0) 
	{
	  double prob= (1.0/ofInterest) * (double)mc->histogram[d];   
	  jointShannonEntropy += -prob*log(prob); 
	}
      
    }
  jointShannonEntropy /= log(2.0); /* rescale log_10 to log_2 */
  if (mc->debugFlag) fprintf(stderr,"Joint entropy: %lg\n",jointShannonEntropy);
  return jointShannonEntropy;
}

double ent_calcJointEntropyFloat( MutualInfoContext* mc,
				   const float* img1,
				   const float* img2,
				   long dx, long dy, long dz, 
				   long stride1, long stride2 )
{
  long a,b,c,d,e,f;
  
  long pixels= dx*dy*dz;
  double min1;
  double min2;
  double range1;
  double range2;
  double boxRange1;
  double boxRange2;
  double jointShannonEntropy= 0.0;
  long outOfRangeCounts= 0;

  if (stride1<1) Abort("calcJointEntropyFloat: nonsense stride1!\n");
  if (stride2<1) Abort("calcJointEntropyFloat: nonsense stride2!\n");

  scanRegionFloat(mc->ec1, img1, pixels, stride1);
  scanRegionFloat(mc->ec2, img1, pixels, stride2);

  min1= ent_getMin(mc->ec1);
  min2= ent_getMin(mc->ec2);
  range1 = ent_getMax(mc->ec1) - ent_getMin(mc->ec1);
  range2 = ent_getMax(mc->ec2) - ent_getMin(mc->ec2);
  if (mc->debugFlag) 
    fprintf(stderr,"min1= %lg, min2= %lg, range1= %lg, range2= %lg\n",
	    min1,min2,range1,range2);
      
  if (range1 == 0.0 && range2 == 0.0) return 0.0;

  /* set up bins */ 
  if (mc->num_of_boxes==0) 
    mc->num_of_boxes= pickNBins(pixels,stride1,NULL);

  boxRange1 = range1 / (mc->num_of_boxes - 1);
  boxRange2 = range2 / (mc->num_of_boxes - 1);
  
  if (mc->num_of_boxes*mc->num_of_boxes > mc->histogramArraySize) {
    if (mc->debugFlag) 
      fprintf(stderr,"Reallocating histogram array from %ld to %ld\n",
	      mc->histogramArraySize,mc->num_of_boxes*mc->num_of_boxes);
    if (mc->histogramArraySize) free(mc->histogram);
    if (!(mc->histogram= 
	  (long*)malloc(mc->num_of_boxes*mc->num_of_boxes*sizeof(long))))
      Abort("ent_calcJointEntropyFloat: unable to allocate %d bytes!\n",
	    mc->num_of_boxes*mc->num_of_boxes*sizeof(long));
    mc->histogramArraySize= mc->num_of_boxes* mc->num_of_boxes;
  }
 
  /* initialize histogram array */  
  for(b = 0; b < mc->num_of_boxes * mc->num_of_boxes; b++)    
     mc->histogram[b] = 0;
     
  /* fill histogram array */
  f= 0;
  c= 0;
  while (c<stride1*pixels) {
    long bin = (long)((img1[c]-min1)*(1.0/boxRange1));
    long bin2 = (long)((img2[f]-min2)*(1.0/boxRange2));
    if (bin>=0 && bin<mc->num_of_boxes
	&& bin2>=0 && bin2<mc->num_of_boxes) {
      long location = (long)(bin*mc->num_of_boxes) + bin2; 
      mc->histogram[location]++;
    }
    else {
      outOfRangeCounts++;
      if (mc->debugFlag) fprintf(stderr,"miss: bins would be %ld %ld\n",bin,bin2);
    }
    c += stride1;
    f += stride2;
  }
  
  if (outOfRangeCounts && mc->verboseFlag)
    Warning(1,"calcJointEntropyFloat: %ld of %ld samples out of binned range!\n",
	    outOfRangeCounts,pixels);

  /* sum histogram entries */
  jointShannonEntropy = 0.0;
  
  /* Calculate the entropy using eqn (2) from the Mutual Information 
   * paper
   */ 
  for(d = 0; d < mc->num_of_boxes*mc->num_of_boxes; d++)
    {
      if(mc->histogram[d] != 0) 
	{
	  double prob= (1.0/pixels) * (double)mc->histogram[d];   
	  jointShannonEntropy += -prob*log(prob); 
	}
    }
  jointShannonEntropy /= log(2.0); /* rescale log_10 to log_2 */
  if (mc->debugFlag) fprintf(stderr,"Joint entropy: %lg\n",jointShannonEntropy);
  return jointShannonEntropy;
}

double ent_calcMaskedMutualInformationFloat( MutualInfoContext* mc,
					      const float* img1, 
					      const float* img2, 
					      const int* mask, 
					      long dx, long dy, long dz,
					      long stride1, long stride2,
					      long maskstride )
{
  double mask1Entropy = ent_calcMaskedImageEntropyFloat(mc->ec1, 
							 img1, mask, dx, dy, 
							 dz, stride1, maskstride);
  double mask2Entropy = ent_calcMaskedImageEntropyFloat(mc->ec2, 
							 img2, mask, dx, dy, 
							 dz, stride2, maskstride);
  double mutualInformation = 
    mask1Entropy + mask2Entropy 
    - ent_calcMaskedJointEntropyFloat(mc, img1, img2, mask, dx, dy, dz, 
				       stride1, stride2, maskstride);

  return mutualInformation;
}

double ent_calcMutualInformationFloat( MutualInfoContext* mc,
					const float* img1, 
					const float* img2, 
					long dx, long dy, long dz,
					long stride1, long stride2 )
{
  double img1Entropy = ent_calcImageEntropyFloat(mc->ec1, 
						  img1, dx, dy, dz, stride1);
  double img2Entropy = ent_calcImageEntropyFloat(mc->ec2, 
						  img2, dx, dy, dz, stride2);
  double mutualInformation = 
    img1Entropy + img2Entropy - ent_calcJointEntropyFloat(mc, img1, img2,
							   dx, dy, dz, 
							   stride1, stride2);

  return mutualInformation;
}

