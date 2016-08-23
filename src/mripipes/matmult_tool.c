/************************************************************
 *      matmult_tool.c                               *
 *                                                          *
 *	Copyright (c) 2004 Pittsburgh Supercomputing Center *
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
 *	History:                                            *
 *		5/04: Written by Joel Welling               *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mri.h>
#include <fmri.h>
#include <stdcrg.h>
#include <kvhash.h>
#include <mripipes.h>
#include "../fmri/lapack.h"

/***********************
 * Notes-
 **********************/

static char rcsid[] = "$Id: matmult_tool.c,v 1.8 2005/08/09 23:11:06 welling Exp $"; 

#define MAX_AT_ONCE (64*1024*1024)

/* Some abbreviations for convenience */
#define LEFT_IN(self) (self->sinkArray[0]->source)
#define RIGHT_IN(self) (self->sinkArray[1]->source)

typedef struct mm_data_struct {
  int complexFlag;
  char* foreachDims;
  int summedOverDim;
  long summed_extent;
  long left_fast_blksize;
  long long left_slow_blksize;
  long right_fast_blksize;
  long long right_slow_blksize;
  long long foreach_blksize;
  long long leftPrevOffset;
  long long rightPrevOffset;
  long leftInbufSize;
  double* leftInbuf;
  long rightInbufSize;
  double* rightInbuf;
  long obufSize;
  long long obufOffset;
  int obufValid;
  double* obuf;
  void (*goMethod)(Tool* self, long long offset);
} MMData;

static void decomposeOffset( Tool* self, const long long outOffset,
			     long long* leftOffset, long long* rightOffset)
{
  MMData* data= (MMData*)self->hook;
  long long offset= outOffset;
  long long foreachOffset;
  long long leftSlowOffset;
  long long rightSlowOffset;
  long long leftFastOffset;
  long long leftFastSlow= data->left_fast_blksize*data->left_slow_blksize;

  /* The output data has structure 
   * left_fast/right_slow/left_slow/foreach_blksize; it
   * is into that array that outOffset indexes.
   *
   * rightFastOffset is always zero, because that is the summed-over
   * dimension.
   */
  foreachOffset= offset/(leftFastSlow*data->right_slow_blksize);
  offset -= foreachOffset*(leftFastSlow*data->right_slow_blksize);
  leftSlowOffset= offset/(data->left_fast_blksize*data->right_slow_blksize);
  offset -= leftSlowOffset*(data->left_fast_blksize*data->right_slow_blksize);
  rightSlowOffset= offset/(data->left_fast_blksize);
  offset -= rightSlowOffset*data->left_fast_blksize;
  leftFastOffset= offset;

  *leftOffset= 
    (foreachOffset*data->left_slow_blksize 
     + leftSlowOffset)*data->left_fast_blksize*data->summed_extent
    + leftFastOffset;
  *rightOffset= 
    (foreachOffset*data->right_slow_blksize 
     + rightSlowOffset)*data->right_fast_blksize*data->summed_extent
    + 0; /* rightFastOffset is always zero */

  if (self->debug) {
    fprintf(stderr,"decomposeOffset: foreach %lld, leftSlow %lld, rightSlow %lld, leftFast %lld\n",
	    foreachOffset,leftSlowOffset,rightSlowOffset,leftFastOffset);
    fprintf(stderr,"decomposeOffset: %lld -> %lld %lld\n",
	    outOffset, *leftOffset, *rightOffset);
  }

  if (leftOffset<0 || rightOffset<0 || foreachOffset>data->foreach_blksize)
    Abort("%s: failed to decompose offset %lld; is it out of range?\n",
	  self->typeName, outOffset);
}

static void mult_general(Tool* self, long long outOffset)
{
  MMData* data= (MMData*)self->hook;
  long long left_base_offset;
  long long left_offset;
  long long right_base_offset;
  long left_fast;
  long summed;

  if (self->debug) fprintf(stderr,"%s: General method used.\n",
			   self->typeName);
  
  decomposeOffset(self, outOffset, &left_base_offset, &right_base_offset);
  if (right_base_offset != data->rightPrevOffset) {
    if (self->debug)
      fprintf(stderr,"%s: Reading %d from right at %lld\n", 
	      self->typeName,data->rightInbufSize,right_base_offset);
    forceGetAllFloat64(RIGHT_IN(self), data->rightInbufSize, 
		       right_base_offset, data->rightInbuf);
    data->rightPrevOffset= right_base_offset;
  }
  for (left_fast=0; left_fast<data->left_fast_blksize; left_fast++)
    data->obuf[left_fast]= 0.0;
  left_offset= left_base_offset;
  for (summed=0; summed<data->summed_extent; summed++) {
    if (left_offset != data->leftPrevOffset) {
      if (self->debug)
	fprintf(stderr,"%s: Reading %d from left at %lld\n",self->typeName,
		data->leftInbufSize, left_offset);
      forceGetAllFloat64(LEFT_IN(self), data->leftInbufSize, left_offset,
			 data->leftInbuf);
      data->leftPrevOffset= left_offset;
    }
    for (left_fast=0; left_fast<data->left_fast_blksize; left_fast++) {
      data->obuf[left_fast] +=
	data->rightInbuf[summed]*data->leftInbuf[left_fast];
    }
    left_offset += data->left_fast_blksize;
  }
  data->obufValid= 1;
  if (self->debug) 
    fprintf(stderr,"%s: Writing block to %lld (size %ld)\n",
	    self->typeName,outOffset,data->obufSize);
}

static void mult_general_preread(Tool* self, long long outOffset)
{
  MMData* data= (MMData*)self->hook;
  long long left_base_offset;
  long long left_offset;
  long long right_base_offset;
  long long right_offset;
  long left_fast;
  long summed;

  if (self->debug) fprintf(stderr,"%s: General_preread method used.\n",
			   self->typeName);
  
  decomposeOffset(self, outOffset, &left_base_offset, &right_offset);
  right_base_offset = right_offset -
    (right_offset % (data->right_slow_blksize*data->summed_extent));
  if (right_base_offset != data->rightPrevOffset) {
    if (self->debug)
      fprintf(stderr,"%s: Reading %d from right at %lld\n", 
	      self->typeName,data->rightInbufSize,right_base_offset);
    forceGetAllFloat64(RIGHT_IN(self), data->rightInbufSize, 
		       right_base_offset, data->rightInbuf);
    data->rightPrevOffset= right_base_offset;
  }
  for (left_fast=0; left_fast<data->left_fast_blksize; left_fast++)
    data->obuf[left_fast]= 0.0;
  left_offset= left_base_offset;
  for (summed=0; summed<data->summed_extent; summed++) {
    if (left_offset != data->leftPrevOffset) {
      if (self->debug)
	fprintf(stderr,"%s: Reading %d from left at %lld\n",self->typeName,
		data->leftInbufSize, left_offset);
      forceGetAllFloat64(LEFT_IN(self), data->leftInbufSize, left_offset,
			 data->leftInbuf);
      data->leftPrevOffset= left_offset;
    }
    for (left_fast=0; left_fast<data->left_fast_blksize; left_fast++) {
      data->obuf[left_fast] +=
	data->rightInbuf[summed+(right_offset-right_base_offset)]
	*data->leftInbuf[left_fast];
    }
    left_offset += data->left_fast_blksize;
  }
  data->obufValid= 1;
  if (self->debug) 
    fprintf(stderr,"%s: Writing block to %lld (size %ld)\n",
	    self->typeName,outOffset,data->obufSize);
}

static void mult_leftsmall(Tool* self, long long outOffset)
{
  MMData* data= (MMData*)self->hook;
  long long left_base_offset;
  long long right_base_offset;
  long left_fast;
  long summed;

  if (self->debug) fprintf(stderr,"%s: Leftsmall method used.\n",
			   self->typeName);
  
  decomposeOffset(self, outOffset, &left_base_offset, &right_base_offset);
  if (left_base_offset != data->leftPrevOffset) {
    if (self->debug)
      fprintf(stderr,"%s: Reading %d from left at %lld\n",self->typeName,
	      data->leftInbufSize, left_base_offset);
    forceGetAllFloat64(LEFT_IN(self), data->leftInbufSize, left_base_offset,
		       data->leftInbuf);
    data->leftPrevOffset= left_base_offset;
  }
  if (right_base_offset != data->rightPrevOffset) {
    if (self->debug)
      fprintf(stderr,"%s: Reading %d from right at %lld\n", 
	      self->typeName,data->rightInbufSize,right_base_offset);
    forceGetAllFloat64(RIGHT_IN(self), data->rightInbufSize, 
		       right_base_offset, data->rightInbuf);
    data->rightPrevOffset= right_base_offset;
  }
  for (left_fast=0; left_fast<data->left_fast_blksize; left_fast++)
    data->obuf[left_fast]= 0.0;
  for (summed=0; summed<data->summed_extent; summed++) {
    for (left_fast=0; left_fast<data->left_fast_blksize; left_fast++) {
      data->obuf[left_fast] +=
	data->rightInbuf[summed]
	*data->leftInbuf[(data->left_fast_blksize)*summed + left_fast];
    }
  }
  data->obufValid= 1;
  if (self->debug) 
    fprintf(stderr,"%s: Writing block to %lld (size %ld)\n",
	    self->typeName,outOffset,data->obufSize);
}

static void mult_leftsmall_preread(Tool* self, long long outOffset)
{
  MMData* data= (MMData*)self->hook;
  long long left_base_offset;
  long long right_offset;
  long long right_base_offset;
  long left_fast;
  long summed;

  if (self->debug) fprintf(stderr,"%s: Leftsmall_preread method used.\n",
			   self->typeName);
  
  decomposeOffset(self, outOffset, &left_base_offset, &right_offset);
  right_base_offset = right_offset -
    (right_offset % (data->right_slow_blksize*data->summed_extent));
  if (left_base_offset != data->leftPrevOffset) {
    if (self->debug)
      fprintf(stderr,"%s: Reading %d from left at %lld\n",self->typeName,
	      data->leftInbufSize, left_base_offset);
    forceGetAllFloat64(LEFT_IN(self), data->leftInbufSize, left_base_offset,
		       data->leftInbuf);
    data->leftPrevOffset= left_base_offset;
  }
  if (right_base_offset != data->rightPrevOffset) {
    if (self->debug)
      fprintf(stderr,"%s: Reading %d from right at %lld\n", 
	      self->typeName,data->rightInbufSize,right_base_offset);
    forceGetAllFloat64(RIGHT_IN(self), data->rightInbufSize, 
		       right_base_offset, data->rightInbuf);
    data->rightPrevOffset= right_base_offset;
  }
  for (left_fast=0; left_fast<data->left_fast_blksize; left_fast++) {
    data->obuf[left_fast]= 0.0;
  }
  for (summed=0; summed<data->summed_extent; summed++) {
    for (left_fast=0; left_fast<data->left_fast_blksize; left_fast++) {
      data->obuf[left_fast] +=
	data->rightInbuf[summed + (right_offset-right_base_offset)]
	*data->leftInbuf[(data->left_fast_blksize)*summed + left_fast];
    }
  }
  data->obufValid= 1;
  if (self->debug) 
    fprintf(stderr,"%s: Writing block to %lld (size %ld)\n",
	    self->typeName,outOffset,data->obufSize);
}

static void mult_leftsmall_accumsmall(Tool* self, long long outOffset)
{
  MMData* data= (MMData*)self->hook;
  double one= 1.0;
  double zero= 0.0;
  int int_one= 1;
  long long left_base_offset;
  long long right_base_offset;
  long long right_offset;
  long long accum_offset;
  int right_slow;
  int rsb= (int)data->right_slow_blksize;
  int lfb= (int)data->left_fast_blksize;
  int se= (int)data->summed_extent;

  if (self->debug) fprintf(stderr,"%s: Leftsmall_Accumsmall method used.\n",
			   self->typeName);
  
  decomposeOffset(self, outOffset, &left_base_offset, &right_base_offset);
  right_base_offset -= 
    (right_base_offset % (data->right_slow_blksize*data->summed_extent));
  if (left_base_offset != data->leftPrevOffset) {
    if (self->debug)
      fprintf(stderr,"%s: Reading %d from left at %lld\n",self->typeName,
	      data->leftInbufSize, left_base_offset);
    forceGetAllFloat64(LEFT_IN(self), data->leftInbufSize, left_base_offset,
		       data->leftInbuf);
    data->leftPrevOffset= left_base_offset;
  }
  right_offset= right_base_offset;
  accum_offset= 0;
  for (right_slow=0; right_slow<data->right_slow_blksize; right_slow++) { 
    if (right_offset != data->rightPrevOffset) {
      if (self->debug)
	fprintf(stderr,"%s: Reading %d from right at %lld\n", 
		self->typeName,data->summed_extent,right_offset);
      forceGetAllFloat64(RIGHT_IN(self), data->rightInbufSize, 
			 right_offset, data->rightInbuf);
      data->rightPrevOffset= right_offset;
    }
    DGEMV( "n", &lfb, &se, 
	   &one, 
	   data->leftInbuf, &lfb, 
	   data->rightInbuf, &int_one,
	   &zero,
	   data->obuf+accum_offset, &int_one );
    accum_offset += data->left_fast_blksize;
    right_offset += data->summed_extent;
  }
  data->obufValid= 1;
  if (self->debug) 
    fprintf(stderr,"%s: Writing block to %lld (size %ld)\n",
	    self->typeName,outOffset,data->obufSize);
}

static void mult_leftsmall_accumsmall_preread(Tool* self, long long outOffset)
{
  MMData* data= (MMData*)self->hook;
  double one= 1.0;
  double zero= 0.0;
  int int_one= 1;
  long long left_base_offset;
  long long right_base_offset;
  int rsb= (int)data->right_slow_blksize;
  int lfb= (int)data->left_fast_blksize;
  int se= (int)data->summed_extent;

  if (self->debug) 
    fprintf(stderr,"%s: Leftsmall_Accumsmall_preread method used.\n",
	    self->typeName);
  
  decomposeOffset(self, outOffset, &left_base_offset, &right_base_offset);
  if (left_base_offset != data->leftPrevOffset) {
    if (self->debug)
      fprintf(stderr,"%s: Reading %d from left at %lld\n", self->typeName,
	      data->leftInbufSize, left_base_offset);
    forceGetAllFloat64(LEFT_IN(self), data->leftInbufSize,
		       left_base_offset, data->leftInbuf);
    data->leftPrevOffset= left_base_offset;
  }
  if (right_base_offset != data->rightPrevOffset) {
    if (self->debug)
      fprintf(stderr,"%s: Reading %d from right at %lld\n", 
	      self->typeName,data->rightInbufSize,right_base_offset);
    forceGetAllFloat64(RIGHT_IN(self), data->rightInbufSize,
		       right_base_offset, data->rightInbuf);
    data->rightPrevOffset= right_base_offset;
  }

  DGEMM( "n", "n", 
	 &lfb, &rsb, &se, 
	 &one, 
	 data->leftInbuf, &lfb, 
	 data->rightInbuf, &se, 
	 &zero, data->obuf, &lfb );
  data->obufValid= 1;
  if (self->debug) 
    fprintf(stderr,"%s: Writing block to %lld (size %ld)\n",
	    self->typeName,outOffset, data->obufSize);
}

static void mult_general_complex(Tool* self, long long outOffset)
{
  MMData* data= (MMData*)self->hook;
  long long left_base_offset;
  long long left_offset= 0;
  long long right_base_offset;
  long left_fast;
  long summed;

  if (self->debug) fprintf(stderr,"%s: General_complex method used.\n",
			   self->typeName);
  
  decomposeOffset(self, outOffset, &left_base_offset, &right_base_offset);
  if (right_base_offset != data->rightPrevOffset) {
    if (self->debug)
      fprintf(stderr,"%s: Reading %d from right at %lld\n", 
	      self->typeName,data->rightInbufSize,right_base_offset);
    forceGetAllFloat64(RIGHT_IN(self), data->rightInbufSize, 
		       right_base_offset, data->rightInbuf);
    data->rightPrevOffset= right_base_offset;
  }
  for (left_fast=0; left_fast<data->left_fast_blksize; left_fast += 2) {
    data->obuf[left_fast+1]= data->obuf[left_fast]= 0.0;
  }
  for (summed=0; summed<data->summed_extent; summed++) {
    double rR= data->rightInbuf[2*summed];
    double rI= data->rightInbuf[2*summed+1];
    if (left_offset != data->leftPrevOffset) {
      if (self->debug)
	fprintf(stderr,"%s: Reading %d from left at %lld\n",self->typeName,
		data->leftInbufSize, left_offset);
      forceGetAllFloat64(LEFT_IN(self), data->leftInbufSize, left_offset,
			 data->leftInbuf);
      data->leftPrevOffset= left_offset;
    }
    for (left_fast=0; left_fast<data->left_fast_blksize; left_fast += 2) {
      double lR= data->leftInbuf[left_fast];
      double lI= data->leftInbuf[left_fast +1];
      data->obuf[left_fast] += rR*lR - rI*lI;
      data->obuf[left_fast+1] += rR*lI + rI*lR;
    }
    left_offset += data->left_fast_blksize;
  }
  data->obufValid= 1;
  if (self->debug) 
    fprintf(stderr,"%s: Writing block to %lld (size %ld)\n",
	    self->typeName,outOffset,data->obufSize);
}

static void mult_general_preread_complex(Tool* self, long long outOffset)
{
  MMData* data= (MMData*)self->hook;
  long long left_base_offset;
  long long left_offset;
  long long right_base_offset;
  long long right_offset;
  long left_fast;
  long summed;

  if (self->debug) fprintf(stderr,"%s: General_preread_complex method used.\n",
			   self->typeName);
  
  decomposeOffset(self, outOffset, &left_base_offset, &right_offset);
  right_base_offset = right_offset -
    (right_offset % (2*data->right_slow_blksize*data->summed_extent));
  if (right_base_offset != data->rightPrevOffset) {
    if (self->debug)
      fprintf(stderr,"%s: Reading %d from right at %lld\n", 
	      self->typeName,data->rightInbufSize,right_base_offset);
    forceGetAllFloat64(RIGHT_IN(self), data->rightInbufSize, 
		       right_base_offset, data->rightInbuf);
    data->rightPrevOffset= right_base_offset;
  }
  left_offset= left_base_offset;
  for (left_fast=0; left_fast<data->left_fast_blksize; left_fast += 2) {
    data->obuf[left_fast+1]= data->obuf[left_fast]= 0.0;
  }
  for (summed=0; summed<data->summed_extent; summed++) {
    long rOff= 2*summed+(right_offset - right_base_offset);
    double rR= data->rightInbuf[rOff];
    double rI= data->rightInbuf[rOff+1];
    if (left_offset != data->leftPrevOffset) {
      if (self->debug)
	fprintf(stderr,"%s: Reading %d from left at %lld\n",self->typeName,
		data->leftInbufSize, left_offset);
      forceGetAllFloat64(LEFT_IN(self), data->leftInbufSize, left_offset,
			 data->leftInbuf);
      data->leftPrevOffset= left_offset;
    }
    for (left_fast=0; left_fast<data->left_fast_blksize; left_fast += 2) {
      double lR= data->leftInbuf[left_fast];
      double lI= data->leftInbuf[left_fast +1];
      data->obuf[left_fast] += rR*lR - rI*lI;
      data->obuf[left_fast+1] += rR*lI + rI*lR;
    }
    left_offset += data->left_fast_blksize;
  }
  data->obufValid= 1;
  if (self->debug) 
    fprintf(stderr,"%s: Writing block to %lld (size %ld)\n",
	    self->typeName,outOffset,data->obufSize);
}

static void mult_leftsmall_complex(Tool* self, long long outOffset)
{
  MMData* data= (MMData*)self->hook;
  long long left_base_offset;
  long long right_base_offset;
  long left_fast;
  long summed;

  if (self->debug) fprintf(stderr,"%s: Leftsmall_complex method used.\n",
			   self->typeName);
  
  decomposeOffset(self, outOffset, &left_base_offset, &right_base_offset);
  if (left_base_offset != data->leftPrevOffset) {
    if (self->debug)
      fprintf(stderr,"%s: Reading %d from left at %lld\n",self->typeName,
	      data->leftInbufSize, left_base_offset);
    forceGetAllFloat64(LEFT_IN(self), data->leftInbufSize, left_base_offset,
		       data->leftInbuf);
    data->leftPrevOffset= left_base_offset;
  }
  if (right_base_offset != data->rightPrevOffset) {
    if (self->debug)
      fprintf(stderr,"%s: Reading %d from right at %lld\n", 
	      self->typeName,data->rightInbufSize,right_base_offset);
    forceGetAllFloat64(RIGHT_IN(self), data->rightInbufSize, 
		       right_base_offset, data->rightInbuf);
    data->rightPrevOffset= right_base_offset;
  }
  for (left_fast=0; left_fast<data->left_fast_blksize; left_fast += 2) {
    data->obuf[left_fast+1]= data->obuf[left_fast]= 0.0;
  }
  for (summed=0; summed<data->summed_extent; summed++) {
    double rR= data->rightInbuf[2*summed];
    double rI= data->rightInbuf[2*summed+1];
    for (left_fast=0; left_fast<data->left_fast_blksize; left_fast += 2) {
      double lR= data->leftInbuf[data->left_fast_blksize*summed + left_fast];
      double lI= data->leftInbuf[data->left_fast_blksize*summed +left_fast +1];
      data->obuf[left_fast] += rR*lR - rI*lI;
      data->obuf[left_fast+1] += rR*lI + rI*lR;
    }
  }
  data->obufValid= 1;
  if (self->debug) 
    fprintf(stderr,"%s: Writing block to %lld (size %ld)\n",
	    self->typeName,outOffset,data->obufSize);
}

static void mult_leftsmall_preread_complex(Tool* self, long long outOffset)
{
  MMData* data= (MMData*)self->hook;
  long long left_base_offset;
  long long right_offset;
  long long right_base_offset;
  long left_fast;
  long summed;

  if (self->debug) fprintf(stderr,"%s: Leftsmall_preread_complex method used.\n",
			   self->typeName);
  
  decomposeOffset(self, outOffset, &left_base_offset, &right_offset);
  right_base_offset = right_offset -
    (right_offset % (2*data->right_slow_blksize*data->summed_extent));
  if (left_base_offset != data->leftPrevOffset) {
    if (self->debug)
      fprintf(stderr,"%s: Reading %d from left at %lld\n",self->typeName,
	      data->leftInbufSize, left_base_offset);
    forceGetAllFloat64(LEFT_IN(self), data->leftInbufSize, left_base_offset,
		       data->leftInbuf);
    data->leftPrevOffset= left_base_offset;
  }
  if (right_base_offset != data->rightPrevOffset) {
    if (self->debug)
      fprintf(stderr,"%s: Reading %d from right at %lld\n", 
	      self->typeName,data->rightInbufSize,right_base_offset);
    forceGetAllFloat64(RIGHT_IN(self), data->rightInbufSize, 
		       right_base_offset, data->rightInbuf);
    data->rightPrevOffset= right_base_offset;
  }
  for (left_fast=0; left_fast<data->left_fast_blksize; left_fast += 2) {
    data->obuf[left_fast+1]= data->obuf[left_fast]= 0.0;
  }
  for (summed=0; summed<data->summed_extent; summed++) {
    long rOff= 2*summed+(right_offset - right_base_offset);
    double rR= data->rightInbuf[rOff];
    double rI= data->rightInbuf[rOff+1];
    for (left_fast=0; left_fast<data->left_fast_blksize; left_fast += 2) {
      double lR= data->leftInbuf[data->left_fast_blksize*summed + left_fast];
      double lI= data->leftInbuf[data->left_fast_blksize*summed +left_fast +1];
      data->obuf[left_fast] += rR*lR - rI*lI;
      data->obuf[left_fast+1] += rR*lI + rI*lR;
    }
  }
  data->obufValid= 1;
  if (self->debug) 
    fprintf(stderr,"%s: Writing block to %lld (size %ld)\n",
	    self->typeName,outOffset,data->obufSize);
}

static void mult_leftsmall_accumsmall_complex(Tool* self, long long outOffset)
{
  MMData* data= (MMData*)self->hook;
  double c_one[2]= {1.0, 0.0};
  double c_zero[2]= {0.0, 0.0};
  int int_one= 1;
  long long left_base_offset;
  long long right_base_offset;
  long long right_offset;
  long long accum_offset;
  int right_slow;
  int rsb= (int)data->right_slow_blksize;
  int hlfb= (int)(data->left_fast_blksize/2);
  int se= (int)data->summed_extent;

  if (self->debug) fprintf(stderr,"%s: Leftsmall_Accumsmall_complex method used.\n",
			   self->typeName);
  
  decomposeOffset(self, outOffset, &left_base_offset, &right_base_offset);
  right_base_offset -= 
    (right_base_offset % (data->right_slow_blksize*data->summed_extent));
  if (left_base_offset != data->leftPrevOffset) {
    if (self->debug)
      fprintf(stderr,"%s: Reading %d from left at %lld\n",self->typeName,
	      data->leftInbufSize, left_base_offset);
    forceGetAllFloat64(LEFT_IN(self), data->leftInbufSize, left_base_offset,
		       data->leftInbuf);
    data->leftPrevOffset= left_base_offset;
  }
  right_offset= right_base_offset;
  accum_offset= 0;
  for (right_slow=0; right_slow<data->right_slow_blksize; right_slow++) { 
    if (right_offset != data->rightPrevOffset) {
      if (self->debug)
	fprintf(stderr,"%s: Reading %d from right at %lld\n", 
		self->typeName,data->rightInbufSize,right_offset);
      forceGetAllFloat64(RIGHT_IN(self), data->rightInbufSize, 
			 right_offset, data->rightInbuf);
      data->rightPrevOffset= right_offset;
    }
    ZGEMV( "n", &hlfb, &se, 
	   c_one, 
	   data->leftInbuf, &hlfb, 
	   data->rightInbuf, &int_one,
	   c_zero,
	   data->obuf+accum_offset, &int_one );
    accum_offset += data->left_fast_blksize;
    right_offset += 2*data->summed_extent;
  }
  data->obufValid= 1;
  if (self->debug) 
    fprintf(stderr,"%s: Writing block to %lld (size %ld)\n",
	    self->typeName,outOffset,data->obufSize);
}

static void mult_leftsmall_accumsmall_preread_complex(Tool* self, 
						      long long outOffset)
{
  MMData* data= (MMData*)self->hook;
  double c_one[2]= {1.0,0.0};
  double c_zero[2]= {0.0,0.0};
  int int_one= 1;
  long long left_base_offset;
  long long right_base_offset;
  int rsb= (int)data->right_slow_blksize;
  int hlfb= (int)(data->left_fast_blksize/2); /* allow for complex */
  int se= (int)data->summed_extent;

  if (self->debug) 
    fprintf(stderr,"%s: Leftsmall_Accumsmall_preread_complex method used.\n",
	    self->typeName);
  
  decomposeOffset(self, outOffset, &left_base_offset, &right_base_offset);
  if (left_base_offset != data->leftPrevOffset) {
    if (self->debug)
      fprintf(stderr,"%s: Reading %d from left at %lld\n", self->typeName,
	      data->leftInbufSize, left_base_offset);
    forceGetAllFloat64(LEFT_IN(self), data->leftInbufSize,
		       left_base_offset, data->leftInbuf);
    data->leftPrevOffset= left_base_offset;
  }
  if (right_base_offset != data->rightPrevOffset) {
    if (self->debug)
      fprintf(stderr,"%s: Reading %d from right at %lld\n", 
	      self->typeName,data->rightInbufSize,right_base_offset);
    forceGetAllFloat64(RIGHT_IN(self), data->rightInbufSize,
		       right_base_offset, data->rightInbuf);
    data->rightPrevOffset= right_base_offset;
  }

  ZGEMM( "n", "n", 
	 &hlfb, &rsb, &se, 
	 c_one, 
	 data->leftInbuf, &hlfb, 
	 data->rightInbuf, &se, 
	 c_zero, data->obuf, &hlfb );
  data->obufValid= 1;
  if (self->debug) 
    fprintf(stderr,"%s: Writing block to %lld (size %ld)\n",
	    self->typeName,outOffset, data->obufSize);
}

static void calcBlockSizes( Tool* self, DataSource* Left, DataSource* Right,
			    DataSource* Out )
{
  MMData* data= (MMData*)self->hook;
  char* dimstr1_orig= strdup(getSourceDims(Left));
  char* dimstr1= dimstr1_orig;
  char* dimstr2_orig= strdup(getSourceDims(Right));
  char* dimstr2= dimstr2_orig;
  
  data->summed_extent= getSourceDimExtent(Left, data->summedOverDim);

  data->foreach_blksize= 1;
  if (*(data->foreachDims)) {
    char* here= strstr(dimstr1,data->foreachDims);
    *here= '\0';
    here= strstr(dimstr2,data->foreachDims);
    *here= '\0';
    here= (char*)data->foreachDims;
    while (*here) data->foreach_blksize *= getSourceDimExtent(Left, *here++);
  }
  
  calcSourceBlockSizes(Left, dimstr1, data->summedOverDim, 
	     &(data->left_fast_blksize), &(data->left_slow_blksize));
  calcSourceBlockSizes(Right, dimstr2, data->summedOverDim, 
	     &(data->right_fast_blksize), &(data->right_slow_blksize));
  if (self->debug) {
    fprintf(stderr,"Left string %s, block sizes: %ld %ld %lld\n",
	    dimstr1,data->left_fast_blksize, data->summed_extent, 
	    data->left_slow_blksize);
    fprintf(stderr,"Right string %s, block sizes: %ld %ld %lld\n",
	    dimstr2,data->right_fast_blksize, data->summed_extent, 
	    data->right_slow_blksize);
    fprintf(stderr,"Foreach string %s, block size %lld\n", 
	    data->foreachDims,data->foreach_blksize );
    fprintf(stderr,"Output dims %s\n",getSourceDims(Out));
  }

  free(dimstr1_orig);
  free(dimstr2_orig);  
}

static void restructureDims( Tool* self, DataSource* Output,
			     DataSource* Factor1, DataSource* Factor2,
			     const int complex, const char* foreach_dims )
{
  /* This routine changes the dimensions of the output to be those
   * appropriate for the product of the inputs.  structure_check()
   * has already verified that the input dataset have the necessary
   * structure for these steps to work.
   */
  char buf[256];
  char key_buf[32];
  char* dimstr1_orig= strdup(getSourceDims(Factor1));
  char* dimstr1= dimstr1_orig;
  char* dimstr2_orig= strdup(getSourceDims(Factor2));
  char* dimstr2= dimstr2_orig;
  int summed_dim; 
  char* work1= strdup(dimstr1);
  char* here= NULL;
  char* there= NULL;
  char* end;

  if (complex) { /* skip the v dimension */
    dimstr1++;
    dimstr2++;
  }

  summed_dim= dimstr2[0];
  here= strchr(work1,summed_dim);
  if (foreach_dims[0]) {
    there= strstr(dimstr1+1,foreach_dims);
    *there= '\0';
    there= strstr(dimstr2+1,foreach_dims);
    *there= '\0';
  }
  /* here now points to the live part of dimstr2, excluding the summed dim,
   * the foreach dims if any, and the v dim if any.
   */

  if (self->verbose) {
    if (complex) {
      Message("Multiplying complex input source<%s>, dims v%s%s and v%s%s.\n",
	      Factor1->name,dimstr1,foreach_dims,dimstr2,foreach_dims);
    }
    else {
      Message("Multiplying input source <%s>, dims %s%s and %s%s.\n",
	      Factor1->name,dimstr1,foreach_dims,dimstr2,foreach_dims);
    }
    Message("Summing over dimension %c, looping over dimensions <%s>\n",
	    summed_dim, foreach_dims);
  }

  /* Break the left dimension string */
  *here= '\0';
  here++; /* now points to left slow dimensions + foreach dims */

  /* Define the new dimension string */
  buf[sizeof(buf)-1]= '\0';
  strncpy(buf,work1,sizeof(buf)-1);
  strncat(buf,dimstr2+1,sizeof(buf)-1);
  strncat(buf,here,sizeof(buf)-1);
  kvDefString(Output->attr,"dimensions",buf);
  if (self->verbose) 
    Message("Output dimension string is <%s>\n",buf);
  
  /* Delete the summed-over extent */
  snprintf(key_buf,sizeof(key_buf),"extent.%c",summed_dim);
  kvDeleteAll(Output->attr,key_buf);

  /* Define the undefined extents */
  for (here= dimstr2+1; *here; here++) {
    snprintf(key_buf,sizeof(key_buf),"extent.%c",*here);
    snprintf(buf,sizeof(buf),"%d",getSourceDimExtent(Factor2,*here));
    kvDefString(Output->attr,key_buf,buf);
  }
  
  free(work1);
  free(dimstr1_orig);
  free(dimstr2_orig);
}

static int structureCheck( Tool* self, DataSource* Left, DataSource* Right,
			   const int complex,
			   int* summed_dim_ptr, char** foreach_dims_ptr)
{
  char* dimstr1_orig= NULL;
  char* dimstr1= NULL;
  char* dimstr2_orig= NULL;
  char* dimstr2= NULL;
  int summed_dim= '\0';
  char* here;
  char* there;
  char* foreach_dims;

  dimstr1_orig= strdup(getSourceDims(Left));
  dimstr1= dimstr1_orig;
  dimstr2_orig= strdup(getSourceDims(Right));
  dimstr2= dimstr2_orig;

  /* If this is to be a complex matrix multiplication, both
   * inputs should have first input dimensions v with extent 2.
   */
  if (complex) {
    if (!kvLookup(Left->attr,"extent.v"))
      Abort("Left input to %s has no extent.v info!\n",self->typeName);
    if (!kvLookup(Right->attr,"extent.v"))
      Abort("Right input to %s has no extent.v info!\n",self->typeName);
    if (dimstr1[0] != 'v' || getSourceDimExtent(Left,'v')!=2 )
      Abort("%s: complex multiplication requested, but first input is not complex!\n",
	    self->typeName);
    if (dimstr2[0] != 'v' || getSourceDimExtent(Right,'v')!=2 )
      Abort("%s: complex multiplication requested, but second input is not complex!\n",
	    self->typeName);
    dimstr1++;
    dimstr2++;
  }

  /* Find any looped-over dimensions and clip them from the dim strings */
  here= dimstr1 + strlen(dimstr1) - 1;
  there= dimstr2 + strlen(dimstr2) - 1;
  while (*here==*there && here>dimstr1 && there>dimstr2) {
    if (getSourceDimExtent(Left,*here) != getSourceDimExtent(Right,*here))
      Abort("%s: extents of looped-over dimension %c do not match!\n",
	    self->typeName, *here);
    here--;
    there--;
  }
  foreach_dims= strdup(here+1);
  *(here+1)= *(there+1)= '\0';

  /* In remaining dim string, there must be exactly one dimension in common */
  for (here=dimstr1; *here; here++) {
    if ((there= strchr(dimstr2, *here)) != NULL) {
      if (summed_dim== '\0') summed_dim= *here;
      else Abort("%s: more than one duplicated dimension in input chunks!\n",
		 self->typeName);
    }
  }

  /* It must be the left-most dimension of the right-hand dataset */
  if (dimstr2[0] != summed_dim)
    Abort("%s: Summed dimension %c must be left-most in input dataset #2!\n",
	  self->typeName, summed_dim);

  /* This dimension must have the same extent in both factors */
  if (getSourceDimExtent(Left,summed_dim) 
      != getSourceDimExtent(Right,summed_dim))
    Abort("%s: The inputs have different extents for summed dimension %c!\n",
	  self->typeName, summed_dim);

  *summed_dim_ptr= summed_dim;
  *foreach_dims_ptr= foreach_dims;
  free(dimstr1_orig);
  free(dimstr2_orig);
  return 1;
}

static int selectMultMethod(Tool* self)
{
  MMData* data= (MMData*)self->hook;
  long accumSize= 0;
  long leftSize= 0;
  long rightSize= 0;

  /* Pick an output buffer size and multiplication method */
  if (data->complexFlag) {

    int preread= 
      ( 2*data->right_slow_blksize*data->summed_extent <= MAX_AT_ONCE );
    if ( preread ) rightSize= 2*data->right_slow_blksize*data->summed_extent;
    else rightSize= 2*data->summed_extent;

    if (data->left_fast_blksize*data->summed_extent <= MAX_AT_ONCE) {
      leftSize= data->left_fast_blksize*data->summed_extent;

      if (data->left_fast_blksize*data->right_slow_blksize <= MAX_AT_ONCE) {
	accumSize= data->left_fast_blksize*data->right_slow_blksize;
	if (preread) data->goMethod= mult_leftsmall_accumsmall_preread_complex;
	else data->goMethod= mult_leftsmall_accumsmall_complex;
      }
      else {
	accumSize= data->left_fast_blksize;
	if (preread) data->goMethod= mult_leftsmall_preread_complex;
	else data->goMethod= mult_leftsmall_complex;
      }

    }
    else {
      accumSize= data->left_fast_blksize;
      leftSize= data->left_fast_blksize;
      if (preread) data->goMethod= mult_general_preread_complex;
      else data->goMethod= mult_general_complex;
    }
  }
  else {

    int preread= 
      ( data->right_slow_blksize*data->summed_extent <= MAX_AT_ONCE );
    if ( preread ) rightSize= data->right_slow_blksize*data->summed_extent;
    else rightSize= data->summed_extent;

    if (data->left_fast_blksize*data->summed_extent <= MAX_AT_ONCE) {
      leftSize= data->left_fast_blksize*data->summed_extent;

      if (data->left_fast_blksize*data->right_slow_blksize <= MAX_AT_ONCE) {
	accumSize= data->left_fast_blksize*data->right_slow_blksize;
	if (preread) data->goMethod= mult_leftsmall_accumsmall_preread;
	else data->goMethod= mult_leftsmall_accumsmall;
      }
      else {
	accumSize= data->left_fast_blksize;
	if (preread) data->goMethod= mult_leftsmall_preread;
	else data->goMethod= mult_leftsmall;
      }

    }
    else {
      accumSize= data->left_fast_blksize;
      leftSize= data->left_fast_blksize;
      if (preread) data->goMethod= mult_general_preread;
      else data->goMethod= mult_general;
    }
  }

  if (!(data->obuf=(double*)malloc(accumSize*sizeof(double))))
    Abort("Unable to allocate %d bytes!\n",accumSize*sizeof(double));
  data->obufSize= accumSize;

  if (!(data->leftInbuf=(double*)malloc(leftSize*sizeof(double))))
    Abort("Unable to allocate %d bytes!\n",leftSize*sizeof(double));
  data->leftInbufSize= leftSize;

  if (!(data->rightInbuf=(double*)malloc(rightSize*sizeof(double))))
    Abort("Unable to allocate %d bytes!\n",rightSize*sizeof(double));
  data->rightInbufSize= rightSize;

  return 1;
}
  
static DataSink* createMatmultLeftSink( Tool* owner )
{
  DataSink* result= createBaseSink(owner);
  return result;
}

static DataSink* createMatmultRightSink( Tool* owner )
{
  DataSink* result= createBaseSink(owner);
  return result;
}

static void fillBuffer(Tool* self, long long offset)
{
  MMData* data= (MMData*)self->hook;
  /* If data->obuf is invalid or contains non-applicable data, 
   * calculate the block base offset which is closest to the 
   * requested offset and fill the obuf from that offset.
   */
  if (!data->obufValid || offset<data->obufOffset 
      || offset>=data->obufOffset+data->obufSize) {
    long long blockBaseOffset= offset - (offset % data->obufSize);
    data->obufValid= 0;
    data->goMethod(self, blockBaseOffset);
    data->obufOffset= blockBaseOffset;
  }
} 

static long getFloat32Chunk(DataSource* self,
				long size, long long offset, float* buf)
{
  MMData* data= (MMData*)(self->owner->hook);
  long shift;
  long n;
  long i;

  fillBuffer(self->owner, offset);

  /* The obuf is now guaranteed to have an applicable segment;
   * copy it to the output.
   */
  shift= (long)(offset - data->obufOffset);
  if (size>data->obufSize-shift) n= data->obufSize-shift;
  else n= size;
  for (i=0; i<n; i++) buf[i]= (float)data->obuf[i+shift];
  return n;
}

static long getFloat64Chunk(DataSource* self,
				long size, long long offset, double* buf)
{
  MMData* data= (MMData*)(self->owner->hook);
  long shift;
  long n;
  long i;

  fillBuffer(self->owner, offset);

  /* The obuf is now guaranteed to have an applicable segment;
   * copy it to the output.
   */
  shift= (long)(offset - data->obufOffset);
  if (size>data->obufSize-shift) n= data->obufSize-shift;
  else n= size;
  for (i=0; i<n; i++) buf[i]= data->obuf[i+shift];
  return n;
}

static DataSource* createMatmultSource( Tool* owner )
{
  DataSource* result= createBaseSource(owner);
  result->pGetFloat64Chunk= getFloat64Chunk;
  result->pGetFloat32Chunk= getFloat32Chunk;

  return result;
}

static void destroySelf(Tool* self)
{
  MMData* data= (MMData*)self->hook;
  if (data->foreachDims) free(data->foreachDims);
  if (data->leftInbufSize) free(data->leftInbuf);
  if (data->rightInbufSize) free(data->rightInbuf);
  if (data->obufSize) free(data->obuf);
  baseToolDestroySelf(self);
}

static int init(Tool* self)
{
  DataSource* mySource= self->sourceArray[0];
  DataSource* leftSource= self->sinkArray[0]->source;
  DataSource* rightSource= self->sinkArray[1]->source;
  MMData* data= (MMData*)self->hook;
  if (!baseToolInit(self)) return 0;
  kvCopyUniqueExceptHashes( mySource->attr, leftSource->attr );
  mySource->pSetName(mySource, leftSource->pGetName(leftSource));
  /* Here we edit the attributes of the output to comply with the 
   * matmult rules.
   */
  if (structureCheck(self,leftSource,rightSource,data->complexFlag,
		     &(data->summedOverDim),&(data->foreachDims))) {
    restructureDims(self,mySource,leftSource,rightSource,
		    data->complexFlag,data->foreachDims);
    calcBlockSizes(self, leftSource, rightSource, mySource);
    if (data->complexFlag) {
      if (data->right_fast_blksize != 2)
	Abort("%s: internal error: right fast block is not 2!\n",
	      self->typeName);
    }
    else {
      if (data->right_fast_blksize != 1)
	Abort("%s: internal error: right fast block is not 1!\n",
	      self->typeName);
    }
  }
  else return 0;
  
  /* Choose buffer sizes and methods */
  if (data->leftInbufSize) { free(data->leftInbuf); data->leftInbufSize= 0; }
  if (data->rightInbufSize) { free(data->rightInbuf); data->rightInbufSize= 0;}
  if (data->obufSize) { free(data->obuf); data->obufSize= 0; }
  if (!selectMultMethod(self)) return 0;

  return 1;
}

Tool* createMatmultTool(Arena* arena) {
  Tool* result= createBaseTool(arena);
  MMData* data= NULL;
  if (!(data=(MMData*)malloc(sizeof(MMData))))
    Abort("Unable to allocate %d bytes!\n",sizeof(MMData));
  result->hook= data;
  data->complexFlag= 0;
  data->foreachDims= NULL;
  data->summedOverDim= '\0';
  data->summed_extent= data->left_fast_blksize=
    data->left_slow_blksize= data->right_fast_blksize=
    data->right_slow_blksize= data->foreach_blksize= 0;
  data->leftInbufSize= 0;
  data->leftInbuf= NULL;
  data->rightInbufSize= 0;
  data->rightInbuf= NULL;
  data->obufSize= 0;
  data->obuf= NULL;
  data->obufValid= 0;
  data->goMethod= NULL;
  data->leftPrevOffset= data->rightPrevOffset= -1;
  result->pInit= init;
  result->pAddSink( result, createMatmultLeftSink(result) );
  result->pAddSink( result, createMatmultRightSink(result) );
  result->pAddSource( result, createMatmultSource(result) );
  result->typeName= "MatMult";
  return result;
}

Tool* createComplexMatmultTool(Arena* arena) {
  Tool* result= createBaseTool(arena);
  MMData* data= NULL;
  if (!(data=(MMData*)malloc(sizeof(MMData))))
    Abort("Unable to allocate %d bytes!\n",sizeof(MMData));
  result->hook= data;
  data->complexFlag= 1;
  data->foreachDims= NULL;
  data->summedOverDim= '\0';
  data->summed_extent= data->left_fast_blksize=
    data->left_slow_blksize= data->right_fast_blksize=
    data->right_slow_blksize= data->foreach_blksize= 0;
  data->leftInbufSize= 0;
  data->leftInbuf= NULL;
  data->rightInbufSize= 0;
  data->rightInbuf= NULL;
  data->obufSize= 0;
  data->obuf= NULL;
  data->obufValid= 0;
  data->goMethod= NULL;
  data->leftPrevOffset= data->rightPrevOffset= -1;
  result->pDestroySelf= destroySelf;
  result->pInit= init;
  result->pAddSink( result, createMatmultLeftSink(result) );
  result->pAddSink( result, createMatmultRightSink(result) );
  result->pAddSource( result, createMatmultSource(result) );
  result->typeName= "MatMult";
  return result;
}


