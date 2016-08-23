/************************************************************
 *      func_2rows_unblocked_tool.c                         *
 *                                                          *
 *	Copyright (c) 2007 Pittsburgh Supercomputing Center *
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

/***********************
 * Notes-
 **********************/

static char rcsid[] = "$Id: special_tool.c,v 1.1 2007/06/21 23:32:43 welling Exp $"; 

/*
 * The first output value will be the RMS error; second will be lag
 */
#define ALG_N_OUTPUTS 2
#define MAX_LAG_MAG 100

static double calcRMSError(const double* left, const double* right, 
			   long nIn, int lag)
{
  int i;
  double diff;
  double sum= 0.0;
  if (lag>=0) {
    for (i=0; i<nIn-lag; i++) {
      diff= left[i]-right[i+lag];
      sum += diff*diff;
    }
    return sqrt(sum)/(nIn-lag);
  }
  else {
    lag= -lag;
    for (i=0; i<nIn-lag; i++) {
      diff= left[i+lag]-right[i];
      sum += diff*diff;
    }
    return sqrt(sum)/(nIn-lag);
  }
  return -666.0; /* never happens; impossible value */
}

static int test(const double* left, const double* right, long nIn,
		double* out, long nOut,
		void* hook)
{
  int lag;
  int bestLag;
  double bestRMS;

  if (nIn<=MAX_LAG_MAG) {
    fprintf(stderr,"Test routine in special_tool: nIn is too small!\n");
    return 0;
  }

  bestLag= 0;
  bestRMS= calcRMSError(left, right, nIn, bestLag);

  for (lag=1; lag<=MAX_LAG_MAG; lag++) {
    double newRMS= calcRMSError(left, right, nIn, lag);
    if (newRMS<bestRMS) {
      bestRMS= newRMS;
      bestLag= lag;
    }
    newRMS= calcRMSError(left, right, nIn, -lag);
    if (newRMS<bestRMS) {
      bestRMS= newRMS;
      bestLag= -lag;
    }
  }

  out[0]= bestRMS;
  out[1]= (double)bestLag;
  return 1;
}

Tool* createSpecialTool(Arena* arena) {
  Tool* result= createFunc2UnblkTool(arena, 
				     (FUNC2UNBLKCBFUNC)test,
				     ALG_N_OUTPUTS,
				     NULL);
  return result;
}



