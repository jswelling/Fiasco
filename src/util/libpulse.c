/************************************************************
 *                                                          *
 *  libpulse.c                                              *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 2003 Department of Statistics,         *
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
 *  Original programming by Joel Welling January 2003       *
 ************************************************************/
/* This library supplies utilities related to pulse         *
 * sequence calculations in k-space.                        */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "pulse.h"
#include "errors.h"
#include "misc.h"

static char rcsid[] = "$Id: libpulse.c,v 1.3 2003/01/31 02:47:57 welling Exp $";

const char* pul_getEpiReadoutTypeName( const EpiReadoutType gradType )
{
  switch (gradType) {
  case EPI_READOUT_SQUARE: return("Square");
  case EPI_READOUT_TRAPEZOIDAL: return("Trapeziodal");
  default: Abort("pul_getEpiReadoutTypeName: invalid type %d!\n",
		 (int)gradType);
  }
  return NULL;
}

static double trapezoidArea(const double ht1, const double ht2, 
			    const double base)
{
  return 0.5*(ht1+ht2)*base;
}

int pul_calcEpiReadoutCoords( const int npts, double* coordOut,
			      const double rampUp, 
			      const double flatTopTime, 
			      const double rampDown,
			      const double sampDelayTime,
			      const double durationADC, 
			      const double lag,
			      const EpiReadoutType gradType )
{
  switch (gradType) {
  case EPI_READOUT_SQUARE:
    {
      double step= 1.0/((double)(npts-1));
      int i;
      for (i=0; i<npts-1; i++) coordOut[i]= i*step;
      coordOut[npts-1]= 1.0; /* be exact about this one */
    }
    break;
  case EPI_READOUT_TRAPEZOIDAL:
    {
      double totalArea= 0.0;
      double dArea;
      double startPt= sampDelayTime;
      double endPt= sampDelayTime+durationADC+lag;
      double step= durationADC/((double)npts);
      double here;
      double prevHere;
      int i= 0;

      /* We implement this by walking forward through time
       * over the interval for which the ADC is on.
       */

      if (lag>=step)
	Abort("libpulse: calcEpiReadoutCoords: lag %f > sample time %f!\n",
	      lag,step);

      /* First sample is taken at time startPt+lag */
      coordOut[i++]= 0.0;
      here= prevHere= startPt+lag;

      while (i<npts) {
	here += step;
	dArea= 0.0;
	if (here<rampUp) {
	  dArea += trapezoidArea( prevHere/rampUp, here/rampUp, 
				  here-prevHere );
	}
	else if (here<rampUp+flatTopTime) {
	  if (prevHere<rampUp) {
	    /* Crossed front boundary */
	    dArea += trapezoidArea( prevHere/rampUp, 1.0, 
				    rampUp-prevHere );
	    prevHere= rampUp;
	  }
	  dArea += trapezoidArea( 1.0, 1.0, here-prevHere );
	}
	else {
	  if (prevHere<rampUp+flatTopTime) {
	    /* Crossed back boundary */
	    dArea += trapezoidArea( 1.0, 1.0, (rampUp+flatTopTime)-prevHere);
	    prevHere= rampUp+flatTopTime;
	  }
	  dArea += trapezoidArea( 1.0-(prevHere-(rampUp+flatTopTime))/rampDown, 
				1.0-(here-(rampUp+flatTopTime))/rampDown, 
				 here-prevHere );
	}
	totalArea += dArea;
	prevHere= here;
	coordOut[i++]= totalArea;
      }
      /* Rescale to range 0.0 to 1.0 */
      for (i=0; i<npts; i++) coordOut[i] /= totalArea;
    }
    break;
  default:
    Abort("pul_calcEpiReadoutCoords: unknown type %d!\n",(int)gradType);
  }
 
  return 1;
}
