/************************************************************
 *                                                          *
 *  pulsecalc.c                                             *
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
/* This supplies some simple tests of libpulse. */

#include <stdio.h>
#include <strings.h>
#include <string.h>
#include <stdlib.h>
#include "pulse.h"
#include "errors.h"
#include "misc.h"

static char rcsid[] = "$Id: pulsecalc.c,v 1.2 2007/03/21 23:54:26 welling Exp $";

int main( int argc, char* argv[] )
{
  if (!strcasecmp(argv[1],"epireadout")) {
    EpiReadoutType type;
    int npts;
    int i;
    double* result;
    double rampUp;
    double rampDown;
    double flatTopTime;
    double sampDelay;
    double durADC;
    double lag;
    double sampTime;
    
    if (argc != 10) {
      fprintf(stderr,
	      "usage: %s epireadout [square | trapezoid] nsamps rampUp rampFlat rampDn sampDelay durADC lag\n",
	      argv[0]);
      exit(-1);
    }

    npts= atoi(argv[3]);
    rampUp= atof(argv[4]);
    flatTopTime= atof(argv[5]);
    rampDown= atof(argv[6]);
    sampDelay= atof(argv[7]);
    durADC= atof(argv[8]);
    sampTime= durADC/((double)npts);
    if (!strcmp(argv[9],"none")) lag= 0.0;
    else if (!strcmp(argv[9],"half")) lag= 0.5*sampTime;
    else if (!strcmp(argv[9],"full")) lag= sampTime;
    else lag= atof(argv[9]);

    if (!(result=(double*)malloc(npts*sizeof(double))))
      Abort("%s: unable to allocate %d bytes!\n",npts*sizeof(double));

    if (!strcasecmp(argv[2],"square")) type= EPI_READOUT_SQUARE;
    else type= EPI_READOUT_TRAPEZOIDAL;

    printf("# Call was calc <%s> type <%s>, npts %d\n",
	   argv[1],argv[2],npts);
    printf("#    up %f flat %f down %f delay %f durationADC %f lag %f\n",
	   rampUp,flatTopTime,rampDown,sampDelay,durADC,lag);
    printf("# type translation <%s>\n",pul_getEpiReadoutTypeName(type));

    if (!pul_calcEpiReadoutCoords(npts, result, 
				  rampUp, flatTopTime, rampDown, 
				  sampDelay, durADC, lag, type))
      fprintf(stderr,"pul_calcEpiReadoutCoords failed!\n");

    printf("# Scaled k-space sample coords follow: \n");
    for (i=0; i<npts; i++) printf("%d %f\n",i,result[i]);

  }
  else fprintf(stderr,"%s: test type <%s> unrecognized!\n",
	       argv[0],argv[1]);
  return 0;
}

