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

#ifndef INCL_PULSE_H
#define INCL_PULSE_H

typedef enum 
{ EPI_READOUT_SQUARE, EPI_READOUT_TRAPEZOIDAL } EpiReadoutType;

extern const char* pul_getEpiReadoutTypeName( const EpiReadoutType type );

/* Routine to calculate relative k-space locations of ramped EPI lines */
extern int pul_calcEpiReadoutCoords( const int npts, double* coordOut,
				     const double rampUp, 
				     const double flatTopTime, 
				     const double rampDown,
				     const double sampDelayTime,
				     const double durationADC, 
				     const double lag,
				     const EpiReadoutType gradType );

#endif

