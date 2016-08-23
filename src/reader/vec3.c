/************************************************************
 *                                                          *
 *  vec3.c                                             *
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
 *  Original programming by Mark Fitzgerald  5-96           *
 *  Modified to exclusively use libmri calls for output,    *
 *       Greg Hood (PSC), 9-98                              *
 *  Modified to read header files, and to use LX2           *
 *       resampling stuff, Joel Welling (PSC/Stats), 5-1999 *
 *  Major restructuring, and the addition of formidable     *
 *       flexibility and intelligence, Joel Welling 5-2002  *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mri.h"
#if (SGI64 || SGI5 || SGIMP)
#include <bstring.h>
#endif

#include "bio.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"
#include "smartreader.h"

static char rcsid[] = "$Id: vec3.c,v 1.4 2005/02/05 02:57:45 welling Exp $";

static char* anatDimNames[]= {"Sag","Cor","Tra"};

int getVec3(KVHash* info, char* name, double* result)
{
  char buf[64];
  int i; 

  if (strlen(name)>60) 
    Abort("%s: multi_reader: internal error; buffer too short!\n",
	  progname);

  for (i=0; i<3; i++) {
    sprintf(buf,"%s.%d",name,i);
    if (!kvLookup(info,buf)) return 0;
    result[i]= kvGetDouble(info,buf);
  }
  return 1;
}

int testVec3(KVHash* info, char* name)
{
  double v[3];
  return(getVec3(info,name,v) != 0);
}

int getVec3_anat(KVHash* info, char* name, double* result)
{
  char buf[128];
  int i; 
  int found= 0;

  if (strlen(name)>120) 
    Abort("%s: multi_reader: internal error; buffer too short!\n",
	  progname);

  for (i=0; i<3; i++) {
    sprintf(buf,"%s%s",name,anatDimNames[i]);
    if (kvLookup(info,buf)) {
      result[i]= kvGetDouble(info,buf);
      found= 1;
    }
    else result[i]= 0.0;
  }
  return (found!=0);
}

void defVec3(KVHash* info, char* name, double* vals)
{
  char buf[64];
  int i;

  if (strlen(name)>60) 
    Abort("%s: multi_reader: internal error; buffer too short!\n",
	  progname);

  for (i=0; i<3; i++) {
    sprintf(buf,"%s.%d",name,i);
    kvDefDouble(info,buf,vals[i]);
  }
}

void subtractVec3( double* result, const double* v1, const double* v2 )
{
  result[0]= v1[0]-v2[0];
  result[1]= v1[1]-v2[1];
  result[2]= v1[2]-v2[2];
}

void multVec3( double* result, const double* v, double factor )
{
  result[0]= v[0]*factor;
  result[1]= v[1]*factor;
  result[2]= v[2]*factor;
}

void copyVec3( double* to, const double* from )
{
  to[0]= from[0];
  to[1]= from[1];
  to[2]= from[2];
}

double dotVec3( const double* v1, const double* v2 )
{
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

double normVec3( const double* v )
{
  double sum= v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
  return sqrt(sum);
}

void normalizeVec3( double* v )
{
  double cpy[3];
  double norm;
  copyVec3( cpy, v );
  norm= normVec3( cpy );
  if (norm==0.0) return;
  else multVec3( v, cpy, 1.0/norm );
}

void crossVec3( double* result, const double* v1, const double* v2 )
{
  result[0]= v1[1]*v2[2]-v1[2]*v2[1];
  result[1]= v1[2]*v2[0]-v1[0]*v2[2];
  result[2]= v1[0]*v2[1]-v1[1]*v2[0];
}

void xplusbyVec3( double* result, const double* v1, const double *v2, 
		  double scale )
{
  result[0]= v1[0]+scale*v2[0];
  result[1]= v1[1]+scale*v2[1];
  result[2]= v1[2]+scale*v2[2];
}

void flipToPositiveHemisphereVec3( double* vec )
{
  /* This routine takes a vector and reflects it to the positive
   * direction if it is in the negative direction.
   */
  static double posVec[]= {1.0,1.0,1.0};
  if (dotVec3(vec,posVec)<0.0) {
    double tmpVec[3];
    multVec3(tmpVec,vec,-1.0);
    copyVec3(vec,tmpVec);
  }
}

