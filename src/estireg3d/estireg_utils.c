/************************************************************
 *                                                          *
 *  estireg_utils.c                                             *
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
 *  Original programming by Joel Welling 3/00               *
 *      10/02: Parallelization, Jenn Bakal                  *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifdef DARWIN
#include <sys/time.h>
#endif
#include <time.h>
#include <math.h>
#include <limits.h>
#include <sys/resource.h>
#include <unistd.h>
#include "../fmri/lapack.h"
#include "mri.h"
#include "fmri.h"
#include "par.h"
#include "stdcrg.h"
#include "misc.h"
#include "array.h"
#include "acct.h"
#include "algorithm.h"
#include "estireg_utils.h"


static char rcsid[] = "$Id: estireg_utils.c,v 1.2 2004/08/11 17:16:31 welling Exp $";

/* Smallest standard deviation contributing to weighting.  This
 * prevents numerical noise weights from overwhelming the signal.
 */
#define STDV_FLOOR 1.0


int checkDatasetDims( char* dsname, MRI_Dataset* ds, char* chunk, 
		      char* dims_required, const char* progname )
{
  char* dimstr;
  char buf[256];
  int i;
  int offset;

  if (strlen(chunk)>200) 
    Abort("%s: chunk name <%s> too long!\n",chunk);

  if( !mri_has( ds, chunk ) ) {
    Message( "%s: dataset <%s> has no \"%s\" chunk!", progname, dsname, chunk);
    return 0;
  }

  sprintf(buf,"%s.dimensions",chunk);
  if ( !mri_has( ds, buf ) ) {
    Message( "%s: dataset <%s> has no dimension string for chunk \"%s\"!\n",
	     progname, dsname, chunk );
    return 0;
  }
  else dimstr= mri_get_string( ds, buf );

  offset= 0;
  for (i=0; i<strlen(dimstr); i++) {
    int thisextent;
    char* tchar;
    sprintf(buf,"%s.extent.%c",chunk,dimstr[i]);
    if (!mri_has(ds,buf)) {
      Message("%s: dataset <%s> has no info for %s!\n",buf);
      return 0;
    }
    thisextent= mri_get_int(ds,buf);
    if ((tchar=strchr(dims_required,dimstr[i])) != NULL) {
      /* Make sure this is farther along dims_required than the last one */
      if (tchar-dims_required >= offset) offset= tchar-dims_required;
      else {
	/* Order of strings doesn't match */
	Message( "%s: dataset <%s> dimension order doesn't match \"%s\"!\n",
		 progname, dsname, dims_required );
	return 0;
      }
    }
    else {
      /* not in required string- extent must be 1 */
      if (thisextent != 1) {
	Message( "%s: dataset <%s> has an inappropriate dimension %c!\n",
		 progname, dsname, dimstr[i] );
	return 0;
      }
    }
  }

  /* Were all of the required dims present? */
  if (offset!=strlen(dims_required)-1) {
      Message( "%s: dataset <%s> doesn't include dimensions \"%s\"!\n",
	       progname, dsname, dims_required );
      return 0;    
  }

  return 1;
}

void loadImage( float* image, MRI_Dataset* ds, int t, int dx, int dy, int dz ) 
{
  int x;
  int y;
  int z;
  int blocksize;
  float* tmp_image;

  /* This routine reads the given image from the file, swaps it to 
   * z-fastest order, and returns it in the buffer supplied.  Dims
   * are assumed to be "xyzt" for the input dataset.
   */
  blocksize= (int)(dx*dy*dz);
  tmp_image= (float*)mri_get_chunk(ds, "images", blocksize,
				   t*blocksize, MRI_FLOAT);
  
  /* Convert to complex-valued for registration */
  for (z=0; z<dz; z++) 
    for(y=0; y<dy; y++)
      for(x=0; x<dx; x++)
	MEM(image,dx,dy,dz,x,y,z)= MEM_BCK(tmp_image,dx,dy,dz,x,y,z);
  
}

void loadImageComplex( FComplex* image, MRI_Dataset* ds, int t,
		       int dx, int dy, int dz ) 
{
  int x;
  int y;
  int z;
  int blocksize;
  float* tmp_image;

  /* This routine reads the given image from the file, swaps it to 
   * z-fastest order, and returns it in the buffer supplied.  Dims
   * are assumed to be "xyzt" for the input dataset.
   */
  blocksize= (int)(dx*dy*dz);
  tmp_image= (float*)mri_get_chunk(ds, "images", blocksize,
				   t*blocksize, MRI_FLOAT);
  
  /* Convert to complex-valued for registration */
  for (z=0; z<dz; z++) 
    for(y=0; y<dy; y++)
      for(x=0; x<dx; x++) {
	MEM(image,dx,dy,dz,x,y,z).real= MEM_BCK(tmp_image,dx,dy,dz,x,y,z);
	MEM(image,dx,dy,dz,x,y,z).imag= 0.0;
      }
}


void buildTimeString( char* time_string, long sLength,
		      struct rusage* start, struct rusage* end )
{
  long s_usec= end->ru_stime.tv_usec - start->ru_stime.tv_usec;
  long s_sec= end->ru_stime.tv_sec - start->ru_stime.tv_sec;
  long u_usec= end->ru_utime.tv_usec - start->ru_utime.tv_usec;
  long u_sec= end->ru_utime.tv_sec - start->ru_utime.tv_sec;

  if (s_usec < 0) {
    s_sec -= 1;
    s_usec += 1000000;
  }
  if (u_usec < 0) {
    u_sec -= 1;
    u_usec += 1000000;
  }
  snprintf(time_string,sLength,"%d.%06du %d.%06ds",u_sec,u_usec,s_sec,s_usec);
  time_string[sLength-1]= '\0';
}

/* Take the inverse of real data in place */
void invertStdvImage( float *image, int dx, int dy, int dz )
{
  long length= dx*dy*dz;
  long i;

  for (i=0; i<length; i++) {
    if (fabs(image[i]) >= STDV_FLOOR) image[i]= 1.0/image[i];
    else image[i]= 0.0;
  }
}

/* Convenience routine */
void copyImage( float *image, float* source, int dx, int dy, int dz )
{
  long length= dx*dy*dz;
  long i;
  for (i=0; i<length; i++) image[i]= source[i];
}

