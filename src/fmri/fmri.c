/************************************************************
 *                                                          *
 *  fmri.c                                                  *
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
 *  Original programming by Mark Fitzgerald  4-96           *
 ************************************************************/
/*************************************************************

  DESCRIPTION OF FMRI.C

  A set of utility functions to support other FIASCO code.

  File contains the following functions:
    get_missing
    get_typesize
    realign_matrix
    efseek
    Modulus
    Phase
    fqsrt

**************************************************************/


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include "bio.h"
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"

static char rcsid[] = "$Id: fmri.c,v 1.14 2004/09/10 17:30:46 welling Exp $";

static void check_validity( MRI_Dataset* DS, long* dzp, long* dtp, 
			    const char* chunk ) 
{
  char buf[256];
  char* dimstr;

  snprintf(buf,sizeof(buf),"%s.dimensions",chunk);
  if (!mri_has(DS,buf))
    Abort("Input dataset has %s info but no %s.dimensions!\n",
	  chunk, buf);
  dimstr= mri_get_string(DS,buf);
  /* There had better be at least one of the z and t dimensions, or
   * this whole notion makes no sense.  If only one is present we 
   * assume an extent of 1 for the other.
   */
  snprintf(buf,sizeof(buf),"%s.extent.z",chunk);
  if (mri_has(DS,buf)) {
    *dzp = mri_get_int( DS, buf );
    snprintf(buf,sizeof(buf),"%s.extent.t",chunk);
    if (mri_has(DS,buf)) {
      *dtp = mri_get_int( DS, buf );
    }
    else *dtp= 1;
  }
  else {
    *dzp= 1;
    snprintf(buf,sizeof(buf),"%s.extent.t",chunk);
    if (mri_has(DS,buf)) {
      *dtp = mri_get_int( DS, buf );
    }
    else {
      Abort("Input dataset images chunk has neither extent.z nor extent.t!\n");
      *dtp= 1; /* not reached */
    }
  }
}

/* Get or create a vector of indicators for missing images */
unsigned char **get_missing( MRI_Dataset* DS )
{
  unsigned char **missing = NULL;
  long t, z, dz, dt, missing_dt, missing_dz;
  int images_present= 0;
  int samples_present= 0;
  int missing_present= 0;
  char* missing_dims= NULL;

  images_present= mri_has( DS,"images" );
  samples_present= mri_has( DS,"samples" );
  missing_present= mri_has( DS,"missing" );

  /* Do we know enough to continue? */
  if (!images_present && !samples_present && !missing_present)
    Abort( "Can't create missing chunk w/o # of slices and images." );

  /* If an images chunk is present, check validity */
  if (images_present) {
    check_validity(DS, &dz, &dt, "images");
  }
  else if (samples_present) {
    check_validity(DS, &dz, &dt, "samples");
  }

  /* If a missing chunk is present, check validity */
  if (missing_present) {
    if (!mri_has(DS,"missing.dimensions"))
      Abort("Input dataset has missing info but no missing.dimensions!\n");
    missing_dims= mri_get_string(DS,"missing.dimensions");
    if (strcmp(missing_dims,"tz") && strcmp(missing_dims,"zt"))
      Abort("Input dataset missing chunk has inappropriate dimensions!\n");
    if (!mri_has(DS,"missing.extent.z"))
      Abort("Input dataset missing chunk has no extent.z!\n");
    if (!mri_has(DS,"missing.extent.t"))
      Abort("Input dataset missing chunk has no extent.t!\n");
    if (!mri_has(DS,"missing.datatype"))
      Abort("Input dataset missing chunk has no datatype!\n");
    if (strcmp( mri_get_string( DS, "missing.datatype" ), "uint8" ))
      Abort("Input dataset missing chunk has wrong datatype!\n");
    missing_dz = mri_get_int( DS, "missing.extent.z" );
    missing_dt = mri_get_int( DS, "missing.extent.t" );
  }

  /* If both are present, are they consistent? */
  if (missing_present) {
    if (images_present) {
      if (dz != missing_dz)
	Abort("Input dataset images.extent.z and missing.extent.z differ!\n");
      if (dt != missing_dt)
	Abort("Input dataset images.extent.z and missing.extent.z differ!\n");
    }
    else {
      dz= missing_dz;
      dt= missing_dt;
    }
  }

  /* Create "missing" vector, initialized to 0 */
  missing = Matrix( dt, dz, unsigned char );

  /* Read in "missing" chunk if it already exists,
   * or create it if it does not.
   */
  if (missing_present) {
    unsigned char* mbuf;
    mbuf= (unsigned char *)mri_get_chunk( DS, "missing", (int) ( dt * dz ), 
					  0, MRI_UNSIGNED_CHAR );
    if (!strcmp(missing_dims,"zt")) {
      memcpy( *missing, mbuf, (long) ( dt * dz ) );
    }
    else {
      /* t varies fastest */
      for (z=0; z<dz; z++)
	for (t=0; t<dt; t++) missing[t][z]= mbuf[z*dt + t];
    }
  }
  else {
    /* Set all as non-missing */
    for( t = 0; t < dt; t++ )
      for( z = 0; z < dz; z++ )
	missing[t][z] = (unsigned char) 0;

    /* Create the chunk if possible */
    if (DS->mode != MRI_READ) {
      /* Set appropriate parameters for chunk and write out */
      mri_create_chunk( DS, "missing" );
      mri_set_string( DS, "missing.datatype", "uint8" );
      mri_set_string( DS, "missing.dimensions", "zt" );
      mri_set_int( DS, "missing.extent.z", 
		   mri_get_int( DS, "images.extent.z" ) );
      mri_set_string( DS, "missing.description.z", "gridded image-space");
      mri_set_int( DS, "missing.extent.t", 
		   mri_get_int( DS, "images.extent.t" ) );
      mri_set_string( DS, "missing.description.t", "gridded image-space");
      mri_set_chunk( DS, "missing", (int) (dz * dt), 0,
		     MRI_UNSIGNED_CHAR, missing[0] );
    }
  }

  return( missing );

}


/* Finds the size of an individual value for a   */
/*   chunk of a dataset                          */
long get_typesize( MRI_Dataset* DS, char* chunkname )
{
  char *datatype_key;
  char *datatype;
  long typesize;

  datatype_key = mri_cat(chunkname, ".datatype");
  if( !mri_has( DS, datatype_key ) )
    Abort( "Can't find type size --- %s key is missing.", datatype_key );
  datatype = mri_get_string( DS, datatype_key );
  if( !strcmp( datatype, "uint8" ) )
    typesize = 1;
  else if( !strcmp( datatype, "int16" ) )
    typesize = 2;
  else if( !strcmp( datatype, "int32" ) )
    typesize = 4;
  else if( !strcmp( datatype, "float32" ) )
    typesize = 4;
  else if( !strcmp( datatype, "float64" ) )
    typesize = 8;
  else
    Abort( "get_typesize: Unrecognized data type --- %s", datatype );
  return( typesize );
}


/* Realign two-dimensional array, after pointer to */
/*   initial array point has been changed          */
/*      (as when mri_get_image is used)            */
void realign_matrix( void** ptr, long n, size_t rowsize )
{
  long i;
  char **A = (char **) ptr;
  for( i = 1; i < n; i++ )
    A[i] = A[0] + i * rowsize;
  return;
}


/* Error-checking fseek function */
void efseek( FILE* fp, long offset, int whence, char* filename )
{
  if( fseek( fp, offset, whence ) )
    Abort( "Error positioning in file %s at byte %ld.",
	   filename, SEEK_CUR );
  return;
}


/* Calculates the modulus of a complex number */
float Modulus( FComplex z )
{
  return( (float) sqrt( ( z.real*z.real ) + ( z.imag*z.imag ) ) );
}

/* Calculates the arg of a complex number */
float Phase( FComplex z )
{
  if( z.real || z.imag )
    return( (float) atan2( z.imag, z.real ) );
  else
    return( (float) 0.0 );
}

static int fqsrt_compare( const void* v1, const void* v2 )
{
  float* f1= (float*)v1;
  float* f2= (float*)v2;
  if (*f1<*f2) return -1;
  else if (*f1==*f2) return 0;
  else return 1;
}

/* A version of quicksort for floating point arrays */
void fqsrt( float* v, long left, long right )
{
  qsort(v, (right-left)+1, sizeof(float), fqsrt_compare);
}


