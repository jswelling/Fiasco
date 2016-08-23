/*
 *	Array allocator module - dynamically allocates 2 and 3 dimensional
 *		arrays so that the elements are contiguous in memory
 *
 *	Copyright (c) 1995  Pittsburgh Supercomputing Center
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
 *
 *	HISTORY
 *		12/95 Written by Greg Hood
 *		1/96  Split off from util.c (ghood)
 */

#include <stdio.h>
#include <stdlib.h>

static char rcsid[] = "$Id: libarray.c,v 1.4 2002/04/04 18:54:42 welling Exp $";

/* Alloc2DFloatArray dynamically allocates a 2-D array
       of floats with dimensions [d1][d2].  The actual
       contents of the array are stored contiguously so
       that block operations may be performed on the
       memory chunk holding the array. */
float **
Alloc2DFloatArray (int d1,
	      int d2)
{
  int i;
  float **a;
  float *contents;

  contents = (float *) malloc(d1*d2*sizeof(float));
  a = (float **) malloc(d1 * sizeof(float *));
  for (i = 0; i < d1; ++i)
    a[i] = contents + i*d2;
  return(a);
}

/* Free2DFloatArray deallocates a 2-D array created
   by Alloc2DFloatArray */
void Free2DFloatArray (float **a)
{
  if (a != NULL)
    {
      free(a[0]);
      free(a);
    }
}

/* Alloc3DFloatArray dynamically allocates a 3-D array
       of floats with dimensions [d1][d2][d3].  The actual
       contents of the array are stored contiguously so
       that block operations may be performed on the
       memory chunk holding the array. */
float ***
Alloc3DFloatArray (int d1,
		   int d2,
		   int d3)
{
  int i, j;
  float ***a;
  float **index;
  float *contents, *offset;

  contents = (float *) malloc(d1*d2*d3*sizeof(float));
  index = (float **) malloc(d1*d2*sizeof(float *));
  a = (float ***) malloc(d1*sizeof(float **));
  for (i = 0; i < d1; ++i)
    {
      a[i] = index + i*d2;
      offset = contents + i*d2*d3;
      for (j = 0; j < d2; ++j)
	a[i][j] = offset + j*d3;
    }
  return(a);
}

/* Free3DFloatArray deallocates a 3-D array created
   by Alloc3DFloatArray */
void Free3DFloatArray (float ***a)
{
  if (a != NULL)
    {
      free(a[0][0]);
      free(a[0]);
      free(a);
    }
}


/* Alloc2DShortArray dynamically allocates a 2-D array
       of short ints with dimensions [d1][d2].  The actual
       contents of the array are stored contiguously so
       that block operations may be performed on the
       memory chunk holding the array. */
short int **
Alloc2DShortArray (int d1,
		   int d2)
{
  int i;
  short int **a;
  short int *contents;

  contents = (short int *) malloc(d1*d2*sizeof(short int));
  a = (short int **) malloc(d1*sizeof(short int *));
  for (i = 0; i < d1; ++i)
    a[i] = contents + i*d2;
  return(a);
}

/* Free2DShortArray deallocates a 2-D array created
   by Alloc2DShortArray */
void Free2DShortArray (short int **a)
{
  if (a != NULL)
    {
      free(a[0]);
      free(a);
    }
}

/* Alloc2DDoubleArray dynamically allocates a 2-D array
       of doubles with dimensions [d1][d2].  The actual
       contents of the array are stored contiguously so
       that block operations may be performed on the
       memory chunk holding the array. */
double **
Alloc2DDoubleArray (int d1,
	      int d2)
{
  int i;
  double **a;
  double *contents;

  contents = (double *) malloc(d1*d2*sizeof(double));
  a = (double **) malloc(d1 * sizeof(double *));
  for (i = 0; i < d1; ++i)
    a[i] = contents + i*d2;
  return(a);
}

/* Free2DDoubleArray deallocates a 2-D array created
   by Alloc2DDoubleArray */
void Free2DDoubleArray (double **a)
{
  if (a != NULL)
    {
      free(a[0]);
      free(a);
    }
}

/* Alloc3DDoubleArray dynamically allocates a 3-D array
       of doubles with dimensions [d1][d2][d3].  The actual
       contents of the array are stored contiguously so
       that block operations may be performed on the
       memory chunk holding the array. */
double ***
Alloc3DDoubleArray (int d1,
		    int d2,
		    int d3)
{
  int i, j;
  double ***a;
  double **index;
  double *contents, *offset;

  contents = (double *) malloc(d1*d2*d3*sizeof(double));
  index = (double **) malloc(d1*d2*sizeof(double *));
  a = (double ***) malloc(d1*sizeof(double **));
  for (i = 0; i < d1; ++i)
    {
      a[i] = index + i*d2;
      offset = contents + i*d2*d3;
      for (j = 0; j < d2; ++j)
	a[i][j] = offset + j*d3;
    }
  return(a);
}

/* Free3DDoubleArray deallocates a 3-D array created
   by Alloc3DDoubleArray */
void Free3DDoubleArray (double ***a)
{
  if (a != NULL)
    {
      free(a[0][0]);
      free(a[0]);
      free(a);
    }
}


