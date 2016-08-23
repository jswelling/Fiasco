/*
 *	Array utility functions which allow the allocation/deallocation
 *		of 2 and 3-dimensional arrays where the elements are
 *		stored contiguously in order to ensure predictable cache
 *		behavior and to have the potential to use optimized
 *		1-dimensional routines for processing and I/O.
 *
 *	Copyright (c) 1996 Pittsburgh Supercomputing Center
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
 *		1/96 Written by Greg Hood (PSC)
 */

extern float **Alloc2DFloatArray (int d1, int d2);
extern void Free2DFloatArray (float **a);
extern float ***Alloc3DFloatArray (int d1, int d2, int d3);
extern void Free3DFloatArray (float ***a);
extern short int **Alloc2DShortArray (int d1, int d2);
extern void Free2DShortArray (short int **a);
extern double **Alloc2DDoubleArray (int d1, int d2);
extern void Free2DDoubleArray (double **a);
extern double ***Alloc3DDoubleArray (int d1, int d2, int d3);
extern void Free3DDoubleArray (double ***a);
