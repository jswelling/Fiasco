/************************************************************
 *                                                          *
 *  estireg_utils.h                                             *
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


/* "$Id: estireg_utils.h,v 1.1 2004/08/05 23:16:02 welling Exp $" */


int checkDatasetDims( char* dsname, MRI_Dataset* ds, char* chunk, 
		      char* dims_required, const char* progname );

void loadImage( float* image, MRI_Dataset* ds, int t, int dx, int dy, int dz );

void loadImageComplex( FComplex* image, MRI_Dataset* ds, int t,
		       int dx, int dy, int dz );

void buildTimeString( char* time_string, long sLength,
			     struct rusage* start, struct rusage* end );

void invertStdvImage( float *image, int dx, int dy, int dz );

void copyImage( float *image, float* source, int dx, int dy, int dz );



