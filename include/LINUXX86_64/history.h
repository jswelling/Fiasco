/************************************************************
 *                                                          *
 *  history.h                                               *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1997 Department of Statistics,         *
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
 *  Original programming by Joel Welling, 10/00             *
 ************************************************************/
/* This module contains routines for manipulating the history info
 * in Pgh MRI files.
 */

const char* hist_get( MRI_Dataset* ds, int i ); /* first entry is i=1! */
void hist_add( MRI_Dataset* ds, const char* s );
void hist_delete_all( MRI_Dataset* ds );
void hist_add_cl( MRI_Dataset* ds, int argc, char** argv );
void hist_dump( MRI_Dataset* ds, FILE* ofile );
void hist_delete_some( MRI_Dataset* ds, int n_to_delete );
void hist_repair( MRI_Dataset* ds );
