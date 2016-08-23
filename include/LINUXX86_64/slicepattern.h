/************************************************************
 *                                                          *
 *  slicepattern.h                                           *
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
 ************************************************************/

/* slp_generateSlicePatternTable() returns a table of slice reordered
 * position by temporal acquisition number, so for example if table[0]==9
 * time, the 0th slice acquired (counting from 0) goes in the spatial
 * position 9.  The Inverted version does the opposite, such that
 * invertedTable[0]==9 means the 9th acquired slice goes in spatial
 * position 0.  The caller owns the memory allocated by these routines
 * for the returned tables.
 */
int* slp_generateSlicePatternTable( int dz, const char* name );
int* slp_invertSlicePatternTable(int dz, const int* directTable);
int* slp_generateInvertedSlicePatternTable( int dz, const char* name );
const char* slp_findSlicePatternNameFromTable(int dz, const int* table);

