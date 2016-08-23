/************************************************************
 *                                                          *
 *  mriu.h                                                  *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 2007 Department of Statistics,         *
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
 *  Original programming by Joel Welling, 7/07              *
 ************************************************************/
/* This module contains supplemental utilities for manipulating 
 * Pgh MRI files.
 */

/* Function prototype to communicate index remapping pattern to
 * mriu_updateLabels().  When called with a given old index, the
 * corresponding new index is returned.  A return value < 0 indicates
 * that the given index no longer exists.
 */
typedef int (*DIMREMAPFUNC)(int oldIndex, void* hook);

/* Utility for modifying label-type tags, of the form 'chunk.label.v.n'
 * where v is a dimension and n is an index.
 */
void mriu_updateLabels(MRI_Dataset* ds, const char* chunk, 
		       char dim, int low, int highplusone, 
		       DIMREMAPFUNC remapFunc, void* hook);
