/************************************************************
 *                                                          *
 *  mr.h                                                    *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1995 Department of Statistics,         *
 *                        Carnegie Mellon University        *
# *                                                          *
# *  This program is distributed in the hope that it will    *
# *  be useful, but WITHOUT ANY WARRANTY; without even the   *
# *  implied warranty of MERCHANTABILITY or FITNESS FOR A    *
# *  PARTICULAR PURPOSE.  Neither Carnegie Mellon University *
# *  nor any of the authors assume any liability for         *
# *  damages, incidental or otherwise, caused by the         *
# *  installation or use of this software.                   *
# *                                                          *
# *  CLINICAL APPLICATIONS ARE NOT RECOMMENDED, AND THIS     *
# *  SOFTWARE HAS NOT BEEN EVALUATED BY THE UNITED STATES    *
# *  FDA FOR ANY CLINICAL USE.                               *
# *                                                          *
 *                                                          *
 *  Original programming by Audris Mockus                   *
 ************************************************************/

#ifndef MR_H
#define MR_H

/* $Id: mr.h,v 1.2 1999/07/07 23:49:11 welling Exp $ */

#include <stdio.h>

typedef struct {
	int type, ndim, veclen, * dimensions;
	unsigned char * data;
} field;

typedef struct {
	int type, ndim, veclen, * dimensions;
	int dataOffset;
	FILE * f;
} fField;

enum {
	CHAR = 0, INT = 1, FLOAT = 2, DOUBLE = 3, SHORT = 4};
#ifdef __cplusplus
extern "C" {
#endif
#define MISSING -1

	FILE *  myfopen (char * fName, char * mode);
	int typeSize (int type);
	char * typeName (int type);
	field allocField (int type, int dim, int vec, int * dims);
	fField fileField (int type, int dim, int vec, int * dims, char * name);
	field mmapField (int type, int dim, int vec, int * dims, FILE * f);
	void  freeField (field XX);
	void  unmapField (field XX);
	field ReadMr (char *filename);
	field ReadMrMmap (char *filename);
	int WriteMr (field data, char *filename);
	int PrintMr (field data, char *filename);
	void Split (int SKIP, char * inFName,
	    char * outFName);
	int STAT (char * inFName, FILE * summfile);
	void SPM (int cond1, int cond2, int l1, int l2,
	    char * inFName1);

	int ImRegister (field * Input, field * RegisterTo,
	    field * iPars, field * oPars,
	    int Loss,
	    int Method, float Power, int Forward);

#ifdef __cplusplus
}
#endif
#endif
