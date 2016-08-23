/************************************************************
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

/*  $Id: read.c,v 1.3 1999/07/07 23:49:01 welling Exp $ */
/*
readField ()
readHead ("final1")
makec("read")
iterateFieldLines("final1")
*/

#include <stdio.h>
#include "../src/mr.h"

static char rcsid[] = "$Id: read.c,v 1.3 1999/07/07 23:49:01 welling Exp $";

static field TheField;
static first = 1;
static char Name [1000];

readFieldHead (char ** name,
		 int * type,
		 int * ndim,
		 int * veclen,
		 int * dims)
{
	int i;
	field * answer = &TheField;
	if (first){
		first = 0;
		TheField = ReadMrMmap (*name);
		strcpy (Name, *name);
	}
	if (strcmp (Name, *name)){
		unmapField (TheField);
		TheField = ReadMrMmap (*name);
		strcpy (Name, *name);
	}
	*type = answer ->type;
	*ndim = answer ->ndim;
	*veclen = answer ->veclen;
	for (i = 0; i < *ndim; i++)
		dims [i] = answer ->dimensions [i];
}

readField (int * dimsL, int * dimsU,
			  float * value)
{
	int i0, i1, i2, i3;
	int d0 = TheField .dimensions [0];
	int d1 = TheField .dimensions [1];
	int d2 = TheField .dimensions [2];
	int d3 = TheField .dimensions [3];
	int dd0 = dimsU [0] - dimsL [0]+1;
	int dd1 = dimsU [1] - dimsL [1]+1;
	int dd2 = dimsU [2] - dimsL [2]+1;
	int dd3 = dimsU [3] - dimsL [3]+1;
	if (first){
		fprintf (stderr, "No field is open\n");
		return -1;
	}
	for (i3 = dimsL [3]; i3 <= dimsU [3]; i3++)
		for (i2 = dimsL [2]; i2 <= dimsU [2]; i2++)
			for (i1 = dimsL [1]; i1 <= dimsU [1]; i1++)
				for (i0 = dimsL [0]; i0 <= dimsU [0]; i0++){
					int ii0 = i0 - dimsL [0];
					int ii1 = i1 - dimsL [1];
					int ii2 = i2 - dimsL [2];
					int ii3 = i3 - dimsL [3];
					int ind1 =
						ii0 + ii1*dd0 + ii2*dd1*dd0 + ii3*dd2*dd1*dd0;
					int ind2 = i0 + i1*d0 + i2*d1*d0 + i3*d2*d1*d0;
					switch (TheField .type){
					case CHAR:
						value [ind1] =
							(float)(((unsigned char*)TheField .data) [ind2]);
						break;
					case INT:
						value [ind1] =
							(float)(((int*)TheField .data) [ind2]);
						break;
					case FLOAT:
						value [ind1] =
							(float)(((float*)TheField .data) [ind2]);
						break;
					case DOUBLE:
						value [ind1] =
							(float)(((double*)TheField .data) [ind2]);
						break;
					case SHORT:
						value [ind1] =
							(float)(((short*)TheField .data) [ind2]);
						break;
					}
				}
}

writeField (char ** name,
		 int * type,
		 int * ndim,
		 int * veclen,
		 int * dims, float * val)
{
	field ans;
	int i;
	ans .type = *type;
	ans .ndim = *ndim;
	ans .veclen = *veclen;
	ans .dimensions  = dims;
	ans .data = (unsigned char *) val;
	WriteMr (ans, *name);
}

/*
readFieldHead (char ** name,
		 int * type,
		 int * ndim,
		 int * veclen,
		 int * dims)
{
	int i;
	field answer = ReadMrMmap (*name);
	*type = answer .type;
	*ndim = answer .ndim;
	*veclen = answer .veclen;
	for (i = 0; i < *ndim; i++)
		dims [i] = answer .dimensions [i];
	unmapField (answer);
}

readField (char ** name,
			  int * dimsL, int * dimsU,
			  float * value)
{
	int i0, i1, i2, i3;
	field answer = ReadMrMmap (*name);
	int d0 = answer .dimensions [0];
	int d1 = answer .dimensions [1];
	int d2 = answer .dimensions [2];
	int d3 = answer .dimensions [3];
	int dd0 = dimsU [0] - dimsL [0]+1;
	int dd1 = dimsU [1] - dimsL [1]+1;
	int dd2 = dimsU [2] - dimsL [2]+1;
	int dd3 = dimsU [3] - dimsL [3]+1;
	for (i3 = dimsL [3]; i3 <= dimsU [3]; i3++)
		for (i2 = dimsL [2]; i2 <= dimsU [2]; i2++)
			for (i1 = dimsL [1]; i1 <= dimsU [1]; i1++)
				for (i0 = dimsL [0]; i0 <= dimsU [0]; i0++){
					int ii0 = i0 - dimsL [0];
					int ii1 = i1 - dimsL [1];
					int ii2 = i2 - dimsL [2];
					int ii3 = i3 - dimsL [3];
					int ind1 =
						ii0 + ii1*dd0 + ii2*dd1*dd0 + ii3*dd2*dd1*dd0;
					int ind2 = i0 + i1*d0 + i2*d1*d0 + i3*d2*d1*d0;
					value [ind1]
						= ((float*)answer .data) [ind2];
				}
	unmapField (answer);
}
*/






