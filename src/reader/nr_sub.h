/************************************************************
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *                                                          *
 *     Copyright (c) 1999 Carnegie Mellon University        *
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
 ***********************************************************/
/* Headers for nr_sub.c */

float **nr_matrix(int nrl, int nrh, int ncl, int nch);
void nr_free_matrix(float**m, int nrl, int nrh, int ncl, int nch);

typedef struct FCOMPLEX {float r,i;} fcomplex;

fcomplex nr_Csub(fcomplex a, fcomplex b);
fcomplex nr_Cmul(fcomplex a, fcomplex b);
fcomplex nr_Complex(float re, float im);
fcomplex nr_Conjg(fcomplex z);
float nr_Cabs(fcomplex z);
fcomplex nr_Csqrt(fcomplex z);

