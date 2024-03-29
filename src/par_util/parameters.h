/************************************************************
 *                                                          *
 *  parameters.h                                            *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 2000 Department of Statistics,         *
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
 *  Original programming by Joel Welling, 6/00              *
 ************************************************************/
/* This package implements methods for dealing with parameter files */

typedef enum { PRM_AFNI3D, PRM_ESTIREG3D, PRM_AIR, PRM_NONE } prm_type;

typedef struct prm_struct {
  prm_type type;
  int dim_x;
  int dim_y;
  int dim_z;
  float vox_x;
  float vox_y;
  float vox_z;
  int dz;
  int dt;
  int npar;
  int datasize;
  double* data;
} PrmStruct;

int prm_needs_dims( prm_type type );
int prm_needs_vox( prm_type type );
int prm_can_load( prm_type type );
PrmStruct* prm_load( prm_type type, int dz, int dt, 
		     int dim_x, int dim_y, int dim_z, 
		     float vox_x, float vox_y, float vox_z,
		     char* fname );
int prm_can_save( prm_type type );
void prm_save( PrmStruct* ps, char* fname );
PrmStruct* prm_translate( PrmStruct* in, prm_type type );
void prm_free( PrmStruct* s );
prm_type prm_name_to_type( char* name );
char* prm_type_to_name( prm_type type );
double prm_access( PrmStruct* in, int z, int t, int p );
void prm_add( PrmStruct* ps, int z, int t, int np, double* data_in );
PrmStruct* prm_copy( PrmStruct* in );
