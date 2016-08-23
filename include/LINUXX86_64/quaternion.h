/************************************************************
 *                                                          *
 *  quaternion.h                                            *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *                                                          *
 *  Copyright (c) 1999 Department of Statistics             *
 *                     Carnegie Mellon University           *
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
 *  Original programming by Joel Welling 2/99               *
 ************************************************************/
/* Header file for quaternions and associated transforms */

typedef struct quaternion_struct {
  double x; 
  double y; 
  double z; 
  double w;
} Quat; /* a quaternion */

typedef double Vec4[4]; /* a homogeneous vector */

typedef double Transform[16]; /* a homogeneous transform */

void trans_identity( Transform t_out );

void trans_copy( Transform t_out, const Transform t_in );

void trans_mult_right( Transform t, const Transform factor );

void trans_mult_left( const Transform factor, Transform t );

void trans_transpose( Transform t );

int trans_inverse( Transform invT, const Transform T );

void trans_vec_mult( const Transform factor, Vec4 v );  

void trans_dump( FILE* ofile, Transform t );

Quat* trans_to_quat( Quat* q_out, Transform t ); /* translations ignored */

void quat_to_trans( Transform t_out, 
		    Quat* q, double dx, double dy, double dz );

Quat* quat_copy( Quat* q_out, Quat* q_in );

/* input quat assumed normalized! */
Quat* quat_nrm_sqrt( Quat* result, const Quat* q_in ); 

Quat* quat_normalize( Quat* q );

Quat* quat_mult_right( Quat* q, Quat* factor );

Quat* quat_mult_left( Quat* factor, Quat* q );

Quat* quat_identity( Quat* q );

Quat* quat_conjugate( Quat* q );

Quat* quat_create( Quat* q, double x, double y, double z, int w_pos );

Quat* quat_from_axis_angle( Quat* q, double x, double y, double z, 
			    double theta );

void quat_to_axis_angle( Quat* q, double* x, double* y, double* z, 
			 double* theta );

Quat* quat_from_euler_RzRyRx( Quat* q, double x_angle, double y_angle,
			      double z_angle );

int quat_to_euler_RzRyRx( Quat* q, double* x_angle, double* y_angle,
			   double* z_angle ); /* returns 1 on success */
