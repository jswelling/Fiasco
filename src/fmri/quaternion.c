/************************************************************
 *                                                          *
 *  quaternion.c                                            *
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
/* Methods for quaternions and associated transforms */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "misc.h"
#include "fmri.h"
#include "lapack.h"

static char rcsid[] = "$Id: quaternion.c,v 1.15 2007/07/06 22:10:26 welling Exp $";

/* This package provides manipulations of unit quaternions.
   Much code was taken from the "Matrix and Quaternions Faq".
   Unfortunately the trans_to_quat algorithm given there seems
   totally screwed up;  this one was rederived using Maple.
 */

/* How long to iterate looking for Euler angle decompositions of quats */
#define EULER_DECOMP_MAX_ITER 100

#define DEBUG 0

static Quat* quat_alloc()
{
  Quat* result;
  if (!(result= (Quat*)malloc(sizeof(Quat))))
      Abort("Unable to allocate %d bytes for quaternion!\n",sizeof(Quat));
  return result;
}

Quat* trans_to_quat( Quat* q_out, Transform t )
{
  double trace;

  if (!q_out) q_out= quat_alloc();

  trace= t[0]+t[5]+t[10]+1.0;
  if (trace>0.5) {
    double s;
    s= 0.5/sqrt(trace);
    q_out->w= 0.25/s;
    q_out->x = ( t[9] - t[6] ) * s;
    q_out->y = ( t[2] - t[8] ) * s;
    q_out->z = ( t[4] - t[1] ) * s;
  }
  else {
    double xd= 1.0 + t[0] - t[5] - t[10];
    double yd= 1.0 + t[5] - t[0] - t[10];
    double zd= 1.0 + t[10] - t[0] - t[5];
    if (xd>1.0) {
      double S= 2.0/sqrt( xd ); /* = 1/X */
      q_out->x = 1.0 / S;
      q_out->y = 0.25*(t[1] + t[4] )*S;
      q_out->z = 0.25*(t[2] + t[8] )*S;
      q_out->w = 0.25*(t[9] - t[6] )*S;
    }
    else if (yd>1.0) {
      double S = 2.0/sqrt( yd ); /* = 1/Y */
      q_out->x = 0.25*(t[1] + t[4])*S;
      q_out->y = 1.0 / S;
      q_out->z = 0.25*(t[6] + t[9])*S;
      q_out->w = 0.25*(t[2] - t[8])*S; 
    }
    else {
      double S = 2.0/sqrt( zd ); /* = 1/Z */
      q_out->x = 0.25*(t[2] + t[8] )*S;
      q_out->y = 0.25*(t[6] + t[9] )*S;
      q_out->z = 1.0 / S;
      q_out->w = 0.25*(t[4] - t[1] ) * S;
    }
  }

  /* quat_normalize(q_out); */

  return q_out;
}

void trans_identity( Transform t_out )
{
  int i;
  for (i=0; i<16; i++) t_out[i]= 0.0;
  for (i=0; i<4; i++) t_out[5*i]= 1.0;
}

void trans_copy( Transform t_out, const Transform t_in )
{
  int i;
  for (i=0; i<16; i++) t_out[i]= t_in[i];
}

void quat_to_trans( Transform t_out, 
		    Quat* q, double dx, double dy, double dz )
{
  double xx, xy, xz, xw, yy, yz, yw, zz, zw;

  xx      = q->x * q->x;
  xy      = q->x * q->y;
  xz      = q->x * q->z;
  xw      = q->x * q->w;
  
  yy      = q->y * q->y;
  yz      = q->y * q->z;
  yw      = q->y * q->w;
  
  zz      = q->z * q->z;
  zw      = q->z * q->w;
  
  t_out[0]  = 1 - 2 * ( yy + zz );
  t_out[1]  =     2 * ( xy - zw );
  t_out[2]  =     2 * ( xz + yw );
  
  t_out[4]  =     2 * ( xy + zw );
  t_out[5]  = 1 - 2 * ( xx + zz );
  t_out[6]  =     2 * ( yz - xw );
  
  t_out[8]  =     2 * ( xz - yw );
  t_out[9]  =     2 * ( yz + xw );
  t_out[10] = 1 - 2 * ( xx + yy );
  
  t_out[3]  = dx;
  t_out[7] = dy;
  t_out[11]= dz;
  t_out[12] = t_out[13] = t_out[14] = 0;
  t_out[15] = 1;
}

Quat* quat_copy( Quat* q_out, Quat* q_in )
{
  if (!q_out) q_out= quat_alloc();

  q_out->x= q_in->x;
  q_out->y= q_in->y;
  q_out->z= q_in->z;
  q_out->w= q_in->w;

  return q_out;
}

double quat_magnitude( Quat* q )
{
  return( sqrt(q->w*q->w+q->x*q->x+ q->y*q->y+q->z*q->z) );
}

Quat* quat_normalize( Quat* q )
{
  double mag= quat_magnitude(q);
  q->x /= mag;
  q->y /= mag;
  q->z /= mag;
  q->w /= mag;

  return q;
}

Quat* quat_mult_right( Quat* q, Quat* f )
{
  Quat t;
  t.w= q->w*f->w - (q->x*f->x + q->y*f->y + q->z*f->z);
  t.x= q->y*f->z - q->z*f->y + q->w*f->x + q->x*f->w;
  t.y= q->z*f->x - q->x*f->z + q->w*f->y + q->y*f->w;
  t.z= q->x*f->y - q->y*f->x + q->w*f->z + q->z*f->w;
  quat_copy(q,&t);
  return q;    
}

Quat* quat_mult_left( Quat* f, Quat* q )
{
  Quat t;
  t.w= f->w*q->w - (f->x*q->x + f->y*q->y + f->z*q->z);
  t.x= f->y*q->z - f->z*q->y + f->w*q->x + f->x*q->w;
  t.y= f->z*q->x - f->x*q->z + f->w*q->y + f->y*q->w;
  t.z= f->x*q->y - f->y*q->x + f->w*q->z + f->z*q->w;
  quat_copy(q,&t);
  return q;    
}

Quat* quat_identity( Quat* q )
{
  if (!q) q= quat_alloc();
  q->w= 1.0;
  q->x= q->y= q->z= 0.0;
  return q;
}

Quat* quat_conjugate( Quat* q )
{
  q->x= -q->x;
  q->y= -q->y;
  q->z= -q->z;
  return q;
}

Quat* quat_create( Quat* q, double x, double y, double z, int w_pos )
{
  if (!q) q= quat_alloc();
  q->x= x;
  q->y= y;
  q->z= z;
  q->w= sqrt(1.0 - (x*x + y*y + z*z));
  if (!w_pos) q->w= -q->w;
  return q;
}

Quat* quat_from_axis_angle( Quat* q, double x, double y, double z, 
			    double theta )
{
  double sin_a;
  double cos_a;

  if (!q) q= quat_alloc();

  sin_a= sin(0.5*theta);
  cos_a= cos(0.5*theta);

  q->x= x * sin_a;
  q->y= y * sin_a;
  q->z= z * sin_a;
  q->w= cos_a;
  quat_normalize( q );

  return q;
}

void quat_to_axis_angle( Quat* q, double* x, double* y, double* z, 
			 double* theta ) {
  Quat qnorm;
  double len;

  quat_copy(&qnorm,q);
  quat_normalize(&qnorm);

  len= sqrt(qnorm.x*qnorm.x + qnorm.y*qnorm.y + qnorm.z*qnorm.z);

  if (len>0.0) {
    double len_inv= 1.0/len;
    *x= qnorm.x*len_inv;
    *y= qnorm.y*len_inv;
    *z= qnorm.z*len_inv;
  }
  else {
    *theta= 0.0;
    *x= *y= 0.0;
    *z= 1.0;
  }

  *theta= 2.0*acos(qnorm.w);
}

Quat* quat_from_euler_RzRyRx( Quat* q, double x_angle, double y_angle,
			      double z_angle )
{
  Quat qx;
  Quat qy;
  Quat qz;

  if (!q) q= quat_alloc();

  quat_from_axis_angle( &qx, 1.0, 0.0, 0.0, x_angle );
  quat_from_axis_angle( &qy, 0.0, 1.0, 0.0, y_angle );
  quat_from_axis_angle( &qz, 0.0, 0.0, 1.0, z_angle );

  quat_copy(q, &qz);
  quat_mult_right(q, &qy);
  quat_mult_right(q, &qx);

  return q;
}

/*
 * This routine finds Euler angles xout, yout, and zout such that the series of
 * axis-aligned rotations RzRyRx for the given angles is equal to the input 
 * quaternion rotation q.  Returns 1 on success, zero on failure.
 */
int quat_to_euler_RzRyRx( Quat* q, double* xout, double* yout,
			  double* zout )
{
  double theta_x=0.0, theta_y=0.0, theta_z=0.0;
  Quat Rx, Ry, Rz, temp;
  double x,y,z,phi;
  int iter= 0;

  if (DEBUG) fprintf(stderr,
		     "quat_to_euler_RzRyRx: decomposing (%g %g %g %g)\n",
		     q->x,q->y,q->z,q->w);

  do {
    (void)quat_from_axis_angle(&Rx,1.0,0.0,0.0,theta_x);
    (void)quat_from_axis_angle(&Ry,0.0,1.0,0.0,theta_y);
    (void)quat_from_axis_angle(&Rz,0.0,0.0,1.0,theta_z);

    /* Calculate RzRyRxRxRyRz */
    (void)quat_copy(&temp,q);
    (void)quat_conjugate(&temp);
    (void)quat_mult_right(&temp, &Rx);
    (void)quat_mult_right(&temp, &Ry);
    (void)quat_mult_right(&temp, &Rz);

    /* Make a new guess */
    quat_to_axis_angle(&temp,&x,&y,&z,&phi);
    /*
      fprintf(stderr,"new guess %g %g %g (radians)-> %g about (%g,%g,%g)\n",
      theta_x,theta_y,theta_z,phi,x,y,z);
    */
    theta_x -= x*phi;
    theta_y -= y*phi;
    theta_z -= z*phi;
  } while (phi!=0.0 && iter++ < EULER_DECOMP_MAX_ITER);

  if (iter>= EULER_DECOMP_MAX_ITER) return 0; /* convergence failure */

  if (DEBUG) {
    fprintf(stderr,"decomp results: %g %g %g after %d iterations\n", 
	    theta_x, theta_y, theta_z, iter);
    (void)quat_from_axis_angle(&Rx,1.0,0.0,0.0,theta_x);
    (void)quat_from_axis_angle(&Ry,0.0,1.0,0.0,theta_y);
    (void)quat_from_axis_angle(&Rz,0.0,0.0,1.0,theta_z);
    (void)quat_identity(&temp);
    (void)quat_mult_right(&temp, &Rx);
    (void)quat_mult_right(&temp, &Ry);
    (void)quat_mult_right(&temp, &Rz);
    fprintf(stderr,"goal: (%g %g %g %g) vs. result (%g %g %g %g)\n",
	    q->x,q->y,q->z,q->w,temp.x,temp.y,temp.z,temp.w);
  }

  *xout= theta_x; 
  *yout= theta_y;
  *zout= theta_z;
  return 1;
}

void trans_mult_right( Transform t, const Transform factor )
{
  int i, row, column;
  Transform newT;

  for (row = 0;row<4;row++)
    for (column= 0;column<4;column++) {
      newT[(4*row)+column] = t[(4*row)] * factor[column];
      newT[(4*row)+column] += t[(4*row)+1] * factor[(4*1)+column];
      newT[(4*row)+column] += t[(4*row)+2] * factor[(4*2)+column];
      newT[(4*row)+column] += t[(4*row)+3] * factor[(4*3)+column];
    }

  for (i=0; i<16; i++) t[i]= newT[i];
}

void trans_mult_left( const Transform factor, Transform t )
{
  int i, row, column;
  Transform newT;

  for (row = 0;row<4;row++)
    for (column= 0;column<4;column++) {
      newT[(4*row)+column] = factor[(4*row)] * t[column];
      newT[(4*row)+column] += factor[(4*row)+1] * t[(4*1)+column];
      newT[(4*row)+column] += factor[(4*row)+2] * t[(4*2)+column];
      newT[(4*row)+column] += factor[(4*row)+3] * t[(4*3)+column];
    }

  for (i=0; i<16; i++) t[i]= newT[i];
}

void trans_transpose( Transform t )
{
  /* Transpose this transform in place */
  register float temp;
  register int i,j;

  for (i=1; i<4; i++)
    for (j=0; j<i; j++) {
        temp= t[4*i +j];
        t[4*i + j]= t[4*j +i];
        t[4*j + i]= temp;
      };
}

int trans_inverse( Transform invT, const Transform t )
{
  int four= 4;
  int sixteen= 16;
  double work[16];
  int pivots[4];
  int info;
  
  trans_copy( invT, t );
  DGETRF( &four, &four, invT, &four, pivots, &info );
  if (info>0) return 0; /* failure */
  else if (info<0)
    Abort("trans_inverse: internal error; arg %d to DGETRF invalid!\n",-info);
  DGETRI( &four, invT, &four, pivots, work, &sixteen, &info );
  if (info==0) return 1; /* success */
  else if (info>0) return 0; /* failure */
  else Abort("trans_inverse: internal error; arg %d to DGETRI invalid!\n",-info);
  return 0; /* to keep compilers happy */
}

void trans_dump( FILE* ofile, Transform t )
{
  fprintf(ofile,"  ( %5f %5f %5f %5f )\n",t[0],t[1],t[2],t[3]);
  fprintf(ofile,"  ( %5f %5f %5f %5f )\n",t[4],t[5],t[6],t[7]);
  fprintf(ofile,"  ( %5f %5f %5f %5f )\n",t[8],t[9],t[10],t[11]);
  fprintf(ofile,"  ( %5f %5f %5f %5f )\n",t[12],t[13],t[14],t[15]);
}

void trans_vec_mult( const Transform factor, Vec4 v )
{
  Vec4 newV;
  int row;

  for (row=0; row<4; row++) {
    newV[row] = factor[4*row]*v[0];
    newV[row] += factor[4*row+1]*v[1];
    newV[row] += factor[4*row+2]*v[2];
    newV[row] += factor[4*row+3]*v[3];
  }
  for (row=0; row<4; row++) v[row]= newV[row];
}


Quat* quat_nrm_sqrt( Quat* result, const Quat* q_in )
{
  if (q_in->w >= 0.0) {
    double k= 1.0/sqrt(2.0*(1.0+q_in->w));
    double rsqr= q_in->x*q_in->x + q_in->y*q_in->y + q_in->z*q_in->z;
    result->w= sqrt(1.0-k*k*rsqr);
    result->x= k*q_in->x;
    result->y= k*q_in->y;
    result->z= k*q_in->z;
    return result;
  }
  else {
    Quat result_lcl;
    Quat tmp;
    Quat tmpsqrt;
    double mag;
    double v[3];

    /* take sqrt of negative of input */
    tmp.x= -q_in->x; tmp.y= -q_in->y; tmp.z= -q_in->z; tmp.w= -q_in->w; 
    quat_nrm_sqrt( &tmpsqrt, &tmp );

    /* We must now shuffle parts to build the sqrt of the original quat. */
    
    /* Assemble result by multiplying that sqrt by (imaginary) i */
    mag= sqrt(tmpsqrt.x*tmpsqrt.x + tmpsqrt.y*tmpsqrt.y + tmpsqrt.z*tmpsqrt.z);
    if (mag==0.0) {
      /* Symmetry condition; pick a default direction */
      v[0]= 1.0;
      v[1]= 0.0;
      v[2]= 0.0;
    }
    else {
      v[0]= tmpsqrt.x/mag;
      v[1]= tmpsqrt.y/mag;
      v[2]= tmpsqrt.z/mag;
    }
    result_lcl.x= tmpsqrt.w*v[0];
    result_lcl.y= tmpsqrt.w*v[1];
    result_lcl.z= tmpsqrt.w*v[2];
    result_lcl.w= -mag;

    if (DEBUG) {
      Quat prod;
      quat_copy(&prod,&result_lcl);
      quat_mult_right(&prod,&result_lcl);
      fprintf(stderr,
	      "quat_nrm_sqrt: result_lcl (%g %g %g %g)\n",
	      result_lcl.x,result_lcl.y,result_lcl.z,result_lcl.w);
      fprintf(stderr,
	      "           goal (%g %g %g %g) vs. result (%g %g %g %g)\n", 
	      q_in->x,q_in->y,q_in->z,q_in->w,
	      prod.x,prod.y,prod.z,prod.w);
    }
    
    *result= result_lcl;
    return result;
  }
}
