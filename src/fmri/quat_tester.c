/************************************************************
 *                                                          *
 *  quat_tester.c                                           *
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "quaternion.h"

#define DEG2RAD (M_PI/180.0)
#define RAD2DEG (180.0/M_PI)

int main( int argc, char* argv[] ) 
{
  Quat q;
  Quat q_conj;
  Quat q_sqrt;
  Transform t;
  double x,y,z,theta;

  if (argc != 1) {
    fprintf(stderr,"Usage: %s\n",argv[0]);
    fprintf(stderr,"    You will be prompted for Euler angles RzRyRx,\n");
    fprintf(stderr,"    or 'q' followed by quaternion components.\n");
    fprintf(stderr,
  "    An empty line causes previously given quats to be be multiplied.\n");
    exit(-1);
  }

  quat_identity(&q);

  while (1) {
    double tx, ty, tz;
    int count;
    char inbuf[81];
    printf("Enter 3 Euler angles X, Y, Z (used as RzRyRx) in degrees\n");
    printf("or 'q' followed by 4 quaternion components.\n");
    fgets(inbuf,80,stdin);
    if (strlen(inbuf)<2) break;
    else if (inbuf[0]=='q') {
      Quat r;
      if ((count=sscanf(inbuf+1,"%lf %lf %lf %lf\n",
			&r.x,&r.y,&r.z,&r.w)) == 4) {
	int neg = (r.w<0.0);
	r.w= sqrt( 1.0-(r.x*r.x + r.y*r.y + r.z*r.z) );
	if (neg) r.w *= -1.0;
      quat_mult_right(&q,&r);
      }
      else printf("Try again.\n");
    }
    else if ((count=sscanf(inbuf,"%lf %lf %lf\n",&tx,&ty,&tz)) == 3) {
      Quat r;
      quat_from_euler_RzRyRx(&r,DEG2RAD*tx,DEG2RAD*ty,DEG2RAD*tz);
      quat_mult_right(&q,&r);
    }
    else printf("Try again.\n");
  }

  printf("Resulting quaternion: (%g %g %g %g)\n",q.x,q.y,q.z,q.w);
  quat_to_axis_angle(&q,&x,&y,&z,&theta);
  printf("Axis-angle representation: %f degrees about (%f,%f,%f)\n",
	 RAD2DEG*theta,x,y,z);
  quat_nrm_sqrt(&q_sqrt,&q);
  printf("Square root: (%g %g %g %g)\n",q_sqrt.x,q_sqrt.y,q_sqrt.z,q_sqrt.w);
  quat_mult_right(&q_sqrt,&q_sqrt);
  printf("Square root squared: (%g %g %g %g)\n",
	 q_sqrt.x,q_sqrt.y,q_sqrt.z,q_sqrt.w);
  printf("Equivalent rotation matrix:\n");
  quat_to_trans(t,&q,0.0,0.0,0.0);
  trans_dump(stdout,t);
  quat_to_euler_RzRyRx( &q, &x, &y, &z );
  printf("Decomposing back to Euler angles (may not match): %g %g %g\n",
	 RAD2DEG*x, RAD2DEG*y, RAD2DEG*z);
  trans_to_quat(&q_conj,t);
  printf("Changing back to quaternion: (%g %g %g %g)\n",q_conj.x,q_conj.y,
	 q_conj.z,q_conj.w);
  quat_conjugate(&q_conj);
  quat_mult_right(&q,&q_conj);
  printf("After multiplication by conjugate: get (%f %f %f %f)\n",
	 q.x,q.y,q.z,q.w);
  quat_to_trans(t,&q,0.0,0.0,0.0);
  trans_dump(stdout,t);

}
