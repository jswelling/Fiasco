/************************************************************
 *                                                          *
 *  rttraj.h                                                *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 2003 Department of Statistics,         *
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
 *  Original programming by Joel Welling March 2003         *
 ************************************************************/
/* This library supplies utilities related to pulse         *
 * sequence calculations in k-space.                        */

/* Older real-time spiral trajectory, as used by "sf11" pulse sequences */
extern int getrttrajvd (int npts, int nl, double ts, double gts,
			double angle, double fsgcm, double opfov,
			int risetime, int densamp, double kxscale,
			double kyscale, float **kx, float **ky,
			int* res_out);

/* Newer real-time spiral trajectory, as used by "splx" and "sprl" 
 * pulse sequences.
 */
extern int getrttrajghg(int opxres, int nl, double dts, double gts, 
			double gamp, double opfov,
                        double gslew, double Tmax, int gtype, 
			float** kx, float**ky );


