/************************************************************
 *                                                          *
 *  vpolygon.h                                               *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1997 Department of Statistics,         *
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
 *  Original programming by Joel Welling, 2/03              *
 ************************************************************/

typedef struct voronoi_wing_struct {
  double vtx[2]; /* location of an associated vertex */
  double theta; /* atan2(deltaY,deltaX) for sorting purposes */
  int keep;     /* zero if one of the forming points of this verts is bogus */
} VWing;

/* a structure to hold a *CONVEX* polyhedron */
typedef struct voroni_poly_struct {
  double ctr[2]; /* A point within the polygon */
  VWing* wings;
  int nWings;
  int wingBufSize;
  double maxDist;
  double maxDistSqr;
} VPoly;

extern void vply_setDebug(const int i);
extern void vply_dump(VPoly* vp, FILE* ofile);
extern VPoly* vply_create(double ctrX, double ctrY, int initialN);
extern void vply_destroy(VPoly* vp);
extern void vply_add(VPoly* vp, double x, double y);
extern void vply_sort(VPoly* vp);
extern double vply_calcArea( VPoly* vp );
extern void vply_writePS(VPoly* vp);
extern void vply_prepPS(char* fname);
extern void vply_finishPS();
extern double vply_getMaxDistance( VPoly* vp );

/* These may free vp! Call it like: vp= vply_clip(vp, clipAgainst) */
extern VPoly* vply_clip(VPoly* vp, const VPoly* clipAgainst);
extern VPoly* vply_clipToCircle(VPoly* vp, 
				const double xCtr, const double yCtr,
				const double radius);
