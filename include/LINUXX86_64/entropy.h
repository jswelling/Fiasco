/************************************************************
 *                                                          *
 *  entropy.h                                            *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1995 Department of Statistics,         *
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
 *  Original programming by Chris Hefferan 6/04             *
 ************************************************************/

/* "$Id: entropy.h,v 1.5 2007/04/12 23:49:23 welling Exp $" */

typedef struct entropy_ctx_struct {
  long* histogram;
  long histogramArraySize;
  int debugFlag;
  int verboseFlag;
  long num_of_boxes;
  double min;
  double max;
  int min_set;
  int max_set;
} EntropyContext;

typedef struct mutual_info_ctx_struct {
  long* histogram;
  long histogramArraySize;
  long num_of_boxes;
  int debugFlag;
  int verboseFlag;
  EntropyContext* ec1;
  EntropyContext* ec2;
} MutualInfoContext;

EntropyContext* ent_createContext(void);
void ent_destroyContext(EntropyContext* ctx);
void ent_setDebug(EntropyContext* ctx, int val);
void ent_setVerbose(EntropyContext* ctx, int val);
int ent_getDebug(const EntropyContext* ctx);
int ent_getVerbose(const EntropyContext* ctx);
void ent_setMin(EntropyContext* ctx, double val);
void ent_setMax(EntropyContext* ctx, double val);
double ent_getMin(const EntropyContext* ctx);
double ent_getMax(const EntropyContext* ctx);
int ent_getMinSet(const EntropyContext* ctx);
int ent_getMaxSet(const EntropyContext* ctx);
void ent_unsetMin(EntropyContext* ctx);
void ent_unsetMax(EntropyContext* ctx);

/* NOTE: unsetNBins(ctx) is accomplished by doing setNBins(0) */
void ent_setNBins(EntropyContext* ctx, const long nbins);
long ent_getNBins(const EntropyContext* ctx);
int ent_getNBinsSet(const EntropyContext* ctx);

MutualInfoContext* ent_createMIContext(void);
void ent_destroyMIContext(MutualInfoContext* ctx);
void ent_setMIDebug(MutualInfoContext* ctx, int val);
void ent_setMIVerbose(MutualInfoContext* ctx, int val);
int ent_getMIDebug(const MutualInfoContext* ctx);
int ent_getMIVerbose(const MutualInfoContext* ctx);
void ent_setMIMin1(MutualInfoContext* ctx, double val);
void ent_setMIMax1(MutualInfoContext* ctx, double val);
double ent_getMIMin1(const MutualInfoContext* ctx);
double ent_getMIMax1(const MutualInfoContext* ctx);
int ent_getMIMin1Set(const MutualInfoContext* ctx);
int ent_getMIMax1Set(const MutualInfoContext* ctx);
void ent_unsetMIMin1(MutualInfoContext* ctx);
void ent_unsetMIMax1(MutualInfoContext* ctx);
void ent_setMIMin2(MutualInfoContext* ctx, double val);
void ent_setMIMax2(MutualInfoContext* ctx, double val);
double ent_getMIMin2(const MutualInfoContext* ctx);
double ent_getMIMax2(const MutualInfoContext* ctx);
int ent_getMIMin2Set(const MutualInfoContext* ctx);
int ent_getMIMax2Set(const MutualInfoContext* ctx);
void ent_unsetMIMin2(MutualInfoContext* ctx);
void ent_unsetMIMax2(MutualInfoContext* ctx);

/* NOTE: unsetNBins(ctx) is accomplished by doing setNBins(0) */
void ent_setMINBins(MutualInfoContext* ctx, const long nbins);
long ent_getMINBins(const MutualInfoContext* ctx);
int ent_getMINBinsSet(const MutualInfoContext* ctx);

double ent_calcMaskedImageEntropyDouble( EntropyContext* ctx,
					 const double* img, const int* mask, 
					 long dx, long dy, long dz, 
					 long stride, long maskstride);
double ent_calcMaskedImageEntropyFloat( EntropyContext* ctx,
					 const float* img, const int* mask, 
					 long dx, long dy, long dz, 
					 long stride, long maskstride);
double ent_calcImageEntropyDouble( EntropyContext* ctx,
				   const double* img, long dx, 
				   long dy, long dz, long stride );
double ent_calcImageEntropyFloat( EntropyContext* ctx,
				   const float* img, long dx, 
				  long dy, long dz, long stride );


double ent_calcMaskedJointEntropyDouble( MutualInfoContext* mc,
					 const double* img1,
					 const double* img2,
					 const int* mask,
					 long dx, long dy, long dz, 
					 long stride1, long stride2,
					 long maskstride );
double ent_calcMaskedJointEntropyFloat( MutualInfoContext* mc,
					const float* img1,
					const float* img2,
					const int* mask,
					long dx, long dy, long dz, 
					long stride1, long stride2,
					long maskstride );
double ent_calcJointEntropyDouble( MutualInfoContext* mc,
				   const double* img1,
				   const double* img2,
				   long dx, long dy, long dz, 
				   long stride1, long stride2 );
double ent_calcJointEntropyFloat( MutualInfoContext* mc,
				   const float* img1,
				   const float* img2,
				   long dx, long dy, long dz, 
				   long stride1, long stride2 );
double ent_calcMaskedMutualInformationDouble( MutualInfoContext* mc,
					      const double* img1, 
					      const double* img2, 
					      const int* mask, 
					      long dx, long dy, long dz,
					      long stride1, long stride2, 
					      long maskstride );
double ent_calcMaskedMutualInformationFloat( MutualInfoContext* mc,
					     const float* img1, 
					     const float* img2, 
					     const int* mask, 
					     long dx, long dy, long dz,
					     long stride1, long stride2,
					     long maskstride );
double ent_calcMutualInformationDouble( MutualInfoContext* mc,
					const double* img1, 
					const double* img2, 
					long dx, long dy, long dz,
					long stride1, long stride2);
double ent_calcMutualInformationFloat( MutualInfoContext* mc,
					const float* img1, 
					const float* img2, 
					long dx, long dy, long dz,
					long stride1, long stride2);

