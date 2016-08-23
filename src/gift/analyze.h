/*
 *	analyze.h - Offsets for ANALYZE header files
 *
 *	Copyright (c) 1995 Pittsburgh Supercomputing Center
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
 *
 *	HISTORY
 *		12/95	- Written by Greg Hood (PSC)
 */

#define ANALYZE_HEADER_SIZE		348

#define ANALYZE_SIZEOF_HDR		0
#define ANALYZE_EXTENTS			32
#define ANALYZE_REGULAR			38
#define ANALYZE_N_DIMS			40
#define ANALYZE_X			42
#define ANALYZE_Y			44
#define ANALYZE_Z			46
#define ANALYZE_T			48
#define ANALYZE_DATATYPE		70
#define ANALYZE_BITPIX			72
#define ANALYZE_VOXEL_X_SIZE		80
#define ANALYZE_VOXEL_Y_SIZE		84
#define ANALYZE_VOXEL_Z_SIZE		88
#define ANALYZE_T_STEP			92
#define ANALYZE_GLMAX			140
#define ANALYZE_GLMIN			144


#define ANALYZE_DATATYPE_UNKNOWN	0
#define ANALYZE_DATATYPE_UINT1		1
#define ANALYZE_DATATYPE_UINT8		2
#define ANALYZE_DATATYPE_INT16		4
#define ANALYZE_DATATYPE_INT32		8
#define ANALYZE_DATATYPE_FLOAT32	16
#define ANALYZE_DATATYPE_COMPLEX64	32
#define ANALYZE_DATATYPE_FLOAT64	64


Boolean AnalyzeCheckFormat ();
void AnalyzeStartReading ();
void AnalyzeReadImage ();
void AnalyzeEndReading ();

void AnalyzeStartWriting ();
void AnalyzeWriteImage ();
void AnalyzeEndWriting ();
