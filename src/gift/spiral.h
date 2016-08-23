/*
 *	spiral.h  - Constants defining byte offsets
 *			within the raw header file
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
 *		12/95 - Written by Greg Hood (PSC)
 */

#define RAW_HDR_SIZE		32	/* total size of the raw
					   header in bytes */

/* all of the following fields are 32-bit integers */
#define RAW_HDR_NCOILS		0	/* number of coils */
#define RAW_HDR_NSLICES		4	/* number of slices */
#define RAW_HDR_NPHASES		8	/* number of phases */
#define RAW_HDR_NPR		12	/* number of projections */
#define RAW_HDR_NDAT		16	/* number of data points per projection */
#define RAW_HDR_RES		20	/* x and y resolution of images in pixels */
#define RAW_HDR_OVER_SAMP	24	/* oversampling ratio */
#define RAW_HDR_GRID_LEN	28	/* grid length */

