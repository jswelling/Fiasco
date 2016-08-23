/*
 *	fiasco.h - Offsets for FIASCO header files
 *
 *	Copyright (c) 1995,1996 Pittsburgh Supercomputing Center
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
 *		1/96	- Prefixed constants with OLD_ (Hood)
 *		4/96	- Removed OLD_ prefixes :-) (Hood)
 */

#define FIASCO_TYPE		0
#define FIASCO_NDIMS		4
#define FIASCO_VEC		8
#define FIASCO_X		12
#define FIASCO_Y		16
#define FIASCO_Z		20
#define FIASCO_T		24

/* the number stored in FIASCO_N_DIMS must be multiplied by 4 and added to
   following offsets to get the true offset */
#define FIASCO_FILENAME		12
#define FIASCO_MISSING		524


/* the following are valid values to store at the FIASCO_TYPE location */
#define FIASCO_MRI_CHAR		0
#define FIASCO_MRI_LONG		1
#define FIASCO_MRI_FLOAT	2
#define FIASCO_MRI_DOUBLE	3
#define FIASCO_MRI_SHORT	4

Boolean FiascoCheckFormat ();
void FiascoStartReading ();
void FiascoReadImage ();
void FiascoEndReading ();

void FiascoStartWriting ();
void FiascoWriteImage ();
void FiascoEndWriting ();
