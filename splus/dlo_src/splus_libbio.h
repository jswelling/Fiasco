/*
 *	Binary I/O utility functions
 *
 *	Copyright (c) 1996 Pittsburgh Supercomputing Center
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
 *		1/96 Written by Greg Hood (PSC)
 */

extern int bio_big_endian_machine;	/* if 1, then the internal representation
					   on this machine is big-endian;
					   the bio.c package will set this according
					   to the machine type */
extern int bio_big_endian_input;	/* may be set by user of this package to
					   control endian-ness of input operations */
extern int bio_big_endian_output;	/* may be set by user of this package to
					   control endian-ness of output operations */

extern int bio_error;		/* the bio package will set this variable to 1
				   if an I/O error occurs during file reading
				   or writing; the user may reset this flag
				   to 0 at any time */

/* Functions that read/write to a location in memory */
int BRdUInt8 (unsigned char *addr);
void BWrUInt8 (unsigned char *addr, int v);
int BRdInt8 (unsigned char *addr);
void BWrInt8 (unsigned char *addr, int v);
int BRdInt16 (unsigned char *addr);
void BWrInt16 (unsigned char *addr, int v);
int BRdInt32 (unsigned char *addr);
void BWrInt32 (unsigned char *addr, int v);
float BRdFloat32 (unsigned char *addr);
void BWrFloat32 (unsigned char *addr, float v);
double BRdFloat64 (unsigned char *addr);
void BWrFloat64 (unsigned char *addr, double v);

/* Functions that read/write arrays to a buffer in memory */
void BRdUInt8Array (unsigned char *buf, unsigned char *a, int n);
void BWrUInt8Array (unsigned char *buf, unsigned char *a, int n);
void BRdInt8Array (unsigned char *buf, signed char *a, int n);
void BWrInt8Array (unsigned char *buf, signed char *a, int n);
void BRdInt16Array (unsigned char *buf, short *a, int n);
void BWrInt16Array (unsigned char *buf, short *a, int n);
void BRdInt32Array (unsigned char *buf, int *a, int n);
void BWrInt32Array (unsigned char *buf, int *a, int n);
void BRdFloat32Array (unsigned char *buf, float *a, int n);
void BWrFloat32Array (unsigned char *buf, float *a, int n);
void BRdFloat64Array (unsigned char *buf, double *a, int n);
void BWrFloat64Array (unsigned char *buf, double *a, int n);

/* Functions that read/write to a file */
int FRdUInt8 (FILE *stream);
void FWrUInt8 (FILE *stream, int v);
int FRdInt8 (FILE *stream);
void FWrInt8 (FILE *stream, int v);
int FRdInt16 (FILE *stream);
void FWrInt16 (FILE *stream, int v);
int FRdInt32 (FILE *stream);
void FWrInt32 (FILE *stream, int v);
float FRdFloat32 (FILE *stream);
void FWrFloat32 (FILE *stream, float v);
double FRdFloat64 (FILE *stream);
void FWrFloat64 (FILE *stream, double v);

/* Functions that read/write arrays to a file */
void FRdUInt8Array (FILE *stream, unsigned char *a, int n);
void FWrUInt8Array (FILE *stream, unsigned char *a, int n);
void FRdInt8Array (FILE *stream, signed char *a, int n);
void FWrInt8Array (FILE *stream, signed char *a, int n);
void FRdInt16Array (FILE *stream, short *a, int n);
void FWrInt16Array (FILE *stream, short *a, int n);
void FRdInt32Array (FILE *stream, int *a, int n);
void FWrInt32Array (FILE *stream, int *a, int n);
void FRdFloat32Array (FILE *stream, float *a, int n);
void FWrFloat32Array (FILE *stream, float *a, int n);
void FRdFloat64Array (FILE *stream, double *a, int n);
void FWrFloat64Array (FILE *stream, double *a, int n);
