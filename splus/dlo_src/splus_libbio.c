/*
 *	Binary I/O utility functions
 *
 *	Copyright (c) 1996  Pittsburgh Supercomputing Center
# *                                                          *
# *  This program is distributed in the hope that it will    *
# *  be useful, but WITHOUT ANY WARRANTY; without even the   *
# *  implied warranty of MERCHANTABILITY or FITNESS FOR A    *
# *  PARTICULAR PURPOSE.  Neither Carnegie Mellon University *
# *  nor any of the authors assume any liability for         *
# *  damages, incidental or otherwise, caused by the         *
# *  installation or use of this software.                   *
# *                                                          *
# *  CLINICAL APPLICATIONS ARE NOT RECOMMENDED, AND THIS     *
# *  SOFTWARE HAS NOT BEEN EVALUATED BY THE UNITED STATES    *
# *  FDA FOR ANY CLINICAL USE.                               *
# *                                                          *
 *
 *	HISTORY
 *		1/96 Written by Greg Hood
 *		3/96 Minor modifications to use with libmri (Hood)
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "splus_libbio.h"

static char rcsid[] = "$Id: splus_libbio.c,v 1.4 2000/08/08 00:50:47 welling Exp $";

/* The following horrible macros and defs work around restrictions on
 * memory allocation in Splus.  See the S-plus Programmer's Manual,
 * section 9.10.1 (Allocating Storage) for details.  Note that
 * calls to realloc() must be explicitly recoded.  Can't seem to
 * include S.h because of a function proto conflict.
 */
char *S_alloc(long n, int size);
char*S_realloc( char* p, long new, long old, int size );
#define malloc( size ) S_alloc(size,1)
#define free( pointer ) /* free is a no-op */

/* Unless otherwise specified, we do all I/O in native-endian format */

#if defined(LITTLE_ENDIAN)
int bio_big_endian_machine = 0;
int bio_big_endian_input = 0;
int bio_big_endian_output = 0;
#else
int bio_big_endian_machine = 1;
int bio_big_endian_input = 1;
int bio_big_endian_output = 1;
#endif

/* This flag will be set to 1 upon detection of any I/O error */
int bio_error = 0;


int
BRdUInt8 (unsigned char *addr)
{
  return(*((unsigned char *) addr));
}

void
BWrUInt8 (unsigned char *addr,
	  int v)
{
  *((unsigned char *) addr) = v;
}

int
BRdInt8 (unsigned char *addr)
{
  return(*((signed char *) addr));
}

void
BWrInt8 (unsigned char *addr,
	 int v)
{
  *((signed char *) addr) = v;
}

int
BRdInt16 (unsigned char *addr)
{
  if (bio_big_endian_input)
    return((((signed char) addr[0]) << 8) |
	   addr[1]);
  else
    return((((signed char) addr[1]) << 8) |
	   addr[0]);
}

void
BWrInt16 (unsigned char *addr,
	  int v)
{
  if (bio_big_endian_output)
    {
      addr[0] = (v >> 8) & 0xff;
      addr[1] = v & 0xff;
    }
  else
    {
      addr[0] = v & 0xff;
      addr[1] = (v >> 8) & 0xff;
    }
}

int
BRdInt32 (unsigned char *addr)
{
  if (bio_big_endian_input)
    return((((signed char) addr[0]) << 24) |
	   (addr[1] << 16) |
	   (addr[2] << 8) |
	   addr[3]);
  else
    return((((signed char) addr[3]) << 24) |
	   (addr[2] << 16) |
	   (addr[1] << 8) |
	   addr[0]);
}

void
BWrInt32 (unsigned char *addr,
	  int v)
{
  if (bio_big_endian_output)
    {
      addr[0] = (v >> 24) & 0xff;
      addr[1] = (v >> 16) & 0xff;
      addr[2] = (v >> 8) & 0xff;
      addr[3] = v & 0xff;
    }
  else
    {
      addr[0] = v & 0xff;
      addr[1] = (v >> 8) & 0xff;
      addr[2] = (v >> 16) & 0xff;
      addr[3] = (v >> 24) & 0xff;
    }
}

float
BRdFloat32 (unsigned char *addr)
{
  union {
    unsigned char temp[4];
    float val;
  } u;

  if (bio_big_endian_machine ^ bio_big_endian_input)
    {
      u.temp[0] = addr[3];
      u.temp[1] = addr[2];
      u.temp[2] = addr[1];
      u.temp[3] = addr[0];
    }
  else
    {
      u.temp[0] = addr[0];
      u.temp[1] = addr[1];
      u.temp[2] = addr[2];
      u.temp[3] = addr[3];
    }
  return (u.val);
}

void
BWrFloat32 (unsigned char *addr,
	    float v)
{
  unsigned char *temp;

  temp = (unsigned char *) &v;
  if (bio_big_endian_machine ^ bio_big_endian_output)
    {
      addr[0] = temp[3];
      addr[1] = temp[2];
      addr[2] = temp[1];
      addr[3] = temp[0];
    }
  else
    {
      addr[0] = temp[0];
      addr[1] = temp[1];
      addr[2] = temp[2];
      addr[3] = temp[3];
    }
}

double
BRdFloat64 (unsigned char *addr)
{
  union {
    unsigned char temp[8];
    double val;
  } u;

  if (bio_big_endian_machine ^ bio_big_endian_input)
    {
      u.temp[0] = addr[7];
      u.temp[1] = addr[6];
      u.temp[2] = addr[5];
      u.temp[3] = addr[4];
      u.temp[4] = addr[3];
      u.temp[5] = addr[2];
      u.temp[6] = addr[1];
      u.temp[7] = addr[0];
    }
  else
    memcpy(u.temp, addr, 8);
  return (u.val);
}

void
BWrFloat64 (unsigned char *addr,
	   double v)
{
  unsigned char *temp;

  temp = (unsigned char *) &v;
  if (bio_big_endian_machine ^ bio_big_endian_output)
    {
      addr[0] = temp[7];
      addr[1] = temp[6];
      addr[2] = temp[5];
      addr[3] = temp[4];
      addr[4] = temp[3];
      addr[5] = temp[2];
      addr[6] = temp[1];
      addr[7] = temp[0];
    }
  else
    memcpy(addr, temp, 8);
}

void
BRdUInt8Array (unsigned char *buf,
	       unsigned char *a,
	       int n)
{
  memcpy(a, buf, n);
}

void
BWrUInt8Array (unsigned char *buf,
	      unsigned char *a,
	      int n)
{
  memcpy(buf, a, n);
}

void
BRdInt8Array (unsigned char *buf,
	      signed char *a,
	      int n)
{
  memcpy(a, buf, n);
}

void
BWrInt8Array (unsigned char *buf,
	      signed char *a,
	      int n)
{
  memcpy(buf, a, n);
}

void
BRdInt16Array (unsigned char *buf,
	       short *a,
	       int n)
{
  int i, len;
  unsigned char *dest, temp;

  if (sizeof(short) == 2)
    {
      /* copy directly into the destination */
      memcpy(a, buf, 2*n);
      if (bio_big_endian_machine ^ bio_big_endian_input)
	{
	  /* swap bytes */
	  dest = (unsigned char *) a;
	  len = 2*n;
	  for (i = 0; i < len; i+=2)
	    {
	      temp = dest[i];
	      dest[i] = dest[i+1];
	      dest[i+1] = temp;
	    }
	}
    }
  else
    for (i = 0; i < n; ++i)
      a[i] = BRdInt16(&buf[2*i]);
}

void
BWrInt16Array (unsigned char *buf,
	       short *a,
	       int n)
{
  int i, len;
  unsigned char temp;

  if (sizeof(short) == 2)
    {
      /* copy directly into the destination */
      memcpy(buf, a, 2*n);
      if (bio_big_endian_machine ^ bio_big_endian_output)
	{
	  /* swap bytes */
	  len = 2*n;
	  for (i = 0; i < len; i += 2)
	    {
	      temp = buf[i];
	      buf[i] = buf[i+1];
	      buf[i+1] = temp;
	    }
	}
    }
  else
    for (i = 0; i < n; ++i)
      BWrInt16(&buf[2*i], a[i]);
}

void
BRdInt32Array (unsigned char *buf,
	       int *a,
	       int n)
{
  int i, len;
  unsigned char *dest, temp;

  if (sizeof(int) == 4)
    {
      /* copy directly into the destination */
      memcpy(a, buf, 4*n);
      if (bio_big_endian_machine ^ bio_big_endian_input)
	{
	  /* reverse byte order */
	  dest = (unsigned char *) a;
	  len = 4*n;
	  for (i = 0; i < len; i+=4)
	    {
	      temp = dest[i];
	      dest[i] = dest[i+3];
	      dest[i+3] = temp;
	      temp = dest[i+2];
	      dest[i+2] = dest[i+1];
	      dest[i+1] = temp;
	    }
	}
    }
  else
    for (i = 0; i < n; ++i)
      a[i] = BRdInt32(&buf[4*i]);
}

void
BWrInt32Array (unsigned char *buf,
	       int *a,
	       int n)
{
  int i, len;
  unsigned char temp;

  if (sizeof(int) == 4)
    {
      /* copy directly into the destination */
      memcpy(buf, a, 4*n);
      if (bio_big_endian_machine ^ bio_big_endian_output)
	{
	  /* reverse byte order */
	  len = 4*n;
	  for (i = 0; i < len; i += 4)
	    {
	      temp = buf[i];
	      buf[i] = buf[i+3];
	      buf[i+3] = temp;
	      temp = buf[i+2];
	      buf[i+2] = buf[i+1];
	      buf[i+1] = temp;
	    }
	}
    }
  else
    for (i = 0; i < n; ++i)
      BWrInt32(&buf[4*i], a[i]);
}

void
BRdFloat32Array (unsigned char *buf,
		 float *a,
		 int n)
{
  int i, len;
  unsigned char *dest, temp;

  if (sizeof(float) == 4)
    {
      /* copy directly into the destination */
      memcpy(a, buf, 4*n);
      if (bio_big_endian_machine ^ bio_big_endian_input)
	{
	  /* reverse byte order */
	  dest = (unsigned char *) a;
	  len = 4*n;
	  for (i = 0; i < len; i+=4)
	    {
	      temp = dest[i];
	      dest[i] = dest[i+3];
	      dest[i+3] = temp;
	      temp = dest[i+2];
	      dest[i+2] = dest[i+1];
	      dest[i+1] = temp;
	    }
	}
    }
  else
    for (i = 0; i < n; ++i)
      a[i] = BRdFloat32(&buf[4*i]);
}

void
BWrFloat32Array (unsigned char *buf,
		 float *a,
		 int n)
{
  int i, len;
  unsigned char temp;

  if (sizeof(int) == 4)
    {
      /* copy directly into the destination */
      memcpy(buf, a, 4*n);
      if (bio_big_endian_machine ^ bio_big_endian_output)
	{
	  /* reverse byte order */
	  len = 4*n;
	  for (i = 0; i < len; i += 4)
	    {
	      temp = buf[i];
	      buf[i] = buf[i+3];
	      buf[i+3] = temp;
	      temp = buf[i+2];
	      buf[i+2] = buf[i+1];
	      buf[i+1] = temp;
	    }
	}
    }
  else
    for (i = 0; i < n; ++i)
      BWrFloat32(&buf[4*i], a[i]);
}

void
BRdFloat64Array (unsigned char *buf,
		 double *a,
		 int n)
{
  int i, len;
  unsigned char *dest, temp;

  if (sizeof(double) == 8)
    {
      /* copy directly into the destination */
      memcpy(a, buf, 8*n);
      if (bio_big_endian_machine ^ bio_big_endian_input)
	{
	  /* reverse byte order */
	  dest = (unsigned char *) a;
	  len = 8*n;
	  for (i = 0; i < len; i+=8)
	    {
	      temp = dest[i];
	      dest[i] = dest[i+7];
	      dest[i+7] = temp;
	      temp = dest[i+1];
	      dest[i+1] = dest[i+6];
	      dest[i+6] = temp;
	      temp = dest[i+2];
	      dest[i+2] = dest[i+5];
	      dest[i+5] = temp;
	      temp = dest[i+3];
	      dest[i+3] = dest[i+4];
	      dest[i+4] = temp;
	    }
	}
    }
  else
    for (i = 0; i < n; ++i)
      a[i] = BRdFloat64(&buf[8*i]);
}

void
BWrFloat64Array (unsigned char *buf,
		 double *a,
		 int n)
{
  int i, len;
  unsigned char temp;

  if (sizeof(double) == 8)
    {
      /* copy directly into the destination */
      memcpy(buf, a, 8*n);
      if (bio_big_endian_machine ^ bio_big_endian_output)
	{
	  /* reverse byte order */
	  len = 8*n;
	  for (i = 0; i < len; i += 8)
	    {
	      temp = buf[i];
	      buf[i] = buf[i+7];
	      buf[i+7] = temp;
	      temp = buf[i+1];
	      buf[i+1] = buf[i+6];
	      buf[i+6] = temp;
	      temp = buf[i+2];
	      buf[i+2] = buf[i+5];
	      buf[i+5] = temp;
	      temp = buf[i+3];
	      buf[i+3] = buf[i+4];
	      buf[i+4] = temp;
	    }
	}
    }
  else
    for (i = 0; i < n; ++i)
      BWrFloat64(&buf[8*i], a[i]);
}

int
FRdUInt8 (FILE *stream)
{
  int c;

  if ((c = fgetc(stream)) == EOF)
    {
      bio_error = 1;
      return(0);
    }
  return(c & 0xff);
}

void
FWrUInt8 (FILE *stream,
	  int v)
{
  if (fputc(v & 0xff, stream) == EOF)
    bio_error = 1;
}

int
FRdInt8 (FILE *stream)
{
  int c;

  if ((c = fgetc(stream)) == EOF)
    {
      bio_error = 1;
      return(0);
    }
  return((int) ((signed char) (c & 0xff)));
}

void
FWrInt8 (FILE *stream,
	 int v)
{
  if (fputc(v & 0xff, stream) == EOF)
    bio_error = 1;
}

int
FRdInt16 (FILE *stream)
{
  unsigned char buf[2];

  if (fread(buf, 2, 1, stream) != 1)
    {
      bio_error = 1;
      return(0);
    }
  return(BRdInt16(buf));
}

void
FWrInt16 (FILE *stream,
	  int v)
{
  unsigned char buf[2];

  BWrInt16(buf, v);
  if (fwrite(buf, 2, 1, stream) != 1)
    bio_error = 1;
}

int
FRdInt32 (FILE *stream)
{
  unsigned char buf[4];

  if (fread(buf, 4, 1, stream) != 1)
    {
      bio_error = 1;
      return(0);
    }
  return(BRdInt32(buf));
}

void
FWrInt32 (FILE *stream,
	  int v)
{
  unsigned char buf[4];

  BWrInt32(buf, v);
  if (fwrite(buf, 4, 1, stream) != 1)
    bio_error = 1;
}

float
FRdFloat32 (FILE *stream)
{
  unsigned char buf[4];

  if (fread(buf, 4, 1, stream) != 1)
    {
      bio_error = 1;
      return(0.0);
    }
  return(BRdFloat32(buf));
}

void
FWrFloat32 (FILE *stream,
	    float v)
{
  unsigned char buf[4];

  BWrFloat32(buf, v);
  if (fwrite(buf, 4, 1, stream) != 1)
    bio_error = 1;
}

double
FRdFloat64 (FILE *stream)
{
  unsigned char buf[8];

  if (fread(buf, 1, 8, stream) != 8)
    {
      bio_error = 1;
      return(0);
    }
  return(BRdFloat64(buf));
}

void
FWrFloat64 (FILE *stream,
	    double v)
{
  unsigned char buf[8];

  BWrFloat64(buf, v);
  if (fwrite(buf, 8, 1, stream) != 1)
    bio_error = 1;
}

void
FRdUInt8Array (FILE *stream,
	       unsigned char *a,
	       int n)
{
  if (fread(a, 1, n, stream) != n)
    bio_error = 1;
}

void
FWrUInt8Array (FILE *stream,
	       unsigned char *a,
	       int n)
{
  if (fwrite(a, 1, n, stream) != n)
    bio_error = 1;
}

void
FRdInt8Array (FILE *stream,
	      signed char *a,
	      int n)
{
  if (fread(a, 1, n, stream) != n)
    bio_error = 1;
}

void
FWrInt8Array (FILE *stream,
	      signed char *a,
	      int n)
{
  if (fwrite(a, 1, n, stream) != n)
    bio_error = 1;
}

void
FRdInt16Array (FILE *stream,
	       short *a,
	       int n)
{
  int i, len;
  unsigned char *buf, temp;

  if (sizeof(short) == 2)
    {
      buf = (unsigned char *) a;
      if (fread(buf, 2, n, stream) != n)
	bio_error = 1;
      if (bio_big_endian_machine ^ bio_big_endian_input)
	{
	  /* swap bytes */
	  len = 2*n;
	  for (i = 0; i < len; i+=2)
	    {
	      temp = buf[i];
	      buf[i] = buf[i+1];
	      buf[i+1] = temp;
	    }
	}
    }
  else
    {
      buf = (unsigned char *) malloc(2*n);
      if (fread(buf, 2, n, stream) != n)
	bio_error = 1;
      BRdInt16Array(buf, a, n);
      free(buf);
    }
}

void
FWrInt16Array (FILE *stream,
	       short *a,
	       int n)
{
  int i, len;
  unsigned char *buf;

  if (sizeof(short) == 2 &&
      !(bio_big_endian_machine ^ bio_big_endian_output))
    {
      /* directly write out to the file */
      if (fwrite(a, 2, n, stream) != n)
	bio_error = 1;
    }
  else
    {
      buf = (unsigned char *) malloc(2*n);
      BWrInt16Array(buf, a, n);
      if (fwrite(buf, 2, n, stream) != n)
	bio_error = 1;
      free(buf);
    }
}

void
FRdInt32Array (FILE *stream,
	       int *a,
	       int n)
{
  int i, len;
  unsigned char *buf, temp;

  if (sizeof(int) == 4)
    {
      buf = (unsigned char *) a;
      if (fread(buf, 4, n, stream) != n)
	bio_error = 1;
      if (bio_big_endian_machine ^ bio_big_endian_input)
	{
	  /* reverse byte order */
	  len = 4*n;
	  for (i = 0; i < len; i+=4)
	    {
	      temp = buf[i];
	      buf[i] = buf[i+3];
	      buf[i+3] = temp;
	      temp = buf[i+1];
	      buf[i+1] = buf[i+2];
	      buf[i+2] = temp;
	    }
	}
    }
  else
    {
      buf = (unsigned char *) malloc(4*n);
      if (fread(buf, 4, n, stream) != n)
	bio_error = 1;
      BRdInt32Array(buf, a, n);
      free(buf);
    }
}

void
FWrInt32Array (FILE *stream,
	       int *a,
	       int n)
{
  unsigned char *buf;

  if (sizeof(int) == 4 &&
      !(bio_big_endian_machine ^ bio_big_endian_output))
    {
      /* directly write out to the file */
      if (fwrite(a, 4, n, stream) != n)
	bio_error = 1;
    }
  else
    {
      buf = (unsigned char *) malloc(4*n);
      BWrInt32Array(buf, a, n);
      if (fwrite(buf, 4, n, stream) != n)
	bio_error = 1;
      free(buf);
    }
}

void
FRdFloat32Array (FILE *stream,
		 float *a,
		 int n)
{
  int i, len;
  unsigned char *buf, temp;

  if (sizeof(float) == 4)
    {
      buf = (unsigned char *) a;
      if (fread(buf, 4, n, stream) != n)
	bio_error = 1;
      if (bio_big_endian_machine ^ bio_big_endian_input)
	{
	  /* reverse byte order */
	  len = 4*n;
	  for (i = 0; i < len; i+=4)
	    {
	      temp = buf[i];
	      buf[i] = buf[i+3];
	      buf[i+3] = temp;
	      temp = buf[i+1];
	      buf[i+1] = buf[i+2];
	      buf[i+2] = temp;
	    }
	}
    }
  else
    {
      buf = (unsigned char *) malloc(4*n);
      if (fread(buf, 4, n, stream) != n)
	bio_error = 1;
      BRdFloat32Array(buf, a, n);
      free(buf);
    }
}

void
FWrFloat32Array (FILE *stream,
		float *a,
		int n)
{
  unsigned char *buf;

  if (sizeof(float) == 4 &&
      !(bio_big_endian_machine ^ bio_big_endian_output))
    {
      /* directly write out to the file */
      if (fwrite(a, 4, n, stream) != n)
	bio_error = 1;
    }
  else
    {
      buf = (unsigned char *) malloc(4*n);
      BWrFloat32Array(buf, a, n);
      if (fwrite(buf, 4, n, stream) != n)
	bio_error = 1;
      free(buf);
    }
}

void
FRdFloat64Array (FILE *stream,
		 double *a,
		 int n)
{
  int i, len;
  unsigned char *buf, temp;

  if (sizeof(double) == 8)
    {
      buf = (unsigned char *) a;
      if (fread(buf, 8, n, stream) != n)
	bio_error = 1;
      if (bio_big_endian_machine ^ bio_big_endian_input)
	{
	  /* reverse byte order */
	  len = 8*n;
	  for (i = 0; i < len; i+=8)
	    {
	      temp = buf[i];
	      buf[i] = buf[i+7];
	      buf[i+7] = temp;
	      temp = buf[i+1];
	      buf[i+1] = buf[i+6];
	      buf[i+6] = temp;
	      temp = buf[i+2];
	      buf[i+2] = buf[i+5];
	      buf[i+5] = temp;
	      temp = buf[i+3];
	      buf[i+3] = buf[i+4];
	      buf[i+4] = temp;
	    }
	}
    }
  else
    {
      buf = (unsigned char *) malloc(8*n);
      if (fread(buf, 8, n, stream) != n)
	bio_error = 1;
      BRdFloat64Array(buf, a, n);
      free(buf);
    }
}

void
FWrFloat64Array (FILE *stream,
		 double *a,
		 int n)
{
  unsigned char *buf;

  if (sizeof(double) == 8 &&
      !(bio_big_endian_machine ^ bio_big_endian_output))
    {
      /* directly write out to the file */
      if (fwrite(a, 8, n, stream) != n)
	bio_error = 1;
    }
  else
    {
      buf = (unsigned char *) malloc(8*n);
      BWrFloat64Array(buf, a, n);
      if (fwrite(buf, 8, n, stream) != n)
	bio_error = 1;
      free(buf);
    }
}
