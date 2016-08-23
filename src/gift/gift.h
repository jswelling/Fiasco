/*
 *	gift.h - common data structures & functions for gift
 *
 *	Copyright (c) 1996  Pittsburgh Supercomputing Center
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
 *		1/96	Written by Greg Hood (PSC)
 */

#include "misc.h"

/* The basic data types are: */
typedef short GiftDataType;                   /* the lowest 8 bits encode the length in bits of the type */
#define GIFT_UINT8		(0x0100 | 8)	/* 1 byte, unsigned integer format */
#define GIFT_INT16		(0x0100 | 16)	/* 2 bytes, big-endian signed integer format */
#define GIFT_INT32		(0x0100 | 32)	/* 4 bytes, big-endian signed integer format */
#define GIFT_FLOAT32		(0x0200 | 32)	/* 4 bytes, big-endian IEEE single-precision format */
#define GIFT_FLOAT64		(0x0200 | 64)	/* 8 bytes, big-endian IEEE double-precision format */

/* The space occuped by each data type: */
#define GiftBitsPerItem(t)	((t) & 0xff)
#define GiftBytesPerItem(t)	(((t) & 0xff) >> 3)

/* Several classes of data are defined: */
typedef short GiftClass;
#define GIFT_I_SPACE		0	/* image-space */
#define GIFT_JX_SPACE		1	/* image-space along X, Fourier-space along Y */
#define GIFT_JY_SPACE		2	/* Fourier-space along X, image-space along Y */
#define GIFT_K_SPACE		3	/* Fourier space */
#define GIFT_P_SPACE		4	/* spiral projection space */

/* Dimensions are identified by the following constants: */
typedef int GiftDimensionType;
#define GIFT_V_DIMENSION	0	/* the vector dimension allowing multiple data items per voxel */
#define GIFT_X_DIMENSION	1	/* the x dimension */
#define GIFT_Y_DIMENSION	2	/* the y dimension */
#define GIFT_Z_DIMENSION	3	/* the z dimension indicating which slice */
#define GIFT_T_DIMENSION	4	/* the time dimension indicating when the image was acquired */

/* Each dimension is described by the following structure: */
typedef struct GiftDimension {
  GiftDimensionType type;	/* what this dimension represents */
  int min;			/* the starting index for this dimension */
  int max;			/* the ending index for this dimension */
  int stride;			/* the interval between successive indices */
  int n;			/* the number of indices along this dimension;
				   must equal (max-min)/stride + 1 */
  float size;			/* the physical length (in mm or seconds) corresponding
				   to one integral step along this dimension */
} GiftDimension;

#define GIFT_MAX_DIMENSIONS		10		/* space is reserved for this many dimensions although
							   we currently support only 5 */
/* GiftHeader is the intermediate header */
typedef struct GiftHeader {
  GiftClass class;		/* the type of image contained in the data file */
  GiftDataType data_type;	/* the type of each individual data item */
  short n_dims;			/* # of dimensions (currently must be 5) */
  short n_image_dims;		/* the number of dimensions within a single image */
  int n_images;			/* the number of images in the data file */
  int n_items_per_image;	/* the number of data items per image */
  GiftDimension dim[GIFT_MAX_DIMENSIONS];	/* the dimensions within the data file */
  char *corrupt;		/* an array indicating which images are considered
                                   corrupt and to be ignored during processing
				   (0 signifies usable, 1 signifies corrupt */
} GiftHeader;

typedef struct FileListElement {
  struct FileListElement *next;
  Filename name;
  int size;
} FileListElement, *FileList;

extern FILE *input;
extern FILE *output;

/* the intermediate format */
extern GiftHeader fh;		/* header */
extern char *image;		/* one image */

extern int slice_start;
extern int slice_end;
extern int slice_stride;
extern int time_start;
extern int time_end;
extern int time_stride;

extern Boolean delete_input;
extern Boolean compress_output;

extern Boolean big_endian_input;
extern Boolean big_endian_output;

void AppendToFileList (FileList *pHead, FileList *pTail, char *name, int n);
Boolean HasExtension (char *name, char *ext);
void ReadImage (int time, int slice);
void WriteImage (int time, int slice);
void Compress (Filename name, int size);
void Uncompress (Filename unc_name, char *name);
int UncompressBatch (Filename unc_name, FileList *fl);
void Delete (Filename name);
