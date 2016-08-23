/*
 *	C Library for Reading and Writing the
 *	Pittsburgh MRI Format
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
 *	History:
 *		3/96: Written by Greg Hood
 */

/* the following constants may modified if necessary,
   but the defaults should work for most sites */
#define MRI_MAX_DIMS		16	/* maximum # of dimensions in
					   a chunk */
#define MRI_MAX_KEY_LENGTH	255	/* the maximum # of characters
					   in a key name */
#define MRI_MAX_VALUE_LENGTH	4095	/* the maximum # of characters
					   in a key's value */
#define MRI_MAX_CHUNKS		256	/* the maximum # of chunks
					   in any one dataset */
#define MRI_MAX_FILENAME_LENGTH	255	/* the maximum # of characters
					   in a filename */
#define MRI_MAX_OPEN_FILES	8	/* the maximum # of files the
					   libmri library is allowed to
					   have open for each dataset
					   at any one time */
#define MRI_SAFE_BUFFER_COUNT	4	/* the number of unretained buffers
					   that a client may read before
					   they starting getting recycled */
#define MRI_MAX_BUFFER_COUNT	8	/* the maximum number of
					   unretained buffers that the libmri
					   library may keep around */
#define MRI_ALIGNMENT_THRESHOLD	65536	/* chunks larger than this number of
					   bytes will be aligned on disk */
#define MRI_ALIGNMENT_BOUNDARY	16384	/* the boundary increment in bytes
					   to which large chunks will be
					   aligned */


/*----------- nothing beyond this point----------------*/
/*----------- should have to be modified---------------*/
/*---------------for installation----------------------*/

#include <stdio.h>

/* the following refer to the data type stored in the files */
typedef int MRI_Datatype;
#define MRI_UINT8	0	/* 8 bit unsigned integer */
#define MRI_INT16	1	/* 16 bit integer */
#define MRI_INT32	2	/* 32 bit integer */
#define MRI_FLOAT32	3	/* 32 bit IEEE floating point */
#define MRI_FLOAT64	4	/* 64 bit IEEE floating point */


/* the following refer to the data type in memory; they correspond
   to the standard C types available for the architecture on which
   this library is compiled */
typedef int MRI_ArrayType;
#define MRI_RAW			0	/* raw bytes, no type conversions */
#define MRI_UNSIGNED_CHAR	1
#define MRI_SHORT		2
#define MRI_INT			3
#define MRI_FLOAT		4
#define MRI_DOUBLE		5

#define MRI_SCALAR		0x00
#define MRI_COMPLEX		0x10
#define MRI_VECTOR		0x20

#define MRI_COMPLEX_SHORT	(MRI_COMPLEX | MRI_SHORT)
#define MRI_COMPLEX_INT		(MRI_COMPLEX | MRI_INT)
#define MRI_COMPLEX_FLOAT	(MRI_COMPLEX | MRI_FLOAT)
#define MRI_COMPLEX_DOUBLE	(MRI_COMPLEX | MRI_DOUBLE)

#define MRI_VECTOR_SHORT	(MRI_VECTOR | MRI_SHORT)
#define MRI_VECTOR_INT		(MRI_VECTOR | MRI_INT)
#define MRI_VECTOR_FLOAT	(MRI_VECTOR | MRI_FLOAT)
#define MRI_VECTOR_DOUBLE	(MRI_VECTOR | MRI_DOUBLE)


typedef int MRI_OpenMode;
#define MRI_READ		0	/* open for reading only */
#define MRI_WRITE		1	/* create a new dataset */
#define MRI_MODIFY		2	/* open for reading and writing */
#define MRI_MODIFY_DATA		3	/* open for reading and writing
					   chunk data -- no other
					   modifications allowed */

typedef int MRI_Error_Handling;
#define MRI_ABORT_ON_ERROR	0	/* print out a message to
					   stderr and then call abort()
					   if any error is detected */
#define MRI_REPORT_ERRORS	1	/* print out a message to
					   stderr, but try to recover */
#define MRI_IGNORE_ERRORS	2	/* don't print anything or abort;
					   just set mri_error to
					   indicate the problem */

/* these constants are used to designate where a chunk will be
   placed in a data file */
typedef int MRI_Order;
#define MRI_EXTERNAL		-2	/* this chunk is located in an
					   external file and cannot be moved */
#define MRI_FIXED_OFFSET	-1	/* place it at the byte offset
					   specified in the ".offset" key */
#define MRI_FIRST		0	/* place it first in the file */
#define MRI_ANYWHERE		100	/* this is a somewhat arbitrary value
					   chosen to be between MRI_FIRST
					   and MRI_LAST */
#define MRI_LAST		10000	/* place it at the end of the file */


#define MRI_UNSPECIFIED		(-2147483647)	/* value returned if
						   key is non-existent */



/*--------- the internals of these structures ----------*/
/*----------- should not concern the user of -----------*/
/*---------------- the libmri library ------------------*/

typedef struct MRI_KeyValue {
  struct MRI_KeyValue *next_in_hash_table;	/* the next key in this
						   hash bucket */
  char *key;			/* the key name */
  char *value;			/* the key's value */
} MRI_KeyValue;

typedef struct MRI_Buffer {
  struct MRI_Buffer *next;	/* the next buffer on the buffers
				   or retained_buffers list */
  int size;			/* the size of the buffer in bytes */
  void *buffer;			/* the buffer area itself */
} MRI_Buffer;

typedef struct MRI_Dataset {
  /* general info */
  char *name;			/* the filename where the dataset's
				   header is stored */
  MRI_OpenMode mode;		/* the mode (read/write/modify) in
				   which this dataset was opened */

  /* file info */
  int n_open_files;		/* the number of files we have open
				   for this dataset */
  struct MRI_File *files;	/* the list of files that we are
				   referencing */

  /* header info */
  int header_size;		/* the number of bytes reserved for
				   the header */
  struct MRI_File *header_file;	/* the file that contains the header */
  int n_keys;			/* the number of keys */
  int hash_table_size;		/* the number of buckets in the hash
				   table */
  MRI_KeyValue **hash_table;	/* a hash table so that we can go quickly
				   from a key's name to it's value */

  /* iteration support */
  struct MRI_KeyValue **iteration_table; /* an alphabetically sorted table
					    of the keys in the dataset */
  int iteration_count;		/* the number of entries in the
				   iteration_table */
  int iteration_index;		/* what index we are on in the
				   iteration_table */

  /* chunk info */
  struct MRI_Chunk *chunks;	/* the chunks associated with this dataset */
  int recompute_positions;	/* if TRUE, we should recompute all chunk
				   positions before doing any reading or
				   writing to any chunk */

  /* user buffering */
  struct MRI_Buffer *buffers;	/* the libmri-managed buffers for reading
				   and writing portions of chunks */
  struct MRI_Buffer *retained_buffers;	/* the buffers that the application
					   program has retained indefinitely
					   for its own use */
  
  /* higher-level image support */
  int std_images;		/* if TRUE, this dataset contains
				   standard format images */
  int std_image_vector_size;	/* extent of the v dimension */
  int std_image_size;		/* the number of voxels in one image */
  int std_n_slices;		/* the number of slices for one time */
} MRI_Dataset;

typedef struct MRI_File {
  struct MRI_File *next;	/* next on file list */
  MRI_Dataset *ds;		/* the dataset this file belongs to */
  char *name;			/* the filename */
  FILE *fp;			/* this will be be non-NULL if we
				   have the file open */
  int writeable;		/* if TRUE, the file pointer (fp)
				   was opened for writing */
  int last_use;			/* we set this field when we read or
				   write a file, so that we can tell
				   which files haven't been used
				   recently, and can be closed if
				   necessary */
  int used;			/* set to TRUE during garbage collection
				   when we determine that this
				   file is still being used */
  int external;			/* if TRUE, this file is strictly
				   not part of the dataset (e.g.
				   a temporary file) */
} MRI_File;

typedef struct MRI_Chunk {
  struct MRI_Chunk *next;	/* next on chunk list */
  MRI_Dataset *ds;		/* the dataset this chunk belongs to */
  char *name;			/* the key name of this chunk */

  /* these fields tell what the chunk
     is supposed to look like */
  MRI_File *file;		/* the file that this chunk should be
				   stored in */
  MRI_Datatype datatype;	/* the datatype of the items in the
				   chunk */
  char *dimensions;		/* the characters representing the
				   dimensions of the data in the chunk */
  int extent[MRI_MAX_DIMS];	/* the size of each dimension */
  int little_endian;		/* if TRUE, the chunk is stored in
				   little-endian format */
  MRI_Order order;		/* specifies where the chunk
				   should go in its file */
  int offset;			/* specifies the absolute byte offset
				   where the chunk should go in its file */
  int size;			/* the size of the chunk in bytes */

  /* the actual_* fields record the state of the chunk
     on disk; these may lag behind the values above,
     in which case the modified field will be TRUE */
  MRI_File *actual_file;
  MRI_Datatype actual_datatype;
  char *actual_dimensions;
  int actual_extent[MRI_MAX_DIMS];
  int actual_little_endian;
  int actual_offset;
  int actual_size;

  int modified;		/* TRUE if the chunk's location
			   or attributes have changed */
  int checked;		/* TRUE if we have checked this
			   chunk for possible repositioning */

  int repositioning;	/* TRUE while we are trying to
			   reposition this chunk */

  int mmapped;		/* TRUE if we have mmap'ped the
			   chunk data into memory */
  char *addr;		/* pointer to the mmap'ped data */

  int ready_to_read;	/* if FALSE, PrepareToRead must be
			   called before reading */
  int ready_to_write;	/* if FALSE, PrepareToWrite must be
			   called before writing */
} MRI_Chunk;


/*-----------------------------------------------------------------------
 *	APPLICATION PROGRAM INTERFACE FUNCTIONS
 *-----------------------------------------------------------------------*/

/*-------- OPENING AND CLOSING DATASETS ---------------------------*/

extern MRI_Dataset *mri_open_dataset (const char *filename, MRI_OpenMode mode);
extern MRI_Dataset *mri_copy_dataset (const char *filename, MRI_Dataset *original);
extern void mri_close_dataset (MRI_Dataset *ds);
extern void mri_destroy_dataset (MRI_Dataset *ds);

/*-------- GETTING KEY VALUES ------------------------------------*/
extern int mri_get_int (MRI_Dataset *ds, const char *key);
extern double mri_get_float (MRI_Dataset *ds, const char *key);
extern char *mri_get_string (MRI_Dataset *ds, const char *key);

/*-------- SETTING KEY VALUES ------------------------------------*/
extern void mri_set_int (MRI_Dataset *ds, const char *key, const int value);
extern void mri_set_float (MRI_Dataset *ds, const char *key, const float value);
extern void mri_set_string (MRI_Dataset *ds, const char *key, const char *value);

/*-------- MISCELLANEOUS KEY FUNCTIONS ---------------------------*/
extern int mri_has (MRI_Dataset *ds, const char *key);
extern void mri_remove (MRI_Dataset *ds, const char *key);
extern void mri_iterate_over_keys (MRI_Dataset *ds);
extern char *mri_next_key (MRI_Dataset *ds);
extern char *mri_cat (const char *s1, const char *s2);

/*-------- CHUNKS ------------------------------------------------*/
extern void mri_create_chunk (MRI_Dataset *ds, const char *key);
extern void mri_update_chunk (MRI_Dataset *ds, const char *key);
extern void *mri_get_chunk (MRI_Dataset *ds, const char *key,
			    int size, int offset, MRI_ArrayType type);
extern void mri_set_chunk (MRI_Dataset *ds, const char *key,
			   int size, int offset, MRI_ArrayType type, void *ptr);

/*--------- BUFFER MANAGEMENT ------------------------------------*/
extern void mri_retain_buffer (MRI_Dataset *ds, void *ptr);
extern void mri_discard_buffer (MRI_Dataset *ds, void *ptr);

/*--------- ERROR HANDLING ---------------------------------------*/
extern char *mri_error;
extern void mri_set_error_handling (MRI_Error_Handling mode);

/*--------- IMAGES -----------------------------------------------*/
extern void *mri_get_image (MRI_Dataset *ds,
			    const int time, const int slice,
			    MRI_ArrayType type);
extern void mri_set_image (MRI_Dataset *ds,
			   const int time, const int slice,
			   MRI_ArrayType type, void *image);
