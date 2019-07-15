/*
 *	C Library for Reading and Writing the
 *	Pittsburgh MRI Format
 *
 *	Copyright (c) 1996,1997,1998,1999 Pittsburgh Supercomputing Center
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
 *		9/99: Revisions to accommodate datasets larger
 *			than 2GB
 */

#if defined(SUN4SOL2) || defined(LINUX) || defined(HPPA) || defined(HPPA20) || defined(ALPHA) || defined(DARWIN) || defined(CYGWIN)
#define NO_FSEEK64
#if !defined(HPPA)
#define HAVE_FSEEKO
#define _GNU_SOURCE
#define _FILE_OFFSET_BITS     64
#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE
#endif
#endif

#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#if ( ! ( DARWIN || CYGWIN ) )
#include <values.h>
#endif
#include <errno.h>
#ifdef AFS
#include <sys/types.h>
#include <sys/wait.h>
#include <errno.h>
#endif
#include "mri.h"
#include "bio.h"

#ifdef DARWIN
#define finite( foo ) isfinite( foo )
#endif

#ifndef LONGLONG_MIN
#define LONGLONG_MIN (-9223372036854775807LL-1LL)
#endif
#ifndef LONGLONG_MAX
#define LONGLONG_MAX 9223372036854775807LL
#endif

#define TRUE		1
#define FALSE		0

#define BUFFER_SIZE	1048575
#define MAX_CHUNK_SIZE  LONGLONG_MAX

typedef struct EmptyBlock {
  long long start;
  long long end;
} EmptyBlock;

typedef struct CopyRequest {
  struct CopyRequest *next;
  struct MRI_File *src_file;
  long long src_offset;
  struct MRI_File *dest_file;
  long long dest_offset;
  long long size;
  MRI_Chunk *chunk;
} CopyRequest;

static char rcsid[] = "$Id: libmri.c,v 1.46 2007/04/26 23:17:23 welling Exp $";

char *mri_error = NULL;

static MRI_Error_Handling error_handling = MRI_ABORT_ON_ERROR;
static unsigned int tmp_num = 0;
static unsigned int file_access = 1;

#ifdef DEBUG
static FILE *logfile = NULL;	/* for debugging */
#endif

/* FORWARD DECLARATIONS */
static MRI_File *GetFile (MRI_Dataset *ds, char *name);
static void CreateNewDataset (MRI_Dataset *ds);
static MRI_Chunk *NewChunk (MRI_Dataset *ds, char *name);
static MRI_File *CreateTempFile (MRI_Dataset *ds);
static void DeallocateDataset (MRI_Dataset *ds);
static int CheckForHooks (MRI_Dataset *ds, MRI_KeyValue *kv);
static char *GetBuffer (MRI_Dataset *ds, long long size);
static void CopyBlock (MRI_File *dest_file, long long dest_offset,
		       MRI_File *src_file, long long src_offset,
		       long long size);
static void ClearBlock (MRI_File *dest_file, long long dest_offset,
			long long size);
static MRI_KeyValue *FindInHashTable (MRI_Dataset *ds, const char *key,
				      int add);
static void RemoveFromHashTable (MRI_Dataset *ds, const char *key);
static int ReadKeyValuePair (MRI_Dataset *ds, FILE *f);
static void ModifyChunk (MRI_Chunk *ch);
static int ReadHeader (MRI_Dataset *ds, FILE *f);
static int ReadString (MRI_Dataset *ds, char *s, int max_len, FILE *f);
static int ReadQuotedString (MRI_Dataset *ds, char *s, int max_len, FILE *f);
static int ReadBackslashedChar (FILE *f);
static int HashFunction (const char *s);
static void WriteHeader (MRI_Dataset *ds, int separator, FILE *f);
static void WriteKeyValuePair (MRI_KeyValue *kv, FILE *f);
static void WriteString (char *s, FILE *f);
static void WriteQuotedString (char *s, FILE *f);
static void ComputeChunkPositions (MRI_Dataset *ds);
static long long ReserveBlock (EmptyBlock *empty_blocks, int *n_empty_blocks,
			       long long offset, long long size);
static int OpenFile (MRI_File *file, int for_writing);
static void CloseFile (MRI_File *file);
static void DestroyFile (MRI_File *file);
static void RepositionChunk (MRI_Chunk *ch, CopyRequest **copy_queue);
static void UpdateChunkAttributes (MRI_Chunk *ch);
static void ConvertChunk (MRI_Chunk *ch, MRI_File *f, long long offset);
static int MRI_TypeLength (MRI_Datatype datatype);
static int MRI_ArrayTypeLength( MRI_ArrayType type, int n );
static void CheckForStdImages (MRI_Dataset *ds);
static int InitializeKeyValue (MRI_Dataset *ds, const char *key,
			       const char *value);
static void SetChunkNotReady (MRI_Chunk *ch);
static int PrepareToRead (MRI_Chunk *ch);
static int PrepareToWrite (MRI_Chunk *ch);
static void CleanFiles (MRI_Dataset *ds);
static int CompareKeys (const void *k1, const void *k2);
static void ChunkSeek (MRI_Chunk *ch, long long offset);
static void DeallocateChunk (MRI_Chunk *chunk);
static int CheckForRemovalHooks (MRI_Dataset *ds, MRI_KeyValue *kv);
static int ConvertDatatype (MRI_Datatype *pdt, char *s);
static void mri_report_error (MRI_Dataset *ds, char *fmt, ...);
static void mri_report_warning (MRI_Dataset *ds, char *fmt, ...);
static void CopyConvFloatLonglong(float* float_buf,
				  long long* longlong_buf,
				  int size,int* error);
static void CopyConvDoubleLonglong(double* float_buf,
				   long long* longlong_buf,
				   int size,int* error);
#ifdef AFS
static int MySystem( const char* command );
static int CheckForAFS( const char* fname );
static void FlushAFS( MRI_Dataset *ds );
#endif

#ifdef DEBUG
static void Log(char *fmt, ...);
#endif
#ifdef NO_FSEEK64
static int mri_fseek (FILE *f, long long offset, int whence);
static long long mri_ftell (FILE *f);
static int mri_ftruncate (int fildes, long long length);
#else
#define mri_fseek		fseek64
#define mri_ftell		ftell64
#define mri_ftruncate(a,b)	ftruncate64((a), (off64_t) (b))
#endif


/* PUBLIC FUNCTIONS */

/*-------- OPENING AND CLOSING DATASETS ---------------------------*/

MRI_Dataset
*mri_open_dataset (const char *name, MRI_OpenMode mode)
{
  MRI_Dataset *ds;
  MRI_KeyValue *kv;
  int i;
  MRI_Chunk *ch;
  long long first_start;
  char file_name[MRI_MAX_FILENAME_LENGTH+1];
  int empty;
  static int bio_initialized = FALSE;

#ifdef DEBUG
  if (logfile == NULL)
    logfile = fopen("libmri.log", "a");
  Log("Opening dataset %s in mode %d\n", name, mode);
#endif

  /* initialize the BIO (binary I/O) package on the first call to
     mri_open_dataset */
  if (!bio_initialized)
    {
      bio_error = 0;
      InitBIO();
      if (bio_error)
	{
	  mri_report_error(NULL, "mri_open_dataset: cannot determine endianness of this computer\n");
	  return(NULL);
	}
      bio_initialized = TRUE;
    }

  /* if there is no .mri extension on the filename,
     add one */
  strcpy(file_name, name);
  if ((i = strlen(file_name)) < 4 ||
      strcmp(&file_name[i-4], ".mri") != 0)
    strcat(file_name, ".mri");

  /* initialize the Dataset structure */
  if ((ds = (MRI_Dataset *) malloc(sizeof(MRI_Dataset))) == NULL ||
      (ds->name = (char *) malloc(strlen(file_name)+1)) == NULL)
    {
      if (ds != NULL)
	free(ds);
      mri_report_error(NULL, "mri_open_dataset: cannot allocate space for MRI_Dataset\n");
      return(NULL);
    }
  strcpy(ds->name, file_name);
  ds->mode = mode;
  ds->n_open_files = 0;
  ds->files = NULL;

  ds->header_size = 512;   /* initially assume the header will be this long */
  ds->header_file = GetFile(ds, "");
  ds->n_keys = 0;
  ds->hash_table_size = 64;
  ds->hash_table = (MRI_KeyValue **) calloc(ds->hash_table_size, sizeof(MRI_KeyValue *));

  ds->iteration_table = NULL;
  ds->iteration_count = 0;
  ds->iteration_index = 0;

  ds->chunks = NULL;
  ds->recompute_positions = FALSE;

  ds->buffers = NULL;
  ds->retained_buffers = NULL;

#ifdef AFS
  ds->some_parts_in_afs= CheckForAFS(ds->name);
  if (ds->some_parts_in_afs) FlushAFS(ds);
#endif

  if (!OpenFile(ds->header_file, ds->mode == MRI_WRITE || ds->mode == MRI_MODIFY))
    {
      mri_report_error(ds, "mri_open_dataset: could not open file %s for %s\n",
		ds->header_file->name,
		(ds->mode == MRI_WRITE || ds->mode == MRI_MODIFY) ? "writing" : "reading");
      DeallocateDataset(ds);
      return(NULL);
    }
  if (mri_fseek(ds->header_file->fp, 0LL, SEEK_END) != 0)
    {
      mri_report_error(ds,
		       "libmri: could not seek to end of file %s\n",
		       ds->header_file->name);
      DeallocateDataset(ds);
      return(NULL);
    }
  empty = mri_ftell(ds->header_file->fp) == 0LL;
  rewind(ds->header_file->fp);
  if (ds->mode == MRI_WRITE || empty)
    CreateNewDataset(ds);
  else if (!ReadHeader(ds, ds->header_file->fp))
    {
      DeallocateDataset(ds);
      return(NULL);
    }
  
  /* If this is a read-only file, we better have found some 
   * key-value pairs!
   */
  if (mode==MRI_READ && ds->n_keys==0) 
    {
      mri_report_error(ds, "mri_open_dataset: %s is unexpectedly empty\n",
		       ds->header_file->name);
      DeallocateDataset(ds);
      return(NULL);
    }


  /* build the associated chunk data structures */
  for (i = 0; i < ds->hash_table_size; ++i)
    for (kv = ds->hash_table[i]; kv != NULL; kv = kv->next_in_hash_table)
      if (strcmp(kv->value, "[chunk]") == 0)
	if (NewChunk(ds, kv->key) == NULL)
	  {
	    mri_report_error(ds, "mri_open_dataset: error in chunk %s\n",
			     kv->key);
	    DeallocateDataset(ds);
	    return(NULL);
	  }

  /* check if there is actually more space reserved for the header */
  first_start = 999999999999999999LL;
  for (ch = ds->chunks; ch != NULL; ch = ch->next)
    if (ch->file == ds->header_file &&
	ch->offset < first_start)
      first_start = ch->offset;
  if (first_start < 999999999999999999LL)
    ds->header_size = first_start;

  return(ds);
}

MRI_Dataset *mri_copy_dataset (const char *filename, MRI_Dataset *original)
{
  MRI_Dataset *ds, *nds;
  int count;
  int len;
  long long total;
  long long n;
  void *ptr;
  MRI_Chunk *ch, *nch;
  MRI_File *f;
  char *key;
  char file_key[MRI_MAX_KEY_LENGTH+1];
  char value[MRI_MAX_VALUE_LENGTH+1];
  char chunk_filename[MRI_MAX_FILENAME_LENGTH+1];

  /* check original dataset */
  if (original == NULL)
    {
      mri_report_error(NULL, "mri_copy_dataset: original dataset in invalid\n");
      return(NULL);
    }
  ds = original;

  /* open new dataset */
  if ((nds = mri_open_dataset(filename, MRI_WRITE)) == NULL)
    {
      mri_report_error(NULL, "mri_copy_dataset: cannot open new dataset\n");
      return(NULL);
    }

  /* copy all the keys over */
  mri_iterate_over_keys(ds);
  while ((key = mri_next_key(ds)) != NULL)
    {
      strcpy(value, mri_get_string(ds, key));

      /* check if the key is a chunk's filename,
	 which may have to be changed to a new name
	 if it is absolute */
      len = strlen(key);
      if (len >= 5 && strcmp(&key[len-5], ".file") == 0 &&
	  value[0] != '.')
	{
	  /* make sure the key belongs to a chunk */
	  strcpy(file_key, key);
	  file_key[len-5] = '\0';
	  for (ch = ds->chunks; ch != NULL; ch = ch->next)
	    if (strcmp(ch->name, file_key) == 0)
	      break;
	  if (ch != NULL)
	    {
	      /* construct the new filename */
	      count = 0;
	      for (f = ds->files; f != NULL && f != ch->file; f = f->next)
		if (!f->external && f->name[0] != '.')
		  ++count;
	      if (count == 0)
		sprintf(value, ".dat");
	      else
		sprintf(value, ".%d.dat", count);
	    }
	}

      /* set the key's value in the new dataset */
      mri_set_string(nds, key, value);
    }

  /* now go through and copy all the chunks */
  for (ch = ds->chunks; ch != NULL; ch = ch->next)
    if (ch->order != MRI_EXTERNAL)
      if (ds->mode == MRI_READ)
	{
	  /* we make a lazy copy of the original chunk */
	  for (nch = nds->chunks; nch != NULL; nch = nch->next)
	    if (strcmp(nch->name, ch->name) == 0)
	      break;
	  if (nch == NULL)
	    {
	      mri_report_error(ds, "libmri: internal error -- chunk name not copied correctly\n");
	      abort();
	    }
	  /* make sure we have an absolute filename (one with a '/' in it)
	     so that the call to GetFile does not try to relativize it to
	     the new dataset name */
	  if (ch->file->name[0] == '/')
	    strcpy(chunk_filename, ch->file->name);
	  else
	    {
	      if (getcwd(chunk_filename, MRI_MAX_FILENAME_LENGTH) == NULL)
		{
		  mri_report_error(ds, "libmri: cannot get current directory name\n");
		  abort();
		}
	      strcat(chunk_filename, "/");
	      strcat(chunk_filename, ch->file->name);
	    }
	  nch->actual_file = GetFile(nds, chunk_filename);
	  nch->actual_file->external = TRUE;
	  nch->actual_datatype = ch->datatype;
	  nch->actual_little_endian = ch->little_endian;
	  nch->actual_offset = ch->offset;
	  nch->actual_size = ch->size;
	  nch->modified = TRUE;
	  nds->recompute_positions = TRUE;
	}
      else
	/* we copy the chunks immediately */
	for (total = 0LL; total < ch->size; total += BUFFER_SIZE)
	  {
	    n = ch->size - total;
	    if (n > BUFFER_SIZE)
	      n = BUFFER_SIZE;
	    ptr = mri_get_chunk(ds, ch->name, n, total, MRI_RAW);
	    mri_set_chunk(nds, ch->name, n, total, MRI_RAW, ptr);
	  }

  return(nds);
}

void
mri_close_dataset (MRI_Dataset *ds)
{
  long long pos;
  MRI_File *temp;
  MRI_Chunk *ch;
  int alone;

#ifdef DEBUG
  Log("Closing dataset %s\n", ds->name);
#endif

  /* if the dataset is read-only or just data-writable,
     we only have to throw away the stuff in memory */
  if (ds->mode == MRI_READ || ds->mode == MRI_MODIFY_DATA)
    {
      DeallocateDataset(ds);
      return;
    }

  /* update the chunk positions */
  if (ds->recompute_positions)
    ComputeChunkPositions(ds);

  /* determine if header is alone in .mri file */
  alone = TRUE;
  for (ch = ds->chunks; ch != NULL; ch = ch->next)
    if (ch->file == ds->header_file)
      {
	alone = FALSE;
	break;
      }

  /* write the header out to a temporary file to check its size */
  temp = CreateTempFile(ds);
  if (!OpenFile(temp, TRUE))
    {
      mri_report_error(ds, "mri_close_dataset: could not open temp file\n");
      return;
    }

  WriteHeader(ds, !alone, temp->fp);
  pos = mri_ftell(temp->fp);
  CloseFile(temp);

  /* will it fit into reserved space? */
  if (alone)
    ds->header_size = pos;
  else if (pos > ds->header_size)
    {
      /* no, we have to expand header space */
      /* we add in a little extra space (10*ds->n_keys)
	 to be sure that we still have enough space
	 in case the printed offsets grow a bit */
      ds->header_size = 1;
      while (ds->header_size < (pos + 10*ds->n_keys))
	ds->header_size *= 2;

      ComputeChunkPositions(ds);

      /* rewrite the header */
      if (!OpenFile(temp, TRUE))
	{
	  mri_report_error(ds, "mri_close_dataset: could not open temp file\n");
	  return;
	}
      WriteHeader(ds, !alone, temp->fp);
      pos = mri_ftell(temp->fp);
      CloseFile(temp);
    }

  /* make sure to write all the chunks out */
  for (ch = ds->chunks; ch != NULL; ch = ch->next)
    if (ch->modified)
      RepositionChunk(ch, NULL);

  /* copy header from temp file to proper location */
  CopyBlock(ds->header_file, 0LL, temp, 0LL, pos);
  /* reset header size to actual size so gap gets zeroed */
  ds->header_size = pos;

  /* throw away temp file */
  DestroyFile(temp);

  /* zero in all the gaps */
  CleanFiles(ds);

#ifdef AFS
  if (ds->some_parts_in_afs) FlushAFS(ds);
#endif

  /* throw away stuff in memory */
  DeallocateDataset(ds);
}

void
mri_destroy_dataset (MRI_Dataset *ds)
{
  MRI_File *f;

  if (ds->mode == MRI_READ || ds->mode == MRI_MODIFY_DATA)
    {
      mri_report_error(ds, "mri_destroy_dataset: cannot destroy read-only dataset!\n");
      return;
    }

  /* unlink files */
  for (f = ds->files; f != NULL; f = f->next)
    if (!f->external)
      {
	CloseFile(f);
#ifdef DEBUG
	Log("GOING TO UNLINK FILE!!! (within mri_destroy_dataset): %s\n",
	    f->name);
#endif
	(void) unlink(f->name);
      }

  /* deallocate */
  DeallocateDataset(ds);
}


/*-------- GETTING KEY VALUES ------------------------------------*/

long long
mri_get_int (MRI_Dataset *ds, const char *key)
{
  long long v;

  if (!mri_has(ds, key))
    {
      mri_report_error(ds, "mri_get_int: non-existent key %s\n",key);
      return(MRI_UNSPECIFIED);
    }
  if (sscanf(mri_get_string(ds, key), "%lld", &v) != 1)
    {
      mri_report_error(ds, "mri_get_int: key value is not an integer\n");
      return(MRI_UNSPECIFIED);
    }
  return(v);
}

double
mri_get_float (MRI_Dataset *ds, const char *key)
{
  double v;

  if (!mri_has(ds, key))
    {
      mri_report_error(ds, "mri_get_float: non-existent key %s\n",key);
      return((double) MRI_UNSPECIFIED);
    }
  if (sscanf(mri_get_string(ds, key), "%lg", &v) != 1)
    {
      mri_report_error(ds, "mri_get_float: key value is not a float\n");
      return((double) MRI_UNSPECIFIED);
    }
  return(v);
}

char *
mri_get_string (MRI_Dataset *ds, const char *key)
{
  MRI_KeyValue *kv;

  kv = FindInHashTable(ds, key, FALSE);
  if (kv == NULL)
    {
      mri_report_error(ds, "mri_get_string: non-existent key %s\n",key);
      return(NULL);
    }
  return(kv->value);
}


/*-------- SETTING KEY VALUES ------------------------------------*/

void
mri_set_int (MRI_Dataset *ds, const char *key, const long long value)
{
  char s[64];

  sprintf(s, "%lld", value);
  mri_set_string(ds, key, s);
}

void
mri_set_float (MRI_Dataset *ds, const char *key, const float value)
{
  char s[64];

  sprintf(s, "%g", value);
  mri_set_string(ds, key, s);
}

void
mri_set_string (MRI_Dataset *ds, const char *key, const char *value)
{
  char *old_value;
  MRI_KeyValue *kv;

  /* check that dataset is writeable */
  if (ds->mode == MRI_READ || ds->mode == MRI_MODIFY_DATA)
    {
      mri_report_error(ds, "libmri: attempt to add key to read-only dataset\n");
      return;
    }

  /* enter key/value into hash table */
  kv = FindInHashTable(ds, key, TRUE);
  old_value = kv->value;
  kv->value = (char *) malloc(strlen(value) + 1);
  strcpy(kv->value, value);
  if (!CheckForHooks(ds, kv))
    {
      free(kv->value);
      kv->value = old_value;
      if (old_value == NULL)
	RemoveFromHashTable(ds, key);
      return;
    }
  if (old_value != NULL)
    free(old_value);
}


/*-------- MISCELLANEOUS KEY FUNCTIONS ---------------------------*/

int
mri_has (MRI_Dataset *ds, const char *key)
{
  return(FindInHashTable(ds, key, FALSE) != NULL);
}

void
mri_remove (MRI_Dataset *ds, const char *key)
{
  MRI_KeyValue *kv;

  /* check that dataset is writeable */
  if (ds->mode == MRI_READ || ds->mode == MRI_MODIFY_DATA)
    {
      mri_report_error(ds, "libmri: attempt to remove key from read-only dataset\n");
      return;
    }

  kv = FindInHashTable(ds, key, FALSE);
  if (kv == NULL)
    return;
  if (!CheckForRemovalHooks(ds, kv))
    return;
  RemoveFromHashTable(ds, kv->key);
}

void
mri_iterate_over_keys (MRI_Dataset *ds)
{
  int n;
  int i;
  MRI_KeyValue *kv;

  if (ds->iteration_table != NULL)
    free(ds->iteration_table);

  /* sort the keys in a table */
  ds->iteration_table = (MRI_KeyValue **) malloc(ds->n_keys * sizeof(MRI_KeyValue *));
  n = 0;
  for (i = 0; i < ds->hash_table_size; ++i)
    for (kv = ds->hash_table[i]; kv != NULL; kv = kv->next_in_hash_table)
      ds->iteration_table[n++] = kv;
  if (n != ds->n_keys)
    {
      mri_report_error(ds, "libmri: internal error -- invalid key count\n");
      abort();
    }
  qsort(ds->iteration_table, n, sizeof(MRI_KeyValue *), CompareKeys);

  ds->iteration_count = n;
  ds->iteration_index = 0;
}

char *
mri_next_key (MRI_Dataset *ds)
{
  if (ds->iteration_table == NULL)
    {
      mri_report_error(ds, "mri_next_key: mri_iterate_over_keys not called first\n");
      return(NULL);
    }
  if (ds->iteration_index >= ds->iteration_count)
    {
      free(ds->iteration_table);
      ds->iteration_table = NULL;
      return(NULL);
    }
  return(ds->iteration_table[ds->iteration_index++]->key);
}

char *
mri_cat (const char *s1, const char *s2)
{
  static char str[MRI_MAX_KEY_LENGTH+1];

  strcpy(str, s1);
  strcat(str, s2);
  return(str);
}


/*-------- CHUNKS ------------------------------------------------*/

void
mri_create_chunk (MRI_Dataset *ds, const char *key)
{
  /* check that dataset is writeable */
  if (ds->mode == MRI_READ || ds->mode == MRI_MODIFY_DATA)
    {
      mri_report_error(ds, "mri_create_chunk: attempt to create chunk in read-only dataset\n");
      return;
    }

  mri_set_string(ds, key, "[chunk]");
}

void
mri_update_chunk (MRI_Dataset *ds, const char *key)
{
  MRI_Chunk *ch;

  for (ch = ds->chunks; ch != NULL; ch = ch->next)
    if (strcmp(ch->name, key) == 0)
      break;
  if (ch == NULL)
    {
      mri_report_error(ds, "mri_update_chunk: no such chunk!\n");
      return;
    }

  if (ds->recompute_positions)
    ComputeChunkPositions(ds);

  if (ch->modified)
    RepositionChunk(ch, NULL);
}

void *
mri_read_chunk (MRI_Dataset *ds, const char *key, long long size,
		long long offset, MRI_ArrayType type, void* buffer)
{
  MRI_Chunk *ch;
  long long i;
  unsigned char *uchar_buf= NULL;
  short *short_buf= NULL;
  int *int_buf= NULL;
  float *float_buf= NULL;
  double *double_buf= NULL;
  long *long_buf= NULL;
  long long *longlong_buf= NULL;
  int conv_error;
  int saved_bio_big_endian_input;
  int saved_bio_error;

  for (ch = ds->chunks; ch != NULL; ch = ch->next)
    if (strcmp(ch->name, key) == 0)
      break;
  if (ch == NULL)
    {
      mri_report_error(ds, "mri_get_chunk: no such chunk named %s\n",
		       key);
      return(NULL);
    }

  /* check that everything is in-bounds */
  if (offset < 0 ||
      (offset+size)*(type == MRI_RAW ? 1 : MRI_TypeLength(ch->datatype)) > ch->size)
    { 
      mri_report_error(ds, "mri_get_chunk: out-of-bounds parameters\n");
      return(NULL);
    }

  if (!ch->ready_to_read &&
      !PrepareToRead(ch))
    return(NULL);
  ch->file->last_use = file_access++;
  saved_bio_big_endian_input = bio_big_endian_input;
  saved_bio_error = bio_error;
  bio_big_endian_input = !ch->little_endian;
  bio_error = FALSE;
  conv_error = FALSE;

  switch (type) {
  case MRI_UNSIGNED_CHAR: 
    uchar_buf= (unsigned char*)buffer;
    break;
  case MRI_SHORT:
    short_buf= (short*)buffer;
    break;
  case MRI_INT:
    int_buf= (int*)buffer;
    break;
  case MRI_LONG:
    long_buf= (long*)buffer;
    break;
  case MRI_LONGLONG:
    longlong_buf= (long long*)buffer;
    break;
  case MRI_FLOAT:
    float_buf= (float*)buffer;
    break;
  case MRI_DOUBLE:
    double_buf= (double*)buffer;
    break;
  }

  if (type == MRI_RAW)
    {
      uchar_buf = (unsigned char*)buffer;
      ChunkSeek(ch, offset);
      FRdUInt8Array(ch->file->fp, uchar_buf, size);
    }
  else switch (ch->datatype)
    {
    case MRI_UINT8:
      if (uchar_buf==NULL) 
	uchar_buf= (unsigned char*)GetBuffer(ds,size*sizeof(char));
      ChunkSeek(ch, offset);
      FRdUInt8Array(ch->file->fp, uchar_buf, size);
      switch (type)
	{
	case MRI_UNSIGNED_CHAR:
	  /* all done */
	  break;

	case MRI_SHORT:
	  for (i = 0; i < size; ++i) short_buf[i] = uchar_buf[i];
	  break;

	case MRI_INT:
	  for (i = 0; i < size; ++i) int_buf[i] = uchar_buf[i];
	  break;

	case MRI_LONG:
	  for (i = 0; i < size; ++i) long_buf[i] = uchar_buf[i];
	  break;

	case MRI_LONGLONG:
	  for (i = 0; i < size; ++i) longlong_buf[i] = uchar_buf[i];
	  break;

	case MRI_FLOAT:
	  for (i = 0; i < size; ++i) float_buf[i] = uchar_buf[i];
	  break;

	case MRI_DOUBLE:
	  for (i = 0; i < size; ++i) double_buf[i] = uchar_buf[i];
	  break;

	default:
	  mri_report_error(ds, "mri_get_chunk: Invalid array type\n");
	  bio_big_endian_input = saved_bio_big_endian_input;
	  bio_error = saved_bio_error;
	  return(NULL);
	}
      break;

    case MRI_INT16:
      if (short_buf==NULL)
	short_buf = (short *)GetBuffer(ds, size*sizeof(short));
      ChunkSeek(ch, 2*offset);
      FRdInt16Array(ch->file->fp, short_buf, size);
      switch (type)
	{
	case MRI_UNSIGNED_CHAR:
	  for (i = 0; i < size; ++i)
	    if (short_buf[i] < 0)
	      {
		uchar_buf[i] = 0;
		conv_error = TRUE;
	      }
	    else if (short_buf[i] > 255)
	      {
		uchar_buf[i] = 255;
		conv_error = TRUE;
	      }
	    else uchar_buf[i] = short_buf[i];
	  break;

	case MRI_SHORT:
	  /* all done */
	  break;

	case MRI_INT:
	  for (i = 0; i < size; ++i) int_buf[i] = short_buf[i];
	  break;

	case MRI_LONG:
	  for (i = 0; i < size; ++i) long_buf[i] = short_buf[i];
	  break;

	case MRI_LONGLONG:
	  for (i = 0; i < size; ++i) longlong_buf[i] = short_buf[i];
	  break;

	case MRI_FLOAT:
	  for (i = 0; i < size; ++i) float_buf[i] = short_buf[i];
	  break;

	case MRI_DOUBLE:
	  for (i = 0; i < size; ++i) double_buf[i] = short_buf[i];
	  break;

	default:
	  mri_report_error(ds, "mri_get_chunk: Invalid array type\n");
	  bio_big_endian_input = saved_bio_big_endian_input;
	  bio_error = saved_bio_error;
	  return(NULL);
	}
      break;

    case MRI_INT32:
      if (int_buf==NULL)
	int_buf = (int *)GetBuffer(ds, size*sizeof(int));
      ChunkSeek(ch, 4*offset);
      FRdInt32Array(ch->file->fp, int_buf, size);
      switch (type)
	{
	case MRI_UNSIGNED_CHAR:
	  for (i = 0; i < size; ++i)
	    if (int_buf[i] < 0)
	      {
		uchar_buf[i] = 0;
		conv_error = TRUE;
	      }
	    else if (int_buf[i] > 255)
	      {
		uchar_buf[i] = 255;
		conv_error = TRUE;
	      }
	    else uchar_buf[i] = int_buf[i];
	  break;

	case MRI_SHORT:
	  for (i = 0; i < size; ++i)
	    if (int_buf[i] < SHRT_MIN)
	      {
		short_buf[i] = SHRT_MIN;
		conv_error = TRUE;
	      }
	    else if (int_buf[i] > SHRT_MAX)
	      {
		short_buf[i] = SHRT_MAX;
		conv_error = TRUE;
	      }
	    else short_buf[i] = int_buf[i];
	  break;

	case MRI_INT:
	  /* all done */
	  break;

	case MRI_LONG:
	  for (i = 0; i < size; ++i) long_buf[i] = int_buf[i];
	  break;

	case MRI_LONGLONG:
	  for (i = 0; i < size; ++i) longlong_buf[i] = int_buf[i];
	  break;

	case MRI_FLOAT:
	  for (i = 0; i < size; ++i) float_buf[i] = int_buf[i];
	  break;

	case MRI_DOUBLE:
	  for (i = 0; i < size; ++i) double_buf[i] = int_buf[i];
	  break;

	default:
	  mri_report_error(ds, "mri_get_chunk: Invalid array type\n");
	  bio_big_endian_input = saved_bio_big_endian_input;
	  bio_error = saved_bio_error;
	  return(NULL);
	}
      break;

    case MRI_INT64:
      if (longlong_buf==NULL)
	longlong_buf = (long long *)GetBuffer(ds, size*sizeof(long long));
      ChunkSeek(ch, 8*offset);
      FRdInt64Array(ch->file->fp, longlong_buf, size);
      switch (type)
	{
	case MRI_UNSIGNED_CHAR:
	  for (i = 0; i < size; ++i)
	    if (longlong_buf[i] < 0)
	      {
		uchar_buf[i] = 0;
		conv_error = TRUE;
	      }
	    else if (longlong_buf[i] > 255)
	      {
		uchar_buf[i] = 255;
		conv_error = TRUE;
	      }
	    else uchar_buf[i] = longlong_buf[i];
	  break;

	case MRI_SHORT:
	  for (i = 0; i < size; ++i)
	    if (longlong_buf[i] < SHRT_MIN)
	      {
		short_buf[i] = SHRT_MIN;
		conv_error = TRUE;
	      }
	    else if (longlong_buf[i] > SHRT_MAX)
	      {
		short_buf[i] = SHRT_MAX;
		conv_error = TRUE;
	      }
	    else short_buf[i] = longlong_buf[i];
	  break;

	case MRI_INT:
	  for (i = 0; i < size; ++i)
	    if (longlong_buf[i] < INT_MIN)
	      {
		int_buf[i] = INT_MIN;
		conv_error = TRUE;
	      }
	    else if (longlong_buf[i] > INT_MAX)
	      {
		int_buf[i] = INT_MAX;
		conv_error = TRUE;
	      }
	    else int_buf[i] = longlong_buf[i];
	  break;

	case MRI_LONG:
	  for (i = 0; i < size; ++i)
	    if (longlong_buf[i] < LONG_MIN)
	      {
		long_buf[i] = LONG_MIN;
		conv_error = TRUE;
	      }
	    else if (longlong_buf[i] > LONG_MAX)
	      {
		long_buf[i] = LONG_MAX;
		conv_error = TRUE;
	      }
	    else long_buf[i] = longlong_buf[i];
	  break;

	case MRI_LONGLONG:
	  /* All done */
	  break;

	case MRI_FLOAT:
	  for (i = 0; i < size; ++i) float_buf[i] = longlong_buf[i];
	  break;

	case MRI_DOUBLE:
	  for (i = 0; i < size; ++i) double_buf[i] = longlong_buf[i];
	  break;

	default:
	  mri_report_error(ds, "mri_get_chunk: Invalid array type\n");
	  bio_big_endian_input = saved_bio_big_endian_input;
	  bio_error = saved_bio_error;
	  return(NULL);
	}
      break;

    case MRI_FLOAT32:
      if (float_buf==NULL)
	float_buf = (float *)GetBuffer(ds, size*sizeof(float));
      ChunkSeek(ch, 4*offset);
      FRdFloat32Array(ch->file->fp, float_buf, size);
      switch (type)
	{
	case MRI_UNSIGNED_CHAR:
	  for (i = 0; i < size; ++i)
	    if (float_buf[i] < 0.0)
	      {
		uchar_buf[i] = 0;
		conv_error = TRUE;
	      }
	    else if (float_buf[i] > 255.0)
	      {
		uchar_buf[i] = 255;
		conv_error = TRUE;
	      }
	    else uchar_buf[i] = (int) float_buf[i];
	  break;

	case MRI_SHORT:
	  for (i = 0; i < size; ++i)
	    if (float_buf[i] < (float) SHRT_MIN)
	      {
		short_buf[i] = SHRT_MIN;
		conv_error = TRUE;
	      }
	    else if (float_buf[i] > (float)SHRT_MAX)
	      {
		short_buf[i] = SHRT_MAX;
		conv_error = TRUE;
	      }
	    else short_buf[i] = (int) float_buf[i];
	  break;

	case MRI_INT:
	  for (i = 0; i < size; ++i)
	    if (float_buf[i] < (float)INT_MIN)
	      {
		int_buf[i] = INT_MIN;
		conv_error = TRUE;
	      }
	    else if (float_buf[i] > (float)INT_MAX)
	      {
		int_buf[i] = INT_MAX;
		conv_error = TRUE;
	      }
	    else int_buf[i] = float_buf[i];
	  break;

	case MRI_LONG:
	  for (i = 0; i < size; ++i)
	    if (float_buf[i] < (float) LONG_MIN)
	      {
		long_buf[i] = LONG_MIN;
		conv_error = TRUE;
	      }
	    else if (float_buf[i] > (float)LONG_MAX)
	      {
		long_buf[i] = LONG_MAX;
		conv_error = TRUE;
	      }
	    else long_buf[i] = float_buf[i];
	  break;

	case MRI_LONGLONG:
	  for (i = 0; i < size; ++i)
	    if (float_buf[i] < (float) LONGLONG_MIN)
	      {
		longlong_buf[i] = LONGLONG_MIN;
		conv_error = TRUE;
	      }
	    else if (float_buf[i] > (float) LONGLONG_MAX)
	      {
		longlong_buf[i] = LONGLONG_MAX;
		conv_error = TRUE;
	      }
	    else longlong_buf[i] = float_buf[i];
	  break;

	case MRI_FLOAT:
	  /* all done */
	  break;

	case MRI_DOUBLE:
	  for (i = 0; i < size; ++i) double_buf[i] = float_buf[i];
	  break;

	default:
	  mri_report_error(ds, "mri_get_chunk: Invalid array type\n");
	  bio_big_endian_input = saved_bio_big_endian_input;
	  bio_error = saved_bio_error;
	  return(NULL);
	}
      break;

    case MRI_FLOAT64:
      if (double_buf==NULL)
	double_buf = (double *)GetBuffer(ds, size*sizeof(double));
      ChunkSeek(ch, 8*offset);
      FRdFloat64Array(ch->file->fp, double_buf, size);
      switch (type)
	{
	case MRI_UNSIGNED_CHAR:
	  for (i = 0; i < size; ++i)
	    if (double_buf[i] < 0.0)
	      {
		uchar_buf[i] = 0;
		conv_error = TRUE;
	      }
	    else if (double_buf[i] > 255.0)
	      {
		uchar_buf[i] = 255;
		conv_error = TRUE;
	      }
	    else uchar_buf[i] = (int) double_buf[i];
	  break;

	case MRI_SHORT:
	  for (i = 0; i < size; ++i)
	    if (double_buf[i] < (double)SHRT_MIN)
	      {
		short_buf[i] = SHRT_MIN;
		conv_error = TRUE;
	      }
	    else if (double_buf[i] > (double)SHRT_MAX)
	      {
		short_buf[i] = SHRT_MAX;
		conv_error = TRUE;
	      }
	    else short_buf[i] = (int) double_buf[i];
	  break;

	case MRI_INT:
	  for (i = 0; i < size; ++i)
	    if (double_buf[i] < (double)INT_MIN)
	      {
		int_buf[i] = INT_MIN;
		conv_error = TRUE;
	      }
	    else if (double_buf[i] > (double)INT_MAX)
	      {
		int_buf[i] = INT_MAX;
		conv_error = TRUE;
	      }
	    else int_buf[i] = double_buf[i];
	  break;

	case MRI_LONG:
	  for (i = 0; i < size; ++i)
	    if (double_buf[i] < (double)LONG_MIN)
	      {
		long_buf[i] = LONG_MIN;
		conv_error = TRUE;
	      }
	    else if (double_buf[i] > (double)LONG_MAX)
	      {
		long_buf[i] = LONG_MAX;
		conv_error = TRUE;
	      }
	    else long_buf[i] = double_buf[i];
	  break;

	case MRI_LONGLONG:
	  for (i = 0; i < size; ++i)
	    if (double_buf[i] < (double) LONGLONG_MIN)
	      {
		longlong_buf[i] = LONGLONG_MIN;
		conv_error = TRUE;
	      }
	    else if (double_buf[i] > (double) LONGLONG_MAX)
	      {
		longlong_buf[i] = LONGLONG_MAX;
		conv_error = TRUE;
	      }
	    else longlong_buf[i] = double_buf[i];
	  break;

	case MRI_FLOAT:
	  for (i = 0; i < size; ++i)
	    if (double_buf[i] < -MAXFLOAT)
	      {
		float_buf[i] = -MAXFLOAT;
		conv_error = TRUE;
	      }
	    else if (double_buf[i] > MAXFLOAT)
	      {
		float_buf[i] = MAXFLOAT;
		conv_error = TRUE;
	      }
	    else float_buf[i] = (float) double_buf[i];
	  break;

	case MRI_DOUBLE:
	  /* all done */
	  break;

	default:
	  mri_report_error(ds, "mri_get_chunk: Invalid array type\n");
	  bio_big_endian_input = saved_bio_big_endian_input;
	  bio_error = saved_bio_error;
	  return(NULL);
	}
      break;

    default:
      mri_report_error(ds, "mri_read_chunk: Invalid chunk datatype\n");
      bio_big_endian_input = saved_bio_big_endian_input;
      bio_error = saved_bio_error;
      return(NULL);
    }
      
  if (bio_error)
    {
      mri_report_error(ds, "mri_read_chunk: could not read chunk data from file %s\n",
		       ch->file->name);
      bio_big_endian_input = saved_bio_big_endian_input;
      bio_error = saved_bio_error;
      return(NULL);
    }
  if (conv_error)
    {
      mri_report_warning(ds, "mri_read_chunk: out-of-range conversions\n");
      bio_big_endian_input = saved_bio_big_endian_input;
      bio_error = saved_bio_error;
    }
  bio_big_endian_input = saved_bio_big_endian_input;
  bio_error = saved_bio_error;
  return(buffer);
}

void *
mri_get_chunk (MRI_Dataset *ds, const char *key, long long size,
	       long long offset, MRI_ArrayType type)
{

  return mri_read_chunk(ds, key, size, offset, type,
			GetBuffer(ds, MRI_ArrayTypeLength(type,size)));
}

#if (HPPA || HPPA20)
#pragma optimize off
#endif
static void CopyConvFloatLonglong(float* float_buf,
				  long long* longlong_buf,
				  int size,int* error)
{
  int i;
  for (i = 0; i < size; ++i)
    if (float_buf[i] < (float)LONGLONG_MIN)
      {
	longlong_buf[i] = LONGLONG_MIN;
	*error = TRUE;
      }
    else if (float_buf[i] > (float)LONGLONG_MAX)
      {
	longlong_buf[i] = LONGLONG_MAX;
	*error = TRUE;
      }
    else longlong_buf[i] = float_buf[i];
}

static void CopyConvDoubleLonglong(double* double_buf,
				   long long* longlong_buf,
				   int size,int* error)
{
  int i;
  for (i = 0; i < size; ++i)
    if (double_buf[i] < (double)LONGLONG_MIN)
      {
	longlong_buf[i] = LONGLONG_MIN;
	*error = TRUE;
      }
    else if (double_buf[i] > (double)LONGLONG_MAX)
      {
	longlong_buf[i] = LONGLONG_MAX;
	*error = TRUE;
      }
    else longlong_buf[i] = double_buf[i];
}
#if (HPPA || HPPA20)
#pragma optimize on
#endif

void
mri_set_chunk (MRI_Dataset *ds, const char *key,
	       long long size, long long offset,
	       MRI_ArrayType type, void *buf)
{
  MRI_Chunk *ch;
  int i;
  unsigned char *uchar_buf;
  short *short_buf;
  int *int_buf;
  long *long_buf;
  long long *longlong_buf;
  float *float_buf;
  double *double_buf;
  int conv_error;
  int saved_bio_big_endian_output;
  int saved_bio_error;

  for (ch = ds->chunks; ch != NULL; ch = ch->next)
    if (strcmp(ch->name, key) == 0)
      break;
  if (ch == NULL)
    {
      mri_report_error(ds, "mri_set_chunk: no such chunk named %s\n", key);
      return;
    }

  /* check that everything is in-bounds */
  if (offset < 0 ||
      (offset+size)*(type == MRI_RAW ? 1 : MRI_TypeLength(ch->datatype)) > ch->size)
    {      
       mri_report_error(ds, "mri_set_chunk: out-of-bounds parameters\n");
      return;
    }

  if (!ch->ready_to_write &&
      !PrepareToWrite(ch))
    return;
  ch->file->last_use = file_access++;
  saved_bio_big_endian_output = bio_big_endian_output;
  saved_bio_error = bio_error;
  bio_big_endian_output = !ch->little_endian;
  bio_error = FALSE;
  conv_error = FALSE;

  switch (type)
    {
    case MRI_RAW:
      ChunkSeek(ch, offset);
      FWrUInt8Array(ch->file->fp, buf, size);
      break;

    case MRI_UNSIGNED_CHAR:
      switch (ch->datatype)
	{
	case MRI_UINT8:
	  ChunkSeek(ch, offset);
	  FWrUInt8Array(ch->file->fp, buf, size);
	  break;
	case MRI_INT16:
	  uchar_buf = (unsigned char *) buf;
	  short_buf = (short *) GetBuffer(ds, size*sizeof(short));
	  for (i = 0; i < size; ++i)
	    short_buf[i] = uchar_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_SHORT, short_buf);
	  break;
	case MRI_INT32:
	  uchar_buf = (unsigned char *) buf;
	  int_buf = (int *) GetBuffer(ds, size*sizeof(int));
	  for (i = 0; i < size; ++i)
	    int_buf[i] = uchar_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_INT, int_buf);
	  break;
	case MRI_INT64:
	  uchar_buf = (unsigned char *) buf;
	  longlong_buf = (long long *) GetBuffer(ds, size*sizeof(long long));
	  for (i = 0; i < size; ++i)
	    longlong_buf[i] = uchar_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_LONGLONG, longlong_buf);
	  break;
	case MRI_FLOAT32:
	  uchar_buf = (unsigned char *) buf;
	  float_buf = (float *) GetBuffer(ds, size*sizeof(float));
	  for (i = 0; i < size; ++i)
	    float_buf[i] = uchar_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_FLOAT, float_buf);
	  break;
	case MRI_FLOAT64:
	  uchar_buf = (unsigned char *) buf;
	  double_buf = (double *) GetBuffer(ds, size*sizeof(double));
	  for (i = 0; i < size; ++i)
	    double_buf[i] = uchar_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_DOUBLE, double_buf);
	  break;
	default:
	  mri_report_error(ds, "mri_set_chunk: Invalid array type\n");
	  break;
	}
      break;

    case MRI_SHORT:
      short_buf = (short *) buf;
      switch (ch->datatype)
	{
	case MRI_UINT8:
	  uchar_buf = (unsigned char *) GetBuffer(ds, size);
	  for (i = 0; i < size; ++i)
	    if (short_buf[i] < 0)
	      {
		uchar_buf[i] = 0;
		conv_error = TRUE;
	      }
	    else if (short_buf[i] > 255)
	      {
		uchar_buf[i] = 255;
		conv_error = TRUE;
	      }
	    else uchar_buf[i] = short_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_UNSIGNED_CHAR, uchar_buf);
	  break;
	case MRI_INT16:
	  ChunkSeek(ch, 2*offset);
	  FWrInt16Array(ch->file->fp, buf, size);
	  break;
	case MRI_INT32:
	  int_buf = (int *) GetBuffer(ds, size*sizeof(int));
	  for (i = 0; i < size; ++i)
	    int_buf[i] = short_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_INT, int_buf);
	  break;
	case MRI_INT64:
	  longlong_buf = (long long *) GetBuffer(ds, size*sizeof(long long));
	  for (i = 0; i < size; ++i)
	    longlong_buf[i] = short_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_LONGLONG, longlong_buf);
	  break;
	case MRI_FLOAT32:
	  float_buf = (float *) GetBuffer(ds, size*sizeof(float));
	  for (i = 0; i < size; ++i)
	    float_buf[i] = short_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_FLOAT, float_buf);
	  break;
	case MRI_FLOAT64:
	  double_buf = (double *) GetBuffer(ds, size*sizeof(double));
	  for (i = 0; i < size; ++i)
	    double_buf[i] = short_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_DOUBLE, double_buf);
	  break;
	default:
	  mri_report_error(ds, "mri_set_chunk: Invalid array type\n");
	  break;
	}
      break;

    case MRI_INT:
      int_buf = (int *) buf;
      switch (ch->datatype)
	{
	case MRI_UINT8:
	  uchar_buf = (unsigned char *) GetBuffer(ds, size);
	  for (i = 0; i < size; ++i)
	    if (int_buf[i] < 0)
	      {
		uchar_buf[i] = 0;
		conv_error = TRUE;
	      }
	    else if (int_buf[i] > 255)
	      {
		uchar_buf[i] = 255;
		conv_error = TRUE;
	      }
	    else uchar_buf[i] = int_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_UNSIGNED_CHAR, uchar_buf);
	  break;
	case MRI_INT16:
	  short_buf = (short *) GetBuffer(ds, size*sizeof(short));
	  for (i = 0; i < size; ++i)
	    if (int_buf[i] < SHRT_MIN)
	      {
		short_buf[i] = SHRT_MIN;
		conv_error = TRUE;
	      }
	    else if (int_buf[i] > SHRT_MAX)
	      {
		short_buf[i] = SHRT_MAX;
		conv_error = TRUE;
	      }
	    else short_buf[i] = int_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_SHORT, short_buf);
	  break;
	case MRI_INT32:
	  ChunkSeek(ch, 4*offset);
	  FWrInt32Array(ch->file->fp, buf, size);
	  break;
	case MRI_INT64:
	  longlong_buf = (long long *) GetBuffer(ds, size*sizeof(long long));
	  for (i = 0; i < size; ++i)
	    longlong_buf[i] = int_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_LONGLONG, longlong_buf);
	  break;
	case MRI_FLOAT32:
	  float_buf = (float *) GetBuffer(ds, size*sizeof(float));
	  for (i = 0; i < size; ++i)
	    float_buf[i] = int_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_FLOAT, float_buf);
	  break;
	case MRI_FLOAT64:
	  double_buf = (double *) GetBuffer(ds, size*sizeof(double));
	  for (i = 0; i < size; ++i)
	    double_buf[i] = int_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_DOUBLE, double_buf);
	  break;
	default:
	  mri_report_error(ds, "mri_set_chunk: Invalid array type\n");
	  break;
	}
      break;

    case MRI_LONG:
      long_buf = (long *) buf;
      switch (ch->datatype)
	{
	case MRI_UINT8:
	  uchar_buf = (unsigned char *) GetBuffer(ds, size);
	  for (i = 0; i < size; ++i)
	    if (long_buf[i] < 0)
	      {
		uchar_buf[i] = 0;
		conv_error = TRUE;
	      }
	    else if (long_buf[i] > 255)
	      {
		uchar_buf[i] = 255;
		conv_error = TRUE;
	      }
	    else uchar_buf[i] = long_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_UNSIGNED_CHAR, uchar_buf);
	  break;
	case MRI_INT16:
	  short_buf = (short *) GetBuffer(ds, size*sizeof(short));
	  for (i = 0; i < size; ++i)
	    if (long_buf[i] < SHRT_MIN)
	      {
		short_buf[i] = SHRT_MIN;
		conv_error = TRUE;
	      }
	    else if (long_buf[i] > SHRT_MAX)
	      {
		short_buf[i] = SHRT_MAX;
		conv_error = TRUE;
	      }
	    else short_buf[i] = long_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_SHORT, short_buf);
	  break;
	case MRI_INT32:
	  int_buf = (int *) GetBuffer(ds, size*sizeof(int));
	  for (i = 0; i < size; ++i)
	    if (long_buf[i] < INT_MIN)
	      {
		int_buf[i] = INT_MIN;
		conv_error = TRUE;
	      }
	    else if (long_buf[i] > INT_MAX)
	      {
		int_buf[i] = INT_MAX;
		conv_error = TRUE;
	      }
	    else int_buf[i] = long_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_INT, int_buf);
	  break;
	case MRI_INT64:
	  longlong_buf = (long long *) GetBuffer(ds, size*sizeof(long long));
	  for (i = 0; i < size; ++i)
	    longlong_buf[i] = long_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_LONGLONG, longlong_buf);
	  break;
	case MRI_FLOAT32:
	  float_buf = (float *) GetBuffer(ds, size*sizeof(float));
	  for (i = 0; i < size; ++i)
	    float_buf[i] = long_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_FLOAT, float_buf);
	  break;
	case MRI_FLOAT64:
	  double_buf = (double *) GetBuffer(ds, size*sizeof(double));
	  for (i = 0; i < size; ++i)
	    double_buf[i] = long_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_DOUBLE, double_buf);
	  break;
	default:
	  mri_report_error(ds, "mri_set_chunk: Invalid array type\n");
	  break;
	}
      break;

    case MRI_LONGLONG:
      longlong_buf = (long long *) buf;
      switch (ch->datatype)
	{
	case MRI_UINT8:
	  uchar_buf = (unsigned char *) GetBuffer(ds, size);
	  for (i = 0; i < size; ++i)
	    if (longlong_buf[i] < 0)
	      {
		uchar_buf[i] = 0;
		conv_error = TRUE;
	      }
	    else if (longlong_buf[i] > 255)
	      {
		uchar_buf[i] = 255;
		conv_error = TRUE;
	      }
	    else uchar_buf[i] = longlong_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_UNSIGNED_CHAR, uchar_buf);
	  break;
	case MRI_INT16:
	  short_buf = (short *) GetBuffer(ds, size*sizeof(short));
	  for (i = 0; i < size; ++i)
	    if (longlong_buf[i] < SHRT_MIN)
	      {
		short_buf[i] = SHRT_MIN;
		conv_error = TRUE;
	      }
	    else if (longlong_buf[i] > SHRT_MAX)
	      {
		short_buf[i] = SHRT_MAX;
		conv_error = TRUE;
	      }
	    else short_buf[i] = longlong_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_SHORT, short_buf);
	  break;
	case MRI_INT32:
	  int_buf = (int *) GetBuffer(ds, size*sizeof(int));
	  for (i = 0; i < size; ++i)
	    if (longlong_buf[i] < INT_MIN)
	      {
		int_buf[i] = INT_MIN;
		conv_error = TRUE;
	      }
	    else if (longlong_buf[i] > INT_MAX)
	      {
		int_buf[i] = INT_MAX;
		conv_error = TRUE;
	      }
	    else int_buf[i] = longlong_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_INT, int_buf);
	  break;
	case MRI_INT64:
	  ChunkSeek(ch, 8*offset);
	  FWrInt64Array(ch->file->fp, buf, size);
	  break;
	case MRI_FLOAT32:
	  float_buf = (float *) GetBuffer(ds, size*sizeof(float));
	  for (i = 0; i < size; ++i)
	    float_buf[i] = longlong_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_FLOAT, float_buf);
	  break;
	case MRI_FLOAT64:
	  double_buf = (double *) GetBuffer(ds, size*sizeof(double));
	  for (i = 0; i < size; ++i)
	    double_buf[i] = longlong_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_DOUBLE, double_buf);
	  break;
	default:
	  mri_report_error(ds, "mri_set_chunk: Invalid array type\n");
	  break;
	}
      break;

    case MRI_FLOAT:
      float_buf = (float *) buf;
      switch (ch->datatype)
	{
	case MRI_UINT8:
	  uchar_buf = (unsigned char *) GetBuffer(ds, size);
	  for (i = 0; i < size; ++i)
	    if (float_buf[i] < 0.0)
	      {
		uchar_buf[i] = 0;
		conv_error = TRUE;
	      }
	    else if (float_buf[i] > 255.0)
	      {
		uchar_buf[i] = 255;
		conv_error = TRUE;
	      }
	    else uchar_buf[i] = (int) float_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_UNSIGNED_CHAR, uchar_buf);
	  break;
	case MRI_INT16:
	  short_buf = (short *) GetBuffer(ds, size*sizeof(short));
	  for (i = 0; i < size; ++i)
	    if (float_buf[i] < (float)SHRT_MIN)
	      {
		short_buf[i] = SHRT_MIN;
		conv_error = TRUE;
	      }
	    else if (float_buf[i] > (float)SHRT_MAX)
	      {
		short_buf[i] = SHRT_MAX;
		conv_error = TRUE;
	      }
	    else short_buf[i] = (int) float_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_SHORT, short_buf);
	  break;
	case MRI_INT32:
	  int_buf = (int *) GetBuffer(ds, size*sizeof(int));
	  for (i = 0; i < size; ++i)
	    if (float_buf[i] < (float)INT_MIN)
	      {
		int_buf[i] = INT_MIN;
		conv_error = TRUE;
	      }
	    else if (float_buf[i] > (float)INT_MAX)
	      {
		int_buf[i] = INT_MAX;
		conv_error = TRUE;
	      }
	    else int_buf[i] = float_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_INT, int_buf);
	  break;
	case MRI_INT64:
	  {
	  longlong_buf = (long long *) GetBuffer(ds, size*sizeof(long long));
	  CopyConvFloatLonglong(float_buf,longlong_buf,size,&conv_error);
	  mri_set_chunk(ds, key, size, offset, MRI_LONGLONG, longlong_buf);
	  }
	  break;
	case MRI_FLOAT32:
	  ChunkSeek(ch, 4*offset);
	  FWrFloat32Array(ch->file->fp, buf, size);
	  break;
	case MRI_FLOAT64:
	  double_buf = (double *) GetBuffer(ds, size*sizeof(double));
	  for (i = 0; i < size; ++i)
	    double_buf[i] = float_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_DOUBLE, double_buf);
	  break;
	default:
	  mri_report_error(ds, "mri_set_chunk: Invalid array type\n");
	  break;
	}
      break;

    case MRI_DOUBLE:
      double_buf = (double *) buf;
      switch (ch->datatype)
	{
	case MRI_UINT8:
	  uchar_buf = (unsigned char *) GetBuffer(ds, size);
	  for (i = 0; i < size; ++i)
	    if (double_buf[i] < 0.0)
	      {
		uchar_buf[i] = 0;
		conv_error = TRUE;
	      }
	    else if (double_buf[i] > 255.0)
	      {
		uchar_buf[i] = 255;
		conv_error = TRUE;
	      }
	    else uchar_buf[i] = (int) double_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_UNSIGNED_CHAR, uchar_buf);
	  break;
	case MRI_INT16:
	  short_buf = (short *) GetBuffer(ds, size*sizeof(short));
	  for (i = 0; i < size; ++i)
	    if (double_buf[i] < (double)SHRT_MIN)
	      {
		short_buf[i] = SHRT_MIN;
		conv_error = TRUE;
	      }
	    else if (double_buf[i] > (double)SHRT_MAX)
	      {
		short_buf[i] = SHRT_MAX;
		conv_error = TRUE;
	      }
	    else short_buf[i] = (int) double_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_SHORT, short_buf);
	  break;
	case MRI_INT32:
	  int_buf = (int *) GetBuffer(ds, size*sizeof(int));
	  for (i = 0; i < size; ++i)
	    if (double_buf[i] < (double)INT_MIN)
	      {
		int_buf[i] = INT_MIN;
		conv_error = TRUE;
	      }
	    else if (double_buf[i] > (double)INT_MAX)
	      {
		int_buf[i] = INT_MAX;
		conv_error = TRUE;
	      }
	    else int_buf[i] = double_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_INT, int_buf);
	  break;
	case MRI_INT64:
	  longlong_buf = (long long *) GetBuffer(ds, size*sizeof(long long));
	  CopyConvDoubleLonglong(double_buf,longlong_buf,size,&conv_error);
	  mri_set_chunk(ds, key, size, offset, MRI_LONGLONG, longlong_buf);
	  break;
	case MRI_FLOAT32:
	  float_buf = (float *) GetBuffer(ds, size*sizeof(float));
	  for (i = 0; i < size; ++i)
	    if (finite(double_buf[i])) {
	      if (double_buf[i] < -MAXFLOAT)
		{
		  float_buf[i] = -MAXFLOAT;
		  conv_error = TRUE;
		}
	      else if (double_buf[i] > MAXFLOAT)
		{
		  float_buf[i] = MAXFLOAT;
		  conv_error = TRUE;
		}
	      else float_buf[i] = (float) double_buf[i];
	    }
	    else {
	      /* We'll let the compiler translate the NaN or Inf */
	      float_buf[i]= (float)double_buf[i];
	    }
	  mri_set_chunk(ds, key, size, offset, MRI_FLOAT, float_buf);
	  break;
	case MRI_FLOAT64:
	  ChunkSeek(ch, 8*offset);
	  FWrFloat64Array(ch->file->fp, buf, size);
	  break;
	default:
	  mri_report_error(ds, "mri_set_chunk: Invalid array type\n");
	  break;
	}
      break;

    default:
      mri_report_error(ds, "mri_set_chunk: invalid array type\n");
      break;
    }

  if (bio_error)
    mri_report_error(ds, "mri_set_chunk: could not write chunk data to file %s\n",
		     ch->file->name);
  if (conv_error)
    mri_report_warning(ds, "mri_set_chunk: out-of-range conversions\n");
  bio_big_endian_output = saved_bio_big_endian_output;
  bio_error = saved_bio_error;
}


/*--------- BUFFER MANAGEMENT ------------------------------------*/

void
mri_retain_buffer (MRI_Dataset *ds, void *ptr)
{
  MRI_Buffer *pb, *b;

  for (pb = NULL, b = ds->buffers; b != NULL; pb = b, b = b->next)
    if (b->buffer == ptr)
      {
	if (pb != NULL)
	  pb->next = b->next;
	else
	  ds->buffers = b->next;
	b->next = ds->retained_buffers;
	ds->retained_buffers = b;
	return;
      }

  for (b = ds->retained_buffers; b != NULL; b = b->next)
    if (b->buffer == ptr)
      return;

  mri_report_error(ds, "mri_retain_buffer: buffer parameter is not valid\n");
}

void 
mri_discard_buffer (MRI_Dataset *ds, void *ptr)
{
  MRI_Buffer *pb, *b;

  for (pb = NULL, b = ds->buffers; b != NULL; pb = b, b = b->next)
    if (b->buffer == ptr)
      {
	if (pb != NULL)
	  pb->next = b->next;
	else
	  ds->buffers = b->next;
	free(b->buffer);
	free(b);
	return;
      }

  for (pb = NULL, b = ds->retained_buffers; b != NULL; pb = b, b = b->next)
    if (b->buffer == ptr)
      {
	if (pb != NULL)
	  pb->next = b->next;
	else
	  ds->retained_buffers = b->next;
	free(b->buffer);
	free(b);
	return;
      }

  mri_report_error(ds, "mri_discard_buffer: buffer parameter is not valid\n");
}


/*--------- ERROR HANDLING ---------------------------------------*/

void mri_set_error_handling (MRI_Error_Handling mode)
{
  if (mode != MRI_ABORT_ON_ERROR &&
      mode != MRI_REPORT_ERRORS &&
      mode != MRI_IGNORE_ERRORS)
    {
      fprintf(stderr, "mri_set_error_handling: invalid mode specified\n");
      abort();
    }
  error_handling = mode;
}

/*--------- IMAGES -----------------------------------------------*/

void *
mri_get_image (MRI_Dataset *ds, const int time, const int slice,
	       MRI_ArrayType type)
{
  if (!ds->std_images)
    {
      mri_report_error(ds, "libmri: non-standard images prevent mri_get_image\n");
      return(NULL);
    }
  switch (type & 0xf0)
    {
    case MRI_SCALAR:
      if (ds->std_image_vector_size != 1)
	{
	  mri_report_error(ds, "mri_get_image: image is not composed of scalars\n");
	  return(NULL);
	}
      break;
    case MRI_COMPLEX:
      if (ds->std_image_vector_size != 2)
	{
	  mri_report_error(ds, "mri_get_image: image is not composed of complex values\n");
	  return(NULL);
	}
      break;
    case MRI_VECTOR:
      break;
    default:
      mri_report_error(ds, "mri_get_image: invalid array type\n");
      return(NULL);
    }
  return(mri_get_chunk(ds, "images",
		       ds->std_image_size,
		       (((long long) time) * ds->std_n_slices + slice) *
		         ds->std_image_size,
		       type & 0xf));
}

void
mri_set_image (MRI_Dataset *ds, const int time, const int slice,
	       MRI_ArrayType type, void *image)
{
  if (!ds->std_images)
    {
      mri_report_error(ds, "libmri: non-standard images prevent mri_set_image\n");
      return;
    }
  switch (type & 0xf0)
    {
    case MRI_SCALAR:
      if (ds->std_image_vector_size != 1)
	{
	  mri_report_error(ds, "mri_set_image: image is not composed of scalars\n");
	  return;
	}
      break;
    case MRI_COMPLEX:
      if (ds->std_image_vector_size != 2)
	{
	  mri_report_error(ds, "mri_set_image: image is not composed of complex values\n");
	  return;
	}
      break;
    case MRI_VECTOR:
      break;
    default:
      mri_report_error(ds, "mri_set_image: invalid array type\n");
      return;
    }
  mri_set_chunk(ds, "images",
		ds->std_image_size,
		(((long long) time) * ds->std_n_slices + slice) *
		  ds->std_image_size,
		type & 0xf, image);
}


/*------------------------INTERNAL FUNCTIONS---------------------------------*/
   
static void
CreateNewDataset (MRI_Dataset *ds)
{
  mri_set_string(ds, "!format", "pgh");
  mri_set_string(ds, "!version", "1.0");
}

static MRI_Chunk *
NewChunk (MRI_Dataset *ds, char *name)
{
  MRI_Chunk *ch;
  char *s;
  long long v;
  int i;
  char fieldname[16];

  ch = (MRI_Chunk *) malloc(sizeof(MRI_Chunk));
  ch->ds = ds;
  ch->name = (char *) malloc(strlen(name)+1);
  strcpy(ch->name, name);
  ch->next = ds->chunks;
  ds->chunks = ch;

  if (!mri_has(ds, mri_cat(name, ".datatype")))
    ch->datatype = MRI_INT16;
  else
    if (!ConvertDatatype(&ch->datatype,  mri_get_string(ds, mri_cat(name, ".datatype"))))
      {
	mri_report_error(ds, "libmri: Invalid chunk datatype.\n");
	ds->chunks = ch->next;
	free(ch->name);
	free(ch);
	return(NULL);
      }

  if (!mri_has(ds, mri_cat(name, ".dimensions")))
    s = "xyzt";
  else
    s = mri_get_string(ds, mri_cat(name, ".dimensions"));
  ch->dimensions = (char *) malloc(strlen(s) + 1);
  strcpy(ch->dimensions, s);
  for (i = 0; i < (int) strlen(s); ++i)
    {
      sprintf(fieldname, ".extent.%c", s[i]);
      if (!mri_has(ds, mri_cat(name, fieldname)))
	ch->extent[i] = 1;
      else
	ch->extent[i] = mri_get_int(ds, mri_cat(name, fieldname));
    }

  if (mri_has(ds, mri_cat(name, ".little_endian")) &&
      mri_get_int(ds, mri_cat(name, ".little_endian")) == 1)
    ch->little_endian = TRUE;
  else
    ch->little_endian = FALSE;

  if (!mri_has(ds, mri_cat(name, ".order")))
    ch->order = 0;
  else
    if (strcmp((s = mri_get_string(ds, mri_cat(name, ".order"))), "fixed_offset") == 0)
      ch->order = MRI_FIXED_OFFSET;
    else if (strcmp(s, "external") == 0)
      ch->order = MRI_EXTERNAL;
    else
      ch->order = mri_get_int(ds, mri_cat(name, ".order"));

  if (!mri_has(ds, mri_cat(name, ".file")))
    s = "";
  else
    s = mri_get_string(ds, mri_cat(name, ".file"));
  ch->file = GetFile(ds, s);

  if (mri_has(ds, mri_cat(name, ".offset")))
    ch->offset = mri_get_int(ds, mri_cat(name, ".offset"));
  else
    ch->offset = 0;

  v = MRI_TypeLength(ch->datatype);
  for (i = 0; i < (int) strlen(ch->dimensions); ++i)
    v *= ch->extent[i];
  if (mri_has(ds, mri_cat(name, ".size")))
    {
      ch->size = mri_get_int(ds, mri_cat(name, ".size"));
      if (ch->size != v)
	{
	  mri_report_warning(ds, "libmri: %s.size field is being corrected to %d\n", name, v);
	  ch->size = v;
	}
    }
  else
    ch->size = v;

  /* set up the actual_* fields */
  ch->actual_file = ch->file;
  ch->actual_datatype = ch->datatype;
  ch->actual_dimensions = malloc(strlen(ch->dimensions)+1);
  strcpy(ch->actual_dimensions, ch->dimensions);
  for (i = 0; i < (int) strlen(ch->dimensions); ++i)
    ch->actual_extent[i] = ch->extent[i];
  ch->actual_little_endian = ch->little_endian;
  ch->actual_offset = ch->offset;
  ch->actual_size = ch->size;

  ch->modified = FALSE;
  ch->checked = FALSE;
  ch->repositioning = FALSE;
  ch->ready_to_read = FALSE;
  ch->ready_to_write = FALSE;

  CheckForStdImages(ds);
  return(ch);
}

static int
ConvertDatatype (MRI_Datatype *pdt, char *s)
{
  if (strcmp(s, "uint8") == 0)
    *pdt = MRI_UINT8;
  else if (strcmp(s, "int16") == 0)
    *pdt = MRI_INT16;
  else if (strcmp(s, "int32") == 0)
    *pdt = MRI_INT32;
  else if (strcmp(s, "float32") == 0)
    *pdt = MRI_FLOAT32;
  else if (strcmp(s, "float64") == 0)
    *pdt = MRI_FLOAT64;
  else if (strcmp(s, "int64") == 0)
    *pdt = MRI_INT64;
  else return(FALSE);
  return(TRUE);
}

static void
CleanFiles (MRI_Dataset *ds)
{
  MRI_Chunk *ch;
  MRI_File *f;
  int i;
  int n_empty_blocks;
  EmptyBlock empty_blocks[MRI_MAX_CHUNKS+2];

  for (f = ds->files; f != NULL; f = f->next)
    if (!f->external)
      {
	/* set up the empty block array */
	if (f == ds->header_file)
	  empty_blocks[0].start = ds->header_size;
	else
	  empty_blocks[0].start = 0;
	empty_blocks[0].end = MAX_CHUNK_SIZE-1;
	n_empty_blocks = 1;
	
	for (ch = ds->chunks; ch != NULL; ch = ch->next)
	  if (ch->file == f)
	    (void) ReserveBlock(empty_blocks, &n_empty_blocks,
				ch->offset, ch->size);

	if (!OpenFile(f, TRUE))
	  abort();
	for (i = 0; i < n_empty_blocks; ++i)
	  if (empty_blocks[i].end < (MAX_CHUNK_SIZE-1))
	    ClearBlock(f, empty_blocks[i].start,
		       empty_blocks[i].end - empty_blocks[i].start + 1);
	  else
	    {
	      fflush(f->fp);
	      mri_ftruncate(fileno(f->fp), empty_blocks[i].start);
	    }
      }
}

static MRI_File *
CreateTempFile (MRI_Dataset *ds)
{
  char file_name[256];
  char *tmp_dir;
  MRI_File *f;
  char dir_name[256];

  if ((tmp_dir = getenv("MRI_TMP_DIR")) == NULL)
    tmp_dir = "/tmp";
  if (tmp_dir[0] != '/')
    {
      if (getcwd(dir_name, 256) == NULL)
	{
	  mri_report_error(ds, "libmri: cannot get current directory name\n");
	  abort();
	}
      sprintf(file_name, "%s/%s/mri%d.%d", dir_name, tmp_dir,
	      (int)getpid(), tmp_num++);
    }
  else
    sprintf(file_name, "%s/mri%d.%d", tmp_dir, (int)getpid(), tmp_num++);
  f = GetFile(ds, file_name);
  f->external = TRUE;
  return(f);
}

static void
DeallocateDataset (MRI_Dataset *ds)
{
  MRI_Chunk *ch, *nch;
  MRI_File *f, *nf;
  MRI_KeyValue *kv, *nkv;
  MRI_Buffer *b, *nb;
  int i;

  /* deallocate the files */
  f = ds->files;
  while (f != NULL)
    {
      nf = f->next;
      CloseFile(f);
      free(f->name);
      free(f);
      f = nf;
    }
  ds->files = NULL;

  /* deallocate the chunks */
  ch = ds->chunks;
  while (ch != NULL)
    {
      nch = ch->next;
      free(ch->dimensions);
      free(ch->actual_dimensions);
      free(ch);
      ch = nch;
    }
  ds->chunks = NULL;

  /* deallocate the buffers */
  b = ds->buffers;
  while (b != NULL)
    {
      nb = b->next;
      free(b->buffer);
      free(b);
      b = nb;
    }
  ds->buffers = NULL;
  b = ds->retained_buffers;
  while (b != NULL)
    {
      nb = b->next;
      free(b->buffer);
      free(b);
      b = nb;
    }
  ds->retained_buffers = NULL;
  
  /* deallocate the key/value pairs */
  for (i = 0; i < ds->hash_table_size; ++i)
    for (kv = ds->hash_table[i]; kv != NULL; )
      {
	nkv = kv->next_in_hash_table;
	free(kv->key);
	free(kv->value);
	free(kv);
	kv = nkv;
      }

  /* deallocate other fields */
  if (ds->iteration_table != NULL)
    free(ds->iteration_table);
  free(ds->name);
  free(ds->hash_table);

  /* deallocate the MRI_Dataset itself */
  free(ds);
}


static int
InitializeKeyValue (MRI_Dataset *ds, const char *key, const char *value)
{
  MRI_KeyValue *kv;

  /* check that the key has not already been read */
  kv = FindInHashTable(ds, key, FALSE);
  if (kv != NULL)
    {
      mri_report_error(ds, "mri_open_dataset: duplicate keys <%s> found in header\n", key);
      return(FALSE);
    }

  /* enter key/value into hash table */
  kv = FindInHashTable(ds, key, TRUE);
  kv->value = (char *) malloc(strlen(value) + 1);
  strcpy(kv->value, value);
  return(TRUE);
}

static int
CheckForHooks (MRI_Dataset *ds, MRI_KeyValue *kv)
{
  char *tail;
  MRI_Datatype new_datatype;
  MRI_Order new_order;
  long long new_offset;
  long long new_extent;
  long long new_size;
  int new_little_endian;
  MRI_Chunk *ch;
  char chunk_name[MRI_MAX_KEY_LENGTH+1];

  /* check if we are creating a chunk */
  if (strcmp(kv->value, "[chunk]") == 0)
    {
      /* check if there is already a chunk with that name */
      for (ch = ds->chunks; ch != NULL; ch = ch->next)
	if (strcmp(ch->name, kv->key) == 0)
	  return(TRUE);

      /* create a brand new chunk */
      ch = NewChunk(ds, kv->key);
      if (ch == NULL)
	return(FALSE);

      ch->actual_file = ds->header_file;
      ch->actual_offset = 0;
      ch->actual_size = 0;
      ModifyChunk(ch);

      if (!mri_has(ds, mri_cat(ch->name, ".little_endian")))
	mri_set_int(ds, mri_cat(ch->name, ".little_endian"),
	            !bio_big_endian_machine);

      return(TRUE);
    }

  /* if the key has no dots in it, there
     is nothing more to check for */
  if (strchr(kv->key, '.') == NULL)
    return(TRUE);

  /* split the key up into chunk name and field name (tail) */
  strcpy(chunk_name, kv->key);
  tail = strrchr(chunk_name, '.');
  *tail++ = '\0';

  /* hack: if the chunk_name ends in ".extent", consider
     that word part of the tail */
  if ((tail - chunk_name) >= 8 &&
      strcmp(tail-8, ".extent") == 0)
    {
      *(tail-1) = '.';
      *(tail-8) = '\0';
      tail -= 7;
    }
  
  /* check if the key has a corresponding chunk */
  for (ch = ds->chunks; ch != NULL; ch = ch->next)
    if (strcmp(ch->name, chunk_name) == 0)
      break;
  if (ch == NULL)
    /* no, it doesn't, so we are through */
    return(TRUE);

  /* handle a change in datatype */
  if (strcmp(tail, "datatype") == 0)
    {
      if (!ConvertDatatype(&new_datatype, kv->value))
	{
	  mri_report_error(ds, "libmri: Invalid datatype specified for chunk.\n");
	  return(FALSE);
	}
      if (new_datatype != ch->datatype)
	{
	  ch->datatype = new_datatype;
	  ModifyChunk(ch);
	}
      return(TRUE);
    }

  /* handle a change in dimensions */
  if (strcmp(tail, "dimensions") == 0)
    {
      if ((int) strlen(kv->value) > MRI_MAX_DIMS)
	{
	  mri_report_error(ds, "libmri: Too many dimensions specified.\n");
	  return(FALSE);
	}
      if (strlen(ch->dimensions) != strlen(kv->value))
	{
	  free(ch->dimensions);
	  ch->dimensions = (char *) malloc(strlen(kv->value)+1);
	}
      strcpy(ch->dimensions, kv->value);
      ModifyChunk(ch);
      return(TRUE);
    }

  /* handle a change in file location */
  if (strcmp(tail, "file") == 0)
    {
      ch->file = GetFile(ds, kv->value);
      ModifyChunk(ch);
      return(TRUE);
    }

  /* handle a change in order within file */
  if (strcmp(tail, "order") == 0)
    {
      if (strcmp(kv->value, "fixed_offset") == 0)
	new_order = MRI_FIXED_OFFSET;
      else if (strcmp(kv->value, "external") == 0)
	new_order = MRI_EXTERNAL;
      else if (sscanf(kv->value, "%d", &new_order) != 1)
	{
	  mri_report_error(ds, "libmri: Invalid chunk order specified.\n");
	  return(FALSE);
	}
      if (new_order != ch->order)
	{
	  ch->order = new_order;
	  ModifyChunk(ch);
	}
      return(TRUE);
    }

  /* handle a change in specified offset */
  if (strcmp(tail, "offset") == 0)
    {
      if (sscanf(kv->value, "%lld", &new_offset) != 1)
	{
	  mri_report_error(ds, "libmri: Invalid chunk offset.\n");
	  return(FALSE);
	}
      if (new_offset != ch->offset)
	{
	  ch->offset = new_offset;
	  ModifyChunk(ch);
	}
      return(TRUE);
    }

  /* handle a change in extent */
  if (strncmp(tail, "extent.", 7) == 0 &&
      strlen(tail) == 8)
    {
      if (sscanf(kv->value, "%lld", &new_extent) != 1)
	{
	  mri_report_error(ds, "libmri: Invalid extent.\n");
	  return(FALSE);
	}
      ModifyChunk(ch);
      return(TRUE);
    }

  /* handle a change in endianness */
  if (strcmp(tail, "little_endian") == 0)
    {
      if (sscanf(kv->value, "%d", &new_little_endian) != 1 ||
	  (new_little_endian != 0 && new_little_endian != 1))
	{
	  mri_report_error(ds, "libmri: Invalid little_endian value.\n");
	  return(FALSE);
	}
      if (new_little_endian != ch->little_endian)
	{
	  ch->little_endian = new_little_endian;
	  ModifyChunk(ch);
	}
      return(TRUE);
    }

  /* handle an attempt to change the size */
  if (strcmp(tail, "size") == 0)
    {
      if (sscanf(kv->value, "%lld", &new_size) != 1 ||
	  new_size != ch->size)
	{
	  mri_report_error(ds, "libmri: Not allowed to set chunk size.\n");
	  return(FALSE);
	}
      return(TRUE);
    }

  /* anything else is OK */
  return(TRUE);
}

static void
ModifyChunk (MRI_Chunk *ch)
{
  int i;
  long long size;
  char fieldname[16];
  char key_name[MRI_MAX_KEY_LENGTH+1];

  /* update the extents array */
  size = MRI_TypeLength(ch->datatype);
  for (i = 0; i < (int) strlen(ch->dimensions); ++i)
    {
      sprintf(fieldname, ".extent.%c", ch->dimensions[i]);
      if (!mri_has(ch->ds, mri_cat(ch->name, fieldname)))
	ch->extent[i] = 1;
      else
	ch->extent[i] = mri_get_int(ch->ds, mri_cat(ch->name, fieldname));
      size *= ch->extent[i];
    }

  /* update the size field */
  ch->size = size;
  sprintf(key_name, "%s.size", ch->name);
  mri_set_int(ch->ds, key_name, ch->size);

  ch->ds->recompute_positions = TRUE;
  ch->modified = TRUE;
  SetChunkNotReady(ch);
  CheckForStdImages(ch->ds);
}

static void
SetChunkNotReady (MRI_Chunk *ch)
{
  ch->ready_to_read = FALSE;
  ch->ready_to_write = FALSE;
}

static void
CheckForStdImages (MRI_Dataset *ds)
{
  MRI_Chunk *ch;

  for (ch = ds->chunks; ch != NULL; ch = ch->next)
    if (strcmp(ch->name, "images") == 0)
      {
	if (strcmp(ch->dimensions, "xyzt") == 0 ||
	    strcmp(ch->dimensions, "xyz") == 0)
	  {
	    ds->std_image_vector_size = 1;
	    ds->std_image_size = ch->extent[0] * ch->extent[1];
	    ds->std_n_slices = ch->extent[2];
	    ds->std_images = TRUE;
	    return;
	  }
	if (strcmp(ch->dimensions, "vxyzt") == 0 ||
	    strcmp(ch->dimensions, "vxyz") == 0)
	  {
	    ds->std_image_vector_size = ch->extent[0];
	    ds->std_image_size = ch->extent[0] * ch->extent[1] * ch->extent[2];
	    ds->std_n_slices = ch->extent[3];
	    ds->std_images = TRUE;
	    return;
	  }
      }
  ds->std_images = FALSE;
}

static int
CheckForRemovalHooks (MRI_Dataset *ds, MRI_KeyValue *kv)
{
  char *tail;
  MRI_Datatype new_datatype;
  MRI_Order new_order;
  long long new_offset;
  int new_little_endian;
  MRI_Chunk *ch;
  char chunk_name[MRI_MAX_KEY_LENGTH+1];

  /* check if we are discarding a chunk */
  if (strcmp(kv->value, "[chunk]") == 0)
    {
      /* check if there is already a chunk with that name */
      for (ch = ds->chunks; ch != NULL; ch = ch->next)
	if (strcmp(ch->name, kv->key) == 0)
	  {
	    DeallocateChunk(ch);
	    break;
	  }
      return(TRUE);
    }

  /* if the key has no dots in it, there
     is nothing more to check for */
  if (strchr(kv->key, '.') == NULL)
    return(TRUE);

  /* split the key up into chunk name and field name (tail) */
  strcpy(chunk_name, kv->key);
  tail = strrchr(chunk_name, '.');
  *tail++ = '\0';

  /* hack: if the chunk_name ends in ".extent", consider
     that word part of the tail */
  if ((tail - chunk_name) >= 8 &&
      strcmp(tail-8, ".extent") == 0)
    {
      *(tail-1) = '.';
      *(tail-8) = '\0';
      tail -= 7;
    }
  
  /* check if the key has a corresponding chunk */
  for (ch = ds->chunks; ch != NULL; ch = ch->next)
    if (strcmp(ch->name, chunk_name) == 0)
      break;
  if (ch == NULL)
    /* no, it doesn't, so we are through */
    return(TRUE);

  /* reset datatype */
  if (strcmp(tail, "datatype") == 0)
    {
      new_datatype = MRI_INT16;
      if (new_datatype != ch->datatype)
	{
	  ch->datatype = new_datatype;
	  ModifyChunk(ch);
	}
      return(TRUE);
    }

  /* reset dimensions */
  if (strcmp(tail, "dimensions") == 0)
    {
      free(ch->dimensions);
      ch->dimensions = (char *) malloc(5);
      strcpy(ch->dimensions, "xyzt");
      ModifyChunk(ch);
      return(TRUE);
    }

  /* reset file location */
  if (strcmp(tail, "file") == 0)
    {
      ch->file = ds->header_file;
      ModifyChunk(ch);
    }

  /* reset order within file */
  if (strcmp(tail, "order") == 0)
    {
      new_order = 0;
      if (new_order != ch->order)
	{
	  ch->order = new_order;
	  ModifyChunk(ch);
	}
      return(TRUE);
    }

  /* reset offset */
  if (strcmp(tail, "offset") == 0)
    {
      new_offset = 0;
      if (new_offset != ch->offset)
	{
	  ch->offset = new_offset;
	  ModifyChunk(ch);
	}
      return(TRUE);
    }

  /* reset an extent */
  if (strncmp(tail, "extent.", 7) == 0 &&
      strlen(tail) == 8)
    {
      ModifyChunk(ch);
      return(TRUE);
    }

  /* reset the endianness */
  if (strcmp(tail, "little_endian") == 0)
    {
      new_little_endian = FALSE;
      if (new_little_endian != ch->little_endian)
	{
	  ch->little_endian = new_little_endian;
	  ModifyChunk(ch);
	}
      return(TRUE);
    }

  /* handle an attempt to remove the size */
  if (strcmp(tail, "size") == 0)
    {
      mri_report_error(ds, "libmri: Not allowed to remove chunk size.\n");
      return(FALSE);
    }

  /* anything else is OK */
  return(TRUE);
}

static void
ChunkSeek (MRI_Chunk *ch, long long offset)
{
  if (mri_fseek(ch->file->fp, ch->offset + offset, SEEK_SET) != 0)
    mri_report_error(ch->ds,
		     "libmri: could not seek to position %lld in file %s\n",
		     ch->offset + offset,
		     ch->file->name);
}

static void
DeallocateChunk (MRI_Chunk *chunk)
{
  MRI_Dataset *ds;
  MRI_Chunk *pch, *ch;

  SetChunkNotReady(chunk);
  ds = chunk->ds;
  for (pch = NULL, ch = ds->chunks; ch != NULL; pch = ch, ch = ch->next)
    if (ch == chunk)
      break;
  if (ch == NULL)
    {
      mri_report_error(ds, "libmri: internal error in deallocating chunk\n");
      abort();
    }
  if (pch != NULL)
    pch->next = ch->next;
  else
    ds->chunks = ch->next;
  free(ch->dimensions);
  free(ch->actual_dimensions);
  free(ch);
  ds->recompute_positions = TRUE;
  CheckForStdImages(ds);
}

static void
mri_report_error (MRI_Dataset *ds, char *fmt, ...)
{
  va_list args;

  mri_error = fmt;
  if (error_handling != MRI_IGNORE_ERRORS)
    {
      va_start(args, fmt);
      vfprintf(stderr, fmt, args);
      va_end(args);
    }
  if (error_handling == MRI_ABORT_ON_ERROR)
    {
      char buf[256];
      time_t t;
      t= time(NULL);
      sprintf(buf,"%s",asctime(localtime(&t)));
      buf[strlen(buf)-1]= '\0'; /* supress trailing blank */
      if (ds) {
	strcat(buf,": fatal error on ");
	strncat(buf,ds->name,sizeof(buf)-1);
	buf[sizeof(buf)-1]= '\0';
      }
      else {
	sprintf(buf,": fatal error on unknown dataset");
      }
      if (errno != 0)
	perror(buf);
      else fprintf(stderr,"%s\n",buf);
      fprintf(stderr, "Aborting...\n");
      abort();
    }
}

static void
mri_report_warning (MRI_Dataset *ds, char *fmt, ...)
{
  va_list args;

  mri_error = fmt;
  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
}


static int
ReadHeader (MRI_Dataset *ds, FILE *f)
{
  int c;

  for (;;)
    {
      /* skip white space */
      while (isspace(c = fgetc(f))) ;
      
      /* have we reached the end of the file? */
      if (c == EOF || c == '\032')
	{
	  ds->header_size = mri_ftell(f);
	  return(TRUE);
	}

      ungetc(c, f);
      if (!ReadKeyValuePair(ds, f))
	return(FALSE);
    }
}

static int
ReadKeyValuePair (MRI_Dataset *ds, FILE *f)
{
  int c;
  char key[MRI_MAX_KEY_LENGTH+1];
  char value[MRI_MAX_VALUE_LENGTH+1];

  /* read the key */
  if (!ReadString(ds, key, MRI_MAX_KEY_LENGTH, f))
    return(FALSE);

  /* skip white space */
  while (isspace(c = fgetc(f))) ;

  /* check for the '=' sign */
  if (c != '=')
    {
      mri_report_error(ds, "mri_open_dataset: = not found in key/value pair");
      return(FALSE);
    }

  /* read in the value */
  if (!ReadString(ds, value, MRI_MAX_VALUE_LENGTH, f))
    return(FALSE);

  if (!InitializeKeyValue(ds, key, value))
    return(FALSE);
  return(TRUE);
}

static int
ReadString (MRI_Dataset *ds, char *s, int max_len, FILE *f)
{
  int c;
  int len = 0;

  /* skip leading white space */
  while (isspace(c = fgetc(f)) && (c != '\n')) ;

  if (c == '"')
    return(ReadQuotedString(ds, s, max_len, f));
  else
    ungetc(c, f);

  /* read characters until we hit a character not in the set
     {32-126, but including tab and excluding '='} */
  while ((c = fgetc(f)) >= 32 && c <= 126 && c != '=' || c == '\t')
    {
      if (len >= max_len)
	{
	  s[len] = '\0';
	  mri_report_error(ds, "mri_open_dataset: string too long\n");
	  return(FALSE);
	}
      s[len++] = c;
    }
  --len;

  /* strip off trailing whitespace */
  while (len >= 0 && isspace(s[len]))
    --len;
  s[len+1] = '\0';

  if (c != EOF)
    ungetc(c, f);
  return(TRUE);
}

static int
ReadQuotedString (MRI_Dataset *ds, char *s, int max_len, FILE *f)
{
  int c;
  int len = 0;

  while ((c = fgetc(f)) != EOF && c != '"')
    {
      if (len >= max_len)
	{
	  s[len] = '\0';
	  mri_report_error(ds, "mri_open_dataset: string too long\n");
	  return(FALSE);
	}
      if (c == '\\')
	c = ReadBackslashedChar(f);
      if (c != 0)
	s[len++] = c;
    }
  s[len] = '\0';
  if (c == EOF)
    {
      mri_report_error(ds, "mri_open_dataset: EOF reached inside quoted string\n");
      return(FALSE);
    }
  return(TRUE);
}

static int
ReadBackslashedChar (FILE *f)
{
  int c;
  int v;

  switch (c = fgetc(f))
    {
    case 'n':
      return('\n');
    case 'r':
      return('\r');
    case 't':
      return('\t');
    case '0':
    case '1':
    case '2':
    case '3':
      v = (c - '0') << 6;
      c = fgetc(f);
      if (c < '0' || c > '7')
	return(0);
      v |= (c - '0') << 3;
      c = fgetc(f);
      if (c < '0' || c > '7')
	return(0);
      v |= (c - '0');
      return(v);
    case '\n':
    case EOF:
      return(0);
    default:
      return(c);
    }
}

static MRI_KeyValue *
FindInHashTable (MRI_Dataset *ds, const char *key,
		 int add)
{
  int hv;
  MRI_KeyValue *kv, *nkv;
  MRI_KeyValue **new;
  int i;

  hv = HashFunction(key);
  if (ds->hash_table != NULL)
    {
      hv &= ds->hash_table_size - 1;
      kv = ds->hash_table[hv];
      while (kv != NULL)
	{
	  if (strcmp(kv->key, key) == 0)
	    return(kv);
	  kv = kv->next_in_hash_table;
	}
    }

  /* not found */
  if (!add)
    return(NULL);

  /* add to hash table */
  if (++ds->n_keys > 2*ds->hash_table_size)
    {
      int new_size= 4*ds->hash_table_size;
      /* quadruple the hash table */
      new = (MRI_KeyValue **) malloc(new_size * sizeof(MRI_KeyValue *));
      for (i = 0; i < new_size; ++i) new[i]= NULL;
      for (i = 0; i < ds->hash_table_size; ++i)
	{
	  kv = ds->hash_table[i];
	  while (kv != NULL)
	    {
	      nkv = kv->next_in_hash_table;
	      hv = HashFunction(kv->key) & (4 * ds->hash_table_size - 1);
	      kv->next_in_hash_table = new[hv];
	      new[hv] = kv;
	      kv = nkv;
	    }
	}
      free(ds->hash_table);
      ds->hash_table = new;
      ds->hash_table_size = new_size;
      hv = HashFunction(key) & (ds->hash_table_size - 1);
    }

  kv = (MRI_KeyValue *) malloc(sizeof(MRI_KeyValue));
  kv->key = (char *) malloc(strlen(key) + 1);
  strcpy(kv->key, key);
  kv->value = NULL;
  kv->next_in_hash_table = ds->hash_table[hv];
  ds->hash_table[hv] = kv;
  return(kv);
}

static void
RemoveFromHashTable (MRI_Dataset *ds, const char *key)
{
  int hv;
  MRI_KeyValue *pkv, *kv;

  hv = HashFunction(key);
  if (ds->hash_table == NULL)
    return;
  hv &= ds->hash_table_size - 1;
  kv = ds->hash_table[hv];
  pkv = NULL;
  while (kv != NULL)
    {
      if (strcmp(kv->key, key) == 0)
	{
	  if (pkv != NULL)
	    pkv->next_in_hash_table = kv->next_in_hash_table;
	  else
	    ds->hash_table[hv] = kv->next_in_hash_table;
	  free(kv->key);
	  if (kv->value != NULL)
	    free(kv->value);
	  free(kv);
	  --ds->n_keys;
	  return;
	}
      pkv = kv;
      kv = kv->next_in_hash_table;
    }
}

static int
HashFunction (const char *s)
{
  int hv = 0;

  while (*s != '\0')
    hv = ((237 * hv) + *s++) & 0x7fffff;
  return(hv);
}

static void
WriteHeader (MRI_Dataset *ds, int separator, FILE *f)
{
  MRI_KeyValue *kv;
  char *key_name;
  
  mri_iterate_over_keys(ds);

  while ((key_name = mri_next_key(ds)) != NULL)
    {
      kv = FindInHashTable(ds, key_name, FALSE);
      if (kv == NULL)
	{
	  mri_report_error(ds, "libmri: internal error - write header iteration\n");
	  abort();
	}
      WriteKeyValuePair(kv, f);
    }

  if (separator)
    {
      /* we have other chunks going into the header file,
	 so write a couple separator characters */
      fputc('\014', f);
      fputc('\032', f);
    }
}

static int
CompareKeys (const void *k1, const void *k2)
{
  return(strcmp((*((MRI_KeyValue **) k1))->key,
		(*((MRI_KeyValue **) k2))->key));
}

static void
WriteKeyValuePair (MRI_KeyValue *kv, FILE *f)
{
  WriteString(kv->key, f);
  fprintf(f, " = ");
  WriteString(kv->value, f);
  fputc('\012', f);
}

static void
WriteString (char *s, FILE *f)
{
  int i;

  /* check if we have to quote the string */
  if (strlen(s)==0) { 
    WriteQuotedString(s,f);
    return;
  }
  for (i = 0; s[i] != '\0'; ++i) {
    if (s[i] < 32 ||
	s[i] > 126 ||
	s[i] == ' ' ||
	s[i] == '=')
      {
	WriteQuotedString(s, f);
	return;
      }
  }
  if (strlen(s)==0) {
    WriteQuotedString(s, f);
    return;
  }
  fputs(s, f);
}

static void
WriteQuotedString (char *s, FILE *f)
{
  char c;

  fputc('"', f);
  while ((c = *s++) != '\0')
    if (c == '"' || c == '\\')
      {
	fputc('\\', f);
	fputc(c, f);
      }
    else if (c >= 32 && c <= 126)
      fputc(c, f);
    else
      {
	fputc('\\', f);
	if (c == '\n')
	  fputc('n', f);
	else
	  fprintf(f, "%.3o", c);
      }
  fputc('"', f);
}

static void
ComputeChunkPositions (MRI_Dataset *ds)
{
  int n_chunks;
  MRI_Chunk *ch, *nch;
  MRI_Chunk *chunks[MRI_MAX_CHUNKS];
  int n, m;
  EmptyBlock empty_blocks[MRI_MAX_CHUNKS+2];
  int n_empty_blocks;
  int first;
  char key_name[MRI_MAX_KEY_LENGTH+1];

  /* mark all chunks as unchecked; also mark files
     as external */
  for (ch = ds->chunks; ch != NULL; ch = ch->next)
    {
      if (ch->order == MRI_EXTERNAL)
	ch->file->external = TRUE;
      ch->checked = FALSE;
    }

  for (ch = ds->chunks; ch != NULL; ch = ch->next)
    {
      /* if we have handled this chunk already or
	 if it is external, then skip it */
      if (ch->checked || ch->order == MRI_EXTERNAL)
	continue;

      /* find all other chunks destined for the same file */
      n_chunks = 0;
      for (nch = ch; nch != NULL; nch = nch->next)
	if (nch->file == ch->file)
	  chunks[n_chunks++] = nch;

      /* set up the empty block array */
      if (ch->file == ds->header_file)
	empty_blocks[0].start = ds->header_size;
      else
	empty_blocks[0].start = 0;
      empty_blocks[0].end = MAX_CHUNK_SIZE-1;
      n_empty_blocks = 1;

      /* reserve space for the fixed offset chunks */
      for (n = 0; n < n_chunks; ++n)
	if (chunks[n]->order == MRI_FIXED_OFFSET)
	  {
	    (void) ReserveBlock(empty_blocks, &n_empty_blocks, chunks[n]->offset, chunks[n]->size);
	    /* remove from table */
	    for (m = n+1; m < n_chunks; ++m)
	      chunks[m-1] = chunks[m];
	    --n_chunks;
	    --n;
	  }

      /* sort them according to order and offset (offset is
         used as a secondary sorting key since it will produce
	 some hysteresis and avoid unnecessary swaps among chunks
	 with equal orders) */
      for (n = 0; n < n_chunks; ++n)
	{
	  first = n;
	  for (m = first+1; m < n_chunks; ++m)
	    if (chunks[n]->order < chunks[first]->order ||
		chunks[n]->order == chunks[first]->order &&
		chunks[n]->offset < chunks[first]->offset)
	      first = m;
	  nch = chunks[n];
	  chunks[n] = chunks[first];
	  chunks[first] = nch;
	}

      /* go through and allocate space for them */
      for (n = 0; n < n_chunks; ++n)
	{
	  chunks[n]->offset = ReserveBlock(empty_blocks, &n_empty_blocks,
					   -1LL, chunks[n]->size);
	  sprintf(key_name, "%s.offset", chunks[n]->name);
	  mri_set_int(ds, key_name, chunks[n]->offset);
	  if (chunks[n]->file != chunks[n]->actual_file ||
	      chunks[n]->offset != chunks[n]->actual_offset)
	    chunks[n]->modified = TRUE;
	  chunks[n]->checked = TRUE;
	}
    }
  ds->recompute_positions = FALSE;
}

static long long
ReserveBlock (EmptyBlock *empty_blocks, int *n_empty_blocks,
	      long long offset, long long size)
{
  int i;
  int m;
  long long start;
  EmptyBlock blk;

  if (offset < 0)
    {
      /* determine the offset */
      for (i = 0; i < *n_empty_blocks; ++i)
	if (size <= (empty_blocks[i].end - empty_blocks[i].start + 1))
	  {
	    /* a potential location; now check if we can do the
	       proper alignment */
	    if (size >= MRI_ALIGNMENT_THRESHOLD)
	      {
		start = ((empty_blocks[i].start + (MRI_ALIGNMENT_BOUNDARY-1)) /
		  MRI_ALIGNMENT_BOUNDARY) * MRI_ALIGNMENT_BOUNDARY;
		if (size > (empty_blocks[i].end - start + 1))
		  continue;
		offset = start;
	      }
	    else
	      offset = empty_blocks[i].start;
	    break;
	  }
      if (i >= *n_empty_blocks)
	{
	  fprintf(stderr, "libmri: cannot allocate filespace for chunk\n");
	  abort();
	}
    }

  /* reserve the space */
  for (i = 0; i < *n_empty_blocks; ++i)
    if (offset <= empty_blocks[i].end &&
	offset + size - 1 >= empty_blocks[i].start)
      {
	/* pull out the empty block */
	blk = empty_blocks[i];
	for (m = i+1; m < *n_empty_blocks; ++m)
	  empty_blocks[m - 1] = empty_blocks[m];
	--(*n_empty_blocks);

	if (blk.start < offset)
	  {
	    /* add the partial block */
	    for (m = *n_empty_blocks; m > i; --m)
	      empty_blocks[m] = empty_blocks[m - 1];
	    empty_blocks[i].start = blk.start;
	    empty_blocks[i].end = offset - 1;
	    ++i;
	    ++(*n_empty_blocks);
	  }

	if (blk.end > offset + size - 1)
	  {
	    /* add the partial block */
	    for (m = *n_empty_blocks; m > i; --m)
	      empty_blocks[m] = empty_blocks[m - 1];
	    empty_blocks[i].start = offset + size;
	    empty_blocks[i].end = blk.end;
	    ++i;
	    ++(*n_empty_blocks);
	  }
	/* compensate for removing one block */
	--i;
      }

  /* return the chosen offset */
  return(offset);
}

static void
MakeChunkFilename (char *name, char *ds_file, char *ch_file)
{
  char *p;

  if (ch_file[0] == '\0')
    {
      /* place the chunk in header's file */
      strcpy(name, ds_file);
    }
  else if (ch_file[0] == '.')
    {
      /* strip any extension off the ds_file and append the new extension to it */
      strcpy(name, ds_file);
      if ((p = strrchr(name, '.')) != NULL)
	*p = '\0';
      strcat(name, ch_file);
    }
  else if (strchr(ch_file, '/') == NULL && strchr(ds_file, '/') != NULL)
    {
      strcpy(name, ds_file);
      p = strrchr(name, '/');
      *(p+1) = '\0';
      strcat(name, ch_file);
    }
  else
    strcpy(name, ch_file);
}

static int
PrepareToWrite (MRI_Chunk *ch)
{
  if (ch->ds->mode == MRI_READ)
    {
      mri_report_error(ch->ds, "libmri: attempt to write to read-only dataset\n");
      return(FALSE);
    }

  if (ch->order == MRI_EXTERNAL)
    {
      mri_report_error(ch->ds, "libmri: attempt to write to external chunk\n");
      return(FALSE);
    }

  if (ch->ds->mode != MRI_MODIFY_DATA)
    {
      if (ch->ds->recompute_positions)
	ComputeChunkPositions(ch->ds);
      
      if (ch->modified)
	RepositionChunk(ch, NULL);
    }
  else
    {
      if (ch->ds->recompute_positions ||
	  ch->modified)
	{
	  fprintf(stderr, "libmri: Chunk moved despite MRI_MODIFY_DATA mode!\n");
	  abort();
	}
    }

  if (!OpenFile(ch->file, TRUE))
    return(FALSE);

  ch->ready_to_read = TRUE;
  ch->ready_to_write = TRUE;
  return(TRUE);
}

static int
PrepareToRead (MRI_Chunk *ch)
{
  if (ch->ds->recompute_positions)
    ComputeChunkPositions(ch->ds);

  if (ch->modified)
    RepositionChunk(ch, NULL);

  if (!OpenFile(ch->file, FALSE))
    return(FALSE);

  ch->ready_to_read = TRUE;
  return(TRUE);
}

static int
OpenFile (MRI_File *file, int for_writing)
{
  MRI_Dataset *ds;
  MRI_File *f;
  MRI_File *oldest;

  ds = file->ds;

  if (for_writing && ds->mode == MRI_READ)
    {
      mri_report_error(ds, "libmri: Attempt to write file in read-only dataset\n", file->name);
      return(FALSE);
    }

  /* check if the file is already open */
  if (file->fp != NULL)
    {
      /* do we already have it open in the right mode? */
      if (file->writeable || !for_writing) {
	file->last_use = file_access++;
	return(TRUE);
      }

      /* no, try to reopen it for writing */
#ifdef DEBUG
      Log("Closing file %s with intent to reopen for writing\n", file->name);
#endif
      fclose(file->fp);
      file->fp = NULL;
    }

  if (ds->n_open_files >= MRI_MAX_OPEN_FILES)
    {
      /* close the least recently used file */
      oldest = NULL;
      for (f = ds->files; f != NULL; f = f->next)
	if (f->fp != NULL &&
	    (oldest == NULL ||
	     f->last_use < oldest->last_use))
	  oldest = f;
      if (oldest == NULL)
	{
	  mri_report_error(ds, "libmri: internal error in OpenFile\n");
	  abort();
	}
      CloseFile(oldest);
    }

  if (!for_writing)
    {
#ifdef DEBUG
      Log("Opening file %s for reading\n", file->name);
#endif
      file->fp = fopen(file->name, "r");
    }
  else if (ds->mode == MRI_WRITE && file->last_use == 0)
    {
#ifdef DEBUG
      Log("Opening file %s for writing (truncating)\n", file->name);
#endif
      /* only truncate if we have not opened this file before */
      file->fp = fopen(file->name, "w+");
    }
  else
    {
#ifdef DEBUG
      Log("Opening file %s for update\n", file->name);
#endif
      file->fp = fopen(file->name, "r+");
      if (file->fp == NULL)
	{
#ifdef DEBUG
	  Log("Could not open file for update, so creating for update\n");
#endif
	  file->fp = fopen(file->name, "w+");
	}
    }
  if (file->fp == NULL)
    {
      mri_report_error(ds, "libmri: Unable to open file %s\n", file->name);
      return(FALSE);
    }
  ++ds->n_open_files;
  file->last_use = file_access++;
  file->writeable = for_writing;
  return(TRUE);
}

static void
CloseFile (MRI_File *file)
{
  MRI_Chunk *ch;

  if (file->fp != NULL)
    {
      /* make sure we don't try to read or write any chunks
	 in that file without reopening the file pointer */
      for (ch = file->ds->chunks; ch != NULL; ch = ch->next)
	if (ch->file == file)
	  SetChunkNotReady(ch);

#ifdef DEBUG
      Log("Closing file %s\n", file->name);
#endif
      fclose(file->fp);
      file->fp = NULL;
      --file->ds->n_open_files;
    }
}

static void
RepositionChunk (MRI_Chunk *ch, CopyRequest **copy_queue)
{
  CopyRequest *queue, *req, *nreq;
  int use_temp_file;
  MRI_File *temp;
  MRI_Chunk *och;
  MRI_File *pf, *f, *nf;

  /* we want to have the chunk be in file ch->file
     at bytes ch->offset through (ch->offset + ch->size - 1);
     it is currently in file ch->actual_file at bytes
     ch->actual_offset through (ch->actual_offset + ch->actual_size - 1) */

  if (ch->order == MRI_EXTERNAL)
    {
      /* we don't have to do anything for chunks
	 located in external files except update
	 their attributes */
      UpdateChunkAttributes(ch);
      return;
    }

  ch->repositioning = TRUE;
  use_temp_file = FALSE;
  queue = NULL;

  /* go through other chunks and see which ones have to be moved */
  for (och = ch->ds->chunks; och != NULL; och = och->next)
    if (ch->file == och->actual_file &&
	ch->offset <= och->actual_offset + och->actual_size + 1 &&
	ch->offset + ch->size + 1 >= och->actual_offset)
      /* there is an overlap */
      if (!och->repositioning)
	/* move the chunk */
	RepositionChunk(och, (copy_queue != NULL) ? copy_queue : &queue);
      else
	if (och != ch)
	  /* we have hit a circularity, so we'll have to use
	     temporary files to resolve the impasse */
	  use_temp_file = TRUE;
	else
	  {
	    /* we are overwriting the chunk itself; check if
	       we can do the copy in-place */
	    if (ch->offset != ch->actual_offset ||
		MRI_TypeLength(ch->datatype) > MRI_TypeLength(ch->actual_datatype))
	      use_temp_file = TRUE;
	  }

  if (use_temp_file)
    {
      temp = CreateTempFile(ch->ds);
      ConvertChunk(ch, temp, 0);
      /* queue for copying later */
      req = (CopyRequest *) malloc(sizeof(CopyRequest));
      req->src_file = temp;
      req->src_offset = 0;
      req->dest_file = ch->file;
      req->dest_offset = ch->offset;
      req->size = ch->size;
      req->chunk = ch;
      if (copy_queue != NULL)
	{
	  req->next = *copy_queue;
	  *copy_queue = req;
	}
      else
	{
	  req->next = queue;
	  queue = req;
	}
    }
  else
    {
      ConvertChunk(ch, ch->file, ch->offset);
      UpdateChunkAttributes(ch);
    }

  /* check if we are back at the top-level call
     to RepositionChunk, and it is time to do
     all the queued copy operations */
  if (copy_queue == NULL)
    {
      for (req = queue; req != NULL; req = nreq)
	{
	  nreq = req->next;
	  CopyBlock(req->dest_file, req->dest_offset, req->src_file, req->src_offset, req->size);
	  /* remove the temporary file holding the chunk */
	  UpdateChunkAttributes(req->chunk);
	  DestroyFile(req->src_file);
	  free(req);
	}

      /* go through the list of files and see if any are unused */
      for (f = ch->ds->files; f != NULL; f = f->next)
	f->used = FALSE;
      ch->ds->header_file->used = TRUE;
      for (och = ch->ds->chunks; och != NULL; och = och->next)
	och->file->used = TRUE;
      pf = NULL;
      f = ch->ds->files;
      while (f != NULL)
	{
	  nf = f->next;
	  if (!f->used && !f->external)
	    {
	      CloseFile(f);
#ifdef DEBUG
	      Log("GOING TO UNLINK FILE!!! (within RepositionChunk): %s\n",
		  f->name);
#endif
	      unlink(f->name);
	      free(f->name);
	      free(f);
	      if (pf != NULL)
		pf->next = nf;
	      else 
		ch->ds->files = nf;
	    }
	  else
	    pf = f;
	  f = nf;
	}
    }
}

static void
UpdateChunkAttributes(MRI_Chunk *ch)
{
  int i;

  ch->actual_file = ch->file;
  ch->actual_datatype = ch->datatype;
  if (strlen(ch->actual_dimensions) != strlen(ch->dimensions))
    {
      free(ch->actual_dimensions);
      ch->actual_dimensions = (char *) malloc(strlen(ch->dimensions) + 1);
    }
  strcpy(ch->actual_dimensions, ch->dimensions);
  for (i = 0; i < (int) strlen(ch->dimensions); ++i)
    ch->actual_extent[i] = ch->extent[i];
  ch->actual_little_endian = ch->little_endian;
  ch->actual_offset = ch->offset;
  ch->actual_size = ch->size;
  ch->modified = FALSE;
  ch->repositioning = FALSE;
}

static void
CopyBlock (MRI_File *dest_file, long long dest_offset,
	   MRI_File *src_file, long long src_offset,
	   long long size)
{
  long long n;
  long long offset;
  char buffer[BUFFER_SIZE];

  /* check if this is a null operation */
  if (src_file == dest_file &&
      src_offset == dest_offset)
    return;

  if (!OpenFile(src_file, FALSE))
    {
      fprintf(stderr, "libmri: Internal copying error (file %s)\n",
	      src_file->name);
      abort();
    }
  if (!OpenFile(dest_file, TRUE))
    {
      fprintf(stderr, "libmri: Internal copying error (file %s)\n",
	      dest_file->name);
      abort();
    }
  offset = 0;
  while (size > 0)
    {
      n = (size > BUFFER_SIZE) ? BUFFER_SIZE : size;
      if (mri_fseek(src_file->fp, src_offset+offset, SEEK_SET) != 0)
	{
	  fprintf(stderr, "libmri: could not seek to position %lld in file %s\n",
		  src_offset+offset, src_file->name);
	  abort();
	}
      if (fread(buffer, (size_t) n, 1, src_file->fp) != 1)
	{
	  fprintf(stderr, "libmri: attempt to read past EOF in %s\n",
			   src_file->name);
	  abort();
	}
      if (mri_fseek(dest_file->fp, dest_offset+offset, SEEK_SET) != 0)
	{
	  fprintf(stderr, "libmri: could not seek to position %lld in file %s\n",
		  dest_offset+offset, dest_file->name);
	  abort();
	}
      if (fwrite(buffer, (size_t) n, 1, dest_file->fp) != 1)
	{
	  fprintf(stderr, "libmri: could not copy to file %s\n",
		  dest_file->name);
	  abort();
	}
      size -= n;
      offset += n;
    }
}

static void
ConvertChunk (MRI_Chunk *ch, MRI_File *f, long long offset)
{
#define N_READ	(BUFFER_SIZE / sizeof(double))
  long long read_size, write_size;
  long long i;
  long long n;
    long long read_offset, write_offset;
    long long count;
    double *dbl;
 
   /* find out max element size */
   union element {
     unsigned char     u;
     short             s;
     int               i;
     long long         l;
     float             f;
     double            d;
   } element;
 
   union pbuffer {
     unsigned char     *u;
     short             *s;
     int               *i;
     long long         *ll;
     float             *f;
     double            *d;
     void      *me;
   } in_buffer, out_buffer;
 
   int saved_bio_big_endian_input;
   int saved_bio_big_endian_output;
   int saved_bio_error;
 
   in_buffer.me = calloc (N_READ, sizeof(element));
   out_buffer.me = calloc (N_READ, sizeof(element));
 
   /* check if a copy operation will suffice */
   if (ch->datatype == ch->actual_datatype &&
       ch->little_endian == ch->actual_little_endian)
    {
      CopyBlock(f, offset, ch->actual_file, ch->actual_offset,
		(ch->size < ch->actual_size) ? ch->size : ch->actual_size);
      n = ch->size - ch->actual_size;
      if (n > 0)
	ClearBlock(f, offset + ch->actual_size, n);
      return;
    }

  /* we have to do some data conversion at this point */
  read_size = MRI_TypeLength(ch->actual_datatype);
  write_size = MRI_TypeLength(ch->datatype);
  saved_bio_big_endian_input = bio_big_endian_input;
  saved_bio_big_endian_output = bio_big_endian_output;
  saved_bio_error = bio_error;
  bio_big_endian_input = !ch->actual_little_endian;
  bio_big_endian_output = !ch->little_endian;

  if (!OpenFile(ch->actual_file, FALSE) ||
      !OpenFile(f, TRUE))
    abort();
      
  dbl = &(out_buffer.d[0]);
  read_offset = 0;
  write_offset = 0;
  count = ch->actual_size / read_size;
  if (count * write_size > ch->size)
    count = ch->size / write_size;
  while (count > 0)
    {
      n = (count > N_READ) ? N_READ : count;
      if (mri_fseek(ch->actual_file->fp,
		ch->actual_offset + read_offset,
		SEEK_SET) != 0)
	mri_report_error(ch->ds,
			 "libmri: could not seek to position %ld in file %s\n",
			 (long) (ch->actual_offset + read_offset),
			 ch->actual_file->name);
      bio_error = FALSE;
      switch (ch->actual_datatype)
	{
	case MRI_UINT8:
	  FRdUInt8Array(ch->actual_file->fp, in_buffer.u, n);
	  for (i = 0; i < n; ++i)
	    dbl[i] = in_buffer.u[i];
	  break;
	case MRI_INT16:
	  FRdInt16Array(ch->actual_file->fp, in_buffer.s, n);
	  for (i = 0; i < n; ++i)
	    dbl[i] = in_buffer.s[i];
	  break;
	case MRI_INT32:
	  FRdInt32Array(ch->actual_file->fp, in_buffer.i, n);
	  for (i = 0; i < n; ++i)
	    dbl[i] = in_buffer.i[i];
	  break;
	case MRI_INT64:
	  FRdInt64Array(ch->actual_file->fp, in_buffer.ll, n);
	  for (i = 0; i < n; ++i)
	    dbl[i] = in_buffer.ll[i];
	  break;
	case MRI_FLOAT32:
	  FRdFloat32Array(ch->actual_file->fp, in_buffer.f, n);
	  for (i = 0; i < n; ++i)
	    dbl[i] = in_buffer.f[i];
	  break;
	case MRI_FLOAT64:
	  FRdFloat64Array(ch->actual_file->fp, in_buffer.d, n);
	  for (i = 0; i < n; ++i)
	    dbl[i] = in_buffer.d[i];
	  break;
	default:
	  break;
	}
      if (bio_error)
	mri_report_error(ch->ds, "libmri: could not read chunk data from file %s\n",
			 ch->actual_file->name);
      read_offset += n * read_size;
       
      if (mri_fseek(f->fp, offset + write_offset, SEEK_SET) != 0)
	mri_report_error(ch->ds,
			 "libmri: could not seek to position %lld in file %s\n",
			 offset + write_offset,
			 f->name);

      bio_error = FALSE;
      switch (ch->datatype)
	{
	case MRI_UINT8:
	  for (i = 0; i < n; ++i)
	    out_buffer.u[i] = (unsigned char) dbl[i];
	  FWrUInt8Array(f->fp, out_buffer.u, n);
	  break;
	case MRI_INT16:
	  for (i = 0; i < n; ++i)
	    out_buffer.s[i] = (short) dbl[i];
	  FWrInt16Array(f->fp, out_buffer.s, n);
	  break;
	case MRI_INT32:
	  for (i = 0; i < n; ++i)
	    out_buffer.i[i] = (int) dbl[i];
	  FWrInt32Array(f->fp, out_buffer.i, n);
	  break;
	case MRI_INT64:
	  for (i = 0; i < n; ++i)
	    out_buffer.ll[i] = (long long) dbl[i];
	  FWrInt64Array(f->fp, out_buffer.ll, n);
	  break;
	case MRI_FLOAT32:
	  for (i = 0; i < n; ++i)
	    out_buffer.f[i] = (float) dbl[i];
	  FWrFloat32Array(f->fp, out_buffer.f, n);
	  break;
	case MRI_FLOAT64:
	  for (i = 0; i < n; ++i)
	    out_buffer.d[i] = dbl[i];
	  FWrFloat64Array(f->fp, out_buffer.d, n);
	  break;
	default:
	  mri_report_error(ch->ds, "libmri: Internal error - invalid chunk datatype\n");
	  abort();
	  break;
	}
      if (bio_error)
	mri_report_error(ch->ds, "libmri: could not write converted chunk data into file %s\n",
			 f->name);
       write_offset += n * write_size;
       count -= n;
     }

   n = ch->size - write_offset;
   if (n > 0)
     ClearBlock(f, offset + write_offset, n);
  bio_big_endian_input = saved_bio_big_endian_input;
  bio_big_endian_output = saved_bio_big_endian_output;
  bio_error = saved_bio_error;

  free (in_buffer.me);
  free (out_buffer.me);
#undef N_READ
}

static void
ClearBlock (MRI_File *dest_file, long long dest_offset, long long size)
{
  char buffer[BUFFER_SIZE];

  if (!OpenFile(dest_file, TRUE))
    {
      fprintf(stderr, "libmri: internal error in ClearBlock\n");
      abort();
    }
  memset(buffer, 0, (size_t)((size > BUFFER_SIZE) ? BUFFER_SIZE : size));
  if (mri_fseek(dest_file->fp, dest_offset, SEEK_SET) != 0)
    {
      fprintf(stderr, "libmri: could not seek to position %lld in file %s\n",
	      dest_offset, dest_file->name);
      abort();
    }
  while (size > 0)
    {
      if (fwrite(buffer,
		 (size_t) ((size > BUFFER_SIZE) ? BUFFER_SIZE : size), 1,
		 dest_file->fp) != 1)
	{
	  fprintf(stderr, "libmri: cannot clear chunk in file %s\n",
		  dest_file->name);
	  abort();
	}
      size -= BUFFER_SIZE;
    }
}

static MRI_File *
GetFile (MRI_Dataset *ds, char *name)
{
  MRI_File *f;
  char filename[MRI_MAX_FILENAME_LENGTH+1];

  MakeChunkFilename(filename, ds->name, name);

  for (f = ds->files; f != NULL; f = f->next)
    if (strcmp(f->name, filename) == 0)
      return(f);
  
  f = (MRI_File *) malloc(sizeof(MRI_File));
  f->next = ds->files;
  f->ds = ds;
  f->name = (char *) malloc(strlen(filename)+1);
  strcpy(f->name, filename);
  f->fp = NULL;
  f->last_use = 0;
  f->used = TRUE;
  f->external = FALSE;
  f->writeable = FALSE;
  ds->files = f;
#ifdef AFS
  ds->some_parts_in_afs |= CheckForAFS(filename);
#endif
  return(f);
}


static void
DestroyFile (MRI_File *file)
{
  MRI_File *pf, *f;
  MRI_Dataset *ds;

  CloseFile(file);
  
  ds = file->ds;
  pf = NULL;
  for (f = ds->files; f != NULL; pf = f, f = f->next)
    if (f == file)
      {
	if (pf != NULL)
	  pf->next = f->next;
	else
	  ds->files = f->next;
#ifdef DEBUG
	Log("GOING TO UNLINK FILE!!! (within DestroyFile): %s\n",
	    f->name);
#endif
	unlink(f->name);
	free(f);
	return;
      }
}


static char *
GetBuffer (MRI_Dataset *ds, long long size)
{
  MRI_Buffer *pb, *b, *nb;
  int count;

  /* first check if there are not recently used buffers
     of the given size */
  pb = NULL;
  b = ds->buffers;
  count = 0;
  while (b != NULL)
    {
      nb = b->next;
      if (++count > MRI_SAFE_BUFFER_COUNT &&
	  b->size == size)
	{
	  pb->next = nb;
	  b->next = ds->buffers;
	  ds->buffers = b;
	  return(b->buffer);
	}
      if (count > MRI_MAX_BUFFER_COUNT &&
	  nb == NULL)
	{
	  /* take the least recently used buffer and reallocate
	     it with the new size */
	  pb->next = NULL;
	  b->next = ds->buffers;
	  ds->buffers = b;
	  b->buffer = realloc(b->buffer, (size_t) size);
	  b->size = size;
	  return(b->buffer);
	}
      pb = b;
      b = nb;
    }
  
  /* allocate a fresh buffer */
  b = (MRI_Buffer *) malloc(sizeof(MRI_Buffer));
  b->buffer = malloc((size_t) size);
  b->size = size;
  b->next = ds->buffers;
  ds->buffers = b;
  return(b->buffer);
}

static int
MRI_TypeLength (MRI_Datatype datatype)
{
  switch (datatype)
    {
    case MRI_UINT8:	return(1);
    case MRI_INT16:	return(2);
    case MRI_INT32:	return(4);
    case MRI_INT64:	return(8);
    case MRI_FLOAT32:	return(4);
    case MRI_FLOAT64:	return(8);
    default:		abort();
    }
  return 0; /* not reached */
}

static int
MRI_ArrayTypeLength (MRI_ArrayType type, int n)
{
  int result;

  switch (type & 0x0f)
    {

    case MRI_RAW: result= 1; break;
    case MRI_UNSIGNED_CHAR: result= 1; break;
    case MRI_SHORT: result= 2; break;
    case MRI_INT: result= 4; break;
    case MRI_FLOAT: result= 4; break;
    case MRI_DOUBLE: result= 8; break;
    case MRI_LONG: result= sizeof(long); break; /* whatever it may be */
    case MRI_LONGLONG: result= 8; break;

    default:		abort();
    }
  result *= n;
  if (type & MRI_COMPLEX) result *= 2;
  return result;
}

#ifdef DEBUG
static void
Log (char *fmt, ...)
{
  va_list args;
  time_t t;
  char ts[64];

  t = time(NULL);
  strcpy(ts, ctime(&t));
  ts[24] = '\0';
  fprintf(logfile, "%s: ", ts);

  va_start(args, fmt);
  vfprintf(logfile, fmt, args);
  va_end(args);
  fflush(logfile);
}
#endif

#ifdef NO_FSEEK64
static int
mri_fseek (FILE *f, long long offset, int whence)
{
#ifdef HAVE_FSEEKO
  return(fseeko(f, (off_t) offset, whence));
#else
  if (offset < 0 || offset >= LONG_MAX)
    {
      mri_report_error(NULL, "mri_fseek: offset out-of-range\n");
      return(-1);
    }
  return(fseek(f, (long) offset, whence));
#endif
}

static long long
mri_ftell (FILE *f)
{
#ifdef HAVE_FSEEKO
  return((long long) ftello(f));
#else
  return((long long) ftell(f));
#endif
}

static int
mri_ftruncate (int fildes, long long length)
{
#ifdef HAVE_FSEEKO
  return(ftruncate(fildes, (off_t) length));
#else

#ifdef LINUX
#define FTRUNCATE_LENGTH_TYPE size_t
#else
#define FTRUNCATE_LENGTH_TYPE off_t
#endif
  if (length < 0 ||
      sizeof(FTRUNCATE_LENGTH_TYPE) == 4 && length >= 2147483647LL)
    {
      mri_report_error(NULL, "mri_ftruncate: length out-of-range\n");
      return(-1);
    }
  return(ftruncate(fildes, (FTRUNCATE_LENGTH_TYPE) length));

#endif
}
#endif

#ifdef AFS

static int MySystem (const char *command) {
  int pid, status;
  extern char **environ;
  extern int errno;
  
  if (command == 0)
    return 1;
  pid = fork();
  if (pid == -1)
    return -1;
  if (pid == 0) {
    char *argv[4];
    argv[0] = "sh";
    argv[1] = "-c";
    argv[2] = (char*)command;
    argv[3] = 0;
    execve("/bin/sh", argv, environ);
    exit(127);
  }
  do {
    if (waitpid(pid, &status, 0) == -1) {
      if (errno != EINTR)
        return -1;
    } else
      return status;
  } while(1);
}


static int CheckForAFS( const char* fname )
{
  char *cmd= "fs whereis";
  char buf[256];
  int retval;
  snprintf(buf,sizeof(buf),"%s %s >/dev/null 2>&1",cmd,(char*)fname);
  retval= MySystem(buf);
  return retval;
}

static void FlushAFS( MRI_Dataset* ds )
{
  MySystem("fs flushvolume >/dev/null 2>&1");
}

#endif

