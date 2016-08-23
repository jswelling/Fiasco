/*
 *	C Library for Reading and Writing the
 *	Pittsburgh MRI Format
 *
 *	Copyright (c) 1996 Pittsburgh Supercomputing Center
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
 *	History:
 *		3/96: Written by Greg Hood
 */
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include "splus_libmri.h"
#include "splus_libbio.h"

#ifdef USE_MMAP
#include <sys/mman.h>
#endif

#define TRUE		1
#define FALSE		0

#define BUFFER_SIZE	65536

typedef struct EmptyBlock {
  unsigned int start;
  unsigned int end;
} EmptyBlock;

typedef struct CopyRequest {
  struct CopyRequest *next;
  struct MRI_File *src_file;
  int src_offset;
  struct MRI_File *dest_file;
  int dest_offset;
  int size;
  MRI_Chunk *chunk;
} CopyRequest;

static char rcsid[] = "$Id: splus_libmri.c,v 1.6 2003/09/25 19:45:42 welling Exp $";

/* The following horrible macros work around restrictions on
 * memory allocation in Splus.  See the S-plus Programmer's Manual,
 * section 9.10.1 (Allocating Storage) for details.  Note that
 * calls to realloc() must be explicitly recoded.  Can't seem to
 * include S.h because of a function proto conflict.
 */
char *S_alloc(long n, int size);
char*S_realloc( char* p, long new, long old, int size );
#define malloc( size ) S_alloc(size,1)
#define free( pointer ) /* free is a no-op */

char *mri_error = NULL;

static MRI_Error_Handling error_handling = MRI_ABORT_ON_ERROR;
static unsigned int tmp_num = 0;
static unsigned int file_access = 1;

/* FORWARD DECLARATIONS */
static MRI_File *GetFile ();
static void CreateNewDataset ();
static MRI_Chunk *NewChunk ();
static MRI_File *CreateTempFile ();
static void DeallocateDataset ();
static int CheckForHooks ();
static char *GetBuffer ();
static void CopyBlock ();
static void ClearBlock ();
static MRI_KeyValue *FindInHashTable ();
static void RemoveFromHashTable ();
static int ReadKeyValuePair ();
static void ModifyChunk ();
static int ReadHeader ();
static int ReadString ();
static int ReadQuotedString ();
static int ReadBackslashedChar ();
static int HashFunction ();
static void WriteHeader ();
static void WriteKeyValuePair ();
static void WriteString ();
static void WriteQuotedString ();
static void ComputeChunkPositions ();
static int ReserveBlock ();
static int OpenFile ();
static void CloseFile ();
static void DestroyFile ();
static void RepositionChunk ();
static void UpdateChunkAttributes ();
static void ConvertChunk ();
static int MRI_TypeLength ();
static void CheckForStdImages ();
static int InitializeKeyValue ();
static void SetChunkNotReady ();
static int PrepareToRead ();
static int PrepareToWrite ();
static void CleanFiles ();
static int CompareKeys ();
static void DeallocateChunk ();
static int CheckForRemovalHooks ();
static int ConvertDatatype ();
static void mri_report_error (MRI_Dataset *ds, char *fmt, ...);

/* PUBLIC FUNCTIONS */

/*-------- OPENING AND CLOSING DATASETS ---------------------------*/

MRI_Dataset
*mri_open_dataset (const char *name, MRI_OpenMode mode)
{
  MRI_Dataset *ds;
  MRI_KeyValue *kv;
  int i;
  MRI_Chunk *ch;
  int first_start;
  char file_name[MRI_MAX_FILENAME_LENGTH+1];

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
      mri_report_error(ds, "mri_load_dataset: cannot allocate space for MRI_Dataset\n");
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
  system("fs flushvolume");
#endif

  if (!OpenFile(ds->header_file, ds->mode == MRI_WRITE || ds->mode == MRI_MODIFY))
    {
      mri_report_error(ds, "mri_open_dataset: could not open file %s for %s\n",
		ds->header_file->name,
		(ds->mode == MRI_WRITE || ds->mode == MRI_MODIFY) ? "writing" : "reading");
      DeallocateDataset(ds);
      return(NULL);
    }
  if (ds->mode == MRI_WRITE)
    CreateNewDataset(ds);
  else if (!ReadHeader(ds, ds->header_file->fp))
    {
      DeallocateDataset(ds);
      return(NULL);
    }

  /* build the associated chunk data structures */
  for (i = 0; i < ds->hash_table_size; ++i)
    for (kv = ds->hash_table[i]; kv != NULL; kv = kv->next_in_hash_table)
      if (strcmp(kv->value, "[chunk]") == 0)
	(void) NewChunk(ds, kv->key);

  /* check if there is actually more space reserved for the header */
  first_start = 999999999;
  for (ch = ds->chunks; ch != NULL; ch = ch->next)
    if (ch->file == ds->header_file &&
	ch->offset < first_start)
      first_start = ch->offset;
  if (first_start < 999999999)
    ds->header_size = first_start;

  return(ds);
}

MRI_Dataset *mri_copy_dataset (const char *filename, MRI_Dataset *original)
{
  MRI_Dataset *ds, *nds;
  int count;
  int total, n;
  void *ptr;
  MRI_Chunk *ch;
  MRI_File *f;
  char *key;
  char file_key[MRI_MAX_KEY_LENGTH+1];
  char value[MRI_MAX_VALUE_LENGTH+1];

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
      n = strlen(key);
      if (n >= 5 && strcmp(&key[n-5], ".file") == 0 &&
	  value[0] != '.')
	{
	  /* make sure the key belongs to a chunk */
	  strcpy(file_key, key);
	  file_key[n-5] = '\0';
	  for (ch = ds->chunks; ch != NULL; ch = ch->next)
	    if (strcmp(ch->name, file_key) == 0)
	      break;
	  if (ch != NULL)
	    {
	      /* construct the new filename */
	      count = 1;
	      for (f = ds->files; f != NULL && f != ch->file; f = f->next)
		if (!f->external && f->name[0] != '.')
		  ++count;
	      sprintf(value, ".%d.dat", count);
	    }
	}

      /* set the key's value in the new dataset */
      mri_set_string(nds, key, value);
    }

  /* now go through and copy all the chunks */
  for (ch = ds->chunks; ch != NULL; ch = ch->next)
    if (ch->order != MRI_EXTERNAL)
      for (total = 0; total < ch->size; total += BUFFER_SIZE)
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
  long pos;
  MRI_File *temp;
  MRI_Chunk *ch;
  int alone;

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
      abort();
    }

  WriteHeader(ds, !alone, temp->fp);
  pos = ftell(temp->fp);
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
	  abort();
	}
      WriteHeader(ds, !alone, temp->fp);
      pos = ftell(temp->fp);
      CloseFile(temp);
    }

  /* make sure to write all the chunks out */
  for (ch = ds->chunks; ch != NULL; ch = ch->next)
    if (ch->modified)
      RepositionChunk(ch, NULL);

  /* copy header from temp file to proper location */
  CopyBlock(ds->header_file, 0, temp, 0, (int) pos);
  /* reset header size to actual size so gap gets zeroed */
  ds->header_size = pos;

  /* throw away temp file */
  DestroyFile(temp);

  /* zero in all the gaps */
  CleanFiles(ds);

#ifdef AFS
  system("fs flushvolume");
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
	(void) unlink(f->name);
      }

  /* deallocate */
  DeallocateDataset(ds);
}


/*-------- GETTING KEY VALUES ------------------------------------*/

int
mri_get_int (MRI_Dataset *ds, const char *key)
{
  int v;

  if (!mri_has(ds, key))
    {
      mri_report_error(ds, "mri_get_int: non-existent key\n");
      return(MRI_UNSPECIFIED);
    }
  if (sscanf(mri_get_string(ds, key), "%d", &v) != 1)
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
      mri_report_error(ds, "mri_get_float: non-existent key\n");
      return((double) MRI_UNSPECIFIED);
    }
  if (sscanf(mri_get_string(ds, key), "%lf", &v) != 1)
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
      mri_report_error(ds, "mri_get_string: non-existent key\n");
      return(NULL);
    }
  return(kv->value);
}


/*-------- SETTING KEY VALUES ------------------------------------*/

void
mri_set_int (MRI_Dataset *ds, const char *key, const int value)
{
  char s[64];

  sprintf(s, "%d", value);
  mri_set_string(ds, key, s);
}

void
mri_set_float (MRI_Dataset *ds, const char *key, const float value)
{
  char s[64];

  sprintf(s, "%f", value);
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
mri_get_chunk (MRI_Dataset *ds, const char *key, int size, int offset, MRI_ArrayType type)
{
  MRI_Chunk *ch;
  void *buf;
  int i;
  unsigned char *uchar_buf;
  short *short_buf;
  int *int_buf;
  float *float_buf;
  double *double_buf;

  for (ch = ds->chunks; ch != NULL; ch = ch->next)
    if (strcmp(ch->name, key) == 0)
      break;
  if (ch == NULL)
    {
      mri_report_error(ds, "mri_get_chunk: no such chunk!\n");
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
  bio_big_endian_input = !ch->little_endian;

  switch (type)
    {
    case MRI_RAW:
      if (ch->mmapped)
	return(ch->addr + offset);
      else
	{
	  buf = GetBuffer(ds, size);
	  fseek(ch->file->fp, (long) (ch->offset + offset), SEEK_SET);
	  FRdUInt8Array(ch->file->fp, buf, size);
	  return(buf);
	}
      /* break; -NOTREACHED- */

    case MRI_UNSIGNED_CHAR:
      if (ch->datatype != MRI_UINT8)
	{
	  mri_report_error(ds, "mri_get_chunk: chunk's datatype will not fit into unsigned char\n");
	  return(NULL);
	}
      if (sizeof(char) == 1 && ch->mmapped)
	return(ch->addr + offset);
      else
	{
	  buf = GetBuffer(ds, size*sizeof(unsigned char));
	  if (ch->mmapped)
	    BRdUInt8Array((unsigned char *) (ch->addr + offset), buf, size);
	  else
	    {
	      fseek(ch->file->fp, (long) (ch->offset + offset), SEEK_SET);
	      FRdUInt8Array(ch->file->fp, buf, size);
	    }
	  return(buf);
	}
      /* break; -NOTREACHED- */

    case MRI_SHORT:
      switch (ch->datatype)
	{
	case MRI_UINT8:
	  short_buf = (short *) GetBuffer(ds, size*sizeof(short));
	  if (sizeof(unsigned char) == 1 && ch->mmapped)
	    uchar_buf = (unsigned char*)(ch->addr + offset);
	  else
	    {
	      uchar_buf = (unsigned char *) GetBuffer(ds, size*sizeof(unsigned char));
	      fseek(ch->file->fp, (long) (ch->offset + offset), SEEK_SET);
	      FRdUInt8Array(ch->file->fp, uchar_buf, size);
	    }
	  for (i = 0; i < size; ++i)
	    short_buf[i] = uchar_buf[i];
	  return(short_buf);
	  /* break; -NOTREACHED- */
	case MRI_INT16:
	  if (sizeof(short) == 2 && (ch->little_endian != bio_big_endian_machine) && ch->mmapped)
	    return(ch->addr + 2*offset);
	  else
	    {
	      buf = GetBuffer(ds, size*sizeof(short));
	      if (ch->mmapped)
		BRdInt16Array((unsigned char *) (ch->addr + 2*offset), buf, size);
	      else
		{
		  fseek(ch->file->fp, (long) (ch->offset + 2*offset), SEEK_SET);
		  FRdInt16Array(ch->file->fp, buf, size);
		}
	      return(buf);
	    }
	  /* break; -NOTREACHED- */
	default:
	  mri_report_error(ds, "mri_get_chunk: chunk's datatype will not fit into short\n");
	  return(NULL);
	}
      /* break; -NOTREACHED- */

    case MRI_INT:
      switch (ch->datatype)
	{
	case MRI_UINT8:
	  int_buf = (int *) GetBuffer(ds, size*sizeof(int));
	  if (sizeof(unsigned char) == 1 && ch->mmapped)
	    uchar_buf = (unsigned char*)(ch->addr + offset);
	  else
	    {
	      uchar_buf = (unsigned char *) GetBuffer(ds, size*sizeof(unsigned char));
	      fseek(ch->file->fp, (long) (ch->offset + offset), SEEK_SET);
	      FRdUInt8Array(ch->file->fp, uchar_buf, size);
	    }
	  for (i = 0; i < size; ++i)
	    int_buf[i] = uchar_buf[i];
	  return(int_buf);
	  /* break; -NOTREACHED- */
	case MRI_INT16:
	  int_buf = (int *) GetBuffer(ds, size*sizeof(int));
	  if (sizeof(short) == 2 && (ch->little_endian != bio_big_endian_machine) && ch->mmapped)
	    short_buf = (short *) (ch->addr + 2*offset);
	  else
	    {
	      short_buf = (short *) GetBuffer(ds, size*sizeof(short));
	      fseek(ch->file->fp, (long) (ch->offset + 2*offset), SEEK_SET);
	      FRdInt16Array(ch->file->fp, short_buf, size);
	    }
	  for (i = 0; i < size; ++i)
	    int_buf[i] = short_buf[i];
	  return(int_buf);
	  /* break; -NOTREACHED- */
	case MRI_INT32:
	  if (sizeof(int) == 4 && (ch->little_endian != bio_big_endian_machine) && ch->mmapped)
	    return(ch->addr + 4*offset);
	  else
	    {
	      buf = GetBuffer(ds, size*sizeof(int));
	      if (ch->mmapped)
		BRdInt32Array((unsigned char *) (ch->addr + 4*offset), buf, size);
	      else
		{
		  fseek(ch->file->fp, (long) (ch->offset + 4*offset), SEEK_SET);
		  FRdInt32Array(ch->file->fp, buf, size);
		}
	      return(buf);
	    }
	  /* break; -NOTREACHED- */
	default:
	  mri_report_error(ds, "mri_get_chunk: chunk's datatype will not fit into int\n");
	  return(NULL);
	}
      /* break; -NOTREACHED- */

    case MRI_FLOAT:
      switch (ch->datatype)
	{
	case MRI_UINT8:
	  float_buf = (float *) GetBuffer(ds, size*sizeof(float));
	  if (sizeof(unsigned char) == 1 && ch->mmapped)
	    uchar_buf = (unsigned char*)(ch->addr + offset);
	  else
	    {
	      uchar_buf = (unsigned char *) GetBuffer(ds, size*sizeof(unsigned char));
	      fseek(ch->file->fp, (long) (ch->offset + offset), SEEK_SET);
	      FRdUInt8Array(ch->file->fp, uchar_buf, size);
	    }
	  for (i = 0; i < size; ++i)
	    float_buf[i] = uchar_buf[i];
	  return(float_buf);
	  /* break; -NOTREACHED- */
	case MRI_INT16:
	  float_buf = (float *) GetBuffer(ds, size*sizeof(float));
	  if (sizeof(short) == 2 && (ch->little_endian != bio_big_endian_machine) && ch->mmapped)
	    short_buf = (short *) (ch->addr + 2*offset);
	  else
	    {
	      short_buf = (short *) GetBuffer(ds, size*sizeof(short));
	      fseek(ch->file->fp, (long) (ch->offset + 2*offset), SEEK_SET);
	      FRdInt16Array(ch->file->fp, short_buf, size);
	    }
	  for (i = 0; i < size; ++i)
	    float_buf[i] = short_buf[i];
	  return(float_buf);
	  /* break; -NOTREACHED- */
	case MRI_INT32:
	  float_buf = (float *) GetBuffer(ds, size*sizeof(float));
	  if (sizeof(int) == 4 && (ch->little_endian != bio_big_endian_machine) && ch->mmapped)
	    int_buf = (int *) (ch->addr + 4*offset);
	  else
	    {
	      int_buf = (int *) GetBuffer(ds, size*sizeof(int));
	      fseek(ch->file->fp, (long) (ch->offset + 4*offset), SEEK_SET);
	      FRdInt32Array(ch->file->fp, int_buf, size);
	    }
	  for (i = 0; i < size; ++i)
	    float_buf[i] = int_buf[i];
	  return(float_buf);
	  /* break; -NOTREACHED- */
	case MRI_FLOAT32:
	  if (sizeof(float) == 4 && (ch->little_endian != bio_big_endian_machine) && ch->mmapped)
	    return(ch->addr + 4*offset);
	  else
	    {
	      buf = GetBuffer(ds, size*sizeof(float));
	      if (ch->mmapped)
		BRdFloat32Array((unsigned char *) (ch->addr + 4*offset), buf, size);
	      else
		{
		  fseek(ch->file->fp, (long) (ch->offset + 4*offset), SEEK_SET);
		  FRdFloat32Array(ch->file->fp, buf, size);
		}
	      return(buf);
	    }
	  /* break; -NOTREACHED- */
	default:
	  mri_report_error(ds, "mri_get_chunk: chunk's datatype will not fit into float\n");
	  return(NULL);
	}
      /* break; -NOTREACHED- */

    case MRI_DOUBLE:
      switch (ch->datatype)
	{
	case MRI_UINT8:
	  double_buf = (double *) GetBuffer(ds, size*sizeof(double));
	  if (sizeof(unsigned char) == 1 && ch->mmapped)
	    uchar_buf = (unsigned char *)(ch->addr + offset);
	  else
	    {
	      uchar_buf = (unsigned char *) GetBuffer(ds, size*sizeof(unsigned char));
	      fseek(ch->file->fp, (long) (ch->offset + offset), SEEK_SET);
	      FRdUInt8Array(ch->file->fp, uchar_buf, size);
	    }
	  for (i = 0; i < size; ++i)
	    double_buf[i] = uchar_buf[i];
	  return(double_buf);
	  /* break; -NOTREACHED- */
	case MRI_INT16:
	  double_buf = (double *) GetBuffer(ds, size*sizeof(double));
	  if (sizeof(short) == 2 && (ch->little_endian != bio_big_endian_machine) && ch->mmapped)
	    short_buf = (short *) (ch->addr + 2*offset);
	  else
	    {
	      short_buf = (short *) GetBuffer(ds, size*sizeof(short));
	      fseek(ch->file->fp, (long) (ch->offset + 2*offset), SEEK_SET);
	      FRdInt16Array(ch->file->fp, short_buf, size);
	    }
	  for (i = 0; i < size; ++i)
	    double_buf[i] = short_buf[i];
	  return(double_buf);
	  /* break; -NOTREACHED- */
	case MRI_INT32:
	  double_buf = (double *) GetBuffer(ds, size*sizeof(double));
	  if (sizeof(int) == 4 && (ch->little_endian != bio_big_endian_machine) && ch->mmapped)
	    int_buf = (int *) (ch->addr + 4*offset);
	  else
	    {
	      int_buf = (int *) GetBuffer(ds, size*sizeof(int));
	      fseek(ch->file->fp, (long) (ch->offset + 4*offset), SEEK_SET);
	      FRdInt32Array(ch->file->fp, int_buf, size);
	    }
	  for (i = 0; i < size; ++i)
	    double_buf[i] = int_buf[i];
	  return(double_buf);
	  /* break; -NOTREACHED- */
	case MRI_FLOAT32:
	  double_buf = (double *) GetBuffer(ds, size*sizeof(double));
	  if (sizeof(float) == 4 && (ch->little_endian != bio_big_endian_machine) && ch->mmapped)
	    float_buf = (float *) (ch->addr + 4*offset);
	  else
	    {
	      float_buf = (float *) GetBuffer(ds, size*sizeof(float));
	      fseek(ch->file->fp, (long) (ch->offset + 4*offset), SEEK_SET);
	      FRdFloat32Array(ch->file->fp, float_buf, size);
	    }
	  for (i = 0; i < size; ++i)
	    double_buf[i] = float_buf[i];
	  return(double_buf);
	  /* break; -NOTREACHED- */
	case MRI_FLOAT64:
	  if (sizeof(double) == 8 && (ch->little_endian != bio_big_endian_machine) && ch->mmapped)
	    return(ch->addr + 8*offset);
	  else
	    {
	      buf = GetBuffer(ds, size*sizeof(double));
	      if (ch->mmapped)
		BRdFloat64Array((unsigned char *) (ch->addr + 8*offset), buf, size);
	      else
		{
		  fseek(ch->file->fp, (long) (ch->offset + 8*offset), SEEK_SET);
		  FRdFloat64Array(ch->file->fp, buf, size);
		}
	      return(buf);
	    }
	  /* break; -NOTREACHED- */
	default:
	  mri_report_error(ds, "mri_get_chunk: chunk's datatype will not fit into double\n");
	  return(NULL);
	}
      /* break; -NOTREACHED- */

    default:
      mri_report_error(ds, "mri_get_chunk: invalid array type\n");
      break;
    }
  return(NULL);
}

void
mri_set_chunk (MRI_Dataset *ds, const char *key,
	       int size, int offset,
	       MRI_ArrayType type, void *buf)
{
  MRI_Chunk *ch;
  char *addr;
  int i;
  unsigned char *uchar_buf;
  short *short_buf;
  int *int_buf;
  float *float_buf;
  double *double_buf;

  for (ch = ds->chunks; ch != NULL; ch = ch->next)
    if (strcmp(ch->name, key) == 0)
      break;
  if (ch == NULL)
    {
      mri_report_error(ds, "mri_set_chunk: no such chunk!\n");
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
  bio_big_endian_output = !ch->little_endian;

  switch (type)
    {
    case MRI_RAW:
      if (ch->mmapped)
	{
	  addr = ch->addr + offset;
	  if (addr != (char*)buf)
	    memcpy(addr, buf, size);
	}
      else
	{
	  fseek(ch->file->fp, (long) (ch->offset + offset), SEEK_SET);
	  FWrUInt8Array(ch->file->fp, buf, size);
	}
      break;

    case MRI_UNSIGNED_CHAR:
      switch (ch->datatype)
	{
	case MRI_UINT8:
	  if (sizeof(unsigned char) == 1 && ch->mmapped)
	    {
	      addr = ch->addr + offset;
	      if (addr != (char*)buf)
		memcpy(addr, buf, size);
	    }
	  else if (ch->mmapped)
	    BWrUInt8Array((unsigned char *) (ch->addr + offset), buf, size);
	  else
	    {
	      fseek(ch->file->fp, (long) (ch->offset + offset), SEEK_SET);
	      FWrUInt8Array(ch->file->fp, buf, size);
	    }
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
	  mri_report_error(ds, "mri_set_chunk: unsigned char will not fit into chunk's datatype\n");
	  break;
	}
      break;

    case MRI_SHORT:
      switch (ch->datatype)
	{
	case MRI_INT16:
	  if (sizeof(short) == 2 && (ch->little_endian != bio_big_endian_machine) && ch->mmapped)
	    {
	      addr = ch->addr + 2*offset;
	      if (addr != (char*)buf)
		memcpy(addr, buf, 2*size);
	    }
	  else if (ch->mmapped)
	    BWrInt16Array((unsigned char *) (ch->addr + 2*offset), buf, size);
	  else
	    {
	      fseek(ch->file->fp, (long) (ch->offset + 2*offset), SEEK_SET);
	      FWrInt16Array(ch->file->fp, buf, size);
	    }
	  break;
	case MRI_INT32:
	  short_buf = (short *) buf;
	  int_buf = (int *) GetBuffer(ds, size*sizeof(int));
	  for (i = 0; i < size; ++i)
	    int_buf[i] = short_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_INT, int_buf);
	  break;
	case MRI_FLOAT32:
	  short_buf = (short *) buf;
	  float_buf = (float *) GetBuffer(ds, size*sizeof(float));
	  for (i = 0; i < size; ++i)
	    float_buf[i] = short_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_FLOAT, float_buf);
	  break;
	case MRI_FLOAT64:
	  short_buf = (short *) buf;
	  double_buf = (double *) GetBuffer(ds, size*sizeof(double));
	  for (i = 0; i < size; ++i)
	    double_buf[i] = short_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_DOUBLE, double_buf);
	  break;
	default:
	  mri_report_error(ds, "mri_set_chunk: short will not fit into chunk's datatype\n");
	  break;
	}
      break;

    case MRI_INT:
      switch (ch->datatype)
	{
	case MRI_INT32:
	  if (sizeof(int) == 4 && (ch->little_endian != bio_big_endian_machine) && ch->mmapped)
	    {
	      addr = ch->addr + 4*offset;
	      if (addr != (char*)buf)
		memcpy(addr, buf, 4*size);
	    }
	  else if (ch->mmapped)
	    BWrInt32Array((unsigned char *) (ch->addr + 4*offset), buf, size);
	  else
	    {
	      fseek(ch->file->fp, (long) (ch->offset + 4*offset), SEEK_SET);
	      FWrInt32Array(ch->file->fp, buf, size);
	    }
	  break;
	case MRI_FLOAT32:
	  int_buf = (int *) buf;
	  float_buf = (float *) GetBuffer(ds, size*sizeof(float));
	  for (i = 0; i < size; ++i)
	    float_buf[i] = int_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_FLOAT, float_buf);
	  break;
	case MRI_FLOAT64:
	  int_buf = (int *) buf;
	  double_buf = (double *) GetBuffer(ds, size*sizeof(double));
	  for (i = 0; i < size; ++i)
	    double_buf[i] = int_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_DOUBLE, double_buf);
	  break;
	default:
	  mri_report_error(ds, "mri_set_chunk: int will not fit into chunk's datatype\n");
	  break;
	}
      break;

    case MRI_FLOAT:
      switch (ch->datatype)
	{
	case MRI_FLOAT32:
	  if (sizeof(float) == 4 && (ch->little_endian != bio_big_endian_machine) && ch->mmapped)
	    {
	      addr = ch->addr + 4*offset;
	      if (addr != (char*)buf)
		memcpy(addr, buf, 4*size);
	    }
	  else if (ch->mmapped)
	    BWrFloat32Array((unsigned char *) (ch->addr + 4*offset), buf, size);
	  else
	    {
	      fseek(ch->file->fp, (long) (ch->offset + 4*offset), SEEK_SET);
	      FWrFloat32Array(ch->file->fp, buf, size);
	    }
	  break;
	case MRI_FLOAT64:
	  float_buf = (float *) buf;
	  double_buf = (double *) GetBuffer(ds, size*sizeof(double));
	  for (i = 0; i < size; ++i)
	    double_buf[i] = float_buf[i];
	  mri_set_chunk(ds, key, size, offset, MRI_DOUBLE, double_buf);
	  break;
	default:
	  mri_report_error(ds, "mri_set_chunk: float will not fit into chunk's datatype\n");
	  break;
	}
      break;

    case MRI_DOUBLE:
      switch (ch->datatype)
	{
	case MRI_FLOAT64:
	  if (sizeof(double) == 8 && (ch->little_endian != bio_big_endian_machine) && ch->mmapped)
	    {
	      addr = ch->addr + 8*offset;
	      if (addr != (char*)buf)
		memcpy(addr, buf, 8*size);
	    }
	  else if (ch->mmapped)
	    BWrFloat64Array((unsigned char *) (ch->addr + 8*offset), buf, size);
	  else
	    {
	      fseek(ch->file->fp, (long) (ch->offset + 8*offset), SEEK_SET);
	      FWrFloat64Array(ch->file->fp, buf, size);
	    }
	  break;
	default:
	  mri_report_error(ds, "mri_set_chunk: double will not fit into chunk's datatype\n");
	  break;
	}
      break;

    default:
      mri_report_error(ds, "mri_set_chunk: invalid array type\n");
      break;
    }
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

  for (b = ds->buffers; b != NULL; b = b->next)
    if (b->buffer == ptr)
      return;

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
mri_get_image (MRI_Dataset *ds, const int time, const int slice, MRI_ArrayType type)
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
		       (time * ds->std_n_slices + slice) * ds->std_image_size,
		       type & 0xf));
}

void
mri_set_image (MRI_Dataset *ds, const int time, const int slice, MRI_ArrayType type, void *image)
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
		(time * ds->std_n_slices + slice) * ds->std_image_size,
		type & 0xf, image);
}


/*------------------------INTERNAL FUNCTIONS----------------------------------*/
   
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
  int v;
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
	abort();
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

  if (mri_has(ds, mri_cat(name, ".size")))
    ch->size = mri_get_int(ds, mri_cat(name, ".size"));
  else
    {
      v = MRI_TypeLength(ch->datatype);
      for (i = 0; i < (int) strlen(ch->dimensions); ++i)
	v *= ch->extent[i];
      ch->size = v;
    }

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
  ch->mmapped = FALSE;
  ch->addr = NULL;
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
	empty_blocks[0].end = 4000000000U;
	n_empty_blocks = 1;
	
	for (ch = ds->chunks; ch != NULL; ch = ch->next)
	  if (ch->file == f)
	    (void) ReserveBlock(empty_blocks, &n_empty_blocks,
				ch->offset, ch->size);

	if (!OpenFile(f, TRUE))
	  abort();
	for (i = 0; i < n_empty_blocks; ++i)
	  if (empty_blocks[i].end < 4000000000U)
	    ClearBlock(f, empty_blocks[i].start,
		       empty_blocks[i].end - empty_blocks[i].start + 1);
	  else
	    {
	      fflush(f->fp);
	      ftruncate(fileno(f->fp), (off_t) (empty_blocks[i].start));
	    }
      }
}

static MRI_File *
CreateTempFile (MRI_Dataset *ds)
{
  char file_name[256];
  char *tmp_dir;
  MRI_File *f;

  if ((tmp_dir = getenv("MRI_TMP_DIR")) == NULL)
    tmp_dir = "/tmp";
  sprintf(file_name, "%s/mri%ld.%d", tmp_dir, getpid(), tmp_num++);
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
      mri_report_error(ds, "mri_open_dataset: duplicate keys found in header\n");
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
  int new_order;
  int new_offset;
  int new_extent;
  int new_size;
  int new_little_endian;
  MRI_Chunk *ch;
  char chunk_name[MRI_MAX_KEY_LENGTH+1];
  char key_name[MRI_MAX_KEY_LENGTH+1];

  /* check if we are creating a chunk */
  if (strcmp(kv->value, "[chunk]") == 0)
    {
      /* check if there is already a chunk with that name */
      for (ch = ds->chunks; ch != NULL; ch = ch->next)
	if (strcmp(ch->name, kv->key) == 0)
	  return(TRUE);

      /* create a brand new chunk */
      ch = NewChunk(ds, kv->key);
      ch->actual_file = ds->header_file;
      ch->actual_offset = 0;
      ch->actual_size = 0;
      ModifyChunk(ch);
      if (!bio_big_endian_machine &&
	  !mri_has(ds, mri_cat(ch->name, ".little_endian")))
	mri_set_int(ds, mri_cat(ch->name, ".little_endian"), 1);
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
      if (sscanf(kv->value, "%d", &new_offset) != 1)
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
      if (sscanf(kv->value, "%d", &new_extent) != 1)
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
      if (sscanf(kv->value, "%d", &new_size) != 1 ||
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
  int v;
  int size;
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
  if (ch->mmapped)
    {
#ifdef USE_MMAP
      if (msync((void *) ch->addr, ch->actual_size, MS_SYNC) < 0)
	{
	  fprintf(stderr, "libmri: could not msync chunk\n");
	  abort();
	}	
      if (munmap((void *) ch->addr, ch->actual_size) < 0)
	{
	  fprintf(stderr, "libmri: could not munmap chunk\n");
	  abort();
	}
#endif
      ch->addr = 0;
      ch->mmapped = FALSE;
    }
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
  int new_order;
  int new_offset;
  int new_extent;
  int new_size;
  int new_little_endian;
  MRI_Chunk *ch;
  char chunk_name[MRI_MAX_KEY_LENGTH+1];
  char key_name[MRI_MAX_KEY_LENGTH+1];

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
      fprintf(stderr, "Aborting...\n");
      abort();
    }
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
	  ds->header_size = ftell(f);
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
      mri_report_error(ds, "mri_load_dataset: = not found in key/value pair");
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
  while (isspace(c = fgetc(f))) ;

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
      /* break; -NOTREACHED- */
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
      /* quadruple the hash table */
      new = (MRI_KeyValue **) malloc(4 * ds->hash_table_size * sizeof(MRI_KeyValue *));
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
      ds->hash_table_size *= 4;
      hv = HashFunction(key);
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
  MRI_KeyValue *pkv, *kv, *nkv;
  MRI_KeyValue **new;
  int i;

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
  int i, n;
  MRI_KeyValue *kv;
  MRI_Chunk *ch;
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
  for (i = 0; s[i] != '\0'; ++i)
    if (s[i] < 32 && s[i] != '\t' ||
	s[i] > 126 ||
	s[i] == '=')
      {
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
  int i, n, m;
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
      empty_blocks[0].end = 4000000000U;
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
	  chunks[n]->offset = ReserveBlock(empty_blocks, &n_empty_blocks, -1, chunks[n]->size);
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

static int
ReserveBlock (EmptyBlock *empty_blocks, int *n_empty_blocks, int offset, int size)
{
  int i;
  int m;
  int start;
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

#ifdef USE_MMAP
  ch->addr = mmap(0, ch->size,
		  (ch->ds->mode == MRI_READ) ? PROT_READ : (PROT_READ | PROT_WRITE),
		  MAP_PRIVATE, fileno(ch->file->fp), ch->offset);
  if (((void *) ch->addr) != ((void *) -1))
    ch->mmapped = TRUE;
#endif

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

#ifdef USE_MMAP
  ch->addr = mmap(0, ch->size,
		  (ch->ds->mode == MRI_READ) ? PROT_READ : (PROT_READ | PROT_WRITE),
		  MAP_PRIVATE, fileno(ch->file->fp), ch->offset);
  if (((void *) ch->addr) != ((void *) -1))
    ch->mmapped = TRUE;
#endif

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
      if (file->writeable || !for_writing)
	return(TRUE);

      /* no, try to reopen it for writing */
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
	abort();
      CloseFile(oldest);
    }

  if (!for_writing)
    file->fp = fopen(file->name, "r");
  else if (ds->mode == MRI_WRITE && file->last_use == 0)
    /* only truncate if we have not opened this file before */
    file->fp = fopen(file->name, "w+");
  else
    {
      file->fp = fopen(file->name, "r+");
      if (file->fp == NULL)
	file->fp = fopen(file->name, "w+");
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
CopyBlock (MRI_File *dest_file, int dest_offset,
	   MRI_File *src_file, int src_offset,
	   int size)
{
  int n;
  int offset;
  char buffer[BUFFER_SIZE];

  /* check if this is a null operation */
  if (src_file == dest_file &&
      src_offset == dest_offset)
    return;

  if (!OpenFile(src_file, FALSE) ||
      !OpenFile(dest_file, TRUE))
    {
      fprintf(stderr, "libmri: Internal copying error.\n");
      abort();
    }
  offset = 0;
  while (size > 0)
    {
      n = (size > BUFFER_SIZE) ? BUFFER_SIZE : size;
      fseek(src_file->fp, (long) (src_offset+offset), SEEK_SET);
      fread(buffer, n, 1, src_file->fp);
      fseek(dest_file->fp, (long) (dest_offset+offset), SEEK_SET);
      fwrite(buffer, n, 1, dest_file->fp);
      size -= n;
      offset += n;
    }
}

static void
ConvertChunk (MRI_Chunk *ch, MRI_File *f, int offset)
{
#define N_READ	(BUFFER_SIZE / sizeof(double))
  int i;
  int n;
  int read_size, write_size;
  int read_offset, write_offset;
  int count;
  double *dbl;
  union {
    unsigned char	u[N_READ];
    short		s[N_READ];
    int			i[N_READ];
    float		f[N_READ];
    double		d[N_READ];
  } in_buffer, out_buffer;

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
      fseek(ch->actual_file->fp, (long) (ch->actual_offset + read_offset), SEEK_SET);
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
      read_offset += n * read_size;
       
       fseek(f->fp, (long) (offset + write_offset), SEEK_SET);
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
       write_offset += n * write_size;
       count -= n;
     }

   n = ch->size - write_offset;
   if (n > 0)
     ClearBlock(f, offset + write_offset, n);
#undef N_READ
}

static void
ClearBlock (MRI_File *dest_file, int dest_offset, int size)
{
  char buffer[BUFFER_SIZE];

  if (!OpenFile(dest_file, TRUE))
    {
      fprintf(stderr, "libmri: internal error in ClearBlock\n");
      abort();
    }
  memset(buffer, 0, (size > BUFFER_SIZE) ? BUFFER_SIZE : size);
  fseek(dest_file->fp, (long) dest_offset, SEEK_SET);
  while (size > 0)
    {
      fwrite(buffer, (size > BUFFER_SIZE) ? BUFFER_SIZE : size, 1, dest_file->fp);
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
	unlink(f->name);
	free(f);
	return;
      }
}


static char *
GetBuffer (MRI_Dataset *ds, int size)
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
	  /* The following subsitutes the Splus equivalent of the
	   * original code:
	   * b->buffer = realloc(b->buffer, size);
	   */
	  b->buffer= S_realloc(b->buffer, size, b->size, 1);
	  /* End substitution */
	  b->size = size;
	  return(b->buffer);
	}
      pb = b;
      b = nb;
    }
  
  /* allocate a fresh buffer */
  b = (MRI_Buffer *) malloc(sizeof(MRI_Buffer));
  b->buffer = malloc(size);
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
    case MRI_FLOAT32:	return(4);
    case MRI_FLOAT64:	return(8);
    default:		abort();
    }
  return 0;
}


