/************************************************************
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *                                                          *
 *     Copyright (c) 1999 Carnegie Mellon University        *
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
 ***********************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define MAGIC_HEADER	689143391
#define MAGIC_TRAILER	1282787955

typedef struct Wrapper {
  size_t size;
  struct Wrapper *prev;
  struct Wrapper *next;
  int header;
} Wrapper;

#define WRAPPER_SIZE	sizeof(struct Wrapper)

static char rcsid[] = "$Id: libmdbg.c,v 1.6 2000/10/05 22:18:06 welling Exp $";

Wrapper *mdbg_list = 0;
int mdbg_count = 0;
int mdbg_start = 0;
int mdbg_end = 2000000000;
int mdbg_freq = 1;

void mdbg_test ();
void my_check_memory ();

void *my_malloc (size_t size)
{
  unsigned char *new;
  Wrapper *w;

  mdbg_test();
  new = (unsigned char *) malloc(size+WRAPPER_SIZE+4);
  w = (Wrapper *) new;
  w->size = size;
  w->prev = 0;
  w->next = mdbg_list;
  w->header = MAGIC_HEADER;
  if (mdbg_list != 0)
      mdbg_list->prev = w;
  mdbg_list = w;

  /* we need to do this byte by byte since the trailer
     is not necessarily word-aligned */
  *(new + (size + WRAPPER_SIZE)) = (MAGIC_TRAILER >> 24) & 0xff;
  *(new + (size + WRAPPER_SIZE + 1)) = (MAGIC_TRAILER >> 16) & 0xff;
  *(new + (size + WRAPPER_SIZE + 2)) = (MAGIC_TRAILER >> 8) & 0xff;
  *(new + (size + WRAPPER_SIZE + 3)) = MAGIC_TRAILER & 0xff;

  return((void *) (new+WRAPPER_SIZE));
}

void *my_strdup( const char* s )
{
  char* sdup;
  sdup= (char*)my_malloc( strlen(s)+1 );
  (void)strcpy( sdup, s );
  return sdup;
}

void my_free (void *ptr)
{
  unsigned char *p;
  Wrapper *w;
  int i;

  mdbg_test();
  p = (unsigned char *) ptr;
  w = (Wrapper *) (p - WRAPPER_SIZE);
  if (w->header != MAGIC_HEADER)
    abort();
  if (*(p + w->size) != ((MAGIC_TRAILER >> 24) & 0xff) ||
      *(p + w->size + 1) != ((MAGIC_TRAILER >> 16) & 0xff) ||
      *(p + w->size + 2) != ((MAGIC_TRAILER >> 8) & 0xff) ||
      *(p + w->size + 3) != (MAGIC_TRAILER & 0xff))
    abort();
  w->header = 0;
  *(p + w->size) = 0;
  *(p + w->size + 1) = 0;
  *(p + w->size + 2) = 0;
  *(p + w->size + 3) = 0;
  /* Zero out the freed memory, in the hopes of tripping up anything
   * that tries to access it later.
   */
  for (i=0; i<w->size; i++) *(p + i)= 0;

  if (w->prev != 0)
    w->prev->next = w->next;
  else
    mdbg_list = w->next;
  if (w->next != 0)
    w->next->prev = w->prev;
  w->prev = 0;
  w->next = 0;
  free((void *) w);
}

void *my_realloc (void *ptr, size_t size)
{
  void *new;
  unsigned char *p;
  Wrapper *w;

  new = my_malloc(size);
  if (ptr == NULL)
    return(new);

  p = (unsigned char *) ptr;
  w = (Wrapper *) (p - WRAPPER_SIZE);
  if (size < w->size)
    memcpy(new, ptr, size);
  else
    memcpy(new, ptr, w->size);

  my_free(ptr);
  return(new);
}

void *my_calloc (size_t nelem, size_t elsize)
{
  void *ptr;

  ptr = my_malloc(nelem * elsize);
  memset(ptr, 0, nelem * elsize);
  return(ptr);
}

void mdbg_test ()
{
  ++mdbg_count;
  if (mdbg_count >= mdbg_start &&
      mdbg_count <= mdbg_end &&
      (mdbg_count % mdbg_freq) == 0)
    my_check_memory();
}

void my_check_memory ()
{
  Wrapper *w;
  unsigned char *p;

  for (w = mdbg_list; w != 0; w = w->next)
    {
      if (w->header != MAGIC_HEADER)
	abort();
      p = ((unsigned char *) w) + WRAPPER_SIZE;
      if (*(p + w->size) != ((MAGIC_TRAILER >> 24) & 0xff) ||
	  *(p + w->size + 1) != ((MAGIC_TRAILER >> 16) & 0xff) ||
	  *(p + w->size + 2) != ((MAGIC_TRAILER >> 8) & 0xff) ||
	  *(p + w->size + 3) != (MAGIC_TRAILER & 0xff))
	abort();
    }
}

