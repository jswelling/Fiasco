/************************************************************
 *                                                          *
 *  kvhash.c                                                *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1995 Department of Statistics,         *
 *                        Carnegie Mellon University        *
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
 *                                                          *
 *  Original programming by Mark Fitzgerald  5-96           *
 *  Modified to exclusively use libmri calls for output,    *
 *       Greg Hood (PSC), 9-98                              *
 *  Modified to read header files, and to use LX2           *
 *       resampling stuff, Joel Welling (PSC/Stats), 5-1999 *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "fmri.h"
#include "misc.h"
#include "kvhash.h"

static char rcsid[] = "$Id: kvhash.c,v 1.6 2003/07/01 21:35:22 welling Exp $";

void kvDumpTableUnique( KVHash* kvh, FILE* ofile, const char* prefix,
			int maxDepth, int depth )
{
  KVIterator* kvi= kvUniqueIteratorFactory(kvh);
  while (kvIteratorHasMorePairs(kvi)) {
    KVPair* p= kvIteratorNextPair(kvi);
    switch (kvType(p)) {
    case KV_STRING:
      fprintf(ofile,"%s.%s: <%s>\n",prefix,p->key,p->v.s);
      break;
    case KV_DOUBLE:
      fprintf(ofile,"%s.%s: %g (double)\n",prefix,p->key,p->v.d);
      break;
    case KV_LONG:
      fprintf(ofile,"%s.%s: %lld (long long)\n",prefix,p->key,p->v.l);
      break;
    case KV_INT:
      fprintf(ofile,"%s.%s: %d (int)\n",prefix,p->key,(int)(p->v.l));
      break;
    case KV_BOOLEAN:
      fprintf(ofile,"%s.%s: %s\n",prefix,p->key,(p->v.l ? "TRUE":"FALSE"));
      break;
    case KV_HASH:
      {
	if (depth<maxDepth) {
	  char* buf;
	  if (!(buf=(char*)malloc(strlen(prefix)+strlen(p->key)+10)))
	    Abort("kvDumpTable: unable to allocate %d bytes!\n",
		  strlen(prefix)+strlen(p->key)+10);
	  sprintf(buf,"%s.%s",prefix,p->key);
	  kvDumpTableUnique(p->v.h, ofile, buf, maxDepth, depth+1);
	  free(buf);
	}
	else fprintf(ofile,"%s.%s: ***hash table***\n",prefix,p->key);
      }
      break;
    }
  }
}

static void stepToNext( KVIterator* kvi )
{
  if (kvi->p) {
    if (kvi->p->next) { kvi->p= kvi->p->next; }
    else if (kvi->id < kvi->hash->size - 1) {
      int i;
      for (i=kvi->id+1; i<kvi->hash->size; i++) {
	if (kvi->hash->table[i]) {
	  kvi->id= i;
	  kvi->p= kvi->hash->table[i];
	  break;
	}
      }
      if (i==kvi->hash->size) {
	kvi->p= NULL;
	kvi->id= kvi->hash->size - 1;
      }
    }
    else kvi->p= NULL;
  }
  /* else do nothing */
}

static void stepToNextUnique( KVIterator* kvi )
{
  if (kvi->p) {
    KVHash* kvh= kvi->hash;
    KVPair* next= kvi->p->next;
    while (next) {
      if (kvLookup(kvh,next->key)==next) {
	kvi->p= next;
	return;
      }
      next= next->next;
    }
    if (kvi->id < kvi->hash->size - 1) {
      /* First pair in a new row is always unique, since that
       * pair was the most recently inserted under that hash code.
       */
      int i;
      for (i=kvi->id+1; i<kvi->hash->size; i++) {
	if (kvi->hash->table[i]) {
	  kvi->id= i;
	  kvi->p= kvi->hash->table[i];
	  break;
	}
      }
      if (i==kvi->hash->size) {
	kvi->p= NULL;
	kvi->id= kvi->hash->size - 1;
      }
    }
    else kvi->p= NULL;
  }
}

KVIterator* kvIteratorFactory( KVHash* kvh )
{
  KVIterator* kvi;
  int i;

  if (!(kvi=(KVIterator*)malloc(sizeof(KVIterator))))
    Abort("kvIteratorFactory: unable to allocate %d bytes!\n",
	  sizeof(KVIterator));
  kvi->hash= kvh;
  kvi->stepToNext= stepToNext;
  kvi->hook= NULL;
  for (i=0; i<kvh->size; i++) {
    if (kvh->table[i]) {
      kvi->id= i;
      kvi->p= kvh->table[i];
      break;
    }
  }
  if (i==kvh->size) kvi->p= NULL;
  return kvi;
}

KVIterator* kvUniqueIteratorFactory( KVHash* kvh )
{
  KVIterator* kvi= kvIteratorFactory(kvh);
  kvi->stepToNext= stepToNextUnique;
  return kvi;
}

static void stepToNextSorted( KVIterator* kvi )
{
  if (kvi->p) {
    kvi->p= ((KVPair**)(kvi->hook))[++(kvi->id)];
  }
  /* else do nothing */
}

static int compareKeys( const void* p1, const void* p2 )
{
  KVPair** pair1= (KVPair**)p1;
  KVPair** pair2= (KVPair**)p2;
  return strcmp((*pair1)->key,(*pair2)->key);
}

KVIterator* kvSortedIteratorFactory( KVHash* kvh )
{
  int count;
  int i;
  KVIterator* unsorted= NULL;
  KVIterator* kvi= kvIteratorFactory(kvh);

  kvi->stepToNext= stepToNextSorted;

  /* Build and sort a table of unique pairs */
  unsorted= kvUniqueIteratorFactory( kvh );
  count= 0;
  while (kvIteratorHasMorePairs(unsorted)) {
    count++;
    kvIteratorNextPair(unsorted);
  }
  kvDestroyIterator(unsorted);
  if (!(kvi->hook= malloc((count+1)*sizeof(KVPair*))))
    Abort("kvSortedIteratorFactory: unable to allocate %d bytes!\n",
	  count*sizeof(KVPair*));
  unsorted= kvUniqueIteratorFactory( kvh );
  i= 0;
  while (kvIteratorHasMorePairs(unsorted)) {
    ((KVPair**)kvi->hook)[i]= kvIteratorNextPair(unsorted);
    i++;
  }
  kvDestroyIterator(unsorted);
  ((KVPair**)kvi->hook)[count]= NULL; /* to mark the end */
  if (count>0)
    qsort(kvi->hook, count, sizeof(KVPair*), compareKeys);

  kvi->id= 0; /* we will use this as offset into sorted table */
  kvi->p= ((KVPair**)kvi->hook)[kvi->id]; /* first element */
  return kvi;
}

void kvDestroyIterator( KVIterator* kvi )
{
  if (kvi->hook) free(kvi->hook);
  free(kvi);
}

int kvIteratorHasMorePairs( KVIterator* kvi ) 
{
  return( kvi->p != NULL );
}

KVPair* kvIteratorNextPair( KVIterator* kvi )
{
  KVPair* result= kvi->p;
  (*(kvi->stepToNext))(kvi);
  return result;
}

const char* kvKey( const KVPair* p )
{
  return (const char*)(p->key);
}

KVType kvType( const KVPair* p) 
{
  return p->type;
}

const char* kvTypeName( KVType t )
{
  switch (t) {
  case KV_STRING: return "string";
  case KV_HASH:   return "hash_table";
  case KV_LONG:   return "long_long";
  case KV_DOUBLE: return "double";
  case KV_BOOLEAN:return "boolean";
  case KV_INT:    return "int";
  }
  Abort("kvTypeName: unknown type %d; internal error!\n",(int)t);
  return("UnknownType");
}

KVHash *kvFactory( int size )
{
  KVHash* result;
  int i;

  if (!(result= (KVHash*)malloc(sizeof(KVHash)))) {
    Abort("kvFactory: unable to allocate %d bytes!\n",sizeof(KVHash));
  }
  if (!(result->table= (KVPair**)malloc(size*sizeof(KVPair*))))
    Abort("kvFactory: unable to allocate %d bytes!\n",size*sizeof(KVPair*));

  result->size= size;
  for (i=0; i<result->size; i++) result->table[i]= NULL;
  return result;
}

static void kvFreePair( KVPair* p )
{
  switch (p->type) {
  case KV_STRING: free(p->v.s); break;
  case KV_HASH: kvDestroy(p->v.h); break;
  case KV_LONG: break;
  case KV_DOUBLE: break;
  case KV_BOOLEAN: break;
  case KV_INT: break;
  }
  free(p->key);
  free(p);
}

void kvDestroy( KVHash* kvh )
{
  int i;

  for (i=0; i<kvh->size; i++) {
    KVPair* thisPair= kvh->table[i];
    while (thisPair) {
      KVPair* nextPair= thisPair->next;
      kvFreePair(thisPair);
      thisPair= nextPair;
    }
  }
  free(kvh->table);
  free(kvh);
}

static int hashKey(KVHash *kvh, const char *string)
/* Returns an integer value between 0 and size-1 to be used as <string>'s
 * index into the hash table
 */
{
  const char *p;
  unsigned h=0,g;
  
  for (p=string;*p != '\0';p++)
    {
      h = (h << 4) + (*p);
      if (g = h&0xf0000000)
        {
          h = h ^ (g >> 24);
          h = h ^ g;
        }
    }
  return(h%(kvh->size));
}

KVPair* kvLookup( KVHash* kvh, const char* string )
{
  KVPair* here= kvh->table[ hashKey(kvh,string) ];
  while (here) {
    if (!strcmp(string,here->key)) break;
    here= here->next;
  }
  return here;
}

KVHash* kvGetHash( KVHash* kvh, const char* key )
{
  KVPair* p= kvLookup(kvh,key);
  if (!p) return NULL;
  if (kvType(p)!=KV_HASH) 
    Abort("kvGetHash: value of <%s> is not a hash table!\n",key);
  return p->v.h;
}

char* kvGetString( KVHash* kvh, const char* key )
{
  KVPair* p= kvLookup(kvh,key);
  if (!p) return NULL;
  if (kvType(p)!=KV_STRING) 
    Abort("kvGetString: value of <%s> is not a string!\n",key);
  return p->v.s;
}

long long kvGetLong( KVHash* kvh, const char* key )
{
  KVPair* p= kvLookup(kvh,key);
  if (!p) Abort("kvGetLong: key <%s> is not defined!\n",key);
  if (kvType(p)!=KV_LONG && kvType(p)!=KV_INT) 
    Abort("kvGetLong: value of <%s> is not a long long!\n",key);
  return p->v.l;
}

int kvGetInt( KVHash* kvh, const char* key )
{
  KVPair* p= kvLookup(kvh,key);
  if (!p) return 0;
  if (kvType(p)!=KV_INT) 
    Abort("kvGetInt: value of <%s> is not an int!\n",key);
  return ((int)(p->v.l));
}

int kvGetBoolean( KVHash* kvh, const char* key )
{
  KVPair* p= kvLookup(kvh,key);
  if (!p) return 0;
  if (kvType(p)!=KV_BOOLEAN) 
    Abort("kvGetBoolean: value of <%s> is not a boolean!\n",key);
  return (p->v.l != 0);
}

double kvGetDouble( KVHash* kvh, const char* key )
{
  KVPair* p= kvLookup(kvh,key);
  if (!p) Abort("kvGetDouble: key <%s> is not defined!\n",key);
  if (kvType(p)!=KV_DOUBLE) 
    Abort("kvGetDouble: value of <%s> is not a double!\n",key);
  return p->v.d;
}

void kvDelete( KVHash* kvh, const char* string )
{
  KVPair* p= kvLookup(kvh,string);
  if (p) {
    KVPair* newPrev= p->prev;
    KVPair* newNext= p->next;
    if (newPrev) newPrev->next= newNext;
    else {
      /* This was the first one in this list */
      kvh->table[hashKey(kvh,string)]= newNext;
    }
    if (newNext) newNext->prev= newPrev;
    kvFreePair(p);
  }
}

void kvDeleteAll( KVHash* kvh, const char* string )
{
  while (kvLookup(kvh,string) != NULL) kvDelete(kvh,string);
}

void kvDefString( KVHash* kvh, const char* key, const char* value )
{
  KVPair* p;
  KVPair* head;
  int id;

  if (!(p= (KVPair*)malloc(sizeof(KVPair))))
    Abort("kvDefString: unable to allocate %d bytes!\n",sizeof(KVPair));
  if (!(p->key= strdup(key)))
    Abort("kvDefString: unable to duplicate a %d-char string!\n",
	  strlen(key));
  if (!(p->v.s= strdup(value)))
    Abort("kvDefString: unable to duplicate a %d-char string!\n",
	  strlen(value));
  p->type= KV_STRING;
  id= hashKey(kvh,key);
  head= kvh->table[ id ];
  if (head) head->prev= p;
  p->next= head;
  p->prev= NULL;
  kvh->table[ id ]= p;
}

void kvDefLong( KVHash* kvh, const char* key, const long long value )
{
  KVPair* p;
  KVPair* head;
  int id;

  if (!(p= (KVPair*)malloc(sizeof(KVPair))))
    Abort("kvDefString: unable to allocate %d bytes!\n",sizeof(KVPair));
  if (!(p->key= strdup(key)))
    Abort("kvDefString: unable to duplicate a %d-char string!\n",
	  strlen(key));
  p->v.l= value;
  p->type= KV_LONG;
  id= hashKey(kvh,key);
  head= kvh->table[ id ];
  if (head) head->prev= p;
  p->next= head;
  p->prev= NULL;
  kvh->table[ id ]= p;
}

void kvDefBoolean( KVHash* kvh, const char* key, const int value )
{
  KVPair* p;
  KVPair* head;
  int id;

  if (!(p= (KVPair*)malloc(sizeof(KVPair))))
    Abort("kvDefString: unable to allocate %d bytes!\n",sizeof(KVPair));
  if (!(p->key= strdup(key)))
    Abort("kvDefString: unable to duplicate a %d-char string!\n",
	  strlen(key));
  p->v.l= value;
  p->type= KV_BOOLEAN;
  id= hashKey(kvh,key);
  head= kvh->table[ id ];
  if (head) head->prev= p;
  p->next= head;
  p->prev= NULL;
  kvh->table[ id ]= p;
}

void kvDefInt( KVHash* kvh, const char* key, const int value )
{
  KVPair* p;
  KVPair* head;
  int id;

  if (!(p= (KVPair*)malloc(sizeof(KVPair))))
    Abort("kvDefString: unable to allocate %d bytes!\n",sizeof(KVPair));
  if (!(p->key= strdup(key)))
    Abort("kvDefString: unable to duplicate a %d-char string!\n",
	  strlen(key));
  p->v.l= value;
  p->type= KV_INT;
  id= hashKey(kvh,key);
  head= kvh->table[ id ];
  if (head) head->prev= p;
  p->next= head;
  p->prev= NULL;
  kvh->table[ id ]= p;
}

void kvDefDouble( KVHash* kvh, const char* key, const double value )
{
  KVPair* p;
  KVPair* head;
  int id;

  if (!(p= (KVPair*)malloc(sizeof(KVPair))))
    Abort("kvDefString: unable to allocate %d bytes!\n",sizeof(KVPair));
  if (!(p->key= strdup(key)))
    Abort("kvDefString: unable to duplicate a %d-char string!\n",
	  strlen(key));
  p->v.d= value;
  p->type= KV_DOUBLE;
  id= hashKey(kvh,key);
  head= kvh->table[ id ];
  if (head) head->prev= p;
  p->next= head;
  p->prev= NULL;
  kvh->table[ id ]= p;
}

void kvDefHash( KVHash* kvh, const char* key, KVHash* value )
{
  KVPair* p;
  KVPair* head;
  int id;

  if (!(p= (KVPair*)malloc(sizeof(KVPair))))
    Abort("kvDefString: unable to allocate %d bytes!\n",sizeof(KVPair));
  if (!(p->key= strdup(key)))
    Abort("kvDefString: unable to duplicate a %d-char string!\n",
	  strlen(key));
  p->v.h= value;
  p->type= KV_HASH;
  id= hashKey(kvh,key);
  head= kvh->table[ id ];
  if (head) head->prev= p;
  p->next= head;
  p->prev= NULL;
  kvh->table[ id ]= p;
}

KVHash *kvCloneUnique(KVHash* original)
{
  /* Create a duplicate of the input hash table, including
   * all unique entries recursively.
   */
  KVHash* result= kvFactory(original->size);
  KVIterator* kvi= kvUniqueIteratorFactory(original);

  while (kvIteratorHasMorePairs(kvi)) {
    KVPair* p= kvIteratorNextPair(kvi);
    switch (p->type) {
    case KV_STRING: kvDefString(result, p->key, p->v.s); break;
    case KV_HASH: kvDefHash(result, p->key, kvCloneUnique(p->v.h)); break;
    case KV_LONG: kvDefLong(result, p->key, p->v.l); break;
    case KV_DOUBLE: kvDefDouble(result, p->key, p->v.d); break;
    case KV_BOOLEAN: kvDefBoolean(result, p->key, p->v.l); break;
    case KV_INT: kvDefInt(result, p->key, (int)(p->v.l)); break;
    }
  }

  kvDestroyIterator(kvi);
  return result;
}

void kvCopyUniqueExceptHashes(KVHash* to, KVHash* from)
{
  /* Copy all unique entries in the "from" hash to the "to" hash,
   * except for entries which are themselves hashes.  We can't
   * copy them because we don't want to clone them and a simple
   * copy would mess up memory management.
   */
  KVIterator* kvi= kvUniqueIteratorFactory(from);

  while (kvIteratorHasMorePairs(kvi)) {
    KVPair* p= kvIteratorNextPair(kvi);
    switch (p->type) {
    case KV_STRING: kvDefString(to, p->key, p->v.s); break;
    case KV_HASH: /* Ignore these */ break;
    case KV_LONG: kvDefLong(to, p->key, p->v.l); break;
    case KV_DOUBLE: kvDefDouble(to, p->key, p->v.d); break;
    case KV_BOOLEAN: kvDefBoolean(to, p->key, p->v.l); break;
    case KV_INT: kvDefInt(to, p->key, (int)(p->v.l)); break;
    }
  }

  kvDestroyIterator(kvi);
}

long kvGetNumUniqueEntries(KVHash* kvh)
{
  long count= 0;
  KVIterator* kvi= kvUniqueIteratorFactory(kvh);
  while (kvIteratorHasMorePairs(kvi)) {
    count++;
    kvIteratorNextPair(kvi);
  }
  return count;
}
