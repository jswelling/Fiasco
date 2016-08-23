/************************************************************
 *                                                          *
 *  kvhash.h                                                *
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
 ************************************************************/

#ifndef INC_KVHASH_H
#define INC_KVHASH_H

/*
 * Structs for simple hash table facility
 */
typedef enum { 
  KV_STRING, KV_LONG, KV_DOUBLE, KV_HASH, KV_BOOLEAN, KV_INT
} KVType;

#define KV_DEFAULT_SIZE 257

struct kv_hash_struct; /* forward definition */

typedef struct kv_pair_struct {
  char* key;
  KVType type;
  union {
    char* s;
    struct kv_hash_struct* h;
    long long l;
    double d;
  } v;
  struct kv_pair_struct* next;
  struct kv_pair_struct* prev;
} KVPair;

typedef struct kv_hash_struct {
  int size; /* number of lists, defined at creation time */
  KVPair** table;
} KVHash;

typedef struct kv_iterator_struct {
  KVHash* hash;
  int id;
  KVPair* p;
  void (*stepToNext)(struct kv_iterator_struct*);
  void* hook;
  int showOnlyUnique;
} KVIterator;

/* Hash table methods */
extern KVHash *kvFactory(int size); /* size should be prime */
extern KVHash *kvCloneUnique(KVHash* original);
extern void kvCopyUniqueExceptHashes(KVHash* to, KVHash* from);
extern void kvDestroy( KVHash* kvh );
extern KVPair *kvLookup( KVHash* kvh, const char* key );
extern KVType kvType( const KVPair* p );
extern const char* kvKey( const KVPair* p );
extern const char* kvTypeName( KVType t );
extern void kvDelete( KVHash* kvh, const char* key );
extern void kvDeleteAll( KVHash* kvh, const char* key );
extern void kvDefString( KVHash* kvh, const char* key, const char* value );
extern void kvDefLong( KVHash* kvh, const char* key, const long long value );
extern void kvDefDouble( KVHash* kvh, const char* key, const double value );
extern void kvDefHash( KVHash* kvh, const char* key, KVHash* value );
extern void kvDefBoolean( KVHash* kvh, const char* key, int value );
extern void kvDefInt( KVHash* kvh, const char* key, int value );
extern char* kvGetString( KVHash* kvh, const char* key );
extern long long kvGetLong( KVHash* kvh, const char* key );
extern double kvGetDouble( KVHash* kvh, const char* key );
extern KVHash* kvGetHash( KVHash* kvh, const char* key );
extern int kvGetBoolean( KVHash* kvh, const char* key );
extern int kvGetInt( KVHash* kvh, const char* key );
extern long kvGetNumUniqueEntries( KVHash* kvh );
extern KVIterator* kvIteratorFactory( KVHash* kvh );
extern KVIterator* kvUniqueIteratorFactory( KVHash* kvh );
extern KVIterator* kvSortedIteratorFactory( KVHash* kvh );
extern void kvDestroyIterator( KVIterator* kvi );
extern KVPair* kvIteratorNextPair( KVIterator* kvi );
extern int kvIteratorHasMorePairs( KVIterator* kvi );
extern void kvDumpTableUnique( KVHash* kvh, FILE* ofile, const char* prefix,
			       int maxDepth, int depth );

#endif /* ifdef INC_KVHASH_H */
