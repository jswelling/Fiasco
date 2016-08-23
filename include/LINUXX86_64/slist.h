/************************************************************
 *                                                          *
 *  slist.h                                                 *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1999 Department of Statistics,         *
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
 *  Original programming by Joel Welling 4/1999             *
 ************************************************************/
/* Interface to slist.c */
/* A simple rooted singly-connected list package.  Note that it does
 * not allocate separate copies of the things stored in it!
 * Null data elements will work.
 */

#ifndef __SLIST_INC__
#define __SLIST_INC__

typedef struct slist_struct SList;

SList* slist_create(void);
void slist_destroy(SList* list, void (*destructor)(void*));
void slist_append(SList* list, void* d_in);
void slist_push(SList* list, void* d_in);
void* slist_pop(SList* list);
void* slist_get(SList* list);
void* slist_getlast(SList* list);
void* slist_next(SList* list);
void slist_totop(SList* list);
int slist_count(SList* list);
int slist_empty(SList* list);
int slist_attop(SList* list);
int slist_atend(SList* list);
void slist_sort(SList* list, int (*compare)(const void*,const void*));

#endif


