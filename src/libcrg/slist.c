/************************************************************
 *                                                          *
 *  slist.c                                                 *
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
/* A simple rooted singly-connected list package.  Note that it does
 * not allocate separate copies of the things stored in it!
 * Null data elements will work.
 */

#include <stdio.h>
#include <stdlib.h>
#include <errors.h>
#include "stdcrg.h"
#include "slist.h"

static char rcsid[] = "$Id: slist.c,v 1.7 2005/02/11 01:23:02 welling Exp $";

typedef struct scell_struct {
  void* data;
  struct scell_struct* next;
} SCell;

struct slist_struct {
  SCell* head;
  SCell* tail;
  SCell* current;
};

int slist_empty(SList* list) {
  return(list->head == NULL);
}

SList* slist_create()
{
  SList* result;
  if (!(result= (SList*)malloc(sizeof(SList)))) {
    Abort("slist_create: unable to allocate %d bytes!\n",sizeof(SList));
  }
  result->head= result->tail= result->current= NULL;
  return result;
}

void slist_destroy(SList* list, void (*destructor)(void*))
{
  void* d;
  while (!slist_empty(list)) {
    d= slist_pop(list);
    if (destructor != NULL) (*destructor)(d);
  }
}

void slist_append(SList* list, void* d_in) {
  SCell* cell;
  if (!(cell= (SCell*)malloc(sizeof(SCell)))) {
    Abort("slist_create: unable to allocate %d bytes!\n",sizeof(SCell));
  }
  cell->data= d_in;
  cell->next= NULL;
  if (slist_empty(list)) {
    list->head= list->current= cell;
  }
  else list->tail->next= cell;
  list->tail= cell;
}

void slist_push(SList* list, void* d_in) {
  SCell* cell;
  if (!(cell= (SCell*)malloc(sizeof(SCell)))) {
    Abort("slist_create: unable to allocate %d bytes!\n",sizeof(SCell));
  }
  cell->data= d_in;
  cell->next= list->head;
  if (slist_empty(list)) {
    list->tail= cell;
    list->current= cell;
  }
  list->head= cell;
}

void* slist_pop(SList* list) {
  if (!slist_empty(list)) {
    void* result= list->head->data;
    SCell* old= list->head;
    list->head= list->head->next;
    if (!list->head) {
      list->tail= list->current= NULL;
    }
    if (list->current==old) list->current= list->current->next;
    free((void*)old);
    return result;
  }
  else return NULL;
}

void* slist_get(SList* list) {
  if (list->current) return list->current->data;
  else return NULL;
}

void* slist_getlast(SList* list) {
  if (list->tail) return list->tail->data;
  else return NULL;
}

void* slist_next(SList* list) {
  if (list->current) {
    void* result= list->current->data;
    list->current= list->current->next;
    return result;
  }
  else return NULL;
}

void slist_totop(SList* list) {
  list->current= list->head;
}

int slist_attop(SList* list) {
  return (list->current==list->head);
}

int slist_atend(SList* list) {
  return (list->current==NULL);
}

int slist_count(SList* list) {
  SCell* here= list->head;
  int count= 0;
  while (here!=NULL) {
    count++;
    here= here->next;
  }
  return count;
}

void slist_sort(SList* l, int (*compare)(const void*,const void*))
{
  void** ptbl= NULL;
  int nItems= slist_count(l);
  int i;

  if (!(ptbl=(void**)malloc(nItems*sizeof(void*))))
    Abort("%s: unable to allocate %d bytes!\n",nItems*sizeof(void*));

  i= 0;
  while (!slist_empty(l))
    ptbl[i++]= slist_pop(l);

  qsort(ptbl, nItems, sizeof(void*), compare);

  for (i=0; i<nItems; i++)
    slist_append(l, ptbl[i]);

  slist_totop(l);

  free(ptbl);
}

