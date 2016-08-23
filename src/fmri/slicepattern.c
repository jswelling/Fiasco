/************************************************************
 *                                                          *
 *  slicepattern.c                                             *
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
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "mri.h"
#if (SGI64 || SGI5 || SGIMP)
#include <bstring.h>
#endif
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

#include "bio.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"
#include "slicepattern.h"

static char rcsid[] = "$Id: slicepattern.c,v 1.2 2006/03/16 23:04:29 welling Exp $";

static const char* knownSlicePatternNames[]= {
  "sequential",
  "reversed_sequential",
  "even/odd",
  "reversed_even/odd",
  "odd/even",
  "reversed_odd/even",
  "halves_low_first",
  "reversed_halves_low_first",
  "halves_high_first",
  "reversed_halves_high_first",
  NULL /* must be last */
};

int* slp_generateSlicePatternTable( int dz, const char* name )
{
  int* indexTable= NULL;
  int i;
  int brk;

  if (!(indexTable=(int*)malloc(dz*sizeof(int))))
    Abort("slp_generateSlicePatternTable: unable to allocate %d bytes!\n",
	  dz*sizeof(int));

  if (!strcmp(name,"sequential")) {
    for (i=0; i<dz; i++) {
      int j= i;
      indexTable[j]=i;
    }
  }
  else if (!strcmp(name,"reversed_sequential")) {
    for (i=0; i<dz; i++) {
      int j= dz-(i+1);
      indexTable[j]=i;
    }
  }
  else if (!strcmp(name,"even/odd")) {
    if (dz%2) brk= (dz-1)/2;
    else brk= (dz-2)/2;
    for (i=0; i<dz; i++) {
      int j= i;
      if (i<=brk) {
	indexTable[j]=2*i;
      }
      else {
	indexTable[j]=2*(i-brk)-1;
      }
    }
  }
  else if (!strcmp(name,"reversed_even/odd")) {
    if (dz%2) brk= (dz-1)/2;
    else brk= (dz-2)/2;
    for (i=0; i<dz; i++) {
      int j= dz-(i+1);
      if (i<=brk) {
	indexTable[j]=2*i;
      }
      else {
	indexTable[j]=2*(i-brk)-1;
      }
    }
  }
  else if (!strcmp(name,"odd/even")) {
    if (dz%2) brk= (dz-3)/2;
    else brk= (dz-2)/2;
    for (i=0; i<dz; i++) {
      int j= i;
      if (i<=brk) {
	indexTable[j]=2*i+1;
      }
      else {
	indexTable[j]=2*(i-(brk+1));
      }
    }
  }
  else if (!strcmp(name,"reversed_odd/even")) {
    if (dz%2) brk= (dz-3)/2;
    else brk= (dz-2)/2;
    for (i=0; i<dz; i++) {
      int j= dz-(i+1);
      if (i<=brk) {
	indexTable[j]=2*i+1;
      }
      else {
	indexTable[j]=2*(i-(brk+1));
      }
    }
  }
  else if (!strcmp(name,"halves_low_first")) {
    if (dz%2) brk= (dz+1)/2;
    else brk= dz/2;
    for (i=0; i<dz; i++) {
      int j= i;
      if (i%2) {
	indexTable[j]= (i-1)/2 + brk;
      }
      else {
	indexTable[j]= i/2;
      }
    }
  }
  else if (!strcmp(name,"reversed_halves_low_first")) {
    if (dz%2) brk= (dz+1)/2;
    else brk= dz/2;
    for (i=0; i<dz; i++) {
      int j= dz-(i+1);
      if (i%2) {
	indexTable[j]= (i-1)/2 + brk;
      }
      else {
	indexTable[j]= i/2;
      }
    }
  }
  else if (!strcmp(name,"halves_high_first")) {
    if (dz%2) brk= (dz-3)/2;
    else brk= (dz-2)/2;
    for (i=0; i<dz; i++) {
      int j= i;
      if (i%2) {
	indexTable[j]= (i-1)/2;
      }
      else {
	indexTable[j]= (i/2)+brk+1;
      }
    }
  }
  else if (!strcmp(name,"reversed_halves_high_first")) {
    if (dz%2) brk= (dz-3)/2;
    else brk= (dz-2)/2;
    for (i=0; i<dz; i++) {
      int j= dz-(i+1);
      if (i%2) {
	indexTable[j]= (i-1)/2;
      }
      else {
	indexTable[j]= (i/2)+brk+1;
      }
    }
  }
  else {
    free(indexTable);
    indexTable= NULL; /* to signal error */
  }

  return indexTable;
}

const char* slp_findSlicePatternNameFromTable(int dz, const int* table)
{
  const char* name= NULL;
  int i= 0;

  i= 0;
  while (name= knownSlicePatternNames[i]) {
    int* stdTable= slp_generateSlicePatternTable(dz, name);
    int j;
    int failed= 0;
    for (j=0; j<dz; j++) {
      if (stdTable[j] != table[j]) {
	failed= 1;
	break;
      }
    }
    free(stdTable);
    if (!failed) break;
    i++;
  }

  return name;
}

int* slp_generateInvertedSlicePatternTable( int dz, const char* name )
{
  int* directTable= slp_generateSlicePatternTable(dz,name);
  int* table= slp_invertSlicePatternTable(dz, directTable);

#ifdef never
  {
   int i;
   for (i=0; i<dz; i++)
     fprintf(stderr,"%d: %d %d dblmap %d\n",i,directTable[i],table[i],directTable[table[i]]);
 }
#endif

  free(directTable);
  return table;
}

int* slp_invertSlicePatternTable(int dz, const int* directTable)
{
  int* table= NULL;
  int i;

  if (!(table=(int*)malloc(dz*sizeof(int)))) 
    Abort("slicepattern: unable to allocate %d bytes!\n",dz*sizeof(int));

  for (i=0; i<dz; i++) {
    table[directTable[i]]= i;
  }

  return table;
}
