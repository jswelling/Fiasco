/************************************************************
 *                                                          *
 *  parsesplit.c                                             *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1997 Department of Statistics,         *
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
 *  Original programming by Joel Welling, 8/98              *
 ************************************************************/
/* This module contains routines for parsing split and condition files */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "fmri.h"
#include "misc.h"
#include "slist.h"

#define INBUF_LENGTH 512
#define FACTOR_NAME_ABBREV_LENGTH 5

typedef struct factor_info_struct {
  char* name;
  KVHash* levels;
  int nLevels; /* for convenience */
} FactorInfo;

static sp_ConditionDef** cond_def_table= NULL;
static int cond_def_table_size= 0;

static int is_NA( char* str )
{
  if (!strcmp(str,"NA") || !strcmp(str,"na") || !strcmp(str,"missing")
      || !strcmp(str,"MISSING")) return 1;
  else return 0;
}

static int cond_rec_compare( const void* r1_in, const void* r2_in )
{
  sp_ConditionDef** r1= (sp_ConditionDef**)r1_in;
  sp_ConditionDef** r2= (sp_ConditionDef**)r2_in;

  if ((*r1)->id < (*r2)->id) return -1;
  else if ((*r1)->id > (*r2)->id) return 1;
  else return 0;
}

static int process_cond_file( FILE* cond_file, 
			      char*** factor_table, int* nfactors, 
			      sp_ConditionDef*** cond_table, int* nconditions )
{
  sp_ConditionDef** cond_def_table= NULL;
  sp_ConditionDef* thisCond= NULL;
  char inbuf[INBUF_LENGTH];
  char condName[INBUF_LENGTH];
  char** factor_lvl_names= NULL;
  int* factor_lvl_counts= NULL;
  int i;
  int j;
  char* tmpp;
  SList* factors= NULL;
  SList* conditions= NULL;
  KVHash* condNameHash= NULL;
  KVIterator* kvi= NULL;
  char* here;
  int linenum= 0;
  int condnum= -1;
  int totalCrosses;
  
  /* Read header line */
  if (!fgets(inbuf,INBUF_LENGTH,cond_file)) {
    Error("parsesplit: failed to read conditions file: %s!\n",strerror(errno));
    return 0;
  }
  inbuf[INBUF_LENGTH-1]= '\0';
  linenum++;
  
  /* Count factors, recording names */
  factors= slist_create();
  conditions= slist_create();
  condNameHash= kvFactory(KV_DEFAULT_SIZE);
  here= strtok_r(inbuf," \t\n",&tmpp);
  while (here) {
    FactorInfo* fi;
    if (!(fi=(FactorInfo*)malloc(sizeof(FactorInfo))))
      Abort("parsesplit: unable to allocate %d bytes!\n",sizeof(FactorInfo));
    fi->name= strdup(here);
    fi->levels= kvFactory(KV_DEFAULT_SIZE);
    fi->nLevels= 0;
    slist_append(factors,fi);
    here= strtok_r(NULL," \t\n",&tmpp);
  }
  *nfactors= slist_count(factors);
  
  /* Read conditions, building factor level and name info */
  condnum= -1;
  while (!feof(cond_file) && !ferror(cond_file)) {
    if (!fgets(inbuf,INBUF_LENGTH,cond_file)) break;
    inbuf[INBUF_LENGTH-1]= '\0';
    linenum++;
    if (ferror(cond_file)) break;
    here= strtok_r(inbuf," \t\n",&tmpp);
    if (here) {
      int thiscond= atoi(here);
      int* levelValues;
      char** levelNames;
      if (thiscond<=condnum) 
	Abort("parsesplit: condition numbers invalid or not sorted!\n");
      condnum= thiscond;
      if (!(levelValues=(int*)malloc(*nfactors*sizeof(int))))
	Abort("parsesplit: unable to allocate %d bytes!\n",
	      *nfactors*sizeof(int));
      if (!(levelNames=(char**)malloc(*nfactors*sizeof(char*))))
	Abort("parsesplit: unable to allocate %d bytes!\n",
	      *nfactors*sizeof(char*));
      if (!(thisCond= (sp_ConditionDef*)malloc(sizeof(sp_ConditionDef))))
	Abort("parsesplit: unable to allocate %d bytes!\n",
	      sizeof(sp_ConditionDef));
      thisCond->factor_lvl= levelNames;
      thisCond->factor_intlvl= levelValues;
      slist_totop(factors);
      condName[0]= '\0';
      i= 0;
      here= strtok_r(NULL," \t\n",&tmpp);
      while (!slist_atend(factors)) {
	FactorInfo* fi= (FactorInfo*)slist_next(factors);
	if (!here) 
	  Abort("parsesplit: condition line %d has too few factor levels!\n",
		linenum);
	if ( !kvLookup(fi->levels,here) ) {
	  kvDefInt(fi->levels,here,fi->nLevels++);
	}
	levelValues[i]= kvGetInt(fi->levels, here);
	levelNames[i]= strdup(here);
	/* This can't overflow; shorter than inbuf */
	strncat(condName, here, FACTOR_NAME_ABBREV_LENGTH); 
	strcat(condName, "-");
	here= strtok_r(NULL," \t\n",&tmpp); 
	i++;
      }
      condName[strlen(condName)-1]= '\0'; /* remove trailing blank */
      if (!kvLookup(condNameHash, condName))
	kvDefInt(condNameHash, condName, condnum);
      thisCond->name= strdup(condName);
      thisCond->id= condnum;
      slist_append(conditions, thisCond);
    }
    else {
      /* skip this empty line */
    }
  }
  if (ferror(cond_file)) {
    Error("parsesplit: Error reading conditions file: %s!\n",strerror(errno));
    return 0;
  }

  *nconditions= slist_count(conditions);


  /* Issue a warning if some possible crosses don't appear in
   * the condition table.
   */
  totalCrosses= 1;
  slist_totop(factors);
  while (!slist_atend(factors)) {
    FactorInfo* fi= (FactorInfo*)(slist_next(factors));
    totalCrosses *= fi->nLevels;
  }
  if (totalCrosses> *nconditions)
    Warning(1,"parsesplit: some factor crosses never occur!\n");

  /* Build the table of factors */
  if (!(*factor_table= (char**)malloc((*nfactors)*sizeof(char*))))
    Abort("parsesplit: failed to allocate %d bytes!\n", 
	  (*nfactors)*sizeof(char*));
  slist_totop(factors);
  i= 0;
  while (!slist_atend(factors)) {
    FactorInfo* fi= (FactorInfo*)slist_pop(factors);
    (*factor_table)[i++]= fi->name;
    kvDestroy(fi->levels);
    free(fi);
  }
  slist_destroy(factors,NULL);
  
  /* build the table of conditions */
  if (!(cond_def_table=  
	(sp_ConditionDef**)malloc(*nconditions*sizeof(sp_ConditionDef*))))
    Abort("parsesplit: unable to allocate %d bytes!\n",
	  *nconditions * sizeof(sp_ConditionDef*));
  *cond_table= cond_def_table;
  slist_totop(conditions);
  i= 0;
  while (!slist_empty(conditions))
    cond_def_table[i++]= (sp_ConditionDef*)slist_pop(conditions);
  slist_destroy(conditions, NULL);

  /* Deal with possibility that the file presented conditions in
     the wrong order, by sorting on id */
  qsort((void*)cond_def_table, *nconditions, sizeof(sp_ConditionDef*), 
	cond_rec_compare);
  for (i=0; i<*nconditions; i++)
    if (cond_def_table[i]->id != i) {
      Error("parsesplit: condition data is not sequential from 0!\n");
      return 0;
    }

  return 1;
}

int sp_parse_conditions( char* fname, 
			 char*** factor_table, int* nfactors,
			 sp_ConditionDef*** cond_table, int* nconditions )
{
  FILE* cond_file;

  if (!(cond_file= fopen(fname,"r"))) {
    Error("Unable to open condition file <%s>!\n",fname);
    return 0;
  }
  if (!process_cond_file(cond_file, factor_table, nfactors,
			 cond_table, nconditions)) {
    Error("Error processing conditions file <%s>!\n",fname);
    return 0;
  }
  if (fclose(cond_file))
    Error("Failed to close conditions file <%s>!\n",fname);

  return 1;
}

static int split_rec_compare( const void* r1_in, const void* r2_in )
{
  sp_SplitRec* r1= (sp_SplitRec*)r1_in;
  sp_SplitRec* r2= (sp_SplitRec*)r2_in;

  if (r1->image < r2->image) return -1;
  else if (r1->image > r2->image) return 1;
  else if (r1->slice < r2->slice) return -1;
  else if (r1->slice > r2->slice) return 1;
  else return 0;
}

static sp_SplitRec* process_split_file( FILE* split_file, 
					int nslices, int nimages,
					int nconditions )
{
  int i;
  sp_SplitRec* split_table= NULL;

  if (!(split_table= 
	(sp_SplitRec*)malloc(nimages*nslices*sizeof(sp_SplitRec))))
    Abort("parsesplit: failed to allocate %d bytes!\n",
	  nimages*nslices*sizeof(sp_SplitRec));

  /* Load the split records */
  i= 0;
  while (!feof(split_file) && !ferror(split_file)) {
    if ( fscanf(split_file,"%d %d %d\n", &(split_table[i].image), 
		&(split_table[i].slice), &(split_table[i].cond)) != 3 ) {
      Error("parsesplit: error parsing split file!\n");
      return NULL;
    }
    if ((split_table[i].image<0) || (split_table[i].image>=nimages)
	|| (split_table[i].slice<0) || (split_table[i].slice>=nslices)
	|| (split_table[i].cond<0) 
	|| (split_table[i].cond>=nconditions)) {
      Error("parsesplit: bad split table line: <%d %d %d>\n",
	    split_table[i].image,split_table[i].slice,
	    split_table[i].cond);
      return NULL;
    }
    i++;
  }

  if (ferror(split_file)) {
    Error("parsesplit: error reading split file!\n");
    return NULL;
  }

  /* Now cleverly sort the records by image and slice */
  qsort((void*)split_table, nimages*nslices, sizeof(sp_SplitRec), 
	split_rec_compare);

  return split_table;
}

int sp_parse_split( char* fname, int nslices, int nimages, int nconditions,
		    sp_SplitRec** split_table )
{
  FILE* split_file;

  if (!(split_file= fopen(fname,"r"))) {
    Error("parsesplit: unable to open split file <%s>!\n",fname);
    return 0;
  }
  if (!(*split_table= 
	process_split_file(split_file, nslices, nimages, nconditions))) {
    Error("parsesplit: error processing split file <%s>!\n",fname);
    return 0;
  }
  if (fclose(split_file))
    Error("parsesplit: failed to close split file <%s>!\n",fname);

  return 1;
}

int sp_match_missing( sp_SplitRec* split_table, unsigned char** missing,
		      int nslices, int nimages )
{
  int t;
  int z;

  for (t=0; t<nimages; t++) 
    for (z=0; z<nslices; z++) {
      int offset= (t*nslices) + z;
      if ((split_table[offset].slice != z) 
	  || (split_table[offset].image != t)) {
	Error("parsesplit: split table is disordered in sp_match_missing!\n");
	return 0;
      }
      if (missing[t][z]) split_table[offset].cond= 0;
      else if (split_table[offset].cond==0) missing[t][z]= 1;
    }
  return 1;
}
