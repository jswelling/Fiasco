/************************************************************
 *                                                          *
 *  mriu.c                                                  *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 2007 Department of Statistics,         *
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
 *  Original programming by Joel Welling, 7/07              *
 ************************************************************/
/* This module contains supplemental utilities for manipulating 
 * Pgh MRI files.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "slist.h"
#include "misc.h"
#include "stdcrg.h"

typedef struct label_pair_struct {
  char* s1;
  char* s2;
} LabelPair;

/* List of known labels requiring remapping when indices changes.
 * For example, the tags associated with the label 'foo' in chunk
 * 'images' follow the pattern 'images.foo.x.n' where x is a dimension
 * and n is an index.
 */
static char* remappableLabels[]= { "label",(char*)NULL };

static LabelPair* createLabelPair( char* s1, char* s2 )
{
  LabelPair* result= NULL;
  if (!(result=(LabelPair*)malloc(sizeof(LabelPair)))) 
    MALLOC_FAILURE(sizeof(LabelPair),byte);
  result->s1= strdup(s1);
  result->s2= strdup(s2);
  return result;
}

static void destroyLabelPair( LabelPair* lp )
{
  free(lp->s1);
  free(lp->s2);
  free(lp);
}

void mriu_updateLabels(MRI_Dataset* ds, const char* chunk, 
		       char dim, int low, int highplusone, 
		       DIMREMAPFUNC remapFunc, void* hook) 
{
  char buf[128];
  char newbuf[128];
  int candidateLabelIndex;
  int i;
  SList* changeList= NULL;

  changeList= slist_create();
  for (candidateLabelIndex= 0; remappableLabels[candidateLabelIndex]!=NULL;
       candidateLabelIndex++) {
    const char* candidate= remappableLabels[candidateLabelIndex];
    for (i=low; i<highplusone; i++) {
      snprintf(buf,sizeof(buf),"%s.%s.%c.%d",chunk,candidate,dim,i);
      /*
      fprintf(stderr,"mriu_updateLabels: Checking <%s>\n",buf);
      */
      if (mri_has(ds, buf)) {
	int newIndex= remapFunc(i,hook);
	if (newIndex>=0) {
	  snprintf(newbuf,sizeof(newbuf),"%s.%s.%c.%d",
		   chunk,candidate,dim,newIndex);
	  slist_append(changeList,createLabelPair(newbuf,
						  mri_get_string(ds,buf)));
	}
	/*
	fprintf(stderr,"mriu_updateLabels: Deleting <%s>:<%s>\n",
		buf,mri_get_string(ds,buf));
	*/
	mri_remove(ds,buf);
      }
    }
  }
  slist_totop(changeList);
  while (!slist_atend(changeList)) {
    LabelPair* lp= (LabelPair*)slist_next(changeList);
    /*
    fprintf(stderr,"mriu_updateLabels: Adding <%s>:<%s>\n",lp->s1,lp->s2);
    */
    mri_set_string(ds,lp->s1,lp->s2);
    destroyLabelPair(lp);
  }
  slist_destroy(changeList,NULL);
}


