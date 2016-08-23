/************************************************************
 *                                                          *
 *  history.c                                               *
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
 *  Original programming by Joel Welling, 10/00             *
 ************************************************************/
/* This module contains routines for manipulating the history info
 * in Pgh MRI files.
 */

#include <string.h>
#include <stdio.h>
#include <time.h>
#include "mri.h"
#include "fmri.h"

const char* hist_get( MRI_Dataset* ds, int i )
{
  char keybuf[64];
  sprintf(keybuf,"history.%d",i);
  if (mri_has(ds,keybuf)) return mri_get_string(ds,keybuf);
  sprintf(keybuf,"history.%03d",i);
  if (mri_has(ds,keybuf)) return mri_get_string(ds,keybuf);
  else return NULL;
}

void hist_add( MRI_Dataset* ds, const char* s )
{
  int i;
  char keybuf[64];

  i= 1;
  while (hist_get(ds,i)) { i++; }; /* walk i up to first free value */
  sprintf(keybuf,"history.%03d",i);
  mri_set_string(ds,keybuf,s);
}

void hist_delete_all( MRI_Dataset* ds )
{
  int i;
  char keybuf[64];

  for (i=1; hist_get( ds, i ); i++) {
    sprintf(keybuf,"history.%d",i);
    if (mri_has(ds,keybuf)) mri_remove(ds,keybuf);
    else {
      sprintf(keybuf,"history.%03d",i);
      mri_remove(ds,keybuf);
    }
  }
}

void hist_delete_some( MRI_Dataset* ds, int n_to_delete )
{
  int i;
  const char* s;
  char keybuf[64];

  for (i=1; i<=n_to_delete && hist_get( ds, i ); i++) {
    sprintf(keybuf,"history.%d",i);
    if (mri_has(ds,keybuf)) mri_remove(ds,keybuf);
    else {
      sprintf(keybuf,"history.%03d",i);
      mri_remove(ds,keybuf);
    }
  }

  /* This last bit relies on the fact that hist_add always starts at 1 */
  for ( /* continue iteration */; (s=hist_get( ds,i ))!=NULL; i++) {
    hist_add(ds, s);
    sprintf(keybuf,"history.%d",i);
    if (mri_has(ds,keybuf)) mri_remove(ds,keybuf);
    else {
      sprintf(keybuf,"history.%03d",i);
      mri_remove(ds,keybuf);
    }
  }
}

void hist_repair( MRI_Dataset* ds )
{
  int i;

  /* This is an automatic function to repair the history tags
   * of a dataset.  
   */

  /* History numbering has changed, so that the first entry is
   * history.001 rather than just history.1.  Look for the old
   * type and repair them.
   */
  for (i=0; i<100; i++) {
    char buf1[64];
    char buf2[64];
    snprintf(buf1,sizeof(buf1),"history.%d",i);
    snprintf(buf2,sizeof(buf1),"history.%03d",i);
    if (mri_has(ds,buf1)) {
      if (mri_has(ds,buf2))
	Abort("hist_repair: corrupted history! Both %s and %s are defined in %s.\n",
	      MRI_DATASET_NAME(ds));
      mri_set_string(ds,buf2,mri_get_string(ds,buf1));
      mri_remove(ds,buf1);
    }
  }
}

void hist_add_cl( MRI_Dataset* ds, int argc, char** argv )
{
  int word;
  char strbuf[512];
  int tot_length;
  time_t t;
  struct tm* tsp;

  /* Let's take this opportunity to fix the history list if necessary. */
  hist_repair(ds);

  time(&t);
  tsp= localtime(&t);
  tot_length= strftime(strbuf,sizeof(strbuf),"[%F %T] ",tsp);
  for (word= 0; word<argc; word++) {
    int length= strlen(argv[word]);
    if (length+tot_length >= sizeof(strbuf)-2) 
      length= sizeof(strbuf)-(tot_length+2);
    strncat(strbuf,argv[word],length);
    tot_length= strlen(strbuf);
    if (word<argc-1 && tot_length<sizeof(strbuf)-1) {
      strcat(strbuf," ");
      tot_length++;
    }
    if (tot_length>=sizeof(strbuf)-1) break;
  }
  strbuf[sizeof(strbuf)-1]= '\0'; /* just for paranoia */

  hist_add(ds, strbuf);
}

void hist_dump( MRI_Dataset* ds, FILE* ofile )
{
  int i; 
  const char* s;

  i= 1;
  while (s=hist_get(ds,i++)) { fputs(s,ofile); fputc('\n',ofile); }
}
