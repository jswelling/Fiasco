/************************************************************
 *                                                          *
 *  wildcard.c                                         *
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
 *   Original programming by Joel Welling 6-02              *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "mri.h"
#if (SGI64 || SGI5 || SGIMP)
#include <bstring.h>
#endif
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>

#include "bio.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"
#include "smartreader.h"

static char rcsid[] = "$Id: wildcard.c,v 1.9 2007/03/22 00:03:47 welling Exp $";
#if ( SUN4SOL2 || DARWIN )
/* This foolish BSD machine does not supply scandir, so we will include
 * code to do so.  This version of scandir for Suns was written by 
 * joehill@users.sourceforge.net and is believed to be open source.
 * See http://archives.seul.org/gdis/gdiscuss/Apr-2001/msg00002.html .
 */

/* 
 * Scan the directory dirname calling select to make a list of selected 
 * directory entries then sort using qsort and compare routine dcomp. 
 * Returns the number of entries and a pointer to a list of pointers to 
 * struct dirent (through namelist). Returns -1 if there were any errors. 
 */ 

/*******************************************************************************
*                                                                              *
*                                   Viewmol                                    *
*                                                                              *
*                              S C A N D I R . C                               *
*                                                                              *
*                 Copyright (c) Joerg-R. Hill, December 2000                   *
*                                                                              *
********************************************************************************
*
* $Id: wildcard.c,v 1.9 2007/03/22 00:03:47 welling Exp $
* $Log: wildcard.c,v $
* Revision 1.9  2007/03/22 00:03:47  welling
* Mods to get things to compile on DARWIN, plus (optional) support for SON
* file format.
*
* Revision 1.8  2005/02/05 02:56:53  welling
* Added '^' as an equivalent to '*', to avoid stupid csh-expansion problems.
*
* Revision 1.7  2005/01/20 23:55:05  welling
* Fixed problem in handling '*' scanning backwards.
*
* Revision 1.6  2003/07/01 20:38:28  welling
* Try to fix a type check warning in the fake scandir included for the
* benefit of SUN4SOL2.
*
* Revision 1.5  2003/04/24 22:17:16  welling
* Improvements to handling of y dimension to better support partial-k; modified
* -multi wildcard support to handle '*' and '?'.
*
* Revision 1.4  2003/01/21 19:18:24  welling
* Mods to wildcard.c to support Suns, which don't seem to provide scandir().
*
* Revision 1.2  2000/12/10 15:37:02  jrh
* Release 2.3
*
* Revision 1.1  1999/05/24 01:29:43  jrh
* Initial revision
*
*/


/* This function is only required for SunOS, all other supported OS
   have this function in their system library */

int scandir(const char *dir, struct dirent ***namelist,
            int (*select)(const struct dirent *),
            int (*compar)(const void*, const void*))
{
  DIR *d;
  struct dirent *entry;
  register int i=0;
  size_t entrysize;

  if ((d=opendir(dir)) == NULL)
     return(-1);

  *namelist=NULL;
  while ((entry=readdir(d)) != NULL)
  {
    if (select == NULL || (select != NULL && (*select)(entry)))
    {
      *namelist=(struct dirent **)realloc((void *)(*namelist),
                 (size_t)((i+1)*sizeof(struct dirent *)));
        if (*namelist == NULL) return(-1);
        entrysize=sizeof(struct dirent)-sizeof(entry->d_name)+strlen(entry->d_name)+1;
        (*namelist)[i]=(struct dirent *)malloc(entrysize);
        if ((*namelist)[i] == NULL) return(-1);
        memcpy((*namelist)[i], entry, entrysize);
        i++;
    }
  }
  if (closedir(d)) return(-1);
  if (i == 0) return(-1);
  if (compar != NULL)
    qsort((void *)(*namelist), (size_t)i, sizeof(struct dirent *), compar);
    
  return(i);
}

int alphasort(const void *a, const void *b)
{
  const struct dirent** ad= (const struct dirent**)a;
  const struct dirent** bd= (const struct dirent**)b;
  return(strcmp((*ad)->d_name, (*bd)->d_name));
}

#endif


/* We want matches to work in both directions, so that a string
 * must match both forwards and backwards.
 */
static int patternMatchBackwards( const char* s, const char* p )
{
  const char* rs= s + strlen(s) - 1;
  const char* rp= p + strlen(p) - 1;
  while (rs > s && rp > p) {
    if (*rs == *rp || *rp == '?' || *rp == '#') {
      rs--;
      rp--;
    }
    else if (*rp == '*' || *rp == '^') {
      rs--;
    }
    else return 0;
  }
  if ((*rs == *rp) || (*rp=='*') || (*rp=='^')) return 1;
  else return 0;
}

static int patternMatch(const char* s_orig, const char* p_orig)
{
  const char* s= s_orig;
  const char* p= p_orig;
  while (*s && *p) {
    if (*s == *p || *p == '?' || *p == '#') {
      s++;
      p++;
    }
    else if (*p == '*' || *p == '^') {
      s++;
    }
    else return 0;
  }
  if (!(*s)  && (!(*p) || (*p=='*') || (*p=='^')))
    return patternMatchBackwards(s_orig,p_orig);
  else return 0;
}

static void splitFilename( const char* fname, 
			   char** dirname, char** rootname )
{
  char* here= NULL;

  /* We need to separate the directory from the root file path */
  *dirname= strdup(fname);
  if ((here=strrchr(*dirname,'/')) != NULL) {
    *rootname= strdup(here+1);
    *here= '\0';
  }
  else {
    *rootname= *dirname;
    *dirname= strdup(".");
  }
}

static char* safe_concat( const char* s1, const char* s2, const char* s3 )
{
  char* result;
  int length= strlen(s1) + strlen(s2) + strlen(s3) + 1;

  if (!(result= (char*)malloc(length)))
    Abort("%s: unable to allocate %d bytes!\n",progname);

  strcpy( result, s1 );
  strcat( result, s2 );
  strcat( result, s3 );
  return result;
}

static int getDirEntries( char* dirname, struct dirent*** namelist )
{
  int numEntries= 0;
  
  if ((numEntries=scandir(dirname, namelist, NULL, alphasort))<0) {
    Abort("%s: Error scanning directory <%s> for multiple files: %s!\n",
	  progname,dirname,strerror(errno));
  }

  return numEntries;
}

char* expandWildcardFirstFilename( char* fname )
{
  char* result= NULL;
  char* dirname= NULL;
  char* rootname= NULL;
  struct dirent** namelist;
  int numEntries;
  int i;

  splitFilename( fname, &dirname, &rootname );
  numEntries= getDirEntries( dirname, &namelist );

  /* scandir's select test is not reentrant, so let's use our own test */
  if (numEntries>0) {
    for (i=0; i<numEntries; i++) {
      if (patternMatch(namelist[i]->d_name, rootname)) {
	result= safe_concat(dirname,"/",namelist[i]->d_name);
	break;
      }
    }
  }
  free(namelist);
  free(dirname);
  free(rootname);
  return result;
}

SList* expandWildcardFilename( char* fname )
{
  SList* result= NULL;
  char* dirname= NULL;
  char* rootname= NULL;
  struct dirent** namelist;
  int numEntries;
  int i;

  splitFilename( fname, &dirname, &rootname );
  numEntries= getDirEntries( dirname, &namelist );

  /* scandir's select test is not reentrant, so let's use our own test */
  if (numEntries>0) {
    for (i=0; i<numEntries; i++) {
      if (patternMatch(namelist[i]->d_name, rootname)) {
	if (!result) result= slist_create();
	slist_append(result,safe_concat(dirname,"/",namelist[i]->d_name));
      }
    }
  }
  if (debug) {
    if (result) fprintf(stderr,"expanding <%s> produced %d matching files\n",
			fname, slist_count(result));
    else fprintf(stderr,"expanding <%s> produced no matching files\n",
		 progname);
  }

  free(namelist);
  free(dirname);
  free(rootname);
  return result;
}

void destroyWildcardFilenameList( SList* list )
{
  while (!slist_empty(list)) {
    free( (char*)slist_pop(list) );
  }
}


