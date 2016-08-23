/************************************************************
 *                                                          *
 *  dicom_uid_dict.c                                          *
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
/* Notes-
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "mri.h"
#if (SGI64 || SGI5 || SGIMP)
#include <bstring.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

#include "bio.h"
#include "fmri.h"
#include "fexceptions.h"
#include "array.h"
#include "stdcrg.h"
#include "misc.h"
#include "smartreader.h"

/* DICOM header files */
#include "dicom_parser.h"
#include "dicom_uid_dict.h"

static char rcsid[] = "$Id: dicom_uid_dict.c,v 1.3 2004/11/02 22:24:25 welling Exp $";

/* For use in recognizing DICOM files */
static DicomDictDataElement RecognitionElement= 
  { 0x8, 0x16, UI, "SOPClassUID" };
#define DICOM_UID_STRING "1.2.840.10008.5.1.4."

/* We will be needing a place to hang a sorted table of UIDs */
static UID* uidDictSortedByUIDString= NULL;

static int uidSearchByNameTest( const void* v1, const void* v2 )
{
  UID* u1= (UID*)v1;
  UID* u2= (UID*)v2;
  return strcmp(u1->name,u2->name);
}

static int uidSearchByUIDStringTest( const void* v1, const void* v2 )
{
  UID* u1= (UID*)v1;
  UID* u2= (UID*)v2;
  return strcmp(u1->uid,u2->uid);
}

const UID* dcm_getUIDByName( const char* name )
{
  UID* result= NULL;
  UID dummy;

  /* We need to assemble a dummy UID because a UID is what is being
   * searched for in the table.
   */
  dummy.name= name;
  dummy.uid= NULL;
  result= (UID*)bsearch(&dummy, uidDict,
		  sizeof(uidDict)/sizeof(UID),
		  sizeof(UID),
		  uidSearchByNameTest);
  return result;
}

const UID* dcm_getUIDByUIDString( const char* uidString )
{
  UID* result= NULL;
  UID dummy;

  if (!uidDictSortedByUIDString) {
    /* First call- build the sorted table */
    if (!(uidDictSortedByUIDString=(UID*)malloc(sizeof(uidDict))))
      Abort("%s: unable to allocate %d bytes!\n",
	    sizeof(uidDict));
    memcpy(uidDictSortedByUIDString, uidDict, sizeof(uidDict));
    qsort(uidDictSortedByUIDString,
	  sizeof(uidDict)/sizeof(UID),
	  sizeof(UID),
	  uidSearchByUIDStringTest);
  }

  /* We need to assemble a dummy UID because a UID is what is being
   * searched for in the table.
   */
  dummy.name= NULL;
  dummy.uid= uidString;
  result= (UID*)bsearch(&dummy, uidDictSortedByUIDString,
		  sizeof(uidDict)/sizeof(UID),
		  sizeof(UID),
		  uidSearchByUIDStringTest);
  return result;
}

