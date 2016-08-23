/************************************************************
 *                                                          *
 *  mri_describe.c                                           *
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
 *  Original programming by Joel Welling, 5/98              *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "slist.h"
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"

#define KEYBUF_SIZE 512

static char rcsid[] = "$Id: mri_describe.c,v 1.5 2005/10/04 17:32:44 welling Exp $";

static const char* progname= NULL;
static int verbose_flag= 0;
static int debug_flag= 0;

static int tryGetVec3( MRI_Dataset* ds, const char* chunk, const char* key,
		       double* vec )
{
  char buf[256];
  int i;
  int success= 1;
  for (i=0; i<3; i++) {
    snprintf(buf,sizeof(buf),"%s.%s.%d",chunk,key,i);
    if (!mri_has(ds,buf)) {
      success=0;
      break;
    }
    vec[i]= mri_get_float(ds,buf);
  }
  return success;
}

static int tryGetVec3XYZ( MRI_Dataset* ds, const char* chunk, const char* key,
			  double* vec )
{
  static char* dims= "xyz";
  char buf[256];
  int i;
  int success= 1;
  for (i=0; i<3; i++) {
    snprintf(buf,sizeof(buf),"%s.%s.%c",chunk,key,dims[i]);
    if (!mri_has(ds,buf)) {
      success=0;
      break;
    }
    vec[i]= mri_get_float(ds,buf);
  }
  return success;
}

static void describeChunk( MRI_Dataset* ds, const char* chunk )
{
  char buf[256];
  char dimbuf[256];
  const char* dimstr;
  const char* datatype;
  const char* runner;
  char* dimrunner;

  buf[sizeof(buf)-1]= '\0';
  dimbuf[sizeof(dimbuf)-1]= '\0';
  snprintf(buf,sizeof(buf)-1,"%s.dimensions",chunk);
  if (!mri_has(ds,buf))
    Abort("%s: This dataset has no dimensions for chunk %s!\n",
	  progname,chunk);
  dimstr= mri_get_string(ds,buf);
  runner= dimstr;
  dimrunner= dimbuf;
  while (*runner) {
    snprintf(buf,sizeof(buf)-1,"%s.extent.%c",chunk,*runner);
    if (!mri_has(ds,buf))
      Abort("%s: This dataset has no extent for dimension %c chunk %s!\n",
	    progname,*runner,chunk);
    snprintf(dimrunner,sizeof(dimbuf)-((dimrunner-dimbuf)+1),
	     "%lld",mri_get_int(ds,buf));
    while (*dimrunner) dimrunner++;
    runner++;
    if (*runner) {
      snprintf(dimrunner,sizeof(dimbuf)-((dimrunner-dimbuf)+1),
	       ":");
      while (*dimrunner) dimrunner++;
    }
  }
  printf("   chunk <%s>\n",chunk);
  printf("      dimensions <%s>\n",dimstr);
  printf("      extents %s\n",dimbuf);
  snprintf(buf,sizeof(buf)-1,"%s.datatype",chunk);
  printf("      datatype %s\n",mri_get_string(ds,buf));
  runner= dimstr;
  while (*runner) {
    snprintf(buf,sizeof(buf)-1,"%s.description.%c",chunk,*runner);
    if (mri_has(ds,buf)) {
      printf("      dim %c: %s\n",*runner,mri_get_string(ds,buf));
    }
    else {
      if (verbose_flag)
      printf("      dim %c: (no description)\n",*runner);
    }
    runner++;
  }
  if (verbose_flag) {
    static char* names[]= {"tlf","trf","tlb","trb","blf","brf","blb","brb",
			   "ctr","slice_norm","slice_tlc","slice_trc",
			   "slice_blc", "slice_brc", NULL};
    static char* namesXYZ[]= {"voxel_spacing",NULL};
    double vec[3];
    int i;
    for (i=0; names[i]; i++) {
      if (tryGetVec3(ds, chunk, names[i], vec)) {
	printf("      %s: (%f, %f, %f)\n",names[i],vec[0],vec[1],vec[2]);
      }
    }
    for (i=0; namesXYZ[i]; i++) {
      if (tryGetVec3XYZ(ds, chunk, namesXYZ[i], vec)) {
	printf("      %s: (%f, %f, %f)\n",namesXYZ[i],vec[0],vec[1],vec[2]);
      }
    }
  }
}


static void describeDS( const char* fname )
{
  MRI_Dataset* ds;
  const char* key= NULL;
  SList* chunkList= slist_create();

  printf("Dataset: %s\n",fname);

  /* Open input dataset */
  ds = mri_open_dataset( fname, MRI_READ );
  mri_iterate_over_keys(ds);
  while ( key= mri_next_key(ds) ) {
    if (!strcmp(mri_get_string(ds,key),"[chunk]"))
      slist_append(chunkList,strdup(key));
  }

  slist_totop(chunkList);
  while (!slist_atend(chunkList)) 
    describeChunk(ds, (const char*)slist_next(chunkList));
  slist_destroy(chunkList,free);
  mri_close_dataset(ds);
}

int main( int argc, char* argv[] ) 
{
  char fname[512];
  SList* fnameList= slist_create();
  int printDelimiters= 0;

  progname= argv[0];

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /*** Parse command line ***/

  cl_scan( argc, argv );

  verbose_flag= cl_present("verbose|v");
  debug_flag= cl_present("debug");

  while (cl_get("","%s",fname)) {
    slist_append(fnameList, strdup(fname));
  }
  
  if (cl_cleanup_check()) {
    int i;
    fprintf(stderr,"%s: invalid argument in command line:\n    ",argv[0]);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(-1);
  }

  /*** End command-line parsing ***/

  if (slist_empty(fnameList)) {
    fprintf(stderr,"%s: you must provide at least one dataset name.\n",
	    progname);
  }
  printDelimiters= (slist_count(fnameList)>1);

  slist_totop(fnameList);
  while (!slist_atend(fnameList)) {
    char* thisFname= (char*)slist_next(fnameList);
    if (printDelimiters) printf("##########################\n");
    describeDS(thisFname);
  }

  slist_destroy(fnameList, free);
  return 0;
}

