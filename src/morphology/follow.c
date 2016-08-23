/************************************************************
 *                                                          *
 *  follow.c                                                *
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
 *  Original programming by Mark Fitzgerald  2-95           *
 *     5-96: Pittsburgh Format, Mark Fitzgerald             *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"

static char rcsid[] = "$Id: follow.c,v 1.1 2003/03/18 23:58:49 welling Exp $";

/* Access for 3D arrays */
#define LOC(matrix,x,y,z) matrix[((((z)*dy)+(y))*dx)+(x)]
#define PLOC(matrix,x,y,z) matrix[((((z+1)*(dy+2))+(y+1))*(dx+2))+(x+1)]

static int debug_flag= 0;
static char* progname;

static int opCount= 0;
static int changes= 0;

static void op(long* before, long* after, long* mask, long* maskAfter,
	       const long dx, const long dy, const long dz)
{
  int i;
  int j;
  int k;
  int changed;

  opCount++;

  changes= 0;
  for (k=0; k<dz; k++)
    for (j=0; j<dy; j++) 
      for (i=0; i<dx; i++) {
	changed= 0;
	if (PLOC(mask,i,j,k) != 0) {
	  PLOC(after,i,j,k)= PLOC(before,i,j,k);
	  PLOC(maskAfter,i,j,k)= 1;
	}

	/* Closest 6 neighbors (one cell away) */
	else if (PLOC(mask,i-1,j,k) != 0) {
	  PLOC(after,i,j,k)= PLOC(before,i-1,j,k); changed++;
	}
	else if (PLOC(mask,i+1,j,k) != 0) {
	  PLOC(after,i,j,k)= PLOC(before,i+1,j,k); changed++;
	}

	else if (PLOC(mask,i,j-1,k) != 0) {
	  PLOC(after,i,j,k)= PLOC(before,i,j-1,k); changed++;
	}
	else if (PLOC(mask,i,j+1,k) != 0) {
	  PLOC(after,i,j,k)= PLOC(before,i,j+1,k); changed++;
	}

	else if (PLOC(mask,i,j,k-1) != 0) {
	  PLOC(after,i,j,k)= PLOC(before,i,j,k-1); changed++;
	}
	else if (PLOC(mask,i,j,k+1) != 0) {
	  PLOC(after,i,j,k)= PLOC(before,i,j,k+1); changed++;
	}

	/* 12 Next-closest neighbors (sqrt(2) cells away) */
	else if (PLOC(mask,i-1,j-1,k) != 0) {
	  PLOC(after,i,j,k)= PLOC(before,i-1,j-1,k); changed++;
	}
	else if (PLOC(mask,i-1,j+1,k) != 0) {
	  PLOC(after,i,j,k)= PLOC(before,i-1,j+1,k); changed++;
	}
	else if (PLOC(mask,i-1,j,k-1) != 0) {
	  PLOC(after,i,j,k)= PLOC(before,i-1,j,k-1); changed++;
	}
	else if (PLOC(mask,i-1,j,k+1) != 0) {
	  PLOC(after,i,j,k)= PLOC(before,i-1,j,k+1); changed++;
	}
	else if (PLOC(mask,i+1,j-1,k) != 0) {
	  PLOC(after,i,j,k)= PLOC(before,i+1,j-1,k); changed++;
	}
	else if (PLOC(mask,i+1,j+1,k) != 0) {
	  PLOC(after,i,j,k)= PLOC(before,i+1,j+1,k); changed++;
	}
	else if (PLOC(mask,i+1,j,k-1) != 0) {
	  PLOC(after,i,j,k)= PLOC(before,i+1,j,k-1); changed++;
	}
	else if (PLOC(mask,i+1,j,k+1) != 0) {
	  PLOC(after,i,j,k)= PLOC(before,i+1,j,k+1); changed++;
	}
	else if (PLOC(mask,i,j-1,k-1) != 0) {
	  PLOC(after,i,j,k)= PLOC(before,i,j-1,k-1); changed++;
	}
	else if (PLOC(mask,i,j-1,k+1) != 0) {
	  PLOC(after,i,j,k)= PLOC(before,i,j-1,k+1); changed++;
	}
	else if (PLOC(mask,i,j+1,k-1) != 0) {
	  PLOC(after,i,j,k)= PLOC(before,i,j+1,k-1); changed++;
	}
	else if (PLOC(mask,i,j+1,k+1) != 0) {
	  PLOC(after,i,j,k)= PLOC(before,i,j+1,k+1); changed++;
	}

	/* 8 farthest neighbors (sqrt(3) cells away) */
	else if (PLOC(mask,i-1,j-1,k-1) != 0) {
	  PLOC(after,i,j,k)= PLOC(before,i-1,j-1,k-1); changed++;
	}
	else if (PLOC(mask,i-1,j-1,k+1) != 0) {
	  PLOC(after,i,j,k)= PLOC(before,i-1,j-1,k+1); changed++;
	}
	else if (PLOC(mask,i-1,j+1,k-1) != 0) {
	  PLOC(after,i,j,k)= PLOC(before,i-1,j+1,k-1); changed++;
	}
	else if (PLOC(mask,i-1,j+1,k+1) != 0) {
	  PLOC(after,i,j,k)= PLOC(before,i-1,j+1,k+1); changed++;
	}
	else if (PLOC(mask,i+1,j-1,k-1) != 0) {
	  PLOC(after,i,j,k)= PLOC(before,i+1,j-1,k-1); changed++;
	}
	else if (PLOC(mask,i+1,j+1,k-1) != 0) {
	  PLOC(after,i,j,k)= PLOC(before,i+1,j+1,k-1); changed++;
	}
	else if (PLOC(mask,i+1,j-1,k+1) != 0) {
	  PLOC(after,i,j,k)= PLOC(before,i+1,j-1,k+1); changed++;
	}
	else if (PLOC(mask,i+1,j+1,k+1) != 0) {
	  PLOC(after,i,j,k)= PLOC(before,i+1,j+1,k+1); changed++;
	}

	if (changed) {
	  PLOC(maskAfter,i,j,k)= 1;
	  changes++;
	}
	else {
	  PLOC(after,i,j,k)= PLOC(before,i,j,k);
	}
      }
}

static int doneTest(long* vol)
{
  if (changes==0) return 1;
  else return 0;
}

static void checkFormat( MRI_Dataset* ds, char* fname )
{
  if( !mri_has( ds, "images" ) )
    Abort( "%s operates only on standard images; %s won't work", 
	   progname, fname );
  if( mri_has( ds, "images.dimensions" ) )
    {
      if (strcmp( mri_get_string(ds,"images.dimensions"),"xyz"))
	Abort("%s: input dataset %s must have dimensions xyz!\n",
	      progname,fname);
    }
  else
    Abort( "%s: %s does not have the images.dimensions key.", 
	   progname,fname);
  if (!mri_has(ds,"images.extent.x"))
    Abort("%s: %s is missing images.extent.x tag!\n",progname,fname);
  if (!mri_has(ds,"images.extent.y"))
    Abort("%s: %s is missing images.extent.y tag!\n",progname,fname);
  if (!mri_has(ds,"images.extent.z"))
    Abort("%s: %s is missing images.extent.z tag!\n",progname,fname);
}

int main( argc, argv ) 
     int argc;
     char **argv;
{

  MRI_Dataset *Input = NULL, *Output = NULL, *Mask= NULL;
  char infile[512], hdrfile[512], maskfile[512];
  long dx, dy, dz;
  long *vol0= NULL, *vol1= NULL, *maskVol0= NULL, *maskVol1= NULL; /* padded */
  long *ioVol= NULL; /* not padded */
  long i, j, k;
  long* buf= NULL;
  long* before;
  long* after;
  long* runner;
  long* maskBefore;
  long* maskAfter;

  progname= argv[0];

  /* Print version number */
  Message( "# %s\n", rcsid );

  /* Check to see if help was requested */
  if( ( argc > 1 ) && !strcmp( argv[1], "-help" ) )
    {
      if( argc == 2 )
	Help( "selecttopic" );
      else
	Help( argv[2] );
    }

  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Get filenames */
  cl_get( "output|out", "%option %s[%]", "follow.mri", hdrfile );
  cl_get( "input|i", "%option %s[%]", "input.mri", infile );
  cl_get( "mask|m", "%option %s[%]", "mask.mri", maskfile );
  debug_flag= cl_present("debug");

  if (cl_cleanup_check()) {
    int i;
    fprintf(stderr,"%s: invalid argument in command line:\n    ",argv[0]);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(-1);
  }

  /*** End command-line parsing ***/

  /* Open input dataset */
  if ( !strcmp( infile, hdrfile ) || !strcmp( infile, maskfile )
      || !strcmp( hdrfile, maskfile ) )
    Abort( "Input, output, and mask files must be distinct." );
  Input = mri_open_dataset( infile, MRI_READ );
  Mask  = mri_open_dataset( maskfile, MRI_READ );

  /* Check that program will function on datasets */
  checkFormat(Input,infile);
  checkFormat(Mask,maskfile);

  dx= mri_get_int(Input,"images.extent.x");
  dy= mri_get_int(Input,"images.extent.y");
  dz= mri_get_int(Input,"images.extent.z");

  if ( mri_get_int(Mask,"images.extent.x") != dx 
       || mri_get_int(Mask,"images.extent.y") != dy 
       || mri_get_int(Mask,"images.extent.z") != dz )
    Abort("%s: input and mask are not commensurate!\n",progname);

  /* Set output dataset */
  Output = mri_copy_dataset( hdrfile, Input );
  hist_add_cl( Output, argc, argv );
  mri_set_string( Output, "images.datatype", "int32" );
  mri_set_string( Output, "images.file", ".dat" );
  mri_set_string( Output, "images.dimensions", "xyz" );

  /* Allocate image and parameter storage */
  if (!(vol0= (long*)malloc((dx+2)*(dy+2)*(dz+2)*sizeof(long)))) 
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,(dx+2)*(dy+2)*(dz+2)*sizeof(long));
  if (!(vol1= (long*)malloc((dx+2)*(dy+2)*(dz+2)*sizeof(long)))) 
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,(dx+2)*(dy+2)*(dz+2)*sizeof(long));
  if (!(maskVol0= (long*)malloc((dx+2)*(dy+2)*(dz+2)*sizeof(long)))) 
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,(dx+2)*(dy+2)*(dz+2)*sizeof(long));
  if (!(maskVol1= (long*)malloc((dx+2)*(dy+2)*(dz+2)*sizeof(long)))) 
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,(dx+2)*(dy+2)*(dz+2)*sizeof(long));
  if (!(ioVol= (long*)malloc(dx*dy*dz*sizeof(long)))) 
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,dx*dy*dz*sizeof(long));

  /* Clear */
  runner= vol0;
  for (k=0; k<dz+2; k++)
    for (j=0; j<dy+2; j++)
      for (i=0; i<dx+2; i++) *runner++ = -1;
  runner= vol1;
  for (k=0; k<dz+2; k++)
    for (j=0; j<dy+2; j++)
      for (i=0; i<dx+2; i++) *runner++ = -1;
  runner= maskVol0;
  for (k=0; k<dz+2; k++)
    for (j=0; j<dy+2; j++)
      for (i=0; i<dx+2; i++) *runner++ = 0;
  runner= maskVol1;
  for (k=0; k<dz+2; k++)
    for (j=0; j<dy+2; j++)
      for (i=0; i<dx+2; i++) *runner++ = 0;
  runner= ioVol;
  for (k=0; k<dz; k++)
    for (j=0; j<dy; j++)
      for (i=0; i<dx; i++) *runner++ = 0;

  /* Load 'em up */
  buf= mri_get_chunk(Input, "images", dx*dy*dz, 0, MRI_LONG);
  for (k=0; k<dz; k++)
    for (j=0; j<dy; j++)
      for (i=0; i<dx; i++) {
	PLOC(vol0,i,j,k)= LOC(buf,i,j,k);
      }
  buf= mri_get_chunk(Mask, "images", dx*dy*dz, 0, MRI_LONG);
  for (k=0; k<dz; k++)
    for (j=0; j<dy; j++)
      for (i=0; i<dx; i++) {
	PLOC(maskVol0,i,j,k)= LOC(buf,i,j,k);
      }
  mri_close_dataset( Mask );
  mri_close_dataset( Input );

  /* Do the thing */
  before= vol0;
  after= vol1;
  maskBefore= maskVol0;
  maskAfter= maskVol1;
  while (1) {
    long* tmp;
    op(before,after,maskBefore,maskAfter,dx,dy,dz);
    fprintf(stderr,"step %d: %d changes\n",opCount,changes);
    if (doneTest(after)) break;
    tmp= after;
    after= before;
    before= tmp;
    tmp= maskAfter;
    maskAfter= maskBefore;
    maskBefore= tmp;
  }
      
  /* Write and close output */
  for (k=0; k<dz; k++)
    for (j=0; j<dy; j++)
      for (i=0; i<dx; i++) {
	LOC(ioVol,i,j,k)= PLOC(after,i,j,k);
      }
  mri_set_chunk( Output, "images", dx*dy*dz, 0, MRI_LONG, ioVol );
  mri_close_dataset( Output );
  
  /* Clean up */
  free(vol0);
  free(vol1);
  free(maskVol0);
  free(maskVol1);
  free(ioVol);

  Message( "#      Morphology operations complete (%d ops).\n",
	   opCount);
  return 0;
}

