/************************************************************
 *                                                          *
 *  pullback.c                                                *
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

#define BLOCKSIZE (1024*1024)
#define KEYBUF_SIZE 512

static char rcsid[] = "$Id: pullback.c,v 1.2 2005/06/17 23:01:26 welling Exp $";

/* Access for 3D arrays */
#define LOC(matrix,x,y,z) matrix[((((z)*dy)+(y))*dx)+(x)]

static int debug_flag= 0;
static int verbose_flag= 0;
static char* progname;

static const char* safe_get_dims(MRI_Dataset* ds, char* chunk)
{
  char key_buf[KEYBUF_SIZE];
  if (snprintf(key_buf,sizeof(key_buf),"%s.dimensions",chunk)
      >= sizeof(key_buf))
    Abort("%s: chunk name is too long!\n",progname);
  if (mri_has(ds,key_buf)) return mri_get_string(ds,key_buf);
  else Abort("%s: input missing tag %s!\n",progname,key_buf);
  return 0; /* not reached */
}

static long safe_get_extent(MRI_Dataset* ds, char* chunk, char dim)
{
  char key_buf[KEYBUF_SIZE];
  if (snprintf(key_buf,sizeof(key_buf),"%s.extent.%c",chunk,dim)
      >= sizeof(key_buf))
    Abort("%s: chunk name is too long!\n",progname);
  if (mri_has(ds,key_buf)) return (long)mri_get_int(ds,key_buf);
  else Abort("%s: input missing tag %s!\n",progname,key_buf);
  return 0; /* not reached */
}

static long long getChunkSize(MRI_Dataset* ds, char* chunk)
{
  long long result = 1;
  const char* dimstr= safe_get_dims(ds,chunk);
  const char* runner= dimstr;
  while (*runner) {
    result *= safe_get_extent(ds,chunk,*runner);
    runner++;
  }
  if (debug_flag) fprintf(stderr,"Dataset %s chunk %s has size %lld\n",
			  ds->name,chunk,result);
  return result;
}

int main( argc, argv ) 
     int argc;
     char **argv;
{

  MRI_Dataset *Input = NULL, *Addr= NULL, *Output = NULL;
  char infile[512], addrfile[512], outfile[512], chunkName[512];
  char key_buf[KEYBUF_SIZE];
  long dx, dy, dz, dq;
  double* inVoxels;
  double* outVol;
  long *addrVoxels;
  long i;
  long* buf= NULL;
  long long infileLength;
  long long baseOffset;
  long voxThisBlock;

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
  cl_get( "chunk|c", "%option %s[%]", "images", chunkName );
  if (!cl_get( "qdim", "%option %d", &dq ))
    Abort("%s: required argument qdim not given\n",progname);
  debug_flag= cl_present("debug");
  verbose_flag= cl_present("verbose|v");

  if (!cl_get( "addr|a", "%option %s", addrfile ))
    fprintf(stderr,"%s: required address file name not given.\n",progname);

  if (!cl_get("", "%s", infile)) {
    fprintf(stderr,"%s: Input file name not given.\n",progname);
    Help( "usage" );
    exit(-1);
  }

  if (!cl_get("", "%s", outfile)) {
    fprintf(stderr,"%s: Output file name not given.\n",progname);
    Help( "usage" );
    exit(-1);
  }


  if (cl_cleanup_check()) {
    fprintf(stderr,"%s: invalid argument in command line:\n    ",progname);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(-1);
  }

  /*** End command-line parsing ***/

  /* Open input dataset */
  if ( !strcmp( infile, outfile ) || !strcmp( infile, addrfile )
      || !strcmp( outfile, addrfile ) )
    Abort( "Input, output, and addr files must be distinct." );
  Input = mri_open_dataset( infile, MRI_READ );
  Addr  = mri_open_dataset( addrfile, MRI_READ );

  /* Check that program will function on datasets */
  infileLength= getChunkSize(Input,chunkName);
  if (getChunkSize(Addr,chunkName)!= infileLength)
    Abort("%s: input file and address file are not commensurate!\n",
	  progname);

  /* Set output dataset */
  Output = mri_copy_dataset( outfile, Input );
  hist_add_cl( Output, argc, argv );
  if (snprintf(key_buf,sizeof(key_buf),"%s.datatype",chunkName)
      >= sizeof(key_buf))
    Abort("%s: chunk name too long!\n",progname);
  mri_set_string( Output, key_buf, "float64" );
  if (snprintf(key_buf,sizeof(key_buf),"%s.file",chunkName)
      >= sizeof(key_buf))
    Abort("%s: chunk name too long!\n",progname);
  mri_set_string( Output, key_buf, ".dat" );
  if (snprintf(key_buf,sizeof(key_buf),"%s.dimensions",chunkName)
      >= sizeof(key_buf))
    Abort("%s: chunk name too long!\n",progname);
  mri_set_string( Output, key_buf, "q" );
  if (snprintf(key_buf,sizeof(key_buf),"%s.extent.q",chunkName)
      >= sizeof(key_buf))
    Abort("%s: chunk name too long!\n",progname);
  mri_set_int( Output, key_buf, dq );

  /* Allocate image and parameter storage */
  if (!(outVol= (double*)malloc(dq*sizeof(double)))) 
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,dq*sizeof(double));

  if (verbose_flag) Message("Mapping %lld voxels into %ld bins\n",
			    infileLength,dq);

  /* Initialize */
  for (i=0; i<dq; i++) outVol[i]= 0.0;

  /* Transfer the voxels */
  baseOffset= 0;
  while (baseOffset<infileLength) {
    voxThisBlock= 
      ((infileLength-baseOffset < BLOCKSIZE)?
       (infileLength-baseOffset):BLOCKSIZE);
    inVoxels= mri_get_chunk(Input,chunkName,voxThisBlock,
			    baseOffset,MRI_DOUBLE);
    addrVoxels= mri_get_chunk(Addr,chunkName,voxThisBlock,
			      baseOffset,MRI_LONG);
    for (i=0; i<voxThisBlock; i++) {
      if (addrVoxels[i]>=0 && addrVoxels[i]<dq)
	outVol[addrVoxels[i]] += inVoxels[i];
      else {
	Warning(1,"%s: address %ld is out of range!\n",progname,addrVoxels[i]);
      }
    }
    baseOffset += voxThisBlock;
  }

  mri_set_chunk( Output, chunkName, dq, 0, MRI_DOUBLE, outVol );
  
  /* Clean up */
  free(outVol);
  mri_close_dataset( Output );
  mri_close_dataset( Addr );
  mri_close_dataset( Input );

  Message( "#      Pullback complete.\n" );
  return 0;
}

