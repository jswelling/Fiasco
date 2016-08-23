/************************************************************
 *                                                          *
 *  mri_permute.c                                               *
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
 *  Original programming by Chris Genovese 11-95            *
 *     6-96: Pittsburgh Format, Mark Fitzgerald             *
 ************************************************************/
/*************************************************************

  DESCRIPTION OF MRI_PERMUTE.C

  mri_permute.c is used to change the storage order of the
    dimensions of (a chunk of) a dataset 

**************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"
#include "tumble.h"

static char rcsid[] = "$Id: mri_permute.c,v 1.6 2005/01/28 18:06:53 welling Exp $";

#define DEFAULT_MEMLIMIT 52428800

int main( int argc, char* argv[] ) 
{

  MRI_Dataset *Input = NULL, *Output = NULL;
  char infile[512], outfile[512];
  char chunkname[512], outorder[512], *inorder = NULL, buf[512];
  long memlimit;
  long num_dim, *idims = NULL, *odims = NULL;
  long *cntr = NULL, *iorder = NULL, *oorder = NULL;
  long *ioffset = NULL, *ooffset = NULL;
  long k, j, match, done;
  long typesize, maxcp, numcp, extramem, maxstor;
  long long size;
  long ip, op, fip, fop;
  char *indata = NULL, *outdata = NULL;
  long default_memlimit= DEFAULT_MEMLIMIT;
  char* here;
  int verboseFlag= 0;

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /* Allow the user to use more memory */
  if ((here=getenv("F_MEMSIZE_HINT")) != NULL) {
    default_memlimit= atol(here);
    if (default_memlimit==0)
      Abort("%s: environment variable F_MEMSIZE_HINT is not a long integer!\n",
	    argv[0]);
  }

  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Deprecate old options */

  if (cl_present( "dataout|d" ))
     Abort ("Option dataout|d has been replaced by infile outfile format.  Please see help file.\n");

  if (cl_present( "headerout|h" ))
     Abort ("Option headerout|h has been replaced by infile outfile format.  Please see help file.\n");

  if (cl_present( "input|i" ))
     Abort ("Option input|i has been replaced by infile outfile format.  Please see help file.\n");
 
  if (cl_present( "m" ))
     Abort ("Option m has been replaced by memlimit|mem.  Please see help file.\n");

  verboseFlag= cl_present("verbose|ver|v");

  /* Get filenames */
  cl_get( "memlimit|mem", "%option %ld[%]", 
	  default_memlimit, &memlimit );
  cl_get( "chunk|chu|c", "%option %s[%]", "images", chunkname );
  cl_get( "order|ord|o", "%option %s[%]", "vtxyz", outorder );

  if(!cl_get("", "%s", infile)) {
    fprintf(stderr, "%s: Input file name not given.\n", argv[0]);
    exit(-1);
  }

  if(!cl_get("", "%s", outfile)) {
    fprintf(stderr, "%s: Output file name not given.\n", argv[0]);
    exit(-1);
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

  /* Print version number */
  if (verboseFlag) Message( "# %s\n", rcsid );

  /* Open input dataset */
  if( !strcmp( infile, outfile ) )
    Abort( "Input and output files must be distinct." );
  Input = mri_open_dataset( infile, MRI_READ );

  /* Check that program will function on data-set */
  if( !mri_has( Input, chunkname ) )
    Abort( "Chunk %s not found in input dataset.", chunkname );

  /* Get input dimension order and make sure that */
  /*   output order has same number of dimensions */
  snprintf( buf, sizeof(buf), "%s.dimensions", chunkname );
  inorder = mri_get_string( Input, buf );
  num_dim = strlen( inorder );
  if( num_dim != strlen( outorder ) )
    Abort( "Number of dimensions in new ordering does not match input dataset. %s -> %s.\n",
	   inorder, outorder );
  if (!strcmp(inorder,outorder))
    Abort( "The new order and the current order are identical: %s\n",
	   inorder );

  /* Set output dataset */
  Output = mri_copy_dataset( outfile, Input );
  hist_add_cl( Output, argc, argv );
  snprintf( buf, sizeof(buf), "%s.file", chunkname );
  mri_set_string( Output, buf, ".dat" );
  snprintf( buf, sizeof(buf), "%s.dimensions", chunkname );
  mri_set_string( Output, buf, outorder );

  /* Allocate space for various dimensions trackers */
  idims = (long*) emalloc( num_dim * sizeof(long) );
  odims = (long*) emalloc( num_dim * sizeof(long) );
  cntr = (long*) emalloc( num_dim * sizeof(long) );
  iorder = (long*) emalloc( num_dim * sizeof(long) );
  oorder = (long*) emalloc( num_dim * sizeof(long) );
  ioffset = (long*) emalloc( num_dim * sizeof(long) );
  ooffset = (long*) emalloc( num_dim * sizeof(long) );

  /* Check to see if all dimension labels in    */
  /*   output order have a match in input order */
  for( k = 0; k < num_dim; k++ )
    {
      match = 0;
      for( j = 0; j < num_dim; j++ )
	if( outorder[j] == inorder[k] )
	  {
	    iorder[j] = k;
	    oorder[k] = j;
	    match = 1;
	    break;
	  }
      if( !match )
	Abort( "Dimensions of new ordering do not match input dataset.  %s -> %s.\n",
	       inorder, outorder );
    }

  /* Use the simplified "tumble" routine if applicable */
  typesize = get_typesize( Input, chunkname );
  if (tumble_test(inorder, outorder) && memlimit>=2*typesize) {
    int retcode;
    retcode= tumble( Input, Output, inorder, outorder, chunkname, memlimit );
    if (retcode != 0) Abort("Tumble operation failed!\n");
  }
  else {
    /* Get extent of each dimension */
    for( k = 0; k < num_dim; k++ )
      {
	snprintf( buf, sizeof(buf), "%s.extent.%c", chunkname, inorder[k] );
	if( !mri_has( Input, buf ) )
	  Abort( "%s key missing from header.\n", buf );
	idims[k] = odims[oorder[k]] = mri_get_int( Input, buf );
      }
    
    /* Calculate total size of dataset */
    size = typesize;
    for( k = 0; k < num_dim; k++ )
      size *= idims[k];
    
    /* Calculate offsets for dimensions */
    ioffset[0] = ooffset[0] = typesize;
    for( k = 1; k < num_dim; k++ )
      ioffset[k] = ioffset[k-1] * idims[k-1];
    for( k = 1; k < num_dim; k++ )
      ooffset[k] = ooffset[k-1] * odims[k-1];
    
    /* Calculate maximum copy size */
    maxcp = 0;
    while( ( maxcp < num_dim ) && ( oorder[maxcp] == maxcp ) )
      maxcp++;
    
    
    /* PERMUTE */
    
    /* If there is enough room to store two copies of the entire dataset   */
    /*   read it all in and store the permuted copy, then write out        */
    if( memlimit > ( size * 2 ) )
      {
	/* Allocate space for two copies */
	outdata = (char*) emalloc( size );
	
	/* Read in all data */
	if (indata) mri_discard_buffer( Input, indata ); /* free old memory */
	indata = (char*) mri_get_chunk( Input, chunkname, size, 0, MRI_RAW );
	
	/* Begin looping through dimensions to re-order data */
	done = 0;
	numcp = num_dim - maxcp;
	for( k = 0; k < numcp; k++ )
	  cntr[k] = 0;
	while( !done )
	  {
	    
	    /* Calculate position in input dataset */
	    ip = 0;
	    for( k = 0; k < numcp; k++ )
	      ip += cntr[k] * ioffset[iorder[k+maxcp]];
	    
	    /* Calculate position in output dataset */
	    op = 0;
	    for( k = 0; k < numcp; k++ )
	      op += cntr[k] * ooffset[k+maxcp];
	    
	    /* Copy input data to output */
	    memcpy( (void*) ( outdata + op ), (void*) ( indata + ip ), 
		    ooffset[maxcp] );
	    
	    /* Update the loop */
	    if( cntr[0] < ( odims[maxcp] - 1 ) )
	      cntr[0]++;
	    else
	      {
		k = 0;
		while( ( k < numcp ) && ( cntr[k] == ( odims[k+maxcp] - 1 ) ) )
		  cntr[k++] = 0;
		if( k == numcp )
		  done = 1;
		else
		  cntr[k]++;
	      }
	  }
	
	/* Write out data */
	mri_set_chunk( Output, chunkname, size, 0, MRI_RAW, outdata );
      }
    
    /* If there is enough room to store one copy of the entire     */
    /*   dataset, read all data in and then write it out piecemeal */
    else if( memlimit >= ( size + ooffset[0] ) )
      {
	/* Calculate maximum write and copy size */
	extramem = memlimit - size;
	maxstor = 0;
	while( ( maxstor < ( num_dim - 1 ) ) && 
	       ( ooffset[maxstor+1] <= extramem ) )
	  maxstor++;
	maxcp = ( maxcp > maxstor )? maxstor: maxcp;
	
	/* Allocate space for the max number of dimension's worth for output */
	outdata = (char*) emalloc( ooffset[maxstor] );
	
	/* Read in all data */
	if (indata) mri_discard_buffer( Input, indata ); /* free old memory */
	indata = (char*) mri_get_chunk( Input, chunkname, size, 0, MRI_RAW );
	
	/* Begin looping through dimensions to re-order data */
	done = op = fop = 0;
	numcp = num_dim - maxcp;
	for( k = 0; k < numcp; k++ )
	  cntr[k] = 0;
	while( !done )
	  {
	    
	    /* Calculate position in input dataset */
	    ip = 0;
	    for( k = 0; k < numcp; k++ )
	      ip += cntr[k] * ioffset[iorder[k+maxcp]];
	    
	    /* Copy input data to output */
	    memcpy( ( outdata + op ), ( indata + ip ), ooffset[maxcp] );
	    op += ooffset[maxcp];
	    
	    /* If outdata has been filled, then write it out */
	    if( op == ooffset[maxstor] )
	      {
		op = 0;
		mri_set_chunk( Output, chunkname, ooffset[maxstor],
			       fop, MRI_RAW, outdata );
		fop += ooffset[maxstor];
	      }
	    
	    /* Update the loop */
	    if( cntr[0] < ( odims[maxcp] - 1 ) )
	      cntr[0]++;
	    else
	      {
		k = 0;
		while( ( k < numcp ) && ( cntr[k] == ( odims[k+maxcp] - 1 ) ) )
		  cntr[k++] = 0;
		if( k == numcp )
		  done = 1;
		else
		  cntr[k]++;
	      }
	  }
	
      }
    
    /* If there is not enough room to store one copy of the */
    /*   entire dataset, read and write data out piecemeal  */
    else if( memlimit > ( ioffset[0] + ooffset[0] ) )
      {
	/* Calculate maximum read size */
	maxstor = 0;
	while( ( maxstor < ( num_dim - 1 ) ) && 
	       ( ioffset[maxstor+1] <= memlimit ) )
	  maxstor++;
	
	/* Begin looping through dimensions to re-order data */
	
	done = fip = fop = 0;
	ip = ioffset[maxstor];
	numcp = num_dim - maxcp;
	for( k = 0; k < numcp; k++ )
	  cntr[k] = 0;
	while( !done )
	  {
	    
	    /* If all of indata has been written out, read in more data */
	    if( ip == ioffset[maxstor] )
	      {
		ip = 0;
		if (indata) 
		  mri_discard_buffer( Input, indata ); /* free old memory */
		indata = (char*) 
		  mri_get_chunk( Input, chunkname, ioffset[maxstor], 
				 fip, MRI_RAW );
		fip += ioffset[maxstor];
	      }
	    
	    /* Calculate position in output */
	    fop = 0;
	    for( k = 0; k < numcp; k++ )
	      fop += cntr[k] * ooffset[oorder[k+maxcp]];
	    
	    /* Write out piece of output */
	    mri_set_chunk( Output, chunkname, ooffset[maxcp],
			   fop, MRI_RAW, ( indata + ip ) );
	    ip += ooffset[maxcp];
	    
	    /* Update the loop */
	    if( cntr[0] < ( idims[maxcp] - 1 ) )
	      cntr[0]++;
	    else
	      {
		k = 0;
		while( ( k < numcp ) && ( cntr[k] == ( idims[k+maxcp] - 1 ) ) )
		  cntr[k++] = 0;
		if( k == numcp )
		  done = 1;
		else
		  cntr[k]++;
	      }
	  }
	
      }
    
    /* Not enough available memory to do anything */
    else
      Abort( "Memory limit too small to work: (%ld)\n.", memlimit );
  }
    
  /* Write and close data-sets */
  mri_close_dataset( Input );
  mri_close_dataset( Output );

  if (verboseFlag) Message( "#      Permutation complete.\n" );

  return 0;
}
      
