/************************************************************
 *                                                          *
 *  image_entropy.c                                            *
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
 *  Original programming by Chris Hefferan 6/04             *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "misc.h"
#include "stdcrg.h"

static char rcsid[] = "$Id: image_entropy.c,v 1.11 2005/09/27 20:22:36 welling Exp $";

#define KEYBUF_SIZE 512

static char* progname;
static int debug= 0;
static int verbose_flag= 0;

static void safe_copy(char* str1, char* str2) {
  strncpy(str1, str2, KEYBUF_SIZE);
  str1[KEYBUF_SIZE-1]= '\0';
}

static void safe_concat(char* str1, char* str2) {
  strncat(str1, str2, (KEYBUF_SIZE-strlen(str1))-1);
}

static int safe_get_extent(MRI_Dataset* ds, char* chunk, char* dim)
{
  char key_buf[KEYBUF_SIZE];
  char dim_buf[4];
  dim_buf[0]= *dim;
  dim_buf[1]= '\0';
  safe_copy(key_buf,chunk);
  safe_concat(key_buf,".extent.");
  safe_concat(key_buf,dim_buf);
  if (mri_has(ds,key_buf)) return mri_get_int(ds,key_buf);
  else Abort("%s: input missing tag %s!\n",progname,key_buf);
  return 0; /* not reached */
}

static int input_valid(MRI_Dataset* imgDS, long* dx, long* dy, long* dz, 
		       long* dt)
{
  char key_buf[KEYBUF_SIZE];
  char* dimstr;
  char* here;

  if (!mri_has(imgDS,"images.dimensions")) {
    Error("%s: weight file has no images.dimensions tag!\n",progname);
    return 0;
  }
  dimstr= mri_get_string(imgDS,"images.dimensions");
  if (!strncmp(dimstr,"xyz",3)) {
    /* This will definitely work. */
  }
  else if (!strncmp(dimstr,"vxyz",4)) {
    if (safe_get_extent(imgDS,"images","v") != 1) {
      Abort("%s: weight dataset must have v extent 1!\n",progname);
      return 0;
    }
  }
  else {
    Error("%s: weight dataset must have dimensions (v)xyz(...)!\n",progname);
    return 0;
  }

  *dx= safe_get_extent(imgDS,"images","x");
  *dy= safe_get_extent(imgDS,"images","y");
  *dz= safe_get_extent(imgDS,"images","z");
  *dt= 1;
  here= strchr(dimstr,'z')+1;
  while (*here) {
    *dt *= safe_get_extent(imgDS,"images",here);
    here++;
  }
  if (debug)
    fprintf(stderr,"Got dims %ld, %ld, %ld, %ld\n",*dx,*dy,*dz,*dt);
  return 1;
}

static FILE* initParFile(char* parfname, const int filtered_flag,
			 const long dx, const long dy, const long dz)
{
  FILE* result;
  time_t tm;
  if (!(result= fopen(parfname,"w"))) {
    Abort("%s: unable to open <%s> for writing!\n",progname,parfname);
  }

  tm= time(NULL);
  fprintf(result,"##Format: order:index_t, type:%s\n",
	  filtered_flag ? "filtered":"raw");
  fprintf(result,"##Format: names:(entropy)\n");
  fprintf(result,"# Generated at %s",asctime(localtime(&tm)));
  fprintf(result,"# dims %ld, %ld, %ld\n",dx, dy, dz);
  fflush(result);
  return result;
}

static void writeOutput( int t, double val, FILE* ofp ) {
  /* Write 'em out */
  fprintf(ofp,"%d %17.10g\n",t,val);
}

int main( int argc, char** argv ) 
{
  char infile[512], parfile[512], maskfile[512];
  FILE *ofp = NULL;
  MRI_Dataset* imgDS= NULL;
  MRI_Dataset* maskDS= NULL;
  long dx=0, dy=0, dz=0, dt=0;
  long t;
  int filtered_flag = 0;
  int mask_present= 0;
  long mask_dt;
  EntropyContext* ec= NULL;
  long num_of_boxes= 0;
  int min_set= 0;
  int max_set= 0;
  double min= 0.0;
  double max= 0.0;

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
  cl_get( "estimates|est|e", "%option %s[%]", "image_entropy.par", parfile );
  debug= cl_present("debug|deb");
  verbose_flag= cl_present("v|verbose");
  cl_get( "nbins","%option %d[0]",&num_of_boxes);
  max_set= cl_get("max","%option %lf",&max);
  min_set= cl_get("min","%option %lf",&min);
  mask_present= cl_get("mask","%option %s", maskfile);

  if(!cl_get("", "%s", infile)) {
    fprintf(stderr, "%s: Input file name not given.\n", argv[0]);
    Help( "usage" );
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

  if (max_set && min_set) {
    if (max<=min) Abort("%s: inconsistent max and min!\n",argv[0]);
  }

  imgDS= mri_open_dataset(infile,MRI_READ);
  if (!input_valid(imgDS, &dx, &dy, &dz, &dt))
    Abort("%s: input file format is inappropriate.\n",argv[0]);
  
  if (mask_present) {
    long mask_dx, mask_dy, mask_dz;
    if (debug) fprintf(stderr,"reading <%s> as a mask\n",maskfile);
    maskDS= mri_open_dataset(maskfile,MRI_READ);
    if (!input_valid(maskDS, &mask_dx, &mask_dy, &mask_dz, &mask_dt))
      Abort("%s: mask file format is inappropriate.\n",argv[0]);
    if (mask_dx != dx || mask_dy != dy || mask_dz != dz)
      Abort("%s: mask dimensions do not match input image dimensions!\n",
	    argv[0]);
  }

  ec= ent_createContext();
  ent_setDebug(ec,debug);
  ent_setVerbose(ec,verbose_flag);
  if (max_set) ent_setMax(ec,max);
  if (min_set) ent_setMin(ec,min);
  if (num_of_boxes != 0) ent_setNBins(ec,num_of_boxes);

  /* Open files */
  if (debug) fprintf(stderr,"writing <%s>\n",parfile);
  ofp= initParFile(parfile, filtered_flag, dx, dy, dz);

  /* Read in all input lines, generating appropriate output */
  if (mask_present){
    for (t=0; t<dt; t++) {
      double* imgBuf= mri_get_chunk(imgDS,"images",dx*dy*dz,
				    t*dx*dy*dz, MRI_DOUBLE);
      int* maskBuf= mri_get_chunk(maskDS,"images",dx*dy*dz,
				  (t%mask_dt)*dx*dy*dz, MRI_INT);
      writeOutput(t, ent_calcMaskedImageEntropyDouble(ec,imgBuf,maskBuf,
						      dx,dy,dz,1,1), 
		  ofp);
    }
  }
  else {
    for (t=0; t<dt; t++) {
      double* imgBuf= mri_get_chunk(imgDS,"images",dx*dy*dz,
				    t*dx*dy*dz, MRI_DOUBLE);
      writeOutput(t, ent_calcImageEntropyDouble(ec,imgBuf,dx,dy,dz,1), ofp);
    }
  }
  
  /* Close files */
  if (fclose(ofp)) {
    perror("Error closing output:");
  }
  mri_close_dataset(imgDS);
  if (mask_present) mri_close_dataset(maskDS);

  Message( "#      Entropy estimates calculated (%d records).\n",dt);

  if (ec) ent_destroyContext(ec);
  return 0;
}

