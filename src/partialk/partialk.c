/************************************************************
 *                                                          *
 *  partialk.c                                              *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1995 Department of Statistics,         *
 *                        Carnegie Mellon University, and   *
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
 *                        Pittsburgh Supercomputing Center  *
 *                                                          *
 *  Original programming by Mark Fitzgerald  2-95           *
 *     5-96: Pittsburgh Format, Mark Fitzgerald             *
 *     9-96: Parallelization, Greg Hood, PSC
 *     12-97: Creation of partialk.c from baseline.c,
 *            Joel Welling, PSC/CMU
 ************************************************************/
/*************************************************************

  DESCRIPTION OF PARTIALK.C

  partialk is used to extend the negative phase direction
    of k-space data which has been collected with all the
    positive phase information but only partial negative
    information.  The result is symmetrized data, with
    equal coverage in the positive and negative phase
    directions.

  partialk   [-reverse even|odd|all|none] -dim nnn 
       -input Input-mri-file -output Output-mri-file

**************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"
#include "par.h"
#include "array.h"
#include "misc.h"
#include "acct.h"

typedef struct Context {
  /* NOTE: any new fields added to this struct should
     be also added to PackContext and UnpackContext */ 
  Filename input_file;		/* the images to work on */
  Filename output_file;		/* the output dataset */

  int dx;			/* # of pixels along the X dimension */
  int dy;			/* # of pixels along the Y dimension */
  int dy_final;                 /* # of pixels along Y after extension */
  int dz;			/* # of slices along the Z dimension */
  int dt;			/* # of images along the T dimension */
  int flip_odd;                 /* # whether or not to flip odd rows for FFT */
  int flip_even;                /* # whether or not to flip even rows fr FFT */
} Context;

typedef struct Task {
  /* NOTE: any new fields added to this struct should
     be also added to PackTask and UnpackTask */ 
  int t;			/* image number to work on */
  int z;			/* slice number to work on */
} Task;

typedef struct Result {
  /* NOTE: any new fields added to this struct should
     be also added to PackResult and UnpackResult */
  int t;
  int z;
} Result;

static char rcsid[] = "$Id: partialk.c,v 1.14 2007/04/19 22:57:52 welling Exp $";

/* GLOBAL VARIABLES FOR MASTER & WORKER */
static char* progname;
Task t;
Context c;
Result r;

/* GLOBAL VARIABLES FOR MASTER */
static MRI_Dataset *mInput;
static MRI_Dataset *mOutput;
static char reverse_string[256];
static int reverse_set= 0;
static int dim_set= 0;

/* GLOBAL VARIABLES FOR WORKER */
static MRI_Dataset *wInput;		/* input dataset */
static MRI_Dataset *wOutput;		/* output dataset */
static FComplex **image = NULL;		/* input image */
static FComplex **corr_image = NULL;	/* corrected image */

/* FORWARD DECLARATIONS */
void MasterTask (int argc, char **argv);
void ReadArguments (int argc, char **argv);
void MasterResult (int task_number);
void WorkerContext ();
void WorkerTask ();
void WorkerFinalize();
void PackContext ();
void UnpackContext ();
void PackTask ();
void UnpackTask ();
void PackResult ();
void UnpackResult ();
static int parse_cl_reverse( const char* s, int* revOdd, int* revEven );
static void parse_mrihdr_reverse( MRI_Dataset* ds, const char* chunk, 
				  int* revOdd, int* revEven );


void flip_rows(FComplex** image);

int
main (int argc,
      char **argv,
      char **envp) 
{
  progname= argv[0];
  par_process(argc, argv, envp,
	      MasterTask, MasterResult,
	      WorkerContext, WorkerTask,
	      WorkerFinalize,
	      PackContext, UnpackContext,
	      PackTask, UnpackTask,
	      PackResult, UnpackResult);
  return(0);
}

/* MASTER PROCEDURES */

static int parse_cl_reverse( const char* s, int* revOdd, int* revEven )
{
  if (!strcmp(s,"none")) {
    *revOdd= 0;
    *revEven= 0;
  }
  else if (!strcmp(s,"odd")) {
    *revOdd= 1;
    *revEven= 0;
  }
  else if (!strcmp(s,"even")) {
    *revOdd= 0;
    *revEven= 1;
  }
  else if (!strcmp(s,"all")) {
    *revOdd= 1;
    *revEven= 1;
  }
  else {
    return 0;
  }
  return 1;
}

static void parse_mrihdr_reverse( MRI_Dataset* ds, const char* chunk, 
				  int* revOdd, int* revEven )
{
  char buf[256];

  buf[sizeof(buf)-1]= '\0';
  snprintf(buf,sizeof(buf)-1,"%s.rowflip",chunk);
  if (mri_has(ds,buf) 
      && mri_get_int(ds,buf) != 0) {
    snprintf(buf,sizeof(buf)-1,"%s.rowflip_pattern",chunk);
    if (mri_has(ds,buf)) {
      char* pattern= mri_get_string(ds,buf);
      if (!strcmp(pattern,"none")) {
	*revOdd= 0;
	*revEven= 0;
      }
      else if (!strcmp(pattern,"odd")) {
	*revOdd= 1;
	*revEven= 0;
      }
      else if (!strcmp(pattern,"even")) {
	*revOdd= 0;
	*revEven= 1;
      }
      else if (!strcmp(pattern,"both")) {
	*revOdd= 1;
	*revEven= 1;
      }
      else Abort("%s: input dataset has invalid %s tag %s!\n",
		 progname,buf,pattern);
    }
    else Abort("%s: input dataset is missing %s info!\n",
	       progname,buf);
  }
  else {
    *revOdd= 0;
    *revEven= 0;
  }
}

void
MasterTask (int argc,
	    char **argv)
{
  int dv;
  char* dimstr;

  verbose = TRUE;
  /* Print version number */
  Report("# %s\n",rcsid);

  /* Check to see if help was requested */
  if( ( argc > 1 ) && !strcmp( argv[1], "-help" ) )
    {
      if( argc == 2 )
	Help( "selecttopic" );
      else
	Help( argv[2] );
    }

  ReadArguments(argc, argv);

  /* Open input dataset */
  if( !strcmp( c.input_file, c.output_file ) )
    Abort( "Input and output files must be distinct." );
  mInput = mri_open_dataset( c.input_file, MRI_READ );

  /* Check that program will function on data-set */
  if( !mri_has( mInput, "images" ) )
    Abort( "%s operates only on standard images.\n", argv[0] );
  if ( !mri_has( mInput, "images.dimensions" ) )
    Abort( "%s has no images.dimensions tag!\n", argv[0] );
  dimstr= mri_get_string( mInput, "images.dimensions" );
  if (dimstr[0] != 'v' || !mri_has( mInput, "images.extent.v" )
      || (mri_get_int( mInput, "images.extent.v" ) != 2 ) )
    Abort( "%s operates only on complex-valued images (v = 2).\n", argv[0] );
  if ( strcmp(dimstr,"vxyzt") && strcmp(dimstr,"vqyzt")
       && strcmp(dimstr,"vxyz") && strcmp(dimstr,"vqyz")
       && strcmp(dimstr,"vxy") && strcmp(dimstr,"vqy") )
    Abort("%s operates only on vxy(z)(t) or vqy(z)(t) images!\n",argv[0]);
  if (strchr(dimstr,'x')) {
    if (mri_has(mInput,"images.extent.x")) 
      c.dx= mri_get_int(mInput,"images.extent.x");
    else Abort("%s: input is missing images.extent.x tag!\n",argv[0]);
  }
  else {
    if (mri_has(mInput,"images.extent.q")) 
      c.dx= mri_get_int(mInput,"images.extent.q");
    else Abort("%s: input is missing images.extent.q tag!\n",argv[0]);
  }
  if (mri_has(mInput,"images.extent.y")) 
    c.dy= mri_get_int(mInput,"images.extent.y");
  else Abort("%s: input is missing images.extent.y tag!\n",argv[0]);
  if (mri_has(mInput,"images.extent.z")) 
    c.dz= mri_get_int(mInput,"images.extent.z");
  else c.dz= 1;
  if (mri_has(mInput,"images.extent.t")) 
    c.dt= mri_get_int(mInput,"images.extent.t");
  else c.dt= 1;

  /* If scan line reversal information is given in the command line,
   * implement it.  Otherwise, infer from the input file.
   */
  if (reverse_set) {
    if (!parse_cl_reverse( reverse_string, &(c.flip_odd), &(c.flip_even) ))
      Abort("%s: invalid -reverse opion <%s> given\n", 
	    progname,reverse_string);
  }
  else parse_mrihdr_reverse(mInput, "images", &(c.flip_odd), &(c.flip_even));

  /* If partialk completion information is given in the command line,
   * implement it.  Otherwise, infer from the input file.
   */
  if (dim_set) {
    /* correct value read directly into context */
  }
  else {
    if (mri_has(mInput,"images.dy_base"))
      c.dy_final= mri_get_int(mInput,"images.dy_base");
    else Abort("%s: information on final Y resolution not provided!\n",
	       progname);
  }

  /* Set output dataset */
  mOutput = mri_copy_dataset( c.output_file, mInput );
  hist_add_cl( mOutput, argc, argv );
  mri_set_string( mOutput, "images.datatype", "float32" );
  mri_set_int( mOutput, "images.extent.y", c.dy_final );
  mri_set_int( mOutput, "images.rowflip", 0 );
  mri_set_int( mOutput, "images.partialk", 0 );
  mri_close_dataset ( mOutput);

  mri_close_dataset ( mInput );

  /* Reality check input dimensions */
  if( ( c.dt <= 0 ) || ( c.dx <= 0 ) || ( c.dy <= 0 ) || ( c.dz <= 0 ) )
    Abort( "%s: images.extent key(s) is non-positive.\n", argv[0] );
  if ( c.dy < (c.dy_final/2) )
    Abort( "%s: insufficient k space coverage in input file\n", argv[0] );
  if ( c.dy > (c.dy_final) )
    Abort( "%s: output y range smaller than input!\n", argv[0] );

  /* set the context for the workers */
  par_set_context();

  /* Loop through images */
  for ( t.t = 0; t.t < c.dt; t.t++ )
    for ( t.z = 0; t.z < c.dz; t.z++ ) 
      par_delegate_task();
  par_finish();
        
  Report( "\n#      Partialk adjustment complete.\n" );
  exit(0);
}

void
ReadArguments (int argc,
	       char **argv)
{
  cl_scan( argc, argv );

  /* Get filenames */
  cl_get( "output|o", "%option %s[%]", "partialk.mri", c.output_file );
  cl_get( "input|i", "%option %s[%]", "input.mri", c.input_file );
  reverse_set= cl_get( "reverse|r", "%option %s", reverse_string );
  dim_set= cl_get( "dim|d", "%option %d", &(c.dy_final) );

  if (cl_cleanup_check()) {
    int i;
    fprintf(stderr,"%s: invalid argument in command line:\n    ",argv[0]);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(-1);
  }
}

void MasterResult (int task_number)
{
  static int results = 0;
  int n;

  if (++results % c.dz != 0)
    return;

  /* Print progress report */
  n = (results-1) / c.dz;
  if (n == 0)
    Report( "      " );

  if (n == (c.dt - 1) || (n+1) % 60 == 0)
    {
      Report( "# %ld\n", (long) ( n + 1 ) );
      if (n != (c.dt - 1))
	Report( "      " );
      else
	results = 0;
    }
  else
    Report("#");
}

/* WORKER PROCEDURES */

void
WorkerFinalize()
{
  if (wInput != NULL) mri_close_dataset(wInput);
  if (wOutput != NULL) mri_close_dataset(wOutput);
  if (image != NULL)
    free(image);
  if (corr_image != NULL)
    FreeMatrix(corr_image);
}

void
WorkerContext ()
{
  static int first = TRUE;
  static int old_dx = 0, old_dy = 0, old_dy_final= 0;
  int x, y;

  if (first)
    {
      wInput = mri_open_dataset(c.input_file, MRI_READ);
      wOutput = mri_open_dataset(c.output_file, MRI_MODIFY_DATA);
      first = FALSE;
    }

  if (c.dx != old_dx || c.dy != old_dy || c.dy_final != old_dy_final)
    {
      if (image != NULL)
	free(image);
      if (corr_image != NULL)
	FreeMatrix(corr_image);

      /* Allocate image storage */
      image = (FComplex **) emalloc( c.dy * sizeof(FComplex *) );
      corr_image = Matrix( c.dy_final, c.dx, FComplex );

      old_dx = c.dx;
      old_dy = c.dy;
      old_dy_final = c.dy_final;
    }
}

void
WorkerTask ()
{
  long x, y, x_flip, y_flip, y_shift;
  long imageSize= 2*c.dx*c.dy;
  long imageSize_final= 2*c.dx*c.dy_final;

  /* Get image and re-align matrix pointers */
  *image = (FComplex *) 
    mri_get_chunk( wInput, "images", imageSize,
		   (t.t*c.dz + t.z)*imageSize, MRI_FLOAT );
  realign_matrix( (void **) image, c.dy, 
		  (long) ( c.dx * sizeof(FComplex) ) );
	
  if (c.flip_odd || c.flip_even) flip_rows( image );

  /*** PARTIALK CORRECTION ***/

  /* Transcribe the known part of the input image to the output.
   * Do this by counting down from the top of k space.
   */

  y_shift= c.dy_final - c.dy;
  for ( y=0; y<c.dy; y++)
    for (x=0; x<c.dx; x++) {
      corr_image[y+y_shift][x].real= image[y][x].real;
      corr_image[y+y_shift][x].imag= image[y][x].imag;
    }

  if (c.dy < c.dy_final) {
    /* Slap in the missing part */
    for ( y=0; y<y_shift; y++) {
      y_flip= c.dy - (y+1);
      for (x=0; x<c.dx; x++) {
	x_flip= c.dx - (x+1);
	corr_image[y][x].real= image[y_flip][x_flip].real;
	corr_image[y][x].imag= -image[y_flip][x_flip].imag;
      }
    }
  }
	
  /*** END PARTIALK CORRECTION ***/
  
  /* Set output image */
  mri_set_chunk( wOutput, "images", imageSize_final, 
		 (t.t*c.dz + t.z)*imageSize_final, MRI_FLOAT, *corr_image );
}

void flip_rows(FComplex** image)
{
  static FComplex* scratch= NULL;
  int x,y;
  /* This routine reverses rows of the image if requested */

  if (!scratch) {
    /* first pass */
    if (!(scratch= (FComplex*)malloc(c.dx*sizeof(FComplex))))
      Abort("Unable to allocate %d float complex pairs!\n",c.dx);
  }

  for (y=0; y<c.dy; y++) {
    if (y%2) {
      if (c.flip_odd) {
	for (x=0; x<c.dx; x++) {
	  scratch[x].real= image[y][x].real;
	  scratch[x].imag= image[y][x].imag;
	}
	for (x=0; x<c.dx; x++) {
	  image[y][c.dx-(x+1)].real= scratch[x].real;
	  image[y][c.dx-(x+1)].imag= scratch[x].imag;
	}
      }
    }
    else {
      if (c.flip_even) {
	for (x=0; x<c.dx; x++) {
	  scratch[x].real= image[y][x].real;
	  scratch[x].imag= image[y][x].imag;
	}
	for (x=0; x<c.dx; x++) {
	  image[y][c.dx-(x+1)].real= scratch[x].real;
	  image[y][c.dx-(x+1)].imag= scratch[x].imag;
	}
      }
    }
  }
}

/* UTILITY FUNCTIONS */

void
PackContext ()
{
  par_pkstr(c.input_file);
  par_pkstr(c.output_file);
  par_pkint(c.dx);
  par_pkint(c.dy);
  par_pkint(c.dy_final);
  par_pkint(c.dz);
  par_pkint(c.dt);
  par_pkint(c.flip_odd);
  par_pkint(c.flip_even);
}

void
UnpackContext ()
{
  par_upkstr(c.input_file);
  par_upkstr(c.output_file);
  c.dx= par_upkint();
  c.dy= par_upkint();
  c.dy_final= par_upkint();
  c.dz= par_upkint();
  c.dt= par_upkint();
  c.flip_odd= par_upkint();
  c.flip_even= par_upkint();
}

void
PackTask ()
{
  par_pkint(t.t);
  par_pkint(t.z);
}

void
UnpackTask ()
{
  t.t= par_upkint();
  t.z= par_upkint();
}

void
PackResult ()
{
  par_pkint(r.t);
  par_pkint(r.z);
}

void
UnpackResult ()
{
  r.t= par_upkint();
  r.z= par_upkint();
}


