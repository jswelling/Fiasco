/************************************************************
 *                                                          *
 *  recon.c                                                 *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1995,1997 Department of Statistics,    *
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
 *			  Pittsburgh Supercomputing Center  *
 *                                                          *
 *  Original programming by Mark Fitzgerald  2-95           *
 *     5-96: Pittsburgh Format, Mark Fitzgerald             *
 *     7-97: Parallelization, Greg Hood, PSC                *
 ************************************************************/
/*************************************************************

  DESCRIPTION OF RECON.C

  recon.c applies the 2-D (inverse) Fourier transform
    to a set of images

  recon.m [-input Input-header-file] [-headerout Output-header-file]
          [-dataout Output-data-file] [-direction Which-transform]
	  [-recon Output-form]

**************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <unistd.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"
#include "par.h"
#include "misc.h"
#include "acct.h"

static char rcsid[] = "$Id: recon.c,v 1.14 2007/04/19 23:04:04 welling Exp $";

typedef struct Context {
  /* NOTE: any new fields added to this struct should
     be also added to PackContext and UnpackContext */ 
  Filename input_file;		/* the images to work on */
  Filename output_header_file;	/* the name of the output header file */
  Filename output_data_file;	/* the name of the output data file */

  int dv;			/* # of values along the V dimension */
  int dx;			/* # of pixels along the X dimension */
  int dy;			/* # of pixels along the Y dimension */
  int dz;			/* # of slices along the Z dimension */
  int dt;			/* # of images along the T dimension */

  int fftsign;			/* if  1, regular FFT
				   if -1, inverse FFT */
  int complex;			/* if 1, output complex images */
  int phase;			/* if complex is 0, this field controls
				   what to output:
				       0 = output modulus
				       1 = output phase */
} Context;

typedef struct Task {
  /* NOTE: any new fields added to this struct should
     be also added to PackTask and UnpackTask */ 
  int z;			/* slice number to work on */
  int t;			/* image number to work on */
} Task;

/* GLOBAL VARIABLES FOR MASTER & WORKER */
Task t;
Context c;

/* GLOBAL VARIABLES FOR MASTER */
static MRI_Dataset *mInput;
static MRI_Dataset *mOutput;

/* GLOBAL VARIABLES FOR WORKER */
static MRI_Dataset *wInput;
static MRI_Dataset *wOutput;
static FComplex **image = NULL;
static float **real_image = NULL;
static float **recon_image = NULL;

/* FORWARD DECLARATIONS */
void MasterTask (int argc, char **argv, char **envp);
void ReadArguments (int argc, char **argv);
void MasterResult (int task_number);
void WorkerContext ();
void WorkerTask ();
void WorkerFinalize();
void PackContext ();
void UnpackContext ();
void PackTask ();
void UnpackTask ();

int
main (int argc,
      char **argv,
      char **envp)
{
  par_process(argc, argv, envp,
	      MasterTask, MasterResult,
	      WorkerContext, WorkerTask,
	      WorkerFinalize,
	      PackContext, UnpackContext,
	      PackTask, UnpackTask,
	      NULL, NULL);
  return(0);
}

/* MASTER PROCEDURES */

void
MasterTask (int argc,
	    char **argv,
	    char **envp)
{
  char cfftsign[512];
  char crecon[512];
  int i;
  int temp_fixed;
  int tt, zz;
  Filename short_name;

  verbose = TRUE;
  Report("# %s\n", rcsid);

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
  cl_get( "dataout|d", "%option %s[%]", ".dat", c.output_data_file );
  cl_get( "headerout|h", "%option %s[%]", "recon.mri", c.output_header_file );
  cl_get( "input|inp", "%option %s[%]", "input.mri", c.input_file );
  cl_get( "direction|dir", "%option %s[%]", "inverse", cfftsign );
  if (!strcasecmp( cfftsign, "forward" ) || !strcmp( cfftsign, "1" )
      || !strcasecmp( cfftsign, "fwd" )) 
    c.fftsign= 1;
  else if (!strcasecmp( cfftsign, "backward" ) || !strcmp( cfftsign, "-1" )
	   || !strcasecmp( cfftsign, "bkwd" ) 
	   || !strcasecmp( cfftsign, "inverse")
	   || !strcasecmp( cfftsign, "inv"))
    c.fftsign = -1;
  else Abort("%s: unrecognized direction string %s!\n",argv[0],cfftsign);
  cl_get( "recon|r", "%option %s[%]", "modulus", crecon );
  c.complex = strcasecmp(crecon, "complex") == 0;
  c.phase = strcasecmp(crecon, "phase") == 0;

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
  if (strcmp(c.input_file, c.output_header_file) == 0)
    Abort( "Input and output files must be distinct." );
  if (c.input_file[0] != '/')
    {
      strcpy(short_name, c.input_file);
      if (getcwd(c.input_file, sizeof(Filename)-strlen(short_name)-3) == NULL)
	Abort("%s: pathname length exceeded", argv[0]);
      strcat(c.input_file, "/");
      strcat(c.input_file, short_name);
    }
  mInput = mri_open_dataset (c.input_file, MRI_READ);

  /* Check that program will function on data-set */
  if( !mri_has( mInput, "images" ) )
    Abort( "%s operates only on standard images.", argv[0] );
  if( mri_has( mInput, "images.dimensions" ) )
    {
      c.dv = ( !strcmp( mri_get_string( mInput, "images.dimensions" ), 
		      "vxyzt" ) )?
	mri_get_int( mInput, "images.extent.v" ): 
	( !strcmp( mri_get_string( mInput, "images.dimensions" ), "xyzt" ) )? 
	1: 0;
      if( ( c.dv < 1 ) || ( c.dv > 2 ) )
	Abort( "%s takes only reals or complex numbers of the form (v)xyzt.",
	       argv[0] );
    }
  else
    Abort( "%s does not have the images.dimensions key.", c.input_file );

  /* Set output dataset */
  mOutput = mri_copy_dataset(c.output_header_file, mInput);
  hist_add_cl( mOutput, argc, argv );
  mri_set_string(mOutput, "images.datatype", "float32");
  mri_set_string(mOutput, "images.file", c.output_data_file);
  mri_set_string(mOutput, "images.dimensions", "vxyzt");
  mri_set_int(mOutput, "images.extent.v", c.complex ? 2 : 1);
  if (c.complex)
    mri_set_string(mOutput, "images.description.v", "complex real/imaginary");
  if (mri_has(mInput,"images.description.x")) {
    const char* oldDesc= mri_get_string(mInput,"images.description.x");
    if (strstr(oldDesc,"ungridded"))
      Warning(1,"%s: Warning: input X dimension is ungridded!\n",argv[0]);
    if (strstr(oldDesc,"image-space")) {
      mri_set_string(mOutput,"images.description.x","gridded k-space");
    }
    else if (strstr(oldDesc,"k-space")) {
      mri_set_string(mOutput,"images.description.x","gridded image-space");
    }
  }
  if (mri_has(mInput,"images.description.y")) {
    const char* oldDesc= mri_get_string(mInput,"images.description.y");
    if (strstr(oldDesc,"ungridded"))
      Warning(1,"%s: Warning: input Y dimension is ungridded!\n",argv[0]);
    if (strstr(oldDesc,"image-space")) {
      mri_set_string(mOutput,"images.description.y","gridded k-space");
    }
    else if (strstr(oldDesc,"k-space")) {
      mri_set_string(mOutput,"images.description.y","gridded image-space");
    }
  }
  mri_close_dataset(mOutput);

  /* Set parameters in local variables */
  if( !mri_has( mInput, "images.extent.t" ) ||
      !mri_has( mInput, "images.extent.x" ) ||
      !mri_has( mInput, "images.extent.y" ) ||
      !mri_has( mInput, "images.extent.z" ) )
    Abort( "images.extent key(s) missing from header." );
  c.dt = mri_get_int( mInput, "images.extent.t" );
  c.dx = mri_get_int( mInput, "images.extent.x" );
  c.dy = mri_get_int( mInput, "images.extent.y" );
  c.dz = mri_get_int( mInput, "images.extent.z" );
  if (c.dt <= 0 || c.dx <= 0 || c.dy <= 0 || c.dz <= 0)
    Abort( "images.extent key(s) is non-positive." );

  mri_close_dataset(mInput);

  /* set context once and for all */
  par_set_context();

  /* go through all times */
  for (t.t = 0; t.t < c.dt; ++t.t)
    /* go through all slices */
    for (t.z = 0; t.z < c.dz; ++t.z)
      par_delegate_task();

  par_finish();

  Report( "#      Image reconstruction complete.\n" );
  exit(0);
}

void MasterResult (int task_number)
{
  static int results = 0;
  int z;
  int t;

  z = results % c.dz;
  t = results / c.dz;
  ++results;

  /* Print progress report */
  if (t == 0 && z == 0)
    Report( "      " );

  if (z == c.dz-1) {
    if ((t == (c.dt - 1)) || ((t+1) % 60 == 0))
      {
	Report( "# %ld\n      ", (long) ( t + 1 ) );
      }
    else
      Report("#");
  }

  if (results==c.dz*c.dt) Report("\n");
}


/* WORKER PROCEDURES */

void
WorkerFinalize()
{
  if (wInput != NULL) mri_close_dataset(wInput);
  if (wOutput != NULL) mri_close_dataset(wOutput);
  if (image != NULL)
    FreeMatrix(image);
  if (real_image != NULL)
    FreeMatrix(real_image);
  if (recon_image != NULL)
    FreeMatrix(recon_image);
}

void
WorkerContext ()
{
  static int first = TRUE;
  static int old_dx = 0, old_dy = 0;

  if (first)
    {
      wInput = mri_open_dataset(c.input_file, MRI_READ);
      wOutput = mri_open_dataset(c.output_header_file, MRI_MODIFY_DATA);
      first = FALSE;
    }

  if (c.dx != old_dx || c.dy != old_dy)
    {
      /* free up old storage */
      if (image != NULL)
	FreeMatrix(image);
      if (real_image != NULL)
	FreeMatrix(real_image);
      if (recon_image != NULL)
	FreeMatrix(recon_image);

      /* Allocate parameter and image storage */
      image = Matrix( c.dy, c.dx, FComplex );
      real_image = Matrix( c.dy, c.dx, float );
      recon_image = Matrix( c.dy, c.dx, float );

      old_dx = c.dx;
      old_dy = c.dy;
    }
}

void
WorkerTask ()
{
  register int x, y;
  FComplex *orig_image;
  float *orig_real_image;

  /* Get image */
  if( c.dv == 2 )
    {
      orig_image = (FComplex *)
	mri_get_image( wInput, t.t, t.z, MRI_COMPLEX_FLOAT );
      memcpy(*image, orig_image, c.dy*c.dx*sizeof(FComplex));
    }
  else
    {
      orig_real_image = (float *) mri_get_image( wInput, t.t, t.z, MRI_FLOAT );
      /* Make real-valued image into complex-valued for FFT */
      for (y = 0; y < c.dy; y++)
	for (x = 0; x < c.dx; x++)
	  {
	    image[y][x].real = *orig_real_image++;
	    image[y][x].imag = 0.0;
	  }
    }

  /* Fourier transform image */
  fft2d( image, (long)(c.dx), (long)(c.dy), (long)(c.fftsign), '2', 
      (long)0, (long)0 );

  if (c.complex)
    /* Complex-valued output, so just write FFT'ed image */
      mri_set_image( wOutput, t.t, t.z, MRI_COMPLEX_FLOAT, *image );
  else
    {
      /* Calculate modulus/phase of image and write out */
      if (c.phase)
	for( y = 0; y < c.dy; y++ )
	  for( x = 0; x < c.dx; x++ )
	    recon_image[y][x] = Phase( image[y][x] );
      else
	for( y = 0; y < c.dy; y++ )
	  for( x = 0; x < c.dx; x++ )
	    recon_image[y][x] = Modulus( image[y][x] );
      mri_set_image( wOutput, t.t, t.z, MRI_FLOAT, *recon_image );
    }
}


/* UTILITY FUNCTIONS */

void
PackContext ()
{
  par_pkstr(c.input_file);
  par_pkstr(c.output_header_file);
  par_pkstr(c.output_data_file);
  par_pkint(c.dv);
  par_pkint(c.dx);
  par_pkint(c.dy);
  par_pkint(c.dz);
  par_pkint(c.dt);
  par_pkint(c.fftsign);
  par_pkint(c.complex);
  par_pkint(c.phase);
}

void
UnpackContext ()
{
  par_upkstr(c.input_file);
  par_upkstr(c.output_header_file);
  par_upkstr(c.output_data_file);
  c.dv= par_upkint();
  c.dx= par_upkint();
  c.dy= par_upkint();
  c.dz= par_upkint();
  c.dt= par_upkint();
  c.fftsign= par_upkint();
  c.complex= par_upkint();
  c.phase= par_upkint();
}

void
PackTask ()
{
  par_pkint(t.z);
  par_pkint(t.t);
}

void
UnpackTask ()
{
  t.z= par_upkint();
  t.t= par_upkint();
}


