/************************************************************
 *
 *    estireg.c - Parallel registration parameter estimator
 *		 for within-slice, 2-D, rigid-body motion
 *
 *  Permission is hereby granted to any individual or
 *  institution for use, copying, or redistribution of
 *  this code and associated documentation, provided
 *  that such code and documentation are not sold for
 *  profit and the following copyright notice is retained
 *  in the code and documentation:
 *     Copyright (c) 1995,1996 Department of Statistics,
 *                           Carnegie Mellon University, and
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
 *                           Pittsburgh Supercomputing Center
 *  Original programming by Mark Fitzgerald  2/95
 *	5/96: Pittsburgh Format, Mark Fitzgerald
 * 	6/96: Parallelization, Greg Hood
 *
 *  Usage:
 *	estireg [-input <input-header-file>]
 *              [-parameters <parameter-file>]
 *	        [-fixed <fixed-image-number>]
 *    where:
 *         -input specifies the dataset (in Pittsburgh format)
 *	   -parameters specifies the parameter file to write
 *	   -fixed specifies the image number to register to
 *             (default is middle image of sequence)
*************************************************************/

/*
 * TO DO:
 *	The FFT routine used here could be replaced with highly optimized
 *	    routines on certain architectures.  These would have to be
 *	    incorporated using #ifdefs so that the code remains portable.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"
#include "par.h"
#include "array.h"
#include "misc.h"
#include "acct.h"


static char rcsid[] = "$Id: estireg.c,v 1.24 2007/04/19 22:31:02 welling Exp $";

void smooth_image( float** image, float** sm_image, int nx, int ny );
void regist( RegPars *regpar, float* mseval );
float mse( float par[3], void* userHook );
void restrt( float par[3], void* userHook );
void bil_shift_rot( RegPars par );

typedef struct Context {
  /* NOTE: any new fields added to this struct should
     be also added to PackContext and UnpackContext */ 
  Filename input_file;		/* the images to work on */
  Filename fixed_file;		/* the image that is considered fixed
				   (if not in input_file) */

  int dv;			/* # of values along the V dimension */
  int dx;			/* # of pixels along the X dimension */
  int dy;			/* # of pixels along the Y dimension */
  int dz;			/* # of slices along the Z dimension */
  int dt;			/* # of images along the T dimension */

  int z;			/* slice number to work on */
  int fixed;			/* # of image to be used as a reference */
} Context;

typedef struct Task {
  /* NOTE: any new fields added to this struct should
     be also added to PackTask and UnpackTask */ 
  int t;			/* image number to work on */
} Task;

typedef struct Result {
  /* NOTE: any new fields added to this struct should
     be also added to PackResult and UnpackResult */
  int t;
  int z;
  RegPars reg_pars;
  float mean_squared_error; /* negative means optimization failed */
} Result;

/* GLOBAL VARIABLES FOR MASTER & SLAVE */
Task t;
Context c;
Result r;

/* GLOBAL VARIABLES FOR MASTER */
static MRI_Dataset *mInput;
static MRI_Dataset *fInput;
static FILE *pOutput;
static Filename parameter_file;	/* the file where registration
				   parameters will be written */
static char fixed_string[128];	/* the string holding the image number
				   which is considered fixed */
RegPars **regpar = NULL;	/* the array holding the registration
				   parameters being computed */
float **mseval = NULL;		/* the array holding the mean square
				   error values for the parameters being
				   computed */
static int min_slice = 0;	/* the minimum slice number in the dataset */

/* GLOBAL VARIABLES FOR SLAVE */
static MRI_Dataset *sInput;
static MRI_Dataset *xInput;
static RegPars adjpar;
static float **fixed_image = NULL, **moved_fixed_image = NULL;
static float **sm_moved_fixed_image = NULL;
static FComplex **comp_fixed_image = NULL, **comp_moved_fixed_image = NULL;
static float **reg_image = NULL, **sm_reg_image = NULL;
static float **moved_sm_reg_image = NULL;
static short **check = NULL;

/* FORWARD DECLARATIONS */
void MasterTask (int argc, char **argv, char **envp);
void ReadEnvironment ();
void ReadArguments (int argc, char **argv);
void MasterResult (int task_number);
void SlaveContext ();
void SlaveTask ();
void SlaveFinalize();
void PackContext ();
void UnpackContext ();
void PackTask ();
void UnpackTask ();
void PackResult ();
void UnpackResult ();

static void writeParFileHeader(FILE* ofile)
{
  fprintf(ofile,"##Format: order:index_tz, type:raw\n");
  fprintf(ofile,"##Format: names:(xshift,yshift,rotation,mse)\n");
}

int
main (int argc,
      char **argv,
      char **envp)
{
  par_process(argc, argv, envp,
	      MasterTask, MasterResult,
	      SlaveContext, SlaveTask,
	      SlaveFinalize,
	      PackContext, UnpackContext,
	      PackTask, UnpackTask,
	      PackResult, UnpackResult);
  return(0);
}

/* MASTER PROCEDURES */

void
MasterTask (int argc,
	    char **argv,
	    char **envp)
{
  unsigned char **missing = NULL;
  int i;
  int temp_fixed= 0;
  int tt, zz;
  Filename short_name;
  int fixed_file_times;

  verbose = TRUE;
  Report("# estireg: version 4.0\n");

  /* Check to see if help was requested */
  if( ( argc > 1 ) && !strcmp( argv[1], "-help" ) )
    {
      if( argc == 2 )
	Help( "selecttopic" );
      else
	Help( argv[2] );
    }

  ReadEnvironment();
  ReadArguments(argc, argv);

  /* Open input dataset */
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
  if( !( mInput -> std_images ) )
    Abort( "%s operates only on standard images.", argv[0] );
  if( !mri_has( mInput, "images.dimensions" ) || 
      !( !strcmp( mri_get_string( mInput, "images.dimensions" ), "vxyzt" ) ||
	 !strcmp( mri_get_string( mInput, "images.dimensions" ), "xyzt" ) ) )
    Abort( "%s operates only on images in order vxyzt or xyzt.", argv[0] );
  if( !mri_has( mInput, "images.extent.v" ))
    c.dv = 1;
  else if ((c.dv = mri_get_int( mInput, "images.extent.v")) != 1)
    Abort( "%s operates only on real-valued images (v = 1).", argv[0] );

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
  if (c.dx<5 || c.dy<5) Abort("X and Y extents must be at least 5!\n");

  if (mri_has( mInput, "images.min.z" ))
    min_slice = mri_get_int ( mInput, "images.min.z" );

  if (c.fixed_file[0] != '\0')
    {
      if (c.fixed_file[0] != '/')
	{
	  strcpy(short_name, c.fixed_file);
	  if (getcwd(c.fixed_file, sizeof(Filename)-strlen(short_name)-3) == NULL)
	    Abort("%s: Fixed file pathname length exceeded", argv[0]);
	  strcat(c.fixed_file, "/");
	  strcat(c.fixed_file, short_name);
	}
      fInput = mri_open_dataset (c.fixed_file, MRI_READ);

      /* Check that program will function on data-set */
      if( !( fInput -> std_images ) )
	Abort( "%s operates only with standard fixed images.", argv[0] );
      if( !mri_has( fInput, "images.dimensions" ) || 
	  strcmp( mri_get_string( fInput, "images.dimensions" ), "vxyzt" ) != 0 &&
	  strcmp( mri_get_string( fInput, "images.dimensions" ), "xyzt" ) != 0 &&
	  strcmp( mri_get_string( fInput, "images.dimensions" ), "vxyz" ) != 0 &&
	  strcmp( mri_get_string( fInput, "images.dimensions" ), "xyz" ) != 0)
	Abort( "%s operates only on fixed images in order (v)xyz(t).", argv[0] );
      if (mri_has( fInput, "images.extent.v") &&
	  mri_get_int( fInput, "images.extent.v") != 1)
	Abort( "%s operates only on real-valued images (v = 1).", argv[0] );
      if (!mri_has( fInput, "images.extent.x") ||
	  mri_get_int(fInput, "images.extent.x") != c.dx)
	Abort( "%s: x dimension of the registered and fixed images must be equal.\n");
      if (!mri_has( fInput, "images.extent.y") ||
	  mri_get_int(fInput, "images.extent.y") != c.dy)
	Abort( "%s: y dimension of the registered and fixed images must be equal.\n");
      if (!mri_has( fInput, "images.extent.z") ||
	  mri_get_int(fInput, "images.extent.z") != c.dz)
	Abort( "%s: z dimension of the registered and fixed images must be equal.\n");
      if (mri_has( fInput, "images.extent.t"))
	fixed_file_times = mri_get_int(fInput, "images.extent.t");
      else
	fixed_file_times = 1;

      mri_close_dataset(fInput);
    }
  else
    fixed_file_times = c.dt;

  /* Determine which image is to be considered fixed */
  if (strcasecmp(fixed_string, "middle" ) == 0 ||
      strcasecmp(fixed_string, "center" ) == 0)
    c.fixed = fixed_file_times / 2;
  else
    {
      c.fixed = atoi( fixed_string );
      c.fixed = ( c.fixed < 0 ) ? 0 : 
	( c.fixed >= fixed_file_times ) ? (fixed_file_times-1) : c.fixed;
    }

  /* Read/Create missing image indicators */
  missing = get_missing(mInput);

  /* Allocate parameter storage */
  regpar = Matrix( c.dt, c.dz, RegPars );
  mseval = Matrix( c.dt, c.dz, float );

  /* open the parameter file for interim output;
     this file will be rewritten later after all
     results are in hand */
  pOutput = efopen(parameter_file, "w");
  writeParFileHeader(pOutput);

  /* go through all slices */
  for (c.z = 0; c.z < c.dz; ++c.z)
    {
      int all_images_missing= 0;
      /* Find first non-missing image occurring */
      /*   from specified fixed image on        */      
      if (c.fixed_file[0] == '\0')
	{
	  for (i = 0;  i < c.dt; ++i)
	    {
	      temp_fixed = (c.fixed + i) % c.dt;
	      if (!missing[temp_fixed][c.z])
		break;
	    }
	  if (i >= c.dt)
	    {
	      Warning(1,"%s: Error: All images for slice %d are missing!\n",
		      argv[0],c.z);
	    }
	  c.fixed = temp_fixed;
	}
      par_set_context();

      /* assign each image to a slave */
      for (t.t = 0; t.t < c.dt; ++t.t)
	if (c.fixed_file[0] == '\0' && t.t == c.fixed || missing[t.t][c.z])
	  {
	    /* By definition, registration is 0 at fixed image */
	    regpar[t.t][c.z].x_shift = regpar[t.t][c.z].y_shift = 
	      regpar[t.t][c.z].rotation = 0.0;
	    mseval[t.t][c.z] = missing[t.t][c.z] ? -1.0 : 0.0;
	  }
	else
	  par_delegate_task();
    }

  par_finish();
  mri_close_dataset(mInput);

  /* Re-write parameters in proper order */
  efclose( pOutput );
  pOutput = efopen( parameter_file, "w" );
  writeParFileHeader(pOutput);
  for( tt = 0; tt < c.dt; tt++ )
    for( zz = 0; zz < c.dz; zz++ )
      fprintf( pOutput, "%5d %3d %14.9g %14.9g %14.9g %14.6g\n",
	       tt, zz + min_slice,
	       regpar[tt][zz].x_shift, regpar[tt][zz].y_shift,
	       regpar[tt][zz].rotation, mseval[tt][zz] );
  efclose( pOutput );


  Message( "#      Registration parameters estimated.\n" );
}

void MasterResult (int task_number)
{
  static int results= 0;

  /* Flush most recent value to parameter file */
  fprintf(pOutput, "%d %d %g %g %g %g\n", r.t, r.z + min_slice, 
	  r.reg_pars.x_shift, r.reg_pars.y_shift,
	  r.reg_pars.rotation, r.mean_squared_error);
  /* We use a terse report rather than the following verbose one 
   * Report("Results of task %d: %ld %ld %g %g %g %g\n", 
   *  task_number, r.t, r.z, r.reg_pars.x_shift, r.reg_pars.y_shift,
   *  r.reg_pars.rotation, r.mean_squared_error);
   */

  ++results;
  if (!(results % 60) || (results == c.dz*(c.dt-1))) 
    Report("# %ld\n", (long)results);
  else 
    Report("#");

  if (results == (c.dt-1)*c.dz) Report("\n");

  regpar[r.t][r.z] = r.reg_pars;
  mseval[r.t][r.z] = r.mean_squared_error;

  /* Print progress report (verbose version) */
  /*  if (r.t != 0)
    Report( "\n      #Slice %ld\n      ", r.z );
  if ((r.t == (c.dt - 1)) || ((r.t + 1) % 60) == 0)
    {
      Report( "# %ld\n      ", (long) (r.t + 1));
      if (r.t != c.dt - 1)
	Report("      ");
    }
  else
    Report( "#" ); 
  */
}

void
ReadEnvironment ()
{
  char *s;

  if ((s = getenv("F_ESTIREG_INP_FILE")) != NULL)
    StringCopy(c.input_file, s, sizeof(Filename));
  else
    StringCopy(c.input_file, "input.mri", sizeof(Filename));
  if ((s = getenv("F_ESTIREG_PAR_FILE")) != NULL)
    StringCopy(parameter_file, s, sizeof(Filename));
  else
    StringCopy(parameter_file, "reg.par", sizeof(Filename));
  if ((s = getenv("F_ESTIREG_FIX_FILE")) != NULL)
    StringCopy(c.fixed_file, s, sizeof(Filename));
  else
    c.fixed_file[0] = '\0';
  if ((s = getenv("F_ESTIREG_FIXED")) != NULL)
    StringCopy(fixed_string, s, sizeof(fixed_string));
  else
    StringCopy(fixed_string, "middle", sizeof(fixed_string));
  if ((s = getenv("F_ESTIREG_HOSTS")) != NULL)
    par_set_hosts(s);
}

void
ReadArguments (int argc,
	       char **argv)
{
  cl_scan(argc, argv);
  cl_get("input|i", "%option %s[%]", c.input_file, c.input_file);
  cl_get("parameters|p", "%option %s[%]", parameter_file, parameter_file);
  cl_get("fixed|f", "%option %s[%]", fixed_string, fixed_string);
  cl_get("fixed_file|align|a", "%option %s[%]", c.fixed_file, c.fixed_file);

  if (cl_cleanup_check()) {
    int i;
    fprintf(stderr,"%s: invalid argument in command line:\n    ",argv[0]);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(-1);
  }
}

/* SLAVE PROCEDURES */

void
SlaveFinalize()
{
  if (xInput != NULL && xInput != sInput) mri_close_dataset(xInput);
  if (sInput != NULL) mri_close_dataset(sInput);
  if (fixed_image != NULL)
    free(fixed_image);
  if (comp_fixed_image != NULL)
    FreeMatrix(comp_fixed_image);
  if (moved_fixed_image != NULL)
    FreeMatrix(moved_fixed_image);
  if (sm_moved_fixed_image != NULL)
    FreeMatrix(sm_moved_fixed_image);
  if (reg_image != NULL)
    free(reg_image);
  if (sm_reg_image != NULL)
    FreeMatrix(sm_reg_image);
  if (moved_sm_reg_image != NULL)
    FreeMatrix(moved_sm_reg_image);
  if (check != NULL)
    FreeMatrix(check);
}

void
SlaveContext ()
{
  static int first = TRUE;
  static int old_dx = 0, old_dy = 0, old_z = -1;
  int x, y;

  if (first)
    {
      sInput = mri_open_dataset(c.input_file, MRI_READ);
      if (c.fixed_file[0] != '\0')
	xInput = mri_open_dataset(c.fixed_file, MRI_READ);
      else
	xInput = sInput;
      first = FALSE;
    }

  if (c.dx != old_dx || c.dy != old_dy)
    {
      /* free up old storage */
      if (fixed_image != NULL)
	free(fixed_image);
      if (comp_fixed_image != NULL)
	FreeMatrix(comp_fixed_image);
      if (moved_fixed_image != NULL)
	FreeMatrix(moved_fixed_image);
      if (sm_moved_fixed_image != NULL)
	FreeMatrix(sm_moved_fixed_image);
      if (reg_image != NULL)
	free(reg_image);
      if (sm_reg_image != NULL)
	FreeMatrix(sm_reg_image);
      if (moved_sm_reg_image != NULL)
	FreeMatrix(moved_sm_reg_image);
      if (check != NULL)
	FreeMatrix(check);

      /* Allocate parameter and image storage */
      fixed_image = (float **) emalloc( c.dy * sizeof(float *) );
      comp_fixed_image = Matrix( c.dy, c.dx, FComplex );
      comp_moved_fixed_image = Matrix( c.dy, c.dx, FComplex );
      moved_fixed_image = Matrix( c.dy, c.dx, float );
      sm_moved_fixed_image = Matrix( c.dy, c.dx, float );
      reg_image = (float **) emalloc( c.dy * sizeof(float *) );
      sm_reg_image = Matrix( c.dy, c.dx, float );
      moved_sm_reg_image = Matrix( c.dy, c.dx, float );
      check = Matrix( c.dy, c.dx, short );

      old_dx = c.dx;
      old_dy = c.dy;
    }

  if (c.z != old_z)
    {
      /* set fixed image for the slice */
      *fixed_image = (float *)
	mri_get_image( xInput, c.fixed, c.z, MRI_FLOAT );
      realign_matrix( (void **) fixed_image, c.dy,
		      (long) ( c.dx * sizeof(float) ) );
      adjpar.rotation = 0.0;
      smooth_image( fixed_image, sm_moved_fixed_image, c.dx, c.dy );

      /* Set Fourier version of fixed image (set imaginary part to 0) */
      for (y = 0; y < c.dy; y++ )
	for (x = 0; x < c.dx; x++ )
	  {
	    comp_fixed_image[y][x].real = fixed_image[y][x];
	    comp_fixed_image[y][x].imag = 0.0;
	  }
      old_z = c.z;
    }
}

void
SlaveTask ()
{
  /* Get register image */
  *reg_image = (float *)
    mri_get_image( sInput, t.t, c.z, MRI_FLOAT);
  realign_matrix( (void **) reg_image, c.dy,
		  (long) ( c.dx * sizeof(float) ) );
  smooth_image( reg_image, sm_reg_image, c.dx, c.dy );
  
  /* start each registration at (0.0, 0.0, 0.0) */
  r.t = t.t;
  r.z = c.z;
  r.reg_pars.x_shift = 0.0;
  r.reg_pars.y_shift = 0.0;
  r.reg_pars.rotation = 0.0;
	      
  /* Estimate parameters */
  regist(&r.reg_pars, &r.mean_squared_error);
}

/* Function to slightly blur images (no blurring done at border of image) */
void smooth_image( float** image, float** sm_image, int nx, int ny )
{
  int x, y;

  for( y = 0; y < ny; y++ )
    for( x = 0; x < nx; x++ )
      {
	if( ( x == 0 ) || ( y == 0 ) || ( x == ( nx - 1 ) ) || ( y == ( ny - 1 ) ) )
	  sm_image[y][x] = image[y][x];
	else
	  sm_image[y][x] = 0.5 * image[y][x] + 
	    0.125 * ( image[y][x-1] + image[y][x+1] +
		      image[y-1][x] + image[y+1][x] );
      }
  return;
}

/* Function to set-up and call Nelder-Mead minimization */
void regist( RegPars *regpar, float* mseval )
{
  int numpar = 3, steps_per_conv_check = 3;
  int max_iter = 400, max_restart = 3;
  int num_iter, num_restart, return_cond;
#ifdef never
  float stopping_val = 0.1;
#endif
  float stopping_val = 0.02;
  RegPars start_vals, scale;
  float tcrit_sqr= 2.353*2.353; /* 95% confidence on 3 degrees of freedom */
  start_vals.x_shift = (*regpar).x_shift;
  start_vals.y_shift = (*regpar).y_shift;
  start_vals.rotation = (*regpar).rotation;
  scale.x_shift = scale.y_shift = 1.0;
  scale.rotation = 2.0;

#ifdef never
  fprintf(stderr,"instance %d: shifts are %g %f\n",par_instance(),start_vals.x_shift,regpar->x_shift);
#endif
  nelmin_t( mse, restrt,
	  &numpar, (float *) &start_vals, (float *) regpar, mseval,
	  &stopping_val, (float *) (&scale), &steps_per_conv_check,
	  &max_iter, &num_iter, &num_restart, &return_cond, &max_restart,
	  &tcrit_sqr, NULL );
  if (return_cond != 0) {
    if (return_cond==2)
      Warning(1,"estireg: Nelder-Mead convergence failed!\n");
    else 
      Warning(1,"estireg: Nelder-Mead call illegal value!\n");
    Warning(1,"         failed at iter = %d of %d, restart= %d of %d\n",
	    num_iter, max_iter, num_restart, max_restart);
    *mseval= -1.0; /* impossible value indicates error */
  }

  return;
}


/* Returns mean squared-error between (adjusted) fixed image  */
/*   and (adjusted) register image, the criterion to be       */
/*   minimized by the Nelder-Mead minization for registration */
float mse( float par[3], void* userHook )
{
  RegPars tmppar;
  double sse;
  long x, y, count;
  register double tmp;

  /* Adjust amount to move register image, to account */
  /*   for movement of fixed image                    */
  tmppar.x_shift = par[0] + adjpar.x_shift;
  tmppar.y_shift = par[1] + adjpar.y_shift;
  tmppar.rotation = par[2] + adjpar.rotation;

  /* Move register image */
  bil_shift_rot( tmppar );

  /* Find mean-squared difference */
  sse = 0;
  count = 0;
  for( y = 2; y < ( c.dy - 2 ); y++ )
    for( x = 2; x < ( c.dx - 2 ); x++ )
      if( check[y][x] )
	{
	  tmp= moved_sm_reg_image[y][x] - sm_moved_fixed_image[y][x];
	  sse += tmp*tmp;
	  count++;
	}
  sse /= (double) count;

  return( (float) sse );
}


void bil_shift_rot( RegPars par )
{
  float rotate, cos_rot, sin_rot;
  float midx, midy;
  long x, y;
  float shiftx, shifty;
  float newx, newy;
  long lowerx, upperx, lowery, uppery;
  float distlx, distly;

  /* Convert rotation units from degrees to radians */
  /*   for trig functions                           */
  rotate = par.rotation * M_PI / 180.0;
  cos_rot = cos( rotate );
  sin_rot = sin( rotate );

  /* Calculate origin of x- and y-dimensions */
  midx = (long) ( c.dx / 2 );
  midy = (long) ( c.dy / 2 );

  for( y = 0; y < c.dy; y++ )
    {
      /* New y-coordinate after shift */
      shifty = (float) ( y - midy ) + par.y_shift;

      for( x = 0; x < c.dx; x++ )
	{
	  /* New x-coordinate after shift */
	  shiftx = (float) ( x - midx ) + par.x_shift;

	  /* New coordinates after rotation */
	  newx = shiftx * cos_rot - shifty * sin_rot + midx;
	  newy = shifty * cos_rot + shiftx * sin_rot + midy;

	  /* Find nearest grid-points and distances thereto */
	  lowerx = newx;
	  if (newx - lowerx < 0.0) lowerx -= 1;
	  lowery = newy;
	  if (newy - lowery < 0.0) lowery -= 1;
	  upperx = lowerx + 1;
	  if ( upperx == c.dx ) upperx = lowerx;
	  uppery = lowery + 1;
	  if ( uppery == c.dy ) uppery = lowery;
	  distlx = newx - lowerx;
	  distly = newy - lowery;

	  if( ( lowerx < 0 ) || ( upperx >= c.dx ) ||
	      ( lowery < 0 ) || ( uppery >= c.dy ) )
	    {
	      /* Moved out of image range */
	      moved_sm_reg_image[y][x] = 0.0;

	      /* Indicate pixel as missing */
	      check[y][x] = 0;
	    }
	  else
	    {
	      /* Bilinear interpolation */
	      moved_sm_reg_image[y][x] = 
		( 1.0 - distly ) * ( ( 1.0 - distlx ) * sm_reg_image[lowery][lowerx] +
				     distlx * sm_reg_image[lowery][upperx] ) +
		distly * ( ( 1.0 - distlx ) * sm_reg_image[uppery][lowerx] +
			   distlx * sm_reg_image[uppery][upperx] );

	      /* Indicate pixel as present */
	      check[y][x] = 1;
	    }

	}
    }

  return;
}


void restrt( float par[3], void* userHook )
{
  long x, y;

  /* Currently, restarts are for rotations only */
  adjpar.x_shift = adjpar.y_shift = 0.0;

  if( adjpar.rotation == -par[2] )
    /* Rotation parameter hasn't changed --- no need */
    /*   to modify fixed image                       */
    return;

  else
    {
      /* Rotate the fixed image, since a new "best"           */
      /*   rotation has been found (or simply an initial one) */
      adjpar.rotation = -par[2];
      fourier_shift_rot( adjpar, comp_fixed_image, comp_moved_fixed_image,
			 c.dx, c.dy, 'i' );

      /* Convert back to real-valued */
      for( y = 0; y < c.dy; y++ )
	for( x = 0; x < c.dx; x++ )
	  moved_fixed_image[y][x] = Modulus( comp_moved_fixed_image[y][x] );

      /* Re-smooth moved image */
      smooth_image( moved_fixed_image, sm_moved_fixed_image, c.dx, c.dy );
      return;
    }

}

/* UTILITY FUNCTIONS */

void
PackContext ()
{
  par_pkstr(c.input_file);
  par_pkstr(c.fixed_file);
  par_pkint(c.dv);
  par_pkint(c.dx);
  par_pkint(c.dy);
  par_pkint(c.dz);
  par_pkint(c.dt);
  par_pkint(c.z);
  par_pkint(c.fixed);
}

void
UnpackContext ()
{
  par_upkstr(c.input_file);
  par_upkstr(c.fixed_file);
  c.dv= par_upkint();
  c.dx= par_upkint();
  c.dy= par_upkint();
  c.dz= par_upkint();
  c.dt= par_upkint();
  c.z= par_upkint();
  c.fixed= par_upkint();
}

void
PackTask ()
{
  par_pkint(t.t);
}

void
UnpackTask ()
{
  t.t= par_upkint();
}

void
PackResult ()
{
  par_pkint(r.t);
  par_pkint(r.z);
  par_pkfloat(r.reg_pars.x_shift);
  par_pkfloat(r.reg_pars.y_shift);
  par_pkfloat(r.reg_pars.rotation);
  par_pkfloat(r.mean_squared_error);
}

void
UnpackResult ()
{
  r.t= par_upkint();
  r.z= par_upkint();
  r.reg_pars.x_shift= par_upkfloat();
  r.reg_pars.y_shift= par_upkfloat();
  r.reg_pars.rotation= par_upkfloat();
  r.mean_squared_error= par_upkfloat();
}

