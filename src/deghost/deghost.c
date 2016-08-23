/************************************************************
 *                                                          *
 *  deghost.c                                               *
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
 *  Original programming by Mark Fitzgerald  6-96           *
 *     7/97: Parallelization, Greg Hood                     *
 *     9/97: Incorporated algorithmic changes by Jana Asher *
 *                                        (Greg Hood, PSC)  *
 ************************************************************/
/*************************************************************

  DESCRIPTION OF DEGHOST

  deghost is used to reduce ghosting effects in images

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
#include "../fmri/lapack.h"

static char rcsid[] = "$Id: deghost.c,v 1.26 2007/03/21 23:46:16 welling Exp $";

static void opp_phase_shift( FComplex** image, long dy, long dx, 
			     float phase, FComplex** corr_image );
static float phase_est();
static void flip_rows(FComplex** image);
static double dbl_eval_phase( double v );
static float eval_phase(float* phase);
static void restrt();
static int parse_cl_reverse( const char* s, int* revOdd, int* revEven );
static void parse_mrihdr_reverse( MRI_Dataset* ds, const char* chunk, 
				  int* revOdd, int* revEven );

typedef struct Context {
  /* NOTE: any new fields added to this struct should
     be also added to PackContext and UnpackContext */ 
  Filename input_file;		/* the images to work on */
  Filename output_file;	        /* the name of the output file */

  int dv;			/* # of values along the V dimension */
  int dx;			/* # of pixels along the X dimension */
  int dy;			/* # of pixels along the Y dimension */
  int dz;			/* # of slices along the Z dimension */
  int dt;			/* # of images along the T dimension */
  int flip_odd;                 /* # whether or not to flip odd rows for FFT */
  int flip_even;                /* # whether or not to flip even rows fr FFT */
} Context;

#define PHASE_EST	0
#define OPP_PHASE_SHIFT	1

typedef enum { MODE_ESTIMATE, MODE_APPLY, MODE_BOTH } Mode;

typedef struct Task {
  /* NOTE: any new fields added to this struct should
     be also added to PackTask and UnpackTask */ 
  int type;		/* type of Task to perform:
			      PHASE_EST = estimate phase shift
			      OPP_PHASE_SHIFT = apply opposite phase shift */
  int z;		/* slice number to work on */
  int t;		/* image number to work on */
  float phase;		/* phase adjustment for that slice */
} Task;

typedef struct Result {
  /* NOTE: any new fields added to this struct should
     be also added to PackResult and UnpackResult */
  int type;
  int z;
  int t;
  float phase;
} Result;

/* GLOBAL VARIABLES FOR MASTER & WORKER */
static char* progname;
Task t;
Context c;
Result r;

/* GLOBAL VARIABLES FOR MASTER */
static MRI_Dataset *mInput= NULL;
static MRI_Dataset *mOutput= NULL;
static float **phase;
static float **smphase;
static Mode mode;

/* GLOBAL VARIABLES FOR WORKER */
static MRI_Dataset *wInput= NULL;
static MRI_Dataset *wOutput= NULL;
static FComplex **tmp_image = NULL;
static FComplex **image = NULL;
static FComplex **corr_image = NULL;

/* FORWARD DECLARATIONS */
void MasterTask (int argc, char **argv, char **envp);
void ReadArguments (int argc, char **argv);
void MasterResult (int task_number);
void WorkerContext ();
void WorkerTask ();
void WorkerFinalize ();
void PackContext ();
void UnpackContext ();
void PackTask ();
void UnpackTask ();
void PackResult ();
void UnpackResult ();

int main (int argc,
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
  return 0;
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
	    char **argv,
	    char **envp)
{
  char infile[512], hdrfile[512], outfile[512];
  char parfile[512], smparfile[512], cphase[512];
  Filename short_name;
  float iphase;
  int fixed_phase_flag= 0;
  float weight, total_weight;
  float distance;
  FILE *par = NULL;
  FILE *smpar = NULL;
  long n;
  int i;
  float bandwidth = 45.0;
  float kpar = 60.0;
  char reverse_string[256];
  int reverse_set= 0;
  char mode_string[256];
  int mode_set= 0;
  Smoother* smoother= NULL;
  unsigned char** missing= NULL;
  char* dimstr;

  verbose = TRUE;
  /* Print version number */
  Message( "# %s\n", rcsid );

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  sm_init();
  sm_set_params( SM_GAUSSIAN, bandwidth, kpar, 0.0, NULL );

  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Deprecate old options */

   if (cl_present( "dataout|d" ))
     Abort ("Option dataout|d has been replaced by infile outfile format.  Please see help file.\n");
   if (cl_present( "headerout|h" ))
     Abort ("Option headerout|h has been replaced by infile outfile format.  Please see help file.\n");
   if (cl_present( "input|i" ))
     Abort ("Option input|i has been replaced by infile outfile format.  Please see help file.\n");
   if (cl_present( "parameters|p" ))       
     Abort ("Option parameters|p has been replaced by estimates|est|e.  Please see help file.\n");
    if (cl_present( "smoothedparameters|s" ))       
     Abort ("Option smoothedparameters|s has been replaced by smoothedestimates|sme.  Please see help file.\n");
 
  /* Get filenames */

  cl_get( "estimates|est|e", "%option %s[%]", "deghost.par", parfile );
  cl_get( "smoothedestimates|sme", "%option %s[%]", "deghost.smoothed.par", 
	  smparfile ); 
  cl_get( "phase|pha", "%option %s[%]", "byslice", cphase );
  reverse_set= cl_get( "reverse|rev|r", "%option %s", reverse_string );
  mode_set= cl_get( "mode|mod", "%option %s", mode_string );

  sm_parse_cl_opts();

  if(!cl_get("", "%s", c.input_file)) {
    fprintf(stderr, "%s: Input file name not given.\n", argv[0]);
    exit(-1);
  }
  if(!cl_get("", "%s", c.output_file)) {
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

  /* Open input dataset */
  if( !strcmp(c.input_file, c.output_file ) )
    Abort( "Input and output files must be distinct." );
  mInput = mri_open_dataset(c.input_file, MRI_READ );

  /* Check that program will function on data-set */
  if( !mri_has( mInput, "images" ) )
    Abort( "%s operates only on standard images.\n", argv[0] );
  if ( !mri_has( mInput, "images.dimensions" ) )
    Abort( "%s has no images.dimensions tag!\n", argv[0] );
  dimstr= mri_get_string( mInput, "images.dimensions" );
  if (dimstr[0] != 'v' || !mri_has( mInput, "images.extent.v" )
      || (mri_get_int( mInput, "images.extent.v" ) != 2 ) )
    Abort( "%s operates only on complex-valued images (v = 2).\n", argv[0] );
  c.dv= 2;
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

  /* Parse the mode string */
  if (mode_set) {
    if (!strcasecmp(mode_string,"estimate")) {
      mode= MODE_ESTIMATE;
    }
    else if (!strcasecmp(mode_string,"apply")) {
      mode= MODE_APPLY;
    }
    else if (!strcasecmp(mode_string,"both")) {
      mode= MODE_BOTH;
    }
    else Abort("%s: unrecognized mode <%s>!\n",progname,mode_string);
  }
  else mode= MODE_BOTH;

  /* Decide how phase is to be determined */
  if( sscanf( cphase, "%f", &iphase ) == 1 ) {
    fixed_phase_flag = 1;
  }
  else if ( strcasecmp(cphase, "byslice") )
    Abort("%s: unrecognized phase string <%s>!\n",cphase);

  /* Read/Create missing image indicators */
  missing = get_missing( mInput );

  /* Set output dataset */
  if (mode != MODE_ESTIMATE) {
    mOutput = mri_copy_dataset( c.output_file, mInput );
    hist_add_cl( mOutput, argc, argv );
    mri_set_int( mOutput, "images.rowflip", 0 );
    mri_set_string( mOutput, "images.datatype", "float32" );
  }

  /* close datasets */
  mri_close_dataset( mInput );
  if (mOutput != NULL) mri_close_dataset( mOutput );

  /* Allocate parameter storage */
  if (mode==MODE_APPLY) {
    phase = NULL;
    smphase = Matrix(c.dz, c.dt, float);
  }
  else {
    phase = Matrix(c.dz, c.dt, float);
    smphase = Matrix(c.dz, c.dt, float);
  }

  /* Set context for workers */
  par_set_context();

  if (mode==MODE_APPLY) {
    /* 
     * Read phases from smoothed parameter file 
     */
    smpar = efopen( smparfile, "r" );
    for (t.z = 0; t.z < c.dz; ++t.z)
      for (t.t = 0; t.t < c.dt; ++t.t)
	fscanf(smpar, "%g\n", &smphase[t.z][t.t]);
    efclose(smpar);
  }
  else {

    /* 
     * Generate the unsmoothed phase shifts 
     */
    if( fixed_phase_flag ) {
      /* User-specified phase --- applies to all slices the same */
      for (t.t = 0; t.t < c.dt; t.t++)
	for (t.z = 0; t.z < c.dz; t.z++)
	  phase[t.z][t.t] = iphase;
    }
    else {
      /* Estimate phase shift */
      t.type = PHASE_EST;
      for (t.t = 0; t.t < c.dt; t.t++)
	for (t.z = 0; t.z < c.dz; t.z++)
	  par_delegate_task();
      while (par_tasks_outstanding() > 0) par_wait(1.0);
    }

    /*
     * Write the unsmoothed phase shifts
     */
    par = efopen( parfile, "w" );
    for (t.z = 0; t.z < c.dz; ++t.z)
      for (t.t = 0; t.t < c.dt; ++t.t)
	fprintf(par, "%g\n", phase[t.z][t.t]);
    efclose( par );
    Report("Finished writing phase estimates to %s\n", parfile);
    
    /*
     * Generate and write the smoothed phase shifts
     */
    smpar = efopen( smparfile, "w" );
    smoother= sm_create_smoother();
    for (t.z = 0; t.z < c.dz; ++t.z) {
      SM_SMOOTH( smoother, phase[t.z], smphase[t.z], c.dt, missing, t.z );
      for (t.t = 0; t.t < c.dt; ++t.t)
	fprintf(smpar, "%g\n", smphase[t.z][t.t]);
    }
    sm_destroy(smoother);
    efclose( smpar );
  }

  /* 
   * Loop through images, applying the smoothed phase shifts
   */
  if (mode != MODE_ESTIMATE) {
    t.type = OPP_PHASE_SHIFT;
    for (t.t = 0; t.t < c.dt; t.t++)
      for (t.z = 0; t.z < c.dz; t.z++) 
	{
	  t.phase = smphase[t.z][t.t];
	  par_delegate_task();
	}
  }
    
  par_finish();

  Message( "#      De-ghosting complete.\n" );
  exit(0);
}

void MasterResult (int task_number)
{
  static int results = 0;
  int n;

  if (r.type == PHASE_EST)
    phase[r.z][r.t] = r.phase;

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

void WorkerFinalize()
{
  if (wInput != NULL) mri_close_dataset(wInput);
  if (wOutput != NULL) mri_close_dataset(wOutput);
  if (image != NULL) free(image);
  if (corr_image != NULL) FreeMatrix(corr_image);
  if (tmp_image != NULL) FreeMatrix(tmp_image);
}

void
WorkerContext ()
{
  static int old_dx = 0, old_dy = 0;

  if (c.dx != old_dx || c.dy != old_dy)
    {
      /* free up old storage */
      if (image != NULL)
	free(image);
      if (corr_image != NULL)
	FreeMatrix(corr_image);
      if (tmp_image != NULL)
	FreeMatrix(tmp_image);

      /* Allocate parameter and image storage */
      image = (FComplex **) emalloc( c.dy * sizeof(FComplex *) );
      corr_image = Matrix( c.dy, c.dx, FComplex );
      tmp_image = Matrix( c.dy, c.dx, FComplex );

      old_dx = c.dx;
      old_dy = c.dy;
    }
}

void
WorkerTask ()
{
  long imageSize= c.dv*c.dx*c.dy;

  r.z = t.z;
  r.t = t.t;
  r.type= t.type;

  switch (t.type)
    {
    case PHASE_EST:

      if (wInput==NULL)
	wInput = mri_open_dataset(c.input_file, MRI_READ);
      *image = 
	(FComplex *)mri_get_chunk(wInput, "images", imageSize,
				  (t.t*c.dz+t.z)*imageSize, MRI_FLOAT);
      realign_matrix( (void **) image, c.dy, 
		      (long) ( c.dx * sizeof(FComplex) ) );
      if (c.flip_odd || c.flip_even) flip_rows( image );
      r.phase = phase_est();
      break;
    case OPP_PHASE_SHIFT:

      if (wInput==NULL)
	wInput = mri_open_dataset(c.input_file, MRI_READ);
      if (wOutput==NULL)
	wOutput = mri_open_dataset(c.output_file, MRI_MODIFY_DATA);

      /* Get image and re-align matrix pointers */
      *image = 
	(FComplex *)mri_get_chunk(wInput, "images", imageSize,
				  (t.t*c.dz+t.z)*imageSize, MRI_FLOAT);
      realign_matrix( (void **) image, c.dy,
			  (long) (c.dx * sizeof(FComplex) ) );

      /* flip rows if needed, de-ghost correction, and flip back */
      if (c.flip_odd || c.flip_even) flip_rows( image );
      opp_phase_shift( image, (long) c.dy, (long) c.dx, t.phase, corr_image );
	  
      /* Set output image */
      mri_set_chunk( wOutput, "images", imageSize,
		     (t.t*c.dz + t.z)*imageSize, MRI_FLOAT, *corr_image );

      r.phase = 0.0;
      break;
    }
}

/* Function which performs opposite phase shifts on odd and even lines */
static void opp_phase_shift( FComplex** image, long dy, long dx, 
			     float phase, FComplex** corr_image )
{
  long y, x;
  float cosph, sinph;

  cosph = cos( phase );
  sinph = sin( phase );

  for( y = 0; y < dy; y++ )
    if( y % 2 )
      for( x = 0; x < dx; x++ )
	{
	  corr_image[y][x].real = cosph * image[y][x].real -
	    sinph * image[y][x].imag;
	  corr_image[y][x].imag = cosph * image[y][x].imag +
	    sinph * image[y][x].real;
	}
    else
      for( x = 0; x < dx; x++ )
	{
	  corr_image[y][x].real = cosph * image[y][x].real +
	    sinph * image[y][x].imag;
	  corr_image[y][x].imag = cosph * image[y][x].imag -
	    sinph * image[y][x].real;
	}
  return;
}

static double dbl_eval_phase( double v )
{
  float fv= (float)v;
  float result= eval_phase( &fv );
  return (double)result;
}

/* Function which estimates the phase shift which best de-ghosts an image */
static float phase_est()
{
  static double machep= 0.0;
  double phase;
  double start_phase= 0.0;
  double start_step= 0.01;
  int max_iter = 400;

  if (machep==0.0) {
    machep= (2.0*DLAMCH("e"));
  }
  if ((1.0+machep) == 1.0) 
    Abort("%s: precision test failed!\n",progname);

  phase= (float)fmin1D( -1.0, 1.0, dbl_eval_phase, 0.0, sqrt(machep), 0);

  return( phase );
}


static float eval_phase( float* phase )
{
  long y, x;
  double sum = 0.0;

  opp_phase_shift( image, (long) c.dy, (long) c.dx, *phase, tmp_image );
  fft2d( tmp_image, (long) c.dx, (long) c.dy, -1, '2', 0, 0 );
  for( y = 0; y < 3; y++ )
    for( x = ( c.dx / 8 ); x < ( 7 * c.dx / 8 ); x++ )
      sum += Modulus( tmp_image[y][x] );
  for( y = ( c.dy - 3 ); y < c.dy; y++ )
    for( x = ( c.dx / 8 ); x < ( 7 * c.dx / 8 ); x++ )
      sum += Modulus( tmp_image[y][x] );

  return( (float) sum );
}

static void flip_rows(FComplex** image)
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

/* Dummy function */
static void restrt()
{
  return;
}

/* UTILITY FUNCTIONS */

void
PackContext ()
{
  par_pkstr(c.input_file);
  par_pkstr(c.output_file);
  par_pkint(c.dv);
  par_pkint(c.dx);
  par_pkint(c.dy);
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
  c.dv= par_upkint();
  c.dx= par_upkint();
  c.dy= par_upkint();
  c.dz= par_upkint();
  c.dt= par_upkint();
  c.flip_odd= par_upkint();
  c.flip_even= par_upkint();
}

void
PackTask ()
{
  par_pkint(t.type);
  par_pkint(t.z);
  par_pkint(t.t);
  par_pkfloat(t.phase);
}

void
UnpackTask ()
{
  t.type= par_upkint();
  t.z= par_upkint();
  t.t= par_upkint();
  t.phase= par_upkfloat();
}

void
PackResult ()
{
  par_pkint(r.type);
  par_pkint(r.z);
  par_pkint(r.t);
  par_pkfloat(r.phase);
}

void
UnpackResult ()
{
  r.type= par_upkint();
  r.z= par_upkint();
  r.t= par_upkint();
  r.phase= par_upkfloat();
}


