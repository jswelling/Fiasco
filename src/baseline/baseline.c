/************************************************************
 *                                                          *
 *  baseline.c                                              *
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
 ************************************************************/

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
  int dz;			/* # of slices along the Z dimension */
  int dt;			/* # of images along the T dimension */

  int reverse_odd;		/* TRUE if we should reverse the odd lines
				   in the image */
  int reverse_even;		/* TRUE if we should reverse the even lines
				   in the image */
  long baseline_area;		/* the area of region over which data is summed */
} Context;

typedef enum task_type_enum { TASK_ESTIMATE, TASK_CORRECT } TaskType;

typedef struct Task {
  /* NOTE: any new fields added to this struct should
     be also added to PackTask and UnpackTask */ 
  int t;			/* image number to work on */
  int z;			/* slice number to work on */
  TaskType type;                /* task TASK_ESTIMATE or TASK_CORRECT */
  float mean_real;              /* real adjustment for slice */
  float mean_imag;              /* imaginary adjustment for slice */
} Task;

typedef struct Result {
  /* NOTE: any new fields added to this struct should
     be also added to PackResult and UnpackResult */
  int t;
  int z;
  float mean_real;
  float mean_imag;
} Result;

static char rcsid[] = "$Id: baseline.c,v 1.17 2007/04/19 22:34:18 welling Exp $";

/* GLOBAL VARIABLES FOR MASTER & WORKER */
static char* progname;
Task t;
Context c;
Result r;

/* GLOBAL VARIABLES FOR MASTER */
static MRI_Dataset *mInput;
static MRI_Dataset *mOutput;
static FILE *pOutput;
static Filename output_data_file; /* the file where the
				     images chunk will
				     reside */
static Filename parameter_file;	/* the file where means
				   will be written */
static Filename smooth_par_file;/* the file where smoothed means
				   will be written */
static char reverse_string[128];
static int reverse_set= 0;
static FComplex **mean = NULL;		/* means of images */
static FComplex **smoothMean = NULL;    /* means after smoothing */

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
void WorkerFinalize ();
void PackContext ();
void UnpackContext ();
void PackTask ();
void UnpackTask ();
void PackResult ();
void UnpackResult ();
static int parse_cl_reverse( const char* s, int* revOdd, int* revEven );
static void parse_mrihdr_reverse( MRI_Dataset* ds, const char* chunk, 
				  int* revOdd, int* revEven );


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
  FILE *fp = NULL;
  float bandwidth = 5.0;
  float kpar = 60.0;
  Smoother* smoother;
  float *scratchReal;
  float *scratchImag;
  float *smoothScratchReal;
  float *smoothScratchImag;
  float *dtblIn[2];
  float *dtblOut[2];
  unsigned char** missing= NULL;
  char* dimstr;

  verbose = TRUE;
  /* Print version number */
  Report("# %s\n",rcsid);

  /* Set up the smoother */
  sm_init();
  sm_set_params( SM_GAUSSIAN, bandwidth, kpar, 0.0, NULL );

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  ReadArguments(argc, argv);

  /* Open input dataset */
  if( !strcmp( c.input_file, c.output_file ) )
    Abort( "Input and output files must be distinct." );
  mInput = mri_open_dataset( c.input_file, MRI_READ );

  /* Check that program will function on data-set */
  if( !mri_has( mInput, "images" ) )
    Abort( "%s operates only on standard images.", argv[0] );
  if( !mri_has( mInput, "images.dimensions" ) )
    Abort("%s: images.dimensions tag is missing!\n",argv[0]);
  dimstr= mri_get_string(mInput, "images.dimensions");
  if (strcmp(dimstr,"vxyzt") && strcmp(dimstr,"vqyzt"))
    Abort( "%s operates only on images in order vxyzt or vqyzt.", argv[0] );
  if( !mri_has( mInput, "images.extent.v" ) ||
      ( ( dv = mri_get_int( mInput, "images.extent.v" ) ) != 2 ) )
    Abort( "%s operates only on complex-valued images (v = 2).", argv[0] );

  /* If scan line reversal information is given in the command line,
   * implement it.  Otherwise, infer from the input file.
   */
  if (reverse_set) {
    if (!parse_cl_reverse( reverse_string, 
			   &(c.reverse_odd), &(c.reverse_even) ))
      Abort("%s: invalid -reverse opion <%s> given\n", 
	    progname,reverse_string);
  }
  else parse_mrihdr_reverse(mInput, "images", 
			    &(c.reverse_odd), &(c.reverse_even));


  /* Read/Create missing image indicators */
  missing = get_missing( mInput );

  /* Set output dataset */
  mOutput = mri_copy_dataset( c.output_file, mInput );
  hist_add_cl( mOutput, argc, argv );
  mri_set_string( mOutput, "images.datatype", "float32" );
  mri_set_int( mOutput, "images.rowflip", 0 );
  mri_close_dataset ( mOutput);

  /* Set parameters in local variables */
  if( !mri_has( mInput, "images.extent.t" ) ||
      (!mri_has( mInput, "images.extent.x" ) 
	&& !mri_has(mInput, "images.extent.q")) ||
      !mri_has( mInput, "images.extent.y" ) ||
      !mri_has( mInput, "images.extent.z" ) )
    Abort( "images.extent key(s) missing from header." );
  c.dt = mri_get_int( mInput, "images.extent.t" );
  if (mri_has(mInput,"images.extent.x"))
    c.dx = mri_get_int( mInput, "images.extent.x" );
  else c.dx= mri_get_int( mInput, "images.extent.q" );
  c.dy = mri_get_int( mInput, "images.extent.y" );
  c.dz = mri_get_int( mInput, "images.extent.z" );
  mri_close_dataset ( mInput );
  if( ( c.dt <= 0 ) || ( c.dx <= 0 ) || ( c.dy <= 0 ) || ( c.dz <= 0 ) )
    Abort( "images.extent key(s) is non-positive." );

  /* Allocate parameter storage */
  mean = Matrix( c.dt, c.dz, FComplex );
  smoothMean = Matrix( c.dt, c.dz, FComplex );
  if (!(scratchReal=(float*)malloc(c.dt*sizeof(float))))
    Abort("unable to allocate %d bytes!\n",c.dt*sizeof(float));
  if (!(scratchImag=(float*)malloc(c.dt*sizeof(float))))
    Abort("unable to allocate %d bytes!\n",c.dt*sizeof(float));
  if (!(smoothScratchReal=(float*)malloc(c.dt*sizeof(float))))
    Abort("unable to allocate %d bytes!\n",c.dt*sizeof(float));
  if (!(smoothScratchImag=(float*)malloc(c.dt*sizeof(float))))
    Abort("unable to allocate %d bytes!\n",c.dt*sizeof(float));
  dtblIn[0]= scratchReal;
  dtblIn[1]= scratchImag;
  dtblOut[0]= smoothScratchReal;
  dtblOut[1]= smoothScratchImag;
  
  /* Calculate area of region over which data is summed */
  c.baseline_area = 
    c.dy * ( (long) ( c.dx / 4 ) + ( c.dx - (long) ( 3 * c.dx / 4 ) ) );

  /* set the context for the workers to calculate means */
  par_set_context();

  /* Loop through images, finding means */
  Report("Estimating baseline corrections\n");
  t.type= TASK_ESTIMATE;
  t.mean_real= 0.0;
  t.mean_imag= 0.0;
  for ( t.t = 0; t.t < c.dt; t.t++ )
    for ( t.z = 0; t.z < c.dz; t.z++ ) 
      par_delegate_task();
  while (par_tasks_outstanding() > 0)
    par_wait(1.0);

  /* Smooth */
  Report("Smoothing\n");
  smoother= sm_create_smoother();

  for (t.z=0; t.z<c.dz; t.z++) {
    for (t.t=0; t.t<c.dt; t.t++) {
      scratchReal[t.t]= mean[t.t][t.z].real;
      scratchImag[t.t]= mean[t.t][t.z].imag;
    }
    SM_SMOOTH_GROUP(smoother, dtblIn, dtblOut, 2, c.dt, missing, t.z);
    for (t.t=0; t.t<c.dt; t.t++) {
      smoothMean[t.t][t.z].real= smoothScratchReal[t.t];
      smoothMean[t.t][t.z].imag= smoothScratchImag[t.t];
    }
  }

  sm_destroy(smoother);
      
  /* Write parameter files */
  Report("Writing smoothed and unsmoothed estimates\n");
  fp = efopen( parameter_file, "w" );
  for ( t.t = 0; t.t < c.dt; t.t++ )
    for ( t.z = 0; t.z < c.dz; t.z++ )
      fprintf( fp, "%15.6f %15.6f\n", 
	       mean[t.t][t.z].real, mean[t.t][t.z].imag );
  efclose( fp );

  fp = efopen( smooth_par_file, "w" );
  for ( t.t = 0; t.t < c.dt; t.t++ )
    for ( t.z = 0; t.z < c.dz; t.z++ )
      fprintf( fp, "%15.6f %15.6f\n", 
	       smoothMean[t.t][t.z].real, smoothMean[t.t][t.z].imag );
  efclose( fp );

  /* Loop through images, applying corrections */
  Report("Implementing estimated baseline corrections\n");
  t.type= TASK_CORRECT;
  for ( t.t = 0; t.t < c.dt; t.t++ )
    for ( t.z = 0; t.z < c.dz; t.z++ ) {
      t.mean_real= smoothMean[t.t][t.z].real;
      t.mean_imag= smoothMean[t.t][t.z].imag;
      par_delegate_task();
    }
  par_finish();

  FreeMatrix(mean);
  FreeMatrix(smoothMean);
  free(scratchReal);
  free(scratchImag);
  free(smoothScratchReal);
  free(smoothScratchImag);
  
  Report( "\n#      Baseline adjustment complete.\n" );
  exit(0);

}

void
ReadArguments (int argc,
	       char **argv)
{
  cl_scan( argc, argv );

  /* Get filenames */

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
 
  
  cl_get( "estimates|est|e", "%option %s[%]", "rawbaseline.par", parameter_file );
  cl_get( "smoothedestimates|sme", "%option %s[%]", "baseline.par", 
	  smooth_par_file );
  reverse_set= cl_get( "reverse|rev|r", "%option %s", reverse_string );

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
}

void MasterResult (int task_number)
{
  static int results = 0;
  int n;

  mean[r.t][r.z].real = r.mean_real;
  mean[r.t][r.z].imag = r.mean_imag;

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
  static int old_dx = 0, old_dy = 0;
  int x, y;

  if (first)
    {
      wInput = mri_open_dataset(c.input_file, MRI_READ);
      wOutput = mri_open_dataset(c.output_file, MRI_MODIFY_DATA);
      first = FALSE;
    }

  if (c.dx != old_dx || c.dy != old_dy)
    {
      if (image != NULL)
	free(image);
      if (corr_image != NULL)
	FreeMatrix(corr_image);

      /* Allocate image storage */
      image = (FComplex **) emalloc( c.dy * sizeof(FComplex *) );
      corr_image = Matrix( c.dy, c.dx, FComplex );

      old_dx = c.dx;
      old_dy = c.dy;
    }
}

void
WorkerTask ()
{
  long x, y;
  double sum_real, sum_imag;
  long imageSize= 2*c.dx*c.dy;

  /* Get image and re-align matrix pointers */
  *image = (FComplex *) 
    mri_get_chunk( wInput, "images", imageSize,
		   (t.t*c.dz + t.z)*imageSize, MRI_FLOAT );
  realign_matrix( (void **) image, c.dy, 
		  (long) ( c.dx * sizeof(FComplex) ) );

  switch (t.type) {
  case TASK_ESTIMATE:
    {
      /* Estimate baseline --- mean of outer quarter-images */
      /*   in x-direction                                   */
      sum_real = sum_imag = 0.0;

      for( y = 0; y < c.dy; y++ )
	{
	  /* Sum over left quarter of image */
	  for( x = 0; x < (long) ( c.dx / 4 ); x++ )
	    {
	      sum_real += (double) image[y][x].real;
	      sum_imag += (double) image[y][x].imag;
	    }
	  
	  /* Sum over right quarter of image */
	  for( x = (long) ( 3 * c.dx / 4 ); x < c.dx; x++ )
	    {
	      sum_real += (double) image[y][x].real;
	      sum_imag += (double) image[y][x].imag;
	    }
	}

      /* Convert sums to means */
      r.t = t.t;
      r.z = t.z;
      r.mean_real = (float) ( sum_real / (double) c.baseline_area );
      r.mean_imag = (float) ( sum_imag / (double) c.baseline_area );
    }
    break;

  case TASK_CORRECT:
    {
      /* Subtract off estimated baseline and   */
      /*   perform line reversals as necessary */

	for( y = 0; y < c.dy; y++ )
	  {
	    
	    if( y % 2 )
	      
	      if ( c.reverse_odd )
		for( x = 0; x < c.dx; x++ )
		  {
		    corr_image[y][x].real = image[y][c.dx-x-1].real - 
		      t.mean_real;
		    corr_image[y][x].imag = image[y][c.dx-x-1].imag - 
		      t.mean_imag;
		  }
	      else 
		for( x = 0; x < c.dx; x++ ) 
		  {
		    corr_image[y][x].real = image[y][x].real - 
		      t.mean_real;
		    corr_image[y][x].imag = image[y][x].imag - 
		      t.mean_imag;
		  }
	    
	    else
	      
	      if ( c.reverse_even )
		for( x = 0; x < c.dx; x++ ) 
		  {
		    corr_image[y][x].real = image[y][c.dx-x-1].real - 
		      t.mean_real;
		    corr_image[y][x].imag = image[y][c.dx-x-1].imag - 
		      t.mean_imag;
		  }
	      else
		for( x = 0; x < c.dx; x++ )
		  {
		    corr_image[y][x].real = image[y][x].real - 
		      t.mean_real;
		    corr_image[y][x].imag = image[y][x].imag - 
		      t.mean_imag;
		  }
	  }
	
	/* Set output image */
	mri_set_chunk( wOutput, "images", imageSize, 
		       (t.t*c.dz + t.z)*imageSize, MRI_FLOAT, *corr_image );
    }
    break;
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
  par_pkint(c.dz);
  par_pkint(c.dt);
  par_pkint(c.reverse_odd);
  par_pkint(c.reverse_even);
  par_pklong(c.baseline_area);
}

void
UnpackContext ()
{
  par_upkstr(c.input_file);
  par_upkstr(c.output_file);
  c.dx= par_upkint();
  c.dy= par_upkint();
  c.dz= par_upkint();
  c.dt= par_upkint();
  c.reverse_odd= par_upkint();
  c.reverse_even= par_upkint();
  c.baseline_area= par_upklong();
}

void
PackTask ()
{
  int task_type;
  task_type= (t.type == TASK_CORRECT) ? 1 : 0;
  
  par_pkint(t.t);
  par_pkint(t.z);
  par_pkint(task_type);
  par_pkfloat(t.mean_real);
  par_pkfloat(t.mean_imag);
}

void
UnpackTask ()
{
  int task_type;
  
  t.t= par_upkint();
  t.z= par_upkint();
  task_type= par_upkint();
  t.mean_real= par_upkfloat();
  t.mean_imag= par_upkfloat();
  
  if (task_type==1) t.type= TASK_CORRECT;
  else t.type= TASK_ESTIMATE;
}

void
PackResult ()
{
  par_pkint(r.t);
  par_pkint(r.z);
  par_pkfloat(r.mean_real);
  par_pkfloat(r.mean_imag);
}

void
UnpackResult ()
{
  r.t= par_upkint();
  r.z= par_upkint();
  r.mean_real= par_upkfloat();
  r.mean_imag= par_upkfloat();
}


