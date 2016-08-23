/************************************************************
 *                                                          *
 *  nav_deghost.c                                           *
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <unistd.h>
#ifdef DARWIN
#include <sys/time.h>
#endif
#include <time.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"
#include "par.h"
#include "misc.h"
#include "acct.h"
#include "../fmri/lapack.h"

static char rcsid[] = "$Id: nav_deghost.c,v 1.11 2007/03/21 23:46:16 welling Exp $";

static void phase_est( double* phaseShift, double* absolutePhase );
static double loc_est();
static void opp_phase_shift( FComplex** image, long dy, long dx, double phase, 
			     FComplex** corr_image );
static void opp_loc_shift( FComplex** image, long dy, long dx, double shift, 
			   FComplex** corr_image );
static void flip_rows(FComplex** image, const int flipOdd, const int flipEven, 
		      const int yLow, const int yHigh);
static int parse_cl_reverse( const char* s, int* revOdd, int* revEven );
static void parse_mrihdr_reverse( MRI_Dataset* ds, const char* chunk, 
				  int* revOdd, int* revEven );
static void check_input_format();

#ifdef never
static double eval_phase(double phase);
static void restrt();
#endif

typedef enum { MODE_ESTIMATE, MODE_APPLY, MODE_BOTH } Mode;

typedef enum { ALG_WT_SUM } Algorithm;

typedef enum { NAV_ESTIMATE, NAV_APPLY } TaskType;

typedef struct Context {
  /* NOTE: any new fields added to this struct should
     be also added to PackContext and UnpackContext */ 
  Filename input_file;		/* The images to work on */
  Filename output_file;         /* The output dataset */
  long dv;			/* # of values along the V dimension */
  long dx;			/* # of pixels along the X dimension */
  long dy;			/* # of pixels along the Y dimension */
  long dz;			/* # of slices along the Z dimension */
  long dt;			/* # of images along the T dimension */
  long dn;                      /* # of navigator lines */
  long yCtr;                    /* y of k-space center */
  int use_yCtr;                 /* Include yCtr line in calculation */
  int flip_odd;                 /* # whether to flip image odd rows */
  int flip_even;                /* # whether to flip image even rows */
  int nav_flip_odd;             /* # whether to flip navigator odd rows */
  int nav_flip_even;            /* # whether to flip navigator even rows */
  Algorithm alg;                /* algorithm to use in making estimates */ 
} Context;

typedef struct Task {
  /* NOTE: any new fields added to this struct should
     be also added to PackTask and UnpackTask */ 
  TaskType type;        /* type of Task to perform */
  long z;		/* slice number to work on */
  long t;		/* image number to work on */
  double phaseShift;    /* used only for TaskType NAV_APPLY */
  double locShift;      /* used only for TaskType NAV_APPLY */
  double absolutePhase; /* used only for TaskType NAV_APPLY */
} Task;

typedef struct Result {
  /* NOTE: any new fields added to this struct should
     be also added to PackResult and UnpackResult */
  TaskType type;
  long z; 
  long t;
  double phaseShift;    /* used only for TaskType NAV_ESTIMATE */
  double locShift;      /* used only for TaskType NAV_ESTIMATE */
  double absolutePhase; /* used only for TaskType NAV_ESTIMATE */
} Result;

/* GLOBAL VARIABLES FOR MASTER & WORKER */
static char* progname;
static int debug= 0;
Task t;
Context c;
Result r;

/* GLOBAL VARIABLES FOR MASTER */
static MRI_Dataset *mInput= NULL;
static MRI_Dataset *mOutput= NULL;
static float **phaseShift;
static float **locShift;
static float **absolutePhase;
static float **smPhaseShift;
static float **smLocShift;
static float **smAbsolutePhase;
static Mode mode;

/* GLOBAL VARIABLES FOR WORKER */
static MRI_Dataset *wInput= NULL;
static MRI_Dataset *wOutput= NULL;
static FComplex **tmp_nav = NULL;
static FComplex **tmp2_nav = NULL;
static FComplex **nav = NULL;
static FComplex **tmp_image = NULL;
static FComplex **tmp2_image = NULL;
static FComplex **image = NULL;

/* FORWARD DECLARATIONS */
static void MasterTask (int argc, char **argv, char **envp);
static void ReadArguments (int argc, char **argv);
static void MasterResult (int task_number);
static void WorkerContext ();
static void WorkerTask ();
static void WorkerFinalize ();
static void PackContext ();
static void UnpackContext ();
static void PackTask ();
static void UnpackTask ();
static void PackResult ();
static void UnpackResult ();

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

static void check_input_format()
{
  char* dimstr;
  char* nav_dimstr;

  /* 
   * Check that program will function on data-set.  
   * Image chunk first. 
   */
  if( !mri_has( mInput, "images" ) )
    Abort( "%s operates only on standard images.\n", progname );
  if ( !mri_has( mInput, "images.dimensions" ) )
    Abort( "%s has no images.dimensions tag!\n", progname );
  dimstr= mri_get_string( mInput, "images.dimensions" );
  if (dimstr[0] != 'v' || !mri_has( mInput, "images.extent.v" )
      || (mri_get_int( mInput, "images.extent.v" ) != 2 ) )
    Abort( "%s operates only on complex-valued images (v = 2).\n", progname );
  c.dv= 2;
  if ( strcmp(dimstr,"vxyzt") && strcmp(dimstr,"vqyzt")
       && strcmp(dimstr,"vxyz") && strcmp(dimstr,"vqyz")
       && strcmp(dimstr,"vxy") && strcmp(dimstr,"vqy") )
    Abort("%s operates only on vxy(z)(t) or vqy(z)(t) images!\n",progname);
  if (strchr(dimstr,'x')) {
    if (mri_has(mInput,"images.extent.x")) 
      c.dx= mri_get_int(mInput,"images.extent.x");
    else Abort("%s: input is missing images.extent.x tag!\n",progname);
  }
  else {
    if (mri_has(mInput,"images.extent.q")) 
      c.dx= mri_get_int(mInput,"images.extent.q");
    else Abort("%s: input is missing images.extent.q tag!\n",progname);
  }
  if (mri_has(mInput,"images.extent.y")) 
    c.dy= mri_get_int(mInput,"images.extent.y");
  else Abort("%s: input is missing images.extent.y tag!\n",progname);
  if (mri_has(mInput,"images.extent.z")) 
    c.dz= mri_get_int(mInput,"images.extent.z");
  else c.dz= 1;
  if (mri_has(mInput,"images.extent.t")) 
    c.dt= mri_get_int(mInput,"images.extent.t");
  else c.dt= 1;
  if( ( c.dt <= 0 ) || ( c.dx <= 0 ) || ( c.dy <= 0 ) || ( c.dz <= 0 ) )
    Abort( "images.extent key(s) is non-positive." );

  /* 
   * Check the navigator chunk dimensions
   */
  if( !mri_has( mInput, "navigator" ) )
    Abort( "%s requires input with an associated navigator chunk.\n", 
	   progname );
  if ( !mri_has( mInput, "navigator.dimensions" ) )
    Abort( "%s has no navigator.dimensions tag!\n", progname );
  nav_dimstr= mri_get_string( mInput, "navigator.dimensions" );
  if (nav_dimstr[0] != 'v' || !mri_has( mInput, "navigator.extent.v" )
      || (mri_get_int( mInput, "navigator.extent.v" ) != 2 ) )
    Abort( "%s operates only on complex-valued navigator (v = 2).\n", 
	   progname );
  if ( strcmp(nav_dimstr,"vxnzt") && strcmp(nav_dimstr,"vqnzt")
       && strcmp(nav_dimstr,"vxnz") && strcmp(nav_dimstr,"vqnz")
       && strcmp(nav_dimstr,"vxn") && strcmp(nav_dimstr,"vqn") )
    Abort("%s operates only on vxn(z)(t) or vqn(z)(t) navigator!\n",progname);
  if (strchr(nav_dimstr,'x')) {
    if (mri_has(mInput,"navigator.extent.x")) {
      if (c.dx != mri_get_int(mInput,"navigator.extent.x"))
	Abort("%s: navigator and image widths incommensurate!\n",progname);
    }
    else Abort("%s: input is missing navigator.extent.x tag!\n",progname);
  }
  else {
    if (mri_has(mInput,"navigator.extent.q")) {
      if (c.dx != mri_get_int(mInput,"navigator.extent.q"))
	Abort("%s: navigator and image widths incommensurate!\n",progname);
    }
    else Abort("%s: input is missing navigator.extent.q tag!\n",progname);
  }
  if (mri_has(mInput,"navigator.extent.n")) 
    c.dn= mri_get_int(mInput,"navigator.extent.n");
  else Abort("%s: input is missing navigator.extent.n tag!\n",progname);
  if (mri_has(mInput,"navigator.extent.z")) {
    if (c.dz != mri_get_int(mInput,"navigator.extent.z"))
      Abort("%s: navigator and image slice counts incommensurate!\n",progname);
  }
  else if (c.dz != 1)
    Abort("%s: input is missing navigator.extent.z tag!\n",progname);
  if (mri_has(mInput,"navigator.extent.t")) {
    if (c.dt != mri_get_int(mInput,"navigator.extent.t"))
      Abort("%s: navigator and image time counts incommensurate!\n",progname);
  }
  else if (c.dt != 1)
    Abort("%s: input is missing navigator.extent.t tag!\n",progname);
}

static FILE* create_parfile( const char* fname, const char* type, 
			     const char* infname )
{
  time_t tm;
  FILE* result= fopen( fname, "w" );
  tm= time(NULL);
  fprintf(result,"##Format: order:index_zt, type:%s\n",type);
  fprintf(result,
     "##Format: names:(NavEchoPhaseShift,NavEchoLocShift,NavEchoAbsPhase)\n");
  fprintf(result,"# Navigator deghost parameters generated %s",
	  asctime(localtime(&tm)));
  fprintf(result,"# Input file: %s\n",infname);
  fflush(result);
  return result;
}

static void read_parfile( const char* fname, float** phases, 
			  float** locs, float** absPhases )
{
  FILE* f;
  long nLoaded= 0;
  long t;
  long z;
  double ph;
  double loc;
  double absPh;
  char line[256];
  long lineCnt= 0;

  if (!(f= fopen(fname,"r")))
    Abort("%s: unable to open <%s> for reading!\n",
	  progname, fname);

  while (!feof(f) && !ferror(f)) {
    if (!fgets(line, sizeof(line), f)) break;
    lineCnt++;
    if (line[0]=='#' || strlen(line)==0) continue;
    if (sscanf(line,"%ld %ld %lg %lg %lg", &z, &t, &ph, &loc, &absPh) != 5)
      Abort("%s: error parsing parameter file <%s> at line %ld!\n",
	    progname, fname, lineCnt);
    if (t<0 || t>=c.dt)
      Abort("%s: t == %ld is out of range in <%s> line %ld!\n",
	    progname, t, fname, lineCnt); 
    if (z<0 || z>=c.dz)
      Abort("%s: z == %ld is out of range in <%s> line %ld!\n",
	    progname, z, fname, lineCnt); 
    phases[z][t]= ph;
    locs[z][t]= loc;
    absPhases[z][t]= absPh;
    nLoaded++;
  }

  (void)fclose(f);

  if (nLoaded != c.dt*c.dz)
    Abort("%s: found only %ld lines in <%s>; expected %ld!\n",
	  progname, nLoaded, fname, c.dt*c.dz);
}

static void
MasterTask (int argc,
	    char **argv,
	    char **envp)
{
  char parfile[512], smparfile[512];
  long n;
  long i;
  float bandwidth = 45.0;
  float kpar = 60.0;
  char reverse_string[256];
  int reverse_set= 0;
  char nav_reverse_string[256];
  int nav_reverse_set= 0;
  int include_flag= 0;
  char mode_string[256];
  int mode_set= 0;
  Smoother* smoother;
  unsigned char** missing= NULL;

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

  sm_init();
  sm_set_params( SM_GAUSSIAN, bandwidth, kpar, 0.0, NULL );

  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Get filenames */
   if (cl_present( "dataout|d" ))
     Abort ("Option dataout|d no longer works.  Please see help file.\n");
   if (cl_present( "headerout|h" ))
     Abort ("Option headerout|h no longer works.  Please see help file.\n");
   if (cl_present( "input|i" ))
     Abort ("Option input|i has been replaced by infile outfile format.  Please see help file.\n");
   if (cl_present( "parameters|p" ))       
     Abort ("Option parameters|p has been replaced by estimates|est|e.  Please see help file.\n");
    if (cl_present( "smoothedparameters|s" ))       
     Abort ("Option smoothedparameters|s has been replaced by smoothedestimates|sme.  Please see help file.\n");

  cl_get( "estimates|est|e", "%option %s[%]", "navdeghost.par", parfile );
  cl_get( "smoothedestimates|sme", "%option %s[%]", "navdeghost.smoothed.par", 
	  smparfile ); 
  reverse_set= cl_get( "reverse|rev|r", "%option %s", reverse_string );
  nav_reverse_set= cl_get( "nav_reverse|nrv", "%option %s", 
			   nav_reverse_string );
  mode_set= cl_get( "mode", "%option %s", mode_string );
  verbose= cl_present("verbose|ver|v");
  debug= cl_present("debug");
  include_flag= cl_present("include|inc");

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

  if( !strcmp(c.input_file, c.output_file ) )
    Abort( "Input and output files must be distinct." );

  if( !strcmp(parfile,smparfile) )
    Abort("Raw and smoothed parameter files must be distinct." );

  /* Open input dataset */
  mInput = mri_open_dataset(c.input_file, MRI_READ );

  /* Test that the program will function on the dataset */
  check_input_format();

  /* If scan line reversal information is given in the command line,
   * implement it.  Otherwise, infer from the input file.
   */
  if (reverse_set) {
    if (!parse_cl_reverse( reverse_string, &(c.flip_odd), &(c.flip_even) ))
      Abort("%s: invalid -reverse opion <%s> given\n", 
	    progname,reverse_string);
  }
  else parse_mrihdr_reverse(mInput, "images", &(c.flip_odd), &(c.flip_even));
  if (nav_reverse_set) {
    if (!parse_cl_reverse( nav_reverse_string, 
			   &(c.nav_flip_odd), &(c.nav_flip_even) ))
      Abort("%s: invalid -navreverse option <%s> given\n", 
	    progname,nav_reverse_string);
  }
  else parse_mrihdr_reverse(mInput, "navigator", 
			    &(c.nav_flip_odd), &(c.nav_flip_even));

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

  c.alg= ALG_WT_SUM;

  if (include_flag) {
    /*
     * If the center line of k-space is known, we can include it in 
     * the deghosting process.  It's equivalent to another navigator.
     */
    if (mri_has(mInput,"images.kspace_ctr.y")) {
      long yCtr= mri_get_int(mInput,"images.kspace_ctr.y");
      if (yCtr>=0 && yCtr<=c.dy) {
	c.use_yCtr= 1;
	c.yCtr= yCtr;
      }
      else {
	c.use_yCtr= 0;
	c.yCtr= 0;
      }
    }
    else {
      c.use_yCtr= 0;
      c.yCtr= 0;
    }
  }
  else {
    c.use_yCtr= 0;
    c.yCtr= 0;
  }

  /* Read/Create missing image indicators */
  missing = get_missing( mInput );

  /* Set output dataset */
  if (mode != MODE_ESTIMATE) {
    mOutput = mri_copy_dataset( c.output_file, mInput );
    hist_add_cl( mOutput, argc, argv );
    mri_set_int( mOutput, "images.rowflip", 0 );
    mri_set_int( mOutput, "navigator.rowflip", 0 );
    mri_set_string( mOutput, "images.datatype", "float32" );
    mri_set_string( mOutput, "images.file", ".dat" );
  }

  /* close datasets */
  mri_close_dataset( mInput );
  if (mOutput != NULL) mri_close_dataset( mOutput );

  /* Allocate parameter storage */
  phaseShift = Matrix(c.dz, c.dt, float);
  locShift = Matrix(c.dz, c.dt, float);
  absolutePhase = Matrix(c.dz, c.dt, float);
  smPhaseShift = Matrix(c.dz, c.dt, float);
  smLocShift = Matrix(c.dz, c.dt, float);
  smAbsolutePhase = Matrix(c.dz, c.dt, float);

  /* Set context for workers */
  par_set_context();

  if (mode==MODE_APPLY) {
    /* 
     * Read estimates from smoothed parameter file 
     */
    read_parfile( smparfile, smPhaseShift, smLocShift, smAbsolutePhase );
  }
  else {
    FILE *par = NULL;
    FILE *smpar = NULL;

    /* Open parameter files */
    par= create_parfile( parfile, "raw", c.input_file );
    smpar= create_parfile( smparfile, "filtered", c.input_file );

    /* 
     * Estimate phase shift 
     */
    t.type= NAV_ESTIMATE;
    for (t.t = 0; t.t < c.dt; t.t++)
      for (t.z = 0; t.z < c.dz; t.z++)
	par_delegate_task();
    while (par_tasks_outstanding() > 0)
      par_wait(1.0);

    /*
     * Write unsmoothed estimates
     */
    for (t.z = 0; t.z < c.dz; ++t.z)
      for (t.t = 0; t.t < c.dt; ++t.t)
	fprintf(par, "%ld %ld %lg %lg %lg\n", 
		t.z, t.t, phaseShift[t.z][t.t], locShift[t.z][t.t],
		absolutePhase[t.z][t.t]);
    efclose( par );
    Report("Finished writing phase estimates to %s\n", parfile);

    /*
     * Generate and write the smoothed estimates
     */
    smoother= sm_create_smoother();
    for (t.z = 0; t.z < c.dz; ++t.z) {
	float* dtbl_in[3];
	float* dtbl_out[3];
	dtbl_in[0]= phaseShift[t.z];
	dtbl_in[1]= locShift[t.z];
	dtbl_in[2]= absolutePhase[t.z];
	dtbl_out[0]= smPhaseShift[t.z];
	dtbl_out[1]= smLocShift[t.z];
	dtbl_out[2]= smAbsolutePhase[t.z];
	SM_SMOOTH_GROUP( smoother, dtbl_in, dtbl_out, 3, c.dt, missing, t.z );
	for (t.t = 0; t.t < c.dt; ++t.t)
	  fprintf(smpar, "%ld %ld %lg %lg %lg\n", 
		  t.z, t.t, smPhaseShift[t.z][t.t], smLocShift[t.z][t.t],
		  smAbsolutePhase[t.z][t.t]);
    }
    efclose( smpar );
    sm_destroy(smoother);
    Report("Finished writing smoothed phase estimates to %s\n", smparfile);
  }

  if (mode != MODE_ESTIMATE) {
    /*
     * Apply the smoothed estimates
     */
    t.type= NAV_APPLY;
    for (t.t = 0; t.t < c.dt; t.t++)
      for (t.z = 0; t.z < c.dz; t.z++) {
	t.phaseShift= smPhaseShift[t.z][t.t];
	t.locShift= smLocShift[t.z][t.t];
	t.absolutePhase= smAbsolutePhase[t.z][t.t];
	par_delegate_task();
      }
    while (par_tasks_outstanding() > 0)
      par_wait(1.0);
  }

  par_finish();

  Message( "#      Navigator de-ghosting estimation complete.\n" );
  exit(0);
}

static void MasterResult (int task_number)
{
  static long results = 0;
  long n;

  switch (r.type) {
  case NAV_ESTIMATE:
    {
      phaseShift[r.z][r.t] = (float)r.phaseShift;
      locShift[r.z][r.t] = (float)r.locShift;
      absolutePhase[r.z][r.t] = (float)r.absolutePhase;
    }
    break;
  case NAV_APPLY:
    {
      /* Do nothing */
    }
    break;
  default: Abort("%s: internal error: unrecognized result type %d\n",
		 progname,(int)r.type);
  }

  if (++results % c.dz == 0) {
    /* Print progress report */
    n = (results-1) / c.dz;
    if (n == 0)
      Report( "      " );
    
    if (n == (c.dt - 1) || (n+1) % 60 == 0) {
      Report( "# %ld\n", (long) ( n + 1 ) );
      if (n != (c.dt - 1))
	Report( "      " );
      else
	results = 0;
    }
    else {
      Report("#");
    }
  }
}

/* WORKER PROCEDURES */

void WorkerFinalize()
{
  if (wInput != NULL) mri_close_dataset(wInput);
  if (wOutput != NULL) mri_close_dataset(wOutput);
  if (image != NULL) free(image);
  if (tmp_image != NULL) FreeMatrix(tmp_image);
  if (tmp2_image != NULL) FreeMatrix(tmp2_image);
  if (nav != NULL) FreeMatrix(nav);
  if (tmp_nav != NULL) FreeMatrix(tmp_nav);
  if (tmp2_nav != NULL) FreeMatrix(tmp2_nav);
}

static void
WorkerContext ()
{
  static long old_dx = 0, old_dy = 0, old_dnEffective= 0;
  long dnEffective = ( c.use_yCtr ? c.dn+1 : c.dn );

  if (c.dx != old_dx || c.dy != old_dy || dnEffective != old_dnEffective)
    {
      /* free up old storage */
      if (image != NULL) free(image);
      if (tmp_image != NULL) FreeMatrix(tmp_image);
      if (tmp2_image != NULL) FreeMatrix(tmp2_image);
      if (nav != NULL) FreeMatrix(nav);
      if (tmp_nav != NULL) FreeMatrix(tmp_nav);
      if (tmp2_nav != NULL) FreeMatrix(tmp2_nav);

      /* Allocate parameter and image storage */
      image = (FComplex**)emalloc( c.dy*sizeof(FComplex*) );
      tmp_image = Matrix( c.dy, c.dx, FComplex );
      tmp2_image = Matrix( c.dy, c.dx, FComplex );
      nav = Matrix( dnEffective, c.dx, FComplex );
      tmp_nav = Matrix( dnEffective, c.dx, FComplex );
      tmp2_nav = Matrix( dnEffective, c.dx, FComplex );

      old_dx = c.dx;
      old_dy = c.dy;
      old_dnEffective = dnEffective;
    }
}

static void
WorkerTask ()
{
  long imageSize= c.dv*c.dx*c.dy;
  long navImageSize= c.dv*c.dx*c.dn;
  long dnEffective = ( c.use_yCtr ? c.dn+1 : c.dn );
  long i;
  long j;

  r.t= t.t;
  r.z= t.z;
  r.type= t.type;

  switch (t.type) {
  case NAV_ESTIMATE: 
    {
      FComplex* navData;
      if (!wInput) wInput= mri_open_dataset(c.input_file, MRI_READ);
      navData= (FComplex*)mri_get_chunk(wInput,"navigator",
					navImageSize,
					(t.t*c.dz+t.z)*navImageSize,
					MRI_FLOAT);
      for (j=0; j<c.dn; j++) 
	for (i=0; i<c.dx; i++)
	  nav[j][i]= navData[j*c.dx + i];
      if (c.nav_flip_odd || c.nav_flip_even) 
	flip_rows( nav, c.nav_flip_odd, c.nav_flip_even, 0, c.dn-1 );
      
      if (c.use_yCtr) {
	FComplex* yCtrData= (FComplex*)mri_get_chunk(wInput,"images",
						     c.dv*c.dx,
						     (t.t*c.dz+t.z)*imageSize,
						     MRI_FLOAT);
	for (i=0; i<c.dx; i++) nav[c.dn][i]= yCtrData[i];
	if (c.flip_odd || c.flip_even) 
	  flip_rows( nav, c.flip_odd, c.flip_even, c.dn, c.dn );
      }
      
#ifdef never
      {
	int x,y;
	for (y=0; y<dnEffective; y++) 
	  for (x=0; x<c.dx; x++)
	    fprintf(stderr,"%d %d %g %g\n",
		    x,y,nav[y][x].real,nav[y][x].imag);
	exit(0);
      }
#endif
      
      phase_est(&r.phaseShift, &r.absolutePhase);
      r.locShift= loc_est();
    }
    break;
  case NAV_APPLY: 
    {
      FComplex* navData;
      double navPhaseShift;
      double navLocShift;
      double imagePhaseShift;
      double imageLocShift;

      navPhaseShift= t.phaseShift;
      navLocShift= t.locShift;
      if (c.dn % 2) {
	/* Odd number of navigator scan lines implies a sign flip
	 * between the navigator and the image.
	 */
	imagePhaseShift= -navPhaseShift;
	imageLocShift= -navLocShift;
      }
      else {
	/* Even number of navigator scan lines implies no sign flip
	 * between the navigator and the image.
	 */
	imagePhaseShift= navPhaseShift;
	imageLocShift= navLocShift;
      }

      if (!wInput) wInput= mri_open_dataset(c.input_file, MRI_READ);
      if (!wOutput) wOutput= mri_open_dataset(c.output_file, MRI_MODIFY_DATA);

      /* Apply the changes to the navigator... */
      navData= (FComplex*)mri_get_chunk(wInput,"navigator",
					navImageSize,
					(t.t*c.dz+t.z)*navImageSize,
					MRI_FLOAT);
      for (j=0; j<c.dn; j++) 
	for (i=0; i<c.dx; i++)
	  nav[j][i]= navData[j*c.dx + i];
      if (c.nav_flip_odd || c.nav_flip_even) 
	flip_rows( nav, c.nav_flip_odd, c.nav_flip_even, 0, c.dn-1 );

      opp_phase_shift( nav, c.dn, c.dx, navPhaseShift, tmp_nav );
      opp_loc_shift( tmp_nav, c.dn, c.dx, navLocShift, tmp2_nav );
      mri_set_chunk( wOutput, "navigator", navImageSize,
		     (t.t*c.dz+t.z)*navImageSize,
		     MRI_FLOAT, *tmp2_nav);

      /* ...and apply the changes to the image chunk */
      *image= (FComplex*)mri_get_chunk(wInput,"images",
				       imageSize,
				       (t.t*c.dz+t.z)*imageSize,
				       MRI_FLOAT);
      realign_matrix( (void**)image, c.dy, (long)(c.dx*sizeof(FComplex)) );
      if (c.flip_odd || c.flip_even) 
	flip_rows( image, c.flip_odd, c.flip_even, 0, c.dy-1 );
      opp_phase_shift( image, c.dy, c.dx, imagePhaseShift, tmp_image );
      opp_loc_shift( tmp_image, c.dy, c.dx, imageLocShift, tmp2_image );
      mri_set_chunk( wOutput, "images", imageSize,
		     (t.t*c.dz+t.z)*imageSize,
		     MRI_FLOAT, *tmp2_image);
    }
    break;
  default: Abort("%s: internal error: unknown task type %d!\n",
		 progname,(int)t.type);
  }
}

static void phase_est( double* phaseShift, double* absolutePhase )
{
  long dnEffective = ( c.use_yCtr ? c.dn+1 : c.dn );
  long i;
  long j;
  double evenPhase= 0.0;
  int evenCount= 0;
  double oddPhase= 0.0;
  int oddCount= 0;
  double thisPhase;
  FComplex sum;

  /* Estimate the phase of the echo.  This works because the echo
   * magnitude is highly peaked.
   */

  for (j=0; j<dnEffective; j++) {
    sum.real= sum.imag= 0.0;
    for (i=0; i<c.dx; i++) {
      double mag= Modulus( nav[j][i] );
      sum.real += mag*nav[j][i].real;
      sum.imag += mag*nav[j][i].imag;
    }
    if (sum.imag != 0.0 || sum.real != 0.0) {
      thisPhase= atan2(sum.imag, sum.real);
      if (thisPhase<0.0) thisPhase += 2.0*M_PI; /* avoids averaging problems */
    }
    else thisPhase= 0.0;
    if( j % 2 ) {
      oddPhase += thisPhase;
      oddCount++;
    }
    else {
      evenPhase += thisPhase;
      evenCount++;
    }
  }
  evenPhase /= evenCount;
  oddPhase /= oddCount;
  *phaseShift= 0.5*(evenPhase-oddPhase);
  *absolutePhase= 0.5*(evenPhase+oddPhase);
}

static double loc_est()
{
  long dnEffective = ( c.use_yCtr ? c.dn+1 : c.dn );
  long i;
  long j;
  double evenSum= 0.0;
  double oddSum= 0.0;
  double evenMagSum= 0.0;
  double oddMagSum= 0.0;
  double result;

  /* Estimate the location of the echo.  This works because the echo
   * magnitude is highly peaked.
   */

  for (j=0; j<dnEffective; j++) {
    if( j % 2 ) {
      for (i=0; i<c.dx; i++) {
	double mag= Modulus( nav[j][i] );
	oddMagSum += mag;
	oddSum += mag*(double)i;
      }
    }
    else {
      for (i=0; i<c.dx; i++) {
	double mag= Modulus( nav[j][i] );
	evenMagSum += mag;
	evenSum += mag*(double)i;
      }
    }
  }
  evenSum /= evenMagSum;
  oddSum /= oddMagSum;
  result= 0.5*(evenSum-oddSum);

  return result;
}

void flip_rows(FComplex** image, const int flipOdd, const int flipEven, 
	       const int yLow, const int yHigh)
{
  static FComplex* scratch= NULL;
  static long dx_old= 0;
  long x, y;

  /* This routine reverses rows of the image if requested */

  if (dx_old != c.dx) {
    if (scratch) free(scratch);
    if (!(scratch= (FComplex*)malloc(c.dx*sizeof(FComplex))))
      Abort("Unable to allocate %d float complex pairs!\n", c.dx);
  }
  
  for (y=yLow; y<=yHigh; y++) {
    if (y%2) {
      if (flipOdd) {
	for (x=0; x<c.dx; x++) scratch[x]= image[y][x];
	for (x=0; x<c.dx; x++) image[y][c.dx-(x+1)]= scratch[x];
      }
    }
    else {
      if (flipEven) {
	for (x=0; x<c.dx; x++) scratch[x]= image[y][x];
	for (x=0; x<c.dx; x++) image[y][c.dx-(x+1)]= scratch[x];
      }
    }
  }
}

/* Function which performs opposite phase shifts on odd and even lines */
static void opp_phase_shift( FComplex** image, long dy, long dx, double phase, 
			     FComplex** corr_image )
{
  long y, x;
  double cosph, sinph;

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
}

/* Function which performs opposite phase shifts on odd and even lines */
static void opp_loc_shift( FComplex** image, long dy, long dx, double shift, 
			   FComplex** corr_image )
{
  long i;
  long j;

  for( j = 0; j < dy; j++ )
    for( i = 0; i < dx; i++ ) {
      corr_image[j][i]= image[j][i];
    }

  fft2d( corr_image, dx, dy, -1, 'x', 0, dy-1 );
  flip_rows( corr_image, 1, 0, 0, dy-1 );
  fshrot3d_set_shift_phases( 0.0, shift, 0.0, 
			     *corr_image, dy, dx, 1,
			     1.0, 1.0, 1.0, 0, 0 );
  flip_rows( corr_image, 1, 0, 0, dy-1 );
  fft2d( corr_image, dx, dy, 1, 'x', 0, dy-1 );
}

/* UTILITY FUNCTIONS */


static void
PackContext ()
{
  par_pkstr(c.input_file);
  par_pkstr(c.output_file);
  par_pklong(c.dv);
  par_pklong(c.dx);
  par_pklong(c.dy);
  par_pklong(c.dz);
  par_pklong(c.dt);
  par_pklong(c.dn);
  par_pklong(c.yCtr);
  par_pkint(c.use_yCtr);
  par_pkint(c.flip_odd);
  par_pkint(c.flip_even);
  par_pkint(c.nav_flip_odd);
  par_pkint(c.nav_flip_even);
  par_pkint((int)c.alg);
}

static void
UnpackContext ()
{
  par_upkstr(c.input_file);
  par_upkstr(c.output_file);
  c.dv= par_upklong();
  c.dx= par_upklong();
  c.dy= par_upklong();
  c.dz= par_upklong();
  c.dt= par_upklong();
  c.dn= par_upklong();
  c.yCtr= par_upklong();
  c.use_yCtr= par_upkint();
  c.flip_odd= par_upkint();
  c.flip_even= par_upkint();
  c.nav_flip_odd= par_upkint();
  c.nav_flip_even= par_upkint();
  c.alg= (Algorithm)par_upkint();
}

static void
PackTask ()
{
  par_pkint((int)t.type);
  par_pklong(t.z);
  par_pklong(t.t);
  par_pkdouble(t.phaseShift);
  par_pkdouble(t.locShift);
  par_pkdouble(t.absolutePhase);
}

static void
UnpackTask ()
{
  t.type= (TaskType)par_upkint();
  t.z= par_upklong();
  t.t= par_upklong();
  t.phaseShift= par_upkdouble();
  t.locShift= par_upkdouble();
  t.absolutePhase= par_upkdouble();
}

static void
PackResult ()
{
  par_pkint((int)r.type);
  par_pklong(r.z);
  par_pklong(r.t);
  par_pkdouble(r.phaseShift);
  par_pkdouble(r.locShift);
  par_pkdouble(r.absolutePhase);
}

static void
UnpackResult ()
{
  r.type= (TaskType)par_upkint();
  r.z= par_upklong();
  r.t= par_upklong();
  r.phaseShift= par_upkdouble();
  r.locShift= par_upkdouble();
  r.absolutePhase= par_upkdouble();
}


