/************************************************************
 *                                                          *
 *  estireg3d.c                                             *
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
 *  Original programming by Joel Welling 3/00               *
 *      10/02: Parallelization, Jenn Bakal                  *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifdef DARWIN
#include <sys/time.h>
#endif
#include <time.h>
#include <math.h>
#include <limits.h>
#include <sys/resource.h>
#include <unistd.h>
#include "../fmri/lapack.h"
#include "mri.h"
#include "fmri.h"
#include "par.h"
#include "stdcrg.h"
#include "misc.h"
#include "array.h"
#include "acct.h"
#include "algorithm.h"
#include "estireg_utils.h"


static char rcsid[] = "$Id: estireg3d.c,v 1.46 2007/04/19 22:32:28 welling Exp $";


/* Notes-
   -We will always take the positive value of w.  This should be OK for 
    typical fMRI displacements.
 */

typedef struct regpar3d_struct {
  Quat q;
  double x;
  double y;
  double z;
  double mse;
} RegPar3D;


typedef struct Context {
  /* NOTE: any new fields added to this struct should
     also be added to PackContext and UnpackContext */ 

  char progname[512];
  int debugLevel;

  Filename input_file;		/* the images to work on */
  float* align_image;           /* the image to align to */
  float* weight_image;          /* weight information    */

  int dx;			/* # of voxels along the X dimension */
  int dy;			/* # of voxels along the Y dimension */
  int dz;			/* # of slices along the Z dimension */
  int dt;			/* # of images along the T dimension */
  double lx;
  double ly;
  double lz;

  Algorithm alg;
} Context;


typedef struct Task {
  /* NOTE: any new fields added to this struct should
     also be added to PackTask and UnpackTask */ 
  int t;			/* image number to work on */
} Task;

typedef struct Result {
  /* NOTE: any new fields added to this struct should
     also be added to PackResult and UnpackResult */
  int t;
  RegPar3D reg_pars;
  float mean_squared_error;     /* negative means optimization failed */
} Result;

typedef  double (*ObjectiveFunctionType)(Algorithm* alg, FComplex* moved, 
					 float* align, float* weight, int* mask, 
					 char* check, int dx, int dy, int dz);

/* GLOBAL VARIABLES FOR MASTER & SLAVE */
Task t;
Context c;
Result r;
unsigned char** alignMissing= NULL; /* will hold missing info for Align ds */
unsigned char** stdvMissing= NULL; /* will hold missing info for Stdv ds */

/* GLOBAL VARIABLES FOR MASTER */
FILE *pf;
int weightUpToDate= 0;

RegPar3D *regpar;       /* array holding the registration parameters
			   being computed */

/* GLOBAL VARIABLES FOR SLAVE */
MRI_Dataset *sInput = NULL;
FComplex* raw_image= NULL;
FComplex* raw_prerotated= NULL;
FComplex* moved_image= NULL;
char* check= NULL;
int* mask= NULL;
int maskUpToDate= 0;
ScalarFunction* targetFunc= NULL;
Quat prerotation;
double preshift[3];
unsigned char** missing= NULL; /* will hold missing info */
#ifdef never
FILE* tp;
#endif


/* FORWARD DECLARATIONS FOR MASTER */
void MasterTask(int argc, char **argv);
void MasterResult(int task_number);
static FILE* initParFile(char* parfname, char* infname, 
			 char* alignfname, char* stdvfname,
			 double xv, double yv, double zv);
static void saveRegParams(FILE* ofile, RegPar3D* p, long t);

/* FORWARD DECLARATIONS FOR SLAVE */
void SlaveContext();
void SlaveTask();
void SlaveFinalize();
static void estimateAlignment(RegPar3D* par);
static double mse(const double* guess, const int npar, void* userHook);
static void restrt(const double* guess, const int npar, void* userHook);
static void restrtMatched(const double* guess, const int npar, void* userHook);

/* FORWARD DECLARATIONS - UTILITY */
void PackContext ();
void UnpackContext ();
void PackTask ();
void UnpackTask ();
void PackResult ();
void UnpackResult ();

/* Weights below alg.weight_floor as a fraction of maximum weight are ignored.
 * Non-zero values here seem to speed mse() but slow convergence by
 * making the optimization surface bumpier.
 */
#define DEFAULT_SEARCH_ALGORITHM \
   "nosmooth,noinplane,norotonly,notransonly,noxonly,opt=praxis,weight=smoothalign,qual=ssqr,inner=shear4,outer=shear4,obj=mse,wtfloor=0.0"

static void regpar3d_copy(RegPar3D* out, RegPar3D* in)
{
  quat_copy(&(out->q),&(in->q));
  out->x= in->x;
  out->y= in->y;
  out->z= in->z;
  out->mse= in->mse;
}

static void regpar3d_invert(RegPar3D* p)
{
  Transform t;
  Transform tInv;
  Vec4 vec;
  if (c.debugLevel)
    fprintf(stderr,"Inverting (%g %g %g %g) %g %g %g\n",
	    p->q.x, p->q.y, p->q.z, p->q.w, p->x, p->y, p->z);
  quat_to_trans(t,&(p->q),p->x,p->y,p->z);
  if (!trans_inverse(tInv,t)) 
    Abort("%s: unable to invert a rotation transform!\n",c.progname);
  trans_to_quat(&(p->q),tInv);
  p->x= tInv[3];
  p->y= tInv[7];
  p->z= tInv[11];
}

static int getCurrentNDim(Algorithm* alg)
{
  if (alg->x_only_flag) {
    /* minimize only in X direction- 1D */
    return 1;
  }
  else if (alg->inplane_flag) {
    /* minimize in-plane only */
    if (alg->rot_only_flag) {
      /* 1 degree of freedom */
      return 1;
    }
    else if (alg->trans_only_flag) {
      /* 2 degrees of freedom */
      return 2;
    }
    else {
      /* 3 degrees of freedom */
      return 3;
    }
  }
  else {
    if (alg->rot_only_flag) {
      /* minimize on 3 degrees of freedom */
      return 3;
    }
    else if (alg->trans_only_flag) {
      /* minimize on 3 degrees of freedom */
      return 3;
    }
    else {
      /* minimize on all 6 degrees of freedom */
      return 6;
    }
  }
}

int main (int argc,
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
	    char **argv)
{
  MRI_Dataset *Input = NULL, *Align = NULL, *Stdv = NULL;
  char infile[512], alignfile[512], stdvfile[512], parfile[512];
  char search_string[512]; 
  double xvoxel, yvoxel, zvoxel;
  int lin_count;
  int shear_count_x, shear_count_y, shear_count_z, phase_count;
  double shears_per_call;
  struct rusage start_rusage, end_rusage;
  char run_time_string[256];
  FILE *pp;
  time_t tm;
  long i;

  c.debugLevel=0;

  strcpy(c.progname, argv[0]);
#ifdef never
  fprintf(stderr, "You're running %s.\n", c.progname);
#endif

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

  algPrep();

  /*** Parse command line ***/

#ifdef never
  fprintf(stderr, "Starting command line parsing.");
#endif

  cl_scan( argc, argv );

  /* Get filenames */
  cl_get( "input|i", "%option %s[%]", "input.mri", c.input_file );
  cl_get( "align|a", "%option %s[%]", "align.mri", alignfile );
  cl_get( "stdv|s", "%option %s[%]", "stdv.mri", stdvfile );
  cl_get( "parameters|p", "%option %s[%]", "reg3d.par", parfile );
  if (!cl_get( "x|xv|xvoxel", "%option %lf", &xvoxel )) {
    fprintf(stderr,"%s: required argument xvoxel omitted.\n", c.progname);
    Help("usage");
    exit(-1);
  }
  if (!cl_get( "y|yv|yvoxel", "%option %lf", &yvoxel )) {
    fprintf(stderr,"%s: required argument yvoxel omitted.\n",c.progname);
    Help("usage");
    exit(-1);
  }
  if (!cl_get( "z|zv|zvoxel", "%option %lf", &zvoxel )) {
    fprintf(stderr,"%s: required argument zvoxel omitted.\n",c.progname);
    Help("usage");
    exit(-1);
  }
  if (cl_present("debug")) c.debugLevel= 1;
  if (cl_present("DEBUG")) c.debugLevel= 2;
  cl_get( "algorithm|alg", "%option %s[%]", DEFAULT_SEARCH_ALGORITHM, 
	  search_string );

  /* Check for args dealing with components held by c.alg */
  algParseOpts(&(c.alg));

  if (cl_cleanup_check()) {
    int i;
    fprintf(stderr,"%s: invalid argument in command line:\n    ",c.progname);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(-1);
  }

  /*** End command-line parsing ***/

  /* The algorithm string initializes most of the components of c.alg */
  if (!algInit(&(c.alg), DEFAULT_SEARCH_ALGORITHM, c.progname, 
	       getCurrentNDim))
    Abort("%s: internal error: invalid default algorithm!\n",c.progname);
  if (!algParseInfoString(&(c.alg), search_string)) {
    fprintf(stderr,"%s: invalid search method specified.\n",c.progname);
    Help("usage");
    exit(-1);
  }
  Message("# algorithm: %s\n",algGetInfoString(&(c.alg)));

  /* Turn on diagnostics if requested, and clear fshrot3d's
   * operation counters.
   */
  if (c.debugLevel>1) {
    algSetDebugLevel(&(c.alg),c.debugLevel);
#ifdef never
    fshrot3d_set_debug(c.debugLevel);
    linrot3d_set_debug(c.debugLevel);
#endif
  }

  fshrot3d_clear_shear_counts();
  linrot3d_clear_counts();

  /* Open input datasets */
  if ( !strcmp( c.input_file, alignfile )
       || (algNeedsStdv(&(c.alg)) 
	   && !strcmp( c.input_file,stdvfile ))
       || (algNeedsStdv(&(c.alg)) 
	   && !strcmp( alignfile, stdvfile )) ) {
    Abort( "Input, alignment, and stdv files must be distinct.\n" );
  }
  Input = mri_open_dataset( c.input_file, MRI_READ );
  Align = mri_open_dataset( alignfile, MRI_READ );
  alignMissing= get_missing(Align);
  if (algNeedsStdv(&(c.alg))) {
    Stdv = mri_open_dataset( stdvfile, MRI_READ );
    stdvMissing= get_missing(Stdv);
  }
  else {
    Stdv= NULL;
    stdvMissing= NULL;
  }

  /* Check that program will function on data-sets */
  if (!checkDatasetDims("input",Input,"images","xyzt",c.progname)
      || !checkDatasetDims("align",Align,"images","xyz",c.progname)
      || (algNeedsStdv(&(c.alg))
	  && !checkDatasetDims("stdv",Stdv,"images","xyz",c.progname))) {
    Abort( "An input dataset was not of the correct type.\n" );
  }

  /* Set parameters in local variables */
  c.dx= mri_get_int( Input, "images.extent.x" );
  c.dy= mri_get_int( Input, "images.extent.y" );
  c.dz= mri_get_int( Input, "images.extent.z" );
  c.dt= mri_get_int( Input, "images.extent.t" );

#ifdef never
  fprintf(stderr, "dx = %d, dy = %d, dz = %d, dt = %d\n", c.dx, c.dy, c.dz, c.dt); 
#endif

  /* Check that relevant dimensions are commensurate */
  if ((mri_get_int( Align, "images.extent.x" ) != c.dx)
      || (algNeedsStdv(&(c.alg))
	  && (mri_get_int( Stdv, "images.extent.x" ) != c.dx))) {
    Abort( "Input dataset x extents are not commensurate.\n" );
  }
  if ((mri_get_int( Align, "images.extent.y") ) != c.dy
      || (algNeedsStdv(&(c.alg))
	  && (mri_get_int( Stdv, "images.extent.y") != c.dy))) {
    Abort( "Input dataset y extents are not commensurate.\n" );
  }
  if ((mri_get_int( Align, "images.extent.z" ) != c.dz)
      || (algNeedsStdv(&(c.alg))
	  && (mri_get_int( Stdv, "images.extent.z" ) != c.dz))) {
    Abort( "Input dataset z extents are not commensurate.\n" );
  }

  /* Allocate image storage */
  if (!(c.align_image= (float*)malloc(c.dx*c.dy*c.dz*sizeof(float)))) 
    Abort("%s: unable to allocate %d bytes!\n",
	  c.progname,c.dx*c.dy*c.dz*sizeof(float));
  if (!(c.weight_image= (float*)malloc(c.dx*c.dy*c.dz*sizeof(float)))) 
    Abort("%s: unable to allocate %d bytes!\n",
	  c.progname,c.dx*c.dy*c.dz*sizeof(float));
  weightUpToDate= 0;
  
  /* Load target image */
  loadImage( c.align_image, Align, 0, c.dx, c.dy, c.dz );
  algMaybeSmoothImage( &(c.alg), c.align_image, c.dx, c.dy, c.dz, 
		       alignMissing, 0 );
  
  /* Prepare weight */
  algBuildWeight(&(c.alg), c.weight_image, c.align_image, Align, Stdv, 
		 c.dx, c.dy, c.dz, alignMissing, stdvMissing);
  weightUpToDate= 1;

  /* Done with those datasets */
  mri_close_dataset( Align );
  if (Stdv) mri_close_dataset( Stdv );

  /* Set up voxel edge dimensions. */
  c.lx= c.dx*xvoxel;
  c.ly= c.dy*yvoxel;
  c.lz= c.dz*zvoxel;

  /* Allocate parameter storage */
  if (!(regpar= (RegPar3D*)malloc(c.dt*sizeof(RegPar3D))))
    Abort("%s: unable to allocate %d bytes!\n",c.dt*sizeof(RegPar3D));

  /* Init parameter file */
  pf= initParFile(parfile, c.input_file, alignfile, stdvfile,
		  xvoxel, yvoxel, zvoxel );

  /* JENN: figure this out later */
  /* Let's keep track of alignment time */
  /* getrusage(RUSAGE_SELF,&start_rusage); */

  par_set_context();

#ifdef never
  fprintf(stderr, "Just set context, about to delegate tasks.\n");
#endif

  /* Loop through images, doing alignment */
  for( t.t = 0; t.t < c.dt; t.t++ ) {
    par_delegate_task();
  }
  par_finish();
#ifdef never
  fprintf(stderr, "Hey are we there yet?\n");     
#endif
  
  /* JENN: figure this out later */
  /*   getrusage(RUSAGE_SELF,&end_rusage);
   *   buildTimeString(run_time_string,sizeof(run_time_string),&start_rusage,&end_rusage);
   *   fprintf(pf,"# %d images in %s\n",c.dt,run_time_string);
   */

  /* fprintf(pf,"# File finished writing %s", asctime(localtime(&tm)));*/


  /* Clean up */
  mri_close_dataset( Input );
  fclose(pf);

#ifdef never
  fprintf(stderr, "Write out 2ndregpar file.\n");
#endif

  pp= initParFile(parfile, c.input_file, alignfile, stdvfile,
		  xvoxel, yvoxel, zvoxel );
  tm=time(NULL);
  fprintf(pp,"# File finished writing %s", asctime(localtime(&tm)));
  fflush(pp);
  for(t.t = 0; t.t<c.dt; t.t++) {
    saveRegParams(pp, &(regpar[t.t]), t.t);
  }
  fclose(pp);
  
  algFinalize(&(c.alg));
  free(c.align_image);
  free(c.weight_image);
  free(regpar);

  /* JENN: figure this out later */
#ifdef JENN
  /* Write out shear and linrot3d counts */
  linrot3d_get_counts(&lin_count);
  Message("# Total transforms by linear interpolation: %d\n",
	  lin_count );
  fshrot3d_get_shear_counts(&shear_count_x, &shear_count_y, &shear_count_z,
			    &phase_count, &shears_per_call);
  Message("# Total shears in x, y, z, rephases: %d %d %d %d (%f per call)\n",
	  shear_count_x,shear_count_y, shear_count_z, phase_count,
	  shears_per_call );
  Message( "# Image registration complete in %s\n", run_time_string );
#endif

  exit(0);
}

/* JENN: what stays in here? */
void MasterResult (int task_number)
{
  static int results= 0;

  ++results;

#ifdef never
  fprintf(stderr, "Master Result, time %d.\n", r.t);
  fprintf(stderr, "result #%d, c.dt=%d\n",results,c.dt);
#endif

  if ((results % 60 == 0 ) || (results == (c.dt))) 
    Message("# %ld\n", (long)results);
  else 
    Message("#");

#ifdef never
  fprintf(stderr, "save reg pars to array.\n");
#endif

  regpar[r.t] = r.reg_pars;

#ifdef never
  fprintf(stderr, "q.x = %11.5lg, q.y = %11.5lg, q.z = %11.5lg, q.w = %11.5lg, x = %11.5lg, y = %11.5lg, z = %11.5lg, mse = %11.5lg \n", regpar[r.t].q.x, regpar[r.t].q.y, regpar[r.t].q.z, regpar[r.t].q.w, regpar[r.t].x, regpar[r.t].y, regpar[r.t].z, regpar[r.t].mse);

  fprintf(stderr, "done saving reg pars to array. \n");
#endif

  saveRegParams( pf, &r.reg_pars, r.t );
}

static FILE* initParFile(char* parfname, char* infname, 
			 char* alignfname, char* stdvfname,
			 double xv, double yv, double zv)
{
  FILE* result;
  time_t tm;
  if (!(result= fopen(parfname,"w"))) {
    Abort("%s: unable to open <%s> for writing!\n",c.progname,parfname);
  }

  tm= time(NULL);
  fprintf(result,"##Format: order:index_t, type:raw\n");
  fprintf(result,"##Format: names:(3d_qbar_x,3d_qbar_y,3d_qbar_z,3d_qbar_w,");
  fprintf(result,"3d_deltabarx,3d_deltabary,3d_deltabarz,mse)\n");
  fprintf(result,"# Alignment parameters generated %s",
	  asctime(localtime(&tm)));
  fprintf(result,"# Input file: %s\n",infname);
  fprintf(result,"# Alignment file: %s\n",alignfname);
  fprintf(result,"# Stdv file: %s\n",stdvfname);
  fprintf(result,"# xdim= %ld, ydim= %ld, zdim= %ld\n",c.dx,c.dy,c.dz);
  fprintf(result,"# voxel size x= %g, y= %g, z= %g\n", xv, yv, zv);
  fprintf(result,"# algorithm: %s\n",algGetInfoString(&(c.alg)));
  fflush(result);
  return result;
}

static void saveRegParams(FILE* ofile, RegPar3D* p_in, long t) {
  RegPar3D p;
  regpar3d_copy(&p, p_in);
  regpar3d_invert(&p);
  fprintf(ofile, 
    "%ld %11.5lg %11.5lg %11.5lg %11.5lg %11.5lg %11.5lg %11.5lg %11.5lg\n",
	  t, p.q.x, p.q.y, p.q.z, p.q.w, p.x, p.y, p.z, p.mse);
  fflush(ofile);
}

/* SLAVE PROCEDURES */

void SlaveFinalize()
{
  if (sInput != NULL) mri_close_dataset(sInput);
  if (raw_image != NULL) free(raw_image);
  if (raw_prerotated != NULL) free(raw_prerotated);
  if (moved_image != NULL) free(moved_image);
  if (check != NULL) free(check);
  if (mask != NULL) free(mask);
  if (targetFunc != NULL) {
    targetFunc->destroySelf(targetFunc);
    targetFunc= NULL;
  }
  algFinalize(&(c.alg));

#ifdef never
  if (tp != NULL) fclose(tp);
#endif
}

void SlaveContext ()
{
  ObjectiveFunctionType objectiveFunction= NULL;

  /* Open copy of input file for slaves to read. */
#ifdef never
  fprintf(stderr, "progname in SC = %s\n", c.progname);
  fprintf(stderr, "input file in SC = %s\n", c.input_file);
#endif
  sInput = mri_open_dataset( c.input_file, MRI_READ );
  missing= get_missing(sInput);

#ifdef never
  fprintf(stderr, "dx = %d, dy = %d, dz = %d\n", c.dx, c.dy, c.dz); 
#endif

  /* Allocate image storage */
  if (!(raw_image= (FComplex*)malloc(c.dx*c.dy*c.dz*sizeof(FComplex)))) 
    Abort("%s: unable to allocate %d bytes!\n",
	  c.progname,c.dx*c.dy*c.dz*sizeof(FComplex));
  if (!(raw_prerotated= (FComplex*)malloc(c.dx*c.dy*c.dz*sizeof(FComplex)))) 
    Abort("%s: unable to allocate %d bytes!\n",
	  c.progname,c.dx*c.dy*c.dz*sizeof(FComplex));
  if (!(moved_image= (FComplex*)malloc(c.dx*c.dy*c.dz*sizeof(FComplex)))) 
    Abort("%s: unable to allocate %d bytes!\n",
	  c.progname,c.dx*c.dy*c.dz*sizeof(FComplex));
  if (!(check= (char*)malloc(c.dx*c.dy*c.dz*sizeof(char))))
    Abort("%s: unable to allocate %d bytes!\n",
	  c.progname,c.dx*c.dy*c.dz*sizeof(char));
  if (algNeedsMask(&(c.alg))) {
    if (!(mask= (int*)malloc(c.dx*c.dy*c.dz*sizeof(int))))
      Abort("%s: unable to allocate %d bytes!\n",
	    c.progname,c.dx*c.dy*c.dz*sizeof(int));
    maskUpToDate= 0; /* since weight has just been updated */
  }
#ifdef never
  if (!(tp= fopen("align_data","w"))) {
    Abort("%s: unable to open align_data for writing!\n",c.progname);
  }
#endif

  /* Based on the algorithm information in the context, build the
   * ScalarFunction object that will be optimized.
   */
  switch (c.alg.objective_method) {
  case OBJECTIVE_MSE: objectiveFunction= algCalcChiSqr; break;
  case OBJECTIVE_MI: objectiveFunction= algCalcMutualInfo; break;
  case OBJECTIVE_JE: objectiveFunction= algCalcJointEntropy; break;
  default: Abort("%s: internal error: unknown objective method!\n",
		 c.progname);
  }
  if (c.alg.inner_search_method == c.alg.outer_search_method) {
    targetFunc= buildSimpleScalarFunction( mse, restrtMatched, 
					   getCurrentNDim(&(c.alg)), 
					   (void*)objectiveFunction );
  }
  else {
    targetFunc= buildSimpleScalarFunction( mse, restrt, 
					   getCurrentNDim(&(c.alg)), 
					   (void*)objectiveFunction );
  }
}

void SlaveTask ()
{
  int i;
  int all_slices_missing;

  if (c.debugLevel)
    Message("Image %ld...\n",t);
  loadImageComplex( raw_image, sInput, t.t, c.dx, c.dy, c.dz );

  r.t = t.t;

  /* Pick a starting point */
  r.reg_pars.q.x= r.reg_pars.q.y= r.reg_pars.q.z= 0.0; 
  r.reg_pars.q.w= 1.0;
  r.reg_pars.mse= -1.0;
  r.reg_pars.x= r.reg_pars.y= r.reg_pars.z= 0.0;

  /* optimize if any slice is not missing */
  all_slices_missing= 1;
  for (i=0; i<c.dz; i++)
    if (!missing[t.t][i]) {
      all_slices_missing= 0;
      break;
    }
  if (!all_slices_missing)
    estimateAlignment( &r.reg_pars);
  
}

static void getParFromDGuess(RegPar3D* par, const double* guess)
{
  if (c.alg.x_only_flag) {
    /* This is a 1D minimization */
      par->q.x= 0.0;
      par->q.y= 0.0;
      par->q.z= 0.0;
      par->q.w= 1.0;
      par->x= c.dx*guess[0];
      par->y= 0.0;
      par->z= 0.0;
  }
  else if (c.alg.inplane_flag) {
    if (c.alg.rot_only_flag) {
      par->q.x= 0.0;
      par->q.y= 0.0;
      par->q.z= guess[0]/(1.0+fabs(guess[0]));
      par->q.w= 
	sqrt(1.0-(par->q.x*par->q.x + par->q.y*par->q.y + par->q.z*par->q.z));
      par->x= 0.0;
      par->y= 0.0;
      par->z= 0.0;
    }
    else if (c.alg.trans_only_flag) {
      par->q.x= 0.0;
      par->q.y= 0.0;
      par->q.z= 0.0;
      par->q.w= 1.0;
      par->x= c.dx*guess[0];
      par->y= c.dy*guess[1];
      par->z= 0.0;
    }
    else {
      par->q.x= 0.0;
      par->q.y= 0.0;
      par->q.z= guess[0]/(1.0+fabs(guess[0]));
      par->q.w= 
	sqrt(1.0-(par->q.x*par->q.x + par->q.y*par->q.y + par->q.z*par->q.z));
      par->x= c.dx*guess[1];
      par->y= c.dy*guess[2];
      par->z= 0.0;
    }
  }
  else {
    if (c.alg.rot_only_flag) {
      par->q.x= guess[0]/(1.0+fabs(guess[0]));
      par->q.y= guess[1]/(1.0+fabs(guess[1]));
      par->q.z= guess[2]/(1.0+fabs(guess[2]));
      par->q.w= 
	sqrt(1.0-(par->q.x*par->q.x + par->q.y*par->q.y + par->q.z*par->q.z));
      par->x= 0.0;
      par->y= 0.0;
      par->z= 0.0;
    }
    else if (c.alg.trans_only_flag) {
      par->q.x= 0.0;
      par->q.y= 0.0;
      par->q.z= 0.0;
      par->q.w= 1.0;
      par->x= c.dx*guess[0];
      par->y= c.dy*guess[1];
      par->z= c.dz*guess[2];
    }
    else {
      par->q.x= guess[0]/(1.0+fabs(guess[0]));;
      par->q.y= guess[1]/(1.0+fabs(guess[1]));;
      par->q.z= guess[2]/(1.0+fabs(guess[2]));;
      par->q.w= 
	sqrt(1.0-(par->q.x*par->q.x + par->q.y*par->q.y + par->q.z*par->q.z));
      par->x= c.dx*guess[3];
      par->y= c.dy*guess[4];
      par->z= c.dz*guess[5];
    }
  }
}

static int getDGuessFromPar( double* dguess, const RegPar3D* par )
{
  if (c.alg.x_only_flag) {
    /* 1 degree of freedom */
      dguess[0]= par->x/c.dx;
  }
  else if (c.alg.inplane_flag) {
    /* minimize in-plane only */
    if (c.alg.rot_only_flag) {
      /* 1 degree of freedom */
      dguess[0]= par->q.z/(1.0-fabs(par->q.z));
    }
    else if (c.alg.trans_only_flag) {
      /* 2 degrees of freedom */
      dguess[0]= par->x/c.dx;
      dguess[1]= par->y/c.dy;
    }
    else {
      /* 3 degrees of freedom */
      dguess[0]= par->q.z/(1.0-fabs(par->q.z));
      dguess[1]= par->x/c.dx;
      dguess[2]= par->y/c.dy;
    }
  }
  else {
    if (c.alg.rot_only_flag) {
      /* minimize on 3 degrees of freedom */
      dguess[0]= par->q.x/(1.0-fabs(par->q.x)); 
      dguess[1]= par->q.y/(1.0-fabs(par->q.y)); 
      dguess[2]= par->q.z/(1.0-fabs(par->q.z)); 
    }
    else if (c.alg.trans_only_flag) {
      /* minimize on 3 degrees of freedom */
      dguess[0]= par->x/c.dx; 
      dguess[1]= par->y/c.dy; 
      dguess[2]= par->z/c.dz; 
    }
    else {
      /* minimize on all 6 degrees of freedom */
      dguess[0]= par->q.x/(1.0-fabs(par->q.x)); 
      dguess[1]= par->q.y/(1.0-fabs(par->q.y)); 
      dguess[2]= par->q.z/(1.0-fabs(par->q.z)); 
      dguess[3]= par->x/c.dx; 
      dguess[4]= par->y/c.dy; 
      dguess[5]= par->z/c.dz; 
    }
  }
  return getCurrentNDim(&(c.alg));
}

static void estimateAlignment( RegPar3D* par)
{
  double dguess[MAX_DOF];
  double mse= -1.0;
  long x,y,z;

  quat_identity(&(prerotation));
  preshift[0]= preshift[1]= preshift[2]= 0.0;
  /* Update the raw_prerotated image to be identical to the raw image,
   * since the initial pre-transformation is the identity.
   */
  for (x=0; x<c.dx; x++)
    for (y=0; y<c.dy; y++)
      for (z=0; z<c.dz; z++) {
	MEM(raw_prerotated,c.dx,c.dy,c.dz,x,y,z).real= 
	  MEM(raw_image,c.dx,c.dy,c.dz,x,y,z).real;
	MEM(raw_prerotated,c.dx,c.dy,c.dz,x,y,z).imag= 
	  MEM(raw_image,c.dx,c.dy,c.dz,x,y,z).imag;
      }

  /* In this geometry, all the "check" data should be 1 */
  for (x=0; x<c.dx; x++)
    for (y=0; y<c.dy; y++)
      for (z=0; z<c.dz; z++) {
	MEM(check,c.dx,c.dy,c.dz,x,y,z)= 1;
      }

  /* If a mask is necessary, compute it from the weights */
  if (!maskUpToDate) {
    algMaybeBuildMask(&(c.alg),c.weight_image,mask,c.dx,c.dy,c.dz);
    maskUpToDate= 1;
  }

  getDGuessFromPar(dguess,par);
  (void)c.alg.opt->go(c.alg.opt, targetFunc, dguess, getCurrentNDim(&(c.alg)), &mse);
  getParFromDGuess(par,dguess);
  par->mse= mse;

  if (c.debugLevel) {
    if (c.alg.x_only_flag) {
      Message("Alignment rslt (X axis only): (%lg %lg %lg %lg) (%lg %lg %lg) -> mse %lg\n",
	      par->q.x, par->q.y, par->q.z, par->q.w, 
	      par->x, par->y, par->z, par->mse);
    }
    else if (c.alg.inplane_flag) {
      Message(
  "Alignment rslt (in-plane): (%lg %lg %lg %lg) (%lg %lg %lg) -> mse %lg\n",
	      par->q.x, par->q.y, par->q.z, par->q.w, 
	      par->x, par->y, par->z, par->mse);
    }
    else {
      Message("Alignment result: (%lg %lg %lg %lg) (%lg %lg %lg) -> mse %lg\n",
	      par->q.x, par->q.y, par->q.z, par->q.w, 
	      par->x, par->y, par->z, par->mse);
    }
  }
  return;
}

/* Returns mean squared-error between (adjusted) fixed image  */
/*   and (adjusted) register image, the criterion to be       */
/*   minimized by the minization op for registration */
static double mse( const double* guess, const int npar, void* userHook )
{
  Quat q;
  Quat tq;
  Quat tqshift;
  double xshift, yshift, zshift;
  double chisqr;
  static double reallyBig= 0.0;
  RegPar3D lclPar;
  ObjectiveFunctionType objectiveFunction= (ObjectiveFunctionType)userHook;

  if (reallyBig==0.0) {
      reallyBig= sqrt(SLAMCH("o"));
  }

  /* Unpack the input info */
  getParFromDGuess(&lclPar,guess);
  q= lclPar.q;
  xshift= lclPar.x;
  yshift= lclPar.y;
  zshift= lclPar.z;

  /* Put up a barrier against search walking out of the home k-space cell */
  if ((xshift>0.5*c.dx) || (xshift<-0.5*c.dx) 
      || (yshift>0.5*c.dy) || (yshift<-0.5*c.dy) 
      || (zshift>0.5*c.dz) || (zshift<-0.5*c.dz)
      || (q.x>=1.0 || q.y >= 1.0 || q.z >= 1.0)
      || (q.x<=-1.0 || q.y <= -1.0 || q.z <= -1.0)) {
    chisqr= reallyBig*(xshift*xshift + yshift*yshift + zshift*zshift
		     + q.x*q.x + q.y*q.y + q.z*q.z);
    if (c.debugLevel>1) {
      Message("mse: (%g %g %g %g %g %g) breaks cell boundary-> chisqr= %g\n",
	      q.x, q.y, q.z, xshift, yshift, zshift, chisqr);
    }
  }
  else {
    /* Compensate for current prerotation.  The first steps compute an
     * adjusted rotation q; the remainder use quaternion algebra to adjust
     * the shift.
     */
    quat_copy(&tq, &(prerotation));
    quat_conjugate(&tq);
    quat_mult_right(&q,&tq);
    quat_normalize(&q);
    tqshift.x= -preshift[0];
    tqshift.y= -preshift[1];
    tqshift.z= -preshift[2];
    tqshift.w= 0.0;
    quat_copy(&tq, &q);
    quat_mult_left(&tq,&tqshift);
    quat_conjugate(&tq);
    quat_mult_right(&tqshift,&tq);
    xshift += tqshift.x;
    yshift += tqshift.y;
    zshift += tqshift.z;

#ifdef never
    fprintf(stderr,"Prerot is (%lg %lg %lg %lg); shift ( %lg %lg %lg )\n",
	    prerotation.x,prerotation.y,
	    prerotation.z,prerotation.w,
	    preshift[0],preshift[1],preshift[2]);
    fprintf(stderr,"Adjusted quat (%lg %lg %lg %lg), shift (%lg %lg %lg)\n",
	    q.x,q.y,q.z,q.w,xshift,yshift,zshift);
#endif
    
    /* Move the raw data into position */
    switch (c.alg.inner_search_method) {
    case SEARCH_TRILIN:
      linear_shift_rot3d( &q, xshift, yshift, zshift, 
			  raw_prerotated, moved_image, 
			  check, c.dx, c.dy, c.dz,
			  c.lx, c.ly, c.lz, 
			  0 );
      break;
    case SEARCH_SHEAR4_FFT:
      fshrot3d_set( FR3D_SHEAR_PATTERN, FR3D_SHEAR_4 );
      fourier_shift_rot3d( &q, xshift, yshift, zshift, 
			   raw_prerotated, moved_image, 
			   c.dx, c.dy, c.dz,
			   c.lx, c.ly, c.lz, 1 );
      break;
    case SEARCH_SHEAR7_FFT:
      fshrot3d_set( FR3D_SHEAR_PATTERN, FR3D_SHEAR_7 );
      fourier_shift_rot3d( &q, xshift, yshift, zshift, 
			   raw_prerotated, moved_image, 
			   c.dx, c.dy, c.dz,
			   c.lx, c.ly, c.lz, 1 );
      break;
    case SEARCH_SHEAR13_FFT:
      fshrot3d_set( FR3D_SHEAR_PATTERN, FR3D_SHEAR_13 );
      fourier_shift_rot3d( &q, xshift, yshift, zshift, 
			   raw_prerotated, moved_image, 
			   c.dx, c.dy, c.dz,
			   c.lx, c.ly, c.lz, 1 );
      break;
    case SEARCH_INVALID:
      Abort("%s: internal error: bad search method in mse!\n",c.progname);
    }
    
    chisqr= objectiveFunction(&(c.alg), moved_image, c.align_image, c.weight_image,
			      mask, check, c.dx, c.dy, c.dz);
    
    if (c.debugLevel>1) {

      Message("mse: trying (%.14g %.14g %.14g %.14g %.14g %.14g) (adjusted) -> chisqr= %.14g\n",
	      q.x, q.y, q.z, xshift, yshift, zshift, chisqr);
    }
  }
  
  return( chisqr );
}


static void restrt( const double* guess, const int npar, void* userHook )
{
  RegPar3D lclPar;
  Quat q;

  /* Unpack the input info */
  getParFromDGuess(&lclPar,guess);

  if (c.alg.outer_search_method == c.alg.inner_search_method) {
    /* No point in rotating the image twice by the same method; it just
     * wastes time and slows convergence.  We'll specify a preshift of
     * 0.0 and a prerotation by 0.0 degrees, essentially making the
     * prerotated image identical to the raw image.
     */
    q.x= q.y= q.z= 0.0;
    q.w= 1.0;
    quat_copy(&prerotation,&q);
    preshift[0]= preshift[1]= preshift[2]= 0.0;
  }
  else {
    quat_copy(&q,&(lclPar.q));
    quat_copy(&prerotation, &q); 
    preshift[0]= lclPar.x;
    preshift[1]= lclPar.y;
    preshift[2]= lclPar.z;
  }
  
  if (c.debugLevel) 
    Message("Restart! prerotation (%lg %lg %lg %lg), preshift (%lg %lg %lg)\n",
	    q.x, q.y, q.z, q.w,
	    preshift[0],preshift[1],preshift[2]);

  switch (c.alg.outer_search_method) {
  case SEARCH_TRILIN:
    linear_shift_rot3d( &q, preshift[0], preshift[1],
			preshift[2],
			raw_image, raw_prerotated, 
			check, 
			c.dx, c.dy, c.dz,
			c.lx, c.ly, c.lz, 0 );
    break;
  case SEARCH_SHEAR4_FFT:
    fshrot3d_set( FR3D_SHEAR_PATTERN, FR3D_SHEAR_4 );
    fourier_shift_rot3d( &q, 
			 preshift[0],preshift[1],
			 preshift[2],
			 raw_image, raw_prerotated, 
			 c.dx, c.dy, c.dz,
			 c.lx, c.ly, c.lz, 1 );
    break;
  case SEARCH_SHEAR7_FFT:
    fshrot3d_set( FR3D_SHEAR_PATTERN, FR3D_SHEAR_7 );
    fourier_shift_rot3d( &q, 
			 preshift[0],preshift[1],
			 preshift[2],
			 raw_image, raw_prerotated, 
			 c.dx, c.dy, c.dz,
			 c.lx, c.ly, c.lz, 1 );
    break;
  case SEARCH_SHEAR13_FFT:
    fshrot3d_set( FR3D_SHEAR_PATTERN, FR3D_SHEAR_13 );
    fourier_shift_rot3d( &q, 
			 preshift[0],preshift[1],
			 preshift[2],
			 raw_image, raw_prerotated, 
			 c.dx, c.dy, c.dz,
			 c.lx, c.ly, c.lz, 1 );
    break;
  case SEARCH_INVALID:
    Abort("%s: internal error: bad search method in restrt!\n",c.progname);
  }

  if (c.debugLevel>1) {
    /* call mse to print out value using newly prerotated volumes */
    (void)mse(guess,npar,userHook);
  }
  algMaybeSmoothImageComplex( &(c.alg), raw_prerotated, c.dx, c.dy, c.dz,
			      missing, t.t );

}

static void restrtMatched( const double* guess, const int npar, 
			   void* userHook )
{
  Quat q;

  /* No point in rotating the image twice by the same method; it just
   * wastes time and slows convergence.  We'll specify a preshift of
   * 0.0 and a prerotation by 0.0 degrees, essentially making the
   * prerotated image identical to the raw image.
   */
  q.x= q.y= q.z= 0.0;
  q.w= 1.0;
  quat_copy(&prerotation,&q);
  preshift[0]= preshift[1]= preshift[2]= 0.0;
  
  if (c.debugLevel) 
    Message("RestartMatched! prerotation (%lg %lg %lg %lg), preshift (%lg %lg %lg)\n",
	    q.x, q.y, q.z, q.w,
	    preshift[0],preshift[1],preshift[2]);

  if (c.debugLevel>1) {
    /* call mse to print out value using newly prerotated volumes */
    (void)mse(guess,npar,userHook);
  }
  algMaybeSmoothImageComplex( &(c.alg), raw_prerotated, c.dx, c.dy, c.dz,
			      missing, t.t );
}

/* UTILITY FUNCTIONS */

void
PackContext ()
{
  char* s;

#ifdef never
  fprintf(stderr, "in Pack: dx = %d, dy = %d, dz = %d\n", c.dx, c.dy, c.dz); 
#endif

  par_pkstr(c.progname);
  par_pkint(c.debugLevel);
  
  par_pkstr(c.input_file);
  par_pkint(c.dx);
  par_pkint(c.dy);
  par_pkint(c.dz);
  par_pkint(c.dt);
  par_pkdouble(c.lx);
  par_pkdouble(c.ly);
  par_pkdouble(c.lz);

  algPack(&(c.alg));

  par_pkfloatarray(c.align_image, c.dx*c.dy*c.dz);
  par_pkfloatarray(c.weight_image, c.dx*c.dy*c.dz);
}

void
UnpackContext ()
{
  char tbuf[256];
  int l;

  par_upkstr(c.progname);
  c.debugLevel= par_upkint();
  
  par_upkstr(c.input_file);
  c.dx= par_upkint();
  c.dy= par_upkint();
  c.dz= par_upkint();
  c.dt= par_upkint();
  c.lx= par_upkdouble();
  c.ly= par_upkdouble();
  c.lz= par_upkdouble();

  if (!(c.align_image= (float*)malloc(c.dx*c.dy*c.dz*sizeof(float))))
    Abort("%s: unable to allocate %d bytes!\n", 
	  c.progname,c.dx*c.dy*c.dz*sizeof(float));
  if (!(c.weight_image= (float*)malloc(c.dx*c.dy*c.dz*sizeof(float)))) 
    Abort("%s: unable to allocate %d bytes!\n",
	  c.progname,c.dx*c.dy*c.dz*sizeof(float));
  
  algUnpack(&(c.alg), c.progname, getCurrentNDim);

  par_upkfloatarray(c.align_image, c.dx*c.dy*c.dz);
  par_upkfloatarray(c.weight_image, c.dx*c.dy*c.dz);
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
  par_pkdoublearray((double*)&r.reg_pars.q, 4);
  par_pkdouble(r.reg_pars.x);
  par_pkdouble(r.reg_pars.y);
  par_pkdouble(r.reg_pars.z);
  par_pkdouble(r.reg_pars.mse);
}

void
UnpackResult ()
{
  
#ifdef never
  fprintf(stderr, "Start of Unpack Result.\n");
#endif
  
  r.t= par_upkint();
  par_upkdoublearray((double*)&r.reg_pars.q, 4);
  r.reg_pars.x= par_upkdouble();
  r.reg_pars.y= par_upkdouble();
  r.reg_pars.z= par_upkdouble();
  r.reg_pars.mse= par_upkdouble();
  
#ifdef never
  fprintf(stderr, "End of Unpack Result.\n");
#endif
}
