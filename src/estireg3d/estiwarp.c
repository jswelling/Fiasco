/************************************************************
 *                                                          *
 *  estiwarp.c                                             *
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
 *      3/04: estireg3d -> estiwarp, Joel Welling           *
 ************************************************************/

/* Notes-
 * -clean out JENN comments!
 */


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

static char rcsid[] = "$Id: estiwarp.c,v 1.10 2007/04/19 22:32:28 welling Exp $";

/* Notes-
 */

typedef struct warppar_struct {
  Transform v;
  double mse;
} WarpPar;


typedef struct Context {
  /* NOTE: any new fields added to this struct should
     also be added to PackContext and UnpackContext */ 

  char progname[512];
  int debugLevel;
  int verbose;

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
  WarpPar warp_par;
} Result;

typedef  double (*ObjectiveFunctionType)(Algorithm* alg, FComplex* moved, 
					 float* align, float* weight, int* mask, 
					 char* check, int dx, int dy, int dz);


/* GLOBAL VARIABLES FOR MASTER & SLAVE */
Task t;
Context c;
Result r;

/* GLOBAL VARIABLES FOR MASTER */
FILE *pf;
int weightUpToDate= 0;
WarpPar *warppar;       /* array holding the registration parameters
			   being computed */
unsigned char** alignMissing= NULL; /* will hold missing info for Align ds */
unsigned char** stdvMissing= NULL; /* will hold missing info for Stdv ds */

/* GLOBAL VARIABLES FOR SLAVE */
MRI_Dataset *sInput = NULL;
FComplex* raw_image;
FComplex* raw_prerotated;
FComplex* moved_image;
char* check;
int* mask= NULL;
int maskUpToDate= 0;
ScalarFunction* targetFunc= NULL;
WarpPar prewarp;
Transform inv_prewarp_trans;
unsigned char** missing= NULL; /* will hold missing info */

/* FORWARD DECLARATIONS FOR MASTER */
void MasterTask(int argc, char **argv);
void MasterResult(int task_number);
static FILE* initParFile(char* parfname, char* infname, 
			 char* alignfname, char* stdvfname,
			 double xv, double yv, double zv);
static void saveWarpParams(FILE* ofile, WarpPar* p, long t);

/* FORWARD DECLARATIONS FOR SLAVE */
void SlaveContext();
void SlaveTask();
void SlaveFinalize();
static void estimateAlignment(WarpPar* par);
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

#define DEFAULT_SEARCH_ALGORITHM \
   "nosmooth,noinplane,opt=praxis,weight=smoothalign,qual=ssqr,inner=trilin,outer=trilin,norotonly,notransonly,noxonly,obj=mse,wtfloor=0.0"

static int getCurrentNDim(Algorithm* alg)
{
  if (alg->x_only_flag) {
    return 2;
  }
  else if (alg->inplane_flag) {
    return 6;
  }
  else {
    return 12;
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
  struct rusage start_rusage, end_rusage;
  char run_time_string[256];
  FILE *pp;
  time_t tm;
  long i;

  c.debugLevel=0;
  c.verbose= 0;

  strcpy(c.progname, argv[0]);

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

  cl_scan( argc, argv );

  /* Get filenames */
  cl_get( "input|i", "%option %s[%]", "input.mri", c.input_file );
  cl_get( "align|a", "%option %s[%]", "align.mri", alignfile );
  cl_get( "stdv|s", "%option %s[%]", "stdv.mri", stdvfile );
  cl_get( "parameters|p", "%option %s[%]", "warp.par", parfile );
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
  if (cl_present("v|verbose")) c.verbose= 1;
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

  /* Turn on diagnostics if requested */
  if (c.debugLevel) {
    algSetDebugLevel(&(c.alg),c.debugLevel);
  }

  /* clear utility's operation counters.*/
  linwarp_clear_counts();

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
  if (!(warppar= (WarpPar*)malloc(c.dt*sizeof(WarpPar))))
    Abort("%s: unable to allocate %d bytes!\n",c.dt*sizeof(WarpPar));

  /* Init parameter file */
  pf= initParFile(parfile, c.input_file, alignfile, stdvfile,
		  xvoxel, yvoxel, zvoxel );

  /* JENN: figure this out later */
  /* Let's keep track of alignment time */
  /* getrusage(RUSAGE_SELF,&start_rusage); */

  par_set_context();

  /* Loop through images, doing alignment */
  for( t.t = 0; t.t < c.dt; t.t++ ) {
    par_delegate_task();
  }
  par_finish();
  
  /* JENN: figure this out later */
  /*   getrusage(RUSAGE_SELF,&end_rusage);
   *   buildTimeString(run_time_string,sizeof(run_time_string),&start_rusage,&end_rusage);
   *   fprintf(pf,"# %d images in %s\n",c.dt,run_time_string);
   */

  /* fprintf(pf,"# File finished writing %s", asctime(localtime(&tm)));*/


  /* Clean up */
  mri_close_dataset( Input );
  fclose(pf);

  pp= initParFile(parfile, c.input_file, alignfile, stdvfile,
		  xvoxel, yvoxel, zvoxel );
  tm=time(NULL);
  fprintf(pp,"# File finished writing %s", asctime(localtime(&tm)));
  fflush(pp);
  for(t.t = 0; t.t<c.dt; t.t++) {
    saveWarpParams(pp, &(warppar[t.t]), t.t);
  }
  fclose(pp);
  
  algFinalize(&(c.alg));
  free(c.align_image);
  free(c.weight_image);
  free(warppar);

  /* JENN: figure this out later */
#ifdef JENN
  /* Write out shear and linrot3d counts */
  linwarp_get_counts(&lin_count);
  Message("# Total transforms by linear interpolation: %d\n",
	  lin_count );
  Message( "# Image registration complete in %s\n", run_time_string );
#endif

  exit(0);

}

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

  warppar[r.t] = r.warp_par;

  saveWarpParams( pf, &r.warp_par, r.t );
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
  fprintf(result,"##Format: names:(a00,a01,a02,a03,a10,a11,a12,a13,");
  fprintf(result,"a20,a21,a22,a23,a30,a31,a32,a33,mse)\n");
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

static void saveWarpParams(FILE* ofile, WarpPar* p, long t) {
  fprintf(ofile, 
    "%ld %11.5lg %11.5lg %11.5lg %11.5lg %11.5lg %11.5lg %11.5lg %11.5lg %11.5lg %11.5lg %11.5lg %11.5lg %11.5lg %11.5lg %11.5lg %11.5lg %11.5lg\n",
	  t, 
	  p->v[0], p->v[1], p->v[2], p->v[3], p->v[4], p->v[5], p->v[6], 
	  p->v[7], p->v[8], p->v[9], p->v[10], p->v[11], p->v[12], p->v[13], 
	  p->v[14],p->v[15], 
	  p->mse);
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
}

void SlaveContext ()
{
  ObjectiveFunctionType objectiveFunction= NULL;

  /* Open copy of input file for slaves to read. */
  sInput = mri_open_dataset( c.input_file, MRI_READ );
  missing= get_missing(sInput);

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

  /* Based on the algorithm information in the context, build the
   * ScalarFunction object that will be optimized.
   */
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
  trans_identity(r.warp_par.v);
  r.warp_par.mse= -1.0;

  /* optimize if any slice is not missing */
  all_slices_missing= 1;
  for (i=0; i<c.dz; i++)
    if (!missing[t.t][i]) {
      all_slices_missing= 0;
      break;
    }
  if (!all_slices_missing)
    estimateAlignment( &r.warp_par);
}

static void getParFromDGuess(WarpPar* par, const double* guess)
{
  int i;
  if (c.alg.x_only_flag) {
    par->v[0]= guess[0];
    par->v[1]= 0.0;
    par->v[2]= 0.0;
    par->v[3]= guess[1];
    par->v[4]= 0.0;
    par->v[5]= 1.0;
    par->v[6]= 0.0;
    par->v[7]= 0.0;
    par->v[8]= 0.0;
    par->v[9]= 0.0;
    par->v[10]= 1.0;
    par->v[11]= 0.0;
    par->v[12]= 0.0;
    par->v[13]= 0.0;
    par->v[14]= 0.0;
    par->v[15]= 1.0;
  }
  else if (c.alg.inplane_flag) {
    par->v[0]= guess[0];
    par->v[1]= guess[1];
    par->v[2]= 0.0;
    par->v[3]= guess[2];
    par->v[4]= guess[3];
    par->v[5]= guess[4];
    par->v[6]= 0.0;
    par->v[7]= guess[5];
    par->v[8]= 0.0;
    par->v[9]= 0.0;
    par->v[10]= 1.0;
    par->v[11]= 0.0;
    par->v[12]= 0.0;
    par->v[13]= 0.0;
    par->v[14]= 0.0;
    par->v[15]= 1.0;
  }
  else {
    par->v[0]= guess[0];
    par->v[1]= guess[1];
    par->v[2]= guess[2];
    par->v[3]= guess[3];
    par->v[4]= guess[4];
    par->v[5]= guess[5];
    par->v[6]= guess[6];
    par->v[7]= guess[7];
    par->v[8]= guess[8];
    par->v[9]= guess[9];
    par->v[10]= guess[10];
    par->v[11]= guess[11];
    par->v[12]= 0.0;
    par->v[13]= 0.0;
    par->v[14]= 0.0;
    par->v[15]= 1.0;
  }
}

static int getDGuessFromPar( double* guess, const WarpPar* par )
{
  if (c.alg.x_only_flag) {
    guess[0]= par->v[0];
    guess[1]= par->v[3];
  }
  else if (c.alg.inplane_flag) {
    guess[0]= par->v[0];
    guess[1]= par->v[1];
    guess[2]= par->v[3];
    guess[3]= par->v[4];
    guess[4]= par->v[5];
    guess[5]= par->v[7];

  }
  else {
    guess[0]= par->v[0];
    guess[1]= par->v[1];
    guess[2]= par->v[2];
    guess[3]= par->v[3];
    guess[4]= par->v[4];
    guess[5]= par->v[5];
    guess[6]= par->v[6];
    guess[7]= par->v[7];
    guess[8]= par->v[8];
    guess[9]= par->v[9];
    guess[10]= par->v[10];
    guess[11]= par->v[11];
  }
  return getCurrentNDim(&(c.alg));
}

static void estimateAlignment( WarpPar* par)
{
  double dguess[MAX_DOF];
  double mse= -1.0;
  long x,y,z;

  trans_identity(prewarp.v);
  /* The inverse of the identity is the identity */
  trans_copy(inv_prewarp_trans,prewarp.v);

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
      Message("Alignment rslt (X axis only): this trans produces mse %lg\n",
	      par->mse);
      trans_dump(stderr,par->v);
    }
    else if (c.alg.inplane_flag) {
      Message("Alignment rslt (in-plane): this trans produces mse %lg\n",
	      par->mse);
      trans_dump(stderr,par->v);
    }
    else {
      Message("Alignment rslt: this trans produces mse %lg\n",
	      par->mse);
      trans_dump(stderr,par->v);
    }
  }
  return;
}

static void map( Transform par, Vec4 pt )
{
  trans_vec_mult( par, pt );
  if (pt[3]==0.0) 
    Abort("%s: map: degenerate homogeneous rescaling!");
  if (pt[3] != 1.0) {
    int i;
    for (i=0; i<3; i++) pt[i] /= pt[3];
    pt[3]= 1.0;
  }
}

/* Returns mean squared-error between (adjusted) fixed image  */
/*   and (adjusted) register image, the criterion to be       */
/*   minimized by the minization op for registration */
static double mse( const double* guess, const int npar, void* userHook )
{
  WarpPar lclPar;
  Vec4 loc;
  double xshift, yshift, zshift;
  double chisqr;
  static double reallyBig= 0.0;
  ObjectiveFunctionType objectiveFunction= (ObjectiveFunctionType)userHook;

  if (reallyBig==0.0) {
      reallyBig= sqrt(SLAMCH("o"));
  }

  /* Unpack the input info */
  getParFromDGuess(&lclPar,guess);

  loc[0]= loc[1]= loc[2]= 0.0;
  loc[3]= 1.0;
  map(lclPar.v,loc);
  xshift= loc[0];
  yshift= loc[1];
  zshift= loc[2];

  /* Put up a barrier against search walking out of the home k-space cell */
  if ((xshift>0.5*c.dx) || (xshift<-0.5*c.dx) 
      || (yshift>0.5*c.dy) || (yshift<-0.5*c.dy) 
      || (zshift>0.5*c.dz) || (zshift<-0.5*c.dz)) {
    chisqr= reallyBig*(xshift*xshift + yshift*yshift + zshift*zshift);
    if (c.debugLevel) {
      Message("mse: this transform breaks cell boundary-> chisqr= %g\n",
	      chisqr);
      trans_dump(stderr,lclPar.v);
    }
  }
  else {

    /* Compensate for current prewarp */
    trans_mult_right(lclPar.v, inv_prewarp_trans);

    /* Move the raw data into position */
    switch (c.alg.inner_search_method) {
    case SEARCH_TRILIN:
      linwarp_warp( lclPar.v, raw_prerotated, moved_image,
		    check, c.dx, c.dy, c.dz, c.lx, c.ly, c.lz, 0 );
      break;
    case SEARCH_INVALID:
      Abort("%s: internal error: bad search method in mse!\n",c.progname);
    }
    
    chisqr= objectiveFunction(&(c.alg), moved_image, c.align_image, c.weight_image,
			      mask, check, c.dx, c.dy, c.dz);
    
    if (c.debugLevel>1) {
      Message("mse: adjusted transform:\n");
      trans_dump(stderr,lclPar.v);
      Message("mse: this trans gives chisqr= %.14g\n",chisqr);
    }
  }
  
  return( chisqr );
}

static void restrt( const double* guess, const int npar, void* userHook )
{
  if (c.alg.outer_search_method == c.alg.inner_search_method) { 
    restrtMatched( guess, npar, userHook );
  }
  else {
    WarpPar lclPar;
    Transform test;

    /* Unpack the input info */
    getParFromDGuess(&lclPar,guess);

    /* Initialize the prewarp transformation and its inverse */
    trans_copy(prewarp.v,lclPar.v);
    if (!trans_inverse(inv_prewarp_trans,prewarp.v)) {
      Message("%s: Fatal error: cannot invert this transform:\n",
	      c.progname);
      trans_dump(stderr,prewarp.v);
      Abort("%s: transform has become non-invertible at time %d!\n",
	    c.progname,t.t);
    }

    if (c.debugLevel) {
      Message("Restart! prewarp transform follows:\n");
      trans_dump(stderr,prewarp.v);
    }
    
    switch (c.alg.outer_search_method) {
    case SEARCH_TRILIN:
      linwarp_warp( lclPar.v, raw_image, raw_prerotated,
		    check, c.dx, c.dy, c.dz, c.lx, c.ly, c.lz, 0 );
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
   
}

static void restrtMatched( const double* guess, const int npar, 
			   void* userHook )
{
  /* No point in moving the image twice by the same method; it just
   * wastes time and slows convergence.  We'll specify a prewarp which
   * is the identity, essentially making the prewarped image
   * identical to the raw image.
   */
  trans_identity(prewarp.v);
  /* The inverse of the identity is the identity */
  trans_copy(inv_prewarp_trans,prewarp.v);

  if (c.debugLevel) {
    Message("RestartMatched! prewarp transform follows:\n");
    trans_dump(stderr,prewarp.v);
  }

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

  par_pkstr(c.progname);
  par_pkint(c.debugLevel);
  par_pkint(c.verbose);
  
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
  c.verbose= par_upkint();
  
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
  par_pkdoublearray((double*)&r.warp_par.v, 16);
  par_pkdouble(r.warp_par.mse);
}

void
UnpackResult ()
{
  
  r.t= par_upkint();
  par_upkdoublearray((double*)&r.warp_par.v, 16);
  r.warp_par.mse= par_upkdouble();
  
}
