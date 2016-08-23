/************************************************************
 *                                                          *
 *  slow_ft.c                                               *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1997 Department of Statistics,         *
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
 *  Original programming by Joel Welling                    *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <math.h>
#include "par.h"
#include "bio.h"
#include "array.h"
#include "acct.h"
#include "mri.h"
#include "fmri.h"
#include "misc.h"
#include "stdcrg.h"
#include "dirichlet.h"
#include "vpolygon.h"
#include "../fmri/lapack.h" /* for DLAMCH() */

#ifdef USE_NFFT
#include <nfft.h>
#endif

static char rcsid[] = "$Id: slow_ft.c,v 1.24 2007/03/22 00:09:18 welling Exp $";

/* Notes-
   -the correction to sampPhase should happen in a separate script?
   -Do we want _specific() versions of nfft_init?
   -There are about a thousand different nfft algorithms or limit vals to
    pick; the current choices are just a first pass.
   -Weights can change for each slice; we're just doing z=0 everywhere I think.
   -I don't want their stupid next_power_of_2() function!
   -The nfft method and the direct method disagree on which way is up.
 */

#define DEFAULT_ALGORITHM "weight=const,nft=direct"

#define LINTERP( A, B, lamda ) ( lamda*B + (1.0-lamda)*A )

/* Weighting algorithms */
typedef enum { WEIGHT_CONST, 
	       WEIGHT_VORONOI,
	       WEIGHT_CIRCVORONOI,
	       WEIGHT_INVALID } WeightMethod;

/* Overall algorithm for the nonuniform FT */
typedef enum { NFT_DIRECT,
	       NFT_NFFT,
	       NFT_INVALID } NftMethod;

typedef struct algorithm_struct {
  WeightMethod weightMethod;
  NftMethod nftMethod;
} Algorithm;

typedef struct context_struct {
  /* NOTE: any new fields added to this struct should
     be also added to PackContext and UnpackContext */ 
  Filename input_fname;
  Filename output_fname;
  Filename lagmap_fname;
  int useLagMap;
  int dp;                       /* # of samples per shot */
  int ds;                       /* # shots per coil */
  int dc;                       /* # of coils */
  int dz;                       /* # of slices */
  int dt;                       /* # of images */
  int dx;                       /* x resolution */
  int dy;                       /* y resolution */
  double xvoxel;                /* x voxel size in mm */
  double yvoxel;                /* y voxel size in mm */
  double samp_nom_fov;          /* nominal ispace FOV size of samples (mm)*/
  double ph_delta_per_sample;   /* "expected" phase shift per sample */
  double sample_lag;            /* fraction of a step of lag back along
				 * the k-space trajectory at the final 
                                 * sample in the FID */
  int verbose;                  /* verbosity level */
  int debug;                    /* debugging level */                   
  Algorithm alg;                /* Algorithm info */

/* variably-sized context info */
  double *****sampLoc;          /* dimensioned as [dz][dc][2][ds][dp] */
  double *****sampWeight;       /* dimensioned as [dz][dc][2][ds][dp] */
  double *lagMap;
} Context;

typedef struct task_struct {
  int z;   /* slice number to work on */
  int t;   /* image number to work on */
} Task;

typedef struct result_struct {
  int z;   /* slice number just completed */
  int t;   /* image number just completed */
} Result;

typedef struct voronoi_wt_struct {
  int z;
} VoronoiWtStruct;

typedef struct voronoi_hook_struct {
  int shot;
  int coil;
  int sample;
  int onHull;
} VoronoiHookStruct;

/* GLOBAL VARIABLES FOR MASTER */
static char worker_hosts[512];
static double phase_delta_scale;
static int xvoxel_set= 0;
static int yvoxel_set= 0;

/* GLOBAL VARIABLES FOR WORKER */
double *image= NULL;
double *samples= NULL;
#ifdef USE_NFFT
static nfft_plan fwd_plan;
static infft_plan inv_plan;
static int nfft_initialized= 0;
#endif

/* GLOBAL VARIABLES FOR MASTER & WORKERS */
Task t;
Context c;
Result r;
MRI_Dataset* Input= NULL;
MRI_Dataset* Output= NULL;
static char* progname= NULL;

/* FORWARD DECLARATIONS */
static void InitMaster ();
static void ReadArguments (const int argc, const char **argv);
static void loadLagMap();
static void PrintMode ();
static void ProcessFile ();
static void LoadSamples();
static void GenerateImage ();
static void SaveImage();
static void CreateOutputDataset (int argc, char** argv);
static void MasterTask (const int argc, const char **argv, const char **envp);
static void MasterResult (int task_number);
static void WorkerContext ();
static void WorkerInitialize();
static void WorkerTask ();
static void WorkerFinalize();
static void PackContext();
static void UnpackContext();
static void PackTask();
static void UnpackTask();
static void PackResult();
static void UnpackResult();

static char* weightMethodName( WeightMethod mthd )
{
  switch (mthd) {
  case WEIGHT_CONST: return "const";
  case WEIGHT_VORONOI: return "voronoi";
  case WEIGHT_CIRCVORONOI: return "circ_voronoi";
  case WEIGHT_INVALID: return "invalid weighting!";
  default: return NULL;
  }
}

static char* nftMethodName( NftMethod mthd )
{
  switch (mthd) {
  case NFT_DIRECT: return "direct";
  case NFT_NFFT: return "nfft";
  case WEIGHT_INVALID: return "invalid NFT algorithm!";
  default: return NULL;
  }
}

static WeightMethod parseWeightMethod(char* s)
{
  if (!strcasecmp(s,"const")) return WEIGHT_CONST;
  else if (!strcasecmp(s,"voronoi")) return WEIGHT_VORONOI;
  else if (!strcasecmp(s,"circle_voronoi")) return WEIGHT_CIRCVORONOI;
  else return WEIGHT_INVALID;
}

static NftMethod parseNftMethod(char* s)
{
  if (!strcasecmp(s,"direct")) return NFT_DIRECT;
  else if (!strcasecmp(s,"nfft")) return NFT_NFFT;
  else return NFT_INVALID;
}

static int parseAlgString(char* search_string) 
{
  char* work;
  char* tok;
  char* tmp;

  work= strdup(search_string); /* make guaranteed-writable copy */
  tok= strtok_r(work, " ,+&:;", &tmp);
  while (tok != NULL) {

    if (!strncasecmp(tok,"weight=",7)) {
      if ((c.alg.weightMethod=parseWeightMethod(tok+7))==WEIGHT_INVALID)
	{ free(work); return 0; }
    }
    else if (!strncasecmp(tok,"nft=",4)) {
      if ((c.alg.nftMethod=parseNftMethod(tok+4))==NFT_INVALID)
	{ free(work); return 0; }
    }
    else { free(work); return 0; }
    
    tok= strtok_r(NULL, " ,+&:;",&tmp);
  }

  free(work);
  return 1;
}

static char* getAlgString()
{
  static char result[512];

  result[0]= '\0';

  snprintf(result,sizeof(result),"weight=%s,nft=%s",
	  weightMethodName(c.alg.weightMethod),
	  nftMethodName(c.alg.nftMethod));

  return result;
}

int
main (int argc,
      char **argv,
      char **envp)
{
  progname= argv[0];
  Acct(PROCESSING);
  par_process(argc, argv, envp,
	      MasterTask, MasterResult,
	      WorkerContext, WorkerTask,
	      WorkerFinalize,
	      PackContext, UnpackContext,
	      PackTask, UnpackTask,
	      PackResult, UnpackResult);
  PrintAcct(argv[0], 0);
  exit(0);
}

/* MASTER PROCEDURES */

void
MasterTask (const int argc,
	    const char **argv,
	    const char **envp)
{
  int i;
  int len;

  InitMaster ();
  ReadArguments (argc, argv);
  if (c.verbose)
    PrintMode ();

  ProcessFile(argc,argv);

  par_finish();
  /* Clean up files and memory */
  if (c.useLagMap) free(c.lagMap);
  if (c.verbose)
    Message("Done!!\n");
}

static void
InitMaster ()
{
  parseAlgString(DEFAULT_ALGORITHM);

  /* initially set the variably-sized arrays to be unallocated */
  c.sampLoc = NULL;
  c.sampWeight = NULL;
  c.lagMap= NULL;
  c.useLagMap= 0;
}

static int safe_get_extent(MRI_Dataset* ds, const char* chunk, char* dim)
{
  char buf[256];
  if (strlen(chunk)>(sizeof(buf)-10))
    Abort("%s: safe_get_extent: absurdly long chunk name!\n",progname);
  snprintf(buf,sizeof(buf)-1,"%s.extent.%c",chunk,*dim);
  buf[sizeof(buf)-1]= '\0';
  if (!mri_has(ds,buf))
    Abort("%s: file unexpectedly has no tag %s!\n",progname,buf);
  return mri_get_int(ds,buf);
}

static int map_valid(MRI_Dataset* imgDS, long* dx, long* dy, long* dz, 
		       long* dt)
{
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
  if (c.debug)
    fprintf(stderr,"Got dims %ld, %ld, %ld, %ld\n",*dx,*dy,*dz,*dt);
  return 1;
}

static void loadLagMap()
{
  long mapDx, mapDy, mapDz, mapDt;
  
  MRI_Dataset* ds= mri_open_dataset(c.lagmap_fname, MRI_READ);

  if (!map_valid(ds, &mapDx, &mapDy, &mapDz, &mapDt))
    Abort("%s: map file must have dimensions (v)xyz(t)!\n",progname);
  if (mapDx != c.dx || mapDy != c.dy || mapDz != c.dz)
    Abort("%s: map dimensions don't match output image!\n",progname);
  if (mapDt != 1)
    Warning(1,"Map contains multiple times; only the first will be used!\n");

  if (!(c.lagMap=(double*)malloc(c.dx*c.dy*c.dz*sizeof(double))))
    Abort("%s: unable to allocate %d bytes!\n",c.dx*c.dy*c.dz*sizeof(double));

  (void)mri_read_chunk(ds, "images", c.dx*c.dy*c.dz, 0, MRI_DOUBLE, c.lagMap);
  mri_close_dataset(ds);
}

static void
ReadArguments (const int argc, const char **argv)
{
  char alg_string[512];

  /* Check to see if help was requested */
  if( ( argc > 1 ) && !strcmp( argv[1], "-help" ) )
    {
      if( argc == 2 )
        Help( "selecttopic" );
      else
        Help( (char*)(argv[2]) );
    }

  /*** Parse command line ***/
  cl_scan( argc, (char**)argv );

  /* first handle the argument-less switches */
  c.verbose= cl_present("v|verbose");
  c.debug= cl_present("debug");

  cl_get("hosts", "%option %s[%]", "", worker_hosts);
  cl_get("alg|algorithm", "%option %s[%]", DEFAULT_ALGORITHM, alg_string);

  cl_get("xres|xrs", "%option %d[64]", &(c.dx));
  cl_get("yres|yrs", "%option %d[64]", &(c.dy));

  if (cl_get("xvoxel|xvx", "%option %lf", &(c.xvoxel)))
    xvoxel_set= 1;
  else xvoxel_set= 0;
  if (cl_get("yvoxel|yvx", "%option %lf", &(c.yvoxel)))
    yvoxel_set= 1;
  else yvoxel_set= 0;

  cl_get("phscale", "%option %lf[1.0]", &phase_delta_scale);
  cl_get("samplag", "%option %lf[0.0]", &(c.sample_lag));

  c.useLagMap= cl_get("lagmap", "%option %s", c.lagmap_fname);

  /* Strip possible quotes off worker_hosts; shell syntax makes
   * them impossible to avoid in some cases.
   */
  if (worker_hosts[0]=='\'' || worker_hosts[0]=='"') {
    worker_hosts[0]= ' ';
    worker_hosts[strlen(worker_hosts)-1]= '\0';
  }

  if (!(cl_get( "", "%s", &(c.input_fname) )))
    Abort("%s:Required input file name not given!\n",argv[0]);
  if (!(cl_get( "", "%s", &(c.output_fname) )))
    Abort("%s:Required output file name not given!\n",argv[0]);

  if (cl_cleanup_check()) {
    int i;
    fprintf(stderr,"%s: invalid argument in command line:\n    ", argv[0]);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(-1);
  }
  /*** End command-line parsing ***/

  vply_setDebug(c.debug);

  if (!parseAlgString(alg_string))
    Abort("%s: invalid algorithm <%s>!\n",argv[0],alg_string);

  if (worker_hosts[0] != '\0')
    par_set_hosts(worker_hosts);
}

static void
PrintMode ()
{
  char *s;
  int v;

  Message("\nProgram slow_ft starting\n");
  Message("\n");

  Message("Input: %s\n", c.input_fname);
  Message("Output: %s\n", c.output_fname);
  if (c.useLagMap) 
    Message("Lag map: %s\n", c.lagmap_fname);
  Message("\n");
  Message("Algorithm: %s\n", getAlgString());
	
  if ((s = getenv("PAR_ENABLE")) != NULL &&
      sscanf(s, "%d", &v) == 1 &&
      v != 0)
    {
      Message("Parallel Operation:\n");
      Message("\tEnabled\n");
      Message("\tHosts: %s\n", worker_hosts);
      Message("\n");
    }
  else if (worker_hosts[0] != '\0')
    {
      Message("Parallel Operation:\n");
      Message("\tDisabled since PAR_ENABLE not set to 1\n");
      Message("\tHosts: %s\n", worker_hosts);
      Message("\n");
    }

  Message("\n");
}

static void free5DDoubleArray( double***** array, int d1, int d2 )
{
  int i1;
  int i2;
  for (i1=0; i1<d1; i1++) {
    for (i2=0; i2<d2; i2++)
      Free3DDoubleArray(array[i1][i2]);
    free( array[i1] );
  }
  free(array);
}

static double***** alloc5DDoubleArray( int d1, int d2, int d3, int d4, int d5 )
{
  int i1;
  int i2;
  double***** result= NULL;

  /* Note that only the last 3 dims are guaranteed to be contiguous! */
  if (!(result= (double*****)malloc(d1*sizeof(double****))))
    Abort("%s: unable to allocate %d bytes!\n",d1*sizeof(double****));
  for (i1=0; i1<d1; i1++) {
    if (!(result[i1]= (double****)malloc(d2*sizeof(double***))))
      Abort("%s: unable to allocate %d bytes!\n",d2*sizeof(double***));
    for (i2=0; i2<d2; i2++) {
      result[i1][i2]= Alloc3DDoubleArray(d3, d4, d5);
    }
  }

  return result;
}

static void
ReadFileHeader (const Filename input_fname, Context *c)
{
  int islice;
  int icoil;
  int ishot;
  int isample;
  double nom_fov; /* nominal image space field of view in cm */
  double nom_res; /* nominal image space resolution */

  if (Input != NULL) {
    Acct(READCLOSE);
    mri_close_dataset(Input);
    Input= NULL;
  }

  Acct(READOPEN);
  /* open the file */
  if (!(Input= mri_open_dataset(input_fname, MRI_READ)))
    Abort("Can't open %s\n", input_fname);
  if (!mri_has(Input,"samples"))
    Abort("Input dataset %s has no 'samples' chunk!\n",input_fname);
  if (!mri_has(Input,"samples.dimensions"))
    Abort("Input dataset %s has no 'samples.dimensions' tag!\n",input_fname);
  if (strcmp(mri_get_string(Input,"samples.dimensions"),"vpsczt"))
    Abort("Input dataset %s sample dimensions are not vpsczt!\n",input_fname);
  if (!mri_has(Input,"sample_kxloc"))
    Abort("Input dataset %s has no 'sample_kxloc' chunk!\n",input_fname);
  if (!mri_has(Input,"sample_kxloc.dimensions"))
    Abort("Input dataset %s has no 'sample_kxloc.dimensions' tag!\n",
	  input_fname);
  if (strcmp(mri_get_string(Input,"sample_kxloc.dimensions"),"pscz"))
    Abort("Input dataset %s kxloc dimensions are not pscz!\n",input_fname);
  if (!mri_has(Input,"sample_kyloc"))
    Abort("Input dataset %s has no 'sample_kyloc' chunk!\n",input_fname);
  if (!mri_has(Input,"sample_kyloc.dimensions"))
    Abort("Input dataset %s has no 'sample_kyloc.dimensions' tag!\n",
	  input_fname);
  if (strcmp(mri_get_string(Input,"sample_kyloc.dimensions"),"pscz"))
    Abort("Input dataset %s kyloc dimensions are not pscz!\n",input_fname);
  
  /* extract the context info */
  c->dz= mri_get_int(Input,"samples.extent.z");
  c->dc= mri_get_int(Input,"samples.extent.c");
  c->dp= mri_get_int(Input,"samples.extent.p");
  c->ds = mri_get_int(Input,"samples.extent.s");
  c->dt = mri_get_int(Input,"samples.extent.t");

  /* We need the nominal FOV of samples to compare with reconstructed FOV */
  if (mri_has(Input,"samples.opfov")) {
    c->samp_nom_fov= 10.0*mri_get_float(Input, "samples.opfov");
  }
  else c->samp_nom_fov= 200.0;

  /* This sets the expected phase step per sample */
  c->ph_delta_per_sample= -0.5*M_PI; /* default to Nyquist limit freq */
  if (mri_has(Input,"samples.fast_rec_lpf") 
      && mri_has(Input,"samples.samp_time")) {
    double fast_rec_lpf= mri_get_float(Input,"samples.fast_rec_lpf");
    double samp_time= mri_get_float(Input,"samples.samp_time");
    double newVal= -2.0*M_PI*1.0e-3*fast_rec_lpf*samp_time;
    if (c->debug)
      fprintf(stderr,
	      "Calculated ph_delta of %f from header info, vs %f default\n",
	      newVal, c->ph_delta_per_sample);
    
    c->ph_delta_per_sample= newVal;
  }
  c->ph_delta_per_sample *= phase_delta_scale; /* scale from command line */

  /* Most of our files have voxel size info now */
  if (!xvoxel_set) {
    if (mri_has(Input,"samples.voxel_spacing.x"))
      c->xvoxel= mri_get_float(Input,"samples.voxel_spacing.x");
    else Abort("%s: x voxel size information is not available!",progname);
  }
  if (!yvoxel_set) {
    if (mri_has(Input,"samples.voxel_spacing.y"))
      c->xvoxel= mri_get_float(Input,"samples.voxel_spacing.y");
    else Abort("%s: y voxel size information is not available!",progname);
  }
  
  /* allocate the sampLoc and sampWeight arrays */
  c->sampLoc= alloc5DDoubleArray(c->dz, c->dc, 2, c->ds, c->dp);
  c->sampWeight= alloc5DDoubleArray(c->dz, c->dc, 2, c->ds, c->dp);

  /* Import trajectory info */
  for (islice=0; islice<c->dz; islice++)
    for (icoil=0; icoil<c->dc; icoil++)
      for (ishot=0; ishot<c->ds; ishot++) {
	long long offset = ((((islice*c->dc)+icoil)*c->ds)+ishot)*c->dp;
	mri_read_chunk(Input, "sample_kxloc", c->dp, offset,
		      MRI_DOUBLE, c->sampLoc[islice][icoil][0][ishot]);
	mri_read_chunk(Input, "sample_kyloc", c->dp, offset,
		      MRI_DOUBLE, c->sampLoc[islice][icoil][1][ishot]);
      }

  /* Initialize weight info */
  for (islice=0; islice<c->dz; islice++)
    for (icoil=0; icoil<c->dc; icoil++)
      for (ishot=0; ishot<c->ds; ishot++)
	for (isample=0; isample<c->dp; isample++) {
	  c->sampWeight[islice][icoil][0][ishot][isample]= 0.0;
	  c->sampWeight[islice][icoil][1][ishot][isample]= 0.0;
	}

  if (c->verbose) {
    Message("Dataset info: %d times, %d slices, %d shots/slice\n",
	   c->dt, c->dz, c->ds);
    Message("              %d coils, %d samples/coil\n", c->dc, c->dp);
  }
}


static float* access_coords( void* data, int i, void* info, void** hookPtr )
{
  static float result[2];
  VoronoiWtStruct* ws= (VoronoiWtStruct*)info;
  VoronoiHookStruct* hs= NULL;
  double***** sampLoc= (double*****)data;
  int ic= i / (c.ds*c.dp);
  int is= (i - ic*(c.ds*c.dp)) / c.dp;
  int ip= i - ((ic*c.ds+is)*c.dp);
  
  if (!(hs= (VoronoiHookStruct*)malloc(sizeof(VoronoiHookStruct))))
    Abort("%s: unable to allocate %d bytes!\n",sizeof(VoronoiHookStruct));

  /* Label this point with info about the sample */
  hs->coil= ic;
  hs->shot= is;
  hs->sample= ip;
  hs->onHull= 0; /* The calling program will eventually mark these */
  *hookPtr= hs;

  result[0]=sampLoc[ws->z][ic][0][is][ip];
  result[1]=sampLoc[ws->z][ic][1][is][ip];

  return result;
}

static void GenerateVoronoiWeights()
{
  dch_Tess* myTess= NULL;
  VoronoiWtStruct ws;
  dch_Pt_list* ptList;
  int iz;
  
  for (iz=0; iz<c.dz; iz++) {
    int ic;
    int is;
    VPoly* convexHull= NULL;
    
    ws.z= iz;
    
    if (c.verbose) {
      Message("Generating Voronoi polygon file polys.ps\n");
      vply_prepPS("polys.ps");
    }
    
    /* Build the Voronoi tesselation */
    if (c.verbose) 
      Message("Creating tess for slice %d\n",ws.z);
    myTess= dch_create_dirichlet_tess( (void*)c.sampLoc, 
				       c.dc*c.ds*c.dp, 2, &ws,
				       access_coords );
    
    /* Crawl around and find the convex hull */
    convexHull= vply_create( myTess->center_pt->coords[0],
			     myTess->center_pt->coords[1],
			     100 );
    ptList= myTess->pt_list;
    while (ptList) {
      dch_Pt* thisPt= ptList->pt;
      dch_Pt_list* neighbors= thisPt->neighbors;
      VoronoiHookStruct* hs= (VoronoiHookStruct*)(thisPt->user_hook);
      
      /* For each non-bogus point, we check to see if any of 
       * its neighbors are bogus.  If so, it's on the convex hull.
       */
      if (thisPt->user_hook != NULL) {
	while (neighbors) {
	  if (neighbors->pt->user_hook == NULL) {
	    vply_add(convexHull, thisPt->coords[0], thisPt->coords[1]);
	    hs->onHull= 1;
#ifdef never
	    fprintf(stderr,"Point %d at %f %f is on hull!\n",
		    thisPt->id,thisPt->coords[0],thisPt->coords[1]);
#endif
	    break;
	  }
	  neighbors= neighbors->next;
	}
      }
      else {
	/* Bogus point; ignore it */
      }
      ptList= ptList->next;
    }
    vply_sort(convexHull);
    if (c.verbose && c.alg.weightMethod==WEIGHT_VORONOI) 
      vply_writePS(convexHull);

    /* Assign weights to all the sample points (all shots and coils) 
     * based on the area of the associated Voronoi polygons.
     */
    ptList= myTess->pt_list;
    while (ptList) {
      double area;
      dch_Pt* thisPt= ptList->pt;
      dch_Vtx_list* vtxList;
      VPoly* vp= NULL;
      if (thisPt->user_hook) {
	/* This is a real point; get its sample info from where it was
	 * attached by access_coords()
	 */
	VoronoiHookStruct* hs= (VoronoiHookStruct*)(thisPt->user_hook);
	int coil= hs->coil;
	int shot= hs->shot;
	int sample= hs->sample;
	
#ifdef never
	fprintf(stderr,
		"Pt at %f %f is sample %d, shot %d, coil %d, dch_pt %d\n",
		thisPt->coords[0],thisPt->coords[1],sample,shot,coil,
		thisPt->id);
#endif

	/* Create a polygon containing all the vertices */
	vp= vply_create(thisPt->coords[0], thisPt->coords[1], 10);
	vtxList= thisPt->verts;
	while (vtxList) {
	  vply_add(vp, vtxList->vtx->coords[0], vtxList->vtx->coords[1]);
	  vtxList= vtxList->next;
	}
	vply_sort(vp);

	/* Clip the polygon; otherwise outside-edge polys will be huge */
	switch (c.alg.weightMethod) {
	case WEIGHT_VORONOI:
	  vp= vply_clip(vp, convexHull);
	  break;
	case WEIGHT_CIRCVORONOI:
	  vp= vply_clipToCircle(vp, 0.0, 0.0, 
				vply_getMaxDistance(convexHull));
	}
	if (c.verbose) vply_writePS(vp);
	area= vply_calcArea(vp);

	/* Assign the area of the polygon area as a point weight.  If
	 * there are multiple samples at this location, divide the 
	 * area equally between them.
	 */
	if (thisPt->twins) {
	  dch_Pt_list* twinList= thisPt->twins;
	  int nTwins= 0;
	  double wtPerTwin;
#ifdef never
	  fprintf(stderr,
		  "Pt %d (shot %d coil %d sample %d) at %f %f has twins!\n",
		  thisPt->id, shot, coil, sample, 
		  thisPt->coords[0],thisPt->coords[1]);
#endif
	  while (twinList) {
	    dch_Pt* twinPt= twinList->pt;
#ifdef never
	    VoronoiHookStruct* twinH= 
	      (VoronoiHookStruct*)(twinPt->user_hook);
	    int twinCoil= twinH->coil;
	    int twinShot= twinH->shot;
	    int twinSample= twinH->sample;
	    fprintf(stderr,"   twin pt %d (shot %d coil %d sample %d)\n",
		    twinPt->id, twinShot, twinCoil, twinSample);
#endif
	    nTwins++;
	    twinList= twinList->next;
	  }
	  wtPerTwin= area/(double)nTwins;
#ifdef never
	  fprintf(stderr,"%d twins in all; each gets %f\n",
		  nTwins,wtPerTwin);
#endif
	  c.sampWeight[iz][coil][0][shot][sample]= wtPerTwin;
	  c.sampWeight[iz][coil][1][shot][sample]= 0.0;
	  twinList= thisPt->twins;
	  while (twinList) {
	    dch_Pt* twinPt= twinList->pt;
	    VoronoiHookStruct* twinH= 
	      (VoronoiHookStruct*)(twinPt->user_hook);
	    int twinCoil= twinH->coil;
	    int twinShot= twinH->shot;
	    int twinSample= twinH->sample;
	    c.sampWeight[iz][twinCoil][0][twinShot][twinSample]= wtPerTwin;
	    c.sampWeight[iz][twinCoil][1][twinShot][twinSample]= 0.0;
	    /* User hook info is not managed by dch_destroy_tesselation() */
	    free(twinH); 
	    twinPt->user_hook= NULL;
	    twinList= twinList->next;
	  }
	}
	else {
	  c.sampWeight[iz][coil][0][shot][sample]= area;
	  c.sampWeight[iz][coil][1][shot][sample]= 0.0;
	}
	vply_destroy(vp);
	
	/* User hook info is not managed by dch_destroy_tesselation() */
	free(hs); 
	thisPt->user_hook= NULL;
      }
      else {
	/* This is a bogus point and should be ignored */
#ifdef never
	fprintf(stderr,"Point %d is a bogus point.\n",ptList->pt->id);
#endif
      }
      ptList= ptList->next;
    }
    dch_destroy_tesselation( myTess );
    if (c.verbose) vply_finishPS();
    vply_destroy(convexHull);
  }
}

static void GenerateWeights()
{

  if (c.alg.weightMethod==WEIGHT_VORONOI
      || c.alg.weightMethod==WEIGHT_CIRCVORONOI) {
    GenerateVoronoiWeights();
  }
  else if (c.alg.weightMethod==WEIGHT_CONST) {
    int iz;
    int ic;
    int is;
    int ip;
    for (iz=0; iz<c.dz; iz++)
      for (ic=0; ic<c.dc; ic++)
	for (is=0; is<c.ds; is++)
	  for (ip=0; ip<c.dp; ip++) {
	    c.sampWeight[iz][ic][0][is][ip]= 1.0;
	    c.sampWeight[iz][ic][1][is][ip]= 0.0;
	  }
  }
  else Abort("%s: internal error: unknown weight method %d!\n",
	     progname,(int)c.alg.weightMethod);
}

static void
ProcessFile (int argc, char** argv)
{
  ReadFileHeader(c.input_fname, &c);
  if (c.useLagMap) loadLagMap();
  GenerateWeights();
  CreateOutputDataset(argc, argv);

  /* Close input so that workers can re-open it */
  if (Input != NULL) {
    Acct(READCLOSE);
    mri_close_dataset(Input);
    Input= NULL;
  }

  par_set_context();

  for (t.t = 0; t.t < c.dt; t.t++) {
    for (t.z = 0; t.z < c.dz; t.z++) {
      par_delegate_task();
    }
  }
}

static void
CreateOutputDataset (int argc, char** argv)
{
  Filename fn;
  char* s;
  char buf[256];
  char* val;
  int i;

  if (Output != NULL) { /* well, that would certainly be a surprise... */
    Acct(WRITECLOSE);
    mri_close_dataset(Output);
    Output= NULL;
  }

  Acct(WRITEOPEN);
  Output = mri_open_dataset(c.output_fname, MRI_WRITE);

  /* Copy in the history, and add the current command */
  for (i=1; hist_get( Input, i ); i++) {
    char keybuf[64];
    snprintf(keybuf,sizeof(keybuf),"history.%d",i);
    mri_set_string(Output, keybuf, hist_get( Input, i ));
  }
  hist_add_cl(Output,argc,argv);

  Acct(WRITING);
  mri_create_chunk(Output, "images");
  mri_set_string(Output, "images.datatype", "float32");
  mri_set_string(Output, "images.dimensions", "vxyzt");
  mri_set_int(Output, "images.extent.v", 2);
  mri_set_string(Output, "images.description.v", "complex real/imaginary");
  mri_set_int(Output, "images.extent.x", c.dx);
  mri_set_string(Output, "images.description.x", "gridded image-space");
  mri_set_int(Output, "images.extent.y", c.dy);
  mri_set_string(Output, "images.description.y", "gridded image-space");
  mri_set_int(Output, "images.extent.z", c.dz);
  mri_set_string(Output, "images.description.z", "gridded image-space");
  mri_set_int(Output, "images.extent.t", c.dt);
  mri_set_string(Output, "images.description.t", "gridded image-space");
  mri_set_string(Output, "images.file", ".dat");

  /* Additional info available from Pfile header */
  mri_set_float( Output, "images.voxel_size.x",c.xvoxel);
  mri_set_float( Output, "images.voxel_size.y",c.yvoxel);
  
  /* Copy "samples.*" tags from input, except those which specifically 
   * deal with the chunk data.  Many are useful scan info.
   */
  mri_iterate_over_keys(Input);
  while ((s=mri_next_key(Input)) != NULL) {
    if (!strncmp(s,"samples.",8)
	&& strcmp(s,"samples.datatype")
	&& strcmp(s,"samples.dimensions")
	&& strcmp(s,"samples.file")
	&& strcmp(s,"samples.order")
	&& strcmp(s,"samples.offset")
	&& strcmp(s,"samples.size")
	&& strcmp(s,"samples.little_endian")
	&& strncmp(s,"samples.extent.",strlen("samples.extent."))
	&& strcmp(val=mri_get_string(Input,s),"[chunk]")) {
      if (snprintf(buf,sizeof(buf),"images.%s",s+strlen("samples."))<0)
	Abort("%s: key <%s> is too long!\n",s);
      mri_set_string(Output,buf,val);
    }
  }

  Acct(WRITECLOSE);
  mri_close_dataset(Output);
  Output= NULL;
}

static void LoadSamples()
{
  long long blocksize= 2*c.dp*c.ds*c.dc;
  long long offset= blocksize*(t.t*c.dz + t.z);
  if (!(mri_read_chunk(Input, "samples", blocksize, offset, 
		       MRI_DOUBLE, samples))) 
    Abort("Unable to read %lld samples starting at %lld from %s!\n",
	  blocksize, offset, c.input_fname);
}

static void SaveImage()
{
  long long blocksize= 2*c.dx*c.dy;
  long long offset= blocksize*(t.t*c.dz + t.z);
  mri_set_chunk(Output, "images", blocksize, offset, MRI_DOUBLE,image);
}

static void GenerateImageDirect()
{
  int ploop;
  int cloop;
  int sloop;
  int i;
  int j;
  double kXScale= (c.dx*c.xvoxel)/c.samp_nom_fov;
  double kYScale= (c.dy*c.yvoxel)/c.samp_nom_fov;

  if (c.debug)
    fprintf(stderr,"kXScale= %f, kYScale= %f\n",kXScale,kYScale);

  for (i=0; i<c.dx; i++)
    for (j=0; j<c.dy; j++) {
      double* rHere= image+2*(i+(j*c.dx));
      double* iHere= rHere+1;
      *rHere= *iHere= 0.0;
    }

  for (cloop=0; cloop<c.dc; cloop++)
    for (sloop=0; sloop<c.ds; sloop++)
      for (ploop=0; ploop<c.dp; ploop++) {
	double sampX= samples[2*(((cloop*c.ds + sloop)*c.dp)+ploop)];
	double sampY= samples[2*(((cloop*c.ds + sloop)*c.dp)+ploop)+1];
	double sampMag;
	double sampPhase;
#ifdef never
	double kX= kXScale*c.sampLoc[t.z][cloop][0][sloop][ploop];
	double kY= kYScale*c.sampLoc[t.z][cloop][1][sloop][ploop];
	double sampWtX= c.sampWeight[t.z][cloop][0][sloop][ploop];
	double sampWtY= c.sampWeight[t.z][cloop][1][sloop][ploop];
#endif
	double kX;
	double kY;
	double sampWtX;
	double sampWtY;
	double sample_lag_here= c.sample_lag*((double)ploop/(double)(c.dp-1));

	if (ploop>0){
	  kX= LINTERP( kXScale*c.sampLoc[t.z][cloop][0][sloop][ploop],
		       kXScale*c.sampLoc[t.z][cloop][0][sloop][ploop-1],
		       sample_lag_here );
	  kY= LINTERP( kYScale*c.sampLoc[t.z][cloop][1][sloop][ploop],
		       kYScale*c.sampLoc[t.z][cloop][1][sloop][ploop-1],
		       sample_lag_here );
	  sampWtX= LINTERP( c.sampWeight[t.z][cloop][0][sloop][ploop],
			    c.sampWeight[t.z][cloop][0][sloop][ploop-1],
			    sample_lag_here );
	  sampWtY= LINTERP( c.sampWeight[t.z][cloop][1][sloop][ploop],
			    c.sampWeight[t.z][cloop][1][sloop][ploop-1],
			    sample_lag_here );
	}
	else {
	  kX= kXScale*c.sampLoc[t.z][cloop][0][sloop][ploop];
	  kY= kYScale*c.sampLoc[t.z][cloop][1][sloop][ploop];
	  sampWtX= c.sampWeight[t.z][cloop][0][sloop][ploop];
	  sampWtY= c.sampWeight[t.z][cloop][1][sloop][ploop];
	}

	/* weight the sample, possibly changing its phase */
	sampX= sampWtX*sampX - sampWtY*sampY;
	sampY= sampWtX*sampY + sampWtY*sampX;

	sampMag= sqrt(sampX*sampX+sampY*sampY);
	sampPhase= atan2(sampY,sampX);
	if (c.debug && sampMag>0.0) 
	  fprintf(stderr,"%d: mag %f, phase %f at %f %f; k's %f %f\n",
		  ploop,sampMag,sampPhase,sampX,sampY,kX,kY);
	if (c.useLagMap) {
	  for (i=0; i<c.dx; i++)
	    for (j=0; j<c.dy; j++) {
	      double* rHere= image+2*(i+(j*c.dx));
	      double* iHere= rHere+1;
	      double mapval= c.lagMap[(((t.z*c.dy) + j)*c.dx) + i];
	      double mag= sampMag;
	      double phase= sampPhase 
		- 2*M_PI*(i-(c.dx/2))*kX/c.dx
		- 2*M_PI*(j-(c.dy/2))*kY/c.dy
		+ c.ph_delta_per_sample
		*LINTERP(ploop,(ploop-1),sample_lag_here*mapval);
	      *rHere += mag*cos(phase);
	      *iHere += mag*sin(phase);
	    }
	}
	else {
	  sampPhase += c.ph_delta_per_sample
	    *LINTERP(ploop,(ploop-1),sample_lag_here);
	  for (i=0; i<c.dx; i++)
	    for (j=0; j<c.dy; j++) {
	      double* rHere= image+2*(i+(j*c.dx));
	      double* iHere= rHere+1;
	      double mag= sampMag;
	      double phase= sampPhase 
		- 2*M_PI*(i-(c.dx/2))*kX/c.dx
		- 2*M_PI*(j-(c.dy/2))*kY/c.dy;
	      *rHere += mag*cos(phase);
	      *iHere += mag*sin(phase);
	    }
	}
      }
}

#ifdef USE_NFFT
static void nfftSetSampleKspaceLocations()
{
  long ploop;
  long cloop;
  long sloop;
#ifdef never
  double kXScale= (c.dx*c.xvoxel)/c.samp_nom_fov;
  double kYScale= (c.dy*c.yvoxel)/c.samp_nom_fov;
#endif
  double kXScale= 1.0/c.samp_nom_fov;
  double kYScale= 1.0/c.samp_nom_fov;
  kXScale*=M_PI;
  kYScale*=M_PI;

  if (c.debug)
    fprintf(stderr,"kXScale= %f, kYScale= %f\n",kXScale,kYScale);

  for (cloop=0; cloop<c.dc; cloop++)
    for (sloop=0; sloop<c.ds; sloop++)
      for (ploop=0; ploop<c.dp; ploop++) {
	double kX;
	double kY;
	double sample_lag_here= c.sample_lag*((double)ploop/(double)(c.dp-1));
	long offset= (cloop*c.ds + sloop)*c.dp + ploop;

	if (ploop>0){
	  kX= LINTERP( kXScale*c.sampLoc[t.z][cloop][0][sloop][ploop],
		       kXScale*c.sampLoc[t.z][cloop][0][sloop][ploop-1],
		       sample_lag_here );
	  kY= LINTERP( kYScale*c.sampLoc[t.z][cloop][1][sloop][ploop],
		       kYScale*c.sampLoc[t.z][cloop][1][sloop][ploop-1],
		       sample_lag_here );
	}
	else {
	  kX= kXScale*c.sampLoc[t.z][cloop][0][sloop][ploop];
	  kY= kYScale*c.sampLoc[t.z][cloop][1][sloop][ploop];
	}
	fwd_plan.x[2*offset]= kX;
	fwd_plan.x[2*offset+1]= kY;
      }
}

static void nfftSetSampleWeights()
{
  long ploop;
  long cloop;
  long sloop;

  for (cloop=0; cloop<c.dc; cloop++)
    for (sloop=0; sloop<c.ds; sloop++)
      for (ploop=0; ploop<c.dp; ploop++) {
	double sampWtX;
	double sampWtY;
	double sample_lag_here= c.sample_lag*((double)ploop/(double)(c.dp-1));
	long offset= (cloop*c.ds + sloop)*c.dp + ploop;

	if (ploop>0){
	  sampWtX= LINTERP( c.sampWeight[t.z][cloop][0][sloop][ploop],
			    c.sampWeight[t.z][cloop][0][sloop][ploop-1],
			    sample_lag_here );
	  sampWtY= LINTERP( c.sampWeight[t.z][cloop][1][sloop][ploop],
			    c.sampWeight[t.z][cloop][1][sloop][ploop-1],
			    sample_lag_here );
	}
	else {
	  sampWtX= c.sampWeight[t.z][cloop][0][sloop][ploop];
	  sampWtY= c.sampWeight[t.z][cloop][1][sloop][ploop];
	}
	inv_plan.w[offset]= sqrt(sampWtX*sampWtX + sampWtY*sampWtY);
      }
}

static void nfftSetSampleDampingFactors()
{
  long k;

  /* There are as many damping factors as points in the reconstructed array */
  for (k=0; k<fwd_plan.N_L; k++) inv_plan.w_hat[k]= 1.0;
}

static void InitializeNFFT()
{
  int my_N[2],my_n[2];
   
  /* initialise fwd_plan */
  my_N[0]=c.dx; my_n[0]=next_power_of_2(c.dx);
  my_N[1]=c.dy; my_n[1]=next_power_of_2(c.dy);

  if (c.useLagMap) 
    Abort("%s: lagmap is not supported when using NFFT reconstruction!\n",
	  progname);

  if (nfft_initialized) {
    infft_finalize(&inv_plan);
    nfft_finalize(&fwd_plan);
    nfft_initialized= 0;
  }

  nfft_init_specific(&fwd_plan, 2, my_N, c.dp*c.ds*c.dc, my_n, 6, 
		     PRE_PHI_HUT| PRE_PSI| PRE_FULL_PSI| MALLOC_X| MALLOC_F_HAT|
		     MALLOC_F| FFTW_INIT| FFT_OUT_OF_PLACE,
		     FFTW_MEASURE| FFTW_DESTROY_INPUT);
#ifdef never
  nfft_init_2d(&fwd_plan,c.dx,c.dy,c.dp*c.ds*c.dc);
#endif
  infft_init_specific(&inv_plan,&fwd_plan,PRECOMPUTE_WEIGHT|CGNR_E);
  nfft_initialized= 1;

  nfftSetSampleKspaceLocations();

  /* precompute psi */
  if (fwd_plan.nfft_flags & PRE_PSI) {
    nfft_precompute_psi(&fwd_plan);
    if (fwd_plan.nfft_flags & PRE_FULL_PSI)
      nfft_full_psi(&fwd_plan,DLAMCH("e"));
  }

  /* initialize weights and damping factors */
  if (inv_plan.infft_flags & PRECOMPUTE_WEIGHT) {
    nfftSetSampleWeights();
  }
  if (inv_plan.infft_flags & PRECOMPUTE_DAMP) {
    nfftSetSampleDampingFactors();
  }
}

static void FinalizeNFFT()
{
  infft_finalize(&inv_plan);
  nfft_finalize(&fwd_plan);
  nfft_initialized= 0;
}

static void GenerateImageNFFT()
{
  long ploop;
  long cloop;
  long sloop;
  long i;
  long j;
  long k;
  long l;
  double oldQuality;

  /* Initialize the 'guess' to zero */
  for (k=0; k<fwd_plan.N_L; k++) {
    inv_plan.f_hat_iter[k][0]= 0.0;
    inv_plan.f_hat_iter[k][1]= 0.0;
  }

  /* Copy in the sample info */
  for (cloop=0; cloop<c.dc; cloop++)
    for (sloop=0; sloop<c.ds; sloop++)
      for (ploop=0; ploop<c.dp; ploop++) {
	double sampX= samples[2*(((cloop*c.ds + sloop)*c.dp)+ploop)];
	double sampY= samples[2*(((cloop*c.ds + sloop)*c.dp)+ploop)+1];
	double sampMag;
	double sampPhase;
	double sample_lag_here= c.sample_lag*((double)ploop/(double)(c.dp-1));
	long offset= (cloop*c.ds + sloop)*c.dp + ploop;

	/* With this algorithm we don't directly weight the sample.  The
	 * weights embedded in the nfft data structures are pure real, so
	 * there is no mechanism to change the sample phase on the fly.
	 */
	sampMag= sqrt(sampX*sampX+sampY*sampY);
	sampPhase= atan2(sampY,sampX);
	sampPhase += c.ph_delta_per_sample
	  *LINTERP(ploop,(ploop-1),sample_lag_here);
	sampX= sampMag*cos(sampPhase);
	sampY= sampMag*sin(sampPhase);
	inv_plan.given_f[offset][0]= sampX;
	inv_plan.given_f[offset][1]= sampY;
      }

  /* Do the solution */
  infft_before_loop(&inv_plan);
#ifdef never
  for (l=0; l<10;l++) {
    fprintf(stderr,"l= %ld: %e\n",l,sqrt(inv_plan.dot_r_iter));
    infft_loop_one_step(&inv_plan);
  }
#endif
  oldQuality= sqrt(inv_plan.dot_r_iter);
  l= 0;
  while (1) {
    double newQuality;
    double ratio;
    infft_loop_one_step(&inv_plan);
    newQuality= sqrt(inv_plan.dot_r_iter);
    ratio= fabs((oldQuality-newQuality)/newQuality);
    fprintf(stderr,"l= %ld: %e %e %e\n",l,oldQuality,newQuality,ratio);
    if (ratio<=0.001) break;
    oldQuality= newQuality;
    l++;
  }
  
  /* Copy out the result */
  for (i=0; i<c.dx; i++)
    for (j=0; j<c.dy; j++) {
      double* rHere= image+2*(j+(i*c.dy));
      double* iHere= rHere+1;
      long offset= (j*c.dx)+i;
      *rHere= inv_plan.f_hat_iter[offset][0];
      *iHere= inv_plan.f_hat_iter[offset][1];
    }

}
#endif

static void GenerateImage()
{
  switch (c.alg.nftMethod) {
  case NFT_DIRECT: 
    {
      GenerateImageDirect();
    }
    break;
#ifdef USE_NFFT
  case NFT_NFFT: 
    {
      GenerateImageNFFT();
    }
    break;
#endif
  default:
    Abort("Unimplemented NFT algorithm %d=<%s>!\n",
	  (int)c.alg.nftMethod,nftMethodName(c.alg.nftMethod));
  }
}

static void
MasterResult (int task_number)
{
  static int results= 0;

  ++results;
  if (!(results % 60) || (results == c.dz*c.dt))
    Message("# %ld\n", (long)results);
  else
    Message("#");
}

static void
WorkerInitialize()
{
  if (Input==NULL) {
    if (!(Input= mri_open_dataset(c.input_fname, MRI_READ)))
      Abort("Can't open %s for reading\n", c.input_fname);
  }
  if (Output==NULL) {
    if (!(Output= mri_open_dataset(c.output_fname, MRI_MODIFY_DATA)))
      Abort("Can't open %s for writing\n", c.output_fname);
  }
  if (!(image=(double*)malloc(2*c.dx*c.dy*sizeof(double))))
    Abort("Can't allocate %d doubles!\n", c.dx*c.dy);
  if (!(samples=(double*)malloc(2*c.dp*c.ds*c.dc*sizeof(double))))
    Abort("Can't allocate %d doubles!\n", 2*c.dp*c.ds*c.dc);

  switch (c.alg.nftMethod) {
  case NFT_DIRECT:
    {
      /* Nothing to do */
    }
    break;
#ifdef USE_NFFT
  case NFT_NFFT:
    {
      InitializeNFFT();
    }
    break;
#endif
  default:
    Abort("Unimplemented NFT algorithm %d=<%s>!\n",
	  (int)c.alg.nftMethod,nftMethodName(c.alg.nftMethod));
  }
}

static void
WorkerContext ()
{
  WorkerInitialize();
}

static void WorkerFinalize()
{
  switch (c.alg.nftMethod) {
  case NFT_DIRECT:
    {
      /* Nothing to do */
    }
    break;
#ifdef USE_NFFT
  case NFT_NFFT:
    {
      FinalizeNFFT();
    }
    break;
#endif
  default:
    Abort("Unimplemented NFT algorithm %d=<%s>!\n",
	  (int)c.alg.nftMethod,nftMethodName(c.alg.nftMethod));
  }

  if (Input != NULL) mri_close_dataset(Input);
  if (Output != NULL) mri_close_dataset(Output);
  if (image != NULL) free(image);
  if (samples != NULL) free(samples);
  if (c.lagMap != NULL) free(c.lagMap);
}

static void
WorkerTask ()
{
  LoadSamples();
  GenerateImage();
  SaveImage();

  /* Let the master know what we've been up to... */
  r.z= t.z;
  r.t= t.t;
}

static void
PackContext ()
{
  int iz;
  int ic;
  par_pkstr(c.input_fname);
  par_pkstr(c.output_fname);
  par_pkint(c.dp);
  par_pkint(c.ds);
  par_pkint(c.dc);
  par_pkint(c.dz);
  par_pkint(c.dt);
  par_pkint(c.dx);
  par_pkint(c.dy);
  par_pkdouble(c.xvoxel);
  par_pkdouble(c.yvoxel);
  par_pkdouble(c.samp_nom_fov);
  par_pkdouble(c.ph_delta_per_sample);
  par_pkdouble(c.sample_lag);
  par_pkint(c.verbose);
  par_pkint(c.debug);
  par_pkint(c.useLagMap);
  par_pkint((int)c.alg.weightMethod);
  par_pkint((int)c.alg.nftMethod);
  for (iz=0; iz<c.dz; iz++)
    for (ic=0; ic<c.dc; ic++) 
      par_pkdoublearray(&c.sampLoc[iz][ic][0][0][0], 2*c.ds*c.dp);
  for (iz=0; iz<c.dz; iz++)
    for (ic=0; ic<c.dc; ic++) 
      par_pkdoublearray(&c.sampWeight[iz][ic][0][0][0], 2*c.ds*c.dp);
  if (c.useLagMap)
    par_pkdoublearray(c.lagMap, c.dx*c.dy*c.dz);
  
}

static void
UnpackContext ()
{
  int iz;
  int ic;
  static int previously_allocated = FALSE;

  par_upkstr(c.input_fname);
  par_upkstr(c.output_fname);
  c.dp= par_upkint();
  c.ds= par_upkint();
  c.dc= par_upkint();
  c.dz= par_upkint();
  c.dt= par_upkint();
  c.dx= par_upkint();
  c.dy= par_upkint();
  c.xvoxel= par_upkdouble();
  c.yvoxel= par_upkdouble();
  c.samp_nom_fov= par_upkdouble();
  c.ph_delta_per_sample= par_upkdouble();
  c.sample_lag= par_upkdouble();  
  c.verbose= par_upkint();
  c.debug= par_upkint();
  c.useLagMap= par_upkint();
  c.alg.weightMethod= (WeightMethod)par_upkint();
  c.alg.nftMethod= (WeightMethod)par_upkint();

  /* deallocate old sampLoc array */
  if (previously_allocated)
    {
      free5DDoubleArray(c.sampLoc, c.dz, c.dc);
      free5DDoubleArray(c.sampWeight,c.dz, c.dc);
      previously_allocated = FALSE;
    }

  /* allocate the sampLoc and kdens arrays */
  c.sampLoc = alloc5DDoubleArray(c.dz, c.dc, 2, c.ds, c.dp);
  c.sampWeight = alloc5DDoubleArray(c.dz, c.dc, 2, c.ds, c.dp);
  previously_allocated = TRUE;

  for (iz=0; iz<c.dz; iz++)
    for (ic=0; ic<c.dc; ic++) 
      par_upkdoublearray(&c.sampLoc[iz][ic][0][0][0], 2*c.ds*c.dp);
  for (iz=0; iz<c.dz; iz++)
    for (ic=0; ic<c.dc; ic++) 
      par_upkdoublearray(&c.sampWeight[iz][ic][0][0][0], 2*c.ds*c.dp);
  if (c.useLagMap) {
    if (c.lagMap) free(c.lagMap);
    if (!(c.lagMap=(double*)malloc(c.dx*c.dy*c.dz*sizeof(double))))
      Abort("%s: unable to allocate %d bytes!\n",
	    c.dx*c.dy*c.dz*sizeof(double));
    par_upkdoublearray(c.lagMap, c.dx*c.dy*c.dz);
  }
}

static void
PackTask ()
{
  par_pkint(t.z);
  par_pkint(t.t);
}

static void
UnpackTask ()
{
  t.z= par_upkint();
  t.t= par_upkint();
}

static void
PackResult ()
{
  par_pkint(r.z);
  par_pkint(r.t);
}

static void
UnpackResult ()
{
  r.z= par_upkint();
  r.t= par_upkint();
}


