/*
 *	srecon.c	Version 1.2
 *
 *    Image Reconstruction for Spiral k-space Acquisitions
 *                           (c) 1992,1993,1994,1995,1996
 *       Copyright by Douglas C. Noll and the University of Pittsburgh and
 *	   the Pittsburgh Supercomputing Center
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
 *         Algorithm is based on O'Sullivan, IEEE Trans. Med. Imag., 
 *         4(4):200-207, 1985 as refined by Jackson, et al., Trans. Med. Imag., 
 *         10(2):473-478, 1991.  Sample density correction as described by
 *         Meyer, et al., Magn. Reson. in Med., 28:202-213, 1992. Time-
 *         segmented inhomogeneity correction as described by
 *         Noll, et al., IEEE Trans. Med. Imag., 10(4):629-637, 1991.
 *
 *  John Pauly of Stanford University wrote code that served as the
 *   basis for this program.  Many changes occurred since then.
 *  Several utility routines came from Numerical Recipes in C,
 *   namely bessi0 and four1.
 *  Versions:
 *     gpr2   - converted to spin-echo proj.
 *     gpr3   - window = 2, added apodization correction
 *     gpr6   - added other convolution windows as well
 *              option for sample density correction
 *     gpr8   - added k-b window
 *     gsp1   - 1st run at spiral recon
 *     gsp2   - fixed set up for multi-slices
 *     gsp3   - added multi-phase, reordered input parm list
 *     gsp4   - added over flag (=2 for 2x over, =1 for no over)
 *            - works, but not as well as Craig's version
 *	      - fixed delay term on readout
 *     gsp4hc - 1st run at homogeneity correction (with field map)
 *     gsp4hc2 - uses a previously generated field map
 *     gsp4hc3 - added fixviews subroutine for phase corrections
 *     gsp4hc4 - added section to read in list of files and uncompress
 *     gsp7  - put shifting, rotating, etc. into gridding subroutine
 *           - brought map making into same program (-m) option
 *           - defaut is no correction, (-h) option will do correction
 *           - compressed input is made an option (-c)
 *     gsp9  - 5.x version
 *     gsp10 - added parameters to rhuser variables in header
 *     gsp11 - made subsequent times through faster than first 
 *     gsp12 - added phase multiplier
 *           - turned off most output unless using the (-v) option
 *           - added image registration (-r)
 *     gsp14 - added multi-coil recon option
 *           - added (-2) option to double translations with (-r)
 *           - added (-o) option to specify starting image number
 *           - added (-m2,-h2) options to do linear term inhomogeneity corr.
 *           - added (-m3,-h3) options for both linear and time-seg corrs.
 *           - added (-l) option for shift correction in obliques
 *    srecon - split reconstruction code off from gsp14 and converted to run
 *	     - in parallel using PVM  (Greg Hood, PSC)
 *	v1.0 - bug fixes; changed to read only one sampim file per slice
 *	     -	(Greg Hood, PSC)
 *	v1.1 - now uses the updated parallel task manager in util/par.c
 *	     -  (Greg Hood, PSC, June 96)
 *	v1.2 - output is now in Pittsburgh File Format
 *	     -  (Greg Hood, PSC, Oct 96)
 *	v1.3 - uses Pittsburgh File Format for raw input
 *	     -  (Greg Hood, PSC, Feb 97)
 */

/*
 * TO DO:
 *	The FFT routine used here could be replaced with highly optimized
 *	    routines on certain architectures.  These would have to be
 *	    incorporated using #ifdefs so that the code remains portable.
 *	A few more globals could possibly be turned into local variables.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "stdcrg.h"
#include "par.h"
#include "bio.h"
#include "array.h"
#include "misc.h"
#include "acct.h"
#include "mri.h"
#include "fmri.h"

#define PI	M_PI

#define SWAP(a,b) {tempr=(a);(a)=(b);(b)=tempr;}
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

typedef struct Context {
  /* NOTE: any new fields added to this struct should
     be also added to PackContext and UnpackContext */ 
  Filename input_directory;	/* where to find input files */
  int big_endian_input;		/* TRUE if input is big-endian */
  Filename tmp_directory;	/* where to place scratch files */
  Filename output_directory;	/* where to place output files */
  Filename output_name;		/* name to give to output dataset */
  Filename output_data;		/* file in which to place output data chunk */
  int big_endian_output;	/* TRUE if output if big-endian */

  int ncoils;			/* # of coils */
  int nslices;			/* # of slices per coil */
  int nimages;			/* # of images (phases) per slice */
  int ndat;			/* # of data items per projection in
				   original spiral data */

  int res;			/* # of pixels along the X and Y dimensions */
  int os_res;			/* # of pixels in the (possibly oversampled)
				   input image ( =res*over_samp) */

  int slice;			/* slice number to work on (0 indicates all slices) */

  float ph_twist;		/* phase twist */

  int im_offset;		/* offset used in output image numbers */

  int start_slice;		/* starting slice number */
  int end_slice;		/* ending slice number */

  int hc_cor;			/* time-segmented correction */

  int lin_map;			/* linear correction maps */
  int gen_map;			/* general correction maps */

  float map_time;
  float samp_time;
  int filter_sz;		/* total width of filter in pixels */

  int over_samp;		/* oversampling ratio */
  float grid_len;		/* grid length (half width of convolution) */
  float gridb;

  int all_coils;		/* TRUE if we should write a separate image
				   for each coil */
  int wr_phase;			/* TRUE if we should write an image phase file */

  int slice_num;		/* slice number (first is 0) */

  float **refim_in;		/* [res][res]  input reference image for this slice */
  float **sampim;		/* [os_res][os_res]  sample info for this slice */
} Context;

typedef struct Task {
  /* NOTE: any new fields added to this struct should
     be also added to PackTask and UnpackTask */ 
  int image_num;	/* image number within a slice (first is 0) (corresponds to phase number in gsp14.c) */
} Task;

static char rcsid[] = "$Id: srecon.c,v 1.9 2005/07/07 20:02:56 welling Exp $";

/* GLOBAL VARIABLES FOR MASTER & SLAVE */
Task t;
Context c;

/* GLOBAL VARIABLES FOR MASTER */
MRI_Dataset *mrds = NULL;	/*			master's raw dataset */

/* GLOBAL VARIABLES FOR SLAVE */
float ***grim = NULL;		/* [2][os_res][os_res]  gridded image data */
float ***im = NULL;		/* [2][os_res][os_res]  image data */
short int **outim = NULL;	/* [res][res]		single coil output image data */
short int **totalim = NULL;	/* [res][res]		combined coil output image data */
short int **imphase = NULL;	/* [res][res]		image phase */
float *wc = NULL;		/* [res]		weighting coefficients */
float **refim_out = NULL;	/* [res][res]		reference image */
float **refim_mag = NULL;	/* [res][res]		reference image magnitude */	
float **fbuf = NULL;		/* [2][os_res] 		FIX what? buffer */
MRI_Dataset *ods = NULL;	/* 			output datset */
MRI_Dataset **cds = NULL;	/* [ncoils]		coil-specific datasets */
MRI_Dataset *srds = NULL;	/* 			slave's raw dataset */


/* FORWARD DECLARATIONS */
void MasterTask (const int argc, const char **argv, const char **envp);
void ReadEnvironment ();
void ReadArguments (const int argc, const char **argv);
void ReadRawHeader ();
void CreateOutputDataset ();
void LoadHomogeneityCorrectionReferenceData (int slice_num);
void SlaveFinalize();
void SlaveContext ();
void SlaveTask ();
void ComputeReferences ();
void ComputeReference (int slice_num, int coil_num, int image_num);
void MakePhase1Filter ();
void MakePhase2Filter ();
void WriteCoilReference (int slice_num, int coil_num);
void ComputeMulticoilReference (int slice_num);
void WriteReference (int slice_num);
void ComputeImages ();
short *CombineCoilImage (int *pmax, int coil_num, int ncoils);
void LoadRawImage (int slice_num, int coil_num, int image_num);
void TransformImage (int slice_num, int coil_num, int image_num);
void BesselKaiserCorrectImage (int *pmax);
void WindowImageData (int middle_sample, int sample_window);
void CopyImageData ();
void ApplyFourierTransform ();
void MakeFinalImage (float ***finalim, float ***im, int segment_num);
void CopyFinalImage (float ***im, float ***finalim);
void WriteSingleCoilImage (int coil_num);
void WriteMulticoilImage (short *totalImage);
void WriteImagePhase (int slice_num, int coil_num, int image_num);
void fft (float *rdat, float *idat, int nn, int isign);
void PackContext();
void UnpackContext();
void PackTask();
void UnpackTask();


int
main (int argc,
      char **argv,
      char **envp)
{
  Acct(PROCESSING);
  par_process(argc, argv, envp,
	      MasterTask, NULL,
	      SlaveContext, SlaveTask,
	      SlaveFinalize,
	      PackContext, UnpackContext,
	      PackTask, UnpackTask,
	      NULL, NULL);
  PrintAcct(argv[0], 0);
  exit(0);
}

/* MASTER PROCEDURES */

void
MasterTask (const int argc,
	    const char **argv,
	    const char **envp)
{
  ReadEnvironment();
  ReadArguments(argc, argv);
  ReadRawHeader();
  CreateOutputDataset();

  /* go through all slices */
  for (c.slice_num = c.start_slice; c.slice_num <= c.end_slice; ++c.slice_num)
    {
      if (c.hc_cor)
	LoadHomogeneityCorrectionReferenceData(c.slice_num);
      par_set_context();

      if (c.gen_map || c.lin_map)
	/* if we are generating reference files, assign an entire slice
	   to each slave; that slave will iterate through all image sets
	   and coils within the slice; this is necessary because the images
	   within a slice must be processed sequentially */
	par_delegate_task();
      else
	/* if we are generating image files, assign each multi-coil image set
	   to a slave; that slave will iterate through all coils */
	for (t.image_num = 0; t.image_num < c.nimages; ++t.image_num)
	  par_delegate_task();
    }

  par_finish();

  /* delete the raw dataset that we just processed */
  mri_destroy_dataset(mrds);

  Report("Done!!\n");
}

void
ReadEnvironment ()
{
  char *fn;
  char *hosts;
  char *dir;
  char *name;

  c.filter_sz = GetEnvInt("F_SRECON_FILTER_SZ");
  c.hc_cor = GetEnvInt("F_SRECON_HC_COR");
  c.lin_map = GetEnvInt("F_SRECON_LIN_MAP");
  c.gen_map = GetEnvInt("F_SRECON_GEN_MAP");
  c.im_offset = GetEnvInt("F_SRECON_IM_OFFSET");
  c.slice = GetEnvInt("F_SRECON_SLICE");
  c.all_coils = GetEnvInt("F_SRECON_ALL_COILS");
  c.wr_phase = GetEnvInt("F_SRECON_WR_PHASE");
  c.map_time = GetEnvFloat("F_SRECON_MAP_TIME");
  c.samp_time = GetEnvFloat("F_SRECON_SAMP_TIME");
  c.ph_twist = GetEnvFloat("F_SRECON_PH_TWIST");

  if ((hosts = getenv("F_SRECON_HOSTS")) != NULL)
    par_set_hosts(hosts);

  strcpy(c.input_directory, ".");
  if ((dir = getenv("F_SRECON_INPUT_DIR")) != NULL)
    StringCopy(c.input_directory, dir, sizeof(Filename));

  strcpy(c.tmp_directory, "/tmp");
  if ((dir = getenv("F_SRECON_TMP_DIR")) != NULL)
    StringCopy(c.tmp_directory, dir, sizeof(Filename));

  strcpy(c.output_directory, ".");
  if ((dir = getenv("F_SRECON_OUT_DIR")) != NULL)
    StringCopy(c.output_directory, dir, sizeof(Filename));

  strcpy(c.output_name, "recon.mri");
  if ((name = getenv("F_SRECON_OUT_NAME")) != NULL)
    StringCopy(c.output_name, name, sizeof(Filename));

  strcpy(c.output_data, ".dat");
  if ((name = getenv("F_SRECON_OUT_DATA")) != NULL)
    StringCopy(c.output_data, name, sizeof(Filename));

  /* initially set the variably-sized arrays to be
     unallocated */
  c.refim_in = NULL;
  c.sampim = NULL;

  /* the endianness of the input and output defaults
     to the native mode of the machine on which we are
     running */
  c.big_endian_input = bio_big_endian_machine;
  c.big_endian_output = bio_big_endian_machine;
}

void
ReadArguments (const int argc,
	       const char **argv)
{
  int a;

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

  if (cl_present("a")) c.all_coils= TRUE; /* -2 flag */
  if (cl_present("bei")) c.big_endian_input= TRUE;
  if (cl_present("beo")) c.big_endian_output= TRUE;
  if (cl_present("lei")) c.big_endian_input= FALSE;
  if (cl_present("leo")) c.big_endian_output= FALSE;
  cl_get("o","%option %d[%]", 0, &(c.im_offset));
  if (cl_present("p")) c.wr_phase= TRUE;
  if (cl_present("v")) verbose= TRUE;

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

void ReadRawHeader ()
{
  Filename mrds_name;

  sprintf(mrds_name, "%s/raw.mri", c.input_directory);
  Acct(READOPEN);
  mrds = mri_open_dataset(mrds_name, MRI_MODIFY);
  Acct(PROCESSING);
  c.ncoils = mri_get_int(mrds, "images.extent.c");
  c.nslices = mri_get_int(mrds, "images.extent.z");
  c.nimages = mri_get_int(mrds, "images.extent.t");

  c.ndat = mri_get_int(mrds, "ndat");
  c.res = mri_get_int(mrds, "res");
  c.over_samp = mri_get_int(mrds, "over_samp");
  c.os_res = c.res * c.over_samp;
  c.grid_len = mri_get_int(mrds, "grid_len");
  c.gridb = PI*c.grid_len;

  /* set the starting and ending slices */
  if (c.slice > 0)
    c.start_slice = c.end_slice = c.slice - 1;
  else
    {
      c.start_slice = 0;
      c.end_slice = c.nslices - 1;
    }

  if (c.hc_cor)
    {
      c.refim_in = Alloc2DFloatArray(c.res, c.res);
      c.sampim = Alloc2DFloatArray(c.os_res, c.os_res);
    }
}

void
CreateOutputDataset ()
{
  Filename fn;
  MRI_Dataset *ds;
  int i;

  Acct(WRITEOPEN);
  sprintf(fn, "%s/%s", c.output_directory, c.output_name);
  ds = mri_open_dataset(fn, MRI_WRITE);
  hist_add(ds,"*srecon*");

  Acct(WRITING);
  mri_create_chunk(ds, "images");
  mri_set_string(ds, "images.datatype", "int16");
  mri_set_string(ds, "images.dimensions", "xyzt");
  mri_set_int(ds, "images.extent.x", c.res);
  mri_set_int(ds, "images.extent.y", c.res);
  mri_set_int(ds, "images.extent.z", c.end_slice - c.start_slice + 1);
  mri_set_int(ds, "images.extent.t", c.nimages);
  if (c.start_slice > 0)
    mri_set_int(ds, "images.min.z", c.start_slice);

  mri_set_string(ds, "images.file", c.output_data);
  Acct(WRITECLOSE);
  mri_close_dataset(ds);

  Acct(PROCESSING);
  /* create the coil-specific datasets */
  if (c.all_coils)
    for (i = 0; i < c.ncoils; ++i)
      {
	sprintf(fn, "%s/%s.c%d", c.output_directory, c.output_name, i);
	Acct(WRITEOPEN);
	ds = mri_open_dataset(fn, MRI_WRITE);
	hist_add(ds,"*srecon coil*");

	Acct(WRITING);
	mri_create_chunk(ds, "images");
	mri_set_string(ds, "images.datatype", "int16");
	mri_set_string(ds, "images.dimensions", "xyzt");
	mri_set_int(ds, "images.extent.x", c.res);
	mri_set_int(ds, "images.extent.y", c.res);
	mri_set_int(ds, "images.extent.z", c.end_slice - c.start_slice + 1);
	mri_set_int(ds, "images.extent.t", c.nimages);
	if (c.start_slice > 0)
	  mri_set_int(ds, "images.min.z", c.start_slice);
	
	if (c.output_data[0] == '.' || c.output_data[0] == '\0')
	  mri_set_string(ds, "images.file", c.output_data);
	else
	  mri_set_string(ds, "images.file", ".dat");
	Acct(WRITECLOSE);
	mri_close_dataset(ds);
	Acct(PROCESSING);
      }
}

void
LoadHomogeneityCorrectionReferenceData (int slice_num)
{
  Filename fn;
  FILE *f;
  int overall_offset;

  Acct(READOPEN);
  sprintf(fn, "%s/ref.r%d.s%.2d", c.input_directory, c.res, slice_num+1);
  if ((f = fopen(fn, "r")) == NULL)
    Abort("Can't open reference file %s!\n", fn);
  Acct(READING);
  bio_error = FALSE;
  bio_big_endian_input = c.big_endian_input;
  FRdFloat32Array(f, &c.refim_in[0][0], c.res*c.res);
  if (bio_error)
    Abort("Can't read reference file %s\n", fn);
  Acct(READCLOSE);
  fclose(f);

  Acct(READING);
  overall_offset = slice_num * c.os_res * c.os_res;
  memcpy(&c.sampim[0][0],
	 mri_get_chunk(mrds, "sampim", c.os_res*c.os_res,
		       overall_offset, MRI_FLOAT),
	 c.os_res*c.os_res*sizeof(float));
  Acct(PROCESSING);
}


/* SLAVE PROCEDURES */

void
SlaveFinalize()
{
  if (srds != NULL)
    {
      Acct(READCLOSE);
      mri_close_dataset(srds);
    }
  if (ods != NULL)
    {
      Acct(WRITECLOSE);
      mri_close_dataset(ods);
    }
  /* deallocate the slave arrays */
  Free3DFloatArray(grim);
  Free3DFloatArray(im);
  Free2DShortArray(outim);
  Free2DShortArray(totalim);
  Free2DShortArray(imphase);
  if (wc != NULL)
    free(wc);
  Free2DFloatArray(refim_out);
  Free2DFloatArray(refim_mag);
  Free2DFloatArray(fbuf);
}

void
SlaveContext ()
{
  static int old_res = 0, old_os_res = 0;
  int hr;
  float wctmp, wcmax;
  Filename fn;
  int i;

  if (c.res != old_res ||
      c.os_res != old_os_res)
    {
      /* (re)open the input dataset */
      if (srds != NULL)
	{
	  Acct(READCLOSE);
	  mri_close_dataset(srds);
	}
      Acct(READOPEN);
      sprintf(fn, "%s/raw.mri", c.input_directory);
      srds = mri_open_dataset(fn, MRI_READ);

      /* (re)open the output dataset */
      if (ods != NULL)
	{
	  Acct(WRITECLOSE);
	  mri_close_dataset(ods);
	}
      Acct(WRITEOPEN);
      sprintf(fn, "%s/%s", c.output_directory, c.output_name);
      ods = mri_open_dataset(fn, MRI_MODIFY_DATA);

      /* (re)open the coil-specific datasets */
      if (c.all_coils)
	{
	  if (cds != NULL)
	    {
	      Acct(WRITECLOSE);
	      for (i = 0; i < c.ncoils; ++i)
		if (cds[i] != NULL)
		  mri_close_dataset(cds[i]);
	      free(cds);
	    }
	  Acct(WRITEOPEN);
	  cds = (MRI_Dataset **) malloc(c.ncoils * sizeof(MRI_Dataset *));
	  for (i = 0; i < c.ncoils; ++i)
	    {
	      sprintf(fn, "%s/%s.c%d", c.output_directory, c.output_name, i);
	      cds[i] = mri_open_dataset(fn, MRI_MODIFY_DATA);
	    }
	}
      Acct(PROCESSING);

      /* reallocate the slave arrays */
      Free3DFloatArray(grim);
      Free3DFloatArray(im);
      Free2DShortArray(outim);
      Free2DShortArray(totalim);
      Free2DShortArray(imphase);
      if (wc != NULL)
	free(wc);
      Free2DFloatArray(refim_out);
      Free2DFloatArray(refim_mag);
      Free2DFloatArray(fbuf);

      grim = Alloc3DFloatArray(2, c.os_res, c.os_res);
      im = Alloc3DFloatArray(2, c.os_res, c.os_res);
      outim = Alloc2DShortArray(c.res, c.res);
      totalim = Alloc2DShortArray(c.res, c.res);
      imphase = Alloc2DShortArray(c.res, c.res);
      wc = (float *) malloc(c.res*sizeof(float));
      refim_out = Alloc2DFloatArray(c.res, c.res);
      refim_mag = Alloc2DFloatArray(c.res, c.res);
      fbuf = Alloc2DFloatArray(2, c.os_res);

      /* precompute the wc array that is used for kaiser-bessel correction */
      hr = c.res / 2;
      wcmax = sinh(sqrt(c.gridb*c.gridb))/sqrt(c.gridb*c.gridb);
      for (i = 0; i < c.res; i++)
	{
	  wctmp = PI*PI*c.grid_len*c.grid_len*(i-hr)*(i-hr)/(c.os_res*c.os_res/4) - c.gridb*c.gridb;
	  if(wctmp == 0.0)
	    wc[i] = wcmax;
	  else if (wctmp < 0.0)
	    wc[i] = sqrt(-wctmp)/sinh(sqrt(-wctmp))*wcmax;
	  else
	    wc[i] = sqrt(wctmp)/sin(sqrt(wctmp))*wcmax;
	}

      old_res = c.res;
      old_os_res = c.os_res;
    }
}

void
SlaveTask ()
{
  if (c.lin_map || c.gen_map)
    ComputeReferences();
  else
    ComputeImages();
}

void
ComputeReferences ()
{
  int coil_num;
  int image_num;

  for (coil_num = 0; coil_num < c.ncoils; ++coil_num)
    for (image_num = 0; image_num < c.nimages; ++image_num)
      ComputeReference(c.slice_num, coil_num, image_num);
}

void
ComputeReference (int slice_num,
		  int coil_num,
		  int image_num)
{
  LoadRawImage(slice_num, coil_num, image_num);
  TransformImage(slice_num, coil_num, image_num);
  
  if (image_num == 0)
    {
      /* we only make a filter for the first image, and do not write
	 a reference file */
      MakePhase1Filter();
      Report("\n");
      return;
    }
  else
    MakePhase2Filter();
  Report("Wr.");
  if (c.ncoils == 1)
    WriteReference(slice_num);
  else
    {
      WriteCoilReference(slice_num, coil_num);
      if (coil_num == c.ncoils - 1)
	{
	  /* if this is the last coil, then combine all of the
	     individual coil references into one */
	  ComputeMulticoilReference(slice_num);
	  WriteReference(slice_num);
	}
    }
  Report("\n");
}

void
MakePhase1Filter ()
{
  int i, j, k, offs;
  int hfsize;	/* half-width of filter-size */
  float sumr,sumi;
  int j_start, j_end;
  int k_start, k_end;

  hfsize = (c.filter_sz - 1)/2;
  offs = (c.over_samp-1)*c.res/2;

  for (i = 0; i < c.res; i++)
  {
    j_start = MAX(0, offs - hfsize);
    j_end = MIN(c.res + offs + hfsize, c.os_res-1);
    for (j = j_start; j <= j_end; j++)
      {
	fbuf[0][j] = 0.0;
	fbuf[1][j] = 0.0;
	k_start = MAX(0, i + offs - hfsize);
	k_end = MIN(i + offs + hfsize, c.os_res-1);
	for (k = k_start; k <= k_end; k++)
	  {
	    fbuf[0][j] += im[0][k][j];
	    fbuf[1][j] += im[1][k][j];
	  }
      }
    for (j = 0; j < c.res; j++)
      {
	sumr = sumi = 0.;
	k_start = MAX(0, j + offs - hfsize);
	k_end = MIN(j + offs + hfsize, c.os_res-1);
	for (k = k_start; k <= k_end; k++)
	  {
	    sumr += fbuf[0][k];
	    sumi += fbuf[1][k];
	  }
	refim_out[i][j] = atan2(sumi,sumr);
      }
  }
}

void
MakePhase2Filter ()
{
  int i, j, k, offs;
  int hfsize;
  float sumr, sumi;
  float phaseoff;
  int j_start, j_end;
  int k_start, k_end;

  /* phaseoff is value to offset freq map if phase twist (ph_twist) != 0 */
  phaseoff = -2.0*PI * (c.map_time / c.samp_time) * c.ph_twist / c.ndat;
  hfsize = (c.filter_sz - 1) / 2;
  offs = (c.over_samp - 1) * c.res/2;

  for (i = 0; i < c.res; i++)
    {
      j_start = MAX(0, offs - hfsize);
      j_end = MIN(c.res + offs + hfsize, c.os_res-1);
      for (j = j_start; j <= j_end; j++)
	{
	  fbuf[0][j] = 0.0;
	  fbuf[1][j] = 0.0;
	  k_start = MAX(0, i + offs - hfsize);
	  k_end = MIN(i + offs + hfsize, c.os_res-1);
	  for (k = k_start; k <= k_end; k++)
	    {
	      fbuf[0][j] += im[0][k][j];
	      fbuf[1][j] += im[1][k][j];
	    }
	}
      for (j = 0; j < c.res; j++)
	{
	  sumr = sumi = 0.0;
	  k_start = MAX(0, j + offs - hfsize);
	  k_end = MIN(j + offs + hfsize, c.os_res - 1);
	  for (k = k_start; k <= k_end; k++)
	    {
	      sumr += fbuf[0][k];
	      sumi += fbuf[1][k];
	    }
	  refim_mag[i][j] = hypot(sumi,sumr);
	  refim_out[i][j] -= (atan2(sumi,sumr) - phaseoff); 
	}
    }
}

void
WriteCoilReference (int slice_num,
		    int coil_num)
			 
{
  Filename fn;
  FILE *f;

  Acct(WRITEOPEN);
  sprintf(fn, "%s/ref.r%d.s%.2d.c%d", c.tmp_directory, c.res, slice_num+1, coil_num+1);
  if ((f = fopen(fn,"w")) == NULL)
    Abort("Can't open %s!\n", fn);
  Acct(WRITING);
  bio_big_endian_output = c.big_endian_output;
  FWrFloat32Array(f, &refim_out[0][0], c.res*c.res);
  Acct(WRITECLOSE);
  fclose(f);

  Acct(WRITEOPEN);
  sprintf(fn, "%s/refm.r%d.s%.2d.c%d", c.tmp_directory, c.res, slice_num+1, coil_num+1);
  if ((f = fopen(fn,"w")) == NULL)
    Abort("Can't open %s!\n", fn);
  Acct(WRITING);
  bio_big_endian_output = c.big_endian_output;
  FWrFloat32Array(f, &refim_mag[0][0], c.res*c.res);
  Acct(WRITECLOSE);
  fclose(f);
  Acct(PROCESSING);
}

void
ComputeMulticoilReference (int slice_num)
{
  int i, j, k;
  Filename fn;
  FILE *f;
  float **maptmp;
  float **maptmpmag;

  maptmp = Alloc2DFloatArray(c.res, c.res);
  maptmpmag = Alloc2DFloatArray(c.res, c.res);

  for (i = 0; i < c.res; i++)
    for (j = 0; j < c.res; j++)
      {
	refim_out[i][j] = 0.0;
	refim_mag[i][j] = 0.0;
      }

  for (k = 0; k < c.ncoils; k++)
    {
      Acct(READOPEN);
      sprintf(fn, "%s/ref.r%d.s%.2d.c%d", c.tmp_directory, c.res, slice_num+1, k+1);
      if ((f = fopen(fn, "r")) == NULL)
	Abort("Can't open reference file %s!\n", fn);
      Acct(READING);
      bio_error = FALSE;
      bio_big_endian_input = c.big_endian_input;
      FRdFloat32Array(f, &maptmp[0][0], c.res*c.res);
      if (bio_error)
	Abort("Can't read reference file %s\n", fn);
      Acct(READCLOSE);
      fclose(f);
      unlink(fn);

      Acct(READOPEN);
      sprintf(fn, "%s/refm.r%d.s%.2d.c%d", c.tmp_directory, c.res, slice_num+1, k+1);
      if ((f = fopen(fn, "r")) == NULL)
	Abort("Can't open reference file %s!\n", fn);
      Acct(READING);
      bio_error = FALSE;
      bio_big_endian_input = c.big_endian_input;
      FRdFloat32Array(f, &maptmpmag[0][0], c.res*c.res);
      if (bio_error)
	Abort("Can't read reference file %s\n", fn);
      Acct(READCLOSE);
      fclose(f);
      unlink(fn);
      Acct(PROCESSING);

      for (i = 0; i < c.res; ++i)
	for (j = 0; j < c.res; ++j)
	  {
	    refim_out[i][j] += maptmp[i][j] * maptmpmag[i][j];
	    refim_mag[i][j] += maptmpmag[i][j];
	  }
    }

  for (i = 0; i < c.res; i++)
    for (j = 0; j < c.res; j++)
      refim_out[i][j] /= (refim_mag[i][j] + 0.001);

  Free2DFloatArray(maptmp);
  Free2DFloatArray(maptmpmag);
}

void
WriteReference (int slice_num)
{
  int i, j;
  Filename fn;
  FILE *f;
  int hr;	/* half the resolution */
  int rm1;	/* resolution minus 1 */	
  float max;
  double tmpxr, tmpyr, tmpxi, tmpyi, tmpc, tmpcr, tmpci;
  float out[3];

  if (c.lin_map)
    {
      hr = c.res / 2;
      rm1 = c.res - 1;

      max = 0.0;
      for (i = 0; i < c.res; i++)
	for (j = 0; j < c.res; j++) 
	  if (hypot((float)(i-hr), (float)(j-hr)) < hr &&
	      refim_mag[i][j] > max)
	    max = refim_mag[i][j];
      max /= 4.0;

      tmpxr = tmpyr = tmpxi = tmpyi = tmpcr = tmpci = 0.0;
      /* find linear terms */
      for (i = 0; i < rm1; i++) 
      for (j = 0; j < rm1; j++) 
	if (refim_mag[i][j] > max)
	  {
	    tmpxr += refim_mag[i][j] * refim_mag[i][j+1]
	             * cos(refim_out[i][j]-refim_out[i][j+1]);
	    tmpxi += refim_mag[i][j] * refim_mag[i][j+1]
	             * sin(refim_out[i][j]-refim_out[i][j+1]);
	    tmpyr += refim_mag[i][j] * refim_mag[i+1][j]
		     * cos(refim_out[i][j]-refim_out[i+1][j]);
	    tmpyi += refim_mag[i][j] * refim_mag[i+1][j]
		     * sin(refim_out[i][j]-refim_out[i+1][j]);
        }

      tmpxr = -atan2(tmpxi,tmpxr);
      tmpyr = -atan2(tmpyi,tmpyr);
      /* uncomment to zero linear terms */
      /* tmpxr = tmpyr = 0.; */ 

      /* find constant term with linears removed */
      for (i = 0; i < rm1; i++) 
	for (j = 0; j < rm1; j++) 
	  if (refim_mag[i][j] > max)
	    {
	      tmpcr += refim_mag[i][j]
		       * cos(refim_out[i][j] - tmpxr*(j-hr) - tmpyr*(i-hr));
	      tmpci += refim_mag[i][j]
		       * sin(refim_out[i][j] - tmpxr*(j-hr) - tmpyr*(i-hr));
	    }

      tmpc = atan2(tmpci,tmpcr);

      out[0] = tmpc * c.samp_time / c.map_time;
      out[1] = tmpxr * c.res * c.samp_time / c.map_time / 2.0 /PI;
      out[2] = tmpyr * c.res * c.samp_time / c.map_time / 2.0 /PI;

      Acct(WRITEOPEN);
      sprintf(fn, "%s/refl.s%.2d", c.output_directory, slice_num+1);
      if ((f = fopen(fn, "w")) == NULL)
	Abort("Can't open %s!\n",fn);
      Acct(WRITING);
      bio_big_endian_output = c.big_endian_output;
      FWrFloat32Array(f, out, 3);
      Acct(WRITECLOSE);
      fclose(f);
      Acct(PROCESSING);

      for (i = 0; i < rm1; i++) 
	for (j = 0; j < rm1; j++)
	  {
	    tmpcr = refim_mag[i][j]
	            * cos(refim_out[i][j] - tmpxr*(j-hr) - tmpyr*(i-hr) - tmpc);
	    tmpci = refim_mag[i][j]
	            * sin(refim_out[i][j] - tmpxr*(j-hr) - tmpyr*(i-hr) - tmpc);
	    refim_out[i][j] = atan2(tmpci,tmpcr);
	  }
    }

  if (c.gen_map)
    {
      Acct(WRITEOPEN);
      sprintf(fn,"%s/ref.r%d.s%.2d", c.output_directory, c.res, slice_num+1);
      if ((f = fopen(fn, "w")) == NULL)
	Abort("Can't open %s!\n",fn);

      Acct(PROCESSING);
      for (i = 0; i < c.res; i++)
	for (j = 0; j < c.res; j++)
	  {
	    if (refim_out[i][j] > PI)
	      refim_out[i][j] -= 2*PI;
	    if (refim_out[i][j] < -PI)
	      refim_out[i][j] += 2*PI;
	    if (refim_out[i][j] > 0.5*PI)
	      refim_out[i][j] = 0.5*PI;
	    if (refim_out[i][j] < -0.5*PI)
	      refim_out[i][j] = -0.5*PI;
	  }
      Acct(WRITING);
      bio_big_endian_output = c.big_endian_output;
      FWrFloat32Array(f, &refim_out[0][0], c.res*c.res);
      Acct(WRITECLOSE);
      fclose(f);
      Acct(PROCESSING);
    }
}

void
ComputeImages ()
{
  int coil_num;
  int max;
  short *totalImage= NULL;

  for (coil_num = 0; coil_num < c.ncoils; ++coil_num)
    {
      LoadRawImage (c.slice_num, coil_num, t.image_num);
      TransformImage (c.slice_num, coil_num, t.image_num);
      BesselKaiserCorrectImage (&max);
      if (c.all_coils)
	WriteSingleCoilImage(coil_num);
      if (c.wr_phase)
	WriteImagePhase(c.slice_num, coil_num, t.image_num);
      totalImage = CombineCoilImage(&max, coil_num, c.ncoils);
      Report("max=%3d", max);
    }
  WriteMulticoilImage(totalImage);
}

void
LoadRawImage (int slice_num,
	      int coil_num,
	      int image_num)
{
  int overall_offset;

  Report("Ld.");

  /* load the grim array */
  Acct(READING);
  overall_offset = ((coil_num * c.nimages + image_num) * c.nslices + slice_num) * 2 * c.os_res * c.os_res;
  memcpy(&grim[0][0][0],
	 mri_get_chunk(srds, "images", 2*c.os_res*c.os_res,
		       overall_offset, MRI_FLOAT),
	 2*c.os_res*c.os_res*sizeof(float));
  Acct(PROCESSING);
}

void
TransformImage (int slice_num,
		int coil_num,
		int image_num)
{
  int nsegs;
  int segment_num;
  int sample_window;
  int middle_sample;
  float ***finalim;

  Report("FT.");
  if (c.hc_cor)
    {
      /* segmented hc recon */
      nsegs = floor(1.0 * c.ndat * c.samp_time / c.map_time);
      sample_window = c.map_time / c.samp_time;
      finalim = Alloc3DFloatArray(2, c.os_res, c.os_res);
      for (segment_num = 0; segment_num <= nsegs; segment_num++)
	{
	  Report("%d.", segment_num+1);
	  middle_sample = segment_num * sample_window;
	  WindowImageData(middle_sample, sample_window);
	  ApplyFourierTransform();
	  MakeFinalImage(finalim, im, segment_num);
	}
      CopyFinalImage(im, finalim);
      Free3DFloatArray(finalim);
    }
  else
    {
      CopyImageData();
      ApplyFourierTransform();
    }
}

void
BesselKaiserCorrectImage (int *pmax)
{
  int i, j;
  int hr;
  int offs;
  float imx, imi, imr, imm;
  Filename fn;
  FILE *f;

  hr = c.res / 2;
  imx = 0.0;
  offs = (c.over_samp-1)*hr;

  for (i = 0; i < c.res; i++)
    for (j = 0; j < c.res; j++)
      if (hypot((float)(i-hr),(float)(j-hr)) < .51*c.res)
	{
	  imr = im[0][-j+c.res+offs-1][-i+c.res+offs-1]*wc[i]*wc[j];
	  imi = im[1][-j+c.res+offs-1][-i+c.res+offs-1]*wc[i]*wc[j];
	  imm = sqrt(imi*imi+imr*imr);
	  if (imm > imx)
	    imx = imm;
	  outim[i][j] = 1024 * imm / (c.os_res*c.os_res);
	}
      else
	outim[i][j] = 0.0;
  *pmax = (int) (imx*1024/(c.os_res*c.os_res));
}

void
WindowImageData (int middle_sample,
		 int sample_window)
{
  float start_sample, end_sample;
  float radius;
  float ww;
  int i, j;

  start_sample = middle_sample - sample_window;
  end_sample = middle_sample + sample_window;

  for (i = 0; i < c.os_res; ++i)
    for (j = 0; j < c.os_res; ++j)
      {
	radius = c.sampim[i][j];
	if (radius > start_sample && radius < end_sample)
	  ww = 0.5 * (1.0 + cos(PI*(radius - middle_sample) / sample_window));
	else
	  ww = 0.0;
	im[0][i][j] = ww*grim[0][i][j];
	im[1][i][j] = ww*grim[1][i][j];
      }
}

void
CopyImageData ()
{
  memcpy(&im[0][0][0], &grim[0][0][0], 2*c.os_res*c.os_res*sizeof(float));
}

void
ApplyFourierTransform ()         /* formerly ft_image */
{
  int i, j;
  int offs, offs2;
  float **w1;			/* [2][os_res]   	temporary array in which the FFT is done  */

  w1 = Alloc2DFloatArray(2, c.os_res);

  offs = c.os_res/2;
  offs2 = (c.over_samp-1)*c.res/2;
  for (i = 0; i < c.os_res; i++)
    {
      for (j = 0; j < offs; j++)
	{
	  w1[0][j] = im[0][i][j+offs];
	  w1[1][j] = im[1][i][j+offs];
	  w1[0][j+offs] = im[0][i][j];
	  w1[1][j+offs] = im[1][i][j];
	}
      fft(w1[0], w1[1], c.os_res, 1);
      for (j = 0; j < offs; j++)
	{
	  im[0][i][j+offs] = w1[0][j];
	  im[1][i][j+offs] = w1[1][j];
	  im[0][i][j] = w1[0][j+offs];
	  im[1][i][j] = w1[1][j+offs];
	}
    }

  for (j = 0; j < c.res; j++)
    {
      for (i = 0; i < offs; i++)
	{
	  w1[0][i] = im[0][i+offs][j+offs2];
	  w1[1][i] = im[1][i+offs][j+offs2];
	  w1[0][i+offs] = im[0][i][j+offs2];
	  w1[1][i+offs] = im[1][i][j+offs2];
	}
      fft(w1[0], w1[1], c.os_res, 1);
      for (i=0; i<offs; i++)
	{
	  im[0][i+offs][j+offs2] = w1[0][i];
	  im[1][i+offs][j+offs2] = w1[1][i];
	  im[0][i][j+offs2] = w1[0][i+offs];
	  im[1][i][j+offs2] = w1[1][i+offs];
	}
    }

  Free2DFloatArray(w1);
}

void
MakeFinalImage (float ***finalim,
		float ***im,
		int segment_num)
{
  int i, j;
  int offs;
  float phase_factor;

  phase_factor = -segment_num;
  offs = (c.over_samp-1)*c.res/2;
  if (segment_num == 0)
    for (i=0; i<c.res; i++)
      {
	memcpy(&finalim[0][i][0],
	       &im[0][i+offs][offs],
	       c.res * sizeof(float));
	memcpy(&finalim[1][i][0],
	       &im[1][i+offs][offs],
	       c.res * sizeof(float));
      }
  else
    for (i = 0; i < c.res; i++)
      for (j = 0; j < c.res; j++)
	{
	  finalim[0][i][j] += 
	    im[0][i+offs][j+offs]*cos(phase_factor*c.refim_in[i][j]) -
	    im[1][i+offs][j+offs]*sin(phase_factor*c.refim_in[i][j]);
	  finalim[1][i][j] += 
	    im[0][i+offs][j+offs]*sin(phase_factor*c.refim_in[i][j]) +
	    im[1][i+offs][j+offs]*cos(phase_factor*c.refim_in[i][j]);
	}
}

void
CopyFinalImage (float ***im, float ***finalim)
{
  int i;
  int offs;

  offs = (c.over_samp-1)*c.res/2;
  for (i = 0; i < c.res; i++)
    {
      memcpy(&im[0][i+offs][offs],
	     &finalim[0][i][0],
	     c.res * sizeof(float));
      memcpy(&im[1][i+offs][offs],
	     &finalim[1][i][0],
	     c.res * sizeof(float));
    }
}

void WriteSingleCoilImage (int coil_num)
{
  Acct(WRITING);
  Report("Wr.");
  mri_set_image(cds[coil_num], t.image_num, c.slice_num - c.start_slice, MRI_SHORT, &outim[0][0]);
  Acct(PROCESSING);
}

void WriteMulticoilImage (short *totalImage)
{
  Acct(WRITING);
  Report("Wr.");
  mri_set_image(ods, t.image_num, c.slice_num - c.start_slice, MRI_SHORT, totalImage);
  Acct(PROCESSING);
  Report("\n");
}

short *CombineCoilImage (int *pmax, int coil_num, int ncoils)
{
  float imt, imc, imm, imx;
  int i, j;
  int hr;

  /* if there is just 1 coil, we just pass back
     a pointer to the outim array to save the
     cost of doing a copy into the totalim array */
  if (ncoils == 1)
    return(&outim[0][0]);

  if (coil_num == 0)
    {
      /* copy the coil image into the total image */
      for (i = 0; i < c.res; ++i)
	for (j = 0; j < c.res; ++j)
	  totalim[i][j] = outim[i][j];
      return(&totalim[0][0]);
    }

  /* combine the image in outim with that in totalim
     and place the result in totalim */
  hr = c.res / 2;
  imx = 0.0;
  for (i = 0; i < c.res; ++i)
    for (j = 0; j < c.res; ++j)
      if (hypot((float)(i-hr),(float)(j-hr)) < .51*c.res)
	{
	  imt = totalim[i][j];
	  imt *= (c.os_res*c.os_res)/1024.0;
	  imc = outim[i][j];
	  imc *= (c.os_res*c.os_res)/1024.0;
	  imm = sqrt(imt*imt + imc*imc);
	  if (imm > imx)
	    imx = imm;
	  totalim[i][j] = 1024*imm/(c.os_res*c.os_res);
	}
  *pmax = (int) (imx * 1024 / (c.os_res * c.os_res));
  return(&totalim[0][0]);
}

void
WriteImagePhase (int slice_num,
		 int coil_num,
		 int image_num)
{
  int i, j;
  int hr;
  int offs, imnum;
  float imi, imr;
  Filename fn;
  FILE *f;

  hr = c.res / 2;
  offs = (c.over_samp-1)*hr;
  for (i = 0; i < c.res; i++)
    for (j = 0; j < c.res; j++)
      if (hypot((float)(i-hr),(float)(j-hr)) < .51*c.res)
	{
	  imr = im[0][-j+c.res+offs-1][-i+c.res+offs-1];
	  imi = im[1][-j+c.res+offs-1][-i+c.res+offs-1];
	  imphase[i][j] = 1000 * atan2(imi,imr);
	}
      else
        imphase[i][j] = 0;

  Acct(WRITEOPEN);
  sprintf(fn, "%s/phs.s%.2d.c%d.i%.3d", c.output_directory, slice_num+1, coil_num+1, image_num + c.im_offset);
  if ((f = fopen(fn, "w")) == NULL)
    Abort("Can't open %s!\n", fn);
  Acct(WRITING);
  bio_big_endian_output = c.big_endian_output;
  FWrInt16Array(f, &imphase[0][0], c.res*c.res);
  Acct(WRITECLOSE);
  fclose(f);
  Acct(PROCESSING);
}


/* UTILITY FUNCTIONS */

void fft (float *rdat,
	  float *idat,
	  int nn,
	  int isign)
{
  int n,mmax,m,j,j1,istep,i,mmaxby2;
  double wtemp,wr,wpr,wpi,wi,theta;
  float tempr,tempi;

  n=nn << 1;
  j=1;
  for (i=0;i<nn;i++)
    {
      j1 = (j-1)/2;
      if (j1 > i)
	{
	  SWAP(rdat[j1],rdat[i]);
	  SWAP(idat[j1],idat[i]);
	}
      m=nn;
      while (m >= 2 && j > m)
	{
	  j -= m;
	  m >>= 1;
	}
      j += m;
    }
  mmax=2;
  while (n > mmax)
    {
      istep=2*mmax;
      mmaxby2 = mmax/2;
      theta=(2.0*PI)/(isign*mmax);
      wtemp=sin(0.5*theta);
      wpr = -2.0*wtemp*wtemp;
      wpi=sin(theta);
      wr=1.0;
      wi=0.0;
      for (m=0;m<mmaxby2;m++)
	{
	  for (i=m;i<nn;i+=mmax)
	    {
	      j=i+mmaxby2;
	      tempr=wr*rdat[j]-wi*idat[j];
	      tempi=wr*idat[j]+wi*rdat[j];
	      rdat[j]=rdat[i]-tempr;
	      idat[j]=idat[i]-tempi;
	      rdat[i] += tempr;
	      idat[i] += tempi;
	    }
	  wr=(wtemp=wr)*wpr-wi*wpi+wr;
	  wi=wi*wpr+wtemp*wpi+wi;
	}
      mmax=istep;
    }
}

void
PackContext ()
{
  par_pkstr(c.input_directory);
  par_pkint(c.big_endian_input);
  par_pkstr(c.tmp_directory);
  par_pkstr(c.output_directory);
  par_pkstr(c.output_name);
  par_pkstr(c.output_data);
  par_pkint(c.big_endian_output);
  par_pkint(c.ncoils);
  par_pkint(c.nslices);
  par_pkint(c.nimages);
  par_pkint(c.ndat);
  par_pkint(c.res);
  par_pkint(c.os_res);
  par_pkint(c.slice);
  par_pkfloat(c.ph_twist);
  par_pkint(c.im_offset);
  par_pkint(c.start_slice);
  par_pkint(c.end_slice);
  par_pkint(c.hc_cor);
  par_pkint(c.lin_map);
  par_pkint(c.gen_map);
  par_pkfloat(c.map_time);
  par_pkfloat(c.samp_time);
  par_pkint(c.filter_sz);
  par_pkint(c.over_samp);
  par_pkfloat(c.grid_len);
  par_pkfloat(c.gridb);
  par_pkint(c.all_coils);
  par_pkint(c.wr_phase);
  par_pkint(c.slice_num);
  if (c.hc_cor) {
    par_pkfloatarray(&c.refim_in[0][0], c.res*c.res);
    par_pkfloatarray(&c.sampim[0][0], c.os_res*c.os_res);
  }
}

void
UnpackContext ()
{
  static int previously_allocated = FALSE;
  
  par_upkstr(c.input_directory);
  c.big_endian_input= par_upkint();
  par_upkstr(c.tmp_directory);
  par_upkstr(c.output_directory);
  par_upkstr(c.output_name);
  par_upkstr(c.output_data);
  c.big_endian_output= par_upkint();
  c.ncoils= par_upkint();
  c.nslices= par_upkint();
  c.nimages= par_upkint();
  c.ndat= par_upkint();
  c.res= par_upkint();
  c.os_res= par_upkint();
  c.slice= par_upkint();
  c.ph_twist= par_upkfloat();
  c.im_offset= par_upkint();
  c.start_slice= par_upkint();
  c.end_slice= par_upkint();
  c.hc_cor= par_upkint();
  c.lin_map= par_upkint();
  c.gen_map= par_upkint();
  c.map_time= par_upkfloat();
  c.samp_time= par_upkfloat();
  c.filter_sz= par_upkint();
  c.over_samp= par_upkint();
  c.grid_len= par_upkfloat();
  c.gridb= par_upkfloat();
  c.all_coils= par_upkint();
  c.wr_phase= par_upkint();
  c.slice_num= par_upkint();
  
  /* deallocate the old arrays */
  if (previously_allocated)
    {
      Free2DFloatArray(c.refim_in);
      Free2DFloatArray(c.sampim);
      previously_allocated = FALSE;
    }
  
  if (c.hc_cor)
    {
      /* allocate the arrays */
      c.refim_in = Alloc2DFloatArray(c.res, c.res);
      c.sampim = Alloc2DFloatArray(c.os_res, c.os_res);
      previously_allocated = TRUE;
      par_upkfloat(&c.refim_in[0][0], c.res*c.res);
      par_upkfloat(&c.sampim[0][0], c.os_res*c.os_res);
    }
}

void
PackTask ()
{
  par_pkint(t.image_num);
}

void
UnpackTask ()
{
  t.image_num= par_upkint();
}

