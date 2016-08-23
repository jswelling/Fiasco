/*
 *	spiral.c	Version 1.0
 *
 *    Image Gridding & Reconstruction for Spiral k-space Acquisitions
 *                           (c) 1992,1993,1994,1995,1996,1997
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
 *     sgrid - split gridding code off from gsp14 and converted to run
 *	     - in parallel using PVM  (Greg Hood, PSC)
 *	v1.0 - bug fixes; changed to write only one sampim file per slice
 *	     -	(Greg Hood, PSC)		
 *	v1.1 - now uses the updated parallel task manager in util/par.c
 *	     -  (Greg Hood, PSC, June 96)
 *	v1.3 - uses Pittsburgh File Format for raw output
 *	     -  (Greg Hood, PSC, Feb 97)
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
 *    spiral - recombined sgrid & srecon into single executable; 
 *	     -  replaced all environment variables with command line switches;
 *	     -  some code reorganization for clarity
 *	     -  all output is now in terms of Pittsburgh File Format 1.0
 *	     -  more documentation
 *	v1.0 -  (Greg Hood, PSC, Oct 97-Jan 98)
 *      v1.1 -  added levels of verbosity control (Greg Hood, PSC, Apr 98)
 *
 */

/*
 * REFERENCE MODE [only with lin_map or gen_map switches]:
 *    INPUTS:
 *	P file containing exactly 2 timesteps
 *    OUTPUTS:
 *	reference dataset (default name: ref.mri) containing the following chunks:
 *		linear [only with lin_map switch] VZ (by default stored in ref.lin)
 *		general [only with gen_map switch] XYZ (by default stored in ref.gen)
 *		raw - [only with write_raw switch] XYVZTC Fourier space
 *			gridded images (by default stored in recon.raw)
 *		sampim - [only with write_raw switch and hc_cor switch] XYZ arrays of
 *			the sample information (by default stored in recon.sam)
 *		magnitude - [only with write_mag switch] XYZTC magnitude of the
 *			Fourier space gridded images (by default stored in recon.mag)
 *              samples - [only with write_samples switch] VPSCZT Fourier space
 *                      table of (non-gridded) samples (by default stored in recon.samp)
 *
 * STANDARD MODE:
 *    INPUTS:
 *	one or more P files
 *	reference dataset - [only with hc_cor or lin_cor switch] (by default
 *				read from ref.mri and subordinate files)
 *		linear [only with lin_cor switch] VZ
 *		general [only with hc_cor switch] XYZ
 *    OUTPUTS:
 *	output dataset (default name: recon.mri) containing the following chunks:
 *		images - reconstructed XYZT images (by default stored in recon.dat)
 *		coil<n> - [only with all_coils switch] XYZT coil-specific reconstructed
 *			images (by default stored in recon.coil<n>)
 *		raw - [only with write_raw switch] XYVZTC Fourier space
 *			gridded images (by default stored in recon.raw)
 *		sampim - [only with write_raw switch and hc_cor switch] XYZ arrays of
 *			the sample information (by default stored in recon.sam)
 *		magnitude - [only with write_mag switch] XYZTC magnitude of the
 *			Fourier space gridded images (by default stored in recon.mag)
 *              samples - [only with write_samples switch] VPSCZT Fourier space
 *                      table of (non-gridded) samples (by default stored in recon.samp)
 *		phase - [only with write_phase switch] XYZTC phase of the
 *			reconstructed images (by default stored in recon.phs)
 *		
 */

/* Notes-
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#ifdef AFS
#include <sys/types.h>
#include <sys/wait.h>
#include <errno.h>
#endif
#include "par.h"
#include "bio.h"
#include "array.h"
#include "acct.h"
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"
#include "rttraj.h"
#include "spiral.h"

static char rcsid[] = "$Id: spiral.c,v 1.27 2007/03/22 00:09:18 welling Exp $";

#define BLOCK_SIZE	(16*65536)	/* maximum # of bytes to copy
					   at one time */

/* Defaults for several task parameters */
#define DEFAULT_OVER 2
#define DEFAULT_GRIDL 1.5

/* GLOBAL VARIABLES FOR MASTER & WORKERS */
Task t;
Context c;
Result r;

/* GLOBAL VARIABLES FOR MASTER */
float **regxs = NULL;	/* [nslices][total_images]    registration x shift */
float **regys = NULL;	/* [nslices][total_images]    registration y shift */
float **regrot = NULL;	/* [nslices][total_images]    registration rotation */
int print_mode = 0;
FilenameList *thisname;
FilenameList *fname_list_root= NULL;
FilenameList *fname_list_tail= NULL;
char spiral_hosts[512];
char reference_name[128];
char id_tag[512];
unsigned char** missing= NULL; /* will hold missing info */

/* FORWARD DECLARATIONS */
static void InitMaster ();
static void ReadArguments (const int argc, const char **argv);
static void PrintMode ();
static void ProcessFirstFile ();
static void ProcessAdditionalFile ();
static void ProcessReferenceFile ();
static void ProcessFirstImageFile ();
static void ProcessImageFile ();
static void SetRegistration ();
static void CreateOutputDataset (int argc, char** argv);
static void WriteMissingInfo ();
static void ComputeCalibrationMap ();
static void LoadRegistrationData ();
static int  CheckReferences(MRI_Dataset*);
static void LoadReferences ();
static void UncompressFile (Filename out, const Filename in);
static void RemoveFile (const Filename name);
#ifdef AFS
static int MySystem( const char* command );
static int CheckForAFS( const char* fname );
static void FlushAFS( const char* fname );
#endif


int
main (int argc,
      char **argv,
      char **envp)
{
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
  return 0;
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
  if (c.verbosity >= VERBOSITY_FILE || print_mode)
    PrintMode ();

  thisname = fname_list_root;
  for (t.file_index = 0; t.file_index < c.nfiles; ++t.file_index)
    {
      Filename tname;
      int uncomp_flag= 0;
      if (((len = strlen(thisname->fname)) > 2)
	  && (!strcasecmp(&((thisname->fname)[len-2]), ".z")))
	{
	  UncompressFile(tname, thisname->fname);
	  uncomp_flag= 1;
	}
      else uncomp_flag= 0;
      strcpy(t.filename, uncomp_flag ? tname : thisname->fname);
      if (t.file_index == 0)
	ProcessFirstFile(argc,argv);
      else
	ProcessAdditionalFile();
      if (uncomp_flag)
	{
          while (par_tasks_outstanding() > 0)
            par_wait(0.1);
	  RemoveFile(tname);
	}
      thisname= thisname->next;
    }

  par_finish();

  /* Re-write the missing chunk if applicable */
  if (missing) WriteMissingInfo();

  /* Clean up files and memory */
  thisname= fname_list_root;
  while (thisname) {
    FilenameList *target = thisname;
    thisname= thisname->next;
    free(target);
  }
  fname_list_root= fname_list_tail= NULL;
  if (c.verbosity >= VERBOSITY_MINIMAL)
    Report("Done!!\n");
}

static void
InitMaster ()
{
  verbose = TRUE;

  /* initially set the variably-sized arrays to be unallocated */
  c.t2k = NULL;
  c.kdens = NULL;
  c.refim_in = NULL;
  c.sampim = NULL;
  c.refl = NULL;
  c.ref_missing= NULL;

  r.slice_sampim = NULL;

  /* Set some additional context state. */
  c.do_recon= TRUE;
}

static void
ReadArguments (const int argc, const char **argv)
{
  int file_count= 0;
  Filename input_file;
  Filename input_fullpath;

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
  /* the endianness of the input defaults to TRUE
     since the P files are generated in big-endian format */
  c.big_endian_input = TRUE;
  /* the endianness of the output defaults to the
     native format of the machine on which we are running */
  c.big_endian_output = bio_big_endian_machine;
  if (cl_present("bei")) c.big_endian_input= TRUE;
  if (cl_present("beo")) c.big_endian_output= TRUE;
  if (cl_present("lei")) c.big_endian_input= FALSE;
  if (cl_present("leo")) c.big_endian_output= FALSE;
  c.reg_2x = cl_present("2") ||
    cl_present("reg_2x");
  c.all_coils = cl_present("a") ||
    cl_present("all_coils");
  c.samp_cor = cl_present("samp_cor") ||
    cl_present("sample_correction");
  c.loc_shift = cl_present("l") || cl_present("loc_shift") ||
    cl_present("location_shift");
  c.write_raw = cl_present("write_raw");
  c.write_mag = cl_present("m") || cl_present("write_mag") ||
    cl_present("write_magnitude");
  c.write_samples = cl_present("write_samples");
  c.write_phase = cl_present("p") ||
    cl_present("write_phase");
  print_mode = cl_present("print_mode");
  c.lin_cor = cl_present("lin_cor") ||
    cl_present("linear_correction");
  c.hc_cor = cl_present("gen_cor") ||
    cl_present("general_correction") ||
    cl_present("hc_cor") ||
    cl_present("homogeneity_correction");
  c.lin_map = cl_present("lin_map") ||
    cl_present("linear_map");
  c.gen_map = cl_present("gen_map") ||
    cl_present("general_map");
  c.output_float = cl_present("float") ||
    cl_present("output_float");

  /* next handle the parameterized switches */
  cl_get("res|resolution", "%option %d[64]", &c.res);
  c.grid_len= DEFAULT_GRIDL;
  cl_get("grid_len|grid_length", "%option %f", &c.grid_len);
  c.over_samp= DEFAULT_OVER;
  cl_get("over_samp|oversampling", "%option %d", &c.over_samp);
  cl_get("slice", "%option %d[0]", &c.slice);
  cl_get("samp_delay|sample_delay", "%option %d[50]", &c.samp_delay);
  cl_get("ph_twist|phase_twist", "%option %f[0.0]", &c.ph_twist);
  cl_get("lr_shift|left_right_shift", "%option %d[0]", &c.lr_shift);
  cl_get("tb_shift|top_bottom_shift", "%option %d[0]", &c.tb_shift);
  cl_get("zoom", "%option %f[1.0]", &c.zoom);
  cl_get("mag_factor|magnitude_factor", "%option %f[0.0]", &c.mag_factor);
  cl_get("ph_factor|phase_factor", "%option %f[0.0]", &c.ph_factor);
  cl_get("filter_sz|filter_size", "%option %d[5]", &c.filter_sz);
  cl_get("r|reg_file|registration_file","%option %s[%]", "", &c.reg_file);
  cl_get("hosts", "%option %s[%]", "", spiral_hosts);
  cl_get("input_dir|input_directory", "%option %s[%]", ".", c.input_directory);
  cl_get("ref_dir|reference_directory", "%option %s[%]", ".", c.reference_directory);
  cl_get("tmp_dir|tmp_directory", "%option %s[%]", ".", c.tmp_directory);
  cl_get("output_dir|output_directory", "%option %s[%]", ".", c.output_directory);
  cl_get("out_name|output_name", "%option %s[%]",
	 (c.lin_map || c.gen_map) ? "ref.mri" : "recon.mri",
	 c.output_name);
  cl_get("out_data|output_data", "%option %s[%]", ".dat", c.output_data);
  cl_get("scale|output_scale", "%option %f[1024.0]", &c.output_scale);
  cl_get("ref|reference", "%option %s[%]", "ref.mri", reference_name);
  cl_get("v|verbosity", "%option %d[1]", &c.verbosity);
  cl_get( "tag", "%option %s[%]", "", id_tag );
  c.mapdel_set= cl_get("mapdel|map_delay","%option %lf", &c.mapdel);

  /* The following are now read from the Pfile header */
  /*
  cl_get("samp_time|sample_time", "%option %f[8.0]", &c.samp_time);
  */

  /* Strip possible quotes off spiral_hosts; shell syntax makes
   * them impossible to avoid in some cases.
   */
  if (spiral_hosts[0]=='\'' || spiral_hosts[0]=='"') {
    spiral_hosts[0]= ' ';
    spiral_hosts[strlen(spiral_hosts)-1]= '\0';
  }

  /* Bad karma from previous incarnations forces us to build a
   * linked list of the file names to be dealt with, and then shortly
   * throw it away.
   */
  file_count= 0;
  while (cl_get( "", "%s", input_file )) {
    file_count++;
    ConstructFilename(input_fullpath, c.input_directory, input_file);
    if (!(thisname = (FilenameList *) malloc(sizeof(FilenameList)))) {
      Abort("Unable to allocate %d bytes for file name parsing!",
	    sizeof(FilenameList));
    }
    strcpy(thisname->fname, input_fullpath);
    thisname->next= NULL;
    if (!fname_list_root) fname_list_root= thisname;
    if (fname_list_tail) fname_list_tail->next= thisname;
    fname_list_tail= thisname;
  }

  if (cl_cleanup_check()) {
    int i;
    fprintf(stderr,"%s: invalid argument in command line:\n    ", argv[0]);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(-1);
  }
  /*** End command-line parsing ***/

  if (file_count==0) {
    fprintf(stderr,"%s: no file argument(s) supplied.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  c.nfiles = file_count;

  if (spiral_hosts[0] != '\0')
    par_set_hosts(spiral_hosts);
}

static void
PrintMode ()
{
  char *s;
  int v;

  Report("\nProgram spiral starting\n");
  Report("\n");

  Report("Input:\n");
  Report("\tInput directory: %s\n", c.input_directory);
  Report("\tTmp directory for uncompression: %s\n", c.tmp_directory);
  Report("\t%s-endian input\n", c.big_endian_input ? "Big" : "Little");
  if (c.slice == 0)
    Report("\tReconstructing: all slices\n");
  else
    Report("\tReconstructing: only slice %d\n", c.slice);
  if (c.lin_map || c.gen_map)
    {
      Report("\tProcessing input as a reference file.\n");
      if (c.lin_map)
	Report("\tProducing a linear map.\n");
      if (c.gen_map)
	Report("\tProducing a general map.\n");
    }
  Report("\n");
	
  Report("Parameters:\n");
  Report("\tResolution: %d\n", c.res);
  Report("\tOversampling: %d\n", c.over_samp);
  Report("\tGrid length: %f\n", c.grid_len);
  Report("\tSample delay: %d\n", c.samp_delay);
  Report("\tSample time: %f\n", c.samp_time);
  if (c.mapdel_set)
    Report("\tMap delay: %f\n", c.mapdel);
  else
    Report("\tMap delay is not set\n");
  if (c.lin_map || c.gen_map)
    Report("\tFilter size: %d\n", c.filter_sz);
  Report("\n");

  Report("Corrections:\n");
  if (c.reg_file[0] != '\0')
    {
      Report("\tRegistration correction applied from file: %s\n", c.reg_file);
      if (c.reg_2x)	
	Report("\tRegistration translations will be doubled\n");
    }
  if (c.lin_cor || c.hc_cor)
    Report("\tReference directory: %s\n", c.reference_directory);
  if (c.lin_cor)
    Report("\tLinear correction applied from reference dataset %s\n", reference_name);
  if (c.hc_cor)
    Report("\tGeneral (homogeneity) correction applied from reference dataset %s\n", reference_name);
  Report("\tSample density correction: %s\n", c.samp_cor ? "yes" : "no");
  Report("\n");

  Report("Adjustments:\n");
  Report("\tLocation shift: %s\n", c.loc_shift ? "yes" : "no");
  Report("\tLeft-right shift: %d\n", c.lr_shift);
  Report("\tTop-bottom shift: %d\n", c.tb_shift);
  Report("\tMagnitude factor: %f\n", c.mag_factor);
  Report("\tPhase factor: %f\n", c.ph_factor);
  Report("\tPhase twist: %f\n", c.ph_twist);
  Report("\tZoom: %f\n", c.zoom);
  Report("\n");

  Report("Output:");
  Report("\tOutput directory: %s\n", c.output_directory);
  Report("\tOutput dataset name: %s\n", c.output_name);
  Report("\tOutput datafile name or extension: %s\n", c.output_data);
  Report("\t%s-endian output\n", c.big_endian_output ? "Big" : "Little");
  Report("\tOutput format: %s\n", c.output_float ? "float32" : "int16");
  Report("\tOutput scaling factor: %f\n", c.output_scale);
  Report("\n");

  if (c.write_raw || c.write_mag || c.write_phase || c.all_coils
      || c.write_samples)
    {
      Report("Auxiliary output options:\n");
      if (c.write_samples)
	Report("\tWriting ungridded K-space samples\n");
      if (c.write_raw)
	Report("\tWriting raw file\n");
      if (c.write_mag)
	Report("\tWriting magnitude file\n");
      if (c.write_phase)
	Report("\tWriting phase file\n");
      if (c.all_coils)
	Report("\tReconstructing all coils\n");
      Report("\n");
    }
  
  if ((s = getenv("PAR_ENABLE")) != NULL &&
      sscanf(s, "%d", &v) == 1 &&
      v != 0)
    {
      Report("Parallel Operation:\n");
      Report("\tEnabled\n");
      Report("\tHosts: %s\n", spiral_hosts);
      Report("\n");
    }
  else if (spiral_hosts[0] != '\0')
    {
      Report("Parallel Operation:\n");
      Report("\tDisabled since PAR_ENABLE not set to 1\n");
      Report("\tHosts: %s\n", spiral_hosts);
      Report("\n");
    }

  Report("\n");
}

static void
ProcessFirstFile (int argc, char** argv)
{
  if (check_filetype(t.filename)==FILE_PGH_MRI) {
    c.input_is_mri= TRUE;
    ReadMRIFileHeader(t.filename, &c);
  }
  else {
    c.input_is_mri= FALSE;
    ReadFileHeader(t.filename, &c, NULL);
  }
  ComputeCalibrationMap();
  CreateOutputDataset(argc, argv);
  LoadReferences();
  LoadRegistrationData();
  par_set_context();
  if (c.gen_map || c.lin_map)
    ProcessReferenceFile();
  else
    ProcessFirstImageFile();
}

static void
ProcessAdditionalFile ()
{
  /* make sure this P file is consistent
     with the first P files that we read */
  if (check_filetype(t.filename)==FILE_PGH_MRI) 
    CheckMRIHeaderInfo(t.filename, &c);
  else
    CheckHeaderInfo(t.filename, &c);
  
  /* dispatch most of the work to the
     workers */
  ProcessImageFile(0);
}

static void
ProcessReferenceFile ()
{
  /* if we are generating reference files, assign an entire slice
     to each worker; that worker will iterate through all image sets
     and coils within the slice; this is necessary because the images
     within a slice must be processed sequentially */
  for (t.slice_num = c.start_slice; t.slice_num <= c.end_slice; ++t.slice_num)
    par_delegate_task();
}

static void
ProcessFirstImageFile ()
{
  /* if we are doing homogeneity correction, process the first image of all slices
     in parallel so that we can obtain sampim information for all slices prior to
     processing subsequent timepoints of each slice */
  if (c.hc_cor)
    {
      t.image_num = 0;
      t.overall_image_num = 0;
      for (t.slice_num = c.start_slice; t.slice_num <= c.end_slice; ++t.slice_num)
	{
	  SetRegistration();
	  par_delegate_task();
	}
      while (par_tasks_outstanding() > 0)
	par_wait(0.1);

      /* reset the context so that the sampim and refim arrays will be
	 distributed to the workers */
      par_set_context();

      /* start with image 1 (the second timepoint), since we have already done
	 the first */
      ProcessImageFile(1);
    }
  else
    ProcessImageFile(0);
}

static void
ProcessImageFile (int start_image_num)
{
  for (t.slice_num = c.start_slice; t.slice_num <= c.end_slice; ++t.slice_num)
    {
      /* if we are generating image files, assign each multi-coil image set
	 to a worker; that worker will iterate through all coils */
      for (t.image_num = start_image_num; t.image_num < c.nimages; ++t.image_num)
	{
	  t.overall_image_num = t.file_index * c.nimages + t.image_num;
	  SetRegistration();
	  par_delegate_task();
	}
    }
}

static void
SetRegistration ()
{
  if (c.reg_file[0] != '\0')
    {
      t.reg_xs = regxs[t.slice_num][t.overall_image_num];
      t.reg_ys = regys[t.slice_num][t.overall_image_num];
      t.reg_rot = regrot[t.slice_num][t.overall_image_num];
    }
  else
    t.reg_xs = t.reg_ys = t.reg_rot = 0.0;
}

static void
CreateOutputDataset (int argc, char** argv)
{
  Filename fn;
  MRI_Dataset *ds= NULL;
  int i;
  char name[32];

  Acct(WRITEOPEN);
  ConstructFilename(fn, c.output_directory, c.output_name);
  ds = mri_open_dataset(fn, MRI_WRITE);

  Acct(WRITING);
  if (!c.lin_map && !c.gen_map)
    {
      mri_create_chunk(ds, "images");
      mri_set_string(ds, "images.datatype", 
		     c.output_float ? "float32" : "int16");
      mri_set_string(ds, "images.dimensions", "xyzt");
      mri_set_int(ds, "images.extent.x", c.res);
      mri_set_string(ds, "images.description.x", "gridded image-space");
      mri_set_int(ds, "images.extent.y", c.res);
      mri_set_string(ds, "images.description.y", "gridded image-space");
      mri_set_int(ds, "images.extent.z", c.end_slice - c.start_slice + 1);
      mri_set_string(ds, "images.description.z", "gridded image-space");
      mri_set_int(ds, "images.extent.t", c.total_images);
      mri_set_string(ds, "images.description.t", "gridded image-space");
      if (c.start_slice > 0)
	mri_set_int(ds, "images.min.z", c.start_slice);
      
      mri_set_string(ds, "images.file", c.output_data);

      /* Additional info available from Pfile header */
      mri_set_float( ds, "images.voxel_size.x",10.0*c.opfov/c.res);
      mri_set_float( ds, "images.voxel_size.y",10.0*c.opfov/c.res);
      mri_set_float( ds, "images.voxel_size.z",c.slthick);
      mri_set_float( ds, "images.voxel_spacing.x",10.0*c.opfov/c.res);
      mri_set_float( ds, "images.voxel_spacing.y",10.0*c.opfov/c.res);
      mri_set_float( ds, "images.voxel_spacing.z",c.spacing+c.slthick);
      mri_set_float( ds, "images.fov.x", 10.0*c.opfov);
      mri_set_float( ds, "images.fov.y", 10.0*c.opfov);
      mri_set_float( ds, "images.fov.z", 
		     (c.nslices*c.slthick) + ((c.nslices-1)*c.spacing));
      mri_set_float( ds, "images.tr", c.tr/1000.0);
      mri_set_float( ds, "images.te", c.te/1000.0);
      if (*(c.date))
	mri_set_string( ds, "images.date", c.date);
      if (*(c.time))
	mri_set_string( ds, "images.time", c.time);
      if (id_tag[0])
	mri_set_string( ds, "images.scan.id", id_tag );

      /* create the coil-specific chunks */
      if (c.all_coils)
	for (i = 0; i < c.ncoils; ++i)
	  {
	    sprintf(name, "coil%d", i);
	    Acct(WRITING);
	    mri_create_chunk(ds, name);
	    mri_set_string(ds, mri_cat(name, ".datatype"),
			   c.output_float ? "float32" : "int16");
	    mri_set_string(ds, mri_cat(name, ".dimensions"), "xyzt");
	    mri_set_int(ds, mri_cat(name, ".extent.x"), c.res);
	    mri_set_string(ds, mri_cat(name,".description.x"), 
				       "gridded image-space");
	    mri_set_int(ds, mri_cat(name, ".extent.y"), c.res);
	    mri_set_string(ds, mri_cat(name,".description.y"), 
				       "gridded image-space");
	    mri_set_int(ds, mri_cat(name, ".extent.z"), c.end_slice - c.start_slice + 1);
	    mri_set_string(ds, mri_cat(name,".description.z"), 
				       "gridded image-space");
	    mri_set_int(ds, mri_cat(name, ".extent.t"), c.total_images);
	    mri_set_string(ds, mri_cat(name,".description.t"), 
				       "gridded image-space");
	    if (c.start_slice > 0)
	      mri_set_int(ds, mri_cat(name, ".min.z"), c.start_slice);
	    
	    if (c.output_data[0] == '.' || c.output_data[0] == '\0')
	      mri_set_string(ds, mri_cat(name, ".file"), c.output_data);
	    else
	      {
		sprintf(fn, ".%s", name);
		mri_set_string(ds, mri_cat(name, ".file"), fn);
	      }
	  }

      if (c.write_phase)
	{
	  mri_create_chunk(ds, "phase");
	  mri_set_string(ds, "phase.datatype", "int16");
	  mri_set_string(ds, "phase.dimensions", "xyztc");
	  mri_set_int(ds, "phase.extent.x", c.res);
	  mri_set_string(ds, "phase.description.x", "gridded image-space");
	  mri_set_int(ds, "phase.extent.y", c.res);
	  mri_set_string(ds, "phase.description.y", "gridded image-space");
	  mri_set_int(ds, "phase.extent.z", c.end_slice - c.start_slice + 1);
	  mri_set_string(ds, "phase.description.z", "gridded image-space");
	  mri_set_int(ds, "phase.extent.t", c.total_images);
	  mri_set_string(ds, "phase.description.t", "gridded image-space");
	  mri_set_int(ds, "phase.extent.c", c.ncoils);
	  mri_set_string(ds, "phase.description.c", "discrete");
	  if (c.start_slice > 0)
	    mri_set_int(ds, "phase.min.z", c.start_slice);
	  mri_set_string(ds, "phase.file", ".phs");
	}
    }

  /* create the reference chunks */
  if (c.lin_map)
    {
      mri_create_chunk(ds, "linear");
      mri_set_string(ds, "linear.datatype", "float32");
      mri_set_string(ds, "linear.dimensions", "vz");
      mri_set_int(ds, "linear.extent.v", 3);
      mri_set_int(ds, "linear.extent.z", c.end_slice - c.start_slice + 1);
      if (c.start_slice > 0)
	mri_set_int(ds, "linear.min.z", c.start_slice);
      mri_set_string(ds, "linear.file", ".lin");
      mri_set_float( ds, "linear.voxel_size.x",10.0*c.opfov/c.res);
      mri_set_float( ds, "linear.voxel_size.y",10.0*c.opfov/c.res);
      mri_set_float( ds, "linear.voxel_size.z",c.slthick);
      mri_set_float( ds, "linear.voxel_spacing.x",10.0*c.opfov/c.res);
      mri_set_float( ds, "linear.voxel_spacing.y",10.0*c.opfov/c.res);
      mri_set_float( ds, "linear.voxel_spacing.z",c.spacing+c.slthick);
      mri_set_float( ds, "linear.fov.x", 10.0*c.opfov);
      mri_set_float( ds, "linear.fov.y", 10.0*c.opfov);
      mri_set_float( ds, "linear.fov.z", 
		     (c.nslices*c.slthick) + ((c.nslices-1)*c.spacing));
      mri_set_float( ds, "linear.tr", c.tr/1000.0);
      mri_set_float( ds, "linear.te", c.te/1000.0);
      if (*(c.date))
	mri_set_string( ds, "linear.scan.date", c.date);
      if (*(c.time))
	mri_set_string( ds, "linear.scan.time", c.time);
      if (id_tag[0])
	mri_set_string( ds, "linear.scan.id", id_tag );
    }
  if (c.gen_map)
    {
      mri_create_chunk(ds, "general");
      mri_set_string(ds, "general.datatype", "float32");
      mri_set_string(ds, "general.dimensions", "xyz");
      mri_set_int(ds, "general.extent.x", c.res);
      mri_set_int(ds, "general.extent.y", c.res);
      mri_set_int(ds, "general.extent.z", c.end_slice - c.start_slice + 1);
      if (c.start_slice > 0)
	mri_set_int(ds, "general.min.z", c.start_slice);
      mri_set_string(ds, "general.file", ".gen");
      mri_set_float( ds, "general.voxel_size.x",10.0*c.opfov/c.res);
      mri_set_float( ds, "general.voxel_size.y",10.0*c.opfov/c.res);
      mri_set_float( ds, "general.voxel_size.z",c.slthick);
      mri_set_float( ds, "general.voxel_spacing.x",10.0*c.opfov/c.res);
      mri_set_float( ds, "general.voxel_spacing.y",10.0*c.opfov/c.res);
      mri_set_float( ds, "general.voxel_spacing.z",c.spacing+c.slthick);
      mri_set_float( ds, "general.fov.x", 10.0*c.opfov);
      mri_set_float( ds, "general.fov.y", 10.0*c.opfov);
      mri_set_float( ds, "general.fov.z", 
		     (c.nslices*c.slthick) + ((c.nslices-1)*c.spacing));
      mri_set_float( ds, "general.tr", c.tr/1000.0);
      mri_set_float( ds, "general.te", c.te/1000.0);
      if (*(c.date))
	mri_set_string( ds, "general.scan.date", c.date);
      if (*(c.time))
	mri_set_string( ds, "general.scan.time", c.time);
      if (id_tag[0])
	mri_set_string( ds, "general.scan.id", id_tag );
    }

  /* create the raw chunks */
  if (c.write_raw)
    {
      Acct(WRITING);
      mri_create_chunk(ds, "raw");
      mri_set_string(ds, "raw.datatype", "float32");
      mri_set_string(ds, "raw.dimensions", "xyvztc");
      mri_set_int(ds, "raw.extent.x", c.os_res);
      mri_set_string(ds, "raw.description.x", "gridded k-space");
      mri_set_int(ds, "raw.extent.y", c.os_res);
      mri_set_string(ds, "raw.description.y", "gridded k-space");
      mri_set_int(ds, "raw.extent.v", 2);
      mri_set_string(ds, "raw.description.v", "complex real/imaginary");
      mri_set_int(ds, "raw.extent.z", c.nslices);
      mri_set_string(ds, "raw.description.z", "gridded image-space");
      mri_set_int(ds, "raw.extent.t", c.total_images);
      mri_set_string(ds, "raw.description.t", "gridded image-space");
      mri_set_int(ds, "raw.extent.c", c.ncoils);
      mri_set_string(ds, "raw.description.c", "discrete");
      mri_set_string(ds, "raw.file", ".raw");

      if (c.hc_cor)
	{
	  mri_create_chunk(ds, "sampim");
	  mri_set_string(ds, "sampim.datatype", "float32");
	  mri_set_string(ds, "sampim.dimensions", "xyz");
	  mri_set_int(ds, "sampim.extent.x", c.os_res);
	  mri_set_string(ds, "sampim.description.x", "gridded k-space");
	  mri_set_int(ds, "sampim.extent.y", c.os_res);
	  mri_set_string(ds, "sampim.description.y", "gridded k-space");
	  mri_set_int(ds, "sampim.extent.z", c.nslices);
	  mri_set_string(ds, "sampim.description.z", "gridded image-space");
	  mri_set_string(ds, "sampim.file", ".sam");
	}
      Acct(PROCESSING);
    }

  /* create the magnitude chunk */
  if (c.write_mag)
    {
      Acct(WRITING);
      mri_create_chunk(ds, "magnitude");
      mri_set_string(ds, "magnitude.datatype", "int16");
      mri_set_string(ds, "magnitude.dimensions", "xyztc");
      mri_set_int(ds, "magnitude.extent.x", c.os_res);
      mri_set_string(ds, "magnitude.description.x", "gridded k-space");
      mri_set_int(ds, "magnitude.extent.y", c.os_res);
      mri_set_string(ds, "magnitude.description.y", "gridded k-space");
      mri_set_int(ds, "magnitude.extent.z", c.nslices);
      mri_set_string(ds, "magnitude.description.z", "gridded image-space");
      mri_set_int(ds, "magnitude.extent.t", c.total_images);
      mri_set_string(ds, "magnitude.description.t", "gridded image-space");
      mri_set_int(ds, "magnitude.extent.c", c.ncoils);
      mri_set_string(ds, "magnitude.description.c", "discrete");
      mri_set_string(ds, "magnitude.file", ".mag");
      Acct(PROCESSING);
    }

  /* create the samples chunk */
  if (c.write_samples)
    {
      int ishot;
      int icoil;
      int islice;

      Acct(WRITING);
      mri_create_chunk(ds, "samples");
      mri_set_string(ds, "samples.datatype", "float32");
      mri_set_string(ds, "samples.dimensions", "vpsczt");
      mri_set_int(ds, "samples.extent.v", 2);
      mri_set_string(ds, "samples.description.v", "complex real/imaginary");
      mri_set_int(ds, "samples.extent.z", c.nslices);
      mri_set_string(ds, "samples.description.z", "gridded image-space");
      mri_set_int(ds, "samples.extent.t", c.total_images);
      mri_set_string(ds, "samples.description.t", "gridded image-space");
      mri_set_int(ds, "samples.extent.c", c.ncoils);
      mri_set_string(ds, "samples.description.c", "discrete");
      mri_set_int(ds, "samples.extent.p", c.ndat);
      mri_set_string(ds, "samples.description.p", "ungridded k-space");
      mri_set_int(ds, "samples.extent.s", c.npr);
      mri_set_string(ds, "samples.description.s", "discrete");

      mri_set_string(ds, "samples.psd", c.psd);
      mri_set_int(ds, "samples.nph1", c.nph1);
      mri_set_int(ds, "samples.nphmult", c.nphmult);
      mri_set_int(ds, "samples.nimages", c.nimages);
      mri_set_int(ds, "samples.ndatfr", c.ndatfr);
      mri_set_int(ds, "samples.chop", c.chop);
      mri_set_int(ds, "samples.nex", c.nex);
      mri_set_int(ds, "samples.densamp", c.densamp);
      mri_set_int(ds, "samples.sliceorder", c.sliceorder);
      mri_set_int(ds, "samples.gtype", c.gtype);
      mri_set_int(ds, "samples.opxres", c.opxres);
      mri_set_float(ds, "samples.slewrate", c.slewrate);
      mri_set_int(ds, "samples.ngap", c.ngap);
      mri_set_int(ds, "samples.concat", c.concat);
      mri_set_float(ds, "samples.fast_rec_lpf", c.fast_rec_lpf);
      if (c.mapdel_set)
	mri_set_float(ds, "samples.mapdel", c.mapdel);
      mri_set_float(ds, "samples.samp_time", c.samp_time);
      mri_set_float(ds, "samples.fsgcm", c.fsgcm);
      mri_set_int(ds, "samples.risetime", c.risetime);
      mri_set_float(ds, "samples.opfov", c.opfov);
      mri_set_float(ds, "samples.ts", c.ts);
      mri_set_float(ds, "samples.gts", c.gts);
      mri_set_string(ds, "samples.time", c.time);
      mri_set_string(ds, "samples.date", c.date);
      mri_set_float(ds, "samples.slthick", c.slthick);
      mri_set_int(ds, "samples.tr", c.tr);
      mri_set_int(ds, "samples.te", c.te);
      mri_set_int(ds, "samples.flip", c.flip);
      mri_set_float(ds, "samples.spacing", c.spacing);
      mri_set_float(ds, "samples.tlc.0", c.tlc[0]);
      mri_set_float(ds, "samples.tlc.1", c.tlc[1]);
      mri_set_float(ds, "samples.tlc.2", c.tlc[2]);
      mri_set_float(ds, "samples.trc.0", c.trc[0]);
      mri_set_float(ds, "samples.trc.1", c.trc[1]);
      mri_set_float(ds, "samples.trc.2", c.trc[2]);
      mri_set_float(ds, "samples.brc.0", c.brc[0]);
      mri_set_float(ds, "samples.brc.1", c.brc[1]);
      mri_set_float(ds, "samples.brc.2", c.brc[2]);
      mri_set_float(ds, "samples.ctr.0", c.ctr[0]);
      mri_set_float(ds, "samples.ctr.1", c.ctr[1]);
      mri_set_float(ds, "samples.ctr.2", c.ctr[2]);
      mri_set_int(ds,"samples.index_rots", c.rotation);
      mri_set_int(ds,"samples.index_tpose", c.transpose);
      if (c.loc_shift) {
	mri_set_float(ds, "samples.pix_shifth", c.pix_shifth);
	mri_set_float(ds, "samples.pix_shiftv", c.pix_shiftv);
      }
      mri_set_string(ds, "samples.file", ".samp");

      mri_create_chunk(ds, "sample_kxloc");
      mri_set_string(ds, "sample_kxloc.datatype", "float32");
      mri_set_string(ds, "sample_kxloc.dimensions", "pscz");
      mri_set_int(ds, "sample_kxloc.extent.z", c.nslices);
      mri_set_string(ds, "sample_kxloc.description.z", "gridded image-space");
      mri_set_int(ds, "sample_kxloc.extent.c", c.ncoils);
      mri_set_string(ds, "sample_kxloc.description.c", "discrete");
      mri_set_int(ds, "sample_kxloc.extent.p", c.ndat);
      mri_set_string(ds, "sample_kxloc.description.p", "ungridded k-space");
      mri_set_int(ds, "sample_kxloc.extent.s", c.npr);
      mri_set_string(ds, "sample_kxloc.description.s", "discrete");
      mri_set_string(ds, "sample_kxloc.file", ".sampxloc");
      mri_create_chunk(ds, "sample_kyloc");
      mri_set_string(ds, "sample_kyloc.datatype", "float32");
      mri_set_string(ds, "sample_kyloc.dimensions", "pscz");
      mri_set_int(ds, "sample_kyloc.extent.z", c.nslices);
      mri_set_string(ds, "sample_kyloc.description.z", "gridded image-space");
      mri_set_int(ds, "sample_kyloc.extent.c", c.ncoils);
      mri_set_string(ds, "sample_kyloc.description.c", "discrete");
      mri_set_int(ds, "sample_kyloc.extent.p", c.ndat);
      mri_set_string(ds, "sample_kyloc.description.p", "ungridded k-space");
      mri_set_int(ds, "sample_kyloc.extent.s", c.npr);
      mri_set_string(ds, "sample_kyloc.description.s", "discrete");
      mri_set_string(ds, "sample_kyloc.file", "c.sampyloc");

      for (islice=0; islice<c.nslices; islice++)
	for (icoil=0; icoil<c.ncoils; icoil++)
	  for (ishot=0; ishot<c.npr; ishot++) {
	    int offset = ((((islice*c.ncoils)+icoil)*c.npr)+ishot)*c.ndat;
	    mri_set_chunk(ds, "sample_kxloc", c.ndat, offset,
			  MRI_FLOAT, c.t2k[0][ishot]);
	    mri_set_chunk(ds, "sample_kyloc", c.ndat, offset,
			  MRI_FLOAT, c.t2k[1][ishot]);
	  }

      Acct(PROCESSING);
    }

  /* Create the missing chunk, and keep a copy of its contents. */
  if (!c.lin_map && !c.gen_map)
    missing= get_missing(ds);
  else {
    mri_create_chunk(ds, "missing");
    mri_set_string(ds, "missing.datatype", "uint8");
    mri_set_string(ds, "missing.dimensions", "zt");
    mri_set_int(ds, "missing.extent.z", c.end_slice - c.start_slice + 1);
    mri_set_string(ds, "missing.description.z", "gridded image-space");
    mri_set_int(ds, "missing.extent.t", 1);
    mri_set_string(ds, "missing.description.t", "gridded image-space");
    missing= NULL; /* the master saves no missing info for reference files */
  }

  /* If we are using a reference map, it was presumably made from
   * the first two images of this dataset.  Thus we are obliged
   * to mark those images missing in the output.
   */
  if ((missing != NULL) && (c.lin_cor || c.hc_cor) && (c.total_images>2)) {
    for (i=0; i<(c.end_slice - c.start_slice + 1); i++)
      missing[0][i]= missing[1][i]= (unsigned char)1;
  }

  if (c.input_is_mri) {
    MRI_Dataset *ids= NULL;
    char* s;
    const char* val;
    int i;

    /* We don't need to flush AFS at this point, since it was
     * flushed very recently.
     */
    Acct(READOPEN);
    if (!(ids= mri_open_dataset(t.filename, MRI_READ)))
      Abort("Can't reopen %s\n", t.filename);

    /* Copy in any history info */
    for (i=1; (val=hist_get(ids,i)) != NULL; i++)
      hist_add(ds, val);

    /* Copy "samples.*" tags from input if available, except those which 
     * specifically deal with the chunk data.  Many are useful scan info.
     */
    Acct(WRITING);
    mri_iterate_over_keys(ids);

    while ((s=mri_next_key(ids)) != NULL) {
      if (!strncmp(s,"samples.",8)
	  && strcmp(s,"samples.datatype")
	  && strcmp(s,"samples.dimensions")
	  && strcmp(s,"samples.file")
	  && strcmp(s,"samples.order")
	  && strcmp(s,"samples.offset")
	  && strcmp(s,"samples.size")
	  && strcmp(s,"samples.little_endian")
	  && strncmp(s,"samples.extent.",strlen("samples.extent."))
	  && strcmp(val=mri_get_string(ids,s),"[chunk]")) {
	char buf[256];
	if (snprintf(buf,sizeof(buf),"images.%s",s+strlen("samples."))<0)
	  Abort("%s: key <%s> is too long!\n",s);
	if (!mri_has(ids,buf))
	  mri_set_string(ds,buf,val);
      }
    }

    /*
     * Also copy 'acq_blocks' info if available.
     */
    if (mri_has(ids,"acq_blocks")
	&& !strcmp(mri_get_string(ids,"acq_blocks"),"[chunk]")) {
      long long size;
      long bytesThisBlock;
      long long offset= 0;
      void* buf= NULL;

      mri_iterate_over_keys(ids);
      
      while ((s=mri_next_key(ids)) != NULL) {
	if (!strncmp(s,"acq_blocks.",11) && 
	    strncmp(mri_get_string(ids,s),"[chunk]",7))
	  mri_set_string(ds,s,mri_get_string(ids,s));
      }

      mri_create_chunk(ds, "acq_blocks");
      size = mri_get_int(ids, "acq_blocks.size");

      offset= 0;
      while (size>0) {
	bytesThisBlock= (size>BLOCK_SIZE)?BLOCK_SIZE:size;
	buf= mri_get_chunk(ids, "acq_blocks", bytesThisBlock, offset, MRI_RAW);
	mri_set_chunk(ds, "acq_blocks", bytesThisBlock, offset, MRI_RAW, buf);
	size -= bytesThisBlock;
	offset += bytesThisBlock;
      }
    }
  
    Acct(READCLOSE);
    mri_close_dataset(ids);
  }

  hist_add_cl(ds,argc,argv);
  Acct(WRITECLOSE);
  mri_close_dataset(ds);
}

static void
WriteMissingInfo ()
{
  Filename fn;
  MRI_Dataset *ds;

  if (missing==NULL) return; 

  Acct(WRITEOPEN);
  ConstructFilename(fn, c.output_directory, c.output_name);
  ds = mri_open_dataset(fn, MRI_MODIFY);

  Acct(WRITING);

  mri_set_chunk( ds, "missing", 
		 (long)(c.total_images*
			(c.end_slice - c.start_slice + 1)), 
		 0, MRI_UNSIGNED_CHAR, *missing );

  Acct(WRITECLOSE);
  mri_close_dataset(ds);
}

static void
ComputeCalibrationMap ()
{
  int i, j, n;
  double kx,ky,kxo,kyo,ggx,ggy,ggm,kksr,kksi;
  FILE *f;
  Filename fn;
  int nsamp;
  int res;

  /* kaiser-bessel functions - initialize terms */
  c.gridb = PI*c.grid_len;
  c.wind = c.grid_len;

  if (!c.input_is_mri) {
    res = getrttrajghg(c.opxres, c.npr, c.ts, c.gts, c.fsgcm, c.opfov,
		       c.slewrate, c.gts*21000, c.gtype, c.t2k[0], c.t2k[1]);
  }
  else {
    /* Trajectory info was loaded when the file header info was read */
  }

  if (c.verbosity >= VERBOSITY_FILE) {
    Report("\tndat = %d, npr = %d\n", c.ndat, c.npr);
    Report("\tts = %f, gts = %f, fsgcm = %f\n", c.ts, c.gts, c.fsgcm);
    Report("\topfov = %f, risetime = %d\n", c.opfov, c.risetime);
    Report("\tacq. matrix = %d, nom. resolution = %f mm\n",
	   c.opxres,10*c.opfov/c.opxres); 
    Report("\tk-space points = %d\n",c.res); 
    Report("\n");
  }

  nsamp = c.ndat;

  if (!c.samp_cor) 
    {
      /* no sample density correction */
      for (i=0; i < c.ndat; i++)
	c.kdens[i] = 1.0; 
      return;
    }

  /* sample density correction */
  /* Using Craig's formula from MRM, 28, 202-213 (1992) */
  c.kdens[0] = 0.0;
  for (i=1; i < nsamp; i++)
    {
      /* might wish to add correction for delay */
      kx = c.t2k[0][0][i];
      ky = c.t2k[1][0][i];
      kxo = c.t2k[0][0][i-1];
      kyo = c.t2k[1][0][i-1];
      kksr = kx*kxo + ky*kyo;
      kksi = ky*kxo - kx*kyo;
      ggx = kx - kxo;
      ggy = ky - kyo;
      if (((ggm = hypot(ggx,ggy)) > 0.0) && (hypot(kksr,kksi) > 0.0) )
	c.kdens[i] = ggm * fabs(sin(atan2(ggy,ggx) - atan2(ky,kx)) * (hypot(kx,ky) - hypot(kxo,kyo))/atan2(kksi,kksr));
      else
	c.kdens[i] = 0.0;
      /* components to density correction:
	 ggm = gradient strength; linear velocity
	 sin(atan2(ggy,ggx) - atan2(ky,kx)) =  outward component of velocity 
	     These first two came from the Meyer paper.
	 (hypot(kx,ky) - hypot(kxo,kyo))/atan2(kksi,kksr) = radial density
	     (distance outward)/(angular distance along arc)
	     This term is different from the previous term in that this
	     measure true line-line distances in the radial direction.  
	     It was added for variable density trajectories (except for the
	     the first 2 or 3 points, it is constant for uniform density
	     spiral trajectories. 
       */
    }
  for (i=nsamp; i < c.ndat; i++)
    c.kdens[i] = 0.0; 
}

static void
LoadRegistrationData ()
{
  FILE *f;
  int ss;
  int pp;
  float mse;
  Filename fn;
  char buf[256];
  int line;

  /* check if a registration file was specified */
  if (c.reg_file[0] != '\0')
    {
      /* allocate the arrays to hold the registration data */
      regxs = Alloc2DFloatArray(c.nslices, c.nimages * c.nfiles);
      regys = Alloc2DFloatArray(c.nslices, c.nimages * c.nfiles);
      regrot = Alloc2DFloatArray(c.nslices, c.nimages * c.nfiles);

      Acct(READOPEN);
      if ((f = fopen(c.reg_file, "r")) == NULL)
	Abort("Can't open %s!\n", c.reg_file);
      Acct(READING);
      line= 0;
      while (!feof(f) && !ferror(f)) {
	float v1, v2, v3, v4;	
	fgets(buf, sizeof(buf), f);	
	if (buf[0]=='#') continue; /* skip comments */
	if (sscanf(buf,"%d %d %f %f %f %f",&pp,&ss,&v1,&v2,&v3,&v4) == 6) {
	  if (ss < 0 || ss >= c.nslices ||
	      pp < 0 || pp >= c.nimages*c.nfiles)
	    Abort("Corrupt registration data: %d (max %d) %d (max %d) near line %d\n",
		  ss, c.nslices, pp, c.nimages, line);
	  regxs[ss][pp]= v1;
	  regys[ss][pp]= v2;
	  regrot[ss][pp]= v3;
	  mse= v4;
	  if (c.verbosity>=VERBOSITY_TIMESTEP) {
	    Report("Registration vals x= %f, y= %f, rot= %f\n",
		    regxs[ss][pp], regys[ss][pp], regrot[ss][pp]);
	  }
	  if (c.reg_2x)
	    {
	      regxs[ss][pp] *= 2.0;
	      regys[ss][pp] *= 2.0;
	    }
	}
	else {
	    Abort("Bad registration data near line %d\n",line);
	}
	line++;
      }
      Acct(READCLOSE);
      fclose(f); 
      Acct(PROCESSING);
    }
}

static int CheckReferences(MRI_Dataset* ds)
{
  if (c.lin_cor) {
    if (!mri_has(ds,"linear")
	|| !mri_has(ds,"linear.dimensions") 
	|| strcmp(mri_get_string(ds,"linear.dimensions"),"vz")
	|| !mri_has(ds,"linear.extent.v") 
	|| (mri_get_int(ds,"linear.extent.v") != 3)
	|| !mri_has(ds,"linear.extent.z") 
	|| (mri_get_int(ds,"linear.extent.z") != c.nslices))
      return 0;
	
  }

  if (c.hc_cor) {
    if (!mri_has(ds,"general")
	|| !mri_has(ds,"general.dimensions") 
	|| strcmp(mri_get_string(ds,"general.dimensions"),"xyz")
	|| !mri_has(ds,"general.extent.x") 
	|| (mri_get_int(ds,"general.extent.x") != c.res)
	|| !mri_has(ds,"general.extent.y") 
	|| (mri_get_int(ds,"general.extent.y") != c.res)
	|| !mri_has(ds,"general.extent.z") 
	|| (mri_get_int(ds,"general.extent.z") != c.nslices))
      return 0;
  }

  if (mri_has(ds, "missing") && 
      !strcmp(mri_get_string(ds,"missing"),"[chunk]")) {
    if (!mri_has(ds,"missing.dimensions") 
	|| strcmp(mri_get_string(ds,"missing.dimensions"),"zt")
	|| !mri_has(ds,"missing.extent.t") 
	|| (mri_get_int(ds,"missing.extent.t") != 1)
	|| !mri_has(ds,"missing.extent.z") 
	|| (mri_get_int(ds,"missing.extent.z") != c.nslices))
      return 0;
  }

  return 1;
}

static void
LoadReferences ()
{
  MRI_Dataset *refds;
  Filename fn;
  int z;

  for (z=0; z<c.nslices; z++) {
    c.refl[z][0] = c.refl[z][1] = c.refl[z][2] = 0.0;
    c.ref_missing[z]= 0;
  }

  if (c.lin_cor || c.hc_cor)
    {
      ConstructFilename(fn, c.reference_directory, reference_name);
      refds = mri_open_dataset(fn, MRI_READ);
      if (!CheckReferences(refds))
	Abort("Reference file is not appropriate for this data!\n");

      if (c.lin_cor) {
	Acct(READING);
	for (z=0; z<c.nslices; z++) {
	  memcpy(c.refl[z],
		 mri_get_chunk(refds, "linear", 3, 3*z, MRI_FLOAT),
		 3*sizeof(float));
	}
	Acct(PROCESSING);
	if (c.reverse) {
	  /* Need to flip the sign of reference on input. */
	  int i;
	  for (z=0; z<c.nslices; z++)
	    for (i=0; i<3; i++)
	      c.refl[z][i] *= -1.0;
	}
      }

      if (c.hc_cor) {
	Acct(READING);
	memcpy(&c.refim_in[0][0][0],
	       mri_get_chunk(refds, "general", c.res*c.res*c.nslices, 0, MRI_FLOAT),
	       c.res*c.res*c.nslices*sizeof(float));
	Acct(PROCESSING);
	if (c.reverse) {
	  /* Need to flip the sign of reference on input. */
	  int i;
	  int j;
	  for (z=0; z<c.nslices; z++)
	    for (i=0; i<c.res; i++)
	      for (j=0; j<c.res; j++)
		c.refim_in[z][i][j] *= -1.0;
	}
      }

      if (mri_has(refds, "missing") && 
	  !strcmp(mri_get_string(refds,"missing"),"[chunk]")) {
	mri_read_chunk(refds, "missing", c.nslices, 0, 
		       MRI_UNSIGNED_CHAR, c.ref_missing);
      }

      mri_close_dataset(refds);
    }
}

void
MasterResult (int task_number)
{
  static int results= 0;
  int i;

  if (r.sampim_valid)
    memcpy(&c.sampim[r.slice_num][0][0],
	   &r.slice_sampim[0][0],
	   c.os_res * c.os_res * sizeof(float));

  if (r.recon_failed) {
    if (r.image_num >= 0) {
      if (missing) missing[r.image_num][r.slice_num]= 1;
    }
    else {
      Error("# construction of reference failed for slice %d!\n",
	    r.slice_num);
    }
  }

  ++results;
  if (c.verbosity <  VERBOSITY_TIMESTEP) { /* be silent if others reporting */
    if (!(results % 60) || (results == c.nslices*c.total_images))
      Report("# %ld\n", (long)results);
    else
      Report("#");
  }
}

static void
UncompressFile (Filename out,
		const Filename in)
{
  static int count = 0;
  char com[2*sizeof(Filename)+128];

  if (c.verbosity >= VERBOSITY_FILE)
    Report("Uncompressing raw file %s\n", in);
  sprintf(out, "%s/sgrid.%d.%.3d", c.tmp_directory, getpid(), count++);
  if (in[strlen(in)-1] == 'Z')
    sprintf(com, "uncompress -c %s > %s", in, out);
  else
    sprintf(com, "zcat %s > %s", in, out);
  system(com);
#ifdef AFS
  if (CheckForAFS(out)) FlushAFS(out);
#endif
}

static void
RemoveFile (const Filename name)
{
  char com[1024];

  sprintf(com, "rm -f %s", name);
  system(com);
}

#ifdef AFS

static int MySystem (const char *command) {
  int pid, status;
  extern char **environ;
  extern int errno;
  
  if (command == 0)
    return 1;
  pid = fork();
  if (pid == -1)
    return -1;
  if (pid == 0) {
    char *argv[4];
    argv[0] = "sh";
    argv[1] = "-c";
    argv[2] = (char*)command;
    argv[3] = 0;
    execve("/bin/sh", argv, environ);
    exit(127);
  }
  do {
    if (waitpid(pid, &status, 0) == -1) {
      if (errno != EINTR)
        return -1;
    } else
      return status;
  } while(1);
}

static int CheckForAFS( const char* fname )
{
  char buf[256];
  int retval;
  snprintf(buf,sizeof(buf),"fs whereis %s >/dev/null 2>&1",(char*)fname);
  retval= MySystem(buf);
  return retval;
}

static void FlushAFS( const char* path )
{
  char buf[256];
  snprintf(buf,sizeof(buf),"fs flushvolume -path %s >/dev/null 2>&1",path);
  MySystem(buf);
}

#endif
