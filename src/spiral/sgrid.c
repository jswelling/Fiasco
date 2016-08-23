/*
 *	sgrid.c		Version 1.1
 *
 *    Image Gridding for Spiral k-space Acquisitions
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
 *     sgrid - split gridding code off from gsp14 and converted to run
 *	     - in parallel using PVM  (Greg Hood, PSC)
 *	v1.0 - bug fixes; changed to write only one sampim file per slice
 *	     -	(Greg Hood, PSC)		
 *	v1.1 - now uses the updated parallel task manager in util/par.c
 *	     -  (Greg Hood, PSC, June 96)
 *	v1.3 - uses Pittsburgh File Format for raw output
 *	     -  (Greg Hood, PSC, Feb 97)
 */

/*
 * TO DO:
 *	Some matrix operations can be optimized here.
 *	The complicated function used in FermiFilter1 could be precomputed
 *          and stored in a table.  FermiFilter1 could then just do a
 *	    table lookup.
 *	A few more globals could possibly be turned into local variables.
 *	Make the terminology more consistent w.r.t. "phases" vs. "images".
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#ifdef AFS
#include <sys/types.h>
#include <sys/wait.h>
#include <errno.h>
#endif
#include "par.h"
#include "bio.h"
#include "array.h"
#include "misc.h"
#include "acct.h"
#include "mri.h"
#include "fmri.h"
#include "rdb.h"
#include "stdcrg.h"
#include "rttraj.h"

#define NWEIGHTS	1024	/* # of times kaiser function is sampled */
#define PI	M_PI

#define SWAP(a,b) {tempr=(a);(a)=(b);(b)=tempr;}


typedef struct Context {
  /* NOTE: any new fields added to this struct should
     be also added to PackContext and UnpackContext */ 
  Filename input_directory;	/* where to find input files */
  int big_endian_input;		/* TRUE if input is big-endian */
  Filename tmp_directory;	/* where to place scratch files */
  Filename output_directory;	/* where to place output files */
  int big_endian_output;	/* TRUE if output if big-endian */

  int nfiles;			/* # of input files */
  int ncoils;			/* # of coils */
  int nslices;			/* # of slices per coil */
  int nimages;			/* # of images in data set (per coil) */
  int nph1;			/* # of baselines per slice */
  int nphmult;			/* # of images (phases) per baseline */
  int npr;			/* # of projections (spirals) */
  int ndat;			/* # of data items per projection */

  int chop;			/* may be either 1 or -1 */
  int risetime;			/* ??? */
  int densamp;			/* ??? */

  double ts;			/* ??? */
  double gts;			/* ??? */ 
  double fsgcm;			/* ??? */
  double opfov;			/* ??? field-of-view */

  float pix_shifth;		/* horizontal pixel shift */
  float pix_shiftv;		/* vertical pixel shift */

  int coil_record_length;	/* # of bytes in one coil record */
  int slice_record_length;	/* # of bytes in one slice record */
  int baseline_length;		/* # of bytes in one baseline record */

  int res;			/* # of pixels along the X and Y dimensions */
  int os_res;			/* # of pixels in the (possibly oversampled)
				   input image ( =res*over_samp) */

  int slice;			/* slice number to work on (0 indicates all slices) */
  int samp_delay;		/* input delay expressed as # of samples */
  int samp_cor;			/* 1 if sample density correction is to be done,
				   0 if no correction is to be done */
  float ph_twist;		/* phase twist */
  int lr_shift;			/* left-right shift in pixels */
  int tb_shift;			/* top-bottom shift in pixels */
  int loc_shift;		/* 1 if location shift should be done
				     to align slice according to information
				     contained in the file header
				   0 if no shift should be done */
  float zoom;			/* zoom factor */
  float mag_factor;		/* magnitude correction factor */
  float ph_factor;		/* phase correction factor */

  int start_slice;		/* starting slice number */
  int end_slice;		/* ending slice number */

  int lin_cor;			/* linear correction */

  Filename reg_file;		/* registration file */
  int reg_2x;			/* double translation in registration file */

  int write_mag;		/* 1 if we should write magnitude files out
				     as well as raw files;
				   0 if we should only write raw files */

  int over_samp;		/* oversampling ratio */
  float grid_len;		/* grid length (half width of convolution) */
  float gridb;
  float wind;			/* window */

  float factxx;			/* image rotation & scaling factors */
  float factxy;
  float factyx;
  float factyy;

/* variably-sized context info */
  float ***t2k;			/* dimensioned as [2][npr][ndat] */
  float *kdens;			/* dimensioned as [ndat] */
} Context;

typedef struct Task {
  /* NOTE: any new fields added to this struct should
     be also added to PackTask and UnpackTask */ 
  int file_index;	/* counts which input file we are
			   currently working on (0 => first) */
  Filename filename;	/* input filename */

  int coil_num;		/* coil number (first is 0) */
  int slice_num;	/* slice number (first is 0) */
  int image_num;	/* image number within a slice (first is 0) (corresponds to phase number in gsp14.c) */

  float reg_xs;		/* the registration x-shift for this image */
  float reg_ys;		/* the registration y-shift for this image */
  float reg_rot;	/* the registration rotation for this image */
} Task;

static char rcsid[] = "$Id: sgrid.c,v 1.17 2005/07/07 20:04:37 welling Exp $";

/* GLOBAL VARIABLE FOR MASTER & SLAVE */
Task t;
Context c;

/* GLOBAL VARIABLES FOR MASTER */
float **regxs = NULL;	/* [nslices][nimages]	registration x shift */
float **regys = NULL;	/* [nslices][nimages]	registration y shift */
float **regrot = NULL;	/* [nslices][nimages]	registration rotation */
float x_size;		/* the physical size of a voxel in mm */
float y_size;
float z_size;
float x_spacing;	/* the physical spacing between adjacent slices */
float y_spacing;
float z_spacing;

/* GLOBAL VARIABLES FOR SLAVE */
MRI_Dataset *rds = NULL;	/* the raw dataset to write into */
float ***prd = NULL;		/* [2][npr][ndat]	projection data (input) */
float ***grim = NULL;		/* [2][os_res][os_res]	gridded image data (output) */
float **ws = NULL;		/* [os_res][os_res]	weighting array  */
short int *in_buf = NULL;	/* [2*ndat]		input buffer for projection data */
float refl[3];			/* reference location */
float weight[NWEIGHTS+1]; 	/* array which holds a sampled kaiser function -- it is used
				   so we can do table-lookup instead of computing the kaiser
				   function over and over */
float maxk;			/* maximum magnitude in the t2k array */


/* FORWARD DECLARATIONS */
void MasterTask (const int argc, const char **argv, const char **envp);
void ReadEnvironment ();
void ProcessFile (const Filename input_file);
void ReadFirstFileHeader (const Filename input_file);
void CheckHeaderInfo (const Filename input_file);
void ComputeCalibrationMap ();
void LoadRegistrationData ();
void CalculateLocation (float *p1, float *p2, float *p3, int rot, int trans);
void CreateOutputDataset ();
void SlaveFinalize();
void SlaveContext ();
void SlaveTask ();
void LoadLinearReferenceData ();
void LoadProjections ();
void Refocus ();
void FixViews ();
void PerformGridding ();
void FermiFilter1 ();
void FermiFilter2 ();
void WriteRaw ();
void WriteSampleInfo (float **sampim);
void WriteMagnitude ();
void UncompressFile (Filename out, const Filename in);
void RemoveFile (const Filename name);
float kaiser (float l, float b, float u);
float bessi0 (float x);
int IntBRdFloat32 (unsigned char *addr);
void PackContext();
void UnpackContext();
void PackTask();
void UnpackTask();
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
  typedef struct fname_list_struct {
    Filename fname;
    struct fname_list_struct* next;
  } FnameList;
  FnameList* thisname;
  FnameList* fname_list_root= NULL;
  FnameList* fname_list_tail= NULL;
  int i;
  int len;
  Filename input_file;
  Filename input_fullpath;
  int file_count= 0;

  ReadEnvironment();

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

  if (cl_present("2")) c.reg_2x= TRUE; /* -2 flag */
  if (cl_present("bei")) c.big_endian_input= TRUE;
  if (cl_present("beo")) c.big_endian_output= TRUE;
  if (cl_present("lei")) c.big_endian_input= FALSE;
  if (cl_present("leo")) c.big_endian_output= FALSE;
  if (cl_present("l")) c.loc_shift= TRUE;
  if (cl_present("m")) c.write_mag= TRUE;
  cl_get("r","%option %s[%]", "", &(c.reg_file));
  if (cl_present("v")) verbose= TRUE;

  /* Bad karma from previous incarnations forces us to build a
   * linked list of the file names to be dealt with, and then shortly
   * throw it away.
   */
  file_count= 0;
  while (cl_get( "", "%s", input_file )) {
    file_count++;
    ConstructFilename(input_fullpath, c.input_directory, input_file);
    if (!(thisname= (FnameList*)malloc(sizeof(FnameList)))) {
      Abort("Unable to allocate %d bytes for file name parsing!",
	    sizeof(FnameList));
    }
    strcpy(thisname->fname, input_fullpath);
    thisname->next= NULL;
    if (!fname_list_root) fname_list_root= thisname;
    if (fname_list_tail) fname_list_tail->next= thisname;
    fname_list_tail= thisname;
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

  if (file_count==0) {
    fprintf(stderr,"%s: no file argument(s) supplied.\n",argv[0]);
    Help("usage");
    exit(-1);
  }

  c.nfiles = file_count;

  thisname= fname_list_root;
  for (t.file_index = 0; t.file_index < c.nfiles; ++t.file_index)
    {
      Filename tname;
      int uncomp_flag= 0;
      if (((len= strlen(thisname->fname)) > 2)
	  && ((!strcmp(&((thisname->fname)[len-2]), ".Z"))
	      || (!strcmp(&((thisname->fname)[len-2]),".z")))) {
	UncompressFile(tname, thisname->fname);
	uncomp_flag= 1;
      }
      else uncomp_flag= 0;
      ProcessFile(uncomp_flag ? tname : thisname->fname);
      if (uncomp_flag)
	{
	  while (par_tasks_outstanding() > 0)
	    par_wait(1.0);
	  RemoveFile(tname);
	}
      thisname= thisname->next;
    }

  par_finish();
  /* Clean up files and memory */
  thisname= fname_list_root;
  while (thisname) {
    FnameList* target= thisname;
    thisname= thisname->next;
    free(target);
  }
  fname_list_root= fname_list_tail= NULL;
  Report("Done!!\n");
}

void
ReadEnvironment ()
{
  char *fn;
  char *hosts;
  char *dir;

  c.res = GetEnvInt("F_SGRID_RESOLUTION");
  c.slice = GetEnvInt("F_SGRID_SLICE");
  c.samp_delay = GetEnvInt("F_SGRID_SAMP_DELAY");
  c.samp_cor = GetEnvInt("F_SGRID_SAMP_COR");
  c.ph_twist = GetEnvFloat("F_SGRID_PH_TWIST");
  c.lr_shift = GetEnvInt("F_SGRID_LR_SHIFT");
  c.tb_shift = GetEnvInt("F_SGRID_TB_SHIFT");
  c.zoom = GetEnvFloat("F_SGRID_ZOOM");
  c.mag_factor = GetEnvFloat("F_SGRID_MAG_FACTOR");
  c.ph_factor = GetEnvFloat("F_SGRID_PH_FACTOR");
  c.lin_cor = GetEnvInt("F_SGRID_LIN_COR");

  if ((fn = getenv("F_SGRID_REG_FILE")) != NULL)
    StringCopy(c.reg_file, fn, MAX_FILENAME_LENGTH);
  else
    c.reg_file[0] = '\0';
  
  c.reg_2x = GetEnvInt("F_SGRID_REG_2X");
  c.over_samp = GetEnvInt("F_SGRID_OVER_SAMP");
  c.grid_len = GetEnvFloat("F_SGRID_GRID_LEN");
  c.loc_shift = GetEnvInt("F_SGRID_LOC_SHIFT");
  c.write_mag = GetEnvInt("F_SGRID_WRITE_MAG");

  if ((hosts = getenv("F_SGRID_HOSTS")) != NULL)
    par_set_hosts(hosts);

  strcpy(c.input_directory, ".");
  if ((dir = getenv("F_SGRID_INPUT_DIR")) != NULL)
    StringCopy(c.input_directory, dir, sizeof(Filename));

  strcpy(c.tmp_directory, "/tmp");
  if ((dir = getenv("F_SGRID_TMP_DIR")) != NULL)
    StringCopy(c.tmp_directory, dir, sizeof(Filename));

  strcpy(c.output_directory, ".");
  if ((dir = getenv("F_SGRID_OUTPUT_DIR")) != NULL)
    StringCopy(c.output_directory, dir, sizeof(Filename));

  /* initially set the variably-sized arrays to be
     unallocated */
  c.t2k = NULL;
  c.kdens = NULL;

  /* the endianness of the input defaults to TRUE
     since the P files are generated in big-endian format */
  c.big_endian_input = TRUE;
  /* the endianness of the output defaults to the
     native format of the machine on which we are running */
  c.big_endian_output = bio_big_endian_machine;
}

void
ProcessFile (const Filename input_file)
{
  static int first_file = TRUE;

  /* read what the master needs to know
     from the the input_file */
  if (first_file)
    {
      ReadFirstFileHeader(input_file);
      CreateOutputDataset();
      ComputeCalibrationMap();
      LoadRegistrationData();
      par_set_context();
      first_file = FALSE;
    }
  else
    CheckHeaderInfo(input_file);

  /* dispatch most of the work to the
     slaves */
  strcpy(t.filename, input_file);
  for (t.coil_num = 0; t.coil_num < c.ncoils; ++t.coil_num)
    for (t.slice_num = c.start_slice; t.slice_num <= c.end_slice; ++t.slice_num)
      for (t.image_num = 0; t.image_num < c.nimages; ++t.image_num)
	{
	  if (c.reg_file[0] != '\0')
	    {
	      t.reg_xs = regxs[t.slice_num][(t.file_index*c.nimages)
					   +t.image_num];
	      t.reg_ys = regys[t.slice_num][(t.file_index*c.nimages)
					   +t.image_num];
	      t.reg_rot = regrot[t.slice_num][(t.file_index*c.nimages)
					     +t.image_num];
	    }
	  else
	    t.reg_xs = t.reg_ys = t.reg_rot = 0.0;
	  par_delegate_task();
	}
}

void
ReadFirstFileHeader (const Filename input_file)
{
  FILE *f;
  unsigned char h[RDB_HDR_SIZE];
  float s_i_offset;
  float tempr2,tempr3;
  float p1[3],p2[3],p3[3];
  float ctr[3];
  float n2[3],n3[3];
  int i, k, n;
  float tmp1, tmp2;
  short rotation, transpose;

  Acct(READOPEN);
  /* open the file */
#ifdef AFS
  if (CheckForAFS(input_file)) FlushAFS(input_file);
#endif
  if ((f = fopen(input_file, "r")) == NULL)
    Abort("Can't open %s!\n", input_file);
  
  /* read the header */
  Acct(READING);
  if (fread(h, RDB_HDR_SIZE, 1, f) != 1)
    Abort("Can't read header of %s\n", input_file);
  Acct(READCLOSE);
  fclose(f);
  Acct(PROCESSING);

  /* extract the context info */
  bio_big_endian_input = c.big_endian_input;
  c.nph1 = IntBRdFloat32(&h[RDB_HDR_USER4]);
  c.nphmult = IntBRdFloat32(&h[RDB_HDR_USER10]);
  if (c.nphmult < 1)
    c.nphmult = 1;
  c.nimages = c.nph1 * c.nphmult;
  c.nslices = BRdInt16(&h[RDB_HDR_NSLICES]) / c.nph1;
  c.ndat = BRdInt16(&h[RDB_HDR_FRAME_SIZE]);
  c.npr = IntBRdFloat32(&h[RDB_HDR_USER5]);
  c.chop = IntBRdFloat32(&h[RDB_HDR_USER7]);
  c.risetime = IntBRdFloat32(&h[RDB_HDR_USER11]);
  c.densamp = IntBRdFloat32(&h[RDB_HDR_USER0]);
  c.fsgcm = BRdFloat32(&h[RDB_HDR_USER12]);
  c.opfov = BRdFloat32(&h[RDB_HDR_USER6]);
  c.ts = BRdFloat32(&h[RDB_HDR_USER2]) * 1e-6;
  c.gts = BRdFloat32(&h[RDB_HDR_USER3]) * 1e-6;
  c.ncoils = BRdInt16(&h[RDB_HDR_DAB_0_STOP_RCV]) -
    BRdInt16(&h[RDB_HDR_DAB_0_START_RCV]) + 1;
  c.coil_record_length = BRdInt32(&h[RDB_HDR_RAW_PASS_SIZE]) / c.ncoils;
  c.baseline_length = c.ndat*4;
  c.slice_record_length = c.nimages*c.npr*c.ndat*4 + c.nph1*c.baseline_length;
  if (((c.nphmult * c.npr) % 2) == 1)
    c.slice_record_length += c.nph1*c.ndat*4;
  c.os_res = c.over_samp * c.res;

  Report("number of time points = %d\n", c.nimages);
  Report("number of slices = %d\n", c.nslices);
  Report("number of samples in readout = %d\n", c.ndat);
  Report("field of view (mm) = %.1f\n", c.opfov);
  Report("number of spirals = %d\n", c.npr);
  Report("time of scan = %s\n",&h[RDB_HDR_SCAN_TIME]);
  Report("date of scan = %s\n",&h[RDB_HDR_SCAN_DATE]);

  x_size = c.opfov / c.res;
  y_size = c.opfov / c.res;
  z_size = BRdFloat32(&h[RDB_HDR_IM_SLTHICK]);

  x_spacing = x_size;
  y_spacing = y_size;
  z_spacing = BRdFloat32(&h[RDB_HDR_IM_SCANSPACING]);

  ctr[0] = BRdFloat32(&h[RDB_HDR_IM_CTR_R]);
  ctr[1] = BRdFloat32(&h[RDB_HDR_IM_CTR_A]);
  ctr[2] = BRdFloat32(&h[RDB_HDR_IM_CTR_S]);

  for (i = 0; i < 3; ++i)
    {
      p1[i] = BRdFloat32(&h[RDB_HDR_GW_POINT1_0 + i*RDB_HDR_GW_POINT_STEP]);
      p2[i] = BRdFloat32(&h[RDB_HDR_GW_POINT2_0 + i*RDB_HDR_GW_POINT_STEP]);
      p3[i] = BRdFloat32(&h[RDB_HDR_GW_POINT3_0 + i*RDB_HDR_GW_POINT_STEP]);
    }
  rotation = BRdInt16(&h[RDB_HDR_ROTATION]);
  transpose = BRdInt16(&h[RDB_HDR_TRANSPOSE]);
  CalculateLocation(p1, p2, p3, rotation, transpose);
  s_i_offset = BRdFloat32(&h[RDB_HDR_IM_TLC_S]) - p1[2];
  ctr[2] -= s_i_offset;

  /* calculate in-plane normal vectors n2,n3 */
  tempr2 = tempr3 = 0.0;
  for (i = 0; i < 3; i++)
    {
      n2[i] = p2[i] - p1[i];
      tempr2 += n2[i]*n2[i];
      n3[i] = p2[i] - p3[i];
      tempr3 += n3[i]*n3[i];
    }
  tempr2 = sqrt(tempr2);
  tempr3 = sqrt(tempr3);
  for (i = 0; i < 3; i++)
    {
      n2[i] /= tempr2;
      n3[i] /= tempr3;
    }
  if (c.loc_shift)
    {
      c.pix_shifth = -(n2[0]*ctr[0] +  n2[1]*ctr[1] +  n2[2]*ctr[2])/tempr2;
      c.pix_shiftv = (n3[0]*ctr[0] +  n3[1]*ctr[1] +  n3[2]*ctr[2])/tempr2;
      Report("  shifth = %f, shiftv = %f (in mm)\n", c.pix_shifth*tempr2, c.pix_shiftv*tempr2);
      for (n = 0; n < c.nslices; n++)
	{
	  for (i = 0; i < 3; i++)
	    {
	      p1[i] = BRdFloat32(&h[RDB_HDR_GW_POINT1_0 + n*RDB_HDR_A_STEP + i*RDB_HDR_GW_POINT_STEP]);
	      p2[i] = BRdFloat32(&h[RDB_HDR_GW_POINT2_0 + n*RDB_HDR_A_STEP + i*RDB_HDR_GW_POINT_STEP]);
	      p3[i] = BRdFloat32(&h[RDB_HDR_GW_POINT3_0 + n*RDB_HDR_A_STEP + i*RDB_HDR_GW_POINT_STEP]);
	    }
	  CalculateLocation(p1, p2, p3, rotation, transpose);
	  Report("slice %d    R/L       A/P       S/I\n",n+1);
	  Report(" TL:  %f, %f, %f\n",p1[0],p1[1],p1[2]+s_i_offset);
	  Report(" TR:  %f, %f, %f\n",p2[0],p2[1],p2[2]+s_i_offset);
	  Report(" BR:  %f, %f, %f\n",p3[0],p3[1],p3[2]+s_i_offset);
	}    
    }
  else
    c.pix_shifth = c.pix_shiftv = 0.0;

  /* compute the rotation & scaling factors */
  c.factxx = -1.0/c.zoom; c.factxy = 0.0; c.factyx = 0.0; c.factyy = -1.0/c.zoom; 
  if (transpose)
    {
      tmp1 = c.factyx; tmp2 = c.factyy;
      c.factyx = c.factxx; c.factyy = c.factxy;
      c.factxx = tmp1; c.factxy = tmp2;
    }
  for (k = 0; k < rotation; k++)
    {
      tmp1 = c.factyx; tmp2 = c.factyy;
      c.factyx = -c.factxx; c.factyy = -c.factxy;
      c.factxx = tmp1; c.factxy = tmp2;
    }

  /* allocate the t2k and kdens arrays */
  c.t2k = Alloc3DFloatArray(2, c.npr, c.ndat);
  c.kdens = (float *) malloc(c.ndat * sizeof(float));

  /* set the starting and ending slices */
  if (c.slice > 0)
    c.start_slice = c.end_slice = c.slice - 1;
  else
    {
      c.start_slice = 0;
      c.end_slice = c.nslices - 1;
    }
}

void
CheckHeaderInfo (const Filename input_file)
{
  FILE *f;
  unsigned char h[RDB_HDR_SIZE];

  Acct(READOPEN);
  /* open the file */
  if ((f = fopen(input_file, "r")) == NULL)
    Abort("Can't open %s!\n", input_file);
  
  /* read the header */
  Acct(READING);
  if (fread(h, RDB_HDR_SIZE, 1, f) != 1)
    Abort("Can't read header of %s\n", input_file);
  Acct(READCLOSE);
  fclose(f);
  Acct(PROCESSING);

  /* check that the context info is the same */
  bio_big_endian_input = c.big_endian_input;
  if (c.nph1 != IntBRdFloat32(&h[RDB_HDR_USER4]) ||
      c.nphmult != IntBRdFloat32(&h[RDB_HDR_USER10]) && (c.nphmult > 1 || IntBRdFloat32(&h[RDB_HDR_USER10]) > 1) ||
      c.nimages != c.nph1 * c.nphmult ||
      c.nslices != BRdInt16(&h[RDB_HDR_NSLICES]) / c.nph1 ||
      c.ndat != BRdInt16(&h[RDB_HDR_FRAME_SIZE]) ||
      c.npr != IntBRdFloat32(&h[RDB_HDR_USER5]) ||
      c.chop != IntBRdFloat32(&h[RDB_HDR_USER7]) ||
      c.risetime != IntBRdFloat32(&h[RDB_HDR_USER11]) ||
      c.densamp != IntBRdFloat32(&h[RDB_HDR_USER0]) ||
      c.fsgcm != BRdFloat32(&h[RDB_HDR_USER12]) ||
      c.opfov != BRdFloat32(&h[RDB_HDR_USER6]) ||
      c.ts != BRdFloat32(&h[RDB_HDR_USER2])*1e-6 ||
      c.gts != BRdFloat32(&h[RDB_HDR_USER3])*1e-6 ||
      c.ncoils != BRdInt16(&h[RDB_HDR_DAB_0_STOP_RCV]) - BRdInt16(&h[RDB_HDR_DAB_0_START_RCV]) + 1 ||
      c.coil_record_length != BRdInt32(&h[RDB_HDR_RAW_PASS_SIZE]) / c.ncoils)
    Abort("The headers of the input files are inconsistent.\n");
}


void
ComputeCalibrationMap ()
{
  int i, j, n;
  float kx,ky,kxo,kyo,ggx,ggy,ggm,kksr,kksi;
  FILE *f;
  Filename fn;
  int nsamp;
  int nsamp_generated;
  int res;

  /* kaiser-bessel functions - initialize terms */
  c.gridb = PI*c.grid_len;
  c.wind = c.grid_len;

  nsamp_generated = 
    getrttrajvd(c.ndat, c.npr, c.ts, c.gts, 0.0, c.fsgcm, c.opfov,
		c.risetime, c.densamp, 1.0, 1.0, c.t2k[0], c.t2k[1], &res);
  Report("acq. matrix = %d, nom. resolution = %f mm\n", res, c.opfov/res);

  nsamp = c.ndat;

  if (c.samp_cor != 1) 
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

void
LoadRegistrationData ()
{
  FILE *f;
  int ss;
  int pp;
  float jj;
  Filename fn;

  /* check if a registration file was specified */
  if (c.reg_file[0] != '\0')
    {
      /* allocate the arrays to hold the registration data */
      regxs = Alloc2DFloatArray(c.nslices, c.nimages * c.nfiles);
      regys = Alloc2DFloatArray(c.nslices, c.nimages * c.nfiles);
      regrot = Alloc2DFloatArray(c.nslices, c.nimages * c.nfiles);

      Acct(READOPEN);
      if ((f = fopen(c.reg_file, "r")) == NULL)
	Abort("Can't open %s!\n", fn);
      Acct(READING);
      while (fscanf(f,"%d %d",&pp,&ss) == 2)
	{
	  if (ss < 0 || ss >= c.nslices ||
	      pp < 0 || pp >= c.nimages*c.nfiles)
	    Abort("Corrupt registration data: %d (max %d) %d (max %d)\n",
		  ss, c.nslices, pp, c.nimages);
	  if (fscanf(f,"%f %f %f %f",
		     &regxs[ss][pp],
		     &regys[ss][pp],
		     &regrot[ss][pp],
		     &jj) != 4)
	    Abort("Bad registration data\n");
	  if (c.reg_2x)
	    {
	      regxs[ss][pp] *= 2.0;
	      regys[ss][pp] *= 2.0;
	    }
	}
      Acct(READCLOSE);
      fclose(f); 
      Acct(PROCESSING);
    }
}

void
CalculateLocation (float *p1,
		   float *p2,
		   float *p3,
		   int rot,
		   int trans)
{
  /* 
    top left     - p1
    bottom left  - p2
    top right    - p3
    bottom right - p4
  */
  float p4[3], tempr;
  int i,n;

  /* calc missing corner (bottom right) */
  for (i=0; i<3; i++)
    p4[i] = p3[i] + (p2[i] - p1[i]);

  if (trans)
    for (i=0; i<3; i++)
      SWAP(p2[i],p3[i]);

  for (n=0; n<rot; n++)
    for (i=0; i<3; i++)
      {
	SWAP(p1[i],p3[i]);
	SWAP(p3[i],p4[i]);
	SWAP(p4[i],p2[i]);
      }

  /* reformat to match image header */
  for (i=0; i<3; i++)
    {
      SWAP(p2[i],p3[i]);
      SWAP(p3[i],p4[i]);
    }
}

void
CreateOutputDataset ()
{
  Filename dsn;
  MRI_Dataset *ds;

  sprintf(dsn, "%s/raw.mri", c.output_directory);
  Acct(WRITEOPEN);
  ds = mri_open_dataset(dsn, MRI_WRITE);
  hist_add(ds,"*sgrid*");
  Acct(WRITING);

  mri_create_chunk(ds, "images");
  mri_set_string(ds, "images.datatype", "float32");
  mri_set_string(ds, "images.dimensions", "xyvztc");
  mri_set_int(ds, "images.extent.x", c.os_res);
  mri_set_int(ds, "images.extent.y", c.os_res);
  mri_set_int(ds, "images.extent.v", 2);
  mri_set_int(ds, "images.extent.z", c.nslices);
  mri_set_int(ds, "images.extent.t", c.nimages * c.nfiles);
  mri_set_int(ds, "images.extent.c", c.ncoils);
  mri_set_float(ds, "images.size.x", x_size);
  mri_set_float(ds, "images.size.y", y_size);
  mri_set_float(ds, "images.size.z", z_size);
  mri_set_float(ds, "images.spacing.x", x_spacing);
  mri_set_float(ds, "images.spacing.y", y_spacing);
  mri_set_float(ds, "images.spacing.z", z_spacing);
  mri_set_string(ds, "images.file", ".dat");

  mri_create_chunk(ds, "sampim");
  mri_set_string(ds, "sampim.datatype", "float32");
  mri_set_string(ds, "sampim.dimensions", "xyz");
  mri_set_int(ds, "sampim.extent.x", c.os_res);
  mri_set_int(ds, "sampim.extent.y", c.os_res);
  mri_set_int(ds, "sampim.extent.z", c.nslices);
  mri_set_string(ds, "sampim.file", ".sam");

  mri_set_int(ds, "npr", c.npr);
  mri_set_int(ds, "ndat", c.ndat);
  mri_set_int(ds, "res", c.res);
  mri_set_int(ds, "over_samp", c.over_samp);
  mri_set_int(ds, "grid_len", c.grid_len);

  Acct(WRITECLOSE);
  mri_close_dataset(ds);
  Acct(PROCESSING);
}


/* SLAVE PROCEDURES */

void
SlaveFinalize()
{
  if (rds != NULL)
    {
      Acct(WRITECLOSE);
      mri_close_dataset(rds);
    }
  /* deallocate the slave arrays */
  Free3DFloatArray(prd);
  Free3DFloatArray(grim);
  Free2DFloatArray(ws);
  if (in_buf != NULL)
    free(in_buf);
}

void
SlaveContext ()
{
  int i;
  Filename rds_name;

  /* reopen the raw dataset */
  if (rds != NULL)
    {
      Acct(WRITECLOSE);
      mri_close_dataset(rds);
    }
  Acct(WRITEOPEN);
  sprintf(rds_name, "%s/raw.mri", c.output_directory);
  /* open the output raw dataset */
  rds = mri_open_dataset(rds_name, MRI_MODIFY_DATA);

  Acct(PROCESSING);
  /* reallocate the slave arrays */
  Free3DFloatArray(prd);
  Free3DFloatArray(grim);
  Free2DFloatArray(ws);
  if (in_buf != NULL)
    free(in_buf);

  prd = Alloc3DFloatArray(2, c.npr, c.ndat);
  grim = Alloc3DFloatArray(2, c.os_res, c.os_res);
  ws = Alloc2DFloatArray(c.os_res, c.os_res);
  in_buf = (short int *) malloc(2*sizeof(short int)*c.ndat);

  /* construct the convolution function table --- kaiser-bessel */
  for (i = 0; i < NWEIGHTS+1; i++)
    weight[i] = kaiser((float)NWEIGHTS, c.gridb, (float)(i-NWEIGHTS/2));
}

void
SlaveTask ()
{
  Filename fn;

  LoadLinearReferenceData();
  Report("Co %d, Sl %2d, Im %2d: Ld.", t.coil_num+1, t.slice_num, t.image_num+1);
  LoadProjections();
  Refocus();
  FixViews();
  Report("Rm.");
  PerformGridding();
  Report("Fl.");
  FermiFilter1();
  FermiFilter2();
  WriteRaw();
  if (c.write_mag)
    WriteMagnitude();
  Report("\n");
}

void
LoadLinearReferenceData ()
{
  Filename fn;
  FILE *f;

  if (c.lin_cor)
    {
      sprintf(fn, "%s/refl.s%.2d", c.input_directory, t.slice_num+1);
      Acct(READOPEN);
      if ((f = fopen(fn, "r")) == NULL)
	Abort("Can't open reference file %s!\n", fn);
      Acct(READING);
      bio_error = FALSE;
      bio_big_endian_input = c.big_endian_input;
      FRdFloat32Array(f, refl, 3);
      if (bio_error)
	Abort("Could not read reference file %s\n", fn);
      Acct(READCLOSE);
      fclose(f);
      Acct(PROCESSING);
    }
  else
    refl[0] = refl[1] = refl[2] = 0.0;
}

void
LoadProjections ()
{
  int n;
  int j;
  int chopper = 1;
  FILE *f;
  int slice_sequence_number;
  long offset;

  Acct(READOPEN);
#ifdef AFS
  if (CheckForAFS(t.filename)) FlushAFS(t.filename);
#endif
  if ((f = fopen(t.filename, "r")) == NULL)
    Abort("Could not open input file %s\n", t.filename);

  Acct(READING);
  /* all odd slices are located before all even slices
     in the P file (assuming a numbering scheme that
     begins with 1); hence we need to compute a sequence
     number in order to pull the right slice out of the
     file; note that our internal slice numbers begin
     with 0) */
  if (t.slice_num % 2 == 0)
    slice_sequence_number = t.slice_num / 2;
  else
    slice_sequence_number = (c.nslices + t.slice_num) / 2;

  offset = RDB_HDR_SIZE +
    t.coil_num * c.coil_record_length +
      slice_sequence_number * c.slice_record_length +
	(t.image_num / c.nphmult + 1) * c.baseline_length +
	  t.image_num * c.npr * c.ndat * 4;
  if (((c.nphmult * c.npr) % 2) == 1)
    offset += (t.image_num / c.nphmult) * c.ndat * 4;
  if (fseek(f, offset, SEEK_SET) != 0)
    Abort("Could not load projections\n");

  bio_error = FALSE;
  bio_big_endian_input = c.big_endian_input;
  for (n = 0; n < c.npr; ++n)
    {
      FRdInt16Array(f, in_buf, 2*c.ndat);
      if (bio_error)
	Abort("Could not load projections\n");
      for (j = 0; j < c.ndat; ++j)
	{
	  prd[0][n][j] = in_buf[2*j] * chopper;
	  prd[1][n][j] = in_buf[2*j+1] * chopper;
	}
      chopper *= c.chop;
    }
  Acct(READCLOSE);
  fclose(f);
  Acct(PROCESSING);
}

void
Refocus ()
{
  int i, j;
  float theta;
  float tr, ti, cs, sn;

  for (j=0; j<c.ndat; j++)
    {
      theta = (-2.0*PI*c.ph_twist/c.ndat - refl[0])*j;
      cs = cos(theta);
      sn = sin(theta);
      for (i=0; i<c.npr; i++)
	{
	  tr = prd[0][i][j];
	  ti = prd[1][i][j];
	  prd[0][i][j] = tr*cs - ti*sn;
	  prd[1][i][j] = tr*sn + ti*cs;
	}
    }
}

void
FixViews ()
{
  int i, j;
  float p0r, p0i, prr, pri, cc, dd, tmp;
  FILE *warn, *f;
  float cs, sn, theta;
  long offset;
  int slice_sequence_number;
  float *tmp_p0, *tmp_p1, *tmp_r0, *tmp_r1;
  short *tmp_s;

  p0r = p0i = 0.0;
  for (i=0; i<c.npr; i++)
    for (j = 0; j <= c.samp_delay; j++)
      { 
	p0r += prd[0][i][j];
	p0i += prd[1][i][j];
      }
  p0r /= c.npr;
  p0i /= c.npr;

  for (i = 0; i < c.npr; i++)
    {
      prr = pri = 0.0;
      for (j = 0; j <= c.samp_delay; j++)
	{ 
	  prr += prd[0][i][j];
	  pri += prd[1][i][j];
	}

      /* phase correction angle */
      tmp = atan2(p0i,p0r) - atan2(pri,prr); 
      if (tmp >  PI)
	tmp -= 2.0*PI;
      if (tmp < -PI)
	tmp += 2.0*PI;
      cc = cos(tmp * c.ph_factor);
      dd = sin(tmp * c.ph_factor);

      /* magnitude correction term   */
      if (prr == 0.0 && pri == 0.0)
	{
	  warn = fopen("WARNING", "a");
	  fprintf(warn, "*** WARNING ****  prr == 0.0 and pri == 0.0 in FixViews\n");
	  fprintf(warn, "****************  i=%d image=%d slice=%d filename=%s\n",
		  i, t.image_num, t.slice_num, t.filename);
	  fprintf(warn, "       tmp = %f cc = %f dd = %f\n", tmp, cc, dd);
	  fprintf(warn, "       c.ph_factor = %f\n", c.ph_factor);
	  fprintf(warn, "       prd values\n");
	  for (j = 0; j <= c.samp_delay; ++j)
	    fprintf(warn, "           prd[0][i][%d] = %f  prd[1][i][%d] = %f\n", j, prd[0][i][j], j, prd[1][i][j]);
	  fprintf(warn, "       current in_buf values\n");
	  for (j = 0; j <= c.samp_delay; ++j)
	    fprintf(warn, "           in_buf[2*%d] = %d     in_buf[2*%d+1] = %d\n",
		    j, in_buf[2*j], j, in_buf[2*j+1]);
	  fprintf(warn, "       c.chop = %d\n", c.chop);
	  fprintf(warn, "       recomputed prd\n");
	  tmp_p0 = (float *) malloc((c.samp_delay + 1) * sizeof(float));
	  tmp_p1 = (float *) malloc((c.samp_delay + 1) * sizeof(float));
	  tmp_r0 = (float *) malloc((c.samp_delay + 1) * sizeof(float));
	  tmp_r1 = (float *) malloc((c.samp_delay + 1) * sizeof(float));
	  for (j = 0; j <= c.samp_delay; ++j)
	    {
	      tmp_p0[j] = in_buf[2*j];
	      tmp_p1[j] = in_buf[2*j+1];
	      theta = (-2.0*PI*c.ph_twist/c.ndat - refl[0])*j;
	      cs = cos(theta);
	      sn = sin(theta);
	      tmp_r0[j] = tmp_p0[j]*cs - tmp_p1[j]*sn;
	      tmp_r1[j] = tmp_p0[j]*sn + tmp_p1[j]*cs;
	      fprintf(warn, "   theta=%f cs=%f sn=%f prd[0][i][%d]=%f prd[1][i][%d]=%f\n",
		      theta, cs, sn, j, tmp_r0[j], j, tmp_r1[j]);
	    }
	  free(tmp_p0);
	  free(tmp_p1);
	  free(tmp_r0);
	  free(tmp_r1);
	  fprintf(warn, "       in_buf reread from disk\n");
	  tmp_s = (short *) malloc((c.samp_delay + 1) * sizeof(short));
	  if ((f = fopen(t.filename, "r")) == NULL)
	    fprintf(warn, "Could not open file %s\n", t.filename);
	  else
	    {
	      if (t.slice_num % 2 == 0)
		slice_sequence_number = t.slice_num / 2;
	      else
		slice_sequence_number = (c.nslices + t.slice_num) / 2;
	      offset = RDB_HDR_SIZE +
		t.coil_num * c.coil_record_length +
		slice_sequence_number * c.slice_record_length +	
		(t.image_num / c.nphmult + 1) * c.baseline_length +
		t.image_num * c.npr * c.ndat * 4 +
		i * c.ndat * 4;
	      if (fseek(f, offset, SEEK_SET) != 0)
		fprintf(warn, "could not seek\n");
	      bio_error = FALSE;
	      bio_big_endian_input = c.big_endian_input;
	      fprintf(warn, "  t.coilnum = %d c.coil_record_length = %d slice_sequence_number=%d \n",
		      t.coil_num, c.coil_record_length, slice_sequence_number);
	      fprintf(warn, "  c.slice_record_length = %d  t.image_num = %d  c.nphmult = %d\n",
		      c.slice_record_length, t.image_num, c.nphmult);
	      fprintf(warn, "  c.baseline_length = %d c.npr = %d c.ndat = %d i = %d\n",
		      c.baseline_length, c.npr, c.ndat, i);
	      fprintf(warn, "  offset = %d\n", (int) offset);
	      FRdInt16Array(f, tmp_s, 2*(c.samp_delay+1));
	      if (bio_error)
		fprintf(warn, "bio_error occurred\n");
	      else
		for (j = 0; j <= c.samp_delay; ++j)
		  fprintf(warn, "      in_buf[2*%d] = %d  in_buf[2*%d+1] = %d\n",
			  j, tmp_s[2*j], j, tmp_s[2*j+1]);
	      fclose(f);
	    }
	  free(tmp_s);

	  fclose(warn);
	  tmp = 0.0;
	}
      else
	tmp = hypot(p0r, p0i) / hypot(prr, pri);
      cc *= tmp*c.mag_factor + (1.0 - c.mag_factor);
      dd *= tmp*c.mag_factor + (1.0 - c.mag_factor);
      for (j = 0; j < c.ndat; j++)
	{ 
	  prr = prd[0][i][j];
	  pri = prd[1][i][j];
	  tmp = cc*prr - dd*pri;
	  prd[1][i][j] = dd*prr + cc*pri;
	  prd[0][i][j] = tmp;
	}
    }
}

void
PerformGridding ()
{
  int i, j;
  int lx, ly;
  float kx, ky;
  float w;
  float pr, pi;
  float mkr;
  float dkx,dky;
  float dwin;
  float w2;
  float wx;
  float tmpd;
  float rot1, rot2;
  float rotfact;
  int maxj;
  float **sampim= NULL;      /* computation of this array may someday be
				moved into srecon */

  for (i=0; i<c.os_res; i++)
    for (j=0; j<c.os_res; j++)
      {
	grim[0][i][j] = 0.0;
	grim[1][i][j] = 0.0;
	ws[i][j] = 0.0;
      }
  if (sampim==NULL)
    {
      sampim = Alloc2DFloatArray(c.os_res, c.os_res);
      for (i=0; i<c.os_res; i++)
	for (j=0; j<c.os_res; j++)
	  sampim[i][j] = 0.0;
    }

  maxk = 0.0;
  dwin = c.wind;
  rotfact = 2.0*PI / c.res / c.over_samp;
  maxj = c.ndat - c.samp_delay - 2;
  for (i = 0; i < c.npr; i++)
    for (j=0; j <= maxj; j++)
      { 
	w2 = c.kdens[j];
	kx = c.t2k[0][i][j];
	ky = c.t2k[1][i][j];
	if ( (mkr = hypot(kx,ky)) > maxk)
	  maxk = mkr;
	pr = prd[0][i][j+c.samp_delay];
	pi = prd[1][i][j+c.samp_delay];
	dkx = (c.factxx*kx + c.factxy*ky)*c.over_samp;
	dky = (c.factyx*kx + c.factyy*ky)*c.over_samp;
	if (c.reg_file[0] != '\0')
	  { 
	    rot1 = cos(rotfact*dkx*(c.pix_shifth*c.res - t.reg_xs + c.lr_shift));
	    rot2 = -sin(rotfact*dkx*(c.pix_shifth*c.res - t.reg_xs + c.lr_shift));
	    tmpd = pr*rot1 + pi*rot2;
	    pi = pi*rot1 - pr*rot2;
	    pr = tmpd;

	    rot1 = cos(rotfact*dky*(c.pix_shiftv*c.res - t.reg_ys + c.tb_shift));
	    rot2 = -sin(rotfact*dky*(c.pix_shiftv*c.res - t.reg_ys + c.tb_shift));
	    tmpd = pr*rot1 + pi*rot2;
	    pi = pi*rot1 - pr*rot2;
	    pr = tmpd;

	    rot1 = cos(PI/180.0 * t.reg_rot);
	    rot2 = -sin(PI/180.0 * t.reg_rot);
	    tmpd = dkx*rot1 - dky*rot2;
	    dky = dky*rot1 + dkx*rot2;
	    dkx = tmpd; 
	  }
	else if (c.lr_shift != 0 || c.tb_shift != 0 || c.loc_shift)
	  { 
	    rot1 = cos(rotfact*dkx*(c.pix_shifth*c.res + c.lr_shift));
	    rot2 = -sin(rotfact*dkx*(c.pix_shifth*c.res + c.lr_shift));
	    tmpd = pr*rot1 + pi*rot2;
	    pi = pi*rot1 - pr*rot2;
	    pr = tmpd;

	    rot1 = cos(rotfact*dky*(c.pix_shiftv*c.res + c.tb_shift));
	    rot2 = -sin(rotfact*dky*(c.pix_shiftv*c.res + c.tb_shift));
	    tmpd = pr*rot1 + pi*rot2;
	    pi = pi*rot1 - pr*rot2;
	    pr = tmpd;
	  }
	dkx += (c.res/2 - refl[2]*j)*c.over_samp; 
	dky += (c.res/2 - refl[1]*j)*c.over_samp;
	for (lx = ceil(dkx-dwin); lx<=floor(dkx+dwin); lx++)
	  {
	    if ((lx<0) || (lx>=c.os_res))
	      continue;
	    wx = weight[ Round((((dkx-lx)/dwin)*NWEIGHTS/2 )) + NWEIGHTS/2] *w2;
	    for (ly = ceil(dky-dwin); ly<=floor(dky+dwin); ly++)
	      {
		if ((ly<0) || (ly>=c.os_res))
		  continue;
		w = wx*weight[ Round((((dky-ly)/dwin)*NWEIGHTS/2 )) + NWEIGHTS/2];
		if (t.image_num == 0 && t.file_index == 0 && t.coil_num == 0)
		  sampim[lx][ly] += w * j;
		grim[0][lx][ly] += pr * w;
		grim[1][lx][ly] += pi * w;
		ws[lx][ly] += w;
	      }
	  }
      }
  
  if (t.image_num == 0 && t.file_index == 0 && t.coil_num == 0)
    {
      for (i=0; i<c.os_res; i++)
	for (j=0; j<c.os_res; j++)
	  sampim[i][j] /= (0.001+ws[i][j]);
      WriteSampleInfo(sampim);
      Free2DFloatArray(sampim);
    }
}

void
FermiFilter1 ()
{
  float ww;
  float rad;
  int i, j;
  int offs;

  rad = (c.res / 2) * c.over_samp;
  if (maxk*c.over_samp/c.zoom < rad)
    rad = maxk*c.over_samp/c.zoom;
  offs = (c.res / 2) * c.over_samp;
  for (i = 0; i < c.os_res; ++i)
    for (j = 0; j < c.os_res; ++j)
      {
	ww = 1.0 / (1.0 + exp((hypot((float)(offs-i), (float)(offs-j)) - rad) *
			      6.4 / rad));
	ws[i][j] = ww / (0.001 + ws[i][j]);
      }
}

void
FermiFilter2 ()
{
  int i, j;

  for (i = 0; i < c.os_res; ++i)
    for (j = 0; j < c.os_res; ++j)
      {
	grim[0][i][j] *= ws[i][j];
	grim[1][i][j] *= ws[i][j];
      }
}

void
WriteRaw ()
{
  int overall_offset;

  overall_offset = (((t.coil_num * c.nfiles + t.file_index) * c.nimages +
    t.image_num) * c.nslices + t.slice_num) * 2 * c.os_res * c.os_res;
  
  Acct(WRITING);
  mri_set_chunk(rds, "images", 2*c.os_res*c.os_res, overall_offset,
		MRI_FLOAT, &grim[0][0][0]);
  Acct(PROCESSING);
}

void
WriteSampleInfo (float **sampim)
{
  int overall_offset;

  overall_offset = t.slice_num * c.os_res * c.os_res;
  
  Acct(WRITING);
  mri_set_chunk(rds, "sampim", c.os_res*c.os_res, overall_offset,
		MRI_FLOAT, &sampim[0][0]);
  Acct(PROCESSING);
}

void WriteMagnitude ()
{
  int i, j;
  float imm, imx;
  Filename fn;
  FILE *f;
  int overall_image_num;	/* the image # within the entire set of input files */
  short int **mag;		/* [os_res][os_res]	magnitude of gridded data */

  /* create a temporary array to hold the magnitudes */
  mag = Alloc2DShortArray(c.os_res, c.os_res);

  /* construct the output buffer from the magnitudes of the values
     in the grim array */
  imx = 0.0;
  for (i = 0; i < c.os_res; i++)
    for (j = 0; j < c.os_res; j++)
      {
	imm = hypot(grim[0][j][i], grim[0][j][i]);
	if (imm > imx)
	  imx = imm;
	mag[i][j] = imm;
      }
  Report("max=%3d", (int) imx);

  /* write the file */
  overall_image_num = t.file_index * c.nimages + t.image_num + 1;
  sprintf(fn, "%s/mag.s%.2d.c%d.i%.3d", c.output_directory, t.slice_num+1, t.coil_num+1, overall_image_num);
  Acct(WRITEOPEN);
  if ((f = fopen(fn, "w")) == NULL)
    Abort("Can't open %s for output!\n", fn);
  Acct(WRITING);
  bio_big_endian_output = c.big_endian_output;
  FWrInt16Array(f, &mag[0][0], c.os_res*c.os_res);
  Acct(WRITECLOSE);
  fclose(f);
  Acct(PROCESSING);
  Free2DShortArray(mag);
}

/* UTILITY FUNCTIONS */


void
UncompressFile (Filename out,
		const Filename in)
{
  char* uncompress_cmd= NULL;
  static int count = 0;
  char com[2*sizeof(Filename)+128];

  if (!(uncompress_cmd= getenv("F_UNCOMPRESS")))
    uncompress_cmd= "uncompress";

  Report("Uncompressing raw file %s\n", in);
  sprintf(out, "%s/sgrid.%d.%.3d", c.tmp_directory, getpid(), count++);
  sprintf(com, "%s -c %s > %s", uncompress_cmd, in, out);
  system(com);
#ifdef AFS
  if (CheckForAFS(out)) FlushAFS(out);
#endif
}

void
RemoveFile (const Filename name)
{
  char com[1024];

  sprintf(com, "rm -f %s", name);
  system(com);
}

float
kaiser (float l,
	float b,
	float u)
{
  float f;

  if (fabs(u) <= 0.5*fabs(l))
    {
      f = bessi0(b*sqrt(1.0-(2.0*u/l)*(2.0*u/l)));
      return f;
    }
  else
    return(0.0);
}

float
bessi0 (float x)
{
  float ax,ans;
  double y;
  
  if ((ax=fabs(x)) < 3.75)
    {
      y=x/3.75;
      y*=y;
      ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
    }
  else
    {
      y=3.75/ax;
      ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
			+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
			+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
			+y*0.392377e-2))))))));
    }
  return(ans);
}

int
IntBRdFloat32 (unsigned char *addr)
{
  return((int) Round(BRdFloat32(addr)));
}


void
PackContext ()
{
  par_pkstr(c.input_directory);
  par_pkint(c.big_endian_input);
  par_pkstr(c.tmp_directory);
  par_pkstr(c.output_directory);
  par_pkint(c.big_endian_output);
  par_pkint(c.nfiles);
  par_pkint(c.ncoils);
  par_pkint(c.nslices);
  par_pkint(c.nimages);
  par_pkint(c.nph1);
  par_pkint(c.nphmult);
  par_pkint(c.npr);
  par_pkint(c.ndat);
  par_pkint(c.chop);
  par_pkint(c.risetime);
  par_pkint(c.densamp);
  par_pkdouble(c.ts);
  par_pkdouble(c.gts);
  par_pkdouble(c.fsgcm);
  par_pkdouble(c.opfov);
  par_pkfloat(c.pix_shifth);
  par_pkfloat(c.pix_shiftv);
  par_pkint(c.coil_record_length);
  par_pkint(c.slice_record_length);
  par_pkint(c.baseline_length);
  par_pkint(c.res);
  par_pkint(c.os_res);
  par_pkint(c.slice);
  par_pkint(c.samp_delay);
  par_pkint(c.samp_cor);
  par_pkfloat(c.ph_twist);
  par_pkint(c.lr_shift);
  par_pkint(c.tb_shift);
  par_pkint(c.loc_shift);
  par_pkfloat(c.zoom);
  par_pkfloat(c.mag_factor);
  par_pkfloat(c.ph_factor);
  par_pkint(c.start_slice);
  par_pkint(c.end_slice);
  par_pkint(c.lin_cor);
  par_pkstr(c.reg_file);
  par_pkint(c.reg_2x);
  par_pkint(c.write_mag);
  par_pkint(c.over_samp);
  par_pkfloat(c.grid_len);
  par_pkfloat(c.gridb);
  par_pkfloat(c.wind);
  par_pkfloat(c.factxx);
  par_pkfloat(c.factxy);
  par_pkfloat(c.factyx);
  par_pkfloat(c.factyy);
  par_pkfloatarray(&c.t2k[0][0][0], 2*c.npr*c.ndat);
  par_pkfloatarray(c.kdens, c.ndat);
}

void
UnpackContext ()
{
  static int previously_allocated = FALSE;

  par_upkstr(c.input_directory);
  c.big_endian_input= par_upkint();
  par_upkstr(c.tmp_directory);
  par_upkstr(c.output_directory);
  c.big_endian_output= par_upkint();
  c.nfiles= par_upkint();
  c.ncoils= par_upkint();
  c.nslices= par_upkint();
  c.nimages= par_upkint();
  c.nph1= par_upkint();
  c.nphmult= par_upkint();
  c.npr= par_upkint();
  c.ndat= par_upkint();
  c.chop= par_upkint();
  c.risetime= par_upkint();
  c.densamp= par_upkint();
  c.ts= par_upkdouble();
  c.gts= par_upkdouble();
  c.fsgcm= par_upkdouble();
  c.opfov= par_upkdouble();
  c.pix_shifth= par_upkfloat();
  c.pix_shiftv= par_upkfloat();
  c.coil_record_length= par_upkint();
  c.slice_record_length= par_upkint();
  c.baseline_length= par_upkint();
  c.res= par_upkint();
  c.os_res= par_upkint();
  c.slice= par_upkint();
  c.samp_delay= par_upkint();
  c.samp_cor= par_upkint();
  c.ph_twist= par_upkfloat();
  c.lr_shift= par_upkint();
  c.tb_shift= par_upkint();
  c.loc_shift= par_upkint();
  c.zoom= par_upkfloat();
  c.mag_factor= par_upkfloat();
  c.ph_factor= par_upkfloat();
  c.start_slice= par_upkint();
  c.end_slice= par_upkint();
  c.lin_cor= par_upkint();
  par_upkstr(c.reg_file);
  c.reg_2x= par_upkint();
  c.write_mag= par_upkint();
  c.over_samp= par_upkint();
  c.grid_len= par_upkfloat();
  c.gridb= par_upkfloat();
  c.wind= par_upkfloat();
  c.factxx= par_upkfloat();
  c.factxy= par_upkfloat();
  c.factyx= par_upkfloat();
  c.factyy= par_upkfloat();
  
  /* deallocate old t2k and kdens arrays */
  if (previously_allocated)
    {
      Free3DFloatArray(c.t2k);
      free(c.kdens);
    }

  /* allocate the t2k and kdens arrays */
  c.t2k = Alloc3DFloatArray(2, c.npr, c.ndat);
  c.kdens = (float *) malloc(c.ndat * sizeof(float));
  previously_allocated = TRUE;

  par_upkfloatarray(&c.t2k[0][0][0], 2*c.npr*c.ndat);
  par_upkfloatarray(c.kdens, c.ndat);
}

void
PackTask ()
{
  par_pkint(t.file_index);
  par_pkstr(t.filename);
  par_pkint(t.coil_num);
  par_pkint(t.slice_num);
  par_pkint(t.image_num);
  par_pkfloat(t.reg_xs);
  par_pkfloat(t.reg_ys);
  par_pkfloat(t.reg_rot);
}

void
UnpackTask ()
{
  t.file_index= par_upkint();
  par_upkstr(t.filename);
  t.coil_num= par_upkint();
  t.slice_num= par_upkint();
  t.image_num= par_upkint();
  t.reg_xs= par_upkfloat();
  t.reg_ys= par_upkfloat();
  t.reg_rot= par_upkfloat();
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
