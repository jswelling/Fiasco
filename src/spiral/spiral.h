/*
 *	spiral.h
 *
 *    Header for "spiral" program, q.v.
 *
 *    Copyright by Douglas C. Noll and the University of Pittsburgh and
 *	  the Pittsburgh Supercomputing Center
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
 *
 *    HISTORY
 *	1/98 - split off from spiral.c (Greg Hood, PSC)
 */

#include <math.h>
#include "misc.h"

#define PI	M_PI

#define SWAP(a,b) {tempr=(a);(a)=(b);(b)=tempr;}
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

/* levels of verbosity */
#define VERBOSITY_NONE		0	/* none */
#define VERBOSITY_MINIMAL	1	/* one line per P file */
#define VERBOSITY_FILE		2	/* summary information per P file */
#define VERBOSITY_SLICE		3       /* summary information per slice */
#define VERBOSITY_TIMESTEP	4       /* summary informtion per timestep
					   per slice */

typedef struct Context {
  /* NOTE: any new fields added to this struct should
     be also added to PackContext and UnpackContext */ 
  Filename input_directory;	/* where to find input files */
  Filename reference_directory;	/* where to find reference files */
  int big_endian_input;		/* TRUE if input is big-endian */
  Filename tmp_directory;	/* where to place temporary files */
  Filename output_directory;	/* where to place output files */
  int big_endian_output;	/* TRUE if output if big-endian */
  Filename output_name;		/* name to give to output dataset */
  Filename output_data;		/* file in which to place output data chunk */
  int input_is_mri;             /* TRUE if we are reading Pgh MRI rather
                                   than a Pfile */

  int nfiles;			/* # of input files */
  int ncoils;			/* # of coils */
  int nslices;			/* # of slices */
  int sliceorder;               /* 0=interleaved, 1=sequential */
  int nimages;			/* # of images in one P file (per coil) */
  int nph1;			/* # of baselines per slice */
  int nphmult;			/* # of images (phases) per baseline */
  int npr;			/* # of projections (spirals) */
  int ndat;			/* # of data items per projection */
  int ndatfr;                   /* # of samples in each acq frame */
  int total_images;		/* # of images over all P files */
  int gtype;                    /* k-space trajectory type (future use) */
  int opxres;                   /* nominal output resolution in image space*/
  int ngap;                     /* ??? for concat readout */
  int concat;                   /* ??? for concat readout */

  int chop;			/* may be either 1 or -1 */
  int risetime;			/* ??? */
  int densamp;			/* ??? */

  int tr;  			/* pulse repetition time (microseconds) */
  int te;               	/* pulse echo time (microseconds) */
  short flip;  		        /* flip angle (degrees) */
  float nex;             	/* number of excitations */
  float tlc[3];			/* top left corner, RAS */
  float trc[3];			/* top right corner, RAS */
  float brc[3];			/* bottom right corner, RAS */
  float ctr[3];                 /* center(?), RAS */
  char psd[36];			/* pulse sequence descriptor */
  float slthick;		/* slice thickness (mm) */
  float spacing;		/* slice spacing (mm) */
  char date[12];		/* date of scan */
  char time[10];		/* time of scan */

  double ts;			/* sample spacing in time (seconds) */
  double gts;			/* gradient time spacing (seconds) */ 
  double fsgcm;			/* design gradient max, G/cm */
  double opfov;			/* field-of-view (cm) */

  double slewrate;              /* design slew rate, (mT/m/ms) */
  double fast_rec_lpf;          /* fast rec cutoff in kHz */
  double mapdel;                /* field mapping offset (us) */
  int mapdel_set;               /* non-zero if mapdel is valid */

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

  int rotation;                 /* rotations of coordinate indices */ 
  int transpose;                /* transpose of coordinate indices */

  int start_slice;		/* starting slice number */
  int end_slice;		/* ending slice number */

  int lin_cor;			/* use linear correction map */
  int lin_map;			/* generate linear correction maps */
  int hc_cor;			/* use general (time-segmented) correction */
  int gen_map;			/* generate general (time-segmented) correction map */

  Filename reg_file;		/* registration file */
  int reg_2x;			/* double translation in registration file */

  int write_raw;		/* 1 if we should write raw files out */
  int write_mag;		/* 1 if we should write magnitude files out */
  int write_samples;            /* 1 if we should write sample files out */
  int write_phase;		/* TRUE if we should write image phase file */
  int do_recon;                 /* TRUE if we should actually reconstruct */

  int over_samp;		/* oversampling ratio */
  float grid_len;		/* grid length (half width of convolution) */
  float gridb;
  float wind;			/* window */

  float factxx;			/* image rotation & scaling factors */
  float factxy;
  float factyx;
  float factyy;

  float samp_time;              /* sample time (usec) */
  int filter_sz;		/* total width of filter in pixels */

  int all_coils;		/* TRUE if we should write a separate image
				   for each coil */
  int output_float;		/* TRUE if we should write the output dataset
				   in float32 format */
  float output_scale;		/* scaling factor if writing in int16 format */
  int verbosity;		/* verbosity level */
  int reverse;                  /* data aquired with reverse spiral path */

/* variably-sized context info */
  float ***t2k;			/* dimensioned as [2][npr][ndat] */
  float *kdens;			/* dimensioned as [ndat] */
  float ***refim_in;		/* [nslices][res][res]  input reference */
  float ***sampim;		/* [nslices][os_res][os_res]  sample info */
  float **refl;		        /* [nslices][3] linear reference data */
  unsigned char *ref_missing;    /* [nslices] reference missing data */
} Context;

typedef struct Task {
  /* NOTE: any new fields added to this struct should
     be also added to PackTask and UnpackTask */
  int file_index;	/* counts which input file we are
			   currently working on (0 => first) */
  Filename filename;	/* input filename */

  int slice_num;	/* slice number (first is 0) */
  int image_num;	/* image number (over time) (first is 0)
			   (corresponds to phase number in gsp14.c) */
  int overall_image_num; /* image number within the entire set of
			    P files */

  float reg_xs;		/* the registration x-shift for this image */
  float reg_ys;		/* the registration y-shift for this image */
  float reg_rot;	/* the registration rotation for this image */
} Task;

typedef struct Result {
  /* NOTE: any new fields added to this struct should
     be also added to PackResult and UnpackResult */
  int image_num;        /* image number */
  int slice_num;	/* slice number */
  int sampim_valid;     /* sampim is valid (used for reference images) */
  int recon_failed;     /* non-zero if recon fails for any reason */
  float **slice_sampim;	/* [os_res][os_res] sample info for the slice */
} Result;

typedef struct FilenameList {
  Filename fname;
  struct FilenameList* next;
} FilenameList;

/* GLOBAL VARIABLES FOR MASTER & WORKERS */
extern Task t;
extern Context c;
extern Result r;

/* EXPORTED FUNCTIONS */
void MasterTask (const int argc, const char **argv, const char **envp);
void MasterResult (int task_number);
void WorkerContext ();
void WorkerTask ();
void WorkerFinalize();
void ReadFileHeader (const Filename input_file, Context *c, 
		     unsigned char *hdr);
void CreateMRIFile (const Filename output_file, int argc, char** argv, 
		    Context *c);
void ReadMRIFileHeader (const Filename input_file, Context *c);
void ReadSliceCoords(Context* c, unsigned char* hdr, int islice,
		     float p1[3], float p2[3], float p3[3], 
		     float* s_i_offset_out);
void CheckHeaderInfo (const Filename input_file, Context *c);
void CheckMRIHeaderInfo (const Filename input_file, Context *c);
void WorkerInitialize();
void WorkerShutdown();
void GenerateReference ();
void GenerateImage ();

void PackContext();
void UnpackContext();
void PackTask();
void UnpackTask();
void PackResult();
void UnpackResult();
