/************************************************************
 *                                                          *
 *  reader.c                                                *
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
 *  Original programming by Mark Fitzgerald  5-96           *
 *  Modified to exclusively use libmri calls for output,    *
 *       Greg Hood (PSC), 9-98                              *
 *  Modified to read header files, and to use LX2           *
 *       resampling stuff, Joel Welling (PSC/Stats), 5-1999 *
 ************************************************************/
/* This routine uses functionality from the "epirecon" code
 * supplied with the GE LX2 scanner system, author 
 * Bryan J. Mock (GE Medical Systems).
 */

/* Notes-
   -Is their homodyne recon really equiv to our partialk?
   -Header provides info necessary to reorient sagital and coronal images,
    but not oblique.
   -What's the blank lines stuff in epirecon?
   -Am I supporting partial FOV scanning?  Warn user either way.
   -I *think* twoshot is in the input with successive slices for
    the same shot consecutively (then moving on to the next shot).
   -What's "dc zipper" effect that baseline data is supposed to get
    rid of for fast receiver case?
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include "mri.h"
#if (SGI64 || SGI5 || SGIMP)
#include <bstring.h>
#endif

#include "bio.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"
#include "nr_sub.h"
#include "rcn.h"

/* epirecon header files */
#include "frozen_header_info.h"

/* Windaq header file */
#include "windaq_header_info.h"

static char rcsid[] = "$Id: reader.c,v 1.31 2007/03/22 00:03:47 welling Exp $";

#define BLANKY_NUM_DEFAULT 2
#define BLANKY_NUM_MAX 10
#define FILE_DEFAULT 0
#define FILE_LX 1
#define FILE_WINDAQ 2


/* This is the range of the FFT of the first image with autoscaling on */
#define AUTOSCALE_RANGE 8192.0

int debug = 0;          /* Global debug flag                         */

static long typesize[6] = { 0, 1, 2, 4, 4, 8 };

static char typeconv[6][8] = { "raw", "uint8", "int16", "int32", "float32", 
			       "float64" };

static char* progname; /* program name */
static int do_phaseref= 0;   /* do phase reference correction */
static float *ref_lp= NULL;      /* Phase correction vectors (linear and */
static float *ref_poff= NULL;    /* const) from ref.dat file */
static int do_bandpass= 0;   /* do bandpass asymmetry correction */
static float bandwidth;            /* recevier bandwidth from header */
static int do_rampsample= 0; /* do ramp sampling correction */
static float **ramp_vrgf_coeff= NULL;    /* vrgf coeffient array */
static float band_aphase[2048];    /* band pass asymmetry variables - phase */
static float band_amag[2048];      /*  & mag */
static int band_nbpc;              /* num points */
static int do_xchop= 0;            /* chop in x dir (phase shift pi/2) */
static int do_ychop= 0;            /* chop in y dir (phase shift pi/2) */
static int do_autoscale= 0;        /* rescale to "optimal" range */
static int autoscale_set= 0;       /* scale has been calculated */
static float autoscale;            /* value to use for autoscaling */
static int filetype=0;             /* regular, LX or Windaq file */
static char plane[4]= "";          /* if pl = L/R > axial or A/P > sag, */
                                   /* S/I > coronal */

typedef struct scan_info_struct {
  float voxel_x, voxel_y; /* in mm */
  float image_x, image_y; /* in voxels */
  float fov_x, fov_y;     /* in mm */
  float slice_thickness, slice_gap; /* in mm */
  float TR, TE; /* in msec */
  int overscan; /* number of scan lines past Ky=0 */
  char* date;
  char* time;
} ScanInfo;

static int scan_LX_data_header(char* readfile, MRI_Datatype *datatype, 
			    long *veclen, long *dx_raw,
			    long *dx, long *dy, long *ds, 
			    long *dz, long *dt, long *vxytz_flag, 
			    long *reorder_flag, long *start_offset,
			    long *skip, long *sliceskip,
			    ScanInfo* scan_info)
{
  FILE *fphead;
  unsigned char header[FRZ_RDB_HEADER_SIZE_BYTES];
  unsigned char* rdbhead;
  unsigned char* acq_tab;
  unsigned char* examhead;
  unsigned char* serieshead;
  unsigned char* imagehead;
  
  /* Read from header.  A two-step translation then happens:
   * from the header fields to the language of Mock's "epirecon",
   * then to the language of this program.
   */
  short exnum;           /* Exam Number */
  int flip;              /* flip angle */
  int xres, yres;        /* acq. resolution */
  int rcxres,rcyres;     /* reconstructed resolution (with fovar=1.0) */
  float fovar;            /* Field of View Aspect Ratio from header */
  int numreps;            /* total reps and rep number (nex in header) */
  int numslice;           /* number of slices */
  int frame;             /* size of each image points calc. from */
                         /* header info */
  int ssp_size;          /* Adding these together (+ header) should */
  int pt_size = 0;        /* For extended Dynamic Range data */
  int rhraw_size;        /* give expected file size in bytes        */
  float hnw;             /* homodyne correction transition width in */
                         /* header */
  float fw, fr;           /* fermi width and radius from header */
  char pl;                /* orientation of slice */
  short swapfp;           /* swapped freq/phase encode direction 1 = */
                          /* yes - determines if images are reoriented */
  short rot, tpose;
  int blanky_num;         /* # of y lines to blank @ top & bottom of image */
  int bviews = 0;         /* # of baseline views */
  int ileaves = 1;        /* # of interleaves if multi-shot EPI */
  int fast_rec = 0;       /* fast_reciever check */
  int vpsht;              /* Ky views per shot */
  int bl_save= 0;         /* 1 for Genesis, 0 for LX */
  int interm_xres;        /* reconstructed x resolution (with fovar != 1.0) */
  float revnum;           /* header revision number */
  int ierror= 0;

  /* This bit is from the GE code io_signa_lx.c, with mods for portability */
  ierror= 0;
  if ((fphead = fopen(readfile,"r"))!=NULL)
    {
      if (fseek(fphead, (long) 0, SEEK_SET)) {
	perror("Error seeking header");
	ierror=1;
      }
      else {
	if (fread(header,sizeof(char),FRZ_RDB_HEADER_SIZE_BYTES, fphead)
	    != FRZ_RDB_HEADER_SIZE_BYTES) {
	  perror("Error reading header");
	  ierror=1;
	}
	else {
	  if (fclose(fphead)) {
	    perror("Error closing header");
	    ierror=1;
	  }
	}
      }
    }
  else {
    perror("Error opening header");
    ierror= 1;
  }
  if (ierror) return 0;

  rdbhead= header+FRZ_RDB_HDR_OFF;
  acq_tab= header+FRZ_RDB_DATAACQ_OFF;
  examhead= header+FRZ_RDB_EXAMDATATYPE_OFF;
  serieshead= header+FRZ_RDB_SERIESDATATYPE_OFF;
  imagehead= header+FRZ_RDB_MRIMAGEDATATYPE_OFF;

  /* We test to see if the data does have the expected endian order-
   revision number should be in the neighborhood of 7.0 */
  revnum= BRdFloat32(rdbhead+FRZ_RDBHEAD_RDB_HDR_RDBM_REV_OFF);
  if (revnum<1.0 || revnum>1000.0) {
    /* Oops, try it the other way! */
    bio_big_endian_input= (bio_big_endian_input ? 0 : 1);
    revnum= BRdFloat32(rdbhead+FRZ_RDBHEAD_RDB_HDR_RDBM_REV_OFF);
  }
  if (debug) fprintf(stderr,"Header indicates %s input\n",
		     bio_big_endian_input ? "bigendian" : "littleendian");

  /* Set acquired matrix size xres,yres, reconstruction size rcxres, */
  /* rcyres, number of slices, frame size of each image, and number of */
  /* reps (nex) from rawheader */
  exnum = BRdInt16(examhead+FRZ_EXAMHEAD_EX_NO_OFF);
  flip = BRdInt32(imagehead+FRZ_IMAGEHEAD_MR_FLIP_OFF);
  xres = BRdInt16(rdbhead+FRZ_RDBHEAD_RDB_HDR_DA_XRES_OFF);
  yres = BRdInt16(rdbhead+FRZ_RDBHEAD_RDB_HDR_DA_YRES_OFF)-1;
  scan_info->slice_thickness =  
    BRdFloat32(imagehead+FRZ_IMAGEHEAD_SLTHICK_OFF);
  scan_info->slice_gap = 
    BRdFloat32(imagehead+FRZ_IMAGEHEAD_SCANSPACING_OFF);
  scan_info->TR = BRdInt32(imagehead+FRZ_IMAGEHEAD_TR_OFF)/1000.;
  scan_info->TE = BRdInt32(imagehead+FRZ_IMAGEHEAD_TE_OFF)/1000.;
  scan_info->fov_x =  BRdFloat32(imagehead+FRZ_IMAGEHEAD_DFOV_OFF);
  scan_info->fov_y =  
    BRdFloat32(imagehead+FRZ_IMAGEHEAD_DFOV_RECT_OFF);
  scan_info->voxel_x= 
    BRdFloat32(imagehead+FRZ_IMAGEHEAD_PIXSIZE_X_OFF);
  scan_info->voxel_y= 
    BRdFloat32(imagehead+FRZ_IMAGEHEAD_PIXSIZE_Y_OFF);
  scan_info->image_x= BRdFloat32(imagehead+FRZ_IMAGEHEAD_DIM_X_OFF);
  scan_info->image_y= BRdFloat32(imagehead+FRZ_IMAGEHEAD_DIM_Y_OFF);
  scan_info->overscan= 
    BRdInt16(rdbhead+FRZ_RDBHEAD_RDB_HDR_HNOVER_OFF);
  scan_info->date= 
    strdup((char*)(rdbhead+FRZ_RDBHEAD_RDB_HDR_SCAN_DATE_OFF));
  scan_info->time= 
    strdup((char*)(rdbhead+FRZ_RDBHEAD_RDB_HDR_SCAN_TIME_OFF));

  fovar = BRdFloat32(rdbhead+FRZ_RDBHEAD_RDB_HDR_PHASE_SCALE_OFF);
  if(fovar == 0.5) scan_info->fov_y /= fovar;
  
  rcxres = BRdInt16(rdbhead+FRZ_RDBHEAD_RDB_HDR_RC_XRES_OFF);
  interm_xres = (int)((float)rcxres/fovar);
  rcyres = BRdInt16(rdbhead+FRZ_RDBHEAD_RDB_HDR_RC_YRES_OFF);
  numreps = BRdFloat32(imagehead+FRZ_IMAGEHEAD_NEX_OFF);
  numslice = BRdInt16(rdbhead+FRZ_RDBHEAD_RDB_HDR_NSLICES_OFF);
  /* allow for ext dyn. range data */
  pt_size = BRdInt16(rdbhead+FRZ_RDBHEAD_RDB_HDR_POINT_SIZE_OFF);
  frame = 
    xres*yres*BRdInt16(rdbhead+FRZ_RDBHEAD_RDB_HDR_POINT_SIZE_OFF);
  ssp_size = BRdInt32(rdbhead+FRZ_RDBHEAD_RDB_HDR_SSPSAVE_OFF);
  rhraw_size = 
    BRdInt32(rdbhead+FRZ_RDBHEAD_RDB_HDR_RAW_PASS_SIZE_OFF);
  hnw = BRdFloat32(rdbhead+FRZ_RDBHEAD_RDB_HDR_NTRAN_OFF);
  fw = BRdFloat32(rdbhead+FRZ_RDBHEAD_RDB_HDR_FERMI_WIDTH_OFF);
  fr = BRdFloat32(rdbhead+FRZ_RDBHEAD_RDB_HDR_FERMI_RADIUS_OFF);
  pl = *((char*)imagehead + FRZ_IMAGEHEAD_LOC_RAS_OFF);
  swapfp = BRdInt16(imagehead+FRZ_IMAGEHEAD_SWAPPF_OFF);
  bandwidth = BRdFloat32(imagehead+FRZ_IMAGEHEAD_VBW_OFF);
  rot = BRdInt16(rdbhead+FRZ_RDBHEAD_RDB_HDR_ROTATION_OFF);
  tpose = BRdInt16(rdbhead+FRZ_RDBHEAD_RDB_HDR_TRANSPOSE_OFF);
  blanky_num = BRdInt16(rdbhead+FRZ_RDBHEAD_RDB_HDR_SLBLANK_OFF);
  bviews = BRdInt16(rdbhead+FRZ_RDBHEAD_RDB_HDR_BASELINE_VIEWS_OFF);
  fast_rec = BRdInt32(rdbhead+FRZ_RDBHEAD_RDB_HDR_FAST_REC_OFF);
  ileaves = BRdInt16(rdbhead+FRZ_RDBHEAD_RDB_HDR_ILEAVES_OFF);

  /* set number of blank lines to 2 if out of bounds */
  if( (blanky_num == 0) || (blanky_num > BLANKY_NUM_MAX) )
    blanky_num = BLANKY_NUM_DEFAULT;

  /* number of views per shot */
  if(ileaves != 0) vpsht= yres/ileaves;
  else vpsht= yres; /* presumably single shot data */

  /* BJM: safety check in case BW field in image header is NULL */
  if(bandwidth == 62.0)
    bandwidth+= 0.5;
  else if (bandwidth == 0.0) /* assume GE epibold.e */
    bandwidth = BRdFloat32(rdbhead+FRZ_RDBHEAD_RDB_HDR_USER16_OFF);
  else if (bandwidth == 0.0) { /* if STILL Zero.. */
    Warning(1,"%s: Warning: Reciever bandwidth was ZERO\n", progname);
    if (do_bandpass) {
      Warning(1,
	      "%s: Not performing band pass asymmetry correction\n",progname);
    }
    do_bandpass= 0;
  }

  /* Determine Image Orientation - use imagehead.plane integer istead */
  /* of loc raster character */
  if(BRdInt16(imagehead+FRZ_IMAGEHEAD_PLANE_OFF) == 2)
    sprintf(plane,"AX");
  else if(BRdInt32(imagehead+FRZ_IMAGEHEAD_PLANE_OFF) == 4)
    sprintf(plane,"SAG");
  else if(BRdInt32(imagehead+FRZ_IMAGEHEAD_PLANE_OFF) == 8)
    sprintf(plane,"COR");
  else {
    Warning(1,"%s: Unable to detect image orientation\n",progname);
    Warning(1,"%s: Assuming AXIAL images; verify orientation!\n",progname);
    sprintf(plane,"AX");
  }

  if (pt_size<4) *datatype= MRI_SHORT;
  else *datatype= MRI_INT;
  *veclen= 2;
  *dx_raw= xres;
  if (do_rampsample) *dx= rcxres;
  else *dx= *dx_raw;
  *dy= vpsht;
  if(ileaves != 0) *ds= ileaves;
  else *ds= 1;
  *dz= numslice;
  *dt= numreps;
  *vxytz_flag= 0;
  *reorder_flag= 1;
  *start_offset= FRZ_POOL_HEADER_SIZE + (2*pt_size*(bviews+bl_save));
  *skip= 0;
  *sliceskip= 0;

   return 1;
}


static int scan_WINDAQ_data_header(char* readfile, MRI_Datatype *datatype, 
			    long *veclen, long *dx_raw,
			    long *dx, long *dy, long *ds, 
			    long *dz, long *dt, long *vxytz_flag, 
			    long *reorder_flag, long *start_offset,
			    long *skip, long *sliceskip)
{
  unsigned char header[WINDAQ_HEADER_SIZE_BYTES];
  FILE *fphead;
  int Windaq_Logo1;
  unsigned long data_bytes;
  int ierror= 0;
  int i, channel[29];
  float slope1[29], intercept1[29];
  double slope2[29], intercept2[29];

  if ((fphead = fopen(readfile,"r"))!=NULL)
    {
      if (fseek(fphead, (long) 0, SEEK_SET)) {
	perror("Error seeking header");
	ierror=1;
      }
      else {
	if (fread(header, sizeof(char), WINDAQ_HEADER_SIZE_BYTES, fphead)
	    != WINDAQ_HEADER_SIZE_BYTES) {
	  perror("Error reading header");
	  ierror=1;
	}
	else {
	  if (fclose(fphead)) {
	    perror("Error closing header");
	    ierror=1;
	  }
	}
      }
    }
  else {
    perror("Error opening header");
    ierror= 1;
  }
  
  if (ierror) return 0;
  
  bio_big_endian_input = (((Windaq_Logo1 = BRdInt8(header+WINDAQ_LOGO_OFF1)) == 1) ? 0 : 1);

  *veclen = (BRdUInt8(header+WINDAQ_CHANNELS_OFF) & 0x1f);
  *dx=*dy=*ds=*dz=*dx_raw=1;

  if (debug) fprintf(stderr, "Header indicates %s input \n", 
	  bio_big_endian_input ? "bigendian" : "littleendian");

  data_bytes=BRdInt32(header+WINDAQ_DATAACQ_SIZE_OFF);
  *dt=((data_bytes)/2)/(*veclen);

  /* to get slope and intercepts for each channel in case we need to 
   * scale the raw data to match Windaq viewing program 
   */

  if (debug) {

    for (i=0; i<*veclen; i++){
      slope1[i]=BRdFloat32(header+(110+36*i));
      intercept1[i]=BRdFloat32(header+(114+36*i));
      slope2[i]=BRdFloat64(header+(118+36*i));
      intercept2[i]=BRdFloat64(header+(126+36*i));
      fprintf(stderr, "slope1: %.2f  int1: %.2f  slope2: %.2f  int2: %.2f \n",
	      slope1[i], intercept1[i], slope2[i], intercept2[i]);
      channel[i]=BRdInt8(header+(142+36*i));
      if (channel[i]>64) fprintf(stderr, "differential channel %d \n", channel[i]-64);
      else fprintf(stderr, "single ended channel %d \n", channel[i]);
    }
  }

  *datatype=MRI_SHORT;

  *vxytz_flag=0;
  *reorder_flag=0;
  *skip=0;
  *sliceskip=0;
  *start_offset=WINDAQ_DATAACQ_OFF;
  

  return 1;
}


static int process_default(char* creorder, char* cdataorder, char* cdatatype, 
			   long *reorder_flag, long *vxytz_flag, MRI_Datatype
			   *datatype_in, MRI_Datatype *datatype_out, long *ds,
			   long *dx_raw, long dx)
{
  /* Convert reorder switch to flag */
  if ( isdigit( creorder[0] ) )
    *reorder_flag = atol( creorder );
  else
    *reorder_flag =
      (( ( creorder[0] == 'T' ) || ( creorder[0] == 't' ) ) ? 1: 0);

  /* Convert data order switch to flag */
  if ( !strcmp(cdataorder, "vxyzt") ) *vxytz_flag= 0;
  else if ( !strcmp(cdataorder, "vxytz") ) *vxytz_flag= 1;
  else Abort("Dataorder argument %s invalid; must be vxyzt or vxytz.\n",
	     cdataorder);

  /* Convert data-type string to number */
  if( !strcmp( cdatatype, "MRI_SHORT" ) ||
      !strcasecmp( cdatatype, "short" ) ||
      !strcasecmp( cdatatype, "s" ) )
    *datatype_in = MRI_SHORT;
  else if( !strcmp( cdatatype, "MRI_FLOAT" ) ||
	   !strcasecmp( cdatatype, "float" ) ||
	   !strcasecmp( cdatatype, "f" ) )
    *datatype_in = MRI_FLOAT;
  else if( !strcmp( cdatatype, "MRI_UNSIGNED_CHAR" ) ||
	   !strcasecmp( cdatatype, "uchar" ) ||
	   !strcasecmp( cdatatype, "u" ) ||
	   !strcasecmp( cdatatype, "c" ) )
    *datatype_in = MRI_UNSIGNED_CHAR;
  else if( !strcmp( cdatatype, "MRI_INT" ) ||
	   !strcasecmp( cdatatype, "long" ) ||
	   !strcasecmp( cdatatype, "l" ) ||
	   !strcasecmp( cdatatype, "int" ) ||
	   !strcasecmp( cdatatype, "i" ) )
    *datatype_in = MRI_INT;
  else if( !strcmp( cdatatype, "MRI_DOUBLE" ) ||
	   !strcasecmp( cdatatype, "double" ) ||
	   !strcasecmp( cdatatype, "d" ) )
    *datatype_in = MRI_DOUBLE;
  else 
    Abort( "Data type unrecognized: %s.", cdatatype );

  if (do_xchop || do_ychop || do_autoscale) {
    *datatype_out= MRI_FLOAT;
    Warning(1,
	   "%s: Warning: changing data to float in order to process flag\n",
	    progname);}
  else *datatype_out= *datatype_in;
  *ds=1; /* count of slice-sequential shots */
  *dx_raw= dx;

  return 0;
}


static int process_LX( char* readfile, MRI_Datatype *datatype_in, long *veclen,
		       long *dx_raw, long *dx, long *dy, long *ds, long *dz, long *dt,
		       long *vxytz_flag, long *reorder_flag, long *start_offset,
		       long *skip, long *sliceskip,
		       ScanInfo* scan_info, MRI_Datatype *datatype_out)
{ 

  if (!(scan_info= (ScanInfo*)malloc(sizeof(ScanInfo))))
    Abort("%s: unable to allocate %d bytes!\n",sizeof(ScanInfo));

  if (!scan_LX_data_header(readfile, datatype_in, 
       	      veclen, dx_raw, dx, dy, ds, dz, dt, 
       	      vxytz_flag, reorder_flag, start_offset,
       	      skip, sliceskip, scan_info))
    Abort("reader: error reading data header file.\n");

  if (scan_info->overscan==0)
    do_ychop= 1; /* Only imaged one side of Ky=0 */
  *datatype_out= MRI_FLOAT; /* since at least filetype = FILE_LX  is set */

  return 0;
}


static int process_WINDAQ( char* readfile, MRI_Datatype *datatype_in, long *veclen,
			   long *dx_raw, long *dx, long *dy, long *ds, long *dz, 
			   long *dt, long *vxytz_flag, long *reorder_flag, 
			   long *start_offset, long *skip, long *sliceskip,
			   MRI_Datatype *datatype_out)
{
  if (!scan_WINDAQ_data_header(readfile, datatype_in, 
       	      veclen, dx_raw, dx, dy, ds, dz, dt, 
       	      vxytz_flag, reorder_flag, start_offset,
       	      skip, sliceskip))
    Abort("reader: error reading data header file.\n");

  if (do_phaseref || do_bandpass || do_rampsample || do_xchop || do_ychop
      || do_autoscale)
    Warning(1, 
	"Windaq file does not allow certain flags to have values. Resetting to zero.\n");

  do_phaseref=0;
  do_bandpass=0;
  do_rampsample=0;
  do_xchop=0;
  do_ychop=0;
  do_autoscale=0;

  fprintf(stderr, "number of channels = %ld \n", *veclen);
  fprintf(stderr, "dx=%ld, dy=%ld, ds=%ld, dz=%ld, dx_raw=%d \n", *dx, *dy, *ds, *dz, *dx_raw);
  fprintf(stderr, "dt=%ld \n", *dt);
  fprintf(stderr, "vxytz_flag=%ld, reorder_flag=%ld \n", *vxytz_flag, *reorder_flag);
  fprintf(stderr, "start_offset=%ld \n", *start_offset);
  fprintf(stderr, "skip=%ld, sliceskip=%ld \n", *skip, *sliceskip);

  *datatype_out = MRI_SHORT;

  return 0;
}


static MRI_Dataset* open_output( char* hdrfile, char* datafile,
				 MRI_Datatype datatype, int vxytz_flag,
				 long veclen, long dx, long dy, long ds,
				 long dz, long dt, ScanInfo* scan_info,
				 char* id_tag, int argc, char** argv)
{
  MRI_Dataset* result;

  /* Initialize header */
  result = mri_open_dataset( hdrfile, MRI_WRITE );
  hist_add_cl(result,argc,argv);
  
  /* Establish header parameters */
  mri_create_chunk( result, "images" );
  mri_set_string( result, "images.file", datafile );
  mri_set_string( result, "images.datatype", typeconv[datatype] );
  if (ds==1) {
    if (vxytz_flag) mri_set_string( result, "images.dimensions", "vxytz" );
    else mri_set_string( result, "images.dimensions", "vxyzt" );
  }
  else {
    if (vxytz_flag) mri_set_string( result, "images.dimensions", "vxystz" );
    else mri_set_string( result, "images.dimensions", "vxyzst" );
    mri_set_int( result, "images.extent.s", ds );
  }
  mri_set_int( result, "images.extent.v", veclen );
  mri_set_int( result, "images.extent.x", dx );
  mri_set_int( result, "images.extent.y", dy );
  mri_set_int( result, "images.extent.z", dz );
  mri_set_int( result, "images.extent.t", dt );
  if (scan_info) {
    float fov_z;
    mri_set_float( result, "images.voxel_size.x",scan_info->voxel_x);
    mri_set_float( result, "images.voxel_size.y",scan_info->voxel_y);
    mri_set_float( result, "images.voxel_size.z",scan_info->slice_thickness);
    mri_set_float( result, "images.voxel_spacing.x",scan_info->voxel_x);
    mri_set_float( result, "images.voxel_spacing.y",scan_info->voxel_y);
    mri_set_float( result, "images.voxel_spacing.z",
		   scan_info->slice_thickness + scan_info->slice_gap);
    mri_set_float( result, "images.fov.x", scan_info->fov_x);
    mri_set_float( result, "images.fov.y", scan_info->fov_y);
    fov_z= (dz*scan_info->slice_thickness)
      + ((dz-1)*scan_info->slice_gap);
    mri_set_float( result, "images.fov.z", fov_z);
    mri_set_float( result, "images.tr", scan_info->TR);
    mri_set_float( result, "images.te", scan_info->TE);
    if (((int)scan_info->image_x) != dx)
      mri_set_float( result, "images.image.x", scan_info->image_x);
    if (((int)scan_info->image_y) != dy)
      mri_set_float( result, "images.image.y", scan_info->image_y);
    if (scan_info->date && *(scan_info->date))
      mri_set_string( result, "images.scan.date", scan_info->date);
    if (scan_info->time && *(scan_info->time))
      mri_set_string( result, "images.scan.time", scan_info->time);
  }
  if (id_tag[0])
    mri_set_string( result, "images.scan.id", id_tag );

  return result;
}


void float_convert(void* temp_out, void* temp_in, long dx_raw, 
		   long dy, long veclen, MRI_Datatype datatype_in)
{
  long length= dx_raw*dy*veclen;
  int i;
  float* fout= (float*)temp_out;

  switch (datatype_in) {
    case MRI_UNSIGNED_CHAR:
      for (i=0; i<length; i++) fout[i]= (float)(*((unsigned char*)temp_in+i));
      break;
    case MRI_SHORT:
      for (i=0; i<length; i++) fout[i]= (float)(*((short*)temp_in+i));
      break;
    case MRI_INT:
      for (i=0; i<length; i++) fout[i]= (float)(*((int*)temp_in+i));
      break;
    case MRI_FLOAT:
      bcopy(temp_in, temp_out, dx_raw*dy*veclen*typesize[datatype_in]);
      break;
    case MRI_DOUBLE:
      for (i=0; i<length; i++) fout[i]= (float)(*((double*)temp_in+i));
      break;
    default:
      Abort("%s: Internal error: unknown datatype %d\n",progname,
	    datatype_in);
  }
}


void baseline_correct( FComplex* image, int dx, int dy )
{
  double res_real= 0.0;
  double res_imag= 0.0;
  long count= 0;
  int i, j;

  for (j=0; j<dy; j++) {
    if (j==dy/4) j += dy/2; /* skip central half of y */
    for (i=0; i<dx; i++) {
      if (i==dx/4) i += dx/2; /* skip central half of x */
      res_real += (double)image[j*dx+i].real;
      res_imag += (double)image[j*dx+i].imag;
      count++;
    }
  }
  res_real /= count;
  res_imag /= count;

  for (j=0; j<dy; j++)
    for (i=0; i<dx; i++) {
      image[j*dx+i].real -= res_real;
      image[j*dx+i].imag -= res_imag;
    }

}


static float calc_scale(FComplex* image, long dx, long dy, float target)
{
  float result;
  float max;
  FComplex* scratch;
  int i;

  if (!(scratch= (FComplex*)malloc(dx*dy*sizeof(FComplex))))
    Abort("%s: unable to allocate %d bytes!\n",dx*dy*sizeof(FComplex));

  for (i=0; i<dx*dy; i++) {
    scratch[i].real= image[i].real;
    scratch[i].imag= image[i].imag;
  }

  fft3d(scratch,dx,dy,1,-1,"xy");

  max= 0.0;
  for (i=0; i<dx*dy; i++) {
    float mod= Modulus(scratch[i]);
    if (mod>max) max= mod;
  }

  result= target/max;

  free(scratch);
  return result;
}


static void apply_scale(FComplex* image, long dx, long dy, float scale)
{
  int i;

  for (i=0; i<dx*dy; i++) {
    image[i].real *= scale;
    image[i].imag *= scale;
  }
}


static void filter_default(void* temp_filtered, void* temp_out, void* temp_in,
			   MRI_Datatype datatype_out, MRI_Datatype datatype_in,
			   long dx_raw, long dx, long dy, long veclen, long z,
			   long t)
{
  if ((datatype_in==datatype_out)
      && !do_xchop && !do_ychop && !do_autoscale)
    bcopy(temp_in, temp_filtered, dx*dy*veclen*typesize[datatype_in]);
  else if ((do_xchop || do_ychop || do_autoscale) && datatype_out==MRI_FLOAT){
    float_convert(temp_out, temp_in, dx_raw, dy, veclen,
		  datatype_in);

    /*
     * Filtering steps follow:
     */

    /* Flip raw data, since filtering routines expect it */
    fliprawy(temp_out,dx_raw,dy);

    /* In this case, filtering is a no-op. */
    bcopy(temp_out, temp_filtered, dx*dy*veclen*sizeof(float));

    /* Baseline correction must be done at this point, or chop 
     * operations will shift baseline error pixel from origin to
     * edge of I-space image.
     */
    baseline_correct(temp_filtered, dx, dy);

    if (do_xchop) choprawx(temp_filtered, dx, dy);
    if (do_ychop) choprawy(temp_filtered, dx, dy);

    /* Implementing routines for autoscaling expect row 
     * flipping to be done, we do it and undo it again after.
     */
    if (do_autoscale) {
      fliprowx(temp_filtered,dx,dy,0,1,NULL,0);

      if (!autoscale_set) {
	  autoscale= calc_scale((FComplex*)temp_filtered,dx,dy,
				AUTOSCALE_RANGE);
	  autoscale_set= 1;
      }
      apply_scale((FComplex*)temp_filtered,dx,dy,autoscale);

      fliprowx(temp_filtered,dx,dy,0,1,NULL,0);
    }


    /* Undo initial y flip */
    fliprawy(temp_filtered,dx,dy);
  }

  else Abort("%s: Output datatype %s not compatible with filtering!\n",
	     progname, typeconv[datatype_out]);  
}


static void filter_LX(void* temp_filtered, void* temp_out, void* temp_in,
			   MRI_Datatype datatype_out, MRI_Datatype datatype_in,
			   long dx_raw, long dx, long dy, long veclen, long z, long t)
{
  if (datatype_out==MRI_FLOAT){
    float_convert(temp_out, temp_in, dx_raw, dy, veclen,
		  datatype_in);

    /*
     * Filtering steps follow:
     */

    /* Flip raw data, since filtering routines expect it */
    fliprawy(temp_out,dx_raw,dy);

    /* Ramp sampling */
    if (do_rampsample) {
      vrgf_filt(temp_filtered,temp_out,ramp_vrgf_coeff,dx,dx_raw,dy);
    }
    else if (dx==dx_raw) {
      bcopy(temp_out,temp_filtered, dx*dy*veclen*sizeof(float));
    }
    else Abort("%s: ramp filtering is implied but missing!\n",progname);

    /* Baseline correction must be done at this point, or chop 
     * operations will shift baseline error pixel from origin to
     * edge of I-space image.
     */
    baseline_correct(temp_filtered, dx, dy);

    if (do_xchop) choprawx(temp_filtered, dx, dy);
    if (do_ychop) choprawy(temp_filtered, dx, dy);

    /* Implementing routines for bandpass, phase reference, and
     * autoscaling expect row flipping to be done, we do it and 
     * undo it again after.
     */
    if (do_bandpass || do_phaseref || do_autoscale) {
      fliprowx(temp_filtered,dx,dy,0,1,NULL,0);

      /* Bandpass asymmetry */
      if (do_bandpass) {
	bp_corr(temp_filtered, dx, dy, bandwidth, band_amag, band_aphase,
		band_nbpc, (int)t, (int)z);
      }

      /* phase reference data correction. */
      if (do_phaseref) {
	float scale= 1.0;
	temp_filtered=
	  epi_correction(temp_filtered,dx,dx,dy,dy,scale,ref_lp,
			 ref_poff);
      }

      /* autoscaling */
      if (do_autoscale) {
	if (!autoscale_set) {
	  autoscale= calc_scale((FComplex*)temp_filtered,dx,dy,
				AUTOSCALE_RANGE);
	  autoscale_set= 1;
	}
	apply_scale((FComplex*)temp_filtered,dx,dy,autoscale);
      }

      fliprowx(temp_filtered,dx,dy,0,1,NULL,0);
    }

    if (!do_phaseref) {
      /* if no ref.dat, shift by pixel in phase-encode direction to
       * match on-line recon for GE LX systems
       */
      pixshifty(temp_filtered, dx, dy, -1);
    }

    /* Undo initial y flip */
    fliprawy(temp_filtered,dx,dy);
  }
  else Abort("%s: Output datatype %s not compatible with filtering!\n",
	     progname, typeconv[datatype_out]);
}


void filter_WINDAQ(void* temp_filtered, void* temp_out, void* temp_in,
		   MRI_Datatype datatype_out, MRI_Datatype datatype_in,
		   long dx, long dy, long veclen)
{
  int i;
  short* in= (short*)temp_in;
  short* filtered= (short*)temp_filtered;

  for (i=0; i<veclen; i++) {
    if (in[i]>>15 == 1) filtered[i] = 0;
    else filtered[i] = in[i]>>4;
    
  }

}



static void move_and_filter( void* temp_filtered, void* temp_out, 
			     void* temp_in, 
			     MRI_Datatype datatype_out, 
			     MRI_Datatype datatype_in,
			     long dx_raw, long dx, long dy, 
			     long veclen, long z, long t )
{
  switch (filetype) {
  case FILE_DEFAULT:
    filter_default(temp_filtered, temp_out, temp_in, datatype_out,
		   datatype_in, dx_raw, dx, dy, veclen, z, t);
    break;
  case FILE_LX:
    filter_LX(temp_filtered, temp_out, temp_in, datatype_out,
	      datatype_in, dx_raw, dx, dy, veclen, z, t);
    break;
  case FILE_WINDAQ:
    filter_WINDAQ(temp_filtered, temp_out, temp_in, datatype_out,
		  datatype_in, dx, dy, veclen);
    break;
  }
}

static int read_binary( FILE* fpi, MRI_Datatype type_in, int n, void* buf, 
			int big_endian_input)
{
  switch (type_in) {
  case MRI_UNSIGNED_CHAR:
    FRdUInt8Array(fpi, (unsigned char*)buf, n);
    break;
  case MRI_SHORT:
    FRdInt16Array(fpi, (short*)buf, n);
    break;
  case MRI_INT:
    FRdInt32Array(fpi, (int*)buf, n);
    break;
  case MRI_FLOAT:
    FRdFloat32Array(fpi, (float*)buf, n);
    break;
  case MRI_DOUBLE:
    FRdFloat64Array(fpi, (double*)buf, n);
    break;
  default:
    Abort("%s: Internal error: unknown datatype %d\n",progname,
	  type_in);
  }
  return( bio_error==0 );
}

static void transfer( FILE* fpi, MRI_Dataset* DS, 
		      MRI_Datatype datatype_in, MRI_Datatype datatype_out,
		      long veclen, long dx_raw, 
		      long dx, long dy, long ds, long dz, long dt, 
		      long vxytz_flag, long reorder_flag, 
		      long start_offset, long skip, long sliceskip, 
		      char* readfile, int big_endian_input )
{
  long image_size_in;
  long image_size_out;
  long filtered_size_out;
  void* temp_in= NULL;
  void* temp_out= NULL;
  void* temp_filtered= NULL;
  long t, z, s, slice_number;
  size_t input_offset_increment;

  /* Allocate temporary space in which to read images;
     we assume that the binary data that we read in is big-endian and
     that it is in exactly the format that a Pittsburgh format chunk
     with MRI_Datatype == datatype would be in; under these assumptions
     we can transfer bytes directly to the chunk using the MRI_RAW mode */
  image_size_in = veclen * dx_raw * dy;
  image_size_out = veclen * dx_raw * dy;
  filtered_size_out = veclen * dx * dy;
  temp_in = emalloc( typesize[datatype_in]*image_size_in );
  temp_out = emalloc( typesize[datatype_out]*image_size_out );
  temp_filtered = emalloc( typesize[datatype_out]*filtered_size_out );
      
  if (vxytz_flag) {

    if (ds != 1) Abort("%s: multishot incompatible with vxytz order!\n");

    input_offset_increment = 
      (((typesize[datatype_in]*image_size_in) + sliceskip) * dt) + skip;
    
    /* reorder_flag is guaranteed false by earlier tests */
    /* Loop through slices */
    for (z = 0; z < dz; z++) {
      /* Move pointer within input files to next slice */
      efseek( fpi, start_offset + z * input_offset_increment,
	      SEEK_SET, readfile );
      for (t=0; t<dt; t++) {
	if (read_binary(fpi, datatype_in, image_size_in, temp_in,
			big_endian_input)
	    != 1)
	  Abort("%s: read failed on slice %d, time %d!\n",
		progname,z,t);
	if (sliceskip)
	  efseek( fpi, sliceskip, SEEK_CUR, readfile );
	
	move_and_filter( temp_filtered, temp_out, temp_in, 
			 datatype_out, datatype_in,
			 dx_raw, dx, dy, veclen, z, t );
	
	mri_set_chunk(DS, "images", filtered_size_out,
		      (z * dt + t) * filtered_size_out,
		      datatype_out, temp_filtered);
      }
    }
  }
  else {
    input_offset_increment = 
      (((typesize[datatype_in]*image_size_in) + sliceskip) * dz) + skip;
    
    /* Loop through time points */
    for (t = 0; t < dt; t++) {
      for (s = 0; s < ds; s++) {
	/* Move pointer within input file to next time point */
	efseek( fpi, start_offset + ((t*ds) + s) * input_offset_increment,
		SEEK_SET, readfile );
	
	/* Read in each slice and write to proper position */
	for (z = 0; z < dz; z++ ) {
	  if (read_binary(fpi, datatype_in, image_size_in, temp_in,
			  big_endian_input)
	      != 1)
	    Abort("%s: read failed on slice %d, time %d!\n",
		  progname,z,t);
	  if (sliceskip)
	    efseek( fpi, sliceskip, SEEK_CUR, readfile );
	  if (reorder_flag)
	    slice_number = ( z < (long) (( dz + 1 ) / 2 )) ?
	      ( z * 2 ) : ( z - (long) (( dz + 1 ) / 2 )) * 2 + 1;
	  else
	    slice_number = z;
	  
	  move_and_filter( temp_filtered, temp_out, temp_in, 
			   datatype_out, datatype_in,
			   dx_raw, dx, dy, veclen, z, t );
	  
	  mri_set_chunk(DS, "images", filtered_size_out,
			(t * dz + slice_number) * filtered_size_out,
			datatype_out, temp_filtered);
	}
      }
    }
  }
  free(temp_in);
  free(temp_out);
  free(temp_filtered);
}


int main(int argc, char* argv[])
{
  char readfile[512], hdrfile[512], datafile[512];
  char cdatatype[512], creorder[512], fulldatafile[512], cdataorder[512];
  char phasereffile[512], bpassdir[512], rampfile[512];
  char id_tag[512];
  int read_header= 0;
  int header_redundant_flag= 0;
  MRI_Datatype datatype_in;
  MRI_Datatype datatype_out;
  long veclen, dx, dy, ds, dz, dt;
  long dx_raw= 0;
  long start_offset, skip, sliceskip, reorder_flag, vxytz_flag;
  MRI_Dataset *DS = NULL;
  FILE *fpi = NULL;
  ScanInfo* scan_info= NULL;
  int big_endian_input= 1;


  /* Print version number */
  Message( "# %s\n", rcsid );

  progname= argv[0];

  /* Check to see if help was requested */
  if( ( argc > 1 ) && !strcmp( argv[1], "-help" ) )
    {
      if( argc == 2 )
	Help( "selecttopic" );
      else
	Help( argv[2] );
    }

  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Get filenames */
  cl_get( "input|i", "%option %s[%]", "input.mri", readfile );
  cl_get( "headerout|h", "%option %s[%]", "output.mri", hdrfile );
  cl_get( "dataout|da", "%option %s[%]", ".dat", datafile );
  cl_get( "tag", "%option %s[%]", "", id_tag );
  do_phaseref= cl_get( "phaseref","%option %s", phasereffile );
  do_bandpass= cl_get( "bandpass","%option %s", bpassdir );
  do_rampsample= cl_get( "rampfile","%option %s", rampfile );
  do_xchop= cl_present( "xchop" );
  do_ychop= cl_present( "ychop" );
  do_autoscale= cl_present( "autoscale" );
  debug= cl_present( "debug" );

  /* Let the user pick endianness if desired */
  if (cl_present("bigendian")) big_endian_input= 1;
  else if (cl_present("littleendian")) big_endian_input= 0;
  else big_endian_input= 1; /* raw fmri data is usually bigendian */

  /* Get flags */
  read_header= cl_present("readheader");

  if (read_header) {
    header_redundant_flag= 
      ( cl_get( "type|t", "%options %s", cdatatype )
	|| cl_get( "veclen|v", "%options %ld", &veclen )
	|| cl_get( "dims|di", "%option %ld %ld %ld %ld", &dx, &dy, &dz, &dt )
	|| cl_get( "offset|o", "%option %ld", &start_offset )
	|| cl_get( "skip|s", "%option %ld", &skip )
	|| cl_get( "sliceskip|ss","%option %ld", &sliceskip )
	|| cl_get( "reorder|r", "%option %s", &creorder )
	|| cl_get( "dataorder|do", "%option %s", "vxyzt", cdataorder ) );
  }
  else {
    /* Permit all flags containing info redundant with header */

    /* Get data-type, vector length, and dimension lengths */
    cl_get( "type|t", "%options %s[%]", "MRI_SHORT", cdatatype );
    cl_get( "veclen|v", "%options %ld[%]", 2, &veclen );
    cl_get( "dims|di", "%option %ld[1] %ld[1] %ld[1] %ld[1]", 
	    &dx, &dy, &dz, &dt );
    
    /* Get byte offsets parameters */
    cl_get( "offset|o", "%option %ld[%]", 0, &start_offset );
    cl_get( "skip|s", "%option %ld[%]", 0, &skip );
    cl_get( "sliceskip|ss","%option %ld[%]", 0, &sliceskip );
    
    /* Get EPI-reordering switch */
    cl_get( "reorder|r", "%option %s[%]", "0", &creorder );

    /* Get data order info */
    cl_get( "dataorder|do", "%option %s[%]", "vxyzt", &cdataorder );
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
  
  /* Test for valid filenames */
  if (!strcmp(readfile, datafile))
    Abort("%s: error: input and data output file names are the same.\n",
	  argv[0]);

  /* Test for compatible options */
  if ((!read_header) && do_bandpass)
    Abort("%s: error: bandpass flag requires readheader flag.\n",argv[0]);
  if ((!read_header) && do_phaseref)
    Abort("%s: error: phaseref flag requires readheader flag.\n",argv[0]);
  if ((!read_header) && do_rampsample)
    Abort("%s: error: rampfile flag requires readheader flag.\n",argv[0]);
  if (header_redundant_flag)
    Abort(
  "%s: error: a flag was given which is incompatible with the readheader.\n",
          argv[0]);

  bio_big_endian_input = big_endian_input;

  if (read_header) {
    filetype=check_filetype(readfile);
    if (filetype == FILE_DEFAULT) {
      /* Failed to find a match */
      Abort("%s: Could not parse the header of <%s>.\n",
	    progname,readfile);
    }
  }
  else filetype=FILE_DEFAULT;

  switch (filetype){
  case FILE_DEFAULT:
    process_default(creorder, cdataorder, cdatatype, 
		    &reorder_flag, &vxytz_flag, &datatype_in,
		    &datatype_out, &ds, &dx_raw, dx);
    break;
  case FILE_LX:
    process_LX(readfile, &datatype_in, 
       	      &veclen, &dx_raw, &dx, &dy, &ds, &dz, &dt, 
       	      &vxytz_flag, &reorder_flag, &start_offset,
       	      &skip, &sliceskip, scan_info, &datatype_out);

    break;
  case FILE_WINDAQ:
    process_WINDAQ(readfile, &datatype_in, 
       	      &veclen, &dx_raw, &dx, &dy, &ds, &dz, &dt, 
       	      &vxytz_flag, &reorder_flag, &start_offset,
       	      &skip, &sliceskip, &datatype_out);
    break;
  default:
    Abort("%s doesn't know how to handle %s-format input!\n",
	  progname,nameof_filetype(filetype));
  }

  /* Since process LX and WINDAQ may change bio_big_endian_input
   * make big_endian_input consistent
   */
  big_endian_input = bio_big_endian_input;

  /* Test for valid image specifications */
  if( (dx < 1) || (dy < 1) || (dz < 1) || (dt < 1) )
    Abort( "All dimension lengths must be specified (and positive) [%ld %ld %ld %ld].",
	   dx, dy, dz, dt );
  if( veclen < 1 )
    Abort( "Vector length must be positive [%ld].", veclen );
  if( (sliceskip < 0) || ( skip < 0 ) || ( start_offset < 0 ) )
    Abort("SliceSkip (%ld), Skip (%ld) and offset (%ld) must be non-negative.",
	   sliceskip, skip, start_offset );
  if( (sliceskip || skip) && reorder_flag )
    Warning( 1, "Reordering is being done with skips in reads." );
  /* Test for compatibility between reorder and dataorder */
  if (reorder_flag && vxytz_flag)
    Abort("Dataorder vtxyz is not compatible with reorder true.\n");
  
  /* Print out important parameters */
  if (id_tag[0])
    Message( "#      This data has the scan id <%s>\n",id_tag);
  if (filetype==FILE_LX) {
    Message( "#      This data is being read from a GE LX Pfile\n");
    Message( "#      Scan plane believed to be %s\n",plane);
  }
  Message( "#      Output Header-file: %s, Data-file: %s,\n",
	   hdrfile, datafile );
  Message( "#      Input Data-type: %s (%s), Output Data-type: %s\n",
	   typeconv[datatype_in], 
	   big_endian_input ? "bigendian" : "littleendian",
	   typeconv[datatype_out] );
  Message( "#      Vector Length: %ld, Dimensions: %ld %ld %ld %ld\n",
	   veclen, dx, dy, dz, dt );
  if (dx_raw!=dx)
    Message( "#      Unresampled X resolution: %d\n",dx_raw);
  if (ds>1) Message( "#      Shots: %ld\n",ds);
  if (scan_info) {
    Message("#      TR= %g ms, TE= %g ms\n",scan_info->TR,scan_info->TE);
    Message("#      Field of View X= %g mm, Y= %g mm\n",
	    scan_info->fov_x, scan_info->fov_y);
    Message("#      Image X= %g voxels, Y= %g voxels\n",
	    scan_info->image_x,scan_info->image_y);
    Message("#      Scan includes %d lines past Ky=0\n",
	    scan_info->overscan);
    Message("#      Voxel X= %g mm, Y= %g mm\n",
	    scan_info->voxel_x,scan_info->voxel_y);
    Message("#      Slice thickness %g mm, gap %g mm\n",
	    scan_info->slice_thickness,scan_info->slice_gap);
    Message("#      Scan acquired %s %s\n",scan_info->date,scan_info->time);
  }
  if (do_xchop) Message( "#      Data will be phase shifted pi/2 in X\n");
  if (do_ychop) Message( "#      Data will be phase shifted pi/2 in Y\n");
  if (do_autoscale)
    Message("#      Will autoscale to approximate range %g\n",
	    AUTOSCALE_RANGE);
  
  /* Read datasets used by corrections as requested */
  if (do_phaseref) {
    if (!(ref_lp= (float*)malloc(dy*ds*sizeof(float))))
      Abort("%s: unable to allocate %d bytes!\n",argv[0],dy*ds*sizeof(float));
    if (!(ref_poff= (float*)malloc(dy*ds*sizeof(float))))
      Abort("%s: unable to allocate %d bytes!\n",argv[0],dy*ds*sizeof(float));

    /* The following may turn off do_phaseref */
    ref_dat_read(ref_lp, ref_poff, dx_raw, dy*ds, &do_phaseref, phasereffile);
    if (!do_phaseref) 
      Warning(1,
  "%s: phase ref data not found or not valid; phase correction skipped!\n",
            argv[0]);

  }
  if (do_rampsample) {
    float* fbuf;

    if (!(fbuf=(float*)malloc(dx_raw*sizeof(float)))) 
      Abort("%s: unable to allocate %d bytes!\n",dx_raw*sizeof(float));

    ramp_vrgf_coeff= Matrix(dx_raw,dx,float);
      
    /* Read ramp sampling data.  The following may turn off do_rampsample. */
    vrgf_read(ramp_vrgf_coeff, fbuf, dx_raw, dx, &do_rampsample, rampfile);
    if (!do_rampsample)
      Abort("%s: ramp sample file <%s> is missing!\n",argv[0],rampfile);

    free(fbuf);
  }
  if (do_bandpass) {
    /* Read bandpass asymmetry correction data.  The following may turn
     * off do_bandpass.
     */
    read_bp(bandwidth,band_amag,band_aphase,&band_nbpc,
	    &do_bandpass,bpassdir);
    if (!do_bandpass) 
      Warning(1,
       "%s: Needed bandpass file missing; no bandpass correction!\n",argv[0]);
  }

  /* Open output dataset.  This will cause libbio to decide what type of
   * machine it's living on (via libmri's call to InitBIO).  We'll then
   * update that information based on command line info.
   */
  DS= open_output( hdrfile, datafile, datatype_out, vxytz_flag, veclen,
		   dx, dy, ds, dz, dt, scan_info, id_tag, argc, argv );
  bio_big_endian_input = big_endian_input;

  /* Open input file for reading */
  fpi = efopen( readfile, "r" );

  /* Transfer input to output */
  transfer( fpi, DS, datatype_in, datatype_out, veclen, dx_raw, 
	    dx, dy, ds, dz, dt, 
	    vxytz_flag, reorder_flag, start_offset, skip, sliceskip, 
	    readfile, big_endian_input );

  /* close input file */
  efclose( fpi );

  /* Close output dataset */
  mri_close_dataset( DS );

  fprintf(stderr,"#      Data converted to standard format.\n");
  return 0;
}

