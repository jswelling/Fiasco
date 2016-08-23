/************************************************************
 *                                                          *
 *	spiral_reader.c	                                    *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 2001 Department of Statistics,         *
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
 *  This program derived from "spiral" reconstruction code  *
 *  by Joel Welling, 9/2001 .                               *
 ************************************************************/

/* Notes-
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include "par.h"
#include "bio.h"
#include "array.h"
#include "acct.h"
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"
#include "spiral.h"

/* Defaults for several task parameters */
#define DEFAULT_OVER 2
#define DEFAULT_GRIDL 1.5

static char rcsid[] = "$Id: spiral_reader.c,v 1.4 2007/03/22 00:09:18 welling Exp $";

/* GLOBAL VARIABLES */
Context c;
Result r;
Task t;

/* VARIABLES FOR THIS MODULE */
int print_mode = 0;
FilenameList *thisname;
FilenameList *fname_list_root= NULL;
FilenameList *fname_list_tail= NULL;
char id_tag[512];

/* FORWARD DECLARATIONS */
static void InitMaster ();
static void ReadArguments (int argc, char **argv);
static void PrintMode ();
static void ProcessFirstFile ();
static void ProcessAdditionalFile ();
static void ProcessImageFile (int);
static void CreateOutputDataset (int argc, char** argv);
static void CalcTrajectory ();
static void UncompressFile (Filename out, const Filename in);
static void RemoveFile (const Filename name);


int
main (int argc,
      char **argv,
      char **envp)
{
  int i;
  int len;

  Acct(PROCESSING);
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
      if (uncomp_flag) RemoveFile(tname);
      thisname= thisname->next;
    }

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
  PrintAcct(argv[0], 0);
  exit(0);
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

  r.slice_sampim = NULL;

  /* Set some flags controlling behavior of worker routines.  Most of
   * these fields are never used, but they should be initialized to
   * something reasonable.
   */
  c.all_coils= TRUE;
  c.samp_cor= FALSE;
  c.loc_shift= TRUE;
  c.write_samples= TRUE;
  c.write_raw= FALSE;
  c.write_mag= FALSE;
  c.write_phase= FALSE;
  c.lin_cor= FALSE;
  c.hc_cor= FALSE;
  c.lin_map= FALSE;
  c.gen_map= FALSE;
  c.output_float= FALSE;
  c.res= 64;
  c.grid_len= DEFAULT_GRIDL;
  c.over_samp= DEFAULT_OVER;
  c.slice= 0;
  c.samp_delay= 50;
  c.ph_twist= 0.0;
  c.lr_shift= 0.0;
  c.tb_shift= 0.0;
  c.zoom= 1.0;
  c.mag_factor= 0.0;
  c.ph_factor= 0.0;
  c.filter_sz= 5.0;
  strcpy((char*)&(c.reg_file),"");
  strcpy((char*)&(c.reference_directory),".");
  strcpy((char*)&(c.output_data),".dat");
  c.output_scale= 1024.0;
  c.input_is_mri= FALSE;
  c.do_recon= FALSE;
}

static void
ReadArguments (int argc, char **argv)
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

  print_mode = cl_present("print_mode");

  /* next handle the parameterized switches */
  cl_get("input_dir|input_directory", "%option %s[%]", ".", c.input_directory);
  cl_get("tmp_dir|tmp_directory", "%option %s[%]", ".", c.tmp_directory);
  cl_get("output_dir|output_directory", "%option %s[%]", ".", c.output_directory);
  cl_get("out_name|output_name", "%option %s[%]",
	 "spiral.mri",c.output_name);
  cl_get("scale|output_scale", "%option %f[1024.0]", &c.output_scale);
  cl_get("v|verbosity", "%option %d[1]", &c.verbosity);
  cl_get( "tag", "%option %s[%]", "", id_tag );

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
}

static void
PrintMode ()
{
  char *s;
  int v;

  Report("\nProgram spiral_reader starting\n");
  Report("\n");

  Report("Input:\n");
  Report("\tInput directory: %s\n", c.input_directory);
  Report("\tTmp directory for uncompression: %s\n", c.tmp_directory);
  Report("\t%s-endian input\n", c.big_endian_input ? "Big" : "Little");
  Report("\n");
	
  Report("Parameters:\n");
  Report("\tResolution: %d\n", c.res);
  Report("\tOversampling: %d\n", c.over_samp);
  Report("\tGrid length: %f\n", c.grid_len);
  Report("\tSample delay: %d\n", c.samp_delay);
  Report("\tSample time: %f\n", c.samp_time);
  Report("\tMap delay: %f\n", c.mapdel);
  if (c.lin_map || c.gen_map)
    Report("\tFilter size: %d\n", c.filter_sz);
  Report("\n");

  Report("Output:");
  Report("\tOutput directory: %s\n", c.output_directory);
  Report("\tOutput dataset name: %s\n", c.output_name);
  Report("\tOutput datafile name or extension: %s\n", c.output_data);
  Report("\t%s-endian output\n", c.big_endian_output ? "Big" : "Little");
  Report("\tOutput format: %s\n", c.output_float ? "float32" : "int16");
  Report("\tOutput scaling factor: %f\n", c.output_scale);
  Report("\n");

  if (c.write_samples)
    Report("Writing ungridded K-space samples.\n");
  Report("\n");
}

static void
ProcessFirstFile (int argc, char** argv)
{
  ReadFileHeader(t.filename, &c, NULL);
  CalcTrajectory();
  CreateOutputDataset(argc, argv);
  WorkerInitialize();
  ProcessImageFile(0);
}

static void
ProcessAdditionalFile ()
{
  /* make sure this P file is consistent
     with the first P files that we read */
  CheckHeaderInfo(t.filename, &c);

  /* dispatch most of the work to the workers */
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
	  GenerateImage();
	}
    }
}

static void
CreateOutputDataset (int argc, char** argv)
{
  Filename fn;
  ConstructFilename(fn, c.output_directory, c.output_name);
  CreateMRIFile( fn, argc, argv, &c );
}

static void
CalcTrajectory ()
{
  int i, j, n;
  double kx,ky,kxo,kyo,ggx,ggy,ggm,kksr,kksi;
  FILE *f;
  Filename fn;
  int nsamp;
  int res;
  extern int getrttrajghg(int, int, double, double, double, double,
			  double, double, int, float**, float**);

  /* kaiser-bessel functions - initialize terms */
  c.gridb = PI*c.grid_len;
  c.wind = c.grid_len;

  res = getrttrajghg(c.opxres, c.npr, c.ts, c.gts, c.fsgcm, c.opfov,
		    c.slewrate, c.gts*21000, c.gtype, c.t2k[0], c.t2k[1]);

  if (c.verbosity >= VERBOSITY_FILE) {
    Report("\tndat = %d, npr = %d\n", c.ndat, c.npr);
    Report("\tts = %f, gts = %f, fsgcm = %f\n", c.ts, c.gts, c.fsgcm);
    Report("\topfov = %f, risetime = %d\n", c.opfov, c.risetime);
    Report("\tacq. matrix = %d, nom. resolution = %f mm\n",
	   c.opxres,10*c.opfov/c.opxres); 
    Report("\tk-space points = %d\n",c.res); 
    Report("\n");
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
  if (strncmp(out, "/afs", 4) == 0)
    {
      sprintf(com, "fs flushvolume -path %s", out);
      system(com);
    }
#endif
}

static void
RemoveFile (const Filename name)
{
  char com[1024];

  sprintf(com, "rm -f %s", name);
  system(com);
}
