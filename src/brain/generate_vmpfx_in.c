/************************************************************
 *                                                          *
 *  generate_vmpfx_in.c                                     *
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
 *  Original programming by Joel Welling 9-97               *
 ************************************************************/
/*************************************************************

  DESCRIPTION OF generate_vmpfx_in.c
  This program takes a prototype file and some Fiasco run 
  information and generates an appropriate input file for
  vmpfx, a component of Chris Genovese' BRAIN toolkit.  The
  output file is written to stdout.

  generate_vmpfx_in -proto proto_file -split split_file -cond cond_file
             -droot input_file_root_name -dpath input_file_path 
             -iai iai_value -fixed fixed_cond_names 
	     -nimage num_images -nslice num_slices_per_image

    -proto prototype file name (default "vmpfx_proto.t")
    -split image-by-image split file as output by intsplit (default "newsplit")
    -cond  condition file as output by intsplit (default "conditions")
    -droot root file name for input dataset (required)
    -dpath path for input dataset (required)
    -iai   inter-acquisition interval value (default 3.0)
    -fixed names of factor levels to hold fixed; this is required
    -nimage number of images (required)
    -nslice slices per image (required)
    -null  flag indicating that only the null model is to be run

**************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "misc.h"
#include "fmri.h"
#include "stdcrg.h"

#define INBUF_LENGTH 512
#define STRING_LENGTH 512
#define INITIAL_MAX_CONDITIONS 32

/*******
 * Experience has shown that vmpfx often becomes unstable if the number
 * of experimental conditions is too high; the limit and mechanism are
 * not well understood.  This is a fairly arbitrary limit for a warning.
 * vmpfx crashes if this problem happens.
 *******/
#define MAX_SAFE_CONDITIONS 25

/*******
 * This is an empirical factor that controls the number of knots used
 * to fit the drift curve to the data.
 *******/
#define IMAGES_PER_KNOT 10

typedef struct cond_def_struct {
  int id;
  int fixed;
  char** factor_lvl;
} ConditionDef;

typedef struct split_rec_struct {
  int image;
  int slice;
  int cond;
  float start_time;
} SplitRec;

static char rcsid[] = "$Id: generate_vmpfx_in.c,v 1.15 2007/04/19 22:37:33 welling Exp $";

static char* progname;

static char inbuf[INBUF_LENGTH];

static char data_root[STRING_LENGTH];
static char data_path[STRING_LENGTH];
static char iai_string[STRING_LENGTH];
static char fixed_conditions[STRING_LENGTH];
static char fixed_cond_out_string[STRING_LENGTH];
static char condition_map_string[STRING_LENGTH];

static int nslices= 0;
static int nimages= 0;
static float iai= 0.0;

static int n_factors= 0;
static char** factor_table;
static int n_conditions= 0;
static int cond_def_table_size= 0;
static ConditionDef** cond_def_table= NULL;
static int n_fixed_conditions= 0;
static int some_images_NA= 0;
static int null_model_only= 0;

static SplitRec* split_table= NULL;

static void check_cond_def_table_size( int n_cond )
{
  if (!cond_def_table) {
    if (!(cond_def_table= (ConditionDef**)malloc(INITIAL_MAX_CONDITIONS
						 *sizeof(ConditionDef*))))
      Abort("%s: unable to allocate %d bytes!\n",progname,
	    INITIAL_MAX_CONDITIONS*sizeof(ConditionDef*));
    cond_def_table_size= INITIAL_MAX_CONDITIONS;
  }

  if (n_cond >= cond_def_table_size) {
    ConditionDef** tmp;
    int i;
    if (!(tmp= (ConditionDef**)malloc(2*cond_def_table_size
				      *sizeof(ConditionDef*))))
      Abort("%s: unable to allocate %d bytes!\n",progname,
	    2*cond_def_table_size*sizeof(ConditionDef*));
    for (i=0; i<cond_def_table_size; i++) 
      tmp[i]=cond_def_table[i];
    free(cond_def_table);
    cond_def_table= tmp;
    cond_def_table_size= 2*cond_def_table_size;
  }
}

static int process_cond_file( FILE* cond_file )
{
  char* tstring;
  int i;
  int j;

  /* Read header line */
  if (!fgets(inbuf,INBUF_LENGTH,cond_file)) {
    perror("process_cond_file read error");
    Abort("%s: failed to read conditions file!\n", progname);
  }

  /* Count factors */
  if (!(tstring= (char*)malloc(strlen(inbuf)+1)))
    Abort("%s: failed to allocate %d bytes!\n", progname, strlen(inbuf)+1);
  strcpy(tstring, inbuf);
  if (strtok(tstring," \t\n")) {
    n_factors= 1;
    while (strtok(NULL," \t\n")) n_factors++;
  }
  free(tstring);

  /* Read factors */
  if (!(factor_table= (char**)malloc(n_factors*sizeof(char*))))
    Abort("%s: failed to allocate %d bytes!\n", progname, 
	  n_factors*sizeof(char*));
  for (i=0; i<n_factors; i++) {
    factor_table[i]= strdup(strtok((i ? NULL : inbuf)," \t\n"));
  }

  /* Read conditions */
  n_conditions= 0;
  while (!feof(cond_file) && !ferror(cond_file)) {
    check_cond_def_table_size(n_conditions);
    if (!(cond_def_table[n_conditions]= 
	  (ConditionDef*)malloc(sizeof(ConditionDef))))
      Abort("%s: failed to allocate %d bytes!\n", progname,
	    sizeof(ConditionDef));
    if (!fgets(inbuf,INBUF_LENGTH,cond_file)) break;
    if (ferror(cond_file)) break;
    tstring= strtok(inbuf," \t\n");
    if (tstring) {
      cond_def_table[n_conditions]->id= atoi(tstring);
      cond_def_table[n_conditions]->fixed= 0;
      if (!(cond_def_table[n_conditions]->factor_lvl= 
	    (char**)malloc(n_factors*sizeof(char*))))
	Abort("%s: failed to allocate %d bytes!\n", progname,
	      n_factors*sizeof(char*));
      for (j=0; j<n_factors; j++) {
	tstring= strtok(NULL," \t\n");
	if (!tstring) 
	  Abort("%s: invalid condition file line <%s>\n",progname,inbuf);
	cond_def_table[n_conditions]->factor_lvl[j]= strdup(tstring);
      }
      n_conditions++;
    }
  }

  if (ferror(cond_file)) {
    perror("Error reading conditions");
    Abort("%s: Error reading conditions file!\n",progname);
  }

  return 1;
}

static int parse_fixed_cond(char* fixed_cond_str)
{
  int ret_code= 1;
  char* this_cond;
  

  /* The missing condition is always fixed by assertion */
  cond_def_table[0]->fixed= 1;
  n_fixed_conditions++;

  this_cond= strtok(fixed_cond_str," ,");
  while (this_cond) {
    int found_this= 0;
    int i;
    for (i=0; i<n_conditions; i++) {
      int j;
      for (j=0; j<n_factors; j++) {
	if (!strcmp(this_cond,cond_def_table[i]->factor_lvl[j])) {
	  found_this= 1;
	  if (!cond_def_table[i]->fixed) {
	    n_fixed_conditions++;
	    cond_def_table[i]->fixed= 1;
	  }
	}
      }
    }
    if (!found_this) {
      Error("%s: fixed condition name <%s> does not match.\n",
	    progname,this_cond);
      ret_code= 0;
    }
    this_cond= strtok(NULL," ,");
  }
  return ret_code;
}

static int check_and_substitute( char* instring, char* test, char* replace )
{
  char* tok_start;
  char* runner;
  int offset;
  int in_length;
  int test_length;

  if ((tok_start= strstr(instring,test)) != NULL) {
    in_length= strlen(instring);
    test_length= strlen(test);

    /* check for overflow */
    if (in_length + strlen(replace) - test_length > INBUF_LENGTH-1)
      Abort("%s: string buffer overflow;  substituted line too long!\n",
	    progname);

    /* shift rest of line */
    /* if string lengths are equal, no shift is needed */
    offset= strlen(replace) - test_length;
    if (offset>0) {
      for (runner= instring+in_length; runner>=tok_start+test_length; 
	   runner--)
	*(runner+offset)= *runner;
    }
    else if (offset<0) {
      for (runner= tok_start+test_length; runner<=instring+in_length; 
	   runner++)
	*(runner+offset)= *runner;
    }

    /* substitute */
    runner= tok_start;
    for (; *replace; replace++)
      *runner++= *replace;

    return 1;
  }

  else return 0;
}

static int split_rec_compare( const void* r1_in, const void* r2_in )
{
  SplitRec* r1= (SplitRec*)r1_in;
  SplitRec* r2= (SplitRec*)r2_in;

  if (r1->start_time < r2->start_time) return -1;
  else if (r1->start_time > r2->start_time) return 1;
  else return 0;
}

static int process_split_file( FILE* split_file )
{
  int i;

  float slice_time= iai/nslices;

  if (!(split_table= (SplitRec*)malloc(nimages*nslices*sizeof(SplitRec))))
    Abort("%s: failed to allocate %d bytes!\n",
	  progname, nimages*nslices*sizeof(SplitRec));

  /* Load the split records, add acquisition start times */
  i= 0;
  while (!feof(split_file) && !ferror(split_file)) {
    if (i>=nimages*nslices) break; /* in case we're fed an over-length file */
    if ( fscanf(split_file,"%d %d %d\n", &(split_table[i].image), 
		&(split_table[i].slice), &(split_table[i].cond)) != 3 ) {
      Error("%s: error parsing split file!\n",progname);
      return 0;
    }
    split_table[i].start_time= 
      ( split_table[i].image * iai ) 
      + ((split_table[i].slice / 2) * slice_time)
      + ((split_table[i].slice % 2) * (0.5*iai));
    i++;
  }

  if (ferror(split_file)) {
    Error("%s: error reading split file!\n",progname);
    return 0;
  }

  /* Now cleverly sort the records by start time */
  qsort((void*)split_table, nimages*nslices, sizeof(SplitRec), 
	split_rec_compare);

  return 1;
}

static int test_and_filter_conditions()
{
  int* cond_guesses;
  int i;
  int icond;
  char* runner= fixed_cond_out_string;

  /* Generate consistent guesses about the conditions in effect
   * for each image.
   */
  if (!(cond_guesses= (int*)malloc(nimages*sizeof(int))))
    Abort("%s: unable to allocate %d bytes!\n",progname,nimages*sizeof(int));
  for (i=0; i<nimages; i++) cond_guesses[i]= 0;
  for (icond=0; icond<nimages*nslices; icond++) {
    if (split_table[icond].cond != 0) {
      int image= split_table[icond].image;

      if (cond_guesses[image]==0) cond_guesses[image]= split_table[icond].cond;
      else {
	/* The following test determines if we have already seen 
	 * a slice in this image with a non-NA cond.  This would
	 * violate the requirements for vmpfx which we choose to
	 * enforce.
	 */
	if (cond_guesses[image] != split_table[icond].cond) return 0;
      }
    }
  }

  /* Go through and refill any slices marked NA (condition 0)
   * with the value appropriate for the rest of that image.
   * This allows some NAs to persist, but for the whole image
   * only.
   */
  for (icond=0; icond<nimages*nslices; icond++) {
    if (split_table[icond].cond==0)
      split_table[icond].cond= cond_guesses[ split_table[icond].image ];
  }

  /* Make a note of whether or not any trace of the NA
   * condition remains.
   */
  some_images_NA= 0;
  for (i=0; i<nimages; i++) 
    if (cond_guesses[i]==0) some_images_NA= 1;

  /* Make sure we have enough fixed conditions to get by */
  if ((some_images_NA && n_fixed_conditions<1)
      || (!some_images_NA && n_fixed_conditions<2))
    Abort("%s: at least one fixed condition must exist in experiment!\n",
	  progname);

  /* Build the output list of fixed conditions, including
   * or excluding the NA condition as appropriate.
   */
  runner= fixed_cond_out_string;
  *runner= 0; /* null terminate */
  if (cond_def_table[0]->fixed && some_images_NA) {
    sprintf(runner,"0 ");
    runner += 2;
  }
  for (i=1; i<n_conditions; i++) {
    if (cond_def_table[i]->fixed) {
      if ((int)(runner-fixed_cond_out_string) > STRING_LENGTH-64)
	Abort("%s: fixed cond string overflow!\n",progname);
      if (runner != fixed_cond_out_string) *runner++= ' ';
      sprintf(runner,"%d", some_images_NA ? i : i-1);
      runner += strlen(runner);
    }
  }

  /* Build a list of Fiasco condition numbers corresponding to
   * the valid range of Brain condition numbers.  Format is
   * a space-separated list of f/(brainCondNumber).(fiascoCondNumber)
   * triples, where f is replaced by v for non-fixed conditions.  
   * It needn't include all fiascoCondNumbers, since the NA 
   * condition may be dropped.  The first brainCondNumber
   * will always be zero.
   */
  runner= condition_map_string;
  *runner= 0; /* null terminate */
  if (some_images_NA) {
    for (i=0; i<n_conditions; i++) {
      if ((int)(runner-condition_map_string) > STRING_LENGTH-64)
	Abort("%s: map cond string overflow!\n",progname);
      sprintf(runner,"%c/%d.%d ",
	      (cond_def_table[i]->fixed ? 'f' : 'v'),
	      i,i);
      runner += strlen(runner);
    }
  }
  else {
    for (i=1; i<n_conditions; i++) {
      if ((int)(runner-condition_map_string) > STRING_LENGTH-64)
	Abort("%s: map cond string overflow!\n",progname);
      sprintf(runner,"%c/%d.%d ",
	      (cond_def_table[i]->fixed ? 'f' : 'v'),
	      i-1,i);
      runner += strlen(runner);
    }
  }

  free(cond_guesses);

  return 1;
}

static int emit_condition_table()
{
  int current_condition= split_table[0].cond;
  int cond_out;
  int current_start;
  int i;
  
  current_start= 0;
  for (i=1; i<nimages*nslices; i++) {
    if (split_table[i].cond != current_condition) {
      cond_out= (some_images_NA) ? 
	current_condition : current_condition-1;
      printf("%d   %d   %d\n", cond_out, 
	     current_start, split_table[i].image-current_start);
      current_start= split_table[i].image;
      current_condition= split_table[i].cond;
    }
  }

  cond_out= (some_images_NA) ? current_condition : current_condition-1;
  printf("%d   %d   %d\n", cond_out, current_start,
	 split_table[(nimages*nslices)-1].image + 1 - current_start);
  return 1;
}

static int process_line( char* instring )
{
  if (!instring || !(*instring)) return 1; /* sometimes we get sent nulls */

  /* Check for conditions forbidden with the -null flag */
  if (null_model_only &&
      (strstr(instring,"****CONDITION_TABLE****")
       || strstr(instring,"****NUM_CONDITIONS****")
       || strstr(instring,"****FIXED_CONDITIONS****")))
    Abort("%s: prototype file is inappropriate for -null flag\n",progname);

  if (strstr(instring,"****CONDITION_TABLE****")) {
    return(emit_condition_table());
  }
  else {
    char tstring[512];

    (void)check_and_substitute(instring,"****DATA_ROOT****",data_root);
    (void)check_and_substitute(instring,"****DATA_PATH****",data_path);
    (void)check_and_substitute(instring,"****IAI****",iai_string);

    sprintf(tstring,"%d",nimages);
    (void)check_and_substitute(instring,"****NUM_IMAGES****",tstring);

    sprintf(tstring,"%d",nslices);
    (void)check_and_substitute(instring,"****NUM_SLICES****",tstring);

    if (some_images_NA) sprintf(tstring,"%d",n_conditions);
    else sprintf(tstring,"%d",n_conditions-1);
    (void)check_and_substitute(instring,"****NUM_CONDITIONS****",tstring);

    (void)check_and_substitute(instring,"****FIXED_CONDITIONS****",
			       fixed_cond_out_string);
    (void)check_and_substitute(instring,"****CONDITION_MAP****",
			       condition_map_string);
    sprintf(tstring,"%d",nimages/IMAGES_PER_KNOT);
    (void)check_and_substitute(instring,"****NUM_KNOTS****",tstring);

    fputs(instring,stdout);
  }
  return 1;
}

int main( int argc, char *argv[] )
{
  char in_fname[512];
  char split_fname[512];
  char cond_fname[512];
  FILE* infile= NULL;
  FILE* cond_file= NULL;
  FILE* split_file= NULL;

  progname= argv[0];

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

  cl_get( "proto", "%option %s[%]", "vmpfx_proto.t", in_fname );
  if (!cl_get( "droot", "%option %s", data_root)) {
    fprintf(stderr,"%s: -droot option required.\n",argv[0]);
    Help( "usage" );
    exit(-1);
  }
  if (!cl_get( "dpath", "%option %s", data_path)) {
    fprintf(stderr,"%s: -dpath option required.\n",argv[0]);
    Help( "usage" );
    exit(-1);
  }
  if (!cl_get( "nimage", "%option %d", &nimages )) {
    fprintf(stderr,"%s: -nimage option required.\n",argv[0]);
    Help( "usage" );
    exit(-1);
  }
  if (!cl_get( "nslice", "%option %d", &nslices )) {
    fprintf(stderr,"%s: -nslice option required.\n",argv[0]);
    Help( "usage" );
    exit(-1);
  }
  cl_get( "iai", "%option %s[%]", "3.0", iai_string);
  null_model_only= cl_present("null");
  if (!null_model_only) {
    cl_get( "split", "%option %s[%]", "split", split_fname);
    cl_get( "cond", "%option %s[%s]", "conditions", cond_fname);
    if (!cl_get( "fixed", "%option %s", fixed_conditions)) {
      fprintf(stderr,"%s: -fixed option required.\n",argv[0]);
      Help( "usage" );
    }
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

  if ((iai=atof(iai_string)) <= 0.0)
    Abort("%s: invalid iai value %s\n",iai_string);

  if (!null_model_only) {

    /* Parse condition file */
    if (!(cond_file= fopen(cond_fname,"r")))
      Abort("%s: unable to open condition file <%s>!\n",
	    argv[0], cond_fname);
    
    if (!process_cond_file(cond_file)) 
      Abort("%s: error processing conditions file <%s>!\n",
	    argv[0],cond_fname);
    
    if (n_conditions>MAX_SAFE_CONDITIONS)
      fprintf(stderr,
	      "%s: *****WARNING***** vmpfx sometimes fails to converge above %d conditions; you have %d!\n",
	      argv[0], MAX_SAFE_CONDITIONS, n_conditions);
    
    if (fclose(cond_file))
      Error("%s: failed to close conditions file!\n",argv[0]);
    
    /* Note fixed conditions ***/
    if (!parse_fixed_cond(fixed_conditions))
      Abort("%s: failed to parse fixed condition list\n",
	    progname);
    
    /* Parse split file */
    if (!(split_file= fopen(split_fname,"r")))
      Abort("%s: unable to open split file <%s>!\n",argv[0],split_fname);
    
    if (!process_split_file(split_file))
      Abort("%s: error processing split file <%s>!\n",argv[0],split_fname);
    
    if (fclose(split_file))
      Error("%s: failed to close split file!\n",argv[0]);
    
    if (!test_and_filter_conditions()) 
      Abort("%s: this experimental design is not compatible with vmpfx!\n",
	    argv[0]);
  }

  /* Translate and substitute prototype file */
  if (!(infile= fopen(in_fname,"r")))
    Abort("%s: unable to open prototype file <%s>\n",argv[0],in_fname);

  while (!feof(infile) && !ferror(infile)) {
    process_line( fgets(inbuf, INBUF_LENGTH, infile) );
  }

  if (ferror(infile))
    Error("%s: error while reading prototype file!\n", argv[0]);

  if (fclose(infile))
    Error("%s: failed to close prototype file!\n", argv[0]);

  return(0);
}
