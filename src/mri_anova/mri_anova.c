/************************************************************
 *                                                          *
 *  mri_anova.c                                             *
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
 *  Original programming by Joel Welling, 8/98              *
 ************************************************************/
/*************************************************************

  DESCRIPTION OF MRI_ANOVA

  mri_anova is a wrapper program for the general linear model
  facility provided in glm.c .  It takes a pgh MRI dataset of 
  type txyz and information about an experimental design, and 
  produces a number of datasets of type xyz containing F statistics 
  for the various design parameters.  vt... datasets are also 
  acceptable, so long as the vector length of v is 1.  The output 
  datasets are of type float32.

**************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include "mri.h"
#include "fmri.h"
#include "misc.h"
#include "stdcrg.h"

/* Notes-
 */

static char rcsid[] = "$Id: mri_anova.c,v 1.17 2007/08/30 20:18:54 welling Exp $";

#define MAX_SOURCES 20
#define CHUNKSIZE 4096
#define FACTOR_NAME_ABBREV_LENGTH 5
#define MAXCHARS BUFSIZ           /* maximum number of chars in lines */

typedef struct source_struct {
  char* name;
  long nfactors;
  long* factors;
} SourceDef;

typedef struct out_file_struct {
  char* source;
  MRI_Dataset* dset;
  double* chunk;
  long long chunk_offset;
  long long file_offset;
} OutFileStruct;

static char* progname= NULL;
static int nfactors= 0;
static char** factor_table;
static int nconditions= 0;
static sp_ConditionDef** cond_def_table= NULL;
static long nsources= 0;
static SourceDef* source_table= NULL; /* MAX_SOURCES in size */
static long nimages= 0;
static long nslices= 0;
static long nparams= 0;
static sp_SplitRec* split_table= NULL;
static OutFileStruct* out_table= NULL;
static OutFileStruct* mse_out= NULL;
static double* factors_unpacked= NULL;
static double* factors= NULL;
static double* factor_means= NULL;
static double* tseries= NULL;
static double* params= NULL;
static Regressor* gbl_r= NULL;

/* This routine generates an approprate dataset name from a 
 * source name.
 */
static char* fname_from_source( char* sourcename )
{
  static char buf[MAXCHARS];
  strcpy(buf,"Fmaps_");
  strncat(buf,sourcename,MAXCHARS-(strlen(buf)+1));
  buf[MAXCHARS-1]= '\0';
  return buf;
}

static void init_out_table( MRI_Dataset* in, long tblsize, int mse_flag,
			    int argc, char** argv )
{
  long i;
  OutFileStruct* out;
  char* key;

  if (!(out_table= (OutFileStruct*)malloc(tblsize*sizeof(OutFileStruct))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname, tblsize*sizeof(OutFileStruct));

  /* Build a prototype output ds in the first table element */
  out= out_table;
  out->source= source_table[0].name;
  if (!(out->chunk= (double*)malloc(CHUNKSIZE*sizeof(double))))
    Abort("%s: unable to allocate %d doubles!\n",progname,CHUNKSIZE);
  out->chunk_offset= 0;
  out->file_offset= 0;
  out->dset= mri_open_dataset( fname_from_source(out->source), MRI_WRITE );
  hist_add_cl( out->dset, argc, argv );
  mri_create_chunk( out->dset, "images" );
  mri_set_string( out->dset, "images.datatype", "float32" );
  mri_set_string( out->dset, "images.dimensions", "xyz" );
  mri_set_string( out->dset, "images.file", ".dat" );
  /* iterate over input dataset tags, discarding those we don't want */
  mri_iterate_over_keys(in);
  while (key= mri_next_key(in)) {
    if (!strcmp(key,"images.datatype")) { /* do nothing */ }
    else if (!strcmp(key,"images.dimensions")) { /* do nothing */ }
    else if (!strcmp(key,"images.file")) { /* do nothing */ }
    else if (!strcmp(key,"images.extent.t")) { /* do nothing */ }
    else if (!strcmp(key,"images.extent.v")) { /* do nothing */ }
    else if (!strcmp(key,"images.size")) { /* do nothing */ }
    else if (!strncmp(key,"missing",strlen("missing")))
      { /* do nothing- don't want to copy missing info */ }
    else {
      /* copy the key to new dataset */
      mri_set_string( out->dset, key, mri_get_string(in, key) );
    }
  }
  
  /* Now use this prototype dataset to create any others needed */
  for (i=1; i<tblsize; i++) {
    out= out_table+i;
    out->source= source_table[i].name;
    if (!(out->chunk= (double*)malloc(CHUNKSIZE*sizeof(double))))
      Abort("%s: unable to allocate %d doubles!\n",progname,CHUNKSIZE);
    out->chunk_offset= 0;
    out->file_offset= 0;
    out->dset = mri_copy_dataset( fname_from_source(out->source), 
				  out_table[0].dset );
  }

  /* If mse_flag is set, set up an output structure for MSE info */
  if (mse_flag) {
    if (!(mse_out= (OutFileStruct*)malloc(sizeof(OutFileStruct))))
      Abort("%s: unable to allocate %d bytes!\n",
	    progname, sizeof(OutFileStruct));
    mse_out->source= "MSE";
    if (!(mse_out->chunk= (double*)malloc(CHUNKSIZE*sizeof(double))))
      Abort("%s: unable to allocate %d doubles!\n",progname,CHUNKSIZE);
    mse_out->chunk_offset= 0;
    mse_out->file_offset= 0;
    mse_out->dset= mri_copy_dataset( "MSE", out_table[0].dset );
  }
}

static void flush_out_table(long tblsize)
{
  long i;
  OutFileStruct* out;

  for (i=0; i<tblsize; i++) {
    out= out_table+i;
    if (out->chunk_offset) {
      mri_set_chunk(out->dset, "images", out->chunk_offset,
		    out->file_offset, MRI_DOUBLE, out->chunk);
      out->file_offset += out->chunk_offset;
      out->chunk_offset= 0;
    }
  }
  if (mse_out) {
    if (mse_out->chunk_offset) {
      mri_set_chunk(mse_out->dset, "images", mse_out->chunk_offset,
		    mse_out->file_offset, MRI_DOUBLE, mse_out->chunk);
      mse_out->file_offset += mse_out->chunk_offset;
      mse_out->chunk_offset= 0;
    }
  }
}

static void close_out_table(long tblsize)
{
  long i;
  OutFileStruct* out;

  for (i=0; i<tblsize; i++) {
    out= out_table+i;
    mri_close_dataset(out->dset);
    free(out->source);
    free(out->chunk);
  }

  if (mse_out) {
    mri_close_dataset(mse_out->dset);
    free(mse_out->source);
    free(mse_out->chunk);
  }

  free(out_table);
  out_table= NULL;
}

static void gen_source_name( SourceDef* source )
{
  long i;

  if (!(source->name= 
	(char*)malloc((FACTOR_NAME_ABBREV_LENGTH+2)*source->nfactors)))
    Abort("%s: unable to allocate %d bytes!\n",progname,
	  (FACTOR_NAME_ABBREV_LENGTH+2)*source->nfactors);
  source->name[0]= '\0';
  for (i=0; i<source->nfactors; i++) {
    strncat(source->name,factor_table[source->factors[i]],
	    FACTOR_NAME_ABBREV_LENGTH);
    if (i!=(source->nfactors-1)) strcat(source->name,"_");
  }
}

static void build_source_table()
{
  long i;
  long j;
  long k;
  long factors_per_source;
  long crosses_this_level;
  long crosses_last_level;
  long runner= 0;
  long last_block= 0;

  if (nfactors>MAX_SOURCES) 
    Abort("%s: %d factors is greater than %d allowed sources!\n",
	  progname,nfactors,MAX_SOURCES);
  if ( !(source_table= (SourceDef*)malloc(MAX_SOURCES*sizeof(SourceDef))) )
    Abort("%s: unable to allocate %d bytes!\n",progname,
	  nsources*sizeof(SourceDef));

  factors_per_source= 1;
  crosses_this_level= nfactors;

  runner= 0;
  for (i=0; i<nfactors; i++) {
    source_table[runner].name= strdup(factor_table[i]);
    source_table[runner].nfactors= factors_per_source;
    if ( !(source_table[runner].factors= 
	   (long*)malloc(factors_per_source*sizeof(long))) )
      Abort("%s: unable to allocate %d longs!\n",progname,factors_per_source);
    source_table[runner].factors[0]= i;
    runner++;
  }
  nsources= runner;

  last_block= 0;
  while (1) {
    factors_per_source += 1;
    if (factors_per_source>nfactors) break;
    crosses_last_level= crosses_this_level;
    crosses_this_level= 
      crosses_last_level*(nfactors+1-factors_per_source)/factors_per_source;
    if (nsources+crosses_this_level > MAX_SOURCES) break;

    for (i=0; i<crosses_last_level; i++) {
      for (j=source_table[last_block+i]
	     .factors[source_table[last_block+i].nfactors-1]+1; 
	   j<nfactors; j++) {
	source_table[runner].nfactors= factors_per_source;
	if ( !(source_table[runner].factors= 
	       (long*)malloc(factors_per_source*sizeof(long))) )
	  Abort("%s: unable to allocate %d longs!\n",
		progname,factors_per_source);
	for (k=0; k<source_table[last_block+i].nfactors; k++)
	  source_table[runner].factors[k]= 
	    source_table[last_block+i].factors[k];
	source_table[runner].factors[factors_per_source-1]= j;
	gen_source_name(source_table+runner);
	runner++;
      }
    }
    last_block += nsources;
    nsources= runner;
  }

}

static void allocate_memory()
{
  if (!(factors_unpacked= (double*)malloc(nimages*nsources*sizeof(double))))
    Abort("%s: unable to allocate %d doubles!\n",nimages*nsources);

  if (!(factors= (double*)malloc(nimages*nsources*sizeof(double))))
    Abort("%s: unable to allocate %d doubles!\n",nimages*nsources);

  if (!(factor_means= (double*)malloc(nsources*sizeof(double))))
    Abort("%s: unable to allocate %d doubles!\n",nimages*nsources);

  if (!(tseries= (double*)malloc(nimages*sizeof(double))))
    Abort("%s: unable to allocate %d doubles!\n",nimages);

  if (!(params= (double*)malloc(nparams*sizeof(double))))
    Abort("%s: unable to allocate %d doubles!\n",nparams);
}

static void build_factor_matrix(long z, double* fac)
{
  long t;
  long s;

  for (t=0; t<nimages; t++) {
    long cond_index= split_table[t*nslices+z].cond;
    for (s=0; s<nfactors; s++) {
      fac[(t*nsources)+s]= 
	cond_def_table[cond_index]->factor_intlvl[s];
    }
  }
}

static void add_crosses_to_factor_matrix(long z, double* fac, 
					 long packed_length)
{
  double val;
  long t;
  long s;
  long i;

  for (s=nfactors; s<nsources; s++) {
    for (t=0; t<packed_length; t++) {
      val= 1.0;
      for (i=0; i<source_table[s].nfactors; i++) 
	val *= fac[(t*nsources)+source_table[s].factors[i]];
      fac[(t*nsources)+s]= val;
    }
  }
}

static long pack_out_missing( double* data_in, unsigned char** missing,
			      long tdim, long nblocks,
			      long z, double* data_out )
{
  long iblock;
  long t;
  double* in;
  double* out;
  long n_valid;

  /* This program works on data for which the block index (corresponding
   * to the factor) varies faster than the tdim index (corresponding to
   * observation).
   */
  n_valid= 0;
  for (iblock=0; iblock<nblocks; iblock++) {
    in= data_in + iblock;
    out= data_out + iblock;
    for (t=0; t<tdim; t++) {
      if (!missing[t][z]) {
	*out= *in;
	out += nblocks;
	n_valid++;
      }
      in += nblocks;
    }
  }

  return n_valid; /* note that this should be doubled for complex! */
}

static void mean_correct_factors( double* factors, double* factor_means,
				  long nfactors, long length )
{
  long ifactor;
  double sum;
  long i;
  double* runner;

  for (ifactor=0; ifactor<nfactors; ifactor++) {
    sum= 0.0;
    runner= factors + ifactor;
    for (i=0; i<length; i++) {
      sum += *runner;
      runner += nfactors;
    }
    sum /= length;
    factor_means[ifactor]= sum;
    runner= factors + ifactor;
    for (i=0; i<length; i++) {
      *runner -= sum;
      runner += nfactors;
    }
  }
}

static double mean_correct( double* tseries, long n )
{
  double mean= 0.0;
  long i;

  for (i=0; i<n; i++) mean += tseries[i];
  mean /= (double)n;
  for (i=0; i<n; i++) tseries[i] -= mean;

  return mean;
}

static double calc_mean_variance( double* tseries, long n, double mean )
{
  double sum= 0.0;
  long i;
  for (i=0; i<n; i++) sum += (tseries[i]-mean)*(tseries[i]-mean);
  return (sum/(double)(n-1));
}

static long data_identical( double* data, long n )
{
  long i;
  double val;

  val= data[0];
  for (i=1; i<n; i++) if (data[i]!=val) return 0;
  return 1;
}

/* for debugging */
static void dump_factor_matrix(double* fac)
{
  long t;
  long s;
  for (t=0; t<nimages; t++) {
    for (s=0; s<nsources; s++) 
      fprintf(stderr,"%f ",fac[(t*nsources)+s]);
    fprintf(stderr,"\n");
  }
}

static void write_voxel_data( long i, long j, long k, double mean, double variance,
			      double* params, long ndata )
{
  long loop;
  double* b;
  double* bvar;
  double* ssr;
  double* ssto;
  double* ortho_measure;
  long total_df;
  double sse;

  b= params;
  bvar= b+nsources;
  ssto= bvar + nsources;
  ssr= ssto + 2;
  ortho_measure= ssr+nsources;

  total_df= ndata-1;

  sse= *ssto;
  for (loop=0; loop<nsources; loop++) sse -= ssr[loop];

  printf("ANOVA for voxel %d %d %d: mean: %g, variance %g\n",
	 i, j, k, mean, variance);
  printf("      Factor orthogonality measure %f\n",*ortho_measure);
  printf("\n");
  printf("   Source               Estimate   Std. Err   df       SS       F\n");
  printf("   ------               --------   --------  ----     -----   -----\n");
  printf("\n");
  for (loop=0; loop<nsources; loop++) {
    double ftest= ((total_df-nsources)*ssr[loop])/sse;
    printf("   %-18s %10g %10g  %4d %10g  %5.3f\n",
	   source_table[loop].name,b[loop],sqrt(bvar[loop]),1,ssr[loop],
	   ftest);
  }
  printf("   %-18s                        %4d %10g\n",
	 "ERROR",total_df-nsources,sse);
  printf("   %-18s                        %4d %10g\n",
	 "CORR. TOTAL",total_df,*ssto);

  printf("\n");
}

static void output_voxel_data( double* params, int zero_flag, long ndata )
{
  long src_id;

  if (zero_flag) {
    for (src_id=0; src_id<nsources; src_id++) {
      out_table[src_id].chunk[out_table[src_id].chunk_offset]= 0.0;
      out_table[src_id].chunk_offset++;
    }
    if (mse_out) {
      mse_out->chunk[mse_out->chunk_offset]= 0.0;
      mse_out->chunk_offset++;
    }
  }
  else {
    double* ssr;
    double* ssto;
    long sse_df;
    double sse;
    
    ssto= params + 2*nsources;
    ssr= ssto + 2;
    sse_df= ndata-(nsources+1);
    
    sse= *ssto;
    for (src_id=0; src_id<nsources; src_id++) sse -= ssr[src_id];
    
    for (src_id=0; src_id<nsources; src_id++) {
      double ftest= (sse_df*ssr[src_id])/sse;
      out_table[src_id].chunk[out_table[src_id].chunk_offset]= ftest;
      out_table[src_id].chunk_offset++;
    }
    if (mse_out) {
      mse_out->chunk[mse_out->chunk_offset]= sse/(double)sse_df;
      mse_out->chunk_offset++;
    }
  }
}

int main( int argc, char** argv ) 
{
  MRI_Dataset *Input = NULL;
  long vec_length;
  long total_voxels;
  long voxels_moved;
  long voxels_this_chunk;
  long i;
  long j;
  long k;
  long t;
  char infile[512], outfile[512];
  char* dimstr= NULL;
  char* crunner;
  double* image_in= NULL;
  int v_flag= 0;
  long v_i;
  long v_j;
  long v_k;
  long xdim;
  long ydim;
  char split_fname[512];
  char cond_fname[512];
  unsigned char** missing= NULL;
  long tseries_length_packed= 0;
  long factors_length_packed= 0;
  double mean;
  double mean_variance;
  int mse_flag= 0;

  progname= argv[0];

  /* Print version number */
  Message( "# %s\n", rcsid );

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Deprecate old options */

  if (cl_present( "v" ))
     Abort ("Option v has been renamed voxel_address|vxa.  Please see help file.\n");
  if (cl_present( "cond" ))
     Abort ("Option cond has been renamed condition|cnd.  Please see help file.\n");

  v_flag= cl_get("voxel_address|vxa","%option %d %d %d",&v_i,&v_j,&v_k);

  cl_get( "split|spl", "%option %s[%]", "split", split_fname);
  cl_get( "condition|cnd", "%option %s[%s]", "conditions", cond_fname);
  mse_flag= cl_present("mse");

  /* Get filename */
  if (!cl_get("", "%s", infile)) {
    fprintf(stderr,"%s: Input file name not given.",argv[0]);
    Help( "usage" );
    exit(-1);
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

  /* Open input dataset */
  Input = mri_open_dataset( infile, MRI_READ );

  /* Check that program will function on data-set */
  if( !mri_has( Input, "images" ) )
    Abort( "%s operates only on standard images.", argv[0] );
  if ( !mri_has( Input, "images.dimensions" ) ) 
    Abort( "%s needs images.dimensions info for input file.", argv[0] );
  dimstr= mri_get_string( Input, "images.dimensions" );
  if (strcmp(dimstr,"txyz")) {
    if ( !strcmp(dimstr,"vtxyz") ) {
      if (!mri_has(Input,"images.extent.v")
	  || (mri_get_int(Input,"images.extent.v") != 1))
	Abort( "%s requires that vtxyz data have vector length 1.",argv[0]);
    }
    else Abort( "%s requires that input data be of type (v)txyz" );
  }

  if ( mri_has( Input, "images.extent.x" ) ) 
    xdim= mri_get_int(Input,"images.extent.x");
  else Abort("%s: input data file missing images.extent.x");
  if ( mri_has( Input, "images.extent.y" ) ) 
    ydim= mri_get_int(Input,"images.extent.y");
  else Abort("%s: input data file missing images.extent.y");
  if ( mri_has( Input, "images.extent.z" ) ) 
    nslices= mri_get_int(Input,"images.extent.z");
  else Abort("%s: input data file missing images.extent.z");
  if ( mri_has( Input, "images.extent.t" ) ) 
    nimages= mri_get_int(Input,"images.extent.t");
  else Abort("%s: input data file missing images.extent.t");

  missing= get_missing(Input);

  /* Parse condition file */
  if (!sp_parse_conditions( cond_fname, &factor_table, &nfactors,
			    &cond_def_table, &nconditions ))
      Abort("%s: fatal error parsing condition file %s!\n",
	    progname,cond_fname);

  /* Parse split file */
  if (!sp_parse_split( split_fname, nslices, nimages, nconditions,
		       &split_table ))
    Abort("%s: error processing split file <%s>!\n",
	  argv[0],split_fname);

  if (!sp_match_missing( split_table, missing, nslices, nimages ))
    Abort("%s: error synchronizing data and split missing info!\n",
	  argv[0]);

  /* Create source table, including crosses of split conditions */
  build_source_table();

  /* Initialize the glm package */
  if (!(gbl_r= glm_create_llsq_regressor()))
    Abort("%s: unable to create GLM regressor!\n",argv[0]);
  glm_set(gbl_r,GLM_COMPLEX,0);
  glm_set(gbl_r,GLM_RESIDUALS,0);
  glm_set(gbl_r,GLM_VARIANCES,1);
  glm_set(gbl_r,GLM_SSQR,1);
  glm_set(gbl_r,GLM_ORTHO,1);
  nparams= glm_n_params(gbl_r,nsources);

  /* Allocate memory */
  allocate_memory();

  if (v_flag) {
    /* Do one voxel and output the results */
    double* data;

    data= mri_get_chunk(Input, "images", 
			nimages, 
			((((v_k*ydim) + v_j)*xdim) + v_i)*nimages,
			MRI_DOUBLE);
    if (!data) Abort("%s: failed to read voxel %d %d %d!\n",
		     progname,v_i,v_j,v_k);

    build_factor_matrix(v_k, factors_unpacked);
    
    factors_length_packed= pack_out_missing(factors_unpacked, missing,
					    nimages, nsources, v_k, 
					    factors);
    mean_correct_factors( factors, factor_means, nsources,
			  factors_length_packed/nsources );
    add_crosses_to_factor_matrix(v_k, factors, factors_length_packed/nsources);

    tseries_length_packed= pack_out_missing(data, missing, nimages,
					    1, v_k, tseries);
    mean= mean_correct( tseries, tseries_length_packed );
    mean_variance= calc_mean_variance( tseries, tseries_length_packed, mean );

    if (data_identical(tseries, tseries_length_packed)) {
      printf("All data for voxel (%d %d %d) identical (value %f)\n",
	     v_i,v_j,v_k,data[0]);
    }
    else {
      if (glm_fit(gbl_r, tseries, factors, NULL, params, 
		  tseries_length_packed, nsources))
	Warning(1,"%s: error in glm_fit: %s!\n", progname, glm_error_msg());
      write_voxel_data(v_i, v_j, v_k, mean, mean_variance, params, 
		       tseries_length_packed);
    }
  }
  else {
    double* data;
    long long read_offset= 0;

    for (k=0; k<nslices; k++) {
      /* Pass anova the condition info */
      
      if (k==0) {
	/* Initialization: Set up output files */
	init_out_table( Input, nsources, mse_flag, argc, argv );
      }

      /* Build the factor matrix for this slice */
      build_factor_matrix(k, factors_unpacked);
      factors_length_packed= pack_out_missing(factors_unpacked, missing,
					      nimages, nsources, k, 
					      factors);
      mean_correct_factors( factors, factor_means, nsources,
			    factors_length_packed/nsources );
      add_crosses_to_factor_matrix(k, factors, factors_length_packed/nsources);

      for (j=0; j<ydim; j++) {
	for (i=0; i<xdim; i++) {
	  data= mri_get_chunk(Input, "images", 
			      nimages, read_offset, MRI_DOUBLE);
	  if (!data) 
	    Abort("%s: failed to read voxel %d %d %d!\n",progname,i,j,k);
	  read_offset += nimages;
	  
	  tseries_length_packed= pack_out_missing(data, missing, nimages,
						  1, k, tseries);
	  mean= mean_correct( tseries, tseries_length_packed );
	  mean_variance= 
	    calc_mean_variance( tseries, tseries_length_packed, mean );
	  
	  if (data_identical(tseries, tseries_length_packed)) {
	    output_voxel_data( NULL, 1, tseries_length_packed );
	  }
	  else {
	    if ( glm_fit(gbl_r, tseries, factors, NULL, params, 
			 tseries_length_packed, nsources) )
	      Warning(1,"%s: error in glm_fit at voxel %d %d %d: %s!\n", 
		      progname, i, j, k, glm_error_msg());
	    output_voxel_data( params, 0, tseries_length_packed );
	  }
	  if (out_table[0].chunk_offset == CHUNKSIZE)
	    flush_out_table( nsources );
	}
      }
      printf("Slice %d complete\n",k);
      fflush(stdout);
    }
    flush_out_table( nsources );
    close_out_table( nsources );
  }
  mri_close_dataset( Input );

  glm_destroy(gbl_r);

  Message( "#      Anova calculations complete.\n" );

  return 0;
}

