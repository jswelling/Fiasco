/************************************************************
 *                                                          *
 *  eventpeaks.c                                           *
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
 *  Original programming by Tom Bagby, 12/98                *
 ************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mri.h"
#include "fmri.h"
#include "misc.h"
#include "stdcrg.h"
#include "polyfit.h"

/* Empirically derived constants */
#define DEFAULT_WINDOW 20
#define DEFAULT_THRESHOLD 200.0
#define KEYBUF_SIZE 512
#define MARK_THRESHOLD 20
#define PEAK -1
#define END_OF_FILE -2
#define BEFORE_FIRST_MARK -99

/* Global parameters */
int window;
float threshold;

int find_next_peak(float **, int, int, int);
void process_data(float *, float **, int);


static void safe_copy(char* str1, char* str2) {
  strncpy(str1, str2, KEYBUF_SIZE);
  str1[KEYBUF_SIZE-1]= '\0';
}


static void safe_concat(char* str1, char* str2) {
  strncat(str1, str2, (KEYBUF_SIZE-strlen(str1))-1);
}


static int safe_get_extent(MRI_Dataset* ds, char* chunk, char* dim)
{
  char key_buf[KEYBUF_SIZE];
  char dim_buf[4];
  int ext;

  dim_buf[0]= *dim;
  dim_buf[1]= '\0';
  safe_copy(key_buf,chunk);
  safe_concat(key_buf,".extent.");
  safe_concat(key_buf,dim_buf);

  if (mri_has(ds,key_buf)) {
    ext=mri_get_int(ds, key_buf);
    return mri_get_int(ds,key_buf);
  }
  else Abort("input file missing tag %s!\n",key_buf);
  return 0; /* not reached */
}


/* process_data:							*/
/* does all the actual data processing taking the sync data in the 	*/
/* input buffer and making our cooresponding ramping function in the 	*/
/* output buffer where n is the size of the buffers.			*/

void process_data(float *in_buffer, float **out_buffer, int n) {
	int i, j, k, l, w2 = window/2;
	int first_peak[2], last_peak[2];

	float *temp_buffer = (float *)malloc(n * sizeof(float));

	/* Set the temp and output buffer contents to 0 */
	for(i = 0; i < n; i++)
		temp_buffer[i] = 0;

	for(i = 0; i < 2; i++)	
	  for(j = 0; j < n; j++)
		out_buffer[i][j] = 0;


	/* Finding peaks by using an approximation for the first    */
	/* derivative.                                              */

	for(i = w2; i < n - w2; i++)
	  temp_buffer[i] = (in_buffer[i+1] - in_buffer[i-1])/2; 

	for(i = w2; i < n - w2; i++) {
		if(temp_buffer[i-1] > temp_buffer[i] && 
		   temp_buffer[i+1] > temp_buffer[i] &&
		   temp_buffer[i] < -threshold) {

		  j=i;

		  while ((out_buffer[i][0]==0) && (out_buffer[i][1]==0)
			 && (j < n - w2)) {
		    if(temp_buffer[j-1] < temp_buffer[j] && 
		       temp_buffer[j+1] < temp_buffer[j] &&
		       temp_buffer[j] > threshold) {
		      if (j-i > MARK_THRESHOLD) out_buffer[i][0] = PEAK;
		      else out_buffer[i][1] = PEAK;
		    }
		  j++;
		  }
		  i=j;
		}
	}



	/* Now linearize everything */


	for (l = 0; l < 2; l++){
	  i = find_next_peak(out_buffer, -1, n, l);
	  first_peak[l]=i;
	  while((j = find_next_peak(out_buffer, i, n, l)) != END_OF_FILE) {
	    for(k = i; k < j; k++) 
	      out_buffer[k][l] = k-i;
	    i=j;
	  }
	  last_peak[l]=i;
	}

	/* Now we have to take care of the sections at the beginning and
	   end of the data that just didn't fit into our window.
	   we copy the slope from the first and last sections where we were
	   able to identify peaks, and fill in the gaps.  */

	for(l = 0; l < 2; l++){
	  for(k = 0; k < first_peak[l]; k++) 
	    out_buffer[k][l] = BEFORE_FIRST_MARK;
	}


	for(l = 0; l < 2; l++){
	  for(j = n - 1; j > last_peak[l]; j--)
	    out_buffer[j][l] = j-i;
	}

}

/* find_next_peak:						*/
/* given the location of the previous peak, finds the next 	*/
/* peak.  if none is to be found, returns END_OF_FILE		*/

int find_next_peak(float **buffer, int start, int size, int v) {
  int i = start + 1;
  while(i < size && buffer[i][v] != PEAK) i++;
  if(i == size) return END_OF_FILE;
  return i;
}

int main(int argc, char *argv[]) {
  MRI_Dataset *input, *output;
  int chunk_size, num_elements, i, j, k;
  float *input_buffer, *output_space;
  float **output_buffer;
  char chunk[512];
  char infile[512];
  char outfile[512];
  char tag[600];
  char t[512]="t";
  char *dimstr;
  int dim_extent;

  /* Check to see if help was requested */
  if( ( argc > 1 ) && !strcmp( argv[1], "-help" ) )
    {
      if( argc == 2 )
	Help( "selecttopic" );
      else
	Help( argv[2] );
    }

  /* First, parse the command line */

  cl_scan(argc, argv);

  /* Then, get the options, input and output filenames */

  cl_get("chunk|chu|c", "%option %s[%]", "images", chunk);
  cl_get("window|win", "%option %d[%]", DEFAULT_WINDOW, &window);
  cl_get("threshold|thr", "%option %f[%]", DEFAULT_THRESHOLD, &threshold); 

  if(!cl_get("", "%s", infile)) {
    fprintf(stderr, "%s: Input file name not given.\n", argv[0]);
    exit(-1);
  }
  if(!cl_get("", "%s", outfile)) {
    fprintf(stderr, "%s: Output file name not given.\n", argv[0]);
    exit(-1);
  }

  if(cl_cleanup_check()) {
    int i;
    fprintf(stderr, "%s: invalid argument in command line:\n", argv[0]);
    for(i = 0; i < argc; i++) fprintf(stderr, "%s ", argv[i]);
    fprintf(stderr, "\n");
    Help("usage");
    exit(-1);
  }

  /* Now that we are done the command-line parsing, lets copy the
     whole dataset and open everything up so we can start messing
     around with it. If the named chunk isn't in the dataset,
     then we abort. We also abort if the input and output filenames
     are the same */

  if(!strcmp(infile, outfile)) {
    fprintf(stderr, "%s: input and output files must be distinct.\n", argv[0]);
    exit(-1);
  }
  input = mri_open_dataset(infile, MRI_READ);
  if(!mri_has(input, chunk)) {
    fprintf(stderr, "%s: requested chunk %s not present.\n", argv[0], chunk);
    exit(-1);
  }


  strcpy(tag, chunk);
  strcat(tag, ".dimensions");
  if (strcmp(dimstr=mri_get_string(input,tag), "t")) {
    while (*dimstr){
      if ((*dimstr != *"t") && (dim_extent=(safe_get_extent(input, chunk, dimstr)>1))) {
	fprintf(stderr,"%s: input chunk should only have dimension 't' greater than 1.\n", argv[0]);
	exit(-1); 
      }
      else dimstr++;
    }
  }

  output = mri_copy_dataset(outfile, input);
  hist_add_cl( output, argc, argv );

  strcpy(tag, chunk);
  strcat(tag, ".datatype");
  mri_set_string( output, tag, "float32");

  /* Set v extent of output to 2 (one for each type of event marker) */
  strcpy(tag, chunk);
  strcat(tag, ".extent.v");
  mri_set_int(output, tag, 2);

  strcpy(tag, chunk);
  strcat(tag, ".dimensions");
  mri_set_string(output, tag, "vt");

  /* Get the length of the chunk */
  strcpy(tag, chunk);
  strcat(tag, ".extent.t");
  chunk_size = mri_get_int(input, tag);


  /* Now we read everything in and get the output buffer.  Since the physio
     data sets should not be ridiculously huge, everything gets read in
     at the same time.  */ 

  input_buffer = mri_get_chunk(input, chunk, chunk_size, 0, MRI_FLOAT);

  output_space = (float *)malloc(2 * chunk_size * sizeof(float));
  output_buffer = (float **)malloc(chunk_size * sizeof(float *));
  for (i=0; i<chunk_size; i++) output_buffer[i]=output_space + (2 * i);

  /* Then we process the data, write out the results and close our
     files.  Then there is much happiness and joy. */

  process_data(input_buffer, output_buffer, chunk_size);
  mri_set_chunk(output, chunk, chunk_size, 0, MRI_FLOAT, output_space);
 
  mri_close_dataset(input);
  mri_close_dataset(output);

  exit(0);

}
