/************************************************************
 *                                                          *
 *  cardiopeaks.c                                           *
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
#include "math.h"

/* Empirically derived constants */
#define DEFAULT_WINDOW 20
#define DEFAULT_THRESHOLD 25.0
#define KEYBUF_SIZE 512
#define DEFAULT_TIMEBASE 0

/* Global parameters */
int window;
float threshold;
float timebase;
int magic_number;

/* the magic number is used to get rid of extra peaks.  it's
   smaller for Becky's data than for Windaq data. */

int find_next_peak(float *, int, int);
void process_data(float *, float *, int);


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
/* does all the actual data processing taking the cardio data in the 	*/
/* input buffer and making our cooresponding ramping function in the 	*/
/* output buffer where n is the size of the buffers.			*/

void process_data(float *in_buffer, float *out_buffer, int n) {
	float coeffs[3], slope;
	int i, i2, j, k, w2 = window/2, rate_window;
	Smoother* in_smoother;

	float *temp_buffer = (float *)malloc(n * sizeof(float));
	float *in_buffer_smoothed = (float *)malloc(n * sizeof(float));
	int *rate_buffer = (int *)malloc(n * sizeof(int));

	/* Set the temp and output buffer contents to 0 */
	for(i = 0; i < n; i++) {
		out_buffer[i] = 0;
		temp_buffer[i] = 0;
		in_buffer_smoothed[i] = 0;
		rate_buffer[i] = 0;
	}

	/* Smooth the in_buffer */

	in_smoother = sm_create_smoother();
	SM_SMOOTH( in_smoother, in_buffer, in_buffer_smoothed,
		   n, NULL, 0);

	/* Make a first pass over the data fitting parabolas */
	for(i = w2; i < n - w2; i++) {
		polyfitlin(in_buffer_smoothed + i - w2, window, coeffs, 3);
		temp_buffer[i] = coeffs[2] * 100;
	}

	/* Now go through and find the peaks from the fitting stage */
	/* at each peak, we seek a local minimum in the original    */
	/* data around it.					    */

	for(i = w2; i < n - w2; i++) 

		if(temp_buffer[i-1] < temp_buffer[i] && 
		   temp_buffer[i+1] < temp_buffer[i] &&
		   temp_buffer[i] > threshold) {
			j = i;

			while(((in_buffer_smoothed[j] > in_buffer_smoothed[j+1])
			       ||
			       (in_buffer_smoothed[j] > in_buffer_smoothed[j-1]))
			      &&
			      (j > 1 && j < n-1))
			  if(in_buffer_smoothed[j+1] < in_buffer_smoothed[j]) 
			    j = j+1;
			  else j = j-1;

			out_buffer[j] = 1;
		}


	/* Calculate and display rate */
	if (timebase > 0) {
	  for (i = 2; i < n; i++) {
	    if (out_buffer[i] == 1) rate_buffer[i]= rate_buffer[i-1]+1;
	    else rate_buffer[i] = rate_buffer[i-1];}
	  rate_window = rint(60/timebase);
	  for (i = rate_window; i < n - rate_window; i+=rate_window)
	    fprintf(stderr, "Heart rate for minute %d is %d. \n", i/rate_window, rate_buffer[i]-rate_buffer[i-rate_window]);
	}
	    


	/* Now linearize everything */

	i2 = -1;
	while((i = find_next_peak(out_buffer, i2, n)) != -1) {
	        i2 = i;
		j = find_next_peak(out_buffer, i, n);
		if(j == -1) break;
		/* This is some magic to take care of extra peaks; i2 makes
		 * sure that i moves to the right place if there is an extra
		 * peak
		 */
		if(((in_buffer[j] > in_buffer[find_next_peak(out_buffer,j,n)]) 
		    && ((find_next_peak(out_buffer, j, n) - j) < magic_number))
		   || (j-i < magic_number)) {
		  i2 = j;
		  j = find_next_peak(out_buffer, j, n);
		}
		slope = 1.0/(j - i - 1);
		for(k = i + 1; k < j; k++) 
			out_buffer[k] = (k - i - 1) * slope;
	}

	/* Now we have to take care of the sections at the beginning and
	   end of the data that just didn't fit into our window for the
	   we copy the slope from the first and last sections where we were
	   able to identify peaks, and fill in the gaps.  */


	i = find_next_peak(out_buffer, -1, n);
	j = find_next_peak(out_buffer, i, n);
	for(k = 0; k < i; k++) out_buffer[k] = out_buffer[j - i + k];

	for(i = n - 1; i >= 0 && out_buffer[i] == 0; i--);
	for(j = i - 1; j >= 0 && out_buffer[j] != 1; j--);
	for(k = i + 1; k < n; k++) out_buffer[k] = out_buffer[j - i + k];

}

/* find_next_peak:						*/
/* given the location of the previous peak, finds the next 	*/
/* peak.  if none is to be found, returns -1			*/

int find_next_peak(float *buffer, int start, int size) {
  int i = start + 1;
  while(i < size && buffer[i] != 1) i++;
  if(i == size) return -1;
  return i;
}

int main(int argc, char *argv[]) {
  MRI_Dataset *input, *output;
  int chunk_size, num_elements, i, j, k;
  float *input_buffer, *output_buffer;
  char chunk[512];
  char infile[512];
  char outfile[512];
  char tag[600];
  char t[512]="t";
  char *dimstr;
  int dim_extent;

  sm_init();

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

  cl_get("chunk|chu|c", "%option %s[%]", "physio", chunk);
  cl_get("window|win", "%option %d[%]", DEFAULT_WINDOW, &window);
  cl_get("threshold|thr", "%option %f[%]", DEFAULT_THRESHOLD, &threshold); 
  cl_get("timebase|tim", "%option %f[%]", DEFAULT_TIMEBASE, &timebase); 

  sm_set_params( SM_TRIANGULAR, 20.0, 0, 0, NULL );

  sm_parse_cl_opts();

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

  if (chunk=="images") magic_number = 200;
  else magic_number = 20;

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

  /* Get the length of the chunk */
  strcpy(tag, chunk);
  strcat(tag, ".extent.t");
  chunk_size = mri_get_int(input, tag);

  /* Now we read everything in and get the output buffer.  Since the physio
     data sets should not be ridiculously huge, everything gets read in
     at the same time.  */ 

  input_buffer = mri_get_chunk(input, chunk, chunk_size, 0, MRI_FLOAT);
  output_buffer = mri_get_chunk(output, chunk, chunk_size, 0, MRI_FLOAT);

  /* Then we process the data, write out the results and close our
     files.  Then there is much happiness and joy. */

  process_data(input_buffer, output_buffer, chunk_size);

  mri_set_chunk(output, chunk, chunk_size, 0, MRI_FLOAT, output_buffer);
 
  mri_close_dataset(input);
  mri_close_dataset(output);

  exit(0);
}
