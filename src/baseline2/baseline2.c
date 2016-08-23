/************************************************************
 *                                                          *
 *  baseline2.c                                             *
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
 *  Original programming by Mark Fitzgerald                 *
 *  Hack by Todd Ogden                                      *
 ************************************************************/
/************************************************************
  DESCRIPTION OF BASELINE2.C
   baseline2.c reads and baseline-corrects right and left images for
   two-shot catch and hold data, combines them into a single image,
   averaging the area of overlap in the middle, and pads the ends with
   zeroes.
   baseline.m -fFile [-r even|odd|all|none][-jX][-oX][-sFile][-wFile]
   -r reverses even, odd, all, or no rows.
   -f reads input data from File
   -j performs a jitter correction with parameter X
   -o specifies that the amount of shift (half the amount of overlap) is
        X complex units (default is 10)
   -s writes output to File rather than baseline.mri
   -w writes baselines to File (default to basadj.par)
*************************************************************/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"

static char rcsid[] = "$Id: baseline2.c,v 1.12 2003/04/22 21:28:05 welling Exp $";

static char* progname;

static int parse_cl_reverse( const char* s, int* revOdd, int* revEven )
{
  if (!strcmp(s,"none")) {
    *revOdd= 0;
    *revEven= 0;
  }
  else if (!strcmp(s,"odd")) {
    *revOdd= 1;
    *revEven= 0;
  }
  else if (!strcmp(s,"even")) {
    *revOdd= 0;
    *revEven= 1;
  }
  else if (!strcmp(s,"all")) {
    *revOdd= 1;
    *revEven= 1;
  }
  else {
    return 0;
  }
  return 1;
}

static void parse_mrihdr_reverse( MRI_Dataset* ds, const char* chunk, 
				  int* revOdd, int* revEven )
{
  char buf[256];

  buf[sizeof(buf)-1]= '\0';
  snprintf(buf,sizeof(buf)-1,"%s.rowflip",chunk);
  if (mri_has(ds,buf) 
      && mri_get_int(ds,buf) != 0) {
    snprintf(buf,sizeof(buf)-1,"%s.rowflip_pattern",chunk);
    if (mri_has(ds,buf)) {
      char* pattern= mri_get_string(ds,buf);
      if (!strcmp(pattern,"none")) {
	*revOdd= 0;
	*revEven= 0;
      }
      else if (!strcmp(pattern,"odd")) {
	*revOdd= 1;
	*revEven= 0;
      }
      else if (!strcmp(pattern,"even")) {
	*revOdd= 0;
	*revEven= 1;
      }
      else if (!strcmp(pattern,"both")) {
	*revOdd= 1;
	*revEven= 1;
      }
      else Abort("%s: input dataset has invalid %s tag %s!\n",
		 progname,buf,pattern);
    }
    else Abort("%s: input dataset is missing %s info!\n",
	       progname,buf);
  }
  else {
    *revOdd= 0;
    *revEven= 0;
  }
}

int main(int argc, char**argv){
  FILE *parf=NULL;
  MRI_Dataset *Input=NULL,*Output=NULL;
  Filename filename, wfilename, parfilename;
  long i,j,k,m,s,t,d0,d1,d2,d3;
  int file_flag=0,jitter_flag=0,wrpar_flag=0,miss_flag,hover=10;
  int flip_odd= 0, flip_even= 0;
  double jitter,count,*fr0=NULL,*fr1=NULL,*fl0=NULL,*fl1=NULL;
  float *ndat=NULL, *fldat=NULL, *nldat=NULL, *frdat=NULL, *nrdat=NULL;
  long newdims[4];
  char* dimstr;
  char* missing_in= NULL;
  char* missing_out= NULL;
  char reverse_string[256];
  int reverse_set= 0;
  
  progname= argv[0];

  fprintf(stderr,"# %s\n",rcsid);

  /* Check to see if help was requested */
  if( ( argc > 1 ) && !strcmp( argv[1], "-help" ) )
    {
      if( argc == 2 )
        Help( "selecttopic" );
      else
        Help( (char*)(argv[2]) );
    }

  /*** Parse command line ***/

  cl_scan( argc, argv );

  if (cl_present( "f" ))
    Abort ("Option f has been replaced by an input output format.  Please see help file.\n");
  if (cl_present( "j" ))
    Abort ("Option j has been replaced by jitter|jit.  Please see help file.\n");
  if (cl_present( "o" ))
    Abort ("Option o has been replaced by shift|shi|s.  Please see help file.\n");
  if (cl_present( "s" ))
    Abort ("Option s has been replaced by an input output format.  Please see help file.\n");
  if (cl_present( "w" ))
    Abort ("Option w has been replaced by estimates|est|e.  Please see help file.\n");

  if (cl_get("jitter|jit","%option %lf", &jitter)) jitter_flag= 1;
  cl_get("shift|shi|s","%option %d", &hover);
  if (cl_get("estimates|est|e","%option %s[%]", "basadj.par", parfilename)) wrpar_flag= 1;
  reverse_set= cl_get( "reverse|rev|r", "%option %s", reverse_string );

  if(!cl_get("", "%s", filename)) {
    fprintf(stderr, "%s: Input file name not given.\n", argv[0]);
    exit(-1);
  }
  if(!cl_get("", "%s", wfilename)) {
    fprintf(stderr, "%s: Output file name not given.\n", argv[0]);
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

#ifdef JENN
  if(!file_flag){
    fprintf(stderr,"%s: Must provide input file.\n", argv[0]);
    Help("usage");
    exit(-1);
  }
#endif
  
  fprintf(stderr,"Overlap of %d\n",2*hover);

  if (!(Input= mri_open_dataset(filename, MRI_READ))) {
    Abort("Can't open file %s for reading.",filename);
  }
  /* If scan line reversal information is given in the command line,
   * implement it.  Otherwise, infer from the input file.
   */
  if (reverse_set) {
    if (!parse_cl_reverse( reverse_string, &(flip_odd), &(flip_even) ))
      Abort("%s: invalid -reverse opion <%s> given\n", 
	    progname,reverse_string);
  }
  else parse_mrihdr_reverse(Input, "images", &(flip_odd), &(flip_even));


  if (!(Output= mri_copy_dataset(wfilename, Input))) {
    Abort("Can't open file %s for writing.",wfilename);
  }
  hist_add_cl( Output, argc, argv );
  mri_set_int( Output, "images.rowflip", 0 );
  
  mri_create_chunk(Output, "missing");

  dimstr= mri_get_string(Input,"images.dimensions");
  if (strcmp(dimstr,"vxyzt"))
    Abort("Only dimension types vxyzt is supported.");

  if (mri_get_int(Input,"images.extent.v") != 2) 
    Abort("Vector Length Error (%d); should be 2.",
	  mri_get_int(Input,"images.extent.v"));

  d0 = mri_get_int(Input,"images.extent.x");
  d1 = mri_get_int(Input,"images.extent.y");
  d2 = mri_get_int(Input,"images.extent.z");
  d3 = mri_get_int(Input,"images.extent.t");

  d3 /= 2;  /*Reader thinks d3 is the total number of scans; it is actually
              half as many, since each image has a left and a right half-
              image*/

  /* Set the new dimensions prior to initializing Output */
  newdims[0] = d0*2;
  newdims[1] = d1;
  newdims[2] = d2; 
  newdims[3] = d3; /* which is half the original d3 . . . */

  mri_set_int(Output,"images.extent.x", newdims[0]);
  mri_set_int(Output,"images.extent.t", newdims[3]);
  mri_set_int(Output,"missing.extent.t", newdims[3]);
  mri_set_int(Output,"missing.extent.z", newdims[2]);
  mri_set_string(Output,"images.datatype", "float32");
  mri_set_string(Output,"missing.dimensions", "zt");
  mri_set_string(Output,"missing.datatype", "uint8");

  if(!(fr0 = (double *) malloc(sizeof(double)*d2*d3)))
    Abort("Memory allocation failed (fr0).");
  if(!(fr1 = (double *) malloc(sizeof(double)*d2*d3)))
    Abort("Memory allocation failed (fr1).");
  if(!(fl0 = (double *) malloc(sizeof(double)*d2*d3)))
    Abort("Memory allocation failed (fl0).");
  if(!(fl1 = (double *) malloc(sizeof(double)*d2*d3)))
    Abort("Memory allocation failed (fl1).");
  if(!(ndat = (float *) malloc(2*sizeof(float)*newdims[0]*newdims[1])))
    Abort("Memory allocation failed (ndat).");
  if(!(nrdat = (float *) malloc(2*sizeof(float)*d0*d1)))
    Abort("Memory allocation failed (nrdat).");
  if(!(nldat = (float *) malloc(2*sizeof(float)*d0*d1)))
    Abort("Memory allocation failed (nldat).");
  if (!(missing_out= (char*)malloc(sizeof(char)*newdims[2]*newdims[3])))
    Abort("Memory allocation failed (missing_out).");
  
  /* Main loop, over slices and time (pairs) */
  for (t=0; t<d3; t++) {
    long t1_offset= 2*t * (2*d0*d1*d2);
    long t2_offset= (2*t+1) * (2*d0*d1*d2);
    long out_offset= t * (2*newdims[0]*newdims[1]*newdims[2]);
    for (s=0; s<d2; s++) {
      
      /* Initialize ndat */
      for(k=0;k<2*newdims[0]*newdims[1];k++)
	ndat[k] = 0.0;
      
      /* Read an image pair */
      frdat= mri_get_image(Input, 2*t, s, MRI_COMPLEX_FLOAT);
      fldat= mri_get_image(Input, 2*t+1, s, MRI_COMPLEX_FLOAT);
      
      /* Compute baseline correction for each half-slice */
      fr0[t*d2+s] = fr1[t*d2+s] = fl0[t*d2+s] = fl1[t*d2+s] = count = 0.0;
      for(j=0;j<d1;j++) {
	if(j==(d1/4)) j += (d1/2);
	for (i=0;i<d0;i++) {
	  if(i==(d0/4)) i += (d0/2);
	  fr0[t*d2+s] += (double) frdat[2*(j*d0+i)];
	  fr1[t*d2+s] += (double) frdat[2*(j*d0+i)+1];
	  fl0[t*d2+s] += (double) fldat[2*(j*d0+i)];
	  fl1[t*d2+s] += (double) fldat[2*(j*d0+i)+1];
	  count += 1.0;
	}
      }
      fr0[t*d2+s] /= count;
      fr1[t*d2+s] /= count;
      fl0[t*d2+s] /= count;
      fl1[t*d2+s] /= count;
      
      /* Subtract off baseline, reversing lines as appropriate */
      for(j=0;j<d1;j++) {
	if(j%2) {
	  if (flip_odd) {
	    for (i=0;i<d0;i++) {
	      nrdat[2*(j*d0+d0-(i+1))] = frdat[2*(j*d0+i)] - fr0[t*d2+s];
	      nrdat[2*(j*d0+d0-(i+1))+1] = frdat[2*(j*d0+i)+1] - fr1[t*d2+s];
	      nldat[2*(j*d0+d0-(i+1))] = fldat[2*(j*d0+i)] - fl0[t*d2+s];
	      nldat[2*(j*d0+d0-(i+1))+1] = fldat[2*(j*d0+i)+1] - fl1[t*d2+s];
	    }
	  }
	  else {
	    for (i=0;i<d0;i++) {
	      nrdat[2*(j*d0+i)] = frdat[2*(j*d0+i)] - fr0[t*d2+s];
	      nrdat[2*(j*d0+i)+1] = frdat[2*(j*d0+i)+1] - fr1[t*d2+s];
	      nldat[2*(j*d0+i)] = fldat[2*(j*d0+i)] - fl0[t*d2+s];
	      nldat[2*(j*d0+i)+1] = fldat[2*(j*d0+i)+1] - fl1[t*d2+s];
	    }
	  }
	}
	else {
	  if (flip_even) {
	    for (i=0;i<d0;i++) {
	      nrdat[2*(j*d0+d0-(i+1))] = frdat[2*(j*d0+i)] - fr0[t*d2+s];
	      nrdat[2*(j*d0+d0-(i+1))+1] = frdat[2*(j*d0+i)+1] - fr1[t*d2+s];
	      nldat[2*(j*d0+d0-(i+1))] = fldat[2*(j*d0+i)] - fl0[t*d2+s];
	      nldat[2*(j*d0+d0-(i+1))+1] = fldat[2*(j*d0+i)+1] - fl1[t*d2+s];
	    }
	  }
	  else {
	    for (i=0;i<d0;i++) {
	      nrdat[2*(j*d0+i)] = frdat[2*(j*d0+i)] - fr0[t*d2+s];
	      nrdat[2*(j*d0+i)+1] = frdat[2*(j*d0+i)+1] - fr1[t*d2+s];
	      nldat[2*(j*d0+i)] = fldat[2*(j*d0+i)] - fl0[t*d2+s];
	      nldat[2*(j*d0+i)+1] = fldat[2*(j*d0+i)+1] - fl1[t*d2+s];
	    }
	  }
	}
      }
      
      /* Check to see if image is missing */
      miss_flag = 1;
      for(j=0;j<d1;j++){
	if(!miss_flag) break;
	for(i=0;i<d0;i++){
	  if(!miss_flag) break;
	  if((fldat[2*(j*d0+i)]!=-1.0) || (fldat[2*(j*d0+i)+1]!=-1.0))
	    miss_flag = 0;
	}
      }
      
      if (miss_flag) {
	for(j=0;j<newdims[1];j++){
	  for(i=0;i<newdims[0];i++){
	    ndat[2*(j*d0+i)] = ndat[2*(j*d0+i)+1] = -1.0;
	  }
	}
      } else {
	/* Put left image into ndat */
	for(j=0;j<d1;j++){
	  for(i=0;i<d0;i++){
	    ndat[2*((2*j*d0)+i+hover)] = nldat[2*(j*d0+i)];
	    ndat[2*((2*j*d0)+i+hover)+1] = nldat[2*(j*d0+i)+1];
	  }
	}
	/* Add right image to ndat */
	for(j=0;j<d1;j++){
	  for(i=0;i<d0;i++){
	    if (i<2*hover) {             /* Overlapped section -- average */
	      ndat[2*((2*j*d0)+d0+i-hover)] =
		0.5*(ndat[2*((2*j*d0)+d0+i-hover)]+nrdat[2*(j*d0+i)]);
	      ndat[2*((2*j*d0)+d0+i-hover)+1] = 
		0.5*(ndat[2*((2*j*d0)+d0+i-hover)+1]+nrdat[2*(j*d0+i)+1]);
	    }
	    else {
	      ndat[2*((2*j*d0)+d0+i-hover)] = nrdat[2*(j*d0+i)];
	      ndat[2*((2*j*d0)+d0+i-hover)+1] = nrdat[2*(j*d0+i)+1];
	    }
	  }
	}
	
	if(jitter_flag){
	  /* Place jitter correction here */
	}
      }
      mri_set_image(Output, t, s, MRI_COMPLEX_FLOAT, ndat);
    }
  }
    
  /* Transcribe "missing" dataset to output; 
   * non-zero means it's missing.
   */
  
  if (mri_has(Input, "missing")){
    missing_in= mri_get_chunk(Input,"missing",2*d2*d3,0,MRI_UNSIGNED_CHAR);
    for (s=0; s<d2; s++)
      for (t=0; t<d3; t++)
	missing_out[t*d2+s]= 
	  (missing_in[(2*t*d2)+s] || missing_in[(2*(t+1)*d2) + s]) ? 1 : 0;
    mri_set_chunk(Output,"missing",newdims[2]*newdims[3],0,
		  MRI_UNSIGNED_CHAR, missing_out);}
  else {
    for (s=0; s<d2; s++)
      for (t=0; t<d3; t++)
	missing_out[t*d2+s]= 0;
    mri_set_chunk(Output,"missing",newdims[2]*newdims[3],0,
		  MRI_UNSIGNED_CHAR, missing_out);}
	

  if(wrpar_flag){
    if(!(parf=fopen(parfilename,"w")))
      Abort("Couldn't open %s for writing.",parfilename);
    for(i=0;i<d3;i++){
      for(s=0;s<d2;s++){
	fprintf(parf,"%15.6lf %15.6lf %15.6lf %15.6lf\n",
           fl0[i*d2+s],fl1[i*d2+s],fr0[i*d2+s],fr1[i*d2+s]);
      }
    }
    fclose(parf);
  }

  free(fr0);
  free(fr1);
  free(fl0);
  free(fl1);
  free(ndat);
  free(nrdat);
  free(nldat);
  free(missing_out);
  mri_close_dataset(Input);
  mri_close_dataset(Output);
  fprintf(stderr,"#      Baseline adjustment complete.\n");
  return 0;
}

 
