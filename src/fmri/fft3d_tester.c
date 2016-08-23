/************************************************************
 *                                                          *
 *  fft3d.c                                                 *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *                                                          *
 *  Copyright (c) 1999 Department of Statistics             *
 *                     Carnegie Mellon University           *
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
 *  Original programming by Joel Welling, 9/98              *
 ************************************************************/

/* This utility reads a file containing Pgh MRI data in "xyz" form,
 * and does various FFT things to it to test "fft3d".  It then
 * writes the resulting data out.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"

int main(int argc, char* argv[]) 
{
  MRI_Dataset *Input= NULL;
  MRI_Dataset *Output= NULL;
  float* obuf;
  float* ibuf;
  char* ops;
  char* tok;
  FComplex* data;
  int nx;
  int ny;
  int nz;
  int i;
  int j; 
  int k;
  int firstOp;

  if (argc != 4) {
    fprintf(stderr,"Usage: %s ops in.mri out.mri\n", argv[0]);
    fprintf(stderr,
	    "   in.mri must be an xyz Pgh MRI dataset in image space,\n");
    fprintf(stderr,
	    "   containing float data.  ops is a comma-separated series\n");
    fprintf(stderr,
	    "   of FFT operations like x,xy,-x,-xy\n");
  }
  else {
    Input= mri_open_dataset( argv[2], MRI_READ );
    if (!Input) {
      fprintf(stderr,"%s: cannot open %s for reading!\n",argv[0],argv[2]);
      exit(-1);
    }
    if (strcmp(mri_get_string(Input,"images.dimensions"),"xyz")) {
      fprintf(stderr,"%s: input data must have dims xyz!\n",argv[0]);
      exit(-1);
    }
    if (strcmp(mri_get_string(Input,"images.datatype"),"float32")) {
      fprintf(stderr,"%s: input data must have type float32!\n",argv[0]);
      exit(-1);
    }

    ops= strdup(argv[1]);

    nx= mri_get_int(Input,"images.extent.x");
    ny= mri_get_int(Input,"images.extent.y");
    nz= mri_get_int(Input,"images.extent.z");

    Output= mri_copy_dataset( argv[3], Input );
    hist_add_cl( Output, argc, argv );
    if (!Output) {
      fprintf(stderr,"%s: error creating %s for writing!\n",argv[0],argv[2]);
      exit(-1);
    }

    if (!(obuf= (float*)malloc(nx*ny*nz*sizeof(float)))) {
      fprintf(stderr,"Unable to allocate %d floats!\n",nx*ny*nz);
      exit(-1);
    }
    if (!(data= (FComplex*)malloc(nx*ny*nz*sizeof(FComplex)))) {
      fprintf(stderr,"Unable to allocate %d FComplex!\n",nx*ny*nz);
      exit(-1);
    }

    ibuf= mri_get_chunk(Input,"images",nx*ny*nz,0,MRI_FLOAT);

    for (i=0; i<nx; i++) {
      for (j=0; j<ny; j++) {
	for (k=0; k<nz; k++) {
	  data[(((i*ny)+j)*nz)+k].real= ibuf[(((k*ny)+j)*nx)+i];
	  data[(((i*ny)+j)*nz)+k].imag= 0.0;
	}
      }
    }

    firstOp= 1;
    while (tok=strtok((firstOp ? ops : NULL),",")) {
      fprintf(stderr,"<%s>\n",tok);
      if (!strcmp(tok,"xy")) fft3d(data,nx,ny,nz,+1,"xy");
      else if (!strcmp(tok,"-xy")) fft3d(data,nx,ny,nz,-1,"xy");
      else if (!strcmp(tok,"yz")) fft3d(data,nx,ny,nz,+1,"yz");
      else if (!strcmp(tok,"-yz")) fft3d(data,nx,ny,nz,-1,"yz");
      else if (!strcmp(tok,"xz")) fft3d(data,nx,ny,nz,+1,"xz");
      else if (!strcmp(tok,"-xz")) fft3d(data,nx,ny,nz,-1,"xz");
      else if (!strcmp(tok,"x")) fft3d(data,nx,ny,nz,+1,"x");
      else if (!strcmp(tok,"-x")) fft3d(data,nx,ny,nz,-1,"x");
      else if (!strcmp(tok,"y")) fft3d(data,nx,ny,nz,+1,"y");
      else if (!strcmp(tok,"-y")) fft3d(data,nx,ny,nz,-1,"y");
      else if (!strcmp(tok,"z")) fft3d(data,nx,ny,nz,+1,"z");
      else if (!strcmp(tok,"-z")) fft3d(data,nx,ny,nz,-1,"z");
      else if (!strcmp(tok,"xyz")) fft3d(data,nx,ny,nz,+1,"xyz");
      else if (!strcmp(tok,"-xyz")) fft3d(data,nx,ny,nz,-1,"xyz");
      else {
	fprintf(stderr,"Unrecognized token <%s>!\n",tok);
	exit(-1);
      }
      firstOp= 0;
    }

    for (i=0; i<nx; i++) {
      for (j=0; j<ny; j++) {
	for (k=0; k<nz; k++) {
	  obuf[(((k*ny)+j)*nx)+i]= data[(((i*ny)+j)*nz)+k].real;
	}
      }
    }

    mri_set_chunk(Output,"images",nx*ny*nz,0,MRI_FLOAT,obuf);

    mri_close_dataset(Input);
    mri_close_dataset(Output);
  }

  return 0;
}
