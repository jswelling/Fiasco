/************************************************************
 *                                                          *
 *  smoother_tester.c                                       *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1998 Department of Statistics,         *
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
/* this module tests smoother */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"
#include "par.h"
#include "misc.h"
#include "acct.h"

void Help(char* foo) { /* fake entry point */ }

int main( int argc, char* argv[] ) {
  char pifname[512];
  char pofname[512];
  sm_type type;
  float band;
  float k;
  float thresh;
  int dt;
  int dz;
  int t;
  int z;
  float** data;
  float** smdata;
  Smoother* smoother;
  FILE* ifile;
  FILE* ofile;

  sm_init();

  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Get filenames */
  if (!cl_get( "pin", "%option %s", pifname)
      || !cl_get( "pout", "%option %s",pofname)
      || !cl_get( "t", "%option %d",&dt)
      || !cl_get( "z", "%option %d",&dz)) {
    fprintf(stderr,
	    "Usage: %s -t dt -z dz -pin file_in -pout file_out [...smoother opts...]\n",
	    argv[0]);
    exit(-1);
  }
  
  sm_parse_cl_opts();
  
  if (cl_cleanup_check()) {
    int i;
    fprintf(stderr,"%s: invalid argument in command line:\n    ",argv[0]);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(-1);
  }
  
  /*** End command-line parsing ***/
  
  sm_get_params( &type, &band, &k, &thresh, NULL );
  fprintf(stderr,"type %d, band %f, k %f, thresh %f\n",
	  (int)type, band, k, thresh);
  
  if (!(ifile= fopen(pifname,"r")))
    Abort("%s: can't open input param file!\n",argv[0]);
  if (!(ofile= fopen(pofname,"w")))
    Abort("%s: can't open output param file!\n",argv[0]);
  
  data= Matrix(dz, dt, float);
  smdata= Matrix(dz, dt, float);
  
  for (z=0; z<dz; z++) {
    for (t=0; t<dt; t++) {
      fscanf(ifile, "%f", &(data[z][t]));
    }
  }
  
  smoother= sm_create_smoother();
  for (z=0; z<dz; z++) {
    SM_SMOOTH( smoother, data[z], smdata[z], dt, NULL, z );
  }
  sm_destroy(smoother);
  
  for (z=0; z<dz; z++) {
    for (t=0; t<dt; t++) {
      fprintf(ofile, "%f\n", smdata[z][t]);
    }
  }

  fclose(ifile);
  fclose(ofile);
}
