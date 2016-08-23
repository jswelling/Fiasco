/************************************************************
 *                                                          *
 *  splus_binary_pipes.c                                    *
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
 *  Original programming by Joel Welling, 8/97              *
 ************************************************************/
#include <stdio.h>

static char rcsid[] = "$Id: splus_binary_pipes.c,v 1.3 1999/07/07 18:47:15 welling Exp $";

static int din; /* input pipe fdes */
static int dout; /* output pipe fdes */
static int initialized= 0;

int bio_set_pipe_fdes( int* din_fdes, int* dout_fdes )
{
  din= *din_fdes;
  dout= *dout_fdes;
  initialized= 1;
  return 0;
}

int bio_read_floats( float* buf, int* n )
{
  int bytesread;

  if (!initialized) {
    fprintf(stderr,"bio_read_floats: not initialized!\n");
    *n= 0;
  }
  else { 
#ifdef never
    bytesread= read(din, (void*)buf, *n * sizeof(float));
    if (bytesread==-1) {
      perror("bio_read_floats: error reading piped floats");
      *n= 0;
    }
    else *n= bytesread/sizeof(float);
#endif
    int i;
    for (i=0; i<*n; i++) buf[i]= 0.001*i;
  }
  return 0;
}

int bio_write_floats( float* buf, int* n )
{
  int bytesread;

  if (!initialized) {
    fprintf(stderr,"bio_write_floats: not initialized!\n");
    *n= 0;
  }
  else { 
    bytesread= write(dout, (void*)buf, *n * sizeof(float));
    if (bytesread==-1) {
      perror("bio_write_floats: error writing piped floats");
      *n= 0;
    }
    else *n= bytesread/sizeof(float);
  }
  return 0;
}


