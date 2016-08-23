/************************************************************
 *                                                          *
 *  windaq_reader.c                                         *
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "mri.h"
#if (SGI64 || SGI5 || SGIMP)
#include <bstring.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

#include "bio.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"
#include "smartreader.h"

/* windaq header files */
#include "windaq_header_info.h"

static char rcsid[] = "$Id: windaq_reader.c,v 1.9 2004/09/10 00:29:48 welling Exp $";

static void windaqRead( FileHandler* self, KVHash* info, 
			long long offset, long n, 
			SRDR_Datatype type, void* buf )
{
  int i;
  short* sbuf= (short*)buf;

  if (type != SRDR_INT16) 
    Abort("%s: windaqRead internal error: data type %s should be int16!\n",
	  progname,srdrTypeName[type]);

  if (kvGetInt(info,"datatype_in") != type)
    Abort("%s: datatype unexpectedly changed from %s to int16!\n",
	  progname, srdrTypeName[type]);

  /* read, then filter in place */
  baseRead( self, info, offset, n, SRDR_INT16, sbuf );
  for (i=0; i<n; i++) {
    if (sbuf[i]>>15 == 1) sbuf[i] = 0;
    else sbuf[i] = sbuf[i]>>4;
  }
}

static int scan_WINDAQ_data_header(KVHash*info, char* readfile)
{
  unsigned char header[WINDAQ_HEADER_SIZE_BYTES];
  FILE *fphead;
  int Windaq_Logo1;
  int veclen;
  unsigned long data_bytes;
  int ierror= 0;

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
  
  bio_big_endian_input = 
    (((Windaq_Logo1 = BRdInt8(header+WINDAQ_LOGO_OFF1)) == 1) ? 0 : 1);
  kvDefBoolean(info,"big_endian_input",bio_big_endian_input);
  if (debug) fprintf(stderr, "Header indicates %s input \n", 
	  bio_big_endian_input ? "bigendian" : "littleendian");

  veclen= (BRdUInt8(header+WINDAQ_CHANNELS_OFF) & 0x1f);
  kvDefString(info,"dimstr","vxyszt");
  kvDefInt(info,"dv",veclen);
  kvDefString(info,"description.v","gridded image-space");
  kvDefInt(info,"dx",1);
  kvDefString(info,"description.x","gridded image-space");
  kvDefInt(info,"dy",1);
  kvDefString(info,"description.y","gridded image-space");
  kvDefInt(info,"ds",1);
  kvDefString(info,"description.s","gridded image-space");
  kvDefInt(info,"dz",1);
  kvDefString(info,"description.z","gridded image-space");

  data_bytes=BRdInt32(header+WINDAQ_DATAACQ_SIZE_OFF);
  kvDefInt(info,"dt", ((data_bytes)/2)/(veclen));
  kvDefString(info,"description.t","gridded image-space");

  /* to get slope and intercepts for each channel in case we need to 
   * scale the raw data to match Windaq viewing program 
   */

  if (debug) {
    int i, channel[29];
    float slope1[29], intercept1[29];
    double slope2[29], intercept2[29];
    for (i=0; i<veclen; i++){
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

  kvDefInt(info,"datatype_in",SRDR_INT16);
  kvDefInt(info,"datatype_out",SRDR_INT16);
  kvDefLong(info,"skip",0);
  kvDefLong(info,"sliceskip",0);
  kvDefLong(info,"start_offset",WINDAQ_DATAACQ_OFF);

  return 1;
}

static void processHeader( FileHandler* self, KVHash* info, SList* cStack )
{
  /* Call the base class method */
  baseProcessHeader( self, info, cStack );

  if (!scan_WINDAQ_data_header(info, self->fileName))
    Abort("%s: error reading data header file.\n",progname);

  if (kvLookup(info,"phaseref")) {
    Warning(1,"Windaq file format does not allow phaseref!\n");
    kvDeleteAll(info,"phaseref");
  }
  if (kvLookup(info,"bandpassdir")) {
    Warning(1,"Windaq file format does not allow bandpass!\n");
    kvDeleteAll(info,"bandpassdir");
  }
  if (kvLookup(info,"rampfile")) {
    Warning(1,"Windaq file format does not allow ramp sample correction!\n");
    kvDeleteAll(info,"rampfile");
  }
  if (kvGetBoolean(info,"xchop")) {
    Warning(1,"Windaq file format does not allow xchop!\n");
    kvDefBoolean(info,"xchop",0);
  }
  if (kvGetBoolean(info,"ychop")) {
    Warning(1,"Windaq file format does not allow ychop!\n");
    kvDefBoolean(info,"ychop",0);
  }
  if (kvGetBoolean(info,"autoscale")) {
    Warning(1,"Windaq file format does not allow autoscale!\n");
    kvDefBoolean(info,"autoscale",0);
  }
  if (kvGetBoolean(info,"reorder")) {
    Warning(1,"Windaq file format does not allow reorder!\n");
    kvDefBoolean(info,"reorder",0);
  }

  if (debug) {
    fprintf(stderr, "number of Windaq channels = %ld \n", 
	    kvGetInt(info,"dv"));
  }
}

FileHandler* windaqFactory(char* fname, KVHash* info)
{
  FileHandler* result= baseFactory(fname);
  result->typeName= strdup( nameof_filetype(FILE_WINDAQ) );
  result->processHeader= processHeader;
  result->read= windaqRead;
  return result;
}

int windaqTester(const char* filename)
{
  return( check_filetype(filename)==FILE_WINDAQ );
}

