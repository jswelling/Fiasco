/************************************************************
 *                                                          *
 *  desmith_reader.c                                         *
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
#include <strings.h>
#include <ctype.h>
#include "mri.h"
#if (SGI64 || SGI5 || SGIMP)
#include <bstring.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <assert.h>

/* Stupid DARWIN! */
#ifdef DARWIN
char* strsep( char** stringp, const char* delim );
#endif

#include "bio.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"
#include "smartreader.h"

#define DESMITH_HEADER_SIZE_BYTES 758
#define MAX_CHANNELS 30
#define MAX_FIELDS 30
#define INITIAL_N_FRAME_LOCS 100
#define TRIGGERED 1
#define UNTRIGGERED 2

static char rcsid[] = "$Id: desmith_reader.c,v 1.9 2009/08/11 05:33:54 welling Exp $";

/* Something to store correct offset and other slice info in */
typedef struct DESmith_data_struct {
  long long offset;
  int slices;
  int slice_count;
  long long *frame_locs;
  long long *timebuf;
  FileHandler* timestampHandler;
  KVHash* timestampInfo;
  long n_frame_locs;
  long long dt;
  int dv;
  int trigger_type;
} DESmithData;

static char* strip( char* s )
{
  char* tail;

  if (!s) return NULL;
  while (isspace(*s)) s++;
  if (*s=='\0') return s;
  tail= s+(strlen(s)-1);
  while (isspace(*tail)) tail--;
  *(tail+1)='\0';
  return s;
}

static void DESmithDestroySelf( FileHandler* self )
{
  DESmithData* data= (DESmithData*)(self->hook);
  free(data->frame_locs);
  /* We do not free timestampHandler or timestampInfo because they
   * are on the caller's handler pair list, and will be freed 
   * when the caller finishes processing them.
   */
  baseDestroySelf(self);
}

static void desmithRead( FileHandler* self, KVHash* info, 
			long long offset, long n, 
			SRDR_Datatype type, void* buf )
{

  DESmithData* data= (DESmithData*)(self->hook);
  int i;
  short* sbuf= (short*)buf;
  long long frameBlocks;
  long long frameLength;
  long long frameStart;
  long long blocks_this_frame= 0;
  long long new_bytes_chomped=0;
  long long starting_bytes_chomped=0;
  long long starting_blocks_chomped= 0;
  long long new_blocks_chomped= 0;
  long long buf_offset=0;
  long long header_offset=data->frame_locs[0];
  int inner_skip_z=34;
  int inner_skip_t=8;
  long timebuf_offset;
  long long outer_start_offset= kvGetLong(info,"start_offset");
  /* outer skip.t and skip.v are presumed to be zero */
  

  if (type != SRDR_INT16) 
    Abort("%s: desmithRead internal error: data type %s should be int16!\n",
	  progname,srdrTypeName[type]);

/*   if (kvGetInt(info,"datatype_in") != type) */
/*     Abort("%s: datatype unexpectedly changed from %s to int16!\n", */
/* 	  progname, srdrTypeName[type]); */

  switch (data->trigger_type) {
  case UNTRIGGERED:
    baseRead(self, info, offset, n, type, sbuf);
    break;

  case TRIGGERED:

    assert(n%data->dv==0); /* assumption is that data comes dv at a time */

    /* Work out the slice the calling routine *thinks* we're dealing with */
    data->slice_count= 
      1+(int)((offset-outer_start_offset)/(srdrTypeSize[SRDR_INT16]*data->dv*data->dt));

    /* Where does it think we are in the block? */
    starting_bytes_chomped= 
      offset - (outer_start_offset
		+ ((data->slice_count - 1)
		   *(srdrTypeSize[SRDR_INT16]*data->dv*data->dt)));
    assert(starting_bytes_chomped % (srdrTypeSize[SRDR_INT16]*data->dv) == 0);
    starting_blocks_chomped= starting_bytes_chomped/(srdrTypeSize[SRDR_INT16]*data->dv);
    new_bytes_chomped= 0;
    frameLength=
      (data->frame_locs[data->slice_count+1]-data->frame_locs[data->slice_count])
      -(inner_skip_z + inner_skip_t);
    frameBlocks= frameLength/(data->dv*srdrTypeSize[SRDR_INT16]);
    frameStart= data->frame_locs[data->slice_count] + inner_skip_z;
    blocks_this_frame= 
      (frameBlocks*data->dv>n)? n/data->dv : frameBlocks;
    timebuf_offset= 
      (data->slice_count-1)*data->dt + starting_blocks_chomped;
#ifdef never
    fprintf(stderr,"n= %ld\n",n);
    fprintf(stderr,"frameLength= %lld\n",frameLength);
    fprintf(stderr,"data->dv= %d\n",data->dv);
    fprintf(stderr,"frameBlocks= %lld\n",frameBlocks);
    fprintf(stderr,"frameLength= %lld\n",frameLength);
    fprintf(stderr,"starting_bytes_chomped= %lld\n",starting_bytes_chomped);
    fprintf(stderr,"starting_blocks_chomped= %lld\n",starting_blocks_chomped);
    fprintf(stderr,"blocks_this_frame= %lld\n",blocks_this_frame);
    fprintf(stderr,"timebuf_offset= %ld\n",timebuf_offset);
#endif
    new_bytes_chomped= 0;
    new_blocks_chomped= 0;
    buf_offset= 0;
    while (new_blocks_chomped<blocks_this_frame){

      long long rd_offset= 
	header_offset+data->frame_locs[data->slice_count]+inner_skip_z
	+starting_bytes_chomped+new_bytes_chomped;
      long long timestamp_offset= rd_offset+data->dv*srdrTypeSize[SRDR_INT16];

      baseRead(self, info, rd_offset, data->dv, SRDR_INT16, sbuf+(buf_offset));
      new_bytes_chomped+=data->dv*srdrTypeSize[SRDR_INT16];
      new_blocks_chomped += 1;
      buf_offset+=data->dv;

      /* After the following read, the new timestamp resides in memory
       * at the location *(data->timebuf + timebuf_offset)
       */
      if (timestamp_offset+srdrTypeSize[SRDR_INT64]<=self->totalLengthBytes)
	baseRead(self, data->timestampInfo, 
		 timestamp_offset, 1, SRDR_INT64, 
		 data->timebuf+timebuf_offset);
      else *(data->timebuf + timebuf_offset)= 0;
      timebuf_offset += 1;
    }
    
    /* pad end of buffer with zeros (since frames are different frameLengths) */
    if (blocks_this_frame*data->dv < n) {
      for (i=blocks_this_frame*data->dv;i<n;i++) sbuf[i]=0;
      for (i=blocks_this_frame; i<n/data->dv; i++) {
	*(data->timebuf + timebuf_offset)= 0;
	timebuf_offset += 1;
      }
    }
    
  }
}

static void checkFrameLocsSize(DESmithData* data, int num_slices)
{
  if (num_slices>=data->n_frame_locs) {
    data->frame_locs= realloc(data->frame_locs,
			      2*data->n_frame_locs*sizeof(long long));
    if (data->frame_locs==NULL)
      Abort("%s: unable to realloc %d bytes!\n",progname,
	    2*data->n_frame_locs*sizeof(long long));
    data->n_frame_locs *= 2; 
 }
}

static int scan_DESMITH_data_header(FileHandler* self, KVHash*info, char* readfile)
{
  char header[DESMITH_HEADER_SIZE_BYTES];
  FILE *fphead;
  int veclen;
  int ierror= 0;
  int i;
  char* header_start, *header_end;
  char* field_strings[MAX_FIELDS];
  char* frame_find, frame_find2;
  int frame_start, frame_start2;
  int dsmith_ftype;
  char** fp;
  char* keys[MAX_CHANNELS][2];
  char* channel_names[MAX_CHANNELS];
  char** chp;
  int channel_field=0;
  char channel_headings[MAX_CHANNELS][512];
  int channels=1;
  DESmithData* data= (DESmithData*)(self->hook); 
  int num_slices=0;
  long long count=0;
  int nseen=0;
  int c;
  int num_fields;
  long long max_frame_length=0;

  if ((fphead = fopen(readfile,"r"))!=NULL)
    {
      if (fseek(fphead, (long) 0, SEEK_SET)) {
	perror("Error seeking header");
	ierror=1;
      }
      else {
	if (fread(header, sizeof(char), DESMITH_HEADER_SIZE_BYTES, fphead)
	    != DESMITH_HEADER_SIZE_BYTES) {
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

  for (i=0; i<DESMITH_HEADER_SIZE_BYTES; i++){ 
    if (header[i]=='\r') header[i]='\n';
  }

  header_start=strstr(header,"<START HEADER>\n");

  /* label data as triggered or untriggered */
  if ((frame_find=strstr(header,"!!FRAME#"))==NULL){
    data->trigger_type=UNTRIGGERED;
/*     Abort("Jenn needs to fix this\n"); */
  }
  else {
    data->trigger_type=TRIGGERED;
  }

  if (debug) fprintf(stderr,"This is file type %d \n",data->trigger_type);

  header_end=strstr(header,"<END HEADER>\n");
  header_end+=strlen("<END HEADER>\n");
  data->frame_locs[0]=(int)(header_end-header);
  if (debug) {
    fprintf(stderr, "Slice #%d marker is at offset %lld.\n",
	    0, data->frame_locs[0]);
  }

  if ((fphead = fopen(readfile,"r"))!=NULL)
    {
      if (fseek(fphead, (long) data->frame_locs[0], SEEK_SET)) {
	perror("Error finding first frame");
	ierror=1;
      }
      else {
	while (!feof(fphead)) {
	  c=fgetc(fphead);
	  switch (nseen) {
	  case 0:
	    if (c=='!') nseen=1;
	    else nseen=0;
	    break;
	  case 1:
	    if (c=='!') nseen=2;
	    else nseen=0;
	    break;
	  case 2:
	    if (c=='F') nseen=3;
	    else nseen=0;
	    break;
	  case 3:
	    if (c=='R') nseen=4;
	    else nseen=0;
	    break;
	  case 4:
	    if (c=='A') nseen=5;
	    else nseen=0;
	    break;
	  case 5:
	    if (c=='M') nseen=6;
	    else nseen=0;
	    break;
	  case 6:
	    if (c=='E') nseen=7;
	    else nseen=0;
	    break;
	  case 7:
	    if (c=='#') {
	      num_slices++;
	      checkFrameLocsSize(data,num_slices);
	      data->frame_locs[num_slices]=count+1;
	      if (debug) {
		fprintf(stderr, "Slice #%d starts at offset %lld.\n",num_slices, count+1);}
	      nseen=0;
	    break;
	    }
	  }
	  count++;
	}      
	/*JENN1*/
	num_slices++;
	checkFrameLocsSize(data,num_slices);
	data->frame_locs[num_slices]=count;
	if (debug) {
	  fprintf(stderr, 
		  "Slice trailer #%d starts at offset %lld.\n",
		  num_slices, count+1);
	}
    
	if (fclose(fphead)) {
	  perror("Error scanning for frames");
	  ierror=1;
	}
      }
    }
  else {
    perror("Error opening file while scanning for frames");
    ierror= 1;
  }

  /*JENN1a*/
  if (!ierror) {
    switch (data->trigger_type) {
    case TRIGGERED:
      data->slices=num_slices-1;
      break;
    case UNTRIGGERED:
      data->slices=num_slices;
    }
  }      


  /* Separate header into parts consisting of field name and value */
  for (fp=field_strings; (*fp=strsep(&header_start, "\n")) != NULL;) {
    if (**fp != '\0')
      if (++fp >= &field_strings[30])
	break;
  }

  /* Separate field_strings into name (keys[][0]) and value (keys [][1]) pairs */
  i=1;
  while ((field_strings[i] != NULL) && (field_strings[i][0] != '<')) {
    keys[i][0]=strsep(&field_strings[i],"=");
    keys[i][1]=field_strings[i];
    i++;
  }
  num_fields = i;

  if (debug) {fprintf(stderr,"There are %d fields\n",num_fields);}

  /* For now, assuming [CHANNEL] is first field after header. */

  i=1;
/*   while (strip(keys[i][0]) != "[CHANNELS]") i++; */
/*   fprintf(stderr, "Channels are at field %d\n",i); */

/*   while ((keys[i][0]) != NULL){ */
/*     fprintf(stderr, "%s\n",keys[i][0]); */
/*   } */

/*   for (i=0; i<num_fields; i++){ */
/*     if (strstr(keys[i][0],"[CHANNELS]")!=NULL) { */
/*       channel_field=i; */
/*       break; */
/*     } */
/*     fprintf(stderr,"%s\n",keys[i][0]); */
/*   } */

  /* Separate the channel field into different channels */
  for (chp=channel_names; (*chp=strsep(&keys[1][1], ",")) != NULL;) {
    if (**chp != '\0')
      if (++chp >= &channel_names[MAX_FIELDS])
	break;
  }

  while(channel_names[channels] != NULL){
    channels++;
  }

  for (i=1; i<channels; i++) {
    sprintf(channel_headings[i],"label.v.%d",i-1);
  }
  
  for (i=1; i<channels; i++) {
    kvDefString(info,channel_headings[i],channel_names[i]);
    if (debug) {
      fprintf(stderr,"i=%d, channel_heading=%s, channel_name=%s\n",i,channel_headings[i],channel_names[i]);
    }
  }

  for (i=2; i<num_fields; i++) {
    kvDefString(info,strip(keys[i][0]),strip(keys[i][1]));
    if (debug) {
      fprintf(stderr,"key %2d  %20s  %10s\n",i,keys[i][0],keys[i][1]);
    }
  }
  veclen=channels-1;
  kvDefString(info,"dimstr","vtz");
  kvDefInt(info,"dv",veclen);
  kvDefString(info,"description.v","channels");
  data->dv=veclen;

  kvDefInt(info,"dz",data->slices);
  kvDefString(info,"description.z","gridded image-space");

  /* need to get dt from largest value between frame_starts */
  switch (data->trigger_type) {
  case TRIGGERED:
    for (i=1; i<=data->slices; i++) {
      if (data->frame_locs[i+1]-data->frame_locs[i]>max_frame_length) {
	max_frame_length=data->frame_locs[i+1]-data->frame_locs[i];
      }
    }
    break;
  case UNTRIGGERED:
    max_frame_length=data->frame_locs[data->slices];
  }
  if (debug) fprintf(stderr,"max_frame_length= %lld\n",max_frame_length);

  /* for TRIGGERED data, 42 because 
     "each scan of the microscope the data begins with the string "!!FRAME#" <8 bytes> 
     followed by an unsigned 16-bit integer <2 bytes> 
       indicating the frame number (Starting with 0) 
     followed by two 128-bit timestamps <16 bytes * 2>"

     veclen+4 because "(each of the data channels followed by 2 32-bit unsigned 
     integer counters...)(2 bytes * 2) " */

  /* 

  /* for UNTRIGGERED data, 
     "the file only contains one continuous chunk of data, with no header 
        before this other than the file header.  
     The data is signed 16-bit integers, one for each channel indicated in the header; 
     however, the additional 64-bits for the timestamps is not included 
       in a rolling mode file."

     Jenn's note:  the bits for the timestamps are still there, just filled with zeros
  */

  switch (data->trigger_type) {
  case TRIGGERED:
    data->dt=(max_frame_length-42)/((veclen+4)*2);
    kvDefLong(info,"skip.t",0);
    break;
  case UNTRIGGERED:
    data->dt=(max_frame_length)/((veclen+4)*2);
    kvDefLong(info,"skip.v",8); 
  }

  kvDefInt(info,"dt",data->dt);

  kvDefString(info,"description.t","gridded image-space");

  kvDefInt(info,"datatype_in",SRDR_INT16);
  kvDefInt(info,"datatype_out",SRDR_INT16);
  kvDefLong(info,"skip.t",0); 
  kvDefLong(info,"sliceskip",0); 
  kvDefLong(info,"start_offset",data->frame_locs[0]);
  data->slice_count=0;

  kvDefBoolean(info,"big_endian_input",1);

  return 1;
}

static void addTimestampChunk( FileHandler* self, KVHash* info, SList* cStack )
{
  /* Add chunks for bandpass directory information */
  DESmithData* data= (DESmithData*)(self->hook); 
  KVHash* defs= kvGetHash(info,"definitions");
  ChunkHandlerPair* chPair= NULL;
  KVHash* subInfo= kvFactory(KV_DEFAULT_SIZE);
  long long* timebuf= NULL;
  
  if (!(timebuf=(long long*)malloc(data->slices*data->dt*sizeof(long long)))) 
    Abort("%s: unable to allocate %d bytes!\n",
	  data->slices*data->dt*sizeof(long long));
  data->timebuf= timebuf;
  
  initInfoHash(subInfo);
  kvDefString(subInfo,"chunkname","timestamp");
  kvDefString(subInfo,"chunkfile",".time");
  kvDefString(subInfo,"dimstr","tz");
  kvDefInt(subInfo,"dt",data->dt);
  kvDefInt(subInfo,"dz",data->slices);
  kvDefLong(subInfo,"start_offset",0);
  kvDefInt(subInfo,"datatype_in",SRDR_INT64);
  kvDefInt(subInfo,"handler_datatype_out",SRDR_INT64);
  
  if (!(chPair= (ChunkHandlerPair*)malloc(sizeof(ChunkHandlerPair))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,sizeof(ChunkHandlerPair));
  chPair->info= subInfo;
  chPair->handler= 
    ramDataHandlerFactory(timebuf, data->dt*data->slices, SRDR_INT64);
  data->timestampHandler= chPair->handler;
  data->timestampInfo= subInfo;
  slist_push(cStack,chPair);
}

static void processHeader( FileHandler* self, KVHash* info, SList* cStack )
{
  DESmithData* data= (DESmithData*)(self->hook); 

  /* Call the base class method */
  baseProcessHeader( self, info, cStack );

  if (!scan_DESMITH_data_header(self, info, self->fileName)) 
    Abort("%s: error reading data header file.\n",progname);

  if (data->trigger_type==TRIGGERED)
    addTimestampChunk(self, info, cStack);

  /* Leaving stuff commented out until Joel looks at file */

/*   if (kvLookup(info,"phaseref")) { */
/*     Warning(1,"Desmith file format does not allow phaseref!\n"); */
/*     kvDeleteAll(info,"phaseref"); */
/*   } */
/*   if (kvLookup(info,"bandpassdir")) { */
/*     Warning(1,"Desmith file format does not allow bandpass!\n"); */
/*     kvDeleteAll(info,"bandpassdir"); */
/*   } */
/*   if (kvLookup(info,"rampfile")) { */
/*     Warning(1,"Desmith file format does not allow ramp sample correction!\n"); */
/*     kvDeleteAll(info,"rampfile"); */
/*   } */
/*   if (kvGetBoolean(info,"xchop")) { */
/*     Warning(1,"Desmith file format does not allow xchop!\n"); */
/*     kvDefBoolean(info,"xchop",0); */
/*   } */
/*   if (kvGetBoolean(info,"ychop")) { */
/*     Warning(1,"Desmith file format does not allow ychop!\n"); */
/*     kvDefBoolean(info,"ychop",0); */
/*   } */
/*   if (kvGetBoolean(info,"autoscale")) { */
/*     Warning(1,"Desmith file format does not allow autoscale!\n"); */
/*     kvDefBoolean(info,"autoscale",0); */
/*   } */
/*   if (kvGetBoolean(info,"reorder")) { */
/*     Warning(1,"Desmith file format does not allow reorder!\n"); */
/*     kvDefBoolean(info,"reorder",0); */
/*   } */

/*   if (debug) { */
/*     fprintf(stderr, "number of Desmith channels = %ld \n",  */
/* 	    kvGetInt(info,"dv")); */

/*   } */
}

FileHandler* desmithFactory(char* fname, KVHash* info)
{
  FileHandler* result= baseFactory(fname);
  DESmithData* data;

  result->typeName=strdup("desmith");
  result->processHeader= processHeader;
  result->read= desmithRead;
  result->destroySelf=DESmithDestroySelf;

  if (!(data= (DESmithData*)malloc(sizeof(DESmithData))))
    Abort("%s: unable to allocate %d bytes!\n",progname,sizeof(DESmithData));

  if (!(data->frame_locs=(long long*)malloc(INITIAL_N_FRAME_LOCS
					    *sizeof(long long))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,INITIAL_N_FRAME_LOCS*sizeof(long long));
  data->n_frame_locs= INITIAL_N_FRAME_LOCS;

  data->offset=0;
  data->timebuf= NULL;
  data->timestampHandler= NULL;
  data->timestampInfo= NULL;

  result->hook= data;

  return result;
}

int desmithTester(const char* filename)
{ char buf[256];
  FILE* fphead= NULL;
  int ierror= 0;

  if ((fphead = fopen(filename,"r"))!=NULL) {
    if (fgets(buf,sizeof(buf),fphead)!=NULL) {
      if (strncasecmp(buf,"<START HEADER>",14)) {
        ierror= 1;
      }
      else {
        if (fclose(fphead)) {
          perror("Error closing header");
          ierror=1;
        }
      }
    }
    else {
      perror("Error reading from header");
      ierror= 1;
    }
  }
  else {
    perror("Error opening header");
    ierror= 1;
  }

  if (ierror) return 0;

  return 1;
}

