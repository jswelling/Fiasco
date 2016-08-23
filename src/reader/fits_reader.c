/************************************************************
 *                                                          *
 *  fits_reader.c                                         *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 2007 Department of Statistics,         *
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
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "mri.h"
#if (SGI64 || SGI5 || SGIMP)
#include <bstring.h>
#endif

#include "bio.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"
#include "smartreader.h"

static char rcsid[] = "$Id: fits_reader.c,v 1.2 2008/10/06 20:32:26 welling Exp $";

/* Notes-
 */

#ifdef USE_FITSIO

#include <fitsio.h>

/* A reader for FITS files */
typedef struct fits_data_struct {
  int isOpen;
  int hduNum;
  int iOwnHduList;
  int ndim;
  long* dimArray; /* Caution: stored in reverse order! */
  SList* hduList;
  fitsfile* fptr;
} FitsData;

static void fitsClose( FileHandler* self )
{
  FitsData* data= (FitsData*)(self->hook);
  if (data->isOpen) {
    int status= 0;
    if (fits_close_file(data->fptr, &status)) 
      fits_report_error(stderr, status); /* but do not abort */
    data->isOpen= 0;
    data->fptr= NULL;
  }
  baseClose(self);
}

#define TEST_WRAP_1( fitsfun, var1, msg )	\
  { \
    int status= 0; \
    FitsData* data= (FitsData*)(self->hook); \
    if (fitsfun(data->fptr, var1, &status)) { \
      fits_report_error(stderr,status); \
      Abort("%s: %s in HDU %d in %s!\n",progname,msg, \
	    data->hduNum,self->fileName); \
    } \
  }

#define TEST_WRAP_2( fitsfun, var1, var2, msg )	\
  { \
    int status= 0; \
    FitsData* data= (FitsData*)(self->hook); \
    if (fitsfun(data->fptr, var1, var2, &status)) {	\
      fits_report_error(stderr,status); \
      Abort("%s: %s in HDU %d in %s!\n",progname,msg, \
	    data->hduNum,self->fileName); \
    } \
  }

#define TEST_WRAP_4( fitsfun, var1, var2, var3, var4, msg )	\
  { \
    int status= 0; \
    FitsData* data= (FitsData*)(self->hook); \
    if (fitsfun(data->fptr, var1, var2, var3, var4, &status)) {	\
      fits_report_error(stderr,status); \
      Abort("%s: %s in HDU %d in %s!\n",progname,msg, \
	    data->hduNum,self->fileName); \
    } \
  }

#define TEST_WRAP_6( fitsfun,var1,var2,var3,var4,var5,var6,msg )	\
  { \
    int status= 0; \
    FitsData* data= (FitsData*)(self->hook); \
    if (fitsfun(data->fptr, var1, var2, var3, var4, \
		var5, var6, &status)) {		    \
      fits_report_error(stderr,status); \
      Abort("%s: %s in HDU %d in %s!\n",progname,msg, \
	    data->hduNum,self->fileName); \
    } \
  }

static void fitsReopen( FileHandler* self )
{
  FitsData* data= (FitsData*)(self->hook);
  /* Do not call baseReopen; it manipulates self->file ! */
  if (!(data->isOpen)) {
    int status= 0;
    if (fits_open_data(&(data->fptr),self->fileName,READONLY,&status)) {
      fits_report_error(stderr,status);
      Abort("%s: fatal fits error reopening %s!\n",progname,self->fileName);
    }
    data->isOpen= 1;
    if (data->hduNum) {
      int hdutype= 0;
      TEST_WRAP_2(fits_movabs_hdu, data->hduNum, &hdutype,
		  "error seeking");
    }
    else {
      /* FITS 'extended filename' may specify hdu */
      data->hduNum= fits_get_hdu_num(data->fptr,&(data->hduNum));
    }
  }
}

static void fitsDestroySelf( FileHandler* self )
{
  FitsData* data= (FitsData*)(self->hook);
  if (data->isOpen) FH_CLOSE(self);
  if (data->iOwnHduList && (data->hduList!=NULL))
    slist_destroy(data->hduList,NULL); 
  if (data->dimArray) free(data->dimArray);
  baseDestroySelf(self);
}

static void fitsZeroRead( FileHandler* self, KVHash* info,
			  long long offset, long n,
			  SRDR_Datatype datatype, void* obuf )
{
  bzero(obuf,n*srdrTypeSize[datatype]);
}

static inline int mapDataType( int fitsDataType )
{
  switch (fitsDataType) {
  case BYTE_IMG: return SRDR_UINT8;
  case SHORT_IMG: return SRDR_INT16;
  case LONG_IMG: return SRDR_INT32;
  case LONGLONG_IMG: return SRDR_INT64;
  case FLOAT_IMG: return SRDR_FLOAT32;
  case DOUBLE_IMG: return SRDR_FLOAT64;
  case USHORT_IMG: return SRDR_UINT16;
  }
  Abort("%s: unknown or unsupported FITS bitpix %d!\n",
	fitsDataType);
}

static inline int inverseMapDataType( int srdrDataType )
{
  /* Returns a FITS datatype, like 'TLONG' */
  switch (srdrDataType) {
  case SRDR_UINT8: return TBYTE;
  case SRDR_INT16: return TSHORT;
  case SRDR_INT32: return TINT;
  case SRDR_INT64: return TLONG;
  case SRDR_FLOAT32: return TFLOAT;
  case SRDR_FLOAT64: return TDOUBLE;
  case SRDR_UINT16: return TUSHORT;
  }
  Abort("%s: internal error: unknown SRDR data type %s!\n",
	srdrTypeName[srdrDataType]);
}

static void fitsRead( FileHandler* self, KVHash* info,
		     long long offset, long n,
		     SRDR_Datatype datatype, void* obuf )
{
  FitsData* data= (FitsData*)(self->hook);
  long* fpix= NULL;
  int i;
  long long count;
  int anynul;
  int fitsDataType= inverseMapDataType(datatype);

  if (!(fpix=(long*)malloc(data->ndim*sizeof(long))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname, data->ndim*sizeof(long));

  count= offset/srdrTypeSize[datatype];
  for (i=0; i<data->ndim; i++) {
    long dimSz= data->dimArray[i];
    fpix[i]= (long)(count%dimSz)+1; /* indices are 1-based */
    count /= dimSz;
    fprintf(stderr,"counting up: %d %ld %lld\n",
	    i,fpix[i],count);
    fprintf(stderr,"that was dimSz=%ld; count now %lld\n",dimSz,count);
    /*
    */
  }

  FH_REOPEN(self);
  fprintf(stderr,"Reading %ld starting at %ld %ld\n",n,fpix[0],fpix[1]);
  /*
  */
  TEST_WRAP_6(fits_read_pix, fitsDataType,
	      fpix, n, NULL, obuf, &anynul,
	      "unable to read data");
  if (debug) fprintf(stderr,"Read %ld of type %s; %s nulls\n",
		     n, srdrTypeName[datatype],
		     (anynul)?"some":"no");
    
  free(fpix);
}

static void transferHduData( FileHandler* self, KVHash* info )
{
  FitsData* data= (FitsData*)(self->hook);
  int bitpix;
  int inputDataType;
  long* dimArray;
  int i;
  char buf[256];
  char* dimstr= NULL;
  int nkeys;

  /* the FileHandler is open at this point */

  snprintf(buf,sizeof(buf),"images_%04d",data->hduNum-1);
  kvDefString(info,"chunkname",buf);
  snprintf(buf,sizeof(buf),".dat_%04d",data->hduNum-1);
  kvDefString(info,"chunkfile",buf);
  kvDefInt(info,"start_offset",0);

  TEST_WRAP_1(fits_get_img_equivtype, &bitpix, 
	      "cannot get bitpix");
  if (debug) fprintf(stderr,"incoming bitpix= %d\n",bitpix);
  inputDataType= mapDataType(bitpix);
  kvDefInt(info,"datatype_in",inputDataType);
  kvDefInt(info,"handler_datatype_out",inputDataType);

  TEST_WRAP_1(fits_get_img_dim, &(data->ndim), 
	      "cannot get number of axes");
  if (debug) fprintf(stderr,"HDU has %d dimensions\n",data->ndim);
  if (data->ndim<0) 
    Abort("%s: unsupported array of %d dimensions!",
	  progname,data->ndim);
  else if (data->ndim==0) {
    /* An unpleasant FITS convention which is incompatible with
     * PghMRI.  We must fake some data.
     */
    kvDefString(info,"dimstr","x");
    kvDefInt(info,"dx",1);
    kvDefString(info,"chunkname","fakedata");
    kvDefString(info,"chunkfile",".dat");
    self->read= fitsZeroRead;
  }
  else {
    if (!(data->dimArray=(long*)malloc(data->ndim*sizeof(long))))
      Abort("%s: unable to allocate %d long!\n",
	    progname,data->ndim);
    TEST_WRAP_2(fits_get_img_size,data->ndim,data->dimArray,
		"cannot get image axis lengths");
    /* Remember, FITS stores data in row-fastest order! */
    switch (data->ndim) {
    case 1: dimstr="x"; break;
    case 2: dimstr="yx"; break;
    case 3: dimstr="zyx"; break;
    case 4: dimstr="wzyx"; break;
    case 5: dimstr="edcba"; break;
    case 6: dimstr="fedcba"; break;
    default: Abort("%s: unsupported array of %d dimensions!",
		   progname,data->ndim);
    }
    kvDefString(info,"dimstr",dimstr);
    for (i=0; i<data->ndim; i++) {
      snprintf(buf,sizeof(buf),"d%c",dimstr[i]);
      kvDefInt(info,buf,data->dimArray[data->ndim-(i+1)]);
    }
  }
  
  TEST_WRAP_2(fits_get_hdrspace,&nkeys,NULL,
	      "unable to get key count");
  if (debug) fprintf(stderr,"There are %d keys\n",nkeys);
  for (i=1; i<=nkeys; i++) {
    char key[84];
    char value[84];
    char comment[84];
    /* I can't find anything that promises that fits_read_keyn will
     * not overwrite the end of the buffer.  These lengths are based 
     * on the notion that the 'card' is 80 bytes.
     */
    TEST_WRAP_4(fits_read_keyn, i, key, value, comment,
		"Unable to read key-value pair");
    key[80]= '\0';
    value[80]= '\0';
    comment[80]= '\0';
    if (debug) fprintf(stderr,"%d <%s> <%s> <%s>\n",i,key,value,comment);
    if (strcasecmp(key,"COMMENT")) kvDefString(info,key,value);
    
  }

}

static void fitsProcessHeader( FileHandler* self, KVHash* info, 
			      SList* chunkStack )
{
  FitsData* data= (FitsData*)(self->hook);
  KVHash* defs= kvGetHash(info,"definitions");
  int status= 0;

  /* We might be dealing with an entire file or an hdu within
   * a file.  If it's an entire file (signified by an impossible
   * hdunum) we will scan it for sub-files, calling 
   * fitsProcessHeader on each.
   */
  if (data->hduNum<1) {
    /* just created by factory call */
    int nFromFname;
    if (fits_parse_extnum(self->fileName, &nFromFname, &status)) {
	fits_report_error(stderr,status);
	Abort("%s: can't parse hdu field in %s!\n",
	      progname,self->fileName);
    }
    if (debug) 
      fprintf(stderr,"Got hdu %d from fname\n",nFromFname);
    if (nFromFname==-99) { /* signal value- this is a raw fname */
      int hdu;
      int nhdus;
      int hdutype= 0;
      SList* hduList= NULL;

      FH_REOPEN(self);

      /* Snag total number of HDUs */
      if (fits_get_num_hdus(data->fptr,&nhdus,&status)) {
	fits_report_error(stderr,status);
	Abort("%s: fits error counting HDUs on %s!\n",
	      progname,self->fileName);
      }

      if (debug) fprintf(stderr,"Got %d hdus\n",nhdus);
      /* Scan for image hdus */
      for (hdu=1; hdu<=nhdus; hdu++) {
	int hdutype;
	if (fits_movabs_hdu(data->fptr, hdu, &hdutype, &status)) {
	  fits_report_error(stderr,status);
	  Abort("%s: cannot move to HDU %d in %s!\n",progname,
		hdu,self->fileName);
	}
	if (hdutype==IMAGE_HDU) {
	  if (!hduList) {
	    hduList= slist_create();
	    data->iOwnHduList= 1;
	    slist_append(hduList,self);
	    data->hduList= hduList;
	    data->hduNum= hdu;
	    if (debug) 
	      fprintf(stderr,"Initialized HDU chain on HDU %d\n",
		      hdu);
	  }
	  else {
	    ChunkHandlerPair* chPair= NULL;
	    KVHash* subInfo= kvFactory(KV_DEFAULT_SIZE);
	    FitsData* kidData= NULL;
	    
	    if (!(chPair= (ChunkHandlerPair*)malloc(sizeof(ChunkHandlerPair))))
	      Abort("%s: unable to allocate %d bytes!\n",
		    progname,sizeof(ChunkHandlerPair));
	    initInfoHash(subInfo);
	    chPair->info= subInfo;
	    if (debug) fprintf(stderr,"Chaining hdu %d\n",hdu);
	    chPair->handler= fitsFactory(self->fileName, subInfo);
	    kidData= (FitsData*)(chPair->handler->hook);
	    kidData->hduNum= hdu;
	    kidData->iOwnHduList= 0;
	    kidData->hduList= hduList;
	    slist_append(chunkStack,chPair);
	    slist_append(data->hduList, chPair->handler);
	  }
	}
	else {
	  if (debug) fprintf(stderr,"skipping non-image HDU %d\n",
			     hdu);
	}
      }
      /* Get back to where we belong */
      TEST_WRAP_2(fits_movabs_hdu, data->hduNum, &hdutype,
		  "error seeking");
    }
    else {
      FH_REOPEN(self);
      data->hduNum= nFromFname;
      data->iOwnHduList= 0;
      data->hduList= NULL;
    }
  }
  else {
    FH_REOPEN(self);
  }

  if (data->iOwnHduList && (data->hduList!=NULL))
    kvDefInt(info,"num_hdus",slist_count(data->hduList));
  else
    kvDefInt(info,"num_hdus",1);
  kvDefInt(info,"current_hdu",data->hduNum);

  transferHduData(self, info);
#ifdef never
  fits_uint_32 width;
  fits_uint_32 height;
  int bit_depth;
  int color_type;
  int interlace_type;
  int compression_type;
  int filter_method;
  int nChannels;
  fits_uint_32 resX;
  fits_uint_32 resY;
  fits_timep timep;
  fits_text* comments= NULL;
  int nComments;
  char str[64];
  char str2[256];
  int i;

  FH_REOPEN(self);
  data->currentOffset= 0;

  if (setjmp(fits_jmpbuf(data->fits_ptr))) {
        fits_destroy_read_struct(&(data->fits_ptr), &(data->info_ptr),
           &(data->end_info));
	FH_CLOSE(self);
        Abort("%s: error within libfits on file %s!\n",progname,self->fileName);
    }

  fits_set_read_fn( data->fits_ptr, self, custom_read_data );
  fits_read_info(data->fits_ptr, data->info_ptr);
  fits_get_IHDR(data->fits_ptr, data->info_ptr, &width, &height,
	       &bit_depth, &color_type, &interlace_type, &compression_type,
	       &filter_method);
  nChannels= fits_get_channels(data->fits_ptr, data->info_ptr);
  if (fits_get_tIME(data->fits_ptr, data->info_ptr, &(timep))
      == FITS_INFO_tIME) {
    bcopy(timep,&(data->modtime),sizeof(fits_time)); 
    snprintf(str,sizeof(str),"%02d/%02d/%04d",
	     (int)data->modtime.month,
	     (int)data->modtime.day,
	     (int)data->modtime.year);
    kvDefString(info,"date",str);
    kvDefString(defs,"date","date modified");
    snprintf(str,sizeof(str),"%02d:%02d:%02d",
	     (int)data->modtime.hour,
	     (int)data->modtime.minute,
	     (int)data->modtime.second);
    kvDefString(info,"time",str);
    kvDefString(defs,"time","time modified");
  }
  else {
    if (debug) fprintf(stderr,"time info in %s is not valid\n",
		       self->fileName);
    bzero(&(data->modtime),sizeof(fits_time));
  }
  if (nChannels==1) kvDefString(info,"dimstr","xy");
  else {
    kvDefString(info,"dimstr","vxy");
    kvDefInt(info,"dv",nChannels);
  }
  kvDefInt(info,"dx",(long)width);
  kvDefInt(info,"dy",(long)height);
  kvDefInt(info,"fits_color_type",color_type);
  kvDefInt(info,"fits_bit_depth",bit_depth);
  kvDefInt(info,"fits_interlace_type",interlace_type);
  kvDefInt(info,"fits_compression_type",compression_type);
  kvDefInt(info,"fits_filter_method",filter_method);

  if (color_type == FITS_COLOR_TYPE_PALETTE)
    fits_set_palette_to_rgb(data->fits_ptr);
  
  /* We can coerce libfits into inflating 1, 2, or 4 bit data to 8,
   * so we only need to worry about the 8 and 16 bit depth cases.
   */
  if (color_type == FITS_COLOR_TYPE_GRAY &&
      bit_depth < 8) fits_set_gray_1_2_4_to_8(data->fits_ptr);
  if (bit_depth==16) {
    kvDefInt(info,"datatype_in",SRDR_UINT16);
    Warning(1,"%s: FITS image pixel data is unsigned shorts; this format is not completely supported!\n",progname);
  }
  else {
    kvDefInt(info,"datatype_in",SRDR_UINT8);
  }
  kvDefInt(info,"handler_datatype_out",kvGetInt(info,"datatype_in"));
 
  nComments= fits_get_text(data->fits_ptr, data->info_ptr, &comments, NULL);
  for (i=0; i<nComments; i++) {
    char* here;
    snprintf(str,sizeof(str),"fits_txt_%s",comments[i].key);
    strncpy(str2,comments[i].text,sizeof(str2));
    here= str2+strlen(str2)-1;
    while (isspace(*here) && here>=str2) *here--= '\0';
    kvDefString(info,str,str2);
  }

  resX= fits_get_x_pixels_per_meter(data->fits_ptr,data->info_ptr);
  resY= fits_get_y_pixels_per_meter(data->fits_ptr,data->info_ptr);
  if (resX != 0) {
    kvDefDouble(info,"voxel_x",1000.0/resX);
    kvDefString(defs,"voxel_x","X voxel size including gap (mm)");
  }
  if (resY != 0) {
    kvDefDouble(info,"voxel_y",1000.0/resY);
    kvDefString(defs,"voxel_y","Y voxel size including gap (mm)");
  }

  fits_read_update_info(data->fits_ptr, data->info_ptr);
#endif

  FH_CLOSE(self);
}

#ifdef never
static int fitsCompare( FileHandler* f1, FileHandler* f2 )
{
  FitsData* data1= (FitsData*)(f1->hook);
  FitsData* data2= (FitsData*)(f2->hook);

  if (0) {} 
  else return (baseCompare(f1,f2));
}
#endif

FileHandler* fitsFactory(char* fname, KVHash* info)
{
  FileHandler* result= baseFactory(fname);
  FitsData* data;
  int status= 0;

  result->typeName= strdup("FITSDataHandler");
  result->file= NULL; /* only manipulate it via fptr */
  result->close= fitsClose;
  result->reopen= fitsReopen;
  result->destroySelf= fitsDestroySelf;
  result->read= fitsRead;
  result->processHeader= fitsProcessHeader;
#ifdef never
  result->compareThisType= fitsCompare;
#endif

  if (!(data= (FitsData*)malloc(sizeof(FitsData))))
    Abort("%s: unable to allocate %d bytes!\n",progname,sizeof(FitsData));
  result->hook= data;

  /* Force reopen method to open the file for us */
  data->isOpen= 0;
  data->hduNum= 0; /* which is an impossible value */
  data->hduList= NULL;
  data->ndim= 0;
  data->dimArray= NULL;

  return result;
}

int fitsTester(const char* filename)
{
  int match= 0;
  fitsfile* fptr= NULL;
  int status= 0;

  if (fits_open_file(&fptr, filename, READONLY, &status)) {
    match= 0;
  }
  else {
    match= 1;
  }
  fits_close_file(fptr, &status);
  return match;
}

#undef TEST_WRAP_1
#undef TEST_WRAP_2
#undef TEST_WRAP_6

#else /* ifdef USE_FITSIO */

FileHandler* fitsFactory(char* fname, KVHash* info)
{
  FileHandler* result= baseFactory("NotARealFile");
  return result;
}

int fitsTester(const char* filename)
{
  return 0;
}

#endif

