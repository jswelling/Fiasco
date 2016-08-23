 /************************************************************
 *                                                          *
 *  tiff_reader.c                                         *
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
#include <sys/mman.h>
#include <stdarg.h>
#include "mri.h"
#if (SGI64 || SGI5 || SGIMP)
#include <bstring.h>
#endif

#include "bio.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"
#include "smartreader.h"

static char rcsid[] = "$Id: tiff_reader.c,v 1.3 2007/06/18 22:46:32 welling Exp $";

/* Notes-
 * -Are the needed defs strings getting set?
 * -The Zeiss LSMInfo structure knows if it's multiple slices or 
 *  multiple times of the same slice, but we're ignoring that info.
 */

#ifdef USE_TIFF

#include <tiff.h>
#include <tiffio.h>

#define TIFFTAG_CZ_LSMINFO 34412
#define TIFF_VARIABLE -1

/* Carl Zeiss' custom tag for laser scanning microscope */
TIFFFieldInfo lsmFieldInfo= { 
  TIFFTAG_CZ_LSMINFO, /* tag */
  TIFF_VARIABLE,  /* read width */
  TIFF_VARIABLE,  /* write width */
  TIFF_BYTE,       /* type */
  FIELD_CUSTOM,   /* field_bit (usually CUSTOM) */
  FALSE,          /* oktochange */
  TRUE,           /* passcount */
  "CZ_LSMINFO" /* name */
};

typedef struct LSM_info_struct {
  uint32 magicNumber;
  int StructureSize;
  int DimensionX;
  int DimensionY;
  int DimensionZ;
  int DimensionChannels;
  int DimensionTime;
  int DataType;
  int ThumbnailX;
  int ThumbnailY;
  double VoxelSizeX;
  double VoxelSizeY;
  double VoxelSizeZ;
  uint32 ScanType;
  uint32 DataType2;
  uint32 OffsetVectorOverlay;
  uint32 OffsetInputLut;
  uint32 offsetOutputLut;
  uint32 OffsetChannelColors;
  double TimeInterval;
  uint32 OffsetChannelDataTypes;
  uint32 OffsetScanInformation;
  uint32 OffsetKsData;
  uint32 OffsetTimeStamps;
  uint32 OffsetEventList;
  uint32 OffsetRoi;
  uint32 OffsetBleachRoi;
  uint32 OffsetNextRecording;
  uint32 Reserved[90];
} LSMInfo;
  
/* Info about a particular IFD */
typedef struct ifd_data_struct {
  long long offset;
  uint32 imageWidth;
  uint32 imageLength;
  uint16 samplesPerPixel;
  uint16 bitsPerSample;
  uint16 planarConfig;
} IFDData;

/* A reader for TIFF files */
typedef struct tiff_data_struct {
  long long currentOffset;
  SList* ifdList;
  TIFF* tif;
  IFDData sampleIFD;
  long long lastVirtualOffsetRead;
  int passesPerDict;
  int scanlinesPerPass;
  int currentPass;
  int currentScanline;
} TiffData;

static TIFFExtendProc _ParentExtender = NULL;

static void _XTIFFDefaultDirectory(TIFF *tif)
{
  if (debug) fprintf(stderr,"_XTIFFDefaultDirectory is happening!\n");
  /* Install the extended Tag field info */
  TIFFMergeFieldInfo(tif, &lsmFieldInfo, 1);
  
  /* Since an XTIFF client module may have overridden
   * the default directory method, we call it now to
   * allow it to set up the rest of its own methods.
   */
  
  if (_ParentExtender) 
    (*_ParentExtender)(tif);
}

static void _XTIFFInitialize(void)
{
  /* We use this routine to tell libtiff about the tags used in the LSM
   * extension to the TIFF standard.
   */
  static int first_time=1;
  if (debug) fprintf(stderr,"_XTIFFInitialize is happening (part 1!)\n");
  
  if (! first_time) return; /* Been there. Done that. */
  first_time = 0;
  
  if (debug) fprintf(stderr,"_XTIFFInitialize is happening (part 2!)\n");
  /* Grab the inherited method and install */
  _ParentExtender = TIFFSetTagExtender(_XTIFFDefaultDirectory);
}

static IFDData* ifd_create( long long offset )
{
  IFDData* result= NULL;
  if (!(result=(IFDData*)malloc(sizeof(IFDData))))
    Abort("%s: unable to allocate %d bytes!\n",progname,sizeof(IFDData));
  result->offset= offset;
  return result;
}

static void ifd_destroy( void* victim )
{
  IFDData* d= (IFDData*)victim;
  free(d);
}

static void tiffDestroySelf( FileHandler* self )
{
  TiffData* data= (TiffData*)(self->hook);
  if (data->tif) {
    TIFFClose(data->tif);
    data->tif= NULL;
  }
  slist_destroy(data->ifdList,ifd_destroy);
  baseDestroySelf(self);
}

static tsize_t custom_ReadProc( thandle_t pSelf, tdata_t buf, tsize_t n )
{
  FileHandler* self= (FileHandler*)pSelf;
  TiffData* data= (TiffData*)(self->hook);
  size_t nRead;
  if (!(self->file)) {
    FH_REOPEN(self);
    if (fseek(self->file, data->currentOffset, SEEK_SET))
      perror("Reopened file and cannot seek!");
  }
  nRead= fread(buf, 1, (size_t)n, self->file);
  if (debug) fprintf(stderr,"Read %d bytes from %s\n",
		     (int)nRead,self->fileName);
  data->currentOffset += nRead;
  return nRead;
}

static tsize_t custom_WriteProc( thandle_t pSelf, tdata_t buf, tsize_t n )
{
  FileHandler* self= (FileHandler*)pSelf;
  TiffData* data= (TiffData*)(self->hook);
  size_t nWritten;
  if (!(self->file)) {
    FH_REOPEN(self);
    if (fseek(self->file, data->currentOffset, SEEK_SET))
      perror("Reopened file and cannot seek!");
  }
  nWritten= fwrite(buf, 1, (size_t)n, self->file);
  if (debug) fprintf(stderr,"Wrote %d bytes to %s\n",nWritten,self->fileName);
  data->currentOffset += nWritten;
  return nWritten;
}

static toff_t custom_SeekProc( thandle_t pSelf, toff_t offset, int whence )
{
  FileHandler* self= (FileHandler*)pSelf;
  TiffData* data= (TiffData*)(self->hook);
  int retcode;

  if (!(self->file)) {
    FH_REOPEN(self);
    if (fseek(self->file, data->currentOffset, SEEK_SET))
      perror("Reopened file and cannot seek!");
  }
  if (retcode= fseek(self->file, offset, whence)) {
    perror("Seek failed!");
    Abort("%s: seek failed on %s\n",progname,self->fileName);
  }
  data->currentOffset= ftell(self->file);
  if (debug) fprintf(stderr,"Seek to %lld in %s\n",
		     data->currentOffset,self->fileName);
  return data->currentOffset;
}

static int custom_CloseProc( thandle_t pSelf )
{
  FileHandler* self= (FileHandler*)pSelf;
  TiffData* data= (TiffData*)(self->hook);
  FH_CLOSE(self);
  if (debug) fprintf(stderr,"Closed %s\n",self->fileName);
  return 0;
}

static toff_t custom_SizeProc( thandle_t pSelf )
{
  FileHandler* self= (FileHandler*)pSelf;
  TiffData* data= (TiffData*)(self->hook);
  if (debug) fprintf(stderr,"Tiff reader queried size of %s: %lld bytes\n",
		     self->fileName, self->totalLengthBytes);
  return (int)self->totalLengthBytes;
}

static int custom_MapProc( thandle_t pSelf, tdata_t *start, toff_t *len )
{
  FileHandler* self= (FileHandler*)pSelf;
  TiffData* data= (TiffData*)(self->hook);
  /* It's not clear what the calling parameters are supposed to mean-
   * but we don't want to enable memory mapping anyway or the lib 
   * will suck the whole file into memory, which could be a disaster
   * if we have many large files.
   */
  if (debug) fprintf(stderr,"Tiff reader requested unimplemented mmap of %s\n",
		     self->fileName);
  return -1;
}

static void custom_UnmapProc( thandle_t pSelf, tdata_t start, toff_t offset )
{
  FileHandler* self= (FileHandler*)pSelf;
  TiffData* data= (TiffData*)(self->hook);
  if (debug) fprintf(stderr,
		     "Tiff reader requested unimplemented munmap of %s\n",
		     self->fileName);
}

static int scanForConsistency( SList* ifdList )
{
  IFDData* proto= NULL;
  slist_totop(ifdList);
  proto= (IFDData*)slist_next(ifdList);
  while (!slist_atend(ifdList)) {
    IFDData* thisIFD= (IFDData*)slist_next(ifdList);
    if ((thisIFD->imageWidth != proto->imageWidth)
	|| (thisIFD->imageLength != proto->imageLength)
	|| (thisIFD->samplesPerPixel != proto->samplesPerPixel)
	|| (thisIFD->bitsPerSample != proto->bitsPerSample)
	|| (thisIFD->planarConfig != proto->planarConfig))
      return 0;
  }
  return 1;
}

static void ifd_print(IFDData* ifd, FILE* ofile)
{
  fprintf(ofile,
	  "IFD at %lld: %d by %d, %hd samples at %hd bits, planar config %hd\n",
	  ifd->offset, (int)ifd->imageWidth, (int)ifd->imageLength,
	  ifd->samplesPerPixel, ifd->bitsPerSample, ifd->planarConfig);
}

static void fillOutIFDInfo( FileHandler* self, IFDData* ifd )
{
  TiffData* data= (TiffData*)(self->hook);
  if (!TIFFGetFieldDefaulted(data->tif, TIFFTAG_IMAGEWIDTH, 
			   &(ifd->imageWidth)))
    ifd->imageWidth= -1;
  if (!TIFFGetFieldDefaulted(data->tif, TIFFTAG_IMAGELENGTH, 
			   &(ifd->imageLength)))
    ifd->imageLength= -1;
  if (!TIFFGetFieldDefaulted(data->tif, TIFFTAG_SAMPLESPERPIXEL, 
			   &(ifd->samplesPerPixel)))
    ifd->samplesPerPixel= -1;
  if (!TIFFGetFieldDefaulted(data->tif, TIFFTAG_BITSPERSAMPLE,
			   &(ifd->bitsPerSample)))
    ifd->samplesPerPixel= -1;
  if (!TIFFGetFieldDefaulted(data->tif, TIFFTAG_PLANARCONFIG,
			   &(ifd->planarConfig)))
    ifd->planarConfig= -1;
}

static int checkForLSMInfo(FileHandler* self, LSMInfo* lsmInfo)
{
  TiffData* data= (TiffData*)(self->hook);
  char* bytes; /* because we declared the data type to be TIFF_BYTE */
  int n;

  if (TIFFGetField(data->tif, TIFFTAG_CZ_LSMINFO, &n, &bytes)) {
    if (debug) fprintf(stderr,"got LSM info! n= %d\n",n);
    bcopy(bytes, lsmInfo, sizeof(LSMInfo) > n ? n : sizeof(LSMInfo));
    return 1;
  }
  else return 0;
}

static void tiffProcessHeader( FileHandler* self, KVHash* info, 
			      SList* chunkStack )
{
  TiffData* data= (TiffData*)(self->hook);
  KVHash* defs= kvGetHash(info,"definitions");
  LSMInfo lsmInfo;
  int dircount= 0;
  uint32 tmp_uint32;
  int dz= 0;

  if (debug)
    fprintf(stderr,"Tiff library version %s\n",TIFFGetVersion());

  /* The 'm' after the mode option 'r' supresses memory mapping of the image */
  data->tif= TIFFClientOpen(self->fileName, "rm", self,
			    custom_ReadProc, custom_WriteProc, custom_SeekProc,
			    custom_CloseProc, custom_SizeProc, 
			    custom_MapProc, custom_UnmapProc);
  if (!(data->tif)) {
    Abort("%s: unexpectedly could not open <%s> as a TIFF file!\n",
	  progname,self->fileName);
  }

  /* The Open operation has pointed us at the first directory */
  dircount= 0;
  do {
    if (debug) {
      fprintf(stderr,"Directory %d:\n",++dircount);
      TIFFPrintDirectory(data->tif,stderr,0);
    }

    if (checkForLSMInfo(self, &lsmInfo)) {
      /* set some LSM-specific tags here */

      /* voxel sizes are in mm, so we must rescale VoxelSize info */
      kvDefDouble(info,"voxel_x",1000.0*lsmInfo.VoxelSizeX);
      kvDefDouble(info,"voxel_y",1000.0*lsmInfo.VoxelSizeY);
      kvDefDouble(info,"voxel_z",1000.0*lsmInfo.VoxelSizeZ);
    }

    if (TIFFGetFieldDefaulted(data->tif, TIFFTAG_SUBFILETYPE, &tmp_uint32)) {
      if (tmp_uint32==FILETYPE_REDUCEDIMAGE) {
	if (debug) fprintf(stderr,"Skipping reduced resolution image\n");
      }
      else if (tmp_uint32==FILETYPE_MASK) {
	if (debug) fprintf(stderr,"Skipping mask image\n");
      }
      else {
	IFDData* ifd= ifd_create(TIFFCurrentDirOffset(data->tif));
	if (TIFFIsTiled(data->tif)) 
	  Abort("%s: %s is a tiled TIFF; tiling is not currently supported.\n",
		progname,self->fileName);
	fillOutIFDInfo( self, ifd );
	if (debug) ifd_print( ifd, stderr );
	slist_append(data->ifdList, ifd);
      }
    }    
  }  while (TIFFReadDirectory(data->tif));

  if (debug)
    Message("There are %d dictionaries total\n",slist_count(data->ifdList));

  if (!scanForConsistency(data->ifdList)) {
    FH_CLOSE(self);
    Abort("%s: %s contains multiple images with inconsistent formats\n",
	  progname,self->fileName);
  }

  dz= slist_count(data->ifdList);
  if (dz==0)
    Abort("%s: no usable images in %s!\n",progname,self->fileName);

  slist_totop(data->ifdList);
  data->sampleIFD= *(IFDData*)slist_get(data->ifdList); /* copy out a sample */

  /* The skip pattern set up here is designed to force reading
   * to occur one scanline at a time, either for all channels or
   * for one channel at a time in the PLANARCONFIG_SEPARATE case.
   */
  if (data->sampleIFD.samplesPerPixel == 1) {
    if (dz==1) {
      kvDefString(info,"dimstr","xy");
    }
    else {
      kvDefString(info,"dimstr","xyz");
      kvDefInt(info,"dz",dz);
    }
    data->passesPerDict= 1;
  }
  else {
    if (data->sampleIFD.planarConfig==PLANARCONFIG_SEPARATE) {
      if (dz==1) {
	kvDefString(info,"dimstr","xyv");
      }
      else {
	kvDefString(info,"dimstr","xyvz");
	kvDefInt(info,"dz",dz);
	kvDefInt(info,"dv",data->sampleIFD.samplesPerPixel);
      }
      data->passesPerDict= data->sampleIFD.samplesPerPixel;
    }
    else {
      if (dz==1) {
	kvDefString(info,"dimstr","vxy");
      }
      else {
	kvDefString(info,"dimstr","vxyz");
	kvDefInt(info,"dz",dz);
	kvDefInt(info,"dv",data->sampleIFD.samplesPerPixel);
      }
      data->passesPerDict= 1;
    }
  }
    
  kvDefInt(info,"dx",data->sampleIFD.imageWidth);
  kvDefInt(info,"dy",data->sampleIFD.imageLength);
  data->scanlinesPerPass= data->sampleIFD.imageLength;
  kvDefInt(info,"skip.x",0);

  switch(data->sampleIFD.bitsPerSample) {
  case 8:
    {
      kvDefInt(info,"datatype_in",SRDR_UINT8);
      kvDefInt(info,"handler_datatype_out",SRDR_UINT8);
    }
    break;
  case 16:
    {
      kvDefInt(info,"datatype_in",SRDR_UINT16);
      kvDefInt(info,"handler_datatype_out",SRDR_UINT16);
    }
    break;
  default:
    Abort("%s: %s has %hd bits per sample (unsupported)!\n",
	  progname, data->sampleIFD.bitsPerSample);
  }

  FH_CLOSE(self); 
  slist_totop(data->ifdList); /* reader will walk the list */
  data->currentPass= data->currentScanline= 0; /* Prep for counting lines */
}

static void tiffRead( FileHandler* self, KVHash* info,
		     long long offset, long n,
		     SRDR_Datatype datatype, void* obuf )
{
  TiffData* data= (TiffData*)(self->hook);
  IFDData* currentIFD= (IFDData*)slist_get(data->ifdList);

  /* This routine relies on moving through the file sequentially! */
  if (data->lastVirtualOffsetRead != offset)
    Abort("%s: internal error in tiffRead: sequential reading rule was violated!\n",
	  progname);
  
  if (data->currentPass==0 && data->currentScanline==0) {
    /* Seek next dict */
    if (debug) fprintf(stderr,"Advancing directory to offset %lld\n",
		       currentIFD->offset);
    if (!TIFFSetSubDirectory(data->tif, currentIFD->offset))
      Abort("%s: error on TIFFSetSubDirectory of %s at offset %lld\n",
	    progname, self->fileName, currentIFD->offset);
    if (debug) fprintf(stderr,"Now at directory at offset %lld\n",
		       currentIFD->offset);
  }

  if (TIFFScanlineSize(data->tif)>n*srdrTypeSize[datatype])
    Abort("%s: internal error: scan line won't fit in buf!\n");
  TIFFReadScanline(data->tif, obuf, data->currentScanline, data->currentPass);
  if (debug) 
    fprintf(stderr,"Read scanline %d of %d, pass %d of %d, dict at %lld\n",
	    data->currentScanline, data->scanlinesPerPass,
	    data->currentPass, data->passesPerDict,
	    currentIFD->offset);

  data->currentScanline++;
  if (data->currentScanline >= data->scanlinesPerPass) {
    data->currentPass++;
    data->currentScanline= 0;
  }
  if (data->currentPass==data->passesPerDict) {
    data->currentPass= 0;
    slist_next(data->ifdList);
  }
  data->lastVirtualOffsetRead= offset+n*srdrTypeSize[datatype];
}

static int tiffCompare( FileHandler* f1, FileHandler* f2 )
{
  TiffData* data1= (TiffData*)(f1->hook);
  TiffData* data2= (TiffData*)(f2->hook);

  return (baseCompare(f1,f2));
}

FileHandler* tiffFactory(char* fname, KVHash* info)
{
  FileHandler* result= baseFactory(fname);
  TiffData* data;

  result->typeName= strdup("TiffDataHandler");

  result->destroySelf= tiffDestroySelf;
  result->read= tiffRead;
  result->processHeader= tiffProcessHeader;
  result->compareThisType= tiffCompare;

  if (!(data= (TiffData*)malloc(sizeof(TiffData))))
    Abort("%s: unable to allocate %d bytes!\n",sizeof(TiffData));
  result->hook= data;
  data->tif= NULL;
  data->currentOffset= 0;
  data->ifdList= slist_create();
  data->passesPerDict= data->scanlinesPerPass= 0;
  data->currentPass= data->currentScanline= 0;
  data->lastVirtualOffsetRead= 0;

  return result;
}

int tiffTester(const char* filename)
{
  TIFF* tif= NULL;

  _XTIFFInitialize();

  tif= TIFFOpen(filename,"r");
  if (tif) {
    TIFFClose(tif);
    return 1;
  }
  else {
    return 0;
  }
}

#else /* ifdef USE_TIFF */

FileHandler* tiffFactory(char* fname, KVHash* info)
{
  FileHandler* result= baseFactory("NotARealFile");
  return result;
}

int tiffTester(const char* filename)
{
  return 0;
}

#endif

