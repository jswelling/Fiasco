/************************************************************
 *                                                          *
 *  smartreader.h                                           *
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
 ************************************************************/

/* Hash table structure and methods */
#include "kvhash.h"

/* List structure and methods */
#include "slist.h"

/*
 * Smartreader's own notion of data storage types
 */
typedef int SRDR_Datatype;
#define SRDR_UINT8	0	/* 8 bit unsigned integer */
#define SRDR_INT16	1	/* 16 bit signed integer */
#define SRDR_UINT16	2	/* 16 bit unsigned integer */
#define SRDR_INT32	3	/* 32 bit integer */
#define SRDR_FLOAT32	4	/* 32 bit IEEE floating point */
#define SRDR_FLOAT64	5	/* 64 bit IEEE floating point */
#define SRDR_INT64	6	/* 64 bit integer */


/* 
 * Structs for file format handlers.
 */

typedef struct file_handler_struct {
  char* fileName;
  char* typeName;
  long long totalLengthBytes;
  FILE *file;
  void (*processHeader)( struct file_handler_struct*, KVHash*, SList* );
  void (*destroySelf)( struct file_handler_struct * );
  void (*read)( struct file_handler_struct *, KVHash*,
		long long offset, long n, 
		SRDR_Datatype datatype_out, void* obuf );
  void (*close)( struct file_handler_struct * );
  void (*reopen)( struct file_handler_struct * ); 
  int (*compareThisType)( struct file_handler_struct *,
			  struct file_handler_struct * ); /* used by sorts */ 
  void* hook;
} FileHandler;

typedef FileHandler* (*FileHandlerFactory)(char*, KVHash*);

typedef struct file_handler_pair_struct {
  int (*tester)(const char*);
  FileHandlerFactory factory;
  const char* name;
} FileHandlerPair;

typedef struct chunk_handler_pair_struct {
  KVHash* info;
  FileHandler* handler;
} ChunkHandlerPair;

typedef struct key_type_pair_struct {
  const char* name;
  int type;
} KeyTypePair;

/* Some easy access macros */
#define FH_PROCESSHEADER( fh, kvh, chunkStack ) \
  (*(fh->processHeader))( fh, kvh, chunkStack )
#define FH_DESTROYSELF( fh ) (*(fh->destroySelf))( fh )
#define FH_READ( fh, kvh, offset, size, dt_out, obuf ) \
  (*(fh->read))( fh, kvh, offset, size, dt_out, obuf )
#define FH_CLOSE( fh ) (*(fh->close))( fh )
#define FH_REOPEN( fh ) (*(fh->reopen))( fh )
#define FH_COMPARE( fh1, fh2 ) (*(fh1->compareThisType))( fh1, fh2 )

/* A global debug flag */
extern int debug;

/* A global verbosity value */
extern int verbose_flg;

/* The program name */
extern char* progname;

/* Some tools for handling wildcarded filenames.
 * The caller owns all the memory returned by these functions.
 */
extern char* expandWildcardFirstFilename( char* fname );
extern SList* expandWildcardFilename( char* fname );
extern void destroyWildcardFilenameList( SList* list );

/* Some tools for table handling */
extern long srdrTypeSize[7];
extern char* srdrTypeName[7];
extern long libmriTypeSize[6];
extern char* libmriTypeName[6];
extern long libmriArrayTypeSize[8];
extern char* libmriArrayTypeName[8];
void loadStringTransTable( KVHash* kvh, char* tbl[][2] );
void loadStringTypeTable( KVHash* kvh, KeyTypePair* tbl );
void initInfoHash(KVHash* info);
int stringTableLookup(const char* s, const char** tbl);

/* 
 * Utilities provided by the smart_utils module 
*/
void smart_dump(KVHash* info);
int legitDimString( const char* str );
int consistency_check(KVHash* info);
void emit_format_summary(KVHash* info, FileHandler* handler);
int smart_reconcile(KVHash* info);

/* 3D math routines */
int getVec3(KVHash* info, char* name, double* result);
int getVec3_anat(KVHash* info, char* name, double* result);
int testVec3(KVHash* info, char* name);
void defVec3(KVHash* info, char* name, double* vals);
void subtractVec3( double* result, const double* v1, const double* v2 );
void multVec3( double* result, const double* v, double factor );
void copyVec3( double* to, const double* from );
double dotVec3( const double* v1, const double* v2 );
double normVec3( const double* v );
void crossVec3( double* result, const double* v1, const double* v2 );
void xplusbyVec3( double* result, const double* v1, const double *v2, 
		  double scale );
void normalizeVec3( double* v );
void flipToPositiveHemisphereVec3( double* vec );

/* A "base class" to supply common functionality.  All of its methods
 * are public so that the "derived classes" can call them.
 */
extern FileHandler* baseFactory();
extern void baseDestroySelf( FileHandler* self );
extern void baseProcessHeader( FileHandler* self, KVHash* info, 
			       SList* cStack );
extern void baseRead( FileHandler* self, KVHash* info,
		      long long offset, long n,
		      SRDR_Datatype datatype, void* obuf );
extern void baseClose( FileHandler* self );
extern void baseReopen( FileHandler* self );
extern int baseCompare( FileHandler* f1, FileHandler* f2 );
extern int bigfile_fseek (FILE *f, long long offset, int whence);

/* General raw data handler */
extern int rawTester(const char* fname);
extern FileHandler* rawFactory(char* fname, KVHash* info);

/* GE LX data handler, and some associated routines.  These
 * exist in multiple versions, depending on the setting of 
 * compile-time flags. */

extern int lxTester_excite(const char* fname);
extern FileHandler* lxFactory_excite(char* fname, KVHash* info);
extern void scanSplxHeader_excite( KVHash* info, unsigned char*, 
				 unsigned char*, unsigned char*, 
				 unsigned char*, unsigned char* );
extern void scanSf11Header_excite( KVHash* info, unsigned char*, 
				 unsigned char*, unsigned char*, 
				 unsigned char*, unsigned char* );
extern void scanEpiboldHeader_excite( KVHash* info, unsigned char*, 
				    unsigned char*, unsigned char*, 
				    unsigned char*, unsigned char* );
extern void scan2dfastHeader_excite( KVHash* info, unsigned char*, 
				    unsigned char*, unsigned char*, 
				    unsigned char*, unsigned char* );
void addSplxTrajChunk_excite(FileHandler* self, KVHash* info, SList* cStack);
void addSf11TrajChunk_excite(FileHandler* self, KVHash* info, SList* cStack);
void addBandpassChunk_excite(FileHandler* self, KVHash* info, SList* cStack);
void addPhaseRefChunk_excite(FileHandler* self, KVHash* info, SList* cStack);
void addRampSampleChunk_excite(FileHandler* self, KVHash* info, SList* cStack);

extern int lxTester_cnv4(const char* fname);
extern FileHandler* lxFactory_cnv4(char* fname, KVHash* info);
extern void scanSplxHeader_cnv4( KVHash* info, unsigned char*, 
				 unsigned char*, unsigned char*, 
				 unsigned char*, unsigned char* );
extern void scanSf11Header_cnv4( KVHash* info, unsigned char*, 
				 unsigned char*, unsigned char*, 
				 unsigned char*, unsigned char* );
extern void scanEpiboldHeader_cnv4( KVHash* info, unsigned char*, 
				    unsigned char*, unsigned char*, 
				    unsigned char*, unsigned char* );
extern void scan2dfastHeader_cnv4( KVHash* info, unsigned char*, 
				    unsigned char*, unsigned char*, 
				    unsigned char*, unsigned char* );
void addSplxTrajChunk_cnv4(FileHandler* self, KVHash* info, SList* cStack);
void addSf11TrajChunk_cnv4(FileHandler* self, KVHash* info, SList* cStack);
void addBandpassChunk_cnv4(FileHandler* self, KVHash* info, SList* cStack);
void addPhaseRefChunk_cnv4(FileHandler* self, KVHash* info, SList* cStack);
void addRampSampleChunk_cnv4(FileHandler* self, KVHash* info, SList* cStack);

extern int lxTester_lx2(const char* fname);
extern FileHandler* lxFactory_lx2(char* fname, KVHash* info);
extern void scanSplxHeader_lx2( KVHash* info, unsigned char*, 
				 unsigned char*, unsigned char*, 
				 unsigned char*, unsigned char* );
extern void scanSf11Header_lx2( KVHash* info, unsigned char*, 
				 unsigned char*, unsigned char*, 
				 unsigned char*, unsigned char* );
extern void scanEpiboldHeader_lx2( KVHash* info, unsigned char*, 
				    unsigned char*, unsigned char*, 
				    unsigned char*, unsigned char* );
extern void scan2dfastHeader_lx2( KVHash* info, unsigned char*, 
				    unsigned char*, unsigned char*, 
				    unsigned char*, unsigned char* );
void addSplxTrajChunk_lx2(FileHandler* self, KVHash* info, SList* cStack);
void addSf11TrajChunk_lx2(FileHandler* self, KVHash* info, SList* cStack);
void addBandpassChunk_lx2(FileHandler* self, KVHash* info, SList* cStack);
void addPhaseRefChunk_lx2(FileHandler* self, KVHash* info, SList* cStack);
void addRampSampleChunk_lx2(FileHandler* self, KVHash* info, SList* cStack);

extern int lxTester_prelx(const char* fname);
extern FileHandler* lxFactory_prelx(char* fname, KVHash* info);
extern void scanSplxHeader_prelx( KVHash* info, unsigned char*, 
				  unsigned char*, unsigned char*, 
				  unsigned char*, unsigned char* );
extern void scanSf11Header_prelx( KVHash* info, unsigned char*, 
				  unsigned char*, unsigned char*, 
				  unsigned char*, unsigned char* );
extern void scanEpiboldHeader_prelx( KVHash* info, unsigned char*, 
				     unsigned char*, unsigned char*, 
				     unsigned char*, unsigned char* );
extern void scan2dfastHeader_prelx( KVHash* info, unsigned char*, 
				    unsigned char*, unsigned char*, 
				    unsigned char*, unsigned char* );
void addSplxTrajChunk_prelx(FileHandler* self, KVHash* info, SList* cStack);
void addSf11TrajChunk_prelx(FileHandler* self, KVHash* info, SList* cStack);
void addBandpassChunk_prelx(FileHandler* self, KVHash* info, SList* cStack);
void addPhaseRefChunk_prelx(FileHandler* self, KVHash* info, SList* cStack);
void addRampSampleChunk_prelx(FileHandler* self, KVHash* info, SList* cStack);

/* WINDAQ data handler */
extern int windaqTester(const char* fname);
extern FileHandler* windaqFactory(char* fname, KVHash* info);

/* DE Smith data handler */
extern int desmithTester(const char* fname);
extern FileHandler* desmithFactory(char* fname, KVHash* info);

/* AFNI data handler */
extern int afniTester(const char* fname);
extern FileHandler* afniFactory(char* fname, KVHash* info);

/* DICOM data handler */
extern int dicomTester(const char* fname);
extern FileHandler* dicomFactory(char* fname, KVHash* info);

/* Siemens k-space format handler */
extern int siemensKspaceTester(const char* fname);
extern FileHandler* siemensKspaceFactory(char* fname, KVHash* info);

/* A fake "reader" for data in ram */
/* Note that it owns the data buffer after it is created! */
FileHandler* ramDataHandlerFactory(void* buf, long length, SRDR_Datatype type);

/* GE LX Image data handler */
extern int lxImageTester(const char* fname);
extern FileHandler* lxImageFactory(char* fname, KVHash* info);

/* NIfTI format data handler */
extern int niftiTester(const char* fname);
extern FileHandler* niftiFactory(char* fname, KVHash* info);

/* ANALYZE format data handler */
extern int analyzeTester(const char* fname);
extern FileHandler* analyzeFactory(char* fname, KVHash* info);

/* Pgh MRI format data handler */
extern int pghmriTester(const char* fname);
extern FileHandler* pghmriFactory(char* fname, KVHash* info);

/* unsigned short data handler */
extern int ushortTester(const char* fname);
extern FileHandler* ushortFactory(char* fname, KVHash* info);

/* PNG format data handler */
extern int pngTester(const char* fname);
extern FileHandler* pngFactory(char* fname, KVHash* info);

/* FITS format data handler */
extern int fitsTester(const char* fname);
extern FileHandler* fitsFactory(char* fname, KVHash* info);

/* TIFF format data handler */
extern int tiffTester(const char* fname);
extern FileHandler* tiffFactory(char* fname, KVHash* info);

/* FIFF format data handler */
extern int fiffTester(const char* fname);
extern FileHandler* fiffFactory(char* fname, KVHash* info);

/* SON format data handler */
extern int sonTester(const char* fname);
extern FileHandler* sonFactory(char* fname, KVHash* info);

/*
 * These special-purpose handlers have nonstandard APIs
 */

/* A reader to handle type conversions */
FileHandler* convertFactory(FileHandler* child, 
			    SRDR_Datatype dt_in, SRDR_Datatype dt_out);

/* A reader to handle multiple files in series */
FileHandler* multiFileHandlerFactory();
void multiFileHandlerAddFile(FileHandler* multi, FileHandler* newChild);

