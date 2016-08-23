/************************************************************
 *                                                          *
 *  smartreader.c                                             *
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
 *  Major restructuring, and the addition of formidable     *
 *       flexibility and intelligence, Joel Welling 5-2002  *
 ************************************************************/

/* Notes-

   "datatype_in" is the nominal datatype of data on disk; the offset
   passed to FileHandler.read() advances by n*srdrTypeSize[datatype_in]
   each time n data objects are read.

   "handler_datatype_out" is the datatype which FileHandler.read()
   delivers to the buffer pointed to by obuf.

   "datatype_out" is the datatype actually delivered to the output
   Pgh MRI file.  The smartreader top-level routine will wrap the
   FileHandler in a converter if handler_datatype_out differs from
   datatype_out.

   General Questions:
   -Should export both voxel_size and voxel_spacing in X and Y
    (with the same values)
   -Original spiral header.c line 253 does a correction to ctr[2]
    which is not implemented here.  What does it mean?
   -All the rotation and transpose stuff in spiral/header.c presumably
    should happen in general LX processing as well.
   -The first chunk created always has autoscale attributes.  That may not
    be what is desired; it makes no sense for some chunks.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include "mri.h"
#if (SGI64 || SGI5 || SGIMP)
#include <bstring.h>
#endif
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

#include "bio.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"
#include "smartreader.h"

static char rcsid[] = "$Id: smartreader.c,v 1.41 2008/10/08 18:52:39 welling Exp $";

/* An array of data format handlers.  Each of the factories in this
 * table must take a parameter list of the form: 
 *
 *    someFactory(char* fileName, KVHash* contextInfo)
 *
 * Another factory exists- ushortFactory- but it is never used
 * directly to parse headers.
 */
static FileHandlerPair handlerTable[]= {
  { pghmriTester, pghmriFactory, (const char*)"Pittsburgh MRI" },
  { lxTester_excite, lxFactory_excite, (const char*)"GE LX (excite)" },
  { lxTester_cnv4, lxFactory_cnv4, (const char*)"GE LX (cnv4)" },
  { lxTester_lx2, lxFactory_lx2, (const char*)"GE LX (lx2)" },
  { lxTester_prelx, lxFactory_prelx, (const char*)"GE LX (pre-lx)" },
  { windaqTester, windaqFactory, (const char*)"Windaq" },
  { desmithTester, desmithFactory, (const char*)"DE Smith Image" },
  { lxImageTester, lxImageFactory, (const char*)"GE LX Image" },
  { afniTester, afniFactory, (const char*)"AFNI" },
  { dicomTester, dicomFactory, (const char*)"DICOM" },
  /* must check for NIfTI before ANALYZE, since NIfTI files qualify
   * as ANALYZE format but not vice versa.
   */
  { niftiTester, niftiFactory, (const char*)"NIfTI" },
  { analyzeTester, analyzeFactory, (const char*)"ANALYZE" },
#ifdef USE_FIFF
  { fiffTester, fiffFactory, (const char*)"FIFF" },
#endif
#ifdef USE_SON
  { sonTester, sonFactory, (const char*)"SON" },
#endif
#ifdef USE_PNG
  { pngTester, pngFactory, (const char*)"PNG" },
#endif
#ifdef USE_TIFF
  { tiffTester, tiffFactory, (const char*)"TIFF" },
#endif
#ifdef USE_FITSIO
  { fitsTester, fitsFactory, (const char*)"FITS" },
#endif
  { siemensKspaceTester, siemensKspaceFactory, (const char*)"Siemens K-space"},
  /* RAW reader comes last, last, since it handles everything */
  { rawTester, rawFactory, (const char*)"raw" } 
};

/* We need to place a limit on the maximum number of "things" in one IO op */
#define MAX_BLOCK (16*1024*1024)

int debug = 0;          /* Global debug flag                         */

int verbose_flg = 0;        /* Global verbosity value (0=off) */

char* progname= NULL; /* program name */

static MRI_Dataset* open_output( char* hdrfile )
{
  MRI_Dataset* result;

  /* Initialize header */
  result = mri_open_dataset( hdrfile, MRI_WRITE );
  /* We're going to wait until later to add the command line to
   * the history file.
   */
  
  return result;
}

static void export_relevant_tags( MRI_Dataset* ds, const char* chunk, 
			       KVHash* info )
{
  /* More accurately, we'll set all tags that aren't irrelevant */
  KVHash* defs= kvGetHash(info,"definitions");
  KVHash* extNames= kvGetHash(info,"external_names");
  KVIterator* kvi= kvUniqueIteratorFactory(info);
  char* dimstr= kvGetString(info,"dimstr");
  int maxHistory= 0;
  char buf[256];
  int i;

  if (strlen(chunk)>sizeof(buf)/2) {
    /* let's be reasonable here! */
    Abort("%s: chunk name <%s> is too long!\n",progname,chunk);
  }

  while (kvIteratorHasMorePairs(kvi)) {
    KVPair* p= kvIteratorNextPair(kvi);
    int i;
    int ignore= 0;
    const char* outname;

    /* Check to see if this is a dimension key; they're handled elsewhere */
    if (strlen(kvKey(p))==2 && kvKey(p)[0]=='d' 
	&& (strchr(dimstr,kvKey(p)[1])!=NULL)) continue;

    /* Check to see if this is a skip key; they're no longer relevant */
    if (strlen(kvKey(p))==6 && !strncmp(kvKey(p),"skip.",5) 
	&& (strchr(dimstr,kvKey(p)[5])!=NULL)) continue;
    
    /* Get the name translation, if any.  Translation to the empty string
     * means to ignore this pair.
     */
    if (kvLookup(extNames,kvKey(p))) {
      /* An empty external name means ignore this key */
      if (strlen(outname= kvGetString(extNames,kvKey(p))) == 0) continue;
    }
    else outname= kvKey(p);
    
    /* If the key is a history value, we want to add it as a history entry
     * in the new file.  We want to preserve history order,however, so
     * for the moment we just note that there is such a history entry.
     * We'll add all the entries in order later.
     *
     * This allows the various readers to set history info if desired.  
     * Otherwise, the key gets appended to the current chunk name.
     */
    if (!strncmp(outname,"history.",strlen("history."))) {
      int histNum= atoi(outname+strlen("history."));
      if (histNum>maxHistory) maxHistory= histNum;
    }
    else {
      strncpy(buf,chunk,sizeof(buf)-2);
      strcat(buf,".");
      strncat(buf,outname,sizeof(buf)-strlen(buf));
      buf[sizeof(buf)-1]= '\0';
      /* There are some ideosyncracies in the actual argument types of
       * the following function calls.
       */
      switch (kvType(p)) {
      case KV_STRING:
	mri_set_string(ds, buf, kvGetString(info,kvKey(p)));
	break;
      case KV_LONG:
	mri_set_int(ds, buf, kvGetLong(info,kvKey(p)));
	break;
      case KV_DOUBLE:
	mri_set_float(ds, buf, (float)kvGetDouble(info,kvKey(p)));
	break;
      case KV_HASH:
	/* Do nothing */
	break;
      case KV_BOOLEAN:
	mri_set_int(ds, buf, (kvGetBoolean(info,kvKey(p)) ? 1 : 0));
	break;
      case KV_INT:
	mri_set_int(ds, buf, kvGetInt(info,kvKey(p)));
	break;
      }
    }
  }
  /* Look back over the input information for history lines, and add
   * them in order.
   */
  for (i=0; i<=maxHistory; i++) {
    sprintf(buf,"history.%d",i);
    if (kvLookup(info,buf)) {
      hist_add(ds, kvGetString(info,buf));
    }
  }
}

static void recursiveTransfer( FileHandler* handler, KVHash* info,
			       MRI_Dataset* ds, const char* chunk, 
			       long long *in_offset, 
			       long long *out_offset,
			       int myDepth, int minDepth,
			       long long sizeAtMinDepth )
{
  /*
   * Note: in_offset is in bytes, but out_offset is in elements!!!
   */
  char buf[64];
  char* dimstr= kvGetString(info,"dimstr");
  int n;
  long long skip;
  
  sprintf(buf,"d%c",dimstr[myDepth]);
  n= kvGetInt(info,buf);
  sprintf(buf,"skip.%c",dimstr[myDepth]);
  skip= ((kvLookup(info,buf)==NULL)? 0 : kvGetLong(info,buf));

  if (myDepth>minDepth) {
    int i;
    for (i=0; i<n; i++) {
      recursiveTransfer(handler, info, ds, chunk, 
			in_offset, out_offset,
			myDepth-1, minDepth, sizeAtMinDepth);
    }
    (*in_offset) += skip;
  }
  else {
    void *buf;
    long long numThisBlock;
    long long numToGo= sizeAtMinDepth;
    SRDR_Datatype dt_in= kvGetInt(info,"datatype_in");
    SRDR_Datatype dt_out= kvGetInt(info,"datatype_out");
    MRI_ArrayType mri_arraytype_out;
    if (!(buf= malloc(MAX_BLOCK*srdrTypeSize[dt_out])))
      Abort( "%s: unable to allocate %d bytes!\n",
	     progname, MAX_BLOCK*srdrTypeSize[dt_out]); 
    switch (dt_out) {
    case SRDR_UINT8: mri_arraytype_out= MRI_UNSIGNED_CHAR; break;
    case SRDR_INT16: mri_arraytype_out= MRI_SHORT; break;
    case SRDR_INT32: mri_arraytype_out= MRI_INT; break;
    case SRDR_FLOAT32: mri_arraytype_out= MRI_FLOAT; break;
    case SRDR_FLOAT64: mri_arraytype_out= MRI_DOUBLE; break;
    case SRDR_INT64: mri_arraytype_out= MRI_LONGLONG; break;
    default:
      Abort("%s: internal error: cannot write type %s to a Pgh MRI file!\n",
	    progname, srdrTypeName[dt_out]);
    }
    while (numToGo>0) {
      numThisBlock=  (numToGo>MAX_BLOCK) ? MAX_BLOCK : numToGo;
      FH_READ( handler, info, *in_offset, numThisBlock, dt_out, buf );
      mri_set_chunk( ds, chunk, numThisBlock, *out_offset, mri_arraytype_out, 
		     buf );

      /* in_offset is in bytes, but out_offset is in elements! */
      (*in_offset) += (numThisBlock*srdrTypeSize[dt_in]);
      (*out_offset) += numThisBlock;
      numToGo -= numThisBlock;
    }
    (*in_offset) += skip;
    free(buf);
  }
}

static void transfer( FileHandler* handler, KVHash* info, MRI_Dataset* ds )
{
  long long sizeAtMinDepth;
  long long in_offset;
  long long out_offset;
  int minDepth;
  int i;
  char* dimstr= kvGetString(info,"dimstr");
  char skipstring[64];
  char extentstring[64];
  char buf[256];
  char* chunk= kvGetString(info,"chunkname");
  MRI_Datatype mri_datatype_out;
  
  /* Note that we want even skip lengths of 0 to terminate this
   * recursion, so that readers can force a break in the read 
   * block size if they need to.
   */
  sizeAtMinDepth= 1;
  for (i=0; i<strlen(dimstr); i++) {
    sprintf(skipstring,"skip.%c",dimstr[i]);
    sprintf(extentstring,"d%c",dimstr[i]);
    sizeAtMinDepth *= kvGetInt(info,extentstring);
    if (kvLookup(info,skipstring)!=NULL) break;
  }
  if (i==strlen(dimstr)) minDepth= i-1;
  else minDepth= i;

  /* Create the requested chunk */
  mri_create_chunk( ds, chunk );
  if (strlen(kvGetString(info,"chunkfile"))) {
    sprintf(buf,"%.200s.file",chunk);
    mri_set_string( ds, buf, kvGetString(info,"chunkfile") );
  }
  else {
    /* Empty chunk file means default behavior, which is to put the
     * chunk data in the .mri file.
     */
  }
  sprintf(buf,"%.200s.datatype",chunk);
  switch (kvGetInt(info,"datatype_out")) {
  case SRDR_UINT8: mri_datatype_out= MRI_UINT8; break;
  case SRDR_INT16: mri_datatype_out= MRI_INT16; break;
  case SRDR_INT32: mri_datatype_out= MRI_INT32; break;
  case SRDR_FLOAT32: mri_datatype_out= MRI_FLOAT32; break;
  case SRDR_FLOAT64: mri_datatype_out= MRI_FLOAT64; break;
  case SRDR_INT64: mri_datatype_out= MRI_INT64; break;
  default:
    Abort("%s: internal error: Pgh MRI files have no equivalent to %s!\n",
	  progname, srdrTypeName[kvGetInt(info,"datatype_out")]);
  }
  mri_set_string( ds, buf, libmriTypeName[mri_datatype_out] );
  sprintf(buf,"%.200s.dimensions",chunk);
  mri_set_string( ds, buf, dimstr );
  for (i=0; i<strlen(dimstr); i++) {
    sprintf(extentstring,"d%c",dimstr[i]);
    sprintf(buf,"%.200s.extent.%c",chunk,dimstr[i]);
    mri_set_int(ds, buf, kvGetInt(info,extentstring));
  }
  
  /* And move the data. */
  in_offset= kvGetLong(info,"start_offset");
  out_offset= 0;
  recursiveTransfer(handler, info, ds, chunk, 
		    &in_offset, &out_offset,
		    strlen(dimstr)-1, minDepth, sizeAtMinDepth);
  
  /* Add whatever other tags seem appropriate */
  export_relevant_tags( ds, chunk, info );

}

static void parse_command_line( KVHash* info, char* readfile, char* hdrfile,
				int argc, char* argv[] )
{
  char string[512];
  long itmp;
  double dtmp;
  KVHash* defs= kvGetHash(info,"definitions");

  cl_scan( argc, argv );

  /* Get filenames */
  cl_get( "input|i", "%option %s[%]", "input.mri", readfile );
  cl_get( "output|out", "%option %s[%]", "output.mri", hdrfile );

  if (cl_get( "tag", "%option %s[%]", "", string )) {
    kvDefString(info,"tag",string);
    kvDefString(defs,"tag","an ID string associated with this data");
  }
  if (cl_get( "phaseref","%option %s", string )) {
    kvDefString(info,"phasereffile",string);
    kvDefString(defs,"phasereffile","filename for phase reference info");
  }
  if (cl_get( "bandpass","%option %s", string )) {
    kvDefString(info,"bandpassdir",string);
    kvDefString(defs,"bandpassdir","directory for scanner band pass info");
  }
  if (cl_get( "rampfile","%option %s", string )) {
    kvDefString(info,"rampfile",string);
    kvDefString(defs,"rampfile","filename for ramp sampling info");
  }
  if (cl_get( "auxfile","%option %s", string )) {
    kvDefString(info,"auxfile",string);
    kvDefString(defs,"auxfile","auxiliary format-specific info file");
  }
  kvDefBoolean(info,"xchop",cl_present( "xchop" ));
  kvDefString(defs,"xchop","phase chop in X?");
  kvDefBoolean(info,"ychop",cl_present( "ychop" ));
  kvDefString(defs,"ychop","phase chop in Y?");
  if (cl_present("autoscale"))
    kvDefBoolean(info,"cl_autoscale",1);
  kvDefString(defs,"autoscale","perform autoscaling?");
  if (cl_get("autoscale_range","%option %lf",&dtmp))
    kvDefDouble(info,"cl_autoscale_range",dtmp);

  debug= cl_present( "debug" );
  verbose_flg= cl_present( "verbose" );

  /* Let the user pick endianness if desired */
  if (cl_present("bigendian")) 
    kvDefBoolean(info,"cl_big_endian_input",1);
  else if (cl_present("littleendian")) 
    kvDefBoolean(info,"cl_big_endian_input",0);

  kvDefBoolean(info,"ignoreheader",cl_present( "ignoreheader" ));
  kvDefString(defs,"ignoreheader","ignore info in the file header?");

  kvDefBoolean(info,"multi",cl_present( "multi" ));
  kvDefString(defs,"multi","read multiple files");

  /* Get data-type, vector length, and dimension lengths */
  if (cl_get( "type|t", "%options %s[%]", "SRDR_INT16", string )) {
    /* Convert data-type string to number */
    if( !strcmp( string, "MRI_SHORT" ) ||
	!strcmp( string, "SRDR_INT16" ) ||
	!strcasecmp( string, "short" ) ||
	!strcasecmp( string, "s" ) )
      kvDefInt(info,"datatype_in",SRDR_INT16);
    else if( !strcmp( string, "SRDR_UINT16" ) ||
	!strcasecmp( string, "ushort" ) ||
	!strcasecmp( string, "us" ) )
      kvDefInt(info,"datatype_in",SRDR_UINT16);
    else if( !strcmp( string, "MRI_FLOAT" ) ||
	     !strcmp( string, "SRDR_FLOAT32" ) ||
	     !strcasecmp( string, "float" ) ||
	     !strcasecmp( string, "f" ) )
      kvDefInt(info,"datatype_in",SRDR_FLOAT32);
    else if( !strcmp( string, "MRI_UNSIGNED_CHAR" ) ||
	     !strcmp( string, "SRDR_UINT8" ) ||
	     !strcasecmp( string, "uchar" ) ||
	     !strcasecmp( string, "u" ) ||
	     !strcasecmp( string, "c" ) )
      kvDefInt(info,"datatype_in",SRDR_UINT8);
    else if( !strcmp( string, "MRI_INT" ) ||
	     !strcmp( string, "SRDR_INT32" ) ||
	     !strcasecmp( string, "long" ) ||
	     !strcasecmp( string, "l" ) ||
	     !strcasecmp( string, "int" ) ||
	     !strcasecmp( string, "i" ) )
      kvDefInt(info,"datatype_in",SRDR_INT32);
    else if( !strcmp( string, "MRI_DOUBLE" ) ||
	     !strcmp( string, "SRDR_FLOAT64" ) ||
	     !strcasecmp( string, "double" ) ||
	     !strcasecmp( string, "d" ) )
      kvDefInt(info,"datatype_in",SRDR_FLOAT64);
    else 
      Abort( "Data type unrecognized: %s.", string );
  }

  /* Get byte offsets parameters */
  if (cl_get( "offset|o", "%option %ld[%]", 0, &itmp ))
    kvDefLong(info,"start_offset",itmp);
  if (cl_get( "skip|s", "%option %ld[%]", 0, &itmp ))
    kvDefLong(info,"cl_skip",itmp);
  if (cl_get( "sliceskip|ss","%option %ld[%]", 0, &itmp ))
    kvDefLong(info,"cl_sliceskip",itmp);
  
  /* Get EPI-reordering switch */
  if (cl_get( "reorder|r", "%option %s", string )) {
    if ( isdigit( string[0] ) )
      kvDefBoolean(info,"cl_reorder",atol( string ));
    else
      kvDefBoolean(info,"cl_reorder",
		   (( ( string[0] == 'T' ) || ( string[0] == 't' ) ) 
		    ? 1: 0));
  }
  
  /* Get data order info */
  if (cl_get( "dimorder|dataorder|do", "%option %s", &string )) {
    if (!legitDimString(string))
      Abort("%s: <%s> is not a legitimate data order!\n",
	    progname,string);
    kvDefString(info,"cl_dim_string",string);
    kvDefString(defs,"cl_dim_string",
		"order of dimensions according to command line");
  }

  if (cl_get( "dims|di", "%option %s", &string )) {
    kvDefString(info,"cl_extent_string",string);
    kvDefString(defs,"cl_extent_string",
		"extents of dimensions according to command line");
  }

  if (cl_get( "skips", "%option %s", &string )) {
    kvDefString(info,"cl_skipstr",string);
    kvDefString(defs,"cl_skipstr",
		"skip lengths of dimensions according to command line");
  }

  /* Allow an arbitrary definition from the command line */
  while (cl_get( "def", "%option %s", &string )) {
    char* key= string;
    char* value= strchr(string,'=');
    if (value==NULL) {
      /* Treat this as a boolean to be set true */
      kvDefBoolean(info,key,1);
    }
    else {
      /* Break the string into a key and a value */
      *value= '\0';
      value++;
      
      if (*value=='\0')
	Abort("%s: empty definition value!\n",progname);

      /* If it's TRUE or FALSE (or any common equivalent), it's boolean.
       * If it consists only of digits, '+', and '-' it's an int.
       * If it consists only of digits, '.', and E or e, it's double.
       */
      if (!strcasecmp(value,"T") || !strcasecmp(value,"TRUE")) {
	kvDefBoolean(info,key,1);
      }
      else if (!strcasecmp(value,"F") || !strcasecmp(value,"FALSE")) {
	kvDefBoolean(info,key,0);
      }
      else if (strspn(value,"0123456789+-")==strlen(value)) {
	kvDefInt(info,key,atoi(value));
      }
      else if (strspn(value,"0123456789+-.eE")==strlen(value)) {
	kvDefDouble(info,key,atof(value));
      }
      else kvDefString(info,key,value);
    }
  }
  
  if (cl_cleanup_check()) {
    int i;
    fprintf(stderr,"%s: invalid argument in command line:\n    ",argv[0]);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(-1);
  }

  kvDefInt(info,"handler_datatype_out", kvGetInt(info,"datatype_in"));
}

static void initStuff(KVHash* info)
{
  initInfoHash(info);

  /* Set up some default values specific to the first chunk */
  kvDefBoolean(info,"big_endian_input",1); /* fmri data usually bigendian */
  kvDefInt(info,"datatype_in",SRDR_INT16); /* by far the most common */
#ifdef never
  kvDefInt(info,"datatype_out",SRDR_FLOAT32);
#endif
  kvDefDouble(info,"autoscale_range",1.0);
  kvDefInt(info,"start_offset",0); 
  kvDefString(info,"chunkname","images");
  kvDefString(info,"chunkfile",".dat");
}

int main(int argc, char* argv[])
{
  char readfile[512], hdrfile[512];
  KVHash* info;
  FileHandlerFactory handlerFactory= NULL;
  MRI_Dataset *outDS = NULL;
  SList* chunkStack= NULL;
  ChunkHandlerPair* chPair= NULL;

  /* Print version number */
  Message( "# %s\n", rcsid );

  progname= argv[0];

  /* Check to see if help was requested */
  if( ( argc > 1 ) && !strcmp( argv[1], "-help" ) )
    {
      if( argc == 2 )
	Help( "selecttopic" );
      else
	Help( argv[2] );
    }

  /* This will serve as a list of things to do */
  chunkStack= slist_create();

  /* (almost) All knowledge will accumulate as key-value pairs */
  info= kvFactory(KV_DEFAULT_SIZE);

  /* Miscellaneous initialization tasks- definition and translation
   * tables, default values, opportunity for libbio to figure out
   * endian-ness of this machine.
   */
  InitBIO();
  initStuff(info);

  /* Parse command line */
  parse_command_line(info, readfile, hdrfile, argc, argv);
  bio_big_endian_input = kvGetBoolean(info,"big_endian_input");

  /* Let's find a handler type that can deal with this file.  If it's
   * a single filename, no problem.  If "multi" is set it may
   * contain wildcards, so we'll pick the first matching file
   * as our prototype for testing.
   */
  if (kvGetBoolean(info,"ignoreheader")) {
    if (kvGetInt(info,"datatype_in")==SRDR_UINT16)
      handlerFactory= ushortFactory;
    else
      handlerFactory= rawFactory;
  }
  else if (kvGetBoolean(info,"multi")) {
    int i;
    char* matchingName= expandWildcardFirstFilename(readfile);
    
    if (!matchingName) 
      Abort("%s: no files with names matching <%s> were found!\n",
	    progname, readfile);
    
    for (i=0; i<sizeof(handlerTable)/sizeof(FileHandlerPair); i++) {
      if ( (*handlerTable[i].tester)(matchingName) ) {
	if (verbose_flg) Message("It's <%s>!\n",handlerTable[i].name);
	handlerFactory= handlerTable[i].factory;
	break;
      }
      else {
	if (verbose_flg) Message("It's not <%s>\n",handlerTable[i].name);
      }
    }
    
    free(matchingName);
  }
  else {
    int i;
    struct stat buf;
    
    if (stat(readfile, &buf)) {
      Abort("%s: could not stat <%s>: %s!\n",
	    progname, readfile, strerror(errno));
    }
    for (i=0; i<sizeof(handlerTable)/sizeof(FileHandlerPair); i++) {
      if ( (*handlerTable[i].tester)(readfile) ) {
	if (verbose_flg) Message("It's <%s>!\n",handlerTable[i].name);
	handlerFactory= handlerTable[i].factory;
	break;
      }
      else {
	if (verbose_flg) Message("It's not <%s>\n",handlerTable[i].name);
      }
    }
  }
  if (handlerFactory==NULL) 
    Abort("%s: could not find a handler type that recognized %s!\n",
	  progname, readfile);

  /* These will be the first pair on our stack of things to do. 
   * If there may be multiple files, we put a wrapper around the
   * handlers we've found to deal with them.
   */
  if (!(chPair= (ChunkHandlerPair*)malloc(sizeof(ChunkHandlerPair))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,sizeof(ChunkHandlerPair));
  if (kvGetBoolean(info,"multi")) {
    SList* fileList= expandWildcardFilename(readfile);
    chPair->handler= multiFileHandlerFactory();
    while (!slist_atend(fileList)) {
      char* fname= (char*)slist_next(fileList);
      multiFileHandlerAddFile(chPair->handler,
			      handlerFactory(fname, info));
    }
    destroyWildcardFilenameList(fileList);
  }
  else {
    chPair->handler= (*handlerFactory)(readfile, info);
  }
  chPair->info= info;
  slist_append(chunkStack,chPair);

  /* We can now loop over the chunks (perhaps discovering more
   * chunks in the process).
   */
  while (!slist_empty(chunkStack)) {
    chPair= (ChunkHandlerPair*)slist_pop(chunkStack);

    FH_PROCESSHEADER(chPair->handler, chPair->info, chunkStack);

    if (!smart_reconcile(chPair->info)){
      smart_dump(chPair->info);
      Abort("Available information is irreconcilable; can't proceed!\n");
    }

    if (!consistency_check(chPair->info)) {
      Error("Consistency test failed!\n");
      smart_dump(chPair->info);
      Abort("%s: File format information is missing or inconsistent; can't proceed!\n",
	    progname);
    }

    /* If we have to do type conversion, we'll put a wrapper around the
     * FileHandler to do so.
     */

    if (kvGetInt(chPair->info,"handler_datatype_out") 
	== kvGetInt(chPair->info,"datatype_out")) {
      /* No action required */
    }
    else {
      if (debug) fprintf(stderr,"Chunk %s needs a converter (%s to %s)\n",
			 kvGetString(chPair->info,"chunkname"),
			 srdrTypeName[kvGetInt(chPair->info,
					   "handler_datatype_out")],
			 srdrTypeName[kvGetInt(chPair->info,
					   "datatype_out")]);
      chPair->handler= 
	convertFactory( chPair->handler, 
			kvGetInt(chPair->info,"handler_datatype_out"),
			kvGetInt(chPair->info,"datatype_out") );
    }
    
    if (verbose_flg) smart_dump(chPair->info);

    if (!outDS) {
      /* This must be the first pass; good time to show some info */
      emit_format_summary(chPair->info, chPair->handler);

      /* Open output dataset.  This will cause libbio to decide what type of
       * machine it's living on (via libmri's call to InitBIO).  We'll then
       * update that information based on command line info.
       */
      outDS= open_output( hdrfile );
      bio_big_endian_input = kvGetBoolean(chPair->info,"big_endian_input");
    }

    transfer( chPair->handler, chPair->info, outDS );

    FH_DESTROYSELF(chPair->handler);
    kvDestroy(chPair->info);
    free(chPair);
  } 

  /* Add this command to the dataset history.  We do this last because
   * history information may have been added by various chunk handlers
   * as well.
   */
  hist_add_cl(outDS, argc, argv);

  /* Close output dataset */
  mri_close_dataset( outDS );

  /* Clean up */
  slist_destroy(chunkStack,NULL);

  if (verbose_flg)
    fprintf(stderr,"#      Data converted to standard format.\n");
  return 0;

}

