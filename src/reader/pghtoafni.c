/************************************************************
 *                                                          *
 *  pghtoafni.c                                             *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 2003 Department of Statistics,         *
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
 *  Derived from smartreader, Joel Welling 6/03
 ************************************************************/

/* Notes-
 */

#include <stdlib.h>
#include <stdarg.h>
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
#include <assert.h>

#include "bio.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"
#include "smartreader.h"
#include "afni_defs.h"

static char rcsid[] = "$Id: pghtoafni.c,v 1.20 2007/07/03 20:04:55 welling Exp $";

#define SINGLE_QUOTE 0x27

/* The following tag name translations occur when internal tags are written
 * to the AFNI header file.  Other name translations may occur as well;
 * see translateNameToAfni() below.  Mapping a name to the empty string
 * will guarantee that it will not be written; specifying a non-empty
 * translation string will guarantee that it will be written.
 */
static char* afniNameTransTable[][2]= {
  {"slice_gap","pghmri.slice_gap"},
  {"slice_thickness","pghmri.slice_thickness"},
  {"afni.IDCODE_STRING", ""}, /* force AFNI to generate a new one */
  {"tlf.0","pghmri.tlf"},
  {"tlf.1","pghmri.tlf"},
  {"tlf.2","pghmri.tlf"},
  {"trf.0","pghmri.trf"},
  {"trf.1","pghmri.trf"},
  {"trf.2","pghmri.trf"},
  {"blf.0","pghmri.blf"},
  {"blf.1","pghmri.blf"},
  {"blf.2","pghmri.blf"},
  {"brf.0","pghmri.brf"},
  {"brf.1","pghmri.brf"},
  {"brf.2","pghmri.brf"},
  {"tlb.0","pghmri.tlb"},
  {"tlb.1","pghmri.tlb"},
  {"tlb.2","pghmri.tlb"},
  {"trb.0","pghmri.trb"},
  {"trb.1","pghmri.trb"},
  {"trb.2","pghmri.trb"},
  {"blb.0","pghmri.blb"},
  {"blb.1","pghmri.blb"},
  {"blb.2","pghmri.blb"},
  {"brb.0","pghmri.brb"},
  {"brb.1","pghmri.brb"},
  {"brb.2","pghmri.brb"},
  {"ctr.0","pghmri.ctr"},
  {"ctr.1","pghmri.ctr"},
  {"ctr.2","pghmri.ctr"},
  {"fov.x","pghmri.fov.x"},
  {"fov.y","pghmri.fov.y"},
  {"fov.z","pghmri.fov.z"},
  {"time","pghmri.time"},
  {"TE","pghmri.TE"},
  {"missing","pghmri.missing"},
  {"missing.dimstr","pghmri.missing.dimensions"},
  {"brick_max",""}, /* If set, these are stored in BRICK_STATS */
  {"brick_min",""}, /* If set, these are stored in BRICK_STATS */
  {NULL, NULL} /* ends list */
};

/* We need to place a limit on the maximum number of "things" in one IO op */
#define MAX_BLOCK (16*1024*1024)

int debug = 0;          /* Global debug flag                         */

int verbose_flg = 0;        /* Global verbosity value (0=off) */

char* progname= NULL; /* program name */

/* This table maps KVTypes to Afni type names */
static char* afniTypeName[] = {
  "string-attribute", 
  "integer-attribute",
  "float-attribute",
  NULL,
  "integer-attribute",
  "integer-attribute"
};
  
static void initStuff(KVHash* info)
{
  initInfoHash(info);

  /* Set up some default values specific to the first chunk */
  kvDefBoolean(info,"big_endian_input",1); /* fmri data usually bigendian */
  kvDefInt(info,"datatype_in",SRDR_INT16);
#ifdef never
  kvDefInt(info,"datatype_out",SRDR_FLOAT32);
#endif
  kvDefDouble(info,"autoscale_range",1.0);
  kvDefInt(info,"start_offset",0); 
  kvDefString(info,"chunkname","images");
  kvDefString(info,"chunkfile",".dat");
}

static void parse_command_line( KVHash* info, 
				char* inputFileName, char* outFileName,
				int argc, char* argv[] )
{
  char string[512];
  int itmp;
  double dtmp;
  KVHash* defs= kvGetHash(info,"definitions");

  cl_scan( argc, argv );

  /* Get filenames */
  if (!cl_get( "", "%s", inputFileName )) {
    fprintf(stderr,"%s: Input file name not given!\n",progname);
    Help("usage");
    exit(-1);
  }
  if (!cl_get( "", "%s", outFileName )) {
    fprintf(stderr,"%s: Output file name not given!\n",progname);
    Help("usage");
    exit(-1);
  }

  debug= cl_present( "debug" );
  verbose_flg= cl_present( "verbose|v" );
  if (cl_present("anat")) kvDefBoolean(info,"cl_hint_anat",1);
  if (cl_present("func")) kvDefBoolean(info,"cl_hint_func",1);
  if (cl_present("head")) kvDefBoolean(info,"cl_hint_head",1);
  if (cl_present("gen")) kvDefBoolean(info,"cl_hint_gen",1);

  kvDefString(info,"cl_input_file_name",inputFileName);

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

static int canTranslate(SList* chunkList)
{
  int found= 0;
  ChunkHandlerPair* imgPair= NULL;
  KVHash* info= NULL;
  const char* dimstr;

  slist_totop(chunkList);
  while (!slist_atend(chunkList)) {
    imgPair= (ChunkHandlerPair*)slist_next(chunkList);
    if (kvLookup(imgPair->info,"chunkname") 
	&& !strcmp(kvGetString(imgPair->info,"chunkname"),"images")) {
      found= 1;
      break;
    }
  }
  slist_totop(chunkList);

  if (!found) {
    Error("%s: input file has no 'images' chunk!\n",progname);
    return 0;
  }

  info= imgPair->info;
  dimstr= kvGetString(info,"dimstr");
  if (strcmp(dimstr,"vxyz") && strcmp(dimstr,"vxyzt")
      && strcmp(dimstr,"xyz") && strcmp(dimstr,"xyzt")
      && strcmp(dimstr,"xyzb") && strcmp(dimstr,"vxyzb")) {
    Error("%s: input file image dimensions <%s> do not translate!\n",
	  progname, dimstr);
    return 0;
  }

  if (*dimstr == 'v') {
    int dv= kvGetInt(info,"dv");
    if (dv<1 || dv>3) {
      Error("%s: input file image v dimension of %d does not translate!\n",
	    progname,dv);
      return 0;
    }
  }

  return 1;
}

static void openOutputFiles(const char* fname, FILE** head, FILE** brik)
{
  char* buf;
  char* here;

  if (!(buf=(char*)malloc(strlen(fname)+6)))
    Abort("%s: unable to allocate %d bytes!\n",progname,strlen(fname)+6);
  strcpy(buf,fname);

  here= strrchr(buf,'.');
  if (here && (!strcmp(here,".HEAD") || !strcmp(here,".BRIK"))) {
    *here= '\0';
  }
  else here= buf + strlen(buf);
      
  strcat(buf,".HEAD");
  if (!(*head=fopen(buf,"w")))
    Abort("%s: cannot open <%s> for writing!\n",progname,buf);
  *here= '\0';
  strcat(buf,".BRIK");
  if (!(*brik=fopen(buf,"w")))
    Abort("%s: cannot open <%s> for writing!\n",progname,buf);

  free(buf);
}

static void closeOutputFiles(FILE* head, FILE* brik)
{
  (void)fclose(head);
  (void)fclose(brik);
}

static void writeOneTag( KVHash* info, FILE* f, KVType type, const char* name,
			 int count, void* buf )
{
  KVHash* wroteThese= kvGetHash(info,"wrote_these_already");

  if (kvGetBoolean(wroteThese,name)) {
    if (debug) fprintf(stderr,"Supressing duplicate tag <%s>\n",name);
    return;
  }
  kvDefBoolean(wroteThese,name,1);

  fprintf(f,"type = %s\n",afniTypeName[(int)type]);
  fprintf(f,"name = %s\n",name);
  fprintf(f,"count = %d\n",count);

  switch (type) {
  case KV_STRING: 
    {
      char* tbuf= (char*)buf;
      int i;
      fputc(SINGLE_QUOTE,f);
      for (i=0; i<count-1; i++) {
#ifdef never
	if (tbuf[i]=='~') fputc('*',f);
	else if (tbuf[i]=='\0') fputc('~',f);
	else fputc(tbuf[i],f);
#endif
	fputc(tbuf[i],f);
      }
      fputc('~',f);
      fputc('\n',f);
    }
   break;
  case KV_LONG:
    {
      long long* tbuf= (long long *)buf;
      int i;
      for (i=0; i<count; i++) {
	fprintf(f,"%lld ",tbuf[i]);
	if (!((i+1)%5)) fputc('\n',f);
      }
      if (count%5) fputc('\n',f);
    }
    break;
  case KV_DOUBLE:
    {
      double* tbuf= (double*)buf;
      int i;
      for (i=0; i<count; i++) {
	fprintf(f,"%lf ",tbuf[i]);
	if (!((i+1)%5)) fputc('\n',f);
      }
      if (count%5) fputc('\n',f);
    }
    break;
  case KV_BOOLEAN:
    {
      int* tbuf= (int*)buf;
      int i;
      for (i=0; i<count; i++) {
	fprintf(f,"%d ",tbuf[i]);
	if (!((i+1)%5)) fputc('\n',f);
      }
      if (count%5) fputc('\n',f);
    }
  case KV_INT:
    {
      int* tbuf= (int*)buf;
      int i;
      for (i=0; i<count; i++) {
	fprintf(f,"%d ",tbuf[i]);
	if (!((i+1)%5)) fputc('\n',f);
      }
      if (count%5) fputc('\n',f);
    }
    break;
  default:
    Abort("%s: internal error: nonsense type in writeOneTag!\n",progname);
  }

  fprintf(f,"\n");
}

static void writeIntTag( KVHash* info, FILE* afniHead, const char* name, 
			 int count, ... )
{
  va_list ap;
  int i;
  static int* buf= NULL;
  int bufSize= 0;

  /* Grow the buffer as needed */
  if (bufSize<count) {
    if (buf) free(buf);
    if (!(buf=(int*)malloc(count*sizeof(int))))
      Abort("%s: unable to allocate %d bytes!\n",
	    progname,count*sizeof(int));
    bufSize= count;
  }

  va_start(ap, count);
  for (i=0; i<count; i++) {
    buf[i]= va_arg(ap,int);
  }
  va_end(ap);
  writeOneTag( info, afniHead, KV_INT, name, count, (void*)buf );
}

static void writeDoubleTag( KVHash* info, FILE* afniHead, const char* name, 
			    int count, ... )
{
  va_list ap;
  int i;
  static double* buf= NULL;
  int bufSize= 0;

  /* Grow the buffer as needed */
  if (bufSize<count) {
    if (buf) free(buf);
    if (!(buf=(double*)malloc(count*sizeof(double))))
      Abort("%s: unable to allocate %d bytes!\n",
	    progname,count*sizeof(double));
    bufSize= count;
  }

  va_start(ap, count);
  for (i=0; i<count; i++) buf[i]= va_arg(ap,double);
  va_end(ap);
  writeOneTag( info, afniHead, KV_DOUBLE, name, count, (void*)buf );
}

static void writeStringTag( KVHash* info, FILE* afniHead, const char* name, 
			    const char* val )
{
  writeOneTag( info, afniHead, KV_STRING, name, strlen(val)+1, (void*)val );
}

static int pickAfniType(KVHash* info)
{
  const char* type1= NULL;
  const char* type2= NULL;
  int head_known= 0;
  int is_head= 0;
  int func_known= 0;
  int is_func= 0;

  /* If the necessary info is available from the command line,
   * we'll let that override.
   */
  if (kvGetBoolean(info,"cl_hint_gen") 
      && !kvGetBoolean(info,"cl_hint_head")) {
    is_head= 0;
    head_known= 1;
  }
  if (kvGetBoolean(info,"cl_hint_head") 
      && !kvGetBoolean(info,"cl_hint_gen")) {
    is_head= 1;
    head_known= 1;
  }

  if (kvGetBoolean(info,"cl_hint_func") 
      && !kvGetBoolean(info,"cl_hint_anat")) {
    is_func= 1;
    func_known= 1;
  }
  if (kvGetBoolean(info,"cl_hint_anat") 
      && !kvGetBoolean(info,"cl_hint_func")) {
    is_func= 0;
    func_known= 1;
  }

  if (!func_known || !head_known) {
    /* Check to see if the input file used to be an AFNI file
     * with the necessary information.
     */
    if (kvLookup(info,"afni.TYPESTRING")) 
      type1= kvGetString(info,"afni.TYPESTRING");
    if (kvLookup(info,"afni.SCENE_DATA.2")) 
      type2= kvGetString(info,"afni.SCENE_DATA.2");
    if (type1) {
      if (type2) {
	if (strcmp(type1,type2))
	  Abort("%s: input file has inconsistent AFNI type string info!\n",
		progname);
      }
    }
    else if (type2) type1= type2;
    if (type1) {
      int oldtype= stringTableLookup(type1, typestring_names);
      if (oldtype<0) 
	Abort("%s: input file has unknown AFNI type string <%s>!\n",
	      progname, type1);
      /* Maybe the command line specified a change to one feature */
      if (func_known) {
	if (is_func) {
	  switch (oldtype) {
	  case THREEDIM_HEAD_ANAT:
	  case THREEDIM_HEAD_FUNC:
	    return THREEDIM_HEAD_FUNC;
	    break;
	  case THREEDIM_GEN_ANAT:
	  case THREEDIM_GEN_FUNC:
	    return THREEDIM_GEN_FUNC;
	    break;
	  }
	}
	else {
	  switch (oldtype) {
	  case THREEDIM_HEAD_ANAT:
	  case THREEDIM_HEAD_FUNC:
	    return THREEDIM_HEAD_ANAT;
	    break;
	  case THREEDIM_GEN_ANAT:
	  case THREEDIM_GEN_FUNC:
	    return THREEDIM_GEN_ANAT;
	    break;
	  }
	}
      }
      else if (head_known) {
	if (is_head) {
	  switch (oldtype) {
	  case THREEDIM_HEAD_ANAT:
	  case THREEDIM_GEN_ANAT:
	    return THREEDIM_HEAD_ANAT;
	    break;
	  case THREEDIM_HEAD_FUNC:
	  case THREEDIM_GEN_FUNC:
	    return THREEDIM_HEAD_FUNC;
	    break;
	  }
	}
	else {
	  switch (oldtype) {
	  case THREEDIM_HEAD_ANAT:
	  case THREEDIM_GEN_ANAT:
	    return THREEDIM_GEN_ANAT;
	    break;
	  case THREEDIM_HEAD_FUNC:
	  case THREEDIM_GEN_FUNC:
	    return THREEDIM_GEN_FUNC;
	    break;
	  }
	}
      }
      else return oldtype;
    }
  }

  /* Make a guess. We'll consider the input file name if that helps.
   */
  if (!func_known && kvLookup(info,"cl_input_file_name")) {
    char* fname= kvGetString(info,"cl_input_file_name");
    if (!strncasecmp(fname,"mean",strlen("mean"))
	|| !strncasecmp(fname,"grandmean",strlen("grandmean"))
	|| !strncasecmp(fname,"struct",strlen("struct"))
	|| !strncasecmp(fname,"strct",strlen("strct"))) {
      is_func= 0;
    }
    else {
      is_func= 1;
    }
    func_known= 1;
  }

  if (!head_known) {
    is_head= 1;
    head_known= 1;
  }

  if (head_known && func_known) {
    if (is_head) {
      if (is_func) return THREEDIM_HEAD_FUNC;
      else return THREEDIM_HEAD_ANAT;
    }
    else {
      if (is_func) return THREEDIM_GEN_FUNC;
      else return THREEDIM_GEN_ANAT;
    }
  }
  else
    return THREEDIM_HEAD_FUNC;
}

static int pickAfniViewType( KVHash* info )
{
  if (kvLookup(info,"afni.SCENE_DATA.0")) {
    const char* s= kvGetString(info,"afni.SCENE_DATA.0");
    int result= stringTableLookup(s, coordsys_names);
    if (result<0) 
      Abort("%s: input file has unknown AFNI coord sys string <%s>!\n",
	    progname,s);
    return result;
  }
  else if (kvLookup(info,"voxel_x") && kvLookup(info,"voxel_y")
	   && kvLookup(info,"voxel_x")) {
    if ((kvGetDouble(info,"voxel_x")==1.0)
	&& (kvGetDouble(info,"voxel_y")==1.0)
	&& (kvGetDouble(info,"voxel_z")==1.0)) {
      if (kvLookup(info,"afni.WARP_DATA.0"))
	return 2; /* "+tlrc" */
      else return 0; /* "+orig" */
    }
    else return 0; /* "+orig" */
  }
  else return 0; /* "+orig" */
}

static int pickAfniFuncType( KVHash* info, const int afniType )
{
  int is_bucket= 0;
  const char* dimstr= kvGetString(info,"dimstr");

  if (dimstr[strlen(dimstr)-1]=='b') {
    if (kvGetInt(info,"db")>1) is_bucket= 1;
    else is_bucket= 0;
  }
  else is_bucket= 0;

  if (afniType==THREEDIM_HEAD_ANAT || afniType==THREEDIM_GEN_ANAT) {
    if (kvLookup(info,"afni.SCENE_DATA.1")) {
      const char* s= kvGetString(info,"afni.SCENE_DATA.1");
      int result= stringTableLookup(s, anat_type_names);
      if (result<0) {
	result= stringTableLookup(s, func_type_names);
	if (result>=0) {
	  Warning(1,"%s: Dataset declared anatomical, used to be functional!\n",
		  progname);
	  if (is_bucket) return ANAT_BUCK_TYPE;
	  else return ANAT_OMRI_TYPE;
	}
	else
	  Abort("%s: input file has unknown AFNI coord sys string <%s>!\n",
		progname,s);
      }
      else return result;
    }
    else {
      if (is_bucket) return ANAT_BUCK_TYPE;
      else return ANAT_OMRI_TYPE;
    }
  }
  else {
    if (kvLookup(info,"afni.SCENE_DATA.1")) {
      const char* s= kvGetString(info,"afni.SCENE_DATA.1");
      int result= stringTableLookup(s, func_type_names);
      if (result<0) {
	result= stringTableLookup(s, anat_type_names);
	if (result>=0) {
	  Warning(1,"%s: Dataset declared functional, used to be anatomical!\n",
		  progname);
	  if (is_bucket) return FUNC_BUCK_TYPE;
	  else return FUNC_FIM_TYPE;
	}
	else
	  Abort("%s: input file has unknown AFNI coord sys string <%s>!\n",
		progname,s);
      }
      else return result;
    }
    else {
      if (is_bucket) return FUNC_BUCK_TYPE;
      else return FUNC_FIM_TYPE;
    }
  }
  return 0; /* NOTREACHED */
}

static int lookupOrientation( const char* s )
{
  int result= stringTableLookup(s, orientation_names);
  if (result<0) 
    Abort("%s: input file has unknown AFNI orientation code <%s>!\n",
	  progname,s);
  return result;
}

static void calcSliceTimeOffsets( KVHash* info, FILE* afniHead )
{
  double TR;
  double TS;
  int dz;
  double* vals;
  int i;

  if (kvLookup(info,"TR")) TR= (double)kvGetInt(info,"TR");
  else TR= 0.0;

  if (kvLookup(info,"dz")) dz= kvGetInt(info,"dz");
  else dz= 1;

  TS= TR/(1000*dz); /* need to convert from usec to msec */

  if (!(vals= (double*)malloc(dz*sizeof(double))))
    Abort("%s: unable to allocate %d bytes!\n", progname, dz*sizeof(double));

  if (kvLookup(info,"reorder_pattern")) {
    const char* s= kvGetString(info,"reorder_pattern");
    if (!strcasecmp(s,"sequential")) {
      for (i=0; i<dz; i++) vals[i]= i*TS;
    }
    else if (!strcasecmp(s,"reversed_sequential")) {
      for (i=0; i<dz; i++) vals[i]= (dz-(i+1))*TS;
    }
    else if (!strcasecmp(s,"even/odd")) {
      double v= 0.0;
      for (i=0; i<dz; i+=2) {
	vals[i]= v;
	v += TS;
      }
      for (i=1; i<dz; i+=2) {
	vals[i]= v;
	v += TS;
      }
    }
    else if (!strcasecmp(s,"odd/even")) {
      double v= 0.0;
      for (i=1; i<dz; i+=2) {
	vals[i]= v;
	v += TS;
      }
      for (i=0; i<dz; i+=2) {
	vals[i]= v;
	v += TS;
      }
    }
    else if (!strcasecmp(s,"reversed_even/odd")) {
      double v= 0.0;
      int maxEven;
      int maxOdd;
      if (dz % 2) { /* odd number of slices */
	maxEven= dz-1;
	maxOdd= dz-2;
      }
      else { /* even number of slices */
	maxEven= dz-2;
	maxOdd= dz-1;
      }
      for (i=maxEven; i>=0; i -= 2) {
	vals[i]= v;
	v += TS;
      }
      for (i=maxOdd; i>=0; i -= 2) {
	vals[i]= v;
	v += TS;
      }
    }
    else if (!strcasecmp(s,"reversed_odd/even")) {
      double v= 0.0;
      int maxEven;
      int maxOdd;
      if (dz % 2) { /* odd number of slices */
	maxEven= dz-1;
	maxOdd= dz-2;
      }
      else { /* even number of slices */
	maxEven= dz-2;
	maxOdd= dz-1;
      }
      for (i=maxOdd; i>=0; i -= 2) {
	vals[i]= v;
	v += TS;
      }
      for (i=maxEven; i>=0; i -= 2) {
	vals[i]= v;
	v += TS;
      }
    }
    else Abort("%s: input file has unknown reorder pattern <%s>!\n",
	       progname, s);
  }
  else {
    /* Assume even/odd reorder, the most common case */
    double v= 0.0;
    for (i=0; i<dz; i+=2) {
      vals[i]= v;
      v += TS;
    }
    for (i=1; i<dz; i+=2) {
      vals[i]= v;
      v += TS;
    }
  }

  writeOneTag(info, afniHead, KV_DOUBLE, "TAXIS_OFFSETS", dz, vals);

  free(vals);
}

static SRDR_Datatype closestTypeToAfni( const SRDR_Datatype type, 
					const int dv )
{
  if (dv==2) {
    return SRDR_FLOAT32; /* we're going for complex */
  }
  else if (dv==3) {
    return SRDR_UINT8; /* my guess for type of rgb data */
  }
  else switch (type) {
  case SRDR_INT16: return SRDR_INT16;
  case SRDR_UINT16: return SRDR_INT16;
  case SRDR_FLOAT64: return SRDR_FLOAT32;
  case SRDR_UINT8: return SRDR_UINT8;
  case SRDR_INT32: return SRDR_INT16;
  case SRDR_FLOAT32: return SRDR_FLOAT32;
  case SRDR_INT64: return SRDR_INT16;
  default: Abort("%s: internal error: unknown type %d in closestAfniType!\n",
		 progname, type);
  }
  return 0; /* not reached */
}

static int afniTypeNum( const SRDR_Datatype type, const int dv )
{
  if (dv==2) return 5; /* complex */
  else if (dv==3) return 6; /* rgb */
  else switch (type) {
  case SRDR_UINT8: return 0; /* byte */
  case SRDR_INT16: return 1; /* short */
  case SRDR_FLOAT32: return 3; /* float */
  default: Abort("%s: unknown type %d in afniTypeNum!\n",
		 progname,type);
  }
  return 0; /* NOTREACHED */
}

static void writeBrickTypes(KVHash* info, FILE* afniHead)
{
  int dt;
  int db;
  int dv;
  int* vals= NULL;
  int type;
  int nBricks;
  int i;

  if (kvLookup(info,"dt")) dt= kvGetInt(info,"dt");
  else dt= 1;
  
  if (kvLookup(info,"db")) db= kvGetInt(info,"db");
  else db= 1;
  
  if (kvLookup(info,"dv")) dv= kvGetInt(info,"dv");
  else dv= 1;

  if (dt==1) {
    if (db==1) nBricks= 1;
    else nBricks= db;
  }
  else {
    if (db != 1) 
      Abort("%s: internal error: multiple buckets at multiple times!\n",
	    progname);
    nBricks= dt;
  }
  
  if (!(vals=(int*)malloc(nBricks*sizeof(int)))) 
    Abort("%s: unable to allocate %d bytes!\n", progname, nBricks*sizeof(int));
  
  type= afniTypeNum( closestTypeToAfni(kvGetInt(info,"handler_datatype_out"),
				       dv ),
		     dv );
  
  for (i=0; i<nBricks; i++) vals[i]= type;
  writeOneTag(info, afniHead, KV_INT, "BRICK_TYPES", nBricks, vals);
  free(vals);
}

static void exportRequiredTags( KVHash* info, FILE* afniHead )
{
  /* This routine is in charge of exporting the info that *all* AFNI
   * files must have.  It bypasses the name translation mechanisms
   * used elsewhere.
   */

  KVHash* defs= kvGetHash(info,"definitions");
  KVHash* extNames= kvGetHash(info,"external_names");
  int dx;
  int dy;
  int dz;
  int dt;
  int dv;
  int db;
  int afniType= pickAfniType(info);
  int afniViewType= pickAfniViewType(info);
  int afniFuncType= pickAfniFuncType(info, afniType);
  double blb[3];

  /* Only a limited range of dimension strings are possible here,
   * due to checking when we made sure this file was translatable.
   */
  if (kvLookup(info,"dv")) dv= kvGetInt(info,"dv");
  else dv= 1;
  dx= kvGetInt(info,"dx");
  dy= kvGetInt(info,"dy");
  if (kvLookup(info,"dz")) dz= kvGetInt(info,"dz");
  else dz= 1;
  if (kvLookup(info,"dt")) dt= kvGetInt(info,"dt");
  else dt= 1;
  if (kvLookup(info,"db")) db= kvGetInt(info,"db");
  else db= 1;
  assert(dx>0);
  assert(dy>0);
  assert(dz>0);
  assert(dt>0);
  assert(dv>=1 && dv<=3);
  assert(db>0);
  if (db>1) {
    if (dt==1) writeIntTag(info,afniHead,"DATASET_RANK",2,3,db);
    else Abort("%s: internal error: both time series and bucket data!\n",
	       progname);
  }
  else 
    writeIntTag(info,afniHead,"DATASET_RANK",2,3,dt);
  writeIntTag(info,afniHead,"DATASET_DIMENSIONS",3,dx,dy,dz);

  writeStringTag(info,afniHead,"TYPESTRING",typestring_names[afniType]);
  writeIntTag(info, afniHead, "SCENE_DATA", 3, afniViewType, afniFuncType, 
	      afniType);

  if (kvLookup(info,"afni.ORIENT_SPECIFIC.0") 
      && kvLookup(info,"afni.ORIENT_SPECIFIC.1") 
      && kvLookup(info,"afni.ORIENT_SPECIFIC.2")) {
    writeIntTag(info,afniHead,"ORIENT_SPECIFIC",3,
		lookupOrientation(kvGetString(info,"afni.ORIENT_SPECIFIC.0")),
		lookupOrientation(kvGetString(info,"afni.ORIENT_SPECIFIC.1")),
		lookupOrientation(kvGetString(info,"afni.ORIENT_SPECIFIC.2")));
  }
  else {
    /* Use Fiasco standard orientation */
    writeIntTag(info,afniHead,"ORIENT_SPECIFIC",3,
		ORI_R2L_TYPE, ORI_A2P_TYPE, ORI_S2I_TYPE);
  }

  if (kvLookup(info,"voxel_x") && kvLookup(info,"voxel_y")
      && kvLookup(info,"voxel_z"))
    writeDoubleTag(info,afniHead, "DELTA", 3, 
		   kvGetDouble(info,"voxel_x"),
		   kvGetDouble(info,"voxel_y"),
		   -kvGetDouble(info,"voxel_z"));
  else Abort("%s: input file does not contain required voxel size info!\n",
	     progname);

  if (getVec3(info,"blb",blb)) {
    writeDoubleTag(info, afniHead, "ORIGIN", 3, blb[0], -blb[1], -blb[2]);
  }
  else {
    /* AFNI considers the origin to be the center of the top left front
     * voxel, and to3d assumes that the volume is centered at (0,0,0)
     * unless otherwise specified.  We'll follow this convention.
     */
    writeDoubleTag(info, afniHead, "ORIGIN", 3,
		   -0.5*(dx-1)*kvGetDouble(info,"voxel_x"),
		   -0.5*(dy-1)*kvGetDouble(info,"voxel_y"),
		   0.5*(dz-1)*kvGetDouble(info,"voxel_z"));
  }

  if (dt>1) {
    double TR= (kvLookup(info,"TR") ? (double)kvGetInt(info,"TR") : 0.0);
    writeIntTag(info, afniHead,"TAXIS_NUMS",3,
		dt, dz, UNITS_MSEC_TYPE);
    writeDoubleTag(info, afniHead, "TAXIS_FLOATS", 5,
		   0.0, TR/1000.0, dt*TR/1000.0, 
		   -blb[2], -kvGetDouble(info,"voxel_z"));
    calcSliceTimeOffsets(info, afniHead);
  }
  else {
    /* At the moment, we only keep the min and max if there is only a 
     * single brick.  AFNI programs will regenerate it if necessary.
     */
    if (kvLookup(info,"brick_min") && kvLookup(info,"brick_max")) {
      writeDoubleTag(info, afniHead, "BRICK_STATS", 2,
		     kvGetDouble(info,"brick_min"),
		     kvGetDouble(info,"brick_max"));
    }
  }

  if (kvLookup(info,"date")) {
    if (kvLookup(info,"time")) {
      char buf[64];
      sprintf(buf,"%s %s",kvGetString(info,"date"),kvGetString(info,"time"));
      writeStringTag(info,afniHead,"IDCODE_DATE",buf);
    }
    else 
      writeStringTag(info,afniHead,"IDCODE_DATE",kvGetString(info,"date"));
  }

  writeBrickTypes(info, afniHead);

  writeStringTag(info, afniHead, "BYTEORDER_STRING",
		 (bio_big_endian_machine ? "MSB_FIRST" : "LSB_FIRST"));
}

static void exportDependentRequiredTags( KVHash* info, FILE* afniHead )
{
  /* Some AFNI tags are required if other tags are present.  This
   * routine is in charge of checking and fulfilling certain left-over
   * requirements of that sort.
   */
  KVHash* defs= kvGetHash(info,"definitions");
  KVHash* extNames= kvGetHash(info,"external_names");
  KVHash* wroteThese= kvGetHash(info,"wrote_these_already");
  int afniType= pickAfniType(info);
  int afniViewType= pickAfniViewType(info);
  int afniFuncType= pickAfniFuncType(info, afniType);

  /* Fake up an empty set of markers fields */
  if ( afniType==THREEDIM_HEAD_ANAT && afniViewType==VIEW_ORIGINAL_TYPE
       && (!kvLookup(info,"dt") || kvGetInt(info,"dt")==1)) {
    if (!kvLookup(wroteThese,"MARKS_XYZ"))
      writeDoubleTag(info, afniHead, "MARKS_XYZ", 30,
		     -999999.0, -999999.0, -999999.0, -999999.0, -999999.0, 
		     -999999.0, -999999.0, -999999.0, -999999.0, -999999.0, 
		     -999999.0, -999999.0, -999999.0, -999999.0, -999999.0, 
		     -999999.0, -999999.0, -999999.0, -999999.0, -999999.0, 
		     -999999.0, -999999.0, -999999.0, -999999.0, -999999.0, 
		     -999999.0, -999999.0, -999999.0, -999999.0, -999999.0);
    if (!kvLookup(wroteThese,"MARKS_FLAGS"))
      writeIntTag(info, afniHead, "MARKS_FLAGS", 8, 1, 1, 0, 0, 0, 0, 0, 0);
    if (!kvLookup(wroteThese,"MARKS_LAB"))
      writeOneTag(info, afniHead, KV_STRING,"MARKS_LAB",200,
		     "AC superior edge\0\0\0\0AC posterior margin\0PC inferior edge\0\0\0\0First mid-sag pt\0\0\0\0Another mid-sag pt\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0");
  }

}

static void transcribeSimplePair( KVHash* info, FILE* afniHead, 
				  const char* afniKey, const KVPair* p )
{
  if (debug) fprintf(stderr,"Transcribing simple key <%s> -> <%s>, type %s\n",
		     kvKey(p),afniKey,kvTypeName(kvType(p)));
  switch (kvType(p)) {
  case KV_STRING: 
    writeStringTag(info, afniHead, afniKey, kvGetString(info, kvKey(p))); 
    break;
  case KV_LONG: 
    writeIntTag(info, afniHead, afniKey, 1, (int)kvGetLong(info, kvKey(p))); 
    break;
  case KV_INT: 
    writeIntTag(info, afniHead, afniKey, 1, kvGetInt(info, kvKey(p))); 
    break;
  case KV_BOOLEAN:  
    writeIntTag(info, afniHead, afniKey, 1, kvGetBoolean(info, kvKey(p))); 
    break;
  case KV_DOUBLE:  
    writeDoubleTag(info, afniHead, afniKey, 1, kvGetDouble(info, kvKey(p))); 
    break;
  default: 
    /* skip it; it's a hash table */
    if (debug) 
      fprintf(stderr,"Unexpectedly hit hash table in transcribeSimplePair!\n");
    break;
  }
}

static void defCompromiseType( KVHash* typeInfo, const char* key, 
			       const KVType type)
{
  /* Sometimes early values in groups will be mistaken for an incorrect type,
   * like '3' being interpreted as an int, a long long, or a double.  This
   * routine picks a compromise between what the group type was thought to be
   * and what it later turns out to be.
   */
  KVPair* p;

  if ((p=kvLookup(typeInfo,key)) != NULL) {
    KVType oldType= (KVType)kvGetInt(typeInfo,key);
    KVType finalType= type;
    switch (oldType) {
    case KV_BOOLEAN: 
      {
	switch (type) {
	case KV_BOOLEAN:
	case KV_INT: 
	  finalType= oldType;
	  break;
	case KV_LONG:
	case KV_DOUBLE:
	case KV_STRING:
	  finalType= type;
	  break;
	default: 
	  Abort("%s: internal error: nonsense type %d in defCompromiseType!\n",
		progname,type);
	}
      }
      break;
    case KV_INT:
      {
	switch (type) {
	case KV_BOOLEAN:
	case KV_INT: 
	  finalType= oldType;
	  break;
	case KV_LONG:
	case KV_DOUBLE:
	case KV_STRING:
	  finalType= type;
	  break;
	default: 
	  Abort("%s: internal error: nonsense type %d in defCompromiseType!\n",
		progname,type);
	}
      }
      break;
    case KV_LONG:
      {
	switch (type) {
	case KV_BOOLEAN:
	case KV_INT: 
	case KV_LONG:
	  finalType= oldType;
	  break;
	case KV_DOUBLE:
	case KV_STRING:
	  finalType= type;
	  break;
	default: 
	  Abort("%s: internal error: nonsense type %d in defCompromiseType!\n",
		progname,type);
	}
      }
      break;
    case KV_DOUBLE:
      {
	switch (type) {
	case KV_BOOLEAN:
	case KV_INT: 
	case KV_LONG:
	case KV_DOUBLE:
	  finalType= oldType;
	  break;
	case KV_STRING:
	  finalType= type;
	  break;
	default: 
	  Abort("%s: internal error: nonsense type %d in defCompromiseType!\n",
		progname,type);
	}
      }
      break;
    case KV_STRING:
      {
	finalType= type;
      }
      break;
    default: 
      Abort("%s: internal error: nonsense old type %d in defCompromiseType!\n",
	    progname, oldType);
    }
    if (debug) 
      fprintf(stderr,"Compromise group type of <%s> is %s %s -> %s\n",
	      key,kvTypeName(oldType),kvTypeName(type),kvTypeName(finalType));
    kvDefInt(typeInfo, key, (int)finalType);
  }
  else {
    if (debug) 
      fprintf(stderr,"First type of group <%s> is %s\n",key,kvTypeName(type));
    kvDefInt(typeInfo, key, (int)type);
  }
}

static const char* coerceToString( KVHash* info, const char* key )
{
  static char buf[64];
  KVPair* p= kvLookup(info,key);
  if (!p) Abort("%s: coerceToString: internal error: key <%s> has no value!\n",
		progname,key);
  switch (kvType(p)) {
  case KV_STRING: 
    return kvGetString(info,key);
  case KV_INT:
    sprintf(buf,"%d",kvGetInt(info,key));
    return buf;
  case KV_LONG:
    sprintf(buf,"%lld",kvGetLong(info,key));
    return buf;
  case KV_DOUBLE:
    sprintf(buf,"%g",kvGetDouble(info,key));
    return buf;
  case KV_BOOLEAN:
    if (kvGetBoolean(info,key)) return "TRUE";
    else return "FALSE";
  default: Abort("%s: type %s cannot be coerced to a string!\n",
		 progname,kvTypeName(kvType(p)));
  }
  return NULL; /* not reached */
}

static double coerceToDouble( KVHash* info, const char* key )
{
  KVPair* p= kvLookup(info,key);
  if (!p) Abort("%s: coerceToDouble: internal error: key <%s> has no value!\n",
		progname,key);
  switch (kvType(p)) {
  case KV_INT:
    return (double)kvGetInt(info,key);
  case KV_LONG:
    return (double)kvGetLong(info,key);
  case KV_DOUBLE:
    return kvGetDouble(info,key);
  case KV_BOOLEAN:
    return (double)kvGetBoolean(info,key);
  case KV_STRING: 
  default: Abort("%s: type %s cannot be coerced to a double!\n",
		 progname,kvTypeName(kvType(p)));
  }
  return 0.0; /* not reached */
}

static long long coerceToLongLong( KVHash* info, const char* key )
{
  KVPair* p= kvLookup(info,key);
  if (!p) Abort("%s: coerceToLongLong: internal error: key <%s> has no value!\n",
		progname,key);
  switch (kvType(p)) {
  case KV_INT:
    return (long long)kvGetInt(info,key);
  case KV_LONG:
    return kvGetLong(info,key);
  case KV_BOOLEAN:
    return (long long)kvGetBoolean(info,key);
  case KV_DOUBLE:
  case KV_STRING: 
  default: Abort("%s: type %s cannot be coerced to a long long!\n",
		 progname,kvTypeName(kvType(p)));
  }
  return 0; /* not reached */
}

static int coerceToInt( KVHash* info, const char* key )
{
  KVPair* p= kvLookup(info,key);
  if (!p) Abort("%s: coerceToInt: internal error: key <%s> has no value!\n",
		progname,key);
  switch (kvType(p)) {
  case KV_INT:
    return kvGetInt(info,key);
  case KV_BOOLEAN:
    return kvGetBoolean(info,key);
  case KV_LONG:
  case KV_DOUBLE:
  case KV_STRING: 
  default: Abort("%s: type %s cannot be coerced to a double!\n",
		 progname,kvTypeName(kvType(p)));
  }
  return 0; /* not reached */
}

const char* translateNameToAfni( KVHash* info, const char* name )
{
  /* Returning NULL from this routine causes the given tag to be dropped. */
  KVHash* afniNames= kvGetHash(info,"afni_names");
  KVPair* p;
  const char* result= NULL;

  if ((p=kvLookup(afniNames,name)) != NULL) {
    if (strlen(result= kvGetString(afniNames,name))) return result;
    else return NULL;
  }
  else {
    /* Not in the name translation table; try some other things */
    char* tail;
    static char* buf= NULL;
    static int bufSize= 0;
    int minBufSize= strlen("HISTORY_NOTE")+1;
    int neededSize= strlen(name)+1;

    /* Grow the scratch buffer if necessary */
    if (neededSize<minBufSize) neededSize= minBufSize;
    if (neededSize > bufSize) {
      if (buf) free(buf);
      bufSize= neededSize;
      if (!(buf=(char*)malloc(neededSize*sizeof(char))))
	Abort("%s: unable to allocate %d bytes!\n",progname,neededSize);
    }

    if (!strncmp(name,"afni.", strlen("afni."))) {
      /* These pairs were originally imported from an AFNI file. */
      strcpy(buf, name+strlen("afni."));
      if (!strcmp(buf,"BRICK_LABS")) {
	/* AFNI requires that the length of BRICK_LABS match the number
	 * of bricks, or times.  Fiasco may have changed this since the
	 * original AFNI file was translated to Pgh MRI.  Drop this tag
	 * if the number of labels isn't appropriate.
	 */
	int count= 0;
	int dt= 1;
	const char* runner= kvGetString(info,name);
	while (runner= strchr(runner,' ')) { 
	  while (isspace(*runner)) runner++;
	  count++;
	}
	if (kvLookup(info,"dt")) dt= kvGetInt(info,"dt");
	if (dt==count+1) result= buf;
	else result= NULL;
      }
      else result= buf;
    }
    else if (!strncmp(name,"history.", strlen("history."))) {
      /* These pairs contain history info */
      strcpy(buf, "HISTORY_NOTE");
      result= buf;
    }
    else return NULL;

    if ((tail=strrchr(buf,'.')) != NULL) {
      long l;
      char* endptr;
      l= strtol(tail+1,&endptr,0);
      while (isspace(*endptr)) endptr++;
      if (*endptr=='\0') {
	/* Found an index; clip it off the tail end of the name */
	*tail= '\0';
      }
      else {
	/* Do nothing; not an index */
      }
    }
  }
    
  return result;
}

static void exportTagGroups( KVHash* info, KVHash* pendingGroups, 
			     KVHash* pendingGroupTypes, 
			     KVHash* pendingGroupInternalNames,
			     FILE* afniHead )
{
  KVHash* wroteThese= kvGetHash(info,"wrote_these_already");
  KVIterator* kvGroupItr= NULL;
  static char* buf= NULL;
  static int bufSize= 0;

  kvGroupItr= kvUniqueIteratorFactory(pendingGroups);
  while (kvIteratorHasMorePairs(kvGroupItr)) {
    KVPair* p= kvIteratorNextPair(kvGroupItr);
    KVType groupType= (KVType)kvGetInt(pendingGroupTypes, kvKey(p));
    const char* lookupName= kvGetString(pendingGroupInternalNames, kvKey(p));
    int min;
    int max;
    
    
    if (kvLookup(wroteThese, kvKey(p))) {
      if (debug) fprintf(stderr,"Supressing duplicate group <%s>\n",kvKey(p));
      continue;
    }

    /* Grow the buffer if necessary */
    if (strlen(lookupName)+64 > bufSize) {
      if (buf) free(buf);
      bufSize= strlen(lookupName)+64;
      if (!(buf=(char*)malloc(bufSize*sizeof(char))))
	Abort("%s: unable to allocate %d bytes!\n",progname,bufSize);
    }

    snprintf(buf,bufSize,"%s.0",lookupName);
    if (kvLookup(info,buf)) min= 0;
    else min= 1;
    max= kvGetInt(pendingGroups,kvKey(p));
    if (debug) 
      fprintf(stderr,"Found group <%s>, index range <%d> to <%d>, type %s\n",
	      kvKey(p), min, max, kvTypeName(groupType));
    switch (groupType) {
    case KV_STRING:
      {
	int i;
	int totalLength= 0;
	char* s= NULL;
	char* here= NULL;

	for (i=min; i<=max; i++) {
	  snprintf(buf,bufSize,"%s.%d",lookupName,i);
	  /* catch special zero-filled numbering of history info */
	  if (!kvLookup(info,buf))
	    snprintf(buf,bufSize,"%s.%03d",lookupName,i);
	  totalLength += strlen(coerceToString(info,buf)) + 2;
	}
	totalLength += 1; /* for safety */
	if (!(s=(char*)malloc(totalLength*sizeof(char))))
	  Abort("%s: unable to allocate %d bytes!\n",
		progname,totalLength*sizeof(char));
	here= s;
	for (i=min; i<=max; i++) {
	  const char* thisVal;
	  snprintf(buf,bufSize,"%s.%d",lookupName,i);
	  /* catch special zero-filled numbering of history info */
	  if (!kvLookup(info,buf))
	    snprintf(buf,bufSize,"%s.%03d",lookupName,i);
	  thisVal= coerceToString(info,buf);
	  if (i==max) {
	    sprintf(here,"%s",thisVal);
	    here += strlen(thisVal);
	  }
	  else {
	    sprintf(here,"%s\\n",thisVal);
	    here += strlen(thisVal)+2;
	  }
	}
	writeStringTag(info, afniHead, kvKey(p), s);
	free(s);
      }
      break;
    case KV_BOOLEAN:
      {
	int* vals;
	int i;
	if (!(vals=(int*)malloc((max+1-min)*sizeof(int))))
	  Abort("%s: unable to allocate %d bytes!\n",
		progname,(max+1-min)*sizeof(int));
	for (i=min; i<=max; i++) {
	  snprintf(buf,bufSize,"%s.%d",lookupName,i);
	  vals[i-min]= kvGetBoolean(info,buf);
	}
	writeOneTag(info, afniHead, KV_BOOLEAN, kvKey(p), max+1-min, vals);
	free(vals);
      }
      break;
    case KV_INT:
      {
	int* vals;
	int i;
	if (!(vals=(int*)malloc((max+1-min)*sizeof(int))))
	  Abort("%s: unable to allocate %d bytes!\n",
		progname,(max+1-min)*sizeof(int));
	for (i=min; i<=max; i++) {
	  snprintf(buf,bufSize,"%s.%d",lookupName,i);
	  vals[i-min]= coerceToInt(info,buf);
	}
	writeOneTag(info, afniHead, KV_INT, kvKey(p), max+1-min, vals);
	free(vals);
      }
      break;
    case KV_LONG:
      {
	long long* vals;
	int i;
	if (!(vals=(long long*)malloc((max+1-min)*sizeof(long long))))
	  Abort("%s: unable to allocate %d bytes!\n",
		progname,(max+1-min)*sizeof(long long));
	for (i=min; i<=max; i++) {
	  snprintf(buf,bufSize,"%s.%d",lookupName,i);
	  vals[i-min]= coerceToLongLong(info,buf);
	}
	writeOneTag(info, afniHead, KV_LONG, kvKey(p), max+1-min, vals);
	free(vals);
      }
      break;
    case KV_DOUBLE:
      {
	double* vals;
	int i;
	if (!(vals=(double*)malloc((max+1-min)*sizeof(double))))
	  Abort("%s: unable to allocate %d bytes!\n",
		progname,(max+1-min)*sizeof(double));
	for (i=min; i<=max; i++) {
	  snprintf(buf,bufSize,"%s.%d",lookupName,i);
	  vals[i-min]= coerceToDouble(info,buf);
	}
	writeOneTag(info, afniHead, KV_DOUBLE, kvKey(p), max+1-min, vals);
	free(vals);
      }
      break;
    default: Abort("%s: internal error: unexpectedly came on type %d!\n",
		   progname,groupType);
    }
  }
  kvDestroyIterator(kvGroupItr);

}

static void exportOtherChunks( KVHash* info, SList* otherChunks, 
			       FILE* afniHead )
{
  slist_totop(otherChunks);
  while (!slist_empty(otherChunks)) {
    ChunkHandlerPair* chPair= (ChunkHandlerPair*)slist_pop(otherChunks);
    KVHash* subInfo= chPair->info;
    FileHandler* handler= chPair->handler;
    const char* chunkName= kvGetString(subInfo,"chunkname");
    const char* outName= translateNameToAfni(info, chunkName);

    if (debug) fprintf(stderr,"exportOtherChunks: exporting <%s> as <%s>\n",
		       chunkName, outName);

    if (outName && strlen(outName)>0) {
      int nItems;
      int type;
      if (!strcmp(chunkName,"missing")) {
	if (strcmp(kvGetString(subInfo,"dimstr"), "zt")
	    && strcmp(kvGetString(subInfo,"dimstr"), "tz")) {
	  Abort("%s: internal error: expected missing chunk to have dimstr 'zt' or 'tz'!\n",
		progname);

	}
	writeStringTag(info, afniHead, 
		       translateNameToAfni(info,"missing.dimstr"),
		       kvGetString(subInfo,"dimstr"));
	nItems= kvGetInt(subInfo,"dz")*kvGetInt(subInfo,"dt");
	type= SRDR_INT32;
      }
      else {
	/* Reality check; get size and type */
	if (strcmp(kvGetString(subInfo,"dimstr"), "a"))
	  Abort("%s: internal error: expected chunk <%s> to have dimstr 'a'!\n",
		progname, chunkName);
	nItems= kvGetInt(subInfo,"da");
	type= kvGetInt(subInfo,"handler_datatype_out");
      }
      
      switch (type) {
      case SRDR_INT32:
	{
	  int* vals;
	  if (!(vals= (int*)malloc(nItems*sizeof(int))))
	    Abort("%s: unable to allocate %d bytes!\n",nItems*sizeof(int));
	  FH_READ( handler, subInfo, 0, nItems, SRDR_INT32, vals );
	  writeOneTag(info, afniHead, KV_INT, outName, nItems, vals);
	  free(vals);
	}
	break;
      case SRDR_FLOAT32:
	{
	  double* vals;
	  if (!(vals= (double*)malloc(nItems*sizeof(double))))
	    Abort("%s: unable to allocate %d bytes!\n",
		  nItems*sizeof(double));
	  FH_READ( handler, subInfo, 0, nItems, SRDR_FLOAT64, vals );
	  writeOneTag(info, afniHead, KV_DOUBLE, outName, nItems, vals);
	  free(vals);
	}
	break;
      default:
	Abort("%s: internal error: unexpected supplementary data type %d!\n",
	      type);
      }
      
      kvDestroy(subInfo);
      FH_DESTROYSELF(handler);
      free(chPair);
    }
  }
}

static void exportTags( KVHash* info, SList* chunkList, FILE* afniHead )
{
  KVHash* pendingGroups= kvFactory(KV_DEFAULT_SIZE);
  KVHash* pendingGroupTypes= kvFactory(KV_DEFAULT_SIZE);
  KVHash* pendingGroupInternalNames= kvFactory(KV_DEFAULT_SIZE);
  KVIterator* kvi= NULL;

  exportRequiredTags(info, afniHead);

  kvi= kvUniqueIteratorFactory(info);
  while (kvIteratorHasMorePairs(kvi)) {
    KVPair* p= kvIteratorNextPair(kvi);
    const char* internalName= kvKey(p);
    const char* afniName= translateNameToAfni(info, internalName);

    /* Keep only the interesting pairs, watching out for
     * indexed groups of tags.
     */
    if (afniName!=NULL) {
      const char* tail= strrchr(internalName,'.');
      int partOfGroup= 0;
      long l;
      char* endptr;

      if (debug) 
	fprintf(stderr,"Handling key <%s> -> <%s>\n",internalName,afniName);

      if (tail) {
	l= strtol(tail+1,&endptr,0);
	while (isspace(*endptr)) endptr++;
	if (*endptr=='\0') partOfGroup= 1;
	else partOfGroup= 0;
      }
      else partOfGroup= 0;

      if (partOfGroup) {
	static char* buf= NULL;
	static int bufSize= 0;

	/* Grow the scratch buffer if necessary */
	if (strlen(internalName)+1 > bufSize) {
	  if (buf) free(buf);
	  bufSize= strlen(internalName)+1;
	  if (!(buf= (char*)malloc(bufSize*sizeof(char))))
	    Abort("%s: unable to allocate %d bytes!\n",progname,bufSize);
	}

	/* Make a version of the string without the index */
	strncpy(buf,internalName,bufSize);
	*(strrchr(buf,'.'))= '\0';

	/* Keep track of max index for group */
	if (!kvLookup(pendingGroups,afniName) 
	    || (kvGetInt(pendingGroups,afniName)<l)) {
	  if (debug) fprintf(stderr,"new max for <%s> is %ld\n", 
			     afniName, l);
	  kvDefInt(pendingGroups,afniName,(int)l);	
	}

	/* Keep track of best type for group */
	defCompromiseType(pendingGroupTypes, afniName, kvType(p));

	/* Keep track of the internal name base string */
	kvDefString(pendingGroupInternalNames,afniName,buf);
      }
      else {
	/* Key does not end in an int, so it's not a group of vals */
	transcribeSimplePair(info, afniHead, afniName, p);
      }
    }
  }
  kvDestroyIterator(kvi);

  exportTagGroups(info, pendingGroups, pendingGroupTypes, 
		  pendingGroupInternalNames, afniHead);

  exportOtherChunks(info, chunkList, afniHead);

  exportDependentRequiredTags(info, afniHead);

  kvDestroy(pendingGroups);
  kvDestroy(pendingGroupTypes);
  kvDestroy(pendingGroupInternalNames);
}


static void kvDefGeneric( KVHash* kvh, const char* key, KVPair* p )
{
  switch (kvType(p)) {
  case KV_STRING: kvDefString(kvh, key, p->v.s);
    break;
  case KV_LONG: kvDefLong(kvh, key, p->v.l);
    break;
  case KV_DOUBLE: kvDefDouble(kvh, key, p->v.d);
    break;
  case KV_HASH: kvDefHash(kvh, key, kvCloneUnique(p->v.h));
    break;
  case KV_BOOLEAN: kvDefBoolean(kvh, key, (p->v.l != 0));
    break;
  case KV_INT: kvDefInt(kvh, key, (int)p->v.l);
    break;
  }
}		  

static KVHash* mergeTags( SList* chunkList )
{
  KVHash* result= kvFactory(KV_DEFAULT_SIZE);
  SList* otherChunks= slist_create();

  while (!slist_empty(chunkList)) {
    ChunkHandlerPair* chPair= (ChunkHandlerPair*)slist_pop(chunkList);
    KVHash* info= chPair->info;
    KVHash* defs= kvGetHash(info,"definitions");
    KVHash* extNames= kvGetHash(info,"external_names");

    if (!strcmp(kvGetString(info,"chunkname"), "images")) {
      KVIterator* kvi= kvUniqueIteratorFactory(info);
      kvDefHash(result,"definitions",
		kvCloneUnique(kvGetHash(info,"definitions")));
      kvDefHash(result,"external_names",
		kvCloneUnique(kvGetHash(info,"external_names")));
      while (kvIteratorHasMorePairs(kvi)) {
	KVPair* p= kvIteratorNextPair(kvi);
	
	kvDefGeneric(result, kvKey(p), p);
      }
      /* chPair->handler is still referenced, for use transferring
       * data in the main chunk.
       */
      kvDestroy(chPair->info);
      free(chPair);
    }
    else {
      /* We need to snag history, "afni.", and misc. chunk
       * information at this point.
       */
      const char* chunkname= kvGetString(info,"chunkname");
      KVIterator* kvi= kvUniqueIteratorFactory(info);
      while (kvIteratorHasMorePairs(kvi)) {
	KVPair* p= kvIteratorNextPair(kvi);
	if (!strncmp(kvKey(p), "history.",strlen("history."))
	    || !strncmp(kvKey(p), "afni.",strlen("afni."))
	    || !strncmp(kvKey(p), "cl_", strlen("cl_")))
	  kvDefGeneric(result, kvKey(p), p);
      }
      slist_append(otherChunks, chPair);
    }
  }  

  /* Slosh the chunks other than "images" back into the now-empty
   * chunk list.
   */
  while (!slist_empty(otherChunks)) 
    slist_append(chunkList, slist_pop(otherChunks));
  slist_destroy(otherChunks,NULL);

  slist_totop(chunkList);

  return result;
}

#define UPDATE_BOUNDS( type ) \
    { \
      type* here= (type*)buf; \
      type* end= here+countThisBlock; \
      if (!(boundsSet)) { \
	*min= *max= (double)*here; \
	*boundsSet= 1; \
	here++; \
	countThisBlock--; \
      } \
      while (here<end) { \
	if (*here<*min) *min=(double)*here; \
	if (*here>*max) *max=(double)*here; \
	here++; \
      } \
    }

static void updateBounds(void* buf, long countThisBlock, int type, 
			 double* min, double* max, int* boundsSet ) 
{
  switch (type) {
  case SRDR_UINT8:
    UPDATE_BOUNDS(unsigned char);
    break;
  case SRDR_INT16:
    UPDATE_BOUNDS(short);
    break;
  case SRDR_INT32:
    UPDATE_BOUNDS(int);
    break;
  case SRDR_FLOAT32:
    UPDATE_BOUNDS(float);
    break;
  case SRDR_FLOAT64:
    UPDATE_BOUNDS(double);
    break;
  default:
    Abort("%s: Internal error: unknown datatype %d\n",progname,type);
  }
}

#undef UPDATE_BOUNDS

static void transfer( FileHandler* handler, KVHash* info, FILE* afniBrik )
{
  long long totalCount;
  long long offset= 0;
  long countThisBlock;
  void* buf= NULL;
  int itype;
  int otype;
  int dv= 1;
  int dx= 1;
  int dy= 1;
  int dz= 1;
  int dt= 1;
  int db= 1;
  double min;
  double max;
  int maintainBoundsFlag= 0;
  int boundsSet= 0;

  if (kvLookup(info,"dv")) dv= kvGetInt(info,"dv");
  if (kvLookup(info,"dx")) dx= kvGetInt(info,"dx");
  if (kvLookup(info,"dy")) dy= kvGetInt(info,"dy");
  if (kvLookup(info,"dz")) dz= kvGetInt(info,"dz");
  if (kvLookup(info,"dt")) dt= kvGetInt(info,"dt");
  if (kvLookup(info,"db")) db= kvGetInt(info,"db");

  if (dt==1) maintainBoundsFlag= 1;
  
  totalCount= dv*dx*dy*dz*dt*db;
  itype= kvGetInt(info,"handler_datatype_out");
  otype= closestTypeToAfni( itype, dv );

  if (!(buf= (void*)malloc(MAX_BLOCK*srdrTypeSize[otype])))
    Abort("%s: unable to allocate %d bytes!\n",progname,
	  MAX_BLOCK*srdrTypeSize[otype]);

  while (totalCount>0) {
    if (MAX_BLOCK<totalCount) countThisBlock= MAX_BLOCK;
    else countThisBlock= totalCount;

    FH_READ( handler, info, offset*srdrTypeSize[itype], countThisBlock,
	     otype, buf );
    if (maintainBoundsFlag)
      updateBounds(buf, countThisBlock, otype, &min, &max, &boundsSet );
    if (fwrite(buf, srdrTypeSize[otype], (size_t) countThisBlock, afniBrik)
	!= (size_t)countThisBlock) {
      perror("Error writing BRIK");
      Abort("%s: error writing output data BRIK!\n",progname);
    }

    totalCount -= countThisBlock;
    offset += countThisBlock;
  }

  if (maintainBoundsFlag) {
    kvDefDouble(info,"brick_max",max);
    kvDefDouble(info,"brick_min",min);
  }
}

static void addCLToHistory( KVHash* info, int argc, char* argv[] )
{
  int i;
  int totLength;
  char buf[64];
  char* clbuf= NULL;
  char* here;
  int thisEntry;

  /* Find the final existing history entry */
  i= 1;
  while (1) {
    snprintf(buf, sizeof(buf), "history.%d", i);
    if (!kvLookup(info,buf)) {
      snprintf(buf, sizeof(buf), "history.%03d", i);
      if (!kvLookup(info,buf)) break;
    }
    i++;
  } 
  thisEntry= i;

  totLength= 1;
  for (i=0; i<argc; i++) totLength += strlen(argv[i])+1;

  if (!(clbuf=(char*)malloc(totLength*sizeof(char))))
    Abort("%s: unable to allocate %d bytes!\n",progname,
	  totLength*sizeof(char));

  here= clbuf;
  sprintf(here,"%s",argv[0]);
  here+= strlen(argv[0]);
  for (i=1; i<argc; i++) {
    sprintf(here," %s",argv[i]);
    here += strlen(argv[i])+1;
  }

  kvDefString(info, buf, clbuf);

  free(clbuf);
}


int main(int argc, char* argv[])
{
  char inputFileName[512], outFileName[512];
  KVHash* info;
  KVHash* wroteThese= NULL;
  KVHash* goodInfo= NULL;
  KVHash* afniNames= NULL;
  FILE* afniHead= NULL;
  FILE* afniBrik= NULL;
  SList* chunkList= NULL;
  SList* processedChunkList= NULL;
  ChunkHandlerPair* chPair= NULL;
  FileHandler* fileHandler= NULL;

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

  /* We'll keep lists of chunks before and after header processing. */
  chunkList= slist_create();
  processedChunkList= slist_create();

  /* (almost) All knowledge will accumulate as key-value pairs */
  info= kvFactory(KV_DEFAULT_SIZE);
  wroteThese= kvFactory(KV_DEFAULT_SIZE);
  afniNames= kvFactory(KV_DEFAULT_SIZE);
  loadStringTransTable(afniNames, afniNameTransTable);

  /* Miscellaneous initialization tasks- definition and translation
   * tables, default values, opportunity for libbio to figure out
   * endian-ness of this machine.
   */
  InitBIO();
  initStuff(info);

  /* Parse command line */
  parse_command_line(info, inputFileName, outFileName, argc, argv);

  /* Is the input really a Pgh MRI file? */
  if (strcmp(inputFileName+strlen(inputFileName)-4, ".mri")) {
    if (strlen(inputFileName)>sizeof(inputFileName)-5)
      Abort("%s: Input file name <%s> is too long!\n",
	    progname,inputFileName);
    strcat(inputFileName, ".mri");
  }
  if (!pghmriTester(inputFileName))
    Abort("%s: %s is not a Pittsburgh MRI file!\n",
	  progname, inputFileName);

  /* The pghmri handler will be the first pair on our chunk list. */
  if (!(chPair= (ChunkHandlerPair*)malloc(sizeof(ChunkHandlerPair))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,sizeof(ChunkHandlerPair));
  chPair->handler= pghmriFactory(inputFileName, info);
  chPair->info= info;
  slist_append(chunkList,chPair);

  /* We can now loop over the chunks (perhaps discovering more
   * chunks in the process).  Make sure to snag the FileHandler
   * for the actual data.
   */
  while (!slist_empty(chunkList)) {
    chPair= (ChunkHandlerPair*)slist_pop(chunkList);
    FH_PROCESSHEADER(chPair->handler, chPair->info, chunkList);
    slist_append(processedChunkList, chPair);
    if (!strcmp(kvGetString(chPair->info,"chunkname"),"images"))
      fileHandler= chPair->handler;
  }
  slist_destroy(chunkList,NULL);

  if (!canTranslate(processedChunkList))
    Abort("%s: %s cannot be translated to AFNI format.\n",
	  progname, inputFileName);

  openOutputFiles(outFileName, &afniHead, &afniBrik);

  /* Unfortunately, useful information is scattered among the
   * several chunks.  Gather it all in one place.  We keep the
   * chunks other than "images" in the chunk list, so relevant
   * info can be extracted from them.
   */
  goodInfo= mergeTags( processedChunkList );
  kvDefHash(goodInfo,"wrote_these_already",wroteThese);
  kvDefHash(goodInfo,"afni_names",afniNames);

  /* Add this command to the dataset history.  We do this last because
   * history information may have been added by various chunk handlers
   * as well.
   */
  addCLToHistory( goodInfo, argc, argv );

  if (verbose_flg) {
    smart_dump(goodInfo);
    slist_totop(processedChunkList);
    while (!slist_atend(processedChunkList)) {
      ChunkHandlerPair* chPair= 
	(ChunkHandlerPair*)slist_next(processedChunkList);
      fprintf(stderr,"Supplementary chunk: %s\n",
	      kvGetString(chPair->info, "chunkname"));
    }
    slist_totop(processedChunkList);
  }

  transfer( fileHandler, goodInfo, afniBrik );

  exportTags(goodInfo, processedChunkList, afniHead);

  /* Close output dataset */
  closeOutputFiles(afniHead, afniBrik);

  /* Clean up */
  slist_destroy(processedChunkList,NULL);
  kvDestroy(goodInfo);

  fprintf(stderr,"#      Data converted to AFNI format.\n");
  return 0;

}

