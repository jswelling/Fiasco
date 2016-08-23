/************************************************************
 *                                                          *
 *  afni_reader.c                                         *
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
 *  Original programming by Joel Welling 7/2002             *
 ************************************************************/

/* This module reads AFNI .HEAD and .BRIK files.  It's based
 * heavily on information from Bob Cox' document "Attributes 
 * in the AFNI Dataset Header".  Some implicit beliefs are:
 *
 *  -The only things corresponding to leading 'v' values in
 *   a FIASCO dimstr are complex (dv=2) and rgb (dv=3).  These
 *   are the only cases where dv is set and dv != 1.
 *
 */

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

#include "afni_defs.h"

static char rcsid[] = "$Id: afni_reader.c,v 1.15 2005/02/18 20:31:28 welling Exp $";

#define SINGLE_QUOTE 0x27

static char* strbuf= NULL;
static int strbuf_length= 0;
static int* intbuf= NULL;
static int intbuf_length= 0;
static float* floatbuf= NULL;
static int floatbuf_length= 0;

typedef struct afni_data_struct {
  float* float_facs;
  int facs_length;
  int flip_z;
  long slicesize;
  long slicesize_bytes;
  long nslices;
  long blocksize;
  long blocksize_bytes;
  long long base_offset;
} AfniData;

static void checkIntValidity( char* name, int* buf, int count,
			      int minVals, int lowerLimit, int upperLimit )
{
  int i;

  /* This set of checks just happens to be common */
  if (count<minVals) 
    Abort("%s: afni %s length should be at least %d!\n",
	  progname,minVals);

  for (i=0; i<minVals; i++)
    if (buf[i]<lowerLimit)
      Abort("%s: afni %s value %d is out of range (%d with min %d)!\n",
	    progname, i, buf[i], lowerLimit);

  if (upperLimit>lowerLimit)
    for (i=0; i<minVals; i++)
      if (buf[i]>upperLimit)
	Abort("%s: afni %s value %d is out of range (%d with max %d)!\n",
	      progname, i, buf[i], upperLimit);

}

static void checkFloatValidity( char* name, float* buf, int count,
				int minVals )
{
  int i;

  /* This set of checks just happens to be common */
  if (count<minVals) 
    Abort("%s: afni %s length should be at least %d!\n",
	  progname,minVals);

}

static int alwaysMakeThisAChunk( const char* name )
{
  if (!strcmp(name,"counts") || !strcmp(name,"counts2") 
      || !strcmp(name,"dof") || !strcmp(name,"dof2") 
      || !strcmp(name,"missing")) 
    return 1;
  else return 0;
}

static int neverMakeThisAChunk( const char* name )
{
  if (!strcmp(name,"TAXIS_OFFSETS"))
    return 1;
  else return 0;
}

static const char* translateName( const char* name )
{
  static char tbuf[THD_MAX_NAME+64];
  if (!strncmp(name,"pghmri.",strlen("pghmri."))) {
    sprintf(tbuf,"%s",name+strlen("pghmri."));
  }
  else sprintf(tbuf,"afni.%s",name);

  return tbuf;
}

static void handleStringAttribute( FileHandler* self, 
				   char* name, char* str, 
				   KVHash* info, SList* cStack)
{
  char tbuf[THD_MAX_NAME+64];
  strcpy(tbuf, translateName(name));
  if (!strcmp(name,"TYPESTRING")) {
    if (kvLookup(info,"afni.SCENE_DATA.2")) {
      /* Enforce a consistency requirement */
      if (strcmp(str,kvGetString(info,"afni.SCENE_DATA.2"))) 
	Abort("%s: afni SCENE_DATA[2] is inconsistent with TYPESTRING!\n",
	      progname);
    }
    kvDefString(info,tbuf,str);
  }
  else if (!strcmp(name,"BYTEORDER_STRING")) {
    if (!strcmp(str,"LSB_FIRST")) {
      kvDefBoolean(info,"big_endian_input",0);
    }
    else if (!strcmp(str,"MSB_FIRST")) {
      kvDefBoolean(info,"big_endian_input",1);
    }
    else Abort("%s: afni BYTEORDER_STRING is unparsable.\n",progname);
  }
  else if (!strcmp(name,"HISTORY_NOTE")) {
    /* Translate this to FIASCO history format */
    char histKey[256];
    char* here= str;
    char* newHere;
    int i= 1;
    while (newHere= strstr(here,"\\n")) {
      *newHere= '\0';
      sprintf(histKey,"history.%d",i++);
      kvDefString(info,histKey, here);
      here= newHere+2; /* skip the 'n' */
    }
    sprintf(histKey,"history.%d",i);
    kvDefString(info,histKey, here);
  }
  else if (!strcmp(name,"BRICK_LABS")) {
    /* Let's just drop them.  The label string can be long enough
     * to confuse libmri, and they don't seem to have much content.
     */
  }
  else if (!strcmp(name,"BRICK_KEYWORDS")) {
    /* Let's just drop them.  The label string can be long enough
     * to confuse libmri, and they don't seem to have much content.
     */
  }
  else kvDefString(info,tbuf,str);
}

static void handleIntAttribute( FileHandler* self,
				char* name, int* buf, int count,
				KVHash* info, SList* cStack)
{
  if (!strcmp(name,"DATASET_RANK")) {
    checkIntValidity( name, buf, count, 2, 0, -1 );
    if (buf[1]<=1) {
      kvDefString(info,"dimstr","xyz");
    }
    else if (buf[1]>1) {
      /* multiple sub-bricks map to multiple times */
      kvDefString(info,"dimstr","xyzt");
      kvDefInt(info,"dt",buf[1]);
      kvDefString(info,"description.t","gridded image-space");
    }
  }
  else if (!strcmp(name,"DATASET_DIMENSIONS")) {
    checkIntValidity( name, buf, count, 3, 1, -1);
    kvDefInt(info,"dx",buf[0]);
    kvDefString(info,"description.x","gridded image-space");
    kvDefInt(info,"dy",buf[1]);
    kvDefString(info,"description.y","gridded image-space");
    kvDefInt(info,"dz",buf[2]);
    kvDefString(info,"description.z","gridded image-space");
  }
  else if (!strcmp(name,"SCENE_DATA")) {
    const char* s;
    checkIntValidity( name, buf, count, 3, 0, -1 );
    if (buf[0]<3) 
      kvDefString(info,"afni.SCENE_DATA.0", coordsys_names[buf[0]]);
    else Abort("%s: invalid afni SCENE_DATA[0] value %d!\n",
	       progname, buf[0]);
  
    if (buf[2]<4) {
      if (kvLookup(info,"afni.TYPESTRING")) {
	/* Enforce a consistency requirement */
	if (strcmp(typestring_names[buf[2]], 
		   kvGetString(info,"afni.TYPESTRING")))
	  Abort("%s: afni SCENE_DATA[2] is inconsistent with TYPESTRING!\n",
		progname);
      }	  
      kvDefString(info,"afni.SCENE_DATA.2",typestring_names[buf[2]]);
    }
    else Abort("%s: invalid afni SCENE_DATA[2] value %d!\n",
	       progname, buf[2]);
    if (buf[2]==0 || buf[2]==2) s= anat_type_names[buf[1]];
    else s= func_type_names[buf[1]];
    kvDefString(info,"afni.SCENE_DATA.1",s);
  }
  else if (!strcmp(name,"ORIENT_SPECIFIC")) {
    checkIntValidity(name, buf, count, 3, 0, 5);
    kvDefString(info,"afni.ORIENT_SPECIFIC.0", orientation_names[buf[0]]);
    kvDefString(info,"afni.ORIENT_SPECIFIC.1", orientation_names[buf[1]]);
    kvDefString(info,"afni.ORIENT_SPECIFIC.2", orientation_names[buf[2]]);
  }
  else if (!strcmp(name,"BRICK_TYPES")) {
    int typeOfFirst= buf[0];
    int i;
    checkIntValidity(name, buf, count, count, 0, 6);
    for (i=1; i<count; i++) {
      if (buf[i] != typeOfFirst)
	Abort("%s: afni_reader: mixed brick types are not supported!\n",
	      progname);
    }
    switch (typeOfFirst) {
    case 0:
      kvDefInt(info, "datatype_in", SRDR_UINT8);
      break;
    case 1:
      kvDefInt(info, "datatype_in", SRDR_INT16);
      break;
    case 2:
      kvDefInt(info, "datatype_in", SRDR_INT32);
      break;
    case 3:
      kvDefInt(info, "datatype_in", SRDR_FLOAT32);
      break;
    case 4:
      kvDefInt(info, "datatype_in", SRDR_FLOAT64);
      break;
    case 5:
      /* AFNI's version of complex is 2 floats */
      kvDefInt(info, "datatype_in", SRDR_FLOAT32);
      kvDefInt(info, "dv", 2);
      break;
    case 6: 
      /* rgb; presumably 3 bytes, but the doc doesn't say */
      kvDefInt(info, "datatype_in", SRDR_UINT8);
      kvDefInt(info, "dv", 3);
      break;
    }
    kvDefInt(info, "handler_datatype_out", kvGetInt(info,"datatype_in"));
  }
  else {
    char tname[THD_MAX_NAME+64];
    if (neverMakeThisAChunk(name) 
	|| (!alwaysMakeThisAChunk(name) && count<=10)) {
      /* Put these in key-value pairs */
      if (count==1) {
	  sprintf(tname,"%s",translateName(name));
	  kvDefInt(info,tname,buf[0]);
      }
      else {
	int i;
	for (i=0; i<count; i++) {
	  sprintf(tname,"%s.%d",translateName(name),i);
	  kvDefInt(info,tname,buf[i]);
	}
      }
    }
    else {
      /* Make these a chunk of their own */
      int* data;
      ChunkHandlerPair* chPair;
      KVHash* subInfo= kvFactory(KV_DEFAULT_SIZE);
      int i;

      if (!(data=(int*)malloc(count*sizeof(int))))
	Abort("%s: unable to allocate %d bytes!\n",
	      progname, count*sizeof(int));
      if (!(chPair= (ChunkHandlerPair*)malloc(sizeof(ChunkHandlerPair))))
	Abort("%s: unable to allocate %d bytes!\n",
	      progname,sizeof(ChunkHandlerPair));
      chPair->info= subInfo;
      chPair->handler= ramDataHandlerFactory(data, count, SRDR_INT32);

      for (i=0; i<count; i++) data[i]= buf[i];
      
      initInfoHash(subInfo);
      strcpy(tname,translateName(name));
      kvDefString(subInfo, "chunkname", tname);
      sprintf(tname,".%s",name);
      kvDefString(subInfo, "chunkfile", tname);
      kvDefString(subInfo, "dimstr", "a");
      kvDefInt(subInfo, "da", count);
      kvDefLong(subInfo, "start_offset", 0);
      kvDefInt(subInfo,"datatype_in",SRDR_INT32);
      kvDefInt(subInfo,"handler_datatype_out",SRDR_INT32);
      slist_push(cStack, chPair);
    }
  }
}

static void handleFloatAttribute( FileHandler* self, 
				  char* name, float* buf, int count,
				  KVHash* info, SList* cStack)
{
  if (!strcmp(name,"ORIGIN")) {
    checkFloatValidity(name, buf, count, 3);
    kvDefDouble(info,"afni.ORIGIN.0",buf[0]);
    kvDefDouble(info,"afni.ORIGIN.1",buf[1]);
    kvDefDouble(info,"afni.ORIGIN.2",buf[2]);
  }
  else if (!strcmp(name,"DELTA")) {
    checkFloatValidity(name,buf,count,3);
    kvDefDouble(info,"afni.DELTA.0", buf[0]);
    kvDefDouble(info,"afni.DELTA.1", buf[1]); 
    kvDefDouble(info,"afni.DELTA.2", buf[2]);
  }
  else if (!strcmp(name,"BRICK_FLOAT_FACS")) {
    int i;
    AfniData* data= (AfniData*)(self->hook);
    if (data->facs_length > 0) {
      Warning(1,"%s: redundant BRICK_FLOAT_FACS information!\n",progname);
      free(data->float_facs);
    }
    for (i=0; i<count; i++) {
      if (buf[i]!=0.0) break;
    }
    if (i==count) {
      /* All the factors are 0.0, so we just turn off scaling */
      data->float_facs= NULL;
      data->facs_length= 0;
    }
    else {
      float* facs= NULL;
      if (!(facs=(float*)malloc(count*sizeof(float))))
	Abort("%s: unable to allocate %d bytes!\n",count*sizeof(float));
      for (i=0; i<count; i++) {
	facs[i]= (buf[i]==0.0) ? 1.0 : buf[i];
      }
      data->facs_length= count;
      data->float_facs= facs;
    }
  }
  else if (!strcmp(name,"BRICK_STATS")) {
    /* FIASCO will not maintain these properly, so it's best to just
     * let this information go to blazes and force AFNI to regenerate
     * it when the dataset is converted back.
     */
#ifdef never
    /* BRICK_STATS is actually a set of min/max limits for bricks; let's
     * make a chunk and fold it to the customary FIASCO order of vt.
     */
    float* data;
    ChunkHandlerPair* chPair;
    KVHash* subInfo= kvFactory(KV_DEFAULT_SIZE);
    int i;
    
    if (!(data=(float*)malloc(count*sizeof(float))))
      Abort("%s: unable to allocate %d bytes!\n",
	    progname, count*sizeof(float));
    if (!(chPair= (ChunkHandlerPair*)malloc(sizeof(ChunkHandlerPair))))
      Abort("%s: unable to allocate %d bytes!\n",
	    progname,sizeof(ChunkHandlerPair));
    chPair->info= subInfo;
    chPair->handler= ramDataHandlerFactory(data, count, SRDR_FLOAT32);
    
    for (i=0; i<count; i++) data[i]= buf[i];
    
    initInfoHash(subInfo);
    kvDefString(subInfo, "chunkname", "afni.BRICK_STATS");
    kvDefString(subInfo, "chunkfile", ".BRICK_STATS");
    kvDefString(subInfo, "dimstr", "vt");
    kvDefInt(subInfo, "dv", 2);
    kvDefString(subInfo, "description.v", "complex real/imaginary");
    kvDefInt(subInfo, "dt", count/2);
    kvDefString(subInfo, "description.t", "gridded image-space");
    kvDefLong(subInfo, "start_offset", 0);
    kvDefInt(subInfo,"datatype_in",SRDR_FLOAT32);
    kvDefInt(subInfo,"datatype_out",SRDR_FLOAT32);
    slist_push(cStack, chPair);
#endif
  }
  else {
    char tname[THD_MAX_NAME+64];
    
    if (neverMakeThisAChunk(name) 
	|| (!alwaysMakeThisAChunk(name) && count<=10)) {
      /* Put these in key-value pairs */
      if (count==1) {
	sprintf(tname,"%s",translateName(name));
	kvDefDouble(info,tname,buf[0]);
      }
      else {
	int i;
	for (i=0; i<count; i++) {
	  sprintf(tname,"%s.%d",translateName(name),i);
	  kvDefDouble(info,tname,buf[i]);
	}
      }
    }
    else {
      /* Make these a chunk of their own */
      float* data;
      ChunkHandlerPair* chPair;
      KVHash* subInfo= kvFactory(KV_DEFAULT_SIZE);
      int i;

      if (!(data=(float*)malloc(count*sizeof(float))))
	Abort("%s: unable to allocate %d bytes!\n",
	      progname, count*sizeof(float));
      if (!(chPair= (ChunkHandlerPair*)malloc(sizeof(ChunkHandlerPair))))
	Abort("%s: unable to allocate %d bytes!\n",
	      progname,sizeof(ChunkHandlerPair));
      chPair->info= subInfo;
      chPair->handler= ramDataHandlerFactory(data, count, SRDR_FLOAT32);

      for (i=0; i<count; i++) data[i]= buf[i];
      
      initInfoHash(subInfo);
      strcpy(tname,translateName(name));
      kvDefString(subInfo, "chunkname", tname);
      sprintf(tname,".%s",name);
      kvDefString(subInfo, "chunkfile", tname);
      kvDefString(subInfo, "dimstr", "a");
      kvDefInt(subInfo, "da", count);
      kvDefLong(subInfo, "start_offset", 0);
      kvDefInt(subInfo,"datatype_in",SRDR_FLOAT32);
      kvDefInt(subInfo,"datatype_out",SRDR_FLOAT32);
      slist_push(cStack, chPair);
    }
  }
}

static void rescueMissing( FileHandler* self, KVHash* info, SList* cstack )
{
  SList* otherChunks= slist_create();
  ChunkHandlerPair* chPair= NULL;
  
  slist_totop(cstack);
  while (!slist_empty(cstack)) {
    ChunkHandlerPair* tmpPair= (ChunkHandlerPair*)slist_pop(cstack);
    if (!strcmp(kvGetString(tmpPair->info,"chunkname"),"missing")) 
      chPair= tmpPair;
    else slist_append(otherChunks, tmpPair);
  }

  if (chPair != NULL) {
    int len= kvGetInt(info,"dt")*kvGetInt(info,"dz");
    int* intBuf= NULL;
    unsigned char* charBuf= NULL;
    KVHash* subInfo= chPair->info;
    FileHandler* handler= chPair->handler;
    int i;

    if (!(intBuf=(int*)malloc(len*sizeof(int))))
      Abort("%s: unable to allocate %d bytes!\n", len*sizeof(int));
    if (!(charBuf=(unsigned char*)malloc(len*sizeof(unsigned char))))
      Abort("%s: unable to allocate %d bytes!\n", len*sizeof(unsigned char));

    kvDefInt(subInfo,"dt",kvGetInt(info,"dt"));
    kvDefString(subInfo,"description.t",kvGetString(info,"description.t"));
    kvDefInt(subInfo,"dz",kvGetInt(info,"dz"));
    kvDefString(subInfo,"description.z",kvGetString(info,"description.z"));
    kvDefString(subInfo,"dimstr",kvGetString(info,"missing.dimensions"));
    kvDefInt(subInfo,"datatype_out",SRDR_UINT8);
    kvDefString(subInfo,"chunkfile","");

    FH_PROCESSHEADER( handler, subInfo, cstack );
    FH_READ( handler, subInfo, 0, len, SRDR_INT32, intBuf );
    FH_CLOSE( handler );
    FH_DESTROYSELF( handler );

    for (i=0; i<len; i++) charBuf[i]= (intBuf[i] != 0);
    chPair->handler= ramDataHandlerFactory(charBuf, len, SRDR_UINT8);
  }
  slist_append(cstack, chPair);

  while (!slist_empty(otherChunks)) {
    slist_append(cstack, slist_pop(otherChunks));
  }
  slist_destroy(otherChunks,NULL);

  kvDelete(info,"missing.dimensions");
}

static void translateSomePairs( FileHandler* self, KVHash* info, 
				SList* cStack )
{
  KVHash* defs= kvGetHash(info,"definitions");
  KVHash* extNames= kvGetHash(info,"external_names");
  KVIterator* kvi= NULL;
  AfniData* data= (AfniData*)(self->hook);

  /* If a typical Fiasco missing chunk dimension string is
   * present, we probably have to rescue a chunk of missing
   * info from mangled form.
   */
  if (kvLookup(info,"missing.dimensions") 
      && (!strcmp(kvGetString(info,"missing.dimensions"),"zt")
	  || !strcmp(kvGetString(info,"missing.dimensions"),"tz")))
    rescueMissing(self,info,cStack);

  /* The TR is the only value in TAXIS_FLOATS typically used in FIASCO */
  if (kvLookup(info,"afni.TAXIS_FLOATS.1")) {
    kvDefInt(info,"TR",
		1000*kvGetDouble(info,"afni.TAXIS_FLOATS.1"));
    kvDefString(defs,"TR","TR (us)");
  }

  /* We may have mistaken a bucketed dataset for a time series.
   * In AFNI datasets, multiple bucketed values are stored in
   * xyzv order.
   */
  if (kvLookup(info,"afni.SCENE_DATA.1")) {
    if (strstr(kvGetString(info,"afni.SCENE_DATA.1"),"_BUCK_")) {
      /* it's a bucket type; call last dimension "b". */
      if (!strcmp(kvGetString(info,"dimstr"),"xyzt")) {
	kvDefString(info,"dimstr","xyzb");
	kvDefInt(info,"db",kvGetInt(info,"dt"));
	kvDeleteAll(info,"dt");
      }
      else if (!strcmp(kvGetString(info,"dimstr"),"vxyzt")) {
	kvDefString(info,"dimstr","vxyzb");
	kvDefInt(info,"db",kvGetInt(info,"dt"));
	kvDeleteAll(info,"dt");
      }
    }
  }

  if (kvLookup(info,"dv") && kvGetInt(info,"dv")!=1) {
    /* Somewhere along the line we added a v dimension. */
    char tbuf[256];
    if (strlen(kvGetString(info,"dimstr")) > sizeof(tbuf)-2)
      Abort("%s: ridiculously long dimension string <%s>!\n",
	    progname, kvGetString(info,"dimstr"));
    sprintf(tbuf,"v%s",kvGetString(info,"dimstr"));
    kvDefString(info,"dimstr",tbuf);
  }

  if (kvLookup(info,"afni.ORIGIN.0") && kvLookup(info,"afni.DELTA.0")) {
    /* Let's calculate the volume corners and center.
     * AFNI sometimes works in left-handed coordinates, so
     * we have to examine the orientation of things.  We must
     * laboriously untangle the coordinate info for each of several
     * known cases.
     */
    double origin[3];
    double vox[3];
    int dx;
    int dy;
    int dz;

    origin[0]= kvGetDouble(info,"afni.ORIGIN.0");
    origin[1]= kvGetDouble(info,"afni.ORIGIN.1");
    origin[2]= kvGetDouble(info,"afni.ORIGIN.2");
    vox[0]= kvGetDouble(info,"afni.DELTA.0");
    vox[1]= kvGetDouble(info,"afni.DELTA.1");
    vox[2]= kvGetDouble(info,"afni.DELTA.2");
    dx= kvGetInt(info,"dx");
    dy= kvGetInt(info,"dy");
    dz= kvGetInt(info,"dz");

    if (!strcmp(kvGetString(info,"afni.ORIENT_SPECIFIC.0"), "R2L")
	&& !strcmp(kvGetString(info,"afni.ORIENT_SPECIFIC.1"), "A2P")) {
      if (!strcmp(kvGetString(info,"afni.ORIENT_SPECIFIC.2"), "S2I")) {
	/* This is the case found in +orig files generated by pghtoafni */
	if (vox[0]<=0.0 || vox[1]<=0.0 || vox[2]>=0.0) {
	  Abort("%s: AFNI ORIENT_SPECIFIC and DELTA tags are incompatible!\n",
		progname);	
#ifdef never
	  Warning(1,"%s: AFNI ORIENT_SPECIFIC and DELTA tags are incompatible!\n",
		progname);	
	  if (vox[0]<0.0) vox[0]= -vox[0];
	  if (vox[1]<0.0) vox[1]= -vox[0];
	  if (vox[2]>0.0) vox[2]= -vox[0];
#endif
	}
	vox[0]= vox[0];
	vox[1]= vox[1];
	vox[2]= -vox[2];
	origin[0]= origin[0];
	origin[1]= -origin[1];
	origin[2]= -origin[2];
	data->flip_z= 0;
      }
      else if (!strcmp(kvGetString(info,"afni.ORIENT_SPECIFIC.2"), "I2S")) {
	/* This is the case found in +tlrc files generated by adwarp */
	if (vox[0]<=0.0 || vox[1]<=0.0 || vox[2]<=0.0) {
	  Abort("%s: AFNI ORIENT_SPECIFIC and DELTA tags are incompatible!\n",
		progname);
#ifdef never
	  Warning(1,"%s: AFNI ORIENT_SPECIFIC and DELTA tags are incompatible!\n",
		progname);
	  if (vox[0]<0.0) vox[0]= -vox[0];
	  if (vox[1]<0.0) vox[1]= -vox[0];
	  if (vox[2]<0.0) vox[2]= -vox[0];
#endif
	}
	vox[0]= vox[0];
	vox[1]= vox[1];
	vox[2]= vox[2];
	origin[0]= origin[0];
	origin[1]= -origin[1];
	origin[2]= -origin[2]-(dz-1)*vox[2];
	data->flip_z= 1;
      }
      else Abort("%s: nonsense AFNI coordinate system %s-%s-%s!\n",
		 progname,kvGetString(info,"afni.ORIENT_SPECIFIC.0"),
		 kvGetString(info,"afni.ORIENT_SPECIFIC.1"),
		 kvGetString(info,"afni.ORIENT_SPECIFIC.2"));
    }
    else Abort("%s: unsupported AFNI coordinate system %s-%s-%s!\n",
	       progname,kvGetString(info,"afni.ORIENT_SPECIFIC.0"),
	       kvGetString(info,"afni.ORIENT_SPECIFIC.1"),
	       kvGetString(info,"afni.ORIENT_SPECIFIC.2"));

    /* We want to get rid of the afni.ORIENT_SPECIFIC values, so that
     * proper values will be generated if this dataset is converted
     * back to AFNI coordinates.  If data->flip_z has been set the
     * data will be reordered on input to smartreader, so the 
     * original afni.ORIENT_SPECIFIC values may not hold.
     */
    kvDefString(extNames,"afni.ORIENT_SPECIFIC.0","");
    kvDefString(extNames,"afni.ORIENT_SPECIFIC.1","");
    kvDefString(extNames,"afni.ORIENT_SPECIFIC.2","");

    /* Use the untangled origin and step to establish corners and
     * voxel info.
     */
    kvDefDouble(info,"voxel_x",vox[0]);
    kvDefDouble(info,"voxel_y",vox[1]);
    kvDefDouble(info,"voxel_z",vox[2]);
    if (!kvLookup(info,"blb.0")) {
      kvDefDouble(info,"blb.0",origin[0]);
      kvDefDouble(info,"blb.1",origin[1]);
      kvDefDouble(info,"blb.2",origin[2]);
      kvDefString(defs,"blb.0", 
		  "bottom left back X (mm), FIASCO internal coords");
      kvDefString(defs,"blb.1", 
		  "bottom left back Y (mm), FIASCO internal coords");
      kvDefString(defs,"blb.2", 
		  "bottom left back Z (mm), FIASCO internal coords");
    }
    if (!kvLookup(info,"trf.0")) {
      kvDefDouble(info,"trf.0",origin[0]+(dx-1)*vox[0]);
      kvDefDouble(info,"trf.1",origin[1]-(dy-1)*vox[1]);
      kvDefDouble(info,"trf.2",origin[2]+(dz-1)*vox[2]);
      kvDefString(defs,"trf.0", 
		  "top right front X (mm), FIASCO internal coords");
      kvDefString(defs,"trf.1", 
		  "top right front Y (mm), FIASCO internal coords");
      kvDefString(defs,"trf.2", 
		  "top right front Z (mm), FIASCO internal coords");
    }
    if (!kvLookup(info,"tlb.0")) {
      kvDefDouble(info,"tlb.0",origin[0]);
      kvDefDouble(info,"tlb.1",origin[1]);
      kvDefDouble(info,"tlb.2",origin[2]+(dz-1)*vox[2]);
      kvDefString(defs,"tlb.0", 
		  "top left back X (mm), FIASCO internal coords");
      kvDefString(defs,"tlb.1", 
		  "top left back Y (mm), FIASCO internal coords");
      kvDefString(defs,"tlb.2", 
		  "top left back Z (mm), FIASCO internal coords");
    }
    if (!kvLookup(info,"blf.0")) {
      kvDefDouble(info,"blf.0",origin[0]);
      kvDefDouble(info,"blf.1",origin[1]-(dy-1)*vox[1]);
      kvDefDouble(info,"blf.2",origin[2]);
      kvDefString(defs,"blf.0", 
		  "bottom left front X (mm), FIASCO internal coords");
      kvDefString(defs,"blf.1", 
		  "bottom left front Y (mm), FIASCO internal coords");
      kvDefString(defs,"blf.2", 
		  "bottom left front Z (mm), FIASCO internal coords");
    }
    if (!kvLookup(info,"ctr.0")) {
      kvDefDouble(info,"ctr.0",origin[0] + ((0.5*dx)-1)*vox[0]);
      kvDefDouble(info,"ctr.1",origin[1] - ((0.5*dy)-1)*vox[1]);
      kvDefDouble(info,"ctr.2",origin[2] + ((0.5*dz)-1)*vox[2]);
      kvDefString(defs,"ctr.0", "center X (mm), FIASCO internal coords");
      kvDefString(defs,"ctr.1", "center Y (mm), FIASCO internal coords");
      kvDefString(defs,"ctr.2", "center Z (mm), FIASCO internal coords");
    }
  }

  if (kvLookup(info,"dz") && kvGetInt(info,"dz")>=3
      && kvLookup(info,"afni.TAXIS_OFFSETS.0")) {
    /* Relative timing data is available for the slices, so we should
     * be able to determine the acquisition pattern. 
     */
    double t0= kvGetDouble(info,"afni.TAXIS_OFFSETS.0");
    double t1= kvGetDouble(info,"afni.TAXIS_OFFSETS.1");
    double t2= kvGetDouble(info,"afni.TAXIS_OFFSETS.2");

    kvDefBoolean(info,"reorder",0); /* AFNI assumes they are in order */
    if (t2>t0) {
      if (t1>t0) {
	if (t2>t1) kvDefString(info,"reorder_pattern","sequential");
	else if (t2<t1) kvDefString(info,"reorder_pattern","even/odd");
	else Abort("%s: unknown slice timing pattern: %f %f %f\n",t0,t1,t2);
      }
      else if (t1<t0) {
	if (t1<t2) kvDefString(info,"reorder_pattern","odd/even");
	else Abort("%s: unknown slice timing pattern: %f %f %f\n",
		   progname,t0,t1,t2);
      }
      else Abort("%s: unknown slice timing pattern: %f %f %f\n",
		 progname,t0,t1,t2);
    }
    else if (t2<t0) { /* reversed top-to-bottom order, like raw Siemens */
      if (t1<t0) {
	if (t2<t1) kvDefString(info,"reorder_pattern","reversed_sequential");
	else if (t2>t1) kvDefString(info,"reorder_pattern",
				    "reversed_even/odd");
	else Abort("%s: unknown slice timing pattern: %f %f %f\n",t0,t1,t2);
      }
      else if (t1>t0) {
	if (t1>t2) kvDefString(info,"reorder_pattern","reversed_odd/even");
	else Abort("%s: unknown slice timing pattern: %f %f %f\n",
		   progname,t0,t1,t2);
      }
      else Abort("%s: unknown slice timing pattern: %f %f %f\n",
		 progname,t0,t1,t2);
    }
    else if (t0 == t1 && t1 == t2) {
      /* Apparently no useful data was saved */
    }
    else Abort("%s: unknown slice timing pattern: %f %f %f\n",
	       progname,t0,t1,t2);
    
  }

  /* Make a note of the slice and image sizes, for convenience 
   * during reading 
   */
  data->slicesize= kvGetInt(info,"dx")*kvGetInt(info,"dy");
  if (kvLookup(info,"dv")) data->slicesize *= kvGetInt(info,"dv");
  data->slicesize_bytes= 
    data->slicesize*srdrTypeSize[kvGetInt(info,"datatype_in")];
  data->blocksize= data->slicesize*kvGetInt(info,"dz");
  data->blocksize_bytes= data->slicesize_bytes*kvGetInt(info,"dz");
  data->nslices= kvGetInt(info,"dz");
  data->base_offset= kvGetLong(info,"start_offset");

  if (data->facs_length != 0) {
    /* Somewhere along the way we've picked up scaling factors.  Make
     * sure there are enough, and otherwise set up to apply them on
     * the fly.
     */

    if (kvLookup(info,"dt") && data->facs_length != kvGetInt(info,"dt"))
      Abort("%s: Number of AFNI scaling factors (%d) doesn't match and time series length (%d)!\n",
	    progname, data->facs_length, kvGetInt(info,"dt"));
    if (kvGetInt(info,"handler_datatype_out") != SRDR_FLOAT64) {
      /* double would certainly be good enough */
      kvDefInt(info,"handler_datatype_out",SRDR_FLOAT32);
    }

    /* Force the reading process to break between timesteps */
    kvDefLong(info,"skip.z",0);
  }  
  if (data->flip_z != 0) {
    /* Force the reading process to break between slices */
    kvDefLong(info,"skip.y",0);
  }

}

static void processHeader( FileHandler* self, KVHash* info, SList* cStack )
{
  FILE *fphead;
  int ierror= 0;

  if ((fphead = fopen(self->fileName,"r"))!=NULL) {
    char name[THD_MAX_NAME];
    char typestr[THD_MAX_NAME];
    int count;
    while (!feof(fphead)) {
      if (fscanf( fphead ," type = %s name = %s count = %d" ,
		  typestr , name , &count ) != 3) break ;
      if (debug) fprintf(stderr,"<%s>, type %s, length %d\n",
			 name, typestr, count);

      
      if (!strcmp(typestr,"string-attribute")) {
	int i;
	int c;
	if (count+32 > strbuf_length) {
	  if (strbuf_length>0) free(strbuf);
	  if (!(strbuf= (char*)malloc(count+32)))
	    Abort("%s: unable to allocate %d bytes!\n",
		  progname, count+32);
	  strbuf_length= count+32;
	}
	while (!feof(fphead) && fgetc(fphead)!=SINGLE_QUOTE)
	  { /* discard characters */ }
	for (i=0; i<count; i++) {
	  strbuf[i]= fgetc(fphead);
#ifdef never
	  if (strbuf[i]=='~') strbuf[i]= ' ';
#endif
	  if (feof(fphead)) break;
	}
	strbuf[count-1]= '\0';
	handleStringAttribute(self, name, strbuf, info, cStack);
      }
      else if (!strcmp(typestr,"integer-attribute")) {
	int i;
	if (count > intbuf_length) {
	  if (intbuf_length>0) free(intbuf);
	  if (!(intbuf= (int*)malloc(count*sizeof(int))))
	    Abort("%s: unable to allocate %d bytes!\n",
		  progname, count*sizeof(int));
	  intbuf_length= count;
	}
	for (i=0; i<count && !feof(fphead); i++) 
	  fscanf(fphead,"%d",&(intbuf[i]));
	if (i<count) for (; i<count; i++) intbuf[i]= 0;
	handleIntAttribute(self, name, intbuf, count, info, cStack);
      }
      else if (!strcmp(typestr,"float-attribute")) {
	int i;
	if (count > floatbuf_length) {
	  if (floatbuf_length>0) free(floatbuf);
	  if (!(floatbuf= (float*)malloc(count*sizeof(float))))
	    Abort("%s: unable to allocate %d bytes!\n",
		  progname, count*sizeof(float));
	  floatbuf_length= count;
	}
	for (i=0; i<count && !feof(fphead); i++) 
	  fscanf(fphead,"%g",&(floatbuf[i]));
	if (i<count) for (; i<count; i++) floatbuf[i]= 0;
	handleFloatAttribute(self, name, floatbuf, count, info, cStack);
      }
      else Abort("%s: unrecognized AFNI attribute type <%s>!\n",
		 progname, typestr);
    }
    if (fclose(fphead)) {
      Abort("%s: Error closing header %s: %s\n",
	    progname, self->fileName, strerror(errno));
    }
  }
  else {
    Abort("%s: Error opening header %s: %s\n",
	  progname, self->fileName, strerror(errno));
  }

  translateSomePairs(self,info,cStack);
}

static void afniClose( FileHandler* self )
{
  /* The file pointer is pointed to the BRIK file, but the process
   * to close it is the same as for the base case.
   */
  baseClose(self);
}

static void afniReopen( FileHandler* self )
{
  /* We actually want to reopen the BRIK file rather than the HEAD. */
  if (!(self->file)) {
    char buf[256];
    char* here;
    strncpy(buf, self->fileName, sizeof(buf));
    buf[sizeof(buf)-1]= '\0';
    if (!(here= strrchr(buf,'.')))
      Abort("%s: Unexpectedly can't find the extension in <%s>!\n",
	    progname, buf);
    if (strcmp(here,".HEAD"))
      Abort("%s: Extension of an AFNI header is not .HEAD!\n",progname);
    strcpy(here,".BRIK"); /* just overwrite the letters */
    if (!(self->file= fopen(buf,"r")))
      Abort("%s: unable to open file <%s> for reading!\n",
	    progname,self->fileName);
  }
}

static void afniRead( FileHandler* self, KVHash* info,
		      long long offset, long n,
		      SRDR_Datatype datatype_out, void* obuf )
{
  AfniData *data= (AfniData*)(self->hook);

  /* We need to point ourselves at the .BRIK file */
  afniReopen(self);

  if (data->facs_length==0 && !data->flip_z) {
    /* no scaling required */
    baseRead( self, info, offset, n, datatype_out, obuf );
  }
  else {
    /* scaling and/or reordering is required.  We will fiendishly do 
     * the scaling in place. obuf is guaranteed to be long enough for 
     * the whole output string, and by previous design we always scale 
     * to floats or doubles, never losing precision.
     */
    long img= (long)((offset-data->base_offset)/(data->blocksize_bytes));
    long slice= (long)((offset-(data->base_offset+img*data->blocksize_bytes))
		       /(data->slicesize_bytes));
    long long reordered_offset;
    long i;
    if (img*data->blocksize_bytes + slice*data->slicesize_bytes 
	+ data->base_offset != offset)
      Abort("%s: afniRead: offset %lld is not on a slice boundary!\n",
	    progname, offset);
    if (n + slice*data->slicesize > data->blocksize)
      Abort("%s: afniRead: tried to read across a block boundary!\n",
	    progname);
    if (data->flip_z) 
      reordered_offset= 
	data->base_offset + img*data->blocksize_bytes
	+ (data->nslices-(slice+1))*data->slicesize_bytes;
    else reordered_offset= offset;
    if (data->facs_length==0) {
      /* No scaling required */
      baseRead( self, info, reordered_offset, n, datatype_out, obuf );
    }
    else {
      if (kvGetInt(info,"handler_datatype_out")==SRDR_FLOAT32) {
	float* fbuf= (float*)obuf;
	switch (kvGetInt(info,"datatype_in")) {
	case SRDR_UINT8:
	  {
	    unsigned char* cbuf= (unsigned char*)obuf;
	    baseRead( self, info, reordered_offset, n, SRDR_UINT8, obuf );
	    for (i=n-1; i>=0; i--) fbuf[i]= data->float_facs[img]*cbuf[i];
	  }
	  break;
	case SRDR_INT16:
	  {
	    short* sbuf= (short*)obuf;
	    baseRead( self, info, reordered_offset, n, SRDR_INT16, obuf );
	    for (i=n-1; i>=0; i--) fbuf[i]= data->float_facs[img]*sbuf[i];
	  }
	  break;
	case SRDR_INT32:
	  {
	    int* ibuf= (int*)obuf;
	    baseRead( self, info, reordered_offset, n, SRDR_INT32, obuf );
	    for (i=n-1; i>=0; i--) fbuf[i]= data->float_facs[img]*ibuf[i];
	  }
	  break;
	case SRDR_FLOAT32:
	  {
	    baseRead( self, info, reordered_offset, n, SRDR_FLOAT32, obuf );
	    for (i=n-1; i>=0; i--) fbuf[i]= data->float_facs[img]*fbuf[i];
	  }
	  break;
	case SRDR_FLOAT64:
	  Abort("%s: afniRead: unexpectedly tried to scale double to float!\n",
		progname);
	  break;
	}
      }
      else if (kvGetInt(info,"handler_datatype_out")==SRDR_FLOAT64) {
	double* dbuf= (double*)obuf;
	switch (kvGetInt(info,"datatype_in")) {
	case SRDR_UINT8:
	  {
	    unsigned char* cbuf= (unsigned char*)obuf;
	    baseRead( self, info, reordered_offset, n, SRDR_UINT8, obuf );
	    for (i=n-1; i>=0; i--) dbuf[i]= data->float_facs[img]*cbuf[i];
	  }
	  break;
	case SRDR_INT16:
	  {
	    short* sbuf= (short*)obuf;
	    baseRead( self, info, reordered_offset, n, SRDR_INT16, obuf );
	    for (i=n-1; i>=0; i--) dbuf[i]= data->float_facs[img]*sbuf[i];
	  }
	  break;
	case SRDR_INT32:
	  {
	    int* ibuf= (int*)obuf;
	    baseRead( self, info, reordered_offset, n, SRDR_INT32, obuf );
	    for (i=n-1; i>=0; i--) dbuf[i]= data->float_facs[img]*ibuf[i];
	  }
	  break;
	case SRDR_FLOAT32:
	  {
	    float* fbuf= (float*)obuf;
	    baseRead( self, info, reordered_offset, n, SRDR_FLOAT32, obuf );
	    for (i=n-1; i>=0; i--) dbuf[i]= data->float_facs[img]*fbuf[i];
	  }
	  break;
	case SRDR_FLOAT64:
	  {
	    baseRead( self, info, reordered_offset, n, SRDR_FLOAT32, obuf );
	    for (i=n-1; i>=0; i--) dbuf[i]= data->float_facs[img]*dbuf[i];
	  }
	  break;
	}
      }
      else Abort("%s: unexpected scaling target %s!\n",
		 progname, srdrTypeName[kvGetInt(info,"handler_datatype_out")]);
    }
  }
}

static void destroySelf( FileHandler* self )
{
  AfniData* data= (AfniData*)(self->hook);

  if (data->facs_length>0) free(data->float_facs);

  baseDestroySelf(self);
}

FileHandler* afniFactory(char* fname, KVHash* info)
{
  AfniData* data= NULL;
  FileHandler* result= baseFactory(fname);

  if (!(data= (AfniData*)malloc(sizeof(AfniData))))
    Abort("%s: unable to allocate %d bytes!\n",sizeof(AfniData));
  result->hook= data;
  data->facs_length= 0;
  data->float_facs= NULL;
  data->blocksize= 1;
  data->blocksize_bytes= 1;
  data->slicesize= 1;
  data->slicesize_bytes= 1;
  data->nslices= 1;
  data->base_offset= 0;
  data->flip_z= 0;

  result->typeName= strdup( "AFNI" );
  result->processHeader= processHeader;
  result->read= afniRead;
  result->destroySelf= destroySelf;
  result->close= afniClose;
  result->reopen= afniReopen;

  return result;
}

int afniTester(const char* filename)
{
  const char* here;
  FILE* f;
  char buf[256];
  int ierror= 0;

  /* The filename should end in ".HEAD" */
  here= filename + strlen(filename) - 5;
  if (strcmp(here,".HEAD")) return 0;

  /* We should be able to find some tags in the header */
  if ((f = fopen(filename,"r"))!=NULL)
    {
      if (fread(buf,sizeof(char),sizeof(buf)-1, f)
	  != sizeof(buf)-1) {
	perror("Error reading header");
	ierror=1;
      }
      else {
	if (fclose(f)) {
	  perror("Error closing header");
	  ierror=1;
	}
      }
    }
  else {
    perror("Error opening header");
    ierror= 1;
  }
  
  buf[sizeof(buf)-1]= '\0';

  if (ierror) return 0;

  if (!strstr(buf,"type")) return 0;
  if (!strstr(buf,"name")) return 0;
  if (!strstr(buf,"count")) return 0;
  if (!strstr(buf," = ")) return 0;
  if (!strstr(buf,"-attribute")) return 0;

  /* looks like an AFNI file to me. */
  return 1;
}

