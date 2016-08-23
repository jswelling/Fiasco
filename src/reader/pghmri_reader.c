/************************************************************
 *                                                          *
 *  pghmri_reader.c                                         *
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

#include "bio.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"
#include "smartreader.h"

static char rcsid[] = "$Id: pghmri_reader.c,v 1.12 2007/07/03 20:04:55 welling Exp $";

/* Notes-
 * -Remember to force first chunk file to be in .dat!
 */

/* Local forward definitions */
static FileHandler* pghmriDummyFactory(char* fname, KVHash* info);

/* This table provides supplemental translations from external 
 * to internal names.  Empty value strings cause the key to be dropped.
 */
static char* suppNameTransTable[][2]= {
  {"datatype", "datatype_in"},
  {"dimensions", "dimstr"},
  {"file", "chunkfile"},
  {"!format", ""},
  {"!version", ""},
  {"offset", ""},
  {"size", ""},
  {"little_endian", ""}, /* libmri will take care of endianness */
  {"pghmri_subtype", ""},
  {"te2", "TE2"},
  {"ti", "TI"},
  {"scan.date", "date"},
  {"scan.time", "time"},
  {"flip_angle","flip"},
  {"normal.x","slice_norm.0"},
  {"normal.y","slice_norm.1"},
  {"normal.z","slice_norm.2"},
  {"cl_extent_string", ""},
  {"cl_dim_string", ""},
  {"cl_reorder", ""},
  {"cl_big_endian_input", ""},
  {"cl_autoscale", ""},
  {NULL, NULL} /* ends list */
};

/* This table provides supplemental value type information. */
static KeyTypePair suppKeyTypeTable[]= {
  {"TE2",KV_INT},
  {"TI",KV_INT},
  {"origin.x",KV_DOUBLE},
  {"origin.y",KV_DOUBLE},
  {"origin.z",KV_DOUBLE},
  {"translate.x",KV_DOUBLE},
  {"translate.y",KV_DOUBLE},
  {"translate.z",KV_DOUBLE},
  {"afni.TYPESTRING",KV_STRING},
  {"afni.ORIGIN.0",KV_DOUBLE},
  {"afni.ORIGIN.1",KV_DOUBLE},
  {"afni.ORIGIN.2",KV_DOUBLE},
  {"afni.DELTA.0",KV_DOUBLE},
  {"afni.DELTA.1",KV_DOUBLE},
  {"afni.DELTA.2",KV_DOUBLE},
  {"afni.IDCODE_STRING",KV_STRING},
  {"afni.IDCODE_DATE",KV_STRING},
  {"afni.BYTEORDER_STRING",KV_STRING},
  {"afni.VOLREG_ROTPARENT_IDCODE",KV_STRING},
  {"afni.VOLREG_ROTPARENT_NAME",KV_STRING},
  {"afni.VOLREG_GRIDPARENT_IDCODE",KV_STRING},
  {"afni.VOLREG_GRIDPARENT_NAME",KV_STRING},
  {"afni.VOLREG_INPUT_IDCODE",KV_STRING},
  {"afni.VOLREG_INPUT_NAME",KV_STRING},
  {"afni.VOLREG_BASE_IDCODE",KV_STRING},
  {"afni.VOLREG_BASE_NAME",KV_STRING},
  {"afni.IDCODE_ANAT_PARENT",KV_STRING},
  {"afni.IDCODE_WARP_PARENT",KV_STRING},
  {"afni.MARKS_LAB",KV_STRING},
  {"afni.MARKS_HELP",KV_STRING},
  {"afni.TAGSET_LABELS",KV_STRING},
  {"afni.DATASET_NAME",KV_STRING},
  {"afni.DATASET_KEYWORDS",KV_STRING},
  {"afni.BRICK_KEYWORDS",KV_STRING},
  {"nifti_sform_code",KV_INT},
  {"nifti_qform_code",KV_INT},
  {"nifti_srow_x_0",KV_DOUBLE},
  {"nifti_srow_x_1",KV_DOUBLE},
  {"nifti_srow_x_2",KV_DOUBLE},
  {"nifti_srow_x_3",KV_DOUBLE},
  {"nifti_srow_y_0",KV_DOUBLE},
  {"nifti_srow_y_1",KV_DOUBLE},
  {"nifti_srow_y_2",KV_DOUBLE},
  {"nifti_srow_y_3",KV_DOUBLE},
  {"nifti_srow_z_0",KV_DOUBLE},
  {"nifti_srow_z_1",KV_DOUBLE},
  {"nifti_srow_z_2",KV_DOUBLE},
  {"nifti_srow_z_3",KV_DOUBLE},
  {NULL, KV_INT} /* ends list */
};

static int wildcard_strcmp( const char* s1, const char* s2 )
{
  while (*s1 && *s2) {
    if ((*s1 == *s2) || (*s2 == '?')) { s1++; s2++; }
    else return 1;
  }
  if (*s1 || *s2) return 1;
  return 0;
}

static void pghmriClose( FileHandler* self )
{
  if (self->hook) mri_close_dataset((MRI_Dataset*)self->hook);
  self->hook= NULL;
}

static void pghmriReopen( FileHandler* self )
{
  self->file= NULL;
  if (!self->hook)
    self->hook= mri_open_dataset( self->fileName, MRI_READ );
}

static void pghmriRead( FileHandler* self, KVHash* info,
		      long long offset, long n,
		      SRDR_Datatype datatype_out, void* obuf )
{
  MRI_Dataset* ds;
  MRI_ArrayType mri_arraytype_in;
  long long scaledOffset;
  
  
  pghmriReopen(self);
  ds= (MRI_Dataset*)self->hook;
  
  switch (datatype_out) {
  case SRDR_UINT8: mri_arraytype_in= MRI_UNSIGNED_CHAR; break;
  case SRDR_INT16: mri_arraytype_in= MRI_SHORT; break;
  case SRDR_INT32: mri_arraytype_in= MRI_INT; break;
  case SRDR_FLOAT32: mri_arraytype_in= MRI_FLOAT; break;
  case SRDR_FLOAT64: mri_arraytype_in= MRI_DOUBLE; break;
  case SRDR_INT64: mri_arraytype_in= MRI_LONGLONG; break;
  default: 
    Abort("%s: internal error: Pgh MRI file reader can't produce type %s!\n",
	  progname,srdrTypeName[datatype_out]);
  }
  scaledOffset= offset/libmriArrayTypeSize[mri_arraytype_in];
  mri_read_chunk(ds, kvGetString(info,"chunkname"), n, scaledOffset, 
		 mri_arraytype_in, obuf);
}

static int isOrphan( const char* key, const char* currentChunk )
{
  if (!currentChunk) return 1; /* if there is no chunk, it must be orphan! */
  if (strncmp(key, currentChunk, strlen(currentChunk)))
    return 1;
  if (key[strlen(currentChunk)] != '.') return 1;

  return 0;
}

char* translateKey( const char* key, const char* currentChunk,
		    KVHash* transHash )
{
  char buf[256];

  /* Any key reaching this point is not an orphan, so it begins with
   * the chunk name followed by '.' .
   */
  if (key[strlen(currentChunk)]=='.')
    key += strlen(currentChunk)+1; /* skip over the '.' */
  else return NULL;

  /* A number of translations are in the hash table provided. */
  if (kvLookup(transHash, key)) key= kvGetString(transHash, key);

  /* Handle some special cases */
  if (!wildcard_strcmp(key,"extent.?")) {
    snprintf(buf,sizeof(buf),"d%c",key[7]);
    key= buf;
  }

  return strdup(key);
}

static KVHash* buildInverseNameMap( KVHash* extNames )
{
  KVHash* result= kvFactory(KV_DEFAULT_SIZE);
  KVIterator* kvi= kvUniqueIteratorFactory( extNames );

  while ( kvIteratorHasMorePairs(kvi) ) {
    KVPair* p= kvIteratorNextPair(kvi);
    const char* key= kvKey( p );
    char* val= kvGetString( extNames, key );
    if (strlen(val)>0) kvDefString(result,val,key);
  }

  return result;
}

static SRDR_Datatype translateDatatype( const char* val )
{
  if (!strcmp(val,"uint8")) return SRDR_UINT8;
  else if (!strcmp(val,"int16")) return SRDR_INT16;
  else if (!strcmp(val,"int32")) return SRDR_INT32;
  else if (!strcmp(val,"float32")) return SRDR_FLOAT32;
  else if (!strcmp(val,"float64")) return SRDR_FLOAT64;
  else Abort("%s: internal error: unknown datatype <%s> in translateDatatype!\n",
	     progname,val);
  return 0; /* not reached */
}

static void defineKVP( KVHash* info, const char* key, const char* val )
{
  KVHash* types= kvGetHash(info,"expected_types");

  if (strstr(key,"datatype")) {
    /* Datatypes are stored as strings in Pgh MRI files, and as 
     * vaguely-related ints internally.
     */
    kvDefInt(info,key,translateDatatype(val));
  }
  else if (kvLookup(types,key)) {
    switch (kvGetInt(types,key)) {
    case KV_STRING:
      kvDefString(info,key,val);
      break;
    case KV_LONG:
      kvDefLong(info,key,atoi(val));
      break;
    case KV_DOUBLE:
      kvDefDouble(info,key,atof(val));
      break;
    case KV_BOOLEAN: 
      kvDefBoolean(info,key,(atoi(val) != 0));
      break;
    case KV_INT:
      kvDefInt(info,key,(int)atoi(val));
      break;
    default:
      Abort("%s: internal error: unknown type %d for key <%s> in defineKVP!\n",
	    progname,kvGetInt(types,key),key);
    }
  }
  else if (!wildcard_strcmp(key,"d?")) {
    kvDefInt(info,key,atoi(val));
  }
  else {
    char* endptr;
    long long lval;
    lval= strtoll(val, &endptr, 0);
    while (isspace(*endptr)) endptr++;
    if (*endptr == '\0') kvDefLong(info,key,lval);
    else {
      double dval= strtod(val, &endptr);
      while (isspace(*endptr)) endptr++;
      if (*endptr == '\0') kvDefDouble(info,key,dval);
      else kvDefString(info,key,val);
    }
  }
}

static void translateValuesBySubtype(KVHash* info)
{
  /* Some adjustment is needed to translate key-value pairs between
   * Pgh MRI files originating from different sources.
   */
  const char* subtype= kvGetString(info,"pghmri_subtype");

  if (!strcasecmp(subtype,"lxconvert")) {
    if (kvLookup(info,"TE")) 
      kvDefInt(info,"TE",1000*kvGetInt(info,"TE"));
    if (kvLookup(info,"TR")) 
      kvDefInt(info,"TR",1000*kvGetInt(info,"TR"));
    if (kvLookup(info,"TE2")) 
      kvDefInt(info,"TE2",1000*kvGetInt(info,"TE2"));
    if (kvLookup(info,"TI")) 
      kvDefInt(info,"TI",1000*kvGetInt(info,"TI"));
    if (!kvLookup(info,"rowflip"))
      kvDefInt(info,"rowflip",0);
    if (!kvLookup(info,"rowflip_pattern"))
      kvDefString(info,"rowflip_pattern","none");

    if (kvLookup(info,"voxel_z") && kvLookup(info,"slice_thickness")) {
      kvDefDouble(info,"slice_gap",
		  kvGetDouble(info,"voxel_z") 
		  - kvGetDouble(info,"slice_thickness"));
    }
    if (kvLookup(info,"origin.x") && kvLookup(info,"origin.y")
	&& kvLookup(info,"origin.z") && kvLookup(info,"fov_x")
	&& kvLookup(info,"fov_y") && kvLookup(info,"fov_z")) {
      double ctr_x= kvGetDouble(info,"origin.x");
      double ctr_y= kvGetDouble(info,"origin.y");
      double ctr_z= kvGetDouble(info,"origin.z");
      double h_fov_x= 0.5*kvGetDouble(info,"fov_x");
      double h_fov_y= 0.5*kvGetDouble(info,"fov_y");
      double h_fov_z= 0.5*kvGetDouble(info,"fov_z");
      kvDefDouble(info,"tlf.0", ctr_x - h_fov_x);
      kvDefDouble(info,"tlf.1", ctr_y + h_fov_y);
      kvDefDouble(info,"tlf.2", ctr_z + h_fov_z);
      kvDefDouble(info,"trf.0", ctr_x + h_fov_x);
      kvDefDouble(info,"trf.1", ctr_y + h_fov_y);
      kvDefDouble(info,"trf.2", ctr_z + h_fov_z);
      kvDefDouble(info,"tlb.0", ctr_x - h_fov_x);
      kvDefDouble(info,"tlb.1", ctr_y - h_fov_y);
      kvDefDouble(info,"tlb.2", ctr_z + h_fov_z);
      kvDefDouble(info,"trb.0", ctr_x + h_fov_x);
      kvDefDouble(info,"trb.1", ctr_y - h_fov_y);
      kvDefDouble(info,"trb.2", ctr_z + h_fov_z);
      kvDefDouble(info,"blf.0", ctr_x - h_fov_x);
      kvDefDouble(info,"blf.1", ctr_y + h_fov_y);
      kvDefDouble(info,"blf.2", ctr_z - h_fov_z);
      kvDefDouble(info,"brf.0", ctr_x + h_fov_x);
      kvDefDouble(info,"brf.1", ctr_y + h_fov_y);
      kvDefDouble(info,"brf.2", ctr_z - h_fov_z);
      kvDefDouble(info,"blb.0", ctr_x - h_fov_x);
      kvDefDouble(info,"blb.1", ctr_y - h_fov_y);
      kvDefDouble(info,"blb.2", ctr_z - h_fov_z);
      kvDefDouble(info,"brb.0", ctr_x + h_fov_x);
      kvDefDouble(info,"brb.1", ctr_y - h_fov_y);
      kvDefDouble(info,"brb.2", ctr_z - h_fov_z);
    }
  }
}

static void finishChunk( KVHash* info, const char* currentChunk )
{
  /* This routine does cleanup stuff once all the tags for a chunk
   * have been read.
   */

  if (debug) fprintf(stderr,"Finishing chunk <%s>!\n",currentChunk);
  if (kvLookup(info,"datatype_in"))
    kvDefInt(info,"handler_datatype_out",kvGetInt(info,"datatype_in"));
  else Abort("%s: chunk %s fails to define datatype!\n",
	     progname,currentChunk);

  if (kvLookup(info,"pghmri_subtype"))
    translateValuesBySubtype(info);
}

static KVHash* startChunk( KVHash* info, const char* chunkName )
{
  KVHash* newInfo= kvFactory(KV_DEFAULT_SIZE);
  if (debug) fprintf(stderr,"Adding new chunk <%s>!\n",chunkName);
  initInfoHash(newInfo);
  kvDefString(newInfo,"chunkname",chunkName);
  kvDefLong(newInfo, "start_offset", 0);
  kvDefString(newInfo, "chunkfile", "");
  return newInfo;
}

static void processHeader( FileHandler* self, KVHash* info, SList* cStack )
{
  MRI_Dataset* ds= NULL;
  char* thisKey;
  char* currentChunk= NULL;
  KVHash* internalNames= NULL;
  KVHash* types= kvGetHash(info,"expected_types");
  int i;

  /* Call the base class method */
  baseProcessHeader( self, info, cStack );

  /* We have to translate some external names to their internal 
   * equivalents, so that they will get properly translated back
   * the other way when the output file is created.
   */
  internalNames= buildInverseNameMap( kvGetHash(info,"external_names") );
  loadStringTransTable(internalNames, suppNameTransTable);
  
  /* We also need to supplement the global type table, to deal with
   * some foreign value types.
   */
  loadStringTypeTable(types, suppKeyTypeTable);

  ds= mri_open_dataset(self->fileName, MRI_READ);

  mri_iterate_over_keys(ds);
  while ((thisKey= mri_next_key(ds)) != NULL) {
    if (!strcmp(mri_get_string(ds,thisKey),"[chunk]")) {
      if (debug) fprintf(stderr,"Found chunk <%s>!\n",thisKey);
      if (currentChunk) {
	/* This is not the first chunk we've dealt with */
	ChunkHandlerPair* chPair= NULL;

	if (!(chPair= (ChunkHandlerPair*)malloc(sizeof(ChunkHandlerPair))))
	  Abort("%s: unable to allocate %d bytes!\n",
		progname,sizeof(ChunkHandlerPair));

	finishChunk(info,currentChunk);
	chPair->info= startChunk(info, thisKey);
	chPair->handler= pghmriDummyFactory(self->fileName,info);
	slist_push(cStack,chPair);

	/* Swap over to the new info table, preserving relevant tags */
	if (kvLookup(info,"pghmri_subtype"))
	  kvDefString(chPair->info,"pghmri_subtype",
		      kvGetString(info,"pghmri_subtype"));
	info= chPair->info;
      }
      else {
	/* First chunk; most of the setup is already done. */
	kvDefString(info,"chunkname",thisKey);
      }
      currentChunk= thisKey;
    }
    else if (!strncmp(thisKey,"history.",strlen("history."))) {
      /* Preserve history info */
      if (debug) fprintf(stderr,"Found history key <%s>\n",thisKey);
      kvDefString(info, thisKey, mri_get_string(ds, thisKey));
    }
    else if (!strncmp(thisKey,"!origin",strlen("!origin"))) {
      if (strstr(mri_get_string(ds,thisKey),"lxconvert")) {
	kvDefString(info,"pghmri_subtype","lxconvert");
      }
    }
    else {
      if (isOrphan(thisKey, currentChunk)) {
	if (debug) fprintf(stderr,"Key <%s> is an orphan.\n",thisKey);
      }
      else {
	char* transKey= translateKey( thisKey, currentChunk, internalNames );
	if (transKey) {
	  if (debug) {
	    const char* s= mri_get_string(ds,thisKey);
	    if (strlen(s)<=20) 
	      fprintf(stderr,"Key <%s> translates to <%s>, value <%.20s>\n",
		      thisKey, transKey, mri_get_string(ds,thisKey));
	    else
	      fprintf(stderr,"Key <%s> translates to <%s>, value <%.20s...>\n",
		      thisKey, transKey, mri_get_string(ds,thisKey));
	  }
	  if (strlen(transKey)) 
	    defineKVP(info, transKey, mri_get_string(ds,thisKey));
	  free(transKey);
	}
	else {
	  if (debug) fprintf(stderr,"Key <%s> translates to NULL\n",
			     thisKey);
	}
      }
    }
  }

  finishChunk(info, currentChunk);
  kvDestroy( internalNames );
  mri_close_dataset(ds);
}

static FileHandler* pghmriDummyFactory(char* fname, KVHash* info)
{
  /* We use this handler type to avoid re-parsing the 
   * Pgh MRI header 
   */

  FileHandler* result= pghmriFactory(fname, info);
  result->processHeader= baseProcessHeader;
  result->typeName= strdup( "Pittsburgh MRI (secondary chunk)" );
  return result;
}

FileHandler* pghmriFactory(char* fname, KVHash* info)
{
  FileHandler* result= baseFactory(fname);
  result->typeName= strdup( "Pittsburgh MRI" );
  result->processHeader= processHeader;
  result->read= pghmriRead;
  result->close= pghmriClose;
  result->reopen= pghmriReopen;
  return result;
}

int pghmriTester(const char* filename)
{
  char buf[256];
  FILE* fphead= NULL;
  int ierror= 0;

  if ((fphead = fopen(filename,"r"))!=NULL) {
    if (fgets(buf,sizeof(buf),fphead)!=NULL) {
      if (strncasecmp(buf,"!format = pgh",13)
	  && strncasecmp(buf,"!format= pgh",12)
	  && strncasecmp(buf,"!format =pgh",12)
	  && strncasecmp(buf,"!format=pgh",11)) {
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

