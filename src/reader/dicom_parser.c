/************************************************************
 *                                                          *
 *  dicom_parser.c                                          *
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
 ************************************************************/
/* Notes-
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
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
#include "fexceptions.h"
#include "array.h"
#include "stdcrg.h"
#include "misc.h"
#include "smartreader.h"

/* DICOM header files */
#include "dicom_parser.h"

/* This include file defines the static dicom element table */
#include "dicom_dict.h"

static char rcsid[] = "$Id: dicom_parser.c,v 1.6 2005/03/09 00:05:27 welling Exp $";

typedef enum { STATE_TOP, STATE_SQ, STATE_ITEM } ParserState;
static const char* stateNameTable[]= { "TOP", "SQ", "ITEM"};

typedef struct state_stack_entry_struct {
  ParserState state;
  long long breakOffset;
} StateStackEntry;

char* dicomElementTypeNames[]= {
  /* These correspond to the elements in the enumerated type, in order */
  "UL", "UI", "SH", "US", "AE", "AT", "LO", "OB", "CS", "DS", "OW", "OX", 
  "DL", "SQ", "PN", "ST", "LT", "TM", "DA", "SS", "IS", "SL", "AS", "FD",
  "DT", "FL", "XS", "UT", "OF", "UNKNOWN", "SPECIAL"
};

/* A few special purpose element types */
#define ITEM_ELEMENT 0xe000
#define ITEM_DELIM_ELEMENT 0xe00d
#define SEQ_DELIM_ITEM_ELEMENT 0xe0dd
static DicomDictDataElement dicomSpecialDict[]= {
  { 0xfffe, ITEM_ELEMENT, SPECIAL, "ItemElement" },
  { 0xfffe, ITEM_DELIM_ELEMENT, SPECIAL, "ItemDelimiterElement" },
  { 0xfffe, SEQ_DELIM_ITEM_ELEMENT, SPECIAL, "sequenceDelimiterItemElement" },
};

/* Forward definitions */
static DicomElementType getElementTypeByName( const char vr[2] );

/* These functions raise EXCEPTION_DICOM */
static long long readDataElement(DicomParser* parser, FILE* f, 
				 DicomDataElement* de, long long offset, 
				 const TransferSyntax* transferSyntax);

/* Method prototypes */
static DCM_METHOD_PROTOTYPE(read_longs);
static DCM_METHOD_PROTOTYPE(read_shorts);
static DCM_METHOD_PROTOTYPE(store_offset);
static DCM_METHOD_PROTOTYPE(read_1_string);
static DCM_METHOD_PROTOTYPE(read_1_date);
static DCM_METHOD_PROTOTYPE(read_1_time);
static DCM_METHOD_PROTOTYPE(read_ascii_floats);
static DCM_METHOD_PROTOTYPE(read_ascii_ints);
static DCM_METHOD_PROTOTYPE(read_hex_ints);
static DCM_METHOD_PROTOTYPE(read_and_map_uid);
static DCM_METHOD_PROTOTYPE(raise_exception);

static StateStackEntry* createStateStackEntry( ParserState state, 
					       long long breakOffset )
{
  StateStackEntry* result= NULL;
  if (!(result=(StateStackEntry*)malloc(sizeof(StateStackEntry))))
    Abort("%s: unable to allocate %d bytes!\n",sizeof(StateStackEntry));
  result->state= state;
  result->breakOffset= breakOffset;
  return result;
}

static void destroyStateStackEntry( void* se )
{
  free(se);
}

static long fixShortSign( long i )
{
  /* This stupid standard stores everything as unsigned shorts, which
   * are not supported on most architectures.  Any negative value
   * that comes in here is actually a 16 bit unsigned positive value.
   */
  return (i>=0) ? i : (long)i + 0xffff + 1;
}

static long long fixLongSign( long long i )
{
  /* This stupid standard stores *some* things as unsigned ints, which
   * are not supported on most architectures.  Any negative value
   * that comes in here is actually a 32 bit unsigned positive value.
   */
  return (i>=0) ? i : (long long)i + 0xffffffff + 1;
}

static DCM_METHOD_PROTOTYPE(raise_exception)
{
  /* This is useful for debugging. */
  if (parser->debug) fprintf(stderr,"Raising exception <%s> <%s>\n",key,def);
  fex_raiseException(EXCEPTION_DICOM,"key <%s> def <%s>",key,def);
}

static DCM_METHOD_PROTOTYPE(store_offset)
{
  if (parser->debug) fprintf(stderr,"Storing payload offset\n");
  kvDefLong(info, key, de->payloadOffset);
  if (def != NULL) {
    KVHash* defs= kvGetHash(info,"definitions");
    if (defs) kvDefString(defs,key,def);
  }
}

static DCM_METHOD_PROTOTYPE(read_1_short)
{
  if (parser->debug) fprintf(stderr,"Reading 1 short\n");
  kvDefInt(info, key, (long)fixShortSign(FRdInt16(f)));
  if (def != NULL) {
    KVHash* defs= kvGetHash(info,"definitions");
    if (defs) kvDefString(defs,key,def);
  }
}

static DCM_METHOD_PROTOTYPE(read_1_string)
{
  if (de->length>0) {
    char* buf;
    if (parser->debug) fprintf(stderr,"Reading 1 string, length %d!\n",de->length);
    if (!(buf=(char*)malloc((de->length+1)*sizeof(char))))
      Abort("%s: unable to allocate %d bytes!\n",
	    progname,(de->length+1)*sizeof(char));
    (void)fgets(buf,de->length+1,f);
    if (buf[de->length-1]==' ') 
      buf[de->length-1]= '\0';
    kvDefString(info,key,buf);
    free(buf);
  }
  else {
    kvDefString(info,key,"");
  }
  if (def != NULL) {
    KVHash* defs= kvGetHash(info,"definitions");
    if (defs) kvDefString(defs,key,def);
  }
}

static DCM_METHOD_PROTOTYPE(read_1_date)
{
  if (de->length>0) {
    char buf[64];
    char tbuf[64];
    if (parser->debug) fprintf(stderr,"Reading one date\n");
    if (de->length<sizeof(buf)-1) {
      (void)fgets(buf,de->length+1,f);
      if (buf[de->length-1]==' ') 
	buf[de->length-1]= '\0';
      if (strlen(buf)==8) {
	sprintf(tbuf,"%c%c/%c%c/%c%c%c%c",
		buf[4],buf[5],
		buf[6],buf[7],
		buf[0],buf[1],buf[2],buf[3]);
      }
      else if (strlen(buf)==7) {
	sprintf(tbuf,"%c%c/0%c/%c%c%c%c",
		buf[4],buf[5],
		buf[6],
		buf[0],buf[1],buf[2],buf[3]);
      }
      else {
	fex_raiseException(EXCEPTION_DICOM,
			   "lost (%x,%x), badly formatted date <%s>\n",
			   de->group, de->element, buf);
      }
    }
    else {
      fex_raiseException(EXCEPTION_DICOM,
			 "lost (%x,%x), payload too long (%ld)",
			 de->group,de->element,de->length);
    }
    kvDefString(info,key,tbuf);
  }
  else {
    kvDefString(info,key,"");
  }
  if (def != NULL) {
    KVHash* defs= kvGetHash(info,"definitions");
    if (defs) kvDefString(defs,key,def);
  }
}

static DCM_METHOD_PROTOTYPE(read_1_time)
{
  if (de->length>0) {
    char buf[64];
    char tbuf[64];
    if (parser->debug) fprintf(stderr,"Reading one date\n");
    if (de->length<sizeof(buf)-1) {
      (void)fgets(buf,de->length+1,f);
      if (buf[de->length-1]==' ') 
	buf[de->length-1]= '\0';
      if (strlen(buf)>=6) {
	sprintf(tbuf,"%c%c:%c%c:%s",
		buf[0],buf[1],
		buf[2],buf[3],
		buf+4);
      }
      else {
	fex_raiseException(EXCEPTION_DICOM,
			   "lost (%x,%x), badly formatted time <%s>",
			   de->group,de->element,buf);
      }
    }
    else {
      fex_raiseException(EXCEPTION_DICOM,
			 "lost (%x,%x), payload too long (%ld)",
			 de->group,de->element,de->length);
    }
    kvDefString(info,key,tbuf);
  }
  else {
    kvDefString(info,key,"");
  }
  if (def != NULL) {
    KVHash* defs= kvGetHash(info,"definitions");
    if (defs) kvDefString(defs,key,def);
  }
}

static DCM_METHOD_PROTOTYPE(read_ascii_floats)
{
  if (de->length>0) {
    char buf[256];
    char* keybuf;
    char* defbuf;
    char* here;
    char* key_here;
    char* def_here;
    char* t1;
    char* t2;
    char* t3;
    if (parser->debug) fprintf(stderr,"Reading multiple floats\n");
    if (de->length<sizeof(buf)-1) {
      (void)fgets(buf,de->length+1,f);
      if (buf[de->length-1]==' ') 
	buf[de->length-1]= '\0';
    }
    else {
      fex_raiseException(EXCEPTION_DICOM,
			 "lost (%x,%x), payload too long (%ld)",
			 de->group,de->element,de->length);
    }
    keybuf= (key!=NULL)?strdup(key):NULL;
    defbuf= (def!=NULL)?strdup(def):NULL;
    here= strtok_r(buf,"\\",&t1);
    key_here= (keybuf!=NULL)?strtok_r(keybuf,"\\",&t2):NULL;
    def_here= (defbuf!=NULL)?strtok_r(defbuf,"\\",&t3):NULL;
    while (here != NULL) {
      if (key_here != NULL && strlen(key_here)>0 && strcmp(key_here," ")) {
	kvDefDouble(info,key_here,atof(here));
	if (def_here != NULL) {
	  KVHash* defs= kvGetHash(info,"definitions");
	  if (defs) kvDefString(defs,key_here,def_here);
	}
      }
      else {
	/* Don't define anything for this one */
      }
      here= strtok_r(NULL,"\\",&t1);
      if (key_here) key_here= strtok_r(NULL,"\\",&t2);
      if (def_here) def_here= strtok_r(NULL,"\\",&t3);
    }
    free(keybuf);
    free(defbuf);
  }
}

static DCM_METHOD_PROTOTYPE(read_hex_ints)
{
  if (de->length>0) {
    char buf[256];
    char* keybuf;
    char* defbuf;
    char* here;
    char* key_here;
    char* def_here;
    char* t1;
    char* t2;
    char* t3;
    if (parser->debug) fprintf(stderr,"Reading multiple hex ints\n");
    if (de->length<sizeof(buf)-1) {
      (void)fgets(buf,de->length+1,f);
      buf[de->length-1]= '\0';
    }
    else {
      fex_raiseException(EXCEPTION_DICOM,
			 "lost (%x,%x), payload too long (%ld)",
			 de->group,de->element,de->length);
    }
    fprintf(stderr,"Buf is <%c%c...>\n",buf[0],buf[1]);
    keybuf= (key!=NULL)?strdup(key):NULL;
    defbuf= (def!=NULL)?strdup(def):NULL;
    here= strtok_r(buf,"\\",&t1);
    key_here= (keybuf!=NULL)?strtok_r(keybuf,"\\",&t2):NULL;
    def_here= (defbuf!=NULL)?strtok_r(defbuf,"\\",&t3):NULL;
    while (here != NULL) {
      if (key_here != NULL && strlen(key_here)>0 && strcmp(key_here," ")) {
	long val= strtol(here,NULL,16);
	if (parser->debug) fprintf(stderr,"<%s> -> 0x%x\n",here,val);
	kvDefInt(info,key_here,val);
	if (def_here != NULL) {
	  KVHash* defs= kvGetHash(info,"definitions");
	  if (defs) kvDefString(defs,key_here,def_here);
	}
      }
      else {
	/* Don't define anything for this one */
      }
      here= strtok_r(NULL,"\\",&t1);
      if (key_here) key_here= strtok_r(NULL,"\\",&t2);
      if (def_here) def_here= strtok_r(NULL,"\\",&t3);
    }
    free(keybuf);
    free(defbuf);
  }
}

static DCM_METHOD_PROTOTYPE(read_ascii_ints)
{
  if (de->length>0) {
    char buf[256];
    char* keybuf;
    char* defbuf;
    char* here;
    char* key_here;
    char* def_here;
    char* t1;
    char* t2;
    char* t3;
    if (parser->debug) fprintf(stderr,"Reading multiple ints\n");
    if (de->length<sizeof(buf)-1) {
      (void)fgets(buf,de->length+1,f);
      if (buf[de->length-1]==' ') 
	buf[de->length-1]= '\0';
    }
    else {
      fex_raiseException(EXCEPTION_DICOM,
			 "lost (%x,%x), payload too long (%ld)",
			 de->group,de->element,de->length);
    }
    keybuf= (key!=NULL)?strdup(key):NULL;
    defbuf= (def!=NULL)?strdup(def):NULL;
    here= strtok_r(buf,"\\",&t1);
    key_here= (keybuf!=NULL)?strtok_r(keybuf,"\\",&t2):NULL;
    def_here= (defbuf!=NULL)?strtok_r(defbuf,"\\",&t3):NULL;
    while (here != NULL) {
      if (key_here != NULL && strlen(key_here)>0 && strcmp(key_here," ")) {
	kvDefInt(info,key_here,atoi(here));
	if (def_here != NULL) {
	  KVHash* defs= kvGetHash(info,"definitions");
	  if (defs) kvDefString(defs,key_here,def_here);
	}
      }
      else {
	/* Don't define anything for this one */
      }
      here= strtok_r(NULL,"\\",&t1);
      if (key_here) key_here= strtok_r(NULL,"\\",&t2);
      if (def_here) def_here= strtok_r(NULL,"\\",&t3);
    }
    free(keybuf);
    free(defbuf);
  }
}

static DCM_METHOD_PROTOTYPE(read_longs)
{
  if (de->length>0) {
    char* keybuf;
    char* defbuf;
    char* key_here;
    char* def_here;
    char* t2;
    char* t3;
    int i;
    if (parser->debug) fprintf(stderr,"Reading multiple binary longs\n");
    keybuf= (key!=NULL)?strdup(key):NULL;
    defbuf= (def!=NULL)?strdup(def):NULL;
    key_here= (keybuf!=NULL)?strtok_r(keybuf,"\\",&t2):NULL;
    def_here= (defbuf!=NULL)?strtok_r(defbuf,"\\",&t3):NULL;
    for (i=0; i<de->length/4; i++) {
      long long val= fixLongSign(FRdInt32(f));
      if (key_here != NULL && strlen(key_here)>0 && strcmp(key_here," ")) {
	kvDefLong(info, key, val);
	if (def_here != NULL) {
	  KVHash* defs= kvGetHash(info,"definitions");
	  if (defs) kvDefString(defs,key_here,def_here);
	}
      }
      else {
	/* Don't define anything for this one */
      }
      if (key_here) key_here= strtok_r(NULL,"\\",&t2);
      if (def_here) def_here= strtok_r(NULL,"\\",&t3);
    }
    free(keybuf);
    free(defbuf);
  }
}

static DCM_METHOD_PROTOTYPE(read_shorts)
{
  if (de->length>0) {
    char* keybuf;
    char* defbuf;
    char* key_here;
    char* def_here;
    char* t2;
    char* t3;
    int i;
    if (parser->debug) fprintf(stderr,"Reading multiple binary shorts\n");
    keybuf= (key!=NULL)?strdup(key):NULL;
    defbuf= (def!=NULL)?strdup(def):NULL;
    key_here= (keybuf!=NULL)?strtok_r(keybuf,"\\",&t2):NULL;
    def_here= (defbuf!=NULL)?strtok_r(defbuf,"\\",&t3):NULL;
    for (i=0; i<de->length/2; i++) {
      int val= fixShortSign(FRdInt16(f));
      if (key_here != NULL && strlen(key_here)>0 && strcmp(key_here," ")) {
	kvDefInt(info,key_here,val);
	if (def_here != NULL) {
	  KVHash* defs= kvGetHash(info,"definitions");
	  if (defs) kvDefString(defs,key_here,def_here);
	}
      }
      else {
	/* Don't define anything for this one */
      }
      if (key_here) key_here= strtok_r(NULL,"\\",&t2);
      if (def_here) def_here= strtok_r(NULL,"\\",&t3);
    }
    free(keybuf);
    free(defbuf);
  }
}

static DCM_METHOD_PROTOTYPE(read_and_map_uid)
{
  if (de->length>0) {
    char* buf;
    const char* name;
    const UID* lookupResult;
    if (parser->debug) 
      fprintf(stderr,
	      "Reading and mapping UID for transfer syntax, length %d!\n",
	      de->length);
    if (!(buf=(char*)malloc((de->length+1)*sizeof(char))))
      Abort("%s: unable to allocate %d bytes!\n",
	    progname,(de->length+1)*sizeof(char));
    (void)fgets(buf,de->length+1,f);
    if (buf[de->length-1]==' ') 
      buf[de->length-1]= '\0';
    lookupResult= dcm_getUIDByUIDString(buf);
    if (lookupResult) {
      if (parser->debug)
	fprintf(stderr,
		"Mapped UID to name <%s> -> <%s>, defined as key <%s>\n",
		buf,lookupResult->name, key);
      kvDefString(info,key,lookupResult->name);
    }
    else {
      if (parser->debug)
	fprintf(stderr,"UID <%s> name is unknown! Raw UID defined as key %s\n",
		buf,key);
      kvDefString(info,key,buf);
    }
    free(buf);
  }
  if (def != NULL) {
    KVHash* defs= kvGetHash(info,"definitions");
    if (defs) kvDefString(defs,key,def);
  }
}

static int dictSearchTest( const void* key, const void* a )
{
  DicomDictDataElement* dictEntry= (DicomDictDataElement*)a;
  DicomDataElement* de= (DicomDataElement*)key;

  if (de->group<dictEntry->group) return -1;
  else if (de->group>dictEntry->group) return 1;
  else {
    if (de->element<dictEntry->element) return -1;
    else if (de->element>dictEntry->element) return 1;
    else return 0;
  }
}

static DicomDictDataElement* dictLookup(const DicomDataElement* de)
{
  DicomDictDataElement* result=
    (DicomDictDataElement*)
    bsearch(de, dicomDict,
	    sizeof(dicomDict)/sizeof(DicomDictDataElement),
	    sizeof(DicomDictDataElement),
	    dictSearchTest);
  if (result==NULL)
    result= (DicomDictDataElement*)
      bsearch(de, dicomSpecialDict,
	      sizeof(dicomSpecialDict)/sizeof(DicomDictDataElement),
	      sizeof(DicomDictDataElement),
	      dictSearchTest);

  return result;
}

static int methodSearchTest( const void* key, const void* a )
{
  DicomElementMethodTableEntry* mthdEntry= (DicomElementMethodTableEntry*)a;
  DicomDataElement* de= (DicomDataElement*)key;

  if (de->group<mthdEntry->group) return -1;
  else if (de->group>mthdEntry->group) return 1;
  else {
    if (de->element<mthdEntry->element) return -1;
    else if (de->element>mthdEntry->element) return 1;
    else return 0;
  }
}

static int methodSortTest( const void* v1, const void* v2 )
{
  DicomElementMethodTableEntry* m1= (DicomElementMethodTableEntry*)v1;
  DicomElementMethodTableEntry* m2= (DicomElementMethodTableEntry*)v2;
  int result;

  if (m1->group<m2->group) result= -1;
  else if (m1->group>m2->group) result= 1;
  else {
    if (m1->element<m2->element) result= -1;
    else if (m1->element>m2->element) result= 1;
    else result= 0;
  }
  return result;
}

static DicomElementMethodTableEntry* 
copyAndSortElementMethodTable( DicomElementMethodTableEntry* table_in,
			       long length)
{
  DicomElementMethodTableEntry* table= NULL;
  if (!(table=(DicomElementMethodTableEntry*)
	malloc(length*sizeof(DicomElementMethodTableEntry))))
    Abort("Unable to allocate %d bytes!\n",
	  length*sizeof(DicomElementMethodTableEntry));
  memcpy(table, table_in, length*sizeof(DicomElementMethodTableEntry));
  qsort(table, length, sizeof(DicomElementMethodTableEntry),
	methodSortTest);
  return table;
}

static DicomElementMethodTableEntry* methodLookup(DicomParser* parser,
						  const DicomDataElement* de)
{
  DicomElementMethodTableEntry* result=
    (DicomElementMethodTableEntry*)
    bsearch(de, parser->elementMethodTable,
	    parser->elementMethodTableSize,
	    sizeof(DicomElementMethodTableEntry),
	    methodSearchTest);
  return result;
}

static void emitDataElementDescription( DicomDataElement* de, FILE* f )
{
  fprintf(stderr,
	  "Element(0x%x, 0x%x), type %s, length %ld, payload offset %lld\n",
	  de->group,de->element, dicomElementTypeNames[(int)de->type],
	  de->length,de->payloadOffset);

  if (bigfile_fseek(f,de->payloadOffset,SEEK_SET))
    fex_raiseException(EXCEPTION_IO,
		       "Could not seek to offset %lld!",de->payloadOffset);
  if (de->dictEntry) {
    switch (de->type) {
    case UI: /* ascii identifier string */
    case SH: /* string */
    case AE: /* Application Entity, probably ascii */
    case AT: /* Some kind of index value, probably ascii */
    case LO: /* ascii organization name */
    case OB: /* no examples available, probably ascii */
    case OW: /* overlay data */
    case OX: /* pixel data */
    case DL: /* no examples available */
    case PN: /* person's name */
    case ST: /* street address? */
    case LT: /* ascii comment */
    case SS: /* no examples available; ascii source string? */
    case SL: /* no examples available, pixel coordinate info, prob ascii */
    case FD: /* no example available, pixel physical coord, prob ascii */
      {
	if (de->length>0) {
	  char* buf;
	  if (!(buf=(char*)malloc((de->length+1)*sizeof(char))))
	    Abort("%s: unable to allocate %d bytes!\n",
		  progname,(de->length+1)*sizeof(char));
	  (void)fgets(buf,de->length+1,f);
	  if (buf[de->length-1]==' ') 
	    buf[de->length-1]= '\0';
	  fprintf(stderr,"Element <%s> type %s: Value is <%s>\n",
		  de->dictEntry->name,
		  dicomElementTypeNames[de->type],buf);
	  free(buf);
	}
	else {
	  fprintf(stderr,"Element <%s> type %s: Value is empty\n",
		  de->dictEntry->name, 
		  dicomElementTypeNames[de->type]);
	}
      }
      break;

    case CS: /* character strings, possibly multiple */
      {
	char* buf;
	if (de->length>0) {
	  char* here;
	  char* t;
	  int i;
	  if (!(buf=(char*)malloc((de->length+1)*sizeof(char))))
	    Abort("%s: unable to allocate %d bytes!\n",
		  progname,(de->length+1)*sizeof(char));
	  (void)fgets(buf,de->length+1,f);
	  if (buf[de->length-1]==' ') 
	    buf[de->length-1]= '\0';
	  fprintf(stderr,
		  "Element name is <%s>, type %s, length %d, val <%s>\n",
		  de->dictEntry->name,
		  dicomElementTypeNames[de->type],
		  de->length,buf);
	  i= 0;
	  here= strtok_r(buf,"\\",&t);
	  while (here != NULL) {
	    fprintf(stderr,"%d: <%s>\n",i,here);
	    here= strtok_r(NULL,"\\",&t);
	    i++;
	  }
	  free(buf);
	}
	else {
	  fprintf(stderr,"Element <%s> type %s: Value is empty\n",
		  de->dictEntry->name, 
		  dicomElementTypeNames[de->type]);
	}
      }
      break;

    case TM: /* ascii time (strange format) */
      {
	if (de->length>0) {
	  char buf[64];
	  char ttime[64];
	  if (de->length<sizeof(buf)-1) {
	    (void)fgets(buf,de->length+1,f);
	    if (buf[de->length-1]==' ') 
	      buf[de->length-1]= '\0';
	    if (strlen(buf)>=6) {
	      sprintf(ttime,"%c%c:%c%c:%s",
		      buf[0],buf[1],
		      buf[2],buf[3],
		      buf+4);
	    }
	    else {
	      Error("%s: dicom_parser: lost <%s>, badly formatted time <%s>\n",
		    progname,de->dictEntry->name,buf);
	    }
	    fprintf(stderr,"Element <%s> type %s: Value is %s\n",
		    de->dictEntry->name,
		    dicomElementTypeNames[de->type],ttime);
	  }
	  else {
	    Error("%s: dicom_parser: lost <%s>, payload too long (%d)!\n",
		  progname,de->dictEntry->name,de->length);
	  }
	}
	else {
	  fprintf(stderr,"Element <%s> type %s: Value is empty\n",
		  de->dictEntry->name, 
		  dicomElementTypeNames[de->type]);
	}
      }
      break;

    case DA: /* ascii date (strange format) */
      {
	if (de->length>0) {
	  char buf[64];
	  char tdate[64];
	  if (de->length<sizeof(buf)-1) {
	    (void)fgets(buf,de->length+1,f);
	    if (buf[de->length-1]==' ') 
	      buf[de->length-1]= '\0';
	    if (strlen(buf)==8) {
	      sprintf(tdate,"%c%c/%c%c/%c%c%c%c",
		      buf[4],buf[5],
		      buf[6],buf[7],
		      buf[0],buf[1],buf[2],buf[3]);
	    }
	    else if (strlen(buf)==7) {
	      sprintf(tdate,"%c%c/0%c/%c%c%c%c",
		      buf[4],buf[5],
		      buf[6],
		      buf[0],buf[1],buf[2],buf[3]);
	    }
	    else {
	      Error("%s: dicom_parser: lost <%s>, badly formatted date <%s>\n",
		    progname,de->dictEntry->name,buf);
	    }
	    fprintf(stderr,"Element <%s> type %s: Value is %s\n",
		    de->dictEntry->name,
		    dicomElementTypeNames[de->type],tdate);
	  }
	  else {
	    Error("%s: dicom_parser: lost <%s>, payload too long (%d)!\n",
		  progname,de->dictEntry->name,de->length);
	  }
	}
	else {
	  fprintf(stderr,"Element <%s> type %s: Value is empty\n",
		  de->dictEntry->name, 
		  dicomElementTypeNames[de->type]);
	}
      }
      break;

    case UL: /* length or offset; binary unsigned long */
      {
	long long val= fixLongSign(FRdInt32(f));
	fprintf(stderr,"Element name is <%s>, type %s, length %d, val %lld\n",
		de->dictEntry->name,
		dicomElementTypeNames[de->type],
		de->length,val);
      }
      break;

    case IS: /* ascii integer */
      {
	char buf[256];
	if (de->length>0) {
	  if (de->length<sizeof(buf)-1) {
	    char* here;
	    char* t;
	    int i;
	    (void)fgets(buf,de->length+1,f);
	    if (buf[de->length-1]==' ') 
	      buf[de->length-1]= '\0';
	    fprintf(stderr,"Element name is <%s>, type %s, length %d, val <%s>\n",
		    de->dictEntry->name,
		    dicomElementTypeNames[de->type],
		    de->length,buf);
	    i= 0;
	    here= strtok_r(buf,"\\ ",&t);
	    while (here != NULL) {
	      long val= atol(here);
	      fprintf(stderr,"%d: %d\n",i,val);
	      here= strtok_r(NULL,"\\ ",&t);
	      i++;
	    }
	  }
	  else {
	    Error("%s: dicom_parser: lost <%s>, payload too long (%d)!\n",
		  progname,de->dictEntry->name,de->length);
	  }
	}
      }
      break;

    case AS: /* ascii integer, age */
      {
	char buf[256];
	long val;
	if (de->length>0) {
	  if (de->length<sizeof(buf)-1) {
	    (void)fgets(buf,de->length+1,f);
	    if (buf[de->length-1]==' ') 
	      buf[de->length-1]= '\0';
	    val= atoi(buf);
	    fprintf(stderr,"Element name is <%s>, type %s, val %d\n",
		    de->dictEntry->name,
		    dicomElementTypeNames[de->type],
		    val);
	  }
	  else {
	    Error("%s: dicom_parser: lost <%s>, payload too long (%d)!\n",
		  progname,de->dictEntry->name,de->length);
	  }
	}
      }
      break;

    case DS: /* ascii float */
      {
	char buf[256];
	if (de->length>0) {
	  if (de->length<sizeof(buf)-1) {
	    char* here;
	    char* t;
	    int i;
	    (void)fgets(buf,de->length+1,f);
	    if (buf[de->length-1]==' ') 
	      buf[de->length-1]= '\0';
	    fprintf(stderr,"Element name is <%s>, type %s, length %d, val <%s>\n",
		    de->dictEntry->name,
		    dicomElementTypeNames[de->type],
		    de->length,buf);
	    i= 0;
	    here= strtok_r(buf,"\\ ",&t);
	    while (here != NULL) {
	      double val= atof(here);
	      fprintf(stderr,"%d: %g\n",i,val);
	      here= strtok_r(NULL,"\\ ",&t);
	      i++;
	    }
	  }
	  else {
	    Error("%s: dicom_parser: lost <%s>, payload too long (%d)!\n",
		  progname,de->dictEntry->name,de->length);
	  }
	}
      }
      break;

    case XS:
    case US: /* unsigned shorts */
      {
	int n;
	int i;
	int val;
	fprintf(stderr,"Element name is <%s>, type %s, length %d\n",
		de->dictEntry->name,
		dicomElementTypeNames[de->type],de->length);	
	n= de->length/2;
	for (i=0; i<n; i++) {
	  val= fixShortSign(FRdInt16(f));
	  fprintf(stderr,"%d: %d\n",i,val);
	}
      }
      break;

    case SQ:
      {
	fprintf(stderr,"Element name is <%s>, type %s\n",
		de->dictEntry->name,
		dicomElementTypeNames[de->type]);
      }
      break;

    default: /* do nothing */
      break;
    }
    
  }
  else fprintf(stderr,"ignoring this unknown element!\n");
}

static void maybeInvokeDataElementMethod(DicomParser* parser,
					 KVHash* info, DicomDataElement* de, 
					 const TransferSyntax* ts, FILE* f)
{
  de->methodEntry= methodLookup(parser,de);
  if (de->methodEntry!=NULL) {
    if (bigfile_fseek(f, de->payloadOffset, SEEK_SET))
      fex_raiseException(EXCEPTION_IO,
			 "Cannot seek to offset %lld",de->payloadOffset);
    if (ts->isLittleEndian)
      bio_big_endian_input= 0;
    else 
      bio_big_endian_input= 1;
    if (de->methodEntry->method != NULL)
      (de->methodEntry->method)(parser,info, de, f, de->methodEntry->key,
				de->methodEntry->def,NULL);
    else {
      /* We have to pick the method based on element type */
      const char* key= de->methodEntry->key;
      const char* def= de->methodEntry->def;
      switch (de->type) {
      case UL:
	read_longs(parser, info, de, f, key, def, NULL);
	break;
      case UI:
	read_and_map_uid(parser, info, de, f, key, def, NULL);
	break;
      case AE:
      case AT:
      case AS:
      case LO:
      case PN:
      case ST:
      case LT:
      case CS:
      case SH:
	read_1_string(parser, info, de, f, key, def, NULL);
	break;
      case US:
      case XS:
	read_shorts(parser, info, de, f, key, def, NULL);
	break;
      case TM:
	read_1_time(parser, info, de, f, key, def, NULL);
	break;
      case DA:
	read_1_date(parser, info, de, f, key, def, NULL);
	break;
      case DS:
	read_ascii_floats(parser, info, de, f, key, def, NULL);
	break;
      case OX:
	store_offset(parser, info, de, f, key, def, NULL);
	break;
      case IS:
	read_ascii_ints(parser, info, de, f, key, def, NULL);
	break;
      case OB:
      case SQ:
	{
	}
	break;
#ifdef never
      case OW:
	{
	}
	break;
      case DL:
	{
	}
	break;
      case SS:
	{
	}
	break;
      case SL:
	{
	}
	break;
      case FD:
	{
	}
	break;
      case DT:
	{
	}
	break;
      case FL:
	{
	}
	break;
      case UT:
	{
	}
	break;
      case OF:
	{
	}
	break;
#endif
      default:
	if (parser->debug) {
	  unsigned char buf[256];
	  short* sbuf= (short*)buf;
	  int i;
	  fprintf(stderr,
		  "This is (0x%04x, 0x%04x), length %ld at %lld\n",
		  de->group, de->element, de->length,de->payloadOffset);
	  FRdUInt8Array (f, buf, sizeof(buf));
	  for (i=0; i<sizeof(buf); i++) {
	    if (!(i%64)) fprintf(stderr,"%03d:  ",i);
	    if (isprint(buf[i])) fputc(buf[i],stderr);
	    else {
	      switch (buf[i]) {
	      case '\0': fprintf(stderr,"~"); break;
	      default: fprintf(stderr,"?"); break;
	      }
	    }
	    if (!((i+1)%64)) fputc('\n',stderr);
	  }
	  for (i=0; i<sizeof(buf); i+=4) {
	    if (!(i%32)) fprintf(stderr,"%03d:  ",i);
	    fprintf(stderr,"%04x %04x ",sbuf[i/2],sbuf[(i/2)+1]);
	    if (!((i+4)%32)) fputc('\n',stderr);
	  }
	}
	fex_raiseException(EXCEPTION_DICOM,
			   "no known handler method for element (%x,%x)!\n",
			   de->group,de->element);
      }
    }
  }
}


static DicomElementType getElementTypeByName( const char vr[2] )
{
  int i;
  for (i=0; i<sizeof(dicomElementTypeNames)/sizeof(char*); i++) {
    if (vr[0]==dicomElementTypeNames[i][0]
	&& vr[1]==dicomElementTypeNames[i][1]) {
      return (DicomElementType)i;
    }
  }
  fex_raiseException(EXCEPTION_DICOM,
		     "encountered unknown DICOM data element type <%c%c>!\n",
		     vr[0],vr[1]);
  return (DicomElementType)0; /* not reached */
}

static long parse_reserved16_length32( FILE* f )
{
  unsigned char reservedBytes[2];
  long result;
  FRdUInt8Array(f,reservedBytes,2);
  if (bio_error) fex_raiseException(EXCEPTION_IO,"premature end of file");
  if (reservedBytes[0]!=0 || reservedBytes[1]!=0) 
    Warning(1,"Unsupported use of reserved bytes in Data Element!\n");
  result= FRdInt32(f);
  if (bio_error) fex_raiseException(EXCEPTION_IO,"premature end of file");
  if (result<0) /* this should be an unsigned int */
    fex_raiseException(EXCEPTION_DICOM,
		       "encountered element with huge length!\n");
  return result;
}

static long long readDataElement(DicomParser* parser, FILE* f, 
				 DicomDataElement* de, long long offset, 
				 const TransferSyntax* transferSyntax)
{
  long long nextOffset;
  char vr[2];

  if (bigfile_fseek(f,offset,SEEK_SET))
    fex_raiseException(EXCEPTION_IO,
		       "Cannot seek to offset %lld!\n",offset);
  if (feof(f)) 
    fex_raiseException(EXCEPTION_IO,"hit EOF at offset %lld",offset);
  if (transferSyntax->isLittleEndian)
    bio_big_endian_input= 0;
  else 
    bio_big_endian_input= 1;
  if (transferSyntax->isExplicitVR) {
    de->group= fixShortSign(FRdInt16(f));
    de->element= fixShortSign(FRdInt16(f));
    if (bio_error) 
      fex_raiseException(EXCEPTION_IO,"premature end of file");
    FRdUInt8Array(f,(unsigned char*)vr,2);
    if (bio_error) 
      fex_raiseException(EXCEPTION_IO,"premature end of file");
    if (parser->debug) fprintf(stderr,"VR is <%c%c>\n",vr[0],vr[1]);
    /* Certain VRs require special handling of length for Explicit reps */
    if (vr[0]=='O') {
      if (vr[1]=='B' || vr[1]=='W') {
	de->length= parse_reserved16_length32( f );
	de->payloadOffset= offset+12;
      }
      else {
	de->length= FRdInt16(f);
	de->payloadOffset= offset+8;
      }
    }
    else if (vr[0]=='S' && vr[1]=='Q') {
      de->length= parse_reserved16_length32( f );
      de->payloadOffset= offset+12;
    }
    else if (vr[0]=='U') {
      if (vr[1]=='N') {
	de->length= parse_reserved16_length32( f );
	de->payloadOffset= offset+12;
      }
      else if (vr[1]=='T') {
	de->length= parse_reserved16_length32( f );
	de->payloadOffset= offset+12;
      }
      else {
	de->length= FRdInt16(f);
	de->payloadOffset= offset+8;
      }
    }
    else {
      de->length= FRdInt16(f);
      de->payloadOffset= offset+8;
    }
    de->dictEntry= dictLookup(de);
    if (de->dictEntry) {
      de->type= de->dictEntry->type;
    }
    else {
      /* Fill in the VR type from the explicit VR data */
      de->type= getElementTypeByName(vr);
    }
  }
  else {
    de->group= fixShortSign(FRdInt16(f));
    de->element= fixShortSign(FRdInt16(f));
    de->length= FRdInt32(f);
    de->payloadOffset= offset+8;
    de->dictEntry= dictLookup(de);
    /************
     * Warning: the order of these classier rules is important!
     ***********/
    if (de->dictEntry) de->type= de->dictEntry->type;
    else if (de->group==0xfffe) de->type= SPECIAL;
    else if (de->length==0xffffffff) de->type= SQ;
    else de->type= UNKNOWN;
  }
  de->methodEntry= NULL;
  
  if (parser->debug) emitDataElementDescription(de,f);

  if (de->length == 0xFFFFFFFF) {
    if (de->type != SQ && de->type != SPECIAL) {
      /* This represents an undefined length; the element should be
       * terminated by a Sequence Delimiter Item.  
       */
      if (de->dictEntry)
	fex_raiseException(EXCEPTION_DICOM,
			   "Element <%s> makes unsupported use of Undefined Length!\n",
			   de->dictEntry->name);
      else 
	fex_raiseException(EXCEPTION_DICOM,
			   "Element <%s> makes unsupported use of Undefined Length!\n",
			   de->dictEntry->name);
    }
  }
  else if (de->length<0) {
    fex_raiseException(EXCEPTION_DICOM,
		       "Element (%x,%x) has negative length at offset %lld",
		       de->group,de->element,offset);
  }

  if (de->type==SQ
      || de->type==SPECIAL && de->element==ITEM_ELEMENT) 
    nextOffset= de->payloadOffset;
  else nextOffset= de->payloadOffset + de->length;

  return nextOffset;
}

void dcm_parseStream(DicomParser* parser, KVHash* info, FILE* f, 
		     long long firstOffset, 
		     const TransferSyntax* transferSyntax,
		     ParserBreakTest breakTest,
		     void* hook)
{
  DicomDataElement de;
  ParserState state= STATE_TOP;
  long long breakOffset= -1;
  SList* stateStack= slist_create();
  long long offset= firstOffset;
  if (parser->debug) 
    fprintf(stderr,"parseDICOMStream starting at offset %lld, syntax <%s>\n",
	    firstOffset,
	    transferSyntax->name);

  if (!dcm_transferSyntaxIsSupported(parser, transferSyntax)) {
    fex_raiseException(EXCEPTION_DICOM,
		       "Transfer syntax <%s> is not supported.",
		       transferSyntax->name);
  }

  /* Since we want de to be on the stack, this stupid language forces
   * us to initialize it by hand.  These values will be seen by 
   * the first breakTest, so they must be valid but innocuous. 
   */
  de.group= 0x0;
  de.element= 0x0;
  de.length= 0;
  de.payloadOffset= 0;
  de.dictEntry= NULL;
  de.methodEntry= NULL;

  while (!feof(f)) {
    long long nextOffset;
    ParserState initialState= state;

    if (ferror(f)) 
      fex_raiseException(EXCEPTION_IO,
			 strerror(errno));
    if (breakTest(info, &de, offset, hook)) break;

    nextOffset= readDataElement(parser,f,&de,offset,transferSyntax);
    maybeInvokeDataElementMethod(parser,info, &de, transferSyntax, f);

    switch (de.type) {
    case SQ: 
      {
	slist_push(stateStack,createStateStackEntry(state, breakOffset));
	state=STATE_SQ;
	if(de.length==0xFFFFFFFF) breakOffset= -1;
	else breakOffset= offset+de.length;
	if (parser->debug) fprintf(stderr,"State pushed onto state stack\n");
      };
      break;
    case SPECIAL:
      {
	/* We are promiscuous here, letting an Item Delim Element
	 * pop a SQ stack element and vice versa.  
	 */
	switch (de.element) {
	case SEQ_DELIM_ITEM_ELEMENT:
	case ITEM_DELIM_ELEMENT:
	  {
	    StateStackEntry* sse= (StateStackEntry*)slist_pop(stateStack);
	    if ((de.element==SEQ_DELIM_ITEM_ELEMENT && state != STATE_SQ)
		|| (de.element==ITEM_DELIM_ELEMENT && state != STATE_ITEM))
	      Warning(1,
		      "DICOM framing error: found <%s> in state %s!\n",
		      de.dictEntry->name, stateNameTable[(int)state]);
	    state= sse->state;
	    breakOffset= sse->breakOffset;
	    destroyStateStackEntry(sse);
	    if (parser->debug) fprintf(stderr,"State stack popped -> %s %lld\n",
			       stateNameTable[(int)state],breakOffset);
	  }
	  break;
	case ITEM_ELEMENT:
	  {
	    slist_push(stateStack,createStateStackEntry(state, breakOffset));
	    state=STATE_ITEM;
	    if(de.length==0xFFFFFFFF) breakOffset= -1;
	    else breakOffset= offset+de.length;
	    if (parser->debug) fprintf(stderr,"State pushed onto state stack\n");
	  }
	  break;
	default:
	  Warning(1,"Skipping over unknown special element (%x,%x) <%s>!\n",
		  de.group, de.element, 
		  ((de.dictEntry!=NULL) ? de.dictEntry->name:"UNKNOWN"));
	}	
      };
      break;
    default:
      {
	while (nextOffset==breakOffset) {
	  StateStackEntry* sse= (StateStackEntry*)slist_pop(stateStack);
	  state= sse->state;
	  breakOffset= sse->breakOffset;
	  destroyStateStackEntry(sse);
	  if (parser->debug) fprintf(stderr,"State stack popped -> %s %lld\n",
			     stateNameTable[(int)state],breakOffset);
	}
      }
      break;
    }

    if (parser->debug && (state != initialState)) 
      fprintf(stderr,"State transition: %s -> %s\n",
	      stateNameTable[(int)initialState],stateNameTable[(int)state]);
    offset= nextOffset;
  }
  slist_destroy(stateStack,destroyStateStackEntry);
}

DicomParser* dcm_createParser(DicomElementMethodTableEntry* table,
			      long tableSize)
{
  DicomParser* result= NULL;
  if (!(result=(DicomParser*)malloc(sizeof(DicomParser))))
    Abort("Unable to allocate %d bytes!\n",sizeof(DicomParser));
  result->debug= 0;
  result->elementMethodTable= copyAndSortElementMethodTable(table, tableSize);
  result->elementMethodTableSize= tableSize;
  return result;
}

void dcm_destroyParser(DicomParser* parser)
{
  if (parser->elementMethodTable) free(parser->elementMethodTable);
  free(parser);
}

void dcm_setDebug(DicomParser* parser, const int val)
{
  parser->debug= val;
}

int dcom_getDebug(DicomParser* parser)
{
  return parser->debug;
}

const char* dcm_getElementTypeName(DicomElementType type)
{
  return dicomElementTypeNames[(int)type];
}
