/************************************************************
 *                                                          *
 *  smart_utils.c                                             *
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
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

static char rcsid[] = "$Id: smart_utils.c,v 1.14 2007/07/03 20:03:19 welling Exp $";

static const char* knownSlicePatternNames[]= {
  "sequential",
  "reversed_sequential",
  "even/odd",
  "reversed_even/odd",
  "odd/even",
  "reversed_odd/even",
  "halves_low_first",
  "reversed_halves_low_first",
  "halves_high_first",
  "reversed_halves_high_first",
  NULL /* must be last */
};

int legitDimString( const char* str )
{
  char* here= (char*)str;
  char* runner;

  while (*here) {
    runner= here+1;
    while (*runner) {
      if (*runner==*here) {
	/* characters at *here and *runner match */
	return 0;
      }
      runner++;
    }
    here++;
  }

  return 1;
}

int consistency_check(KVHash* info)
{
  char* c;
  char* dimstr;
  char buf[256];

  if (!(dimstr=kvGetString(info,"dimstr"))) {
    Error("No information available on data dimension order!\n");
    return 0;
  }
  if (!legitDimString(dimstr)) {
    Error("The dimension string <%s> is inconsistent!\n");
    return 0;
  }
  for (c=dimstr; *c; c++) {
    sprintf(buf,"d%c",*c);
    if (!kvLookup(info,buf)) {
      Error("Missing extent for dimension %c!\n",*c);
      return 0;
    }
    if (kvGetInt(info,buf)<1) {
      Error("Extent of dimension %c is not positive!\n",*c);
      return 0;
    }
    sprintf(buf,"skip.%c",*c);
    if (kvLookup(info,buf))
      if (kvGetLong(info,buf)<0) {
	Error("Skip of dimension %c is negative!\n",*c);
	return 0;
      }
  }

  if (!kvLookup(info,"start_offset")) {
    Error("Initial offset in data is not defined!\n");
    return 0;
  }

  if (kvLookup(info,"start_offset") && kvGetLong(info,"start_offset")<0) {
    Error("%s: start_offset (%lld) must be non-negative!\n",
	  progname,kvGetLong(info,"start_offset"));
    return 0;
  }
  
  if (!kvLookup(info,"chunkname")) {
    Error("%s: Chunk name for output dataset is unknown\n");
    return 0;
  }

  if (!kvLookup(info,"chunkfile")) {
    Error("%s: Chunk file for output dataset is unknown\n");
    return 0;
  }

  return 1;
}

void emit_format_summary(KVHash* info, FileHandler* handler)
{
  char* c;

  /* Print out important parameters */
  if (kvLookup(info,"tag") && (c= kvGetString(info,"tag")) && c[0])
    Message( "#      This data has the scan id <%s>\n",c);
  Message( "#      This data is being read from a %s\n",handler->typeName);
  Message( "#      Input Data-type: %s, Output Data-type: %s\n",
	   srdrTypeName[kvGetInt(info,"datatype_in")], 
	   srdrTypeName[kvGetInt(info,"datatype_out")] );
  if (kvLookup(info,"pulse_seq"))
    Message( "#      Scan was acquired with pulse sequence <%s>\n",
	     kvGetString(info,"pulse_seq"));
  if (kvLookup(info,"plane"))
    Message( "#      Scan plane believed to be %s\n",
	     kvGetString(info,"plane"));
  if (c=kvGetString(info,"dimstr")) {
    char extentbuf[256];
    char* here= extentbuf;
    Message("#      Dataset dimensions are <%s>\n",c);
    extentbuf[255]= '\0';
    for (c=kvGetString(info,"dimstr"); *c; c++) {
      char buf[32];
      sprintf(buf,"d%c",*c);
      snprintf(here,extentbuf+(sizeof(extentbuf)-1)-here,
	       "%d:",kvGetInt(info,buf));
      while (*here) here++;
      if (here-extentbuf >= sizeof(extentbuf)-1) {
	strcpy(extentbuf,"??????????????");
	break;
      }
    }
    *(--here)= '\0'; /* peel off trailing ':' */
    Message("#      extents in order are %s\n",extentbuf);
  }
  if (kvLookup(info,"dx_resampled")!=NULL)
    Message( "#      Expected resampled X resolution: %d\n",
	     kvGetInt(info,"dx_resampled"));
  if (kvLookup(info,"TR") && kvLookup(info,"TE"))
    Message("#      TR= %d us, TE= %d us\n",
	    kvGetInt(info,"TR"), kvGetInt(info,"TE"));
  if (kvLookup(info,"fov_x") && kvLookup(info,"fov_y"))
    Message("#      Field of View X= %g mm, Y= %g mm\n",
	    kvGetDouble(info,"fov_x"), kvGetDouble(info,"fov_y"));
  if (kvLookup(info,"image_x") && kvLookup(info,"image_y"))
    Message("#      Image X= %g voxels, Y= %g voxels\n",
	    kvGetDouble(info,"image_x"), kvGetDouble(info,"image_y"));
  if (kvLookup(info,"overscan")) 
    Message("#      Scan includes %d lines past Ky=0\n",
	    kvGetInt(info,"overscan"));
  if (kvLookup(info,"voxel_x") && kvLookup(info,"voxel_y"))
    Message("#      Voxel X= %g mm, Y= %g mm",
	    kvGetDouble(info,"voxel_x"), kvGetDouble(info,"voxel_y"));
  if (kvLookup(info,"slice_thickness") && kvLookup(info,"slice_gap"))
    Message("\n#      Slice thickness %g mm, gap %g mm\n",
	    kvGetDouble(info,"slice_thickness"), 
	    kvGetDouble(info,"slice_gap"));
  else if (kvLookup(info,"voxel_z"))
    Message(", Z= %g mm\n",kvGetDouble(info,"voxel_z"));
  else Message("\n");
  if (kvLookup(info,"date") && kvLookup(info,"time"))
    Message("#      Scan acquired %s %s\n",
	    kvGetString(info,"date"), kvGetString(info,"time"));
  if (kvGetBoolean(info,"xchop"))
    Message( "#      Data should be phase chopped in X\n");
  if (kvGetBoolean(info,"ychop"))
    Message( "#      Data should be phase chopped in Y\n");
  if (kvGetBoolean(info,"autoscale"))
    Message("#      Will autoscale to approximate range %g\n",
	    kvGetDouble(info,"autoscale_range"));
}

void smart_dump(KVHash* info)
{
  KVIterator* kvi= kvSortedIteratorFactory(info);
  KVHash* defs= NULL;

  if (kvLookup(info,"definitions")) defs= kvGetHash(info,"definitions");
  else defs= kvFactory(3); /* a tiny little fake hash table */

  fprintf(stderr,"Available information is:\n");
  while (kvIteratorHasMorePairs(kvi)) {
    KVPair* p= kvIteratorNextPair(kvi);
    char buf[64];
    char* def= "";
    if (kvLookup(defs,p->key)) def= kvGetString(defs,p->key);
    else if (strstr(p->key,"datatype_")) {
      sprintf(buf,"%d means %s",
	      kvGetInt(info,p->key),srdrTypeName[kvGetInt(info,p->key)]);
      def= buf;
    }
    switch (kvType(p)) {
    case KV_STRING:
      fprintf(stderr,"%21s: %-9s\t%s\n",p->key,p->v.s,def);
      break;
    case KV_DOUBLE:
      fprintf(stderr,"%21s: %-9g\t%s\n",p->key,p->v.d,def);
      break;
    case KV_LONG:
      fprintf(stderr,"%21s: %-9lld\t%s\n",p->key,p->v.l,def);
      break;
    case KV_INT:
      fprintf(stderr,"%21s: %-9d\t%s\n",p->key,(int)(p->v.l),def);
      break;
    case KV_BOOLEAN:
      fprintf(stderr,"%21s: %-9s\t%s\n",p->key,(p->v.l ? "TRUE":"FALSE"),def);
      break;
    case KV_HASH:
      {
        fprintf(stderr,"%21s: **hash**\n",p->key);
#ifdef never
	if (strcmp(p->key,"definitions") && strcmp(p->key,"external_names"))
	  kvDumpTableUnique(p->v.h, stderr, p->key, 0, 0);
#endif
      }
      break;
    }
  }
  fprintf(stderr,"\n");
}

static int parseDelimitedLengthString( KVHash* info, const char* dimstr,
				       const char* prefix, const char* instr,
				       const int allowElided, 
				       const int useLong )
{
  char* dupstr= strdup(instr);
  char* here= dupstr;
  char* str;
  int ns;
  char buf[64];
  int i;

  for (i = 0; i < strlen(dimstr); ++i) {
    char* runner= here;
    while (*runner) {
      if (*runner==',' || *runner==':' || *runner==';' || *runner=='*' || 
	  *runner=='x' || *runner=='-' || *runner==' ' || *runner=='\t' ) 
	break;
      runner++;
    }
    if (*runner) {
      *runner= '\0';
      str= here;
      here= runner+1;
    }
    else {
      if (runner != here) { str= here; here= runner; }
      else str= NULL;
    }
    
    if (str == NULL || str[0]=='\0') {
      if (allowElided) continue; 
      free(dupstr);
      return 0;
    }
    else {
      if (sscanf(str, "%d", &ns) != 1) {
	free(dupstr);
	return 0;
      }
      snprintf(buf,sizeof(buf),"%s%c",prefix,dimstr[i]);
      if (useLong) kvDefLong(info,buf,ns);
      else kvDefInt(info,buf,ns);
    }
  }
  free(dupstr);
  return 1;
}

static int parseSkipString( KVHash* info, 
			    const char* dimstr, const char* skipstr )
{
  return parseDelimitedLengthString( info, dimstr, "skip.", skipstr, 1, 1 );
}

static void assign_skips(KVHash* info)
{
  char* dimstr= NULL;
  char* skipstr= NULL;
  char* cl_skipstr= NULL;
  char* c;
  char buf[256];

  /* From the information given, we make up key-value pairs with
   * keys of the form "skip.x" where x is one of the dimensions
   * in dimstr.  The reading rule is that we skip that many bytes
   * every time we've read in a full block up through x.  For
   * example, if dimstr is "xyzt" and skip.y is 4, 4 bytes will
   * be skipped after each xy slice.
   *
   * We take care to define no skips (even 0) for dimensions 
   * which need no skips, because readers may use a skip of 0
   * to force a break in reading.
   */
  if (!(dimstr= kvGetString(info,"dimstr")))
    return;
  if (kvLookup(info,"skipstr")!=NULL) 
    skipstr= kvGetString(info,"skipstr");
  if (kvLookup(info,"cl_skipstr")!=NULL) 
    skipstr= kvGetString(info,"cl_skipstr");

  /* Parse the archaic 'skip' and 'sliceskip' tags */
  if (strlen(dimstr)>=2) {
    c= dimstr + strlen(dimstr)-2; /* 2nd from last dimension */
    sprintf(buf,"skip.%c",*c);
    if (kvLookup(info,buf)==NULL) {
      if (kvLookup(info,"skip")!=NULL && kvGetLong(info,"skip")!=0) 
	kvDefLong(info, buf, (int)kvGetLong(info,"skip"));
    }
    if (kvLookup(info,"cl_skip")!=NULL && kvGetLong(info,"cl_skip")!=0) 
      kvDefLong(info, buf, (int)kvGetLong(info,"cl_skip"));

    if (c>dimstr) { /* could be only 2 dimensions */
      c--; /* next-to-last dimension */
      sprintf(buf,"skip.%c",*c);
      if (kvLookup(info,buf)==NULL) {
	if (kvLookup(info,"sliceskip")!=NULL 
	    && kvGetLong(info,"sliceskip")!=0) 
	  kvDefLong(info, buf, (int)kvGetLong(info,"sliceskip"));
      }
      if (kvLookup(info,"cl_sliceskip")!=NULL 
	  && kvGetLong(info,"cl_sliceskip")!=0) 
	kvDefLong(info, buf, (int)kvGetLong(info,"cl_sliceskip"));
    }
  }

  /* Parse skipstr, the internally-generated skip string */
  if (kvLookup(info,"skipstr"))
    parseSkipString(info,dimstr,kvGetString(info,"skipstr"));

  /* Parse cl_skipstr, the command line skip string.  It has the final
   * say in assigning skips.
   */
  if (kvLookup(info,"cl_skipstr"))
    parseSkipString(info,dimstr,kvGetString(info,"cl_skipstr"));
}

static int parseExtentString( KVHash* info, 
			      const char* dimstr, const char* extstr )
{
  return parseDelimitedLengthString( info, dimstr, "d", extstr, 0, 0 );
}


int smart_reconcile(KVHash* info)
{
  /* This routine resolves potential conflicts between the multiple keys
   * which may share meanings, including command line settings.
   */

  /* If there is a dimension string from the command line, let it
     override any provided by the file handler. */
  if (kvLookup(info,"cl_dim_string") != NULL)
    kvDefString(info,"dimstr",kvGetString(info,"cl_dim_string"));

  /* If there is an extent string from the command line, we need to 
     apply it to the given dimensions.  */
  if (kvLookup(info,"dimstr") && kvLookup(info,"cl_extent_string") != NULL) {
    if (!parseExtentString(info, kvGetString(info,"dimstr"),
			   kvGetString(info,"cl_extent_string"))) {
      Error("%s: The dimension order <%s> is not consistent with the extents <%s>!\n",
	    progname, ((kvLookup(info,"dimstr")!=NULL)?
		       kvGetString(info,"dimstr") : "UNKNOWN"),
	    kvGetString(info,"cl_extent_string"));
      return 0;
    }
  }

  /* Parse any skip information into its final form.  This may be command
   * line flags or a skip delimiter string set by the handler.
   */
  assign_skips(info);

  /* For the user's convenience, we let the command line setting for
   * autoscale override that from the input file.
   */
  if (kvLookup(info,"cl_autoscale"))
    kvDefBoolean(info,"autoscale",kvGetBoolean(info,"cl_autoscale"));
  if (kvLookup(info,"cl_autoscale_range"))
    kvDefDouble(info,"autoscale_range",kvGetDouble(info,"cl_autoscale_range"));

  /* Sometimes the headers seem to lie about the byte order of 
   * associated data; let the command line override.
   */
  if (kvLookup(info,"cl_big_endian_input"))
    kvDefBoolean(info,"big_endian_input",
		 kvGetBoolean(info,"cl_big_endian_input"));

  /* reordering switch on the command line overrides anything in header */
  if (kvLookup(info,"cl_reorder") != NULL)
    kvDefInt(info,"reorder", (kvGetBoolean(info,"cl_reorder")?1:0));

  /* For automatically-generated chunks there may be no predefined
   * datatype_out.  In this case, do our best to deal with what the
   * handler supplies.
   */
  if (!kvLookup(info,"datatype_out")) {
    int handler_datatype_out= kvGetInt(info,"handler_datatype_out");
    switch (handler_datatype_out) {
    case SRDR_UINT16:
      kvDefInt(info,"datatype_out", SRDR_INT32);
      break;
    default:
      kvDefInt(info,"datatype_out", handler_datatype_out);
    }
    if (verbose_flg) 
      fprintf(stderr,"reconciling chunk %s to type %s\n",
	      kvGetString(info,"chunkname"), 
	      srdrTypeName[kvGetInt(info,"datatype_out")]);
  }

  return 1;
}

