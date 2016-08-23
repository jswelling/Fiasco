/***************************************************************
 *                                                             *
 *  summary.c                                                  *
 *                                                             *
 *  Permission is hereby granted to any individual or          *
 *  institution for use, copying, or redistribution of         *
 *  this code and associated documentation, provided           *
 *  that such code and documentation are not sold for          *
 *  profit and the following copyright notice is retained      *
 *  in the code and documentation:                             *
 *     Copyright (c) 1995 Department of Statistics,            *
 *                        Carnegie Mellon University           *
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
 *                                                             *
 *  Original programming by Mark Fitzgerald and Bill Eddy 2-95 *
 *     6-96: Pittsburgh Format, Mark Fitzgerald                *
 ***************************************************************/
/*************************************************************

  DESCRIPTION OF SUMMARY.C

  summary.c calculates summary statistics for the various
    parameter outputs in FIASCO

  summary.m [-headerinput Input-header-file] [-list Parameter-list-file]
            [-split split-file] [-cond cond-file]
	    [-out Summary-output-file] [-fixed Fixed-image-number]

  Note: summary.m looks for the following environment
        variables for the purpose of printing:
        F_CREDATE  F_PRTDATE  F_HEADER  F_DESCRIPTION

**************************************************************/

#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"

static char rcsid[] = "$Id: summary.c,v 1.39 2007/03/22 00:09:45 welling Exp $";

#define PAGE_TOP 70000
#define LINE_STEP 1200
#define LINES_PER_PAGE 56
#define HEADER_LINES 6
#define MAX_PS_SLICES_PER_BLOCK 8

#define MAX_PARS_PER_FILE 10

typedef enum { PAR_TYPE_UNSET, PAR_TYPE_RAW, PAR_TYPE_FILTERED } ParType;
typedef enum { 
  PAR_ORDER_UNSET,
  PAR_ORDER_Z_FASTEST,
  PAR_ORDER_T_FASTEST,
  PAR_ORDER_Z_ONLY,
  PAR_ORDER_T_ONLY,
  PAR_ORDER_INDEX_T,
  PAR_ORDER_INDEX_Z,
  PAR_ORDER_INDEX_TZ,
  PAR_ORDER_INDEX_ZT,
  PAR_ORDER_ZST /* needed for some multishot data */
} ParOrder;

typedef enum { OUT_PS, OUT_HTML } OutMode;

typedef struct missing_delta_struct {
  char* step;
  int* count;
  struct missing_delta_struct* next;
} MissingDelta;

typedef struct format_info_struct {
  char* fname;
  char* parname;
  ParType type;
  ParOrder order;
  int dz;
  int dt;
  int nFields;
  int nShots;
  char* names[MAX_PARS_PER_FILE];
  struct format_info_struct *next;
} FormatInfo;

static long dt, dz;
static float ***pars = NULL;
static double **mean = NULL, **stdv = NULL, **min = NULL, **max = NULL;
static unsigned char **missing = NULL;
static char* progname;
static int* counts= NULL;
static int* totals= NULL;
static int* missing_totals= NULL;
static int nconditions= 0;
static sp_ConditionDef** cond_table= NULL;
static FormatInfo* formatsHead= NULL;
static FormatInfo* formatsTail= NULL;
static OutMode outMode= OUT_PS;
static long psLocation;
static long psIncrement;

static void stats();

static char* trimString( char* s, int n )
{
  static char buf[256];
  int nToCopy;
  nToCopy= strlen(s);
  if (nToCopy>n) nToCopy= n;
  if (nToCopy>255) nToCopy= 255;
  strncpy(buf,s,nToCopy);
  buf[nToCopy]= '\0';
  return buf;
}

static void readParLine(FILE* ifp, char* fname, int nskip, int n, 
			int t, int z, int parOffset)
{
  char buf[512];
  int i;
  int foundData= 0;
  char* where= NULL;
  char* tok;

  while (!foundData) {
    if (feof(ifp)) 
      Abort("Out of data reading image %d, slice %d in file %s!\n",
	    t,z,fname);
    /* Read, ignoring comments */
    if (fgets(buf, sizeof(buf), ifp) && !ferror(ifp)
	&& strlen(buf)>0 && buf[0] != '#') {
      foundData= 1;
    }
    else {
      if (ferror(ifp)) {
	perror("Error reading parameters");
	Abort("Error occured reading image %d, slice %d from file %s!\n",
	      t,z,fname);
      }
    }
  }
  
  where= NULL;
  for (i=0; i<n+nskip; i++) {
    tok= strtok_r(((i==0) ? buf : 0)," \t;,",&where);
    if (!tok) 
      Abort("Line too short reading image %d, slice %d from file %s!\n",
	    t,z,fname);
    if (i>=nskip) pars[parOffset+i-nskip][t][z]= atof(tok);
  }
}

static FormatInfo* createFormatInfo(char* parname, char* fname, 
				    int dz, int dt) {
  FormatInfo* result= (FormatInfo*)malloc(sizeof(FormatInfo));
  if (!result) Abort("Unable to allocate %d bytes!\n");

  result->parname= strdup(parname);
  result->fname= strdup(fname);
  result->type= PAR_TYPE_UNSET;
  result->order= PAR_ORDER_UNSET;
  result->dz= dz;
  result->dt= dt;
  result->nFields= 0;
  result->nShots= 1;
  result->next= NULL;
  return result;
}

static void destroyFormatInfo( FormatInfo* fi ) {
  int i;
  for (i=0; i<fi->nFields; i++) free(fi->names[i]);
  free(fi->fname);
  free(fi->parname);
  free(fi);
}

static void addFieldName( FormatInfo* fi, char* name )
{
  if (fi->nFields >= MAX_PARS_PER_FILE - 1) 
    Abort("A parameter file has more than %d fields!\n",MAX_PARS_PER_FILE);
  fi->names[ fi->nFields ]= strdup(name);
  fi->nFields++;
}

static char* getParName( FormatInfo* fi )
{
  return fi->parname;
}

static char* getTrimmedParName( FormatInfo* fi, int n )
{
  static char buf[256];
  int nToCopy;
  char* str= getParName(fi);
  nToCopy= strlen(str);
  if (nToCopy>n) nToCopy= n;
  if (nToCopy>255) nToCopy= 255;
  strncpy(buf,str,nToCopy);
  buf[nToCopy]= '\0';
  return buf;
}

static char* getFileName( FormatInfo* fi )
{
  return fi->fname;
}

static char* getFieldName( FormatInfo* fi, int i )
{
  if (i<fi->nFields) return fi->names[i];
  else {
    fprintf(stderr,"internal error: field name index out of bounds!\n");
    return NULL;
  }
}

static char* getTrimmedFieldName( FormatInfo* fi, int i, int n )
{
  static char buf[256];
  int nToCopy;
  char* str= getFieldName(fi,i);
  nToCopy= strlen(str);
  if (nToCopy>n) nToCopy= n;
  if (nToCopy>255) nToCopy= 255;
  strncpy(buf,str,nToCopy);
  buf[nToCopy]= '\0';
  return buf;
}

static int formatValid( FormatInfo* fi )
{
  if ((fi->order != PAR_ORDER_UNSET)
      && (fi->dz>0) && (fi->dt>0) && (fi->nFields>0)) return 1;
  else return 0;
}

static void updateFormat( FormatInfo* fi )
{
  FILE* f;
  char line[512];

  f= efopen( fi->fname, "r");
  while (!feof(f) && !ferror(f)) {
    if (fgets(line, sizeof(line), f)) {
      line[strlen(line)-1]= '\0'; /* get rid of trailing linefeed */
      if (!strncasecmp(line,"##Format:",9)) { 
	/* Found something interesting! */
	char* where;
	char* tok;

	tok= strtok_r(line+9," \t,;",&where);
	while (tok) {
	  if (*tok == ')') {
	    Error("Error scanning %s: misplaced ')'!\n",fi->fname);
	  }
	    
	  if (!strncasecmp(tok,"order:",6)) {
	    tok += 6;
	    if (!strncasecmp(tok,"z_fastest",9)) 
	      fi->order= PAR_ORDER_Z_FASTEST;
	    else if (!strncasecmp(tok,"t_fastest",9)) 
	      fi->order= PAR_ORDER_T_FASTEST;
	    else if (!strncasecmp(tok,"z_only",6)) 
	      fi->order= PAR_ORDER_Z_ONLY;
	    else if (!strncasecmp(tok,"t_only",6)) 
	      fi->order= PAR_ORDER_T_ONLY;
	    else if (!strncasecmp(tok,"index_tz",8)) 
	      fi->order= PAR_ORDER_INDEX_TZ;
	    else if (!strncasecmp(tok,"index_zt",8)) 
	      fi->order= PAR_ORDER_INDEX_ZT;
	    else if (!strncasecmp(tok,"index_t",7)) 
	      fi->order= PAR_ORDER_INDEX_T;
	    else if (!strncasecmp(tok,"index_z",7)) 
	      fi->order= PAR_ORDER_INDEX_Z;
	    else 
	      Error("Error scanning %s: invalid order %s!\n",
		    fi->fname,tok);
	  }
	  else if (!strncasecmp(tok,"type:",5)) {
	    tok += 5;
	    if (!strncasecmp(tok,"raw",3)) fi->type= PAR_TYPE_RAW;
	    else if (!strncasecmp(tok,"filtered",8)) 
	      fi->type= PAR_TYPE_FILTERED;
	    else Error("Error scanning %s: invalid type %s\n",
		       fi->fname,tok);
	  }
	  else if (!strncasecmp(tok,"names:(",7)) {
	    tok += 7;
	    fi->nFields= 0;
	    while ((tok != NULL) && strncmp(tok,")",1)) {
	      if (strlen(tok)>0) {
		char* names_end;
		if (names_end=strchr(tok,')')) {
		  *names_end= '\0';
		  addFieldName(fi,tok);
		  break;
		}
		else addFieldName(fi,tok);
	      }
	      tok= strtok_r(NULL," \t,;",&where);
	    }
	  }
	  tok= strtok_r(NULL," \t,;",&where);
	}
      }
    }
    if (ferror(f)) {
      Error("Error reading from <%s>!\n",fi->fname);
      break;
    }
  }
  fclose(f);
  
}

static void load( FormatInfo* fi )
     /* This loads to the global buffer "pars" */
{
  FILE* fp;
  int t;
  int z;
  int i;

  if (!formatValid(fi)) Abort("Invalid format info for file %s!\n",
			      getFileName(fi));

  fp= efopen( getFileName(fi), "r" );

  switch (fi->order) {
  case PAR_ORDER_UNSET:
    {
      Abort("Invalid format info for file %s!\n",getFileName(fi));
    }
    break;
  case PAR_ORDER_Z_FASTEST:
    {
      for (t=0; t<fi->dt; t++) 
	for (z=0; z<fi->dz; z++) {
	  readParLine(fp, getFileName(fi), 0, fi->nFields, t, z, 0);
	}
    }
    break;
  case PAR_ORDER_T_FASTEST:
    {
      for (z=0; z<fi->dz; z++) {
	for (t=0; t<fi->dt; t++) 
	  readParLine(fp, getFileName(fi), 0, fi->nFields, t, z, 0);
      }
    }
    break;
  case PAR_ORDER_Z_ONLY:
    {
      for (z=0; z<fi->dz; z++) {
	readParLine(fp, getFileName(fi), 0, fi->nFields, 0, z, 0);
	for (i=0; i<fi->nFields; i++) 
	  for( t = 1; t < fi->dt; t++ )
	    pars[i][t][z]= pars[i][0][z];
      }
    }
    break;
  case PAR_ORDER_T_ONLY:
    {
      for (t=0; t<fi->dt; t++) {
	readParLine(fp, getFileName(fi), 0, fi->nFields, t, 0, 0);
	for (i=0; i<fi->nFields; i++) 
	  for( z = 1; z < fi->dz; z++ )
	    pars[i][t][z]= pars[i][t][0];
      }
    }
    break;
  case PAR_ORDER_INDEX_T:
    {
      for (t=0; t<fi->dt; t++) {
	readParLine(fp, getFileName(fi), 1, fi->nFields, t, 0, 0);
	for (i=0; i<fi->nFields; i++) 
	  for( z = 1; z < fi->dz; z++ )
	    pars[i][t][z]= pars[i][t][0];
      }
    }
    break;
  case PAR_ORDER_INDEX_Z:
    {
      for (z=0; z<fi->dz; z++) {
	readParLine(fp, getFileName(fi), 1, fi->nFields, 0, z, 0);
	for (i=0; i<fi->nFields; i++) 
	  for( t = 1; t < fi->dt; t++ )
	    pars[i][t][z]= pars[i][0][z];
      }
    }
    break;
  case PAR_ORDER_INDEX_TZ:
    {
      for (t=0; t<fi->dt; t++) 
	for (z=0; z<fi->dz; z++)
	  readParLine(fp, getFileName(fi), 2, fi->nFields, t, z, 0);
    }
    break;
  case PAR_ORDER_INDEX_ZT:
    {
      for (z=0; z<fi->dz; z++) {
	for (t=0; t<fi->dt; t++) 
	  readParLine(fp, getFileName(fi), 2, fi->nFields, t, z, 0);
      }
    }
    break;
  case PAR_ORDER_ZST: /* useful for some multishot data */
    {
      int shotFields= fi->nFields/2;
      int shot;
      for (t=0; t<fi->dt; t++) {
	for (shot= 0; shot<fi->nShots; shot++) {
	  for (z=0; z<fi->dz; z++) {
	    readParLine(fp, getFileName(fi), 2, shotFields, 
			t, z, shot*shotFields);
	  }
	}
      }
    }
    break;

  }

  fclose(fp);
}

static void init_postscript( FILE* ofp )
{
  fprintf( ofp, "%%!PS-Adobe-1.0 \n" );
  fprintf( ofp, "1 100 div dup scale \n" );
  fprintf( ofp, "/BodyF /Courier findfont 1100 scalefont def \n" );
  fprintf( ofp, "/CW BodyF setfont ( ) stringwidth pop def \n" );
  fprintf( ofp, "/L         { CW mul add exch moveto show } def \n" );
  fprintf( ofp, "/T         { exch moveto show } def \n" );
  fprintf( ofp, "/StartPage { /SavedPage save def \n" );
  fprintf( ofp, "  BodyF setfont 0 setgray } def \n" );
  fprintf( ofp, "/EndPage   {showpage SavedPage restore } def \n" );
  fprintf( ofp, "/outsidecircletext \n" );
  fprintf( ofp, " { circtextdict begin \n" );
  fprintf( ofp, "  /radius exch def /centerangle exch def \n" );
  fprintf( ofp, "  /ptsize exch def /str exch def \n" );
  fprintf( ofp, "  /xradius radius ptsize 4 div add def \n" );
  fprintf( ofp, "  gsave \n" );
  fprintf( ofp, "   centerangle str findhalfangle add rotate \n" );
  fprintf( ofp, "   str \n" );
  fprintf( ofp, "    {/charcode exch def \n" );
  fprintf( ofp, "     (A) dup 0 charcode put outsideplacechar \n" );
  fprintf( ofp, "   } forall \n" );
  fprintf( ofp, "   grestore \n" );
  fprintf( ofp, "   end \n" );
  fprintf( ofp, "   } def \n" );
  fprintf( ofp, "/circtextdict 16 dict def \n" );
  fprintf( ofp, "circtextdict begin \n" );
  fprintf( ofp, " /findhalfangle \n" );
  fprintf( ofp, "   {stringwidth pop 2 div \n" );
  fprintf( ofp, "    2 xradius mul pi mul div 360 mul \n" );
  fprintf( ofp, "   } def \n" );
  fprintf( ofp, "/outsideplacechar \n" );
  fprintf( ofp, "  {/char exch def \n" );
  fprintf( ofp, "   /halfangle char findhalfangle def \n" );
  fprintf( ofp, "   gsave \n" );
  fprintf( ofp, "    halfangle neg rotate \n" );
  fprintf( ofp, "    radius 0 translate \n" );
  fprintf( ofp, "    -90 rotate \n" );
  fprintf( ofp, "    char stringwidth pop 2 div neg 0 moveto \n" );
  fprintf( ofp, "    char show \n" );
  fprintf( ofp, "   grestore \n" );
  fprintf( ofp, "   halfangle 2 mul neg rotate \n" );
  fprintf( ofp, "  } def \n" );
  fprintf( ofp, " /pi 3.1415923 def \n" );
  fprintf( ofp, " end \n\n" );
}

static void init_output( FILE* ofp )
{
  switch (outMode) {
  case OUT_PS: 
    {
      init_postscript(ofp);
    }
    break;
  case OUT_HTML:
    {
      char* hdr;

      fprintf(ofp,"<!DOCTYPE HTML PUBLIC \"-//IETF//DTD HTML//EN\">\n");
      fprintf(ofp,"<html>\n");
      fprintf(ofp,"<head>\n");
      if (hdr=getenv("F_HEADER")) {
	fprintf(ofp,"<title>Fiasco Summary for %s</title>\n",hdr);
	fprintf(ofp,"</head>\n");
	fprintf(ofp,"<body>\n");
	fprintf(ofp,"<center><h1>Fiasco Summary for %s</h1></center>\n",hdr);
      }
      else {
	fprintf(ofp,"<title>Fiasco Summary</title>\n");
	fprintf(ofp,"</head>\n");
	fprintf(ofp,"<body>\n");
	fprintf(ofp,"<center><h1>Fiasco Summary</h1></center>\n");
      }
    }
    break;
  }
}

static void end_output( FILE* ofp )
{
  switch (outMode) {
  case OUT_PS: 
    {
      fprintf( ofp, "EndPage \n" );
    }
    break;
  case OUT_HTML:
    {
      char* credate= NULL;
      credate= getenv("F_CREDATE");
      if (credate) 
	fprintf(ofp,"<!-- Created: %s -->\n",credate);
      fprintf( ofp,"<br><font size=-1>Automatically generated by %s</font>\n",
	       rcsid );
      fprintf( ofp,"</body>\n");
      fprintf( ofp,"</html>\n");
    }
    break;
  }
}

static void step_line(int n)
{
  switch (outMode) {
  case OUT_PS:
    {
      psLocation -= n * psIncrement;
    }
    break;
  case OUT_HTML:
    {
    }
    break;
  }
}

static void out_line( FILE* ofp, int size, char* psFormat, char* htmlFormat,
			... )
{
  va_list args;

  va_start(args, htmlFormat);
  switch (outMode) {
  case OUT_PS:
    { 
      if (!strncmp(psFormat,"%%%%",4)) {
	/* This is meant to be a Postscript comment */
	vfprintf(ofp, psFormat, args);
	fprintf( ofp," %ld 500 T\n", psLocation);
      }
      else {
	fprintf( ofp,"(");
	vfprintf(ofp, psFormat, args);
	fprintf( ofp,") %ld 500 T\n", psLocation);
      }
    }
    break;
  case OUT_HTML:
    {
      if (size != 0) fprintf(ofp,"<font size=%c%d>",(size>0)?'+':'-',size);
      vfprintf(ofp, htmlFormat, args);
      if (size != 0) fprintf(ofp,"</font>");
      fprintf(ofp,"\n");
    }
    break;
  }
  va_end(args);
  step_line(1);
}

static void emit_header( FILE* ofp )
{
  char* credate;
  char* prtdate;

  credate= getenv("F_CREDATE");
  prtdate= getenv("F_PRTDATE");
  if (credate) 
    out_line( ofp, 1, "Created: %s", 
	      "<center>Created: %s</center>", credate );
  else
    out_line( ofp, 1, "%%%%CREATIONMSG%%%%", 
	      "<center>Creation date unknown</center>");
  if (prtdate) 
    out_line( ofp, 1, "Printed: %s", 
	      "<center>Printed: %s</center>", prtdate );
  else
    out_line( ofp, 1, "%%%%PRINTMSG%%%%", 
	      "<center>Print date unknown</center>");
  out_line( ofp, 1, "Study Name: %s",
	    "<center>Study Name: %s</center>",
	    getenv("F_HEADER") );
  out_line( ofp, 1, "Study Description: %s",
	    "<center>Study Description: %s</center>",
	    getenv( "F_DESCRIPTION" ) );
  out_line( ofp, 1, "User: %s", 
	    "<center>User: %s</center>", 
	    getenv( "LOGNAME" ) );
  out_line( ofp, 1, "Path: %s", 
	    "<center>Path: %s</center>", 
	    getenv( "PWD" ) );
  out_line( ofp, 1, "Subject Age: %s", 
	    "<center>Subject Age: %s</center>", 
	    getenv( "F_SUBJ_AGE" ) );
  out_line( ofp, 1, "Subject Gender: %s", 
	    "<center>Subject Gender: %s</center>", 
	    getenv( "F_SUBJ_SEX" ) );
  out_line( ofp, 1, "Subject Diagnosis 1: %s", 
	    "<center>Subject Diagnosis 1: %s</center>", 
	    getenv( "F_SUBJ_1DIAG" ) );
  out_line( ofp, 1, "Subject Diagnosis 2: %s", 
	    "<center>Subject Diagnosis 2: %s</center>", 
	    getenv( "F_SUBJ_2DIAG" ) );
  out_line( ofp, 1, "Subject Diagnosis 3: %s", 
	    "<center>Subject Diagnosis 3: %s</center><p>", 
	    getenv( "F_SUBJ_3DIAG" ) );
}

static void start_table( FILE* ofp, int size, char* titleFormat, ...)
{
  va_list args;
  static int tableNum= 0;

  va_start(args, titleFormat);
  switch (outMode) {
  case OUT_PS: 
    {
      fprintf( ofp,"(");
      vfprintf(ofp, titleFormat, args);
      fprintf( ofp,") %ld 500 T\n", psLocation);
      step_line(2);
    }
    break;
  case OUT_HTML:
    {
      char buf[512];
      fprintf(ofp,"<a name=\"TABLE%d\">\n",tableNum++);
      strcpy(buf,"<caption><b>");
      if (size != 0) 
	sprintf(buf,"%s<font size=%c%d>",buf,(size>0)?'+':'-',size);
      strncat(buf,titleFormat,450-strlen(buf));
      if (size != 0) sprintf(buf,"%s</font>",buf);
      strcat(buf,"</b></caption>\n");
      fprintf(ofp,"<table border cellpadding=5 width=100%% bgcolor=white>\n");
      vfprintf(ofp,buf,args);
      
    }
    break;
  }
  va_end(args);
}

static void end_table( FILE* ofp )
{
  switch (outMode) {
  case OUT_PS: 
    {
      step_line(1);
    }
    break;
  case OUT_HTML:
    {
      fprintf(ofp,"</table><p>\n");
    }
    break;
  }
}

static void start_page( FILE* ofp )
{
  static int pageNum= 0;
      
  switch (outMode) {
  case OUT_PS: 
    {
      if (pageNum>0) fprintf( ofp, "EndPage \n" );
      pageNum++;
      fprintf( ofp, "%%%%Page: %ld %ld \n", pageNum, pageNum );
      fprintf( ofp, "(%ld) (summary.ps) StartPage \n", pageNum);
      fprintf( ofp, "/Helvetica findfont 1500 scalefont setfont \n" );
      fprintf( ofp, "(%s) 5000 5000 1 L\n", getenv( "FIASCO_VERSION" ) );
      fprintf( ofp, "/Helvetica findfont 2000 scalefont setfont \n" );
      fprintf( ofp, "6600 5000 translate \n" );
      fprintf( ofp, "(FIASCO) \n" );
      fprintf( ofp, "2000 90 2000 outsidecircletext \n" );
      fprintf( ofp, "/Courier findfont 1100 scalefont setfont \n" );
      psLocation = PAGE_TOP;
      psIncrement = LINE_STEP;
      emit_header(ofp);
    }
    break;
  case OUT_HTML:
    {
      /* Set some valid values */
      psLocation = PAGE_TOP;
      psIncrement = LINE_STEP;
      if (pageNum==0) emit_header(ofp);
      pageNum++;
      fprintf(ofp,"</a>\n"); /* closes block setting name of table */
    }
    break;
  }
}

static int page_full(int lines_to_come)
{
  switch (outMode) {
  case OUT_PS:
    {
      if (lines_to_come>LINES_PER_PAGE) {
	return 1; /* no hope, but start on fresh page */
      }
      if ((PAGE_TOP - psLocation) 
	  >= ((LINES_PER_PAGE-lines_to_come)*LINE_STEP)) return 1;
      else return 0;
    }
    /* break; -NOTREACHED- */
  case OUT_HTML:
    {
      return 0;
    }
    /* break; -NOTREACHED- */
  }
  return 0;
}

static void emit_contents( FILE* ofp, int env_flag )
{
  /* Write a table of contents if the medium supports it */
  switch (outMode) {
  case OUT_PS:
    {
      /* I don't know of a way to do this */
    }
    break;
  case OUT_HTML:
    {
      int i;
      int tblNum;
      fprintf(ofp,"<h2>Contents:</h2>\n");
      fprintf(ofp,"<ul>\n");
      tblNum= 1;
      for (i=0; i<dz; i++) {
	fprintf(ofp,"<li><a href=\"#TABLE%d\">Slice %d</a>\n",tblNum++,i);
      }
      if (env_flag) {
	fprintf(ofp,"<li><a href=\"#TABLE%d\">Fiasco Environment</a>\n",
		tblNum++);
      }
      fprintf(ofp,"</ul>\n");
      fprintf(ofp,"<p>\n");
    }
    break;
  }
}

static void fake_slice_counts( int nslices, int nimages )
{
  int z;
  int t;
  int cond;

  nconditions= 2; /* NA and one for everything else */

  if (!(counts= (int*)malloc(nconditions*nslices*sizeof(int))))
    Abort("%s: unable to allocate %d ints!\n", progname, 
	  nconditions*nslices);
  if (!(missing_totals= (int*)malloc(nslices*sizeof(int))))
    Abort("%s: unable to allocate %d ints!\n", progname, nslices);
  if (!(totals= (int*)malloc(nslices*sizeof(int))))
    Abort("%s: unable to allocate %d ints!\n", progname, nslices);
  if (!(cond_table= 
	(sp_ConditionDef**)malloc(nconditions*sizeof(sp_ConditionDef*))))
    Abort("%s: unable to allocate %d bytes!\n", progname,
	  nconditions*sizeof(sp_ConditionDef*));

  for (cond=0; cond<nconditions; cond++) {
    if (!(cond_table[cond]=(sp_ConditionDef*)malloc(sizeof(sp_ConditionDef))))
      Abort("%s: unable to allocate %d bytes!\n", progname,
	    nconditions*sizeof(sp_ConditionDef*));
    cond_table[cond]->id= 0;
    cond_table[cond]->name= ((cond==1) ? "All Conditions" : "NA");
    cond_table[cond]->factor_lvl= NULL;
    cond_table[cond]->factor_intlvl= NULL;
  }

  for (z=0; z<nslices; z++) {
    missing_totals[z]= 0;
    for (t=0; t<nimages; t++) if (missing[t][z]) missing_totals[z]++;
  }

  /* All the counts go into our "all conditions" condition. */
  for (cond=0; cond<nconditions; cond++) {
    for (z=0; z<nslices; z++) counts[(cond*nslices)+z]= 0;
  }
  for (t=0; t<nimages; t++)
    for (z=0; z<nslices; z++) {
      if (!missing[t][z]) counts[ (1*nslices)+z ]++;
    }

  for (z=0; z<nslices; z++) {
    totals[z]= 0;
    for (cond=0; cond<nconditions; cond++) 
      totals[z] += counts[(cond*nslices)+z];
  }
  for (z=0; z<nslices; z++) totals[z] += missing_totals[z];
}

static void do_slice_counts(char* condfile, char* splitfile,
			    int nslices, int nimages)
{
  char** factor_table;
  int nfactors;
  sp_SplitRec* split_table;
  int z;
  int t;
  int cond;

  if (!sp_parse_conditions( condfile, &factor_table, &nfactors,
                            &cond_table, &nconditions ))
    Abort("%s: fatal error parsing condition file %s!\n",
	  progname,condfile);

  /* Parse split file */
  if (!sp_parse_split( splitfile, nslices, nimages, nconditions,
                       &split_table ))
    Abort("%s: error processing split file <%s>!\n",
          progname,splitfile);
  
  if (!(counts= (int*)malloc(nconditions*nslices*sizeof(int))))
    Abort("%s: unable to allocate %d ints!\n", progname, 
	  nconditions*nslices);
  if (!(missing_totals= (int*)malloc(nslices*sizeof(int))))
    Abort("%s: unable to allocate %d ints!\n", progname, nslices);
  if (!(totals= (int*)malloc(nslices*sizeof(int))))
    Abort("%s: unable to allocate %d ints!\n", progname, nslices);

  for (z=0; z<nslices; z++) {
    missing_totals[z]= 0;
    for (t=0; t<nimages; t++) if (missing[t][z]) missing_totals[z]++;
  }

  for (cond=0; cond<nconditions; cond++) {
    for (z=0; z<nslices; z++) counts[(cond*nslices)+z]= 0;
  }
  for (t=0; t<nimages; t++)
    for (z=0; z<nslices; z++) {
      sp_SplitRec* s= split_table+((t*nslices)+z);
      if (!missing[t][z]) counts[ (s->cond*nslices) + s->slice ]++;
    }

  for (z=0; z<nslices; z++) {
    totals[z]= 0;
    for (cond=0; cond<nconditions; cond++) 
      totals[z] += counts[(cond*nslices)+z];
  }
  for (z=0; z<nslices; z++) totals[z] += missing_totals[z];

}

static void emit_count_table( FILE* ofp, MissingDelta* missing_deltas )
{
  char* buf;
  int buf_size;
  int z;
  int cond;
  int slices_this_block;
  int slices_so_far;
  int max_slices_per_block= MAX_PS_SLICES_PER_BLOCK;

  switch (outMode) {
  case OUT_PS:
    max_slices_per_block= MAX_PS_SLICES_PER_BLOCK;
    break;
  case OUT_HTML:
    max_slices_per_block= dz;
    break;
  }

  /* Output buffer can get pretty long for big HTML tables */
  buf_size= 22*max_slices_per_block;
  if (buf_size<512) buf_size= 512;
  if (!(buf=(char*)malloc(buf_size*sizeof(char))))
    Abort("%s: unable to allocate %d bytes!\n",progname,buf_size);

  start_table(ofp, 2, "Image Counts by Condition and Slice:");
  step_line(1);

  slices_so_far= 0;
  while (slices_so_far < dz) {
    slices_this_block= (dz-slices_so_far > max_slices_per_block) ?
      max_slices_per_block : dz - slices_so_far;

    switch (outMode) {
    case OUT_PS:
      {
	sprintf(buf, "                      ");
	for (z=0; z<slices_this_block; z++) 
	  sprintf(buf,"%s -%2d- ",buf,z+slices_so_far);
      }
      break;
    case OUT_HTML:
      {
	sprintf(buf,"<td></td><td></td>");
	for (z=0; z<slices_this_block; z++)
	  sprintf(buf,"%s <td><b>%2d</b></td> ",buf,z+slices_so_far);
      }
      break;
    }
    out_line(ofp, 0, "%s", "<tr>%s</tr>", buf);
    step_line(1);
    for (cond=0; cond<nconditions; cond++) {
      if (page_full(0)) {
	start_page(ofp);
	step_line(1);
      }
      switch (outMode) {
      case OUT_PS:
	{
	  sprintf(buf, "%2d %-18s", cond, cond_table[cond]->name);
	  for (z=0; z<slices_this_block; z++)
	    sprintf(buf, "%s  %4d", buf, counts[(cond*dz)+(z+slices_so_far)]);
	}
	break;
      case OUT_HTML:
	{
	  sprintf(buf, "<td><b>%2d</b></td> <td><b>%s</b></td> ", 
		  cond, cond_table[cond]->name);
	  for (z=0; z<slices_this_block; z++)
	    sprintf(buf,"%s <td>%4d</td>",
		    buf, counts[(cond*dz)+(z+slices_so_far)]);
	}
	break;
      }
      out_line( ofp, 0, "%s", "<tr>%s</tr>", buf );
    }
    if (missing_deltas) { /* do step-by-step missing info */
      MissingDelta* tmp= missing_deltas;
      while (tmp != NULL) {
	if (page_full(2)) {
	  start_page(ofp);
	  step_line(2);
	}
	switch (outMode) {
	case OUT_PS:
	  {
	    sprintf(buf,"   drop by %-10s ",trimString(tmp->step,10));
	    for (z=0; z<slices_this_block; z++)
	      sprintf(buf, "%s %4d ", buf, tmp->count[z+slices_so_far]);
	  }
	  break;
	case OUT_HTML:
	  {
	    sprintf(buf,"<td></td><td><b>drop by %-10s</b></td>",tmp->step);
	    for (z=0; z<slices_this_block; z++)
	      sprintf(buf, "%s <td>%4d</td>", 
		      buf, tmp->count[z+slices_so_far]);
	  }
	  break;
	}
	out_line(ofp, 0, "%s", "<tr>%s</tr>", buf);
	tmp= tmp->next;
      }
    }
    else { /* only one missing entry */
      if (page_full(2)) {
	start_page(ofp);
	step_line(2);
      }
      switch (outMode) {
      case OUT_PS:
	{
	  sprintf(buf,"   dropped by Fiasco  ");
	  for (z=0; z<slices_this_block; z++)
	    sprintf(buf, "%s %4d ", buf, missing_totals[z+slices_so_far]);
	}
	break;
      case OUT_HTML:
	{
	  sprintf(buf,"<td></td><td><b>dropped by Fiasco</b></td> ");
	  for (z=0; z<slices_this_block; z++)
	    sprintf(buf, "%s <tr>%4d</tr>", 
		    buf, missing_totals[z+slices_so_far]);
	}
	break;
      }
      out_line(ofp, 0, "%s", "<tr>%s</tr>", buf);
    }
    switch (outMode) {
    case OUT_PS:
      {
	sprintf(buf,"   totals             ");
	for (z=0; z<slices_this_block; z++)
	  sprintf(buf, "%s %4d ", buf, totals[z+slices_so_far]);
      }
      break;
    case OUT_HTML:
      {
	sprintf(buf,"<td></td><td><b>totals</b></td>");
	for (z=0; z<slices_this_block; z++)
	  sprintf(buf, "%s <td>%4d</td>", buf, totals[z+slices_so_far]);
      }
      break;
    }
    out_line(ofp, 0, "%s", "<tr>%s</tr>", buf);

    slices_so_far += slices_this_block;
    if (slices_so_far<dz) step_line(2);
  }

  end_table(ofp);

  free(buf);
}

static void emit_count_sum_table( FILE* ofp, MissingDelta* missing_deltas )
{
  char* buf;
  int buf_size;
  int z;
  int cond;

  /* Output buffer can get pretty long for big HTML tables */
  buf_size= 512;
  if (!(buf=(char*)malloc(buf_size*sizeof(char))))
    Abort("%s: unable to allocate %d bytes!\n",progname,buf_size);

  start_table(ofp, 2, "Image Counts by Condition, Summed Over Slices:");
  step_line(1);

  switch (outMode) {
  case OUT_PS:
    {
      sprintf(buf, "                      ");
      sprintf(buf,"%s -sum- ",buf);
    }
    break;
  case OUT_HTML:
    {
      sprintf(buf,"<td></td><td></td>");
      sprintf(buf,"%s <td><b>sum</b></td> ",buf);
    }
    break;
  }
  out_line(ofp, 0, "%s", "<tr>%s</tr>", buf);
  step_line(1);
  
  for (cond=0; cond<nconditions; cond++) {
    if (page_full(0)) {
      start_page(ofp);
      step_line(1);
    }
    switch (outMode) {
    case OUT_PS:
      {
	int sum= 0;
	for (z=0; z<dz; z++) sum += counts[(cond*dz)+z];
	sprintf(buf, "%2d %-18s", cond, cond_table[cond]->name);
	sprintf(buf, "%s  %4d", buf, sum);
      }
      break;
    case OUT_HTML:
      {
	int sum= 0;
	for (z=0; z<dz; z++) sum += counts[(cond*dz)+z];
	sprintf(buf, "<td><b>%2d</b></td> <td><b>%s</b></td> ", 
		cond, cond_table[cond]->name);
	sprintf(buf,"%s <td>%4d</td>",
		buf, sum);
      }
      break;
    }
    out_line( ofp, 0, "%s", "<tr>%s</tr>", buf );
  }
  if (missing_deltas) { /* do step-by-step missing info */
    MissingDelta* tmp= missing_deltas;
    while (tmp != NULL) {
      if (page_full(2)) {
	start_page(ofp);
	step_line(2);
      }
      switch (outMode) {
      case OUT_PS:
	{
	  int sum= 0;
	  for (z=0; z<dz; z++) sum += tmp->count[z];
	  sprintf(buf,"   drop by %-10s ",trimString(tmp->step,10));
	  sprintf(buf, "%s %4d ", buf, sum);
	}
	break;
      case OUT_HTML:
	{
	  int sum= 0;
	  for (z=0; z<dz; z++) sum += tmp->count[z];
	  sprintf(buf,"<td></td><td><b>drop by %-10s</b></td>",tmp->step);
	  sprintf(buf, "%s <td>%4d</td>", buf, sum);
	}
	break;
      }
      out_line(ofp, 0, "%s", "<tr>%s</tr>", buf);
      tmp= tmp->next;
    }
  }
  else { /* only one missing entry */
    if (page_full(2)) {
      start_page(ofp);
      step_line(2);
    }
    switch (outMode) {
    case OUT_PS:
      {
	int sum= 0;
	for (z=0; z<dz; z++) sum += missing_totals[z];
	sprintf(buf,"   dropped by Fiasco  ");
	sprintf(buf, "%s %4d ", buf, sum);
      }
      break;
    case OUT_HTML:
      {
	int sum= 0;
	for (z=0; z<dz; z++) sum += missing_totals[z];
	sprintf(buf,"<td></td><td><b>dropped by Fiasco</b></td> ");
	sprintf(buf, "%s <tr>%4d</tr>", buf, sum);
      }
      break;
    }
    out_line(ofp, 0, "%s", "<tr>%s</tr>", buf);
  }
  switch (outMode) {
  case OUT_PS:
    {
      int sum= 0;
      for (z=0; z<dz; z++) sum += totals[z];
      sprintf(buf,"   total              ");
      sprintf(buf, "%s %4d ", buf, sum);
    }
    break;
  case OUT_HTML:
    {
      int sum= 0;
      for (z=0; z<dz; z++) sum += totals[z];
      sprintf(buf,"<td></td><td><b>total</b></td>");
      sprintf(buf, "%s <td>%4d</td>", buf, sum);
    }
    break;
  }
  out_line(ofp, 0, "%s", "<tr>%s</tr>", buf);

  end_table(ofp);

  free(buf);
}

static MissingDelta* create_MissingDelta()
{
  MissingDelta* result;
  int i;

  if (!(result= (MissingDelta*)malloc(sizeof(MissingDelta))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname, sizeof(MissingDelta));
  result->step= NULL;
  result->next= NULL;
  if (!(result->count= (int*)malloc(dz*sizeof(int))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname, sizeof(MissingDelta));
  for (i=0; i<dz; i++) result->count[i]= 0;
  return result;
}

static void destroy_MissingDelta( MissingDelta* d ) 
{
  if (!d) return;
  if (d->next != NULL) {
    destroy_MissingDelta( d->next );
    d->next= NULL;
  }
  if (d->step) {
    free(d->step);
    d->step= NULL;
  }
  if (d->count) {
    free(d->count);
    d->count= NULL;
  }
  free(d);
}

static int is_an_int( const char* str )
{
  for (; *str; str++) if (!isdigit(*str)) return 0;
  return 1;
}

static MissingDelta* parse_missing_trail( char* fname )
{
  MissingDelta* result= NULL;
  MissingDelta* prev= NULL;
  MissingDelta* new= NULL;
  FILE* fp= NULL;
  char buf[512];
  int i;
  int icheck, n;
  int keep_flag= 0;
  char* tstring;

  fp= fopen(fname,"r");
  if (!fp)
    Abort("%s: missing trail input file <%s> not found!\n",
	  progname, fname);

  while (!feof(fp)) {
    if (!fgets(buf,sizeof(buf),fp)) break;
    if (buf[0]=='#') continue; /* ignore comments */
    if (buf[strlen(buf)-1]=='\n') buf[strlen(buf)-1]= '\0';
    new= create_MissingDelta();
    /* We drop leading "epi." from names like "epi.deghost" to save space */
    tstring= strchr(buf,'.');
    if (tstring && !is_an_int(tstring+1)) new->step= strdup(tstring+1);
    else new->step= strdup(buf);
    for (i=0; i<dz; i++) {
      do { 
	if (!fgets(buf,sizeof(buf),fp)) break; 
	if (buf[strlen(buf)-1]=='\n') buf[strlen(buf)-1]= '\0';
      } while (buf[0]=='#'); /* ignore comments */
      if (sscanf(buf,"%d %d\n",&icheck,&n) != 2) 
	Abort("%s: last step of missing trail incomplete!\n",progname);
      if (icheck != i)
	Abort("%s: bad entry in missing file for step %s!\n",
	      progname,new->step);
      new->count[i]= n;
    }
    if (!prev) {
      prev= new;
      keep_flag= 0;
      for (i=0; i<dz; i++) if (new->count[i]) { keep_flag= 1; break; }
      if (keep_flag) {
	MissingDelta* tmp= create_MissingDelta();
	tmp->step= strdup(new->step);
	for (i=0; i<dz; i++) tmp->count[i]= new->count[i];
	tmp->next= result;
	result= tmp;
      }
    }
    else {
      keep_flag= 0;
      for (i=0; i<dz; i++) {
	if (new->count[i] != prev->count[i]) { keep_flag= 1; break; }
      }
      if (keep_flag) {
	MissingDelta* tmp= create_MissingDelta();
	tmp->step= strdup(new->step);
	for (i=0; i<dz; i++) 
	  tmp->count[i]= new->count[i] - prev->count[i];
	tmp->next= result;
	result= tmp;
	destroy_MissingDelta(prev);
	prev= new;
	new=NULL;
      }
      else {
	destroy_MissingDelta(new);
	new= NULL;
      }
    }
    if (ferror(fp)) break;
  }
  if (ferror(fp)) {
    perror(progname);
    Abort("%s: error reading missing file %s!\n",progname,fname);
  }
    
  /* Check to make sure totals match expectation */
  if (prev && missing_totals) {
    for (i=0; i<dz; i++) 
      if (prev->count[i] != missing_totals[i]) {
	Abort("%s: missing file not compatible with input MRI file!\n",
	      progname);
      }
  }

  /* Do some cleanup */
  fclose(fp);
  if (prev) destroy_MissingDelta(prev);
  if (new) destroy_MissingDelta(new);

  return result;
}


/* Function to calculate summary statistics on a set of parameters */
static void stats( numpars, offset )
     long numpars, offset;
{
  long p, z, t, count;

  /* Loop through number of parameters over which to calculate */
  /*   summaries and then over slices                          */
  for( p = 0; p < numpars; p++ )
    for( z = 0; z < dz; z++ )
      {
	/* Initialize values */
	mean[z][offset+p] = stdv[z][offset+p] = 0.0;
	t = 0;
	while( t<dt && missing[t][z] )
	  t++;
	if( t == dt )
	  {
	    Warning( 1, "All images for slice %ld are missing.", z );
	    continue;
	  }
	min[z][offset+p] = max[z][offset+p] = pars[p][t][z];
	

	/* Loop through images, summing up and checking for new extremes */
	count = 0;
	for( t = 0; t < dt; t++ )
	  {
	    if( !missing[t][z] )
	      {
		count++;
		mean[z][offset+p] += pars[p][t][z];
		stdv[z][offset+p] += ( pars[p][t][z] * pars[p][t][z] );
		min[z][offset+p] = ( pars[p][t][z] < min[z][offset+p] )?
		  pars[p][t][z]: min[z][offset+p];
		max[z][offset+p] = ( pars[p][t][z] > max[z][offset+p] )?
		  pars[p][t][z]: max[z][offset+p];
	      }
	  }

	/* Convert sums to stats */
	if( count )
	  {
	    mean[z][offset+p] /= (double) count;
	  }
	if( count > 1 )
	  {  
	    stdv[z][offset+p] = ( stdv[z][offset+p] - (double) count * 
	                          pow( mean[z][offset+p], 2.0 ) ) /
	      (double) ( count - 1 );
	    stdv[z][offset+p] = ( stdv[z][offset+p] > 0 )?
	      sqrt( stdv[z][offset+p] ): 0.0;
	  }
	else
	  stdv[z][offset+p] = 0;

      }
  return;
}

static void write_env( FILE* ofp, char* envp[] ) 
{
  int i;
  char buf[512];
  char* word2;
  
  start_page(ofp);
  
  step_line(1);
  start_table(ofp,2,"Fiasco Environment");
  step_line(1);

  for (i=0; envp[i]; i++) {
    if (!strncmp(envp[i],"F_",2)) {
      strncpy(buf,envp[i],sizeof(buf));
      buf[sizeof(buf)-1]= '\0';
      if (word2=strchr(buf,'=')) {
	*word2= '\0';
	word2++;
      }
      else word2="";
      out_line( ofp, 0, "%s=<%s>", "<tr><td>%s</td> <td>%s</td></tr>", 
		buf, word2);
      if (page_full(1)) {
	start_page(ofp);
	step_line(1);
	start_table(ofp,2,"Fiasco Environment (cont)");
	step_line(1);
      }
    }
  }

  end_table(ofp);
}

static void scanParList(char* listfile)
{
  FILE* lfp= NULL;
  char parfile[512], scanline[512], whole_parname[512];
  char* parname;
  long m;
  FormatInfo *thisFI= NULL;


  /* Open parameter list file */
  lfp = efopen( listfile, "r" );

  /* Read parameter list to find and classify the parameter subfiles.
   * We make the assumption that anything with a name like 
   * "epi.baseline" can be classified according to the second
   * part of the name, that is, as "baseline".  Note that this
   * gets screwed up if the last character of parname is '.'
   */
  formatsHead= formatsTail= NULL;
  while( !feof( lfp ) )
    {
      fscanf( lfp, "%510[^\n]%*[\n]", scanline );
      m = sscanf( scanline, "%s %s", whole_parname, parfile );
      if( m == 1 )
	{
	  Warning( 1, "%s missing filename.\n", parname );
	  continue;
	}
      if( m !=2 )
	continue;
      if ((parname= strchr(whole_parname,'.'))==NULL)
	parname= whole_parname;
      else parname++; /* move off the . */

      thisFI= createFormatInfo(parname,parfile,dz,dt);

      /* See what we know about this parameter name */
      if (!strcmp( whole_parname, "ts.deghost" )
	  || !strcmp( whole_parname, "ts.deghost_raw" )) {
	thisFI->order= PAR_ORDER_ZST;
	thisFI->nShots= 2;
	addFieldName(thisFI,"shot 1");
	addFieldName(thisFI,"shot 2");
      }
      else if( !strcmp( parname, "meanc" ) || !strcmp( parname, "phadj" ) ||
	       !strcmp( parname, "outlier" ) ||
	       !strcmp( parname, "displace" ) ||
	       !strcmp( parname, "smdisplace" ) ||
	       !strcmp( parname, "deghost_raw" ) || 
	       !strcmp( parname, "deghost" ) ||
	       !strcmp( parname,"deghost2") ) {
	thisFI->order= PAR_ORDER_Z_FASTEST;
	addFieldName(thisFI,"");
      }
      else if( !strcmp( parname, "baseline" ) ) {
	  thisFI->order= PAR_ORDER_Z_FASTEST;
	  addFieldName(thisFI,"Real");
	  addFieldName(thisFI,"Imaginary");
	}
      else if( !strcmp( parname, "baseline_raw" ) ) {
	  thisFI->order= PAR_ORDER_Z_FASTEST;
	  addFieldName(thisFI,"Real");
	  addFieldName(thisFI,"Imaginary");
	}
      else if( !strcmp( parname, "baseline2" ) ) {
	thisFI->order= PAR_ORDER_Z_FASTEST;
	addFieldName(thisFI,"LeftReal");
	addFieldName(thisFI,"LeftImag");
	addFieldName(thisFI,"RightReal");
	addFieldName(thisFI,"RightImag");
      }
      else if( !strcmp( parname, "estireg" ) ||
	       !strcmp( parname, "parsm" ) ) {
	thisFI->order= PAR_ORDER_INDEX_TZ;
	addFieldName(thisFI,"X-shift");
	addFieldName(thisFI,"Y-shift");
	addFieldName(thisFI,"Rotation");
	addFieldName(thisFI,"MSE");
      }
      else if( !strcmp( parname, "estireg3d" ) ||
	       !strcmp( parname, "parsm3d" ) ) {
	thisFI->order= PAR_ORDER_INDEX_T;
	addFieldName(thisFI,"quat X");
	addFieldName(thisFI,"quat Y");
	addFieldName(thisFI,"quat Z");
	addFieldName(thisFI,"quat W");
	addFieldName(thisFI,"X-shift");
	addFieldName(thisFI,"Y-shift");
	addFieldName(thisFI,"Z-shift");
	addFieldName(thisFI,"MSE");
      }
      else if( !strcmp( parname, "displace3d_raw" ) ||
	       !strcmp( parname, "displace3d" ) ) {
	thisFI->order= PAR_ORDER_INDEX_T;
	addFieldName(thisFI,"3D Rot");
	addFieldName(thisFI,"3D Transl");
	addFieldName(thisFI,"3D Disp");
      }

      /* Let any header info in the file override known format info */
      updateFormat(thisFI);

      if (formatValid(thisFI)) {
	if (formatsHead==NULL) formatsHead= thisFI;
	if (formatsTail!=NULL) formatsTail->next= thisFI;
	formatsTail= thisFI;
      }
      else {
	Error("Warning: format for file %s is unknown!\n",getFileName(thisFI));
	destroyFormatInfo(thisFI);
      }
      
    }

  fclose(lfp);
}

static void emit_stats(FILE* ofp)
{
  int lines_per_block;
  long z, parnum;
  FormatInfo *thisFI= NULL;
  
  start_page(ofp);

  /* Loop through slices */
  lines_per_block= 0; /* count 'em first time through */
  for( z = 0; z < dz; z++ ) {
    long start_loc= psLocation;
    
    /* Set up new page as needed */
    if (page_full(lines_per_block)) {
      start_page(ofp);
    }
    
    /* Initialize for slice printing */
    step_line(2);
    start_table(ofp, 2, "Slice %ld", z);
    step_line(1);
    out_line(ofp, 0, 
	     "                                 Mean    St.Dev.   Min      Max",
	     "<tr><td></td><td></td><td><b>Mean</b></td><td><b>St.Dev.</b></td><td><b>Min</b></td><td><b>Max</b></td></tr>");
    step_line(1);
    
    /* Loop through parameter files and print summary statistics */
    parnum = 0;
    thisFI= formatsHead;
    while (thisFI) {
      int i;
      for (i=0; i<thisFI->nFields; i++) {
	out_line(ofp,0,
		 "%-14s    %-12s %8.2g %8.2g %8.2g %8.2g",
		 "<tr><td><b>%s</b></td><td><b>%s</b></td><td>%8.2g</td><td>%8.2g</td><td>%8.2g</td><td>%8.2g</td></tr>",
		 getTrimmedParName(thisFI,14), getFieldName(thisFI,i),
		 mean[z][parnum], stdv[z][parnum],
		 min[z][parnum], max[z][parnum], psLocation );
	parnum++;
      }
      thisFI= thisFI->next;
    }
    
    end_table(ofp);
    
    if (!lines_per_block) {
      /* first time through, we count lines for pagination purposes */
      lines_per_block= (start_loc - psLocation)/psIncrement;
    }
  }
}

static void calcStats(long fixed_image, unsigned char* tmpmiss)
{
  FormatInfo *thisFI= NULL;
  long z, parnum;

  /* CALCULATE STATS */
  /* Go back to beginning of parameter list and loop through the list.
   * We make the assumption that anything with a name like 
   * "epi.baseline" can be classified according to the second
   * part of the name, that is, as "baseline".  Note that this
   * gets screwed up if the last character of parname is '.'
   */
  parnum = 0;
  thisFI= formatsHead;
  while (thisFI) {
    int skipFixed= 0;

    load(thisFI);

    if ( !strcmp( getParName(thisFI), "estireg" )
	 || !strcmp( getParName(thisFI), "parsm" )
	 || !strcmp( getParName(thisFI), "estireg3d" )
	 ||  !strcmp( getParName(thisFI), "parsm3d" )
	 ||  !strcmp( getParName(thisFI), "displace3d_raw" )
	 || !strcmp( getParName(thisFI), "displace3d" ) ) 
      skipFixed= 1;
    else skipFixed= 0;

    if (skipFixed) {
      if( ( fixed_image >= 0 ) && ( fixed_image < dt ) )
	for( z = 0; z < dz; z++ )
	  {
	    tmpmiss[z] = missing[fixed_image][z];
	    missing[fixed_image][z] = (unsigned char) 1;
	  }
      stats( thisFI->nFields, parnum );
      if( ( fixed_image >= 0 ) && ( fixed_image < dt ) )
	for( z = 0; z < dz; z++ )
	  {
	    missing[fixed_image][z] = tmpmiss[z];
	  }
    }
    else {
      stats( thisFI->nFields, parnum );
    }
    parnum += thisFI->nFields;
    thisFI= thisFI->next;
  }

  /* FINISHED CALCULATING STATS */
}

int main( int argc, char**argv, char* envp[] ) 
{
  MRI_Dataset *Input = NULL;
  FILE *ifp = NULL, *ofp = NULL;
  char infile[512], listfile[512], missfile[512], outfile[512];
  char splitfile[512], condfile[512], fixed_image_string[512];
  long fixed_image, t, z, m;
  long num_pars, parnum, *numnonmiss = NULL;
  unsigned char *tmpmiss = NULL;
  int missing_trail_flag= 0;
  int env_flag= 0;
  MissingDelta* missing_deltas= NULL;
  FormatInfo *thisFI= NULL;

  progname= argv[0];

  /* Print version number */
  Message( "# %s\n", rcsid );

  /* Check to see if help was requested */
  if( ( argc > 1 ) && !strcmp( argv[1], "-help" ) )
    {
      if( argc == 2 )
	Help( "selecttopic" );
      else
	Help( argv[2] );
    }

  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Get filenames */
  cl_get( "headerinput|h", "%option %s[%]", "input.mri", infile );
  cl_get( "list|l", "%option %s[%]", "parlist", listfile );
  cl_get( "out|o", "%option %s[%]", "summary.ps", outfile );
  cl_get( "fixed|f", "%option %s[%]", "-1", fixed_image_string );
  cl_get( "split|f", "%option %s[%]", "newsplit", splitfile );
  cl_get( "cond|c", "%option %s[%]", "condition", condfile );
  env_flag= cl_present("env");
  missing_trail_flag= cl_get( "missing|m", "%option %s", missfile );
  if (cl_present("html")) {
    outMode= OUT_HTML;
    if (cl_present("ps")) {
      fprintf(stderr,"%s: the -ps and -html flags are mutually exclusive.\n",
	      argv[0]);
      exit(-1);
    }
  }
  else outMode= OUT_PS;

  if (cl_cleanup_check()) {
    int i;
    fprintf(stderr,"%s: invalid argument in command line:\n    ",argv[0]);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(-1);
  }
  
  /*** End command-line parsing ***/

  /* Open input dataset */
  Input = mri_open_dataset( infile, MRI_READ );

  /* Check that program will function on data-set */
  if( !mri_has( Input, "images" ) )
    Abort( "%s must have images key.", infile );
  if( !mri_has( Input, "images.extent.t" ) ||
      !mri_has( Input, "images.extent.z" ) )
    Abort( "%s must have images.extent keys for t- and z- dimensions.",
	   infile );

  /* Set parameters in local variables */
  dt = mri_get_int( Input, "images.extent.t" );
  dz = mri_get_int( Input, "images.extent.z" );
  if( dz <= 0 )
    Abort( "images.extent.z is non-positive." );
  if( dt <= 1 )
    Abort( "images.extent.t is less than 2." );

  /* Find a fixed image based on the given string.  It could be an
   * integer, or a word like "middle", "center", or "mean", or possibly
   * a filename.
   */
  if (strcasecmp(fixed_image_string, "middle" ) == 0 ||
      strcasecmp(fixed_image_string, "center" ) == 0)
    fixed_image = dt / 2;
  else if (!strcasecmp(fixed_image_string, "mean")) {
    fixed_image= -1; /* use no fixed image */
  }
  else {
    char* runner= fixed_image_string;
    if (!isdigit(*runner) && !(index("+-",*runner))) {
      /* This is presumably a filename; ignore it */
      fixed_image= -1; /* use no fixed image */
    }
    else {
      runner++;
      for ( ; *runner; runner++ )
	if (!isdigit(*runner) && !index("+-",*runner)) {
	  /* First character is a digit but it's not an integer */
	  Abort("%s: unrecognized fixed image <%s>!\n",
		argv[0],fixed_image_string);
	}
      fixed_image= atoi(fixed_image_string);
      if (fixed_image >= dt) fixed_image= dt-1;
    }
  }

  /* Read/Create missing image indicators */
  missing = get_missing( Input );

  /* Done with header file */
  mri_close_dataset( Input );

  /* Generate the table of slice counts by condition */
  if (!access(condfile,R_OK) && !access(splitfile,R_OK))
    do_slice_counts(condfile, splitfile, dz, dt);
  else {
    Message("%s: condition or split file not given or unreadable!\n", 
	    argv[0]);
    fake_slice_counts(dz, dt);
  }
      
  /* Parse the trail of missing counts by step if it is provided */
  if (missing_trail_flag) 
    missing_deltas= parse_missing_trail(missfile);
  else missing_deltas= NULL;

  /* Scan the parameter list file, building a linked list of params 
   * to be dealt with.
   */
  scanParList(listfile);

  /* Count pars */
  thisFI= formatsHead;
  num_pars= 0;
  while (thisFI) {
    num_pars += thisFI->nFields;
    thisFI= thisFI->next;
  }

  /* Allocate parameter and summary storage */
  pars = (float ***) emalloc( MAX_PARS_PER_FILE * sizeof(float**) );
  for( m = 0; m < MAX_PARS_PER_FILE; m++ )
    pars[m] = Matrix( dt, dz, float );
  mean = Matrix( dz, num_pars, double );
  stdv = Matrix( dz, num_pars, double );
  min = Matrix( dz, num_pars, double );
  max = Matrix( dz, num_pars, double );
  tmpmiss = (unsigned char *) emalloc( dz * sizeof(unsigned char) );
  numnonmiss = (long *) emalloc( dz * sizeof(long) );

  /* Calculate number of non-missing images per slice */
  for( z = 0; z < dz; z++ )
    {
      numnonmiss[z] = dt;
      for( t = 0; t < dt; t++ )
	if( missing[t][z] )
	  numnonmiss[z]--;
    }

  /* OK, we have all the info necessary to calculate the statistics */
  calcStats(fixed_image, tmpmiss);

  /* Open output file and set up */
  ofp = efopen( outfile, "w" );
  if (!ofp) Abort("%s: unable to open %s for writing!\n",argv[0],outfile);

  /* Write out all the stuff we've calculated */
  init_output(ofp);
  start_page(ofp);
  step_line(2); /* skip 2 lines */
  out_line(ofp, 2,
	   "Number of images is %ld; Number of slices is %ld",
	   "<center><font size=+1>Number of images is %ld; Number of slices is %ld</font></center><p>",
	   dt,dz);
  step_line(1);
  emit_count_table( ofp, missing_deltas );
  emit_count_sum_table( ofp, missing_deltas );
  emit_contents(ofp, env_flag);
  emit_stats(ofp);
  if (env_flag) {
    write_env( ofp, envp );
  }

  /* Close up and quit */
  end_output(ofp);
  fclose( ofp );

  Message( "#      Summary statistics calculation complete.\n" );

  return 0;
}

