/*
 *	Miscellaneous utility functions
 *
 *	Copyright (c) 1996  Pittsburgh Supercomputing Center
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
 *
 *	HISTORY
 *		12/95 Written by Greg Hood (PSC)
 */
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>
#include <string.h>
#ifdef PVM
#include "pvm3.h"
#endif
#include "errors.h"
#include "misc.h"

#define reportFile	miscReportFile

static char rcsid[] = "$Id: libmisc.c,v 1.10 2005/03/04 00:39:26 welling Exp $";

Boolean verbose = FALSE;
FILE *reportFile = NULL;

/* Report is used for normal messages to the user which are only displayed
   if the program is in verbose mode */
void
Misc_Report (char *fmt, ...)
{
  va_list args;

  if (!verbose)
    return;

  va_start(args, fmt);
  vprintf(fmt, args);
  va_end(args);
  fflush(stdout);
}

/* Message is used for normal messages to the user. */
void
Misc_Message (char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);
  vprintf(fmt, args);
  va_end(args);
  fflush(stdout);
}

/* Error is used for recoverable errors which the user should be notified
   about */
void
Misc_Error (char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
}

/* Warnings are here for compatibility; they are like errors. */
void
Misc_Warning (int severity, char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
}


/* Abort is used for unrecoverable errors which the user should be notified
   about */
void
Misc_Abort (char *fmt, ...)
{
  va_list args;
  FILE *f;
  Filename fn;
  char hn[128];

  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);

  if (getenv("F_DEBUG")) {
    /* we also print the error into a file in the home directory */
    char* cwd= getcwd(NULL,0);
    gethostname(hn, 128);
    sprintf(fn, "~/abort.%s.%d.log", hn, getpid());
    f = fopen(fn, "w");
    fprintf(f,"Error while running on %s in %s\n",hn,cwd);
    va_start(args, fmt);
    vfprintf(f, fmt, args);
    va_end(args);
    fclose(f);
    free(cwd);
  }

#ifdef T3D
  /* on the T3D we force all PEs to exit, not just the calling PE */
  globalexit(1);
#else
  abort();
#endif
}

void Misc_SysError ( char* file, int line, char *fmt, ... )
{
  va_list args;
  FILE *f;
  Filename fn;
  char hn[128];

  va_start(args, fmt);
  if (line>0) 
    fprintf(stderr, "Function %s at line %d: ", file, line);
  else
    fprintf(stderr, "Function %s: ", file);
  vfprintf(stderr, fmt, args);
  va_end(args);

  if (getenv("F_DEBUG")) {
    /* we also print the error into a file in the home directory */
    char* cwd= getcwd(NULL,0);
    gethostname(hn, 128);
    sprintf(fn, "/%s/abort.%s.%d.log", getenv("HOME"), hn, getpid());
    f = fopen(fn, "w");
    fprintf(f,"Error while running on %s in %s\n",hn,cwd);
    if (line>0)
      fprintf(f, "Function %s at line %d: ", file, line);
    else
      fprintf(f, "Function %s: ", file);
    va_start(args, fmt);
    vfprintf(f, fmt, args);
    va_end(args);
    fclose(f);
    free(cwd);
  }

#ifdef T3D
  /* on the T3D we force all PEs to exit, not just the calling PE */
  globalexit(1);
#else
  abort();
#endif
}

Boolean
Misc_StringCopy (char *to,
	    const char *from,
	    const int maxChars)
{
  if (strlen(from) < maxChars)
    {
      strcpy(to, from);
      return(TRUE);
    }
  else
    {
      strncpy(to, from, maxChars-1);
      to[maxChars-1] = '\0';
      return(FALSE);
    }
}


/* GetEnvInt returns the integer value of an environment variable */
int
Misc_GetEnvInt (const char *var_name)
{
  char *val;
  int i;

  if ((val = getenv(var_name)) == NULL)
    return(0);
  if (sscanf(val, "%d", &i) != 1)
    Abort("Environment variable %s must be set to an integral value.\n",
	  var_name);
  return(i);
}

/* GetEnvFloat return the float value of an environment variable */
float
Misc_GetEnvFloat (const char *var_name)
{
  char *val;
  float f;

  if ((val = getenv(var_name)) == NULL)
    return(0);
  if (sscanf(val, "%f", &f) != 1)
    Abort("Environment variable %s must be set to a numerical value.\n",
	  var_name);
  return(f);
}

/* ConstructFilename returns a pathname for a file given a
   directory and a possibly unprefixed filename */
void
Misc_ConstructFilename (Filename out,
		   const char *dir,
		   const char *in)
{
  int len;

  if (in[0] == '/' || in[0] == '.')
    /* the filename is already prefixed so just copy it */
    strcpy(out, in);
  else
    {
      len = strlen(dir);
      if (len+1+strlen(in) >= sizeof(Filename))
	Abort("Filename too long.\n");
      strcpy(out, dir);
      out[len] = '/';
      strcpy(&out[len+1], in);
    }
}

int
Misc_Round (float f)
{
  int i;

  /* this is to make the behavior of Round conform
     to that of rint with the default rounding mode */
  i = floor(f + (1.0/2.0));
  if ((((float) i) - f) == (1.0/2.0))
    i &= ~1;
  return(i);
}
