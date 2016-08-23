/*
 *	printmrifield.c - Example MRI library program
 *
 *	This program prints out one field from an MRI dataset.  It is
 *	useful in shell scripts to query the value of particular fields.
 *	If a field does not exist it does not print anything but returns
 *	with an exit value of 1; otherwise 0.
 *
 *	Copyright (c) 1997 Pittsburgh Supercomputing Center
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
 *	HISTORY:
 *		2/97 - Written by Greg Hood (PSC)
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"

static char rcsid[] = "$Id: mri_printfield.c,v 1.4 2003/08/07 20:43:16 bakalj Exp $";

static int nofail_flag= 0;
static int wildcard_flag= 0;
static int verbose_flag= 0;

static void emitKeyValue( MRI_Dataset* ds, char* key )
{
  char* s = mri_get_string(ds, key);
  if (s == NULL) {
    if (nofail_flag) printf("\n");
    else exit(1);
  }
  else {
    if (verbose_flag) printf("<%s> = <%s>\n", key, s);
    else printf("%s\n", s);
  }
}

/* We want matches to work in both directions, so that a string
 * must match both forwards and backwards.
 */
static int patternMatchBackwards( const char* s, const char* p )
{
  const char* rs= s + strlen(s) - 1;
  const char* rp= p + strlen(p) - 1;
  while (rs > s && rp > p) {
    if (*rs == *rp || *rp == '?') {
      rs--;
      rp--;
    }
    else if (*rp == '*') {
      s--;
    }
    else return 0;
  }
  if ((*rs == *rp) || (*rp=='*')) return 1;
  else return 0;
}

static int patternMatch(const char* s_orig, const char* p_orig)
{
  const char* s= s_orig;
  const char* p= p_orig;
  while (*s && *p) {
    if (*s == *p || *p == '?') {
      s++;
      p++;
    }
    else if (*p == '*') {
      s++;
    }
    else return 0;
  }
  if (!(*s)  && (!(*p) || (*p=='*'))) 
    return patternMatchBackwards(s_orig,p_orig);
  else return 0;
}

int main (argc, argv)
     int argc;
     char **argv;
{
  MRI_Dataset *ds;
  char infile[512];
  char field[512];

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /*** Parse command line ***/
  cl_scan( argc, argv );

  /* Deprecate old options */

  if (cl_present( "input|i" ))
     Abort ("Option input|i has been replaced by infile outfile format.  Please see help file.\n");
  if (cl_present( "n" ))
     Abort ("Option n has been expanded to nofail|n.  Please see help file.\n");
  if (cl_present( "f" ))
     Abort ("Option f has been expanded to field|fld .  Please see help file.\n");
  if (cl_present( "w" ))
     Abort ("Option w has been expanded to wildcard|wld.  Please see help file.\n");

  /* Get input params */
  cl_get( "field|fld", "%option %s[%]", "images.dimensions", field );
  if (cl_present("nofail|nof")) nofail_flag= 1;
  if (cl_present("verbose|ver|v")) verbose_flag= 1;
  if (cl_present("wildcard|wld")) wildcard_flag= 1;

  if(!cl_get("", "%s", infile)) {
    fprintf(stderr, "%s: Input file name not given.\n", argv[0]);
    exit(-1);
  }

  if (cl_cleanup_check()) {
    int i;
    fprintf(stderr,"%s: invalid argument in command line:\n    ",argv[0]);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(-1);
  }
  /*** End command-line parsing ***/

  mri_set_error_handling(MRI_IGNORE_ERRORS);
  ds = mri_open_dataset(infile, MRI_READ);
  if (ds == NULL)
    {
      fprintf(stderr, "printmrifield: cannot open dataset %s\n", infile);
      exit(1);
    }
  if (wildcard_flag) {
    char* key;
    mri_iterate_over_keys(ds);
    while ((key= mri_next_key(ds)) != NULL) {
      if (patternMatch(key,field)) emitKeyValue(ds,key);
    }
  }
  else {
    emitKeyValue(ds,field);
  }
  mri_close_dataset(ds);
  return(0);
}
