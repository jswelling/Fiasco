/*
 *	mri_setfield.c
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
 *		2/97 - printmrifield Written by Greg Hood (PSC)
 *              9/02 - transformed into mri_setfield by Joel Welling (PSC)
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "mri.h"
#include "fmri.h"
#include "misc.h"
#include "stdcrg.h"

static char rcsid[] = "$Id: mri_setfield.c,v 1.5 2004/02/07 01:10:52 welling Exp $";

int main (argc, argv)
     int argc;
     char **argv;
{
  MRI_Dataset *ds;
  char infile[512];
  char field[512];
  char value[512];
  int delete_flag= 0;
  int fld_present= 0;
  int val_present= 0;
  int all012_present= 0;
  int allxyz_present= 0;
  char *s;

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /*** Parse command line ***/
  cl_scan( argc, argv );

  /* Deprecate old options */

  if (cl_present( "input|i" ))
    Abort ("Option input|i has been replaced by infile outfile format.  Please see help file.\n");
  if (cl_present( "f" ))
     Abort ("Option f has been expanded to field|fld.  Please see help file.\n");
  if (cl_present( "v" ))
     Abort ("Option v has been expanded to value|val.  Please see help file.\n");
  

  /* Get input params */
  fld_present= cl_get( "field|fld", "%option %s", field );
  val_present= cl_get( "value|val", "%option %s", value );
  if (cl_present("delete|del")) delete_flag= 1;
  all012_present= cl_present("all012");
  allxyz_present= cl_present("allxyz");

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

  if (!fld_present)
    Abort("%s: required field name not given!\n",argv[0]);
  if (!val_present && !delete_flag)
    Abort("%s: must either give a value or specify -delete!\n",argv[0]);

  ds = mri_open_dataset(infile, MRI_MODIFY);
  if (ds == NULL) 
    Abort("%s: cannot open dataset %s for modification!\n", infile);
  if (delete_flag) {
    if (all012_present) {
      char fullField[512];
      int i;
      for (i=0; i<3; i++) {
	snprintf(fullField,sizeof(fullField)-1,"%s.%d",field,i);
	fullField[sizeof(fullField)-1]= '\0';
	mri_remove(ds,fullField);
      }    
    }
    else if (allxyz_present) {
      char fullField[512];
      int i;
      for (i=0; i<3; i++) {
	snprintf(fullField,sizeof(fullField)-1,"%s.%c",field,'x'+i);
	fullField[sizeof(fullField)-1]= '\0';
	mri_remove(ds,fullField);
      }    
    }
    else {
      mri_remove(ds,field);
    }
  }
  else {
    if (all012_present || allxyz_present) {
      char fullField[512];
      char oneVal[512];
      char* tmp;
      int i;
      int start= 0;
      int end= strlen(value)-1;
      /* Strip things like quotes and spaces.  This is necessary
       * because of annoying behavior on the part of clproc.
       */
      while ((ispunct(value[start]) && value[start]!='-') 
	     || isspace(value[start])) start++;
      while (ispunct(value[end]) || isspace(value[start])) { 
	value[end]= '\0'; end--; 
      }

      for (i=0; i<3; i++) {
	char* tok= strtok_r((i==0)?value+start:NULL,":,",&tmp);
	if (!tok) Abort("%s: can't get token for value %d!\n",argv[0],i);
	if (all012_present) 
	  snprintf(fullField,sizeof(fullField)-1,"%s.%d",field,i);
	else if (allxyz_present)
	  snprintf(fullField,sizeof(fullField)-1,"%s.%c",field,'x'+i);
	fullField[sizeof(fullField)-1]= '\0';
	mri_set_string(ds,fullField,tok);
      } 
    }
    else {
      mri_set_string(ds,field,value);
    }
  }
  hist_add_cl(ds, argc, argv);
  mri_close_dataset(ds);
  return(0);
}
