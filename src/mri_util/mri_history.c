/************************************************************
 *                                                          *
 *  mri_history.c                                           *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1997 Department of Statistics,         *
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
 *  Original programming by Joel Welling, 5/98              *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"

#define KEYBUF_SIZE 512

static char rcsid[] = "$Id: mri_history.c,v 1.7 2005/06/28 23:37:18 welling Exp $";

int main( int argc, char* argv[] ) 
{
  char fname[512];
  char comment[512];
  int add_flag;
  int del_flag;
  int pdel_flag;
  int num_to_del;
  MRI_Dataset* ds;

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Deprecate old options */

  if (cl_present( "a" ))
     Abort ("Option a(add) has been expanded to add.  Please see help file.\n");

  add_flag= cl_get("add", "%option %s",comment);
  pdel_flag= cl_get("partialdelete|pdl","%option %d",&num_to_del);
  del_flag= cl_present("delete|del");
  if (!cl_get("", "%s", fname)) {
    fprintf(stderr,"%s: File name not given.\n",argv[0]);
    Help( "usage" );
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

  /* Open input dataset */
  ds = mri_open_dataset( fname, 
			 ((add_flag || del_flag || pdel_flag) 
			  ? MRI_MODIFY : MRI_READ) );

  if (del_flag || add_flag || pdel_flag) {
    if (del_flag) {
      hist_delete_all( ds );
      if (!add_flag) {
	/* Leave some sort of mark about what happened */
	hist_add_cl( ds, argc, argv );
      }
    }
    else if (pdel_flag) {
      hist_delete_some( ds, num_to_del );
      if (!add_flag) {
	/* Leave some sort of mark about what happened */
	hist_add_cl( ds, argc, argv );
      }
    }
    else if (add_flag) {
      /* we must add the timestamp, since hist_add() does not do so. */
      char strbuf[256];
      time_t t;
      struct tm* tsp;
      long tot_length;
      time(&t);
      tsp= localtime(&t);
      tot_length= strftime(strbuf,sizeof(strbuf),"added [%F %T] ",tsp);
      strncat(strbuf,comment,sizeof(strbuf)-(tot_length+1));
      strbuf[sizeof(strbuf)-1]= '\0'; /* just for paranoia */
      hist_add(ds, strbuf);
    }
  }
  else {
    if (hist_get(ds,1)) hist_dump(ds, stdout);
    else fprintf(stdout,"No history info present in %s.\n",fname);
  }
  mri_close_dataset(ds);

  return 0;
}

