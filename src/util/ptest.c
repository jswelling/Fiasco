/************************************************************
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *                                                          *
 *     Copyright (c) 1999 Carnegie Mellon University        *
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
 ***********************************************************/
/*
 *	ptest.c - test libpar.c
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include "par.h"
#include "bio.h"
#include "array.h"
#include "misc.h"
#include "acct.h"

typedef struct Context {
  /* NOTE: any new fields added to this struct should
     be also added to PackContext and UnpackContext */ 
  char context_name[80];
} Context;

typedef struct Task {
  /* NOTE: any new fields added to this struct should
     be also added to PackTask and UnpackTask */ 
  char task_name[80];
} Task;

static char rcsid[] = "$Id: ptest.c,v 1.5 2000/10/06 00:32:22 welling Exp $";

/* GLOBAL VARIABLE FOR MASTER & SLAVE */
Task t;
Context c;

/* GLOBAL VARIABLES FOR MASTER */
FILE *f;

/* GLOBAL VARIABLES FOR SLAVE */
FILE *sf;

/* FORWARD DECLARATIONS */
void MasterTask (int argc, const char **argv);
void ReadEnvironment ();
int ReadArguments ();
void SlaveContext ();
void SlaveTask ();

#ifdef PVM
void PackContext();
void UnpackContext();
void PackTask();
void UnpackTask();
#endif


int
main (int argc,
      char **argv,
      char **envp)
{
#ifdef PVM
  par_process(argc, argv, envp,
	      MasterTask, NULL,
	      SlaveContext, SlaveTask,
	      PackContext, UnpackContext,
	      PackTask, UnpackTask,
	      NULL, NULL);
#else
  par_process(argc, argv, envp,
	      MasterTask, NULL,
	      SlaveContext, SlaveTask,
	      NULL, NULL,
	      NULL, NULL,
	      NULL, NULL);
#endif
  exit(0);
}

/* MASTER PROCEDURES */

void
MasterTask (int argc,
	    const char **argv)
{
  int i;
  int n;
  char name[80];

  sprintf(name, "master%d", getpid());
  f = fopen(name, "w");
  fprintf(f, "Master started.\n");
  fflush(f);
  Report("Master starting...\n");
  sprintf(c.context_name, "context0");
  par_set_context();
  for (n = 0; n < 20; ++n)
    {
      sprintf(t.task_name, "task%d", n);
      par_delegate_task();
    }

  par_finish();
  fprintf(f, "Master finished.\n");
  fflush(f);
  fclose(f);
  Report("Done!!\n");
}

/* SLAVE PROCEDURES */

void
SlaveContext ()
{
  char name[80];

  sprintf(name, "slave%d", getpid());
  sf = fopen(name, "w");
  fprintf(sf, "Got new context: %s\n", c.context_name);
  fflush(f);
}

void
SlaveTask ()
{
  fprintf(sf, "Got a task: %s\n", t.task_name);
  fflush(f);
  sleep(4);
  fprintf(sf, "Finished task: %s\n", t.task_name);
  fflush(f);
}

/* UTILITY FUNCTIONS */

/* PVM-RELATED FUNCTIONS */

#ifdef PVM
void
PackContext ()
{
  if (pvm_pkstr(c.context_name) < 0)
    Abort("Could not pack context.\n");
}

void
UnpackContext ()
{
  if (pvm_upkstr(c.context_name) < 0)
    Abort("Could not unpack context.\n");
}

void
PackTask ()
{
  if (pvm_pkstr(t.task_name) < 0)
    Abort("Could not pack task\n");
}

void
UnpackTask ()
{
  if (pvm_upkstr(t.task_name) < 0)
    Abort("Could not unpack task\n");
}
#endif

