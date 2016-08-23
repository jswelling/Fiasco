/*
 *	Very simple elapsed time accounting functions
 *		(to enable these, compile with -DACCT)
 *
 *	Copyright (c) 1995  Pittsburgh Supercomputing Center
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
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <unistd.h>
#include "acct.h"
#include "misc.h"

static char rcsid[] = "$Id: libacct.c,v 1.3 1999/07/07 20:14:45 welling Exp $";

int acct_bins[8] = {0, 0, 0, 0, 0, 0, 0, 0};

void
Acct (int n)
{
#ifdef ACCT
#ifdef T3D
  static int lastAcct = -1;
  static double last_time;
  double this_time;
  int elapsedMicroseconds;

  this_time = ((long) rtclock()) / ((double) CLK_TCK);
  if (lastAcct >= 0)
    {
      elapsedMicroseconds = (this_time - last_time) * 1000000.0;
      acct_bins[lastAcct] += elapsedMicroseconds;
    }
  lastAcct = n;
  last_time = this_time;
#else
  static int lastAcct = -1;
  static struct timeval last_time;
  struct timeval this_time;
  int elapsedMicroseconds;

  gettimeofday(&this_time, NULL);
  if (lastAcct >= 0)
    {
      elapsedMicroseconds = (this_time.tv_sec - last_time.tv_sec) * 1000000 +
	(this_time.tv_usec - last_time.tv_usec);
      acct_bins[lastAcct] += elapsedMicroseconds;
    }
  lastAcct = n;
  last_time = this_time;
#endif
#endif
}


void
PrintAcct (char *name, int tid)
{
#ifdef ACCT
  float totalTime;
  float totalRead;
  float totalWrite;
  FILE *f;
  Filename fn;
  char hn[128];
  
  Acct(PROCESSING);		/* this makes sure the times are updated to
				   reflect everything up to this routine */
  sprintf(fn, "%s.%d.%d.acct", name, getpid(), tid);
  f = fopen(fn, "w");
  totalRead = acct_bins[READOPEN] +
    acct_bins[READING] +
    acct_bins[READCLOSE];
  totalWrite = acct_bins[WRITEOPEN] +
    acct_bins[WRITING] +
    acct_bins[WRITECLOSE];
  totalTime = acct_bins[PROCESSING] +
    totalRead +
    totalWrite +
    acct_bins[MSGWAITING];

  gethostname(hn, 128);
  fprintf(f, "\nHOST: %s\n", hn);
  fprintf(f, "Time spent: %.2f sec.:\n\n", totalTime / 1000000.0);
  fprintf(f, "\tProcessing :                 %.2f sec.   %.2f%%\n",
	  acct_bins[PROCESSING] / 1000000.0,
	  100.0 * acct_bins[PROCESSING] / totalTime);
  fprintf(f, "\tReading (total) :            %.2f sec.   %.2f%%\n",
	  totalRead / 1000000.0,
	  100.0 * totalRead / totalTime);
  fprintf(f, "\tWriting (total) :            %.2f sec.   %.2f%%\n",
	  totalWrite / 1000000.0,
	  100.0 * totalWrite / totalTime);
  fprintf(f, "\tMsg waiting :                %.2f sec    %.2f%%\n",
	  acct_bins[MSGWAITING] / 1000000.0,
	  100.0 * acct_bins[MSGWAITING] / totalTime);

  fprintf(f, "\tFile opening for read :      %.2f sec    %.2f%%\n",
	  acct_bins[READOPEN] / 1000000.0,
	  100.0 * acct_bins[READOPEN] / totalTime);
  fprintf(f, "\tFile reading :               %.2f sec    %.2f%%\n",
	  acct_bins[READING] / 1000000.0,
	  100.0 * acct_bins[READING] / totalTime);
  fprintf(f, "\tFile closing after read:     %.2f sec    %.2f%%\n",
	  acct_bins[READCLOSE] / 1000000.0,
	  100.0 * acct_bins[READCLOSE] / totalTime);

  fprintf(f, "\tFile opening for write :     %.2f sec    %.2f%%\n",
	  acct_bins[WRITEOPEN] / 1000000.0,
	  100.0 * acct_bins[WRITEOPEN] / totalTime);
  fprintf(f, "\tFile writing :               %.2f sec    %.2f%%\n",
	  acct_bins[WRITING] / 1000000.0,
	  100.0 * acct_bins[WRITING] / totalTime);
  fprintf(f, "\tFile closing after write:    %.2f sec    %.2f%%\n",
	  acct_bins[WRITECLOSE] / 1000000.0,
	  100.0 * acct_bins[WRITECLOSE] / totalTime);

  fclose(f);
#endif
}
