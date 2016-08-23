/*
 *	worker.c
 *
 *    Interface to worker routines in the parallel spiral code.  The
 *    routines themselves have moved elsewhere to facilitate reuse.
 *
 *    Copyright (c) 1998 by Douglas C. Noll and the University of Pittsburgh and
 *	  the Pittsburgh Supercomputing Center 
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
 *    HISTORY
 *	1/98 - split off from spiral.c (Greg Hood, PSC)
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mri.h"
#include "fmri.h"
#include "par.h"
#include "bio.h"
#include "array.h"
#include "acct.h"
#include "stdcrg.h"
#include "spiral.h"

static char rcsid[] = "$Id: worker.c,v 1.10 2003/04/11 23:45:40 welling Exp $";

void
WorkerFinalize()
{
  WorkerShutdown();
}

void
WorkerContext ()
{
  static int first_context = TRUE;

  if (first_context)
    {
      verbose = TRUE;
      r.slice_sampim = NULL;
      first_context = FALSE;
    }

  WorkerInitialize();
}

void
WorkerTask ()
{
  if (c.lin_map || c.gen_map)
    GenerateReference();
  else
    GenerateImage();
}

