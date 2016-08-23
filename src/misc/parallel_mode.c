/************************************************************
 *                                                          *
 *  parallel_mode.c                                        *
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
 *  Original programming by Joel Welling 7-03               *
 ************************************************************/
/* All this does is return a string indicating the type of 
 * parallel support compiled in.
 */

#include <stdio.h>

static char rcsid[] = "$Id: parallel_mode.c,v 1.1 2003/07/15 18:10:02 welling Exp $";

int main(int argc, char* argv[])
{
#ifdef PVM
  printf("PVM\n");
#else
#ifdef MPI
  printf("MPI\n");
#else
  printf("NONE\n");
#endif /* MPI */
#endif /* PVM */

  return 0;
}

