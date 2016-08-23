/************************************************************
 *                                                          *
 *  pfileorder.c                                            *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1997 Department of Psychology and      *
 *                        Department of Statistics,         *
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
 *  Original programming by Kate Fissell  1997              *
 *     12/97: modified for Fiasco, Joel Welling             *
 ************************************************************/
/****************************************************************************/
/*                                                                          */
/* File:       pfileorder.c                                                 */
/*                                                                          */
/* Desc.:       Take a reference file name and an arbitrarily ordered list  */
/*              of filenames including the reference filename.  Print out   */
/*              the list of filenames (not including the reference filename)*/
/*              in numerical order, starting with the first filename >      */
/*              than the reference file name and continuing in ascending    */
/*              order and wrapping around to include filenames < the        */
/*              reference file name at the end.  Used to create the         */
/*              correct list of pfiles as input to the sgrid program.       */
/*                                                                          */
/* Args:        reference_filename all_filenames...                         */
/*              eg pfileorder reffile `ls P*`                               */
/*                                                                          */
/* Ret:         prints list to stdout, on error no output.                  */
/*                                                                          */
/* NOTE:        This code is designed to handle a single wrap case,         */
/*              where the correct ordering of files might be:               */
/*              P30000 P40000 P50000 P60000 P00000 P10000                   */
/*              it does _NOT_ handle a possible double wrap case where      */
/*              the correct ordering is:                                    */
/*              P30000 P40000 P50000 P60000 P00000 P10000 P35000            */
/*                                                                          */
/*                                                                          */
/* Kate Fissell 4/9/97                                                      */
/*                                                                          */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"
 
static char rcsid[] = "$Id: pfileorder.c,v 1.5 2003/02/07 21:25:24 welling Exp $";

int compare();
 
int main(int argc, char* argv[])
{
  int i;
  int refind;
 
  /* Check to see if help was requested */
  if( ( argc > 1 ) && !strcmp( argv[1], "-help" ) )
    {
      if( argc == 2 )
	Help( "selecttopic" );
      else
	Help( argv[2] );
    }

  /* must have min 2 args */
  if (argc < 3)
    exit(1);
 
  /* make sure input is in order */
  qsort(&argv[2], argc-2, sizeof (char *), compare);
 
 
  /* find reffile */
  for (refind=-1,i=2; i<argc; i++)
    if (!strcmp(argv[1],argv[i])) {
      refind = i;
      break;
    }
  if (refind == -1) {
    Abort("%s: reference file name not repeated!\n",argv[0]);
  }
 
  /* print after reffile to end, then beginning up to reffile, not reffile */
  for (i=refind; i<argc; i++)
    fprintf(stdout, "%s ",argv[i]);
  for (i=2; i<refind; i++)
    fprintf(stdout, "%s ",argv[i]);
  
  return(0);
}
 
 
/* numerical comparison of 2 strings, formatted as 1 char followed 
 * by numerals
*/
int compare(s1,s2)
char **s1, **s2;
{
  float f1,f2;
 
  f1=atof(*s1+1);
  f2=atof(*s2+1);
 
  if (f1<f2) return -1;
  if (f1==f2) return 0;
  return 1;
}

