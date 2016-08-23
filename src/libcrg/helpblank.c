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

static char rcsid[] = "$Id: helpblank.c,v 1.2 1999/07/07 23:16:42 welling Exp $";

#include <stdio.h>
#include <stdlib.h>

/*
 * Below are dummy versions of Help_init() and Help_listtopics().
 * If not defined elsewhere (i.e., in the output of buildhelp),
 * these will be used.  If they are defined elsewhere, these
 * will not be linked.
 *
 * Note: this works because these functions are compiled into a library.
 *
 */
 
void   Help_init( void )
{
    fprintf( stderr, "Help system not available; consult other documentation for usage information.\n" );
    exit( 0 );
}

void   Help_listtopics( void )
{
}
