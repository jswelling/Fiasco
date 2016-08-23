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

static char rcsid[] = "$Id: tmperror.c,v 1.4 1999/07/07 23:12:41 welling Exp $";

#include   <stdio.h>
#include   <stdlib.h>
#include   <stddef.h>
#include   <stdarg.h>
#include   <string.h>

#include  "stdcrg.h"

#include   "tmperror.h"

/*

  These are completely temporary stubs that are replaced by the
  ErrorStream facility in version 1.

  Got it?  TEMPORARY, NOT REAL, JUST FOR NOW, A SOUPED UP INTEGER STATUS CODE

  Since it would be a pain to install the error streams in the beta version,
  they are represented here as a pointer to an int.  This can be constructed by
  es_new() and deconstructed by es_delete() or simply given a pointer
  to an existing integer.  This is a subset of the interface that covers
  what we need for now.

  ErrorStreams are used throughout version 1, including by the optimizer.
  The actual ErrorStream object is more complicated and es_new and
  es_delete are always used in version 1.

 */
   

ErrorStream     es_new( void )
{
    return( (int*)calloc( 1, sizeof(int) ) );
}

/*  ATTN: Only use on one constructed by es_new */

void            es_delete( ErrorStream* pes )
{
    if( *pes )
        Free( *pes );
    
    *pes = 0;
}

int             es_any_errors( ErrorStream es )
{
    return( *es > 0 );
}

void            es_write_mesg( ErrorStream es, es_level_type level, const char* format, ... )
{
    va_list     ap;

    *es = (int)level + 1;

    va_start( ap, format );
    vfprintf( stderr, format, ap );
    va_end( ap );

    fprintf( stderr, "\n" );

    if( level == ES_FATAL )
        exit( -1 );
}

void            es_write_error( ErrorStream es, es_level_type level, int code, const char* format, ... )
{
    va_list     ap;

    *es = code;

    va_start( ap, format );
    vfprintf( stderr, format, ap );
    va_end( ap );

    fprintf( stderr, "\n" );

    if( level == ES_FATAL )
        exit( code );
}

void            es_set_action( ErrorStream es, es_level_type level, es_action_type action )
{
}

void            es_set_stream( ErrorStream es, es_level_type level, FILE* fstream )
{
}

void            es_set_prefix( ErrorStream es, es_level_type level, char* prefix )
{
}


int             es_clear( ErrorStream es )
{
    int         last;

    last = *es > 0;
    *es = 0;

    return( last );
}

