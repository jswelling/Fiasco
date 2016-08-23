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

static char rcsid[] = "$Id: sorting.c,v 1.2 1999/07/07 23:18:51 welling Exp $";

#include <string.h>

/*
 * NUMERICALLY is used with the library routine qsort to sort values
 * numerically.   Assumes type double, float, long, int, short as suffix
 * is either d, f, l, i, s.
 *
 */
 
int  Numerically_d( const void *a, const void *b )
{
        
    if( *(double *)a < *(double *)b )
        return( -1 );
    else if( *(double *)a > *(double *)b )
        return( 1 );
    else
        return( 0 );
}

int  Numerically_f( const void *a, const void *b )
{
        
    if( *(float *)a < *(float *)b )
        return( -1 );
    else if( *(float *)a > *(float *)b )
        return( 1 );
    else
        return( 0 );
}

int  Numerically_l( const void *a, const void *b )
{
        
    if( *(long *)a < *(long *)b )
        return( -1 );
    else if( *(long *)a > *(long *)b )
        return( 1 );
    else
        return( 0 );
}

int  Numerically_i( const void *a, const void *b )
{
        
    if( *(int *)a < *(int *)b )
        return( -1 );
    else if( *(int *)a > *(int *)b )
        return( 1 );
    else
        return( 0 );
}

int  Numerically_s( const void *a, const void *b )
{
        
    if( *(short *)a < *(short *)b )
        return( -1 );
    else if( *(short *)a > *(short *)b )
        return( 1 );
    else
        return( 0 );
}


/*
 * ALPHABETICALLY is used with the library routine qsort to sort values
 * alphabetically.
 *
 */
 
int  Alphabetically( const void *a, const void *b )
{
   return( strcmp( (char *)a, (char *)b ) );
}

