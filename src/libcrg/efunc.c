/************************************************************
 *                                                          *
 *  efunc.c                                                 *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1995 Christopher R. Genovese,          *
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
 *                                                          *
 *  This is intended to provide effective message routines  *
 *  and a flexible interface to standard library functions  *
 *  with convenient exception handling.                     *
 *                                                          *
 *  These functions are user configurable on a number of    *
 *  options.  See below for details.                        *
 *                                                          *
 *  Original Author:  Chris Genovese                        *
 *  Last Modified:    December 1995                         *
 *  Modified By:      Chris Genovese                        *
 *                                                          *
 ************************************************************/

#include   <stdio.h>
#include   <stdlib.h>
#include   <stdarg.h>
#include <string.h>
#include   <strings.h>
#include   <ctype.h>

#include   <errors.h>
#include   "stdcrg.h"

#include   <errno.h>
#include   <sys/stat.h>

/*
 * A basic error handling interface intended to make 
 * error control more convenient and readable.  Behavior
 * of the routines is configurable.  Also included are
 * robust versions of some common system calls for memory
 * allocation, stream handling, and block i/o.  Such functions
 * are named by putting an 'e' at the beginning of the 
 * standard name.
 *
 */
  

    /* Global Variables     */

enum efunc_options
{
    EFUNC_NOFATAL_ERR,              /* Errors are not fatal              */
    EFUNC_ALLOC_TRACK,              /* Allocs are tracked in output      */
    EFUNC_AUTO_NOCLOB,              /* Files not overwritten if exists   */
    EFUNC_FAIL_ON_EOF,              /* Read quota does fail on EOF       */
    EFUNC_AUTO_STDIN,               /* Name auto-interpreted as stdin    */
    EFUNC_AUTO_STDOUT,              /* Name auto-interpreted as stdout   */
    EFUNC_MESG_PREFIX,              /* Prefix string before all mesgs    */
    EFUNC_ERROR_FILE,               /* File name or ptr for error output */
    EFUNC_NUM_OPTIONS
};

#define   EFUNC_BUF_SIZE   1024   


static char       efunc_extra_buf[EFUNC_BUF_SIZE];

static char      *efunc_auto_stdin  = NULL;          /* File name implying auto-stdin (e.g., "-")  */
static char      *efunc_auto_stdout = NULL;          /* File name implying auto-stdout (e.g., "-") */
static char      *efunc_mesg_prefix = "";            /* Message Prefix                             */

static char      *efunc_opt_txt[] = {
                                        "no fatal error", 
                                        "alloc track", 
                                        "auto noclobber", 
                                        "fail on eof", 
                                        "auto stdin", 
                                        "auto stdout", 
                                        "mesg prefix", 
                                        "log file"
                                    };


static int        efunc_opt_arr[] = { 0, 0, 0, 0, 0, 0, 0, 0 };  /* All options default off */


static int       initialized=0;

static FILE      *efunc_fp = NULL;

extern int        errno;


    /* Source Code          */


/*
 * The following is generic interface for indicating  errors, warnings,
 * and internal errors in programs.  They can be replaced by more specialized
 * functions of the same name for specific applications.
 *
 * These functions have roughly the following interpretations:

     Function   Message Type                      Intended Audience
     --------   ------------                      -----------------
     Error      General errors                    User
     Warning    Warnings and non-fatal errors     User
     Message    Message without error             User
     Internal   Unexpected error deep in code     Installer/Developer
     SysError   Unexpected error in system call   Installer/Developer
 *   
 *
 * Note that \n need not be included as part of the string;
 * it is automatically appended.
 *
 * Also: File and line information can easily be supplied to internal
 * by using __FILE__ and __LINE__.  It might be useful to define
 * a macro WHERE  to be  __FILE__,__LINE__ for this purpose.
 *
 * The first argument to warning is intended to support easy
 * suppression of the message.  For instance, messages at different
 * levels can be provided by
 *
         Warning( DIAGNOSTIC, "General Message unless output suppressed" );
         Warning( DIAGNOSTIC == 1, "Message for 1" );
         Warning( DIAGNOSTIC == 2, "Message for 2" );
         Warning( DIAGNOSTIC == 3, "Message for 3" );
         Warning( DIAGNOSTIC == 4, "Message for 4" );

 * where DIAGNOSTIC is a C-preprocessor macro, potentially set on the command
 * line at compile time.  This is not ideal as a function is still called even
 * when the condition is not true.  A CPP front end does not support variable
 * length argument lists, so conditional compilation is necessary to prevent
 * the function call.
 *
 */

static void initialize()
{
    efunc_fp= stderr;
    initialized= 1;
}

void     Stdcrg_Error( char *fmt, ... )
{
    va_list ap;

    if (!initialized) initialize();

    if( efunc_fp != NULL )
    {
	va_start( ap, fmt );
	fprintf( efunc_fp, "%sError: ", efunc_mesg_prefix );    
	vfprintf( efunc_fp, fmt, ap );
	fprintf( efunc_fp, "\n" );

	va_end( ap );
    }
    
    efunc_exit( 1 );
}

void     Stdcrg_Warning( int diag, char *fmt, ... )
{
    va_list ap;

    if (!initialized) initialize();

    if( diag && (efunc_fp != NULL) )
    {
        va_start( ap, fmt );
        fprintf( efunc_fp, "%sWarning: ", efunc_mesg_prefix );
        vfprintf( efunc_fp, fmt, ap );
        fprintf( efunc_fp, "\n" );
    }

    va_end( ap );
}

void     Stdcrg_Internal( char *file, int line, char *fmt, ... )
{
    va_list ap;

    if (!initialized) initialize();

    if( efunc_fp != NULL )
    {
	va_start( ap, fmt );

	if( line > 0 )
	    fprintf( efunc_fp, "%sError in %s at line %d: ", efunc_mesg_prefix, file, line );
	else
	    fprintf( efunc_fp, "%sError in %s: ", efunc_mesg_prefix, file );

	vfprintf( efunc_fp, fmt, ap );
	fprintf( efunc_fp, "\n" );

	va_end( ap );
    }
    
    efunc_exit( 1 );
}

void     Stdcrg_SysError( char *file, int line, char *fmt, ... )
{
    char    *p;
    va_list ap;

    if (!initialized) initialize();

    if( efunc_fp != NULL )
    {
	va_start( ap, fmt );

	if( line > 0 )
	    sprintf( efunc_extra_buf,
		     "%sError in %s at line %d with errno %d: ", efunc_mesg_prefix, file, line, errno );
	else
	    sprintf( efunc_extra_buf,
		     "%sError in %s with errno %d: ", efunc_mesg_prefix, file, errno );

	vsprintf( efunc_extra_buf + strlen(efunc_extra_buf), fmt, ap );

	p = strerror( errno );

	if( p != NULL )
	    fprintf( efunc_fp, "%s (%s)\n", efunc_extra_buf, p );
	else
	    fprintf( efunc_fp, "%s\n", efunc_extra_buf );

	va_end( ap );
    }
    
    efunc_exit( 1 );
}


void     Stdcrg_Usage( char *fmt, ... )
{
    va_list ap;

    if (!initialized) initialize();

    if( efunc_fp != NULL )
    {
	va_start( ap, fmt );

	fprintf( efunc_fp, "%sUsage: ", efunc_mesg_prefix );
	vfprintf( efunc_fp, fmt, ap );
	fprintf( efunc_fp, "\n" );

	va_end( ap );
    }
    
    efunc_exit( 0 );
}

void     Stdcrg_Message( char *fmt, ... )
{
    va_list ap;

    if (!initialized) initialize();

    if( efunc_fp != NULL )
    {
	va_start( ap, fmt );

	fprintf( efunc_fp, "%s", efunc_mesg_prefix );  /* Message() does not add a \n to string */
	vfprintf( efunc_fp, fmt, ap );

	va_end( ap );
    }
}


void     efunc_exit( int code )
{
    if (!initialized) initialize();

    if( !efunc_opt_arr[EFUNC_NOFATAL_ERR] )
        exit( code );
}


/*
 * EFUNC_CONFIG can be used to configure the behavior of the error-handling
 * system.  String arguments are of the form "option = value" with additional
 * arguments as required.  The option texts are given in the array efunc_opt_txt[];
 * only partial strings are required, but no effort to disambiguate is made by
 * the function, so uniqueness is the caller's responsibility.  The table
 * in the commentary above efunc_opt_txt[] shows the acceptable choices for
 * value and additional arguments in each case.  
 *
 * A single call can include multiple settings separated by ";".  Thus:
 *
 *   "option1 = value1; option2 = value2; ..."
 *
 * is acceptable.  Whitespace before and after ';' and '=' is ignored.
 *
 */

int      efunc_config( char *str, ... )
{
    int        len, status = 1;

    char      *opt, *val, *p, *pcpy, *fn;

    va_list   ap;


    va_start( ap, str );

    if (!initialized) initialize();

    pcpy = p = strdup( str );     /* Make private copy that we can modify */
    
    p = strtok( p, ";" );

    while( p != NULL )
    {
        /* Parse Token string for form  "option = value"       */
        /* White space is ignored before, after and around '=' */

        while( *p && isspace(*p) )
            p++;

        if( *p == '\0' )                /* Empty Token; ignore.  Get next and move on */
        {
            p = strtok( NULL, ";" );
            continue;
        }

        opt = p;

        while( *p && *p != '=' )
            p++;

        if( *p == '\0' )                /* No equal sign.  Return Error Indication   */
	{
	    status = 0;
	    break;
	}

        val = p + 1;                    /* Save entry point to value    */

        for( p--; isspace( *p ); p-- )  /* Set ending to opt, no spaces */
            ;
        *(p + 1) = '\0';

        while( *val && isspace(*val) )  /* Search forward for value beg */
            val++;                   

        if( *val == '\0' )              /* No option value. Return Error Indication  */
	{
	    status = 0;
	    break;
	}

        for( p = val; *p; p++ )         /* Eliminate Ending space in value string.   */
            ;

        for( p--; isspace(*p); p-- )    /* Loop Automatically Sentinelled            */
            ;

        len = strlen( opt );            /* Partial strings fine, but no disambiguation done */


        /* Match tokens found and process options */

        if( !strncasecmp( opt, efunc_opt_txt[EFUNC_NOFATAL_ERR], len ) )
        {
            if( !strcasecmp( val, "true" ) || !strcasecmp( val, "1" ) )
                efunc_opt_arr[EFUNC_NOFATAL_ERR] = 1;
            else if( !strcasecmp( val, "false" ) || !strcasecmp( val, "0" ) )
                efunc_opt_arr[EFUNC_NOFATAL_ERR] = 0;
            else if( !strcasecmp( val, "%d" ) )
                efunc_opt_arr[EFUNC_NOFATAL_ERR] = va_arg( ap, int );
            else
	    {
		status = 0;
		break;
	    }
        }
        else if( !strncasecmp( opt, efunc_opt_txt[EFUNC_ALLOC_TRACK], len ) )
        {
            if( !strcasecmp( val, "true" ) || !strcasecmp( val, "1" ) )
                efunc_opt_arr[EFUNC_ALLOC_TRACK] = 1;
            else if( !strcasecmp( val, "false" ) || !strcasecmp( val, "0" ) )
                efunc_opt_arr[EFUNC_ALLOC_TRACK] = 0;
            else if( !strcasecmp( val, "%d" ) )
                efunc_opt_arr[EFUNC_ALLOC_TRACK] = va_arg( ap, int);
            else
	    {
		status = 0;
		break;
	    }
        }
        else if( !strncasecmp( opt, efunc_opt_txt[EFUNC_AUTO_NOCLOB], len ) )
        {
            if( !strcasecmp( val, "return" ) || !strcasecmp( val, "2" ) )
                efunc_opt_arr[EFUNC_AUTO_NOCLOB] = 2;
            else if( !strcasecmp( val, "fatal" ) || !strcasecmp( val, "on" ) || !strcasecmp( val, "1" ) )
                efunc_opt_arr[EFUNC_AUTO_NOCLOB] = 1;
            else if( !strcasecmp( val, "off" ) || !strcasecmp( val, "0" ) )
                efunc_opt_arr[EFUNC_AUTO_NOCLOB] = 0;
            else if( !strcasecmp( val, "%d" ) )
                efunc_opt_arr[EFUNC_AUTO_NOCLOB] = va_arg( ap, int);
            else
	    {
		status = 0;
		break;
	    }
        }
        else if( !strncasecmp( opt, efunc_opt_txt[EFUNC_FAIL_ON_EOF], len ) )
        {
            if( !strcasecmp( val, "true" ) || !strcasecmp( val, "1" ) )
                efunc_opt_arr[EFUNC_FAIL_ON_EOF] = 1;
            else if( !strcasecmp( val, "false" ) || !strcasecmp( val, "0" ) )
                efunc_opt_arr[EFUNC_FAIL_ON_EOF] = 0;
            else if( !strcasecmp( val, "%d" ) )
                efunc_opt_arr[EFUNC_FAIL_ON_EOF] = va_arg( ap, int);
            else
	    {
		status = 0;
		break;
	    }
        }
        else if( !strncasecmp( opt, efunc_opt_txt[EFUNC_AUTO_STDIN], len ) )
        {
            if( efunc_opt_arr[EFUNC_AUTO_STDIN] )
                free( efunc_auto_stdin );

            efunc_auto_stdin = strdup( val );
            efunc_opt_arr[EFUNC_AUTO_STDIN] = 1;
        }
        else if( !strncasecmp( opt, efunc_opt_txt[EFUNC_AUTO_STDOUT], len ) )
        {
            if( efunc_opt_arr[EFUNC_AUTO_STDOUT] )
                free( efunc_auto_stdout );

            efunc_auto_stdout = strdup( val );
            efunc_opt_arr[EFUNC_AUTO_STDOUT] = 1;
        }
        else if( !strncasecmp( opt, efunc_opt_txt[EFUNC_MESG_PREFIX], len ) )
        {
            if( efunc_opt_arr[EFUNC_MESG_PREFIX] )
                free( efunc_mesg_prefix );

            efunc_mesg_prefix = strdup( val );
            efunc_opt_arr[EFUNC_MESG_PREFIX] = 1;
        }
        else if( !strncasecmp( opt, efunc_opt_txt[EFUNC_ERROR_FILE], len ) )
        {
            if( !strcasecmp( val, "%p" ) )        /* Opened (but not written) file or NULL to disable */
                efunc_fp = va_arg( ap, FILE * ); 
	    else if( !strcasecmp( val, "%s" ) )      
	    {
		fn = va_arg( ap, char * );        /* Given name, re-open into file */

		efunc_fp = freopen( fn, "w", efunc_fp );       

		if( !efunc_fp )
		    SysError( "freopen", -1, "Cannot open file %s\n", fn );
	    }
	    else if( !strcasecmp( val, ">%s" ) )      
	    {
		fn = va_arg( ap, char * );        /* Given name, re-open append onto file */

		efunc_fp = freopen( fn, "a", efunc_fp );       

		if( !efunc_fp )
		    SysError( "freopen", -1, "Cannot open file %s\n", fn );
	    }
            else                                  /* Open new file, overwriting if necessary */
	    {
                efunc_fp = freopen( val, "w", efunc_fp );

		if( !efunc_fp )
		    SysError( "freopen", -1, "Cannot open file %s\n", val );
	    }

	    if( efunc_fp )
		setvbuf( efunc_fp, NULL, _IOLBF, BUFSIZ );    /* Line buffering for log file */
        }
        else
	{
	    status = 0;
	    break;
	}

        p = strtok( NULL, ";" );
    }

    va_end( ap );

    Free( pcpy );

    return( status );
}


/***  Memory Allocation Routines  ***/

void *emalloc( size_t n )
{
    void  *p;

    if (!initialized) initialize();

    if( !( p = malloc( n ) ) )
        SysError( "emalloc", NOLINE, "%d bytes requested\n", n );

    if( efunc_opt_arr[EFUNC_ALLOC_TRACK] )
        fprintf( efunc_fp, "%sTracking emalloc: %12d bytes requested\n", efunc_mesg_prefix, (int)n );

    return( p );
}

void *ecalloc( size_t n, size_t size )
{
    void  *p;

    if (!initialized) initialize();

    if( !( p = calloc( n, size ) ) )
        SysError( "ecalloc", -1, "%d elements of size %d requested\n", 
		  n, size );
 
    if( efunc_opt_arr[EFUNC_ALLOC_TRACK] )
        fprintf( efunc_fp, "%sTracking ecalloc: %12d bytes requested (%d elements of size %d)\n",
                            efunc_mesg_prefix, (int)(n*size), (int)n, (int)size );

    return( p );
}

void *erealloc( void *ptr, size_t n )
{
    void     *p;

    if (!initialized) initialize();

    if( !( p = realloc( ptr, n ) ) )
        SysError( "erealloc", -1, "Unable to resize block of %d bytes.\n", 
		  n );

    return( p );
}

void   **matrix( long n, long p, size_t unit )
{
    int    i;
    char   **A;
  
    if (!initialized) initialize();

    A = (char **)ecalloc( n,     sizeof(void *) );
    *A = (char *)ecalloc( n * p, unit );
  
    for( i = 1; i < n; i++ )
	A[i] = A[0]  +  i * p * unit;
  
    return( (void **)A );
}


/***  Stream Handling Routines ***/

static struct stat   statbuf;

FILE *efopen( char *fn, char *mode )
{
    int    i;
    FILE *fp;

    if (!initialized) initialize();

    if( mode[0] == 'r' && efunc_opt_arr[EFUNC_AUTO_STDIN] && !strcmp( fn, efunc_auto_stdin ) )
    {
        if( strchr( mode, '+' ) )
            Error( "Update requested for opening stdin (mode = %s).", mode );
        
        return( stdin );
    }
    else if( (mode[0] == 'w' || mode[0] == 'a') && efunc_opt_arr[EFUNC_AUTO_STDOUT]
                                                && !strcmp( fn, efunc_auto_stdout ) )
    {
        return( stdout );
    }
    
    i = stat( fn, &statbuf );

    if( efunc_opt_arr[EFUNC_AUTO_NOCLOB] && (mode[0] == 'w') )
    {
        /* Check if file exists but flag other problems that appear */

        if( i && (efunc_opt_arr[EFUNC_AUTO_NOCLOB] == 2) )
            return( NULL );
        else if( i )
            Error( "Attempt to efopen existing file (%s) for writing when auto noclob set.", fn );
        else if( errno != ENOENT )
            SysError( "efopen", -1, "Cannot open file %s\n", fn );
    }

    if (mode[0]=='r' && i) 
      SysError( "efopen", -1, "File %s not found!\n", fn );

    if( !(fp = fopen( fn, mode )) ) {
      if (mode[0]=='r') {
        SysError( "efopen", -1, "Cannot read file %s\n", fn );
      }
      else {
        SysError( "efopen", -1, "Cannot open file %s\n", fn );
      }
    }

    
    return( fp );
}

int efclose( FILE *fp )
{
    if (!initialized) initialize();

    if( fclose( fp ) )
        SysError( "efclose", -1, "Cannot close file\n" );

    return( 0 );
}


/***  Block I/O Routines ***/

size_t efread( void *ptr, size_t size, size_t nitems, FILE *stream )
{
    int   r;

    if (!initialized) initialize();

    r = fread( ptr, size, nitems, stream );

    if( r != nitems && (efunc_opt_arr[EFUNC_FAIL_ON_EOF] || !feof( stream )) )
        Internal( "efread", -1, "Read only %d out of %d items.", r, nitems );

    return( r );
}


size_t efwrite( void *ptr, size_t size, size_t nitems, FILE *stream )
{
    int     w;

    if (!initialized) initialize();

    w = fwrite( ptr, size, nitems, stream );

    if( w != nitems )
        Internal( "efwrite", -1, "Wrote only %d out of %d items.", w, nitems );

    return( w );
}

