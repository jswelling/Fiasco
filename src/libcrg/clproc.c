#include   <stdio.h>
#include   <stdlib.h>
#include   <stdarg.h>
#include <string.h>
#include   <strings.h>
#include   <ctype.h>

#include   <errors.h>

#include   <time.h>

/************************************************************
 *                                                          *
 *  clproc.c                                                *
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
 *                                                          *
 *  Original Author:  Chris Genovese                        *
 *  Last Modified:    March 1996                            *
 *  Modified By:      Chris Genovese                        *
 *                                                          *
 ************************************************************/

/*** Function Prototypes ***/

#include "stdcrg.h"

/***
 * Implemented functions:
void     cl_scan( int argc, char **argv );                    // Read Command Line and Prepare   //
int      cl_get( char *name, char *fmt, ... );                // Get Argument or Flag Values     //
void     cl_configure( char *name, char *val );                  // Set default interpretations     //
void     cl_unmatched( void );                                // Find unmatched args; not implemented             //
int      cl_grep( char *pattern, char **arg );                // Find arg beginning with pattern- not implemented //
void     cl_mark( char *prefix, char *tag, char *showtime );  // Echo Information to Log File    //
int      cl_cmd_line( int *, char *** );                      // Copy command line structure     //
void     cl_cleanup( void );                                  // Clean-up Internal CL Info       //
int cl_cleanup_check( void );                                 // Clean-up Cl info, w/ error code //
int cl_present( char* flag_op );                              // check presence of flag //
****/

static int     get_next_arg( char *names, int class, int start, char **ptr, ... );
static char   *cnvt_code( char *s, int *suppress, int *size, int *type );
static char   *cnvt_def( char *s );
static int     assoc_type( int ftype, int size );
static void    classify_args( void );

/*** Macros and Constants ***/

enum types
{
    VOID_T, CHAR_T, SHORT_T, INT_T, LONG_T, FLOAT_T, DOUBLE_T, CHARPTR_T,
    UNSIGNED_CHAR_T, UNSIGNED_SHORT_T, UNSIGNED_INT_T, UNSIGNED_LONG_T
};

enum sizes     { STANDARD_S, CHAR_S, SHORT_S, LONG_S, LONG_DOUBLE_S };
enum ftypes    { VOID_F, CHAR_F, INTEGER_F, UNSIGNED_F, REAL_F, STRING_F };

enum class     { FLAG, MULTI, OPTION, NONDELIMITED, ONLY_NONDELIMITED };
enum options   { TAGSET, FLAG_MATCH, FLAG_TYPE, MULTI_MATCH, OPTION_SEP, CASE_FOLD, NUM_OPTIONS };
enum match     { MATCH_FALSE, MATCH_TRUE };
enum argtype   { DELIM=1, NONDELIM, BOTH };


/*** Global Variables ***/

static   int     cl_argc = 0;
static   char  **cl_argv = NULL;
static   int    *cl_argseen;
static   int     cl_lastarg = -1;

static   int     tagset_len = 0;
static   char   **tagset = NULL;

static   int     case_fold = 0;
static   char  **option_sep_set = NULL;
static   int     opt_sep_set_len = 0;
static   int    *flag_match = NULL;
static   int    *multi_match = NULL;
static   int    *local_flag_match = NULL;

    /* Scratch buffers, not much space is needed */

static   const int cl_bufsize = 256;

static   void   *cl_buf  = NULL;
static   void   *cl_vbuf = NULL;
static   char   *cl_cbuf = NULL;

    /* Data structures with arg classifications */

static   int    *cl_delimind = NULL;
static   int    *cl_delimtag = NULL;
static   int    *cl_nondelimind = NULL;
static   int    *cl_argtype = NULL;

static   char   **cl_delimarg = NULL;
static   char   **cl_nondelimarg = NULL;

    /* Option names for cl_config */

static char *cl_options[] = {
                              "tag-set",
                              "flag-match",
                              "flag-type",
                              "multi-match",
                              "option-sep",
                              "case-fold"
                            };

static int debug= 0;

/*
 *
 * cl_scan peruses the command line structure for later processing.  Typically,
 * this is invoked with cl_scan( argc, argv ).  (Eventually (but not yet), I may implement
 * a unique structure for each scan; I don't see an application for this at this time.)
 *
 * cl_config allows various options and defaults to be changed, including:
 *
 *    1. Tag Set                               (Default -)
 *    2. Flag tags for true or false settings  (Default -=T,=F)
 *    3. Flag type                             (Default int)
 *    4. Multi tags for true false             (Default -=T,=F)
 *    5. Option separator Set                  (Default "", implies space or no-space both acceptable)
 *    6. Case sensitive?                       (Default 1)
 *
 * Note that the tag set specifies the set of strings that can begin a
 * non-delimited argument.  Two common examples of tag set are "-|+" and
 * "-|-no-". The flag tags then associate to each possibility a value.  The
 * flag type is the type of the 0,1 variable for a flag value.  For most
 * applications, these will be sufficient.  For others, special choices can
 * be set; for example, Flag tags (+=T,-=F,=F), Tags "-|+", and Option
 * Separator "=|" is another reasonable regime.  This function is invoked
 * with
 *
 *      cl_config( name, value )
 *
 * where names is an option name (see below) and value is a string as above.  Some
 * of these options can be set locally in format strings.  See cl_get().
 *
 * cl_get gets off the command line using the given interpretation.
 * There is a quite flexible syntax as described below.  Values are
 * read off a scanf-like format string.  This allows options, flags,
 * and other arguments. This is invoked with
 *
 *    cl_get( name, formatstring, pointers... ).
 *
 * The name is of the form "name1|..." which specifies synonymous
 * names for the option (e.g., "h|help"), for, say, long and short
 * forms.
 *
 * The format string takes the form  "%class[optional qualifiers] ...". See below
 *
 * The pointers are an arbitrary list of pointers, one per % specification in
 * the format string.
 *
 * Format Classes
 * --------------
 *
 *   %flag             Stand-alone flags (on/off switches)         
 *   %multi            Single-character flags that can appear together with a bunch of others
 *   %option           An option that takes one or more arguments
 *   scanf controls    Corresponds to values to be read into pointers (%d, %f,...)
 *                     Some scanf() controls are not supported, see below.
 *
 * 
 * Qualifiers
 * ----------
 * Qualifiers can follow any of the class specifications including scanf controls.
 * Qualifiers are strings surrounded by [ ]'s with no space between the opening
 * bracket and the class name.  Qualifiers allow the specification of default
 * values and local option settings.  See examples and details below.
 *
 * If the same argument is given multiple times, only the next on the command line
 * is used.  This allows for repeated arguments; in that case, the number of the appearence
 * is made available for assignment.
 *
 * Some examples:
 *
 * 1. cl_get( "f|full", "%flag[-=F,+=T,=T|-=T,-no-=F,=T]", &full_flag );
 *
 * This is sensitive to a standalone flag of the form -f or -full.
 * The tags in the first part indicates that -f means set full_flag to
 * 0, and +f means set full_flag to 1.  The empty string in the tag
 * description indicates a default value if the flag is not present.
 * The second part corresponds to the -full flag.  This indicates that
 * -full and -no-full correspond to 1 and 0 respectively, again with a
 * default value of true.
 *
 * If the two empty specifications disagree, the first one is used.
 *
 * 2. cl_get( "f|full", "%flag", &full_flag )
 *
 * This is as above but uses the default settings directly.
 *
 * 3. cl_get( "f|full", "%multi", &full_flag )
 *
 * This is as above but allows the -f to be a multi while -full is a standalone flag.
 *
 * 4. cl_get( "a", "%multi[-=T,=F]", &a_flag )
 *
 * This specifies -a as part of a multi char flag.  For instance "foo -abcx" or "foo -xvaf"
 *
 * 5. cl_get( "s|sep", "%option %d %lf[0.0]", &sep_int_val, &sep_double_val )
 *
 * specifies an option -sep or -s that takes two arguments one an integer
 * and the other a double with a default value of 0.  Default values are parsed
 * specially, so if the argument after the integer does not match a floating point number
 * the default value is used.
 *
 * 6. cl_get( "", "%d %lf %s %s", &arg_int, &arg_double, &arg_str1, &arg_str2 )
 *
 * Gets the next non-delimited arguments in order, requiring that they match
 * the given specification.  Default values are allowed and give more flexibility
 * in the matching. 
 *
 * Note: Order of getting the arguments is important in that flags and options should
 * be read before non-delimited arguments.  This does not pose a real burden.
 *
 *
 * cl_next is a utility routine designed to give more direct control over argument
 * parsing.  This gets the argument following that most recently cl_get'd, marking
 * it as seen.  This is useful for dealing with contingencies where another part
 * of an option is to be given only when a certain condition is met.
 *
 * cl_mark and cl_cleanup are utility routines.
 *
 *
 * VARARGS Note: these routines are set up to use the stdlib.h variadic
 *               argument list formalism of Standard C, rather than
 *               varargs.h or Kernigan and Richie C (K&R C).  Specifically,
 *               it is assumed that floats are promoted to doubles and
 *               (signed or unsigned) chars or shorts are promoted to
 *               signed ints when passed in the variable part or a 
 *               parameter list to cl_get().
 *
 */

int cl_cleanup_check()
{
  int args_remain_unparsed= 0;
  int i;

  if (!(cl_argc && cl_argv)) /* not yet init'd */
    return 0;

  for (i=1; i<cl_argc; i++) {
    if (!cl_argseen[i]) {
      args_remain_unparsed= 1;
      break;
    }
  }

  cl_cleanup();
  return args_remain_unparsed;
}

int cl_present( char* flag_op )
{
  int flag= 0;
  cl_get( flag_op,"%flag", &flag);
  return flag;
}

void     cl_scan( int argc, char **argv )
{
    int     i;

    /* Re-set values if scanned before */

    if (getenv("F_CLPROC_DEBUG")) debug= 1;

    if( cl_argv )
    {
	for( i = 0; i < cl_argc; i++ )
	    Free( cl_argv[i] );

	Free( cl_argv );
	Free( cl_argseen );
	Free( cl_buf );
	Free( cl_vbuf );
	Free( cl_cbuf );

	Free( cl_delimind );
	Free( cl_nondelimind );
	Free( cl_delimtag );
	Free( cl_delimarg );
	Free( cl_nondelimarg );
	Free( cl_argtype );
    }

    if( tagset )
    {
	for( i = 0; i <= tagset_len; i++ )
	    Free( tagset[i] );

	Free( tagset );
	Free( flag_match );
	Free( local_flag_match );
	Free( multi_match );
    }

    if( option_sep_set )
    {
	for( i = 0; i < opt_sep_set_len; i++ )
	    Free( option_sep_set[i] );

	Free( option_sep_set );
    }

    /* Copy Command Line Structure  */
    
    cl_argc = argc;
    cl_argv = Malloc( argc, char * );
    cl_argseen = Calloc( argc, int );

    for( i = 0; i < cl_argc; i++ )
	cl_argv[i] = strdup( argv[i] );

    /* Prepare Scratch Buffers      */

    cl_buf  = (void *)Malloc( cl_bufsize, char );
    cl_vbuf = (void *)Malloc( cl_bufsize, char );
    cl_cbuf = Malloc( cl_bufsize, char );

    /* Precomputed Classifications  */

    cl_delimind = Calloc( cl_argc, int );
    cl_delimtag = Calloc( cl_argc, int );
    cl_nondelimind = Calloc( cl_argc, int );
    cl_argtype = Calloc( cl_argc, int );

    cl_delimarg = Calloc( cl_argc, char * );
    cl_nondelimarg = Calloc( cl_argc, char * );
    
    /* Set-up Tags and Flag Matches */

    tagset_len = 1;
    tagset = Malloc( tagset_len + 1, char * );
    tagset[0] = strdup( "-" );
    tagset[1] = strdup( "" );

    flag_match = Calloc( tagset_len + 1, int );      
    flag_match[0] = MATCH_TRUE;
    flag_match[1] = MATCH_FALSE;        /* Extra entry is value for missing flag */
	
    local_flag_match = Malloc( tagset_len + 1, int );

    multi_match = Calloc( tagset_len + 1, int );
    multi_match[0] = MATCH_TRUE;
    multi_match[1] = MATCH_FALSE;       /* Extra entry is value for missing flag */

    /* Set-up Option Separator List */

    opt_sep_set_len = 1;

    option_sep_set = Malloc( 1, char * );
    option_sep_set[0] = strdup( "" );       /* Space or nothing both ok */

    /* Scan and Classify Argument List */

    classify_args();
}


int      cl_get( char *name, char *fmt, ... )
{
    char      ch;

    int       i, j, nxtarg, ac = 0, start = 0, def_matched, len;
    int       suppress, size, ftype, type;

    char      rfmt[32], rfmtB[32];
    char      *p, *q, *arg;
    char      *def_value;

    va_list   ap;
    
    /* Holding variables */

    char      defv_c;
    short     defv_s;
    int       defv_i;
    long      defv_l;
    float     defv_f;
    double    defv_d;
    char     *defv_p;

    unsigned char      defv_uc;  /* These aren't really necessary but what the heck */
    unsigned short     defv_us;
    unsigned int       defv_ui;
    unsigned long      defv_ul;


    va_start( ap, fmt );

    rfmt[0] = '%';         /* Basic Format String */
    rfmt[1] = '\0';


    if( name[0] == '\0' )  /* Empty Name implies non-delimited argument */
    {
	p = strchr(fmt, '%');

	for( ac = 0; p; p = strchr(q, '%') )
	{
	    suppress = 0;         /* Reset for next argument */
	    size = STANDARD_S;
            ftype = VOID_F;
	    rfmt[1] = '\0';
	    def_matched = 0;

	    nxtarg = get_next_arg( "", NONDELIMITED, start, &arg );     /* Get Next Argument */

	    if( !(q = cnvt_code( p, &suppress, &size, &ftype )) )
		return( -ac );

	    def_value = cnvt_def( q + 1 );

	    if( suppress )
	    {
		/* Skip format string assume no spaces in default string */

		for( q = p; *q && !isspace(*q); q++ )
		    ;
		
		continue; /* There are no assignments or defaults in arg list */
	    }
	    else
		p++;

	    strncat( rfmt, p, 1 + q - p );     /* Set conversion format to use */
	    type = assoc_type( ftype, size );  /* Find underlying read type    */
	    
	    switch( type )
	    {
	      case CHAR_T:

		if (def_value) {
		  if( *def_value == '%' )          /* Default passed as arg */
		    defv_c = va_arg( ap, int );
		  else
		    sscanf( def_value, rfmt, &defv_c );
		}

		if (debug) {
		  if (def_value)
		    printf("#clproc:%s:undelimited%d::char:default=%c\n",
			   cl_argv[0],ac,defv_c);
		  else
		    printf("#clproc:%s:undelimited%d::char:default=(null)\n",
			   cl_argv[0],ac);
		}

		if( !arg && !def_value )
		    return( -ac );
		else if( !def_value )
		{
		    j = sscanf( arg, rfmt, va_arg( ap, char * ) );

		    if( j != 1 )
			return( -ac );
		}
		else
		{
		    strcpy( rfmtB, rfmt );
		    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
		    /* Check if matching argument found */

		    if( arg )
			j = sscanf( arg, rfmtB, (char *)cl_vbuf, (char *)cl_buf ); 

		    if( arg && j == 1 )
			sscanf( arg, rfmt, va_arg( ap, char * ) );    /* Use Command Line Value  */
		    else
		    {
			*(va_arg( ap, char * )) = defv_c;             /* Use Given Default Value */
			def_matched = 1;
		    }
		}
		
		break;
		
	      case SHORT_T:

		if (def_value) {
		  if( *def_value == '%' )          /* Default passed as arg */
		    defv_s = va_arg( ap, int );  
		  else
		    sscanf( def_value, rfmt, &defv_s );
		}

		if (debug) {
		  if (def_value)
		    printf("#clproc:%s:undelimited%d::short:default=%d\n",
			   cl_argv[0],ac,(int)defv_s);
		  else
		    printf("#clproc:%s:undelimited%d::short:default=(null)\n",
			   cl_argv[0],ac);
		}

		if( !arg && !def_value )
		    return( -ac );
		else if( !def_value )
		{
		    j = sscanf( arg, rfmt, va_arg( ap, short * ) );

		    if( j != 1 )
			return( -ac );
		}
		else
		{
		    strcpy( rfmtB, rfmt );
		    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
		    /* Check if matching argument found */

		    if( arg )
			j = sscanf( arg, rfmtB, (short *)cl_vbuf, (char *)cl_buf ); 

		    if( arg && j == 1 )
			sscanf( arg, rfmt, va_arg( ap, short * ) );    /* Use Command Line Value  */
		    else
		    {
			*(va_arg( ap, short * )) = defv_s;             /* Use Given Default Value */
			def_matched = 1;
		    }
		}
		
		break;
		
	      case INT_T:

		if (def_value) {
		  if( *def_value == '%' )          /* Default passed as arg */
		    defv_i = va_arg( ap, int );
		  else
		    sscanf( def_value, rfmt, &defv_i );
		}

		if (debug) {
		  if (def_value)
		    printf("#clproc:%s:undelimited%d::int:default=%d\n",
			   cl_argv[0],ac,defv_i);
		  else
		    printf("#clproc:%s:undelimited%d::int:default=(null)\n",
			   cl_argv[0],ac);
		}

		if( !arg && !def_value )
		    return( -ac );
		else if( !def_value )
		{
		    j = sscanf( arg, rfmt, va_arg( ap, int * ) );

		    if( j != 1 )
			return( -ac );
		}
		else
		{
		    strcpy( rfmtB, rfmt );
		    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
		    /* Check if matching argument found */

		    if( arg )
			j = sscanf( arg, rfmtB, (int *)cl_vbuf, (char *)cl_buf ); 

		    if( arg && j == 1 )
			sscanf( arg, rfmt, va_arg( ap, int * ) );    /* Use Command Line Value  */
		    else
		    {
			*(va_arg( ap, int * )) = defv_i;             /* Use Given Default Value */
			def_matched = 1;
		    }
		}
		
		break;

	      case LONG_T:

		if (def_value) {
		    if( *def_value == '%' )        /* Default passed as arg */
			defv_l = va_arg( ap, long );
		    else
			sscanf( def_value, rfmt, &defv_l );
		}

		if (debug) {
		  if (def_value)
		    printf("#clproc:%s:undelimited%d::long:default=%ld\n",
			   cl_argv[0],ac,defv_l);
		  else
		    printf("#clproc:%s:undelimited%d::long:default=(null)\n",
			   cl_argv[0],ac);
		}

		if( !arg && !def_value )
		    return( -ac );
		else if( !def_value )
		{
		    j = sscanf( arg, rfmt, va_arg( ap, long * ) );

		    if( j != 1 )
			return( -ac );
		}
		else
		{
		    strcpy( rfmtB, rfmt );
		    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
		    /* Check if matching argument found */

		    if( arg )
			j = sscanf( arg, rfmtB, (long *)cl_vbuf, (char *)cl_buf ); 

		    if( arg && j == 1 )
			sscanf( arg, rfmt, va_arg( ap, long * ) );    /* Use Command Line Value  */
		    else
		    {
			*(va_arg( ap, long * )) = defv_l;             /* Use Given Default Value */
			def_matched = 1;
		    }
		}
		
		break;

	      case FLOAT_T:

		if (def_value) {
		  if( *def_value == '%' )          /* Default passed as arg */
		    defv_f = va_arg( ap, double );   /* CAUTION: Value Converted!! */
		  else
		    sscanf( def_value, rfmt, &defv_f );
		}

		if (debug) {
		  if (def_value)
		    printf("#clproc:%s:undelimited%d::float:default=%g\n",
			   cl_argv[0],ac,defv_f);
		  else
		    printf("#clproc:%s:undelimited%d::float:default=(null)\n",
			   cl_argv[0],ac);
		}

		if( !arg && !def_value )
		    return( -ac );
		else if( !def_value )
		{
		    j = sscanf( arg, rfmt, va_arg( ap, float * ) );

		    if( j != 1 )
			return( -ac );
		}
		else
		{
		    strcpy( rfmtB, rfmt );
		    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
		    /* Check if matching argument found */

		    if( arg )
			j = sscanf( arg, rfmtB, (float *)cl_vbuf, (char *)cl_buf ); 

		    if( arg && j == 1 )
			sscanf( arg, rfmt, va_arg( ap, float * ) );    /* Use Command Line Value  */
		    else
		    {
			*(va_arg( ap, float * )) = defv_f;             /* Use Given Default Value */
			def_matched = 1;
		    }
		}
		
		break;

	      case DOUBLE_T:

		if (def_value) {
		    if( *def_value == '%' )         /* Default passed as arg */
			defv_d = va_arg( ap, double );
		    else
			sscanf( def_value, rfmt, &defv_d );
		}

		if (debug) {
		  if (def_value)
		    printf("#clproc:%s:undelimited%d::double:default=%g\n",
			   cl_argv[0],ac,defv_d);
		  else
		    printf("#clproc:%s:undelimited%d::double:default=(null)\n",
			   cl_argv[0],ac);
		}

		if( !arg && !def_value )
		    return( -ac );
		else if( !def_value )
		{
		    j = sscanf( arg, rfmt, va_arg( ap, double * ) );

		    if( j != 1 )
			return( -ac );
		}
		else
		{
		    strcpy( rfmtB, rfmt );
		    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
		    /* Check if matching argument found */

		    if( arg )
			j = sscanf( arg, rfmtB, (double *)cl_vbuf, (char *)cl_buf ); 

		    if( arg && j == 1 )
			sscanf( arg, rfmt, va_arg( ap, double * ) );    /* Use Command Line Value  */
		    else
		    {
			*(va_arg( ap, double * )) = defv_d;             /* Use Given Default Value */
			def_matched = 1;
		    }
		}
		
		break;

	      case CHARPTR_T:
		
		if (def_value) {
		  if( *def_value == '%' )          /* Default passed as arg */
		    defv_p = strdup( va_arg( ap, char * ) );
		  else
		    {
		      defv_p = strdup( def_value );
		      
		      for( i = 0; defv_p[i] && defv_p[i] != ']'; i++ );
		      defv_p[i] = '\0';
		    }
		}
		
		if (debug) {
		  if (def_value)
		    printf("#clproc:%s:undelimited%d::string:default=<%s>\n",
			   cl_argv[0],ac,defv_p);
		  else
		    printf("#clproc:%s:undelimited%d::string:default=<(null)>\n",
			   cl_argv[0],ac);
		}

		if( !arg && !def_value )
		    return( -ac );
		else if( !def_value )
		{
		    if (!strcmp(rfmt,"%s")) 
		      strcpy( va_arg( ap, char* ), arg );
		    else {
		      j = sscanf( arg, rfmt, va_arg( ap, char * ) );
		      if( j != 1 )
			return( -ac );
		    }
		}
		else
		{
		    /* Check if matching argument found */

		    if( arg )
			j = sscanf( arg, rfmt, (int *)cl_vbuf, (char *)cl_buf ); 

		    if( arg && j == 1 )
			strcpy( va_arg(ap, char*), arg );    /* Use Command Line Value  */
		    else
		    {
			strcpy( va_arg( ap, char * ), defv_p );       /* Use Given Default Value */
			def_matched = 1;
		    }

		    Free( defv_p );
		}

		break;

	      case UNSIGNED_CHAR_T:

		if (def_value) {
		    if( *def_value == '%' )         /* Default passed as arg */
			defv_uc = va_arg( ap, int );
		    else
			sscanf( def_value, rfmt, &defv_uc );
		}

		if (debug) {
		  if (def_value)
		    printf("#clproc:%s:undelimited%d::uchar:default=%d\n",
			   cl_argv[0],ac,(int)defv_uc);
		  else
		    printf("#clproc:%s:undelimited%d::uchar:default=(null)\n",
			   cl_argv[0],ac);
		}

		if( !arg && !def_value )
		    return( -ac );
		else if( !def_value )
		{
		    j = sscanf( arg, rfmt, va_arg( ap, unsigned char * ) );

		    if( j != 1 )
			return( -ac );
		}
		else
		{
		    /* Check if matching argument found */

		    if( arg )
			j = sscanf( arg, rfmt, (unsigned char *)cl_vbuf, (char *)cl_buf ); 

		    if( arg && j == 1 )
			sscanf( arg, rfmt, va_arg( ap, unsigned char * ) );   /* Use Command Line Value  */
		    else
		    {
			*(va_arg( ap, unsigned char * )) = defv_uc;           /* Use Given Default Value */
			def_matched = 1;
		    }
		}
		
		break;
		
	      case UNSIGNED_SHORT_T:

		if (def_value) {
		    if( *def_value == '%' )       /* Default passed as arg */
			defv_us = va_arg( ap, int );
		    else
			sscanf( def_value, rfmt, &defv_us );
		}

		if (debug) {
		  if (def_value)
		    printf("#clproc:%s:undelimited%d::ushort:default=%d\n",
			   cl_argv[0],ac,(int)defv_us);
		  else
		    printf("#clproc:%s:undelimited%d::ushort:default=(null)\n",
			   cl_argv[0],ac);
		}

		if( !arg && !def_value )
		    return( -ac );
		else if( !def_value )
		{
		    j = sscanf( arg, rfmt, va_arg( ap, unsigned short * ) );

		    if( j != 1 )
			return( -ac );
		}
		else
		{
		    strcpy( rfmtB, rfmt );
		    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
		    /* Check if matching argument found */

		    if( arg )
			j = sscanf( arg, rfmtB, (unsigned short *)cl_vbuf, (char *)cl_buf ); 

		    if( arg && j == 1 )
			sscanf( arg, rfmt, va_arg( ap, unsigned short * ) );   /* Use Command Line Value  */
		    else
		    {
			*(va_arg( ap, unsigned short * )) = defv_us;           /* Use Given Default Value */
			def_matched = 1;
		    }
		}
		
		break;
		
	      case UNSIGNED_INT_T:

		if (def_value) {
		  if( *def_value == '%' )     /* Default passed as arg */
		    defv_ui = va_arg( ap, unsigned int );
		  else
		    sscanf( def_value, rfmt, &defv_ui );
		}

		if (debug) {
		  if (def_value)
		    printf("#clproc:%s:undelimited%d::uint:default=%d\n",
			   cl_argv[0],ac,defv_ui);
		  else
		    printf("#clproc:%s:undelimited%d::uint:default=(null)\n",
			   cl_argv[0],ac);
		}

		if( !arg && !def_value )
		    return( -ac );
		else if( !def_value )
		{
		    j = sscanf( arg, rfmt, va_arg( ap, unsigned int * ) );

		    if( j != 1 )
			return( -ac );
		}
		else
		{
		    strcpy( rfmtB, rfmt );
		    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
		    /* Check if matching argument found */

		    if( arg )
			j = sscanf( arg, rfmtB, (unsigned int *)cl_vbuf, (char *)cl_buf ); 

		    if( arg && j == 1 )
			sscanf( arg, rfmt, va_arg( ap, unsigned int * ) );  /* Use Command Line Value  */
		    else
		    {
			*(va_arg( ap, unsigned int * )) = defv_ui;          /* Use Given Default Value */
			def_matched = 1;
		    }
		}
		
		break;

	      case UNSIGNED_LONG_T:

		if (def_value) {
		  if( *def_value == '%' )          /* Default passed as arg */
		    defv_ul = va_arg( ap, unsigned long );
		  else
		    sscanf( def_value, rfmt, &defv_ul );
		}

		if (debug) {
		  if (def_value)
		    printf("#clproc:%s:undelimited%d::ulong:default=%ld\n",
			   cl_argv[0],ac,defv_ul);
		  else
		    printf("#clproc:%s:undelimited%d::ulong:default=(null)\n",
			   cl_argv[0],ac);
		}

		if( !arg && !def_value )
		    return( -ac );
		else if( !def_value )
		{
		    j = sscanf( arg, rfmt, va_arg( ap, unsigned long * ) );

		    if( j != 1 )
			return( -ac );
		}
		else
		{
		    strcpy( rfmtB, rfmt );
		    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
		    /* Check if matching argument found */

		    if( arg )
			j = sscanf( arg, rfmtB, (unsigned long *)cl_vbuf, (char *)cl_buf ); 

		    if( arg && j == 1 )
			sscanf( arg, rfmt, va_arg( ap, unsigned long * ) );  /* Use Command Line Value  */
		    else
		    {
			*(va_arg( ap, unsigned long * )) = defv_ul;          /* Use Given Default Value */
			def_matched = 1;
		    }
		}
		
		break;

	      case VOID_T:
	      default:

		return( -ac );
	    }

	    /* Only mark an argument as seen if it was in fact processed */

	    if( !def_matched )
	    {
		cl_argseen[nxtarg] = 1;
		start = cl_lastarg = nxtarg;
	    }

	    /* Skip rest of conversion format, assume no spaces in default string */

	    for( q++; *q && !isspace(*q); q++ )
		;

	    ac++;    /* Defaults add to the count */
	}

	return( ac );
    }
    else if( !strncmp( fmt, "%flag", 5 )  )
    {
      if (debug) {
	printf("#clproc:%s:flag:%s::fmt=<%s>\n",cl_argv[0],name,fmt+5);
      }

	if( fmt[5] == '[' ) /* Get local tags */
	{
	    for( i = 0; i <= tagset_len; i++ )
		local_flag_match[i] = flag_match[i];

	    for( p = fmt + 6; p && *p;  ) /* Extract Tags */
	    {
		if( *p == '=' )
		{
		    j = sscanf( p, "=%c%[], \t]%n", &ch, (char*)cl_vbuf, &len ); /* = not allowed in tag */

		    if( j != 2 )
			break;

		    if( toupper(ch) == 'T' )
			local_flag_match[tagset_len] = MATCH_TRUE;
		    else if( toupper(ch) == 'F' )
			local_flag_match[tagset_len] = MATCH_FALSE;
		}
		else
		{                                       /* = not allowed in tag */
		    j = sscanf( p, "%[^=]=%c%[], \t]%n", (char*)cl_buf, &ch, (char*)cl_vbuf, &len );

		    if( j != 3 )
			break;
		
		    for( i = 0; i < tagset_len; i++ )
		    {
			if( !strcmp( tagset[i], cl_buf ) )
			    break;
		    }
		
		    if( toupper(ch) == 'T' )
			local_flag_match[i] = MATCH_TRUE;
		    else if( toupper(ch) == 'F' )
			local_flag_match[i] = MATCH_FALSE;
		}

		p += len;
	    }

	    nxtarg = get_next_arg( name, FLAG, 1, &arg, &j );     /* Get Flag Argument */

	    if( !arg )
		*(va_arg( ap, int * )) = local_flag_match[tagset_len];
	    else
	    {
		*(va_arg( ap, int * )) = local_flag_match[j];
		cl_argseen[nxtarg] = 1;
		cl_lastarg = nxtarg;
	    }
	}
	else
	{
	    nxtarg = get_next_arg( name, FLAG, 1, &arg, &j );     /* Get Flag Argument */

	    if( !arg )
		*(va_arg( ap, int * )) = flag_match[tagset_len];
	    else
	    {
		*(va_arg( ap, int * )) = flag_match[j];
		cl_argseen[nxtarg] = 1;
		cl_lastarg = nxtarg;
	    }
	}
	
	return( ac );
    }
    else if( !strncmp( fmt, "%multi", 6 ) )
    {
      if (debug) {
	printf("#clproc:%s:flag_multi:%s::fmt=<%s>\n",cl_argv[0],name,fmt+6);
      }

	nxtarg = get_next_arg( name, MULTI, 1, &arg, &j, &p );     /* Get Multi flag Argument */

	if( !arg )
	    *(va_arg( ap, int * )) = multi_match[tagset_len];
	else
	{
	    *(va_arg( ap, int * )) = multi_match[j];

	    while( *arg )            /* Eliminate flag from list */
	    {
		arg[0] = arg[1];
		arg++;
	    }

	    if( p[0] == '\0' )
		cl_argseen[nxtarg] = 1;

	    cl_lastarg = nxtarg;
	}
	
	return( ac );
    }
    else if( !strncmp( fmt, "%option", 7) )
    {
	nxtarg = get_next_arg( name, OPTION, 1, &arg, strchr(fmt + 7, '%'), &p );  /* Get Flag Argument */

	if( arg )
	{
	    cl_argseen[nxtarg] = 1;              /* We've seen the option */
	    cl_lastarg = nxtarg;

	    q = strchr( p, '|' );                /* Skip the name that matched */
	    len = ( q ) ? q - p : strlen( p );
	    arg += len;

	    if( *arg == '\0' )       /* No attached arg, move to next argument */
	    {
		start = nxtarg + 1;
		nxtarg = get_next_arg( "", ONLY_NONDELIMITED, start, &arg );
	    }
	    else                     /* Get attached argument */
	    {
		/* Check against option separators */

		for( j = 0; j < opt_sep_set_len; j++ )
		{
		    len = strlen( option_sep_set[j] );

		    if( !strncmp( arg, option_sep_set[j], len) )
			break;
		}

		if( j < opt_sep_set_len )
		    arg += len;
	    }
	}
	
	/* Skip to args, if any */

	for( p = fmt; *p && !isspace(*p); p++ )
	    ;

	if( p ) 
	    p = strchr(p, '%'); 

	/* Process Args */

	for( ac = 0; p; p = strchr(q, '%') )
	{
	    suppress = 0;         /* Reset for next argument */
	    size = STANDARD_S;
            ftype = VOID_F;
	    rfmt[1] = '\0';
	    def_matched = 0;

	    if( !(q = cnvt_code( p, &suppress, &size, &ftype )) )
		return( -ac );

	    def_value = cnvt_def( q + 1 );

	    if( suppress )
	    {
		/* Skip format string, assume no spaces in default string */

		for( q = p; *q && !isspace(*q); q++ )
		    ;
		
		start = nxtarg + 1;
		nxtarg = get_next_arg( "", ONLY_NONDELIMITED, start, &arg );

		continue; /* There are no assignments or defaults in arg list */
	    }
	    else
		p++;

	    strncat( rfmt, p, 1 + q - p );     /* Set conversion format to use */
	    type = assoc_type( ftype, size );  /* Find underlying read type    */
	    
	    switch( type )
	    {
	      case CHAR_T:

		if (def_value) {
		  if (*def_value=='%') 
		    defv_c= va_arg( ap, int ); /* Default passed as arg */
		  else sscanf( def_value, rfmt, &defv_c );
		}

		if (debug) {
		  if (def_value)
		    printf("#clproc:%s:option%d:%s:char:default=%c\n",
			   cl_argv[0],ac,name,defv_c);
		  else
		    printf("#clproc:%s:option%d:%s:char:default=(null)\n",
			   cl_argv[0],ac,name);
		}

		if( !arg && !def_value )
		    return( -ac );
		else if( !def_value )
		{
		    j = sscanf( arg, rfmt, va_arg( ap, char * ) );

		    if( j != 1 )
			return( -ac );
		}
		else
		{
		    strcpy( rfmtB, rfmt );
		    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
		    /* Check if matching argument found */

		    if( arg )
			j = sscanf( arg, rfmtB, (char *)cl_vbuf, (char *)cl_buf );

		    if( arg && j == 1 )
			sscanf( arg, rfmt, va_arg( ap, char * ) );    /* Use Command Line Value  */
		    else
		    {
			*(va_arg( ap, char * )) = defv_c;             /* Use Given Default Value */
			def_matched = 1;
		    }
		}
		
		break;
		
	      case SHORT_T:

		if (def_value) {
		  if (*def_value=='%') 
		    defv_s= va_arg( ap, int ); /* Default passed as arg */
		  else sscanf( def_value, rfmt, &defv_s );
		}

		if (debug) {
		  if (def_value)
		    printf("#clproc:%s:option%d:%s:short:default=%d\n",
			   cl_argv[0],ac,name,(int)defv_s);
		  else
		    printf("#clproc:%s:option%d:%s:short:default=(null)\n",
			   cl_argv[0],ac,name);
		}

		if( !arg && !def_value )
		    return( -ac );
		else if( !def_value )
		{
		    j = sscanf( arg, rfmt, va_arg( ap, short * ) );

		    if( j != 1 )
			return( -ac );
		}
		else
		{

		    strcpy( rfmtB, rfmt );
		    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
		    /* Check if matching argument found */

		    if( arg )
			j = sscanf( arg, rfmtB, (short *)cl_vbuf, (char *)cl_buf ); 

		    if( arg && j == 1 )
			sscanf( arg, rfmt, va_arg( ap, short * ) );    /* Use Command Line Value  */
		    else
		    {
			*(va_arg( ap, short * )) = defv_s;             /* Use Given Default Value */
			def_matched = 1;
		    }
		}
		
		break;
		
	      case INT_T:

		if (def_value) {
		  if (*def_value=='%') 
		    defv_i= va_arg( ap, int ); /* Default passed as argument */
		  else sscanf( def_value, rfmt, &defv_i );
		}

		if (debug) {
		  if (def_value)
		    printf("#clproc:%s:option%d:%s:int:default=%d\n",
			   cl_argv[0],ac,name,defv_i);
		  else
		    printf("#clproc:%s:option%d:%s:int:default=(null)\n",
			   cl_argv[0],ac,name);
		}

		if( !arg && !def_value )
		    return( -ac );
		else if( !def_value )
		{
		    j = sscanf( arg, rfmt, va_arg( ap, int * ) );

		    if( j != 1 )
			return( -ac );
		}
		else
		{
		    strcpy( rfmtB, rfmt );
		    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
		    /* Check if matching argument found */

		    if( arg )
			j = sscanf( arg, rfmtB, (int *)cl_vbuf, (char *)cl_buf ); 

		    if( arg && j == 1 )
			sscanf( arg, rfmt, va_arg( ap, int * ) );    /* Use Command Line Value  */
		    else
		    {
			*(va_arg( ap, int * )) = defv_i;             /* Use Given Default Value */
			def_matched = 1;
		    }
		}
		
		break;

	      case LONG_T:

		if (def_value) {
		  if (*def_value=='%') 
		    defv_l= va_arg( ap, long ); /* Default passed as arg */
		  else sscanf( def_value, rfmt, &defv_l );
		}

		if (debug) {
		  if (def_value)
		    printf("#clproc:%s:option%d:%s:long:default=%ld\n",
			   cl_argv[0],ac,name,defv_l);
		  else
		    printf("#clproc:%s:option%d:%s:long:default=(null)\n",
			   cl_argv[0],ac,name);
		}

		if( !arg && !def_value )
		    return( -ac );
		else if( !def_value )
		{
		    j = sscanf( arg, rfmt, va_arg( ap, long * ) );

		    if( j != 1 )
			return( -ac );
		}
		else
		{
		    strcpy( rfmtB, rfmt );
		    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
		    /* Check if matching argument found */

		    if( arg )
			j = sscanf( arg, rfmtB, (long *)cl_vbuf, (char *)cl_buf ); 

		    if( arg && j == 1 )
			sscanf( arg, rfmt, va_arg( ap, long * ) );    /* Use Command Line Value  */
		    else
		    {
			*(va_arg( ap, long * )) = defv_l;             /* Use Given Default Value */
			def_matched = 1;
		    }
		}
		
		break;

	      case FLOAT_T:

		if (def_value) {
		  if (*def_value=='%') 
		    defv_f= va_arg( ap, double ); /* Default passed as arg */
		  else sscanf( def_value, rfmt, &defv_f );
		}

		if (debug) {
		  if (def_value)
		    printf("#clproc:%s:option%d:%s:float:default=%f\n",
			   cl_argv[0],ac,name,defv_f);
		  else
		    printf("#clproc:%s:option%d:%s:float:default=(null)\n",
			   cl_argv[0],ac,name);
		}

		if( !arg && !def_value )
		    return( -ac );
		else if( !def_value )
		{
		    j = sscanf( arg, rfmt, va_arg( ap, float * ) );

		    if( j != 1 )
			return( -ac );
		}
		else
		{
		    strcpy( rfmtB, rfmt );
		    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
		    /* Check if matching argument found */

		    if( arg )
			j = sscanf( arg, rfmtB, (float *)cl_vbuf, (char *)cl_buf ); 

		    if( arg && j == 1 ) {
		      sscanf( arg, rfmt, va_arg( ap, float * ) );    /* Use Command Line Value  */
                    }
		    else
		    {
			*(va_arg( ap, float * )) = defv_f;             /* Use Given Default Value */
			def_matched = 1;
		    }
		}
		
		break;

	      case DOUBLE_T:

		if (def_value) {
		  if (*def_value=='%') 
		    defv_d= va_arg( ap, double ); /* Default passed as arg */
		  else sscanf( def_value, rfmt, &defv_d );
		}

		if (debug) {
		  if (def_value)
		    printf("#clproc:%s:option%d:%s:double:default=%f\n",
			   cl_argv[0],ac,name,defv_d);
		  else
		    printf("#clproc:%s:option%d:%s:double:default=(null)\n",
			   cl_argv[0],ac,name);
		}

		if( !arg && !def_value )
		    return( -ac );
		else if( !def_value )
		{
		    j = sscanf( arg, rfmt, va_arg( ap, double * ) );

		    if( j != 1 )
			return( -ac );
		}
		else
		{
		    strcpy( rfmtB, rfmt );
		    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
		    /* Check if matching argument found */

		    if( arg )
			j = sscanf( arg, rfmtB, (double *)cl_vbuf, (char *)cl_buf ); 

		    if( arg && j == 1 )
			sscanf( arg, rfmt, va_arg( ap, double * ) );    /* Use Command Line Value  */
		    else
		    {
			*(va_arg( ap, double * )) = defv_d;             /* Use Given Default Value */
			def_matched = 1;
		    }
		}
		
		break;

	      case CHARPTR_T:
		
		if (def_value) {
		  if (*def_value=='%') 
		    defv_p= 
		      strdup(va_arg( ap, char* )); /* Default passed as arg */
		  else {
		    defv_p = strdup( def_value );
		    for( i = 0; defv_p[i] && defv_p[i] != ']'; i++ );
		    defv_p[i] = '\0';
		  }
		}

		if (debug) {
		  if (def_value)
		    printf("#clproc:%s:option%d:%s:string:default=<%s>\n",
			   cl_argv[0],ac,name,defv_p);
		  else
		    printf("#clproc:%s:option%d:%s:string:default=<(null)>\n",
			   cl_argv[0],ac,name);
		}

		if( !arg && !def_value )
		    return( -ac );
		else if( !def_value )
		{
		    if (!strcmp(rfmt,"%s")) 
		      strcpy( va_arg( ap, char* ), arg );
		    else {
		      j = sscanf( arg, rfmt, va_arg( ap, char * ) );
		      if( j != 1 )
			return( -ac );
		    }
		}
		else
		{
		    /* Check if matching argument found */

		    if( arg )
			j = sscanf( arg, rfmt, (int *)cl_vbuf, (char *)cl_buf ); 
		    if( arg && j == 1 ) {
			strcpy( va_arg(ap, char*), arg );    /* Use Command Line Value  */
		    }
		    else
		    {
			strcpy( va_arg( ap, char * ), defv_p );       /* Use Given Default Value */
			def_matched = 1;
		    }

		    Free( defv_p );
		}

		break;

	      case UNSIGNED_CHAR_T:

		if (def_value) {
		  if (*def_value=='%') 
		    defv_uc= va_arg( ap, int ); /* Default passed as arg */
		  else sscanf( def_value, rfmt, &defv_uc );
		}

		if (debug) {
		  if (def_value)
		    printf("#clproc:%s:option%d:%s:uchar:default=%d\n",
			   cl_argv[0],ac,name,(int)defv_uc);
		  else
		    printf("#clproc:%s:option%d:%s:uchar:default=(null)\n",
			   cl_argv[0],ac,name);
		}

		if( !arg && !def_value )
		    return( -ac );
		else if( !def_value )
		{
		    j = sscanf( arg, rfmt, va_arg( ap, unsigned char * ) );

		    if( j != 1 )
			return( -ac );
		}
		else
		{
		    /* Check if matching argument found */

		    if( arg )
			j = sscanf( arg, rfmt, (unsigned char *)cl_vbuf, (char *)cl_buf ); 

		    if( arg && j == 1 )
			sscanf( arg, rfmt, va_arg( ap, unsigned char * ) );   /* Use Command Line Value  */
		    else
		    {
			*(va_arg( ap, unsigned char * )) = defv_uc;           /* Use Given Default Value */
			def_matched = 1;
		    }
		}
		
		break;
		
	      case UNSIGNED_SHORT_T:

		if (def_value) {
		  if (*def_value=='%') 
		    defv_us= va_arg( ap, int ); /* Default passed as arg */
		  else sscanf( def_value, rfmt, &defv_us );
		}

		if (debug) {
		  if (def_value)
		    printf("#clproc:%s:option%d:%s:ushort:default=%d\n",
			   cl_argv[0],ac,name,(int)defv_us);
		  else
		    printf("#clproc:%s:option%d:%s:ushort:default=(null)\n",
			   cl_argv[0],ac,name);
		}

		if( !arg && !def_value )
		    return( -ac );
		else if( !def_value )
		{
		    j = sscanf( arg, rfmt, va_arg( ap, unsigned short * ) );

		    if( j != 1 )
			return( -ac );
		}
		else
		{

		    strcpy( rfmtB, rfmt );
		    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
		    /* Check if matching argument found */

		    if( arg )
			j = sscanf( arg, rfmtB, (unsigned short *)cl_vbuf, (char *)cl_buf ); 

		    if( arg && j == 1 )
			sscanf( arg, rfmt, va_arg( ap, unsigned short * ) );   /* Use Command Line Value  */
		    else
		    {
			*(va_arg( ap, unsigned short * )) = defv_us;           /* Use Given Default Value */
			def_matched = 1;
		    }
		}
		
		break;
		
	      case UNSIGNED_INT_T:

		if (def_value) {
		  if (*def_value=='%') 
		    defv_ui= va_arg( ap, unsigned int ); /* Default passed as arg */
		  else sscanf( def_value, rfmt, &defv_ui );
		}

		if (debug) {
		  if (def_value)
		    printf("#clproc:%s:option%d:%s:uint:default=%d\n",
			   cl_argv[0],ac,name,defv_ui);
		  else
		    printf("#clproc:%s:option%d:%s:uint:default=(null)\n",
			   cl_argv[0],ac,name);
		}

		if( !arg && !def_value )
		    return( -ac );
		else if( !def_value )
		{
		    j = sscanf( arg, rfmt, va_arg( ap, unsigned int * ) );

		    if( j != 1 )
			return( -ac );
		}
		else
		{
		    strcpy( rfmtB, rfmt );
		    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
		    /* Check if matching argument found */

		    if( arg )
			j = sscanf( arg, rfmtB, (unsigned int *)cl_vbuf, (char *)cl_buf ); 

		    if( arg && j == 1 )
			sscanf( arg, rfmt, va_arg( ap, unsigned int * ) );  /* Use Command Line Value  */
		    else
		    {
			*(va_arg( ap, unsigned int * )) = defv_ui;          /* Use Given Default Value */
			def_matched = 1;
		    }
		}
		
		break;

	      case UNSIGNED_LONG_T:

		if (def_value) {
		  if (*def_value=='%') 
		    defv_ul= va_arg( ap, unsigned long ); /* Default passed as arg */
		  else sscanf( def_value, rfmt, &defv_ul );
		}

		if (debug) {
		  if (def_value)
		    printf("#clproc:%s:option%d:%s:ulong:default=%ld\n",
			   cl_argv[0],ac,name,defv_ul);
		  else
		    printf("#clproc:%s:option%d:%s:ulong:default=(null)\n",
			   cl_argv[0],ac,name);
		}

		if( !arg && !def_value )
		    return( -ac );
		else if( !def_value )
		{
		    j = sscanf( arg, rfmt, va_arg( ap, unsigned long * ) );

		    if( j != 1 )
			return( -ac );
		}
		else
		{
		    if( *def_value == '%' )          /* Default passed as an argument */
			defv_ul = va_arg( ap, unsigned long );
		    else
			sscanf( def_value, rfmt, &defv_ul );

		    strcpy( rfmtB, rfmt );
		    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
		    /* Check if matching argument found */

		    if( arg )
			j = sscanf( arg, rfmtB, (unsigned long *)cl_vbuf, (char *)cl_buf ); 

		    if( arg && j == 1 )
			sscanf( arg, rfmt, va_arg( ap, unsigned long * ) );  /* Use Command Line Value  */
		    else
		    {
			*(va_arg( ap, unsigned long * )) = defv_ul;          /* Use Given Default Value */
			def_matched = 1;
		    }
		}
		
		break;

	      case VOID_T:
	      default:

		return( -ac );
	    }

	    /* Only mark an argument as seen if it was in fact processed */

	    if( !def_matched )
	    {
		cl_argseen[nxtarg] = 1;
		start = cl_lastarg = nxtarg;
	    }

	    /* Skip rest of conversion format, assume no spaces in default string */

	    for( q++; *q && !isspace(*q); q++ )
		;

	    /* Get Next Argument */

	    start = nxtarg + 1;
	    nxtarg = get_next_arg( "", ONLY_NONDELIMITED, start, &arg );

	    ac++;    /* Defaults add to the count */
	}

	return( ac );
    }
    else
	return( 0 );

    va_end( ap );
}

void     cl_configure( char *name, char *value )
{
    /* Re-Scan and Classify Argument List in case something has changed */

    classify_args();
}


void     cl_mark( char *prefix, char *tag, char *showtime  )
{
    int      i;
    time_t   exectime;

    Message( "%s%s", prefix, tag );

    for( i = 0; i < cl_argc; i++ )
	Message( "%s ", cl_argv[i] );

    Message( "\n" );

    if( showtime && (time( &exectime ) != -1) )
	Message( "%s%s %s\n", prefix, showtime, ctime( &exectime ) );
}

int      cl_cmd_line( int *pargc, char ***pargv )
{
    int      i;

    if( cl_argc && cl_argv )
    {
	*pargc = cl_argc;
	*pargv = Calloc( cl_argc, char * );

	for( i = 0; i < cl_argc; i++ )
	    (*pargv)[i] = strdup( cl_argv[i] );

	return( 1 );
    }

    return( 0 );
}

void     cl_cleanup( void )
{
    int     i;

    if (debug) printf("#clproc:%s:cleanup\n",cl_argv[0]);
    if( cl_argc && cl_argv )
    {
	for( i = 0; i < cl_argc; i++ )
	    Free( cl_argv[i] );

	Free( cl_argv );

	cl_argc = 0;
	cl_argv = NULL;

	if( tagset != NULL )
	{
	    for( i = 0; i <= tagset_len; i++ )
		Free( tagset[i] );

	    Free( tagset );
	    Free( flag_match );
	    Free( local_flag_match );
	}
    }
}

/*** Private Functions ***/

static int     get_next_arg( char *names, int class, int start, char **ptr, ... )
{
    char      ch;
    int       i, j, k, namelen;
    int       supp, size, type;
    char      *p, *q, *fmt;
    int       (*cmp)( const char *, const char * );
    int       (*cmpn)( const char *, const char *, size_t );
    va_list   ap;

    va_start( ap, ptr );
    
    if( class == NONDELIMITED )
    {
	for( i = Max(start,1); i < cl_argc; i++ )
	{
	    if( cl_argseen[i] || !(cl_argtype[i] & NONDELIM) )
		continue;

	    *ptr = cl_argv[i];
	    return( i );
	}
    }
    else if( class == ONLY_NONDELIMITED )   /* Make sure arg at start is not seen or delimited */
    {
	if( start >= 0 && start < cl_argc && !cl_argseen[start] && (cl_argtype[start] & NONDELIM) )
	{
	    *ptr = cl_argv[start];
	    return( start );
	}
    }
    else if( class == MULTI )
    {
	for( k = 0; cl_delimind[k] >= 0; k++ )
	{
	    i = cl_delimind[k];
	    
	    if( cl_argseen[i] || i < start )
		continue;

	    for( p = names - 1; p; p = strchr( p, '|' ) )
	    {
		ch = *++p;            /* Skip delimiter or, in first iteration, move to beginning */

		if( case_fold )
		    ch = tolower(ch);

		if( *ptr = strchr( cl_delimarg[k], ch ) )  /* Find match */
		    break;
	    }

	    if( p )
	    {
		*(va_arg(ap, int *)) = cl_delimtag[k];    /* Record tag set element for later use */
		*(va_arg(ap, char **)) = cl_delimarg[k];  /* Record flag string for later use     */

		return( i );
	    }
	}
    }
    else if( class == FLAG )
    {
	if( case_fold )
	    cmp = strcasecmp;
	else
	    cmp = strcmp;

	for( k = 0; cl_delimind[k] >= 0; k++ )
	{
	    i = cl_delimind[k];
	    
	    if( cl_argseen[i] || i < start )
		continue;

	    /* If so, check part after tag against each name */
            /* If there is a match, return it.               */

	    for( p = names; p;  )
	    {
		q = strchr( p, '|' );
		namelen = ( q ) ? q - p : strlen( p );

		strncpy( cl_cbuf, p, namelen );
		cl_cbuf[namelen] = '\0';
		
		if( !((*cmp)( cl_cbuf, cl_delimarg[k] )) )
		    break;
		
		p = ( q ) ? q + 1 : NULL;
	    }

	    if( p )
	    {
		*(va_arg(ap, int *)) = cl_delimtag[k];    /* Record tag set element for later use */

		*ptr = cl_delimarg[k];
		return( i );
	    }
	}
    }
    else  /* option -- does not handle no separator well */
    {
	if( case_fold )
	{
	    cmp = strcasecmp;
	    cmpn = strncasecmp;
	}
	else
	{
	    cmp = strcmp;
	    cmpn = strncmp;
	}

	fmt = va_arg( ap, char * );

	for( k = 0; cl_delimind[k] >= 0; k++ )
	{
	    i = cl_delimind[k];
	    
	    if( cl_argseen[i] || i < start )
		continue;

	    /* If so, check part after tag against each name */
            /* If there is a match, return it.               */

	    for( p = names; p;  )
	    {
		q = strchr( p, '|' );
		namelen = ( q ) ? q - p : strlen( p );

		strncpy( cl_cbuf, p, namelen );
		cl_cbuf[namelen] = '\0';
		
		if( !((*cmp)( cl_cbuf, cl_delimarg[k] )) )    /* Exact Match */
		    break;
		else if( fmt != NULL && !((*cmpn)( cl_cbuf, cl_delimarg[k], namelen )) )
		{
		    /* Currently only checking first conversion string */
		    /* Eventually deal with the possibility that this  */
		    /* one is missing and check the next one, etc.     */

		    /* If conversion is %s, then don't accept.         */

		    j = cnvt_code( fmt, &supp, &size, &type ) - fmt;

		    if( fmt[j] != 's' )
		    {
			strcat( cl_cbuf, "%" );

			if( size != STANDARD_S )
			    sprintf( cl_cbuf + namelen + 1, "%c%c", fmt[j-1], fmt[j] );
			else
			    sprintf( cl_cbuf + namelen + 1, "%c", fmt[j] );

			if( sscanf( cl_delimarg[k], cl_cbuf, (char *)cl_vbuf ) == 1 )
			    break;
		    }
		}
		
		p = ( q ) ? q + 1 : NULL;
	    }

	    if( p )
	    {
		*(va_arg(ap, char **)) = p;      /* Record flag string for later use     */

		*ptr = cl_delimarg[k];
		return( i );
	    }
	}
    }

    va_end( ap );

    *ptr = NULL;
    return( -1 );
}

static char   *cnvt_code( char *p, int *suppress, int *size, int *type )
{
    if( p[0] != '%' )
	return( NULL );

    p++;
    
    *suppress = 0;
    *size = STANDARD_S;
    *type = VOID_F;
    
    if( *p == '*' )
    {
	*suppress = 1;
	p++;
    }

    if( *p == 'l' )
    {
	*size = LONG_S;
	p++;
    }
    else if( *p == 'h' )
    {
	*size = SHORT_S;
	p++;
    }
    else if( *p == 'L' )
    {
	*size = LONG_DOUBLE_S;
	p++;
    }

    switch( *p )
    {
      case 'd':
      case 'i':

	*type = INTEGER_F;
	break;
		
      case 'u':
      case 'o':
      case 'x':
      case 'X':

	*type = UNSIGNED_F;
	break;
		
      case 'e':
      case 'E':
      case 'f':
      case 'g':
      case 'G':

	*type = REAL_F;
	break;

      case 'c':

	*type = CHAR_F;
	break;
	
      case 's':

	*type = STRING_F;
	break;

      default:

	return( NULL );
    }

    return( p );
}

static char     *cnvt_def( char *s )
{
    if( *s != '[' )
	return( NULL );
    else
    {
	do
	{
	    s++;
	    if( *s == '\0' || *s == ']' )
		return( NULL );
	}
	while( isspace(*s) );

	return( s );
    }
}

static int      assoc_type( int ftype, int size )
{
    switch( ftype )
    {
      case INTEGER_F:

	if( size == STANDARD_S )
	    return( INT_T );
	else if( size == LONG_S )
	    return( LONG_T );
	else if( size == SHORT_S )
	    return( SHORT_T );
	else
	    return( VOID_T );
	

      case UNSIGNED_F:

	if( size == STANDARD_S )
	    return( UNSIGNED_INT_T );
	else if( size == LONG_S )
	    return( UNSIGNED_LONG_T );
	else if( size == SHORT_S )
	    return( UNSIGNED_SHORT_T );
	else if( size == CHAR_S )
	    return( UNSIGNED_CHAR_T );
	else
	    return( VOID_T );
	
	
      case REAL_F:
	
	if( size == STANDARD_S )
	    return( FLOAT_T );
	else if( size == LONG_S )
	    return( DOUBLE_T );
	else
	    return( VOID_T );
	

      case CHAR_F:

	if( size == STANDARD_S )
	    return( CHAR_T );
	else
	    return( VOID_T );


      case STRING_F:

	if( size == STANDARD_S )
	    return( CHARPTR_T );
	else
	    return( VOID_T );


      default:
	return( VOID_T );
    }
}

static int isThisTokenANumber(const char* s)
{
  const char* runner= s;
  while (*runner && !isspace(*runner)) {
    if (!isdigit(*runner) && !strchr("-.",*runner)) return 0;
    runner++;
  }
  return 1;
}

static void    classify_args( void )
{
    int       i, j, k, len, namelen;
    int       delim_ind = 0, nondelim_ind = 0;
    char      *p, *q;
    
    for( i = 1; i < cl_argc; i++ )   /* Do not include command name */
    {
	/* Check if it matches an element of tag set */

	for( j = 0; j < tagset_len; j++ )
	{
	    len = strlen( tagset[j] );
	    
	    if( !strncmp( cl_argv[i], tagset[j], len) )  /* Tag matches always case sensitive */
		break;
	}

	/* If so, add to delimited list, record which tag */
	/* Otherwise, add to non-delimited list           */
	/* Note that when the tagset contains "-", then   */
	/* negative numbers are put in both lists, to be  */
	/* taken off first come first served.             */

  	    /* Don't be fooled by negative #'s or by a string "-" */

	if( j == tagset_len || !strcmp( cl_argv[i], "-" ) || isThisTokenANumber(cl_argv[i]) ) 
	{
	    cl_nondelimind[nondelim_ind] = i;
	    cl_nondelimarg[nondelim_ind] = cl_argv[i];
	    cl_argtype[i] |= NONDELIM;

	    nondelim_ind++;
	}
	else  /*  if( j < tagset_len ) */
	{
	    cl_delimtag[delim_ind] = j;
	    cl_delimind[delim_ind] = i;
	    cl_delimarg[delim_ind] = cl_argv[i] + len;
	    cl_argtype[i] |= DELIM;
	    delim_ind++;
	}
    }

    /* Set Sentinels, there is space because the command name is not used */

    cl_delimind[delim_ind] = -1;
    cl_nondelimind[nondelim_ind] = -1;
}
