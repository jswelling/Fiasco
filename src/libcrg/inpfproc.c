#include   <stdio.h>
#include   <stdlib.h>
#include   <stdarg.h>
#include <string.h>
#include   <strings.h>
#include   <ctype.h>

#include   <errors.h>
#include   "stdcrg.h"

/************************************************************
 *                                                          *
 *  inpfproc.c                                              *
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
 *  Last Modified:    May 1996                              *
 *  Modified By:      Chris Genovese                        *
 *                                                          *
 ************************************************************/

/*** Data Structures  ***/

struct hashnode
{
    unsigned int      level;               /* Level of input file                           */
    
    char             *tag;                 /* Tag string                                    */
    char             *value;               /* Value string separated by \0 with \0\0 at end */
    char             *info;                /* File and line information for tag             */

    unsigned int      count;               /* How many fields are represented in value      */
    unsigned int      where;               /* How many values have been read and matched    */
    unsigned int      matched;             /* How many times node has been fully matched    */
    unsigned short    set;                 /* How many times node has been set              */

    void             (*command)(char *);

    struct hashnode  *next;
};

typedef struct hashnode HASHNODE;

/*** Macros and Constants ***/

#define          REALLOC_EXPAND      2
#define          TABLE_DEFAULT      32

enum types
{
    VOID_T, CHAR_T, SHORT_T, INT_T, LONG_T, FLOAT_T, DOUBLE_T, CHARPTR_T, CHARPTRPTR_T,
    UNSIGNED_CHAR_T, UNSIGNED_SHORT_T, UNSIGNED_INT_T, UNSIGNED_LONG_T
};

enum sizes     { STANDARD_S, CHAR_S, SHORT_S, LONG_S, LONG_DOUBLE_S };
enum ftypes    { VOID_F, CHAR_F, INTEGER_F, UNSIGNED_F, REAL_F, STRING_F, VECTOR_F, MATRIX_F };
enum retval    { FAILURE, SUCCESS };
enum cmd_names { APPEND, ARG, COND, DATA, END, ENDINPUT, INCLUDE, REQUIRE, NUM_COMMANDS };
enum cross_eof { END_AT_EOF, ACROSS_EOF, ONE_LINE };

enum special_chars
{
    comment_char  = '#',       /* Comment Character                 */
    command_open  = '<',       /* Opening delimiter for commands    */
    command_close = '>',       /* Closing delimiter for commands    */
    quote_char    = '\"',      /* Character for protecting strings  */
    var_ref_char  = '$',       /* Reference tags in data $<tagname> */
    escape_char   = '\\'       /* Escape Sequence Lead-in           */
};

/*** Function Prototypes ***/

    /* Public Functions  */

int      inpf_scan( char *fn );                    /* Read Input File and Prepare   */
int      inpf_scanx( FILE *, int, char * );        /* Expert version                */

int      inpf_get( char *tag, char *fmt, ... );    /* Get Argument or Flag Values   */
void     inpf_arg( char *tag, char *label );       /* Override with CL argument     */

char    *inpf_data( char *tag );                   /* Data string for tag           */
int      inpf_unmatched( char *tag );              /* Has tag been matched?         */
char    *inpf_taginfo( char * );                   /* File and Line info for Tag    */

char    *inpf_failures( void );                    /* Any unmatched tags?           */
int      inpf_errors( void );                      /* Error Reporting (check)       */
int      inpf_gerror( int *, char *** );           /* Error Reporting (store)       */
int      inpf_perror( void );                      /* Error Reporting (print)       */

void     inpf_config( char *opt, char *val );      /* Set default interpretations   */
void     inpf_cleanup( void );                     /* Clean-up Internal Info        */

FILE    *inpf_insert( char *, FILE * );            /* Copy input file               */


    /* Private Functions */

static void      inpf_initialize( void );            /* Initialize Internal Objects   */
static void      set_error( int, char *, ... );      /* Register processing error     */
static int       check_error( void );                /* Check error status, no change */

static char   	*get_tag( char ** );                 /* Input parsing                 */
static int       get_dat( char **, char **, char **, int * );

static void   	 do_append( char * );                /* Built-in commands             */
static void   	 do_arg( char * );
static void      do_cond( char * );
static void   	 do_data( char * );
static void   	 do_include( char * );
static void   	 do_require( char * );

static HASHNODE *lookup( char * );                   /* Hash Table Support            */
static HASHNODE *install( char *, char * );
static int       hash( char * );

                                                     /* General Parsing               */
static char     *gobble_ws( char *s );
static char     *gobble_line( char *s );
static char     *skip_quoted( char *s, char *stopset, char *stop );
static char     *skip_non_ws( char *s, char *stopset, char *stop );

                                                     /* Conversion Code Parsing       */
static char     *cnvt_code( char *s, int *suppress, int *size, int *type );
static char     *cnvt_def( char *s, int ftype, char **q );
static int       parse_format( int, char *, char **, char *, HASHNODE *, va_list * );
static int       parse_vector( int, char *, char **, char *, HASHNODE *, va_list * );
static int       parse_matrix( int, char *, char **, char *, HASHNODE *, va_list * );
static int       single_format( int, char *, char **, char *, HASHNODE *, void * );
static int       assoc_type( int ftype, int size );
static size_t    TypeSize( int type );


#ifdef USE_INLINE
#  define INLINE inline
#else
#  define INLINE
#endif

static INLINE char*    buf_add_string( char* bp, char* sp, size_t len, char** bufp, int* buflenp );
static INLINE char*    buf_add_char( char* bp, char c, char** bufp, int* buflenp );
static INLINE char*    buf_add_char( char* bp, char c, char** bufp, int* buflenp );
static INLINE char*    buf_inc_pntr( char* bp, char** bufp, int* buflenp );
static INLINE char*    buf_check( size_t inc, char* bp, char** bufp, int* buflenp );


static char     *get_line( int );
static int       get_token( char **psp, char **pbp, char **bufp, int *buflenp, int code );
static void      scan_to_end( char **, char **, int * );


/*** Global Variables ***/

    /* Working Buffers        */

static   int     line_bufsize = 1024;          /* Input Lines     */
static   char   *line_buf = NULL;

static   int     err_bufsize = 768;            /* Error Messages  */
static   char   *err_buf = NULL;

static   int     tag_bufsize = 256;            /* Tag Names       */
static   char   *tag_buf = NULL;

static   int     data_bufsize = 4096;          /* Data Strings    */
static   char   *data_buf = NULL;

static   int     def_bufsize = 1024;           /* Default Values  */
static   void   *def_buf = NULL;    
static   void   *def_vbuf = NULL;


    /* Hash Table and Support */

static int        inpf_hash_length = 199; 
static HASHNODE **inpf_hash_table = NULL;


    /* Error Handling Support */

static   int     errc = 0;
static   char  **errv = NULL;
static   int     errlen = 16;

static   char   *LastFailure = NULL;

    /* File Pointer Stack     */

static int        fp_stackind = -1;
static int        fp_stacklen = 16;
static FILE     **fp_stack = NULL;
static int       *Line_num = NULL;
static char     **Names = NULL;

    /* Command Line Support   */

static int        cmd_line_set = 0;
static int        inpf_argc = 0;
static char     **inpf_argv = NULL;

    /* Option names and default values */

static char *inpf_options[] = {
                                "case-fold",
				"tag-sep",
                                "warnings"
                              };

static   int     case_fold = 0;
static   char   *tag_sep = "";   /* White space always ok, if tag_sep != "" allow it surrounded by ws */
static   int     warnings = 0;


    /* Built-in Command Names          */

static char *inpf_commands[] = {
                                 "<append\0",
                                 "<arg\0",
				 "<cond\0",
                                 "<data\0",
                                 "<end\0",
                                 "<endinput\0",
                                 "<include\0",
				 "<require\0"
                               };


/*
 *  general comments here
 *
 */

int      inpf_scan( char *filename )
{
    int        status;
    FILE      *fp;

    /* Prepare file for reading */

    fp = fopen( filename, "r" );

    if( !fp )
    {
	inpf_initialize();
	set_error( fp_stacklen, "Cannot open file %s", filename );
	return( FAILURE );
    }

    status = inpf_scanx( fp, 0, filename );

    efclose( fp );
    return( status );
}

int      inpf_scanx( FILE *fp, int lines_so_far, char *fname )
{
    int      j, status = SUCCESS;

    char     *sp, *tp, *tag;

    HASHNODE  *hp;

    /* Initialization           */

    inpf_initialize();

    /* Prepare file for reading */

    fp_stackind++;
    fp_stack[fp_stackind] = fp;
    Line_num[fp_stackind] = lines_so_far;
    Names[fp_stackind]    = strdup( fname );

    /* Read file one line at a time */

    for( ;; )
    {
	if( !get_line( ACROSS_EOF ) )
	    break;
	
	/* Parse first token, behavior depends on first non-ws character */

	sp = gobble_ws( line_buf );                     
	
	if( *sp == '\0' || *sp == comment_char )        
	    continue;                              /* Skip until end of line  */
	else if( *sp == command_open )
	{
	    tp = sp;                               /* Search for comand close */
    
	    while( *tp && (*tp != '\n') && (*tp != command_close) && (*tp != comment_char) )
		tp++;

	    if( *tp != command_close )
	    {
		if( *tp == comment_char )
		    set_error( fp_stackind, "Illegal comment character in command designator" );
		else
		    set_error( fp_stackind, "Runaway command designator, no closing >" );

		status = FAILURE;
		continue;
	    }

	    *tp++ = '\0';
	    
	    if( !strcmp( sp, inpf_commands[ENDINPUT] ) )
	    {
		if( fp_stackind > 0 )
		{
		    fclose( fp_stack[fp_stackind] );
		    fp_stack[fp_stackind--] = NULL;
		}
		else
		    break;
	    }
	    else if( !strcmp( sp, inpf_commands[END] ) )
	    {
		set_error( fp_stackind, "Encountered %c%s%c command not matched by %c%s%c",
			   command_open, inpf_commands[END], command_close, command_open,
			   inpf_commands[DATA], command_close );

		status = FAILURE;
		continue;
	    }
	    else
	    {
		hp = lookup( sp );         /* Verify that command is correct */

		if( hp )
		{
		    hp->command( tp );     /* Run Command */

		    if( check_error() )
		    {
			status = FAILURE;
			continue;
		    }
		}
		else
		{
		    set_error( fp_stackind, "Unrecognized command %s", sp );
		    status = FAILURE;
		    continue;
		}
	    }
	}
	else if( *sp == command_close )
	{
	    set_error( fp_stackind, "Line should not begin with %c character", command_close );
	    status = FAILURE;
	    continue;
	}
	else
	{
	    /* Parse tag name (e.g., strip surrounding quotes) */

	    tag = get_tag( &sp );

	    if( check_error() )
	    {
		status = FAILURE;
		continue;
	    }

	    /* Deal with separator */

	    sp = gobble_ws( sp );

	    if( tag_sep[0] != '\0' )      /* Look for optional separator */
	    {
		j = strlen(tag_sep);
		
                if( !strncmp( tag_sep, sp, j ) )
		{
		    sp += j;
		    sp = gobble_ws( sp );
		}
	    }
	    
	    /* Parse (e.g., strip quotes, comments) and read data values */

	    tp = data_buf;
	    
	    get_dat( &sp, &tp, &data_buf, &data_bufsize );

	    if( check_error() )
	    {
		status = FAILURE;
		continue;
	    }

	    install( tag, data_buf );
	}
    }

    fp_stackind = -1;       /* Empty Stack */
    return( status );
}


/* Find command line over-ride for tag      */
/* No fancy processing, but optional prefix */

void     inpf_arg( char *tag, char *label )
{
    int        i, len, cnt;

    char      *p, *dat;
    HASHNODE  *hp;
    

    if( cmd_line_set )
    {
	p = (label) ? label : tag;
	len = strlen( p );

	for( i = 0; i < inpf_argc; i++ )
	{
	    if( !strncmp( p, inpf_argv[i], len ) && inpf_argv[i][len] == '=' )
	    {
		p = inpf_argv[i] + len + 1;       /* Throw away p, now safe to use data_buf */
		hp = install( tag, NULL );

		dat = data_buf;

		cnt = get_dat( &p, &dat, &data_buf, &data_bufsize );

		hp->value = Calloc( dat - data_buf + 1, char );
		memcpy( hp->value, data_buf, dat - data_buf );
		hp->count = cnt;

		break;
	    }
	}
    }
}


int      inpf_get( char *tag, char *fmt, ... )
{
    int       i;
    int       suppress, size, ftype, type;
    int       status = SUCCESS;

    char      rfmtA[32];
    char      *p, *q, *val;
    char      *def_value;
    char      *fmtcpy;

    double    hold[8];
    
    HASHNODE  *hp;

    va_list   ap;

    /* Search for given Tag */
    
    hp = lookup( tag );

    if( !hp )
    {
	hp = install( tag, "\0" );
	hp->set = 0;
    }
    
    if( fmt == NULL || fmt[0] == '\0' )    /* No format, just search and mark argument as matched */
    {
	hp->matched++;
	return( SUCCESS );
    }


    /* Process format string using stored value */

    va_start( ap, fmt );

    rfmtA[0] = '%';         /* Basic Format String */
    rfmtA[1] = '\0';

    fmtcpy = strdup( fmt );
    p = strchr(fmtcpy, '%');   /* Skip initial chaff  */


    /* Start after last datum previously matched */

    for( val = hp->value, i = 0; i < hp->where; i++ )
    {
	if( val[0] == '\0' )
	    break;

	val += strlen(val) + 1;
    }

    if( i < hp->where )
	hp->where = i;
    
    /* Scan format string storing requested values or their defaults */

    for( ; p; p = strchr(q, '%') )
    {
	suppress = '\0';         /* Reset for next argument */
	size = STANDARD_S;
	ftype = VOID_F;
	rfmtA[1] = '\0';

	if( !(q = cnvt_code( p, &suppress, &size, &ftype )) )
	{
	    set_error( fp_stackind, "In tag %s, failed to match data field %d", tag, hp->where + 1 );
	    status = FAILURE;
	    break;
	}

	def_value = cnvt_def( q + 1, ftype, &q ); 
	p++;               

	/*
	 * Where are we:                                                   
	 *   p now points to the null-padded format string without '%'     
	 *   def_value if non-null points to the null padded default string
	 *   q now points to the character beyond the format string        
	 *
	 * Next, we parse the format and assign the values unless suppress
	 * has been indicated by a '*' in the format string after '%'
	 *
	 */

	if( suppress )
	{
	    p++;             /* Move past the suppress marker */

	    if( ftype == VECTOR_F || ftype == MATRIX_F )
	    {
		set_error( fp_stackind, "Cannot suppress matrix or vector format types" );
		status = FAILURE;
		break;
	    }
		
	    type = assoc_type( ftype, size );   /* Find underlying read type    */
	    strcat( rfmtA, p );                 /* Set conversion format to use */

	    status = single_format( type, def_value, &val, rfmtA, hp, (void *)hold );

	    continue;
	}
	else if( ftype == VECTOR_F )
	    status = parse_vector( size, def_value, &val, rfmtA, hp, &ap );
	else if( ftype == MATRIX_F )
	    status = parse_matrix( size, def_value, &val, rfmtA, hp, &ap );
	else
	{
	    type = assoc_type( ftype, size );   /* Find underlying read type    */
	    strcat( rfmtA, p );                 /* Set conversion format to use */

	    status = parse_format( type, def_value, &val, rfmtA, hp, &ap );
	}
	
	if( status == FAILURE )
	    break;
    }

    va_end( ap );

    if( status == SUCCESS )
	hp->matched++;
    else
	LastFailure = hp->tag;

    Free( fmtcpy );

    return( status );
}

char    *inpf_data( char *tag )
{
    HASHNODE       *hp;

    hp = lookup( tag );

    if( hp )
	return( hp->value );
    else
	return( NULL );
}

char    *inpf_taginfo( char *tag )
{
    HASHNODE       *hp;

    hp = lookup( tag );

    if( hp )
	return( hp->info );
    else
	return( NULL );
}

int      inpf_unmatched( char *tag )
{
    HASHNODE       *hp;

    hp = lookup( tag );

    if( !hp )
	return( -1 );             /* Not Set             */
    else if( !(hp->matched) )
	return( 1 );              /* Set but not matched */
    else
	return( 0 );              /* Matched             */
}


char    *inpf_failures( void )
{
    return( LastFailure );
}

int      inpf_errors( void )
{
    return( errc );
}

int      inpf_gerror( int *ecp, char ***evp )
{
    int     i;

    if( ecp )
    {
	*ecp   = errc;
	(*evp) = Calloc( errc, char * );

	for( i = 0; i < errc; i++ )
	    (*evp)[i] = strdup( errv[i] );
    }
    
    return( errc );
}

int      inpf_perror( void )
{
    int     i;

    for( i = 0; i < errc; i++ )
	Message( "%s\n", errv[i] );
    
    return( errc );
}


void     inpf_configure( char *name, char *value )
{
}

static int initialized = 0;

void     inpf_cleanup( void )
{
    int        i;
    HASHNODE  *node, *hold;
    
    if( !initialized )
	return;
    
    /* Clean Up Buffers */

    Free( line_buf );
    Free( err_buf );
    Free( tag_buf );
    Free( data_buf );

    Free( def_buf );
    Free( def_vbuf );

    line_buf = NULL;
    err_buf  = NULL;
    tag_buf  = NULL;
    data_buf = NULL;

    def_buf  = NULL;
    def_vbuf = NULL;

    /* Clean Up Hash Table */

    for( i = 0; i < inpf_hash_length; i++ )
    {
	node = inpf_hash_table[i];

	while( node )
	{
	    hold = node;
	    node = node->next;

	    if( hold->tag )
		Free( hold->tag );

	    if( hold->value )
		Free( hold->value );

	    Free( hold );
	}
    }

    Free( inpf_hash_table );
    inpf_hash_table = NULL;

    /* Clean Up File Pointer Stack */

    for( i = 0; i < fp_stacklen; i++ )
    {
	if( Names[i] )
	    Free( Names[i] );
    }

    Free( fp_stack );
    Free( Line_num );
    Free( Names );

    fp_stack = NULL;
    Line_num = NULL;
    Names = NULL;
    fp_stackind = -1;

    /* Clean Up Error Structure    */

    for( i = 0; i < errlen; i++ )
	if( errv[i] )
	    Free( errv[i] );

    Free( errv );
    errv = NULL;
    errc = 0;

    LastFailure = NULL;


    /* Clean up command line copy  */

    if( cmd_line_set )
    {
	for( i = 0; i < inpf_argc; i++ )
	    if( inpf_argv[i] )
		Free( inpf_argv[i] );

	Free( inpf_argv );
	inpf_argv = NULL;
	inpf_argc = 0;
    }


    /* Mark as clean               */

    initialized = 0;
}

static void     inpf_initialize( void )
{
    int           i;
    char          num[16];
    char          *scr;
    HASHNODE      *hp;
    
    if( initialized )
	return;

    /* Initialize Buffers */

    line_buf = Calloc( line_bufsize, char );
    err_buf  = Calloc( err_bufsize,  char );
    tag_buf  = Calloc( tag_bufsize,  char );
    data_buf = Calloc( data_bufsize, char );

    def_buf  = (void *)Malloc( def_bufsize, char );
    def_vbuf = (void *)Malloc( def_bufsize, char );

    /* Get private copy of command line if scanned */

    cmd_line_set = cl_cmd_line( &inpf_argc, &inpf_argv );
    
    /* Initialize Hash Table */

    inpf_hash_table = Calloc( inpf_hash_length, HASHNODE * );

    hp = install( inpf_commands[APPEND], NULL );
    hp->command = do_append;
    
    hp = install( inpf_commands[ARG], NULL );
    hp->command = do_arg;
    
    hp = install( inpf_commands[COND], NULL );
    hp->command = do_cond;
    
    hp = install( inpf_commands[DATA], NULL );
    hp->command = do_data;
    
    hp = install( inpf_commands[INCLUDE], NULL );
    hp->command = do_include;
    
    hp = install( inpf_commands[REQUIRE], NULL );
    hp->command = do_require;
    
    if( cmd_line_set )
    {
	/* Insert command line arguments as special tags */

 	    /* NOTE: argv[0] gets stripped of leading path!  */
	    /*       Since this is usually what is desired.  */
	    /*       "00" is the full argv[0] string.        */

	scr = strrchr( inpf_argv[0], '/' );

	if( !scr )
	    scr = inpf_argv[0];
	else
	    scr++;

	sprintf( num,      "%d", 0 );
	sprintf( data_buf, "%s%c%c", scr, '\0', '\0' );
	    
	install( num, data_buf );

	sprintf( num,      "%d%d", 0, 0 );
	sprintf( data_buf, "%s%c%c", inpf_argv[0], '\0', '\0' );
	    
	install( num, data_buf );


	    /* Now, do the rest */

	for( i = 1; i < inpf_argc; i++ )
	{
	    sprintf( num,      "%d", i );
	    sprintf( data_buf, "%s%c%c", inpf_argv[i], '\0', '\0' );
	    
	    install( num, data_buf );
	}
    }
    
    /* Initialize File Pointer Stack */

    fp_stackind = -1;
    fp_stack = Calloc( fp_stacklen, FILE * );
    Line_num = Calloc( fp_stacklen, int );
    Names    = Calloc( fp_stacklen, char * );

    /* Initialize Error Structure */

    errc = 0;
    errv = Calloc( errlen, char * );
    LastFailure = NULL;

    initialized = 1;
}

/*
 * Copy an input file to an output stopping at EOF or <endinput>.
 * Currently, this is just a hack to get it working; can be more
 * careful.
 *
 */

FILE   *inpf_insert( char *name, FILE *ofp )
{
    int      c, i;
    FILE    *fp;
    char     check[32];
    
    if( !ofp || ferror(ofp) )
	return( NULL );

    fp = efopen( name, "r" );
    
    while( (c = getc(fp)) != EOF )
    {
	if( c == '<' )
	{
	    for( i = 0; c != EOF; i++ )
	    {
		check[i] = c;
		check[i+1] = '\0';

		if( i == 9 )
		    break;

		c = getc(fp);
	    }

	    if( i == 9 && !strcmp( check, "<endinput>" ) )
		return( fp );

	    fprintf( ofp, "%s", check );      /* Add it to stream */
	}
	else
	    putc( c, ofp );
    }

    return( fp );
}


/*** Private Functions ***/

static void    set_error( int stackind, char *fmt, ... )
{
    int            len;
    va_list        ap;

    va_start( ap, fmt );
    
    if( stackind < 0 )                   /* No line info         */
	sprintf( err_buf, "Input Error, File %s: ", Names[-stackind - 1] );
    else if( stackind < fp_stacklen )    /* File and line info   */
	sprintf( err_buf, "Input Error, File %s, Line %d: ", Names[stackind], Line_num[stackind] );
    else                                 /* No file or line info */
	sprintf( err_buf, "Input Error: " );

    len = strlen( err_buf );
    vsprintf( err_buf + len, fmt, ap );

    if( errc + 1 < errlen )
	errv[errc++] = strdup( err_buf );
    else
    {
	errlen *= REALLOC_EXPAND;
	errv = Realloc( errv, errlen, char * );
	errv[errc++] = strdup( err_buf );
    }

    va_end( ap );
}

static int     check_error( void )
{
    if( errc )
	return( 1 );
    else
	return( 0 );
}


static char   *get_tag( char **psp )
{
    char      stop;
    char      *tag, *sp;
    char      stopset[5];

    int       len;

    sp = *psp;

    sp = gobble_ws( sp );     /* Skip Leading Spaces */

    if( *sp == quote_char )
    {
	tag = sp + 1;
    
	sp = skip_quoted(sp, "\n", &stop);

	if( stop != quote_char )
	    set_error( fp_stackind, "Runaway tag --- no closing quotation mark" );
    }
    else
    {
	tag = sp;

	stopset[0] = comment_char;
	stopset[1] = command_open;
	stopset[2] = command_close;
	stopset[3] = tag_sep[0];
	stopset[4] = '\0';

	for( ;; )
	{
	    sp = skip_non_ws(sp, stopset, &stop);

	    switch( stop )
	    {
	      case command_open:

		set_error( fp_stackind, "Illegal command open character in tag" );
		break;
		
	      case command_close:

		set_error( fp_stackind, "Illegal command close character in tag" );
		break;
		
	      case comment_char:

		*sp = '\0';       /* End tag at comment  */

		while( *sp  )     /* Skip to end of line */
		    sp++;
		
		break;
	    }

	    if( isspace( stop ) || stop == '\0' )
		break;

	    if( tag_sep[0] != '\0' && stop != tag_sep[0] )
		break;
	    else if( tag_sep[0] != '\0' )  /* Skip Separator */
	    {
		/* This small shortcut opens up a bug.  For instance, if tag_sep = "=", */
		/* then a line of the form "tag = =" will lose the second =.  This is   */
		/* rather small, so I'll leave it for now.  The solution is to have a   */
		/* pushback character but I don't feel like doing that now.  NOTE!      */

		len = strlen( tag_sep );

		if( !strncmp( tag_sep, sp, len ) )
		{
		    *sp = '\0';
		    sp += len;
		    break;
		}
	    }
	}
    }

    if( *sp != '\0' )
    {
	*psp = sp + 1;
	*sp = '\0';
    }
    else
	*psp = sp;

    /* Save tag name in case of subsequent refreshes of input buffer */
    /* Ignore any characters that exceed tag_bufsize                 */

    strncpy( tag_buf, tag, tag_bufsize - 1 );
    
    return( tag_buf );
}

static int     get_dat( char **psp, char **pdp, char **datap, int *datalenp )
{
    int        tc, len, cnt;
    char       *sp, *dp;

    sp = gobble_ws( *psp );   /* Skip Leading Spaces        */
    dp = *pdp;

    for( len = cnt = 0; ; cnt++ )
    {
	tc = get_token( &sp, &dp, datap, datalenp, ONE_LINE );
	len += tc;

	if( !tc )
	    break;
    }

    *dp++ = '\0';     /* Cap it off with extra null */

    *psp = sp;        /* Adjust input pointers      */
    *pdp = dp;

    return( cnt );
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
	*suppress = '*';
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

      case 'V':

	*type = VECTOR_F;
	break;

      case 'M':

	*type = MATRIX_F;
	break;

      default:

	return( NULL );
    }

    return( p );
}

static char     *cnvt_def( char *s, int ftype, char **pp )
{
    int          i, j, plevel;

    *pp = s;
    
    if( *s != '[' )
	return( NULL );
    else
    {
	*s++ = '\0';    /* Pad format string */

	if( ftype != STRING_F )  /* Skip leading white space */
	{
	    while( isspace(*s) && *s != ']' && *s )
		s++;

	    if( *s == '\0' )
	    {
		*pp = s;
		set_error( fp_stackind, "Runaway default specifier in format string" );
		return( NULL );
	    }
	    else if( *s == ']' )
	    {
		*s = '\0';
		*pp = s + 1;
		return( NULL );
	    }
	}

	/* Now process default string, s points to its beginning */

	for( i = j = plevel = 1; s[i]; i++, j++ )
	{
	    if( s[i] == '\\' && (s[i+1] == '[' || s[i+1] == ']' || s[i+1] == '\\') )
		i++;
	    else if( s[i] == '[' )
		plevel++;
	    else if( s[i] == ']' )
	    {
		if( !(--plevel) )
		    break;
	    }

	    s[j] = s[i];
	}

	if( s[i] == '\0' )
	{
	    set_error( fp_stackind, "Runaway default specifier in format string" );
	    *pp = s + i;
	}
	else
	    *pp = s + i + 1;
	
	s[i] = s[j] = '\0';
	
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
	else if( size == LONG_S )
	    return( CHARPTRPTR_T );
	else
	    return( VOID_T );


      default:
	return( VOID_T );
    }
}

static size_t      TypeSize( int type )
{
    switch( type )
    {
      case CHAR_T:

	return( sizeof(char) );

      case SHORT_T:

	return( sizeof(short) );

      case INT_T:

	return( sizeof(int) );

      case LONG_T:

	return( sizeof(long) );

      case FLOAT_T:

	return( sizeof(float) );

      case DOUBLE_T:

	return( sizeof(double) );

      case CHARPTR_T:

	return( sizeof(char *) );

      case CHARPTRPTR_T:

	return( sizeof(char **) );

      case UNSIGNED_CHAR_T:

	return( sizeof(unsigned char) );

      case UNSIGNED_SHORT_T:

	return( sizeof(unsigned short) );

      case UNSIGNED_INT_T:

	return( sizeof(unsigned int) );

      case UNSIGNED_LONG_T:

	return( sizeof(unsigned long) );

      default:

	return( (size_t)1 );
    }
}


static int  single_format( int type, char *def_value, char **parg, char *rfmtA, HASHNODE *hp, void *v )
{
    int       i = 1, j, def_matched = 0;
    
    char      rfmtB[32];
    char      *arg;

        /* Holding Variables     */

    char               defv_c;
    short              defv_s;
    int                defv_i;
    long               defv_l;
    float              defv_f;
    double             defv_d;
    char              *defv_p;

    unsigned char      defv_uc;    /* These aren't really necessary but what the heck */
    unsigned short     defv_us;
    unsigned int       defv_ui;
    unsigned long      defv_ul;

        /* End Holding Variables */


    arg = *parg;

    /* Process given format string by type */

    switch( type )
    {
      case CHAR_T:

	if( !(*arg) && !def_value )
	    return( FAILURE );
	else if( !def_value )
	{
	    j = sscanf( arg, rfmtA, (char *)v );

	    if( j != 1 )
		return( FAILURE );
	}
	else
	{
	    i = sscanf( def_value, rfmtA, &defv_c );

	    strcpy( rfmtB, rfmtA );
	    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
	    /* Check if matching argument found */

	    if( *arg )
		j = sscanf( arg, rfmtB, (char *)def_vbuf, (char *)def_buf ); 

	    if( *arg && j == 1 )                    /* Use Given Value   */
		sscanf( arg, rfmtA, (char *)v );    
	    else                                    /* Use Default Value */
	    {
		if( !i )
		    return( FAILURE );                  /* Default string does not match */

		*((char *)v) = defv_c;              
		def_matched = 1;
	    }
	}
		
	break;
		
      case SHORT_T:

	if( !(*arg) && !def_value )
	    return( FAILURE );
	else if( !def_value )
	{
	    j = sscanf( arg, rfmtA, (short *)v );

	    if( j != 1 )
		return( FAILURE );
	}
	else
	{
	    i = sscanf( def_value, rfmtA, &defv_s );

	    strcpy( rfmtB, rfmtA );
	    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
	    /* Check if matching argument found */

	    if( *arg )
		j = sscanf( arg, rfmtB, (short *)def_vbuf, (char *)def_buf ); 

	    if( *arg && j == 1 )
		sscanf( arg, rfmtA, (short *)v );    /* Use Given Value   */
	    else
	    {
		if( !i )
		    return( FAILURE );                  /* Default string does not match */

		*((short *)v) = defv_s;              /* Use Default Value */
		def_matched = 1;
	    }
	}
		
	break;
		
      case INT_T:

	if( !(*arg) && !def_value )
	    return( FAILURE );
	else if( !def_value )
	{
	    j = sscanf( arg, rfmtA, (int *)v );

	    if( j != 1 )
		return( FAILURE );
	}
	else
	{
	    i = sscanf( def_value, rfmtA, &defv_i );

	    strcpy( rfmtB, rfmtA );
	    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
	    /* Check if matching argument found */

	    if( *arg )
		j = sscanf( arg, rfmtB, (int *)def_vbuf, (char *)def_buf ); 

	    if( *arg && j == 1 )
		sscanf( arg, rfmtA, (int *)v );    /* Use Given Value   */
	    else
	    {
		if( !i )
		    return( FAILURE );                  /* Default string does not match */

		*((int *)v) = defv_i;              /* Use Default Value */
		def_matched = 1;
	    }
	}
		
	break;

      case LONG_T:

	if( !(*arg) && !def_value )
	    return( FAILURE );
	else if( !def_value )
	{
	    j = sscanf( arg, rfmtA, (long *)v );

	    if( j != 1 )
		return( FAILURE );
	}
	else
	{
	    i = sscanf( def_value, rfmtA, &defv_l );

	    strcpy( rfmtB, rfmtA );
	    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
	    /* Check if matching argument found */

	    if( *arg )
		j = sscanf( arg, rfmtB, (long *)def_vbuf, (char *)def_buf ); 

	    if( *arg && j == 1 )
		sscanf( arg, rfmtA, (long *)v );    /* Use Given Value   */
	    else
	    {
		if( !i )
		    return( FAILURE );                  /* Default string does not match */

		*((long *)v) = defv_l;              /* Use Default Value */
		def_matched = 1;
	    }
	}
		
	break;

      case FLOAT_T:

	if( !(*arg) && !def_value )
	    return( FAILURE );
	else if( !def_value )
	{
	    j = sscanf( arg, rfmtA, (float *)v );

	    if( j != 1 )
		return( FAILURE );
	}
	else
	{
	    i = sscanf( def_value, rfmtA, &defv_f );

	    strcpy( rfmtB, rfmtA );
	    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
	    /* Check if matching argument found */

	    if( *arg )
		j = sscanf( arg, rfmtB, (float *)def_vbuf, (char *)def_buf ); 

	    if( *arg && j == 1 )
		sscanf( arg, rfmtA, (float *)v );    /* Use Given Value   */
	    else
	    {
		if( !i )
		    return( FAILURE );                  /* Default string does not match */

		*((float *)v) = defv_f;              /* Use Default Value */
		def_matched = 1;
	    }
	}
		
	break;

      case DOUBLE_T:

	if( !(*arg) && !def_value )
	    return( FAILURE );
	else if( !def_value )
	{
	    j = sscanf( arg, rfmtA, (double *)v );

	    if( j != 1 )
		return( FAILURE );
	}
	else
	{
	    i = sscanf( def_value, rfmtA, &defv_d );

	    strcpy( rfmtB, rfmtA );
	    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
	    /* Check if matching argument found */

	    if( *arg )
		j = sscanf( arg, rfmtB, (double *)def_vbuf, (char *)def_buf ); 

	    if( *arg && j == 1 )
		sscanf( arg, rfmtA, (double *)v );    /* Use Given Value   */
	    else
	    {
		if( !i )
		    return( FAILURE );                  /* Default string does not match */

		*((double *)v) = defv_d;              /* Use Default Value */
		def_matched = 1;
	    }
	}
		
	break;

      case CHARPTR_T:
		
	if( !(*arg) && !def_value )
	    return( FAILURE );
	else if( !def_value )
	    strcpy( (char *)v, arg );
	else
	{
	    defv_p = strdup( def_value );
		    
	    /* Check if matching argument found */

	    j = sscanf( arg, rfmtA, (int *)def_vbuf, (char *)def_buf ); 

	    if( j == 1 )
		strcpy( (char *)v, arg );           /* Use Given Value   */
	    else
	    {
		strcpy( (char *)v, defv_p );        /* Use Default Value */
		def_matched = 1;
	    }

	    Free( defv_p );
	}

	break;

      case CHARPTRPTR_T:
		
	if( !(*arg) && !def_value )
	    return( FAILURE );
	else if( !def_value )
	{
	    *((char **)v) = strdup(arg);
	}
	else
	{
	    defv_p = strdup( def_value );
		    
	    /* Check if matching argument found */

	    j = sscanf( arg, rfmtA, (int *)def_vbuf, (char *)def_buf ); 

	    if( j == 1 )
	    {
		*((char **)v) = strdup(arg);        /* Use Given Value   */
		Free( defv_p );
	    }
	    else
	    {
		*((char **)v) = defv_p;             /* Use Default Value */
		def_matched = 1;
	    }
	}

	break;

      case UNSIGNED_CHAR_T:

	if( !(*arg) && !def_value )
	    return( FAILURE );
	else if( !def_value )
	{
	    j = sscanf( arg, rfmtA, (unsigned char *)v );

	    if( j != 1 )
		return( FAILURE );
	}
	else
	{
	    i = sscanf( def_value, rfmtA, &defv_uc );

	    /* Check if matching argument found */

	    j = sscanf( arg, rfmtA, (unsigned char *)def_vbuf, (char *)def_buf ); 

	    if( j == 1 )
		sscanf( arg, rfmtA, (unsigned char *)v );   /* Use Given Value   */
	    else
	    {
		if( !i )
		    return( FAILURE );                  /* Default string does not match */

		*((unsigned char *)v) = defv_uc;            /* Use Default Value */
		def_matched = 1;
	    }
	}
		
	break;
		
      case UNSIGNED_SHORT_T:

	if( !(*arg) && !def_value )
	    return( FAILURE );
	else if( !def_value )
	{
	    j = sscanf( arg, rfmtA, (unsigned short *)v );

	    if( j != 1 )
		return( FAILURE );
	}
	else
	{
	    i = sscanf( def_value, rfmtA, &defv_us );

	    strcpy( rfmtB, rfmtA );
	    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
	    /* Check if matching argument found */

	    if( *arg )
		j = sscanf( arg, rfmtB, (unsigned short *)def_vbuf, (char *)def_buf ); 

	    if( *arg && j == 1 )
		sscanf( arg, rfmtA, (unsigned short *)v );   /* Use Given Value   */
	    else
	    {
		if( !i )
		    return( FAILURE );                  /* Default string does not match */

		*((unsigned short *)v) = defv_us;            /* Use Default Value */
		def_matched = 1;
	    }
	}
		
	break;
		
      case UNSIGNED_INT_T:

	if( !(*arg) && !def_value )
	    return( FAILURE );
	else if( !def_value )
	{
	    j = sscanf( arg, rfmtA, (unsigned int *)v );

	    if( j != 1 )
		return( FAILURE );
	}
	else
	{
	    i = sscanf( def_value, rfmtA, &defv_ui );

	    strcpy( rfmtB, rfmtA );
	    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
	    /* Check if matching argument found */

	    if( *arg )
		j = sscanf( arg, rfmtB, (unsigned int *)def_vbuf, (char *)def_buf ); 

	    if( *arg && j == 1 )
		sscanf( arg, rfmtA, (unsigned int *)v );  /* Use Given Value   */
	    else
	    {
		if( !i )
		    return( FAILURE );                  /* Default string does not match */

		*((unsigned int *)v) = defv_ui;           /* Use Default Value */
		def_matched = 1;
	    }
	}
		
	break;

      case UNSIGNED_LONG_T:

	if( !(*arg) && !def_value )
	    return( FAILURE );
	else if( !def_value )
	{
	    j = sscanf( arg, rfmtA, (unsigned long *)v );

	    if( j != 1 )
		return( FAILURE );
	}
	else
	{
	    i = sscanf( def_value, rfmtA, &defv_ul );

	    strcpy( rfmtB, rfmtA );
	    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
	    /* Check if matching argument found */

	    if( *arg )
		j = sscanf( arg, rfmtB, (unsigned long *)def_vbuf, (char *)def_buf ); 

	    if( *arg && j == 1 )
		sscanf( arg, rfmtA, (unsigned long *)v );  /* Use Given Value   */
	    else
	    {
		if( !i )
		    return( FAILURE );                  /* Default string does not match */

		*((unsigned long *)v) = defv_ul;           /* Use Default Value */
		def_matched = 1;
	    }
	}
		
	break;

      case VOID_T:
      default:

	return( FAILURE );
    }

    /* Only mark a datum as seen if it was in fact processed */

    if( !def_matched && arg[0] != '\0' )
    {
	hp->where++;
	(*parg) += strlen(arg) + 1;   /* Move to next value at char after null */
    }

    return( SUCCESS );
}


static int  parse_format( int type, char *def_value, char **parg, char *rfmtA, HASHNODE *hp, va_list *pap )
{
    int       i = 1, j, def_matched = 0;
    
    char      rfmtB[32];
    char     *arg;

        /* Holding Variables     */

    char               defv_c;
    short              defv_s;
    int                defv_i;
    long               defv_l;
    float              defv_f;
    double             defv_d;
    char              *defv_p;

    unsigned char      defv_uc;    /* These aren't really necessary but what the heck */
    unsigned short     defv_us;
    unsigned int       defv_ui;
    unsigned long      defv_ul;

        /* End Holding Variables */


    arg = *parg;

    /* Process given format string by type */

    switch( type )
    {
      case CHAR_T:

	if( !(*arg) && !def_value )
	    return( FAILURE );
	else if( !def_value )
	{
	    j = sscanf( arg, rfmtA, va_arg( (*pap), char * ) );

	    if( j != 1 )
		return( FAILURE );
	}
	else
	{
	    if( *def_value == '%' )          /* Default passed as an argument */
#if ( LINUX || DARWIN )
		defv_c = (char)va_arg( (*pap), int );
#else
		defv_c = va_arg( (*pap), char );
#endif
	    else
		i = sscanf( def_value, rfmtA, &defv_c );

	    strcpy( rfmtB, rfmtA );
	    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
	    /* Check if matching argument found */

	    if( *arg )
		j = sscanf( arg, rfmtB, (char *)def_vbuf, (char *)def_buf ); 

	    if( *arg && j == 1 )
		sscanf( arg, rfmtA, va_arg( (*pap), char * ) );    /* Use Given Value   */
	    else
	    {
		if( !i )
		    return( FAILURE );                  /* Default string does not match */

		*(va_arg( (*pap), char * )) = defv_c;              /* Use Default Value */
		def_matched = 1;
	    }
	}
		
	break;
		
      case SHORT_T:

	if( !(*arg) && !def_value )
	    return( FAILURE );
	else if( !def_value )
	{
	    j = sscanf( arg, rfmtA, va_arg( (*pap), short * ) );

	    if( j != 1 )
		return( FAILURE );
	}
	else
	{
	    if( *def_value == '%' )          /* Default passed as an argument */
#if ( LINUX || DARWIN )
		defv_s = (int)va_arg( (*pap), int );  /* CHECK CONVERSION HERE!!! */
#else
		defv_s = va_arg( (*pap), short );  /* CHECK CONVERSION HERE!!! */
#endif
	    else
		i = sscanf( def_value, rfmtA, &defv_s );

	    strcpy( rfmtB, rfmtA );
	    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
	    /* Check if matching argument found */

	    if( *arg )
		j = sscanf( arg, rfmtB, (short *)def_vbuf, (char *)def_buf ); 

	    if( *arg && j == 1 )
		sscanf( arg, rfmtA, va_arg( (*pap), short * ) );    /* Use Given Value   */
	    else
	    {
		if( !i )
		    return( FAILURE );                  /* Default string does not match */

		*(va_arg( (*pap), short * )) = defv_s;              /* Use Default Value */
		def_matched = 1;
	    }
	}
		
	break;
		
      case INT_T:

	if( !(*arg) && !def_value )
	    return( FAILURE );
	else if( !def_value )
	{
	    j = sscanf( arg, rfmtA, va_arg( (*pap), int * ) );

	    if( j != 1 )
		return( FAILURE );
	}
	else
	{
	    if( *def_value == '%' )          /* Default passed as an argument */
		defv_i = va_arg( (*pap), int );
	    else
		i = sscanf( def_value, rfmtA, &defv_i );

	    strcpy( rfmtB, rfmtA );
	    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
	    /* Check if matching argument found */

	    if( *arg )
		j = sscanf( arg, rfmtB, (int *)def_vbuf, (char *)def_buf ); 

	    if( *arg && j == 1 )
		sscanf( arg, rfmtA, va_arg( (*pap), int * ) );    /* Use Given Value   */
	    else
	    {
		if( !i )
		    return( FAILURE );                  /* Default string does not match */

		*(va_arg( (*pap), int * )) = defv_i;              /* Use Default Value */
		def_matched = 1;
	    }
	}
		
	break;

      case LONG_T:

	if( !(*arg) && !def_value )
	    return( FAILURE );
	else if( !def_value )
	{
	    j = sscanf( arg, rfmtA, va_arg( (*pap), long * ) );

	    if( j != 1 )
		return( FAILURE );
	}
	else
	{
	    if( *def_value == '%' )          /* Default passed as an argument */
		defv_l = va_arg( (*pap), long );
	    else
		i = sscanf( def_value, rfmtA, &defv_l );

	    strcpy( rfmtB, rfmtA );
	    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
	    /* Check if matching argument found */

	    if( *arg )
		j = sscanf( arg, rfmtB, (long *)def_vbuf, (char *)def_buf ); 

	    if( *arg && j == 1 )
		sscanf( arg, rfmtA, va_arg( (*pap), long * ) );    /* Use Given Value   */
	    else
	    {
		if( !i )
		    return( FAILURE );                  /* Default string does not match */

		*(va_arg( (*pap), long * )) = defv_l;              /* Use Default Value */
		def_matched = 1;
	    }
	}
		
	break;

      case FLOAT_T:

	if( !(*arg) && !def_value )
	    return( FAILURE );
	else if( !def_value )
	{
	    j = sscanf( arg, rfmtA, va_arg( (*pap), float * ) );

	    if( j != 1 )
		return( FAILURE );
	}
	else
	{
	    if( *def_value == '%' )          /* Default passed as an argument */
		defv_f = va_arg( (*pap), double );   /* CAUTION: Value Converted!! */
	    else
		i = sscanf( def_value, rfmtA, &defv_f );

	    strcpy( rfmtB, rfmtA );
	    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
	    /* Check if matching argument found */

	    if( *arg )
		j = sscanf( arg, rfmtB, (float *)def_vbuf, (char *)def_buf ); 

	    if( *arg && j == 1 )
		sscanf( arg, rfmtA, va_arg( (*pap), float * ) );    /* Use Given Value   */
	    else
	    {
		if( !i )
		    return( FAILURE );                  /* Default string does not match */

		*(va_arg( (*pap), float * )) = defv_f;              /* Use Default Value */
		def_matched = 1;
	    }
	}
		
	break;

      case DOUBLE_T:

	if( !(*arg) && !def_value )
	    return( FAILURE );
	else if( !def_value )
	{
	    j = sscanf( arg, rfmtA, va_arg( (*pap), double * ) );

	    if( j != 1 )
		return( FAILURE );
	}
	else
	{
	    if( *def_value == '%' )          /* Default passed as an argument */
		defv_d = va_arg( (*pap), double );
	    else
		i = sscanf( def_value, rfmtA, &defv_d );

	    strcpy( rfmtB, rfmtA );
	    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
	    /* Check if matching argument found */

	    if( *arg )
		j = sscanf( arg, rfmtB, (double *)def_vbuf, (char *)def_buf ); 

	    if( *arg && j == 1 )
		sscanf( arg, rfmtA, va_arg( (*pap), double * ) );    /* Use Given Value   */
	    else
	    {
		if( !i )
		    return( FAILURE );                  /* Default string does not match */

		*(va_arg( (*pap), double * )) = defv_d;              /* Use Default Value */
		def_matched = 1;
	    }
	}
		
	break;

      case CHARPTR_T:
		
	if( !(*arg) && !def_value )
	    return( FAILURE );
	else if( !def_value )
	    strcpy( va_arg( (*pap), char * ), arg );
	else
	{
	    if( *def_value == '%' )          /* Default passed as an argument */
		defv_p = strdup( va_arg( (*pap), char * ) );
	    else
		defv_p = strdup( def_value );
		    
	    /* Check if matching argument found */

	    j = sscanf( arg, rfmtA, (int *)def_vbuf, (char *)def_buf ); 

	    if( j == 1 )
		strcpy( va_arg( (*pap), char * ), arg );           /* Use Given Value   */
	    else
	    {
		strcpy( va_arg( (*pap), char * ), defv_p );        /* Use Default Value */
		def_matched = 1;
	    }

	    Free( defv_p );
	}

	break;

      case CHARPTRPTR_T:
		
	if( !(*arg) && !def_value )
	    return( FAILURE );
	else if( !def_value )
	{
	    *(va_arg( (*pap), char ** )) = strdup(arg);
	}
	else
	{
	    if( *def_value == '%' )          /* Default passed as an argument */
		defv_p = strdup( va_arg( (*pap), char * ) );
	    else
		defv_p = strdup( def_value );
		    
	    /* Check if matching argument found */

	    j = sscanf( arg, rfmtA, (int *)def_vbuf, (char *)def_buf ); 

	    if( j == 1 )
	    {
		*(va_arg( (*pap), char ** )) = strdup(arg);        /* Use Given Value   */
		Free( defv_p );
	    }
	    else
	    {
		*(va_arg( (*pap), char ** )) = defv_p;             /* Use Default Value */
		def_matched = 1;
	    }
	}

	break;

      case UNSIGNED_CHAR_T:

	if( !(*arg) && !def_value )
	    return( FAILURE );
	else if( !def_value )
	{
	    j = sscanf( arg, rfmtA, va_arg( (*pap), unsigned char * ) );

	    if( j != 1 )
		return( FAILURE );
	}
	else
	{
	    if( *def_value == '%' )          /* Default passed as an argument */
#if ( LINUX || DARWIN )
		defv_uc = (unsigned char)va_arg( (*pap), int );
#else
		defv_uc = va_arg( (*pap), unsigned char );
#endif
	    else
		i = sscanf( def_value, rfmtA, &defv_uc );

	    /* Check if matching argument found */

	    j = sscanf( arg, rfmtA, (unsigned char *)def_vbuf, (char *)def_buf ); 

	    if( j == 1 )
		sscanf( arg, rfmtA, va_arg( (*pap), unsigned char * ) );   /* Use Given Value   */
	    else
	    {
		if( !i )
		    return( FAILURE );                  /* Default string does not match */

		*(va_arg( (*pap), unsigned char * )) = defv_uc;            /* Use Default Value */
		def_matched = 1;
	    }
	}
		
	break;
		
      case UNSIGNED_SHORT_T:

	if( !(*arg) && !def_value )
	    return( FAILURE );
	else if( !def_value )
	{
	    j = sscanf( arg, rfmtA, va_arg( (*pap), unsigned short * ) );

	    if( j != 1 )
		return( FAILURE );
	}
	else
	{
	    if( *def_value == '%' )          /* Default passed as an argument */
#if ( LINUX || DARWIN )
		defv_us = (unsigned short)va_arg( (*pap), int );  /* CHECK CONVERSION HERE!!! */
#else
		defv_us = va_arg( (*pap), unsigned short );  /* CHECK CONVERSION HERE!!! */
#endif
	    else
		i = sscanf( def_value, rfmtA, &defv_us );

	    strcpy( rfmtB, rfmtA );
	    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
	    /* Check if matching argument found */

	    if( *arg )
		j = sscanf( arg, rfmtB, (unsigned short *)def_vbuf, (char *)def_buf ); 

	    if( *arg && j == 1 )
		sscanf( arg, rfmtA, va_arg( (*pap), unsigned short * ) );   /* Use Given Value   */
	    else
	    {
		if( !i )
		    return( FAILURE );                  /* Default string does not match */

		*(va_arg( (*pap), unsigned short * )) = defv_us;            /* Use Default Value */
		def_matched = 1;
	    }
	}
		
	break;
		
      case UNSIGNED_INT_T:

	if( !(*arg) && !def_value )
	    return( FAILURE );
	else if( !def_value )
	{
	    j = sscanf( arg, rfmtA, va_arg( (*pap), unsigned int * ) );

	    if( j != 1 )
		return( FAILURE );
	}
	else
	{
	    if( *def_value == '%' )          /* Default passed as an argument */
		defv_ui = va_arg( (*pap), unsigned int );
	    else
		i = sscanf( def_value, rfmtA, &defv_ui );

	    strcpy( rfmtB, rfmtA );
	    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
	    /* Check if matching argument found */

	    if( *arg )
		j = sscanf( arg, rfmtB, (unsigned int *)def_vbuf, (char *)def_buf ); 

	    if( *arg && j == 1 )
		sscanf( arg, rfmtA, va_arg( (*pap), unsigned int * ) );  /* Use Given Value   */
	    else
	    {
		if( !i )
		    return( FAILURE );                  /* Default string does not match */

		*(va_arg( (*pap), unsigned int * )) = defv_ui;           /* Use Default Value */
		def_matched = 1;
	    }
	}
		
	break;

      case UNSIGNED_LONG_T:

	if( !(*arg) && !def_value )
	    return( FAILURE );
	else if( !def_value )
	{
	    j = sscanf( arg, rfmtA, va_arg( (*pap), unsigned long * ) );

	    if( j != 1 )
		return( FAILURE );
	}
	else
	{
	    if( *def_value == '%' )          /* Default passed as an argument */
		defv_ul = va_arg( (*pap), unsigned long );
	    else
		i = sscanf( def_value, rfmtA, &defv_ul );

	    strcpy( rfmtB, rfmtA );
	    strcat( rfmtB, " %s" ); /* This catches the extra stuff */
		    
	    /* Check if matching argument found */

	    if( *arg )
		j = sscanf( arg, rfmtB, (unsigned long *)def_vbuf, (char *)def_buf ); 

	    if( *arg && j == 1 )
		sscanf( arg, rfmtA, va_arg( (*pap), unsigned long * ) );  /* Use Given Value   */
	    else
	    {
		if( !i )
		    return( FAILURE );                  /* Default string does not match */

		*(va_arg( (*pap), unsigned long * )) = defv_ul;           /* Use Default Value */
		def_matched = 1;
	    }
	}
		
	break;

      case VOID_T:
      default:

	return( FAILURE );
    }

    /* Only mark a datum as seen if it was in fact processed */

    if( !def_matched && arg[0] != '\0' )
    {
	hp->where++;
	(*parg) += strlen(arg) + 1;   /* Move to next value at char after null */
    }

    return( SUCCESS );
}

static int parse_vector( int vsize, char *def_value, char **pval, char *rfmtA, HASHNODE *hp, va_list *pap )
{
    char          code;

    int           i, j, numel = 0;
    int           suppress, size, ftype, type;
    int           status = SUCCESS;

    size_t        typesize;

    char          *sp, *q, *last;
    char          rfmtB[32];

    void          *v;


    /*
     * First, determine how many elements the vector has.
     *
     * There are several possibilities, given by the presence (or absence)
     * of an appropriate code before the formats, with the layout being
     * either [code:format] or [format].
     *
     *  Code:
     *
     *   =     Next argument on list is a pointer to an integer; store the number of elements
     *         in that integer
     *   %     Next argument on list is an integer giving the number of elements
     *
     *   -     Previous format string is a required integer, use that value for num elements.
     *
     *  number Use given number
     *
     * If no code is present, then all of the remaining elements are used until a failure
     * to match occurs.  Thus, with no code, defaults are never used.
     *
     */

    if( !def_value )
    {
	set_error( fp_stackind, "Empty format list for Vector (V) format string" );
	return( FAILURE );
    }

    if( def_value[0] == ':' )
    {
	code = '\0';
	sp = def_value + 1;
    }
    else if( def_value[0] == '%' && (isspace(def_value[1]) || def_value[1] == ':') )
    {
	code = '%';
	
	for( sp = def_value + 1; *sp && isspace(*sp); sp++ )
	    ;

	if( *sp == '\0' || *sp == ']' )
	{
	    set_error( fp_stackind, "Empty format list for Vector (V) format string" );
	    return( FAILURE );
	}
	else if( *sp != ':' )
	    set_error( fp_stackind, "No closing ':' in vector (V) format specification" );
	else
	    sp++;
    }
    else if( def_value[0] == '-' || def_value[0] == '=' )
    {
	code = def_value[0];

	for( sp = def_value + 1; *sp && isspace(*sp); sp++ )
	    ;

	if( *sp == '\0' || *sp == ']' )
	{
	    set_error( fp_stackind, "Empty format list for Vector (V) format string" );
	    return( FAILURE );
	}
	else if( *sp != ':' )
	    set_error( fp_stackind, "No closing ':' in vector (V) format specification" );
	else
	    sp++;
    }
    else if( isdigit( def_value[0] ) )
    {
	code = '1';
	sscanf( def_value, "%d", &numel );

	for( sp = def_value + 1; *sp && (isdigit(*sp) || isspace(*sp)); sp++ )
	    ;

	if( *sp == '\0' || *sp == ']' )
	{
	    set_error( fp_stackind, "Empty format list for Vector (V) format string" );
	    return( FAILURE );
	}
	else if( *sp != ':' )
	    set_error( fp_stackind, "No closing ':' in vector (V) format specification" );
	else
	    sp++;
    }
    else
    {
	code = '\0';
	sp = def_value;
    }

    sp = strchr( sp, '%' );    /* Find leading % */

    if( !sp )
    {
	set_error( fp_stackind, "Empty format list for Vector (V) format string" );
	return( FAILURE );
    }

    /* Now, we've encoded how to find the number of elements and sp points to the format */


    /* Parse the format string */

    suppress = '\0';
    size = STANDARD_S;
    ftype = VOID_F;
    rfmtA[1] = '\0';

    if( !(q = cnvt_code( sp, &suppress, &size, &ftype )) )
    {
	set_error( fp_stackind, "Cannot parse component format string in vector (V) format." );
	return( FAILURE );
    }

    def_value = cnvt_def( q + 1, ftype, &q ); 
    sp++;               

	/*
	 * Where are we:                                                   
	 *   sp now points to the null-padded format string without '%'     
	 *   def_value if non-null points to the null padded default string
	 *   q now points to the character beyond the format string        
	 *
	 * Next, we parse the format and assign the values unless suppress
	 * has been indicated by a '*' in the format string after '%'
	 *
	 */

    if( suppress )
    {
	set_error( fp_stackind, "Cannot suppress matrix or vector format types" );
	return( FAILURE );
    }
    else if( ftype == VECTOR_F || ftype == MATRIX_F )
    {
	set_error( fp_stackind, "Cannot nest matrix (M) or vector (V) format codes" );
	return( FAILURE );
    }

    type = assoc_type( ftype, size );   /* Find underlying read type    */
    strcat( rfmtA, sp );                /* Set conversion format to use */

    strcpy( rfmtB, rfmtA );
    strcat( rfmtB, " %s" );             /* This catches the extra stuff */

    switch( code )    /* Find and assign number of elements depending on code  */
    {
      case '%':

	numel = va_arg( (*pap), int );
	break;

      case '-':       /* Use last stored field value (assume required) */

	if( !hp->where )
	{
	    set_error( fp_stackind, "No previous count to use in vector (V) format with - indicator" );
	    return( FAILURE );
	}
	
	last = *pval - ((hp->where + 1 == hp->count) ? 2 : 1);

	while( *last )
	    last--;

	last++;

	j = sscanf( last, "%d", &numel );

	if( j != 1 )
	{
	    set_error( fp_stackind, "Could not match previous count to use in vector (V) format with - indicator" );
	    return( FAILURE );
	}

	break;

      case '=':       /* Use all remaining elements that match format */
      case '\0':      /* No defaults allowed                          */

	for( numel = 0, last = *pval; last[0] != '\0'; numel++ )
	{
	    j = sscanf( last, rfmtB, def_buf, def_vbuf );

	    if( j != 1 )
		break;
	    
	    last += strlen( last ) + 1;
	}

	if( code == '=' )
	    *(va_arg( (*pap), int * )) = numel;

	def_value = NULL;
	
	break;
    }

    typesize = TypeSize(type);

    if( vsize == LONG_S )
    {
	v = emalloc( numel * typesize );
	*(va_arg( (*pap), void **)) = v;
    }
    else
	v = va_arg( (*pap), void * );
	
    for( i = 0, status = SUCCESS; i < numel && *pval != '\0'; i++ )
    {
	status = single_format( type, def_value, pval, rfmtA, hp, (void *)((char *)v + i*typesize) );

	if( status == FAILURE )
	    break;
    }

    return( status );
}

static int parse_matrix( int msize, char *def_value, char **pval, char *fmt, HASHNODE *hp, va_list *pap )
{
    char          code;

    int           i, j, k;
    int	          nrow = 0, ncol = 0;
    int           suppress, size, ftype, type;
    int           status = SUCCESS;

    size_t        tabsize;

    char         *sp, *q, *last;
    char          rfmtA[32], rfmtB[32];

    size_t       *tsiz;
    char        **cols, **defs;
    void        **ptrs;

    /*
     * First, determine how to find out the number of rows of the matrix
     *
     * There are several possibilities, given by the presence (or absence)
     * of an appropriate code before the formats, with the layout being
     * either [code:format] or [format].
     *
     *  Code:
     *
     *   =     Next argument on list is a pointer to an integer; store the number of elements
     *         in that integer
     *   %     Next argument on list is an integer giving the number of elements
     *
     *   -     Previous format string is a required integer, use that value for num elements.
     *
     *  number Use given number
     *
     * If no code is present or '=' is given, then all of the remaining elements are used
     * until a failure to match occurs.  Thus, with no code or '=', defaults are never used.
     *
     */

    if( !def_value )
    {
	set_error( fp_stackind, "Empty format list for Vector (V) format string" );
	return( FAILURE );
    }

    if( def_value[0] == ':' )
    {
	code = '\0';
	sp = def_value + 1;
    }
    else if( def_value[0] == '%' && (isspace(def_value[1]) || def_value[1] == ':') )
    {
	code = '%';
	
	for( sp = def_value + 1; *sp && isspace(*sp); sp++ )
	    ;

	if( *sp == '\0' || *sp == ']' )
	{
	    set_error( fp_stackind, "Empty format list for Vector (V) format string" );
	    return( FAILURE );
	}
	else if( *sp != ':' )
	    set_error( fp_stackind, "No closing ':' in vector (V) format specification" );
	else
	    sp++;
    }
    else if( def_value[0] == '-' || def_value[0] == '=' )
    {
	code = def_value[0];

	for( sp = def_value + 1; *sp && isspace(*sp); sp++ )
	    ;

	if( *sp == '\0' || *sp == ']' )
	{
	    set_error( fp_stackind, "Empty format list for Vector (V) format string" );
	    return( FAILURE );
	}
	else if( *sp != ':' )
	    set_error( fp_stackind, "No closing ':' in vector (V) format specification" );
	else
	    sp++;
    }
    else if( isdigit( def_value[0] ) )
    {
	code = '1';
	sscanf( def_value, "%d", &nrow );

	for( sp = def_value + 1; *sp && (isdigit(*sp) || isspace(*sp)); sp++ )
	    ;

	if( *sp == '\0' || *sp == ']' )
	{
	    set_error( fp_stackind, "Empty format list for Vector (V) format string" );
	    return( FAILURE );
	}
	else if( *sp != ':' )
	    set_error( fp_stackind, "No closing ':' in vector (V) format specification" );
	else
	    sp++;
    }
    else
    {
	code = '\0';
	sp = def_value;
    }

    sp = strchr( sp, '%' );    /* Find leading % */

    if( !sp )
    {
	set_error( fp_stackind, "Missing format list in matrix (M) format" );
	return( FAILURE );
    }
    
    /* Now, we've encoded how to find the number of elements and sp points to the format */
    /* Begin by counting the number of columns and recording this information in a table */

    cols = Calloc( TABLE_DEFAULT, char * );
    defs = Calloc( TABLE_DEFAULT, char * );
    tsiz = Calloc( TABLE_DEFAULT, size_t );
    tabsize = TABLE_DEFAULT;

    for( q = sp, ncol = 0; q; q = strchr( q, '%' ), ncol++ )
    {
	cols[ncol] = q;

	if( !(q = cnvt_code( q, &suppress, &size, &ftype )) )
	{
	    set_error( fp_stackind, "Cannot parse component format string in matrix (M) format." );
	    Free( cols );
	    Free( defs );
	    return( FAILURE );
	}

	if( suppress )
	{
	    set_error( fp_stackind, "Cannot suppress matrix or vector format types" );
	    Free( cols );
	    Free( defs );
	    return( FAILURE );
	}
	else if( ftype == VECTOR_F || ftype == MATRIX_F )
	{
	    set_error( fp_stackind, "Cannot nest matrix (M) or vector (V) format codes" );
	    Free( cols );
	    Free( defs );
	    return( FAILURE );
	}

	tsiz[ncol] = TypeSize( assoc_type( ftype, size ) ); 

	defs[ncol] = cnvt_def( q + 1, ftype, &q ); 

	if( ncol + 1 > tabsize )
	{
	    tabsize *= REALLOC_EXPAND;
	    
	    cols = Realloc( cols, tabsize, char * );
	    defs = Realloc( defs, tabsize, char * );
	    tsiz = Realloc( tsiz, tabsize, size_t );
	}
    }

    /*
     * Now, examine the first row to see which match 
     * Use these matchings throughout, and use them to set number of elements
     *
     */

    switch( code )    /* Find and assign number of elements depending on code  */
    {
      case '%':

	nrow = va_arg( (*pap), int );
	break;

      case '-':       /* Use last stored field value (assume required) */

	if( !hp->where )
	{
	    set_error( fp_stackind, "No previous count to use in matrix (M) format with - indicator" );
	    Free( cols );
	    Free( defs );
	    Free( tsiz );
	    return( FAILURE );
	}
	
	last = *pval - ((hp->where + 1 == hp->count) ? 2 : 1);

	while( *last )
	    last--;

	last++;

	j = sscanf( last, "%d", &nrow );

	if( j != 1 )
	{
	    set_error( fp_stackind, "Could not match previous count to use in matrix (M) format with - indicator" );
	    Free( cols );
	    Free( defs );
	    Free( tsiz );
	    return( FAILURE );
	}

	break;

      case '=':       /* Use all remaining elements that match format */
      case '\0':      /* No defaults allowed                          */

	for( i = 0, k = 0, last = *pval; *last; k++ )
	{
	    sprintf( rfmtB, "%s %%s", cols[i] );
	    j = sscanf( last, rfmtB, def_buf, def_vbuf );
	
	    if( j != 1 && defs[i] == NULL )  /* No match and no default */
	    {
		Free( cols );
		Free( defs );
		Free( tsiz );
		return( FAILURE );
	    }
	    else if( j != 1 )
	    {
		j = sscanf( defs[i], rfmtB, def_buf, def_vbuf );

		if( j != 1 )                /* No match and default doesn't match either */
		{
		    Free( cols );
		    Free( defs );
		    Free( tsiz );
		    return( FAILURE );
		}
	    }
	    else
		last += strlen(last) + 1;
	
	    i = (i + 1) % ncol;
	}

	if( i )   /* Did not complete last row */
	{
	    k += ncol - i;

	    for( ; i < ncol; i++ )
	    {
		sprintf( rfmtB, "%s %%s", cols[i] );
		j = sscanf( defs[i], rfmtB, def_buf, def_vbuf );

		if( j != 1 )                /* No match and default doesn't match either */
		{
		    Free( cols );
		    Free( defs );
		    Free( tsiz );
		    return( FAILURE );
		}
	    }
	}

	nrow = k/ncol;

	if( code == '=' )
	    *(va_arg( (*pap), int * )) = nrow;

	break;
    }

    /* Get pointers to columns, treat generically */

    ptrs = Calloc( tabsize, void * );

    for( i = 0; i < ncol; i++ )
    {
	if( msize == LONG_S )
	{
	    ptrs[i] = emalloc( nrow * tsiz[i] );
	    *(va_arg( (*pap), void **)) = ptrs[i];
	}
	else
	    ptrs[i] = va_arg( (*pap), void * );
    }
    
    /* Parse the format strings */

    for( i = 0; i < nrow && *pval != '\0'; i++ )
    {
	for( j = 0; j < ncol && *pval != '\0'; j++ )
	{
	    cnvt_code( cols[j], &suppress, &size, &ftype );   /* We know from above this is non-null */

	    type = assoc_type( ftype, size );   /* Find underlying read type    */

	    strcpy( rfmtA, cols[j] );           /* Set conversion format to use */
	    sprintf( rfmtB, "%s %%s", rfmtA );
	
	    status = single_format( type, defs[j], pval, rfmtA, hp, (void *)((char *)ptrs[j] + i*tsiz[j]) );

	    if( status == FAILURE )
	    {
		Free( cols );
		Free( defs );
		Free( tsiz );
		Free( ptrs );
		
		return( FAILURE );
	    }
	}
    }
    
    Free( cols );
    Free( defs );
    Free( tsiz );
    Free( ptrs );
		
    return( status );
}



static char  *gobble_ws( char *s )
{
    while( *s && isspace(*s) )
	s++;

    return( s );
}


static char  *gobble_line( char *s )
{
    while( *s )
	s++;

    return( s );
}

static char  *skip_quoted( char *s, char *stopset, char *stop )
{
    int       len, sonq = 1;
    char      buf[32];
    char     *p, *r;

    if( !strchr( stopset, '\"' ) )
    {
	len = strlen( stopset );
	sonq = 0;

	if( len > 28 )
	    r = Malloc( len + 3, char );
	else
	    r = buf;

	strcpy( r, stopset );
	strcat( r, "\"" );
    }
    else
	r = stopset;

    p = s + 1;
    len = strlen( s );

    for( ;; )
    {
	p = strpbrk( p, r );

	if( !p || !(*p == '\"' && p[-1] == '\\' && !sonq) )
	    break;
    }

    if( p )
	*stop = *p;
    else
    {
	p = s + len;
	*stop = '\0';
    }

    if( r != buf && r != stopset )
	Free( r );

    return( p );
}

static char  *skip_non_ws( char *s, char *stopset, char *stop )
{
    int      len;

    char     buf[32];
    char     *p, *r;
    
    len = strlen( stopset );

    if( len > 21 )
	r = Malloc( len + 10, char );
    else
	r = buf;

    strcpy( r, stopset );
    strcat( r, " \t\n\r\f\v\b\a" );

    p = s + 1;
    len = strlen( s );

    p = strpbrk( p, r );

    if( p )
	*stop = *p;
    else
    {
	p = s + len;
	*stop = '\0';
    }

    if( r != buf )
	Free( r );

    return( p );
}


/* scan until an <end> token is reached or EOF, discarding tokens */

static void      scan_to_end( char **psp, char **pbuf, int *pbuflen )
{
    int          tc;
    char         *q;

    q = *pbuf;
    
    for( ;; )
    {
	tc = get_token( psp, &q, pbuf, pbuflen, END_AT_EOF );

	if( !tc || !strcmp( q - tc, "<end>" ) )
	    break;
    }
}


/*
 * Built-in Commands
 *
 * Each command receives the position in the current line directly
 * after the invoking command construct <name>.
 *
 */


static void   	 do_append( char *tp )
{
    int           len, cnt, off = 0;
    char         *tag, *old, *p;
    HASHNODE     *hp;

    tag = get_tag( &tp );

    if( check_error() )
    {
	set_error( fp_stackind, "Could not read tag after <append> command" );
	return;
    }

    p = data_buf;
    cnt = get_dat( &tp, &p, &data_buf, &data_bufsize );
    off = p - data_buf;

    if( check_error() )
    {
	set_error( fp_stackind, "Could not read data after <append> %s command", tag );
	return;
    }
    
    hp = lookup( tag );

    if( hp )
    {
	old = hp->value;
	
	for( p = old, len = strlen( p ); len > 0; len = strlen( p ) )
	    p += len + 1;

	len = p - old;

	hp->value = Calloc( len + off + 2, char );
	memcpy( hp->value,       old,      len );
	memcpy( hp->value + len, data_buf, off );
	hp->count += cnt;

	Free( old );
    }
    else
    {
	hp = install( tag, NULL );
	hp->value = Calloc( off + 2, char );
	hp->count = cnt;

	memcpy( hp->value, data_buf, off );
    }
}

static void   	 do_arg( char *tp )
{
    int        any;
    char      *tag, *p;

    tag = get_tag( &tp );

    if( check_error() )
    {
	set_error( fp_stackind, "Could not read tag after <arg> command" );
	return;
    }

    p = data_buf;
    any = get_token( &tp, &p, &data_buf, &data_bufsize, ONE_LINE );

    if( any )
	inpf_arg( tag, data_buf );
    else
	inpf_arg( tag, NULL );
}

static void      do_cond( char *tp )
{
    int           len, toks, odd_end = 0;
    char         *tag, *var, *p, *q;
    
    tag = get_tag( &tp );

    if( check_error() )
    {
	set_error( fp_stackind, "Could not read tag after <cond> command" );
	scan_to_end( &tp, &data_buf, &data_bufsize );
	return;
    }

    p = data_buf;
    toks = get_token( &tp, &p, &data_buf, &data_bufsize, ONE_LINE );

    if( !toks || check_error() )
    {
	set_error( fp_stackind, "Could not read switch after <cond> %s command", tag );
	scan_to_end( &tp, &data_buf, &data_bufsize );
	return;
    }

    var = data_buf;

    for( ;; )
    {
	q = p;
	len = get_token( &tp, &q, &data_buf, &data_bufsize, END_AT_EOF );

	if( !len )
	{
	    odd_end = 1;
	    break;
	}
	else if( !strcmp( p, var ) )
	{
	    q = p;
	    
	    get_dat( &tp, &q, &data_buf, &data_bufsize );

	    if( check_error() )
	    {
		set_error( fp_stackind, "Could not read data after case %s in <cond> command", var );
		return;
	    }

	    install( tag, p );
	    scan_to_end( &tp, &data_buf, &data_bufsize );
	    break;
	}
	else if( !strcmp( p, "<default>" ) )
	{
	    q = p;
	    
	    get_dat( &tp, &q, &data_buf, &data_bufsize );

	    if( check_error() )
	    {
		set_error( fp_stackind, "Could not read data after case %s in <cond> command", var );
		return;
	    }

	    install( tag, p );
	    scan_to_end( &tp, &data_buf, &data_bufsize );
	    break;
	}
	else if( !strcmp( p, "<end>" ) )
	    break;
	else
	{
	    tp = get_line( END_AT_EOF );

	    if( !tp )
	    {
		odd_end = 1;
		break;
	    }
	}
    }

    if( odd_end )
	set_error( fp_stackind, "Unexpected end-of-file encountered during <cond> command" );

    return;
}

static void   	 do_data( char *tp )
{
    int        tc, len, cnt;

    char       *dp, *tag;
    HASHNODE   *hp;


    tag = get_tag( &tp );

    if( check_error() )
    {
	set_error( fp_stackind, "Could not read tag after <data> command" );
	scan_to_end( &tp, &data_buf, &data_bufsize );
	return;
    }

    hp = install( tag, NULL );     /* Prepare to install new data */

    /* Loop until line beginning with <end> processing and appending data */

    dp = data_buf;

    for( cnt = len = 0;; cnt++ )
    {
	tc = get_token( &tp, &dp, &data_buf, &data_bufsize, END_AT_EOF );
	len += tc;

	if( !tc )
	{
	    set_error( fp_stackind, "Runaway <data> command, no matching <end>" );
	    break;
	}
	else if( !strcmp( dp - tc, "<end>" ) )
	{
	    dp -= tc;
	    len -= tc;
	    break;
	}
    }

    *dp++ = '\0';     /* Cap it off with extra null */
    len++;

    if( len > data_bufsize )
	set_error( fp_stackind, "Exceeded Data Buffer Size" );

    hp->value = Malloc( len + 1, char );
    memcpy( hp->value, data_buf, len );
    hp->count = cnt;
}

static void   	 do_include( char *tp )
{
    int          any;
    char        *p;
    FILE        *fp;

    p = data_buf;
    any = get_token( &tp, &p, &data_buf, &data_bufsize, ONE_LINE );

    if( any )
    {
	fp = fopen( data_buf, "r" );

	if( !fp )
	{
	    set_error( fp_stackind, "Cannot open file %s", data_buf );
	    return;
	}

	fp_stack[++fp_stackind] = fp;
	Names[fp_stackind] = strdup( data_buf );
	Line_num[fp_stackind] = 0;
    }
}

static void   	 do_require( char *tp )
{
    char         *tag;

    tag = get_tag( &tp );

    if( !lookup( tag ) )
	set_error( fp_stackind, "Required tag %s not set in input file", tag );
}


/* Hash Table Support */

static int       hashvalue;

static HASHNODE *lookup( char *tag )
{
    HASHNODE    	*cmp;
    int                 (*sc)(const char *, const char *);

    sc = (case_fold) ? strcasecmp : strcmp;

    hashvalue = hash( tag );
    
    for( cmp = inpf_hash_table[hashvalue]; cmp ; cmp = cmp->next )
        if( !(*sc)( tag, cmp->tag) )
            return(cmp);

    return(NULL);
}

static HASHNODE *install( char *tag, char *data )
{
    int                  cnt, len;
    char                *p;
    HASHNODE            *new;

    new = lookup( tag );

    if( !new )
    {
	new = Malloc( 1, HASHNODE );

	new->tag     = strdup( tag );
	new->where   = 0;
	new->matched = 0;
	new->command = NULL;
	new->set     = 1;
	new->count   = 0;

	/*
	if( fp_stackind >= 0 )
	{
	    sprintf( err_buf, "File %s, Line %d", Names[fp_stackind], Line_num[fp_stackind] );
	    new->info = strdup( err_buf );
	}
	else
	    new->info = Calloc( 2, char );
	    */

	if( !data )
	    new->value = NULL;
	else if( *data == '\0' )
	    new->value = Calloc( 3, char );  /* Empty Data Array \0\0 */
	else
	{
	    for( cnt = 0, p = data, len = strlen( p ); len > 0; len = strlen( p ), cnt++ )
		p += len + 1;

	    new->value = Calloc( p - data + 1, char );
	    memcpy( new->value, data, p - data );
	    new->count = cnt;
	}

	new->next = inpf_hash_table[hashvalue];
	inpf_hash_table[hashvalue] = new;
    }
    else
    {
	if( new->value )
	    Free( new->value );

	new->where = 0;
	new->count = 0;

	if( !data )
	    new->value = NULL;
	else if( *data == '\0' )
	    new->value = strdup( "" );
	else
	{
	    for( cnt = 0, p = data, len = strlen( p ); len > 0; len = strlen( p ), cnt++ )
		p += len + 1;

	    new->value = Calloc( p - data + 1, char );
	    memcpy( new->value, data, p - data );
	    new->count = cnt;
	}
    }

    return( new );
}

static int       hash( char *hashstr )
{
    unsigned int	g;

    for( hashvalue = 0; *hashstr; )
    {
        hashvalue = (hashvalue << 4) + *hashstr++;

        if( (g = (hashvalue & 0xf0000000LU)) )
        {
            hashvalue = hashvalue ^ (g >> 24);
            hashvalue = hashvalue ^ g;
        }
    } 

    hashvalue %= inpf_hash_length;

    return( hashvalue );
}


/** File Input Management **/


/*
 * GET_LINE is the generic input function.  It reads a line from
 * the appropriate input file, managing the file stack and adjusting
 * the line numbers as necessary
 *
 */

static char    *get_line( int state )
{
    int         len;
    
    while( !fgets( line_buf, line_bufsize, fp_stack[fp_stackind] ) )
    {
	if( state == END_AT_EOF && feof(fp_stack[fp_stackind]) )
	    return( NULL );

	if( fp_stackind > 0 )
	{
	    fclose( fp_stack[fp_stackind] );
	    fp_stack[fp_stackind--] = NULL;
	}
	else
	    return( NULL );
    }

    Line_num[fp_stackind]++;
    len = strlen( line_buf );

    while( (len == line_bufsize - 1) && (line_buf[len] != '\n') )
    {
	/* Increase bufsize and read in rest of the line */

	line_bufsize *= REALLOC_EXPAND;
	line_buf = Realloc( line_buf, line_bufsize, char );

	if( !fgets( line_buf + len, line_bufsize - len, fp_stack[fp_stackind] ) )
	    break;
	    
	len = strlen( line_buf );
    }

    return( line_buf );
}


/*
 * GET_TOKEN is the lexical analyzer/pre-processor.  It decomposes the input 
 * stream into a sequence of data tokens, after processing strings, tag references,
 * escape sequences, and so forth.
 *
 * This is the basic function used by get_dat and do_data for reading data tokens
 * after tags.
 *
 */

enum states
{
    DONE, END_OF_LINE, SKIPPING_WS, SKIPPING_COMMENT, PROCESSING_REFERENCE,
    PROCESSING_ESCAPE,  PROCESSING_STRING, PROCESSING_TOK 
};

#define        tpush(i)                (state[++ind] = (i))
#define        tpop()                  ((void)state[ind--])          


/*
 * buf_add_*() add characters to the buffer at the given location bp.
 * These first check that there is enough space in the buffer, enlarging
 * it if not.  Then, they add the specified characters.  The new location
 * in the buffer is returned.
 *
 * buf_inc_pntr() increments the buffer pointer after checking the length.
 * buf_check() does the checking and reallocating.
 *
 */

static INLINE
char*           buf_add_string( char* bp, char* sp, size_t len, char** bufp, int* buflenp )
{
    bp = buf_check( len, bp, bufp, buflenp );

    memcpy( bp, sp, len );
    bp += len;

    return( bp );
}

static INLINE
char*           buf_add_char( char* bp, char c, char** bufp, int* buflenp )
{
    bp = buf_check( 1, bp, bufp, buflenp );

    *bp++ = c;
    
    return( bp );
}

static INLINE
char*           buf_inc_pntr( char* bp, char** bufp, int* buflenp )
{
    bp = buf_check( 1, bp, bufp, buflenp );

    return( ++bp );
}

static INLINE
char*           buf_check( size_t inc, char* bp, char** bufp, int* buflenp )
{
    static int realloc_expand = 1;

    unsigned long        diff;

    diff = bp - *bufp;            /* Current offset in buffer */

    if( diff + inc > *buflenp )
    {
        char*                newp;
        unsigned long        newsize;

        newsize = *buflenp + (*buflenp)/realloc_expand;
    
        newp = (char*) realloc( *bufp, newsize * sizeof(char) );

        while( !newp && (diff + inc <= newsize) )   /* Too much memory requested; try smaller size */
        {
            realloc_expand++;
            
            newsize = *buflenp + (*buflenp)/realloc_expand;
            newp = (char*) realloc( *bufp, newsize * sizeof(char) );
        }

        if( !newp )
        {
            set_error( fp_stackind, "Cannot allocate sufficient memory for token, %d bytes requested.", newsize );
            inpf_perror( );
            exit( 1 );  /* ATTN: What to do here?? */
        }
        else
        {
            *bufp = newp;
            *buflenp = newsize;
            bp = *bufp + diff;
        }
    }

    return( bp );
}

/*
 * get_token() Reads a token from the string *psp and stores it in the
 * buffer *bufp (with length *buflenp) starting at the position *pbp.
 * The buffer is enlarged if necessary to fit the token.
 *
 * The input argument code determines the boundaries of the input string;
 * see get_line().
 *
 */

static int     get_token( char** psp, char** pbp, char** bufp, int* buflenp, int code )
{
    int         ind = -1;
    int         i, j, len = 0;

    char*       bp;
    char*       sp;
    char*       q;
    
    int         state[8];

    HASHNODE*   hp;


    bp = *pbp;
    sp = *psp;

    if( *sp == '\n' || *sp == '\0' )
	tpush(END_OF_LINE);
    else if( isspace(*sp) )
	tpush(SKIPPING_WS);
    else if( *sp == escape_char && isspace(sp[1]) )
    {
	tpush(SKIPPING_WS);
	tpush(PROCESSING_ESCAPE);
    }
    else if( *sp == comment_char )
	tpush(SKIPPING_COMMENT);
    else if( *sp == var_ref_char )
    {
	tpush(PROCESSING_TOK);
	tpush(PROCESSING_REFERENCE);
    }
    else if( *sp == quote_char )
    {
	sp++;
	tpush(PROCESSING_TOK);
	tpush(PROCESSING_STRING);
    }
    else
	tpush(PROCESSING_TOK);

    for( ;; )
    {
	if( state[ind] == DONE )
	    break;

	if( state[ind] == END_OF_LINE )
	{
	    tpop();
	    
	    /* If we have no token continue, otherw. */
	    /* return with token                     */

	    if( ind >= 0 || code == ONE_LINE )
		tpush(DONE);
	    else                                     /* Continue */
	    {
		sp = get_line( code );

		if( !sp || *sp == '\0' )
		    tpush(DONE);
		else if( *sp == '\n' )
		    tpush(END_OF_LINE);
		else if( isspace(*sp) )
		    tpush(SKIPPING_WS);
		else if( *sp == escape_char && isspace(sp[1]) )
		{
		    tpush(SKIPPING_WS);
		    tpush(PROCESSING_ESCAPE);
		}
		else if( *sp == comment_char )
		    tpush(SKIPPING_COMMENT);
		else if( *sp == var_ref_char )
		{
		    tpush(PROCESSING_TOK);
		    tpush(PROCESSING_REFERENCE);
		}
		else if( *sp == quote_char )
		{
		    sp++;
		    tpush(PROCESSING_TOK);
		    tpush(PROCESSING_STRING);
		}
		else
		    tpush(PROCESSING_TOK);
	    }
	}
	else if( state[ind] == SKIPPING_WS )
	{
	    tpop();
	    
	    while( *sp && *sp != '\n' && isspace(*sp) )
		sp++;

	    if( *sp == '\0' || *sp == '\n' )
		tpush(END_OF_LINE);
	    else if( *sp == comment_char )
		tpush(SKIPPING_COMMENT);
	    else if( *sp == escape_char && isspace(sp[1]) )
	    {
		tpush(SKIPPING_WS);
		tpush(PROCESSING_ESCAPE);
	    }
	    else if( *sp == var_ref_char )
	    {
		tpush(PROCESSING_TOK);
		tpush(PROCESSING_REFERENCE);
	    }
	    else if( *sp == quote_char )
	    {
		sp++;
		tpush(PROCESSING_TOK);
		tpush(PROCESSING_STRING);
	    }
	    else
		tpush(PROCESSING_TOK);
	}
	else if( state[ind] == SKIPPING_COMMENT )
	{
	    tpop();
	    
	    while( *sp && *sp != '\n' )
		sp++;

	    tpush(END_OF_LINE);
	}
	else if( state[ind] == PROCESSING_REFERENCE )
	{
	    tpop();
	    
	    if( sp[1] == command_open )
	    {
		sp += 2;

		for( q = sp; *q && *q != '\n' && *q != command_close; )
		    q++;

		if( *q == command_close )
		{
		    *q = '\0';

		    hp = lookup(sp);

		    if( hp )
		    {
			sp = hp->value;
			j = strlen(sp);

                        bp = buf_add_string( bp, sp, j, bufp, buflenp );
			sp += j + 1;
                        
			for( i = 1; i < hp->count; i++ )
			{
			    if( state[ind] == PROCESSING_STRING )
				bp = buf_add_char( bp, ' ', bufp, buflenp );
			    else
				bp = buf_inc_pntr( bp, bufp, buflenp );     /* ATTN: What is this for? */
			    
			    j = strlen(sp);

                            bp = buf_add_string( bp, sp, j, bufp, buflenp );
                            sp += j + 1;
			}
		    }

		    sp = q + 1;
		}
		else
		{
		    set_error( fp_stackind, "Runaway variable reference, no closing >" );
		    tpush(END_OF_LINE);

		    sp = q;
		}
	    }
	    else
                bp = buf_add_char( bp, *sp++, bufp, buflenp );
	}
	else if( state[ind] == PROCESSING_ESCAPE )
	{
	    tpop();
	    sp++;               /* Skip Initial Escape Character */

	    if( *sp && !isspace(*sp) )
	    {
		switch( *sp )
		{
		  case 'n':

		    *sp = '\n';
		    break;
		
		  case 't':

		    *sp = '\t';
		    break;
		
		  case 'v':

		    *sp = '\v';
		    break;
		
		  case 'b':

		    *sp = '\b';
		    break;
		
		  case 'r':

		    *sp = '\r';
		    break;
		
		  case 'f':

		    *sp = '\f';
		    break;
		
		  case 'a':

		    *sp = '\a';
		    break;
		}

                bp = buf_add_char( bp, *sp++, bufp, buflenp ); /* Put it in buffer to avoid downstream processing */
	    }
	}
	else if( state[ind] == PROCESSING_STRING )
	{
	    for( ;; )
	    {
		if( *sp == '\0' || *sp == '\n' )
		{
		    tpop();
		    tpush(END_OF_LINE);
		    set_error( fp_stackind, "Runaway string, missing \" inserted\n" );
		    break;
		}
		else if( *sp == var_ref_char )
		{
		    tpush(PROCESSING_REFERENCE);
		    break;
		}
		else if( *sp == quote_char )
		{
		    sp++;
		    tpop();
		    break;
		}
		else if( *sp == escape_char && (sp[1] == var_ref_char || sp[1] == quote_char) )
		{
		    sp++;
                    bp = buf_add_char( bp, *sp++, bufp, buflenp );
		}
		else if( *sp == escape_char )
		{
		    tpush(PROCESSING_ESCAPE);
		    break;
		}
		else
                    bp = buf_add_char( bp, *sp++, bufp, buflenp );
	    }
	}
	else
	{
	    for( ;; )
	    {
		if( *sp == '\0' || *sp == '\n' )
		{
		    tpush(END_OF_LINE);
		    break;
		}
		else if( isspace( *sp ) || *sp == comment_char )
		{
		    tpop();
		    tpush(DONE);
		    break;
		}
		else if( *sp == var_ref_char )
		{
		    tpush(PROCESSING_REFERENCE);
		    break;
		}
		else if( *sp == quote_char )
		{
		    sp++;
		    tpush(PROCESSING_STRING);
		    break;
		}
		else if( *sp == escape_char &&
			 (sp[1] == var_ref_char || sp[1] == quote_char || *sp == comment_char) )
		{
		    sp++;
                    bp = buf_add_char( bp, *sp++, bufp, buflenp );
		}
		else if( *sp == escape_char )
		{
		    tpush(PROCESSING_ESCAPE);
		    break;
		}
		else
                    bp = buf_add_char( bp, *sp++, bufp, buflenp );
	    }
	}
    }

    len = bp - *pbp;
    
    if( len )           /* If we've gotten something, pad it with a null byte */
    {
        bp = buf_add_char( bp, '\0', bufp, buflenp );
	len++;
    }
    
    *psp = sp;
    *pbp = bp;

    return( len );
}





/*
 * Local Variables:        
 * mode: c
 * eval: (make-local-variable 'compile-command)
 * compile-command: "make -k inpfproc.o"
 * End:
 *
 */

