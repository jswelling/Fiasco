/************************************************************
 *                                                          *
 *  stdcrg.h                                                *
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
 *  This is intended as a set of standard definitions       *
 *  for general use.  At the moment, it is not fully        *
 *  organized.                                              *
 *                                                          *
 *  Original Author:  Chris Genovese                        *
 *  Last Modified:    December 1995                         *
 *  Modified By:      Chris Genovese                        *
 *                                                          *
 ************************************************************/

#if !defined(DIAGNOSTIC)
#    define  DIAGNOSTIC   0
#endif

#if !defined(DEBUG)
#    define  DEBUG        0
#endif


/************************************************************
 *                                                          *
 *                  Error Functions                         *
 *                                                          *
 ************************************************************/
 

void       Stdcrg_Error( char *, ... );
void       Stdcrg_Warning( int, char *, ... );
void       Stdcrg_Internal( char *, int, char *, ... );
void       Stdcrg_SysError( char *, int, char *, ... );
void       Stdcrg_Usage( char *, ... );
void       Stdcrg_Message( char *, ... );

void       efunc_exit( int );
int        efunc_config( char *, ... );

FILE      *efopen( char *, char * );
int        efclose( FILE * );

size_t     efread( void *, size_t, size_t, FILE * );
size_t     efwrite( void *, size_t, size_t, FILE * );

void      *emalloc( size_t );
void      *ecalloc( size_t, size_t );
void      *erealloc( void *, size_t );
void     **matrix( long, long, size_t );


#define    HEREIAM      __FILE__,__LINE__
#define    NOLINE       -1


    /* More user-friendly interfaces to allocation routines */
    /*   Key:  n = number, t = type, p = pointer            */

#define    Malloc(n,t)       ((t*) emalloc( (n) * sizeof(t) ))
#define    Calloc(n,t)       ((t*) ecalloc( (n),  sizeof(t) ))
#define    Realloc(p,n,t)    ((t*) erealloc( (p), (n) * sizeof(t) ))
#define    Free(p)           (free(p))

#define    Matrix(n,p,t)     ((t **)matrix( (n), (p), sizeof(t) ))
#ifndef FreeMatrix
#define FreeMatrix(p) { free(*(p)); free((p)); }
#endif

/************************************************************
 *                                                          *
 *                  Help System                             *
 *                                                          *
 ************************************************************/
 
int      testHelp( int* argc, char** argv ); /* convenience function */
void     Help( char * );        /* Entry port to help system: interactive or lookup. */
void     Help_init(void);       /* Generated function that initializes the help structure */
void     Help_listtopics(void); /* Generated function to give pretty list of topics       */
void       Help_interp( char * );  /* core of Help(), but without table initialization */


#define    HELP_MAX_TOPICS       1024
#define    HELP_DEFAULT_PAGER    "more"
#define    HELP_COMPLETION_CHAR  '*'


/************************************************************
 *                                                          *
 *                  Command Line Options                    *
 *                                                          *
 ************************************************************/
 
void     cl_scan( int argc, char **argv );                    /* Read Command Line and Prepare   */
int      cl_get( char *name, char *fmt, ... );                /* Get Argument or Flag Values     */
void     cl_configure( char *name, char *val );                  /* Set default interpretations     */
/*
void     cl_unmatched( void );                                // Find unmatched args- not implemented //
int      cl_grep( char *pattern, char **arg );                // Find arg beginning with pattern //
*/
void     cl_mark( char *prefix, char *tag, char *showtime );  /* Echo Information to Log File    */
int      cl_cmd_line( int *, char *** );                      /* Copy command line structure     */
void     cl_cleanup( void );                                  /* Clean-up Internal CL Info       */
int      cl_present( char* flag_op );                         /* Check for presence of flag-type operand */
int      cl_cleanup_check( void );                            /* Clean-up CL Info, w/ error code */

/************************************************************
 *                                                          *
 *                  Input File Processing                   *
 *                                                          *
 ************************************************************/
 
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



/************************************************************
 *                                                          *
 *                  Binary Archive Management               *
 *                                                          *
 ************************************************************/



 


/************************************************************
 *                                                          *
 *                  Sorting and Searching                   *
 *                                                          *
 ************************************************************/


    /* Comparison functions for qsort etc. */

int  Numerically_d( const void *, const void * );
int  Numerically_f( const void *, const void * );
int  Numerically_l( const void *, const void * );
int  Numerically_i( const void *, const void * );
int  Numerically_s( const void *, const void * );
int  Alphabetically( const void *, const void * );




/************************************************************
 *                                                          *
 *                  Mathematical Operations                 *
 *                                                          *
 ************************************************************/
 

    /* Max and Min  */

#define     Max(a,b)       (((a) < (b)) ? (b) : (a))
#define     Min(a,b)       (((a) > (b)) ? (b) : (a))

#define     dMax(a,b)      (0.5*((a)+(b) + fabs((b) - (a))))
#define     dMin(a,b)      (0.5*((a)+(b) - fabs((b) - (a))))

#ifdef __GNUC__
#define     fMax(a,b)      (0.5*((a)+(b) + (float)fabs((b) - (a))))
#define     fMin(a,b)      (0.5*((a)+(b) - (float)fabs((b) - (a))))
#else
#define     fMax(a,b)      (0.5*((a)+(b) + fabsf((b) - (a))))
#define     fMin(a,b)      (0.5*((a)+(b) - fabsf((b) - (a))))
#endif

#define     Maxabs(a,b)    ((fabs((double)(a)) < fabs((double)(b))) ? (b) : (a))
#define     Minabs(a,b)    ((fabs((double)(a)) > fabs((double)(b))) ? (b) : (a))

    /* Small Powers */

#define     Square(a)      ((a)*(a))
#define     Cube(a)        ((a)*(a)*(a))
#define     Fourth(a)      ((a)*(a)*(a)*(a))
#define     Fifth(a)       ((a)*(a)*(a)*(a)*(a))

    /* Common Mathematical Functions -- work with both float and double */

#define     Log(x)         (log((double)(x)))
#define     Exp(x)         (exp((double)(x)))  
#define     Cos(x)         (cos((double)(x)))
#define     Sin(x)         (sin((double)(x)))
#define     Sqrt(x)        (sqrt((double)(x)))

 
    /* Important Constants  */

#define  M_SQRT2PI        2.5066282746310005024

#ifdef never
    /* For passing to Fortran routines */

static  int        iOne = 1; 
static  int        iZero = 0;
static  int        imOne = -1;

static  float      fOne = 1.0;
static  float      fZero = 0.0;
static  float      fmOne = -1.0;

static  double     One = 1.0;
static  double     Zero = 0.0;
static  double     mOne = -1.0;
#endif



/************************************************************
 *                                                          *
 *                  Miscellaneous                           *
 *                                                          *
 ************************************************************/
 
#define FILE_MAX   256      /* Replace with enum limits containing various limits */


