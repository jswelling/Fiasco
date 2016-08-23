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

static char rcsid[] = "$Id: vmpfx.c,v 1.26 2007/03/21 23:45:49 welling Exp $";

/*
 *  vmpfx.c
 * 
 *  Voxel-wise, fixed effects, maximum posterior estimation
 *  using basic time-course model.  As this develops, it will
 *  incorporate dips and other noise models (e.g., not just
 *  iid N(0,\sigma^2)).
 *
 *  Original Author:    Christopher Genovese
 *  Last Modification:  May 1996 (Created)
 *  Modified By:        Christopher Genovese
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include  <math.h>
#include  "mri.h"        /* Pittsburgh MRI format definition */
#include  "fmri.h"       /* More Pittsburgh MRI stuff        */
#include  "misc.h"
#include  "stdcrg.h"
#include  "minqnb.h"     /* quasi-newton bounded optimization */
#include  "tmperror.h"   /* temporary version of error streams */

#ifdef DARWIN
#define finite( foo ) isfinite( foo )
#endif

#if defined(NOT_FIASCO)
#   undef Error             
#   undef Abort
#   define Abort Error
#endif /* Not Fiasco compilation */

#include  "stdcrg.h"

/*
 * Old Pittsburgh format is being phased out, and this clunky interface
 * with it.  For now, "old-pgh" format uses these functions.
 *
 */

#if defined(USE_MRIFUNC)
#   include  "mrifunc.h"
#else
    enum mri_types { CHAR_T, LONG_T, FLOAT_T, DOUBLE_T, SHORT_T, VOID_T, RECORD_T };
#endif

#include  <signal.h>
#include  <time.h>
#include  <errno.h>
#include  <sys/stat.h>
#include  <assert.h>

#include  "lapack.h"

#if defined(USE_IMSL)
#   include  "imsl.h"
#endif

/*** Macros and Constants      ***/

#if !defined(SMOOTH_BELLS)
#    define  SMOOTH_BELLS     1
#endif

enum shape_pars { LAG_ON, ATTACK, LAG_OFF, DECAY, RISE, FALL, DIP_HGT, DIP_SKEW, TOTAL_SHAPE_PARAMS };

enum drift_lims { DEGREE_MAX = 4, KNOT_MAX = 4096 };
enum field_lims { FIELD_DIM = 4, MAX_DIAGNOSTICS = 32 };

enum files      { INPUT, DATA, INIT, OUTPUT, BINARY, FILES };

enum products   
{
    ESTIMATES, STANDARD_ERRORS, COVARIANCES,
    RESIDUALS, FITTED_VALUES, DIAGNOSTICS,
    NULL_ESTIMATES, NULL_STANDARD_ERRORS, NULL_COVARIANCES,
    NULL_RESIDUALS, NULL_FITTED_VALUES, NULL_DIAGNOSTICS,
    NUM_PRODUCTS, SET_IT, WRITE_IT
};

enum model_specs { NULL_MODEL=1, MAXIMAL_MODELS=2, FULL_MODEL=4, ALL_MODELS=8, NUM_MODEL_SPECS=4 };
enum model_probs { NULL_WEIGHTED, INDEPENDENT, ARBITRARY, NUM_MODEL_PROBS };

enum inits
{
    BASELINE, DRIFT_COEFS, DRIFT_KNOTS, SHAPE, RESP, O_SIGMASQ, NUMINITS
};

enum covariance_methods { COV_SCD_DIRECT=0, COV_SCD_RICHARDSON=1, COV_FCD_DIRECT=2, COV_FCD_RICHARDSON=3,
                          COV_ANALYTIC=4, COV_FCD_MASK=2 }; /* ATTN: CHANGE MASKS AND CODES */

static ErrorStream errorStream;

static double const Tolerance  = 1.0e-8;
static double const MTolerance = 1.0e-8;
static double const DTolerance = 1.0e-5;
static double const STolerance = 1.0e-4;

static double const Dscale     = 1.0e-3;
static double const EpsLog     = 1.0e-256;  /* Small number for Log near 0 */
static double const BigNumber  = 1.0e256;

static double const LN_SQRT_2PI = 0.9189385332046726;

#define FILE_MAX                 256
#define DEFAULT_CHUNKSIZE        ((unsigned long)(1 << 30))


static const char *Type_String[] = { "char", "long", "float", "double", "short", "void", "record", "" };

    /*
     * The End_Hdr_String is placed at the end of the output file header to distinguish the
     * content that has to be updated with each module.  This is a hack which will be improved
     * on later.
     *
     */

static const char  End_Hdr_String[] = "####################END HEADER####################\n";

    /*
     * Constants for Fortran Functions
     *
     *  These are not declared const because their address needs to be passed
     *  to call-by-reference Fortran functions, but they are used as constants
     *  and are not to be changed.
     *
     */

#ifndef USE_FORTRAN_CONSTS

static  int        iOne = 1;
static  int        iZero = 0;
static  int        imOne = -1;

static  float      fOne = 1.0;
static  float      fZero = 0.0;
static  float      fmOne = -1.0;

static  double     One = 1.0;
static  double     Zero = 0.0;
static  double     mOne = -1.0;

#endif /* USE_FORTRAN_CONSTS */

/*** Types and Data Structures ***/

typedef struct output_struct
{
    char          *buf;
    long           bufsiz;
    long           offset;
    FILE          *stream;
} CookieJar;


    /*
     * A structure to manage access to various input file types.  This is
     * rather clumsy because I haven't gone back and absorbed old init_data
     * param stuff into this structure.
     *       
     *       contributed by Joel Welling
     *
     */

typedef struct data_accessor_struct
{
  void   (*setup)(struct data_accessor_struct* this,
                  FILE* fp, long n_vox_skip, int times_per_voxel, int type_size);
  void   (*read)(struct data_accessor_struct* this,
                 void* buf, FILE* fp, int times_per_voxel, int type_size);
  void   (*close)(struct data_accessor_struct* this, FILE* fp);
  void*    state;
} DataAccessor;

typedef struct mri_accessor_state_struct
{
  MRI_Dataset* ds;
  long         offset;
} MriDataState;

/*** Function Prototypes ***/

    /* Initialization */

void    init_fixed( int, int * );
void    init_params( void );
FILE   *init_init( int *, char ** );
void    init_init2( void );
FILE   *init_data( char *, char *, long *, long *, int *, int *, int **, DataAccessor * );
void    init_design( void );
void    init_products( int *, char ***, int **, int **, long **, long **, long **, long *, long );
void    init_basis( void );

void    init_par_indices( void );
void    init_bounds( void );
void    init_obj_scale( int, double *, double *, int *, double *, int );
void    init_standard_errors( double *, int *, double *, int );
void    init_detrend( int, int, int, double *, int );
void    init_binom_coefs( void );

    /* Voxelwise Fitting */

int    voxel_init( int, double *, int, int, int, double *, double *, double *, int );
void    voxel_fit( int, double *, int, int, double *, double *, int *, int *, double *,
                   double , int *, double *, int, int, void (*)( int *, double *, double * ),
                   void (*)( int *, double *, double * ) );

void    set_inits( int, double *, int, int, double *, int *, double *, double *, int );

void    parse_model( char *, long * );

    /* Posterior Computation */

void    lnpost_fixed( int *npar, double *p, double *val );
void    Dlnpost_fixed( int *npar, double *p, double *deriv );
void    DDlnpost_fixed( int neff, int n, double *p, double *obsinfo, int *skip, double *work, int lwork );

void    lnpost_variable( int *npar, double *p, double *val );
void    Dlnpost_variable( int *npar, double *p, double *deriv );
void    lnpost_null_fixed( int *npar, double *p, double *val );
void    Dlnpost_null_fixed( int *npar, double *p, double *deriv );
void    DDlnpost_null_fixed( int, int, double *, double *, int *, double *, int );
void    lnpost_null_variable( int *npar, double *p, double *val );
void    Dlnpost_null_variable( int *npar, double *p, double *deriv );

void    lnpost_shape( int *, double *, double * );
void    Dlnpost_shape( int *, double *, double * );
void    set_lnpshape( double *, double * );


void   profile_active( int, double *, double *, int, double *, int * );
void   profile_active2( int, double *, double *, int, double *, int * );

#if !defined(NO_INLINE_DRIFTPROF)
#   define  profile_drift( nt, d, basis, coefs, prof, work, lwork ) \
    dgemv( &NoTrans, &(nt), &(d), &One, (basis)[0], &(nt), (coefs), &iOne, &Zero, (prof), &iOne )
#else
void    profile_drift( int, int, double **, double *, double *, double *, int );
#endif

void    profiles( double *p, double *dr_prof, double *ac_prof, double *work, int lwork ); /* DEFUNCT */

void    poly_bell4( long, double *, int *, double *, double, double, double, double );
void    poly_bell6( long, double *, int *, double *, double, double, double, double );
void    poly_bell8( long, double *, int *, double *, double, double, double, double );

void    shape_deriv_4( double *, double **, double *, double, double, double * );
void    shape_deriv_6( double *, double **, double *, double, double, double * );
void    shape_deriv_8( double *, double **, double *, double, double, double * );

void    make_active_mats( int, int, double *, double *, double *, double *,
                            double *, double *, double *, int * );

#if defined(NON_ADDITIVE_PASTE)
double  smooth_paste( double x, double y, double s );
void    Dsmooth_paste( double x, double y, double s, double *grad );
#endif

    /* Covariance Computation */

void    compute_cov( int neff, int n, double *p, double *obsinfo, int method,
                     int *skip, double *work, int worksize, int *rank, double *lndet, int *info,
                     int null_model );

void    hessian_scd( int n, int nbig, double *p, double *plb, double *pub, double *I, int method,
                     double *hinit, int *skip, double maxmin, double errf, double tol,
                     double *work, int worksize, void (*func)(int *, double *, double *), int *info );

void    hessian_fcd( int n, int nbig, double *p, double *plb, double *pub, double *I, int method,
                     double *hinit, int *skip, double maxmin, double errf, double tol,
                     double *work, int worksize, void (*func)(int *, double *, double *), int *info );

void    reduce_hessian( int neff, double *obsinfo, int *skip, int firstskip );
void    expand_hessian( int neff, int ns, double *obsinfo,
                        int *skip, double *diag, int firstskip, double *work );

    /* Basis Management */

void    make_drift_basis( int, int, int, double **, double *, double *, double *, int );
void    make_drift_qforms( int, int, int, double *, double *, double *, double *,
                           double *, double *, int );

    /* Optimization */
    /* defined in minqnb.h */

    /* Data Access */

void    defaultDataSetup( DataAccessor*, FILE*, long, int, int );
void    defaultDataRead( DataAccessor* this, void* buf, FILE* fp, int times_per_voxel, int type_size );
void    defaultDataClose( DataAccessor* this, FILE* fp );
void    mriDataSetup( DataAccessor* this, FILE* fp, long n_vox_skip, int times_per_voxel, int type_size );
void    mriDataRead( DataAccessor* this, void* buf, FILE* fp, int times_per_voxel, int type_size );
void    mriDataClose(DataAccessor* this, FILE* fp);

    /* Output Routines */

void    output_cookie( int, ... );
void    bake_cookie( CookieJar*, void *, size_t );
void    eat_cookies( CookieJar* );

static void  output_null_record( const char* mesg, int voxel, long type,
                                 int nout, int tdim, int null_model, int nnul,
                                 double* pars, double* timec, double* imat );

    /* Signal handlers   */

void    set_current_voxel( int v );
void    handler( int sig );


/*** Global Variables ***/

int       V;               /* Number of total data voxels   */
int       Veff;            /* Number of voxels processed    */
int       T;               /* Number of Images              */
int       B;               /* Number of Blocks              */
int       K;               /* Number of Conditions          */
int       Keff;            /* Effective # of Conditions     */
int       D;               /* Total Size of Drift Basis     */
int       Dg;              /* Polynomial Degree of Drift    */
int       Dk;              /* # interior knots, drift term  */

int       NParams;         /* Total # of Parameters         */
int       NeffParams;      /* Effective # of Parameters     */
int       NoutParams;      /* # of Parameters Output        */
int       NnulParams;      /* # of Parameters in Null Model */
int       UniqueStims;     /* # of Unique Stimulus Lengths  */
int       Shape_Params;    /* # Response Curve Parameters   */

double   *Param;           /* Parameter Vector              */
double   *Param_LB;        /* Lower bounds on parameters    */
double   *Param_UB;        /* Upper bounds on parameters    */
double   *saveParam;       /* Scratch Parameter Vector      */

int       Iter_Max;        /* Maximum Iterations per Voxel  */
int       Eval_Max;        /* Maximum Evaluations per Voxel */
double    Pscale;          /* Scaling for Pasting Function  */

int      *BlockStart;      /* Index -> Block Beginning Map  */
int      *BlockEnd;        /* Index -> Block Ending Map     */
int      *WhichStim;       /* Block -> Unique Stim Mapping  */
int      *Cond;            /* Block -> Condition Mapping    */
int      *Constrained;     /* 0: Free, 1= Resp Fixed at 0.0 */
int      *CondMap;         /* T vector condition by image   */
int      *CondImages;      /* # images in each condition    */
int     **CondMatrix;      /* K x T condition dummy matrix  */
int     **CondBlock;
double   *Stims;           /* Unique Stimulus Lengths       */
double   *StimulusStart;   /* Block -> Stimulus Beginning   */
double   *StimulusLen;     /* Block -> Stimulus Length      */
double   *StimShift;       /* Shift for Unaligned Blocks    */

double   *Bells;           /* Bell Functions for each stim  */
double  **DBells;          /* dBell/dShape functions        */
int      *Bell_Length;     /* Length of Bell Support        */

int       RampLen;         /* Length of Spline Ramps        */
double   *AttackRamp;      /* Spline Ramps for making Bells */
double   *DecayRamp;
double   *RiseRamp;
double   *FallRamp;

double   *DAttackRamp;     /* Derivatives of Spline Ramps   */
double   *DDecayRamp;
double   *DRiseRamp;
double   *DFallRamp;

int      *Image_Dims=NULL; /* Spatial Dimensions of Image */

    /* Drift Basis Management */

double  *Pnorm;
double  *Pcurv;
double  *Pcomb;
double  *Ptot;
double  *basisR;
double  *basisR_inv;
double  **Basis;           /* Drift Basis Functions         */


    /* Data and Intermediary  */

double    TR;              /* Repeat Time (Input)           */
double    SSR;             /* Residual sum of squares       */

double   *Y;               /* Data for current voxel        */
double   *R;               /* Residual for current voxel    */
double   *Rsave;

double   *Drift_prof;      /* Drift Terms, current voxel    */ 
double   *Active_prof;     /* Activation Profile, cur. vox. */

int       Knots_Fixed;     /* Indicator that knots fixed    */
double   *Fixed_Knots;     /* Values for fixed knots        */

double   *Work;            /* Work Space                    */
int       WorkSize;        /* Size of Work Space            */
int      *IWork;           /* Integer Work Space            */
int       IWorkSize;       /* Size of Integer Work Space    */

double   *Save = NULL;     /* Scratch buffer                */

int       minimizerWorkSize; /* Scratch space size specific to minimizer */ 
double   *minimizerWork;     /* Scratch space specific to minimizer */
int       minimizerIWorkSize; /*int scratch space size specific to minimizer*/
int      *minimizerIWork;    /* int scratch space specific to minimizer */

    /* Standard Error Compt'n */

int       SEworksize;

double    Scd_Tolerance;

double   *Hinit;           /* Initial Values for Hessian    */
double   *ObsInfo;         /* Obs. Info Matrix (analogue)   */
double   *SEpsc;           /* Initial Richardson Scale      */
double   *SEwork;

    /* Parameter Index Arrays */

int       Mu;              /* Baseline         */
int       o_SigmaSq;       /* Dispersion       */
int      *Drift;           /* Drift            */
int      *Coefs;               /* Coefficients */
int      *Knots;               /* Knot Posit'n */
int      *Resp;            /* Responsiveness   */
int      *Shape;           /* Shape            */

int      *RespEff;         /* Eff'v Resp Pars  */
int       RespInd;         /* First Eff'v Resp */
int      *Shape_Use;       /* Which Par's used */

    /* Prior Hyper-Parameters */

double    Mu_a, Mu_b;                  /* Baseline: t_1 Center, Spread */
double    Drift_a;                     /* Drift Curvature Penalty      */
double    Drift_b;                     /* Knot Prior Unif Order Stat   */
double    Drift_c;                     /* df for drift smooth (0=free) */
double    SmoothPar = 1.0;             /* drift smooth par. (1/lambda) */
double    Resp_a, Resp_b;              /* Responsiveness: Unif-Exp     */
double    *Shape_a;                    /* Shape:  Gamma distribution   */
double    *Shape_b;
double    SigmaSq_a, SigmaSq_b;        /* Variance: LogNormal          */
double    o_SigmaSq_bsq;

double    Smooth_Target = 1.0;         /* Target coefficient of P      */
double    Drift_df;                    /* Effective df in all cases    */

    /* Initialization Backfitting */

double   *Xty;
double   *XtS0;
double   *S0ty;
double   *StS;
double   *St1;
double   *Sty;
double   *S;


    /* Posterior Odds Computation */

int       Null_Model = 1;

    /* Input File Parsing */

int       Design_Images;
int       Design_Aligned;

double    IAI;              /* Image Acquisition Interval (seconds) */

int       Num_Slices;
int      *Slice_Order;
double    Slice_Offset;

char      Drift_scaled_str[32];
int       Drift_scaled = 0;

int       Skip = 0;

int       Num_Fixed = 0;
int      *Fixed;

int       Num_Init = 0;
int       Init_len;
int       Init_size;
int       Init_type;
char     *Init_Type_String;
int      *Init;
int      *Init_offsets;
char    **Init_List;
int       Init_Max_Iter;
double    Init_Tolerance;

int       Num_Products;
int      *Produce;
int      *Product_Type;
char    **Products;

long      Record_Size;
long     *Units;
long     *UnitSize;
long     *Type;

int       Num_Contents = 0;
char    **Contents;
unsigned long  *Contents_Chunk;
unsigned long  *Contents_Address;

int       Num_QNB_Int_Opts = 0;
char    **QNB_Int_Options_str;
double   *QNB_Int_Options_val;
int       Num_QNB_Double_Opts = 0;
char    **QNB_Double_Options_str;
double   *QNB_Double_Options_val;

int       Output_exists = -1;
int       Append_out;
int      *Out;
char     *Append;
char     *Output_fn;
char     *Output_ext;
char     *Res_fn = NULL;
char     *Log_fn = NULL;

unsigned long      Output_chunksize;

    /* External */

extern    int   errno;

    /* Configurable Functions */

void    (*lnpost)( int *, double *, double * );
void    (*Dlnpost)( int *, double *, double * );
void    (*DDlnpost)( int, int, double *, double *, int *, double *, int );
void    (*poly_bell)( long, double *, int *, double *, double, double, double, double );
void    (*shape_deriv)( double *, double **, double *, double, double, double * );

void    (*lnpost_null)( int *, double *, double * );
void    (*Dlnpost_null)( int *, double *, double * );
void    (*DDlnpost_null)( int, int, double *, double *, int *, double *, int );


#if DIAGNOSTIC >= 20
FILE         *debugfp = NULL; /* ATTN: TEMP */
#endif
    

/*** Source Code ***/

#ifndef USE_MRIFUNC
static size_t mriTypeSize( long type )
{
    switch( type )
    {
      case CHAR_T:

        return( sizeof(char) );

      case SHORT_T:

        return( sizeof(short) );

      case LONG_T:

        return( sizeof(long) );

      case FLOAT_T:

        return( sizeof(float) );

      case DOUBLE_T:

        return( sizeof(double) );

      case VOID_T:
      case RECORD_T:

	fprintf(stderr,"mriTypeSize: sizeof(void) requested!\n");
	exit(-1);
	return 0; /* not reached */
    }
    return 0; /* not reached */
}

/*
 * General conversion between two vectors of arbitrary types and the same length.
 * Usual C conventions used in carrying out the conversion.
 *
 * Here,  y <- x.
 *
 */

void   mriTypeConvert( int T, int type_y, void *y, int type_x, void *x )
{
    int     k;
    
    switch( type_x ) 
    {
      case SHORT_T:

        switch( type_y )
        {
          case SHORT_T:

            for( k = 0; k < T; k++ )
                ((short *)y)[k] = ((short*)x)[k];

            break;

          case FLOAT_T:

            for( k = 0; k < T; k++ )
                ((float *)y)[k] = ((short*)x)[k];

            break;

          case LONG_T:

            for( k = 0; k < T; k++ )
                ((long *)y)[k] = ((short*)x)[k];

            break;

          case CHAR_T:

            for( k = 0; k < T; k++ )
                ((char *)y)[k] = ((short*)x)[k];

            break;

          case DOUBLE_T:

            for( k = 0; k < T; k++ )
                ((double *)y)[k] = ((short*)x)[k];

            break;
        }
        
        break;

      case FLOAT_T:

        switch( type_y )
        {
          case SHORT_T:

            for( k = 0; k < T; k++ )
                ((short *)y)[k] = ((float *)x)[k];

            break;

          case FLOAT_T:

            for( k = 0; k < T; k++ )
                ((float *)y)[k] = ((float *)x)[k];

            break;

          case LONG_T:

            for( k = 0; k < T; k++ )
                ((long *)y)[k] = ((float *)x)[k];

            break;

          case CHAR_T:

            for( k = 0; k < T; k++ )
                ((char *)y)[k] = ((float *)x)[k];

            break;

          case DOUBLE_T:

            for( k = 0; k < T; k++ )
                ((double *)y)[k] = ((float *)x)[k];

            break;
        }
        
        break;
                    
      case LONG_T:

        switch( type_y )
        {
          case SHORT_T:

            for( k = 0; k < T; k++ )
                ((short *)y)[k] = ((long*)x)[k];

            break;

          case FLOAT_T:

            for( k = 0; k < T; k++ )
                ((float *)y)[k] = ((long*)x)[k];

            break;

          case LONG_T:

            for( k = 0; k < T; k++ )
                ((long *)y)[k] = ((long*)x)[k];

            break;

          case CHAR_T:

            for( k = 0; k < T; k++ )
                ((char *)y)[k] = ((long*)x)[k];

            break;

          case DOUBLE_T:

            for( k = 0; k < T; k++ )
                ((double *)y)[k] = ((long*)x)[k];

            break;
        }
        
        break;

      case CHAR_T:

        switch( type_y )
        {
          case SHORT_T:

            for( k = 0; k < T; k++ )
                ((short *)y)[k] = ((char *)x)[k];

            break;

          case FLOAT_T:

            for( k = 0; k < T; k++ )
                ((float *)y)[k] = ((char *)x)[k];

            break;

          case LONG_T:

            for( k = 0; k < T; k++ )
                ((long *)y)[k] = ((char *)x)[k];

            break;

          case CHAR_T:

            for( k = 0; k < T; k++ )
                ((char *)y)[k] = ((char *)x)[k];

            break;

          case DOUBLE_T:

            for( k = 0; k < T; k++ )
                ((double *)y)[k] = ((char *)x)[k];

            break;
        }
        
        break;

      case DOUBLE_T:

        switch( type_y )
        {
          case SHORT_T:

            for( k = 0; k < T; k++ )
                ((short *)y)[k] = ((double*)x)[k];

            break;

          case FLOAT_T:

            for( k = 0; k < T; k++ )
                ((float *)y)[k] = ((double*)x)[k];

            break;

          case LONG_T:

            for( k = 0; k < T; k++ )
                ((long *)y)[k] = ((double*)x)[k];

            break;

          case CHAR_T:

            for( k = 0; k < T; k++ )
                ((char *)y)[k] = ((double*)x)[k];

            break;

          case DOUBLE_T:

            for( k = 0; k < T; k++ )
                ((double *)y)[k] = ((double*)x)[k];

            break;
        }
        
        break;
    }
}

#endif

int main( int argc, char **argv )
{
    int      i, j, k, m, s, t, v;
    int      iter, nfunc, info = 0;
    int      cov_method;                /* Use Richardson Method for SE computation */
    int      fastopt = 1;               /* Use Quasi-Newton Optimization            */
    int      cov_rank, cov0_rank;
    int      first_voxel, last_voxel;

    int     *atbounds;
    
    size_t   size, size_r, tsz;

    long     data_hdr = 0;
    long     type = -1, type_r = -1;
    long     addr = 0;
    long     output_bufsize;
    long     memlimit;
    
    unsigned long  mpsize;

    int      idiagnostics[MAX_DIAGNOSTICS];
    float    rdiagnostics[MAX_DIAGNOSTICS];
    
    double   u, val, val0, fsc;
    double   cov_lndet, cov0_lndet;

    char     ifn[FILE_MAX];
    char     scratch[1024];
    char     *name, *path;

    long     model;
    char    *model_string;
    double   model_prior;

    char     *data_fn, *data_fmt;
    DataAccessor inputDataAccessor;

    void     *voxdat, *initbuf;
    char     *output_buf;

    double   *psc, *initdat;

    time_t   exectime;

    FILE     *fp[FILES];

    struct stat  statbuf;

    int init_retcode;

    /** Initialize error stream **/
    errorStream= es_new();

    /** Offer help if requested **/

    if( argc == 1 )
        Help( "Usage" );
    else if( !strcmp( argv[1], "-help" ) )
    {
        if( argc == 2 )
            Help( "selecttopic" );
        else
            Help( argv[2] );
    }

    time( &exectime );       /* Record current time */

    /** Process Command Line and Input File **/

    cl_scan( argc, argv );
    cl_get( "", "%s", ifn );              /* Input File Name */

    if( !inpf_scan( ifn ) )
    {
        Message( "Errors encountered in input file %s\n", ifn );
        inpf_perror( );
        exit( 0 );
    }

        /* Output Specification     */
    
    inpf_get( "output.file",      "%ls",       &Output_fn );
    inpf_get( "output.logfile",   "%ls",       &Log_fn );
    inpf_get( "output.append",    "%ls[true]", &Append );

    inpf_arg( "output.file", "output" );


            /* Prepare output files  */

    if( inpf_unmatched( "output.file" ) )
    {
        char template[]= "mpfxoutXXXXXX";
	(void)mkstemp( template );
	Output_fn = strdup(template);
        Output_exists = 0;
    }
    else
    {
        if( !stat( Output_fn, &statbuf ) )
            Output_exists = 1;
        else
            Output_exists = 0;
    }

    if( Output_exists )
        Abort( "%s currently set to avoid overwriting, %s exists\n", argv[0], Output_fn );

    Append_out = (!strcasecmp( Append, "true" ) && Output_exists);

    Output_ext = strrchr( Output_fn, '.' );    /* Find extension if any */

    if( Output_ext )
        *Output_ext = '\0';

    if( !Output_exists )               /* Default Results Storage */
    {
        Res_fn = Calloc( strlen(Output_fn) + 5, char ); 
        sprintf( Res_fn, "%s.bin", Output_fn );
    }
    
    if( inpf_unmatched( "output.logfile" ) )   /* Log File        */
    {
        Log_fn = Calloc( strlen(Output_fn) + 5, char ); 
        sprintf( Log_fn, "%s.log", Output_fn );
    }

    if( Output_ext )
        *Output_ext = '.';
    
        /* Record Input Configuration in Log File */

    if( strcmp( Log_fn, "-" ) )     /* If -, use stderr, otherwise use given name */
    {
        if( Append_out || !stat( Log_fn, &statbuf ) )
            efunc_config( "log file = >%s", Log_fn );
        else
            efunc_config( "log file = %s", Log_fn );
    }
    

    Message( "##\n## Module %s invoked at %s", argv[0], ctime( &exectime ) );
    Message( "## Command Line:   " );

    for( i = 0; i < argc; i++ )
        Message( "%s ", argv[i] );
    Message( "\n##\n\n\n" );

    Message( "##\n## Input Configuration:\n##\n\n" );
    fp[INPUT] = inpf_insert( ifn, stderr );         /* Insert input file for documentation */
    efclose( fp[INPUT] );

    Message( "\n\n##\n## Execution Log:\n##\n\n" );    /* Prepare to record messages          */


        /* Experiment Specification */

    inpf_get( "conditions",    "%d",                &K );
    inpf_get( "design",        "%lM[=:%d %lf %lf]", &B, &Cond, &StimulusStart, &StimulusLen );
    inpf_get( "design.images", "%d[1]",             &Design_Images );    
    inpf_get( "design.aligned","%d[1]",             &Design_Aligned );
    inpf_get( "IAI",           "%lf[1]",            &IAI );

    inpf_get( "fixed",         "%lV[=:%d]",         &Num_Fixed, &Fixed );
    inpf_get( "slice.offset",  "%lf[0]",            &Slice_Offset );
    inpf_get( "slice.order",   "%lV[=:%d]",         &Num_Slices, &Slice_Order );


        /* Parameter Specification  */

    inpf_arg( "drift.degree",  "degree" );

    inpf_get( "ramp.length",   "%d[0]", &RampLen );
    inpf_get( "drift.degree",  "%d[3]", &Dg ); 
    inpf_get( "skip",          "%d[0]", &Skip );
    inpf_get( "shape.params",  "%d[4]", &Shape_Params );


#if defined(NON_ADDITIVE_PASTE)

    inpf_get( "paste.form",    "%d[0]",           &Pform );  /* 0 additive, 1 smooth_paste */
    inpf_get( "paste.scale",   "%lf[0.001]",      &Pscale );

    if( Pform )
        Abort( "Non-additive paste is temporarily unavailable." );

#endif

        /* Prior Specification      */

    inpf_get( "prior.baseline",     "%lf %lf[1]",            &Mu_a, &Mu_b );
    inpf_get( "prior.drift",        "%lf[1] %lf[1] %lf[0]",  &Drift_a, &Drift_b, &Drift_c );
    inpf_get( "prior.resp",         "%lf[2.5e-3] %lf[50]",   &Resp_a, &Resp_b );
    inpf_get( "prior.shape",        "%lM[%:%lf[2] %lf[2]]",  TOTAL_SHAPE_PARAMS, &Shape_a, &Shape_b );
    inpf_get( "prior.sigmasq",      "%lf[2.5] %lf[1.25]",    &SigmaSq_a, &SigmaSq_b );

    inpf_get( "prior.model",        "%ls[null-weighted] %lf[0.0]", &model_string, &model_prior );

    for( s = 0; s < TOTAL_SHAPE_PARAMS; s++ )
        if( Shape_a[s] == 1.0 )
            Abort( "Currently, first parameter for Shape prior should be > 1." );  /* ATTN: Revise this */

    inpf_get( "prior.drift.scaled", "%s[false]",             Drift_scaled_str ); /* Deprecated */

    if( !strcmp( "true", Drift_scaled_str ) )     
        Drift_scaled = 1;
    else
        Drift_scaled = 0;

        /* Input Specification      */

    inpf_get( "data.file",      "%ls",       &data_fn );
    inpf_get( "data.format",    "%ls[list]", &data_fmt ); /* list, pgh, pgh-old, raw */

    inpf_arg( "data.file",      "data" );
    inpf_arg( "data.format",    "format" );

    if( !strcasecmp( data_fmt, "raw" ) )
    {
        inpf_get( "data.voxels",  "%d",         &V );
        inpf_get( "data.dim.t",   "%d",         &T );
        inpf_get( "data.type",    "%ld[%]",     FLOAT_T, &type );
        inpf_get( "data.header",  "%ld[0]",     &data_hdr );
        inpf_get( "data.missing", "%s[none]",   scratch );     /* Deal with this later */
        inpf_get( "data.image",   "%s[no]",     scratch );

        if( !strcasecmp( scratch, "yes" ) )
        {
            Image_Dims = Calloc( FIELD_DIM, int );
            Image_Dims[0] = T;
            
            inpf_get( "data.dim.x",  "%d\n",    &Image_Dims[1] );
            inpf_get( "data.dim.y",  "%d\n",    &Image_Dims[2] );
            inpf_get( "data.dim.z",  "%d\n",    &Image_Dims[3] );

            if( V != Image_Dims[1] * Image_Dims[2] * Image_Dims[3] )
                Warning( DIAGNOSTIC >= 1, "Apparent inconsistency between # of voxels (%d) and "
                         "product of spatial dimensions from input file." );
        }
        else
            Image_Dims = NULL;

        if( inpf_unmatched( "data.voxels" ) )
            Abort( "In raw data format, number of voxels (tag data.voxels) required in input file." );

        if( inpf_unmatched( "data.dim.t" ) )
            Abort( "In raw data format, number of time points (tag data.dim.t) required in input file." );
    }

        /* Algorithm Specification  */

    inpf_get( "vmpfx.iterations",  "%d[5000]",       &Iter_Max );
    inpf_get( "vmpfx.evaluations", "%d[%]",          10*Iter_Max, &Eval_Max );
    inpf_get( "vmpfx.memlimit",    "%d[%]",          1048576, &memlimit );     /* Output Buffer Size: 1 Mb default */

    inpf_get( "vmpfx.qnb.options.int",  "%lM[=:%ls %d]",  &Num_QNB_Int_Opts, 
	      &QNB_Int_Options_str, &QNB_Int_Options_val );
    inpf_get( "vmpfx.qnb.options.real", "%lM[=:%ls %lf]", 
	      &Num_QNB_Double_Opts, &QNB_Double_Options_str, 
	      &QNB_Double_Options_val );

    inpf_get( "vmpfx.scd.tol",     "%lf[%]",         DTolerance, &Scd_Tolerance );

    inpf_get( "vmpfx.cov.method",  "%s[analytic]",   scratch );

    if( !strcasecmp( scratch, "analytic" ) )
        cov_method = COV_ANALYTIC;
    else if( !strcasecmp( scratch, "scd.richardson" ) )
        cov_method = COV_SCD_RICHARDSON;
    else if( !strcasecmp( scratch, "scd.direct" ) )
        cov_method = COV_SCD_DIRECT;
    else if( !strcasecmp( scratch, "fcd.richardson" ) )
        cov_method = COV_FCD_RICHARDSON;
    else if( !strcasecmp( scratch, "fcd.direct" ) )
        cov_method = COV_FCD_DIRECT;
    
    if( cov_method & COV_FCD_MASK )
        Abort( "Gradient First Central Differences not yet supported for Hessian computation." );

    inpf_get( "vmpfx.knots",       "%d[6]",          &Dk );
    inpf_get( "vmpfx.fix.knots",   "%s[true]",       scratch );

    if( !strcasecmp( scratch, "true" ) )
    {
        Knots_Fixed = 1;
        Fixed_Knots = Malloc( Dk, double );

        if( inpf_unmatched( "vmpfx.knot.values" ) < 0 )
        {
            for( i = 0; i < Dk; i++ ) 
                Fixed_Knots[i] = (i + 1.0)/(Dk + 1.0);
        }
        else
        {
            i = inpf_get( "vmpfx.knot.values", "%V[%:%lf]", Dk, Fixed_Knots );

            if( !i )
            {
                Warning( DIAGNOSTIC >= 1, "Could not read full set of fixed knot values." );
                
                for( i = 0; i < Dk; i++ ) 
                    Fixed_Knots[i] = (i + 1.0)/(Dk + 1.0);
            }
        }
    }
    else
    {
        Knots_Fixed = 0;                          /* Set to -1 for transient fixing later ...    */

        Fixed_Knots = Malloc( Dk, double );       /* ... so we record default knots just in case */

        for( i = 0; i < Dk; i++ ) 
            Fixed_Knots[i] = (i + 1.0)/(Dk + 1.0);
    }

        /* Initialization Specification  */

    inpf_get( "vmpfx.init",       "%lV[=:%ls]", &Num_Init, &Init_List );
    inpf_get( "vmpfx.init.type",  "%ls[float]", &Init_Type_String );
    inpf_get( "vmpfx.init.iter",  "%d[50]",      &Init_Max_Iter );
    inpf_get( "vmpfx.init.tol",   "%lf[1e-4]",  &Init_Tolerance );

    if( inpf_unmatched( "vmpfx.init" ) )
        Num_Init = 0;    /* Auto Initialization for all parameters */

    fp[INIT] = init_init( &Num_Init, Init_List );

        /* Output Product Specification  */

    inpf_get( "vmpfx.produce",    "%lV[=:%ls[estimates]]", &Num_Products, &Products );
    inpf_get( "vmpfx.model",      "%ls[full]", &model_string );

    parse_model( model_string, &model );

    inpf_get( "vmpfx.null.model", "%s[false]", scratch );       /* Deprecated: Use vmpfx.model instead */

    if( !strcasecmp( scratch, "true" ) )
        model |= NULL_MODEL;

    if( model & NULL_MODEL ) 
        Null_Model = 1;
    else
        Null_Model = 0;

    Free( model_string );  /* Use Later */
    model_string = NULL;


    inpf_arg( "vmpfx.first.voxel", "first" );
    inpf_arg( "vmpfx.last.voxel",  "last" );

    inpf_get( "vmpfx.first.voxel", "%d[0]",  &first_voxel );
    inpf_get( "vmpfx.last.voxel",  "%d[-1]", &last_voxel );

        /* Check and Report any input errors    */

    if( inpf_errors( ) )
    {
        Message( "Errors encountered in input file:\n" );
        inpf_perror( );
        exit( 0 );
    }

        /* Process unmatched parameters */

    if( inpf_unmatched( "data.file" ) )
        Abort( "Data file required but none supplied in input file" );


        /* Clean up input structures */

    inpf_cleanup();
    cl_cleanup();


        /* Set Auxilliary Values     */

    o_SigmaSq_bsq = 1.0/(SigmaSq_b * SigmaSq_b);


        /* Validate Input Parameters */

    if( IAI <= 0.0 )
        Abort( "IAI must be a positive number (IAI = %lf).", IAI );

    if( Dg > DEGREE_MAX )
        Abort( "Drifts of degree > %d not supported (Degree = %d).", DEGREE_MAX, Dg );

    if( Dk > KNOT_MAX )
        Abort( "No more than %d knots are currently supported (Knots = %d).", KNOT_MAX, Dk );


    /** Initialization and Set-Up **/

        /*
         * Initialize Parameter Configuration
         *
         * In total, there are parameters for Baseline, Drift, Knots, Resp,
         * Shape, and Noise.  However, not all of these are used in any
         * specific run.  The parameter vector is re-ordered so that the
         * active variables are in the beginning, and then for each component
         * we set an index vector to record the indices of the component
         * parameters.
         *
         * The basic structure of the parameter vector is as follows.
         *
         * [ Baseline  Drift  Noise  Shape  Resp | Resp Drift Shape ]
         * 
         *         active components              inactive components
         *
         * This ordering is maintained throughout.  Depending on the user's
         * specification different parameters within Drift, Resp, and Shape
         * will be fixed.  The Baseline and Drift components are contiguous to
         * ease some computations.  The Resp components are at the boundary to
         * facilitate rapid changes of model.  The Shape parameters will be
         * folded so all of them should typically be used.  This is all set below.
         *
         * Parameter Descriptions:
         *
         *   NParams              Contains the total length of the parameter vector
         *   NeffParams           Contains the number of active parameters
         *   Dg, Dk, D            Drift Degree, Number of Knots, Number of Drift Basis Functions
         *   K, Keff              Total and active number of conditions
         *   Shape_Params         Number of shape parameters used
         *   
         *
         */    

    init_fixed( K, &Keff );  /* Set up constrained conditions */

        /* How many parameters do we have? */

            /* Total Number of Basis Functions not including constant             */

    D = Dg + Dk;

            /* Check for null activation boundary case */

    if( !Keff )
        Shape_Params = 0;

    if( !Shape_Params )
    {
        for( k = 0; k < K; k++ )
            Constrained[k] = 1;
        
        Keff = 0;
    }

            /* Total Pars: Baseline +  Drift (+ Knots) + Resp + Shape + Noise     */
            /* Effective:  Eliminate fixed resp's, knots, and unused shape params */

    NParams     = 1 + D + Dk + K + TOTAL_SHAPE_PARAMS + 1;
    NeffParams  = 1 + D + Keff + Shape_Params + 1;
    NoutParams  = 1 + D + K + Shape_Params + 1;
    NnulParams  = 2 + D;

    if( !Knots_Fixed )
    {
        NeffParams += Dk;
        NoutParams += Dk;
        NnulParams += Dk;
    }

        /* Set Parameter Indices */

    init_params( );

        /*
         * Prepare Data
         *
         *
         */

    fp[DATA] = init_data( data_fn, data_fmt, &data_hdr, &type, &V, &T, &Image_Dims, &inputDataAccessor );

    if( (last_voxel < 0) || (last_voxel >= V - 1) )
        last_voxel = V - 1;

    if( (first_voxel < 0) || (first_voxel >= V) )
        Abort( "Improper vmpfx.first.voxel provided (%d): %s\n",
               first_voxel, (first_voxel < 0) ? "Negative Index" : "Index larger than total count" );

    Veff = last_voxel - first_voxel + 1;    /* Number of voxels to be processed */

    
        /*
         * Initialize Design and Bell Configuration
         *
         * Constructs tables that specify the design.  This is used
         * in the construction of bell functions and activation profiles.
         * All stimulus times are converted into units of seconds
         * if given in terms of images.  Which tables are useful depends
         * on whether or not the design is aligned.  See the comments
         * in init_design() for more details.
         *
         */    

    init_design();

        /*
         * Additional Memory Allocation and Data Preparation
         *
         * The given product list is converted into the information
         * needed to output the components.  This depends on the 
         * parameterization used.
         *
         */

    init_products( &Num_Products, &Products, &Produce, &Product_Type,
                   &Units, &UnitSize, &Type, &Record_Size, type );

    Y           = Calloc( T+1, double );
    R           = Calloc( T+1, double );
    Rsave       = Calloc( T+1, double );
    Active_prof = Calloc( T+1, double );
    Drift_prof  = Calloc( T+1, double );

    Bells  = Calloc( UniqueStims * 2 * T, double );
    DBells = Calloc( TOTAL_SHAPE_PARAMS, double * );
    Bell_Length = Calloc( UniqueStims + 1, int );

    DBells[0] = Calloc( UniqueStims * TOTAL_SHAPE_PARAMS * 2 * T, double );

    for( i = 1; i < TOTAL_SHAPE_PARAMS; i++ )
        DBells[i] = DBells[0] + i * UniqueStims * 2 * T;

    atbounds = Calloc( NParams + 2, int );

    if( !RampLen )
    {
        RampLen = 8192;     /* NOTE: Put more context dependent value in later */
    }
    
    AttackRamp  = Calloc( RampLen + 20, double );
    DecayRamp   = Calloc( RampLen + 20, double );
    DAttackRamp = Calloc( RampLen + 20, double );
    DDecayRamp  = Calloc( RampLen + 20, double );

    if( Shape_Params > 4 )
    {
        RiseRamp  = Calloc( RampLen + 20, double );
        FallRamp  = Calloc( RampLen + 20, double );
        DRiseRamp = Calloc( RampLen + 20, double );
        DFallRamp = Calloc( RampLen + 20, double );
    }
    
    /*
     * TEMPORARY!
     *
     * Read in Attack and Decay Ramps from fixed file.
     * Eventually code these directly, give the function
     * call here.
     *
     */

    {
        FILE  *sfp;

        sfp = fopen( "ramps.txt", "r" );

        if( sfp == NULL )
        {

          char* fiasco_path;
          fiasco_path= getenv("FIASCO");
          if (!fiasco_path) Abort("FIASCO environment variable not set!\n");
            strcpy( scratch, fiasco_path );
            strcat( scratch, "/../../src/brain/");
            strcat( scratch, "ramps.txt");
           
            sfp = fopen( scratch, "r" );

            if( sfp == NULL )
                Abort( "Cannot open ramp file ramps.txt" );

        }

        fscanf( sfp, "%d", &RampLen );

        for( i = 0; i <= RampLen; i++ )
            fscanf( sfp, "%lf %lf %lf %lf", AttackRamp+i, DecayRamp+i, DAttackRamp+i, DDecayRamp+i );

        for( i = RampLen + 1; i < RampLen + 20; i++ )    /* Pad the attack */
            AttackRamp[i] = 1.0;
    }

    /* END TEMPORARY */

    WorkSize = Max( (2 * (T + 1) * (NParams + 5)), 4*(NParams * (NParams + 1) + 4*(T+1)) );
    IWorkSize = 2*(NParams + 2);

    Work  = Malloc( WorkSize,  double );
    IWork = Calloc( IWorkSize, int );
    Save =  Calloc( IWorkSize, double );

    minimizerWorkSize= 5000; /* no easy way to calculate from here */
    minimizerIWorkSize= 500;

    minimizerWork= Calloc( minimizerWorkSize, double );
    minimizerIWork= Calloc( minimizerIWorkSize, int );

        /*
         * Initialize Drift Basis
         *
         *
         */

    init_basis( );

    
        /*
         * Initialize Matrices for Initial Fitting
         *
         */

    init_init2( );


        /*
         * Initialize Objective Scale
         *
         */

    psc = Calloc( NParams + 1, double );
    fsc = 0.0;

    init_obj_scale( fastopt, psc, &fsc, IWork, Work, WorkSize );

        /*
         * Initialize Standard Error Computation
         *
         * (if applicable)
         *
         */

    if( Null_Model || Produce[STANDARD_ERRORS] || Produce[COVARIANCES] )
    {
        SEpsc = Calloc( NParams + 1, double );
        
        init_standard_errors( SEpsc, IWork, Work, WorkSize );
    }

    
    /** Scan previous contents of results file (if any) **/

    if( Output_exists )
    {
        if( !inpf_scan( Output_fn ) )
        {
            Message( "Errors encountered in results descriptor file %s\n", Output_fn );
            inpf_perror( );
            exit( 0 );
        }

        inpf_get( "file.name",     "%ls",     &name );
        inpf_get( "file.contents", "%lM[=:%ls %lu %lu[0]]", 
                  &Num_Contents, &Contents, &Contents_Address, &Contents_Chunk );
        inpf_get( "chunksize",     "%lu[%]",  DEFAULT_CHUNKSIZE, &Output_chunksize );

        /* The first two fields above are mandatory; check */

        if( inpf_unmatched( "file.name" ) )
            Abort( "Binary file name (tag file.name) missing in descriptor file %s", Output_fn );

        if( inpf_unmatched( "file.contents" ) )
            Abort( "Binary file contents (tag file.contents) missing in descriptor file %s", Output_fn );

        /* Maintain same path for descriptor and binary file */

        if( (*name != '/') && (path = strrchr(Output_fn,'/')) ) 
        {
            Res_fn = Calloc( path - Output_fn + strlen(name) + 2, char );
            strncpy( Res_fn, Output_fn, path - Output_fn + 1 );
            Res_fn[path - Output_fn + 1] = '\0';
            strcat( Res_fn, name );

            Free( name );
        }
        else
        {
            Res_fn = name;
            name = NULL;
        }

        inpf_cleanup( );
    }
    else
    {
        Num_Contents = 0;
        Contents = NULL;
        Contents_Chunk = NULL;
        Contents_Address = NULL;

        Output_chunksize = DEFAULT_CHUNKSIZE;
    }

    /** Begin Output File **/

        /* File Contents: Volumes */

    mpsize = (unsigned long)(Veff*Record_Size);

    if( !Output_exists || !Append_out )
    {
        fp[OUTPUT] = efopen( Output_fn, "w" );

        path = strrchr( Output_fn, '/' );

        if( path && !strncmp( Output_fn, Res_fn, path - Output_fn + 1 ) )
            fprintf( fp[OUTPUT], "file.name         %s\n",  &Res_fn[path - Output_fn + 1] );
        else
            fprintf( fp[OUTPUT], "file.name         %s\n",  Res_fn );

        fprintf( fp[OUTPUT], "\n" );

        fprintf( fp[OUTPUT], "<data> file.contents\n" );
        fprintf( fp[OUTPUT], "#      component              length      chunks\n" );
        fprintf( fp[OUTPUT], "       max.posterior      %10lu      %5lu\n",
                  mpsize % Output_chunksize, mpsize/Output_chunksize );
        fprintf( fp[OUTPUT], "<end>\n" );
        fprintf( fp[OUTPUT], "\n" );

        fprintf( fp[OUTPUT], "chunksize    %lu\n", Output_chunksize );

        fprintf( fp[OUTPUT], "\n%s\n", End_Hdr_String );
    }
    else        
    {
        /* Not yet decided on whether to allow overwrites, see temp code above */
    }

        /* File Contents: Records in max.posterior volume */
    
    fprintf( fp[OUTPUT], "\n#\n# Maximum Posterior Fit\n#\n\n" );

    fprintf( fp[OUTPUT], "<data> max.posterior.contents\n" );
    fprintf( fp[OUTPUT], "#      name                 units    size    type\n" );
    fprintf( fp[OUTPUT], "       mp                   %4d     %lu    record\n", Veff, Record_Size );
    fprintf( fp[OUTPUT], "<end>\n" );
    fprintf( fp[OUTPUT], "\n" );
    
    if( Num_Products > 0 )
    {                         
        fprintf( fp[OUTPUT], "<data> max.posterior/mp.contents\n" );
        fprintf( fp[OUTPUT], "#      name                     units    size  type\n" );

        for( i = 0; i < Num_Products; i++ )
        {
            j = Product_Type[i];
            fprintf( fp[OUTPUT], "       %-24s  %4ld    %4ld  %s\n",
                     Products[i], Units[j], UnitSize[j], Type_String[Type[j]] );
        }

        fprintf( fp[OUTPUT], "<end>\n" );
        fprintf( fp[OUTPUT], "\n" );
    }

    if( Produce[ESTIMATES] )
    {
        tsz = sizeof(double);

        fprintf( fp[OUTPUT], "<data> max.posterior/mp/estimates.contents\n" );
        fprintf( fp[OUTPUT], "#     name            units    size    type\n" );
        fprintf( fp[OUTPUT], "      baseline         %4d    %4d    %s\n", 1,  (int)tsz, "double" );
        fprintf( fp[OUTPUT], "      drift.coef       %4d    %4d    %s\n", D,  (int)tsz, "double" );

            if( !Knots_Fixed )
                fprintf( fp[OUTPUT], "      drift.knots      %4d    %4d    %s\n", Dk, (int)tsz, "double" );

        fprintf( fp[OUTPUT], "      noise.precision  %4d    %4d    %s\n", 1,  (int)tsz, "double" );
        fprintf( fp[OUTPUT], "      shape            %4d    %4d    %s\n", Shape_Params, (int)tsz, "double");
        fprintf( fp[OUTPUT], "      responsiveness   %4d    %4d    %s\n", K,  (int)tsz, "double" );
        fprintf( fp[OUTPUT], "<end>\n" );
        fprintf( fp[OUTPUT], "\n" );

        if( Null_Model )
        {
            fprintf( fp[OUTPUT], "<data> max.posterior/mp/null.estimates.contents\n" );
            fprintf( fp[OUTPUT], "#     name            units    size    type\n" );
            fprintf( fp[OUTPUT], "      baseline         %4d    %4d    %s\n", 1,  (int)tsz, "double" );
            fprintf( fp[OUTPUT], "      drift.coef       %4d    %4d    %s\n", D,  (int)tsz, "double" );

            if( !Knots_Fixed )
                fprintf( fp[OUTPUT], "      drift.knots      %4d    %4d    %s\n", Dk, (int)tsz, "double" );

            fprintf( fp[OUTPUT], "      noise.precision  %4d    %4d    %s\n", 1,  (int)tsz, "double" );
            fprintf( fp[OUTPUT], "<end>\n" );
            fprintf( fp[OUTPUT], "\n" );
        }
    }

    if( Produce[STANDARD_ERRORS] || Produce[COVARIANCES] )
    {
        tsz = sizeof(double);

        fprintf( fp[OUTPUT], "<data> max.posterior/mp/standard.errors.contents\n" );
        fprintf( fp[OUTPUT], "#     name            units    size    type\n" );
        fprintf( fp[OUTPUT], "      baseline         %4d    %4d    %s\n", 1,  (int)tsz, "double" );
        fprintf( fp[OUTPUT], "      drift.coef       %4d    %4d    %s\n", D,  (int)tsz, "double" );

            if( !Knots_Fixed )
                fprintf( fp[OUTPUT], "      drift.knots      %4d    %4d    %s\n", Dk, (int)tsz, "double" );

        fprintf( fp[OUTPUT], "      noise.precision  %4d    %4d    %s\n", 1,  (int)tsz, "double" );
        fprintf( fp[OUTPUT], "      shape            %4d    %4d    %s\n", Shape_Params, (int)tsz, "double");
        fprintf( fp[OUTPUT], "      responsiveness   %4d    %4d    %s\n", K,  (int)tsz, "double" );
        fprintf( fp[OUTPUT], "<end>\n" );
        fprintf( fp[OUTPUT], "\n" );

        if( Null_Model )
        {
            fprintf( fp[OUTPUT], "<data> max.posterior/mp/null.standard.errors.contents\n" );
            fprintf( fp[OUTPUT], "#     name            units    size    type\n" );
            fprintf( fp[OUTPUT], "      baseline         %4d    %4d    %s\n", 1,  (int)tsz, "double" );
            fprintf( fp[OUTPUT], "      drift.coef       %4d    %4d    %s\n", D,  (int)tsz, "double" );

            if( !Knots_Fixed )
                fprintf( fp[OUTPUT], "      drift.knots      %4d    %4d    %s\n", Dk, (int)tsz, "double" );

            fprintf( fp[OUTPUT], "      noise.precision  %4d    %4d    %s\n", 1,  (int)tsz, "double" );
            fprintf( fp[OUTPUT], "<end>\n" );
            fprintf( fp[OUTPUT], "\n" );
        }
    }

    if( Produce[DIAGNOSTICS] )
    {
        tsz = sizeof(float);

        fprintf( fp[OUTPUT], "<data> max.posterior/mp/diagnostics.contents\n" );
        fprintf( fp[OUTPUT], "#     name            units    size    type\n" );
        fprintf( fp[OUTPUT], "      posterior.max    %4d    %4d    %s\n", 1, (int)tsz, "float"  );

        if( Null_Model || Produce[STANDARD_ERRORS] || Produce[COVARIANCES] )
        {
            fprintf( fp[OUTPUT], "      cov.rank         %4d    %4d    %s\n", 1, (int)tsz,       "float"  );
            fprintf( fp[OUTPUT], "      cov.lndet        %4d    %4d    %s\n", 1, (int)tsz,       "float"  );
        }
        
        fprintf( fp[OUTPUT], "      num.iterations   %4d    %4d    %s\n", 1, (int)(sizeof(int)), "integer" );
        fprintf( fp[OUTPUT], "      num.evaluations  %4d    %4d    %s\n", 1, (int)(sizeof(int)), "integer" );
        
        fprintf( fp[OUTPUT], "<end>\n" );
        fprintf( fp[OUTPUT], "\n" );

        if( Null_Model )
        {
            fprintf( fp[OUTPUT], "<data> max.posterior/mp/null.diagnostics.contents\n" );
            fprintf( fp[OUTPUT], "#     name            units    size    type\n" );
            fprintf( fp[OUTPUT], "      posterior.max    %4d    %4d    %s\n", 1,   (int)tsz, "float"  );
            fprintf( fp[OUTPUT], "      cov.rank         %4d    %4d    %s\n", 1,   (int)tsz, "float"  );
            fprintf( fp[OUTPUT], "      cov.lndet        %4d    %4d    %s\n", 1,   (int)tsz, "float"  );
            fprintf( fp[OUTPUT], "      ln.bayes.factor  %4d    %4d    %s\n", 1,   (int)tsz, "float"  );
            fprintf( fp[OUTPUT], "      num.iterations   %4d    %4d    %s\n", 1, (int)(sizeof(int)), "integer" );
            fprintf( fp[OUTPUT], "      num.evaluations  %4d    %4d    %s\n", 1, (int)(sizeof(int)), "integer" );

            fprintf( fp[OUTPUT], "<end>\n" );
            fprintf( fp[OUTPUT], "\n" );
        }
    }


    /** Prepare Binary Results File **/

    if( Append_out )
        fp[BINARY] = efopen( Res_fn, "ab" );
    else
        fp[BINARY] = efopen( Res_fn, "wb" );

    output_bufsize = Record_Size * (memlimit/Record_Size);  /* Fit as many records as possible */
    output_buf = Malloc( output_bufsize, char );


    /** Select Algorithm and Configurable Functions **/

        /* Set Configurable Functions */

    if( Knots_Fixed )
    {
        lnpost = lnpost_fixed;
        Dlnpost = Dlnpost_fixed;
        DDlnpost = DDlnpost_fixed;
        lnpost_null = lnpost_null_fixed;
        Dlnpost_null = Dlnpost_null_fixed;
        DDlnpost_null = DDlnpost_null_fixed;
        fastopt = 1;
    }
    else
    {
        Abort( "Variable Knots not yet supported in vmpfx." );
        
        lnpost = lnpost_variable;
        Dlnpost = Dlnpost_variable;
        lnpost_null = lnpost_null_variable;
        Dlnpost_null = Dlnpost_null_variable;
        fastopt = 0;

        DDlnpost = NULL;
        DDlnpost_null = NULL;
        
        if( cov_method == COV_ANALYTIC )
            Abort( "Analytic Hessian Not Supported for Variable Knots." );
    }

    switch( Shape_Params )
    {
      case 0:
      case 2:
      case 3:
      case 4:

        poly_bell = poly_bell4;
        shape_deriv = shape_deriv_4;
        break;

      case 5:
      case 6:
            
        poly_bell = poly_bell6;
        shape_deriv = shape_deriv_6;
        break;

      case 7:
      case 8:

        poly_bell = poly_bell8;
        shape_deriv = shape_deriv_8;
        break;
    }

    /* Optimization package initialization */
    qnb_initialize(NeffParams, NParams);


        /* Set Optimization Parameters */

    if( fastopt )
    {
        qnb_config_int( "max iterations", Iter_Max );
        qnb_config_int( "max evals",      Eval_Max );
        qnb_config_int( "max grad evals", Eval_Max );
        qnb_config_int( "hessian type",   1 );
	qnb_config_int( "corrections", 5 );

        qnb_config_double( "gradient tolerance", 1.0e-4 );  /* ATTN */
        
        for( i = 0; i < Num_QNB_Int_Opts; i++ )
            qnb_config_int( QNB_Int_Options_str[i], 
			    QNB_Int_Options_val[i] );

        for( i = 0; i < Num_QNB_Double_Opts; i++ )
            qnb_config_double( QNB_Double_Options_str[i], 
			       QNB_Double_Options_val[i] );

    }
    else
    {
        min_config_nmb_int( "auto clean", 0 );
        min_config_nmb_double( "init change", 0.001 );
        min_config_nmb_double( "threshold", 5.0 );
    }

#   if defined(USE_IMSL)

    erset( &iZero, &imOne, &iZero );     /* Set IMSL Error Handling to not stop */

#endif


    /***  Main Processing  ***/

    /*
     * For each voxel:
     *
     *    o Convert data to double
     *    o Initialize parameters (auto and user specified)
     *    o Compute Maximum Posterior estimates
     *    o Derive and store desired output quantities
     *
     */

    size = mriTypeSize( type );
    voxdat = emalloc( T * sizeof(double) );
    initbuf = emalloc( Init_len * Init_size );
    initdat = Malloc( Init_len, double );

    if( first_voxel > 0 )    /* Skip Initial Time Series */
        (*(inputDataAccessor.setup))( &inputDataAccessor, fp[DATA], first_voxel, T, size );

    output_cookie( SET_IT, output_buf, output_bufsize, fp[BINARY] );   /* Prepare Output Buffer    */


    for( v = first_voxel; v <= last_voxel; v++ )
    {
#if defined( USE_SIGNALS )
        set_current_voxel( v );      /* To describe where we are on an abort */
#endif

        /* Input next time course and convert */

        (*(inputDataAccessor.read))(&inputDataAccessor, voxdat, fp[DATA], T, size);
        mriTypeConvert( T, DOUBLE_T, Y, type, voxdat );

        /* Check for zero or constant time course */
        
        {
            double     sumysq = 0.0;
            double     sumy   = 0.0;
            
            for( i = 0; i < T; i++ )
            {
                sumysq += Y[i]*Y[i];
                sumy   += Y[i];
            }

            if( sumysq < Tolerance || (sumysq - sumy*sumy/T) < Tolerance ) /* All 0 or sample variance 0 */
            {
                output_null_record( "Time course constant", v, type, NoutParams, T, Null_Model, NnulParams,
                                    saveParam, voxdat, ObsInfo );
                continue;  /* On to the next voxel */
            }
        }
        
        /*
         * ATTN: convert the above to a generic function get_voxel_ts()
         *       that can be sensitive to the data format.
         *
         *       Also if Image data, set slice offset here based on order and general offset
         *
         */

        /*
         * TEMPORARY HACK:  Dealing with Missing Data
         *
         *   Replace -1's with average of neighbors
         *
         */

        if( Y[0] == -1 )
            Y[0] = Y[1];

        if( Y[T-1] == -1 )
            Y[T-1] = Y[T-2];

        for( t = 1; t < T-1; t++ )
            if( Y[t] == -1 )
                Y[t] = 0.5 * (Y[t-1] + Y[t+1]);

        /* END TEMPORARY HACK */

        /* Initialize Parameter Vector        */

        if( Init_len )
        {
            efread( initbuf, Init_size, Init_len, fp[INIT] );
            mriTypeConvert( Init_len, DOUBLE_T, initdat, Init_type, initbuf );
            set_inits( T, Y, NParams, NeffParams, Param, Init, initdat, Work, WorkSize );
        }
        else
            init_retcode= voxel_init( T, Y, NParams, NeffParams, 0, Param, psc, Work, WorkSize );

        /* Check that initialized values are well defined */

        for( i = 0; i < NeffParams; i++ )
            if( !finite( Param[i] ) )
                break;

        if( !init_retcode || i < NeffParams )
        {
	    Warning(1,"Voxel initialization failed for voxel %d!\n",v);
            output_null_record( "NaN's encountered in initialized estimates", v, type, NoutParams, T, Null_Model, NnulParams,
                                saveParam, voxdat, ObsInfo );
            continue;  /* On to the next voxel */
        }
        

        /* Fit Parameters */

        voxel_fit( T, Y, NParams, NeffParams, Param, &val, &iter, &nfunc, psc, fsc,
                   atbounds, Save, IWorkSize, fastopt, lnpost, Dlnpost );

        /* Construct Covariance Matrix */

        if( Null_Model || Produce[STANDARD_ERRORS] || Produce[COVARIANCES] )
        {
            /* Save Residuals          */

            if( (Produce[RESIDUALS] || Produce[FITTED_VALUES]) && (cov_method != COV_ANALYTIC) )
            {
                for( t = 0; t < T; t++ )
                    Rsave[t] = R[t];
            }
                
            /* Compute Estimated Covariance Matrix */

            compute_cov( NeffParams, NParams, Param, ObsInfo, cov_method, atbounds,
                         SEwork, SEworksize, &cov_rank, &cov_lndet, &info, 0 );

            if( info )
                Warning( DIAGNOSTIC >= 1, "Could not compute covariances for voxel %d (info = %d).",
                         v, info );

            /* Restore Residuals       */

            if( (Produce[RESIDUALS] || Produce[FITTED_VALUES]) && (cov_method != COV_ANALYTIC) )
            {
                for( t = 0; t < T; t++ )
                    R[t] = Rsave[t];
            }
        }

        /* Output */

#if defined(MONITOR) && (MONITOR > 0)

        if( !( (v - first_voxel) % MONITOR) )  /* ATTN: TEMP TO MONITOR PROGRESS */
        {
            time( &exectime );
            Message( "Voxel %4d    %20.15lf in %d iterations at %s.\n\n", v, -val, iter, ctime( &exectime ) );
        }

#endif

#if DIAGNOSTIC >= 5
            Message( "Voxel %4d    %20.15lf in %d iterations.\n\n", v, -val, iter );
#endif

        for( i = 0; i < Num_Products; i++ )
        {
            if( Product_Type[i] == ESTIMATES )
            {
                for( m = 0; m < NoutParams; m++ )
                    saveParam[Out[m]] = Param[m];

                output_cookie( ESTIMATES, NoutParams, saveParam );
            }
            else if( Product_Type[i] == STANDARD_ERRORS )
            {
                /* Record SE Values in Parameter File */

                for( m = 0; m < NeffParams; m++ )
                {
                    u = ObsInfo[m + (m*(m+1))/2];
                    saveParam[Out[m]] = (u > 0.0) ? sqrt( u ) : (2.0 * fabs(Param[m]));
                }
                
                for( ; m < NoutParams; m++ )
                    saveParam[Out[m]] = 0.0;

                output_cookie( STANDARD_ERRORS, NoutParams, saveParam );
            }
            else if( Product_Type[i] == RESIDUALS )
            {
                if( type == DOUBLE_T || type == FLOAT_T )
                    type_r = type;
                else
                    type_r = FLOAT_T;

                size_r = mriTypeSize( type_r );
            
                mriTypeConvert( T, type_r, voxdat, DOUBLE_T, R );  /* voxel_fit() must set R at end */
            
                output_cookie( RESIDUALS, T, size_r, voxdat );
            }
            else if( Product_Type[i] == FITTED_VALUES )
            {
                if( type == DOUBLE_T || type == FLOAT_T )
                    type_r = type;
                else
                    type_r = FLOAT_T;

                size_r = mriTypeSize( type_r );

                for( t = 0; t < T; t++ )
                    Work[t] = Y[t] - R[t];                      /* voxel_fit() must set R at end */

                mriTypeConvert( T, type_r, voxdat, DOUBLE_T, Work );
            
                output_cookie( FITTED_VALUES, T, size_r, voxdat );
            }
            else if( Product_Type[i] == COVARIANCES )
            {
                /* Output Covariance matrix in packed form */

                output_cookie( COVARIANCES, (NoutParams * (NoutParams+1))/2, ObsInfo );
            }
            else if( Product_Type[i] == DIAGNOSTICS )
            {
                rdiagnostics[0] = -val;                 /* Maximized Log Posterior  */

                if( Null_Model || Produce[STANDARD_ERRORS] || Produce[COVARIANCES] )
                {
                    rdiagnostics[1] = cov_rank;         /* Rank                     */
                    rdiagnostics[2] = cov_lndet;        /* Log Determinant          */
                }

                idiagnostics[0] = iter;    /* Number of Iterations     */
                idiagnostics[1] = nfunc;   /* Number of Function Evals */

                j = Null_Model | Produce[STANDARD_ERRORS] | Produce[COVARIANCES];

                output_cookie( DIAGNOSTICS, 1 + j*2, 2, rdiagnostics, idiagnostics );
            }
        }
        
        /* ATTN: TEMP */
#if DIAGNOSTIC >= 7
        if( Produce[STANDARD_ERRORS] || Produce[COVARIANCES] )
            for( i = 0; i < NoutParams; i++ )
                fprintf( stderr, "%20.10lg   %20.10lg\n", Param[i], sqrt(ObsInfo[i + (i*(i+1))/2]) );
        else
            for( i = 0; i < NoutParams; i++ )
                fprintf( stderr, "%20.10lg\n", Param[i] );
#endif

        /* END TEMP */


        /* Null Model Fit and Output */

        if( Null_Model )
        {
            /* ATTN: In collapsed case, no need to refit */

            init_retcode= voxel_init( T, Y, NParams, 2 + D, 1, Param, NULL, Work, WorkSize );

            for( k = 0; k < K; k++ )              /* Cosmetic Alone */
                Param[Resp[k]] = 0.0;

            for( i = 0; i < Shape_Params; i++ )
                Param[Shape[i]] = 0.0;

            voxel_fit( T, Y, NParams, 2 + D, Param, &val0, &iter, &nfunc, psc, fsc,
                      atbounds, Save, IWorkSize, fastopt, lnpost_null, Dlnpost_null );

            /* Prepare Output of Null Results */

                /* Save Residuals          */

            if( (Produce[RESIDUALS] || Produce[FITTED_VALUES]) && (cov_method != COV_ANALYTIC) )
            {
                for( t = 0; t < T; t++ )
                    Rsave[t] = R[t];
            }
                
                /* Compute Estimated Covariance Matrix */

            compute_cov( 2 + D, NParams, Param, ObsInfo, cov_method, atbounds,
                         SEwork, SEworksize, &cov0_rank, &cov0_lndet, &info, 1 );

            if( info )
                Warning( DIAGNOSTIC >= 1, "Could not compute covariances for voxel %d (info = %d).",
                         v, info );

                /* Restore Residuals       */

            if( (Produce[RESIDUALS] || Produce[FITTED_VALUES]) && (cov_method != COV_ANALYTIC) )
            {
                for( t = 0; t < T; t++ )
                    R[t] = Rsave[t];
            }

            /* Output */

#if DIAGNOSTIC >= 6
            Message( "\nNull Model:  Voxel %4d    %8.3lf in %d iterations.\n\n", v, -val0, iter );
#endif

            /* ATTN: TEMP */

#if DIAGNOSTIC >= 8
            if( Produce[STANDARD_ERRORS] || Produce[COVARIANCES] )
                for( i = 0; i < D + 2; i++ )
                    fprintf( stderr, "%20.10lg   %20.10lg\n", Param[i], sqrt(ObsInfo[i + (i*(i+1))/2]) );
            else
                for( i = 0; i < D + 2; i++ )
                    fprintf( stderr, "%20.10lg\n", Param[i] );
#endif

            /* END TEMP */

            for( i = 0; i < Num_Products; i++ )
            {
                if( Product_Type[i] == ESTIMATES )
                {
                    for( m = 0; m < NnulParams; m++ )
                        saveParam[Out[m]] = Param[m];

                    output_cookie( ESTIMATES, NnulParams, saveParam );
                }
                else if( Product_Type[i] == STANDARD_ERRORS )
                {
                    /* Record SE Values in Parameter File */

                    for( m = 0; m < NnulParams; m++ )
                    {
                        u = ObsInfo[m + (m*(m+1))/2];
                        saveParam[Out[m]] = (u > 0.0) ? sqrt( u ) : (2.0 * fabs(Param[m]));
                    }
                
                    output_cookie( STANDARD_ERRORS, NnulParams, saveParam );
                }
                else if( Product_Type[i] == RESIDUALS )
                {
                    if( type == DOUBLE_T || type == FLOAT_T )
                        type_r = type;
                    else
                        type_r = FLOAT_T;

                    size_r = mriTypeSize( type_r );
            
                    mriTypeConvert( T, type_r, voxdat, DOUBLE_T, R );  /* voxel_fit() must set R at end */
            
                    output_cookie( RESIDUALS, T, size_r, voxdat );
                }
                else if( Product_Type[i] == FITTED_VALUES )
                {
                    if( type == DOUBLE_T || type == FLOAT_T )
                        type_r = type;
                    else
                        type_r = FLOAT_T;

                    size_r = mriTypeSize( type_r );

                    for( t = 0; t < T; t++ )
                        Work[t] = Y[t] - R[t];                      /* voxel_fit() must set R at end */

                    mriTypeConvert( T, type_r, voxdat, DOUBLE_T, Work );
            
                    output_cookie( FITTED_VALUES, T, size_r, voxdat );
                }
                else if( Product_Type[i] == COVARIANCES )
                {
                    /* Output Covariance matrix in packed form */

                    output_cookie( COVARIANCES, (NnulParams * (NnulParams + 1))/2, ObsInfo );
                }
                else if( Product_Type[i] == DIAGNOSTICS )
                {
                    rdiagnostics[0] = -val0;                                        /* Maximized Log Posterior  */
                    rdiagnostics[1] = cov0_rank;                                    /* Rank                     */
                    rdiagnostics[2] = cov0_lndet;                                   /* Log Determinant          */

                    /* Carefull of case when all resp 0 */
		    rdiagnostics[3] = -((val - val0) - 0.5*(cov_lndet - cov0_lndet)   /* Approximate Ln B-Factor  */
					- LN_SQRT_2PI * (cov_rank - cov0_rank));


                    idiagnostics[0] = iter;    /* Number of Iterations     */
                    idiagnostics[1] = nfunc;   /* Number of Function Evals */

                    /* ATTN: TEMP     */
#if DIAGNOSTIC >= 9
                    fprintf( stderr, "\nApproximate Ln(Bayes Factor) = %lg\n\n", rdiagnostics[3] );
#endif
                    /* END TEMP */

                    output_cookie( DIAGNOSTICS, 4, 2, rdiagnostics, idiagnostics );
                }
            }
        }


#if DIAGNOSTIC >= 7
            Message( "---------\n\n" );
#endif
    }

    output_cookie( WRITE_IT );


    /* Done: Close files and clean up...   */

    time( &exectime );       /* Record finish time (not perfect but useful info) */

    Message( "\n##\n## Done at %s##\n\n\n", ctime( &exectime ) );

    (*(inputDataAccessor.close))(&inputDataAccessor, fp[DATA]);
    qnb_cleanup();
    efclose( fp[OUTPUT] );
    efclose( fp[BINARY] );
    es_delete( &errorStream );

    exit( 0 );
}

int   voxel_init( int nt, double *y, int npar, int neffpar, int nullmodel,
                   double *p, double *psc, double *work, int lwork )
{
    int           d, i, j, k, m, s, t;
    int           iter, anyresp, info = 0;

    double        u, v, u2, v2, mulast, smp;
    double        value, fsc = 1.0;
    double        sum_y, sum_ysq;
    double        change_tol = 1.0e10;

    double       *ZtZ, *Rlast;
    double       *ap, *dp;

    int           minimizer_evals;
    double        minimizer_grad_tol;

    /*
     * Find good initial Point when no user initialization and fixed basis
     *
     * The strategy is to use component-wise backfitting to get a very good
     * starting point for the optimization.  If the model is null, then the
     * shape and responsiveness components are not fit.  A very simple initial
     * fit is used to get preliminary estimates for the noise level and the
     * baseline.  Given shape, a linear fit is used for drift and responsiveness,
     * shape is optimized nonlinearly, Smooth is set based on the noise level and
     * the target d.f., and the residual mean and (corrected) variance are used
     * to get mu and sigma^2.
     *
     * NOTE: If a null model is being fit, then the preliminary fit is all that
     *       need be computed.  
     *
     * NOTE: If StimulusStart[j+1] < StimulusEnd[j] for any j, the preliminary
     *       initialization will not work optimally.
     *
     * When user initialization occurs, more careful separation of the
     * possibilities is done in set_inits.  This is supposed to be
     * an automatic method in the absence of other info, which will
     * be the typical default.
     *
     */

    /** Set Initial Constant Quantities **/

    for( d = 0; d < D; d++ )
        Xty[d] = 0.0;

    for( k = 0; k < K; k++ )
    {
        S0ty[k] = 0.0;

        for( d = 0; d < D; d++ )
            XtS0[d + k*D] = 0.0;
    }

    for( t = 0, sum_y = sum_ysq = 0.0; t < T; t++ )
    {
        v = y[t];
        
        sum_y += v;
        sum_ysq += v * v;

        for( d = 0; d < D; d++ )
            Xty[d] += Basis[d][t] * v;

        for( k = 0; k < K; k++ )
        {
            u = CondMatrix[k][t];

            S0ty[k] += u * v;

            for( d = 0; d < D; d++ )
                XtS0[d + k*D] += u * Basis[d][t];
        }
    }

    /** Preliminary Fit **/

    smp = Smooth_Target;

    if( !Shape_Params || !Keff || nullmodel )
    {
        work[0] = sum_y;
        
        for( d = 0; d < D; d++ )
            work[d+1] = Xty[d];
        
        m = D + 1;
        ZtZ = work + m;

        lwork -= m*(m+1);

        ZtZ[0] = T;

        for( i = 1; i < m; i++ )
            ZtZ[i] = ZtZ[i*m] = 0.0;

        for( i = 1; i < m; i++ )
        {
            for( j = 1; j < i; j++ )
                ZtZ[i + j*m] = ZtZ[j + i*m] = smp * Ptot[j + i*m];

            ZtZ[i + i*m] = 1.0 + smp * Ptot[i + i*m];
        }

        dsytrf( &Upper, &m, ZtZ, &m, IWork, work + m*(m+1), &lwork, &info );

        if( !info )
            dsytrs( &Upper, &m, &iOne, ZtZ, &m, IWork, work, &m, &info );

        if( info ) {
            Warning(1, "Could not solve regression equations during initialization (info = %d)", info );
	    return 0;
	}
                
        lwork += m*(m+1);
        
        /* Now assign coefficients */

        p[Mu] = work[0];

        for( d = 0; d < D; d++ )
            p[Coefs[d]] = work[1 + d];

        /*
         * Set RSS and estimate r
         *
         * With orthogonal spline basis  (ATTN: B-splines)
         *
         *   || y - mu 1 - X delta ||^2 = y'y + mu^2 T + delta'delta - 2 mu sum_y - 2 delta'X'y
         */

        u = sum_ysq + p[Mu]*p[Mu]*T - 2.0 * p[Mu] * sum_y;

        for( d = 0; d < D; d++ )
            u += p[Coefs[d]]*p[Coefs[d]] - 2.0 * p[Coefs[d]] * Xty[d];

        /* p[o_SigmaSq] = (nt - Drift_df - 1)/u;         ATTN */
        p[o_SigmaSq] = (0.5 * nt + SigmaSq_a - 1.0)/(0.5 * u + SigmaSq_b);

        SmoothPar = (smp * p[o_SigmaSq]);           /* 1/lambda = target * r */

        return 1;  /* In null model, all done */
    }

        /* Otherwise, non-null model... */

    for( d = 0; d < D; d++ )
        work[d] = Xty[d];

    for( k = 0; k < K; k++ )
        work[D + k] = S0ty[k];
        
    m = D + K;
    ZtZ = work + m;
    lwork -= m*(m + 1);

    for( i = 0; i < D; i++ )      /* X'X + smp P  diagonal block */
    {
        for( j = 0; j < i; j++ )
            ZtZ[i + j*m] = ZtZ[j + i*m] = smp * Ptot[(j+1) + (i+1)*(D+1)];

        ZtZ[i + i*m] = 1.0 + smp * Ptot[(i+1) + (i+1)*(D+1)];
    }

    for( i = D; i < m; i++ )      /* S0'S0 diagonal block */
    {
        for( j = D; j < i; j++ )
            ZtZ[i + j*m] = ZtZ[j + i*m] = 0.0;

        ZtZ[i + i*m] = CondImages[i-D];
    }

    for( i = D; i < m; i++ )      /* X'S0 off-diagonal block with symmetry */
        for( j = 0; j < D; j++ )
            ZtZ[i + j*m] = ZtZ[j + i*m] = XtS0[ j + (i - D)*D ];

        
    dsytrf( &Upper, &m, ZtZ, &m, IWork, work + m*(m+1), &lwork, &info );

    if( !info )
        dsytrs( &Upper, &m, &iOne, ZtZ, &m, IWork, work, &m, &info );

    if( info ) {
      Warning(1, "Could not solve regression equations during initialization (info = %d)", info );
      return 0;
    }

    lwork += m*(m + 1);
        
    /*
     * Now assign baseline, drift, and resp coefficients
     *
     *  o Drift coefficients as computed
     *  o Baseline is obtained as a weighted average of constrained conditions
     *    and any that would be negative based on the averages.  The latter
     *    have responsiveness set to 0.
     *  o The rest of the resp parameters are set to the proportional change from Mu
     * 
     */

    for( d = 0; d < D; d++ )
        p[Coefs[d]] = work[d];

    if( Keff < K )
    {
        for( k = 0, u = v = 0.0; k < K; k++ )  /* Average a priori zero perturbations */
        {
            if( Constrained[k] )
            {
                u += CondImages[k] * work[D + k];
                v += CondImages[k];
            }
        }

        for( k = 0, u2 = u, v2 = v; k < K; k++ )  /* Add in those that are negative */
        {
            if( !Constrained[k] && work[D+k] < u/v )
            {
                u2 += CondImages[k] * work[D + k];
                v2 += CondImages[k];
            }
        }

        p[Mu] = u2/v2;

        for( k = 0; k < K; k++ )               /* Now set the parameters */
        {
            if( Constrained[k] || (work[D+k] < u/v) )
                p[Resp[k]] = 0.0;
            else
                p[Resp[k]] = (work[D + k] - p[Mu])/p[Mu];
        }
    }
    else
    {
        for( k = 0, i = -1, v = BigNumber; k < K; k++ )
            if( work[D + k] < v )
            {
                v = work[D + k];
                i = k;
            }

        p[Mu] = v;

        for( k = 0; k < K; k++ )
        {
            if( Constrained[k] || k == i )
                p[Resp[k]] = 0.0;
            else
                p[Resp[k]] = (work[D + k] - p[Mu])/p[Mu];
        }
    }
        

    /*
     * Set RSS and estimate r
     *
     * With orthogonal spline basis  (ATTN: B-splines)
     *
     *   || y - mu 1 - X delta - S gamma mu||^2 =
     *        y'y + mu^2 T + delta'delta + mu^2 gamma S'S gamma
     *          - 2 mu sum_y - 2 delta'X'y - 2 mu^2 1'S gamma - 2 y'S gamma mu - 2 delta X'S gamma mu
     *
     */

    u = sum_ysq + p[Mu]*p[Mu]*T - 2.0 * p[Mu] * sum_y;

    for( k = 0; k < K; k++ )
        u += p[Mu] * p[Mu] * p[Resp[k]] * p[Resp[k]] * CondImages[k] - 
            2.0 * p[Mu] * p[Resp[k]] * (S0ty[k] - p[Mu] * CondImages[k]); 

    for( d = 0; d < D; d++ )
        u += p[Coefs[d]]*p[Coefs[d]] - 2.0 * p[Coefs[d]] * Xty[d];

    for( d = 0; d < D; d++ )
    {
        v = 2.0 * p[Coefs[d]] * p[Mu];

        for( k = 0; k < K; k++ )
            u += v * p[Resp[k]] * XtS0[d + k*D];
    }

    p[o_SigmaSq] = (nt - Drift_df - Keff)/u;

    SmoothPar = (smp * p[o_SigmaSq]);           /* 1/lambda = target * r */


    /* Shape Parameters: use prior Mode */

    for( s = 0; s < TOTAL_SHAPE_PARAMS; s++ )
        if( Shape_Use[s] )
            p[Shape[s]] = (Shape_a[s] - 1.0)/Shape_b[s];    /* ATTN: Previously used Prior Mean */


    /* Drift Knots --- Knots_Fixed is true here */

    for( i = 0; i < Dk; i++ ) 
        p[Knots[i]] = Fixed_Knots[i];


    /** Backfitting Iterations **/

    /*
     * There are two choices for how to do the backfitting.  The first
     * is to use the Bayes estimators derived from the complete
     * conditional distributions and the second is to do a more direct
     * component analyses that maximize the likelihood/posterior.  The former is
     * substantially more work and will likely provide only small
     * gain.  Therefore, we do the latter here essentially.  If it
     * should turn out that more accuracy is needed at initialization
     * then we will implement the other one.  We do expect this to get
     * very close to the optimum.  Drift, Responsiveness, and Baseline
     * are basically computed using linear models (with truncation for
     * resp).  The shape parameters actually maximize a restricted
     * posterior (hopefully with not too much trouble), and the noise
     * is estimated based on the residuals, with Smoothness updated
     * appropriately.
     *
     * The partial residuals are used and the total residuals compared
     * each iteration to see if there is a sufficiently small change to
     * stop.  Otherwise, go for iteration limit.
     *
     */

    ap = work;
    dp = work + nt + 1;

    work  += 2*(nt+1);
    lwork -= 2*(nt+1);

    ZtZ = work;
    
    work  += D * D;
    lwork -= D * D;

    Rlast = work;
    work += nt;
    lwork -= nt;

      /* Set dp */

    profile_drift( nt, D, Basis, &p[Drift[0]], dp, work, lwork );

      /* ATTN: TEMP */

#if DIAGNOSTIC >= 20

    if( debugfp == NULL )
        debugfp = fopen( "/tmp/debug", "wb" );

    fprintf( stderr, "Preliminary Initialization\n\n" );

    lnpost_fixed( &neffpar, p, &v2 );
    fprintf( stderr, "Outcome: RSS = %lg, Post = %lg\n", 2.0*SSR, v2 );

    /* for( i = 0; i < neffpar; i++ )
        fprintf( stderr, "%lg\n", p[i] ); */

    profile_active( nt, p, ap, 0, Bells, Bell_Length );

    for( t = 0; t < nt; t++ )
        R[t] = p[Mu] + dp[t] + ap[t];

    fwrite( dp, sizeof(double), nt, debugfp );
    fwrite( R, sizeof(double), nt, debugfp );

    Dlnpost_fixed( &neffpar, p, Save );

    fprintf( stderr, "--------------------------\n\n" );

#endif

      /* END TEMP */

      /* Begin Iterations */

    for( i = 0; i < Init_Max_Iter; i++ )
    {
        /* Shape */

        if( i == 0 )
        {
            for( t = 0; t < nt; t++ )
            {
                R[t] = y[t] - p[Mu] - dp[t];
                Rlast[t] = R[t];  /* Initialize Rlast, ignoring ap term */
            }
        }
        else
        {
            for( t = 0; t < nt; t++ )
                R[t] = y[t] - p[Mu] - dp[t];
        }
        
        /* ATTN: If all resp 0, set to prior maximum.  Otherwise, maximize partial. */

        for( k = 0, anyresp = 0; k < Keff; k++ )      /* Find if all Resp pars are 0 */
            anyresp += (p[RespInd + k] > 0.0) ? 1 : 0;

        if( anyresp )
        {
            /* Set copy of parameters */

            for( j = 0; j < NParams; j++ )
                Save[j] = p[j];

            for( s = 0; s < Shape_Params; s++ )
                work[s] = p[Shape[s]];

            set_lnpshape( R, Save );

            /* Optimize */

            qnb_minimize( Shape_Params, Shape_Params, work, 
			  &Param_LB[Shape[0]], &Param_UB[Shape[0]],
			  &psc[Shape[0]], fsc, 
			  minimizerWork, minimizerWorkSize,
			  minimizerIWork, minimizerIWorkSize,
			  &value, &iter, 
			  &minimizer_evals, &minimizer_grad_tol,
                          lnpost_shape, Dlnpost_shape, errorStream );

            /* Store and Restore Changed Values */

            for( s = 0; s < Shape_Params; s++ )
                p[Shape[s]] = work[s];

            set_lnpshape( NULL, NULL );
        }
        else
        {
            /* Only prior impact on objective: use prior mode */
            
            for( s = 0; s < Shape_Params; s++ )
                p[Shape[s]] = (Shape_a[s] - 1.0)/Shape_b[s];
        }
        
        /* Responsiveness */

            /* Make S, S'R, S'S etc.; residuals are the same */

        for( s = 0; s < UniqueStims; s++ )  /* Set Bells for each unique stimulus length */
            (*poly_bell)(T, &Bells[s*2*T], &Bell_Length[s], p, Stims[s], Slice_Offset, IAI, StimShift[s]);

        make_active_mats( nt, Keff, p, S, R, StS, Sty, St1, Bells, Bell_Length );

            /* Solve system, adjust, and re-normalize */

        dsytrf( &Upper, &Keff, StS, &Keff, IWork, work, &lwork, &info );

        if( !info )
        {
            for( k = 0, u = 1.0/p[Mu]; k < Keff; k++ )
                Sty[k] = u * Sty[k] -  u * u * Resp_b/p[o_SigmaSq];
            
            dsytrs( &Upper, &Keff, &iOne, StS, &Keff, IWork, Sty, &Keff, 
		    &info );
        }

        if( info ) {
            Warning(1, "Could not solve regression equations during initialization (info = %d)", info );
	    return 0;
	}

        for( k = 0; k < Keff; k++ )  /* Assign Normalized and Adjusted Coefficients */
            if( Sty[k] < 0.0 )
                p[RespInd + k] = 0.0;
            else
                p[RespInd + k] = Sty[k];

            /* Compute Resulting Activation Profile */

        profile_active( nt, p, ap, 0, Bells, Bell_Length );

        /* Drift */

        for( d = 0; d < D; d++ )
        {
            p[Coefs[d]] = 0.0;

            for( j = 0; j < d; j++ )
                ZtZ[j + d*D] = ZtZ[d + j*D] = smp * Ptot[ (j+1) + (d+1)*(D+1) ];

            ZtZ[d + d*D] = 1.0 + smp * Ptot[ (d+1) + (d+1)*(D+1) ];
        }

        for( t = 0; t < nt; t++ )
        {
            R[t] = y[t] - p[Mu] - ap[t];

            for( d = 0; d < D; d++ )
                p[Coefs[d]] += Basis[d][t] * R[t];
        }

        dsytrf( &Upper, &D, ZtZ, &D, IWork, work, &lwork, &info );

        if( !info )
            dsytrs( &Upper, &D, &iOne, ZtZ, &D, IWork, &p[Coefs[0]], &D, &info );

        if( info ) {
            Warning(1, "Could not solve regression equations during initialization (info = %d)", info );
	    return 0;
	}

            /* Form profile */

        profile_drift( nt, D, Basis, &p[Drift[0]], dp, work, lwork );

        /* Baseline */

        mulast = p[Mu];

        for( t = 0, u2 = u = 0.0, v2 = 1.0/mulast; t < nt; t++ )
        {
            u += (y[t] - dp[t]) * (1.0 + ap[t]*v2);
            u2 += (1.0 + ap[t]*v2) * (1.0 + ap[t]*v2);
        }

        v2 = mulast - Mu_a;
        v = (2.0 * Mu_b * v2/(1.0  +  Mu_b * v2 * v2))/p[o_SigmaSq];

        p[Mu] = (u - v)/u2;    /* u/u2 */

        if( p[Mu] < 1 )      /* ATTN: TEMPORARY STABILIZER, prevents mu = 0 */
        {
            Warning( DIAGNOSTIC >= 1, "During initialization, baseline approaches 0." );
            p[Mu] = mulast;
            change_tol = 1.0e10;
            break;
        }

        /*
        p[Mu] = u/nt;

        for( k = 0; k < Keff; k++ )
            p[RespInd + k] *= mulast/p[Mu];
         */

        /* Noise and Smoothness */

        /* u = v - T * (p[Mu] - mulast)*(p[Mu] - mulast);  */ /* RSS when Mu = residual average above */

            /* NOTE: Also update and check Rlast for termination condition */

        for( t = 0, u = 0.0, u2 = 0.0; t < T; t++ )
        {
            ap[t] *= p[Mu]/mulast;                  /* Update actn profile */

            v = y[t] - p[Mu] - dp[t] - ap[t];       /* New Residual Value  */
            u += v*v;                               /* Accumulate RSS      */

            u2 += (Rlast[t] - v) * (Rlast[t] - v);  /* Change in Residual  */
            Rlast[t] = v;                           /* Store for next loop */
        }

        p[o_SigmaSq] = (0.5 * T + SigmaSq_a - 1.0)/(0.5 * u + SigmaSq_b);

        SmoothPar = (smp * p[o_SigmaSq]);           /* 1/lambda = target * r */

        change_tol = sqrt(u2/(1.0 + u));

        /* ATTN: TEMP */
#if DIAGNOSTIC >= 20

        fprintf( stderr, "Initialization Iteration %d\n\n", i );
        
        lnpost_fixed( &neffpar, p, &v2 );
        fprintf( stderr, "Outcome: RSS = %lg, Like = %lg, Post = %lg, ct = %lg\n",
                 u, 0.5 * u * p[o_SigmaSq] - 0.5 * T * Log(p[o_SigmaSq]), v2, change_tol );

        for( t = 0; t < nt; t++ )
            R[t] = p[Mu] + dp[t] + ap[t];

        fwrite( dp, sizeof(double), nt, debugfp );
        fwrite( R, sizeof(double), nt, debugfp );

        Dlnpost_fixed( &neffpar, p, Save );

        fprintf( stderr, "--------------------------\n\n" );

#endif

        /* END TEMP */

        if( change_tol <= Init_Tolerance )
            break;  /* Done */
    }

    if( i == Init_Max_Iter )
        change_tol = 1.0e10;


    /* Polish Initial Value */

    for( t = 0; t < nt; t++ )
        R[t] = y[t] - p[Mu] - dp[t];

        /* Responsiveness */

            /* Make S, S'R, S'S etc.; residuals are the same */

    for( s = 0; s < UniqueStims; s++ )  /* Set Bells for each unique stimulus length */
        (*poly_bell)(T, &Bells[s*2*T], &Bell_Length[s], p, Stims[s], Slice_Offset, IAI, StimShift[s]);

    make_active_mats( nt, Keff, p, S, R, StS, Sty, St1, Bells, Bell_Length );

            /* Solve system, adjust, and re-normalize */

    dsytrf( &Upper, &Keff, StS, &Keff, IWork, work, &lwork, &info );

    if( !info )
    {
        for( k = 0, u = 1.0/p[Mu]; k < Keff; k++ )
            Sty[k] = u * Sty[k] -  u * u * Resp_b/p[o_SigmaSq];
            
        dsytrs( &Upper, &Keff, &iOne, StS, &Keff, IWork, Sty, &Keff, &info );
    }

    if( info ) {
      Warning(1, "Could not solve regression equations during initialization (info = %d)", info );
      return 0;
    }

    for( k = 0; k < Keff; k++ )  /* Assign Normalized and Adjusted Coefficients */
        if( Sty[k] < 0.0 )
            p[RespInd + k] = 0.0;
        else
            p[RespInd + k] = Sty[k];

    profile_active( nt, p, ap, 0, Bells, Bell_Length );

        /* Noise and Smoothness */

    for( t = 0, u = 0.0, u2 = 0.0; t < T; t++ )
    {
        v = R[t] = y[t] - p[Mu] - dp[t] - ap[t];
        u += v*v;

        u2 += (Rlast[t] - v) * (Rlast[t] - v);
    }

    p[o_SigmaSq] = (0.5 * T + SigmaSq_a - 1.0)/(0.5 * u + SigmaSq_b);
    SmoothPar = (smp * p[o_SigmaSq]);           /* 1/lambda = target * r */

    if( sqrt( u2/(1.0 + u) ) > 1.25 * change_tol )
        Warning( DIAGNOSTIC >= 2, "Parameter refinement caused non-trivial change" );

    /* ATTN: TEMP */

#if DIAGNOSTIC >= 20

    fprintf( stderr, "Initialization Refinement:\n\n" );
        

    lnpost_fixed( &neffpar, p, &v2 );
    fprintf( stderr, "Outcome: RSS = %lg, Like = %lg, Post = %lg, ct = %lg\n",
             u, 0.5 * u * p[o_SigmaSq] - 0.5 * T * Log(p[o_SigmaSq]), v2, sqrt(u2/(1.0 + u)) );

#if DIAGNOSTIC >= 20 
    fprintf( stderr, "Varying r...\n" );
    for( i = 0, u = p[o_SigmaSq], v = u/2.0; i <= 1000; i++, v += u/1000.0 )
    {
        p[o_SigmaSq] = v;
        lnpost_fixed( &neffpar, p, &v2 );
        fprintf( stderr, "%lg %.16lg\n", v, v2 );
    }
    p[o_SigmaSq] = u;
    fprintf( stderr, "...Done\n" );
    
    fprintf( stderr, "Varying RespInd+1...\n" );
    for( i = 0, u = p[RespInd + 1], v = u/2.0; i <= 1000; i++, v += u/1000.0 )
    {
        p[RespInd + 1] = v;
        lnpost_fixed( &neffpar, p, &v2 );
        fprintf( stderr, "%lg %.16lg\n", v, v2 );
    }
    p[RespInd + 1] = u;
    fprintf( stderr, "...Done\n" );

    fprintf( stderr, "Varying Shape[0]...\n" );
    for( i = 0, u = p[Shape[0]], v = u/2.0; i <= 1000; i++, v += u/1000.0 )
    {
        p[Shape[0]] = v;
        lnpost_fixed( &neffpar, p, &v2 );
        fprintf( stderr, "%lg %.16lg\n", v, v2 );
    }
    p[Shape[0]] = u;
    fprintf( stderr, "...Done\n" );
#endif

    for( t = 0; t < nt; t++ )
        R[t] = p[Mu] + dp[t] + ap[t];

    fwrite( dp, sizeof(double), nt, debugfp );
    fwrite( R, sizeof(double), nt, debugfp );

    Dlnpost_fixed( &neffpar, p, Save );

    fprintf( stderr, "--------------------------\n\n" );

#endif

    /* END TEMP */

    return 1; /* successful initialization */
}

void   set_inits( int nt, double *vox, int npar, int neffpar, double *p, int *init,
                   double *initdat, double *work, int lwork )
{
    int           d, i, j, k, m, s, t;
    int           lag, nc, info = 0;

    double        u,  v;
    double       *tau, *scr;
    

    /*
     * Find decent initial Point
     *
     * The strategy here is to do a standard analysis for those
     * parameters that are not directly initialized by the user.
     *
     * When the baseline, responsiveness, drift, and precision are not
     * user-initialized we do a careful default procedure that uses
     * the initial knots and shape parameters.  If some of these are
     * user-initialized, we do a subset of this analysis.
     *
     * This initialization procedure is given below in the comments.
     *
     * Note that a negative value of Knots_Fixed implies that the knots are
     * fixed for this minimization but are not in general.  In this case,
     * we need to build the basis components.  In any case, we must build
     * the drift basis corresping to the initial knots.
     * 
     */

    if( Knots_Fixed <= 0 )         /* Knots Free or fixed in a transient way */
    {
        if( init[DRIFT_KNOTS] )
        {
            j = Init_offsets[DRIFT_KNOTS];
        
            for( i = 0; i < Dk; i++ ) 
                p[Knots[i]] = initdat[j + i];
        }
        else                         /* Prior Mean for knot positions      */
        {
            for( i = 0; i < Dk; i++ ) 
                p[Knots[i]] = (i + 1.0)/(Dk + 1.0);
        }

        if( Knots_Fixed )
        {
            qsort( &p[Knots[0]], Dk, sizeof(double), Numerically_d );

            make_drift_basis( T, Dg, Dk, &Basis[-1], &p[Knots[0]], basisR, Work, WorkSize );
            make_drift_qforms( T, Dg, Dk, &p[Knots[0]], Pnorm, Pcurv, NULL, &p[Knots[0]], Work, WorkSize );

            /* Set Pcomb = a Pnorm + b Pcurv */

            m = 1 + Dg + Dk;

            for( i = 0; i < m; i++ )       
            {
                for( j = 0; j < i; j++ )
                    Pcomb[i + j*m] = Pcomb[j + i*m] = (Drift_a * Pnorm[i + j*m] + Drift_b * Pcurv[i + j*m]);

                Pcomb[i + i*m] = (Drift_a * Pnorm[i + i*m] + Drift_b * Pcurv[i + i*m]);
            }
        }
        else
            make_drift_basis( T, Dg, Dk, &Basis[-1], &p[Knots[0]], basisR, Work, WorkSize );
    }
    else
    {
        for( i = 0; i < Dk; i++ ) 
            p[Knots[i]] = Fixed_Knots[i];
    }

    
    if( init[SHAPE] )
    {
        j = Init_offsets[SHAPE];

        for( s = 0; s < TOTAL_SHAPE_PARAMS; s++ )
            if( Shape_Use[s] )
                p[Shape[s]] = initdat[j + s];
    }
    else                         /* Prior Mode for shape parameters */
    {
        for( s = 0; s < TOTAL_SHAPE_PARAMS; s++ )
            if( Shape_Use[s] )
                p[Shape[s]] = (Shape_a[s] - 1.0)/Shape_b[s];
    }


    if( init[DRIFT_COEFS] )
    {
        j = Init_offsets[DRIFT_COEFS];
        
        for( d = 0; d < D; d++ )
            p[Coefs[d]]= initdat[j + d];
    }


    if( init[RESP] )
    {
        j = Init_offsets[RESP];
        
        for( k = 0; k < K; k++ )
            p[Resp[k]] = (Constrained[k]) ? 0.0 : initdat[j + k];
    }

    if( !init[BASELINE] || !init[DRIFT_COEFS] || !init[RESP] || !init[O_SIGMASQ] )
    {
        if( !init[DRIFT_COEFS] && !init[RESP] )
        {
            /*
             * Regress lagged stepped functions (constant on condition blocks,
             * equal within conditions) and trend basis on data.  The minimum
             * condition mean becomes the initial baseline, the responsiveness is
             * initialized by (condition mean - initial baseline)/initial baseline,
             * and the precision by the residual sum of squares.  The drift coefficients
             * give the initial drift profile.  Note that there is a question of
             * whether to regularize here, at the moment I have opted not to because
             * the number of knots will typically be sufficiently small and regularly
             * spaced initially.
             *
             * NOTE:  If StimulusStart[j+1] < StimulusEnd[j] for any j, this will
             *        not work optimally.
             *
             */

            lag = ceil( p[Shape[LAG_ON]]/IAI );

            for( t = 0; t < nt; t++ )
            {
                work[t] = vox[t];                            /* Data                    */

                for( j = 0; j < K; j++ )                     /* Lagged condition blocks */
                    work[t + (j + 1)*nt] = CondMatrix[j][Max(t - lag,0)];

                for( j = 0; j < D; j++ )                     /* Drift basis (no const)  */
                    work[t + (j + K + 1)*nt] = Basis[j][t];
            }

            nc = K + D;

            tau = work + (nc + 1)*nt;
            scr = work + (nc + 1)*(nt + 1);
            lwork -= (nc + 1)*(nt + 1);

                /* QR Decomposition */

            dgeqrf( &nt, &nc, work + nt, &nt, tau, scr, &lwork, &info );

            if( info )
                Warning( DIAGNOSTIC >= 1, "Could not compute QR decomposition in voxel_init (info = %d)",
                         info );
            

                /* Compute Q'y */
            
            dormqr( &Left, &Trans, &nt, &iOne, &nc, work + nt, &nt, tau, work, &nt, scr, &lwork, &info );

            if( info )
                Warning( DIAGNOSTIC >= 1, "Could not compute Q'y in voxel_init (info = %d)", info );
            
            /* Compute RSS by sum squares of work[nc..(nt-1)] */

            for( i = nc, u = 0.0; i < nt; i++ )
                u += work[i] * work[i];

            /* Compute Coefs by solving triangular system: (Q'y) = R work[0..(nc-1)] */
            
            dtrtrs( &Upper, &NoTrans, &NoUnit, &nc, &iOne, work + nt, &nt, work, &nc, &info );

            if( info )
                Warning( DIAGNOSTIC >= 1, "Could not compute coefs in voxel_init (using dtrtrs, info=%d)",
                         info );

            /* Now assign coefficients */

            for( k = 0, i = -1, v = BigNumber; k < K; k++ )
                if( work[k] < v )
                {
                    v = work[k];
                    i = k;
                }

            if( i < 0 || v <= 0.0 )
                Warning( DIAGNOSTIC >= 1, "Unexpected problem computing initial baseline, min = %lg", v );

            p[Mu] = v;

            for( k = 0; k < K; k++ )
            {
                if( Constrained[k] )
                    p[Resp[k]] = 0.0;
                else if( k == i )
                    p[Resp[k]] = 0.0;
                else
                    p[Resp[k]] = (work[k] - v)/v;
            }
            
            u /= (nt - nc);

            p[o_SigmaSq] = 1.0/u;

            for( d = 0; d < D; d++ )
                p[Coefs[d]] = work[K + d];
        }
        else if( !init[RESP] )
        {
            Abort( "This intitialization pattern not yet supported." );
        }
        else if( !init[DRIFT_COEFS] )
        {
            Abort( "This intitialization pattern not yet supported." );
        }
        else
        {
            Abort( "This intitialization pattern not yet supported." );
        }
    }
    
    /* If the baseline and precision numbers are user-initialized, over-write them */

    if( init[BASELINE] )
        p[Mu] = initdat[Init_offsets[BASELINE]];
    
    if( init[O_SIGMASQ] )
        p[o_SigmaSq] = initdat[Init_offsets[O_SIGMASQ]];
}

void   voxel_fit( int nt, double *vox, int npar, int neffpar, double *p, double *val,
                  int *iter, int *nfunc, double *psc, double fsc, int *skip,
                  double *guess, int lguess, int method, void (*func)( int *, double *, double * ),
                  void (*Dfunc)( int *, double *, double * ) )
{
    int           i, j, k, s, info = 0;
    int minimizer_evals;
    double minimizer_grad_tol;
    
    /* Optimize         */

    if( method == 1 )  /* Smooth Optimization---Fast  */
    {

#if !defined( FIX_SIGMA )

        qnb_minimize( neffpar, npar, p, Param_LB, Param_UB,
		      psc, fsc, 
		      minimizerWork, minimizerWorkSize,
		      minimizerIWork, minimizerIWorkSize,
		      val, iter, &minimizer_evals, &minimizer_grad_tol,
                      func, Dfunc, errorStream );

#else

#define SIGMA_REFINE       3    /* ATTN: SIGMA_REFINE is TEMPORARY and PROVISIONAL */

        for( i = 0; i < SIGMA_REFINE; i++ ) 
        {
	    qnb_minimize( neffpar, npar, p, Param_LB, Param_UB,
			  psc, fsc, 
			  minimizerWork, minimizerWorkSize,
			  minimizerIWork, minimizerIWorkSize,
			  val, iter, &minimizer_evals, &minimizer_grad_tol,
			  func, Dfunc, errorStream );

            for( t = 0, u = 0.0; t < T; t++ )
                u += R[t] * R[t];

            v = p[o_SigmaSq];

            p[o_SigmaSq] = (0.5 * T + SigmaSq_a - 1.0)/(0.5 * u + SigmaSq_b);

            if( fabs(v - p[o_SigmaSq]) < Init_Tolerance * 0.5 * (v + p[o_SigmaSq]) )
                break;
        }

#endif
    }
    else               /* General Optimization---Slow */
    {
        *iter = Iter_Max;

        minimize_nmb( neffpar, npar, p, val, iter, Tolerance, Param_LB, Param_UB,
                      psc, fsc, guess, lguess, func, Dfunc, &info );

        if( info )
            Message( "  %s\n", min_errors_nmb(info) );
    }

    *nfunc = guess[0];
    
    /* Synchronize internal quantities, esp. R */

    (*func)( &neffpar, p, val );

    /* ATTN: TEMP */

#if DIAGNOSTIC >= 20

    fprintf( stderr, "Final Result\n\n" );
    fprintf( stderr, "Outcome: RSS = %lg, Post = %lg, Iter = %d\n\n\n", 2.0*SSR, *val, *iter );

    profile_drift( T, D, Basis, &p[Drift[0]], Drift_prof, Work, WorkSize );
    profile_active( T, p, Active_prof, 1, Bells, Bell_Length );

    for( i = 0; i < T; i++ )
        R[i] = p[Mu] + Drift_prof[i] + Active_prof[i];

    fwrite( Drift_prof, sizeof(double), T, debugfp );
    fwrite( R, sizeof(double), T, debugfp );
    
    Dlnpost_fixed( &neffpar, p, Save );

#if DIAGNOSTIC >= 20
    {
        int    k;
        double u,  v;
        double u2, v2;

        fprintf( stderr, "Varying r...\n" );
        for( i = 0, u = p[o_SigmaSq], v = u/2.0; i <= 1000; i++, v += u/1000.0 )
        {
            p[o_SigmaSq] = v;
            lnpost_fixed( &neffpar, p, &v2 );
            fprintf( stderr, "%lg %.16lg\n", v, v2 );
        }
        p[o_SigmaSq] = u;
        fprintf( stderr, "...Done\n" );
    
        fprintf( stderr, "Varying RespInd+1...\n" );
        for( i = 0, u = p[RespInd + 1], v = u/2.0; i <= 1000; i++, v += u/1000.0 )
        {
            p[RespInd + 1] = v;
            lnpost_fixed( &neffpar, p, &v2 );
            fprintf( stderr, "%lg %.16lg\n", v, v2 );
        }
        p[RespInd + 1] = u;
        fprintf( stderr, "...Done\n" );

        fprintf( stderr, "Varying Shape[0]...\n" );
        for( i = 0, u = p[Shape[0]], v = u/2.0; i <= 1000; i++, v += u/1000.0 )
        {
            p[Shape[0]] = v;
            lnpost_fixed( &neffpar, p, &v2 );
            fprintf( stderr, "%lg %.16lg\n", v, v2 );
        }
        p[Shape[0]] = u;
        fprintf( stderr, "...Done\n" );
    }
#endif

    fprintf( stderr, "--------------------------\n\n" );

#endif

    /* END TEMP */

    (*func)( &neffpar, p, val );

    /* Indicate which parameters are at their bounds */

    for( i = 0; i < neffpar; i++ )
    {
        if( p[i] <= Param_LB[i] )
            skip[i] = -1;
        else if( p[i] >= Param_UB[i] )
            skip[i] = 1;
        else
            skip[i] = 0;
    }

    /* ATTN: If Shape_a == 1, then the hessian is 0 when all gamma's are zero */
    /* For the moment, mark the corresponding parameters to be skipped        */
    /* Generally avoid Shape_a == 1, evidently.  This will be revised.        */

    for( s = 0; s < Shape_Params; s++ )
    {
        if( Shape_a[s] == 1.0 && Keff )
        {
            for( k = j = 1; k < Keff; k++ )
                j *= skip[RespInd + k] * skip[RespInd + k];

            if( j ) /* All Resp's marked to skip */
                skip[Shape[s]] = -1;
        }
    }
}


void     compute_cov( int neff, int n, double *p, double *obsinfo, int method,
                      int *skip, double *work, int worksize, int *rank, double *lndet, int *pinfo,
                      int null_model )
{
    int       i, j, info = 0;
    int       ns, firstskip = -1, bigskip = 0;
    int       inda, indb;

    double    u;
    
    /** Compute Hessian at Optimum **/

    if( method == COV_ANALYTIC )              /* Analytic Form */
    {
        if( null_model )
            (*DDlnpost_null)( neff, n, p, obsinfo, skip, work, worksize );
        else
        {
            /* Set initial scale for shape params (while they are done numerically) */

            for( i = 0; i < Shape_Params; i++ )
                Hinit[Shape[i]] = SEpsc[Shape[i]];

            /* Compute Hessian */
            
            (*DDlnpost)( neff, n, p, obsinfo, skip, work, worksize );
        }
    }
    else if( method & COV_FCD_MASK )          /* First (Central) Differences using Gradient  */
    {
        /* Set initial scale */

        for( i = 0; i < n; i++ )
            Hinit[i] = SEpsc[i];

        /* First Central Difference Hessian Approximation */

        if( null_model )
            hessian_fcd( neff, n, p, Param_LB, Param_UB, obsinfo, method ^ COV_FCD_MASK, Hinit,
                         skip, -1.0, 0.0, Scd_Tolerance, work, worksize, Dlnpost_null, &info );
        else
            hessian_fcd( neff, n, p, Param_LB, Param_UB, obsinfo, method ^ COV_FCD_MASK, Hinit,
                         skip, -1.0, 0.0, Scd_Tolerance, work, worksize, Dlnpost, &info );
    }
    else                            /* Second (Central) Differences using Gradient */
    {
        /* Set initial scale */

        for( i = 0; i < n; i++ )
            Hinit[i] = SEpsc[i];


        /* Second Central Difference Hessian Approximation */

        if( null_model )
            hessian_scd( neff, n, p, Param_LB, Param_UB, obsinfo, method, Hinit, skip, -1.0, 0.0,
                         Scd_Tolerance, work, worksize, lnpost_null, &info );
        else
            hessian_scd( neff, n, p, Param_LB, Param_UB, obsinfo, method, Hinit, skip, -1.0, 0.0,
                         Scd_Tolerance, work, worksize, lnpost, &info );
    }

    /** Invert Hessian to Estimate Covariance Matrix **/

    ns = neff;  /* Size of info matrix with skipped vars removed */

    /* Check for skipped variables and reduce matrix size */

    for( i = 0; i < neff; i++ )
    {
        if( skip[i] )
        {
            ns--;

            if( firstskip < 0 )
            {
                firstskip = i;

#if DIAGNOSTIC>=2
                Warning( DIAGNOSTIC >= 2, "  Variables skipped in Hessian computation" );
                Message( "           Indices: " );
#endif    
            }

            if( skip[i] > 1 )
                bigskip++;

#if DIAGNOSTIC>=2
            Message( "%d ", i );
#endif
        }
    }

#if DIAGNOSTIC >= 1
    if( ns < neff )
        Message( "\n" );
#endif

    /* Check Hessian Status */

    if( info > 0 )
        Warning( DIAGNOSTIC >= 1, "Could not compute hessian accurately (code = %d).", info );
    else if( info < 0 )
    {
        Warning( DIAGNOSTIC >= 1,"  Possible inaccuracy in computed hessian (code=%d, index=%d,%d).\n",
                 (-info)%16, ((-info)/16)%neff, ((-info)/16)/neff );

        if( bigskip )
            Warning( DIAGNOSTIC >= 1, "  %d variables removed in hessian computation.\n", bigskip );
    }

    /*
     * Analytic Hessian computes all entries
     * If variables are being skipped, need to reduce matrix before inverting.
     * Later expansion will be the same as in the numerical case.
     *
     */

    if( (method == COV_ANALYTIC) && (ns < neff) )
        reduce_hessian( neff, obsinfo, skip, firstskip );


    /* Invert information matrix if possible */

    dsptrf( &Upper, &ns, obsinfo, IWork, &info );

    if( info  )
        Warning( DIAGNOSTIC >= 1, "Could not decompose hessian (info = %d).", info );

    dsptri( &Upper, &ns, obsinfo, IWork, work, &info );
    
    if( info )
    {
        Warning( DIAGNOSTIC >= 1, "Could not invert hessian (info = %d).", info );

        for( i = 0, j = (ns * (ns + 1))/2; i < j; i++ )
            obsinfo[i] = 0.0;

        *pinfo = info;
        *rank = 0;         /* Flag a big problem */
        *lndet = 0.0;

        return;
    }

        /* Compute Log Determinant of Covariance Matrix  */

    if( Null_Model )  /* ATTN: Comparing_Models */
    {
        for( i = 0, j = (ns * (ns + 1))/2; i < j; i++ )
            work[i] = obsinfo[i];

        dpptrf( &Upper, &ns, work, &info );

        if( info )
        {
            Warning( DIAGNOSTIC >= 1, "Covariance has non-positive leading %d minor, could not factor.",
                     info );

            u = 0.0;    /* ATTN: Better way of flagging this?  This shouldn't happen though. */
        }
        else
            for( i = 0, u = 0.0; i < ns; i++ )
                u += 2.0 * Log( fabs(work[i + (i*(i+1))/2]) );  /* fabs for roundoff, 2 for squaring */
    }

    if( ns < neff )
    {
        /*
         * Here, some variables are on or very close to their bounds
         * so the numerical second derivatives will not be accurate
         * because the standard errors will be on the order of the estimates.
         * We thus remove these variables from the hessian computation and
         * fold them back in afterwards.  Similarly, if there is a catastrophic
         * inaccuracy in the computation of a derivative, that variable has skip
         * set to a large positive value; these also are skipped.
         *
         * Expand hessian matrix to put in skipped variables.
         * Unfortunately, this necessitates some copying.
         *
         */

        expand_hessian( neff, ns, obsinfo, skip, p, firstskip, work );
    }

    *pinfo = info;
    *rank = ns;
    *lndet = u;

    return;
}


void     reduce_hessian( int neff, double *obsinfo, int *skip, int firstskip )
{
    int        i, j, k;
    int        inda, indb, indlimit;

    /*
     * Trusting user to give accurate first skip index.
     * A value of 0 is always safe as a default.
     * If the value is inappopriate or the given index is not skipped,
     * the index is recomputed.  (That's why 0 is safe either way.)
     *
     * The reason for trusting the user is to take advantage of prior
     * computation (i.e., performance).
     *
     */

    if( firstskip < 0 || !skip[firstskip] )  /* Recompute first skip index */
    {
        for( firstskip = 0; firstskip < neff; firstskip++ )
            if( skip[firstskip] )
                break;

        if( firstskip == neff )
            return;
    }

    /* Contract matrix to remove skipped variables */

    inda = (firstskip * (firstskip + 1))/2;
    indlimit = (neff * (neff + 1))/2;

    /* Find next non-skipped column */

    for( j = firstskip + 1; j < neff; j++ )
        if( !skip[j] )
            break;

    if( j == neff )      /* All done; zero out remainder of matrix */
    {
        while( inda < indlimit )
            obsinfo[inda++] = 0.0;

        return;
    }
        
    /*
     * Copy values from indb to inda skipping as we go
     *
     * Here, i will be the column currently being copied into
     *       j will be the column currently being copied out of
     *
     *       k is a loop index variable
     *
     */

    for( i = firstskip; i < neff; i++ )
    {
        indb = (j * (j + 1))/2;   /* First row of next non-skip column */

        /* inda is advanced the right number of times in this loop     */
        /* because all of the columns  (and rows) between [i and j)    */
        /* are skipped.  Hence, there are exactly i non-skipped values */
        /* in the jth column of the matrix.                            */

        for( k = 0; k <= j; k++, indb++ )  
            if( !skip[k] )
                obsinfo[inda++] = obsinfo[indb];        

        for( j++; j < neff; j++ )
            if( !skip[j] )
                break;

        if( j == neff )      /* All done; zero out remainder of matrix */
        {
            while( inda < indlimit )
                obsinfo[inda++] = 0.0;

            return;
        }
    }
}


void    expand_hessian( int neff, int ns, double *obsinfo,
                        int *skip, double *diag, int firstskip, double *work )
{
    int        i, j;
    int        inda, indb;

    if( firstskip < 0 )
        firstskip = 0;     /* Zero is acceptable default as it works with the whole matrix */

    /* Expand matrix to put in skipped variables */

    inda = (firstskip * (firstskip + 1))/2;

    for( i = inda, j = (ns * (ns + 1))/2; i < j; i++ )
        work[i] = obsinfo[i];

    for( i = firstskip, indb = inda; i < neff; i++ )
    {
        if( skip[i] )
        {
            for( j = 0; j < i; j++ )
                obsinfo[inda++] = 0.0;

            obsinfo[inda++] = 2.0 * fabs(diag[i]);
        }
        else
        {
            for( j = 0; j <= i; j++ )
            {
                if( skip[j] )
                    obsinfo[inda++] = 0.0;
                else
                    obsinfo[inda++] = work[indb++];
            }
        }
    }
}


/*** Posterior Computation Code ***/


/*
 * LNPOST_FIXED evaluates the NEGATIVE of the log-likelihood at the given point.
 * Assumes that the drift knots (and thus the drift basis) are fixed.
 *
 * Side Effects: Sets the values of R and SSR to the residual vector and
 *               residual sum of squares respectively.
 *
 */

void   lnpost_fixed( int *npar, double *p, double *val )
{
    int      i, j, m;
    int      d, k, s, t;
    int      lwork, info = 0;

    double   c, u, v;
    double   lnr, ll;
    double   dscaling;

    /* Compute and store profiles  */

    profile_drift( T, D, Basis, &p[Drift[0]], Drift_prof, Work, WorkSize );  /* ATTN: Work */
    profile_active( T, p, Active_prof, 1, Bells, Bell_Length );

    /* Construct Residuals and RSS */

    for( t = 0, SSR = 0.0; t < T; t++ )
    {
        R[t] = u = Y[t] - p[Mu] - Drift_prof[t] - Active_prof[t];
        SSR += u*u;
    }

    SSR *= 0.5;


    /* Likelihood Term */

    lnr = Log(p[o_SigmaSq]);

    ll = SSR * p[o_SigmaSq] - 0.5 * T * lnr;
    
    /* Prior Terms */
                                                          /* Baseline */

    ll += Log( 1.0 + Mu_b * (p[Mu] - Mu_a)*(p[Mu] - Mu_a) );

                                                          /* Drift    */

    dscaling = (Drift_scaled) ? (1.0/(p[Mu] * p[Mu])) : 1.0;

    /* Compute 0.5 delta' Ptot delta = 0.5 delta' R^-1' Pcomb R^-1 delta */

    for( d = 0, m = D + 1, u = 0.0; d < D; d++ )       
    {
        v = p[Coefs[d]];
        
        for( j = 0; j < d; j++ )
            u += v * p[Coefs[j]] * Ptot[(j + 1) + (d + 1)*m];

        u += 0.5 * v * v * Ptot[(d + 1) * (m + 1)];
    }

    ll += dscaling * SmoothPar * u; 
    /* ll += dscaling * (Smooth_Target * p[o_SigmaSq]) * u; *//* ATTN */
    

    for( k = 0; k < Keff; k++ )                           /* Responsiveness */
        ll += Resp_b * p[RespInd + k];
 
    for( s = 0; s < Shape_Params; s++ )                   /* Shape */
    {
        ll += Shape_b[s]*p[Shape[s]] - (Shape_a[s] - 1.0) * Log(p[Shape[s]]);
    }

                                                          /* o_SigmaSq */

    ll += SigmaSq_b * p[o_SigmaSq] - (SigmaSq_a - 1.0) * lnr;
       /* 0.5 * (lnr + SigmaSq_a)*(lnr + SigmaSq_a)*o_SigmaSq_bsq - lnr; */

    *val = ll;

    /* ATTN: TEMP */
#if DIAGNOSTIC >= 10
    fprintf( stderr, "\n" );
    fprintf( stderr, "* ll = %.16lg, mu = %.8lg, r = %.8lg\n",
             ll, p[Mu], 100.0*p[o_SigmaSq] );
    fprintf( stderr, "*              resp = (%.8lg, %.8lg, %.8lg)\n",
             100.0*p[RespInd + 1], 100.0*p[RespInd + 2], 100.0*p[RespInd + 3] );
    fprintf( stderr, "*              shape = (%.8lg, %.8lg, %.8lg, %.8lg)\n",
             p[Shape[0]], p[Shape[1]], p[Shape[2]], p[Shape[3]] );
#endif
    /* END TEMP */
}

void   Dlnpost_fixed( int *npar, double *p, double *deriv )
{
    int      b, d, i, j, k, l, m, s, t;
    int      info;

    double   v, u, r;
    double   hgt, dscaling;

    double   *df;

    /* Compute Profiles and Profile Derivatives */

    profile_drift( T, D, Basis, &p[Drift[0]], Drift_prof, Work, WorkSize );
    profile_active( T, p, Active_prof, 1, Bells, Bell_Length );
    (*shape_deriv)( p, DBells, Stims, Slice_Offset, IAI, StimShift );

    /* Construct Residuals and RSS */

    for( t = 0, SSR = 0.0; t < T; t++ )
    {
        R[t] = u = Y[t] - p[Mu] - Drift_prof[t] - Active_prof[t];
        SSR += u*u;
    }

    SSR *= 0.5;

    /* Set derivative components one at a time */

    r = p[o_SigmaSq];

        /* Baseline */

    for( t = 0, u = 1.0/p[Mu], v = 0.0; t < T; t++ )
        v += R[t] * (1.0 + u * Active_prof[t]);

    u = p[Mu] - Mu_a;                     

    deriv[Mu] = -v * r  +                                 /* Likelihood Term */
                (2.0 * Mu_b * u/(1.0  +  Mu_b * u * u));  /* Prior Term      */


        /* Drift */

            /* Prepare Prior Term First:  Compute Ptot delta = (R^-1' P R^-1) delta */

    dscaling = (Drift_scaled) ? (1.0/(p[Mu] * p[Mu])) : 1.0;
    m = D + 1;

    dgemv( &NoTrans, &D, &D, &One, Ptot + m + 1, &m, &p[Coefs[0]], &iOne, &Zero, Work, &iOne );

            /* Now combine with Likelihood Term */

    for( d = 0; d < D; d++ )
    {
        for( t = 0, df = Basis[d], v = 0.0; t < T; t++ )
            v += R[t] * df[t];

        deriv[Coefs[d]] = -v * r  +  dscaling * SmoothPar * Work[d]; /* Likelihood + Prior Term */
        /* deriv[Coefs[d]] = -v * r  +  dscaling * (Smooth_Target * p[o_SigmaSq]) * Work[d];*/  /* ATTN */

    }

        /* Responsiveness */

    for( k = 0, u = p[Mu]; k < Keff; k++ )
    {
        for( i = CondBlock[RespEff[k]][0], j = 1, v = 0.0; j <= i; j++ )
        {
            b = CondBlock[RespEff[k]][j];
            s = WhichStim[b];
            d = 2 * T * s;        /* 2T is col length in Bells */
            
            for( t = BlockStart[b], m = 0; m < Bell_Length[s] && t < T; t++, m++ )
                v += R[t] * Bells[m + d];
        }

        deriv[RespInd + k] = -v * r * u + Resp_b;  /* Likelihood + Prior Term */  
    }


        /* Shape */    

    for( s = 0; s < Shape_Params; s++ )
        deriv[Shape[s]] = (Shape_b[s] - (Shape_a[s] - 1.0)/p[Shape[s]]);   /* Prior Term      */

    for( b = 0, u = p[Mu]; b < B; b++ )
    {
        k = BlockStart[b];
        j = WhichStim[b];
        i = Cond[b];

        if( !Constrained[i] && (p[Resp[i]] > 0.0) )
        {
            l = Min( Bell_Length[j], T - k - 1 );

            for( s = 0; s < Shape_Params; s++ )
            {
                for( i = 0, df = DBells[s], v = 0.0; i <= l; i++ )
                    v += R[i + k] * df[i + j*2*T];

                deriv[Shape[s]] += -v * r * u * p[Resp[Cond[b]]];    /* Likelihood Term */
            }
        }
    }
            

        /* o_SigmaSq */

#if !defined( FIX_SIGMA )

    deriv[o_SigmaSq] = (SSR - 0.5 * T / r) +                            /* Likelihood Term */
                       (SigmaSq_b - (SigmaSq_a - 1.0)/r);               /* Prior Term      */
                    /* ((Log(r) + SigmaSq_a) * o_SigmaSq_bsq - 1.0)/r; */

#else

    deriv[o_SigmaSq] = 0.0; 

#endif

    /* ATTN: TEMP */

#if DIAGNOSTIC >= 10
    {
        int    k;
        double u2, v2;

        fprintf( stderr, "D(Mu) = %lg\n", deriv[Mu] );
        fprintf( stderr, "D(r) = %lg\n", (SSR - 0.5 * T / r) + (SigmaSq_b - (SigmaSq_a - 1.0)/r) );/*deriv[o_SigmaSq] );*/

        for( k = 1, v2 = fabs(deriv[Coefs[0]]); k < D; k++ )
        {
            u2 = fabs(deriv[Coefs[k]]);

            if( u2 > v2 )
                v2 = u2;
        }

        fprintf( stderr, "D(d0) = %lg\n", v2 );
    
        for( k = 0; k < Shape_Params; k++ )
            fprintf( stderr, "D(Shape[%d]) = %lg\n", k, deriv[Shape[k]] );
    
        for( k = 0; k < Keff; k++ )
            fprintf( stderr, "D(RespEff[%d]) = %lg\n", k, deriv[RespInd + k] );
    }
#endif

    /* END TEMP */

}


/*
 * DDLNPOST_FIXED
 *   Computes Analytic Second-Derivatives for fixed knot posterior
 *   Works with any model specification
 *
 * The Hessian is stored in LAPACK symmetric packed format. 
 * All entries can be computed analytically except for the shape
 * sub-matrix.  Actually, this can as well but it is a real pain;
 * for now, use numerical approach on that small block.  This should
 * still be reasonably efficient and not scale with the problem size.
 *
 * WARNING:  This function takes several liberties with the ordering
 *           of the parameter vector.  The HESSIAN macro below assumes
 *           that the first argument is <= the second.  It could check
 *           but at a performance hit.  Since the orderings of the
 *           parameter blocks are stable  (Mu Drift o_SigmaSq Resp Shape),
 *           this should be ok, but be careful.
 *
 */
 
#define  HESSIAN(i,j)    (obsinfo[(i) + ((j)*((j)+1))/2])  /* Temporary, i <= j */

void   DDlnpost_fixed( int neff, int n, double *p, double *obsinfo, int *skip, double *work, int lwork )
{
    int      i, j, m, t;
    int      info, stim;
    int      active, anyresp = 0, bsplines = 0;

    int      b, b2;    /* Indices for Block, Drift, Resp, Shape */
    int      d, d2;
    int      k, k2;
    int      s, s2;

    double   u, v, r;
    double   dscale, omu;
    double   res_t, ap_t;

    double   *pr;
    double   *respd, *shaped;
    double   *shape_hess;


    /** Prepare Basic Quantities **/

    active = Shape_Params * Keff;    /* Non-zero if active components in model */

    if( active )
    {
        for( k = 0, anyresp = 0; k < Keff; k++ )         /* Find if all Resp pars are 0  */
            anyresp += (p[RespInd + k] > 0.0) ? 1 : 0;   /* ATTN: Non-negativity assumed */
    }
    

    v = p[Mu];
    r = p[o_SigmaSq];
    m = D + 1;

    dscale = (Drift_scaled) ? (1.0/(p[Mu] * p[Mu])) : 1.0;
    dscale *= SmoothPar; /*(r * Smooth_Target);*/   /* SmoothPar; */ /* ATTN */


    /** Initialize Hessian **/

        /* Set self-second-derivatives contributed by priors */
        /* and zero all of the off-diagonal terms            */

    u = p[Mu] - Mu_a;
    
    HESSIAN(Mu,Mu) = 2.0 * Mu_b * (1.0 - Mu_b* u * u)/( (1.0 + Mu_b * u * u) * (1.0 + Mu_b * u * u) );

    for( i = Mu+1; i < neff; i++ )   /* Zero out the rest of the row */
        HESSIAN(Mu,i) = 0.0;


    for( d = 0; d < D; d++ )
    {
        i = Coefs[d];
        
        for(j = Coefs[D-1] + 1; j < neff; j++ )    /* Zero out rest of the row */ /* ATTN: Order Dependence */
            HESSIAN(i,j) = 0.0;
    }

    HESSIAN(o_SigmaSq,o_SigmaSq) = (0.5 * T  + SigmaSq_a - 1.0) / (r * r);

    for( i = o_SigmaSq+1; i < neff; i++ )  /* Zero out rest of the row */
        HESSIAN(o_SigmaSq,i) = 0.0;

    if( active )
    {
        for( k = 0; k < Keff; k++ )
            for( j = RespInd; j < neff; j++ )
                HESSIAN( RespInd, j ) = 0.0;    /* All these prior terms are zero + rest of row */

        if( !anyresp )
        {
            /* Shape prior term is (Shape_a[s] - 1.0)/(p[Shape[s]]*p[Shape[s]]) for Shape_a != 1 */
            /* and 0 for Shape_a == 1.  Without Resp, no impact on likelihood, so when a > 1 the */
            /* maximum is at (a - 1)/b, and the only the diagonal elements in the Shape part are */
            /* non-zero.  Hence, the value is b^2/(a-1).                                         */

            for( s = 0; s < Shape_Params; s++ )
            {
                i = Shape[s];

                if( Shape_a[s] == 1.0 )     /* ATTN: Shape Prior Settings; Shape should be skipped in this case */
                    HESSIAN( i, i ) = 0.0;
                else
                    HESSIAN( i, i ) = (Shape_b[s]*Shape_b[s])/(Shape_a[s] - 1.0);

                for( j = i +1; j < neff; j++ )    /* Zero out rest of row */
                    HESSIAN( i, j ) = 0.0;
            }
        }
        else
        {
            /* Prepare for adding in shape component; zero out entire part of matrix */

            for( s = 0; s < Shape_Params; s++ )
            {
                i = Shape[s];
            
                for( j = i; j < neff; j++ )
                    HESSIAN( i, j ) = 0.0;
            }
        }
    }

    /* ATTN: If all of the responsiveness parameters are zero (currently this means on the boundary), */
    /* then this assumes that they will be skipped in the above calculations and ignores them here.   */
    /* This speeds the computations involving the other parameters, particularly the shape parameters */
    /* which can now be set analytically rather than numerically.                                     */
    /* The more general alternative is to use anyresp within active but for now go for the speed.     */

    if( active && !anyresp )
        active = 0;         /* Remaining Computations ignore Shape and Resp */

    /** Compute Drift Profiles **/ 

    profile_drift( T, D, Basis, &p[Drift[0]], Drift_prof, work, lwork );

    /** Compute Active Profile and Derivatives **/

    if( active )
    {
        profile_active( T, p, Active_prof, 1, Bells, Bell_Length );
        (*shape_deriv)( p, DBells, Stims, Slice_Offset, IAI, StimShift );

        /* Prepare some temporary quantities for later */

        pr = work;
        work += T;
        lwork -= T;

        respd = work;
        work += Keff;
        lwork -= Keff;

        shaped = work;
        work += Shape_Params;
        lwork -= Shape_Params;

        shape_hess = work;
        work += Shape_Params * Shape_Params;
        lwork -= Shape_Params * Shape_Params;
    }

    /** Compute Hessian **/

        /* Set All Components including Likelihood Diagonal Contribution */

    if( active )
    {
        for( t = 0, b = 0, omu = 1.0/p[Mu]; t < T; t++ )
        {
            /* Construct Residuals, Partial Residuals, and other useful quantities */

            pr[t] = Y[t] - p[Mu] - Drift_prof[t];     /* Partial Residual at Shape */
            res_t = pr[t] - Active_prof[t];           /* Full Residual             */
            ap_t  = 1.0 + omu * Active_prof[t];       /* Normalized Act Profile    */

            /*
             * Find Resp and Shape blocks that impact this time
             *
             * Here, b represents the earliest block whose bell
             * can overlap the current time, and b2 represents
             * the blocks >= b over which we search for a contribution.
             * This search continues until the starting time is > t.
             *
             * Over this range of blocks, we accumulate the values of the
             * Bells or DBells for each Resp and Shape parameter.  This,
             * by linearity, is the term that contributes to the derivatives.
             *
             */

            while( BlockStart[b] + Bell_Length[WhichStim[b]] < t && b < B )
                b++;

            for( k = 0; k < Keff; k++ )
                respd[k] = 0.0;

            for( s = 0; s < Shape_Params; s++ )
                shaped[s] = 0.0;

            for( b2 = b; b2 < B && BlockStart[b2] <= t; b2++ )
            {
                stim = WhichStim[b2];

                if( !Constrained[Cond[b2]] )
                {
                    k = Resp[Cond[b2]] - RespInd;  /* Index from 0..Keff-1 */
                
                    respd[k] += v * Bells[t - BlockStart[b2] + 2*T*stim];
                }
                
                for( s = 0; s < Shape_Params; s++ )
                    shaped[s] += v * p[Resp[Cond[b2]]] * DBells[s][t - BlockStart[b2] + 2*T*stim];
            }
                     
            /* Set Hessian values */

            HESSIAN(Mu,Mu) += r * ap_t * ap_t;

            for( d = 0; d < D; d++ )
            {
                i = Coefs[d];
                u = Basis[d][t];

                HESSIAN(Mu,i) += r * u * ap_t;
                HESSIAN(i,o_SigmaSq) += -res_t * u;

                for( k = 0; k < Keff; k++ )
                    HESSIAN(i,RespInd+k) += r * u * respd[k];

                for( s = 0; s < Shape_Params; s++ )
                    HESSIAN(i,Shape[s]) += r * u * shaped[s];
            }
            
                    
            HESSIAN(Mu,o_SigmaSq) += -res_t * ap_t;

            for( k = 0; k < Keff; k++ )
            {
                HESSIAN(Mu,RespInd+k) += r * respd[k] * (ap_t - res_t/v);

                HESSIAN(o_SigmaSq, RespInd+k) += -res_t * respd[k];
                
                for( k2 = 0; k2 <= k; k2++ )  /* Second term in self-derivative is 0 */
                    HESSIAN(RespInd+k2,RespInd+k) += r * respd[k] * respd[k2];

                for( s = 0; s < Shape_Params; s++ )
                    HESSIAN(Shape[s],RespInd+k) += r * respd[k] * shaped[s];
            }
            
            for( s = 0; s < Shape_Params; s++ )
            {
                HESSIAN(Mu,Shape[s]) += r * shaped[s] * (ap_t - res_t/v);

                HESSIAN(o_SigmaSq, Shape[s]) += -res_t * shaped[s];
            }
        }

        /* Compute Coef sub-matrix directly */

        if( !bsplines )
        {
            for( d = 0; d < D; d++ )
            {
                i = Coefs[d];

                for( d2 = 0; d2 < d; d2++ )
                    HESSIAN( Coefs[d2], i ) = dscale * Ptot[(d2+1) + (d+1)*m];
                
                HESSIAN( i, i ) = r + dscale * Ptot[(d+1)*(m+1)];

            }
        }
        else
        {
            /* B-spline matrices are Banded with bandwidth Dg + 1 */

            Abort( "B-spline Hessian not yet supported" );
        }

        /* Compute the Shape sub-matrix numerically */

        /* ATTN: Eventually do this analytically    */
        /* Start by adding in the first derivative  */
        /* parts of the shape derivatives above     */

        set_lnpshape( pr, p );

        for( s = 0; s < Shape_Params; s++ )
            shaped[s] = p[Shape[s]];

        hessian_scd( Shape_Params, Shape_Params, shaped, &Param_LB[Shape[0]], &Param_UB[Shape[0]],
                     shape_hess, COV_SCD_RICHARDSON, &Hinit[Shape[0]], &skip[Shape[0]], -1.0, 0.0,
                     Scd_Tolerance, work, lwork, lnpost_shape, &info );

        set_lnpshape( NULL, NULL );

        for( s = 0; s < Shape_Params; s++ )   /* Copy packed shape hessian into Hessian matrix */
        {
            i = Shape[s];

            for( s2 = 0; s2 < s; s2++ )
                HESSIAN(Shape[s2],i) = shape_hess[s2 + (s*(s+1))/2];

            HESSIAN( i, i ) = shape_hess[s + (s*(s+1))/2];
        }
    }
    else
    {
        HESSIAN(Mu,Mu) += T * r;

        for( t = 0; t < T; t++ )
        {
            /* Compute Residuals and accumulate RSS */

            res_t = Y[t] - p[Mu] - Drift_prof[t];

            /* Compute Hessian terms */

            HESSIAN(Mu,o_SigmaSq) += -res_t;

            for( d = 0; d < D; d++ )
                HESSIAN(Coefs[d],o_SigmaSq) += -res_t * Basis[d][t];
        }
        
        if( !bsplines )
        {
            for( d = 0; d < D; d++ )
            {
                HESSIAN( Mu, Coefs[d] ) = 0.0;

                for( d2 = 0; d2 < d; d2++ )
                    HESSIAN( Coefs[d2], Coefs[d] ) = dscale * Ptot[(d2+1) + (d+1)*m];
                
                HESSIAN( Coefs[d], Coefs[d] ) = r + dscale * Ptot[(d+1)*(m+1)];

            }
        }
        else  /* ATTN BSPLINES */
        {
            /* B-spline matrices are Banded with bandwidth Dg + 1 */

            assert( 0 );
        }
    }
}

#undef HESSIAN

void   lnpost_variable( int *npar, double *p, double *val )
{
}

void   Dlnpost_variable( int *npar, double *p, double *deriv )
{
}


void   lnpost_null_fixed( int *npar, double *p, double *val )
{
    int      i, j, m;
    int      d, k, s, t;
    int      lwork, info = 0;

    double   c, u, v;
    double   lnr, ll;
    double   dscaling;

    /* Compute and Store Drift Profile */

    profile_drift( T, D, Basis, &p[Drift[0]], Drift_prof, Work, WorkSize );

    /* Construct Residuals and RSS */

    for( t = 0, SSR = 0.0; t < T; t++ )
    {
        R[t] = u = Y[t] - p[Mu] - Drift_prof[t];
        SSR += u*u;
    }

    SSR *= 0.5;


    /* Likelihood Term */

    lnr = Log(p[o_SigmaSq]);

    ll = SSR * p[o_SigmaSq] - 0.5 * T * lnr;
    
    /* Prior Terms */
                                                          /* Baseline */

    ll += Log( 1.0 + Mu_b * (p[Mu] - Mu_a)*(p[Mu] - Mu_a) );

                                                          /* Drift    */

    dscaling = (Drift_scaled) ? (1.0/(p[Mu] * p[Mu])) : 1.0;

    /* Compute 0.5 delta' Ptot delta = 0.5 delta' R^-1' Pcomb R^-1 delta */

    for( d = 0, m = D + 1, u = 0.0; d < D; d++ )       
    {
        v = p[Coefs[d]];
        
        for( j = 0; j < d; j++ )
            u += v * p[Coefs[j]] * Ptot[(j + 1) + (d + 1)*m];

        u += 0.5 * v * v * Ptot[(d + 1) * (m + 1)];
    }

    ll += dscaling * SmoothPar * u;

                                                          /* o_SigmaSq */

    ll += SigmaSq_b * p[o_SigmaSq] - (SigmaSq_a - 1.0) * lnr;
       /* 0.5 * (lnr + SigmaSq_a)*(lnr + SigmaSq_a)*o_SigmaSq_bsq - lnr; */

    *val = ll;
}

void   Dlnpost_null_fixed( int *npar, double *p, double *deriv )
{
    int      b, d, i, j, k, l, m, s, t;
    int      info;

    double   v, u, r;
    double   hgt, dscaling;

    double   *df;

    /* Compute and Store Drift Profile */

    profile_drift( T, D, Basis, &p[Drift[0]], Drift_prof, Work, WorkSize );

    /* Construct Residuals and RSS */

    for( t = 0, SSR = 0.0; t < T; t++ )
    {
        R[t] = u = Y[t] - p[Mu] - Drift_prof[t];
        SSR += u*u;
    }

    SSR *= 0.5;

    /* Set derivative components one at a time */

    r = p[o_SigmaSq];

        /* Baseline */

    for( t = 0, v = 0.0; t < T; t++ )
        v += R[t];

    u = p[Mu] - Mu_a;                     

    deriv[Mu] = -v * r  +                                 /* Likelihood Term */
                (2.0 * Mu_b * u/(1.0  +  Mu_b * u * u));  /* Prior Term      */


        /* Drift */

            /* Prepare Prior Term First:  Compute Ptot delta = (R^-1' P R^-1) delta */

    dscaling = (Drift_scaled) ? (1.0/(p[Mu] * p[Mu])) : 1.0;
    m = D + 1;

    dgemv( &NoTrans, &D, &D, &One, Ptot + m + 1, &m, &p[Coefs[0]], &iOne, &Zero, Work, &iOne );

            /* Now combine with Likelihood Term */

    for( d = 0; d < D; d++ )
    {
        for( t = 0, df = Basis[d], v = 0.0; t < T; t++ )
            v += R[t] * df[t];

        deriv[Coefs[d]] = -v * r  +  dscaling * SmoothPar * Work[d];  /* Likelihood + Prior Term */

    }

        /* o_SigmaSq */

    deriv[o_SigmaSq] = (SSR - 0.5 * T / r) +                         /* Likelihood Term */
                       (SigmaSq_b - (SigmaSq_a - 1.0)/r);            /* Prior Term      */
                    /* ((Log(r) + SigmaSq_a) * o_SigmaSq_bsq - 1.0)/r; */
}

void   DDlnpost_null_fixed( int neff, int n, double *p, double *obsinfo, int *skip,
                            double *work, int lwork )
{
    int  k;

    k = Keff;
    Keff = 0;

    DDlnpost_fixed( neff, n, p, obsinfo, skip, work, lwork );

    Keff = k;
}

void   lnpost_null_variable( int *npar, double *p, double *val )
{
}

void   Dlnpost_null_variable( int *npar, double *p, double *deriv )
{
}


/*
 * LNPOST_SHAPE evaluates the NEGATIVE of the log-likelihood at the given point
 * varying only the shape parameters.
 *
 * DLNPOST_SHAPE evaluates the derivative of this function with respect to the
 * shape parameters at the given point.
 *
 * The static variables lnpshape_p and lnpshape_y hold pointers to the
 * current settings of the other (fixed) parameters and the data used in the posterior
 * computation.  These are set by the function   set_lnpshape( double *, double * ) below.
 *
 * It also assumes that the drift knots (and thus the drift basis) are fixed.
 *
 * SIDE EFFECTS
 *
 * The shape components of the vector pointed to by lnpshape_p are changed
 * but restored before the function exits;
 *
 * The variables Bells, Bell_Length, DBells and Active_prof are overwritten
 * by these functions.
 *
 */

static double     *lnpshape_p = NULL;
static double     *lnpshape_y = NULL;

void   lnpost_shape( int *nshape, double *shape, double *val )
{
    int      i, j, m;
    int      d, k, s, t;
    int      lwork, info = 0;

    double   c, u, v;
    double   lnr, ll;
    double   ssr;
    double   save[TOTAL_SHAPE_PARAMS];

    /* Compute and store profiles  */

    for( s = 0; s < Shape_Params; s++ )
    {
        save[s] = lnpshape_p[Shape[s]];
        lnpshape_p[Shape[s]] = shape[s];
    }

    profile_active( T, lnpshape_p, Active_prof, 1, Bells, Bell_Length );

    /* Construct Residuals and RSS */

    for( t = 0, ssr = 0.0; t < T; t++ )
    {
        u = lnpshape_y[t] - Active_prof[t];
        ssr += u*u;
    }

    ssr *= 0.5;


    /* Likelihood Term */

    lnr = Log(lnpshape_p[o_SigmaSq]);

    ll = ssr * lnpshape_p[o_SigmaSq] - 0.5 * T * lnr;
    
    /* Prior Terms */

    for( s = 0; s < Shape_Params; s++ )                   /* Shape */
    {
        ll += Shape_b[s]*shape[s] - (Shape_a[s] - 1.0) * Log(shape[s]);
    }

    /* Restore state variables */

    for( s = 0; s < Shape_Params; s++ )
        lnpshape_p[Shape[s]] = save[s];

    /* Record objective value */

    *val = ll;
}

void   Dlnpost_shape( int *nshape, double *shape, double *deriv )
{
    int      b, d, i, j, k, l, m, s, t;
    int      info;

    double   v, u, r, ssr, hgt;

    double   *df;
    double   save[TOTAL_SHAPE_PARAMS];

    /* Compute Profiles and Profile Derivatives */

    for( s = 0; s < Shape_Params; s++ )
    {
        save[s] = lnpshape_p[Shape[s]];
        lnpshape_p[Shape[s]] = shape[s];
    }

    profile_active( T, lnpshape_p, Active_prof, 1, Bells, Bell_Length );
    (*shape_deriv)( lnpshape_p, DBells, Stims, Slice_Offset, IAI, StimShift );

    /* Construct Residuals and RSS */

    for( t = 0, ssr = 0.0; t < T; t++ )
    {
       u = lnpshape_y[t] - Active_prof[t];
       ssr += u*u;
    }

    ssr *= 0.5;

    /* Set derivative components one at a time */

    r = lnpshape_p[o_SigmaSq];

        /* Shape */    

    for( s = 0; s < Shape_Params; s++ )
        deriv[s] = (Shape_b[s] - (Shape_a[s] - 1.0)/shape[s]);   /* Prior Term      */

    for( b = 0, u = lnpshape_p[Mu]; b < B; b++ )
    {
        k = BlockStart[b];
        j = WhichStim[b];
        i = Cond[b];

        hgt = lnpshape_p[Resp[i]];

        if( !Constrained[i] && (hgt > 0.0) )
        {
            l = Min( Bell_Length[j], T - k - 1 );

            for( s = 0; s < Shape_Params; s++ )
            {
                for( i = 0, df = DBells[s], v = 0.0; i <= l; i++ )
                    v += (lnpshape_y[i + k] - Active_prof[i + k]) * df[i + j*2*T];

                deriv[s] += -v * r * u * hgt;    /* Likelihood Term */
            }
        }
    }

    /* Restore state variables */

    for( s = 0; s < Shape_Params; s++ )
        lnpshape_p[Shape[s]] = save[s];


    ssr *= 1.0;    /* ATTN: TEMP  (for debugger) */
}

void    set_lnpshape( double *data, double *pars )
{
    lnpshape_y = data;
    lnpshape_p = pars;
}


/* Profiles: Activation and Drift */

void   profile_active( int nt, double *p, double *ac_prof, int newbells, double *bells, int *bell_length )
{
    int       b, i, j, k, l;
    double    len, hgt;

    /* Set Bells for each unique stimulus length */

    if( newbells )
    {
        for( l = 0; l < UniqueStims; l++ )
            (*poly_bell)( T, &bells[l*2*T], &bell_length[l], p, Stims[l], Slice_Offset, IAI, StimShift[l] );
    }
    

    /* Set first block then initialize profile */

    i = 0;

    if( !Constrained[Cond[0]] && ((hgt = p[Mu] * p[Resp[Cond[0]]]) > 0.0) )
    {
        j = WhichStim[0];
        l = Min( bell_length[j], T - 1 );

        for( ; i <= l; i++ )
            ac_prof[i] = hgt * bells[i + j*2*T];
    }
    
    for( ; i < T; i++ )
        ac_prof[i] = 0.0;

    /* Set remaining blocks */

    for( b = 1; b < B; b++ )
    {
        k = BlockStart[b];
        j = WhichStim[b];
        i = Cond[b];

        if( !Constrained[i] && ((hgt = p[Mu] * p[Resp[i]]) > 0.0) )
        {
            l = Min( bell_length[j], T - k - 1 );

            for( i = 0; i <= l; i++ )
                ac_prof[i + k] += hgt * bells[i + j*2*T];
        }
    }
}


void   profile_active2( int nt, double *p, double *apnew, int newbells, double *bells, int *bell_length )
{
    int       b, i, j, k, l, m, n, t;
    int       next;

    double    hgt;

    /* Set Bells for each unique stimulus length */

    if( newbells )
        for( l = 0; l < UniqueStims; l++ )
            (*poly_bell)( T, &bells[l*2*T], &bell_length[l], p, Stims[l], Slice_Offset, IAI, StimShift[l] );

    /* Set first bell or zero out first block if constrained */

    next = BlockStart[1];

    if( !Constrained[Cond[0]] && ((hgt = p[Mu] * p[Resp[Cond[0]]]) > 0.0) )
    {
        j = WhichStim[0];
        l = Min( bell_length[j], T );
        n = Max( l, next );
        
        for( t = 0; t < n; t++ )
            apnew[t] = hgt * bells[t + j*2*T];
    }
    else
    {
        for( t = 0; t < next; t++ )
            apnew[t] = 0.0;
    }

    /* Set remaining blocks */

    for( b = 1; b < B; b++ )
    {
        k = next;
        next = BlockStart[b+1];
        j = WhichStim[b];
        i = Cond[b];

        if( !Constrained[i] && ((hgt = p[Mu] * p[Resp[i]]) > 0.0) )
        {
            l = Min( bell_length[j], T - k );
            n = Min( l, t - k );

            for( m = 0, i = k; m < n; m++, i++ )
                apnew[i] += hgt * bells[m + j*2*T];

            for( ; m < l; m++, t++ )
                apnew[t] = hgt * bells[m + j*2*T];
        }
        else
        {
            for( ; t < next; t++ )
                apnew[t] = 0.0;
        }
    }

    /* ATTN: Any follow-up required here??? */
}

#if !defined(NO_INLINE_DRIFTPROF)

    /*
     * Assume that Basis has been just set to contain basis without constant term.
     *
     * Now, construct profile.  (Do not include constant term)
     *
     */

#define  profile_drift( nt, d, basis, coefs, prof, work, lwork ) \
    dgemv( &NoTrans, &(nt), &(d), &One, (basis)[0], &(nt), (coefs), &iOne, &Zero, (prof), &iOne )

#else

void     profile_drift( int nt, int d, double **basis, double *coefs, double *prof, double *work, int lwork )
{
    
    /*
     * Assume that Basis has been just set to contain basis without constant term.
     *
     * Now, construct profile.  (Do not include constant term)
     *
     */

    dgemv( &NoTrans, &nt, &d, &One, basis[0], &nt, coefs, &iOne, &Zero, prof, &iOne );
}

#endif


void   profiles( double *p, double *dr_prof, double *ac_prof, double *work, int lwork )/* DEFUNCT */
{
    int       i, j, k, l;
    int       b;

    double    len, hgt;

    /*
     * Drift Profile
     *
     * Need to construct basis for current drift knots and second derivative (wrt time) profile
     *
     * Assumes that Mu and Drift[0] are contiguous.
     *
     */

    profile_drift( T, D, Basis, &p[Drift[0]], dr_prof, work, lwork );

    /* Generate Activation Profile */

      /* Set Bells for each unique stimulus length */
      /* Do not compute bell for fixed blocks      */

    for( i = 0; i < UniqueStims; i++ )
        (*poly_bell)( T, &Bells[i*2*T], &Bell_Length[i], p, Stims[i], Slice_Offset, IAI, StimShift[i] );
    
      /* Combine over blocks */

          /* Set first block then initialize profile */

    i = 0;

    if( !Constrained[Cond[0]] && ((hgt = p[Mu] * p[Resp[Cond[0]]]) > 0.0) )
    {
        j = WhichStim[0];
        l = Bell_Length[j];

        for( ; i <= l; i++ )
            ac_prof[i] = hgt * Bells[i + j*2*T];
    }
    
    for( ; i < T; i++ )
        ac_prof[i] = 0.0;

          /* Set remaining blocks */

    for( b = 1; b < B; b++ )
    {
        k = BlockStart[b];
        j = WhichStim[b];
        i = Cond[b];

        if( !Constrained[i] && ((hgt = p[Mu] * p[Resp[i]]) > 0.0) )
        {
            l = Min( Bell_Length[j], T - k - 1 );

            for( i = 0; i <= l; i++ )
                ac_prof[i + k] += hgt * Bells[i + j*2*T];
        }
    }
}



/*
 * POLY_BELL: Aligned version; 4, 6, and 8 parameter models
 *
 * Computes a bell function corresponding to a given stimulus length,
 * slight offset, and inter-acquisition interval.  Assumes that the
 * condition beginnings are aligned with the image acquisitions.
 *
 * The bells are pre-computed on a fine grid in [0,1] but must be interpolated
 * onto the image acquisition grid.  ...
 *
 * Assuming that the image acquisition grid is aligned with the stimulus presentation
 * times (i.e., every stimulus presentation lies on an image acquisition time).
 * When this assumption does not hold, unaligned profiles should be used
 * (see profiles_unaligned).  Given that the grids are aligned, we assume that
 * the fine grid is chosen to fit some multiple into a single inter-image space.
 * 
 * NOTE: ABOVE COMMENT APPEARS OBSOLETE, HANDLES ALIGNED AND UNALIGNED CASES
 *
 * The pre-computed functions are AttackRamp and DecayRamp.  Both are vectors
 * of length RampLen + 1; they sample an attack and decay ramp between 0 and 1
 * over a fine regular grid over [0,1].  (Actually a bit longer and padded with
 * 0's and 1's as appropriate to prevent overrun.)
 *
 * Eventually this will include a dip and varying rise/fall rates.
 *
 * WARNING:  Care must be taken that the rounding schemes used here are
 *           portable.  Use of floor() throughout would greatly increase
 *           the expense, but this approach should be handled carefully.
 *
 *
 */

void  poly_bell4( long t, double *b, int *blen, double *p, double stimlen,
                  double offset, double IAI, double shift )
{
    int       i, ja, jd;
    int       len = 0;

    double    u, v;
    double    a1, a2, d1, d2;

    double    lag_on, lag_off;
    double    attack, decay;

    double    awgt, dwgt;    /* Interpolation Weights       */
    double    astart, aend;  /* Beginning and End of Attack */
    double    dstart, dend;  /* Beginning and End of Decay  */

    /* Store shape parameters for later use */

    lag_on = p[Shape[LAG_ON]] + shift;
    lag_off = p[Shape[LAG_OFF]];
    attack = p[Shape[ATTACK]];
    decay = p[Shape[DECAY]];

    /*
     * Compute b = AB * DB
     *
     * Compute AB and DB piecewise so that we loop over the vector only
     * once, multiplying them together as we go.
     *
     */

        /* Mark important boundaries */

    astart = Min( (lag_on)/IAI, t );
    aend   = Min( (lag_on + attack)/IAI, t );
    dstart = Min( (stimlen + lag_off)/IAI, t );
    dend   = Min( (stimlen + lag_off + decay)/IAI, t );

        /* Pre-compute values that will be needed */

    a1 = RampLen * (offset - lag_on)/attack;
    a2 = RampLen * IAI/attack;
    d1 = RampLen * (offset - stimlen - lag_off)/decay;
    d2 = RampLen * IAI/decay;

        /* Create the bell function */

    for( i = 0; i < astart; i++ )
        b[i] = 0.0;

    if( dstart < aend )
    {
        for( ; i < dstart; i++ )
        {
            u = a1 + i * a2;
            ja = u;           /* Want effect of floor(u) here */
            awgt = u - ja;
         
            b[i] = (1.0 - awgt) * AttackRamp[ja]  +  awgt * AttackRamp[ja + 1];
        }

        if( aend < dend )
            for( ; i <= aend; i++ )
            {
                u = a1 + i * a2;
                ja = u;           /* Want effect of floor(u) here */
                awgt = u - ja;
         
                v = d1 + i * d2;
                jd = v;           /* Want effect of floor(v) here */
                dwgt = v - jd;

                b[i] = ((1.0 - awgt) * AttackRamp[ja]  +  awgt * AttackRamp[ja + 1]) *
                       ((1.0 - dwgt) * DecayRamp[jd]   +  dwgt * DecayRamp[jd + 1]);
            }
        else
            for( ; i <= dend; i++ )
            {
                u = a1 + i * a2;
                ja = u;           /* Want effect of floor(u) here */
                awgt = u - ja;
         
                v = d1 + i * d2;
                jd = v;           /* Want effect of floor(v) here */
                dwgt = v - jd;

                b[i] = ((1.0 - awgt) * AttackRamp[ja]  +  awgt * AttackRamp[ja + 1]) *
                       ((1.0 - dwgt) * DecayRamp[jd]   +  dwgt * DecayRamp[jd + 1]);
            }
    }
    else
    {
        for( ; i < aend; i++ )
        {
            u = a1 + i * a2;
            ja = u;           /* Want effect of floor(u) here */
            awgt = u - ja;
         
            b[i] = (1.0 - awgt) * AttackRamp[ja]  +  awgt * AttackRamp[ja + 1];
        }

        for( ; i < dstart; i++ )
            b[i] = 1.0;
    }
    

    for( ; i < dend; i++ )
    {
        v = d1 + i * d2;
        jd = v;           /* Want effect of floor(v) here */
        dwgt = v - jd;

        b[i] = (1.0 - dwgt) * DecayRamp[jd]  +  dwgt * DecayRamp[jd + 1];
    }

    b[i] = 0.0;
    *blen = i;  /* No need for Min( i, t ) since Min's done above for dend */
}


void  poly_bell6( long t, double *b, int *blen, double *p, double stimlen,
                  double offset, double IAI, double shift )
{
}

void  poly_bell8( long t, double *b, int *blen, double *p, double stimlen,
                  double offset, double IAI, double shift )
{
}


/*
 * SHAPE_DERIV:  Aligned version; 4, 6, or 8 parameter models
 *
 * Computes derivative bells corresponding to each shape parameter
 * for every unique stimulus length.
 *
 *
 */

void      shape_deriv_4( double *p, double **dbells, double *stims,
                         double offset, double IAI, double *shift )
{
    int       i, ja, jd, k;
    int       len = 0, ind;

    double    u, v;
    double    a1, a2, d1, d2;

    double    stimlen;

    double    lag_on, lag_off;
    double    attack, decay;
    double    o_na, o_nasq, o_nd, o_ndsq;

    double    awgt, dwgt;    /* Interpolation Weights       */
    double    astart, aend;  /* Beginning and End of Attack */
    double    dstart, dend;  /* Beginning and End of Decay  */

    /* Store shape parameters for later use */

    lag_off = p[Shape[LAG_OFF]];
    attack = p[Shape[ATTACK]];
    decay = p[Shape[DECAY]];

    o_na = -1.0/attack;
    o_nasq = -IAI/(attack * attack);
    o_nd = -1.0/decay;
    o_ndsq = -IAI/(decay * decay);


    /*
     * Compute AB and DB and their derivatives piecewise so
     * that we loop over the vector only once, multiplying them together as we go.
     * The derivatives for each of the parameters is computed in the inner loop
     * to avoid repeating the overhead of re-computing the interpolation.
     *
     * See poly_bell for discussion of this algorithm.
     *
     */

    for( k = 0; k < UniqueStims; k++ )
    {
        lag_on = p[Shape[LAG_ON]] + shift[k];
        stimlen = stims[k];
        ind = 2*k*T;

        /* Mark important boundaries */

        astart = Min( (lag_on)/IAI, T );
        aend   = Min( (lag_on + attack)/IAI, T );
        dstart = Min( (stimlen + lag_off)/IAI, T );
        dend   = Min( (stimlen + lag_off + decay)/IAI, T );

        /* Pre-compute values that will be needed */

        a1 = RampLen * (offset - lag_on)/attack;
        a2 = RampLen * IAI/attack;
        d1 = RampLen * (offset - stimlen - lag_off)/decay;
        d2 = RampLen * IAI/decay;

        /* Create the bell function */

        for( i = 0; i < astart; i++ )
        {
            dbells[LAG_ON][i + ind] = 0.0;
            dbells[ATTACK][i + ind] = 0.0;
            dbells[LAG_OFF][i + ind] = 0.0;
            dbells[DECAY][i + ind] = 0.0;
        }
        
        if( dstart < aend )
        {
            for( ; i < dstart; i++ )
            {
                u = a1 + i * a2;
                ja = u;           /* Want effect of floor(u) here */
                awgt = u - ja;

                dbells[LAG_ON][i + ind] = o_na * ((1.0 - awgt)*DAttackRamp[ja] + awgt*DAttackRamp[ja + 1]);
                dbells[ATTACK][i + ind] = o_nasq * (i - astart + offset) *
                                                 ((1.0 - awgt)*DAttackRamp[ja] + awgt*DAttackRamp[ja + 1]);
                dbells[LAG_OFF][i + ind] = 0.0;
                dbells[DECAY][i + ind] = 0.0;
            }

            if( aend < dend )
                for( ; i <= aend; i++ )
                {
                    u = a1 + i * a2;
                    ja = u;           /* Want effect of floor(u) here */
                    awgt = u - ja;
         
                    v = d1 + i * d2;
                    jd = v;           /* Want effect of floor(v) here */
                    dwgt = v - jd;

                    dbells[LAG_ON][i + ind] = o_na *
                                                 ((1.0 - awgt)*DAttackRamp[ja] + awgt*DAttackRamp[ja + 1]) *
                                                 ((1.0 - dwgt)*DecayRamp[jd]   + dwgt*DecayRamp[jd + 1]);
                    dbells[ATTACK][i + ind] = o_nasq * (i - astart + offset) *
                                                 ((1.0 - awgt)*DAttackRamp[ja] + awgt*DAttackRamp[ja + 1]) *
                                                 ((1.0 - dwgt)*DecayRamp[jd]   + dwgt*DecayRamp[jd + 1]);
                    dbells[LAG_OFF][i + ind] = o_nd *
                                                 ((1.0 - awgt)*AttackRamp[ja] + awgt*AttackRamp[ja + 1]) *
                                                 ((1.0 - dwgt)*DDecayRamp[jd]   + dwgt*DDecayRamp[jd + 1]);
                    dbells[DECAY][i + ind] = o_ndsq * (i - dstart + offset) *
                                                 ((1.0 - awgt)*AttackRamp[ja] + awgt*AttackRamp[ja + 1]) *
                                                 ((1.0 - dwgt)*DDecayRamp[jd]   + dwgt*DDecayRamp[jd + 1]);
                }
            else
                for( ; i <= dend; i++ )
                {
                    u = a1 + i * a2;
                    ja = u;           /* Want effect of floor(u) here */
                    awgt = u - ja;
         
                    v = d1 + i * d2;
                    jd = v;           /* Want effect of floor(v) here */
                    dwgt = v - jd;

                    dbells[LAG_ON][i + ind] = o_na *
                                                 ((1.0 - awgt)*DAttackRamp[ja] + awgt*DAttackRamp[ja + 1]) *
                                                 ((1.0 - dwgt)*DecayRamp[jd]   + dwgt*DecayRamp[jd + 1]);
                    dbells[ATTACK][i + ind] = o_nasq * (i - astart + offset) *
                                                 ((1.0 - awgt)*DAttackRamp[ja] + awgt*DAttackRamp[ja + 1]) *
                                                 ((1.0 - dwgt)*DecayRamp[jd]   + dwgt*DecayRamp[jd + 1]);
                    dbells[LAG_OFF][i + ind] = o_nd *
                                                 ((1.0 - awgt)*AttackRamp[ja] + awgt*AttackRamp[ja + 1]) *
                                                 ((1.0 - dwgt)*DDecayRamp[jd]   + dwgt*DDecayRamp[jd + 1]);
                    dbells[DECAY][i + ind] = o_ndsq * (i - dstart + offset) *
                                                 ((1.0 - awgt)*AttackRamp[ja] + awgt*AttackRamp[ja + 1]) *
                                                 ((1.0 - dwgt)*DDecayRamp[jd]   + dwgt*DDecayRamp[jd + 1]);
                }
        }
        else
        {
            for( ; i < aend; i++ )
            {
                u = a1 + i * a2;
                ja = u;           /* Want effect of floor(u) here */
                awgt = u - ja;
         
                dbells[LAG_ON][i + ind] = o_na * ((1.0 - awgt)*DAttackRamp[ja] + awgt*DAttackRamp[ja + 1]);
                dbells[ATTACK][i + ind] = o_nasq * (i - astart + offset) *
                                                 ((1.0 - awgt)*DAttackRamp[ja] + awgt*DAttackRamp[ja + 1]);
                dbells[LAG_OFF][i + ind] = 0.0;
                dbells[DECAY][i + ind] = 0.0;
            }

            for( ; i < dstart; i++ )
            {
                dbells[LAG_ON][i + ind] = 0.0;
                dbells[ATTACK][i + ind] = 0.0;
                dbells[LAG_OFF][i + ind] = 0.0;
                dbells[DECAY][i + ind] = 0.0;
            }
        }
    

        for( ; i < dend; i++ )
        {
            v = d1 + i * d2;
            jd = v;           /* Want effect of floor(v) here */
            dwgt = v - jd;

            dbells[LAG_ON][i + ind] = 0.0;
            dbells[ATTACK][i + ind] = 0.0;
            dbells[LAG_OFF][i + ind] = o_nd * ((1.0 - dwgt)*DDecayRamp[jd] + dwgt*DDecayRamp[jd + 1]);
            dbells[DECAY][i + ind] = o_ndsq * (i - dstart + offset) *
                                                 ((1.0 - dwgt)*DDecayRamp[jd] + dwgt*DDecayRamp[jd + 1]);
        }

        dbells[LAG_ON][i + ind] = 0.0;
        dbells[ATTACK][i + ind] = 0.0;
        dbells[LAG_OFF][i + ind] = 0.0;
        dbells[DECAY][i + ind] = 0.0;
    }
}

void      shape_deriv_6( double *p, double **dbells, double *stims, double offset, double IAI, double *shift )
{
}

void      shape_deriv_8( double *p, double **dbells, double *stims, double offset, double IAI, double *shift )
{
}

void      make_active_mats( int T, int Keff, double *p, double *S, double *y, double *StS,
                            double *Sty, double *St1, double *bells, int *bell_len )
{
    int        b, d, i, j, k, kless, ktrue, l, m, t, w;
    int        blk_cnt, last;

    double     u, v, sum, doty;
    
    int       *cb;
    
    for( k = 0; k < Keff; k++ )
    {
        ktrue = RespEff[k];
        cb = CondBlock[ktrue];

        /* Initialize Loop */
        
        blk_cnt = cb[0]; 
        last = 0;
        sum = 0.0;
        doty = 0.0;

        for( kless = 0; kless <= k; kless++ )
            StS[ kless + k*Keff ] = 0.0;
        
        for( i = 0; i < blk_cnt; i++ )
        {
            b = cb[i+1];
            j = BlockStart[b];
            w = WhichStim[b];
            l = Min( bell_len[w], T - j );

            t = Min( last, j );
            
            for( ; t < j; t++ )
                S[t + k*T] = 0.0;

            /* Previous bell overlaps current block, add new contribution  */
            /* All of the terms except StS are linear in S, so the overlap */
            /* just adds on wlog.  In the case of StS, however, we need to */
            /* be more careful about the combination.  Here, we **ASSUME** */
            /* that the overlap is of only two pieces of bell.  Hence, we  */
            /* record the previous and new terms (call them x and y).  The */
            /* contributions to the other matrices use y (since x already  */
            /* is included), but diag(StS) which has x^2 needs y^2 + 2 x y */
            /* to be correct.  The case of two overlaps is uncommon enough,*/
            /* to have 3, the bell length would have to be quite long;     */
            /* will bet for now that this will not happen, but must look   */
            /* into it at a later time.                                    */
            /*                                                             */
            /* NOTE: I think the conditional at the end of this loop       */
            /*       fixes the problem by keeping last equal to the        */
            /*       last point ever set.  As long as that is > t, we      */
            /*       always are in overlap mode.  ATTN: Check this!        */
            /* N.B.  last is 1 ahead of what has last been set.            */

            for( m = 0; (m < last - j) && (m < l); m++, t++ )   
            {
                u = bells[m + w*2*T];
                v = S[t + k*T];
                S[t + k*T] += u;

                for( kless = 0; kless < k; kless++ )
                    StS[ kless + k*Keff ] += u * S[ t + kless*T ];

                StS[k + k*Keff] += u * u + 2.0 * u * v;
                sum += u;
                doty += u * y[t];
            }

            /* Continue over region of no overlap */

            for( ; m < l; m++, t++ )
            {
                u = S[t + k*T] = bells[m + w*2*T];

                for( kless = 0; kless < k; kless++ )
                    StS[ kless + k*Keff ] += u * S[ t + kless*T ];

                StS[k + k*Keff] += u * u;
                sum += u;
                doty += u * y[t];
            }

            if( last < t )
                last = t;
        }

        for( t = last; t < T; t++ )  /* Fill remainder of profile with zeros   */
            S[t + k*T] = 0.0;

        for( kless = 0; kless < k; kless++ )                 /* StS symmetric  */
            StS[ k + kless*Keff ] = StS[ kless + k*Keff ];

        St1[k] = sum;                                        /* Record (S'1)_k */
        Sty[k] = doty;                                       /* Record (S'y)_k */
    }
}


/*
 * SMOOTH_PASTE is a C^2 approximation to the function Max(x,y).
 * The quality of the approximation (and thus the magnitude of the derivatives)
 * is controlled by the parameter s.  The goal of this function is to come
 * up with a reasonable approximation to Max() without getting too close; setting
 * s too small will result in such large derivatives that the function might
 * as well be non-smooth.
 *
 */

double  smooth_paste( double x, double y, double s )
{
    double   u, v, w, z;

    u = (double)x/(double)s;
    v = (double)y/(double)s;
    w = u - v;
    z = 0.5*(w + 1.0);
    
    if( w >= 1.0 )
        return( (double)x );
    else if( w <= -1 )
        return( (double)y );
    else
    {
        w = 10.0 * Cube(z) - 15.0 * Fourth(z) + 6.0 * Fifth(z);
        
        return( s * (w * u + (1.0 - w) * v) );
    }
}

/* dSMOOTH_PASTE computes the gradient of the smooth_paste function */

void   Dsmooth_paste( double x, double y, double s, double *grad )
{
    double   a, d, r, u, v, w, z;

    r = 1.0/(double)s; 
    u = x * r;
    v = y * r;
    w = u - v;
    z = 0.5*(w + 1.0);
    
    if( w >= 1.0 )
    {
        grad[0] = r;
        grad[1] = 0.0;
    }
    else if( w <= -1 )
    {
        grad[0] = 0.0;
        grad[1] = r;
    }        
    else
    {
        a = 10.0 * Cube(z) - 15.0 * Fourth(z) + 6.0 * Fifth(z);
        d = 30.0 * Square(z) - 60.0 * Cube(z) + 30 * Fourth(z);   /* da/dd */
        
        grad[0] =  d * 0.5 * w + a * r;
        grad[1] =  d * 0.5 * w + (1.0 - a) * r;
    }
}

/*** Set-Up Routines ***/

/*
 * PARSE_MODEL determines which subsets of the conditions
 * are to be used in a given fit.  Eventually a full specification
 * will be allowed, but in the current version, the possibilities are
 *
 *    name          #sub-models      which conditions included
 *    ----------    -----------      -------------------------
 *    full              1             all non-fixed conditions
 *    full+null         2             all non-fixed conditions; none
 *    null+full         2             all non-fixed conditions; none
 *    all-full       2^Keff - 1       all subsets of non-fixed except full model
 *    null              1             none
 *    all             2^Keff          all subsets of non-fixed     
 *
 */

void    parse_model( char *modstr, long *modind )
{
    if( !strcasecmp( modstr, "full" ) )
        *modind = FULL_MODEL;
    else if( !strcasecmp( modstr, "full+null" ) || !strcasecmp( modstr, "null+full" ) )
        *modind = FULL_MODEL | NULL_MODEL;
    else if( !strcasecmp( modstr, "null" ) )
        *modind = NULL_MODEL;
    else if( !strcasecmp( modstr, "maximal" ) )
        *modind = MAXIMAL_MODELS | NULL_MODEL;
    else if( !strcasecmp( modstr, "all-full" ) )
        *modind = ALL_MODELS ^ FULL_MODEL;
    else if( !strcasecmp( modstr, "all" ) )
        *modind = ALL_MODELS;
}



/*** Initialization Routines ***/

void    init_design( void )
{
    int       b, i, j, k, t;
    int       *condcount;
    double    u;
    
        /* Initialize Condition and Stimulus Tables */

    Stims       = Malloc( B, double );
    WhichStim   = Malloc( B, int );
    BlockStart  = Malloc( B, int );
    BlockEnd    = Malloc( B, int );
    StimShift   = Calloc( B, double );

    CondBlock   = Calloc( K, int * );
    CondImages  = Calloc( K, int );
    CondMap     = Calloc( T, int );
    CondMatrix  = Matrix( K, T, int );   /* Initializes entries to zero */


    if( Design_Images )    /* Put all times in units of seconds */
    {
        for( b = 0; b < B; b++ )
        {
            StimulusStart[b] *= IAI;
            StimulusLen[b]   *= IAI; 
        }
    }
    
        /* Set Block Addresses and Uniquify the Stimulus Lengths */

    for( b = 0, u = 0.0; b < B; b++ )          
    {
        BlockStart[b] = StimulusStart[b]/IAI;                   /* Inaccurate in unaligned case */
        BlockEnd[b] = (StimulusStart[b] + StimulusLen[b])/IAI;  /* Not Used                     */
    }

    if( Design_Aligned )
    {
        for( b = 0, j = 0; b < B; b++ )          
        {
            for( i = 0; i < j; i++ )
                if( Stims[i] == StimulusLen[b] )
                    break;

            if( i < j )
            {
                WhichStim[b] = i;
                continue;
            }
            else
            {
                Stims[j] = StimulusLen[b];
                WhichStim[b] = j;
                StimShift[j] = 0.0;  /* No shift in aligned case */
                j++;
            }
        }
                
        UniqueStims = j;
    }
    else
    {
        /* NOTE:  Could make this smaller in *some* cases by mimicking the     */
        /* above code using both the length and the shift as cues in selecting */
        /* the unique stims.  Probably not worthwhile in general though.       */
        
        UniqueStims = B;

        for( b = 0; b < B; b++ )
        {
            WhichStim[b] = b;
            Stims[b] = StimulusLen[b];
            StimShift[b] = StimulusStart[b] - BlockStart[b]*IAI;
        }
    }
    
    
    /* Each entry of CondBlock is a vector  (L, block1, ...,blockL )              */
    /* where blocki is the ith block in associated with that condition.           */
    /* Go through a little trouble here to avoid waste of memory.                 */

    condcount = Calloc( K, int );

    for( b = 0; b < B; b++ )
        condcount[Cond[b]]++;

    for( k = 0; k < K; k++ )
    {
        CondBlock[k] = Calloc( 1 + condcount[k], int );
        CondBlock[k][0] = condcount[k];
        condcount[k] = 0;
    }

    for( b = 0; b < B; b++ )
    {
        k = Cond[b];
        
        CondBlock[k][++condcount[k]] = b;
    }

    Free( condcount );
    

    /* CondMap and CondMatrix are not well-defined when the ends of the condition */
    /* blocks can occur after the beginning of the subsequent block.              */
    
    for( t = 0, b = 0, u = 0.0; t < T; t++ )
    {
        k = Cond[b];

        CondMap[t] = k;
        CondMatrix[k][t] = 1;
        CondImages[k]++;

        u += IAI;

        if( b < B - 1 && u >= StimulusStart[b+1] )
            b++;
    }
}

void    init_fixed( int ktot, int *keff )
{
    int       i, k;

    *keff = ktot;

    /* Fix responsiveness for specified conditions */

    Constrained = Calloc( K, int );

    for( i = 0; i < Num_Fixed; i++ )
    {
        k = Fixed[i];
            
        if( k >= K || k < 0 )
        {
            Warning( DIAGNOSTIC >= 1, "Condition %d to be fixed must fall in range 0..%d, ignoring...",
                     K, k );
            continue;
        }
                
        *keff -= (!Constrained[k]) ? 1 : 0;
        Constrained[k] = 1;
    }
}

void     init_params( void )
{
    init_par_indices( );
    init_bounds( );
}

void     init_par_indices( void )
{
    int           d, j, k, keff, s;
    
    /* Set up Parameters */
 
    Param     = Calloc( NParams,       double );
    saveParam = Calloc( 2*NParams + 2, double );
    Out       = Malloc( NParams,       int );

    Drift = Calloc( D + Dk, int );
    Resp  = Calloc( K, int );
    Shape = Calloc( TOTAL_SHAPE_PARAMS, int );

    Coefs = Drift;        /* Aliases into drift structure */
    Knots = Drift + D;

    RespEff = Calloc( K, int );

    /* Set coefficient indices */

    Mu = 0;                               /* Baseline */

    for( d = 0, j = 1; d < D; d++ )       /* Drift Coefficients */
        Coefs[d] = j++;

    if( !Knots_Fixed )                    /* Drift Knots if active */
    {
        for( d = 0; d < Dk; d++ )
            Knots[d] = j++;      
    }


    o_SigmaSq = j++;                      /* Noise  */
    

                                          /* Shape */

    Shape_Use = Calloc( TOTAL_SHAPE_PARAMS, int );

    switch( Shape_Params )  /* Select Shape Parameterization */
    {
      case 0:     /* No bells at all ==> no resp also put in place */
        
        Shape[LAG_ON] = j++;         /* Keep these before Resp anyway, neither used */
        Shape[ATTACK] = j++;
        Shape[LAG_OFF] = j++;
        Shape[DECAY] = j++;

        Shape[RISE] = NParams - 4;
        Shape[FALL] = NParams - 3;
        Shape[DIP_HGT] = NParams - 2;
        Shape[DIP_SKEW] = NParams - 1;

        break;

      case 2:     /* LAG_ON(=LAG_OFF), ATTACK(=DECAY), RISE=FALL=0, DIP_HGT=0, DIP_SKEW=0 */

        Shape[LAG_ON] = j++;
        Shape[ATTACK] = j++;
        Shape[LAG_OFF] = Shape[LAG_ON];
        Shape[DECAY] = Shape[ATTACK];

        Shape[RISE] = NParams - 4;
        Shape[FALL] = NParams - 3;
        Shape[DIP_HGT] = NParams - 2;
        Shape[DIP_SKEW] = NParams - 1;

        Shape_Use[LAG_ON] = 1;
        Shape_Use[ATTACK] = 1;

        break;

      case 3:     /* LAG_ON(=LAG_OFF), ATTACK, DECAY, RISE=FALL=0, DIP_HGT=0, DIP_SKEW=0 */

        Shape[LAG_ON] = j++;
        Shape[ATTACK] = j++;
        Shape[DECAY] = j++;
        Shape[LAG_OFF] = Shape[LAG_ON];

        Shape[RISE] = NParams - 4;
        Shape[FALL] = NParams - 3;
        Shape[DIP_HGT] = NParams - 2;
        Shape[DIP_SKEW] = NParams - 1;

        Shape_Use[LAG_ON] = 1;
        Shape_Use[ATTACK] = 1;
        Shape_Use[DECAY] = 1;

        break;

      case 4:     /* LAG_ON, ATTACK, LAG_OFF, DECAY, RISE=FALL=0, DIP_HGT=0, DIP_SKEW=0 */

        Shape[LAG_ON] = j++;
        Shape[ATTACK] = j++;
        Shape[LAG_OFF] = j++;
        Shape[DECAY] = j++;

        Shape[RISE] = NParams - 4;
        Shape[FALL] = NParams - 3;
        Shape[DIP_HGT] = NParams - 2;
        Shape[DIP_SKEW] = NParams - 1;

        Shape_Use[LAG_ON] = 1;
        Shape_Use[ATTACK] = 1;
        Shape_Use[LAG_OFF] = 1;
        Shape_Use[DECAY] = 1;

        break;

      case 5:     /* LAG_ON, ATTACK, LAG_OFF, DECAY, RISE(=FALL), DIP_HGT=0, DIP_SKEW=0 */

        Shape[LAG_ON] = j++;
        Shape[ATTACK] = j++;
        Shape[LAG_OFF] = j++;
        Shape[DECAY] = j++;
        Shape[RISE] = j++;
        Shape[FALL] = Shape[RISE];
        
        Shape[DIP_HGT] = NParams - 2;
        Shape[DIP_SKEW] = NParams - 1;

        Shape_Use[LAG_ON] = 1;
        Shape_Use[ATTACK] = 1;
        Shape_Use[LAG_OFF] = 1;
        Shape_Use[DECAY] = 1;
        Shape_Use[RISE] = 1;

        break;

      case 6:     /* LAG_ON, ATTACK, LAG_OFF, DECAY, RISE, FALL, DIP_HGT=0, DIP_SKEW=0 */

        Shape[LAG_ON] = j++;
        Shape[ATTACK] = j++;
        Shape[LAG_OFF] = j++;
        Shape[DECAY] = j++;
        Shape[RISE] = j++;
        Shape[FALL] = j++;
        
        Shape[DIP_HGT] = NParams - 2;
        Shape[DIP_SKEW] = NParams - 1;

        Shape_Use[LAG_ON] = 1;
        Shape_Use[ATTACK] = 1;
        Shape_Use[LAG_OFF] = 1;
        Shape_Use[DECAY] = 1;
        Shape_Use[RISE] = 1;
        Shape_Use[FALL] = 1;

        break;

      case 7:     /* LAG_ON, ATTACK, LAG_OFF, DECAY, RISE, FALL, DIP_HGT, DIP_SKEW=0 */

        Shape[LAG_ON] = j++;
        Shape[ATTACK] = j++;
        Shape[LAG_OFF] = j++;
        Shape[DECAY] = j++;
        Shape[RISE] = j++;
        Shape[FALL] = j++;
        Shape[DIP_HGT] = j++;

        Shape[DIP_SKEW] = NParams - 1;

        Shape_Use[LAG_ON] = 1;
        Shape_Use[ATTACK] = 1;
        Shape_Use[LAG_OFF] = 1;
        Shape_Use[DECAY] = 1;
        Shape_Use[RISE] = 1;
        Shape_Use[FALL] = 1;
        Shape_Use[DIP_HGT] = 1;

        break;

      case 8:     /* LAG_ON, ATTACK, LAG_OFF, DECAY, RISE, FALL, DIP_HGT, DIP_SKEW */

        Shape[LAG_ON] = j++;
        Shape[ATTACK] = j++;
        Shape[LAG_OFF] = j++;
        Shape[DECAY] = j++;
        Shape[RISE] = j++;
        Shape[FALL] = j++;
        Shape[DIP_HGT] = j++;
        Shape[DIP_SKEW] = j++;

        Shape_Use[LAG_ON] = 1;
        Shape_Use[ATTACK] = 1;
        Shape_Use[LAG_OFF] = 1;
        Shape_Use[DECAY] = 1;
        Shape_Use[RISE] = 1;
        Shape_Use[FALL] = 1;
        Shape_Use[DIP_HGT] = 1;
        Shape_Use[DIP_SKEW] = 1;

        break;

      default:

        Abort( "Illegal number (%d) of shape parameters (tag shape.params)", Shape_Params );
    }

    RespInd = j;                     /* Mark index of first effective responsiveness */

    for( k = 0, keff = 0; k < K; k++ )         /* Unconstrained conditions */
        if( !Constrained[k] )
        {
            Resp[k] = j++;
            RespEff[keff++] = k;
        }

    for( k = 0; k < K; k++ )         /* Constrained Conditions: index > NeffParams */
        if( Constrained[k] )         /* These are not touched during optimization  */
            Resp[k] = j++;


    if( Knots_Fixed )
    {
        for( d = 0; d < Dk; d++ )
            Knots[d] = j++;      
    }


    /* Set Output Order */

    Out[Mu] = 0;

    for( d = 0, j = 1; d < D; d++ )
        Out[Coefs[d]] = j++;

    if( !Knots_Fixed )
        for( d = 0; d < Dk; d++ )
            Out[Knots[d]] = j++;      

    Out[o_SigmaSq] = j++;

    for( s = 0; s < Shape_Params; s++ )
        Out[Shape[s]] = j++;

    for( k = 0; k < K; k++ )
        Out[Resp[k]] = j++;
}

void    init_bounds( void )
{
    int        i;
    
    Param_LB = Calloc( NParams, double );
    Param_UB = Calloc( NParams, double );

    Param_LB[Mu] = 0.0;
    Param_UB[Mu] = BigNumber;

    for( i = 0; i < D; i++ )
    {
        Param_LB[Coefs[i]] = -BigNumber;
        Param_UB[Coefs[i]] =  BigNumber;
    }

    for( i = 0; i < Dk; i++ )
    {
        Param_LB[Knots[i]] = 0.0;
        Param_UB[Knots[i]] = 1.0;
    }

    for( i = 0; i < TOTAL_SHAPE_PARAMS; i++ )
    {
        Param_LB[Shape[i]] = 1.0e-9; /* 0.0; */
        Param_UB[Shape[i]] = BigNumber;
    }

    Param_LB[Shape[RISE]] = -1.0;
    Param_UB[Shape[RISE]] =  1.0;
    Param_LB[Shape[FALL]] = -1.0;
    Param_UB[Shape[FALL]] =  1.0;
        
    for( i = 0; i < K; i++ )
    {
        Param_LB[Resp[i]] = 0.0;
        Param_UB[Resp[i]] = BigNumber;
    }

    Param_LB[o_SigmaSq] = 1.0e-16; /* 0.0; */
    Param_UB[o_SigmaSq] = BigNumber;
}

void     init_basis( void )
{
    int      i, j, m, info = 0;

    /*
     * Drift Basis Matrix
     *
     * The drift basis is represented by a D vectors
     * of length T (stored as a D x T matrix in col-major order)
     * that are the samples of functions on the equally
     * spaced grid of image acquisition times.  The canonical
     * basis does not include the constant term since the
     * baseline is treated separately.  However, for the purposes
     * of orthogonalization, it helps to have Basis contain the
     * constant vector; it is shifted to the -1 position to keep
     * it in the background.  Thus, passing Basis[-1] passes the
     * full (constant included) matrix, while passing Basis
     * gives the reduced matrix.
     *
     */

    Basis = Matrix( (D + 1), T, double ) + 1;

    init_binom_coefs( );
    init_detrend( T, D, Dk, NULL, 0 );

    /*
     * Drift Prior Matrices
     *
     * The drift prior matrices are the quadratic forms defined
     * by S1 and S2 in terms of the current drift coefficients
     * that correspond to the two-norm and Sobolev(2,2) semi-norm
     * of the drift profile.
     *
     * Note that Si = (R^-1)' Pi (R^-1) where Pi is the corresponding
     * matrix for the power basis and R is taken from the QR-decomposition
     * that defines the basis matrix above.  Hence,
     *   a S1 + b S2 = (R^-1)' (a P1 + b P2) (R^-1)
     * so
     *   Det( a S1 + b S2 ) = Det(R^-1)^2 Det( a P1 + b P2 )
     * and
     *   d' (a S1 + b S2) d = dr'(a P1 + b P2) dr
     * where dr = (R^-1) d.  Since R is triangular, this makes computation
     * of the quadratic form and determinants easier.
     *
     * The components of the Pi's can be computed analytically as follows:
     *
     * The matrices are given the name Pnorm for the L2 matrix and Pcurv
     * for the L2,2 matrix.
     *
     */
    
    Pnorm = Calloc( (D+1) * (D+1), double );
    Pcurv = Calloc( (D+1) * (D+1), double );

    Pcomb = Calloc( (D+1) * (D+1), double );     /* Combined Power Basis     */
    Ptot  = Calloc( (D+1) * (D+1), double );     /* Combined, orth'l basis   */

    basisR = Calloc( T * (D + 1), double );      /* Storage for QR structure */
    basisR_inv = Calloc( (D+1) * (D+1), double );/* R^-1 from QR structure   */

    if( Knots_Fixed )         /* Prepare Drift Computations once and for all */
    {
        qsort( Fixed_Knots, Dk, sizeof(double), Numerically_d );

        make_drift_basis( T, Dg, Dk, &Basis[-1], Fixed_Knots, basisR, Work, WorkSize );

        if( Drift_c >= D )
        {
            Warning( DIAGNOSTIC >= 1, "Requested d.f. greater than drift order, no smoothing used" );
            Drift_a = Drift_b = Drift_c = 0.0;
        }

        if( Drift_c == 0.0 )
            Drift_df = D;
        else
            Drift_df = Drift_c;

        if( Drift_a > 0.0 || Drift_b > 0.0 )
        {
            make_drift_qforms( T, Dg, Dk, Fixed_Knots, Pnorm, Pcurv, NULL, Fixed_Knots, Work, WorkSize );

            /* Set Pcomb = a Pnorm + b Pcurv */

            m = 1 + Dg + Dk;

            for( i = 0; i < m; i++ )       
            {
                for( j = 0; j < i; j++ )
                    Pcomb[i + j*m] = Pcomb[j + i*m] = (Drift_a*Pnorm[i + j*m] + Drift_b*Pcurv[i + j*m]);

                Pcomb[i + i*m] = (Drift_a * Pnorm[i + i*m] + Drift_b * Pcurv[i + i*m]);
            }

            /* Set Ptot = t(R^-1) . Pcomb . R^-1 */

            for( i = 0; i < m; i++ )       
            {
                for( j = 0; j < i; j++ )
                {
                    basisR_inv[j + i*m] = basisR[j + i*T];
                    Ptot[i + j*m] = Ptot[j + i*m] = Pcomb[i + j*m];
                }

                basisR_inv[i + i*m] = basisR[i + i*T];
                Ptot[i + i*m] = Pcomb[i + i*m];
            }

            dtrtri( &Upper, &NoUnit, &m, basisR_inv, &m, &info );

            if( info )
                Abort( "Could not invert R during initialization (info = %d).", info );

            dtrmm( &Right, &Upper, &NoTrans, &NoUnit, &m, &m, &One, basisR_inv, &m, Ptot, &m );
            dtrmm( &Left,  &Upper, &Trans,   &NoUnit, &m, &m, &One, basisR_inv, &m, Ptot, &m );
        }

        if( Drift_c > 0.0 )
        {
            for( i = 0; i < D; i++ )       
            {
                for( j = 0; j < i; j++ )
                    Work[j + i*D] = Ptot[(j+1) + (i+1)*(D+1)];

                Work[i + i*D] = Ptot[(i+1) + (i+1)*(D+1)];
            }

            find_drift_df( D, Work, Drift_c, &Smooth_Target, NULL, NULL, Work + D*D,
                           Work + D*(D+1), WorkSize - D*(D+1), 0.01 );
        }
        else
            Smooth_Target = 1.0;
    }
}

void     init_obj_scale( int method, double *psc, double *pfsc, int *iwork, double *work, int lwork )
{
    int     i, d, k, m, s, info = 0;

    if( method == 1 )   /* Fast Optimization:  Currently lbfgs qnb routine */
    {
        /* Default Settings: No information */

        *pfsc = 1.0;

        for( i = 0; i < NParams; i++ )
            psc[i] = 1.0;

        /*
         * Performance of the IMSL qnb algorithm improves
         * when the small parameters  o_SigmaSq, Resp, and Shape
         * are scaled.  Without the scalings below, the trial
         * steps searching to build a trust region initially
         * are orders of magnitude too large.  Not only does this
         * waste cycles but it effectively obliterates any
         * initial information passed to the algorithm.
         *
         * Here, we scale them roughly by the number of decimal places
         * we expect in each case.  This can be improved but seems to
         * do a good job overall.
         *
         * ATTN: DO A BETTER SCALING JOB HERE 
         *
         */
        
        for( k = 0; k < Keff; k++ )
            psc[RespInd + k] = 1000.0;

        psc[o_SigmaSq] = 1000.0;

        for( k = 0; k < Shape_Params; k++ )
            psc[Shape[k]] = 10.0;
    }
    else  /* Amoeba */
    {
        /* Initialize the scaling properties of the objective function      */
        /* Basically, use standard deviations of priors for each parameter  */
        /* Some ad hoc changes are made to keep the init reasonable         */

        psc[Mu] = 4.0 * sqrt(Mu_b);  /* In spirit only, the sd of t_1 does not exist */

        for( k = 0; k < Dk; k++ )
            psc[Knots[k]] = 0.2/(Dk + 1.0);

        for( s = 0; s < TOTAL_SHAPE_PARAMS; s++ )
            psc[Shape[s]] = sqrt(Shape_a[s])/Shape_b[s];

        for( k = 0; k < K; k++ )
            psc[Resp[k]] = (Constrained[k]) ? 0.0 : Max(0.5*Resp_a,0.25/Resp_b);

        psc[o_SigmaSq] = 0.25 * exp( -(SigmaSq_a + 0.5 * SigmaSq_b * SigmaSq_b) )/
            sqrt( exp(SigmaSq_b * SigmaSq_b) - 1.0 );

        /*
          if( Drift_a > 0.0 && Drift_b > 0.0 )
          {
          for( d = 0; d < (D+1)*(D+1); d++ )
          work[d] = Pcomb[d];

          lwork -= (D+1)*(D+1);
          m = D + 1;
        
          dsytrf( &Upper, &m, work, &m, iwork, work + (D+1)*(D+1), &lwork, &info );

          if( !info )
          dsytri( &Upper, &m, work, &m, iwork, work + (D+1)*(D+1), &info );

          if( info )
          {
          for( d = 0; d < D; d++ )
          psc[Coefs[d]] = sqrt( (double)T );

          Warning( DIAGNOSTIC >= 1,
                  "Could not compute inverse drift prior covariance matrix during init (%d)", info );
          }
          else
          for( d = 0; d < D; d++ )
          psc[Coefs[d]] = sqrt(work[d + d*D]) * ((Drift_scaled) ? Mu_a : 1.0);
          }
          else
          for( d = 0; d < D; d++ )
          psc[Coefs[d]] = sqrt( (double)T );
          */

        /* Temporary because above is wrong, using Pcomb instead of R^-1... */

        for( d = 0; d < D; d++ )
            psc[Coefs[d]] = sqrt( (double)T );

        *pfsc = 0.01 * (T/2.0);   /* One percent of the mean, unormalized, true likelihood component */
    }
}

void     init_standard_errors( double *psc, int *iwork, double *work, int lwork )
{
    int     d, k, s, info = 0;

    Hinit   = Calloc( NParams, double );
    ObsInfo = Calloc( (NParams * (NParams + 1))/2, double );

    SEworksize  = 10 * NParams + (NParams * (NParams + 1))/2;
    SEwork  = Calloc( SEworksize, double );
    
    /* Initialize the scaling for richardson extrapolation */

    psc[Mu] = 0.1 * sqrt(Mu_b);  /* In spirit only, the sd of t_1 does not exist */

    for( k = 0; k < Dk; k++ )
        psc[Knots[k]] = 0.2/(Dk + 1.0);

    for( s = 0; s < TOTAL_SHAPE_PARAMS; s++ )
        psc[Shape[s]] = 0.1 * sqrt(Shape_a[s])/Shape_b[s];

    for( k = 0; k < K; k++ )
        psc[Resp[k]] = 0.5 * ((Constrained[k]) ? 0.0 : Max(0.5*Resp_a,0.25/Resp_b));

    psc[o_SigmaSq] = exp( -(SigmaSq_a + 0.5 * SigmaSq_b * SigmaSq_b) )/
                      sqrt( exp(SigmaSq_b * SigmaSq_b) - 1.0 );


    /*
    if( Drift_a > 0.0 && Drift_b > 0.0 )
    {
        for( d = 0; d < D*D; d++ )
            work[d] = Pcomb[d];

        lwork -= D*D;
        
        dsytrf( &Upper, &D, work, &D, iwork, work + D*D, &lwork, &info );

        if( !info )
            dsytri( &Upper, &D, work, &D, iwork, work + D*D, &info );

        if( info )
        {
            for( d = 0; d < D; d++ )
                psc[Coefs[d]] = sqrt( (double)T );

            Warning( DIAGNOSTIC >= 1,
                     "Could not compute inverse drift prior covariance matrix during init (%d)", info );
        }
        else
            for( d = 0; d < D; d++ )
                psc[Coefs[d]] = 0.1 * sqrt(work[d + d*D]) * ((Drift_scaled) ? Mu_a : 1.0);
    }
    else
        for( d = 0; d < D; d++ )
            psc[Coefs[d]] = 0.1 * sqrt( (double)T );
    */

    for( d = 0; d < D; d++ )
        psc[Coefs[d]] = 0.1 * sqrt( (double)T );
}

FILE    *init_data( char *fn, char *fmt, long *phdr, long *ptype, int *pV, int *pT, int **pIdim, DataAccessor *dataAccessor )
{
    int      i, v;
    
    long     ndim, vlen;
    long     dims[FIELD_DIM];

    char    *name, *scratch;
    char    *path, *fullname;
    char    *typename, *missing;

    FILE    *fp;
    
    /*
     * Prepare data file for reading in time-series by time-series.
     * Several possible formats are supported, but in general we assume
     * the data field has been permuted into txyz order, from fastest to
     * slowest varying.
     *
     * The basic format is the voxel "list" format, which consists of
     * an optional input file giving the essential information, and a
     * binary file consisting of consecutive time series for each
     * voxel in the list.
     *
     * Other supported formats include the Pittsburgh MRI format
     * (both new and old) and a raw binary format in which the required
     * information must be included in the current input file.
     *
     */

    if( !strcasecmp( fmt, "list" ) )
    {
        if( !inpf_scan( fn ) )
        {
            Message( "Errors encountered in input file %s\n", fn );
            inpf_perror( );
            exit( 0 );
        }

        inpf_get( "file",    "%ls",        &name );
        inpf_get( "voxels",  "%d",         pV );
        inpf_get( "dim.t",   "%d",         pT );
        inpf_get( "type",    "%ls[float]", &typename );
        inpf_get( "header",  "%d[0]",      phdr );
        inpf_get( "missing", "%ls[none]",  &missing );     /* Deal with this later */
        inpf_get( "image",   "%ls[no]",    &scratch );

        for( i = 0, *ptype = -1; *Type_String[i]; i++ )
        {
            if( !strcasecmp( typename, Type_String[i] ) )
            {
                *ptype = i;  /* Assumes that MRI types go from 0..something */
                break;
            }
        }

        if( *ptype == -1 )
            Abort( "Improper type (%s) given in data file (tag type)", typename );
        
        if( !strcasecmp( scratch, "yes" ) )
        {
            (*pIdim) = Calloc( FIELD_DIM, int );
            (*pIdim)[0] = *pT;
            
            inpf_get( "dim.x",  "%d\n", &((*pIdim)[1]) );
            inpf_get( "dim.y",  "%d\n", &((*pIdim)[2]) );
            inpf_get( "dim.z",  "%d\n", &((*pIdim)[3]) );

            if( *pV != (*pIdim)[1] * (*pIdim)[2] * (*pIdim)[3] )
                Warning( DIAGNOSTIC >= 1, "Apparent inconsistency between # of voxels (%d) and "
                            "product of spatial dimensions, input file %s.", fn );
        }
        else
            (*pIdim) = NULL;

        if( inpf_unmatched( "file" ) )
            Abort( "Data file (tag file) missing from input file %s", fn );

        if( inpf_unmatched( "voxels" ) )
            Abort( "Number of voxels (tag voxels) missing from input file %s.", fn );

        if( inpf_unmatched( "dim.t" ) )
            Abort( "Number of time points (tag dim.t) missing from input file %s.", fn );

        /* Open data file -- but relatively to the path of the given input file      */
        /* If input file contains a path component and data file is relative, use it */
        /* The part of the name after last slash is file name, before is the path    */

        if( (*name != '/') && (path = strrchr(fn,'/')) ) 
        {
            fullname = Calloc( path - fn + strlen(name) + 2, char );
            strncpy( fullname, fn, path - fn + 1 );
            fullname[path - fn + 1] = '\0';
            strcat( fullname, name );

            fp = efopen( fullname, "rb" );
            Free( fullname );
        }
        else
            fp = efopen( name, "rb" );

        if( *phdr )
            fseek( fp, *phdr, SEEK_SET );

        Free( name );
        Free( scratch );
        Free( missing );

        inpf_cleanup( );

	dataAccessor->setup = defaultDataSetup;
	dataAccessor->read  = defaultDataRead;
	dataAccessor->close = defaultDataClose;
	dataAccessor->state = NULL;
    }
    else if( !strcasecmp( fmt, "pgh" ) )
    {
        /* This code contributed by Joel Welling */

        MRI_Dataset*       Input;
        char*              images_dimensions = NULL;
        char               tag_name_buf[64];
        int                count;
        int                ndims;
        int                i;

        /* Open input dataset */

        Input = mri_open_dataset( fn, MRI_READ );

        if (!Input)
            Abort( "Error opening input dataset %f.\n", fn );

        /*
         * Check that program will function on data-set.  We require
         * that the first dimension be t. It is also permitted for
         * the first dimension to be v, if the dimensionality of
         * v is 1 (trivial case).
         *
         */

        if( !mri_has( Input, "images" ) )
            Abort( "Input file lacks standard images.\n" );

        if (!mri_has( Input, "images.dimensions" ))
            Abort("Required images.dimensions element absent from input data.\n");
    
        images_dimensions = mri_get_string( Input, "images.dimensions" );

        if( !((*images_dimensions == 't') || ((*images_dimensions == 'v') && (*(images_dimensions+1) == 't')
                                              && mri_has(Input, "images.extent.v")
                                              && (mri_get_int(Input, "images.extent.v") == 1))) )
            Abort( "Input data error: t must vary fastest.\n");

        /* Ignore vectors of length 1 in remainder of calc */
    
        if( *images_dimensions == 'v')
            images_dimensions++;

        ndims = strlen(images_dimensions);
        (*pIdim) = Malloc( ndims, int );

        if( mri_has( Input, "images.extent.t" ) )
            (*pIdim)[0] = *pT= mri_get_int(Input, "images.extent.t");
        else
            Abort( "Input data error: images.extent.t not present.\n");

        count = 1;

        for( i = 1; i < ndims; i++ )
        {
            sprintf(tag_name_buf,"images.extent.%c",images_dimensions[i]);

            if( mri_has(Input,tag_name_buf) )
            {
                (*pIdim)[i]= mri_get_int(Input, tag_name_buf);
                count *= (*pIdim)[i];
            }
            else
                Abort( "input data error: %s not present.\n", tag_name_buf );
        }

        *pV = count;
        *phdr = 0;              /* since libmri will handle offsets   */
        *ptype = DOUBLE_T;      /* since libmri will convert on input */

        dataAccessor->setup = mriDataSetup;
        dataAccessor->read = mriDataRead;
        dataAccessor->close = mriDataClose;
        dataAccessor->state = NULL;
        if( !(dataAccessor->state = (MriDataState*)malloc(sizeof(MriDataState))) )
            Abort("Error attempting to allocate %d bytes!\n", sizeof(MriDataState));

        ((MriDataState*)(dataAccessor->state))->ds = Input;
        ((MriDataState*)(dataAccessor->state))->offset = 0;
    }
    else if( !strcasecmp( fmt, "pgh-old" ) )
    {
#if defined(USE_MRIFUNC)

        fp = mriReadHeader( fn, ptype, &ndim, &vlen, dims );

        (*pIdim) = Malloc( FIELD_DIM, int );

        for( i = 0, v = 1; i < ndim; i++ )
        {
            (*pIdim)[i] = dims[i];
            v *= dims[i];
        }

        for( i = ndim; i < FIELD_DIM; i++ )
            (*pIdim)[i] = 1;

        *pT = dims[0];
        *pV = v/dims[0];
        *phdr = sizeof(long) * (3 + ndim);

	dataAccessor->setup = defaultDataSetup;
	dataAccessor->read  = defaultDataRead;
	dataAccessor->close = defaultDataClose;
	dataAccessor->state = NULL;

#else

        Abort( "Format old-pgh not supported in this compilation, use pgh" );

#endif
    }
    else if( !strcasecmp( fmt, "raw" ) )
    {
        fp = efopen( fn, "rb" );

        if( *phdr )
            fseek( fp, *phdr, SEEK_SET );

	dataAccessor->setup = defaultDataSetup;
	dataAccessor->read  = defaultDataRead;
	dataAccessor->close = defaultDataClose;
	dataAccessor->state = NULL;
    }
    else
        Abort( "Unrecognized data format %s. Check value for data.format tag in input file.", fmt );

    return( fp );
}

FILE    *init_init( int *pnum, char **initv )
{
    int       i, num, off;
    FILE     *fp;

    Init_len = 0;
    Init_size = 0;

    num = *pnum;

    Init         = Calloc( NUMINITS, int );    /* Indicates which fields are user initialized     */
    Init_offsets = Calloc( NUMINITS, int );    /* Indicates offsets of each field in user records */

    if( !num || (num == 1 && !strcasecmp( initv[0], "auto" )) )
    {
        *pnum = 0;
        return( NULL );
    }

    fp = efopen( initv[0], "rb" );

    for( i = 0, off = 0; i < num; i++ )
    {
        if( !strcasecmp( initv[i], "baseline" ) )
        {
            Init[BASELINE] = 1;
            Init_offsets[BASELINE] = off;

            off++;
        }
        else if( !strcasecmp( initv[i], "coefficients" ) )
        {
            Init[DRIFT_COEFS] = 1;
            Init_offsets[DRIFT_COEFS] = off;

            off += D;
        }
        else if( !strcasecmp( initv[i], "knots" ) )
        {
            Init[DRIFT_KNOTS] = 1;
            Init_offsets[DRIFT_KNOTS] = off;

            off += Dk;
        }
        else if( !strcasecmp( initv[i], "shape" ) )
        {
            Init[SHAPE] = 1;
            Init_offsets[SHAPE] = off;

            off += TOTAL_SHAPE_PARAMS;
        }
        else if( !strcasecmp( initv[i], "responsiveness" ) )
        {
            Init[RESP] = 1;
            Init_offsets[RESP] = off;

            off += K;
        }
        else if( !strcasecmp( initv[i], "noise.precision" ) )
        {
            Init[O_SIGMASQ] = 1;
            Init_offsets[O_SIGMASQ] = off;

            off++;
        }
    }

    Init_len = off;

    if( !strcasecmp( Init_Type_String, "float" ) )
    {
        Init_type = FLOAT_T;
        Init_size = mriTypeSize( Init_type );
    }
    else if( !strcasecmp( Init_Type_String, "short" ) )
    {
        Init_type = SHORT_T;
        Init_size = mriTypeSize( Init_type );
    }
    else if( !strcasecmp( Init_Type_String, "int" ) )
    {
        Init_type = LONG_T;
        Init_size = mriTypeSize( Init_type );
    }
    else if( !strcasecmp( Init_Type_String, "long" ) )
    {
        Init_type = LONG_T;
        Init_size = mriTypeSize( Init_type );
    }
    else if( !strcasecmp( Init_Type_String, "double" ) )
    {
        Init_type = DOUBLE_T;
        Init_size = mriTypeSize( Init_type );
    }
    else if( !strcasecmp( Init_Type_String, "char" ) )
    {
        Init_type = CHAR_T;
        Init_size = mriTypeSize( Init_type );
    }
    else
        Abort( "Illegal initialization type %s; check tag vmpfx.init.type in input file.",
               Init_Type_String );

    *pnum = num - 1;
    return( fp );
}

void   init_init2( void )
{
        Xty  = Calloc( D,           double );
        XtS0 = Calloc( K * D,       double );
        S0ty = Calloc( K,           double );
        S    = Calloc( Keff * T,    double );
        StS  = Calloc( Keff * Keff, double );
        St1  = Calloc( Keff,        double );
        Sty  = Calloc( Keff,        double );
}

void   init_products( int *pnum, char ***prodv, int **produce, int **prtype,
                      long **units, long **unitsize, long **type, long *rsize, long dtype )
{
    int       i, j = -1, num, np;
    long      addr = 0, type_r = FLOAT_T;
    size_t    size_r;

    num = *pnum;                   /* This will double if Null_Model is true      */
    np = NoutParams;               /* Do not keep copying knots if they are fixed */

    if( Null_Model )
    {
        *pnum *= 2;
        (*prodv) = Realloc( (*prodv), (*pnum), char * );
    }

    (*prtype)   = Calloc( (1 + Null_Model) * num, int );        /* Type of the ith product */

    (*produce)  = Calloc( NUM_PRODUCTS, int );  
    (*units)    = Calloc( NUM_PRODUCTS, long );
    (*unitsize) = Calloc( NUM_PRODUCTS, long );
    (*type)     = Calloc( NUM_PRODUCTS, long );
    

    if( dtype == DOUBLE_T )    /* Residuals and Fitted values are float unless specified as double */
        type_r = dtype;

    size_r = mriTypeSize( type_r );

    /* Set record map for each product type */

    for( i = 0; i < num; i++ )
    {
        if( !strcmp( (*prodv)[i], "estimates" ) )
        {
            (*prtype)[i] = ESTIMATES;
            (*produce)[ESTIMATES] = 1;
            (*units)[ESTIMATES] = np;
            (*unitsize)[ESTIMATES] = sizeof(double);
            (*type)[ESTIMATES] = DOUBLE_T;

            if( Null_Model )
            {
                (*prodv)[i+num] = strdup( "null.estimates" );
                (*prtype)[i+num] = NULL_ESTIMATES;
                (*produce)[NULL_ESTIMATES] = 1;
                (*units)[NULL_ESTIMATES] = NnulParams;
                (*unitsize)[NULL_ESTIMATES] = sizeof(double);
                (*type)[NULL_ESTIMATES] = DOUBLE_T;
            }
        }
        else if( !strcmp( (*prodv)[i], "standard.errors" ) )
        {
            (*prtype)[i] = STANDARD_ERRORS;
            (*produce)[STANDARD_ERRORS] = 1;
            (*units)[STANDARD_ERRORS] = np;
            (*unitsize)[STANDARD_ERRORS] = sizeof(double);
            (*type)[STANDARD_ERRORS] = DOUBLE_T;

            if( Null_Model )
            {
                (*prodv)[i+num] = strdup( "null.standard.errors" );
                (*prtype)[i+num] = NULL_STANDARD_ERRORS;
                (*produce)[NULL_STANDARD_ERRORS] = 1;
                (*units)[NULL_STANDARD_ERRORS] = NnulParams;
                (*unitsize)[NULL_STANDARD_ERRORS] = sizeof(double);
                (*type)[NULL_STANDARD_ERRORS] = DOUBLE_T;
            }
        }
        else if( !strcmp( (*prodv)[i], "covariances" ) )
        {
            (*prtype)[i] = COVARIANCES;
            (*produce)[COVARIANCES] = 1;
            (*units)[COVARIANCES] = ((np * (np + 1))/2);
            (*unitsize)[COVARIANCES] = sizeof(double);
            (*type)[COVARIANCES] = DOUBLE_T;

            if( Null_Model )
            {
                (*prodv)[i+num] = strdup( "null.covariances" );
                (*prtype)[i+num] = NULL_COVARIANCES;
                (*produce)[NULL_COVARIANCES] = 1;
                (*units)[NULL_COVARIANCES] = ((NnulParams * NnulParams)/2);
                (*unitsize)[NULL_COVARIANCES] = sizeof(double);
                (*type)[NULL_COVARIANCES] = DOUBLE_T;
            }
        }
        else if( !strcmp( (*prodv)[i], "residuals" ) )
        {
            (*prtype)[i] = RESIDUALS;
            (*produce)[RESIDUALS] = 1;
            (*units)[RESIDUALS] = T;
            (*unitsize)[RESIDUALS] = size_r;
            (*type)[RESIDUALS] = type_r;

            if( Null_Model )
            {
                (*prodv)[i+num] = strdup( "null.residuals" );
                (*prtype)[i+num] = NULL_RESIDUALS;
                (*produce)[NULL_RESIDUALS] = 1;
                (*units)[NULL_RESIDUALS] = T;
                (*unitsize)[NULL_RESIDUALS] = size_r;
                (*type)[NULL_RESIDUALS] = type_r;
            }
        }
        else if( !strcmp( (*prodv)[i], "fitted.values" ) )
        {
            (*prtype)[i] = FITTED_VALUES;
            (*produce)[FITTED_VALUES] = 1;
            (*units)[FITTED_VALUES] = T;
            (*unitsize)[FITTED_VALUES] = size_r;
            (*type)[FITTED_VALUES] = type_r;

            if( Null_Model )
            {
                (*prodv)[i+num] = strdup( "null.fitted.values" );
                (*prtype)[i+num] = NULL_FITTED_VALUES;
                (*produce)[NULL_FITTED_VALUES] = 1;
                (*units)[NULL_FITTED_VALUES] = T;
                (*unitsize)[NULL_FITTED_VALUES] = size_r;
                (*type)[NULL_FITTED_VALUES] = type_r;
            }
        }
        else if( !strcmp( (*prodv)[i], "diagnostics" ) )
        {
            j = Null_Model | (*produce)[STANDARD_ERRORS] | (*produce)[COVARIANCES];

            (*prtype)[i] = DIAGNOSTICS;
            (*produce)[DIAGNOSTICS] = 1;
            (*units)[DIAGNOSTICS] = 1;
            (*unitsize)[DIAGNOSTICS] = (1 + j*2) * sizeof(float) + 2*sizeof(int);
            (*type)[DIAGNOSTICS] = RECORD_T;

            if( Null_Model )
            {
                (*prodv)[i+num] = strdup( "null.diagnostics" );
                (*prtype)[i+num] = NULL_DIAGNOSTICS;
                (*produce)[NULL_DIAGNOSTICS] = 1;
                (*units)[NULL_DIAGNOSTICS] = 1;
                (*unitsize)[NULL_DIAGNOSTICS] = 4 * sizeof(float) + 2*sizeof(int);
                (*type)[NULL_DIAGNOSTICS] = RECORD_T;
            }
        }
        else
            Abort( "Unrecognized product type in vmpfx: %s", (*prodv)[i] );

        addr += ((*units)[(*prtype)[i]]) * ((*unitsize)[(*prtype)[i]]);

        if( Null_Model )
            addr += ((*units)[(*prtype)[i+num]]) * ((*unitsize)[(*prtype)[i+num]]);
    }

    /* Adjust diagnostic and record sizes for "later" setting of standard_errors or covariances  */
    /* This is stored but if j == 0 above because standard.errors or covariances set afterwards, */
    /* then this was not set above.  Note that j is only set in the DIAGNOSTICS section          */
    /* ATTN: Rename j to something more specialized...                                           */

    if( !Null_Model && !j && ((*produce)[STANDARD_ERRORS] || (*produce)[COVARIANCES]) )
    {
        (*unitsize)[DIAGNOSTICS] += 2 * sizeof(float);
        addr += 2 * sizeof(float);
    }

    *rsize = addr;
}

/* Functions to support DataAccessors, transitional */

void defaultDataSetup( DataAccessor* this, FILE* fp, long n_vox_skip, int times_per_voxel, int type_size )
{
    fseek( fp, n_vox_skip * times_per_voxel * type_size, SEEK_SET );
}

void defaultDataRead( DataAccessor* this, void* buf, FILE* fp, int times_per_voxel, int type_size )
{
    efread( buf, type_size, times_per_voxel, fp );
}

void defaultDataClose( DataAccessor* this, FILE* fp )
{
    efclose( fp );
}

void mriDataSetup( DataAccessor* this, FILE* fp, long n_vox_skip, int times_per_voxel, int type_size )
{
    ((MriDataState*)(this->state))->offset += n_vox_skip*times_per_voxel;
}

void mriDataRead( DataAccessor* this, void* buf, FILE* fp, int times_per_voxel, int type_size )
{
    double* mribuf;
    double* runner;
    double* outbuf = (double*)buf;

    mribuf = mri_get_chunk( ((MriDataState*)(this->state))->ds, "images", times_per_voxel,
                            ((MriDataState*)(this->state))->offset, MRI_DOUBLE );

    ((MriDataState*)(this->state))->offset += times_per_voxel;

    for( runner = mribuf; runner < mribuf + times_per_voxel; runner++ )
        *outbuf++ = *runner;
}

void mriDataClose(DataAccessor* this, FILE* fp)
{
    mri_close_dataset( ((MriDataState*)(this->state))->ds );
}

/* Output Handling */

static CookieJar      cookies;


void      output_cookie( int code, ... )
{
    int           i, j;

    void         *x;
    size_t        size;
    va_list       ap;


    va_start( ap, code );

    switch( code )
    {
      case SET_IT:

        cookies.buf    = va_arg( ap, char * );
        cookies.bufsiz = va_arg( ap, long );
        cookies.stream = va_arg( ap, FILE * );
        cookies.offset = 0;

        break;

      case WRITE_IT:

        eat_cookies( &cookies );

        break;

      case ESTIMATES:

        i = va_arg( ap, int );
        x = va_arg( ap, double * );

        size = i * sizeof(double);

        bake_cookie( &cookies, (void *)x, size );

        break;

      case STANDARD_ERRORS:

        i = va_arg( ap, int );
        x = va_arg( ap, double * );

        size = i * sizeof(double);
        
        bake_cookie( &cookies, (void *)x, size );

        break;

      case RESIDUALS:

        i    = va_arg( ap, int );
        size = va_arg( ap, size_t );
        x    = va_arg( ap, void * );

        size *= i;
        
        bake_cookie( &cookies, (void *)x, size );

        break;

      case FITTED_VALUES:

        i    = va_arg( ap, int );
        size = va_arg( ap, size_t );
        x    = va_arg( ap, void * );

        size *= i;
        
        bake_cookie( &cookies, (void *)x, size );

        break;

      case COVARIANCES:

        i    = va_arg( ap, int );
        x    = va_arg( ap, double * );

        size = (i * (i + 1)/2) * sizeof(double);
        
        bake_cookie( &cookies, (void *)x, size );

        break;

      case DIAGNOSTICS:

        i    = va_arg( ap, int );
        j    = va_arg( ap, int );

            /* Real Diagnostics    */

        x    = va_arg( ap, float * );
        size = i * sizeof(float);
        
        bake_cookie( &cookies, (void *)x, size );

            /* Integer Diagnostics */

        x    = va_arg( ap, int * );
        size = j * sizeof(int);
        
        bake_cookie( &cookies, (void *)x, size );

        break;

      default:    /* Shouldn't happen */

        Warning( DIAGNOSTIC >= 1,
                 "Unrecognized product code (%d) at output, result addresses may be unreliable.", code );
        break;
    }

    va_end( ap );
}

void      bake_cookie( CookieJar *c, void *x, size_t size )
{
    if( size > c->bufsiz )
        Abort( "Output buffer is not large enough for output records (%d < %d).", c->bufsiz, size );
    else if( c->offset + size > c->bufsiz )
        eat_cookies( c );
        
    memcpy( (void *)(c->buf + c->offset), x, size );
    c->offset += size;
}

void      eat_cookies( CookieJar *c )
{
    fwrite( (void *)c->buf, 1, c->offset, c->stream );

    if( ferror(c->stream) )
        Abort( "Error condition encountered when outputing cookies to binary file." );

    c->offset = 0;
}

/*
 * Outputs a record of zeros as simply as possible.
 * Assumes nnul < nout and that pars, timec, and imat
 * are large enough spaces.  These spaces are overwritten
 * with 0's.
 *
 * This is a temporary fix.
 *
 */

static
void            output_null_record( const char* mesg, int voxel, long type,
                                    int nout, int tdim, int null_model, int nnul,
                                    double* pars, double* timec, double* imat )
{
    int         i, j;
    size_t      size;
    int         idiagnostics[MAX_DIAGNOSTICS];
    float       rdiagnostics[MAX_DIAGNOSTICS];
    static double* imat_lcl= NULL;
    static int imat_lcl_size= 0;

    if( type != DOUBLE_T && type != FLOAT_T )
        type = FLOAT_T;

    size = mriTypeSize( type );

    Warning( DIAGNOSTIC >= 1, "Null time-course encountered: %s, voxel %d", mesg, voxel );

    /* imat_lcl is a hack to work around an initialization bug- this
     * routine was getting called with imat uninitialized.  
     */
    if (imat_lcl_size < (nout*(nout+1))/2) {
      if (imat_lcl) Free(imat_lcl);
      imat_lcl_size= (nout*(nout+1))/2;
      imat_lcl= Calloc( imat_lcl_size, double );
    }

    /* Set Output Objects to all 0 */

    memset( pars,  0, nout*sizeof(double) );
    memset( timec, 0, tdim*size );
    memset( imat_lcl,  0, sizeof(double)*(nout*(nout+1))/2 );
    if (imat) memcpy( imat, imat_lcl, sizeof(double)*(nout*(nout+1))/2 );

    memset( idiagnostics, 0, sizeof(idiagnostics) );
    memset( rdiagnostics, 0, sizeof(rdiagnostics) );

    /* Output Null Record for full model */

    for( i = 0; i < Num_Products; i++ )
    {
        if( Product_Type[i] == ESTIMATES )
            output_cookie( ESTIMATES, nout, pars );
        else if( Product_Type[i] == STANDARD_ERRORS )
            output_cookie( STANDARD_ERRORS, nout, pars );
        else if( Product_Type[i] == RESIDUALS )
            output_cookie( RESIDUALS, tdim, size, timec );
        else if( Product_Type[i] == FITTED_VALUES )
            output_cookie( FITTED_VALUES, tdim, size, timec );
        else if( Product_Type[i] == COVARIANCES )
            output_cookie( COVARIANCES, (nout * (nout+1))/2, imat_lcl );
        else if( Product_Type[i] == DIAGNOSTICS )
        {
            j = Null_Model | Produce[STANDARD_ERRORS] | Produce[COVARIANCES];

            output_cookie( DIAGNOSTICS, 1 + j*2, 2, rdiagnostics, idiagnostics );
        }
    }

    /* Output Null Record for null model */

    if( null_model )
    {
        for( i = 0; i < Num_Products; i++ )
        {
            if( Product_Type[i] == ESTIMATES )
                output_cookie( ESTIMATES, nnul, pars );
            else if( Product_Type[i] == STANDARD_ERRORS )
                output_cookie( STANDARD_ERRORS, nnul, pars );
            else if( Product_Type[i] == RESIDUALS )
                output_cookie( RESIDUALS, tdim, size, timec );
            else if( Product_Type[i] == FITTED_VALUES )
                output_cookie( FITTED_VALUES, tdim, size, timec );
            else if( Product_Type[i] == COVARIANCES )
                output_cookie( COVARIANCES, (nnul * (nnul+1))/2, imat_lcl );
            else if( Product_Type[i] == DIAGNOSTICS )
                output_cookie( DIAGNOSTICS, 4, 2, rdiagnostics, idiagnostics );
        }
    }
}



/* Signal Handling */

static  int     current_voxel = -1;

void  set_current_voxel( int v )
{
    current_voxel = v;
}

void  handler( int sig )
{
    time_t   exectime;
    
    if( sig == SIGQUIT )
    {
        time( &exectime );
        
        Message( "Execution Aborted at %s while processing voxel %d.\n", ctime(&exectime), current_voxel );

        eat_cookies( &cookies );  /* Flush output buffer */
        
        exit( 0 );
    }
}

/*
 * Local Variables:        
 * mode: c
 * eval: (make-local-variable 'compile-command)
 * compile-command: "make -k vmpfx"
 * End:
 *
 */
