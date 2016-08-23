/************************************************************
 *                                                          *
 *  blind_minimize.c                                        *
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
 *  Original programming by Joel Welling 3-03               *
 ************************************************************/

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <time.h>
#ifdef DARWIN
#include <sys/time.h>
#endif
#include <math.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/resource.h>
#include "../fmri/lapack.h"
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"

static char rcsid[] = "$Id: blind_minimize.c,v 1.9 2007/03/21 23:54:26 welling Exp $";

/* Notes-
 */

#define DEFAULT_SEARCH_ALGORITHM "opt=praxis"

/* This structure defines the algorithm to be used.
 */
typedef struct alg_context_struct {
  Optimizer* opt;
  ScalarFunction* sf;
} Algorithm;

typedef struct par_struct {
  double* vals;
  int ndim;
  double mse;
} Par;

/* Globals.  Some are necessary because the optimization routines
 * expect to access data this way.
 */
static char* progname= NULL;
static int debug= 0;
static Algorithm gblAlg;
static Par gblPar;
static double gblScale;
static char script[512];

static void par_init( Par* par )
{
  par->ndim= 0;
  par->vals= NULL;
  par->mse= 0.0;
}

static void par_parse( Par* par, const char* str )
{
  int ndim= 0;
  int i;
  char* here;
  char* tmp;
  char* tok;
  double* vals= NULL;
  ndim= 0;
  here= strdup(str);
  while (strtok_r((ndim ? NULL : here), " :,", &tmp)) ndim++;
  if (ndim==0)
    Abort("%s: no tokens in start string <%s>!\n",progname,str);

  par->ndim= ndim;
  par->mse= 0.0;
  if (!(vals=(double*)malloc(ndim*sizeof(double))))
    Abort("%s: cannot allocate %d bytes!\n",ndim*sizeof(double));
  free(here);
  here= strdup(str);
  i= 0;
  while ((tok=strtok_r((i ? NULL : here), " :,", &tmp)) != NULL ) {
    vals[i]= atof(tok);
    i++;
  }
  par->vals= vals;
  free(here);

}

static void par_print( const char* hdr, Par* par, FILE* ofile )
{
  int i;
  fprintf(ofile,"%s: ",hdr);
  for (i=0;i<par->ndim;i++) 
    fprintf(ofile,"%lg ",par->vals[i]);
  fprintf(ofile,"-> %18.10lg\n",par->mse);
}

static Optimizer* parseOptMethod(char* s)
{
  /* Two-sided 95% T cutoffs for DF up to 20 */
  static float tcrit[]= { 0.0, 12.706176098809, 4.30263617020413,
			  3.18244545308107, 2.77644502408949, 2.57058182376418,
			  2.44680441622004, 2.36457539633475, 2.30597941403506,
			  2.26214359559273, 2.2281309162355, 2.20098027417758,
			  2.17880969240071, 2.16036657035061, 2.14478525936996,
			  2.13144854260425, 2.11990457994238, 2.10981505255424,
			  2.10092165066567, 2.09302376167972, 2.08596322490106 
  };

  if (gblPar.ndim> (sizeof(tcrit)/sizeof(float))-1) 
    Abort("%s: tcrit for %d degrees of freedom is not compiled in!\n",
	  progname, gblPar.ndim);

  if (!strcasecmp(s,"nelmin")) 
    return createNelminOptimizer(0.02,1.0,0);
  else if (!strcasecmp(s,"nelmin_t")) 
    return createNelminTOptimizer(0.02,1.0,tcrit[gblPar.ndim],0);
  else if (!strcasecmp(s,"praxis")) 
    return createPraxisOptimizer(0.000001,1.0);
  else if (!strcasecmp(s,"none")) return createNoneOptimizer();
  else return NULL;
}

static int parseAlgString(char* search_string) 
{
  char* work;
  char* tok;
  double tol= 0.0;
  int tolSet= 0;
  double scale= 0.0;
  int scaleSet= 0;

  work= strdup(search_string); /* make guaranteed-writable copy */
  tok= strtok(work, " ,+&:;");
  while (tok != NULL) {

    if (!strncasecmp(tok,"opt=",4)) {
      if (gblAlg.opt != NULL) {
	gblAlg.opt->destroySelf(gblAlg.opt);
      }
      if ((gblAlg.opt=parseOptMethod(tok+4))==NULL)
	{ free(work); return 0; }
    }
    else if (!strncasecmp(tok,"optol=",6)) {
      tol= atof(tok+6);
      tolSet= 1;
    }
    else if (!strncasecmp(tok,"opscale=",8)) {
      scale= atof(tok+8);
      scaleSet= 1;
    }
    else { free(work); return 0; }
    
    tok= strtok(NULL, " ,+&:;");
  }
  if (gblAlg.opt) {
    if (tolSet) gblAlg.opt->setTol(gblAlg.opt,tol);
    if (scaleSet) gblAlg.opt->setScale(gblAlg.opt,scale);
  }

  free(work);
  return 1;
}

static char* getAlgInfoString()
{
  char scratch[128];
  static char result[512];
  int offset= 0;

  result[0]= '\0';

  sprintf(scratch,"opt=%s,",gblAlg.opt->getMethodName(gblAlg.opt));
  strncat(result,scratch,sizeof(result)-(offset+1));
  offset += strlen(scratch);
  if (offset>=sizeof(result)-1) return result;  
  sprintf(scratch,"optol=%f,",gblAlg.opt->getTol(gblAlg.opt));
  strncat(result,scratch,sizeof(result)-(offset+1));
  offset += strlen(scratch);
  if (offset>=sizeof(result)-1) return result;  
  sprintf(scratch,"opscale=%f,",gblAlg.opt->getScale(gblAlg.opt));
  strncat(result,scratch,sizeof(result)-(offset+1));
  offset += strlen(scratch);
  if (offset>=sizeof(result)-1) return result;  

  return result;
}

static void spawn_and_listen(char** cmd, char* result, 
			     int resultLength)
{
  int myPipe[2];
  pid_t kidPid;

#ifdef never
  fprintf(stderr,"Sending command <%s>\n",cmd[0]);
#endif

  if (pipe(myPipe)) {
    perror("pipe failed");
    Abort("%s: unable to create a pipe!\n",progname);
  }

  result[resultLength-1]= '\0';
  resultLength -= 1;

  if ((kidPid=fork()) == 0) {
    close(myPipe[0]);
    dup2(myPipe[1],STDOUT_FILENO);
    execvp(cmd[0], cmd);
    /* Should never reach this point */
    perror("Child execvp failed");
    Abort("%s: could not execute script <%s>\n",progname,cmd[0]);
  }
  else {
    int status;
    int len;
    
    close(myPipe[1]);
    if(waitpid( kidPid, &status, 0 )<0) {
      perror("waitpid failed");
      Abort("%s: something bad happend to the child process!\n",progname);
    }
    if (WEXITSTATUS(status) != 0) {
      Abort("%s: child exited with error status %d!\n",progname,status);
    }
    len= read(myPipe[0], result, resultLength);
    if (len<0) {
      perror("Read failed");
      Abort("%s: unable to read from child process!\n",progname);
    }
    result[len]= '\0';
    
    if (debug) {
      char* here;
      if ((here=strchr(result,'\n')) != NULL) *here= '\0';
#ifdef never
      fprintf(stderr,"Child process responded <%s>\n",result);
#endif
    }
    close(myPipe[0]);
  }
}

/* This routine calculates the current weight function based on the global
 * task info.  This value is ultimately returned by mse().
 */
static double calcChiSqr(Par* par)
{
  char tbuf[128];
  char** cmd;
  int i;

  if (!(cmd= (char**)malloc((par->ndim+3)*sizeof(char*))))
    Abort("%s: Unable to allocate %d bytes!\n",(par->ndim+3)*sizeof(char*));

  cmd[0]= strdup(script);
  for (i=0; i<par->ndim; i++) {
    snprintf(tbuf, sizeof(tbuf), "%lf", par->vals[i]);
    cmd[i+1]= strdup(tbuf);
  }
  cmd[par->ndim+1]= NULL;

  spawn_and_listen( cmd, tbuf, sizeof(tbuf) );

  for (i=0; i<=par->ndim; i++) free(cmd[i]);
  free(cmd);

  return( gblScale*atof(tbuf) );
}

/* Returns mean squared-error between (adjusted) fixed image  */
/*   and (adjusted) register image, the criterion to be       */
/*   minimized by the minization op for registration */
static double mse( const double* guess, const int npar, void* userHook )
{
  int i;
  Par par;
  
  if (npar != gblPar.ndim)
    Abort("%s: mse found %d parameters, not %d!\n",progname,
	  npar, gblPar.ndim);
  par.ndim= npar;
  par.vals= (double*)guess;
  par.mse= calcChiSqr(&par);
    
  if (debug) {
    par_print("mse", &par, stderr);
  }
  
  return( par.mse );
}

static void restrt( const double* guess, const int npar, void* userHook )
{
  if (debug) {
    Par par;
    int i;
    if (npar != gblPar.ndim)
      Abort("%s: restrt found %d parameters, not %d!\n",progname,
	    npar, gblPar.ndim);
    par.ndim= npar;
    par.vals= (double*)guess;
    par.mse= 0.0;
    par_print("restrt", &par, stderr);
  }
}

static void buildTimeString( char* time_string, 
			     struct rusage* start, struct rusage* end )
{
  long s_usec= end->ru_stime.tv_usec - start->ru_stime.tv_usec;
  long s_sec= end->ru_stime.tv_sec - start->ru_stime.tv_sec;
  long u_usec= end->ru_utime.tv_usec - start->ru_utime.tv_usec;
  long u_sec= end->ru_utime.tv_sec - start->ru_utime.tv_sec;

  if (s_usec < 0) {
    s_sec -= 1;
    s_usec += 1000000;
  }
  if (u_usec < 0) {
    u_sec -= 1;
    u_usec += 1000000;
  }
  sprintf(time_string,"%d.%06du %d.%06ds",u_sec,u_usec,s_sec,s_usec);
}

int main( int argc, char* argv[] ) 
{
  char startString[512];
  char search_string[512]; 
  struct rusage start_rusage, end_rusage;
  char run_time_string[256];
  int script_given= 0;
  int startString_given= 0;

  progname= argv[0];
  par_init( &gblPar );

  /* Print version number */
  Message( "# %s\n", rcsid );

  /* Check to see if help was requested */
  if( ( argc > 1 ) && !strcmp( argv[1], "-help" ) )
    {
      if( argc == 2 )
	Help( "selecttopic" );
      else
	Help( argv[2] );
    }

  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Get filenames */
  script_given= cl_get( "script", "%option %s", script );
  startString_given= cl_get( "start", "%option %s", startString );
  cl_get("scale", "%option %lf[%]", 1.0, &gblScale);
  if (cl_present("debug")) debug= 1;
  cl_get( "algorithm|alg", "%option %s[%]", DEFAULT_SEARCH_ALGORITHM, 
	  search_string );

  if (cl_cleanup_check()) {
    int i;
    fprintf(stderr,"%s: invalid argument in command line:\n    ",argv[0]);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(-1);
  }

  if (!script_given) {
    fprintf(stderr,"%s: required script name not given!\n",argv[0]);
    Help( "usage" );
    exit(-1);
  }
  if (!startString_given) {
    fprintf(stderr,"%s: required start point not given!\n",argv[0]);
    Help( "usage" );
    exit(-1);
  }

  /*** End command-line parsing ***/

  /* Parse the starting point */
  par_parse(&gblPar, startString);

  /* Parse the search method string to define the search pattern.  This
   * and the construction of the smoother define the global algorithm. 
   */
  if (!parseAlgString(DEFAULT_SEARCH_ALGORITHM)) /* set defaults */
    Abort("%s: internal error: invalid default algorithm!\n",argv[0]);
  if (!parseAlgString(search_string)) {
    fprintf(stderr,"%s: invalid search method specified.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  if (debug) 
    gblAlg.opt->setDebugLevel(gblAlg.opt,10);
  Message("# algorithm: %s\n",getAlgInfoString());
  Message("# optimizer details: %s\n",gblAlg.opt->getStringRep(gblAlg.opt));
  gblAlg.sf= buildSimpleScalarFunction(mse,restrt,gblPar.ndim,NULL);

  /* Let's keep track of optimization time */
  getrusage(RUSAGE_CHILDREN,&start_rusage);

  /* Run the optimizaiton */
  (void)(gblAlg.opt->go(gblAlg.opt, gblAlg.sf, 
			gblPar.vals, gblPar.ndim, &(gblPar.mse)));

  getrusage(RUSAGE_CHILDREN,&end_rusage);
  buildTimeString(run_time_string,&start_rusage,&end_rusage);

  par_print("Final",&gblPar,stdout);

  Message("# Optimization complete in %s\n",run_time_string);

  return 0;
}

