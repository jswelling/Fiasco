/************************************************************
 *                                                          *
 *  rpn_engine.c                                          *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1997 Department of Statistics,         *
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
 *  Original programming by Joel Welling, 8/97              *
 ************************************************************/

/* Notes-
 */

/* We define the following macro to help old versions of
 * Linux find signbit().
 */
#if ( defined(LINUX) || defined(LINUXI386) || defined(LINUXX86_64) )
#define _ISOC99_SOURCE
#endif

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <time.h>
#include "mri.h"
#include "fmri.h"
#include "misc.h"
#include "stdcrg.h"
#include "dcdflib.h"
#include "lapack.h"

#ifdef DARWIN
#define finite( foo ) isfinite( foo )
#endif

#define MAX_INPUT_FILES 20
#define MAX_STACK 300
#define MAX_SCRIPT_CHARS 512
#define INITIAL_CODE_BLOCK_SIZE 100

static char rcsid[] = "$Id: rpn_engine.c,v 1.11 2007/03/21 23:50:20 welling Exp $";

static ParsePair opnames[]= {
  {"+",OP_PLUS}, {"-",OP_MINUS}, {"*",OP_MULT}, {"/",OP_DIV},
  {"sqrt",OP_SQRT}, {"ln",OP_LN}, {"exp",OP_EXP}, {"abs",OP_ABS},
  {"dup",OP_DUP}, {"swap",OP_SWAP}, {"pop",OP_POP}, 
  {">",OP_GT}, {"<",OP_LT}, {"==",OP_EQ}, {">=",OP_GE},
  {"<=",OP_LE}, {"!=",OP_NE}, {"and",OP_AND}, {"or",OP_OR},
  {"if_keep",OP_IF_KEEP}, {"cbeta",OP_CBETA}, {"inv_cbeta", OP_INV_CBETA},
  {"cbinom",OP_CBINOM}, {"inv_cbinom",OP_INV_CBINOM},
  {"cchisqr",OP_CCHISQR}, {"inv_cchisqr",OP_INV_CCHISQR},
  {"cf",OP_CF}, {"inv_cf",OP_INV_CF}, 
  {"ct",OP_CT}, {"inv_ct",OP_INV_CT},
  {"cgamma",OP_CGAMMA}, {"inv_cgamma",OP_INV_CGAMMA},
  {"cnormal",OP_CNORMAL}, {"inv_cnormal",OP_INV_CNORMAL},
  {"cpoisson",OP_CPOISSON}, {"inv_cpoisson",OP_INV_CPOISSON},
  {"fcbeta",OP_FCBETA}, {"inv_fcbeta", OP_INV_FCBETA},
  {"fcbinom",OP_FCBINOM}, {"inv_fcbinom",OP_INV_FCBINOM},
  {"fcchisqr",OP_FCCHISQR}, {"inv_fcchisqr",OP_INV_FCCHISQR},
  {"fcf",OP_FCF}, {"inv_fcf",OP_INV_FCF}, 
  {"fct",OP_FCT}, {"inv_fct",OP_INV_FCT},
  {"fcgamma",OP_FCGAMMA}, {"inv_fcgamma",OP_INV_FCGAMMA},
  {"fcnormal",OP_FCNORMAL}, {"inv_fcnormal",OP_INV_FCNORMAL},
  {"fcpoisson",OP_FCPOISSON}, {"inv_fcpoisson",OP_INV_FCPOISSON},
  {"rand",OP_RANDOM}, 
  {"sin",OP_SIN}, {"cos",OP_COS}, {"tan",OP_TAN},
  {"asin",OP_ASIN}, {"acos",OP_ACOS}, {"atan",OP_ATAN},
  {"min",OP_MIN}, {"max",OP_MAX},
  {"missing",OP_MISSING}, {"is_finite",OP_IS_FINITE},
  {"round",OP_ROUND}, {"floor",OP_FLOOR}, {"ceiling",OP_CEILING},
  {"%",OP_MODULO},
  {"cx_mag",OP_CX_MAG}, {"cx_phase",OP_CX_PHASE}, 
  {"cx_+",OP_CX_PLUS}, {"cx_-",OP_CX_MINUS},
  {"cx_*",OP_CX_MULT}, {"cx_/",OP_CX_DIV}, {"cx_conj",OP_CX_CONJ},
  {"cx_tophase",OP_CX_TOPHASEREP}, {"cx_fmphase",OP_CX_FMPHASEREP},
  {"cx_if_keep",OP_CX_IF_KEEP}, {"signbit",OP_SIGNBIT},
  {NULL,(Op)0} /* must end with null! */
};

static MacroPair macronames[]= {
  {"foldp","dup,1,swap,-,-1,*,swap,dup,0.5,>=,if_keep"},
  {"inv_foldp","dup,dup,1,+,swap,signbit,if_keep"},
  {NULL, NULL} /* must end with null! */  
};

static Program* compile( RpnEngine* re, char* string, Clock* clock );
static void destroyProgram(Program* prog);

static void setErrorStr(RpnEngine* re, const char* fmt, ... )
{
  char buf[256];
  va_list args;
  va_start( args, fmt );
  if (re->errorString) free(re->errorString);
  vsnprintf(buf,sizeof(buf)-1,fmt,args);
  buf[sizeof(buf)-1]= '\0';
  re->errorString= strdup(buf);
  va_end(args);
}

void rpnSetOutputFlag(RpnEngine* re, int flag)
{
  re->outfileFlag= flag;
}

void rpnSetVerbose(RpnEngine* re, int flag)
{
  re->verboseFlag= flag;
}

void rpnSetDebug(RpnEngine* re, int flag)
{
  re->debugFlag= flag;
}

void rpnSetComplex(RpnEngine* re, int flag)
{
  re->complexFlag= flag;
}

static Clock* clock_init(RpnEngine* re, int whichInput)
{
  Clock* clock= NULL;
  char* crunner;
  int i;
  long long l;

  if (!(clock=(Clock*)malloc(sizeof(Clock))))
    Abort("rpn_engine: unable to allocate %d bytes!\n",sizeof(Clock));

  clock->string= strdup( re->getDimensionsCB(whichInput,re->usrHook) );

  if (!(clock->limits= (int*)malloc(strlen(clock->string)*sizeof(int))))
    Abort("rpn_engine: unable to allocate %d bytes!\n",
	  strlen(clock->string)*sizeof(int));
  if (!(clock->strides= (long long*)malloc(strlen(clock->string)
					   *sizeof(long long))))
    Abort("rpn_engine: unable to allocate %d bytes!\n",
	  strlen(clock->string)*sizeof(long long));

  /* Figure out how many instances of this vector */
  i=0;
  l= 1;
  clock->t_offset= -1;
  clock->z_offset= -1;
  for (crunner= clock->string; *crunner; crunner++) {
    if (re->complexFlag) {
      if (crunner==clock->string && *crunner=='v' /* first dimension is v */
	  && re->getDimExtentCB(whichInput,*crunner,re->usrHook)==2)
	clock->limits[i]= 1; /* count complex pair as 1 */
      else clock->limits[i]= re->getDimExtentCB(whichInput,*crunner,
						re->usrHook);

    }
    else {
      clock->limits[i]= re->getDimExtentCB(whichInput,*crunner,re->usrHook);
    }
    clock->strides[i]= l;
    l *= clock->limits[i];
    if (*crunner=='t') clock->t_offset= i;
    else if (*crunner=='z') clock->z_offset= i;
    i++;    
  }
  if (clock->t_offset>=0 && clock->z_offset>=0) clock->missing_enable= 1;
  else clock->missing_enable= 0;

  return clock;
}

static int check_dims_acceptable(RpnEngine* re, int whichInput)
{
  char* dimstr;
  int offset;
  int extent;
  char* runner;

  dimstr= strdup( re->getDimensionsCB(whichInput,re->usrHook) );

  /* Rules to enforce:
   * 1) In complex mode, the only 'v' must be the first dimension and
   *    it must have an extent of 2.
   */

  if (re->complexFlag) {
    int vdim= 0;
    char* vloc= strchr(dimstr,'v');
    if (vloc) {
      if (vloc != dimstr) {
	setErrorStr(re, "In complex mode, v dimension is not first in input %d!",
		    whichInput);
	free(dimstr);
	return 0;
      }
      vdim= re->getDimExtentCB(whichInput,'v',re->usrHook);
      if (vdim != 1 && vdim != 2) {
	setErrorStr(re, 
		    "In complex mode, input %d has a v extent greater than 2!",
		    whichInput);
	free(dimstr);
	return 0;
      }
    }
  }
  
  free(dimstr);
  return 1;
}

static int check_dims_compatible(RpnEngine* re, int whichInput, 
				 Clock* masterClock)
{
  char* dimstr;
  int offset;
  int extent;
  char* runner;
  int* limits= NULL;
  int clock_offset;
  char* clock_runner;
  int allow_trivial_only;

  dimstr= strdup( re->getDimensionsCB(whichInput,re->usrHook) );

  if (!(limits= (int*)malloc(strlen(dimstr)*sizeof(int))))
    Abort("rpn_engine: unable to allocate %d bytes!\n",
	  strlen(dimstr)*sizeof(int));
  for (offset=0; dimstr[offset]; offset++) {
    limits[offset]= re->getDimExtentCB(whichInput,dimstr[offset],
				       re->usrHook);
    if (re->complexFlag 
	&& offset==0 && dimstr[offset]=='v' && limits[offset]==2) {
      limits[offset]= 1; /* see rule (6) below */
    }
  }

  /* Rules to enforce:
   * 1) If a dimension is present in both, it must have the same extent in
   * both or it must be trivial in the test string.
   * 2) If a dimension exists in the test string but not in the clock string,
   * it must be trivial.
   * 3) If a dimension is trivial in the test string but not in the clock
   * string, it can be followed only by trivial dimensions in the 
   * test string.
   * 4) Nontrivial dimensions must appear in the same order in the clock and
   * test strings.
   * 5) If a dimension exists in the clock string but not in the test string,
   * and it is non-trivial in the clock string, it must be followed only by 
   * trivial dimensions in the test string.
   * 6) In complex mode, if the first dim is v and its extent is 2, 
   * treat its length as 1.
   *
   * Scan the test string, checking these rules. 
   */
  allow_trivial_only= 0;
  clock_offset= 0;
  for (offset=0; dimstr[offset]; offset++) {
    extent= limits[offset];
    if (extent != 1 && allow_trivial_only) {
      /* rule 3 */
      setErrorStr(re, "non-trivial dimension %c follows trivial mismatched dimension!",
		  dimstr[offset]);
      free(limits);
      free(dimstr);
      return 0;
    }
    clock_runner= strchr(masterClock->string,dimstr[offset]);
    if (!clock_runner) {
      /* This dimension doesn't occur in the clock */
      if (extent != 1) {
	/* rule 1 */
	setErrorStr(re, "non-trivial unmatched dimension %c!", dimstr[offset]);
	free(limits);
	free(dimstr);
	return 0;
      }
      else {
	/* OK by rule 1 */
      }
    }
    else {
      if (clock_runner < masterClock->string+clock_offset) {
	/* Dimension order mismatch, rule 4 */
	setErrorStr(re, "dimension %c is out of order!",dimstr[offset]);
	free(limits);
	free(dimstr);
	return 0;
      }
      clock_offset= clock_runner - masterClock->string;
      if (extent != masterClock->limits[clock_offset]) {
	/* Dimension order OK, but extents differ */
	if (extent != 1) {
	  /* Fails by rule 2 */
	  setErrorStr(re, "non-trivial mismatched dimension %c!",dimstr[offset]);
	  free(limits);
	  free(dimstr);
	  return 0;
	}
	else {
	  allow_trivial_only= 1; /* rule 3 */
	}
      }
    }
  }

  /* Now scan the clock string, checking any rules remaining. */
  allow_trivial_only= 0;
  offset= -1;
  for (clock_runner=masterClock->string; *clock_runner; clock_runner++) {
    if (masterClock->limits[clock_runner-masterClock->string] != 1) {
      char* here= strchr(dimstr,*clock_runner);
      if (!here) {
	int i;
	for (i=offset+1; dimstr[i]; i++) {
	  if (limits[i] != 1) {
	    /* rule 5 */
	    setErrorStr(re, "nontrivial dimension %c follows unmatched nontrivial dimension!",
			dimstr[i]);
	    free(limits);
	    free(dimstr);
	    return 0;
	  }
	}
      }
      else {
	offset= here-dimstr;
	allow_trivial_only= 1;
      }
    }
  }

  free(limits);
  free(dimstr);
  return 1;
}

static Program* parse_macro( RpnEngine* re, char* tok, Clock* clock )
{
  MacroPair* macroPair= macronames;
  while (macroPair->string) {
    if (!strcasecmp(tok,macroPair->string)) {
      char* codeCopy;
      Program* result;
      if (re->debugFlag)
	fprintf(stderr,"macro expansion <%s> -> <%s>\n",tok,macroPair->code);
      codeCopy= strdup(macroPair->code); /* writable copy */
      result= compile( re, codeCopy, clock );
      free(codeCopy);
      return result;
    }
    macroPair++;
  }
  return NULL;
}

static void parse_instruction( RpnEngine* re,
			       Instruction* inst, char* tok, Clock* clock )
{
  ParsePair* parsePair;
  int matched;

  matched= 0;
  parsePair= opnames;
  while (parsePair->string) {
    if (!strcasecmp(tok,parsePair->string)) {
      inst->op= parsePair->op;
      inst->param.l= 0;
      if (re->debugFlag) 
	fprintf(stderr,"%s (%d)\n",parsePair->string,(int)parsePair->op);
      matched= 1;
      break;
    }
    parsePair++;
  }
  
  if (!matched) {
    if (!(strncasecmp(tok,"if_print_",9))) {
      char* tail= tok + 9; /* after if_print_ */
      long count= atoi(tail);
      if (count<1) 
	Abort("rpn_engine: invalid ""if_print_"" expression!\n");
      inst->op= OP_IF_PRINT;
      inst->param.l= count;
      if (re->debugFlag) 
	fprintf(stderr,"OP_IF_PRINT %d\n",(int)(inst->param.l));
      matched= 1;
    }
  }
  
  if (!matched) {
    if (!(strncasecmp(tok,"rot_",4))) {
      char* tail= tok + 4; /* after rot_ */
      long count= atoi(tail);
      inst->op= OP_ROT;
      inst->param.l= count;
      if (re->debugFlag) fprintf(stderr,"OP_ROT %d\n",(int)(inst->param.l));
      matched= 1;
    }
  }
  
  if (!matched) {
    if (!(strncasecmp(tok,"switch_",7))) {
      char* tail= tok + 7; /* after rot_ */
      long count= atoi(tail);
      inst->op= OP_SWITCH;
      inst->param.l= count;
      if (re->debugFlag) fprintf(stderr,"OP_SWITCH %d\n",(int)(inst->param.l));
      matched= 1;
    }
  }
  
  if (!matched) {
    if (!(strcmp(tok,"nan"))) {
      inst->op= OP_CONST;
      inst->param.f= log(-1.0);
      if (re->debugFlag) fprintf(stderr,"OP_CONST %g\n",inst->param.f);
      matched= 1;
    }
  }

  if (!matched) {
    if (!(strcmp(tok,"inf"))) {
      inst->op= OP_CONST;
      inst->param.f= -log(0.0);
      if (re->debugFlag) fprintf(stderr,"OP_CONST %g\n",inst->param.f);
      matched= 1;
    }
  }

  if (!matched) {
    if (!(strcmp(tok,"ninf"))) {
      inst->op= OP_CONST;
      inst->param.f= log(0.0);
      if (re->debugFlag) fprintf(stderr,"OP_CONST %g\n",inst->param.f);
      matched= 1;
    }
  }

  if (!matched) {
    if (!(strcmp(tok,"pi"))) {
      inst->op= OP_CONST;
      inst->param.f= M_PI;
      if (re->debugFlag) fprintf(stderr,"OP_CONST %g\n",inst->param.f);
      matched= 1;
    }
  }

  if (!matched) {
    if (!(strcmp(tok,"eps"))) {
      inst->op= OP_CONST;
      inst->param.f= DLAMCH("e");
      if (re->debugFlag) fprintf(stderr,"OP_CONST %g\n",inst->param.f);
      matched= 1;
    }
  }

  if (!matched) {
    if (!(strcmp(tok,"sfmin"))) {
      inst->op= OP_CONST;
      inst->param.f= DLAMCH("sfmin");
      if (re->debugFlag) fprintf(stderr,"OP_CONST %g\n",inst->param.f);
      matched= 1;
    }
  }

  if (!matched) {
    if (!(strcmp(tok,"rmax"))) {
      inst->op= OP_CONST;
      inst->param.f= DLAMCH("o");
      if (re->debugFlag) fprintf(stderr,"OP_CONST %g\n",inst->param.f);
      matched= 1;
    }
  }

  if (!matched) {
    if (*tok=='$') {
      int which_in= atoi(tok+1);
      if (which_in==0) {
	char* here;
	/* Might be error, might be variable */
	if (strlen(tok)==2 && ((here=strchr(clock->string,tok[1]))!=NULL)) { 
	  /* clock reference */
	  inst->op= OP_CLOCK;
	  inst->param.l= (long)(here-clock->string);
	  if (re->debugFlag) fprintf(stderr,"OP_CLOCK %d\n",(int)inst->param.l);
	}
	else if (strlen(tok)==5 && !strcmp(tok+2,"dim")
		 && ((here=strchr(clock->string,tok[1]))!=NULL)) { 
	  /* dim reference */
	  inst->op= OP_DIM;
	  inst->param.l= (long)(here-clock->string);
	  if (re->debugFlag) fprintf(stderr,"OP_DIM %d\n",(int)inst->param.l);
	}
	else {
	  Abort("rpn_engine: expression <%s> is invalid!\n",tok);
	}
      }
      else if ((which_in<0) || (which_in > re->nInputs)) {
	Abort("rpn_engine: expression <%s> accesses invalid input file!\n",
	     tok);
      }
      else {
	/* It's a load request */
	inst->op= OP_LOAD;
	inst->param.l= which_in - 1;
	if (re->debugFlag) fprintf(stderr,"OP_LOAD %d\n",(int)inst->param.l);
      }
    }
    else {
      /* This one must be last! */
      /* Try to handle this as a double */
      inst->op= OP_CONST;
      if (sscanf(tok,"%lg",&(inst->param.f)) != 1)
	Abort("rpn_engine: unable to recognize the token <%s>\n",tok);
      if (re->debugFlag) fprintf(stderr,"OP_CONST %f\n",inst->param.f);
    }
  }
}

static Program* compile( RpnEngine* re, char* string, Clock* clock )
{
  Program* result;
  Instruction* code;
  Instruction* runner;
  long buf_length;
  char* tok;
  char* ptr= NULL;

  if (!(result= (Program*)malloc(sizeof(Program))))
    Abort("rpn_engine: compile: unable to allocate %d bytes!\n",
	  sizeof(Program));

  if (!(result->code= code= 
	(Instruction*)malloc(INITIAL_CODE_BLOCK_SIZE*sizeof(Instruction))))
    Abort("rpn_engine: compile: unable to allocate %d bytes!\n",
	  INITIAL_CODE_BLOCK_SIZE*sizeof(Instruction));
  buf_length= INITIAL_CODE_BLOCK_SIZE;

  tok= strtok_r(string," ,",&ptr);
  runner= code;
  while (tok) {
    long nInstructionsThisTok;
    Program* macroCode= parse_macro( re, tok, clock );
    if (macroCode) nInstructionsThisTok= macroCode->code_length;
    else nInstructionsThisTok= 1;

    /* Grow the instruction buffer if necessary */
    while ((runner - code)+nInstructionsThisTok > buf_length) {
      if (!(result->code= code= 
	    (Instruction*)realloc(code, 2*buf_length*sizeof(Instruction))))
	Abort("rpn_engine: compile: unable to allocate %d bytes!\n",
	      2*buf_length*sizeof(Instruction));
      buf_length *= 2;
    }
    
    if (macroCode) {
      long i;
      for (i=0; i<macroCode->code_length; i++) 
	*runner++= macroCode->code[i];
      destroyProgram(macroCode);
    }
    else {
      parse_instruction(re, runner, tok, clock);
      runner++;
    }
    tok= strtok_r(NULL," ,",&ptr);
  }

  result->code_length= runner-code;
  return result;
}

static Program* compile_file( RpnEngine* re, const char* fname, Clock* clock )
{
  Program* result;
  Instruction* code;
  Instruction* runner;
  long buf_length;
  char* tok;
  FILE* infile;
  char inbuf[512];
  char* ptr;

  if (!strcmp(fname,"-")) infile= stdin;
  else if (!(infile= fopen(fname,"r"))) {
    perror("Error opening script");
    Abort("rpn_engine: unable to open script file %s\n",fname);
  }

  if (!(result= (Program*)malloc(sizeof(Program))))
    Abort("rpn_engine: compile: unable to allocate %d bytes!\n",
	  sizeof(Program));

  if (!(result->code= code= 
	(Instruction*)malloc(INITIAL_CODE_BLOCK_SIZE*sizeof(Instruction))))
    Abort("rpn_engine: compile: unable to allocate %d bytes!\n",
	  INITIAL_CODE_BLOCK_SIZE*sizeof(Instruction));
  buf_length= INITIAL_CODE_BLOCK_SIZE;

  runner= code;
  while (!feof(infile)) {
    if (!fgets(inbuf,sizeof(inbuf),infile)) break;
    if (inbuf[strlen(inbuf)-1]=='\n') inbuf[strlen(inbuf)-1]= '\0';
    tok= strtok_r(inbuf," ,",&ptr);
    while (tok) {
      long nInstructionsThisTok;
      Program* macroCode= parse_macro( re, tok, clock );
      if (macroCode) nInstructionsThisTok= macroCode->code_length;
      else nInstructionsThisTok= 1;

      /* Grow the instruction buffer if necessary */
      while ((runner - code)+nInstructionsThisTok > buf_length) {
	if (!(result->code= code= 
	      (Instruction*)realloc(code, 2*buf_length*sizeof(Instruction))))
	  Abort("rpn_engine: compile_file: unable to allocate %d bytes!\n",
		2*buf_length*sizeof(Instruction));
	buf_length *= 2;
      }
      
      if (macroCode) {
	long i;
	for (i=0; i<macroCode->code_length; i++) 
	  *runner++= macroCode->code[i];
	destroyProgram(macroCode);
      }
      else {
	parse_instruction(re, runner, tok, clock);
	runner++;
      }
      tok= strtok_r(NULL," ,",&ptr);
    }
  }

  result->code_length= runner-code;
  if (infile != stdin)
    if (fclose(infile)) perror("Error closing script (ignored)");
  return result;
}

#define BAILOUT(msg) { setErrorStr(re, msg); return NULL; }
#define BAILOUT1(msg,msg1) { setErrorStr(re, msg,msg1); return NULL; }
#define BAILOUT2(msg,msg1,msg2) { setErrorStr(re, msg,msg1,msg2); return NULL; }
#define BAILOUT3(msg,msg1,msg2,msg3) \
{ setErrorStr(re, msg,msg1,msg2,msg3); return NULL; }

double* rpnRun( RpnEngine* re, long length, long long offset )
{
  double* stack_first= &(re->stackBlock[0][0]);
  double* stack_last= (double*)(&(re->stackBlock[MAX_STACK-1][0]));
  double* stack_top;
  Instruction* runner;
  int i;
  Program* prog= re->program;
  Clock* clock= re->clock;

  if (length>RPN_CHUNKSIZE) 
    BAILOUT2("Internal error; rpnRun request of %d values vs. chunksize %d!",
	     length,RPN_CHUNKSIZE);

  stack_top= stack_first - RPN_CHUNKSIZE;

  for (runner= prog->code; (runner - prog->code) < prog->code_length; 
       runner++) {
    switch (runner->op) {
    case OP_CLOCK:
      {
	/* fprintf(stderr,"exec: OP_CLOCK %d\n",(int)runner->param.l); */
	if (stack_top>=stack_last) 
	  BAILOUT("stack overflow on push of clock!\n");
	stack_top += RPN_CHUNKSIZE;
	for (i=0; i<length; i++) {
	  int loc= (((offset + i)/clock->strides[runner->param.l])
	    % clock->limits[runner->param.l]);
	  *(stack_top+i)= (double)(loc);
	}
      }
    break;
    case OP_DIM:
      {
	/* fprintf(stderr,"exec: OP_DIM %d\n",(int)runner->param.l); */
	if (stack_top>=stack_last) 
	  BAILOUT("stack overflow on push of dimension!\n");
	stack_top += RPN_CHUNKSIZE;
	for (i=0; i<length; i++) 
	  *(stack_top+i)= (double)(clock->limits[runner->param.l]);
      }
    break;
    case OP_LOAD: 
      {
	/* fprintf(stderr,"exec: OP_LOAD %d\n",(int)runner->param.l); */
	if (stack_top>=stack_last) 
	  BAILOUT("stack overflow on load!\n");
	if (re->complexFlag) {
	  stack_top += 2*RPN_CHUNKSIZE;
	  re->getInputComplexCB((int)runner->param.l, length, offset,
				stack_top-RPN_CHUNKSIZE, stack_top,
				re->usrHook);
	}
	else {
	  stack_top += RPN_CHUNKSIZE;
	  re->getInputCB((int)runner->param.l, length, offset, stack_top,
			 re->usrHook);
	}
      }
    break;
    case OP_PLUS:
      {
	int i;
	/* fprintf(stderr,"exec: OP_PLUS\n"); */
	if (stack_top <= stack_first)
	  BAILOUT("stack underflow on plus!\n");
	stack_top -= RPN_CHUNKSIZE;
	for (i=0; i<length; i++) *(stack_top+i) += *(stack_top+RPN_CHUNKSIZE+i);
      }
    break;
    case OP_MINUS:
      {
	int i;
	/* fprintf(stderr,"exec: OP_MINUS\n"); */
	if (stack_top <= stack_first)
	  BAILOUT("stack underflow on minus!\n");
	stack_top -= RPN_CHUNKSIZE;
	for (i=0; i<length; i++) *(stack_top+i) -= *(stack_top+RPN_CHUNKSIZE+i);
      }
    break;
    case OP_MULT:
      {
	int i;
	/* fprintf(stderr,"exec: OP_MULT\n"); */
	if (stack_top <= stack_first)
	  BAILOUT("stack underflow on mult!\n");
	stack_top -= RPN_CHUNKSIZE;
	for (i=0; i<length; i++) *(stack_top+i) *= *(stack_top+RPN_CHUNKSIZE+i);
      }
    break;
    case OP_DIV:
      {
	int i;
	/* fprintf(stderr,"exec: OP_DIV\n"); */
	if (stack_top <= stack_first)
	  BAILOUT("stack underflow on div!\n");
	stack_top -= RPN_CHUNKSIZE;
	for (i=0; i<length; i++) *(stack_top+i) /= *(stack_top+RPN_CHUNKSIZE+i);
      }
    break;
    case OP_SQRT:
      {
	int i;
	/* fprintf(stderr,"exec: OP_SQRT\n"); */
	if (stack_top<stack_first) 
	  BAILOUT("stack underflow on sqrt!\n");
	for (i=0; i<length; i++) *(stack_top+i) = sqrt(*(stack_top+i));
      }
    break;
    case OP_LN:
      {
	int i;
	/* fprintf(stderr,"exec: OP_LN\n"); */
	if (stack_top<stack_first) 
	  BAILOUT("stack underflow on ln!\n");
	for (i=0; i<length; i++) *(stack_top+i) = log(*(stack_top+i));
      }
    break;
    case OP_EXP:
      {
	int i;
	/* fprintf(stderr,"exec: OP_EXP\n"); */
	if (stack_top<stack_first) 
	  BAILOUT("stack underflow on exp!\n");
	for (i=0; i<length; i++) *(stack_top+i) = exp(*(stack_top+i));
      }
    break;
    case OP_SIN:
      {
	int i;
	/* fprintf(stderr,"exec: OP_SIN\n"); */
	if (stack_top<stack_first) 
	  BAILOUT("stack underflow on sin!\n");
	for (i=0; i<length; i++) *(stack_top+i) = sin(*(stack_top+i));
      }
    break;
    case OP_COS:
      {
	int i;
	/* fprintf(stderr,"exec: OP_COS\n"); */
	if (stack_top<stack_first) 
	  BAILOUT("stack underflow on cos!\n");
	for (i=0; i<length; i++) *(stack_top+i) = cos(*(stack_top+i));
      }
    break;
    case OP_TAN:
      {
	int i;
	/* fprintf(stderr,"exec: OP_TAN\n"); */
	if (stack_top<stack_first) 
	  BAILOUT("stack underflow on tan!\n");
	for (i=0; i<length; i++) *(stack_top+i) = tan(*(stack_top+i));
      }
    break;
    case OP_ASIN:
      {
	int i;
	/* fprintf(stderr,"exec: OP_ASIN\n"); */
	if (stack_top<stack_first) 
	  BAILOUT("stack underflow on asin!\n");
	for (i=0; i<length; i++) *(stack_top+i) = asin(*(stack_top+i));
      }
    break;
    case OP_ACOS:
      {
	int i;
	/* fprintf(stderr,"exec: OP_ACOS\n"); */
	if (stack_top<stack_first) 
	  BAILOUT("stack underflow on acos!\n");
	for (i=0; i<length; i++) *(stack_top+i) = acos(*(stack_top+i));
      }
    break;
    case OP_ATAN:
      {
	int i;
	/* fprintf(stderr,"exec: OP_ATAN\n"); */
	if (stack_top<stack_first) 
	  BAILOUT("stack underflow on atan!\n");
	for (i=0; i<length; i++) *(stack_top+i) = atan(*(stack_top+i));
      }
    break;
    case OP_ABS:
      {
	int i;
	/* fprintf(stderr,"exec: OP_ABS\n"); */
	if (stack_top<stack_first) 
	  BAILOUT("stack underflow on abs!\n");
	for (i=0; i<length; i++) *(stack_top+i) = fabs(*(stack_top+i));
      }
    break;
    case OP_DUP:
      {
	int i;
	/* fprintf(stderr,"exec: OP_DUP\n"); */
	if (stack_top>=stack_last) 
	  BAILOUT("stack overflow on dup!\n");
	for (i=0; i<length; i++) 
	  *(stack_top+RPN_CHUNKSIZE+i) = *(stack_top+i);
	stack_top += RPN_CHUNKSIZE;
      }
    break;
    case OP_SWAP:
      {
	int i;
	/* fprintf(stderr,"exec: OP_SWAP\n"); */
	if (stack_top < stack_first+RPN_CHUNKSIZE)
	  BAILOUT("stack underflow on swap!\n");
	for (i=0; i<length; i++) {
	  double tmp= *(stack_top+i);
	  *(stack_top+i) = *(stack_top-RPN_CHUNKSIZE+i);
	  *(stack_top-RPN_CHUNKSIZE+i)= tmp;
	}
      }
    break;
    case OP_POP:
      {
	int i;
	/* fprintf(stderr,"exec: OP_POP\n"); */
	if (stack_top <= stack_first)
	  BAILOUT("stack underflow on pop!\n");
	stack_top -= RPN_CHUNKSIZE;
      }
    break;
    case OP_GT:
      {
	int i;
	/* fprintf(stderr,"exec: OP_GT\n"); */
	if (stack_top < stack_first+RPN_CHUNKSIZE)
	  BAILOUT("stack underflow on compare!\n");
	stack_top -= RPN_CHUNKSIZE;
	for (i=0; i<length; i++) {
	  if (*(stack_top+RPN_CHUNKSIZE+i) > *(stack_top+i))
	    *(stack_top+i)= 1.0;
	  else *(stack_top+i)= 0.0;
	}
      }
    break;
    case OP_LT:
      {
	int i;
	/* fprintf(stderr,"exec: OP_LT\n"); */
	if (stack_top < stack_first+RPN_CHUNKSIZE)
	  BAILOUT("stack underflow on compare!\n");
	stack_top -= RPN_CHUNKSIZE;
	for (i=0; i<length; i++) {
	  if (*(stack_top+RPN_CHUNKSIZE+i) < *(stack_top+i))
	    *(stack_top+i)= 1.0;
	  else *(stack_top+i)= 0.0;
	}
      }
    break;
    case OP_EQ:
      {
	int i;
	/* fprintf(stderr,"exec: OP_EQ\n"); */
	if (stack_top < stack_first+RPN_CHUNKSIZE)
	  BAILOUT("stack underflow on compare!\n");
	stack_top -= RPN_CHUNKSIZE;
	for (i=0; i<length; i++) {
	  if (*(stack_top+RPN_CHUNKSIZE+i) == *(stack_top+i))
	    *(stack_top+i)= 1.0;
	  else *(stack_top+i)= 0.0;
	}
      }
    break;
    case OP_GE:
      {
	int i;
	/* fprintf(stderr,"exec: OP_GE\n"); */
	if (stack_top < stack_first+RPN_CHUNKSIZE)
	  BAILOUT("stack underflow on compare!\n");
	stack_top -= RPN_CHUNKSIZE;
	for (i=0; i<length; i++) {
	  if (*(stack_top+RPN_CHUNKSIZE+i) >= *(stack_top+i))
	    *(stack_top+i)= 1.0;
	  else *(stack_top+i)= 0.0;
	}
      }
    break;
    case OP_LE:
      {
	int i;
	/* fprintf(stderr,"exec: OP_LE\n"); */
	if (stack_top < stack_first+RPN_CHUNKSIZE)
	  BAILOUT("stack underflow on compare!\n");
	stack_top -= RPN_CHUNKSIZE;
	for (i=0; i<length; i++) {
	  if (*(stack_top+RPN_CHUNKSIZE+i) <= *(stack_top+i))
	    *(stack_top+i)= 1.0;
	  else *(stack_top+i)= 0.0;
	}
      }
    break;
    case OP_NE:
      {
	int i;
	/* fprintf(stderr,"exec: OP_NE\n"); */
	if (stack_top < stack_first+RPN_CHUNKSIZE)
	  BAILOUT("stack underflow on compare!\n");	
	stack_top -= RPN_CHUNKSIZE;
	for (i=0; i<length; i++) {
	  if (*(stack_top+RPN_CHUNKSIZE+i) != *(stack_top+i))
	    *(stack_top+i)= 1.0;
	  else *(stack_top+i)= 0.0;
	}
      }
    break;
    case OP_IF_KEEP:
      {
	int i;
	/* fprintf(stderr,"exec: OP_IF_KEEP\n"); */
	if (stack_top < stack_first+(2*RPN_CHUNKSIZE))
	  BAILOUT("stack underflow on if_keep!\n");
	stack_top -= (2*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  if (*(stack_top+(2*RPN_CHUNKSIZE)+i) != 0.0) {
	    *(stack_top+i)= *(stack_top+RPN_CHUNKSIZE+i);
	  }
	  /* no else needed, since right value is in place */
	}
      }
    break;
    case OP_IF_PRINT:
      {
	int i;
	int j;
	long n= runner->param.l;
	/*fprintf(stderr,"exec: OP_IF_PRINT %ld\n", n);*/
	if (stack_top < stack_first+(n*RPN_CHUNKSIZE))
	  BAILOUT1("stack underflow on if_print %ld!\n", n);
	stack_top -= ((n+1)*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  if (*(stack_top+((n+1)*RPN_CHUNKSIZE)+i) != 0.0) {
	    for (j=1; j<=n; j++) 
	      printf("%.15g   ",*(stack_top+(j*RPN_CHUNKSIZE)+i));
	    printf("\n");
	  }
	}
      }
    break;
    case OP_CONST:
      {
	int i;
	double val;
	/*fprintf(stderr,"exec: OP_CONST %f\n",runner->param.f);*/
	val= runner->param.f;
	if (stack_top>=stack_last) 
	  BAILOUT("stack overflow on load const!\n");
	stack_top += RPN_CHUNKSIZE;
	for (i=0; i<length; i++) *(stack_top+i)= val;
      }
    break;

    case OP_CT:
      {
	int i;
	int one= 1;
	int status;
	double bound;
	double val;
	double qval;
	/* fprintf(stderr,"exec: OP_CT\n"); */
	if (stack_top < stack_first+RPN_CHUNKSIZE)
	  BAILOUT("stack underflow on ct!\n");
	stack_top -= (1*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  /* We encounter normalized maps with T=0 and no counts in
	   * the unsampled regions.  This hack gives an appropriate
	   * value at these points.
	   */
	  if (*(stack_top+i) == 0.0) qval= 0.5;
	  else {
	    cdf_t(&one, &val, &qval, stack_top+i, stack_top+RPN_CHUNKSIZE+i,
		  &status, &bound);
	    if (status != 0)
	      BAILOUT2("can't get P from T= %f, df= %f\n",
		       *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i));
	  }
	  *(stack_top+i)= val;
	}
      }
    break;
    case OP_INV_CT:
      {
	int i;
	int two= 2;
	int status;
	double bound;
	double val;
	double qval;
	/* fprintf(stderr,"exec: OP_INV_CT\n"); */
	if (stack_top < stack_first+RPN_CHUNKSIZE)
	  BAILOUT("stack underflow on inv_ct!\n");
	stack_top -= (1*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  if (*(stack_top+i) != 0.0) {
	    if (*(stack_top+i) != 1.0) {
	      qval= 1.0- *(stack_top+i);
	      cdf_t(&two, stack_top+i, &qval, &val, stack_top+RPN_CHUNKSIZE+i,
		    &status, &bound);
	      if (status != 0)
		BAILOUT2("can't get T from P= %f, df= %f\n",
			 *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i));
	      *(stack_top+i)= val;
	    }
	    else {
	      *(stack_top+i)= DLAMCH("o");
	    }
	  }
	  else {
	    *(stack_top+i)= -DLAMCH("o");
	  }
	}
      }
    break;

    case OP_CPOISSON:
      {
	int i;
	int one= 1;
	int status;
	double bound;
	double val;
	double qval;
	/* fprintf(stderr,"exec: OP_CPOISSON\n"); */
	if (stack_top < stack_first+RPN_CHUNKSIZE)
	  BAILOUT("stack underflow on cpoisson!\n");
	stack_top -= (1*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  cdf_poi( &one, &val, &qval, stack_top+i, stack_top+RPN_CHUNKSIZE+i,
		&status, &bound);
	  if (status != 0)
	    BAILOUT2("can't get P from X= %f, mean= %f\n",
		  *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i));
	  *(stack_top+i)= val;
	}
      }
    break;
    case OP_INV_CPOISSON:
      {
	int i;
	int two= 2;
	int status;
	double bound;
	double val;
	double qval;
	/* fprintf(stderr,"exec: OP_INV_CPOISSON\n"); */
	if (stack_top < stack_first+RPN_CHUNKSIZE)
	  BAILOUT("stack underflow on inv_cpoisson!\n");
	stack_top -= (1*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  if (*(stack_top+i) != 0.0) {
	    if (*(stack_top+i) != 1.0) {
	      qval= 1.0- *(stack_top+i);
	      cdf_poi(&two, stack_top+i, &qval, &val, 
		      stack_top+RPN_CHUNKSIZE+i, &status, &bound);
	      if (status != 0)
		BAILOUT2("can't get X from P= %f, mean= %f\n",
			 *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i));
	      *(stack_top+i)= val;
	    }
	    else {
	      *(stack_top+i)= DLAMCH("o");
	    }
	  }
	  else {
	    *(stack_top+i)= 0.0;
	  }
	}
      }
    break;

    case OP_CF:
      {
	int i;
	int one= 1;
	int status;
	double bound;
	double val;
	double qval;
	/* fprintf(stderr,"exec: OP_CF\n"); */
	if (stack_top < stack_first+(2*RPN_CHUNKSIZE))
	  BAILOUT("stack underflow on cf!\n");
	stack_top -= (2*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  cdf_f(&one, &val, &qval, stack_top+i, stack_top+RPN_CHUNKSIZE+i,
		stack_top+(2*RPN_CHUNKSIZE)+i, &status, &bound);
	  if (status != 0)
	    BAILOUT3("can't get P from F= %f, dfn= %f, dfd= %f\n",
		  *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i),
		  *(stack_top+(2*RPN_CHUNKSIZE)+i));
	  *(stack_top+i)= val;
	}
      }
    break;
    case OP_INV_CF:
      {
	int i;
	int two= 2;
	int status;
	double bound;
	double val;
	double qval;
	/* fprintf(stderr,"exec: OP_INV_CF\n"); */
	if (stack_top < stack_first+(2*RPN_CHUNKSIZE))
	  BAILOUT("stack underflow on inv_cf!\n");
	stack_top -= (2*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  if (*(stack_top+i) != 1.0) {
	    qval= 1.0- *(stack_top+i);
	    cdf_f(&two, stack_top+i, &qval, &val, stack_top+RPN_CHUNKSIZE+i,
		  stack_top+(2*RPN_CHUNKSIZE)+i, &status, &bound);
	    if (status != 0)
	      BAILOUT3("can't get F from P= %f, dfn= %f, dfd= %f\n",
		       *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i),
		       *(stack_top+(2*RPN_CHUNKSIZE)+i));
	    *(stack_top+i)= val;
	  }
	  else *(stack_top+i)= DLAMCH("o");
	}
      }
    break;

    case OP_CCHISQR:
      {
	int i;
	int one= 1;
	int status;
	double bound;
	double val;
	double qval;
	/* fprintf(stderr,"exec: OP_CCHISQR\n"); */
	if (stack_top < stack_first+RPN_CHUNKSIZE)
	  BAILOUT("stack underflow on cchisqr!\n");
	stack_top -= (1*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  cdf_chi(&one, &val, &qval, stack_top+i, stack_top+RPN_CHUNKSIZE+i,
		&status, &bound);
	  if (status != 0)
	    BAILOUT2(" can't get P from chisqr= %f, df= %f\n",
		     *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i));
	  *(stack_top+i)= val;
	}
      }
    break;
    case OP_INV_CCHISQR:
      {
	int i;
	int two= 2;
	int status;
	double bound;
	double val;
	double qval;
	/* fprintf(stderr,"exec: OP_INV_CCHISQR\n"); */
	if (stack_top < stack_first+RPN_CHUNKSIZE)
	  BAILOUT("stack underflow on inv_cchisqr!\n");
	stack_top -= (1*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  if (*(stack_top+i)!=1.0) {
	    qval= 1.0- *(stack_top+i);
	    cdf_chi(&two, stack_top+i, &qval, &val, stack_top+RPN_CHUNKSIZE+i,
		    &status, &bound);
	    if (status != 0)
	      BAILOUT2(" can't get chisqr from P= %f, df= %f\n",
		       *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i));
	    *(stack_top+i)= val;
	  }
	  else *(stack_top+i)= DLAMCH("o");
	}
      }
    break;
    
    case OP_CBETA:
      {
	int i;
	int one= 1;
	int status;
	double bound;
	double val;
	double qval;
	double b;
	/* fprintf(stderr,"exec: OP_CBETA\n"); */
	if (stack_top < stack_first+(2*RPN_CHUNKSIZE))
	  BAILOUT("stack underflow on cbeta!\n");
	stack_top -= (2*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  b= 1.0- *(stack_top+i);
	  cdf_bet(&one, &val, &qval, stack_top+i, &b, stack_top+RPN_CHUNKSIZE+i,
		stack_top+(2*RPN_CHUNKSIZE)+i, &status, &bound);
	  if (status != 0)
	    BAILOUT3("can't get P from Beta= %f, A= %f, B= %f\n",
		     *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i),
		     *(stack_top+(2*RPN_CHUNKSIZE)+i));
	  *(stack_top+i)= val;
	}
      }
    break;
    case OP_INV_CBETA:
      {
	int i;
	int two= 2;
	int status;
	double bound;
	double val;
	double qval;
	double b;
	/* fprintf(stderr,"exec: OP_INV_CBETA\n"); */
	if (stack_top < stack_first+(2*RPN_CHUNKSIZE))
	  BAILOUT("stack underflow on inv_cbeta!\n");
	stack_top -= (2*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  qval= 1.0- *(stack_top+i);
	  cdf_bet(&two, stack_top+i, &qval, &val, &b, stack_top+RPN_CHUNKSIZE+i,
		stack_top+(2*RPN_CHUNKSIZE)+i, &status, &bound);
	  if (status != 0)
	    BAILOUT3(" can't get Beta from P= %f, A= %f, B= %f\n",
		     *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i),
		     *(stack_top+(2*RPN_CHUNKSIZE)+i));
	  *(stack_top+i)= val;
	}
      }
    break;

    case OP_CBINOM:
      {
	int i;
	int one= 1;
	int status;
	double bound;
	double val;
	double qval;
	double ompr;
	/* fprintf(stderr,"exec: OP_CBINOM\n"); */
	if (stack_top < stack_first+(2*RPN_CHUNKSIZE))
	  BAILOUT("stack underflow on cbinom!\n");
	stack_top -= (2*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  ompr= 1.0- *(stack_top+(2*RPN_CHUNKSIZE)+i);
	  cdf_bin(&one, &val, &qval, stack_top+i, stack_top+RPN_CHUNKSIZE+i,
		stack_top+(2*RPN_CHUNKSIZE)+i, &ompr, &status, &bound);
	  if (status != 0)
	    BAILOUT3(" can't get P from S= %f, trials= %f, prob= %f\n",
		     *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i),
		     *(stack_top+(2*RPN_CHUNKSIZE)+i));
	  *(stack_top+i)= val;
	}
      }
    break;
    case OP_INV_CBINOM:
      {
	int i;
	int two= 2;
	int status;
	double bound;
	double val;
	double qval;
	double ompr;
	/* fprintf(stderr,"exec: OP_INV_CBINOM\n"); */
	if (stack_top < stack_first+(2*RPN_CHUNKSIZE))
	  BAILOUT("stack underflow on inv_cbinom!\n");
	stack_top -= (2*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  if (*(stack_top+i) != 0.0) {
	    qval= 1.0- *(stack_top+i);
	    ompr= 1.0- *(stack_top+(2*RPN_CHUNKSIZE)+i);
	    cdf_bin(&two, stack_top+i, &qval, &val, stack_top+RPN_CHUNKSIZE+i,
		    stack_top+(2*RPN_CHUNKSIZE)+i, &ompr, &status, &bound);
	    if (status != 0)
	      BAILOUT3(" can't get S from P= %f, trials= %f, prob= %f\n",
		       *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i),
		       *(stack_top+(2*RPN_CHUNKSIZE)+i));
	    *(stack_top+i)= val;
	  }
	  else *(stack_top+i)= 0.0;
	}
      }
    break;

    case OP_CGAMMA:
      {
	int i;
	int one= 1;
	int status;
	double bound;
	double val;
	double qval;
	/* fprintf(stderr,"exec: OP_CGAMMA\n"); */
	if (stack_top < stack_first+(2*RPN_CHUNKSIZE))
	  BAILOUT("stack underflow on cgamma!\n");
	stack_top -= (2*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  cdf_gam(&one, &val, &qval, stack_top+i, stack_top+RPN_CHUNKSIZE+i,
		stack_top+(2*RPN_CHUNKSIZE)+i, &status, &bound);
	  if (status != 0)
	    BAILOUT3(" can't get P from X= %f, shape= %f, scale= %f\n",
		     *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i),
		     *(stack_top+(2*RPN_CHUNKSIZE)+i));
	  *(stack_top+i)= val;
	}
      }
    break;
    case OP_INV_CGAMMA:
      {
	int i;
	int two= 2;
	int status;
	double bound;
	double val;
	double qval;
	/* fprintf(stderr,"exec: OP_INV_CGAMMA\n"); */
	if (stack_top < stack_first+(2*RPN_CHUNKSIZE))
	  BAILOUT("stack underflow on inv_cgamma!\n");
	stack_top -= (2*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  qval= 1.0- *(stack_top+i);
	  cdf_gam(&two, stack_top+i, &qval, &val, stack_top+RPN_CHUNKSIZE+i,
		stack_top+(2*RPN_CHUNKSIZE)+i, &status, &bound);
	  if (status != 0)
	    BAILOUT3(" can't get X from P= %f, shape= %f, scale= %f\n",
		     *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i),
		     *(stack_top+(2*RPN_CHUNKSIZE)+i));
	  *(stack_top+i)= val;
	}
      }
    break;

    case OP_CNORMAL:
      {
	int i;
	int one= 1;
	int status;
	double bound;
	double val;
	double qval;
	/* fprintf(stderr,"exec: OP_CNORMAL\n"); */
	if (stack_top < stack_first+(2*RPN_CHUNKSIZE))
	  BAILOUT("stack underflow on cnormal!\n");
	stack_top -= (2*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  cdf_nor(&one, &val, &qval, stack_top+i, stack_top+RPN_CHUNKSIZE+i,
		stack_top+(2*RPN_CHUNKSIZE)+i, &status, &bound);
	  if (status != 0)
	    BAILOUT3(" can't get P from X= %f, mean= %f, stdv= %f\n",
		     *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i),
		     *(stack_top+(2*RPN_CHUNKSIZE)+i));
	  *(stack_top+i)= val;
	}
      }
    break;
    case OP_INV_CNORMAL:
      {
	int i;
	int two= 2;
	int status;
	double bound;
	double val;
	double qval;
	/* fprintf(stderr,"exec: OP_INV_CNORMAL\n"); */
	if (stack_top < stack_first+(2*RPN_CHUNKSIZE))
	  BAILOUT("stack underflow on inv_cnormal!\n");
	stack_top -= (2*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  if (*(stack_top+i)!=0.0) {
	    if (*(stack_top+i)!=1.0) {
	      qval= 1.0- *(stack_top+i);
	      cdf_nor(&two, stack_top+i, &qval, &val, stack_top+RPN_CHUNKSIZE+i,
		      stack_top+(2*RPN_CHUNKSIZE)+i, &status, &bound);
	      if (status != 0)
		BAILOUT3(" can't get X from P= %f, mean= %f, stdv= %f\n",
			 *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i),
			 *(stack_top+(2*RPN_CHUNKSIZE)+i));
	      *(stack_top+i)= val;
	    }
	    else {
	      *(stack_top+i)= DLAMCH("o");
	    }
	  }
	  else {
	    *(stack_top+i)= -DLAMCH("o");
	  }
	}
      }
    break;

    case OP_FCT:
      {
	int i;
	int one= 1;
	int status;
	double bound;
	double val;
	double qval;
	/* fprintf(stderr,"exec: OP_FCT\n"); */
	if (stack_top < stack_first+RPN_CHUNKSIZE)
	  BAILOUT("stack underflow on ct!\n");
	stack_top -= (1*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  /* We encounter normalized maps with T=0 and no counts in
	   * the unsampled regions.  This hack gives an appropriate
	   * value at these points.
	   */
	  if (*(stack_top+i) == 0.0) val= qval= 0.5;
	  else {
	    cdf_t(&one, &val, &qval, stack_top+i, stack_top+RPN_CHUNKSIZE+i,
		  &status, &bound);
	    if (status != 0)
	      BAILOUT2("can't get P from T= %f, df= %f\n",
		       *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i));
	  }
	  *(stack_top+i)= (qval>=val)?val:-qval;
	}
      }
    break;
    case OP_INV_FCT:
      {
	int i;
	int two= 2;
	int status;
	double bound;
	double val;
	double qval;
	double x;
	/* fprintf(stderr,"exec: OP_INV_FCT\n"); */
	if (stack_top < stack_first+RPN_CHUNKSIZE)
	  BAILOUT("stack underflow on inv_ct!\n");
	stack_top -= (1*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  if (*(stack_top+i)!=0.0) {
	    if (*(stack_top+i)>=0.0) {
	      val= *(stack_top+i);
	      qval= 1.0-val;
	    }
	    else {
	      qval= -(*(stack_top+i));
	      val= 1.0-qval;
	    }
	    cdf_t(&two, &val, &qval, &x, 
		  stack_top+RPN_CHUNKSIZE+i, &status, &bound);
	    if (status != 0)
	      BAILOUT2("can't get T from P= %f, df= %f\n",
		       *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i));
	    *(stack_top+i)= x;
	  }
	  else {
	    /* These values *should* be NInf and Inf, but... */
#if ( defined(SGI5) || defined(SGI6) || defined(SGI64) || defined(SGIMP64) )
	    if (_signbit(*(stack_top+i))) {
#else
	    if (signbit(*(stack_top+i))) {
#endif
	      *(stack_top+i)= DLAMCH("o");
	    }
	    else {
	      *(stack_top+i)= -DLAMCH("o");
	    }
	  }
	}
      }
    break;

    case OP_FCPOISSON:
      {
	int i;
	int one= 1;
	int status;
	double bound;
	double val;
	double qval;
	/* fprintf(stderr,"exec: OP_FCPOISSON\n"); */
	if (stack_top < stack_first+RPN_CHUNKSIZE)
	  BAILOUT("stack underflow on cpoisson!\n");
	stack_top -= (1*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  cdf_poi( &one, &val, &qval, stack_top+i, stack_top+RPN_CHUNKSIZE+i,
		&status, &bound);
	  if (status != 0)
	    BAILOUT2("can't get P from X= %f, mean= %f\n",
		  *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i));
	  *(stack_top+i)= (qval>=val)?val:-qval;
	}
      }
    break;
    case OP_INV_FCPOISSON:
      {
	int i;
	int two= 2;
	int status;
	double bound;
	double val;
	double qval;
	double x;
	/* fprintf(stderr,"exec: OP_INV_FCPOISSON\n"); */
	if (stack_top < stack_first+RPN_CHUNKSIZE)
	  BAILOUT("stack underflow on inv_cpoisson!\n");
	stack_top -= (1*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  if (*(stack_top+i)!=0.0) {
	    if (*(stack_top+i)>=0.0) {
	      val= *(stack_top+i);
	      qval= 1.0-val;
	    }
	    else {
	      qval= -(*(stack_top+i));
	      val= 1.0-qval;
	    }
	    cdf_poi(&two, &val, &qval, &x, stack_top+RPN_CHUNKSIZE+i,
		    &status, &bound);
	    if (status != 0)
	      BAILOUT2("can't get X from P= %f, mean= %f\n",
		       *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i));
	    *(stack_top+i)= x;
	  }
	  else {
	    /* These values *should* be 0.0 and Inf, but... */
#if ( defined(SGI5) || defined(SGI6) || defined(SGI64) || defined(SGIMP64) )
	    if (_signbit(*(stack_top+i))) {
#else
	    if (signbit(*(stack_top+i))) {
#endif
	      *(stack_top+i)= DLAMCH("o");
	    }
	    else {
	      *(stack_top+i)= 0.0;
	    }
	  }
	}
      }
    break;

    case OP_FCF:
      {
	int i;
	int one= 1;
	int status;
	double bound;
	double val;
	double qval;
	/* fprintf(stderr,"exec: OP_FCF\n"); */
	if (stack_top < stack_first+(2*RPN_CHUNKSIZE))
	  BAILOUT("stack underflow on cf!\n");
	stack_top -= (2*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  cdf_f(&one, &val, &qval, stack_top+i, stack_top+RPN_CHUNKSIZE+i,
		stack_top+(2*RPN_CHUNKSIZE)+i, &status, &bound);
	  if (status != 0)
	    BAILOUT3("can't get P from F= %f, dfn= %f, dfd= %f\n",
		  *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i),
		  *(stack_top+(2*RPN_CHUNKSIZE)+i));
	  *(stack_top+i)= (qval>=val)?val:-qval;
	}
      }
    break;
    case OP_INV_FCF:
      {
	int i;
	int two= 2;
	int status;
	double bound;
	double val;
	double qval;
	double x;
	/* fprintf(stderr,"exec: OP_INV_FCF\n"); */
	if (stack_top < stack_first+(2*RPN_CHUNKSIZE))
	  BAILOUT("stack underflow on inv_cf!\n");
	stack_top -= (2*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  if (*(stack_top+i)!=0.0) {
	    if (*(stack_top+i)>=0.0) {
	      val= *(stack_top+i);
	      qval= 1.0-val;
	    }
	    else {
	      qval= -(*(stack_top+i));
	      val= 1.0-qval;
	    }
	    cdf_f(&two, &val, &qval, &x, stack_top+RPN_CHUNKSIZE+i,
		  stack_top+(2*RPN_CHUNKSIZE)+i, &status, &bound);
	    if (status != 0)
	      BAILOUT3("can't get F from P= %f, dfn= %f, dfd= %f\n",
		       *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i),
		       *(stack_top+(2*RPN_CHUNKSIZE)+i));
	    *(stack_top+i)= x;
	  }
	  else {
	    /* These values *should* be 0.0 and Inf, but... */
#if ( defined(SGI5) || defined(SGI6) || defined(SGI64) || defined(SGIMP64) )
	    if (_signbit(*(stack_top+i))) {
#else
	    if (signbit(*(stack_top+i))) {
#endif
	      *(stack_top+i)= DLAMCH("o");
	    }
	    else {
	      *(stack_top+i)= 0.0;
	    }
	  }
	}
      }
      break;

    case OP_FCCHISQR:
      {
	int i;
	int one= 1;
	int status;
	double bound;
	double val;
	double qval;
	/* fprintf(stderr,"exec: OP_FCCHISQR\n"); */
	if (stack_top < stack_first+RPN_CHUNKSIZE)
	  BAILOUT("stack underflow on cchisqr!\n");
	stack_top -= (1*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  cdf_chi(&one, &val, &qval, stack_top+i, stack_top+RPN_CHUNKSIZE+i,
		&status, &bound);
	  if (status != 0)
	    BAILOUT2(" can't get P from chisqr= %f, df= %f\n",
		     *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i));
	  *(stack_top+i)= (qval>=val)?val:-qval;
	}
      }
    break;
    case OP_INV_FCCHISQR:
      {
	int i;
	int two= 2;
	int status;
	double bound;
	double val;
	double qval;
	double x;
	/* fprintf(stderr,"exec: OP_INV_FCCHISQR\n"); */
	if (stack_top < stack_first+RPN_CHUNKSIZE)
	  BAILOUT("stack underflow on inv_cchisqr!\n");
	stack_top -= (1*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  if (*(stack_top+i)!=0.0) {
	    if (*(stack_top+i)>=0.0) {
	      val= *(stack_top+i);
	      qval= 1.0-val;
	    }
	    else {
	      qval= -(*(stack_top+i));
	      val= 1.0-qval;
	    }
	    cdf_chi(&two, &val, &qval, &x, stack_top+RPN_CHUNKSIZE+i,
		    &status, &bound);
	    if (status != 0)
	      BAILOUT2(" can't get chisqr from P= %f, df= %f\n",
		       *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i));
	    *(stack_top+i)= x;
	  }
	  else {
	    /* These values *should* be 0.0 and Inf, but... */
#if ( defined(SGI5) || defined(SGI6) || defined(SGI64) || defined(SGIMP64) )
	    if (_signbit(*(stack_top+i))) {
#else
	    if (signbit(*(stack_top+i))) {
#endif
	      *(stack_top+i)= DLAMCH("o");
	    }
	    else {
	      *(stack_top+i)= 0.0;
	    }
	  }
	}
      }
    break;
    
    case OP_FCBETA:
      {
	int i;
	int one= 1;
	int status;
	double bound;
	double val;
	double qval;
	double x;
	double y;
	/* fprintf(stderr,"exec: OP_FCBETA\n"); */
	if (stack_top < stack_first+(2*RPN_CHUNKSIZE))
	  BAILOUT("stack underflow on cbeta!\n");
	stack_top -= (2*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  if (*(stack_top+i)>=0.0) {
	    x= *(stack_top+i);
	    y= 1.0-x;
	  }
	  else {
	    y= -(*(stack_top+i));
	    x= 1.0-y;
	  }
	  if (x==0.0) *(stack_top+i)= 0.0;
	  else if (y==0.0) *(stack_top+i)= -0.0;
	  else {
	    cdf_bet(&one, &val, &qval, &x, &y, stack_top+RPN_CHUNKSIZE+i,
		    stack_top+(2*RPN_CHUNKSIZE)+i, &status, &bound);
	    if (status != 0)
	      BAILOUT3("can't get P from Beta= %f, A= %f, B= %f\n",
		       *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i),
		       *(stack_top+(2*RPN_CHUNKSIZE)+i));
	    *(stack_top+i)= (qval>=val)?val:-qval;
	  }
	}
      }
    break;
    case OP_INV_FCBETA:
      {
	int i;
	int two= 2;
	int status;
	double bound;
	double val;
	double qval;
	double x;
	double y;
	/* fprintf(stderr,"exec: OP_INV_FCBETA\n"); */
	if (stack_top < stack_first+(2*RPN_CHUNKSIZE))
	  BAILOUT("stack underflow on inv_cbeta!\n");
	stack_top -= (2*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  if (*(stack_top+i)!=0.0) {
	    if (*(stack_top+i)>=0.0) {
	      val= *(stack_top+i);
	      qval= 1.0-val;
	    }
	    else {
	      qval= -(*(stack_top+i));
	      val= 1.0-qval;
	    }
	    cdf_bet(&two, &val, &qval, &x, &y, stack_top+RPN_CHUNKSIZE+i,
		    stack_top+(2*RPN_CHUNKSIZE)+i, &status, &bound);
	    if (status != 0)
	      BAILOUT3(" can't get Beta from P= %f, A= %f, B= %f\n",
		       *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i),
		       *(stack_top+(2*RPN_CHUNKSIZE)+i));
	    *(stack_top+i)= (y>=x)?x:-y;
	  }
	  else {
	    /* These values *should* be 0.0 and 1.0, but... */
#if ( defined(SGI5) || defined(SGI6) || defined(SGI64) || defined(SGIMP64) )
	    if (_signbit(*(stack_top+i))) {
#else
	    if (signbit(*(stack_top+i))) {
#endif
	      *(stack_top+i)= 1.0;
	    }
	    else {
	      *(stack_top+i)= 0.0;
	    }
	  }
	}
      }
    break;

    case OP_FCBINOM:
      {
	int i;
	int one= 1;
	int status;
	double bound;
	double val;
	double qval;
	double ompr;
	/* fprintf(stderr,"exec: OP_FCBINOM\n"); */
	if (stack_top < stack_first+(2*RPN_CHUNKSIZE))
	  BAILOUT("stack underflow on cbinom!\n");
	stack_top -= (2*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  ompr= 1.0- *(stack_top+(2*RPN_CHUNKSIZE)+i);
	  cdf_bin(&one, &val, &qval, stack_top+i, stack_top+RPN_CHUNKSIZE+i,
		stack_top+(2*RPN_CHUNKSIZE)+i, &ompr, &status, &bound);
	  if (status != 0) {
	    BAILOUT3(" can't get P from S= %f, trials= %f, prob= %f\n",
		     *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i),
		     *(stack_top+(2*RPN_CHUNKSIZE)+i));
	  }
	  *(stack_top+i)= (qval>=val)?val:-qval;
	}
      }
    break;
    case OP_INV_FCBINOM:
      {
	int i;
	int two= 2;
	int status;
	double bound;
	double val;
	double qval;
	double x;
	double ompr;
	/* fprintf(stderr,"exec: OP_INV_FCBINOM\n"); */
	if (stack_top < stack_first+(2*RPN_CHUNKSIZE))
	  BAILOUT("stack underflow on inv_cbinom!\n");
	stack_top -= (2*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  if (*(stack_top+i)!=0.0) {
	    if (*(stack_top+i)>=0.0) {
	      val= *(stack_top+i);
	      qval= 1.0-val;
	    }
	    else {
	      qval= -(*(stack_top+i));
	      val= 1.0-qval;
	    }
	    ompr= 1.0- *(stack_top+(2*RPN_CHUNKSIZE)+i);
	    cdf_bin(&two, &val, &qval, &x, stack_top+RPN_CHUNKSIZE+i,
		    stack_top+(2*RPN_CHUNKSIZE)+i, &ompr, &status, &bound);
	    if (status != 0)
	      BAILOUT3(" can't get S from P= %f, trials= %f, prob= %f\n",
		       *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i),
		       *(stack_top+(2*RPN_CHUNKSIZE)+i));
	    *(stack_top+i)= x;
	  }
	  else {
#if ( defined(SGI5) || defined(SGI6) || defined(SGI64) || defined(SGIMP64) )
	    if (_signbit(*(stack_top+i))) {
#else
	    if (signbit(*(stack_top+i))) {
#endif
	      *(stack_top+i)= *(stack_top+RPN_CHUNKSIZE+i);
	    }
	    else {
	      *(stack_top+i)= 0.0;
	    }
	  }
	}
      }
    break;

    case OP_FCGAMMA:
      {
	int i;
	int one= 1;
	int status;
	double bound;
	double val;
	double qval;
	/* fprintf(stderr,"exec: OP_FCGAMMA\n"); */
	if (stack_top < stack_first+(2*RPN_CHUNKSIZE))
	  BAILOUT("stack underflow on cgamma!\n");
	stack_top -= (2*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  cdf_gam(&one, &val, &qval, stack_top+i, stack_top+RPN_CHUNKSIZE+i,
		stack_top+(2*RPN_CHUNKSIZE)+i, &status, &bound);
	  if (status != 0)
	    BAILOUT3(" can't get P from X= %f, shape= %f, scale= %f\n",
		     *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i),
		     *(stack_top+(2*RPN_CHUNKSIZE)+i));
	  *(stack_top+i)= (qval>=val)?val:-qval;
	}
      }
    break;
    case OP_INV_FCGAMMA:
      {
	int i;
	int two= 2;
	int status;
	double bound;
	double val;
	double qval;
	double x;
	/* fprintf(stderr,"exec: OP_INV_FCGAMMA\n"); */
	if (stack_top < stack_first+(2*RPN_CHUNKSIZE))
	  BAILOUT("stack underflow on inv_cgamma!\n");
	stack_top -= (2*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  if (*(stack_top+i)!=0.0) {
	    if (*(stack_top+i)>=0.0) {
	      val= *(stack_top+i);
	      qval= 1.0-val;
	    }
	    else {
	      qval= -(*(stack_top+i));
	      val= 1.0-qval;
	    }
	    cdf_gam(&two, &val, &qval, &x, stack_top+RPN_CHUNKSIZE+i,
		    stack_top+(2*RPN_CHUNKSIZE)+i, &status, &bound);
	    if (status != 0)
	      BAILOUT3(" can't get X from P= %f, shape= %f, scale= %f\n",
		       *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i),
		       *(stack_top+(2*RPN_CHUNKSIZE)+i));
	    *(stack_top+i)= x;
	  }
	  else {
#if ( defined(SGI5) || defined(SGI6) || defined(SGI64) || defined(SGIMP64) )
	    if (_signbit(*(stack_top+i))) {
#else
	    if (signbit(*(stack_top+i))) {
#endif
	      *(stack_top+i)= DLAMCH("o");
	    }
	    else {
	      *(stack_top+i)= 0.0;
	    }
	  }
	}
      }
    break;

    case OP_FCNORMAL:
      {
	int i;
	int one= 1;
	int status;
	double bound;
	double val;
	double qval;
	/* fprintf(stderr,"exec: OP_FCNORMAL\n"); */
	if (stack_top < stack_first+(2*RPN_CHUNKSIZE))
	  BAILOUT("stack underflow on cnormal!\n");
	stack_top -= (2*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  cdf_nor(&one, &val, &qval, stack_top+i, stack_top+RPN_CHUNKSIZE+i,
		stack_top+(2*RPN_CHUNKSIZE)+i, &status, &bound);
	  if (status != 0)
	    BAILOUT3(" can't get P from X= %f, mean= %f, stdv= %f\n",
		     *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i),
		     *(stack_top+(2*RPN_CHUNKSIZE)+i));
	  *(stack_top+i)= (qval>=val)?val:-qval;
	}
      }
    break;
    case OP_INV_FCNORMAL:
      {
	int i;
	int two= 2;
	int status;
	double bound;
	double val;
	double qval;
	double x;
	/* fprintf(stderr,"exec: OP_INV_FCNORMAL\n"); */
	if (stack_top < stack_first+(2*RPN_CHUNKSIZE))
	  BAILOUT("stack underflow on inv_cnormal!\n");
	stack_top -= (2*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  if (*(stack_top+i) != 0.0) {
	    if (*(stack_top+i)>0.0) {
	      val= *(stack_top+i);
	      qval= 1.0-val;
	    }
	    else {
	      qval= -(*(stack_top+i));
	      val= 1.0-qval;
	    }
	    cdf_nor(&two, &val, &qval, &x, stack_top+RPN_CHUNKSIZE+i,
		    stack_top+(2*RPN_CHUNKSIZE)+i, &status, &bound);
	    if (status != 0)
	      BAILOUT3(" can't get X from P= %f, mean= %f, stdv= %f\n",
		       *(stack_top+i),*(stack_top+RPN_CHUNKSIZE+i),
		       *(stack_top+(2*RPN_CHUNKSIZE)+i));
	    *(stack_top+i)= x;
	  }
	  else {
	    /* These values *should* be NInf and Inf, but... */
#if ( defined(SGI5) || defined(SGI6) || defined(SGI64) || defined(SGIMP64) )
	    if (_signbit(*(stack_top+i))) {
#else
	    if (signbit(*(stack_top+i))) {
#endif
	      *(stack_top+i)= DLAMCH("o");
	    }
	    else {
	      *(stack_top+i)= -DLAMCH("o");
	    }
	  }
	}
      }
    break;

    case OP_RANDOM:
      {
	int i;
	double val;
	/* fprintf(stderr,"exec: OP_RANDOM\n"); */
	if (stack_top>=stack_last) 
	  BAILOUT("stack overflow on load random!\n");
	stack_top += RPN_CHUNKSIZE;
	for (i=0; i<length; i++) {
	  do { val= drand48(); } while (val == 0.0);
	  *(stack_top+i)= val;
	}

      }
    break;

    case OP_MIN:
      {
	int i;
	/* fprintf(stderr,"exec: OP_MIN\n"); */
	if (stack_top < stack_first+RPN_CHUNKSIZE)
	  BAILOUT("stack underflow on min!\n");	
	stack_top -= RPN_CHUNKSIZE;
	for (i=0; i<length; i++) {
	  if (*(stack_top+RPN_CHUNKSIZE+i) < *(stack_top+i))
	    *(stack_top+i)= *(stack_top+RPN_CHUNKSIZE+i);
	}
      }
      break;

    case OP_MAX:
      {
	int i;
	/* fprintf(stderr,"exec: OP_MAX\n"); */
	if (stack_top < stack_first+RPN_CHUNKSIZE)
	  BAILOUT("stack underflow on max!\n");	
	stack_top -= RPN_CHUNKSIZE;
	for (i=0; i<length; i++) {
	  if (*(stack_top+RPN_CHUNKSIZE+i) > *(stack_top+i))
	    *(stack_top+i)= *(stack_top+RPN_CHUNKSIZE+i);
	}
      }
      break;

    case OP_ROT:
      {
	int i;
	long depth= abs(runner->param.l);
	/* fprintf(stderr,"exec: OP_ROT %ld\n",runner->param.l); */
	if (stack_top < stack_first+((depth-1)*RPN_CHUNKSIZE))
	  BAILOUT1("stack underflow on rot_%ld!\n",runner->param.l);
	if (runner->param.l>1) {
	  /* rotate forward */
	  for (i=0; i<length; i++) {
	    int j;
	    double val= *(stack_top-((depth-1)*RPN_CHUNKSIZE)+i);
	    for (j=depth-1; j>0; j--)
	      *(stack_top-(j*RPN_CHUNKSIZE)+i)=
		*(stack_top-((j-1)*RPN_CHUNKSIZE)+i);
	    *(stack_top+i)= val;
	  }
	}
	if (runner->param.l<-1) {
	  /* rotate backward */
	  for (i=0; i<length; i++) {
	    int j;
	    double val= *(stack_top+i);
	    for (j=1; j<depth; j++) 
	      *(stack_top-((j-1)*RPN_CHUNKSIZE)+i)=
		*(stack_top-(j*RPN_CHUNKSIZE)+i);
	    *(stack_top-((depth-1)*RPN_CHUNKSIZE)+i)= val;
	  }
	}
      }
      break;

    case OP_MISSING:
      {
	/* fprintf(stderr,"exec: OP_MISSING\n"); */
	if (stack_top>=stack_last) 
	  BAILOUT("stack overflow on load missing!\n");
	if (clock->missing_enable) {
	  stack_top += RPN_CHUNKSIZE;
	  for (i=0; i<length; i++) {
	    int t= (((offset + i)/clock->strides[clock->t_offset])
		    % clock->limits[clock->t_offset]);
	    int z= (((offset + i)/clock->strides[clock->z_offset])
		    % clock->limits[clock->z_offset]);
	    *(stack_top+i)= (re->getMissingCB(z,t,re->usrHook) ? 1.0 : 0.0);
	  }
	}
	else BAILOUT("reference to missing with an inappropriate chunk!");
      }
      break;

    case OP_IS_FINITE:
      {
	int i;
	/* fprintf(stderr,"exec: OP_IS_FINITE\n"); */
	if (stack_top < stack_first)
	  BAILOUT("stack underflow on is_finite!\n");
	for (i=0; i<length; i++) {
	  if (finite(*(stack_top+i)))
	    *(stack_top+i)= 1.0;
	  else *(stack_top+i)= 0.0;
	}
      }
      break;

    case OP_ROUND:
      {
	int i;
	/* fprintf(stderr,"exec: OP_ROUND\n"); */
	if (stack_top < stack_first)
	  BAILOUT("stack underflow on round!\n");
	for (i=0; i<length; i++) *(stack_top+i) = rint(*(stack_top+i));
      }
      break;

    case OP_FLOOR:
      {
	int i;
	/* fprintf(stderr,"exec: OP_FLOOR\n"); */
	if (stack_top < stack_first)
	  BAILOUT("stack underflow on floor!\n");
	for (i=0; i<length; i++) *(stack_top+i) = floor(*(stack_top+i));
      }
      break;

    case OP_CEILING:
      {
	int i;
	/* fprintf(stderr,"exec: OP_CEILING\n"); */
	if (stack_top < stack_first)
	  BAILOUT("stack underflow on ceiling!\n");
	for (i=0; i<length; i++) *(stack_top+i) = ceil(*(stack_top+i));
      }
      break;

    case OP_MODULO:
      {
	int i;
	/* fprintf(stderr,"exec: OP_MODULO\n"); */
	if (stack_top <= stack_first)
	  BAILOUT("stack underflow on modulo!\n");
	stack_top -= RPN_CHUNKSIZE;
	for (i=0; i<length; i++) 
	  *(stack_top+i) = 
	    (long)rint(*(stack_top+i)) % (long)rint(*(stack_top+RPN_CHUNKSIZE+i));
      }
    break;

    case OP_SWITCH:
      {
	int i;
	long depth= abs(runner->param.l);
	/*fprintf(stderr,"exec: OP_SWITCH %ld\n",runner->param.l);*/
	if (stack_top < stack_first+(depth*RPN_CHUNKSIZE))
	  BAILOUT1("stack underflow on switch_%ld!\n",
		   runner->param.l);
	for (i=0; i<length; i++) {
	  int which= rint( *(stack_top+i) );
	  if (which<0)
	    BAILOUT3(" invalid index %g = %ld on switch_%ld!\n",
		    *(stack_top+i),which,depth);
	  if (which>=depth) 
	    BAILOUT3("stack underflow on switch_%ld with index %g = %ld!\n",
		     depth,*(stack_top+i),which);
	  *(stack_top + i - depth*RPN_CHUNKSIZE)=
	    *(stack_top + i - (which+1)*RPN_CHUNKSIZE);
	}
	stack_top -= depth*RPN_CHUNKSIZE;
      }
      break;

    case OP_CX_MAG:
      {
	int i;
	/* fprintf(stderr,"exec: OP_CX_MAG\n"); */
	if (stack_top < stack_first+RPN_CHUNKSIZE)
	  BAILOUT("stack underflow on cx_mag!\n");
	stack_top -= RPN_CHUNKSIZE;
	for (i=0; i<length; i++) {
	  double r1= *(stack_top+i);
	  double i1= *(stack_top+i+RPN_CHUNKSIZE);
	  *(stack_top+i)= sqrt( r1*r1 + i1*i1 );
	}
      }
    break;

    case OP_CX_PHASE:
      {
	int i;
	/* fprintf(stderr,"exec: OP_CX_PHASE\n"); */
	if (stack_top < stack_first+RPN_CHUNKSIZE)
	  BAILOUT("stack underflow on cx_mag!\n");
	stack_top -= RPN_CHUNKSIZE;
	for (i=0; i<length; i++) {
	  double r1= *(stack_top+i);
	  double i1= *(stack_top+i+RPN_CHUNKSIZE);
	  *(stack_top+i)= atan2(i1,r1);
	}
      }
    break;

    case OP_CX_PLUS:
      {
	int i;
	/* fprintf(stderr,"exec: OP_CX_PLUS\n"); */
	if (stack_top < stack_first+3*RPN_CHUNKSIZE)
	  BAILOUT("stack underflow on cx_+!\n");
	stack_top -= 2*RPN_CHUNKSIZE;
	for (i=0; i<length; i++) {
	  *(stack_top+i) += *(stack_top+2*RPN_CHUNKSIZE+i);
	  *(stack_top+i-RPN_CHUNKSIZE) += *(stack_top+RPN_CHUNKSIZE+i);
	}
      }
    break;

    case OP_CX_MINUS:
      {
	int i;
	/* fprintf(stderr,"exec: OP_CX_MINUS\n"); */
	if (stack_top < stack_first+3*RPN_CHUNKSIZE)
	  BAILOUT("stack underflow on cx_-!\n");
	stack_top -= 2*RPN_CHUNKSIZE;
	for (i=0; i<length; i++) {
	  *(stack_top+i) -= *(stack_top+2*RPN_CHUNKSIZE+i);
	  *(stack_top+i-RPN_CHUNKSIZE) -= *(stack_top+RPN_CHUNKSIZE+i);
	}
      }
    break;

    case OP_CX_MULT:
      {
	int i;
	/* fprintf(stderr,"exec: OP_CX_MULT\n"); */
	if (stack_top < stack_first+3*RPN_CHUNKSIZE)
	  BAILOUT("stack underflow on cx_*!\n");
	stack_top -= 2*RPN_CHUNKSIZE;
	for (i=0; i<length; i++) {
	  double r1= *(stack_top+i-RPN_CHUNKSIZE);
	  double i1= *(stack_top+i);
	  double r2= *(stack_top+i+RPN_CHUNKSIZE);
	  double i2= *(stack_top+i+2*RPN_CHUNKSIZE);
	  *(stack_top+i-RPN_CHUNKSIZE)= r1*r2 - i1*i2;
	  *(stack_top+i)= r1*i2 + r2*i1;
	}
      }
    break;

    case OP_CX_DIV:
      {
	int i;
	/* fprintf(stderr,"exec: OP_CX_DIV\n"); */
	if (stack_top < stack_first+3*RPN_CHUNKSIZE)
	  BAILOUT("stack underflow on cx_/!\n");
	stack_top -= 2*RPN_CHUNKSIZE;
	for (i=0; i<length; i++) {
	  double r1= *(stack_top+i-RPN_CHUNKSIZE);
	  double i1= *(stack_top+i);
	  double r2= *(stack_top+i+RPN_CHUNKSIZE);
	  double i2= *(stack_top+i+2*RPN_CHUNKSIZE);
	  double denom= r2*r2 + i2*i2;
	  *(stack_top+i-RPN_CHUNKSIZE)= (r1*r2 + i1*i2)*(1.0/denom);
	  *(stack_top+i)= (i1*r2 - r1*i2)*(1.0/denom);
	}
      }
    break;

    case OP_CX_CONJ:
      {
	int i;
	/* fprintf(stderr,"exec: OP_CX_CONJ\n"); */
	if (stack_top < stack_first+RPN_CHUNKSIZE)
	  BAILOUT("stack underflow on cx_conj!\n");
	for (i=0; i<length; i++) {
	  *(stack_top+i) *= -1.0; /* flip sign of imaginary part */
	}
      }
    break;

    case OP_CX_TOPHASEREP:
      {
	int i;
	/* fprintf(stderr,"exec: OP_CX_TOPHASEREP\n"); */
	if (stack_top < stack_first+RPN_CHUNKSIZE)
	  BAILOUT("stack underflow on cx_tophaserep!\n");
	for (i=0; i<length; i++) {
	  double r1= *(stack_top+i-RPN_CHUNKSIZE);
	  double i1= *(stack_top+i);
	  *(stack_top+i-RPN_CHUNKSIZE)= sqrt( r1*r1 + i1*i1 );
	  *(stack_top+i)= atan2(i1,r1);
	}
      }
    break;

    case OP_CX_FMPHASEREP:
      {
	int i;
	/* fprintf(stderr,"exec: OP_CX_FMPHASEREP\n"); */
	if (stack_top < stack_first+RPN_CHUNKSIZE)
	  BAILOUT("stack underflow on cx_fmphaserep!\n");
	for (i=0; i<length; i++) {
	  double r= *(stack_top+i-RPN_CHUNKSIZE);
	  double phi= *(stack_top+i);
	  *(stack_top+i-RPN_CHUNKSIZE)= r*cos(phi);
	  *(stack_top+i)= r*sin(phi);
	}
      }
    break;

    case OP_CX_IF_KEEP:
      {
	int i;
	/* fprintf(stderr,"exec: OP_CX_IF_KEEP\n"); */
	if (stack_top < stack_first+(4*RPN_CHUNKSIZE))
	  BAILOUT("stack underflow on cx_if_keep!\n");
	stack_top -= (3*RPN_CHUNKSIZE);
	for (i=0; i<length; i++) {
	  if (*(stack_top+(3*RPN_CHUNKSIZE)+i) != 0.0) {
	    *(stack_top+i)= *(stack_top+2*RPN_CHUNKSIZE+i);
	    *(stack_top-RPN_CHUNKSIZE+i)= *(stack_top+RPN_CHUNKSIZE+i);
	  }
	  /* no else needed, since right value is in place */
	}
      }
    break;

    case OP_SIGNBIT:
      {
	int i;
	/* fprintf(stderr,"exec: OP_SIGNBIT\n"); */
	if (stack_top<stack_first) 
	  BAILOUT("stack underflow on sqrt!\n");
	for (i=0; i<length; i++) 
#if ( defined(SGI5) || defined(SGI6) || defined(SGI64) || defined(SGIMP64) )
	  *(stack_top+i) = (_signbit(*(stack_top+i)) ? 1:0);
#else
	  *(stack_top+i) = (signbit(*(stack_top+i)) ? 1:0);
#endif
      }
    break;

    case OP_NOOP: /* do nothing */
      /* fprintf(stderr,"exec: OP_NOOP\n"); */
      break;
    }

  }

  if (re->complexFlag) {
    if (re->outfileFlag) {
      /* need 2 stack elements for complex output */
      if (stack_top < stack_first+RPN_CHUNKSIZE)
	BAILOUT("rpnRun: complex mode stack underflow on exit!");
    }
    else {
      if (stack_top < stack_first-RPN_CHUNKSIZE)
	BAILOUT("rpnRun: complex mode stack underflow on exit!");
    }
  }
  else {
    if (re->outfileFlag) {
      if (stack_top < stack_first)
	BAILOUT("rpnRun: stack underflow on exit!");
    }
    else {
      if (stack_top < stack_first-RPN_CHUNKSIZE)
	BAILOUT("rpnRun: stack underflow on exit!");
    }
  }
  return( stack_top );
}

#undef BAILOUT
#undef BAILOUT1
#undef BAILOUT2
#undef BAILOUT3


static void destroyProgram(Program* prog)
{
  free(prog->code);
  free(prog);
}

static void destroyClock(Clock* clock)
{
  if (clock->limits) free(clock->limits);
  if (clock->strides) free(clock->strides);
  if (clock->string) free(clock->string);
}

RpnEngine* createRpnEngine(const int nInputs_in, void* usrHook_in,
			   const char* (*getDimensionsCB_in)(const int which,
							     void* hook),
			   const long (*getDimExtentCB_in)(const int which,
							   const char dim,
							   void* hook),
			   void (*inputCB_in)(const int which, const long n,
					      const long long offset,
					      double* buf, void* hook),
			   void (*inputComplexCB_in)(const int which, 
						     const long n,
						     const long long offset,
						     double* buf1, 
						     double* buf2,
						     void* hook),
			   int (*missingCB_in)(const long z, const long t,
					       void* hook))
{
  long i;
  RpnEngine* result= NULL;
  if (!(result=(RpnEngine*)malloc(sizeof(RpnEngine))))
    Abort("createRpnEngine: unable to allocate %d bytes!\n",
	  sizeof(RpnEngine));
  result->nInputs= nInputs_in;
  result->usrHook= usrHook_in;
  result->getDimensionsCB= getDimensionsCB_in;
  result->getDimExtentCB= getDimExtentCB_in;
  result->getInputCB= inputCB_in;
  result->getInputComplexCB= inputComplexCB_in;
  result->getMissingCB= missingCB_in;
  result->complexFlag= 0;
  result->outfileFlag= 0;
  result->verboseFlag= 0;
  result->debugFlag= 0;
  result->clock= NULL;
  result->program= NULL;
  result->errorString= NULL;
  setErrorStr(result, "There is no error");
  if (!(result->stackBlock=(double**)malloc(MAX_STACK*sizeof(double*))))
    Abort("createRpnEngine: unable to allocate %d bytes!\n",
	  MAX_STACK*sizeof(double*));
  if (!(result->stackBlock[0]=
	(double*)malloc(MAX_STACK*RPN_CHUNKSIZE*sizeof(double))))
    Abort("createRpnEngine: unable to allocate %d bytes!\n",
	  MAX_STACK*RPN_CHUNKSIZE*sizeof(double));
  for (i=1; i<MAX_STACK; i++) 
    result->stackBlock[i]= result->stackBlock[0]+i*RPN_CHUNKSIZE;
  return result;
}

void rpnDestroyEngine( RpnEngine* re )
{
  if (re->stackBlock) {
    if (re->stackBlock[0]) free(re->stackBlock[0]);
    free(re->stackBlock);
  }
  if (re->clock) destroyClock(re->clock);
  if (re->program) destroyProgram(re->program);
  free(re);
}

const char* rpnGetErrorString(RpnEngine* re)
{
  return re->errorString;
}

int rpnCompile( RpnEngine* re, const char* script )
{
  char* myScript= strdup(script);
  if (!(re->program= compile(re, myScript, re->clock))) {
    setErrorStr(re,"error compiling expressiong <%s>!\n",script);
    free(myScript);
    return 0;
  }
  free(myScript);
  return 1;
}

int rpnCompileFile( RpnEngine* re, const char* fname )
{
  if (!(re->program= compile_file(re, fname, re->clock))) {
    setErrorStr(re,"error compiling script file <%s>!\n",fname);
    return 0;
  }
  return 1;
}

int rpnInit( RpnEngine* re )
{
  int i;
  
  /* The structure of the first file sets the structure of the output.
   * Is there any reason we can't use this file?
   */
  if (!check_dims_acceptable(re,0))
    return 0;

  /* Initialize the clock */
  re->clock= clock_init(re,0);

  /* Scan other inputs for consistency with the clock.  This catches
   * user errors involving files permuted relative to each other.
   */
  for (i=1; i<re->nInputs; i++) {
    if (!(check_dims_acceptable(re,i))) {
      char* s= strdup(rpnGetErrorString(re));
      setErrorStr(re, "problem with input %d: %s",i+1,s);
      free(s);
      return 0;
    }
    if (!(check_dims_compatible(re,i,re->clock))) {
      char* s= strdup(rpnGetErrorString(re));
      setErrorStr(re, "mismatch between inputs 1 and %d: %s",i+1,s);
      free(s);
      return 0;
    }
  }

  return 1;
}

