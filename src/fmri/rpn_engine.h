/************************************************************
 *                                                          *
 *  rpn_engine.h                                          *
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

/* $Id: rpn_engine.h,v 1.4 2004/12/08 19:38:03 welling Exp $ */

#define RPN_CHUNKSIZE 4096

typedef enum { OP_LOAD, OP_PLUS, OP_MINUS, OP_MULT, OP_DIV, OP_SQRT,
	       OP_DUP, OP_SWAP, OP_POP, OP_CONST, OP_LN, OP_EXP, OP_ABS,
	       OP_NOOP, OP_GT, OP_LT, OP_EQ, OP_GE, OP_LE, OP_NE, OP_AND,
               OP_OR, OP_IF_KEEP, OP_CBETA, OP_INV_CBETA, OP_CBINOM, 
	       OP_INV_CBINOM, OP_CCHISQR, OP_INV_CCHISQR, OP_CF, OP_INV_CF,
               OP_CT, OP_INV_CT, OP_CGAMMA, OP_INV_CGAMMA, OP_CNORMAL, 
	       OP_INV_CNORMAL, OP_CPOISSON, OP_INV_CPOISSON,
               OP_CLOCK, OP_DIM, OP_IF_PRINT, OP_RANDOM, OP_SIN,
               OP_COS, OP_TAN, OP_ASIN, OP_ACOS, OP_ATAN, OP_MIN, 
	       OP_MAX, OP_ROT, OP_MISSING, OP_IS_FINITE, OP_ROUND,
               OP_FLOOR, OP_CEILING, OP_MODULO, OP_SWITCH,
               OP_CX_MAG, OP_CX_PHASE, OP_CX_PLUS, OP_CX_MINUS, OP_CX_MULT, 
	       OP_CX_DIV, OP_CX_CONJ, OP_CX_TOPHASEREP, OP_CX_FMPHASEREP,
               OP_CX_IF_KEEP, OP_FCBETA, OP_INV_FCBETA, OP_FCBINOM, 
	       OP_INV_FCBINOM, OP_FCCHISQR, OP_INV_FCCHISQR, OP_FCF, 
	       OP_INV_FCF, OP_FCT, OP_INV_FCT, OP_FCGAMMA, OP_INV_FCGAMMA, 
	       OP_FCNORMAL, OP_INV_FCNORMAL, OP_FCPOISSON, OP_INV_FCPOISSON,
	       OP_SIGNBIT,
} Op;

typedef struct instruction_struct {
  Op op;
  union {
    long l;
    double f;
    } param;
} Instruction;

typedef struct program_struct {
  Instruction* code;
  long code_length;
} Program;

typedef struct parsepair_struct {
  char* string;
  Op op;
} ParsePair;

typedef struct macropair_struct {
  char* string;
  char* code;
} MacroPair;

typedef struct clock_struct {
  int* limits;
  long long* strides;
  char* string;
  int t_offset;
  int z_offset;
  int missing_enable;
} Clock;

typedef struct rpn_engine_struct {
  double** stackBlock;
  Clock* clock;
  Program* program;
  int nInputs;
  void* usrHook;
  const char* (*getDimensionsCB)(const int which, void* usrHook);
  const long (*getDimExtentCB)(const int which, const char dim, void* usrHook);
  void (*getInputCB)(const int which, const long n, 
		     const long long offset, double* buf, void* usrHook );
  void (*getInputComplexCB)(const int which, const long n, 
			    const long long offset, double* buf1, double* buf2,
			    void* usrHook);
  int (*getMissingCB)( const long z, const long t, void* usrHook );
  int complexFlag;  /* read complex from input files */
  int outfileFlag; /* will we be producing output, or just printing? */
  int verboseFlag;
  int debugFlag;
  char* errorString;
} RpnEngine;
	
RpnEngine* createRpnEngine(const int nInputs_in, void* usrHook_in,
			   const char* (*getDimensionsCB_in)(const int which,
							     void* hook),
			   const long (*getDimExtentCB_in)(const int which,
							   const char dim,
							   void* hook),
			   void (*inputCB_in)(const int which, const long n,
					      const long long offset,
					      double* buf,void* hook),
			   void (*inputComplexCB_in)(const int which, 
						     const long n,
						     const long long offset,
						     double* buf1, 
						     double* buf2, void* hook),
			   int (*missingCB_in)(const long z, const long t,
					       void* hook));
void rpnDestroyEngine( RpnEngine* re );
int rpnInit( RpnEngine* re);
int rpnCompile( RpnEngine* re, const char* script );
int rpnCompileFile( RpnEngine* re, const char* fname );
double* rpnRun( RpnEngine* re, long n, long long offset );
const char* rpnGetErrorString(RpnEngine* re);
void rpnSetOutputFlag(RpnEngine* re, int flag);
void rpnSetVerbose(RpnEngine* re, int flag);
void rpnSetDebug(RpnEngine* re, int flag);
void rpnSetComplex(RpnEngine* re, int flag);
