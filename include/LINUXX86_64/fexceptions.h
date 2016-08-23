/************************************************************
 *                                                          *
 *  fexceptions.h                                           *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 2004 Department of Statistics,         *
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
 *  Original programming by Joel Welling 10/04              *
 ************************************************************/

/* "$Id: fexceptions.h,v 1.6 2007/04/12 23:49:23 welling Exp $" */

/* Since the macros do setjmp, we must include it here. */
#include <setjmp.h>

/* Note that this must match __fex_exceptionTypeNames[] */
typedef enum { 
  EXCEPTION_BASE, EXCEPTION_IO, EXCEPTION_DICOM, EXCEPTION_MEM
} ExceptionType;

typedef struct exception_struct {
  char* str;
  int type;
} Exception;

typedef struct exception_context_struct {
  sigjmp_buf env;
  struct exception_context_struct* parentExceptionContext;
  Exception* parentException;
} ExceptionContext;

const char* fex_getExceptionTypeName(Exception* e);
const char* fex_getExceptionString(Exception* e);
ExceptionType fex_getExceptionType(Exception* e);
void fex_raiseException(ExceptionType type, const char* fmt, ...);

/* Call these only through the macros */
Exception* __fex_getCurrentException(void);
void __fex_destroyException(Exception* e);
void __fex_windExceptionContext( ExceptionContext* eCtx );
void __fex_unwindExceptionContext(void);
void __fex_rethrowThisException(void);

/***************************************************************
 *
 * An example of the use of these macros would be:
 *
 *
 *   void func(int i)
 *   {
 *     ....
 *     if (somethingWentWrong)
 *        fex_raiseException(EXCEPTION_IO,"Something went wrong!");
 *     ....
 *   }
 *   
 *   int main(int argc, char* argv[])
 *   {
 *     FEX_TRY
 *       ({
 *       ...
 *       func();
 *       ...
 *       })
 *       FEX_CATCH(EXCEPTION_BASE, e,
 *       {
 *         printf("I just got exception %s <%s>!\n",
 *   	     fex_getExceptionTypeName(e), fex_getExceptionString(e));
 *       });
 *       FEX_CATCH(EXCEPTION_IO, e,
 *       {
 *         printf("I just got IO exception %s <%s>!\n",
 *   	     fex_getExceptionTypeName(e), fex_getExceptionString(e));
 *       });
 *     FEX_END_TRY;
 *   
 *     fprintf(stderr,"All is now well!\n");
 *   }
 *
 **************************************************************/

#define FEX_TRY(code) \
    ExceptionContext eCtx; \
    __fex_windExceptionContext(&eCtx); \
    if (!sigsetjmp(eCtx.env,1)) { \
      code; \
      __fex_unwindExceptionContext(); } \
    else { \
      if (!__fex_getCurrentException()) \
         Abort("I can't find the exception that was just raised!\n"); \
      switch (__fex_getCurrentException()->type) {

#define FEX_CATCH(etype, evar, code) \
    case etype: { \
       Exception* evar= __fex_getCurrentException();\
       __fex_unwindExceptionContext(); \
       code; \
       __fex_destroyException(evar); \
    }; break;

#define FEX_END_TRY \
    default: { __fex_rethrowThisException(); } } }

