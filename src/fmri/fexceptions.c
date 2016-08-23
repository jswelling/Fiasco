/************************************************************
 *                                                          *
 *  fexceptions.c                                            *
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

static char rcsid[] = "$Id: fexceptions.c,v 1.6 2006/09/28 18:45:42 welling Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <fmri.h>
#include <fexceptions.h>

/* Note that this must match the ExceptionType enumerated type */
static const char* __fex_exceptionTypeNames[]= {
  "BaseException",
  "IOException",
  "DicomException",
  "MemoryException"
};

static Exception* __fex_currentException= NULL;
static ExceptionContext* __fex_currentExceptionContext= NULL;

Exception* __fex_getCurrentException()
{
  return __fex_currentException;
}

const char* fex_getExceptionTypeName(Exception* e) {
  return __fex_exceptionTypeNames[(int)(e->type)];
}

const char* fex_getExceptionString(Exception* e)
{
  return e->str;
}

ExceptionType fex_getExceptionType(Exception* e)
{
  return (ExceptionType)e->type;
}

void __fex_destroyException(Exception* e)
{
  free(e->str);
  free(e);
}

void fex_raiseException(ExceptionType type, const char* fmt, ...)
{
  Exception* e;
  va_list args;

  if (!(e=(Exception*)malloc(sizeof(Exception))))
    Abort("Unable to allocate %d bytes while raising exception <%s>!\n",
	  sizeof(Exception), fmt);
  e->type= type;
  if (fmt!=NULL) {
    char buf[256];
    va_start(args,fmt);
    vsnprintf(buf, sizeof(buf), fmt, args);
    va_end(args);
    if (!(e->str= strdup(buf)))
      Abort("Unable to allocate %d bytes while raising exception <%s>!\n",
	    strlen(buf),buf);
  }
  __fex_currentException= e;
  if (__fex_currentExceptionContext)
    siglongjmp(__fex_currentExceptionContext->env, 1);
  else 
    Abort("Uncaught exception type %s <%s>!\n",
	  fex_getExceptionTypeName(e), e->str);
}

void __fex_windExceptionContext( ExceptionContext* eCtx )
{
  /* Keep this exception context on the global linked list */
  eCtx->parentExceptionContext= __fex_currentExceptionContext;
  eCtx->parentException= __fex_currentException;
  __fex_currentExceptionContext= eCtx;
}

void __fex_unwindExceptionContext()
{
  /* This context returned without an exception being thrown */
  ExceptionContext *eCtx= __fex_currentExceptionContext;
  __fex_currentException= eCtx->parentException;
  __fex_currentExceptionContext= eCtx->parentExceptionContext;
}

void __fex_rethrowThisException()
{
  /* Pop this exception context off the global linked list;
   * re-throw up to the next context.
   */
  ExceptionContext *eCtx= __fex_currentExceptionContext;
  Exception *e= __fex_currentException;
  __fex_currentExceptionContext= eCtx->parentExceptionContext;
  if (!__fex_currentExceptionContext)
    Abort("Uncaught exception type %s <%s>!\n",
	  fex_getExceptionTypeName(e), e->str);
  siglongjmp(__fex_currentExceptionContext->env, 1);  
}

