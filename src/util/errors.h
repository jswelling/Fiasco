/************************************************************
 *                                                          *
 *  errors.h                                                *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1998 Department of Statistics,         *
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
 ************************************************************/
/* The following defines direct error messages, etc. to 
 * the entry points in the misc library.  libcrg supplies an
 * alternate version.
 */
#define Message  Misc_Message
#define Usage    Misc_Message
#define Report   Misc_Report                     
#define Warning  Misc_Warning                     
#define Error    Misc_Error                      
#define Internal Misc_Abort
#define SysError Misc_SysError
#define Abort    Misc_Abort 

void Misc_Message (char *fmt, ...);
void Misc_Report (char *fmt, ...);
void Misc_Error (char *fmt, ...);
void Misc_Abort (char *fmt, ...);
void Misc_Warning (int severity, char *fmt, ...);
void Misc_SysError (char *file, int line, char *fmt, ... );

