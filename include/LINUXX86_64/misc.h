/*
 *	Miscellaneous utility functions
 *
 *	Copyright (c) 1996 Pittsburgh Supercomputing Center
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
 *
 *	HISTORY
 *		1/96 Written by Greg Hood (PSC)
 *		8/96 Redefined all names since they were conflicting
 *			with other libraries; if you want to use another
 *			library's definition for a symbol after including
 *			misc.h, then simply #undef that symbol immediately
 *			after the #include "misc.h";  Greg Hood (PSC)
 */

#ifndef INCL_MISC_H
#define INCL_MISC_H 1

#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE 1
#endif

#define MAX_HOSTNAME_LENGTH		128
#define MAX_FILENAME_LENGTH		256
#define Boolean				Misc_Boolean
#define Hostname                   	Misc_Hostname                   
#define Filename                   	Misc_Filename                   
#define verbose                    	misc_verbose                    
#define StringCopy                 	Misc_StringCopy                 
#define GetEnvInt                  	Misc_GetEnvInt                  
#define GetEnvFloat                	Misc_GetEnvFloat                
#define ConstructFilename          	Misc_ConstructFilename          
#define Round                      	Misc_Round

typedef int Boolean;

typedef char Hostname[MAX_HOSTNAME_LENGTH];
typedef char Filename[MAX_FILENAME_LENGTH];

extern Boolean verbose;

Boolean StringCopy (char *to, const char *from, const int maxChars);

int GetEnvInt (const char *var_name);
float GetEnvFloat (const char *var_name);

void ConstructFilename (Filename out, const char *dir, const char *in);

int Round (float f);

#endif /* ifndef INCL_MISC_H */

