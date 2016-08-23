/************************************************************
 *                                                          *
 *  filetypes.c                                                *
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
 *  Original programming by Mark Fitzgerald  5-96           *
 *  Modified to exclusively use libmri calls for output,    *
 *       Greg Hood (PSC), 9-98                              *
 *  Modified to read header files, and to use LX2           *
 *       resampling stuff, Joel Welling (PSC/Stats), 5-1999 *
 ************************************************************/
/* This routine tries to guess the type of a given file by looking 
 * at its header.
 */
/* NOTE- it makes the assumption that the location and value of the
 *       GE LX rdb header logo don't change between revisions of
 *       rdbm.h .  So far, that seems to be true.
 */

#include <sys/stat.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include "mri.h"
#if (SGI64 || SGI5 || SGIMP)
#include <bstring.h>
#endif

#include "bio.h"
#include "fmri.h"

/* epirecon header files */
#include "frozen_header_info.h"

/* Windaq header file */
#include "windaq_header_info.h"

static char rcsid[] = "$Id: filetypes.c,v 1.5 2007/03/21 23:50:20 welling Exp $";

#define FILE_DEFAULT 0
#define FILE_LX 1
#define FILE_WINDAQ 2

int check_filetype(const char* readfile) 
{
  FILE *fphead;
  char* actual_name;
  unsigned char *header;
  int logotype;
  struct stat file_info;
#define MYMAX(a,b) ((a>b) ? a : b)
  int read_header_length = MYMAX(FRZ_RDBHEAD_RDB_HDR_LOGO_OFF + 
				 FRZ_RDBHEAD_RDB_HDR_LOGO_SIZE, 
				 WINDAQ_LOGO_OFF1 + WINDAQ_LOGO_SIZE);
#undef MYMAX

  /* Does the file exist? */
  if (stat(readfile,&file_info)) {
    /* Well, nothing by that name.  Try adding ".mri". */
    if (!(actual_name=(char*)malloc(strlen(readfile)+10)))
      Abort("Unable to allocate %d bytes!\n",strlen(readfile)+10);
    strcpy(actual_name,readfile);
    strcat(actual_name,".mri");
    if (stat(actual_name,&file_info)) {
      /* Nope, that one's not there either. */
      Abort("Input file <%s> not found!\n",readfile);
    }
    else {
      /* Well, if it has a .mri extension, it must be .mri */
      logotype= FILE_PGH_MRI;
    }
  }
  else actual_name= strdup(readfile);

  if (read_header_length > file_info.st_size) 
    read_header_length= file_info.st_size;

  if (!(header= (unsigned char*)malloc(read_header_length))) 
    Abort("%s: unable to allocate %d bytes!\n",read_header_length);

  logotype= FILE_DEFAULT;
  if ((fphead = fopen(actual_name,"r"))!=NULL)
    {
      if (fseek(fphead, (long) 0, SEEK_SET)) {
	perror("Error seeking header");
	Abort("Error reading data header file <%s>.\n",actual_name);
      }
      else {
	if (fread(header, sizeof(char), read_header_length, fphead)
	    != read_header_length) {
	  perror("Error reading header");
	  Abort("Error reading data header file <%s>.\n",actual_name);
	}
	else {
	  if (read_header_length>=FRZ_RDBHEAD_RDB_HDR_LOGO_OFF+64) {
	    if (!(strncmp((char*)header+FRZ_RDBHEAD_RDB_HDR_LOGO_OFF, 
			  FRZ_RDBHEAD_RDB_VALID_LOGO, 
			  FRZ_RDBHEAD_RDB_HDR_LOGO_SIZE))) 
	      logotype = FILE_LX;
	    else if (!(strncmp((char*)header+FRZ_RDBHEAD_RDB_HDR_LOGO_OFF, 
			  FRZ_RDBHEAD_RDB_INVALID_LOGO, 
			  FRZ_RDBHEAD_RDB_HDR_LOGO_SIZE))) 
	      logotype = FILE_LX;
	  }
	  if (read_header_length>=WINDAQ_LOGO_OFF2) {
	    unsigned char Windaq_Logo1, Windaq_Logo2;
	    Windaq_Logo1 = BRdInt8(header+WINDAQ_LOGO_OFF1);
	    Windaq_Logo2 = BRdInt8(header+WINDAQ_LOGO_OFF2);
	    
	    if (((Windaq_Logo1 == 1) && (Windaq_Logo2 == 128)) || 
		((Windaq_Logo2 == 1) && (Windaq_Logo1 == 128)))  
	      logotype = FILE_WINDAQ;
	  }
	  if (read_header_length>=63) {
	    /* We need to check for initial string "!format = pgh", but
	     * with great tolerance for diversity.
	     */
	    char string[64];
	    char* here;
	    strncpy( string, (char*)header, 63 );
	    string[63]= '\0';
	    here= string;
	    for (; isspace(*here); here++); /* skip spaces */
	    if (!strncasecmp(here,"!format",strlen("!format"))) {
	      here += strlen("!format");
	      for (; isspace(*here); here++); /* skip spaces */
	      if (*here=='=') {
		here += 1;
		for (; isspace(*here); here++); /* skip spaces */
		if (!strncasecmp(here,"pgh",strlen("pgh"))) {
		  here += strlen("pgh");
		  logotype= FILE_PGH_MRI;
		}
	      }
	    }
	  }

	  if (fclose(fphead)) {
	    perror("Error closing header (ignored)");
	  }
	}
      }
    }
  else {
    perror("Error opening header");
    Abort("Error reading data header file <%s>.\n",readfile);
  }
  free(header);
  free(actual_name);

  return logotype;
} 

const char* nameof_filetype(const int type)
{
  switch (type) {
  case FILE_DEFAULT:  return "default";
  case FILE_WINDAQ:   return "windaq";
  case FILE_LX:       return "GE LX";
  case FILE_PGH_MRI:  return "Pittsburgh MRI";
  default:            return "***UNKNOWN***";
  }
}

