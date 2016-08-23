/************************************************************
 *                                                          *
 *  filetypes.h                                                *
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
/* Header file for file type identification */

#define FILE_DEFAULT 0
#define FILE_LX 1
#define FILE_WINDAQ 2
#define FILE_PGH_MRI 3

extern int check_filetype(const char* readfile);
extern const char* nameof_filetype(const int type);

