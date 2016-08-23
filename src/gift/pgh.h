/*
   pgh.h - Routines for reading/writing PGH files
              in the gift program
   
   Copyright (c) 1996 Pittsburgh Supercomputing Center
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

   HISTORY:
   4/96 - written by Greg Hood (PSC)
*/

Boolean PghCheckFormat ();
void PghStartReading ();
void PghReadImage ();
void PghEndReading ();

void PghStartWriting ();
void PghWriteImage ();
void PghEndWriting ();
