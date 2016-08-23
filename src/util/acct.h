/*
 *	Very simple elapsed time accounting functions
 *		(to enable these, compile with -DACCT)
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
 */
 

#define PROCESSING	0
#define READOPEN	1
#define READING		2
#define READCLOSE	3
#define WRITEOPEN	4
#define WRITING		5
#define WRITECLOSE	6
#define MSGWAITING	7

void Acct (int n);
void PrintAcct ();
