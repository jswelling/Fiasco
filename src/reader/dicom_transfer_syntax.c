/************************************************************
 *                                                          *
 *  dicom_transfer_syntax.c                                          *
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
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "mri.h"
#if (SGI64 || SGI5 || SGIMP)
#include <bstring.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

#include "bio.h"
#include "fmri.h"
#include "fexceptions.h"
#include "array.h"
#include "stdcrg.h"
#include "misc.h"
#include "smartreader.h"

/* DICOM header files */
#include "dicom_parser.h"
#include "dicom_transfer_syntax.h"

static char rcsid[] = "$Id: dicom_transfer_syntax.c,v 1.3 2005/03/09 00:05:27 welling Exp $";

int dcm_transferSyntaxIsSupported(DicomParser* parser, 
				  const TransferSyntax* ts)
{

  if (!(ts->isEncapsulated) && !(ts->isBZip2ed) && !(ts->isDeflated))
    return 1;
  else return 0;
}

const TransferSyntax* dcm_guessTransferSyntax(DicomParser* parser,
					      FILE* fp, long long offset)
{
  /* This draws heavily on the example in the Medical Image Format FAQ */
  const TransferSyntax* result;
  char b[8];
  int bigendian = 0;
  int explicitvr = 0;
  bigfile_fseek(fp, offset, SEEK_SET);
  if (fread(b,1,8,fp) != 8) 
    fex_raiseException(EXCEPTION_IO,
		       "Could not read 8 bytes trying to infer DICOM syntax!");
  /* examine probable group number ... assume <= 0x00ff */
  if (b[0] < b[1]) bigendian=1;
  else if (b[0] == 0 && b[1] == 0) {
    /* blech ... group number is zero
     * no point in looking at element number
     * as it will probably be zero too (group length)
     * try the 32 bit value length of implicit vr
     */
      if (b[4] < b[7]) bigendian=1;
    }
  else {
    bigendian= 0; /* littleendian */
  }
  if (isupper(b[4]) && isupper(b[5])) explicitvr=1;
  /* else unrecognized ... assume default */
  if (bigendian)
    if (explicitvr)
      result= dcm_getTransferSyntaxByName(parser,"ExplicitVRBigEndian");
    else
      result= dcm_getTransferSyntaxByName(parser,"ImplicitVRBigEndian");
  else
    if (explicitvr)
      result= dcm_getTransferSyntaxByName(parser,"ExplicitVRLittleEndian");
    else
      result= dcm_getTransferSyntaxByName(parser,"ImplicitVRLittleEndian");

  if (parser->debug) 
    fprintf(stderr,
	    "guessTransferSyntax guessed <%s> based on [%x%x%x%x%x%x%x%x]\n",
	    result->name,b[0],b[1],b[2],b[3],b[4],b[5],b[6],b[7]);
  return result;
}

const TransferSyntax* dcm_getTransferSyntaxByUIDString(DicomParser* parser,
						       const char* uid)
{
  int i;
  for (i=0; i<sizeof(transferSyntaxTable)/sizeof(TransferSyntax); i++) {
    if (!strcmp(transferSyntaxTable[i].UID,uid))
      return &(transferSyntaxTable[i]);
  }
  fex_raiseException(EXCEPTION_DICOM,
		     "Unknown transfer syntax uid <%s>",uid);
  return NULL; /* not reached */
}

const TransferSyntax* dcm_getTransferSyntaxByName(DicomParser* parser,
						  const char* name)
{
  int i;
  for (i=0; i<sizeof(transferSyntaxTable)/sizeof(TransferSyntax); i++) {
    if (!strcmp(name,transferSyntaxTable[i].name)) {
      if (parser->debug) fprintf(stderr,"Syntax <%s> has index %d!\n",name,i);
      return &(transferSyntaxTable[i]);
    }
  }
  fex_raiseException(EXCEPTION_DICOM,
		     "unknown DICOM transfer syntax name <%s>",name);
  return NULL; /* not reached */
}


