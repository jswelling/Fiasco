/************************************************************
 *                                                          *
 *  tumble.c                                                *
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
 *  Author Joel Welling 2/98                                *
 ************************************************************/
/* This code provides the ability to do a simple permute on a Pgh 1.0
 * MRI file.  There must be only two points of folding, and the first
 * block must be at the beginning of both the input and output dim strings.
 * For example, converting "vxyzt" to "vtxyz" will work, since the three
 * dimension substrings "v", "xyz" and "t" remain contiguous and in order.
 * This is much more restrictive than the general case (which might
 * include "vxyzt" to "tzyxv", but is much faster and is the only
 * case currently used by Fiasco.
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"
#include "tumble.h"

static char rcsid[] = "$Id: tumble.c,v 1.3 2005/01/28 18:06:53 welling Exp $";

char* substring(char* string, int offset, int length) {
  static char buf[256];
  int i;
  for (i=0; i<length; i++) {
    if (!string[i+offset]) { buf[i++]='*'; break; }
    buf[i]= string[i+offset];
  }
  buf[i]= '\0';
  return buf;
}

static int find_blocks(char* instr, char* outstr, int* blk1_length, 
		       int* blk1_offset, int* blk1_out_offset,
		       int* blk2_length, int* blk2_offset, 
		       int* blk2_out_offset, int* blk3_length, 
		       int* blk3_offset, int* blk3_out_offset,
		       int* blk4_length, int* blk4_offset,
		       int* blk4_out_offset)
{
  /* For instr "abcdefgh" and outstr "abefcdgh", blk1 is "ab", blk2 is "cd",
   * blk3 is "ef", and blk4 is "gh".
   *
   * We count on the calling program to make sure all letters in instr
   * appear exactly once in outstr.
   */

  int inlen= strlen(instr);
  int outlen= strlen(outstr);
  int i;

  if (inlen != outlen) return 0; /* dim string lengths don't match */

  for (i=0; i<=inlen && instr[i]==outstr[i]; i++) {};
  *blk1_length= i; /* may be zero */
  if (*blk1_length==inlen) 
    return 0; /* No tumble necessary */

  *blk1_offset= *blk1_out_offset= 0;
  *blk2_offset= *blk1_length;

  for (i=inlen; i>0 && instr[i-1]==outstr[i-1]; i--) {};
  *blk4_length= inlen-i;
  *blk4_offset= *blk4_out_offset= inlen - *blk4_length;


  for (i=*blk1_length; i<*blk4_out_offset; i++) 
    if (outstr[i]==instr[*blk2_offset]) break;
  if (i<=outlen) *blk2_out_offset= i;
  else return 0; /* can't find block 2 in output */

  *blk2_length= *blk4_out_offset - *blk2_out_offset;
  for (i=0; i<*blk4_out_offset-*blk2_out_offset; i++) {
    if (instr[i+*blk2_offset] != outstr[i+*blk2_out_offset]) 
      return 0; /* block 2 is not contiguous in output */
  }

  *blk3_offset= *blk2_offset + *blk2_length;
  *blk3_out_offset= *blk1_length;
  *blk3_length= inlen - (*blk1_length + *blk2_length + *blk4_length);

  for (i=0; i<*blk3_length; i++) {
    if (instr[i+*blk3_offset] != outstr[i+*blk3_out_offset])
	return 0; /* block 3 is not contiguous in output */
  }

  return 1;
}
		       

int tumble_test(char* instr, char* outstr)
{
  int blk1_length= 0;
  int blk1_offset= 0;
  int blk1_out_offset= 0;
  int blk2_length= 0;
  int blk2_offset= 0;
  int blk2_out_offset= 0;
  int blk3_length= 0;
  int blk3_offset= 0;
  int blk3_out_offset= 0;
  int blk4_length= 0;
  int blk4_offset= 0;
  int blk4_out_offset= 0;

  if (!find_blocks(instr, outstr, 
		   &blk1_length, &blk1_offset, &blk1_out_offset,
		   &blk2_length, &blk2_offset, &blk2_out_offset,
		   &blk3_length, &blk3_offset, &blk3_out_offset,
		   &blk4_length, &blk4_offset, &blk4_out_offset)) {
    return 0;
  }

  return 1;
}

static long long count_elements(MRI_Dataset* in, char* chunk, char* instr, 
				int offset, int length)
{
  if (length==0) return 1;
  else {
    int i;
    char buf[256];
    long long result= 1;
    for (i=0; i<length; i++) {
      sprintf(buf,"%s.extent.%c", chunk, instr[offset+i]);
      if (mri_has(in,buf)) 
	result *= mri_get_int(in,buf);
      else Abort("permute: input file missing tag %s!\n",buf);
    }
    return result;
  }
}

int tumble(MRI_Dataset* in, MRI_Dataset* out, char* indims, char* outdims,
	   char* chunkname, long memlimit)
{
  int blk1_length;
  int blk1_offset;
  int blk1_out_offset;
  int blk2_length;
  int blk2_offset;
  int blk2_out_offset;
  int blk3_length;
  int blk3_offset;
  int blk3_out_offset;
  int blk4_length;
  int blk4_offset;
  int blk4_out_offset;
  long long blk2_stride;
  long long blk3_stride;
  long long blk1_bytes;
  long long blk1_elements;
  long long blk2_elements;
  long long blk3_elements;
  long long blk4_elements;
  long long tot_size;
  long long buf_size;
  char buf[256]; /* for composing tag names */
  char* obuf;
  char* ibuf= NULL;
  int blk2_loop;
  int blk2_min;
  int blk2_max;
  int blk3_loop;
  int blk3_min;
  int blk3_max;
  int blk4_loop;
  int blk4_min;
  int blk4_max;
  long long iframe;
  long long oframe;
  long long bytes_this_ichunk;
  long long bytes_this_ochunk;
  long long ioff;
  long long ooff;
  long long bytes_copied_from;
  long long bytes_copied_to;
  long long hits= 0;

  /* to avoid a barely conceivable error due to very long chunk names */
  if (strlen(chunkname)>=240) Abort("permute: chunk name too long!\n");

  if (!find_blocks(indims, outdims, 
		   &blk1_length, &blk1_offset, &blk1_out_offset,
		   &blk2_length, &blk2_offset, &blk2_out_offset,
		   &blk3_length, &blk3_offset, &blk3_out_offset,
		   &blk4_length, &blk4_offset, &blk4_out_offset)) {
    Abort("permute internal error: unexpectedly failed tumble test!\n");
    return 0;
  }

  blk1_elements= count_elements(in,chunkname,indims,blk1_offset,blk1_length);
  blk2_elements= count_elements(in,chunkname,indims,blk2_offset,blk2_length);
  blk3_elements= count_elements(in,chunkname,indims,blk3_offset,blk3_length);
  blk4_elements= count_elements(in,chunkname,indims,blk4_offset,blk4_length);
  blk1_bytes= blk1_elements * get_typesize(in, chunkname);
  blk2_stride= blk1_bytes * blk2_elements;
  blk3_stride= blk1_bytes * blk2_elements * blk3_elements;

  sprintf(buf,"%s.size",chunkname);
  if (mri_has(in,buf)) tot_size= mri_get_int(in,buf);
  else Abort("permute: input file missing tag %s\n",buf);

  buf_size= memlimit/2;
  if (buf_size>tot_size) buf_size= tot_size;
  buf_size = buf_size - (buf_size % blk1_bytes); /* avoid trailing bytes */

  if (!(obuf= (char*)malloc(buf_size))) {
    Abort("permute: unable to allocate %d byte output buffer!\n",buf_size);
  }

  oframe= 0;
  iframe= 0;
  bytes_copied_to= bytes_copied_from= 0;
  bytes_this_ichunk= bytes_this_ochunk= 
    (tot_size > buf_size) ? buf_size : tot_size;

  do {

    /* copy everything appropriate in ichunk to ochunk */
    blk4_min= oframe/blk3_stride;
    blk4_max= blk4_elements;
    blk4_loop= blk4_min;
    do {
      blk3_min= 0;
      blk3_max= blk3_elements;
      for (blk3_loop=blk3_min; blk3_loop<blk3_max; blk3_loop++) {
	blk2_min= 0;
	blk2_max= blk2_elements;
	for (blk2_loop= blk2_min; blk2_loop<blk2_max; blk2_loop++) {
	  
	  ioff= blk2_loop*blk1_bytes + blk3_loop*blk2_stride
	    + blk4_loop*blk3_stride;
	  ooff= blk3_loop*blk1_bytes + blk2_loop*blk3_elements*blk1_bytes
	    + blk4_loop*blk3_stride;

	  if (ooff-oframe < 0) continue;

	  if ((!ibuf) || ioff-iframe<0 || ioff-iframe>=bytes_this_ichunk) {
	    /* read the frame in which ioff resides */
	    iframe= ioff - (ioff % buf_size);
	    bytes_this_ichunk= 
	      (tot_size-iframe > buf_size) ? buf_size : (tot_size-iframe);
	    hits= 0;
	    if (ibuf) mri_discard_buffer(in,ibuf); /* free old memory */
	    ibuf= mri_get_chunk(in, chunkname, bytes_this_ichunk, 
				iframe, MRI_RAW);
	    bytes_copied_from= 0;
	  }
	  
	  if ((ooff-oframe >= 0)
	      && (ooff-oframe < bytes_this_ochunk)) {
	    /* Accelerate a couple of common special cases */
	    if (blk1_bytes==sizeof(long)) {
	      *((long*)(obuf + (ooff-oframe)))=
		*((long*)(ibuf  + (ioff-iframe)));
	    }
	    else if (blk1_bytes==sizeof(short)) {
	      *((short*)(obuf + (ooff-oframe)))=
		*((short*)(ibuf  + (ioff-iframe)));
	    }
	    else
	      memcpy( obuf + (ooff-oframe), ibuf  + (ioff-iframe), 
		      blk1_bytes );
	    bytes_copied_from += blk1_bytes;
	    bytes_copied_to += blk1_bytes;
	    hits++;
	    if (bytes_copied_to==bytes_this_ochunk) {
	      mri_set_chunk(out, chunkname, 
			    bytes_this_ochunk, oframe, MRI_RAW, obuf);
	      oframe += bytes_this_ochunk;
	      bytes_copied_to= 0;
	      bytes_this_ochunk= ((tot_size-oframe> buf_size) ? 
				  buf_size : tot_size-oframe);
	      if (buf_size<tot_size)
		Message("finished %d%%\n",
			(int)(100.0*((float)oframe)/(float)tot_size));
	    }
	  }
	  if (ooff-oframe >= bytes_this_ochunk) break;
	}
      }
      blk4_loop += 1;
    } while (blk4_loop*blk3_stride < oframe+bytes_this_ochunk);

  } while (oframe<tot_size);

  free(obuf);
  
  return 0;
}
