/*
 *	spiral_sl_write.c - routines for writing out SPIRAL_SL files
 *
 *	Copyright (c) 1996  Pittsburgh Supercomputing Center
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
 *		1/96	Written by Greg Hood (PSC)
 */

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include "fmri.h"
#include "bio.h"
#include "gift.h"
#include "spiral_sl.h"

static char rcsid[] = "$Id: spiral_sl_write.c,v 1.5 2003/02/07 21:25:46 welling Exp $";

Filename spiral_sl_write_basename;
Filename spiral_sl_write_name;
int spiral_sl_write_file_size;

void
SpiralSLStartWriting (Filename basename)
{
  if (fh.dim[0].n != 1)
    Abort("SPIRAL_SL format cannot handle more than 1 item per voxel\n");
  if (fh.dim[1].n != fh.dim[2].n)
    Abort("SPIRAL_SL format requires equal X and Y resolutions\n");
  if (fh.data_type != GIFT_INT16)
    Abort("SPIRAL_SL format requires 16-bit integers\n");
  strcpy(spiral_sl_write_basename, basename);
  spiral_sl_write_file_size = 2 * fh.n_items_per_image;
  if (spiral_sl_write_basename[0] != '\0' &&
      spiral_sl_write_basename[strlen(spiral_sl_write_basename)-1] == '/')
    mkdir(spiral_sl_write_basename, 0755);
}

void
SpiralSLWriteImage (int time,
		    int slice)
{
  sprintf(spiral_sl_write_name,
	  "%ssl%1d.%.3d",
	  spiral_sl_write_basename,
	  (slice - fh.dim[3].min) / fh.dim[3].stride + 1,
	  (time - fh.dim[4].min) / fh.dim[4].stride + 1);
  if ((output = fopen(spiral_sl_write_name, "w")) == NULL)
    Abort("Can't open %s for writing.\n", spiral_sl_write_name);
  WriteImage(time, slice);
  fclose(output);
  if (compress_output)
    Compress(spiral_sl_write_name, spiral_sl_write_file_size);
}

void
SpiralSLEndWriting ()
{
  /* nothing to do */
}
