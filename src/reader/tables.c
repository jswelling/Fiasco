/************************************************************
 *                                                          *
 *  tables.c                                             *
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
 *  Original programming by Mark Fitzgerald  5-96           *
 *  Modified to exclusively use libmri calls for output,    *
 *       Greg Hood (PSC), 9-98                              *
 *  Modified to read header files, and to use LX2           *
 *       resampling stuff, Joel Welling (PSC/Stats), 5-1999 *
 *  Major restructuring, and the addition of formidable     *
 *       flexibility and intelligence, Joel Welling 5-2002  *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include "mri.h"
#if (SGI64 || SGI5 || SGIMP)
#include <bstring.h>
#endif
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

#include "bio.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"
#include "smartreader.h"

static char rcsid[] = "$Id: tables.c,v 1.8 2007/07/07 15:45:28 welling Exp $";

/* The following are tag definitions. */
static char* nameDefTable[][2]= {
  {"dv", "values per sample"},
  {"ds", "number of shots per slice"},
  {"dz", "number of slices"},
  {"dt", "number of time points"},
  {"dx_resampled", "intended x grid dimension after resampling"},
  {"dx", "x grid dimension"},
  {"dy", "y grid dimension"},
  {"dq", "q grid dimension"},
  {"dp", "number of samples per shot"},
  {"start_offset", "starting byte offset in data file"},
  {"autoscale_range", "range for autoscaling"},
  {"dimstr", "order of data dimensions (left fastest)"},
  {"reorder", "slices require reordering"},
  {"reorder_pattern", "slice acqisition order"},
  {"voxel_x", "X voxel size (mm)"},
  {"voxel_y", "Y voxel size (mm)"},
  {"voxel_z", "Z voxel size including gap (mm)"},
  {"slice_gap", "Unscanned space btwn adjacent slices (mm)"},
  {"slice_thickness", "Scanned thickness of a slice (mm)"},
  {"blf.0", "scan volume bottom left front X (mm)"},
  {"blf.1", "scan volume bottom left front Y (mm)"},
  {"blf.2", "scan volume bottom left front Z (mm)"},
  {"tlf.0", "scan volume top left front X (mm)"},
  {"tlf.1", "scan volume top left front Y (mm)"},
  {"tlf.2", "scan volume top left front Z (mm)"},
  {"brf.0", "scan volume bottom right front X (mm)"},
  {"brf.1", "scan volume bottom right front Y (mm)"},
  {"brf.2", "scan volume bottom right front Z (mm)"},
  {"trf.0", "scan volume top right front X (mm)"},
  {"trf.1", "scan volume top right front Y (mm)"},
  {"trf.2", "scan volume top right front Z (mm)"},
  {"blb.0", "scan volume bottom left back X (mm)"},
  {"blb.1", "scan volume bottom left back Y (mm)"},
  {"blb.2", "scan volume bottom left back Z (mm)"},
  {"tlb.0", "scan volume top left back X (mm)"},
  {"tlb.1", "scan volume top left back Y (mm)"},
  {"tlb.2", "scan volume top left back Z (mm)"},
  {"brb.0", "scan volume bottom right back X (mm)"},
  {"brb.1", "scan volume bottom right back Y (mm)"},
  {"brb.2", "scan volume bottom right back Z (mm)"},
  {"trb.0", "scan volume top right back X (mm)"},
  {"trb.1", "scan volume top right back Y (mm)"},
  {"trb.2", "scan volume top right back Z (mm)"},
  {"flip", "flip angle (degrees)"},
  {NULL, NULL}
};

/* The following tag name translations happen when key value pairs 
 * are written to the Pgh MRI file.  The file handlers may add other 
 * translations. A translation to the empty string "" causes the 
 * given key to not be written to the MRI file.
 */
static char* nameTransTable[][2]= {
  {"skip", ""},
  {"sliceskip", ""},
  {"start_offset", ""},
  {"chunkname", ""},
  {"datatype_in", ""},
  {"datatype_out", ""},
  {"handler_datatype_out", ""},
  {"big_endian_input", ""},
  {"big_endian_output", ""},
  {"definitions", ""},
  {"external_names", ""},
  {"dimstr", ""},
  {"ignoreheader", ""},
  {"cl_extent_string", ""},
  {"cl_dim_string", ""},
  {"cl_reorder", ""},
  {"cl_big_endian_input", ""},
  {"cl_autoscale", ""},
  {"cl_autoscale_range", ""},
  {"voxel_x", "voxel_spacing.x"},
  {"voxel_y", "voxel_spacing.y"},
  {"voxel_z", "voxel_spacing.z"},
  {"slice_thickness", "voxel_size.z"},
  {"fov_x", "fov.x"},
  {"fov_y", "fov.y"},
  {"fov_z", "fov.z"},
  {"TR", "tr"},
  {"TE", "te"},
  {"tag", "scan.id"},
  {"chunkfile", ""},
  {"multi", ""},
  {NULL, NULL} /* ends list */
};

/* Some key strings have expected value types. In addition to these,
 * extents ("d?") are KV_INT and skips ("skip.?") are KV_LONG.
 */
static KeyTypePair keyTypeTable[]= {
  {"datatype", KV_INT},
  {"datatype_out", KV_INT},
  {"handler_datatype_out", KV_INT},
  {"datatype_in", KV_INT},
  {"definitions", KV_HASH},
  {"external_names", KV_HASH},
  {"expected_types", KV_HASH},
  {"dimstr", KV_STRING},
  {"cl_dim_string", KV_STRING},
  {"cl_extent_string", KV_STRING},
  {"chunkname", KV_STRING},
  {"chunkfile", KV_STRING},
  {"start_offset", KV_LONG},
  {"skip", KV_LONG},
  {"sliceskip", KV_LONG},
  {"tag", KV_STRING},
  {"big_endian_input", KV_BOOLEAN},
  {"big_endian_output", KV_BOOLEAN},
  {"cl_big_endian_input", KV_BOOLEAN},
  {"cl_big_endian_output", KV_BOOLEAN},
  {"pulse_seq", KV_STRING},
  {"plane", KV_STRING},
  {"dx_resampled", KV_INT},
  {"TR", KV_INT},
  {"TE", KV_INT},
  {"image_x", KV_DOUBLE},
  {"image_y", KV_DOUBLE},
  {"fov_x", KV_DOUBLE},
  {"fov_y", KV_DOUBLE},
  {"fov_z", KV_DOUBLE},
  {"voxel_x", KV_DOUBLE},
  {"voxel_y", KV_DOUBLE},
  {"voxel_z", KV_DOUBLE},
  {"overscan", KV_INT},
  {"slice_thickness", KV_DOUBLE},
  {"slice_gap", KV_DOUBLE},
  {"date", KV_STRING},
  {"time", KV_STRING},
  {"xchop", KV_BOOLEAN},
  {"ychop", KV_BOOLEAN},
  {"autoscale", KV_BOOLEAN}, 
  {"autoscale_range", KV_DOUBLE}, 
  {"reorder", KV_BOOLEAN},
  {"cl_reorder", KV_BOOLEAN},
  {"ignoreheader", KV_BOOLEAN},
  {"multi", KV_BOOLEAN},
  {"flip", KV_DOUBLE},
  {"tlf.0", KV_DOUBLE},
  {"tlf.1", KV_DOUBLE},
  {"tlf.2", KV_DOUBLE},
  {"trf.0", KV_DOUBLE},
  {"trf.1", KV_DOUBLE},
  {"trf.2", KV_DOUBLE},
  {"blf.0", KV_DOUBLE},
  {"blf.1", KV_DOUBLE},
  {"blf.2", KV_DOUBLE},
  {"brf.0", KV_DOUBLE},
  {"brf.1", KV_DOUBLE},
  {"brf.2", KV_DOUBLE},
  {"tlb.0", KV_DOUBLE},
  {"tlb.1", KV_DOUBLE},
  {"tlb.2", KV_DOUBLE},
  {"trb.0", KV_DOUBLE},
  {"trb.1", KV_DOUBLE},
  {"trb.2", KV_DOUBLE},
  {"blb.0", KV_DOUBLE},
  {"blb.1", KV_DOUBLE},
  {"blb.2", KV_DOUBLE},
  {"brb.0", KV_DOUBLE},
  {"brb.1", KV_DOUBLE},
  {"brb.2", KV_DOUBLE},
  {"nifti_qform_code",KV_INT},
  {NULL, KV_INT} /* ends list */
};

/* Maps SRDR_* types to their sizes and names:
 * SRDR_UINT8, SRDR_INT16, SRDR_UINT16, SRDR_INT32, SRDR_FLOAT32,
 * SRDR_FLOAT64, SRDR_INT64 in that order.
 */
long srdrTypeSize[7] = { 1, 2, 2, 4, 4, 8, 8 };
char* srdrTypeName[7] = { "uint8", "int16", "uint16", "int32", "float32", 
		      "float64", "int64" };

/* Maps libmri MRI_Datatype types to their sizes and names:
 * MRI_UINT8, MRI_INT16, MRI_INT32, MRI_FLOAT32, MRI_FLOAT64, MRI_INT64
 * in that order.
 */
long libmriTypeSize[6] = { 1, 2, 4, 4, 8, 8 };
char* libmriTypeName[6] = { "uint8", "int16", "int32", "float32", 
		      "float64", "int64" };

/* Maps libmri MRI_ArrayType types to their sizes and names:
 * MRI_RAW, MRI_UNSIGNED_CHAR, MRI_SHORT, MRI_INT, MRI_FLOAT, 
 * MRI_DOUBLE, MRI_LONG, MRI_LONGLONG in that order.
 */
long libmriArrayTypeSize[8] = { 0, 1, 2, 4, 4, 8, 4, 8 };
char* libmriArrayTypeName[8] = { "raw", "unsigned_char", "short", "int", 
				 "float", "double", "long", "longlong" };

void loadStringTransTable( KVHash* kvh, char* tbl[][2] )
{
  int i;
  for (i=0; tbl[i][0]!=NULL; i++) kvDefString(kvh, tbl[i][0], tbl[i][1]);
}

void loadStringTypeTable( KVHash* kvh, KeyTypePair* tbl )
{
  int i;
  for (i=0; tbl[i].name!=NULL; i++)     
    kvDefInt(kvh, tbl[i].name, (int)tbl[i].type);
}

void initInfoHash(KVHash* info) 
{
  int i;
  KVHash* defs= kvFactory(KV_DEFAULT_SIZE);
  KVHash* extNames= kvFactory(KV_DEFAULT_SIZE);
  KVHash* types= kvFactory(KV_DEFAULT_SIZE);

  /* We'll want sub-hashes for definitions and name translations */
  kvDefHash(info,"definitions",defs);
  kvDefHash(info,"external_names",extNames);
  kvDefHash(info,"expected_types",types);

  /* Fill the name translation hash table */
  loadStringTransTable( extNames, nameTransTable );

  /* Define some name definitions */
  loadStringTransTable( defs, nameDefTable );

  /* Expected key value types (not exhaustive) */
  loadStringTypeTable( types, keyTypeTable );
}

int stringTableLookup(const char* s, const char** tbl)
{
  int i= 0;
  while (tbl[i] != NULL) {
    if (!strcasecmp(s,tbl[i])) return i;
    i++;
  }
  return -1;
}

