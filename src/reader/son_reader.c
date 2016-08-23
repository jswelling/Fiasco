/************************************************************
 *                                                          *
 *  son_reader.c                                         *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 2006 Department of Statistics,         *
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
#include <time.h>
#include <assert.h>
#include "mri.h"
#if (SGI64 || SGI5 || SGIMP)
#include <bstring.h>
#endif

#include "bio.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"
#include "smartreader.h"

static char rcsid[] = "$Id: son_reader.c,v 1.1 2007/03/22 00:05:30 welling Exp $";

/* Notes-
 * -Should fdouble be 'long double' (size 12 bytes)?
 */

#ifdef USE_SON

/* The following serves the function of the archaic 'machine.h' file
 * provided with the SON includes.
 */
typedef short BOOLEAN;
typedef char* LPSTR;
typedef const char* LPCSTR;
typedef unsigned short WORD;
typedef unsigned long DWORD;
typedef unsigned char BYTE;
typedef double fdouble;
typedef long Coord;
typedef void* THandle;
#include <sonintl.h>

/* A reader for SON files */
typedef struct son_data_struct {
  time_t unixTimeSecs;
} SonData;

#ifdef never
typedef struct son_block_info {
  int type;
  int parentType;
  long long start_offset;
} SonBlockInfo;

typedef struct son_data_buffer_info {
  long long offset;
  long size; /* in bytes */
  int type;
} SonDataBufferInfo;

static SonDataBufferInfo* allocFDBInfo( long long offset, long size, 
					 int type )
{
  SonDataBufferInfo* result= NULL;
  if (!(result=(SonDataBufferInfo*)malloc(sizeof(SonDataBufferInfo))))
    Abort("%s:%d: unable to allocate %d bytes!\n",
	  __FILE__,__LINE__,sizeof(SonDataBufferInfo));
  result->offset= offset;
  result->size= size;
  result->type= type;
  return result;
}

static void freeFDBInfo( void* v )
{
  free(v);
}

static const char* sonTypeToString( int type )
{
  switch (type) {
  case SONT_VOID: return "VOID";
  case SONT_BYTE: return "BYTE";
  case SONT_SHORT: return "SHORT";
  case SONT_INT: return "INT";
  case SONT_FLOAT: return "FLOAT";
  case SONT_DOUBLE: return "DOUBLE";
  case SONT_JULIAN: return "JULIAN"; 
  case SONT_USHORT: return "USHORT";
  case SONT_UINT: return "UINT";
  case SONT_STRING: return "STRING";
    /* case SONT_ASCII: return "ASCII"; */ /* redundant */
  case SONT_DAU_PACK13: return "DAU_PACK13";
  case SONT_DAU_PACK14: return "DAU_PACK14";
  case SONT_DAU_PACK16: return "DAU_PACK16";
  case SONT_OLD_PACK: return "OLD_PACK";
  case SONT_CH_INFO_STRUCT: return "CH_INFO_STRUCT";
  case SONT_ID_STRUCT: return "ID_STRUCT";
  case SONT_DIR_ENTRY_STRUCT: return "DIR_ENTRY_STRUCT";
  case SONT_DIG_POINT_STRUCT: return "DIG_POINT_STRUCT";
  case SONT_CH_POS_STRUCT: return "CH_POS_STRUCT";
  case SONT_COORD_TRANS_STRUCT: return "COORD_TRANS_STRUCT";
  case SONT_DIG_STRING_STRUCT: return "DIG_STRING_STRUCT";
  default: return "***UNKNOWN***";
  }
}

static const char* sonBlockToString(int kind)
{
  if (kind==-1) return "TOP-OF-STACK";
  else 
    switch (kind) {
    case SONB_ROOT: return "ROOT";
    case SONB_MEAS: return "MEAS";
    case SONB_MEAS_INFO: return "MEAS_INFO";
    case SONB_RAW_DATA: return "RAW_DATA";
    case SONB_PROCESSED_DATA: return "PROCESSED_DATA";
    case SONB_EVOKED: return "EVOKED";
      /* case SONB_MCG_AVE: return "MCG_AVE"; */ /* redundant */
    case SONB_ASPECT: return "ASPECT";
    case SONB_SUBJECT: return "SUBJECT";
    case SONB_ISOTRAK: return "ISOTRAK";
    case SONB_HPI_MEAS: return "HPI_MEAS";
    case SONB_HPI_RESULT: return "HPI_RESULT";
    case SONB_HPI_COIL: return "HPI_COIL";
    case SONB_PROJECT: return "PROJECT";
    case SONB_CONTINUOUS_DATA: return "CONTINUOUS_DATA";
    case SONB_VOID: return "VOID";
    case SONB_EVENTS: return "EVENTS";
    case SONB_INDEX: return "INDEX";
    case SONB_DACQ_PARS: return "DACQ_PARS";
    case SONB_REF: return "REF";
    case SONB_SMSH_RAW_DATA: return "SMSH_RAW_DATA";
    case SONB_SMSH_ASPECT: return "SMSH_ASPECT";
    case SONB_MRI: return "MRI";
    case SONB_MRI_SET: return "MRI_SET";
    case SONB_MRI_SLICE: return "MRI_SLICE";
    case SONB_MRI_SCENERY: return "MRI_SCENERY";
    case SONB_MRI_SCENE: return "MRI_SCENE";
    case SONB_MRI_SEG: return "MRI_SEG";
    case SONB_MRI_SEG_REGION: return "MRI_SEG_REGION";
    case SONB_SPHERE: return "SPHERE";
    case SONB_BEM: return "BEM";
    case SONB_BEM_SURF: return "BEM_SURF";
    case SONB_XFIT_PROJ: return "XFIT_PROJ";
    case SONB_XFIT_PROJ_ITEM: return "XFIT_PROJ_ITEM";
    case SONB_XFIT_AUX: return "XFIT_AUX";
    case SONB_VOL_INFO: return "VOL_INFO";
    case SONB_DATA_CORRECTION: return "DATA_CORRECTION";
    case SONB_CHANNEL_DECOUPLER: return "CHANNEL_DECOUPLER";
    case SONB_SSS_INFO: return "SSS_INFO";
    case SONB_SSS_CAL_ADJUST: return "SSS_CAL_ADJUST";
    case SONB_SSS_ST_INFO: return "SSS_ST_INFO";
    case SONB_SSS_BASES: return "SSS_BASES";
    case SONB_SMARTSHIELD: return "SMARTSHIELD";
    case SONB_PROCESSING_HISTORY: return "PROCESSING_HISTORY";
    case SONB_PROCESSING_RECORD: return "PROCESSING_RECORD";
    default: return "***UNKNOWN***";
    }
}

static const char* sonKindToString(int v)
{
  switch (v) {
  case SON_NEW_FILE: return "NEW_FILE";
  case SON_CLOSE_FILE: return "CLOSE_FILE";
  case SON_DISCARD_FILE: return "DISCARD_FILE";
  case SON_ERROR_MESSAGE: return "ERROR_MESSAGE";
  case SON_SUSPEND_READING: return "SUSPEND_READING";
  case SON_FATAL_ERROR_MESSAGE: return "FATAL_ERROR_MESSAGE";
  case SON_CONNECTION_CHECK: return "CONNECTION_CHECK";
  case SON_SUSPEND_FILING: return "SUSPEND_FILING";
  case SON_RESUME_FILING: return "RESUME_FILING";
  case SON_RAW_PREBASE: return "RAW_PREBASE";
  case SON_RAW_PICK_LIST: return "RAW_PICK_LIST";
  case SON_ECHO: return "ECHO";
  case SON_RESUME_READING: return "RESUME_READING";
  case SON_DACQ_SYSTEM_TYPE: return "DACQ_SYSTEM_TYPE";
  case SON_SELECT_RAW_CH: return "SELECT_RAW_CH";
  case SON_PLAYBACK_MODE: return "PLAYBACK_MODE";
  case SON_CONTINUE_FILE: return "CONTINUE_FILE";
  case SON_JITTER_MAX: return "JITTER_MAX";
  case SON_DECIMATION_FACTOR: return "DECIMATION_FACTOR";
#ifdef _DATA_SERVER
  case SON_MEM_DATA_BUFFER: return "MEM_DATA_BUFFER";
#endif
  case SON_FILE_ID: return "FILE_ID";
  case SON_DIR_POINTER: return "DIR_POINTER";
  case SON_DIR: return "DIR";
  case SON_BLOCK_ID: return "BLOCK_ID";
  case SON_BLOCK_START: return "BLOCK_START";
  case SON_BLOCK_END: return "BLOCK_END";
  case SON_FREE_LIST: return "FREE_LIST";
  case SON_FREE_BLOCK: return "FREE_BLOCK";
  case SON_NOP: return "NOP";
  case SON_PARENT_FILE_ID: return "PARENT_FILE_ID";
  case SON_PARENT_BLOCK_ID: return "PARENT_BLOCK_ID";
  case SON_BLOCK_NAME: return "BLOCK_NAME";
  case SON_BLOCK_VERSION: return "BLOCK_VERSION";
  case SON_CREATOR: return "CREATOR";
  case SON_MODIFIER: return "MODIFIER";
  case SON_REF_ROLE: return "REF_ROLE";
  case SON_REF_FILE_ID: return "REF_FILE_ID";
  case SON_REF_FILE_NUM: return "REF_FILE_NUM";
  case SON_REF_FILE_NAME: return "REF_FILE_NAME";
  case SON_REF_BLOCK_ID: return "REF_BLOCK_ID";
  case SON_DACQ_PARS: return "DACQ_PARS";
  case SON_DACQ_STIM: return "DACQ_STIM";
  case SON_NCHAN: return "NCHAN";
  case SON_SFREQ: return "SFREQ";
  case SON_DATA_PACK: return "DATA_PACK";
  case SON_CH_INFO: return "CH_INFO";
  case SON_MEAS_DATE: return "MEAS_DATE";
  case SON_SUBJECT: return "SUBJECT";
  case SON_COMMENT: return "COMMENT";
  case SON_NAVE: return "NAVE";
  case SON_FIRST_SAMPLE: return "FIRST_SAMPLE";
  case SON_LAST_SAMPLE: return "LAST_SAMPLE";
  case SON_ASPECT_KIND: return "ASPECT_KIND";
  case SON_REF_EVENT: return "REF_EVENT";
  case SON_EXPERIMENTER: return "EXPERIMENTER";
  case SON_DIG_POINT: return "DIG_POINT";
  case SON_CH_POS_VEC: return "CH_POS_VEC";
  case SON_HPI_SLOPES: return "HPI_SLOPES";
  case SON_HPI_NCOIL: return "HPI_NCOIL";
  case SON_REQ_EVENT: return "REQ_EVENT";
  case SON_REQ_LIMIT: return "REQ_LIMIT";
  case SON_LOWPASS: return "LOWPASS";
  case SON_BAD_CHS: return "BAD_CHS";
  case SON_ARTEF_REMOVAL: return "ARTEF_REMOVAL";
  case SON_COORD_TRANS: return "COORD_TRANS";
  case SON_HIGHPASS: return "HIGHPASS";
  case SON_CH_CALS_VEC: return "CH_CALS_VEC";
  case SON_HPI_BAD_CHS: return "HPI_BAD_CHS";
  case SON_HPI_CORR_COEFF: return "HPI_CORR_COEFF";
  case SON_EVENT_COMMENT: return "EVENT_COMMENT";
  case SON_NO_SAMPLES: return "NO_SAMPLES";
  case SON_FIRST_TIME: return "FIRST_TIME";
  case SON_SUBAVE_SIZE: return "SUBAVE_SIZE";
  case SON_SUBAVE_FIRST: return "SUBAVE_FIRST";
  case SON_NAME: return "NAME";
    /* case SON_DESCRIPTION: return "DESCRIPTION"; */ /* redundant */
  case SON_DIG_STRING: return "DIG_STRING";
  case SON_LINE_FREQ: return "LINE_FREQ";
  case SON_HPI_COIL_FREQ: return "HPI_COIL_FREQ";
  case SON_HPI_COIL_MOMENTS: return "HPI_COIL_MOMENTS";
  case SON_HPI_FIT_GOODNESS: return "HPI_FIT_GOODNESS";
  case SON_HPI_FIT_ACCEPT: return "HPI_FIT_ACCEPT";
  case SON_HPI_FIT_GOOD_LIMIT: return "HPI_FIT_GOOD_LIMIT";
  case SON_HPI_FIT_DIST_LIMIT: return "HPI_FIT_DIST_LIMIT";
  case SON_HPI_COIL_NO: return "HPI_COIL_NO";
  case SON_HPI_COILS_USED: return "HPI_COILS_USED";
  case SON_HPI_DIGITIZATION_ORDER: return "HPI_DIGITIZATION_ORDER";
  case SON_CH_SCAN_NO: return "CH_SCAN_NO";
  case SON_CH_LOGICAL_NO: return "CH_LOGICAL_NO";
  case SON_CH_KIND: return "CH_KIND";
  case SON_CH_RANGE: return "CH_RANGE";
  case SON_CH_CAL: return "CH_CAL";
  case SON_CH_POS: return "CH_POS";
  case SON_CH_UNIT: return "CH_UNIT";
  case SON_CH_UNIT_MUL: return "CH_UNIT_MUL";
  case SON_CH_DACQ_NAME: return "CH_DACQ_NAME";
  case SON_SSS_FRAME: return "SSS_FRAME";
  case SON_SSS_JOB: return "SSS_JOB";
  case SON_SSS_ORIGIN: return "SSS_ORIGIN";
  case SON_SSS_ORD_IN: return "SSS_ORD_IN";
  case SON_SSS_ORD_OUT: return "SSS_ORD_OUT";
  case SON_SSS_NMAG: return "SSS_NMAG";
  case SON_SSS_COMPONENTS: return "SSS_COMPONENTS";
  case SON_SSS_CAL_CHANS: return "SSS_CAL_CHANS";
  case SON_SSS_CAL_CORRS: return "SSS_CAL_CORRS";
  case SON_SSS_ST_CORR: return "SSS_ST_CORR";
  case SON_SSS_BASE_IN: return "SSS_BASE_IN";
  case SON_SSS_BASE_OUT: return "SSS_BASE_OUT";
  case SON_SSS_BASE_VIRT: return "SSS_BASE_VIRT";
  case SON_SSS_NORM: return "SSS_NORM";
  case SON_DATA_BUFFER: return "DATA_BUFFER";
  case SON_DATA_SKIP: return "DATA_SKIP";
  case SON_EPOCH: return "EPOCH";
  case SON_DATA_SKIP_SAMP: return "DATA_SKIP_SAMP";
  case SON_SUBJ_ID: return "SUBJ_ID";
  case SON_SUBJ_FIRST_NAME: return "SUBJ_FIRST_NAME";
  case SON_SUBJ_MIDDLE_NAME: return "SUBJ_MIDDLE_NAME";
  case SON_SUBJ_LAST_NAME: return "SUBJ_LAST_NAME";
  case SON_SUBJ_BIRTH_DAY: return "SUBJ_BIRTH_DAY";
  case SON_SUBJ_SEX: return "SUBJ_SEX";
  case SON_SUBJ_HAND: return "SUBJ_HAND";
  case SON_SUBJ_WEIGHT: return "SUBJ_WEIGHT";
  case SON_SUBJ_HEIGHT: return "SUBJ_HEIGHT";
  case SON_SUBJ_COMMENT: return "SUBJ_COMMENT";
  case SON_SUBJ_HIS_ID: return "SUBJ_HIS_ID";
  case SON_PROJ_ID: return "PROJ_ID";
  case SON_PROJ_NAME: return "PROJ_NAME";
  case SON_PROJ_AIM: return "PROJ_AIM";
  case SON_PROJ_PERSONS: return "PROJ_PERSONS";
  case SON_PROJ_COMMENT: return "PROJ_COMMENT";
  case SON_EVENT_CHANNELS: return "EVENT_CHANNELS";
  case SON_EVENT_LIST: return "EVENT_LIST";
  case SON_SQUID_BIAS: return "SQUID_BIAS";
  case SON_SQUID_OFFSET: return "SQUID_OFFSET";
  case SON_SQUID_GATE: return "SQUID_GATE";
  case SON_DECOUPLER_MATRIX: return "DECOUPLER_MATRIX";
  case SON_SPARSE_CH_NAME_LIST: return "SPARSE_CH_NAME_LIST";
  case SON_REF_PATH: return "REF_PATH";
    /* case SON_MRI_SOURCE_PATH: return "MRI_SOURCE_PATH"; */ /* redundant */
  case SON_MRI_SOURCE_FORMAT: return "MRI_SOURCE_FORMAT";
  case SON_MRI_PIXEL_ENCODING: return "MRI_PIXEL_ENCODING";
  case SON_MRI_PIXEL_DATA_OFFSET: return "MRI_PIXEL_DATA_OFFSET";
  case SON_MRI_PIXEL_SCALE: return "MRI_PIXEL_SCALE";
  case SON_MRI_PIXEL_DATA: return "MRI_PIXEL_DATA";
  case SON_MRI_PIXEL_OVERLAY_ENCODING: return "MRI_PIXEL_OVERLAY_ENCODING";
  case SON_MRI_PIXEL_OVERLAY_DATA: return "MRI_PIXEL_OVERLAY_DATA";
  case SON_MRI_BOUNDING_BOX: return "MRI_BOUNDING_BOX";
  case SON_MRI_WIDTH: return "MRI_WIDTH";
  case SON_MRI_WIDTH_M: return "MRI_WIDTH_M";
  case SON_MRI_HEIGHT: return "MRI_HEIGHT";
  case SON_MRI_HEIGHT_M: return "MRI_HEIGHT_M";
  case SON_MRI_DEPTH: return "MRI_DEPTH";
  case SON_MRI_DEPTH_M: return "MRI_DEPTH_M";
  case SON_MRI_THICKNESS: return "MRI_THICKNESS";
  case SON_MRI_SCENE_AIM: return "MRI_SCENE_AIM";
  case SON_MRI_ORIG_SOURCE_PATH: return "MRI_ORIG_SOURCE_PATH";
  case SON_MRI_ORIG_SOURCE_FORMAT: return "MRI_ORIG_SOURCE_FORMAT";
  case SON_MRI_ORIG_PIXEL_ENCODING: return "MRI_ORIG_PIXEL_ENCODING";
  case SON_MRI_ORIG_PIXEL_DATA_OFFSET: return "MRI_ORIG_PIXEL_DATA_OFFSET";
  case SON_MRI_VOXEL_DATA: return "MRI_VOXEL_DATA";
  case SON_MRI_VOXEL_ENCODING: return "MRI_VOXEL_ENCODING";
  case SON_MRI_MRILAB_SETUP: return "MRI_MRILAB_SETUP";
  case SON_MRI_SEG_REGION_ID: return "MRI_SEG_REGION_ID";
  case SON_CONDUCTOR_MODEL_KIND: return "CONDUCTOR_MODEL_KIND";
  case SON_SPHERE_ORIGIN: return "SPHERE_ORIGIN";
  case SON_SPHERE_COORD_FRAME: return "SPHERE_COORD_FRAME";
  case SON_SPHERE_LAYERS: return "SPHERE_LAYERS";
  case SON_BEM_SURF_ID: return "BEM_SURF_ID";
  case SON_BEM_SURF_NAME: return "BEM_SURF_NAME";
  case SON_BEM_SURF_NNODE: return "BEM_SURF_NNODE";
  case SON_BEM_SURF_NTRI: return "BEM_SURF_NTRI";
  case SON_BEM_SURF_NODES: return "BEM_SURF_NODES";
  case SON_BEM_SURF_TRIANGLES: return "BEM_SURF_TRIANGLES";
  case SON_BEM_SURF_NORMALS: return "BEM_SURF_NORMALS";
  case SON_BEM_SURF_CURVS: return "BEM_SURF_CURVS";
  case SON_BEM_SURF_CURV_VALUES: return "BEM_SURF_CURV_VALUES";
  case SON_BEM_POT_SOLUTION: return "BEM_POT_SOLUTION";
  case SON_BEM_APPROX: return "BEM_APPROX";
  case SON_BEM_COORD_FRAME: return "BEM_COORD_FRAME";
  case SON_BEM_SIGMA: return "BEM_SIGMA";
  case SON_SOURCE_DIPOLE: return "SOURCE_DIPOLE";
  case SON_XFIT_LEAD_PRODUCTS: return "XFIT_LEAD_PRODUCTS";
  case SON_XFIT_MAP_PRODUCTS: return "XFIT_MAP_PRODUCTS";
  case SON_XFIT_GRAD_MAP_PRODUCTS: return "XFIT_GRAD_MAP_PRODUCTS";
  case SON_XFIT_VOL_INTEGRATION: return "XFIT_VOL_INTEGRATION";
  case SON_XFIT_INTEGRATION_RADIUS: return "XFIT_INTEGRATION_RADIUS";
  case SON_XFIT_CONDUCTOR_MODEL_NAME: return "XFIT_CONDUCTOR_MODEL_NAME";
  case SON_XFIT_CONDUCTOR_MODEL_TRANS_NAME: return "XFIT_CONDUCTOR_MODEL_TRANS_NAME";
  case SON_PROJ_ITEM_KIND: return "PROJ_ITEM_KIND";
  case SON_PROJ_ITEM_TIME: return "PROJ_ITEM_TIME";
    /* case SON_PROJ_ITEM_DIPOLE: return "PROJ_ITEM_DIPOLE"; */ /* redundant */
  case SON_PROJ_ITEM_IGN_CHS: return "PROJ_ITEM_IGN_CHS";
  case SON_PROJ_ITEM_NVEC: return "PROJ_ITEM_NVEC";
  case SON_PROJ_ITEM_VECTORS: return "PROJ_ITEM_VECTORS";
    /* case SON_PROJ_ITEM_COMMENT: return "PROJ_ITEM_COMMENT"; */ /* redundant */
    /* case SON_PROJ_ITEM_DESCRIPTION: return "PROJ_ITEM_DESCRIPTION"; */ /* redundant */
  case SON_PROJ_ITEM_DEFINITION: return "PROJ_ITEM_DEFINITION";
    /* case SON_PROJ_ITEM_CH_NAME_LIST: return "PROJ_ITEM_CH_NAME_LIST"; */ /* redundant */
    /*  case SON_XFIT_PROJ_ITEM_KIND: return "XFIT_PROJ_ITEM_KIND"; */ /* redundant */
    /* case SON_XFIT_PROJ_ITEM_TIME: return "XFIT_PROJ_ITEM_TIME"; */ /* redundant */
    /* case SON_XFIT_PROJ_ITEM_DIPOLE: return "XFIT_PROJ_ITEM_DIPOLE"; */ /* redundant */
    /* case SON_XFIT_PROJ_ITEM_IGN_CHS: return "XFIT_PROJ_ITEM_IGN_CHS"; */ /*redundant */
    /* case SON_XFIT_PROJ_ITEM_NVEC: return "XFIT_PROJ_ITEM_NVEC"; */ /*redundant*/
    /* case SON_XFIT_PROJ_ITEM_VECTORS: return "XFIT_PROJ_ITEM_VECTORS"; */ /*redundant*/
    /* case SON_XFIT_PROJ_ITEM_COMMENT: return "XFIT_PROJ_ITEM_COMMENT"; */ /*redundant*/
  case SON_XPLOTTER_LAYOUT: return "XPLOTTER_LAYOUT";
  case SON_VOL_ID: return "VOL_ID";
  case SON_VOL_NAME: return "VOL_NAME";
  case SON_VOL_OWNER_ID: return "VOL_OWNER_ID";
  case SON_VOL_OWNER_NAME: return "VOL_OWNER_NAME";
  case SON_VOL_OWNER_REAL_NAME: return "VOL_OWNER_REAL_NAME";
  case SON_VOL_TYPE: return "VOL_TYPE";
  case SON_VOL_HOST: return "VOL_HOST";
  case SON_VOL_REAL_ROOT: return "VOL_REAL_ROOT";
  case SON_VOL_SYMBOLIC_ROOT: return "VOL_SYMBOLIC_ROOT";
  case SON_VOL_MOUNT_POINT: return "VOL_MOUNT_POINT";
  case SON_VOL_BLOCKS: return "VOL_BLOCKS";
  case SON_VOL_FREE_BLOCKS: return "VOL_FREE_BLOCKS";
  case SON_VOL_AVAIL_BLOCKS: return "VOL_AVAIL_BLOCKS";
  case SON_VOL_BLOCK_SIZE: return "VOL_BLOCK_SIZE";
  case SON_VOL_DIRECTORY: return "VOL_DIRECTORY";
  case SON_INDEX_KIND: return "INDEX_KIND";
  case SON_INDEX: return "INDEX";
  default: return "*** Unknown SON kind ***";
  }
}

static const char* sonAspectToString( int kind )
{
  switch (kind) {
    /* Normal average of epochs */
  case SONV_ASPECT_AVERAGE: return "AVERAGE";
    /* Std. error of mean */
  case SONV_ASPECT_STD_ERR: return "STD_ERR";
    /* Single epoch cut out from the continuous data */
  case SONV_ASPECT_SINGLE: return "SINGLE";
    /* ??? */
  case SONV_ASPECT_SUBAVERAGE: return "SUBAVERAGE";
    /* Alternating subaverage */
  case SONV_ASPECT_ALTAVERAGE: return "ALTAVERAGE";
    /* A sample cut out by graph */
  case SONV_ASPECT_SAMPLE: return "SAMPLE";
    /* Power density spectrum */
  case SONV_ASPECT_POWER_DENSITY: return "POWER_DENSITY";
    /* Dipole amplitude curve */
  case SONV_ASPECT_DIPOLE_WAVE: return "DIPOLE_WAVE";
  default: return "*** Unknown SON aspect ***";
  }
}

static sonTagRec* readSonTag(FILE* f)
{
  sonTagRec* result= NULL;

  if (!(result=(sonTagRec*)malloc(sizeof(sonTagRec))))
    Abort("%s:%d: unable to allocate %d bytes!\n",__FILE__,__LINE__,
	  sizeof(sonTagRec));

  result->kind= FRdInt32(f);
  result->type= FRdInt32(f);
  result->size= FRdInt32(f);
  result->next= FRdInt32(f);
  if (result->size<0 || result->size>(1<<30)) {
    if (debug) fprintf(stderr,"Found tag with implausible data size %ld\n",
		       result->size);
    free(result);
    return NULL;
  }
  if (!(result->data=(son_data_t*)malloc(result->size)))
    Abort("%s:%d: unable to allocate %d bytes!\n",__FILE__,__LINE__,
	  result->size);
  if (debug)
    fprintf(stderr,"kind %s, type %s, size %d\n",
	    sonKindToString(result->kind),
	    sonTypeToString(result->type),
	    result->size);
  switch (result->type) {
  case SONT_VOID: 
    {
      FRdUInt8Array(f,(unsigned char*)result->data,result->size/1);
    }
    break;
  case SONT_BYTE: 
    {
      FRdUInt8Array(f,(unsigned char*)result->data,result->size/1);
    }
    break;
  case SONT_SHORT: 
    {
      FRdInt16Array(f,(short*)result->data,result->size/2);
    }
    break;
  case SONT_INT: 
    {
      FRdInt32Array(f,(int*)result->data,result->size/4);
    }
    break;
  case SONT_FLOAT: 
    {
      FRdFloat32Array(f,(float*)result->data,result->size/4);
    }
    break;
  case SONT_DOUBLE: 
    {
      FRdFloat64Array(f,(double*)result->data,result->size/8);
    }
    break;
  case SONT_JULIAN: 
    {
      FRdInt32Array(f,(int*)result->data,result->size/4);
    }
    break;
  case SONT_USHORT: 
    {
      FRdInt16Array(f,(short*)result->data,result->size/2);
    }
    break;
  case SONT_UINT: 
    {
      FRdInt32Array(f,(int*)result->data,result->size/4);
    }
    break;
  case SONT_STRING: 
    {
      char* buf;
      if (!(buf=(char*)malloc(result->size + 1)))
	Abort("%s:%d: unable to allocate %d bytes!\n",__FILE__,__LINE__,
	      result->size+1);
      /* Allocated space does not include room for trailing null! */
      FRdUInt8Array(f,(char*)result->data,result->size/1);
      strncpy(buf,(char*)result->data,result->size);
      buf[result->size]= '\0';
      free(result->data);
      result->data= buf;
    }
    break;
  case SONT_DAU_PACK13:
    {
      FRdInt16Array(f,(short*)result->data,result->size/2);
    }
    break;
  case SONT_DAU_PACK14:
    {
      FRdInt16Array(f,(short*)result->data,result->size/2);
    }
    break;
  case SONT_DAU_PACK16:
    {
      FRdInt16Array(f,(short*)result->data,result->size/2);
    }
    break;
  case SONT_OLD_PACK:
    {
      int nsamp= (result->size - 8)/2;
      short* shortPtr= (short*)(((float*)result->data)+2);
      FRdFloat32Array(f,(float*)result->data,2);
      FRdInt16Array(f,shortPtr,nsamp);
    }
    break;
  case SONT_CH_INFO_STRUCT:
    {
      sonChInfoRec* rec= (sonChInfoRec*)result->data;
      rec->scanNo= FRdInt32(f);
      rec->logNo= FRdInt32(f);
      rec->kind= FRdInt32(f);
      rec->range= FRdFloat32(f);
      rec->cal= FRdFloat32(f);
      rec->chpos.coil_type= FRdInt32(f);
      FRdFloat32Array(f,&(rec->chpos.r0[0]),3);
      FRdFloat32Array(f,&(rec->chpos.ex[0]),3);
      FRdFloat32Array(f,&(rec->chpos.ey[0]),3);
      FRdFloat32Array(f,&(rec->chpos.ez[0]),3);
      rec->unit= FRdInt32(f);
      rec->unit_mul= FRdInt32(f);
      FRdUInt8Array(f,rec->ch_name,16);
    }
    break;
  case SONT_ID_STRUCT:
    {
      sonIdRec* rec= (sonIdRec*)result->data;
      rec->version= FRdInt32(f);
      FRdInt32Array(f,rec->machid,2);
      rec->time.secs= FRdInt32(f);
      rec->time.usecs= FRdInt32(f);
    }
    break;
  case SONT_DIR_ENTRY_STRUCT:
    {
      sonDirEntryRec* rec= (sonDirEntryRec*)result->data;
      rec->kind= FRdInt32(f);
      rec->type= FRdInt32(f);
      rec->size= FRdInt32(f);
      rec->pos= FRdInt32(f);
    }
    break;
  case SONT_DIG_POINT_STRUCT:
    {
      sonDigPointRec* rec= (sonDigPointRec*)result->data;
      rec->kind= FRdInt32(f);
      rec->ident= FRdInt32(f);
      FRdFloat32Array(f,rec->r,3);
    }
    break;
  case SONT_CH_POS_STRUCT:
    {
      int n= result->size/44;
      int i;
      sonChPosRec* rec= (sonChPosRec*)result->data;
      for (i=0;i<n;i++) {
	rec->coil_type= FRdInt32(f);
	FRdFloat32Array(f,rec->r0,3);
	FRdFloat32Array(f,rec->ex,3);
	FRdFloat32Array(f,rec->ey,3);
	FRdFloat32Array(f,rec->ez,3);
	rec += 1;
      }
    }
    break;
  case SONT_COORD_TRANS_STRUCT:
    {
      int n= result->size/80;
      int i;
      sonCoordTransRec* rec= (sonCoordTransRec*)result->data;
      for (i=0; i<n; i++) {
	rec->from= FRdInt32(f);
	rec->to= FRdInt32(f);
	FRdFloat32Array(f,(float*)&(rec->rot[0]),9);
	FRdFloat32Array(f,rec->move,3);
	FRdFloat32Array(f,(float*)&(rec->invrot[0]),9);
	FRdFloat32Array(f,rec->invmove,3);
	rec += 1;
      }
    }
    break;
  case SONT_DIG_STRING_STRUCT:
    {
      float* buf= NULL;
      float** bufp= NULL;
      sonDigStringRec* dsr= NULL;
      int i;

      if (!(dsr=(sonDigStringRec*)malloc(sizeof(sonDigStringRec))))
	Abort("%s:%d: unable to allocate %d bytes!\n",
	      __FILE__,__LINE__,sizeof(sonDigStringRec));
      if (result->data) free(result->data);
      result->data= dsr;
      dsr->kind= FRdInt32(f);
      dsr->ident= FRdInt32(f);
      dsr->np= FRdInt32(f);
      dsr->rr= NULL;
      if (!(bufp=(float**)malloc(dsr->np*sizeof(float*))))
	Abort("%s:%d: unable to allocate %d bytes!\n",
	      __FILE__,__LINE__,dsr->np*sizeof(float*));
      if (!(buf=(float*)malloc(dsr->np*3*sizeof(float))))
	Abort("%s:%d: unable to allocate %d bytes!\n",
	      dsr->np*3*sizeof(float));
      FRdFloat32Array(f,buf,dsr->np*3);
      for (i=0; i<dsr->np; i++)
	bufp[i]= buf + 3*i;
      dsr->rr= bufp;
    }
    break;
  default:
    if (debug) 
      fprintf(stderr,
	      "***Encountered unknown tag type %d, kind %d, size %d***!\n",
	      result->type,result->kind,result->size);

  }
  return result;
}

static SRDR_Datatype sonTypeMap( int type )
{
  switch (type) {
  case SONT_BYTE: return SRDR_UINT8;
  case SONT_SHORT: return SRDR_INT16;
  case SONT_INT: return SRDR_INT32;
  case SONT_FLOAT: return SRDR_FLOAT32;
  case SONT_DOUBLE: return SRDR_FLOAT64;
  case SONT_USHORT: return SRDR_UINT16;
  case SONT_DAU_PACK13: return SRDR_INT16;
  case SONT_DAU_PACK14: return SRDR_INT16;
  case SONT_DAU_PACK16: return SRDR_INT16;
  case SONT_OLD_PACK: return SRDR_FLOAT32;
  default: Abort("%s:%d: SON type %s is untranslatable!\n",
		 __FILE__,__LINE__,sonTypeToString(type));
  }
}

static int sonNDataItems( int type, long nbytes )
{
  switch (type) {
  case SONT_BYTE: return nbytes;
  case SONT_SHORT: return nbytes/2;
  case SONT_INT: return nbytes/4;
  case SONT_FLOAT: return nbytes/4;
  case SONT_DOUBLE: return nbytes/8;
  case SONT_USHORT: return nbytes/2;
  case SONT_DAU_PACK13: return nbytes/2;
  case SONT_DAU_PACK14: return nbytes/2;
  case SONT_DAU_PACK16: return nbytes/2;
  case SONT_OLD_PACK: return (nbytes-8)/2;
  default: Abort("%s:%d: internal error: no counting method for type %s!\n",
		 __FILE__,__LINE__,sonTypeToString(type));
  }
}

static void freeTag( sonTagRec* rec )
{
  if (rec->data) {
    if (rec->type==SONT_DIG_STRING_STRUCT) {
      sonDigStringRec* d= (sonDigStringRec*)(rec->data);
      if (d->rr) {
	if (d->rr[0]) free(d->rr[0]);
	free(d->rr);
      }
    }
    free(rec->data);
  }
  free(rec);
}

static long long getNextTagOffset( const sonTagRec* rec, long long offset )
{
  if (rec->next!=SONV_NEXT_SEQ) {
    if (rec->next==SONV_NEXT_NONE) {
      Abort("%s:%d: internal error: ran off end of linked list!\n",
	    __FILE__,__LINE__);
    }
    else {
      offset= rec->next;
      if (debug) fprintf(stderr,"Next sends us to %lld\n",offset);
    }
  }
  else {
    offset += 4*sizeof(son_int_t)+rec->size;
    if (debug) fprintf(stderr,"Following tag in order: %lld\n",offset);
  }
  return offset;
}

static void checkBlockListConsistency( SList* blockList )
{
  SonDataBufferInfo* firstBlock= NULL;
  long blockNum= 0;

  slist_totop(blockList);
  while (!slist_atend(blockList)) {
    SonDataBufferInfo* block= (SonDataBufferInfo*)slist_next(blockList);
    if (!firstBlock) {
      firstBlock= block;
      blockNum= 0;
    }
    else {
      if (block->size != firstBlock->size)
	Abort("%s: Data block %ld is the wrong size!\n",__FILE__,blockNum);
      if (block->type != firstBlock->type)
	Abort("%s: Data block %ld is the wrong type!\n",__FILE__,blockNum);
    }
    blockNum++;
  }
}

static void scanBlockList( FileHandler* self, KVHash* info )
{
  SonData* data= (SonData*)(self->hook);
  KVHash* defs= kvGetHash(info,"definitions");
  SList* rawList= data->rawBufferList;
  SList* epochList= data->epochBufferList;
  SRDR_Datatype srdrType;
  int nRaw= slist_count(rawList);
  int nEpoch= slist_count(epochList);

  if (debug)
    fprintf(stderr,"Raw blocks: %d  Epoch blocks: %d\n",nRaw,nEpoch);

  if (nRaw && nEpoch)
    Abort("%s:%d: file %s contains both raw and epoch data- not supported!\n",
	  __FILE__,__LINE__,self->fileName);

  if (nRaw) {
    SonDataBufferInfo* firstBlock= NULL;
    checkBlockListConsistency(rawList);
    slist_totop(rawList);
    firstBlock= (SonDataBufferInfo*)slist_get(rawList);
    kvDefString(info,"dimstr","cbt");
    kvDefInt(info,"dc",kvGetInt(info,"nchannels"));
    srdrType= sonTypeMap(firstBlock->type);
    kvDefInt(info,"datatype_in",srdrType);
    kvDefInt(info,"handler_datatype_out",srdrType);
    kvDefInt(info,"db",
	     sonNDataItems(firstBlock->type,firstBlock->size)/kvGetInt(info,"dc"));
    kvDefLong(info,"skip.b",0);
    kvDefInt(info,"dt",nRaw);
  }
  else if (nEpoch) {
    SonDataBufferInfo* firstBlock= NULL;
    if (nEpoch != kvGetInt(info,"nchannels"))
      Abort("%s:%d: number of epoch blocks (%d) does not match number of channels!\n",
	    __FILE__,__LINE__,nEpoch);
    checkBlockListConsistency(epochList);
    slist_totop(epochList);
    firstBlock= (SonDataBufferInfo*)slist_get(epochList);
    kvDefString(info,"dimstr","bc");
    kvDefInt(info,"dc",nEpoch);
    srdrType= sonTypeMap(firstBlock->type);
    kvDefInt(info,"datatype_in",srdrType);
    kvDefInt(info,"handler_datatype_out",srdrType);
    kvDefInt(info,"db",sonNDataItems(firstBlock->type,firstBlock->size));
    kvDefLong(info,"skip.b",0);
  }
  else {
    Abort("%s:%d: input file %s contains no data!\n",
	  __FILE__,__LINE__,self->fileName);
  }
}
#endif /* ifdef never */

static void sonProcessHeader( FileHandler* self, KVHash* info, 
			      SList* chunkStack )
{
  SonData* data= (SonData*)(self->hook);
  KVHash* defs= kvGetHash(info,"definitions");
  long long nextOffset= 0;
  int old_bio_big_endian_input= bio_big_endian_input;

#ifdef never
  sonTagRec* rec= NULL;
  int thisIsTheLastTag= 0;
  FILE* f;
  SList* blockStack= slist_create();
  int currentBlockType= -1;
  int numMeasBlocks= 0;

  bio_big_endian_input= 1;
  kvDefBoolean(info,"big_endian_input",1);

  if (!(f=fopen(self->fileName,"r")))
    Abort("%s:%d: unexpectedly unable to open <%s>!\n",
	  __FILE__,__LINE__,self->fileName);

  thisIsTheLastTag= 0;
  while (!thisIsTheLastTag) {
    long long offset= nextOffset;
    bigfile_fseek(f,offset,SEEK_SET);
    rec= readSonTag(f);
    if (rec->next==SONV_NEXT_NONE) {
      thisIsTheLastTag= 1;
    }
    else {
      nextOffset= getNextTagOffset(rec,offset);
    }

    /* Code to handle block start and end */
    if (rec->kind==SON_BLOCK_START) {
      int ival= *(int*)rec->data;
      SonBlockInfo* thisBlockInfo= NULL;
      if (!(thisBlockInfo=(SonBlockInfo*)malloc(sizeof(SonBlockInfo))))
	Abort("%s:%d: unable to allocate %d bytes!\n",
	      __FILE__,__LINE__,sizeof(SonBlockInfo));
      thisBlockInfo->type= ival;
      thisBlockInfo->parentType= currentBlockType;
      thisBlockInfo->start_offset= offset;
      currentBlockType= ival;
      if (debug) fprintf(stderr,"Push block type to %s\n",
			 sonBlockToString(currentBlockType));
      if (currentBlockType==SONB_MEAS && ++numMeasBlocks>1)
	Abort("%s:s: SON file contains multiple MEAS blocks!\n",
	      __FILE__,__LINE__);
      slist_push(blockStack,thisBlockInfo);
    }
    else if (rec->kind==SON_BLOCK_END) {
      int ival= *(int*)rec->data;
      SonBlockInfo* thisBlockInfo= slist_pop(blockStack);
      if (ival != thisBlockInfo->type)
	fprintf(stderr,
		"*** block framing error! (found type %d, expected %d ***\n",
		ival,thisBlockInfo->type);
      currentBlockType= thisBlockInfo->parentType;
      if (debug) fprintf(stderr,"Pop block type to %s\n",
			 sonBlockToString(currentBlockType));
      free(thisBlockInfo);
    }

    switch (currentBlockType) 
      {
      case SONB_SUBJECT:
	switch (rec->kind) 
	  {
	  case SON_SUBJ_ID: 
	    {
	      kvDefInt(info,"subject_id",*(int*)rec->data);
	      kvDefString(defs,"subject_id","subject identifier");
	    }
	    break;
	  case SON_SUBJ_SEX: 
	    {
	      int ival= *(int*)rec->data;
	      if (ival==1)
		kvDefString(info,"sex","male");
	      else if (ival==2)
		kvDefString(info,"sex","female");
	      kvDefString(defs,"sex","Subject gender");
	    }
	    break;
	  case SON_SUBJ_WEIGHT: 
	    {
	      kvDefDouble(info,"weight",*(float*)rec->data);
	      kvDefString(defs,"weight","subject weight in kilos");
	    }
	    break;
	  case SON_SUBJ_HEIGHT: 
	    {
	      kvDefDouble(info,"height",*(float*)rec->data);
	      kvDefString(defs,"height","subject height in meters");
	    }
	    break;
	  case SON_SUBJ_COMMENT:
	    {
	      kvDefString(info,"subject_comment",(char*)rec->data);
	      kvDefString(defs,"subject_commend","misc subject info");
	    }
	    break;
	  }
	break;
      case SONB_DACQ_PARS:
	switch (rec->kind) {
	case SON_BLOCK_START:
	case SON_BLOCK_END:
	  break; /* avoid printing debugging info */
	default:
	  if (debug && rec->type==SONT_STRING)
	    fprintf(stderr,"Val is <%s>\n",(char*)rec->data);
	  break;
	}
	break;
      case SONB_MEAS_INFO:
	switch (rec->kind) 
	  {
	  case SON_EXPERIMENTER: 
	    {
	      kvDefString(info,"experimenter",(char*)rec->data);
	      kvDefString(defs,"experimenter","experimenter");
	    }
	    break;
	  case SON_COMMENT:
	    {
	      kvDefString(info,"comment",(char*)rec->data);
	      kvDefString(defs,"comment","meas info comment");
	    }
	    break;
	  case SON_PROJ_ID: 
	    {
	      kvDefInt(info,"project_id",*(int*)rec->data);
	      kvDefString(defs,"project_id","project identifier");
	    }
	    break;
	  case SON_PROJ_NAME:
	    {
	      kvDefString(info,"project_name",(char*)rec->data);
	      kvDefString(defs,"project_name","project name");
	    }
	    break;
	  case SON_PROJ_COMMENT:
	    {
	      kvDefString(info,"project_comment",(char*)rec->data);
	      kvDefString(defs,"project_commend","misc project info");
	    }
	    break;
	  case SON_PROJ_PERSONS:
	    {
	      kvDefString(info,"project_persons",(char*)rec->data);
	      kvDefString(defs,"project_persons","people in charge");
	    }
	    break;
	  case SON_NCHAN:
	    {
	      kvDefInt(info,"nchannels",*(int*)rec->data);
	      kvDefString(defs,"nchannels","number of channels collected");
	    }
	    break;
	  case SON_CH_INFO:
	    {
	      char buf[64];
	      sonChInfoRec* ch= (sonChInfoRec*)(rec->data);
	      
	      snprintf(buf,64,"channel.%3d",ch->scanNo);
	      kvDefString(info,buf,ch->ch_name);
	    }
	    break;
	  case SON_MEAS_DATE:
	    {
	      struct tm tmstr;
	      char buf[64];
	      time_t t= *(int*)rec->data;
	      data->unixTimeSecs= t;
	      localtime_r(&t,&tmstr);
	      snprintf(buf,sizeof(buf),"%d:%02d:%02d",
		       tmstr.tm_hour,tmstr.tm_min,tmstr.tm_sec);
	      kvDefString(info,"time",buf);
	      kvDefString(defs,"time","acquisition time");
	      snprintf(buf,sizeof(buf),"%d/%d/%d",
		       tmstr.tm_mon+1,tmstr.tm_mday,tmstr.tm_year+1900);
	      kvDefString(info,"date",buf);
	      kvDefString(defs,"date","acquisition date (mm/dd/yyyy)");
	    }
	    break;
	  case SON_SFREQ:
	    {
	      kvDefDouble(info,"sfreq",*(float*)(rec->data));
	      kvDefString(defs,"sfreq","sampling frequency (hz)");
	    }
	    break;
	  case SON_LOWPASS:
	    {
	      kvDefDouble(info,"lowpass",*(float*)(rec->data));
	      kvDefString(defs,"lowpass","analog lowpass thresh (hz)");
	    }
	    break;
	  case SON_HIGHPASS:
	    {
	      kvDefDouble(info,"highpass",*(float*)(rec->data));
	      kvDefString(defs,"highpass","analog highpass thresh (hz)");
	    }
	    break;
	  case SON_BAD_CHS:
	    /* This block apparently notes which channels are bad, but
	     * I don't have any example of what value denotes badness.
	     */
	    {
	      int i;
	      int n= sonNDataItems(rec->type,rec->size);
	      int* buf= (int*)rec->data;
	      for (i=0; i<n; i++) fprintf(stderr,"bad? %d: %d\n",i,buf[i]);
	    }
	    break;
	  case SON_DATA_PACK:
	    /* This information is redundant with the tag type for the
	     * actual data records.
	     */
	    break;
	  case SON_BLOCK_START:
	  case SON_BLOCK_END:
	    break; /* avoid printing debugging info */
	  default:
	    fprintf(stderr,"Found %s in block type %s, type %s, size %d!\n",
		    sonKindToString(rec->kind),
		    sonBlockToString(currentBlockType),
		    sonTypeToString(rec->type),rec->size);
	  }
	break;
      case SONB_RAW_DATA:
	switch (rec->kind) {
	case SON_DATA_BUFFER:
	  {
	    slist_append(data->rawBufferList,
			 allocFDBInfo( offset + SONC_DATA_OFFSET,
				       rec->size,
				       rec->type ));
	  }
	  break;
	  case SON_BLOCK_START:
	  case SON_BLOCK_END:
	    break; /* avoid printing debugging info */
	  default:
	    if (debug)
	      fprintf(stderr,"Found %s in block type %s!\n",
		      sonKindToString(rec->kind),
		      sonBlockToString(currentBlockType));
	}
	break;
      case SONB_PROCESSED_DATA:
	switch (rec->kind) {
	case SON_BLOCK_START:
	case SON_BLOCK_END:
	  break; /* avoid printing debugging info */
	default:
	  fprintf(stderr,"Found %s in block type %s!\n",
		  sonKindToString(rec->kind),
		  sonBlockToString(currentBlockType));
	}
	break;
      case SONB_EVOKED:
	switch (rec->kind) {
	case SON_BLOCK_START:
	case SON_BLOCK_END:
	  break; /* avoid printing debugging info */
	case SON_COMMENT:
	case SON_EVENT_COMMENT:
	  if (debug) fprintf(stderr,"%s in block type %s: <%s>\n",
			     sonKindToString(rec->kind),
			     sonBlockToString(currentBlockType),
			     (char*)rec->data);		  
	  break;
	default:
	  { /* Do nothing */ }
	}
	break;
      case SONB_CONTINUOUS_DATA:
	switch (rec->kind) {
	  case SON_BLOCK_START:
	  case SON_BLOCK_END:
	    break; /* avoid printing debugging info */
	default:
	  { /* Do nothing */ }
	}
	break;
      case SONB_ASPECT:
	switch (rec->kind) {
	case SON_BLOCK_START:
	case SON_BLOCK_END:
	  break; /* avoid printing debugging info */
	case SON_NAVE:
	  {
	    kvDefInt(info,"n_ave",*(int*)rec->data);
	    kvDefString(defs,"n_ave","Number of averages");
	  }
	  break;
	case SON_ASPECT_KIND:
	  {
	    kvDefString(info,"aspect",sonAspectToString(*(int*)rec->data));
	    kvDefString(defs,"aspect","how epochs were subsampled");
	  }
	  break;
	case SON_EPOCH:
	  {
	    slist_append(data->epochBufferList,
			 allocFDBInfo( offset + SONC_DATA_OFFSET,
				       rec->size,
				       rec->type ));
	  }
	  break;
	default:
	  fprintf(stderr,"Found %s in block type %s; type is %s, size %d!\n",
		  sonKindToString(rec->kind),
		  sonBlockToString(currentBlockType),
		  sonTypeToString(rec->type),rec->size);
	}
	break;
      }
    freeTag(rec);
  }

  if (fclose(f)) perror("Error closing input file");

  if (!slist_empty(blockStack))
    fprintf(stderr,"*** Unmatched block_start tag! ***\n");
  slist_destroy(blockStack,free);

  scanBlockList(self,info);
#endif /* ifdef never */

  bio_big_endian_input= old_bio_big_endian_input;
}

static void sonRead( FileHandler* self, KVHash* info,
		     long long offset, long n,
		     SRDR_Datatype datatype, void* obuf )
{
  SonData* data= (SonData*)(self->hook);

#ifdef never
  SList* blockList= NULL;
  SonDataBufferInfo* block= NULL;
  SRDR_Datatype srdrType;
  const char* dimstr= kvGetString(info,"dimstr");

  if (!strcmp(dimstr,"cbt")) {
    blockList= data->rawBufferList;
  }
  else if (!strcmp(dimstr,"bc")) {
    blockList= data->epochBufferList;
  }
  else {
    Abort("%s:%d: internal error: dim string <%s> unrecognized on read!\n",
	  __FILE__,__LINE__,dimstr);
  }

  if (offset==0) slist_totop(blockList);
  block= (SonDataBufferInfo*)slist_next(blockList);
  if (!block) 
    Abort("%s:%d: ran out of data!\n",__FILE__,__LINE__);
  srdrType= sonTypeMap(block->type);
  if (srdrType!=datatype)
    Abort("%s:%d: Son type %s does not match %s!\n",
	  __FILE__,__LINE__,sonTypeToString(block->type),
	  srdrTypeName[block->type]);
  if (n != sonNDataItems(block->type,block->size))
    Abort("%s:%d: data block n data item not what was expected! (%ld vs %ld)\n",
	  __FILE__,__LINE__,sonNDataItems(block->type,block->size),n);
  if (block->type==SONT_OLD_PACK) {
    float scale;
    float shift;
    float* fbuf= (float*)obuf;
    short* trickBuf;
    int i;
    
    baseReopen(self);
    bigfile_fseek(self->file, offset, SEEK_SET);
    scale= FRdFloat32(self->file);
    shift= FRdFloat32(self->file);
    /* We will reuse (the upper half of) obuf for input to avoid
     * allocating a new buffer.
     */
    trickBuf= ((short*)(fbuf+n)-n);
    for (i=n-1; i>=0; i--)
      fbuf[i]= (scale*trickBuf[i])+shift;
  }
  else baseRead(self, info, block->offset, n, datatype, obuf);

#endif /* ifdef never */
  int srdrType= 0;
  if (debug) 
    fprintf(stderr,
	    "Read %d values sz %d, offset %lld\n",
	    n,srdrTypeSize[srdrType],offset);
}

static int sonCompare( FileHandler* f1, FileHandler* f2 )
{
  /* Do comparison based on Unix time, secs since Unix epoch, which
   * we hopefully found stored in the MEAS_INFO.DATE field.
   */
  SonData* data1= (SonData*)(f1->hook);
  SonData* data2= (SonData*)(f2->hook);

  if (data1->unixTimeSecs<data2->unixTimeSecs) return -1;
  else if (data1->unixTimeSecs<data2->unixTimeSecs) return 1;
  else return 0;
}

static void sonDestroySelf( FileHandler* self )
{
  SonData* data= (SonData*)(self->hook);
  baseDestroySelf(self);
}

FileHandler* sonFactory(char* fname, KVHash* info)
{
  FileHandler* result= baseFactory(fname);
  SonData* data;

  result->typeName= strdup("SonDataHandler");

  result->destroySelf= sonDestroySelf;
  result->read= sonRead;
  result->processHeader= sonProcessHeader;
  result->compareThisType= sonCompare;

  if (!(data= (SonData*)malloc(sizeof(SonData))))
    Abort("%s: unable to allocate %d bytes!\n",sizeof(SonData));
  
  data->unixTimeSecs= 0;
  result->hook= data;

  return result;
}

int sonTester(const char* filename)
{
  long long offset= 0;
  FILE* f= NULL;
  int old_bio_big_endian_input;
  int match= 1;

  /* This code makes some assumptions about the sizes of things; 
   *  test them here.
   */
  assert(sizeof(short)==2);
  assert(sizeof(long)==4);
  assert(sizeof(double)==8);

  fprintf(stderr,"A long double has size %d\n",sizeof(long double));

#ifdef never
  sonTagRec* rec= NULL;

  if ((f = fopen(filename,"r"))==NULL) return 0;

  /* SON files are always bigendian */
  old_bio_big_endian_input= bio_big_endian_input;
  bio_big_endian_input= 1;

  /* If this is a SON file, the first three tags should be
   * SON_FILE_ID, SON_DIR_POINTER, and SON_FREE_LIST.
   */
  bigfile_fseek(f,offset,SEEK_SET);
  rec= readSonTag(f);
  if ((rec!=NULL) && (rec->type==SONT_ID_STRUCT) && (rec->kind==SON_FILE_ID)
      && (rec->next != SONV_NEXT_NONE)) {
    freeTag(rec);
    offset= getNextTagOffset(rec,offset);
    bigfile_fseek(f,offset,SEEK_SET);
    rec= readSonTag(f);
    if ((rec->type==SONT_INT) && (rec->kind==SON_DIR_POINTER)
	&& (rec->next != SONV_NEXT_NONE)) {
      freeTag(rec);
      offset= getNextTagOffset(rec,offset);
      bigfile_fseek(f,offset,SEEK_SET);
      rec= readSonTag(f);
      if ((rec->type==SONT_INT) && (rec->kind==SON_FREE_LIST)) {
	/* Do nothing; freeing happens below */
      }
      else match= 0;
    }
    else match= 0;
  }
  else match= 0;
  if (rec) freeTag(rec);

  /* Restore bio state */
  bio_big_endian_input= old_bio_big_endian_input;

  if (fclose(f)) {
    perror("Error closing file");
    return 0;
  }
#endif

  return match;
}

#else /* ifdef USE_SON */

FileHandler* sonFactory(char* fname, KVHash* info)
{
  FileHandler* result= baseFactory("NotARealFile");
  return result;
}

int sonTester(const char* filename)
{
  return 0;
}

#endif

