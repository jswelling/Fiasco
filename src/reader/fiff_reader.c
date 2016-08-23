/************************************************************
 *                                                          *
 *  fiff_reader.c                                           *
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
#include "mri.h"
#if (SGI64 || SGI5 || SGIMP)
#include <bstring.h>
#endif

#include "bio.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"
#include "smartreader.h"

static char rcsid[] = "$Id: fiff_reader.c,v 1.9 2006/11/27 19:54:31 welling Exp $";

/* Notes-
 * -maybe don't bother to read the data during scanning phase?
 * -We could add the kind, range, and calibration for each channel as
 *  additional per-channel key-value pairs, but it seems a bit much.
 * -FIFF_DACQ_PARS in FIFFB_DACQ_PARS is a huge ascii table, but I
 *  don't know what to do with it.
 */

#ifdef USE_FIFF

#include "fiff_types.h"
#include "fiff_file.h"

/* A reader for FIFF files */
typedef struct fiff_data_struct {
  SList* rawBufferList;
  SList* epochBufferList;
  SList* continuousBufferList;
  time_t unixTimeSecs;
} FiffData;

typedef struct fiff_block_info {
  int type;
  int parentType;
  long long start_offset;
} FiffBlockInfo;

typedef struct fiff_data_buffer_info {
  long long offset;
  long size; /* in bytes */
  int type;
} FiffDataBufferInfo;

static FiffDataBufferInfo* allocFDBInfo( long long offset, long size, 
					 int type )
{
  FiffDataBufferInfo* result= NULL;
  if (!(result=(FiffDataBufferInfo*)malloc(sizeof(FiffDataBufferInfo))))
    Abort("%s:%d: unable to allocate %d bytes!\n",
	  __FILE__,__LINE__,sizeof(FiffDataBufferInfo));
  result->offset= offset;
  result->size= size;
  result->type= type;
  return result;
}

static void freeFDBInfo( void* v )
{
  free(v);
}

static void fiffDestroySelf( FileHandler* self )
{
  FiffData* data= (FiffData*)(self->hook);
  slist_destroy(data->rawBufferList,freeFDBInfo);
  slist_destroy(data->epochBufferList,freeFDBInfo);
  slist_destroy(data->continuousBufferList,freeFDBInfo);
  baseDestroySelf(self);
}

static const char* fiffTypeToString( int type )
{
  switch (type) {
  case FIFFT_VOID: return "VOID";
  case FIFFT_BYTE: return "BYTE";
  case FIFFT_SHORT: return "SHORT";
  case FIFFT_INT: return "INT";
  case FIFFT_FLOAT: return "FLOAT";
  case FIFFT_DOUBLE: return "DOUBLE";
  case FIFFT_JULIAN: return "JULIAN"; 
  case FIFFT_USHORT: return "USHORT";
  case FIFFT_UINT: return "UINT";
  case FIFFT_STRING: return "STRING";
    /* case FIFFT_ASCII: return "ASCII"; */ /* redundant */
  case FIFFT_DAU_PACK13: return "DAU_PACK13";
  case FIFFT_DAU_PACK14: return "DAU_PACK14";
  case FIFFT_DAU_PACK16: return "DAU_PACK16";
  case FIFFT_OLD_PACK: return "OLD_PACK";
  case FIFFT_CH_INFO_STRUCT: return "CH_INFO_STRUCT";
  case FIFFT_ID_STRUCT: return "ID_STRUCT";
  case FIFFT_DIR_ENTRY_STRUCT: return "DIR_ENTRY_STRUCT";
  case FIFFT_DIG_POINT_STRUCT: return "DIG_POINT_STRUCT";
  case FIFFT_CH_POS_STRUCT: return "CH_POS_STRUCT";
  case FIFFT_COORD_TRANS_STRUCT: return "COORD_TRANS_STRUCT";
  case FIFFT_DIG_STRING_STRUCT: return "DIG_STRING_STRUCT";
  default: return "***UNKNOWN***";
  }
}

static const char* fiffBlockToString(int kind)
{
  if (kind==-1) return "TOP-OF-STACK";
  else 
    switch (kind) {
    case FIFFB_ROOT: return "ROOT";
    case FIFFB_MEAS: return "MEAS";
    case FIFFB_MEAS_INFO: return "MEAS_INFO";
    case FIFFB_RAW_DATA: return "RAW_DATA";
    case FIFFB_PROCESSED_DATA: return "PROCESSED_DATA";
    case FIFFB_EVOKED: return "EVOKED";
      /* case FIFFB_MCG_AVE: return "MCG_AVE"; */ /* redundant */
    case FIFFB_ASPECT: return "ASPECT";
    case FIFFB_SUBJECT: return "SUBJECT";
    case FIFFB_ISOTRAK: return "ISOTRAK";
    case FIFFB_HPI_MEAS: return "HPI_MEAS";
    case FIFFB_HPI_RESULT: return "HPI_RESULT";
    case FIFFB_HPI_COIL: return "HPI_COIL";
    case FIFFB_PROJECT: return "PROJECT";
    case FIFFB_CONTINUOUS_DATA: return "CONTINUOUS_DATA";
    case FIFFB_VOID: return "VOID";
    case FIFFB_EVENTS: return "EVENTS";
    case FIFFB_INDEX: return "INDEX";
    case FIFFB_DACQ_PARS: return "DACQ_PARS";
    case FIFFB_REF: return "REF";
    case FIFFB_SMSH_RAW_DATA: return "SMSH_RAW_DATA";
    case FIFFB_SMSH_ASPECT: return "SMSH_ASPECT";
    case FIFFB_MRI: return "MRI";
    case FIFFB_MRI_SET: return "MRI_SET";
    case FIFFB_MRI_SLICE: return "MRI_SLICE";
    case FIFFB_MRI_SCENERY: return "MRI_SCENERY";
    case FIFFB_MRI_SCENE: return "MRI_SCENE";
    case FIFFB_MRI_SEG: return "MRI_SEG";
    case FIFFB_MRI_SEG_REGION: return "MRI_SEG_REGION";
    case FIFFB_SPHERE: return "SPHERE";
    case FIFFB_BEM: return "BEM";
    case FIFFB_BEM_SURF: return "BEM_SURF";
    case FIFFB_XFIT_PROJ: return "XFIT_PROJ";
    case FIFFB_XFIT_PROJ_ITEM: return "XFIT_PROJ_ITEM";
    case FIFFB_XFIT_AUX: return "XFIT_AUX";
    case FIFFB_VOL_INFO: return "VOL_INFO";
    case FIFFB_DATA_CORRECTION: return "DATA_CORRECTION";
    case FIFFB_CHANNEL_DECOUPLER: return "CHANNEL_DECOUPLER";
    case FIFFB_SSS_INFO: return "SSS_INFO";
    case FIFFB_SSS_CAL_ADJUST: return "SSS_CAL_ADJUST";
    case FIFFB_SSS_ST_INFO: return "SSS_ST_INFO";
    case FIFFB_SSS_BASES: return "SSS_BASES";
    case FIFFB_SMARTSHIELD: return "SMARTSHIELD";
    case FIFFB_PROCESSING_HISTORY: return "PROCESSING_HISTORY";
    case FIFFB_PROCESSING_RECORD: return "PROCESSING_RECORD";
    default: return "***UNKNOWN***";
    }
}

static const char* fiffKindToString(int v)
{
  switch (v) {
  case FIFF_NEW_FILE: return "NEW_FILE";
  case FIFF_CLOSE_FILE: return "CLOSE_FILE";
  case FIFF_DISCARD_FILE: return "DISCARD_FILE";
  case FIFF_ERROR_MESSAGE: return "ERROR_MESSAGE";
  case FIFF_SUSPEND_READING: return "SUSPEND_READING";
  case FIFF_FATAL_ERROR_MESSAGE: return "FATAL_ERROR_MESSAGE";
  case FIFF_CONNECTION_CHECK: return "CONNECTION_CHECK";
  case FIFF_SUSPEND_FILING: return "SUSPEND_FILING";
  case FIFF_RESUME_FILING: return "RESUME_FILING";
  case FIFF_RAW_PREBASE: return "RAW_PREBASE";
  case FIFF_RAW_PICK_LIST: return "RAW_PICK_LIST";
  case FIFF_ECHO: return "ECHO";
  case FIFF_RESUME_READING: return "RESUME_READING";
  case FIFF_DACQ_SYSTEM_TYPE: return "DACQ_SYSTEM_TYPE";
  case FIFF_SELECT_RAW_CH: return "SELECT_RAW_CH";
  case FIFF_PLAYBACK_MODE: return "PLAYBACK_MODE";
  case FIFF_CONTINUE_FILE: return "CONTINUE_FILE";
  case FIFF_JITTER_MAX: return "JITTER_MAX";
  case FIFF_DECIMATION_FACTOR: return "DECIMATION_FACTOR";
#ifdef _DATA_SERVER
  case FIFF_MEM_DATA_BUFFER: return "MEM_DATA_BUFFER";
#endif
  case FIFF_FILE_ID: return "FILE_ID";
  case FIFF_DIR_POINTER: return "DIR_POINTER";
  case FIFF_DIR: return "DIR";
  case FIFF_BLOCK_ID: return "BLOCK_ID";
  case FIFF_BLOCK_START: return "BLOCK_START";
  case FIFF_BLOCK_END: return "BLOCK_END";
  case FIFF_FREE_LIST: return "FREE_LIST";
  case FIFF_FREE_BLOCK: return "FREE_BLOCK";
  case FIFF_NOP: return "NOP";
  case FIFF_PARENT_FILE_ID: return "PARENT_FILE_ID";
  case FIFF_PARENT_BLOCK_ID: return "PARENT_BLOCK_ID";
  case FIFF_BLOCK_NAME: return "BLOCK_NAME";
  case FIFF_BLOCK_VERSION: return "BLOCK_VERSION";
  case FIFF_CREATOR: return "CREATOR";
  case FIFF_MODIFIER: return "MODIFIER";
  case FIFF_REF_ROLE: return "REF_ROLE";
  case FIFF_REF_FILE_ID: return "REF_FILE_ID";
  case FIFF_REF_FILE_NUM: return "REF_FILE_NUM";
  case FIFF_REF_FILE_NAME: return "REF_FILE_NAME";
  case FIFF_REF_BLOCK_ID: return "REF_BLOCK_ID";
  case FIFF_DACQ_PARS: return "DACQ_PARS";
  case FIFF_DACQ_STIM: return "DACQ_STIM";
  case FIFF_NCHAN: return "NCHAN";
  case FIFF_SFREQ: return "SFREQ";
  case FIFF_DATA_PACK: return "DATA_PACK";
  case FIFF_CH_INFO: return "CH_INFO";
  case FIFF_MEAS_DATE: return "MEAS_DATE";
  case FIFF_SUBJECT: return "SUBJECT";
  case FIFF_COMMENT: return "COMMENT";
  case FIFF_NAVE: return "NAVE";
  case FIFF_FIRST_SAMPLE: return "FIRST_SAMPLE";
  case FIFF_LAST_SAMPLE: return "LAST_SAMPLE";
  case FIFF_ASPECT_KIND: return "ASPECT_KIND";
  case FIFF_REF_EVENT: return "REF_EVENT";
  case FIFF_EXPERIMENTER: return "EXPERIMENTER";
  case FIFF_DIG_POINT: return "DIG_POINT";
  case FIFF_CH_POS_VEC: return "CH_POS_VEC";
  case FIFF_HPI_SLOPES: return "HPI_SLOPES";
  case FIFF_HPI_NCOIL: return "HPI_NCOIL";
  case FIFF_REQ_EVENT: return "REQ_EVENT";
  case FIFF_REQ_LIMIT: return "REQ_LIMIT";
  case FIFF_LOWPASS: return "LOWPASS";
  case FIFF_BAD_CHS: return "BAD_CHS";
  case FIFF_ARTEF_REMOVAL: return "ARTEF_REMOVAL";
  case FIFF_COORD_TRANS: return "COORD_TRANS";
  case FIFF_HIGHPASS: return "HIGHPASS";
  case FIFF_CH_CALS_VEC: return "CH_CALS_VEC";
  case FIFF_HPI_BAD_CHS: return "HPI_BAD_CHS";
  case FIFF_HPI_CORR_COEFF: return "HPI_CORR_COEFF";
  case FIFF_EVENT_COMMENT: return "EVENT_COMMENT";
  case FIFF_NO_SAMPLES: return "NO_SAMPLES";
  case FIFF_FIRST_TIME: return "FIRST_TIME";
  case FIFF_SUBAVE_SIZE: return "SUBAVE_SIZE";
  case FIFF_SUBAVE_FIRST: return "SUBAVE_FIRST";
  case FIFF_NAME: return "NAME";
    /* case FIFF_DESCRIPTION: return "DESCRIPTION"; */ /* redundant */
  case FIFF_DIG_STRING: return "DIG_STRING";
  case FIFF_LINE_FREQ: return "LINE_FREQ";
  case FIFF_HPI_COIL_FREQ: return "HPI_COIL_FREQ";
  case FIFF_HPI_COIL_MOMENTS: return "HPI_COIL_MOMENTS";
  case FIFF_HPI_FIT_GOODNESS: return "HPI_FIT_GOODNESS";
  case FIFF_HPI_FIT_ACCEPT: return "HPI_FIT_ACCEPT";
  case FIFF_HPI_FIT_GOOD_LIMIT: return "HPI_FIT_GOOD_LIMIT";
  case FIFF_HPI_FIT_DIST_LIMIT: return "HPI_FIT_DIST_LIMIT";
  case FIFF_HPI_COIL_NO: return "HPI_COIL_NO";
  case FIFF_HPI_COILS_USED: return "HPI_COILS_USED";
  case FIFF_HPI_DIGITIZATION_ORDER: return "HPI_DIGITIZATION_ORDER";
  case FIFF_CH_SCAN_NO: return "CH_SCAN_NO";
  case FIFF_CH_LOGICAL_NO: return "CH_LOGICAL_NO";
  case FIFF_CH_KIND: return "CH_KIND";
  case FIFF_CH_RANGE: return "CH_RANGE";
  case FIFF_CH_CAL: return "CH_CAL";
  case FIFF_CH_POS: return "CH_POS";
  case FIFF_CH_UNIT: return "CH_UNIT";
  case FIFF_CH_UNIT_MUL: return "CH_UNIT_MUL";
  case FIFF_CH_DACQ_NAME: return "CH_DACQ_NAME";
  case FIFF_SSS_FRAME: return "SSS_FRAME";
  case FIFF_SSS_JOB: return "SSS_JOB";
  case FIFF_SSS_ORIGIN: return "SSS_ORIGIN";
  case FIFF_SSS_ORD_IN: return "SSS_ORD_IN";
  case FIFF_SSS_ORD_OUT: return "SSS_ORD_OUT";
  case FIFF_SSS_NMAG: return "SSS_NMAG";
  case FIFF_SSS_COMPONENTS: return "SSS_COMPONENTS";
  case FIFF_SSS_CAL_CHANS: return "SSS_CAL_CHANS";
  case FIFF_SSS_CAL_CORRS: return "SSS_CAL_CORRS";
  case FIFF_SSS_ST_CORR: return "SSS_ST_CORR";
  case FIFF_SSS_BASE_IN: return "SSS_BASE_IN";
  case FIFF_SSS_BASE_OUT: return "SSS_BASE_OUT";
  case FIFF_SSS_BASE_VIRT: return "SSS_BASE_VIRT";
  case FIFF_SSS_NORM: return "SSS_NORM";
  case FIFF_DATA_BUFFER: return "DATA_BUFFER";
  case FIFF_DATA_SKIP: return "DATA_SKIP";
  case FIFF_EPOCH: return "EPOCH";
  case FIFF_DATA_SKIP_SAMP: return "DATA_SKIP_SAMP";
  case FIFF_SUBJ_ID: return "SUBJ_ID";
  case FIFF_SUBJ_FIRST_NAME: return "SUBJ_FIRST_NAME";
  case FIFF_SUBJ_MIDDLE_NAME: return "SUBJ_MIDDLE_NAME";
  case FIFF_SUBJ_LAST_NAME: return "SUBJ_LAST_NAME";
  case FIFF_SUBJ_BIRTH_DAY: return "SUBJ_BIRTH_DAY";
  case FIFF_SUBJ_SEX: return "SUBJ_SEX";
  case FIFF_SUBJ_HAND: return "SUBJ_HAND";
  case FIFF_SUBJ_WEIGHT: return "SUBJ_WEIGHT";
  case FIFF_SUBJ_HEIGHT: return "SUBJ_HEIGHT";
  case FIFF_SUBJ_COMMENT: return "SUBJ_COMMENT";
  case FIFF_SUBJ_HIS_ID: return "SUBJ_HIS_ID";
  case FIFF_PROJ_ID: return "PROJ_ID";
  case FIFF_PROJ_NAME: return "PROJ_NAME";
  case FIFF_PROJ_AIM: return "PROJ_AIM";
  case FIFF_PROJ_PERSONS: return "PROJ_PERSONS";
  case FIFF_PROJ_COMMENT: return "PROJ_COMMENT";
  case FIFF_EVENT_CHANNELS: return "EVENT_CHANNELS";
  case FIFF_EVENT_LIST: return "EVENT_LIST";
  case FIFF_SQUID_BIAS: return "SQUID_BIAS";
  case FIFF_SQUID_OFFSET: return "SQUID_OFFSET";
  case FIFF_SQUID_GATE: return "SQUID_GATE";
  case FIFF_DECOUPLER_MATRIX: return "DECOUPLER_MATRIX";
  case FIFF_SPARSE_CH_NAME_LIST: return "SPARSE_CH_NAME_LIST";
  case FIFF_REF_PATH: return "REF_PATH";
    /* case FIFF_MRI_SOURCE_PATH: return "MRI_SOURCE_PATH"; */ /* redundant */
  case FIFF_MRI_SOURCE_FORMAT: return "MRI_SOURCE_FORMAT";
  case FIFF_MRI_PIXEL_ENCODING: return "MRI_PIXEL_ENCODING";
  case FIFF_MRI_PIXEL_DATA_OFFSET: return "MRI_PIXEL_DATA_OFFSET";
  case FIFF_MRI_PIXEL_SCALE: return "MRI_PIXEL_SCALE";
  case FIFF_MRI_PIXEL_DATA: return "MRI_PIXEL_DATA";
  case FIFF_MRI_PIXEL_OVERLAY_ENCODING: return "MRI_PIXEL_OVERLAY_ENCODING";
  case FIFF_MRI_PIXEL_OVERLAY_DATA: return "MRI_PIXEL_OVERLAY_DATA";
  case FIFF_MRI_BOUNDING_BOX: return "MRI_BOUNDING_BOX";
  case FIFF_MRI_WIDTH: return "MRI_WIDTH";
  case FIFF_MRI_WIDTH_M: return "MRI_WIDTH_M";
  case FIFF_MRI_HEIGHT: return "MRI_HEIGHT";
  case FIFF_MRI_HEIGHT_M: return "MRI_HEIGHT_M";
  case FIFF_MRI_DEPTH: return "MRI_DEPTH";
  case FIFF_MRI_DEPTH_M: return "MRI_DEPTH_M";
  case FIFF_MRI_THICKNESS: return "MRI_THICKNESS";
  case FIFF_MRI_SCENE_AIM: return "MRI_SCENE_AIM";
  case FIFF_MRI_ORIG_SOURCE_PATH: return "MRI_ORIG_SOURCE_PATH";
  case FIFF_MRI_ORIG_SOURCE_FORMAT: return "MRI_ORIG_SOURCE_FORMAT";
  case FIFF_MRI_ORIG_PIXEL_ENCODING: return "MRI_ORIG_PIXEL_ENCODING";
  case FIFF_MRI_ORIG_PIXEL_DATA_OFFSET: return "MRI_ORIG_PIXEL_DATA_OFFSET";
  case FIFF_MRI_VOXEL_DATA: return "MRI_VOXEL_DATA";
  case FIFF_MRI_VOXEL_ENCODING: return "MRI_VOXEL_ENCODING";
  case FIFF_MRI_MRILAB_SETUP: return "MRI_MRILAB_SETUP";
  case FIFF_MRI_SEG_REGION_ID: return "MRI_SEG_REGION_ID";
  case FIFF_CONDUCTOR_MODEL_KIND: return "CONDUCTOR_MODEL_KIND";
  case FIFF_SPHERE_ORIGIN: return "SPHERE_ORIGIN";
  case FIFF_SPHERE_COORD_FRAME: return "SPHERE_COORD_FRAME";
  case FIFF_SPHERE_LAYERS: return "SPHERE_LAYERS";
  case FIFF_BEM_SURF_ID: return "BEM_SURF_ID";
  case FIFF_BEM_SURF_NAME: return "BEM_SURF_NAME";
  case FIFF_BEM_SURF_NNODE: return "BEM_SURF_NNODE";
  case FIFF_BEM_SURF_NTRI: return "BEM_SURF_NTRI";
  case FIFF_BEM_SURF_NODES: return "BEM_SURF_NODES";
  case FIFF_BEM_SURF_TRIANGLES: return "BEM_SURF_TRIANGLES";
  case FIFF_BEM_SURF_NORMALS: return "BEM_SURF_NORMALS";
  case FIFF_BEM_SURF_CURVS: return "BEM_SURF_CURVS";
  case FIFF_BEM_SURF_CURV_VALUES: return "BEM_SURF_CURV_VALUES";
  case FIFF_BEM_POT_SOLUTION: return "BEM_POT_SOLUTION";
  case FIFF_BEM_APPROX: return "BEM_APPROX";
  case FIFF_BEM_COORD_FRAME: return "BEM_COORD_FRAME";
  case FIFF_BEM_SIGMA: return "BEM_SIGMA";
  case FIFF_SOURCE_DIPOLE: return "SOURCE_DIPOLE";
  case FIFF_XFIT_LEAD_PRODUCTS: return "XFIT_LEAD_PRODUCTS";
  case FIFF_XFIT_MAP_PRODUCTS: return "XFIT_MAP_PRODUCTS";
  case FIFF_XFIT_GRAD_MAP_PRODUCTS: return "XFIT_GRAD_MAP_PRODUCTS";
  case FIFF_XFIT_VOL_INTEGRATION: return "XFIT_VOL_INTEGRATION";
  case FIFF_XFIT_INTEGRATION_RADIUS: return "XFIT_INTEGRATION_RADIUS";
  case FIFF_XFIT_CONDUCTOR_MODEL_NAME: return "XFIT_CONDUCTOR_MODEL_NAME";
  case FIFF_XFIT_CONDUCTOR_MODEL_TRANS_NAME: return "XFIT_CONDUCTOR_MODEL_TRANS_NAME";
  case FIFF_PROJ_ITEM_KIND: return "PROJ_ITEM_KIND";
  case FIFF_PROJ_ITEM_TIME: return "PROJ_ITEM_TIME";
    /* case FIFF_PROJ_ITEM_DIPOLE: return "PROJ_ITEM_DIPOLE"; */ /* redundant */
  case FIFF_PROJ_ITEM_IGN_CHS: return "PROJ_ITEM_IGN_CHS";
  case FIFF_PROJ_ITEM_NVEC: return "PROJ_ITEM_NVEC";
  case FIFF_PROJ_ITEM_VECTORS: return "PROJ_ITEM_VECTORS";
    /* case FIFF_PROJ_ITEM_COMMENT: return "PROJ_ITEM_COMMENT"; */ /* redundant */
    /* case FIFF_PROJ_ITEM_DESCRIPTION: return "PROJ_ITEM_DESCRIPTION"; */ /* redundant */
  case FIFF_PROJ_ITEM_DEFINITION: return "PROJ_ITEM_DEFINITION";
    /* case FIFF_PROJ_ITEM_CH_NAME_LIST: return "PROJ_ITEM_CH_NAME_LIST"; */ /* redundant */
    /*  case FIFF_XFIT_PROJ_ITEM_KIND: return "XFIT_PROJ_ITEM_KIND"; */ /* redundant */
    /* case FIFF_XFIT_PROJ_ITEM_TIME: return "XFIT_PROJ_ITEM_TIME"; */ /* redundant */
    /* case FIFF_XFIT_PROJ_ITEM_DIPOLE: return "XFIT_PROJ_ITEM_DIPOLE"; */ /* redundant */
    /* case FIFF_XFIT_PROJ_ITEM_IGN_CHS: return "XFIT_PROJ_ITEM_IGN_CHS"; */ /*redundant */
    /* case FIFF_XFIT_PROJ_ITEM_NVEC: return "XFIT_PROJ_ITEM_NVEC"; */ /*redundant*/
    /* case FIFF_XFIT_PROJ_ITEM_VECTORS: return "XFIT_PROJ_ITEM_VECTORS"; */ /*redundant*/
    /* case FIFF_XFIT_PROJ_ITEM_COMMENT: return "XFIT_PROJ_ITEM_COMMENT"; */ /*redundant*/
  case FIFF_XPLOTTER_LAYOUT: return "XPLOTTER_LAYOUT";
  case FIFF_VOL_ID: return "VOL_ID";
  case FIFF_VOL_NAME: return "VOL_NAME";
  case FIFF_VOL_OWNER_ID: return "VOL_OWNER_ID";
  case FIFF_VOL_OWNER_NAME: return "VOL_OWNER_NAME";
  case FIFF_VOL_OWNER_REAL_NAME: return "VOL_OWNER_REAL_NAME";
  case FIFF_VOL_TYPE: return "VOL_TYPE";
  case FIFF_VOL_HOST: return "VOL_HOST";
  case FIFF_VOL_REAL_ROOT: return "VOL_REAL_ROOT";
  case FIFF_VOL_SYMBOLIC_ROOT: return "VOL_SYMBOLIC_ROOT";
  case FIFF_VOL_MOUNT_POINT: return "VOL_MOUNT_POINT";
  case FIFF_VOL_BLOCKS: return "VOL_BLOCKS";
  case FIFF_VOL_FREE_BLOCKS: return "VOL_FREE_BLOCKS";
  case FIFF_VOL_AVAIL_BLOCKS: return "VOL_AVAIL_BLOCKS";
  case FIFF_VOL_BLOCK_SIZE: return "VOL_BLOCK_SIZE";
  case FIFF_VOL_DIRECTORY: return "VOL_DIRECTORY";
  case FIFF_INDEX_KIND: return "INDEX_KIND";
  case FIFF_INDEX: return "INDEX";
  default: return "*** Unknown FIFF kind ***";
  }
}

static const char* fiffAspectToString( int kind )
{
  switch (kind) {
    /* Normal average of epochs */
  case FIFFV_ASPECT_AVERAGE: return "AVERAGE";
    /* Std. error of mean */
  case FIFFV_ASPECT_STD_ERR: return "STD_ERR";
    /* Single epoch cut out from the continuous data */
  case FIFFV_ASPECT_SINGLE: return "SINGLE";
    /* ??? */
  case FIFFV_ASPECT_SUBAVERAGE: return "SUBAVERAGE";
    /* Alternating subaverage */
  case FIFFV_ASPECT_ALTAVERAGE: return "ALTAVERAGE";
    /* A sample cut out by graph */
  case FIFFV_ASPECT_SAMPLE: return "SAMPLE";
    /* Power density spectrum */
  case FIFFV_ASPECT_POWER_DENSITY: return "POWER_DENSITY";
    /* Dipole amplitude curve */
  case FIFFV_ASPECT_DIPOLE_WAVE: return "DIPOLE_WAVE";
  default: return "*** Unknown FIFF aspect ***";
  }
}

static fiffTagRec* readFiffTag(FILE* f)
{
  fiffTagRec* result= NULL;

  if (!(result=(fiffTagRec*)malloc(sizeof(fiffTagRec))))
    Abort("%s:%d: unable to allocate %d bytes!\n",__FILE__,__LINE__,
	  sizeof(fiffTagRec));

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
  if (!(result->data=(fiff_data_t*)malloc(result->size)))
    Abort("%s:%d: unable to allocate %d bytes!\n",__FILE__,__LINE__,
	  result->size);
  if (debug)
    fprintf(stderr,"kind %s, type %s, size %d\n",
	    fiffKindToString(result->kind),
	    fiffTypeToString(result->type),
	    result->size);
  switch (result->type) {
  case FIFFT_VOID: 
    {
      FRdUInt8Array(f,(unsigned char*)result->data,result->size/1);
    }
    break;
  case FIFFT_BYTE: 
    {
      FRdUInt8Array(f,(unsigned char*)result->data,result->size/1);
    }
    break;
  case FIFFT_SHORT: 
    {
      FRdInt16Array(f,(short*)result->data,result->size/2);
    }
    break;
  case FIFFT_INT: 
    {
      FRdInt32Array(f,(int*)result->data,result->size/4);
    }
    break;
  case FIFFT_FLOAT: 
    {
      FRdFloat32Array(f,(float*)result->data,result->size/4);
    }
    break;
  case FIFFT_DOUBLE: 
    {
      FRdFloat64Array(f,(double*)result->data,result->size/8);
    }
    break;
  case FIFFT_JULIAN: 
    {
      FRdInt32Array(f,(int*)result->data,result->size/4);
    }
    break;
  case FIFFT_USHORT: 
    {
      FRdInt16Array(f,(short*)result->data,result->size/2);
    }
    break;
  case FIFFT_UINT: 
    {
      FRdInt32Array(f,(int*)result->data,result->size/4);
    }
    break;
  case FIFFT_STRING: 
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
  case FIFFT_DAU_PACK13:
    {
      FRdInt16Array(f,(short*)result->data,result->size/2);
    }
    break;
  case FIFFT_DAU_PACK14:
    {
      FRdInt16Array(f,(short*)result->data,result->size/2);
    }
    break;
  case FIFFT_DAU_PACK16:
    {
      FRdInt16Array(f,(short*)result->data,result->size/2);
    }
    break;
  case FIFFT_OLD_PACK:
    {
      int nsamp= (result->size - 8)/2;
      short* shortPtr= (short*)(((float*)result->data)+2);
      FRdFloat32Array(f,(float*)result->data,2);
      FRdInt16Array(f,shortPtr,nsamp);
    }
    break;
  case FIFFT_CH_INFO_STRUCT:
    {
      fiffChInfoRec* rec= (fiffChInfoRec*)result->data;
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
  case FIFFT_ID_STRUCT:
    {
      fiffIdRec* rec= (fiffIdRec*)result->data;
      rec->version= FRdInt32(f);
      FRdInt32Array(f,rec->machid,2);
      rec->time.secs= FRdInt32(f);
      rec->time.usecs= FRdInt32(f);
    }
    break;
  case FIFFT_DIR_ENTRY_STRUCT:
    {
      fiffDirEntryRec* rec= (fiffDirEntryRec*)result->data;
      rec->kind= FRdInt32(f);
      rec->type= FRdInt32(f);
      rec->size= FRdInt32(f);
      rec->pos= FRdInt32(f);
    }
    break;
  case FIFFT_DIG_POINT_STRUCT:
    {
      fiffDigPointRec* rec= (fiffDigPointRec*)result->data;
      rec->kind= FRdInt32(f);
      rec->ident= FRdInt32(f);
      FRdFloat32Array(f,rec->r,3);
    }
    break;
  case FIFFT_CH_POS_STRUCT:
    {
      int n= result->size/44;
      int i;
      fiffChPosRec* rec= (fiffChPosRec*)result->data;
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
  case FIFFT_COORD_TRANS_STRUCT:
    {
      int n= result->size/80;
      int i;
      fiffCoordTransRec* rec= (fiffCoordTransRec*)result->data;
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
  case FIFFT_DIG_STRING_STRUCT:
    {
      float* buf= NULL;
      float** bufp= NULL;
      fiffDigStringRec* dsr= NULL;
      int i;

      if (!(dsr=(fiffDigStringRec*)malloc(sizeof(fiffDigStringRec))))
	Abort("%s:%d: unable to allocate %d bytes!\n",
	      __FILE__,__LINE__,sizeof(fiffDigStringRec));
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

static SRDR_Datatype fiffTypeMap( int type )
{
  switch (type) {
  case FIFFT_BYTE: return SRDR_UINT8;
  case FIFFT_SHORT: return SRDR_INT16;
  case FIFFT_INT: return SRDR_INT32;
  case FIFFT_FLOAT: return SRDR_FLOAT32;
  case FIFFT_DOUBLE: return SRDR_FLOAT64;
  case FIFFT_USHORT: return SRDR_UINT16;
  case FIFFT_DAU_PACK13: return SRDR_INT16;
  case FIFFT_DAU_PACK14: return SRDR_INT16;
  case FIFFT_DAU_PACK16: return SRDR_INT16;
  case FIFFT_OLD_PACK: return SRDR_FLOAT32;
  default: Abort("%s:%d: FIFF type %s is untranslatable!\n",
		 __FILE__,__LINE__,fiffTypeToString(type));
  }
}

static int fiffNDataItems( int type, long nbytes )
{
  switch (type) {
  case FIFFT_BYTE: return nbytes;
  case FIFFT_SHORT: return nbytes/2;
  case FIFFT_INT: return nbytes/4;
  case FIFFT_FLOAT: return nbytes/4;
  case FIFFT_DOUBLE: return nbytes/8;
  case FIFFT_USHORT: return nbytes/2;
  case FIFFT_DAU_PACK13: return nbytes/2;
  case FIFFT_DAU_PACK14: return nbytes/2;
  case FIFFT_DAU_PACK16: return nbytes/2;
  case FIFFT_OLD_PACK: return (nbytes-8)/2;
  default: Abort("%s:%d: internal error: no counting method for type %s!\n",
		 __FILE__,__LINE__,fiffTypeToString(type));
  }
}

static void freeTag( fiffTagRec* rec )
{
  if (rec->data) {
    if (rec->type==FIFFT_DIG_STRING_STRUCT) {
      fiffDigStringRec* d= (fiffDigStringRec*)(rec->data);
      if (d->rr) {
	if (d->rr[0]) free(d->rr[0]);
	free(d->rr);
      }
    }
    free(rec->data);
  }
  free(rec);
}

static long long getNextTagOffset( const fiffTagRec* rec, long long offset )
{
  if (rec->next!=FIFFV_NEXT_SEQ) {
    if (rec->next==FIFFV_NEXT_NONE) {
      Abort("%s:%d: internal error: ran off end of linked list!\n",
	    __FILE__,__LINE__);
    }
    else {
      offset= rec->next;
      if (debug) fprintf(stderr,"Next sends us to %lld\n",offset);
    }
  }
  else {
    offset += 4*sizeof(fiff_int_t)+rec->size;
    if (debug) fprintf(stderr,"Following tag in order: %lld\n",offset);
  }
  return offset;
}

static void checkBlockListConsistency( SList* blockList )
{
  FiffDataBufferInfo* firstBlock= NULL;
  long blockNum= 0;

  slist_totop(blockList);
  while (!slist_atend(blockList)) {
    FiffDataBufferInfo* block= (FiffDataBufferInfo*)slist_next(blockList);
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
  FiffData* data= (FiffData*)(self->hook);
  KVHash* defs= kvGetHash(info,"definitions");
  SList* rawList= data->rawBufferList;
  SList* epochList= data->epochBufferList;
  SList* continuousList= data->continuousBufferList;
  SRDR_Datatype srdrType;
  int nRaw= slist_count(rawList);
  int nEpoch= slist_count(epochList);
  int nContinuous= slist_count(continuousList);
  int nTypesPresent= (nRaw?1:0) + (nEpoch?1:0) + (nContinuous?1:0);

  if (debug)
    fprintf(stderr,"Raw blocks: %d  Epoch blocks: %d\n",nRaw,nEpoch);

  if (nTypesPresent>1)
  if (nRaw && nEpoch)
    Abort("%s:%d: file %s contains more than one of raw, epoch and continuous data- not supported!\n",
	  __FILE__,__LINE__,self->fileName);

  if (nRaw) {
    FiffDataBufferInfo* firstBlock= NULL;
    checkBlockListConsistency(rawList);
    slist_totop(rawList);
    firstBlock= (FiffDataBufferInfo*)slist_get(rawList);
    kvDefString(info,"dimstr","cbt");
    kvDefInt(info,"dc",kvGetInt(info,"nchannels"));
    srdrType= fiffTypeMap(firstBlock->type);
    kvDefInt(info,"datatype_in",srdrType);
    kvDefInt(info,"handler_datatype_out",srdrType);
    kvDefInt(info,"db",
	     fiffNDataItems(firstBlock->type,firstBlock->size)/kvGetInt(info,"dc"));
    kvDefLong(info,"skip.b",0);
    kvDefInt(info,"dt",nRaw);
  }
  else if (nEpoch) {
    FiffDataBufferInfo* firstBlock= NULL;
    if (nEpoch != kvGetInt(info,"nchannels"))
      Abort("%s:%d: number of epoch blocks (%d) does not match number of channels!\n",
	    __FILE__,__LINE__,nEpoch);
    checkBlockListConsistency(epochList);
    slist_totop(epochList);
    firstBlock= (FiffDataBufferInfo*)slist_get(epochList);
    kvDefString(info,"dimstr","bc");
    kvDefInt(info,"dc",nEpoch);
    srdrType= fiffTypeMap(firstBlock->type);
    kvDefInt(info,"datatype_in",srdrType);
    kvDefInt(info,"handler_datatype_out",srdrType);
    kvDefInt(info,"db",fiffNDataItems(firstBlock->type,firstBlock->size));
    kvDefLong(info,"skip.b",0);
  }
  else if (nContinuous) {
    FiffDataBufferInfo* firstBlock= NULL;
#ifdef never
    if (nContinuous != kvGetInt(info,"nchannels"))
      Abort("%s:%d: number of continuous blocks (%d) does not match number of channels!\n",
	    __FILE__,__LINE__,nContinuous);
#endif
    checkBlockListConsistency(continuousList);
    slist_totop(continuousList);
    firstBlock= (FiffDataBufferInfo*)slist_get(continuousList);
    kvDefString(info,"dimstr","ct");
    kvDefInt(info,"dt",nContinuous);
    srdrType= fiffTypeMap(firstBlock->type);
    kvDefInt(info,"datatype_in",srdrType);
    kvDefInt(info,"handler_datatype_out",srdrType);
    kvDefInt(info,"dc",fiffNDataItems(firstBlock->type,firstBlock->size));
    kvDefLong(info,"skip.c",0);
  }
  else {
    Abort("%s:%d: input file %s contains no data!\n",
	  __FILE__,__LINE__,self->fileName);
  }
}

static void fiffProcessHeader( FileHandler* self, KVHash* info, 
			      SList* chunkStack )
{
  FiffData* data= (FiffData*)(self->hook);
  KVHash* defs= kvGetHash(info,"definitions");
  long long nextOffset= 0;
  fiffTagRec* rec= NULL;
  int old_bio_big_endian_input= bio_big_endian_input;
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
    rec= readFiffTag(f);
    if (rec->next==FIFFV_NEXT_NONE) {
      thisIsTheLastTag= 1;
    }
    else {
      nextOffset= getNextTagOffset(rec,offset);
    }

    /* Code to handle block start and end */
    if (rec->kind==FIFF_BLOCK_START) {
      int ival= *(int*)rec->data;
      FiffBlockInfo* thisBlockInfo= NULL;
      if (!(thisBlockInfo=(FiffBlockInfo*)malloc(sizeof(FiffBlockInfo))))
	Abort("%s:%d: unable to allocate %d bytes!\n",
	      __FILE__,__LINE__,sizeof(FiffBlockInfo));
      thisBlockInfo->type= ival;
      thisBlockInfo->parentType= currentBlockType;
      thisBlockInfo->start_offset= offset;
      currentBlockType= ival;
      if (debug) fprintf(stderr,"Push block type to %s\n",
			 fiffBlockToString(currentBlockType));
      if (currentBlockType==FIFFB_MEAS && ++numMeasBlocks>1)
	Abort("%s:s: FIFF file contains multiple MEAS blocks!\n",
	      __FILE__,__LINE__);
      slist_push(blockStack,thisBlockInfo);
    }
    else if (rec->kind==FIFF_BLOCK_END) {
      int ival= *(int*)rec->data;
      FiffBlockInfo* thisBlockInfo= slist_pop(blockStack);
      if (ival != thisBlockInfo->type)
	fprintf(stderr,
		"*** block framing error! (found type %d, expected %d ***\n",
		ival,thisBlockInfo->type);
      currentBlockType= thisBlockInfo->parentType;
      if (debug) fprintf(stderr,"Pop block type to %s\n",
			 fiffBlockToString(currentBlockType));
      free(thisBlockInfo);
    }

    switch (currentBlockType) 
      {
      case FIFFB_SUBJECT:
	switch (rec->kind) 
	  {
	  case FIFF_SUBJ_ID: 
	    {
	      kvDefInt(info,"subject_id",*(int*)rec->data);
	      kvDefString(defs,"subject_id","subject identifier");
	    }
	    break;
	  case FIFF_SUBJ_SEX: 
	    {
	      int ival= *(int*)rec->data;
	      if (ival==1)
		kvDefString(info,"sex","male");
	      else if (ival==2)
		kvDefString(info,"sex","female");
	      kvDefString(defs,"sex","Subject gender");
	    }
	    break;
	  case FIFF_SUBJ_WEIGHT: 
	    {
	      kvDefDouble(info,"weight",*(float*)rec->data);
	      kvDefString(defs,"weight","subject weight in kilos");
	    }
	    break;
	  case FIFF_SUBJ_HEIGHT: 
	    {
	      kvDefDouble(info,"height",*(float*)rec->data);
	      kvDefString(defs,"height","subject height in meters");
	    }
	    break;
	  case FIFF_SUBJ_COMMENT:
	    {
	      kvDefString(info,"subject_comment",(char*)rec->data);
	      kvDefString(defs,"subject_commend","misc subject info");
	    }
	    break;
	  }
	break;
      case FIFFB_DACQ_PARS:
	switch (rec->kind) {
	case FIFF_BLOCK_START:
	case FIFF_BLOCK_END:
	  break; /* avoid printing debugging info */
	default:
	  if (debug && rec->type==FIFFT_STRING)
	    fprintf(stderr,"Val is <%s>\n",(char*)rec->data);
	  break;
	}
	break;
      case FIFFB_MEAS_INFO:
	switch (rec->kind) 
	  {
	  case FIFF_EXPERIMENTER: 
	    {
	      kvDefString(info,"experimenter",(char*)rec->data);
	      kvDefString(defs,"experimenter","experimenter");
	    }
	    break;
	  case FIFF_COMMENT:
	    {
	      kvDefString(info,"comment",(char*)rec->data);
	      kvDefString(defs,"comment","meas info comment");
	    }
	    break;
	  case FIFF_PROJ_ID: 
	    {
	      kvDefInt(info,"project_id",*(int*)rec->data);
	      kvDefString(defs,"project_id","project identifier");
	    }
	    break;
	  case FIFF_PROJ_NAME:
	    {
	      kvDefString(info,"project_name",(char*)rec->data);
	      kvDefString(defs,"project_name","project name");
	    }
	    break;
	  case FIFF_PROJ_COMMENT:
	    {
	      kvDefString(info,"project_comment",(char*)rec->data);
	      kvDefString(defs,"project_commend","misc project info");
	    }
	    break;
	  case FIFF_PROJ_PERSONS:
	    {
	      kvDefString(info,"project_persons",(char*)rec->data);
	      kvDefString(defs,"project_persons","people in charge");
	    }
	    break;
	  case FIFF_NCHAN:
	    {
	      kvDefInt(info,"nchannels",*(int*)rec->data);
	      kvDefString(defs,"nchannels","number of channels collected");
	    }
	    break;
	  case FIFF_CH_INFO:
	    {
	      char buf[64];
	      fiffChInfoRec* ch= (fiffChInfoRec*)(rec->data);
	      
	      snprintf(buf,64,"channel.%3d",ch->scanNo);
	      kvDefString(info,buf,ch->ch_name);
	    }
	    break;
	  case FIFF_MEAS_DATE:
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
	  case FIFF_SFREQ:
	    {
	      kvDefDouble(info,"sfreq",*(float*)(rec->data));
	      kvDefString(defs,"sfreq","sampling frequency (hz)");
	    }
	    break;
	  case FIFF_LOWPASS:
	    {
	      kvDefDouble(info,"lowpass",*(float*)(rec->data));
	      kvDefString(defs,"lowpass","analog lowpass thresh (hz)");
	    }
	    break;
	  case FIFF_HIGHPASS:
	    {
	      kvDefDouble(info,"highpass",*(float*)(rec->data));
	      kvDefString(defs,"highpass","analog highpass thresh (hz)");
	    }
	    break;
	  case FIFF_BAD_CHS:
	    /* This block apparently notes which channels are bad, but
	     * I don't have any example of what value denotes badness.
	     */
#ifdef never
	    {
	      int i;
	      int n= fiffNDataItems(rec->type,rec->size);
	      int* buf= (int*)rec->data;
	      for (i=0; i<n; i++) fprintf(stderr,"bad? %d: %d\n",i,buf[i]);
	    }
#endif
	    break;
	  case FIFF_DATA_PACK:
	    /* This information is redundant with the tag type for the
	     * actual data records.
	     */
	    break;
	  case FIFF_BLOCK_START:
	  case FIFF_BLOCK_END:
	    break; /* avoid printing debugging info */
	  default:
	    fprintf(stderr,"Found %s in block type %s, type %s, size %d!\n",
		    fiffKindToString(rec->kind),
		    fiffBlockToString(currentBlockType),
		    fiffTypeToString(rec->type),rec->size);
	  }
	break;
      case FIFFB_RAW_DATA:
	switch (rec->kind) {
	case FIFF_DATA_BUFFER:
	  {
	    slist_append(data->rawBufferList,
			 allocFDBInfo( offset + FIFFC_DATA_OFFSET,
				       rec->size,
				       rec->type ));
	  }
	  break;
	  case FIFF_BLOCK_START:
	  case FIFF_BLOCK_END:
	    break; /* avoid printing debugging info */
	  default:
	    if (debug)
	      fprintf(stderr,"Found %s in block type %s!\n",
		      fiffKindToString(rec->kind),
		      fiffBlockToString(currentBlockType));
	}
	break;
      case FIFFB_PROCESSED_DATA:
	switch (rec->kind) {
	case FIFF_BLOCK_START:
	case FIFF_BLOCK_END:
	  break; /* avoid printing debugging info */
	default:
	  fprintf(stderr,"Found %s in block type %s!\n",
		  fiffKindToString(rec->kind),
		  fiffBlockToString(currentBlockType));
	}
	break;
      case FIFFB_EVOKED:
	switch (rec->kind) {
	case FIFF_BLOCK_START:
	case FIFF_BLOCK_END:
	  break; /* avoid printing debugging info */
	case FIFF_COMMENT:
	case FIFF_EVENT_COMMENT:
	  if (debug) fprintf(stderr,"%s in block type %s: <%s>\n",
			     fiffKindToString(rec->kind),
			     fiffBlockToString(currentBlockType),
			     (char*)rec->data);		  
	  break;
	default:
	  { /* Do nothing */ }
	}
	break;
      case FIFFB_CONTINUOUS_DATA:
	switch (rec->kind) {
	case FIFF_DATA_BUFFER:
	  {
	    slist_append(data->continuousBufferList,
			 allocFDBInfo( offset + FIFFC_DATA_OFFSET,
				       rec->size,
				       rec->type ));
	  }
	  break;
	  case FIFF_BLOCK_START:
	  case FIFF_BLOCK_END:
	    break; /* avoid printing debugging info */
	default:
	  { /* Do nothing */ }
	}
	break;
      case FIFFB_ASPECT:
	switch (rec->kind) {
	case FIFF_BLOCK_START:
	case FIFF_BLOCK_END:
	  break; /* avoid printing debugging info */
	case FIFF_NAVE:
	  {
	    kvDefInt(info,"n_ave",*(int*)rec->data);
	    kvDefString(defs,"n_ave","Number of averages");
	  }
	  break;
	case FIFF_ASPECT_KIND:
	  {
	    kvDefString(info,"aspect",fiffAspectToString(*(int*)rec->data));
	    kvDefString(defs,"aspect","how epochs were subsampled");
	  }
	  break;
	case FIFF_EPOCH:
	  {
	    slist_append(data->epochBufferList,
			 allocFDBInfo( offset + FIFFC_DATA_OFFSET,
				       rec->size,
				       rec->type ));
	  }
	  break;
	default:
	  fprintf(stderr,"Found %s in block type %s; type is %s, size %d!\n",
		  fiffKindToString(rec->kind),
		  fiffBlockToString(currentBlockType),
		  fiffTypeToString(rec->type),rec->size);
	}
	break;
      case FIFFB_SSS_INFO:
	switch (rec->kind) {
	case FIFF_BLOCK_START:
	case FIFF_BLOCK_END:
	  break; /* avoid printing debugging info */
	case FIFF_SSS_JOB:
	  {
	    kvDefInt(info,"sss_job",*(int*)rec->data);
	    kvDefString(defs,"sss_job","????");
	  }
	  break;
	case FIFF_SSS_FRAME:
	  {
	    kvDefInt(info,"sss_frame",*(int*)rec->data);
	    kvDefString(defs,"sss_frame","????");
	  }
	  break;
	case FIFF_SSS_ORD_IN:
	  {
	    kvDefInt(info,"sss_ord_in",*(int*)rec->data);
	    kvDefString(defs,"sss_ord_in","????");
	  }
	  break;
	case FIFF_SSS_ORD_OUT:
	  {
	    kvDefInt(info,"sss_ord_out",*(int*)rec->data);
	    kvDefString(defs,"sss_ord_out","????");
	  }
	  break;
	case FIFF_SSS_NMAG:
	  {
	    kvDefInt(info,"sss_nmag",*(int*)rec->data);
	    kvDefString(defs,"sss_nmag","????");
	  }
	  break;
	case FIFF_SSS_ORIGIN:
	  {
	    float* vals= (float*)rec->data;
	    kvDefDouble(info,"sss_origin.0",vals[0]);
	    kvDefDouble(info,"sss_origin.1",vals[1]);
	    kvDefDouble(info,"sss_origin.2",vals[1]);
	    kvDefString(defs,"sss_origin.0","????");
	    kvDefString(defs,"sss_origin.1","????");
	    kvDefString(defs,"sss_origin.2","????");
	  }
	  break;
	case FIFF_SSS_COMPONENTS:
	  {
	    int* vals= (int*)rec->data;
	    int nvals= rec->size/sizeof(int);
	    int i;
	    char buf[256];
	    if (rec->size % sizeof(int)) 
	      Abort("%s:%d: internal error: size of SSS_COMPONENTS %d is wrong for ints!\n",
		    __FILE__,__LINE__,rec->size);
	    for (i=0; i<nvals; i++) {
	      snprintf(buf,sizeof(buf),"ss_components.%03d",i);
	      kvDefInt(info,buf,vals[i]);
	    }
	  }
	  break;
	default:
	  fprintf(stderr,"Found %s in block type %s; type is %s, size %d!\n",
		  fiffKindToString(rec->kind),
		  fiffBlockToString(currentBlockType),
		  fiffTypeToString(rec->type),rec->size);
	}

      }
    freeTag(rec);
  }

  if (fclose(f)) perror("Error closing input file");

  if (!slist_empty(blockStack))
    fprintf(stderr,"*** Unmatched block_start tag! ***\n");
  slist_destroy(blockStack,free);

  scanBlockList(self,info);

  bio_big_endian_input= old_bio_big_endian_input;
}

static void fiffRead( FileHandler* self, KVHash* info,
		     long long offset, long n,
		     SRDR_Datatype datatype, void* obuf )
{
  FiffData* data= (FiffData*)(self->hook);
  SList* blockList= NULL;
  FiffDataBufferInfo* block= NULL;
  SRDR_Datatype srdrType;
  const char* dimstr= kvGetString(info,"dimstr");

  if (!strcmp(dimstr,"cbt")) {
    blockList= data->rawBufferList;
  }
  else if (!strcmp(dimstr,"bc")) {
    blockList= data->epochBufferList;
  }
  else if (!strcmp(dimstr,"ct")) {
    blockList= data->continuousBufferList;
  }
  else {
    Abort("%s:%d: internal error: dim string <%s> unrecognized on read!\n",
	  __FILE__,__LINE__,dimstr);
  }

  if (offset==0) slist_totop(blockList);
  block= (FiffDataBufferInfo*)slist_next(blockList);
  if (!block) 
    Abort("%s:%d: ran out of data!\n",__FILE__,__LINE__);
  srdrType= fiffTypeMap(block->type);
  if (srdrType!=datatype)
    Abort("%s:%d: Fiff type %s does not match %s!\n",
	  __FILE__,__LINE__,fiffTypeToString(block->type),
	  srdrTypeName[block->type]);
  if (n != fiffNDataItems(block->type,block->size))
    Abort("%s:%d: data block n data item not what was expected! (%ld vs %ld)\n",
	  __FILE__,__LINE__,fiffNDataItems(block->type,block->size),n);
  if (block->type==FIFFT_OLD_PACK) {
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

  if (debug) 
    fprintf(stderr,
	    "Read %d values sz %d, real offset %lld, virt offset %lld\n",
	    n,srdrTypeSize[srdrType],block->offset,offset);
}

static int fiffCompare( FileHandler* f1, FileHandler* f2 )
{
  /* Do comparison based on Unix time, secs since Unix epoch, which
   * we hopefully found stored in the MEAS_INFO.DATE field.
   */
  FiffData* data1= (FiffData*)(f1->hook);
  FiffData* data2= (FiffData*)(f2->hook);

  if (data1->unixTimeSecs<data2->unixTimeSecs) return -1;
  else if (data1->unixTimeSecs<data2->unixTimeSecs) return 1;
  else return 0;
}

FileHandler* fiffFactory(char* fname, KVHash* info)
{
  FileHandler* result= baseFactory(fname);
  FiffData* data;

  result->typeName= strdup("FiffDataHandler");

  result->destroySelf= fiffDestroySelf;
  result->read= fiffRead;
  result->processHeader= fiffProcessHeader;
  result->compareThisType= fiffCompare;

  if (!(data= (FiffData*)malloc(sizeof(FiffData))))
    Abort("%s: unable to allocate %d bytes!\n",sizeof(FiffData));
  
  data->rawBufferList= slist_create();
  data->epochBufferList= slist_create();
  data->continuousBufferList= slist_create();
  data->unixTimeSecs= 0;
  result->hook= data;

  return result;
}

int fiffTester(const char* filename)
{
  long long offset= 0;
  fiffTagRec* rec= NULL;
  FILE* f= NULL;
  int old_bio_big_endian_input;
  int match= 1;

  if ((f = fopen(filename,"r"))==NULL) return 0;

  /* FIFF files are always bigendian */
  old_bio_big_endian_input= bio_big_endian_input;
  bio_big_endian_input= 1;

  /* If this is a FIFF file, the first three tags should be
   * FIFF_FILE_ID, FIFF_DIR_POINTER, and FIFF_FREE_LIST.
   */
  bigfile_fseek(f,offset,SEEK_SET);
  rec= readFiffTag(f);
  if ((rec!=NULL) && (rec->type==FIFFT_ID_STRUCT) && (rec->kind==FIFF_FILE_ID)
      && (rec->next != FIFFV_NEXT_NONE)) {
    freeTag(rec);
    offset= getNextTagOffset(rec,offset);
    bigfile_fseek(f,offset,SEEK_SET);
    rec= readFiffTag(f);
    if ((rec->type==FIFFT_INT) && (rec->kind==FIFF_DIR_POINTER)
	&& (rec->next != FIFFV_NEXT_NONE)) {
      freeTag(rec);
      offset= getNextTagOffset(rec,offset);
      bigfile_fseek(f,offset,SEEK_SET);
      rec= readFiffTag(f);
      if ((rec->type==FIFFT_INT) && (rec->kind==FIFF_FREE_LIST)) {
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

  return match;
}

#else /* ifdef USE_FIFF */

FileHandler* fiffFactory(char* fname, KVHash* info)
{
  FileHandler* result= baseFactory("NotARealFile");
  return result;
}

int fiffTester(const char* filename)
{
  return 0;
}

#endif

