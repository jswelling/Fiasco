/************************************************************
 *                                                          *
 *  afni_defs.h                                             *
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
 *  Original programming by Joel Welling 7/2002             *
 ************************************************************/

/* From the AFNI header files and documentation */
#define THD_MAX_NAME      256
#define ANAT_SPGR_TYPE   0
#define ANAT_FSE_TYPE    1
#define ANAT_EPI_TYPE    2
#define ANAT_MRAN_TYPE   3
#define ANAT_CT_TYPE     4
#define ANAT_SPECT_TYPE  5
#define ANAT_PET_TYPE    6
#define ANAT_MRA_TYPE    7
#define ANAT_BMAP_TYPE   8
#define ANAT_DIFF_TYPE   9
#define ANAT_OMRI_TYPE   10
#define ANAT_BUCK_TYPE   11
#define FUNC_FIM_TYPE   0  /* 1 value           */
#define FUNC_THR_TYPE   1  /* obsolete          */
#define FUNC_COR_TYPE   2  /* fico: correlation */
#define FUNC_TT_TYPE    3  /* fitt: t-statistic */
#define FUNC_FT_TYPE    4  /* fift: F-statistic */
#define FUNC_ZT_TYPE    5  /* fizt: z-score     */
#define FUNC_CT_TYPE    6  /* fict: Chi squared */
#define FUNC_BT_TYPE    7  /* fibt: Beta stat   */
#define FUNC_BN_TYPE    8  /* fibn: Binomial    */
#define FUNC_GT_TYPE    9  /* figt: Gamma       */
#define FUNC_PT_TYPE    10 /* fipt: Poisson     */
#define FUNC_BUCK_TYPE  11 /* fbuc: bucket      */
#define THREEDIM_HEAD_ANAT 0
#define THREEDIM_HEAD_FUNC 1
#define THREEDIM_GEN_ANAT 2
#define THREEDIM_GEN_FUNC 3
#define UNITS_MSEC_TYPE  77001  
#define UNITS_SEC_TYPE   77002 
#define UNITS_HZ_TYPE    77003
#define ORI_R2L_TYPE  0  /* Right to Left         */
#define ORI_L2R_TYPE  1  /* Left to Right         */
#define ORI_P2A_TYPE  2  /* Posterior to Anterior */
#define ORI_A2P_TYPE  3  /* Anterior to Posterior */
#define ORI_I2S_TYPE  4  /* Inferior to Superior  */
#define ORI_S2I_TYPE  5  /* Superior to Inferior  */
#define MARKS_MAXNUM  10
#define MARKS_MAXLAB  20
#define MARKS_MAXHELP 256
#define MARKS_MAXFLAG 8
#define VIEW_ORIGINAL_TYPE    0
#define VIEW_ACPCALIGNED_TYPE 1
#define VIEW_TALAIRACH_TYPE   2
#define VIEW_REGISTERED_TYPE  3

static const char* typestring_names[]= { "3DIM_HEAD_ANAT", "3DIM_HEAD_FUNC",
				   "3DIM_GEN_ANAT", "3DIM_GEN_FUNC", 
				   (char*)NULL };
static const char* coordsys_names[]= { 
  "+orig", "+acpc", "+tlrc", (char*)NULL };
static const char* orientation_names[]= {
  "R2L", "L2R", "P2A", "A2P", "I2S", "S2I", (char*)NULL };
static const char* anat_type_names[]= {
  "ANAT_SPGR_TYPE","ANAT_FSE_TYPE","ANAT_EPI_TYPE","ANAT_MRAN_TYPE",
  "ANAT_CT_TYPE","ANAT_SPECT_TYPE","ANAT_PET_TYPE","ANAT_MRA_TYPE",
  "ANAT_BMP_TYPE","ANAT_DIFF_TYPE","ANAT_OMRI_TYPE","ANAT_BUCK_TYPE",
  (char*)NULL
};
static const char* func_type_names[]= {
  "FUNC_FIM_TYPE","FUNC_THR_TYPE","FUNC_COR_TYPE","FUNC_TT_TYPE",
  "FUNC_FT_TYPE","FUNC_ZT_TYPE","FUNC_CT_TYPE","FUNC_BT_TYPE",
  "FUNC_BN_TYPE","FUNC_GT_TYPE","FUNC_PT_TYPE","FUNC_BUCK_TYPE",
  (char*)NULL
};
