/************************************************************
 *                                                          *
 *  parsesplit.h                                            *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1997 Department of Statistics,         *
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
 *  Original programming by Joel Welling, 8/98              *
 ************************************************************/
/* This module contains routines for parsing split and condition files */

typedef struct sp_cond_def_struct {
  int id;
  char* name;
  char** factor_lvl;
  int* factor_intlvl;
} sp_ConditionDef;

typedef struct sp_split_rec_struct {
  int image;
  int slice;
  int cond;
} sp_SplitRec;

int sp_parse_conditions( char* fname, 
			 char*** factor_table, int* nfactors,
			 sp_ConditionDef*** cond_table, int* nconditions );

int sp_parse_split( char* fname, int nslices, int nimages, int nconditions,
		    sp_SplitRec** split_table );

int sp_match_missing( sp_SplitRec* split_table, unsigned char** missing,
		      int nslices, int nimages );

