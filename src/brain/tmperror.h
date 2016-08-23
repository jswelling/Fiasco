/************************************************************
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *                                                          *
 *     Copyright (c) 1999 Carnegie Mellon University        *
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
 ***********************************************************/

#ifndef __TMPERROR_LOADED__
#define __TMPERROR_LOADED__

typedef int* ErrorStream;

enum es_level_t { ES_WARN_DEBUG, ES_WARN_INFO, ES_WARN_A, ES_WARN_B, ES_WARN_C, ES_WARN_D,
                  ES_ABORTIVE,   ES_FATAL,     ES_NUM_LEVELS,        ES_EVERY=ES_NUM_LEVELS };

typedef enum es_level_t es_level_type;

enum es_actions { ES_GOBBLE, ES_PRINT, ES_STORE, ES_ECHO,  ES_QUIT,  ES_EXIT };

typedef enum es_actions es_action_type;


ErrorStream     es_new( void );
void            es_delete( ErrorStream* pes );
int             es_any_errors( ErrorStream es );
void            es_write_mesg( ErrorStream es, es_level_type level, const char* format, ... );
void            es_write_error( ErrorStream es, es_level_type level, int code, const char* format, ... );

void            es_set_action( ErrorStream es, es_level_type level, es_action_type action );
void            es_set_stream( ErrorStream es, es_level_type level, FILE* fstream );
void            es_set_prefix( ErrorStream es, es_level_type level, char* prefix );

int             es_clear( ErrorStream es );

#endif /* __TMPERROR_LOADED__ */
