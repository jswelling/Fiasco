/************************************************************
 *                                                          *
 *  slave_splus.h                                           *
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
 *  Original programming by Joel Welling, 8/97              *
 ************************************************************/

/* spawn returns 0 on failure */

#ifdef LINUX
#include <sys/types.h>
#endif

typedef struct slave_splus_struct {
  pid_t pid;
  FILE* in;
  FILE* out;
  int din;
  int slave_side_din;
  int dout;
  int slave_side_dout;
} SlaveSplus;

extern SlaveSplus* spawn_slave_splus( char* exename, char* s_dlo, 
				      char* s_bio_init );
/* returns null on failure */

extern void kill_slave_splus( SlaveSplus* target );
