/************************************************************
 *                                                          *
 *  slave_splus.c                                           *
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
 *  Original programming by Joel Welling, 8/97              *
 ************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/signal.h>
#include <sys/wait.h>
#include <string.h>
#include "slave_splus.h"

#define INBUF_LENGTH 256
#define READ_TIMEOUT_SEC 20

static char rcsid[] = "$Id: slave_splus.c,v 1.7 2003/09/25 19:45:42 welling Exp $";

static const char* s_BioSetFdes_return_string= "[1] 1";

static char* fgets_with_timeout( char *s, int n, FILE *stream, int seconds )
{
  unsigned int prev_alarm;
  char* p;

  prev_alarm= alarm(seconds);
  p= fgets(s, n, stream);
  alarm(prev_alarm);
  return p;
}

SlaveSplus* spawn_slave_splus( char* exename, char* s_dlo, char* s_bio_init )
{
  pid_t pid;
  int to_slave[2];
  int from_slave[2];
  int to_slave_data[2];
  int from_slave_data[2];
  int i;
  char buf[INBUF_LENGTH];
  char *fake_argv[2];

  if (pipe(to_slave)) {
    perror("Error creating pipe to slave");
    return NULL;
  }
  if (pipe(from_slave)) {
    perror("Error creating pipe from slave");
    return NULL;
  }
  if (pipe(to_slave_data)) {
    perror("Error creating data pipe to slave");
    return NULL;
  }
  if (pipe(from_slave_data)) {
    perror("Error creating data pipe from slave");
    return NULL;
  }

  fflush(stdout);
  if ( (pid= fork()) == 0 ) { /* This is the spawned process */
    /* Close stdin and replace it with pipe input */
    close(0);
    if (dup(to_slave[0]) == -1) perror("dup error in child stdin");
    close(to_slave[1]);

    /* Close stdout and replace it with pipe output */
    close(1);
    if (dup(from_slave[1]) == -1) perror("dup error in child stdout");
    close(from_slave[0]);

    /* Close the uninteresting sides of the data pipes */
    close(to_slave_data[1]);
    close(from_slave_data[0]);

    /* Become the new executable */
    fake_argv[0]= exename;
    fake_argv[1]= NULL;
    execvp( exename, fake_argv );
    fprintf(stderr,"Couldn't execute \"%s\"!\n",exename);
    perror("spawn error");
    return NULL;  /* Should never reach this point */
  }
  else if (pid == -1) {
    perror("spawn_slave_splus: fork failed");
    return NULL;
  }
  else { /* This is the parent */

    SlaveSplus* result= NULL;

    if (!(result= (SlaveSplus*)malloc(sizeof(SlaveSplus)))) {
      fprintf(stderr,"spawn_slave_splus: unable to allocate %d bytes!\n",
	      (int)(sizeof(SlaveSplus)));
      exit(-1);
    }
    result->pid= pid;
    if (!(result->in= fdopen(to_slave[1],"w"))) {
      perror("Error for parent opening slave input stream");
      return NULL;
    }
    close(to_slave[0]);
    if (!(result->out= fdopen(from_slave[0],"r"))) {
      perror("Error for parent opening slave output stream");
      return NULL;
    }
    close(from_slave[1]);

    /* Keep info re data pipes */
    result->din= to_slave_data[1];
    result->slave_side_din= to_slave_data[0];
    close(to_slave_data[0]);
    result->dout= from_slave_data[0];
    result->slave_side_dout= from_slave_data[1];
    close(from_slave_data[1]);

    /* Carry out initial transactions */
    for (i=0; i<4; i++) {
      if (!fgets_with_timeout(buf,INBUF_LENGTH,result->out,
			      READ_TIMEOUT_SEC)) {
	if (ferror(result->out)) 
	  perror("Error reading slave splus startup info");
	else fprintf(stderr,"EOF reading slave splus startup info!\n");
	return NULL;
      }
      if (strlen(buf)>0) buf[strlen(buf)-1]= 0;
      fprintf(stderr,"Slave splus: <%s>\n",buf);
    }
    
    fprintf(result->in,"dyn.load2(\"%s\")\n",s_dlo);
    fflush(result->in);
    
    fprintf(result->in,"source(\"%s\")\n",s_bio_init);
    fflush(result->in);
    
    fprintf(result->in,"BioSetFdes( %d , %d )\n",
	    result->slave_side_din, result->slave_side_dout);
    fflush(result->in);
    if (!fgets_with_timeout(buf,INBUF_LENGTH,result->out,
			    READ_TIMEOUT_SEC)) {
      if (ferror(result->out)) 
	perror("Error reading slave splus initialization info");
      else fprintf(stderr,"EOF reading slave splus initialization info!\n");
      return NULL;
    }
    if (strlen(buf)>0) buf[strlen(buf)-1]= 0;
    if (strcmp(buf,s_BioSetFdes_return_string)) {
      fprintf(stderr,"got <%s> rather than <%s> from BioSetFdes\n",
	      buf,s_BioSetFdes_return_string);
      fprintf(stderr,"Slave splus apparent initialization error.\n");
      return NULL;
    }
    else fprintf(stderr,"Slave splus initialization appears correct.\n");
    
    return result;
  }

#if ( !( SGI64 || SGI5 || SGI6 || SGIMP64 ) )
  return NULL; /* never happens, but some compilers complain without it */
#endif
}

void kill_slave_splus( SlaveSplus* target ) {
  int wait_stat;
  fprintf(target->in,"q()\n");
  fflush(target->in);
  if ( waitpid(-1,&wait_stat,0) != target->pid ) 
    fprintf(stderr,"kill_slave_splus: Trouble waiting on Splus exit!\n");
  if (WIFEXITED(wait_stat)) {
    if (WEXITSTATUS(wait_stat) != 0) 
      fprintf(stderr,
	      "kill_slave_splus: Slave Splus process exited with status %d\n",
	      WEXITSTATUS(wait_stat));
  }
  else {
    fprintf(stderr,
	    "kill_slave_splus: Slave Splus (pid %d) failed to exit!\n",
	    target->pid);
  }
  if (fclose(target->in)) 
    perror("kill_slave_splus: Error closing Splus input");
  if (fclose(target->out)) 
    perror("kill_slave_splus: Error closing Splus output");
  if (close(target->din)) 
    perror("kill_slave_splus: Error closing Splus data input");
  if (close(target->dout)) 
    perror("kill_slave_splus: Error closing Splus data output");
  free(target);
}


