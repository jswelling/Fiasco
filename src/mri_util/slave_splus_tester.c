#include <stdio.h>
#include "slave_splus.h"

static char rcsid[] = "$Id: slave_splus_tester.c,v 1.2 1997/11/26 00:39:03 welling Exp $";

static const char *exename= "/usr/statlocal/bin/Splus";
static const char* s_dlo= "splus_binary_pipes.o";
static const char* s_bio_init= "bio_init.S";

main(int argc, char* argv[])
{
  SlaveSplus* sp= NULL;

  if (!(sp= spawn_slave_splus(exename,s_dlo,s_bio_init))) {
    fprintf(stderr,"%s: spawn failed!\n",argv[0]);
  }

  fprintf(stderr,"Bang!\n");
  kill_slave_splus(sp);
  fprintf(stderr,"Rtn from kill\n");
  
}
