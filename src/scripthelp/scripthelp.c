/************************************************************
 *                                                          *
 *  shellhelp.c                                     *
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
 *  Original programming by Joel Welling, 5/98              *
 ************************************************************/
/*************************************************************

  DESCRIPTION OF SCRIPTHELP

  This tool provides a help interface for use by shell scripts.

  shellhelp scriptname [topic]

**************************************************************/

#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include "mri.h"
#include "fmri.h"
#include "misc.h"
#include "stdcrg.h"
#include "helphelp.h"

static char rcsid[] = "$Id: scripthelp.c,v 1.3 2003/02/07 21:32:13 welling Exp $";

extern   HelpTable *Help_Table; /* Hook to the global help table */
extern   char*      Help_Progname; /* ... and Help's notion of current prog */

typedef struct topic_pair_struct {
  char* name;
  HelpTable* tbl;
} TopicPair;

static char* progname;
static SList* topics=NULL; /* a list of TopicPairs */

static void parseHelp( HelpTable* t )
{
  int i;
  if (topics) 
    slist_destroy(topics,(void (*)(void *))help_tables_destroy);
  topics= slist_create();

  for (i=0; i<t->topicCount; i++) {
    char* name= t->topicNames[i];
    SList* subtopics= t->topicTable[i];
    SList* subtexts= t->textTable[i];
    slist_totop(subtopics);
    slist_totop(subtexts);
    if (slist_get(subtopics)) {
      TopicPair* p;
      if (!(p= (TopicPair*)malloc(sizeof(TopicPair))))
	Abort("parseHelp: unable to allocate %d bytes!\n",sizeof(TopicPair));
      p->name= strdup(name);
      p->tbl= help_tables_create();
      do {
	help_add_to_tables((char*)slist_next(subtopics),
			   (char*)slist_next(subtexts),
			   p->tbl);
      } while (!slist_atend(subtopics));
      help_add_to_tables("helphelp", NULL, p->tbl);
      help_add_to_tables("alltopics", NULL, p->tbl);
      help_add_to_tables("selecttopic", NULL, p->tbl);
      help_add_to_tables("listtopics", NULL, p->tbl);
      help_add_to_tables("htmldoc", NULL, p->tbl);
      help_add_to_tables("quithelp", NULL, p->tbl);
      slist_append(topics,p);
    }
  }
}

static int prepHelp( char* name )
{
  HelpTable* t;

  slist_totop( topics );
  do {
    TopicPair* p= (TopicPair*)slist_next(topics);
    if (!strcmp(p->name,name)) {
      help_unload_topics( p->tbl );
      return 1;
    }
  } while (!slist_atend(topics));
  return 0;
}

int main( int argc, char** argv ) 
{
  HelpTable* tblForAll;
  char* topic;

  progname= argv[0];

  if (argc > 1) {
    /* we want to grab the root name of any command passed in. */
    if (topic= strrchr(argv[1],'/')) {
      topic++; /* skip past the '/' */
    }
    else topic= argv[1]; 
  }
  else topic= "selecttopic";

  Help_init(); /* Set up the overall help data structure for all topics */
  tblForAll= Help_Table; /* Snatch the help table from the global structures */
  parseHelp(tblForAll);

  if (argc == 1) {
    help_unload_topics( tblForAll );
    Help_Progname= "FiascoScripts";
    Help_interp( "selecttopic" );
  }
  else if (argc == 2) {
    if (prepHelp( topic )) {
      Help_Progname= topic;
      Help_interp( "selecttopic" );
    }
    else {
      help_unload_topics( tblForAll );
      Help_Progname= "FiascoScripts";
      Help_interp( topic );
    }
  }
  else if (argc == 3) {
    if (!strcmp(argv[1],"-help")) {
      help_unload_topics( tblForAll );
      Help_Progname= "FiascoScripts";
      Help_interp( argv[2] );
    }
    else {
      if (!prepHelp( topic )) {
	Error("Requested script %s has no help information.\n",topic);
	Help("Usage");
	exit(-1);
      }
      Help_Progname= topic;
      Help_interp( argv[2] );
    }
  }
  else {
    Help("Usage");
    exit(-1);
  }

  return 0;
}

