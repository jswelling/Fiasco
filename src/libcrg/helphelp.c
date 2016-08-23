/************************************************************
 *                                                          *
 *  helphelp.c                                              *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1999 Department of Statistics,         *
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
 *  Original programming by Joel Welling 4/1999             *
 ************************************************************/
/* Utility functions used both by buildhelp and the actual running
 * help system.
 */

static char rcsid[] = "$Id: helphelp.c,v 1.13 2006/07/10 22:41:59 welling Exp $";

#include   <stdio.h>
#include   <stdlib.h>
#include   <stdarg.h>
#include   <string.h>
#include   <ctype.h>
#if ( defined(HPPA) || defined(DARWIN) || defined(CYGWIN) )
#include   <sys/ioctl.h>
#else
#include   <stropts.h>
#endif

#ifdef LINUX
#include   <asm/termios.h>
#else
#include   <termios.h>
#endif
#include   <errors.h>
#include   "stdcrg.h"
#include   "slist.h"
#include   "helphelp.h"

extern   int      Help_Num_Topics;
extern   char    *Help_Topics[];
extern   char    *Help_Text[];
extern   HelpTable *Help_Table;

static int slist_contains( SList* list, char* s ) {
  if (slist_empty(list)) return 0;
  slist_totop(list);
  do {
    char* here;
    here= (char*)slist_next(list);
    if (!strcmp(here,s)) return 1;
  } while (!slist_atend(list));
  return 0;
}

HelpTable* help_tables_create() 
{
  HelpTable* result= (HelpTable*)malloc(sizeof(HelpTable));
  if (!result) {
    Abort("help_create_tables: unable to allocate %d bytes!\n",
	  sizeof(HelpTable));
  }
  result->topicCount= 0;
  result->knownSubtopics= slist_create();
  return result;
}

void help_tables_destroy( HelpTable* t )
{
  int i;
  for (i=0; i<t->topicCount; i++) {
    free( t->topicNames[i] );
    slist_destroy( t->topicTable[i], free );
    slist_destroy( t->textTable[i], free );
  }
  slist_destroy( t->knownSubtopics, free );
  free( t );
}

void help_dump_tables(FILE* ofile, HelpTable* t) {
  int i;
  fprintf(ofile,"Topic count: %d\n",t->topicCount);
  for (i=0; i<t->topicCount; i++) {
    char* name= t->topicNames[i];
    SList* subtopics= t->topicTable[i];
    SList* subtexts= t->textTable[i];
    fprintf(ofile,"  %s:\n ",name);
    slist_totop(subtopics);
    slist_totop(subtexts);
    do {
      fprintf(ofile," ---%s---\n",(char*)slist_next(subtopics));
      fprintf(ofile,"%s\n",(char*)slist_next(subtexts));
      fprintf(ofile," --------\n");
    } while (!slist_atend(subtopics));
  }
}

void help_add_to_tables(char* topic, char* text, HelpTable* t) {
  int i;
  char* topicstr= strdup(topic); /* so we can edit it */
  char* subtopic= strrchr(topicstr,':');
  if (!subtopic) {
    /* This is top-level, not a subtopic. It gets notated as an
     * entry in t->topicNames with an associated null entry in
     * the appropriate t->topicTable slot and text in the 
     * corresponding entry in t->textTable.
     */
    free(topicstr);
    for (i=0; i<t->topicCount; i++) {
      if (!strcmp(t->topicNames[i],topic)) {
	Abort("help_add_to_tables: top-level topic name <%s> is not unique!\n",
	      topic);
      }
    }
    if (t->topicCount>=HELP_MAX_TOPICS) {
      Abort("help_add_to_tables: %d separate topics is too many!\n",
	    t->topicCount+1);
    }
    /* Bare topic "Usage" gets mapped to "Overview", for historical
     * reasons.
     */
    if (!strcmp(topic,"Usage"))
      t->topicNames[t->topicCount]= strdup("Overview");
    else t->topicNames[t->topicCount]= strdup(topic);
    t->topicTable[t->topicCount]= slist_create();
    slist_append(t->topicTable[t->topicCount], NULL);
    t->textTable[t->topicCount]= slist_create();
    if (!text) slist_append(t->textTable[t->topicCount], NULL);
    else slist_append(t->textTable[t->topicCount], strdup(text));
    t->topicCount++;
  }
  else {
    int found= 0;
    /* This is a topic-subtopic pair.  We must check to make sure
     * that no top-level topic of the same name exists, find the
     * right subtopic list, and add this text.  We also check to
     * make sure that no identical subtopic exists for another
     * topic, since that would prevent access by subtopic name.
     */
    if (*(subtopic+1)=='\0') {
      Abort("help_add_to_tables: invalid topic name <%s>\n",
	    topic);
    }
    *subtopic= '\0'; /* separate topic and subtopic */
    subtopic++;      /* Get off the : */
    for (i=0; i<t->topicCount; i++) {
      if (!strcmp(t->topicNames[i],topicstr)) {
	/* Found the right slot */
	found= 1;
	slist_totop(t->topicTable[i]);
	if (slist_get(t->topicTable[i])==NULL) {
	  Abort(
	   "help_add_to_tables: topic name <%s> appears both with and without subtopics!\n",
	   topicstr);
	}
	slist_append(t->knownSubtopics, strdup(subtopic));
	slist_append(t->topicTable[i],strdup(subtopic));
	if (!text) slist_append(t->textTable[i],NULL);
	else slist_append(t->textTable[i],strdup(text));
	break;
      }
    }
    if (!found) {
      /* New topic, first subtopic */
      if (t->topicCount>=HELP_MAX_TOPICS) {
	Abort("help_add_to_tables: %d separate topics is too many!\n",
	      t->topicCount+1);
      }
      slist_append(t->knownSubtopics, strdup(subtopic));
      t->topicNames[t->topicCount]= strdup(topicstr);
      t->topicTable[t->topicCount]= slist_create();
      slist_append(t->topicTable[t->topicCount], strdup(subtopic));
      t->textTable[t->topicCount]= slist_create();
      if (!text) slist_append(t->textTable[t->topicCount], NULL);
      else slist_append(t->textTable[t->topicCount], strdup(text));
      t->topicCount++;
    }
    free(topicstr);
  }
}

HelpTable* help_load_topics(int n, char* topic[], char* text[])
/* 
 * This routine creates a new HelpTable from the contents of the given tables.
 */
{
  int i;
  HelpTable* t= help_tables_create();

  for (i=0; i<n; i++) help_add_to_tables(topic[i],text[i],t);

  return t;
}

void help_unload_topics( HelpTable* t )
/*
 * This routine unpacks the given HelpTable into the global table structure.
 */
{
  int n= 0;
  int i;

  /* Out with the old... */
  for (i=0; i<Help_Num_Topics; i++) {
    if (Help_Topics[i]) free (Help_Topics[i]);
    if (Help_Text[i]) free (Help_Text[i]);
    Help_Topics[i]= Help_Text[i]= NULL;
  }
  Help_Num_Topics= 0;

  if (!t) return; /* in case we just want to free old help info */

  for (i=0; i<t->topicCount; i++) {
    char* name= t->topicNames[i];
    SList* subtopics= t->topicTable[i];
    SList* subtexts= t->textTable[i];
    slist_totop(subtopics);
    slist_totop(subtexts);
    if (slist_get(subtopics)) {
      do {
	char* thisSubtopic= (char*)slist_next(subtopics);
	char* thisSubtext= (char*)slist_next(subtexts);
	if (n >= HELP_MAX_TOPICS)
	  Abort("help_unload_tables: more than %d topics!\n",HELP_MAX_TOPICS);
	if (!(Help_Topics[n]= 
	      (char*)malloc( strlen(name)+strlen(thisSubtopic)+2 )))
	  Abort("help_unload_topics: unable to allocate %d bytes!\n",
		strlen(name)+strlen(thisSubtopic)+2);
	sprintf(Help_Topics[n],"%s:%s",name,thisSubtopic);
	Help_Text[n]= strdup(thisSubtext);
	n++;
      } while (!slist_atend(subtopics));
    }
    else {
      if (n >= HELP_MAX_TOPICS)
	Abort("help_unload_tables: more than %d topics!\n",HELP_MAX_TOPICS);
      Help_Topics[n]= strdup(name);
      if (slist_get(subtexts)) {
	Help_Text[n]= strdup( (char*)slist_next(subtexts) );
      }
      else {
	Help_Text[n]= NULL;
      }
      n++;
    }
  }
  Help_Num_Topics= n;
  Help_Table= t;
}

static int line_width()
{
  struct winsize ws;
  ioctl(1,TIOCGWINSZ,&ws);
  return (int)(ws.ws_col);
}

void help_emit_topics(HelpTable* t, FILE* ofile)
{
  int i;
  int width= line_width();
  int where= 0;
  int indent= 0;
  for (i=0; i<t->topicCount; i++) {
    char* name= t->topicNames[i];
    SList* subtopics= t->topicTable[i];
    slist_totop(subtopics);
    if (slist_get(subtopics)) {
      fprintf(ofile,"  %s:   ",name);
      where= indent= strlen(name)+6;
      do {
	char* s= (char*)slist_next(subtopics);
	int delta= strlen(s) + 2;
	if (where+delta >= width) {
	  fprintf(ofile,"\n      ");
	  where= 6;
	}
	fprintf(ofile,"%s  ",s);
	where += strlen(s)+2;
      } while (!slist_atend(subtopics));
      fprintf(ofile,"\n");
    }
    else fprintf(ofile,"  %s\n",name);
  }  
}
