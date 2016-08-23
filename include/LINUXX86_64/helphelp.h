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
 *  Original programming by Joel Welling 4/1999             *
 ************************************************************/
/* Header for helphelp.c
 */

#include "slist.h"

typedef struct help_table_struct {
  char* topicNames[HELP_MAX_TOPICS];
  SList* topicTable[HELP_MAX_TOPICS];
  SList* textTable[HELP_MAX_TOPICS];
  SList* knownSubtopics;
  int topicCount;
} HelpTable;

HelpTable* help_tables_create();

void help_tables_destroy(HelpTable* t);

void help_dump_tables(FILE* ofile, HelpTable* t);

void help_add_to_tables(char* topic, char* text, HelpTable* t);

HelpTable* help_load_topics(int n, char* topic[], char* text[] );

void help_unload_topics(HelpTable* t);

void help_emit_topics(HelpTable* t, FILE* ofile);

