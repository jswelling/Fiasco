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

static char rcsid[] = "$Id: help.c,v 1.17 2007/03/21 23:53:28 welling Exp $";

#include   <stdio.h>
#include   <stdlib.h>
#include   <stdarg.h>
#include <string.h>
#include   <strings.h>
#include   <ctype.h>

#include  <errors.h>
#include  "stdcrg.h"
#include  "helphelp.h"

/*
 * HELP engages an automated help system to provide the user of a program
 * with context sensitive help.  
 *
 * Associated with each program using Help() is a text file containing help
 * text organized by selected topics.  The program buildhelp contains documentation
 * on the format of this file; buildhelp is invoked during compilation to convert
 * the help text file into a loadable object file containing the functions 
 * Help_init(), which initializes the internal data structures, and Help_listtopics()
 * and Help_help(), which display useful information.
 *
 * Help can be engaged with a list of topic names or one of the reserved names
 * ALLTOPICS or SELECTTOPIC, described below.  In the former case, the topic
 * list is a space-separated list of topic names.
 
 * In either invoked or interactive mode, Help() requires only enough of a string
 * to disambiguously select a topic.  All names are compared without regard to case,
 * and white space before or after topic names are ignored.  When completed, Help()
 * aborts execution with a call to exit(0).
 *
 * Reserved Names:
 *
     o  ALLTOPICS      Print the entire help file in order with burst lines indicating topic
     o  SELECTTOPIC    Engages an interactive help system where user is prompted for
                       successive topic names.  If there is only one topic, this is equivalent
                       to ALLTOPICS
     o  LISTTOPICS     Lists available topics, only available in interactive mode.
     o  HELPHELP       Provides help text concerning the help system itself.
     o  HTMLDOC        Produces an html document of help info.
     o  QUITHELP       Quits help, only available in interactive mode.
 *
 *
 * In most cases, these can be abbreviated as "all", "select", "list", "quit", "help", "html" and often
 * it will be possible to abbreviate them as "a", "s", "l", "e", "ht" and "he".
 *
 */

int      Help_Num_Topics = 0; 

short    Help_Matching_Topics[HELP_MAX_TOPICS];

char    *Help_Progname= NULL;
char    *Help_Creation_Time= NULL;
char    *Help_Topics[HELP_MAX_TOPICS];
char    *Help_Text[HELP_MAX_TOPICS];
char    *Help_Intro;
HelpTable* Help_Table;

int      Help_match( char *, int *, int );

static char* Help_helphelp_text= "\n\
  The FIASCO Help System provides interactive and non-interactive\n\
  help for FIASCO tools and generates the documentation web pages\n\
  for the tools.  For a brief summary of a tool named \"thistool\" \n\
  try the command:\n\
\n\
    thistool -help usage\n\
\n\
  For more detailed interactive help try:\n\
\n\
    thistool -help \n\
\n\
  To exit interactive help type \"quit\".\n";


static void emit_html(const char* text)
{
  const char* runner;
  printf( "<pre>\n");
  for (runner= text; *runner; runner++) {
    switch (*runner) {
    case '>': fputs("&gt;",stdout); break;
    case '<': fputs("&lt;",stdout); break;
    case '&': fputs("&amp;",stdout); break;
    case '\\': if (*(runner+1)) { runner++; putchar(*runner); }; break;
    default: putchar(*runner); break;
    }
  }
  printf( "</pre>\n");
}

static void emit_plaintext(FILE* f, const char* text)
{
  typedef enum { 
    OUT_MODE, TAG_MODE, ESC_MODE, TAG_ESC_MODE
  } Mode;
  Mode mode= OUT_MODE;
  const char* runner= text;
  while (*runner) {
    switch (mode) {
    case OUT_MODE:
      {
	if (*runner == '\\') mode= ESC_MODE;
	else {
	  fputc(*runner, f);
	}
      }
      break;
    case TAG_MODE:
      {
	if (*runner == '\\') mode= TAG_ESC_MODE;
      }
      break;
    case ESC_MODE:
      {
	if (*runner == '<') mode= TAG_MODE;
	else {
	  mode= OUT_MODE;
	  fputc(*runner,f);
	}
      }
      break;
    case TAG_ESC_MODE:
      {
	if (*runner=='>') mode= OUT_MODE;
	else mode= TAG_MODE;
      }
      break;
    default:
      Abort("%s: emit_plaintext: internal error on string <%s>\n",text);
    }
    runner++;
  }
  fputc('\n',f);
}

void Help_interp( char* topic )
{
    int      i, j, matches, ind, len, comp, useout = 0, interactive = 0;
    char     *p, *pager;
    char     line[1024];
    char     topic_cpy[1024];
    FILE     *pp;

    /* Make a local copy of topic, in case it's stored in 
     * read-only memory so that strtok chokes on it.
     */
    strncpy(topic_cpy, topic, sizeof(topic_cpy));
    topic_cpy[sizeof(topic_cpy)-1]= '\0';

    /* Grab first topic name; if none, exit below. */
    p = strtok( topic_cpy, " " );  

    while( p != NULL )
    {
        matches = Help_match( p, &len, HELP_COMPLETION_CHAR );  /* Use NULL in 3rd arg for no completion */

	/* May have to substitute "Overview" or "Introduction" for "Usage" for 
	 * historical reasons.
	 */
	if (!matches && (!strcasecmp(p,"usage")))
	  matches = Help_match( "Overview", &len, HELP_COMPLETION_CHAR );  /* Use NULL in 3rd arg for no completion */
	if (!matches && (!strcasecmp(p,"usage")))
	  matches = Help_match( "Introduction", &len, HELP_COMPLETION_CHAR );  /* Use NULL in 3rd arg for no completion */

        if( !matches )
            fprintf( stderr, "Help Error: Topic %s not found in help data base.  Try again.\n", p );
        else if( matches > 1 )
        {
            printf( "Given topic \"%s\" is ambiguous.\nPossible completions are: \n", p );

            for( i = 0; i < ind; i++ )
                printf( "\t %s\n", Help_Topics[Help_Matching_Topics[i]] );
        }
        else if( matches == 1 )
        {
            ind = Help_Matching_Topics[0];
                
            if( Help_Text[ind] == NULL ) /* Special Action Required */
            {
	        if ( !strncasecmp( p, "htmldoc", len ) ) {

		  printf("<!DOCTYPE HTML PUBLIC \"-//IETF//DTD HTML//EN\">\n");
		  printf("<html>\n");
		  printf("  <head>\n");
		  printf("    <title>%s</title>\n",Help_Progname);
		  printf("  </head>\n");
		  printf("  <body>\n");
		  printf("    <h1>%s</h1>\n",Help_Progname);

		  for (i=0; i<Help_Table->topicCount; i++) {
		    char* name= Help_Table->topicNames[i];
		    SList* subtopics= Help_Table->topicTable[i];
		    SList* subtexts= Help_Table->textTable[i];
		    slist_totop(subtopics);
		    slist_totop(subtexts);
		    if (slist_get(subtexts)) {
		      if (slist_get(subtopics)) {
			printf( "<h2>%s</h2>\n", name );
			do {
			  printf("<h3>%s</h3>\n", 
				 (char*)slist_next(subtopics));
			  emit_html((char*)slist_next(subtexts));
			} while (!slist_atend(subtopics));
		      }
		      else {
			printf( "<h2>%s</h2>\n", name );
			emit_html((char*)slist_get(subtexts));
		      }
		    }
		  }

		  printf("    <hr>\n");
		  printf("    (automatically generated by %s version of %s)\n",
			 Help_Progname,Help_Creation_Time);
		  printf("  </body>\n");
		  printf("</html>\n");

	        }
                else if( !strncasecmp( p, "alltopics", len ) )
                {
                    for( i = 0; i < Help_Num_Topics; i++ )
                    {
                        if( Help_Text[i] == NULL || !strcasecmp("helphelp", Help_Topics[i]) )
                            continue;
                    
                        printf( "Topic: %s\n", Help_Topics[i] );
			emit_plaintext(stdout, Help_Text[i]);
                    }

                    break;
                }
                else if( !strncasecmp( p, "selecttopic", len ) )
                {
                    interactive = 1;
                    break;
                }
                else if( !strncasecmp( p, "listtopics", len ) )
                    Help_listtopics();
                else if( !strncasecmp( p, "quithelp", len ) )
                    fprintf( stderr, "Help Error: QUITHELP command only allowed in interactive mode.\n" );
		else if ( !strncasecmp( p, "helphelp", len ) ) {
		  printf( "Topic: %s\n", p );
		  printf( "%s\n", Help_helphelp_text );
		}
                else
                    Error( "Unrecognized help command %s.", p ); /* SHOULD NOT HAPPEN */
            }
            else
            {
                printf( "Topic: %s\n", Help_Topics[ind] );
		emit_plaintext(stdout, Help_Text[ind]);
            }
        }
        else if( matches < 0 )    /* Topic name completion */
        {
            matches *= -1;

            for( i = 0; i < matches; i++ )
            {
                ind = Help_Matching_Topics[i];

                if( Help_Text[ind] != NULL ) /* Commands are ignored during completion */
                {
                    printf( "Topic: %s\n", Help_Topics[ind] );
		    emit_plaintext(stdout, Help_Text[ind]);
                }
            }
        }
        
        p = strtok( NULL, " " );   /* Grab next topic name; if none, exit below. */
    }

    if( interactive )
    {
        pager = getenv( "PAGER" );

        printf( "%s\n", Help_Intro );     /* Intro burst */
        Help_listtopics();
        printf( "\n" );

        for( ;; )
        {
            printf( "Help? " );

            if( fgets( line, 1024, stdin ) == NULL )
                break;

            p = strtok( line, " \n\t" );      /* Grab first topic name; if none, continue. */

            if( p == NULL )
                continue;
            
            matches = Help_match( p, &len, HELP_COMPLETION_CHAR );

            if( !matches )
            {
                printf( "Topic %s could not be found in help data base.  Try again.\n", p );
                continue;
            }
            else if( matches > 1 )
            {
                printf( "Given topic \"%s\" is ambiguous.\nPossible completions are: \n", p );

                for( i = 0; i < matches; i++ )
                    printf( "\t %s\n", Help_Topics[Help_Matching_Topics[i]] );
            }
            else if( matches == 1 )
            {
                ind = Help_Matching_Topics[0];
                
                if( Help_Text[ind] == NULL ) /* Special Action Required */
                {
		    if (!strncasecmp( p, "htmldoc", len )) {
		      fprintf(stderr," HTMLDOC command not allowed in interactive mode.\n");
		    }
		    else if( !strncasecmp( p, "alltopics", len ) )
		      {
			fprintf( stderr, "ALLTOPICS command not allowed in interactive mode.\n" );
			continue;
		      }
                    else if( !strncasecmp( p, "selecttopic", len ) )
		      {
			fprintf( stderr, "SELECTTOPIC command not allowed in interactive mode.\n" );
			continue;
		      }
                    else if( !strncasecmp( p, "listtopics", len ) )
		      {
			Help_listtopics();
			continue;
		      }
                    else if( !strncasecmp( p, "quithelp", len ) )
		      break;
		    else if ( !strncasecmp( p, "helphelp", len ) ) {
		      printf( "Topic: %s\n", p );
		      printf( "%s\n", Help_helphelp_text );
		    }
                    else
		      Error( "Unrecognized help command %s.", p ); /* SHOULD NOT HAPPEN */
                }
                else
                {
#ifndef T3D
                    if( pager != NULL )
                        pp = popen( pager, "w" );
                    else
                        pp = popen( HELP_DEFAULT_PAGER, "w" );

                    if( pp == NULL )
                    {
                        Warning( 1, "Pager could not be opened, using stderr instead." );
                        pp = stdout;
                        useout = 1;
                    }
#else
		  pp = stdout;
		  useout = 1;
#endif

                    fprintf( pp, "Topic: %s\n", Help_Topics[ind] );
		    emit_plaintext( pp, Help_Text[ind] );

#ifndef T3D
                    if( !useout && (127 == pclose( pp )) )
                        Warning( 1, "Pager failed because popen() could not execute shell." );
#endif
                }
            }
            else if( matches < 0 )    /* Topic name completion */
            {
                matches *= -1;

#ifndef T3D
                if( pager != NULL )
                    pp = popen( pager, "w" );
                else
                    pp = popen( HELP_DEFAULT_PAGER, "w" );

                if( pp == NULL )
                {
                    Warning( 1, "Pager could not be opened, using stderr instead." );
                    pp = stdout;
                    useout = 1;
                }
#else
		pp = stdout;
		useout = 1;
#endif

                for( i = 0; i < matches; i++ )
                {
                    ind = Help_Matching_Topics[i];

                    if( Help_Text[ind] != NULL )        /* Commands are ignored during completion */
                    {
                        fprintf( pp, "Topic: %s\n", Help_Topics[ind] );
			emit_plaintext(pp, Help_Text[ind]);
                    }
                }

#ifndef T3D
                if( !useout && (127 == pclose( pp )) )
                    Warning( 1, "Pager failed because popen() could not execute shell." );
#endif
            }
        }
    }

    exit(0);
}

int testHelp( int* argc, char* argv[] )
{
  /* Check to see if help was requested */
  if( ( *argc > 1 ) && !strcmp( argv[1], "-help" ) )
    {
      if( *argc == 2 )
	Help( "selecttopic" );
      else
	Help( argv[2] );
      return 1;
    }
  else return 0;
}

void   Help( char *topic )
{
    Help_init();               /* Initialize Help System */
    Help_interp( topic );      /* Handle the topic */
}

/*
 * HELP_MATCH searches for a unambiguous match of the (partial) topic name given
 * to the list of available help topics.  Since it is unlikely that there will
 * be a very long list of available topics, a simple search is used instead of
 * a fancy search or hashing algorithm.  Matching topic names are stored by index
 * in Help_Matching_Topics and the number of matching topics is returned.  If a
 * completion character is supplied and the given topic name ends with the completion
 * character, then the return value is flagged by negation.  This allows a completion
 * character to be used to list a set of topics that match an initial substring; for
 * instance, to allow for subtopics to be listed as a whole. If no completion character
 * is used, then an exact match takes precedence over a partial match.
 *
 */

int    Help_match( char *p, int *lenp, int compl_ch )
{

  int      i, len, matches = 0, complete = 0;
  
  *lenp = len = strlen(p);

  if( compl_ch && (p[len-1] == compl_ch) )
    {
      complete = 1;
      p[len-1] = '\0';
      *lenp = --len;
    }

  for( i = 0; i < Help_Num_Topics; i++ ) {
    char* target= Help_Topics[i];
    char* subtopic= strrchr(target,':');
    if (subtopic) {
      target= subtopic+1; /* advance to after : */
      if (*target == '\0') continue; /* should have been caught earlier */
    }
    
    if( strncasecmp( p, target, len ) )    /* Partial Match */
      continue;
    
    if( !complete && !strcasecmp( p, target ) )  /* Exact Match   */
      {
	matches = 1;
	Help_Matching_Topics[0] = i;
	break;
      }
    
    Help_Matching_Topics[matches++] = i;
  }

  if( matches > 1 && complete )
    matches *= -1;

  return( matches );
}


/*
 * Below are dummy versions of Help_init() and Help_listtopics().
 * If not defined elsewhere (i.e., in the output of buildhelp),
 * these will be used.  If they are defined elsewhere, these
 * will not be linked.
 *
 * Note: this works because these functions are compiled into a library.
 *
 */
 
/*  --this technique will not work on the Alphas -- ghood@psc.edu
 void   Help_init( void )
{
    fprintf( stderr, "Help system not available; consult other documentation for usage information.\n" );
    exit( 0 );
}

void   Help_listtopics( void )
{
}
*/

