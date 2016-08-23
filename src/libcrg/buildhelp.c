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

static char rcsid[] = "$Id: buildhelp.c,v 1.11 2000/11/02 20:06:47 welling Exp $";

#include   <stdio.h>
#include   <stdlib.h>
#include   <stdarg.h>
#include   <string.h>
#include   <ctype.h>

#include   <errors.h>
#include   "stdcrg.h"
#include   "slist.h"
#include   "helphelp.h"
#include   <sys/stat.h>
#include   <time.h>

/*  Usage:  buildhelp [-d] [-i indent] [-t topic] helpfile progname                           */
/*     -d        makes generated file progname_help.c a driver file (i.e., includes a main()) */
/*     -i indent Indent each line of help text by indent spaces (indent >= 0).  The default   */
/*               is 1 which thus includes the leading space on each line of help text.        */
/*     -t topic  makes topic the argument to Help() in the driver (default "selecttopic")     */
/*               This is only relevant if -d is supplied.                                     */
/*                                                                                            */
/*  Topic names are flagged by a * in the first column and are white-space delimited          */

#define      MAX_LINE      1024
#define      WIDTH           80

static char  *reserved[] = { "helphelp", "alltopics", "selecttopic", "listtopics", "htmldoc", "quithelp" };
static int   num_reserved = 6;
    
typedef struct string_struct {
  char* buf;
  int len;
  int last;
} String;

static char* exename= NULL;

static char      *myHelp_Intro, *myHelp_Help;

static HelpTable* hT= NULL;

static String* string_create()
{
  String* result;
  if (!(result= (String*)malloc(sizeof(String)))) {
    Abort("string_create: unable to allocate %d bytes!\n",sizeof(String));
  }
  if (!(result->buf= (char*)malloc(MAX_LINE))) {
    Abort("string_create: unable to allocate %d bytes!\n",MAX_LINE);
  }
  result->buf[0]= '\0';
  result->last= 0;
  result->len= MAX_LINE;
  return result;
}

static void string_destroy(String* s)
{
  free(s->buf);
  free(s);
}

static void string_add( String* s, int c )
{
  if (s->last>=s->len-1) {
    char* newbuf= (char*)malloc(2*s->len);
    if (!newbuf) {
      Abort("string_add: unable to allocate %d bytes!\n",2*s->len);
    }
    strcpy(newbuf,s->buf);
    free(s->buf);
    s->buf= newbuf;
    s->len= 2*s->len;
  }
  s->buf[s->last++]= (char)c;
  s->buf[s->last]= '\0';
}

static void string_clear(String* s) 
{
  s->last= 0;
  s->buf[s->last]= '\0';
}

static String* string_dup(char* s)
{
  String* result;
  if (!(result= (String*)malloc(sizeof(String)))) {
    Abort("create_string: unable to allocate %d bytes!\n",sizeof(String));
  }
  if (!(result->buf= strdup(s))) {
    Abort("create_string: unable to allocate %d bytes!\n",strlen(s)+1);
  }
  result->len= strlen(s);
  return result;
}

int main( int argc, char **argv )
{
    char     *line, *p;
    char     *helpfile, *progname;
    char     *dtopic = NULL;
    String*  tbuf= NULL;
    
    int      i, j, k, c;
    int      argnum = 0;
    int      driver = 0, indent = 1;
    int      linelen, linenum;
    int      stdin_flg= 0;

    FILE     *hfp, *ofp;

    struct stat  sbuf;

    exename= argv[0];

    hT= help_tables_create();

    if( (argc == 2) && !strcmp( argv[1], "-help" ) )
        Help( "selecttopic" );
    else if( argc < 3 )
        Usage( "%s [-d] [-i indent] [-t topic] helpfile progname\n"
               "       %s [-help [topic]]\n", argv[0], argv[0] );
    else if( !strcmp( argv[1], "-help" ) )
        Help( argv[2] );

    for( i = 1; i < argc; i++ )
    {
        if( argv[i][0] == '-' )
        {
	  if (argv[i][1] == '\0') {
	    if (i==argc-2) { /* read standard input */
	      stdin_flg= 1;
	      helpfile= "standard input";
	      argnum++;
	    }
	    else {
	      Abort("%s: - for stdin is out of place\n",argv[0]);
	    }
	  }
	  else {
            switch( argv[i][1] )
            {
              case 'd':

                driver = 1;
                break;

              case 'i':

                if( argv[i][2] == '\0' )
                    indent = atoi( argv[++i] );
                else
                    indent = atoi( &argv[i][2] );

                break;

              case 't':

                if( argv[i][2] == '\0' )
                    dtopic = argv[++i];
                else
                    dtopic = &argv[i][2];

                break;

              default:

                Warning( 1, "Unrecognized flag %s, continuing...", argv[i] );
                break;
            }
	  }
        }
        else
        {
            if( !argnum )
            {
                helpfile = argv[i];
                argnum++;
            }
            else if( argnum == 1 )
            {
                progname = argv[i];
                argnum++;
            }
        }
    }

    /* Make sure we have at least enough space to store entire help file */
    /* This is conservative, but might actually be needed in some cases  */
    
    if( !stdin_flg && (i = stat( helpfile, &sbuf )) )
        SysError( "stat", 0, "unable to stat file %s", helpfile );

    if (stdin_flg)
      line = (char *) ecalloc( MAX_LINE, sizeof(char) );
    else
      line = (char *) ecalloc( MAX_LINE + sbuf.st_size, sizeof(char) );

    /* Open Required Files: helpfile and  prognameHelp.c */

    sprintf( line, "%sHelp.c", progname );

    ofp = efopen( line, "w" );
    hfp = (stdin_flg ? stdin : efopen( helpfile, "r" ));

    /* Output header information */

    fprintf( ofp, "#include <stdio.h>\n" );
    fprintf( ofp, "#include <stdlib.h>\n" );
    fprintf( ofp, "#include <string.h>\n" );
    fprintf( ofp, "#include \"stdcrg.h\"\n" );
    fprintf( ofp, "#include \"helphelp.h\"\n" );
    fprintf( ofp, "\n" );
    fprintf( ofp, "static   int      Help_Initialized_ = 0;\n" );
    fprintf( ofp, "\n" );
    fprintf( ofp, "extern   int      Help_Num_Topics;\n" );
    fprintf( ofp, "extern   char    *Help_Progname;\n");
    fprintf( ofp, "extern   char    *Help_Creation_Time;\n");
    fprintf( ofp, "extern   char    *Help_Intro;\n" );
    fprintf( ofp, "extern   char    *Help_Topics[];\n" );
    fprintf( ofp, "extern   char    *Help_Text[];\n" );
    fprintf( ofp, "extern   HelpTable *Help_Table;\n" );
    fprintf( ofp, "\n" );
    fprintf( ofp, "void     Help_init( void );\n" );    
    fprintf( ofp, "void     Help_listtopics( void );\n" );    
    fprintf( ofp, "\n\n\n" );

    /* If making a driver routine, output appropriate main() */

    if( driver )
    {
        fprintf( ofp, "void     main( void )\n" );    
        fprintf( ofp, "{\n" );

        if( dtopic != NULL )
            fprintf( ofp, "    Help( \"%s\" );\n", dtopic );
        else
            fprintf( ofp, "    Help( \"%s\" );\n", "selecttopic" );

        fprintf( ofp, "}\n" );
        fprintf( ofp, "\n\n" );
    }

    /* Set message for help on help and introductory burst  */
    myHelp_Help = strdup(
                     "\\nHelp can be invoked either interactively or via a function call.\\n\\n"
                     "The interactive help system accepts commands and topic names at the prompt.\\n"
                     "The topic names are not sensitive to case, and only as many characters are\\n"
                     "required as are necessary to uniquely select a topic.  If a topic name is\\n"
                     "ambiguous, a list of possible completions will be supplied.  Topic names\\n"
                     "should generally not end in * because that character is used as a completion\\n"
                     "character in the interactive system.  If a partial topic name ends in *,\\n"
                     "then the help text associated with all topics whose initial substrings match\\n"
                     "the given string (not including the *) is output.  This is particularly useful\\n"
                     "for managing subtopics (e.g., by having topic names foo, foo:A, foo:B, foo:C\\n"
                     "that can be accessed separately or together with the name foo*).\\n"
                     "Note that if no completion character is used, then exact matches take precedence\\n"
                     "over partial matches.  Giving the topic name \\\"helphelp\\\" results in this\\n"
                     "message.\\n\\n"
                     "The available commands are \\\"listtopics\\\" and \\\"quithelp\\\" which list the\\n"
                     "available topics and quit the help system, respectively.  Commands are also\\n"
                     "case insensitive, and they require only as many characters as necessary to\\n"
                     "uniquely define them.\\n\\n"
                     "Invoking help by calling the function Help() requires an argument string\\n"
                     "containing one or more space separated topic names or commands.  The available\\n"
                     "commands are \\\"alltopics\\\", \\\"selecttopic\\\", \\\"listtopics\\\",  and \\\"htmldoc\\\".\\n"
                     "The first outputs the text for all topics, the second engages interactive mode,\\n"
                     "and the third provides a formatted list of available topics (not including\\n"
                     "the built-in topic \\\"helphelp\\\".  Topic completion is also available in\\n"
                     "invokation mode; it is used as described above.\\n\\n"
                      );
    
    sprintf( line, "Interactive Help System for %s.  Select topics at prompt.\\n"
                   "Use quithelp to quit system.\\n", progname );
    myHelp_Intro = strdup( line );

    /* Process help file: build list of topics and text     */

    tbuf= string_create();

    c = fgetc( hfp );
    linenum = 1;

    /* Skip to first asterisk */
    while( c != '*' && c != EOF )
    {
        c = fgetc( hfp );

        if( c == '\n' )
            linenum++;
    }

    while( c != EOF )
    {
        /* Get Topic Name */
            
        for( i = 0; c != '\n' && c != EOF; i++ )
            c = line[i] = fgetc( hfp );

        line[i] = '\0';

        if( c == EOF )
        {
            Warning( 1, "Unexpected end of file amidst topic definition in %s, line %d.\n", 
                    helpfile, linenum );
            break;
        }

        linenum++;
        p = strtok( line, " \t\n" );

        if( p == NULL )
            Abort( "Empty topic name in definition in %s, line %d.\n", helpfile, linenum-1 );

        /*
         * Get Text
         *
         *  Grab all text until next '*' at the beginning of a line.
         *  Only a few characters get special treatment:
         *
             \n            Converted to '\\' 'n' for printing, then followed by appropriate indentation
                           and search for stars
             \t\v\b\r\f\a  Converted to '\\' x form for printing.
             \             Converted to '\\' '\\' form for printing.
             ",'           Converted to '\\' x form for printing.              
             EOF           Finish up
         *
         */

	string_clear(tbuf);
        for( i = 0; c != EOF; )
        {
            c = fgetc( hfp );

            if( c == EOF )
                break;
            else if( strchr( "\n\t\v\b\r\f\a\\\"\'", c ) )    /* Special character */
            {
                if( c == '\\' )
                {
                    string_add(tbuf, '\\');
                    string_add(tbuf, c);
                }
                else if( c == '\"' )
                {
                    string_add(tbuf, '\\');
                    string_add(tbuf, c);
                }
                else if( c == '\'' )
                {
                    string_add(tbuf, '\\');
                    string_add(tbuf, c);
                }
                else if( c == '\t' )
                {
                    string_add(tbuf, '\\');
                    string_add(tbuf, 't');
                }
                else if( c == '\v' )
                {
                    string_add(tbuf, '\\');
                    string_add(tbuf, 'v');
                }
                else if( c == '\b' )
                {
                    string_add(tbuf, '\\');
                    string_add(tbuf, 'b');
                }
                else if( c == '\r' )
                {
                    string_add(tbuf, '\\');
                    string_add(tbuf, 'r');
                }
                else if( c == '\f' )
                {
                    string_add(tbuf, '\\');
                    string_add(tbuf, 'f');
                }
                else if( c == '\a' )
                {
                    string_add(tbuf, '\\');
                    string_add(tbuf, 'a');
                }
                else if( c == '\n' )
                {
                    while( c == '\n' )
                    {
                        string_add(tbuf, '\\');
                        string_add(tbuf, 'n');
                        linenum++;
                    
                        c = fgetc( hfp );       

                        if( c == EOF || c == '*' ) /* End of topic and/or file. Continue  */
                            break; 
                        else if( c == ' ' || c == '\t' ) /* Initial TAB treated just like space */
                        {
                            for( j = 0; j < indent; j++ )
                                string_add(tbuf, ' ');
                        }
                        else if( c == '\n' ) /* Let blank lines slide               */
                            continue;
                        else
                            Abort( "Line %d of %s should begin with a space or *.", linenum, helpfile );
                    }

                    if( c == EOF || c == '*' )
                        break;
                }
            }
            else
                string_add(tbuf, c);
        }
	string_add(tbuf,'\0');

	/* c should be either EOF or '*' (at beginning of line) here */
	help_add_to_tables(p, tbuf->buf, hT);

    }

    string_destroy(tbuf);

    /* Add reserved names */

    for( i = 0; i < num_reserved; i++ )
    {
      help_add_to_tables(reserved[i], NULL, hT);
    }

    if (!stdin_flg) efclose( hfp );

    /* Set Help_init (including Help_Intro)                 */

    fprintf( ofp, "void     Help_init( void )\n" );    
    fprintf( ofp, "{\n" );
    fprintf( ofp, "    if( !Help_Initialized_ )\n" );
    fprintf( ofp, "    {\n" );
    fprintf( ofp, "        Help_Initialized_ = 1;\n" );
    {
      char* progname_tail;
      progname_tail= progname + strlen(progname) - 1;
      while ( progname_tail > progname ) {
	if (*(progname_tail-1) == '/') break;
	progname_tail--;
      }
      fprintf( ofp, "        Help_Progname = \"%s\";\n", progname_tail );
    }
    {
      time_t tval;
      char* t_cptr;
      char t_cptr_no_newline[32];
      tval= time(NULL);
      t_cptr= ctime(&tval);
      strncpy(t_cptr_no_newline, t_cptr, 26); /* 26 is std length */
      t_cptr_no_newline[24]= '\0'; /* remove newline */
      fprintf( ofp, "        Help_Creation_Time = \"%s\";\n",
	       t_cptr_no_newline );
    }
    fprintf( ofp, "\n" );

    fprintf( ofp, "        Help_Intro = strdup( \"%s\" );\n", myHelp_Intro );
    fprintf( ofp, "\n" );

    j= 0;
    for( i = 0; i < hT->topicCount; i++ )
    {
      slist_totop(hT->topicTable[i]);
      slist_totop(hT->textTable[i]);
      do {
	char* s= slist_next(hT->topicTable[i]);
	char* t= slist_next(hT->textTable[i]);
	if (s) 
	  fprintf( ofp, "        Help_Topics[%d] = strdup( \"%s:%s\" );\n", 
		   j, hT->topicNames[i],s );
	else 
	  fprintf( ofp, "        Help_Topics[%d] = strdup( \"%s\" );\n", 
		   j, hT->topicNames[i] );
	if (t)
	  fprintf( ofp, "        Help_Text[%d]   = strdup( \"%s\" );\n", 
		   j, t );
	else fprintf( ofp, "        Help_Text[%d]   = NULL;\n", j );
	j++;
      }
      while (!slist_atend(hT->topicTable[i]));
    }
    fprintf( ofp, "\n" );
    fprintf( ofp, "        Help_Num_Topics = %d;\n", j );
    
    fprintf( ofp, "        Help_Table= help_load_topics(Help_Num_Topics,\n" );
    fprintf( ofp, "                                     Help_Topics,\n" );
    fprintf( ofp, "                                     Help_Text);\n" );

    fprintf( ofp, "    }\n" );
    fprintf( ofp, "}\n" );    
    fprintf( ofp, "\n\n" );

    /* Set Help_listtopics() with prettified listing */

    fprintf( ofp, "void     Help_listtopics( void )\n" );    
    fprintf( ofp, "{\n" );
    fprintf( ofp, "    help_emit_topics( Help_Table, stdout );\n");

    fprintf( ofp, "}\n" );    
    fprintf( ofp, "\n\n" );

    /* Done */

    efclose( ofp );
    exit(0); /* Suns seem to want this */
    return 0;
}

