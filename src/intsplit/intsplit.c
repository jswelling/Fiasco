/************************************************************
 *                                                          *
 *  intsplit.c                                              *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1995 Department of Statistics,         *
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
 *  Original programming by Mark Fitzgerald  6-96           *
 ************************************************************/
/*************************************************************

  DESCRIPTION OF INTSPLIT.C

  intsplit takes a file with information on experimental conditions
  and converts it to an image-by-image condition specification for use
  by anova and spm. 

**************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"
#include "slist.h"

static char rcsid[] = "$Id: intsplit.c,v 1.20 2007/03/21 23:52:35 welling Exp $";

#define SAFE_MALLOC( val, type, count ) \
  if (!(val=(type*)malloc(count*sizeof(type)))) \
    Abort("%s:%s: unable to allocate %d bytes!\n",\
          __FILE__,__LINE__,count*sizeof(type))

typedef struct factor_info_struct {
  char* name;
  KVHash* levels;
  int nLevels; /* for convenience */
} FactorInfo;

typedef struct condition_info_struct {
  KVHash* conditions;
  long nConditions;
} ConditionInfo;

char *progname;
static SList *factors= NULL;
static ConditionInfo* condInfo= NULL;
#ifdef never
static KVHash* conditions;
#endif

static void to_from(const char* str, long* first, long* last, long max, 
		    long linenum);

/* Returns factor number, defining new factor of necessary */
static int facinfo_getFactorLevelID( FactorInfo* fi, const char* s )
{
  int thisFactorLevelID;
  if ( kvLookup(fi->levels,s) ) 
    thisFactorLevelID= kvGetInt(fi->levels, s);
  else {
    thisFactorLevelID= fi->nLevels;
    kvDefInt(fi->levels,s,fi->nLevels++);
  }
  return thisFactorLevelID;
}

/* Creates FactorInfo struct */
static FactorInfo* facinfo_create(const char* name)
{
  FactorInfo* result= NULL;

  SAFE_MALLOC(result,FactorInfo,1);
  result->name= strdup(name);
  result->levels= kvFactory(7); /* a nice small hash table */
  result->nLevels= 0;

  (void)facinfo_getFactorLevelID( result,"NA" );

  return result;
}

/* Returns condition number, defining new condition of necessary */
static int condinfo_getCond( ConditionInfo* ci, const char* s )
{
  int thisConID;
  if ( kvLookup(ci->conditions,s) ) 
    thisConID= kvGetInt(ci->conditions, s);
  else {
    thisConID= ci->nConditions;
    kvDefInt(ci->conditions,s,ci->nConditions++);
  }
  return thisConID;
}

/* Creates ConditionInfo struct */
static ConditionInfo* condinfo_create(int nFactors)
{
  ConditionInfo* result= NULL;
  char condName[512];
  int i;

  SAFE_MALLOC(result,ConditionInfo,1);
  result->conditions= kvFactory(KV_DEFAULT_SIZE);
  result->nConditions= 0;

  strcpy( condName, "NA" );
  for ( i=1; i<nFactors; i++) strncat(condName," NA",sizeof(condName)-1);
  condName[sizeof(condName)-1]= '\0';
  (void)condinfo_getCond(result,condName);

  return result;
}

int main( int argc, char* argv[] ) 
{

  MRI_Dataset *Input = NULL;
  char codefile[512], infile[512], sinfile[512], soutfile[512];
  char condName[512];
  long dt=0, dz=0;
  FILE *ifp = NULL, *ofp = NULL;
  char scanline[512];
  unsigned char **missing = NULL;
  long **condnums = NULL;
  int thisConID;
  int total_conds;
  char cbynumimages[512];
  int num_factors, bynumimages;
  long firstslice, lastslice, firstimage, lastimage;
  int linenum, linelen, m, t, z;
  int i;
  int n_unspecified= 0;
  int isMissing;
  char* here;
  KVIterator* kvi;
  char* tmpp;
  const char** sortbuf;

  progname= argv[0];

  /* Print version number */
  Message( "# %s\n", rcsid );

  /* Check to see if help was requested */
  if( ( argc > 1 ) && !strcmp( argv[1], "-help" ) )
    {
      if( argc == 2 )
	Help( "selecttopic" );
      else
	Help( argv[2] );
    }


  /*** Parse command line ***/

  cl_scan( argc, argv );

  /* Deprecate old options */
  if (cl_present( "c" ))
     Abort ("Option c no longer exists.  Please see help file.\n");
  if (cl_present( "input|i" ))
     Abort ("Option input|i has been replaced by infile outfile format.  Please see help file.\n");
  if (cl_present( "spliti" ))
     Abort ("Option spliti has been replaced by splitin|spi.  Please see help file.\n");
  if (cl_present( "splito" ))
     Abort ("Option splito has been replaced by splitout|spo.  Please see help file.\n");
  if (cl_present( "max_cpf|m" ))
     Abort ("Option max_cpf|m no longer exists.  Please see help file.\n");

  /* Get filenames */
  cl_get( "code", "%option %s[%]", "codes", codefile );
  cl_get( "splitin|spi", "%option %s[%]", "split", sinfile );
  cl_get( "splitout|spo", "%option %s[%]", "newsplit", soutfile );

  if(!cl_get("", "%s", infile)) {
    fprintf(stderr, "%s: Input file name not given.\n", argv[0]);
    exit(-1);
  }
  if (cl_cleanup_check()) {
    int i;
    fprintf(stderr,"%s: invalid argument in command line:\n    ",argv[0]);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(-1);
  }

  /*** End command-line parsing ***/

  /* Get number of images and slices and missing info from header */
  Input = mri_open_dataset( infile, MRI_READ );
  if( mri_has( Input, "images.extent.t" ) &&
      mri_has( Input, "images.extent.z" ) ) {
    dt = mri_get_int( Input, "images.extent.t" );
    dz = mri_get_int( Input, "images.extent.z" );
  }
  else if (mri_has( Input, "samples.extent.t" ) &&
      mri_has( Input, "samples.extent.z" )) {
    dt = mri_get_int( Input, "samples.extent.t" );
    dz = mri_get_int( Input, "samples.extent.z" );
  }
  else
    Abort( "Header-file [%s] does not have extent keys for z and t.\n", 
	   infile );
  missing = get_missing( Input );
  mri_close_dataset( Input );

  /* Initialize data structures */
  factors= slist_create();
  condnums = Matrix( dt, dz, long );
  if (!condnums || !(condnums[0]))
    Abort("%s: unable to allocate %d by %d array of condition numbers!\n",
	  argv[0],dt,dz);
  total_conds= 0;

  /* Open input files */
  ifp = efopen( sinfile, "r" );

  /* Read and assess set-up (first) line of input file */
  fscanf( ifp, "%510[^\n]%*[\n]", scanline );
  m = sscanf( scanline, "%ld %510s", &num_factors, cbynumimages );
  if( m <= 0 )
    Abort( "Couldn't read set-up line of %s.\n", sinfile );
  if ( num_factors <= 0 )
    Abort( "Number of factors must be greater than 0! [%ld]\n.", num_factors );
  if( ( m > 1 ) && !strcasecmp( cbynumimages, "bynumimages" ) )
    {
      bynumimages = 1;
      firstslice = 0;
      lastslice = dz;
    }
  else
    bynumimages = 0;

  condInfo= condinfo_create(num_factors);

  /* Read factor name (second) line of input file */
  fscanf( ifp, "%510[^\n]%*[\n]", scanline );
  here= strtok_r(scanline," \t",&tmpp);
  while (here) {
    FactorInfo* fi= facinfo_create(here);
    slist_append(factors, fi);
    here= strtok_r(NULL," \t",&tmpp);
  }
  
  /* Initialize experimental conditions to something impossible */
  for( t = 0; t < dt; t++ )
    for( z = 0; z < dz; z++ )
      condnums[t][z] = -1;    

  /* Read and process the rest of input file */
  linenum = lastimage = 0;

  while( !feof( ifp ) ) {
    /* Read a line */
    fscanf( ifp, "%510[^\n]%*[\n]", scanline );

    linenum++;
    isMissing= 0;
    condName[0]= '\0';
    here= strtok_r(scanline," \t",&tmpp);
    if (!here) continue; /* blank line */

    /* Build the condition name, noting missing */
    slist_totop(factors);
    while (!slist_atend(factors)) {
      FactorInfo* fi= (FactorInfo*)slist_next(factors);
      if (!here) 
	Abort( "%s: Line #%ld of %s is too short!\n", 
               progname, linenum, sinfile);
      if (!strcasecmp( here, "missing" ) || !strcmp( here, "NA" ))
	isMissing= 1;
      (void)facinfo_getFactorLevelID(fi,here); /* define factor if necessary */
      if (strlen(here)+strlen(condName)>sizeof(condName)-1)
	Abort("%s: condition string too long on split line %d of %s!\n",
	      progname, linenum, sinfile);
      strncat(condName,here,sizeof(condName)-1);
      strcat(condName," ");
      here= strtok_r(NULL," \t",&tmpp);
    }
    condName[strlen(condName)-1]= '\0'; /* trim off trailing blank */

    /* Determine condition number, defining the condition if necessary */
    if (isMissing) thisConID= 0;
    else thisConID= condinfo_getCond(condInfo,condName);
    
    /* Determine image and slice numbers for condition */
    if( bynumimages ) {
      if (!here)
	Abort( "%s: Line #%ld of %s is too short!\n", 
	       progname, linenum, sinfile);
      if (*here=='#') /* comment field */
	Abort( "%s: Line #%ld of %s is too short (excluding comment)!\n",
	       progname, linenum, sinfile);
      firstimage = lastimage;
      if( firstimage >= dt )
	{
	  Abort("Ran out of images before condition file ended.\n");
	}
      /* Test for a common user error */
      if (strchr(here,'-')) 
	Abort("%s: Found a range <%s> where an integer was expected!\n",
	      progname, here);
      lastimage += atol( here );
      if( lastimage > dt ) lastimage = dt;
    }
    else {
      if (!here)
	Abort( "%s: Line #%ld of %s is too short!\n", 
	       progname, linenum, sinfile);
      if (*here=='#') /* comment field */
	Abort( "%s: Line #%ld of %s is too short (excluding comment)!\n",
	       progname, linenum, sinfile);
      to_from( here, &firstimage, &lastimage, dt, linenum );
      here= strtok_r(NULL," \t",&tmpp);
      if( here && *here != '#')
	to_from( here, &firstslice, &lastslice, dz, linenum );
      else {
	firstslice = 0;
	lastslice = dz;
      }
    }
    
    /* Set experimental condition to appropriate image/slice numbers */
    for( t = firstimage; t < lastimage; t++ )
      for( z = firstslice; z < lastslice; z++ ) {
	condnums[t][z] = ( missing[t][z] )? 0: thisConID;
      }
  }

  /* Done reading */
  efclose( ifp );

  /* Count number of unspecified conditions.  Set them to missing,
   * and issue a warning if appropriate.
   */
  n_unspecified= 0;
  for (t=0; t<dt; t++)
    for (z=0; z<dz; z++) {
      if (condnums[t][z]==-1) {
	n_unspecified++;
	condnums[t][z]= 0; /* missing */
      }
    }
  if (n_unspecified>0) 
    Warning(1,"%s: split file left the conditions of %d slices unspecified!\n",
	    argv[0], n_unspecified);

  /* Write the output condition file */
  ofp = efopen( soutfile, "w" );
  for( t = 0; t < dt; t++ )
    for( z = 0; z < dz; z++ )
      fprintf( ofp, "%7ld %7ld %8ld\n", t, z, condnums[t][z] );
  efclose( ofp );

  /* Write the codings for the experimental conditions, in numerical order.
   */
  SAFE_MALLOC(sortbuf,const char*,condInfo->nConditions);
  kvi= kvUniqueIteratorFactory(condInfo->conditions);
  while (kvIteratorHasMorePairs(kvi)) {
    KVPair* p= kvIteratorNextPair(kvi);
    sortbuf[kvGetInt(condInfo->conditions,kvKey(p))]= kvKey(p);
  }
  kvDestroyIterator(kvi);
  ofp = efopen( codefile, "w" );
  fprintf( ofp, "          " );
  slist_totop( factors );
  while (!slist_atend(factors)) {
    FactorInfo* fi= (FactorInfo*)slist_next(factors);
    fprintf( ofp, "%s ", fi->name );
  }
  fprintf( ofp, "\n" );
  for( m = 0; m < condInfo->nConditions; m++ )
    fprintf( ofp, "%8ld  %s\n", m, sortbuf[m] );
  efclose( ofp );
  free(sortbuf);

#ifdef never
  /* A little diagnostic output about factor levels */
  slist_totop(factors);
  while (!slist_atend(factors)) {
    FactorInfo* fi= (FactorInfo*)slist_next(factors);
    fprintf(stderr,"Factor <%s> has the following levels:\n",fi->name);
    kvDumpTableUnique( fi->levels, stdout, fi->name, 1, 0 );
  }
#endif

  Message( "#      Experimental conditions file created.\n" );	  
  return 0;
}
    

/* Function to convert image/slice numbering string to     */
/*   integers indicating first and last number of sequence */
static void to_from(const char* str, long* first, long* last, long max,
		    long linenum)
{
  char cfirst[64], clast[64];
  long ll;

  /* Read string, checking for dash */
  ll = sscanf( str, "%64[^-]%*[-]%64s", cfirst, clast );
  
  if( ll < 1 )
    Abort( "to_from subroutine: empty string passed." );

  if( ll == 1 )
    {
      /* No dash, just a single image/slice number */
      if( !strcasecmp( cfirst, "all" ) )
	{
	  *first = 0;
	  *last = max;
	}
      else
	{
	  *first = atol( cfirst );
	  *last = *first + 1;
	  if( *last > max )
	    {
	      Warning( 1,
		       "Value %d out of bounds; line near %d ignored.\n",
		       *first, linenum );
	      *last = *first;
	    }
	}
    }
  else
    {
      /* A dash is present --- a sequence is indicated */
      *first = atol( cfirst );
      *last = atol( clast ) + 1;
      if( *last > max )
	{
	  Warning( 1,
		   "Last number of sequence <%s> out of bounds near line %d: set to boundary %d.\n",
		   str, linenum, max );
	  *last = max;
	}
    }

  if( *first > *last )
    {
      Warning( 1,
	       "First number was larger than last number in sequence: ignoring near line %d.\n",
	       linenum);
      *first = 0;
      *last = 0;
    }
	       
  if( *first < 0 )
    {
      Warning( 1,
	       "First number of a sequence was not positive: set to zero near line %d.\n",
	       linenum);
      *first = 0;
    }

  return;

}

