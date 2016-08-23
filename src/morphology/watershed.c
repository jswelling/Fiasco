/************************************************************
 *                                                          *
 *  watershed.c                                             *
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
 *  Original programming by Mark Fitzgerald  2-95           *
 *     5-96: Pittsburgh Format, Mark Fitzgerald             *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <assert.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"
#include "slist.h"

static char rcsid[] = "$Id: watershed.c,v 1.7 2007/03/21 23:58:25 welling Exp $";

/* Default algorithm */
#define DEFAULT_ALG_STRING "nomarkboundaries,nofatboundaries,floodmaxima"

/* A label to designate unlabeled voxels.  It *must* be zero. */
#define LBL_EMPTY 0
/* A label for boundary voxels.  */
#define LBL_BOUNDARY -1

/* Initial size of region table */
#define INITIAL_RGN_TBL_SIZE 64

/* Access for 3D arrays.  If you change these, remember to change
 * indicesFromPOffset() below!
 */
#define LOC(matrix,x,y,z) matrix[((((z)*dy)+(y))*dx)+(x)]
#define PLOC(matrix,x,y,z) matrix[((((z+1)*(dy+2))+(y+1))*(dx+2))+(x+1)]

typedef struct algorithm_struct {
  int markBoundaries;
  int fatBoundaries;
  int alwaysFlood;
  int breakTies;
} Algorithm;

typedef struct live_voxel_struct {
  double val;
  long offset;
} LiveVoxel;

typedef struct region_struct {
  double max_val;
  long max_offset;
  long i, j, k;
  long nvoxels;
} Region;

static int debug_flag= 0;
static int verbose_flag= 0;
static char* progname;
static Region* regionTable= NULL;
static long regionTableSize= 0;

static void growRegionTable()
{
  if (regionTableSize==0) {
    if (!(regionTable= (Region*)malloc(INITIAL_RGN_TBL_SIZE*sizeof(Region))))
      Abort("%s: unable to allocate %d bytes!\n",progname,
	    INITIAL_RGN_TBL_SIZE*sizeof(Region));
    regionTableSize= INITIAL_RGN_TBL_SIZE;
  }
  else {
    if (!(regionTable= (Region*)realloc(regionTable,
					2*regionTableSize*sizeof(Region))))
      Abort("%s: unable to realloc %d bytes!\n",progname,
	    2*regionTableSize*sizeof(Region));
    regionTableSize *= 2;
  }
  if (debug_flag) fprintf(stderr,"Region table size is now %ld\n",
			  regionTableSize);
}

static void indicesFromPOffset(long offset, long* i, long* j, long* k,
			       long dx, long dy, long dz)
{
  long remainder= offset;
  long chunk= remainder/(dx+2);
  double* dummy;
  *i= remainder-((chunk*(dx+2))+1);
  remainder= chunk;
  chunk= remainder/(dy+2);
  *j= remainder-((chunk*(dy+2))+1);
  *k= chunk-1;
  
  assert(offset==&PLOC(dummy,*i,*j,*k)-dummy);
}

static void checkFormat( MRI_Dataset* ds, char* fname )
{
  if( !mri_has( ds, "images" ) )
    Abort( "%s operates only on standard images; %s won't work", 
	   progname, fname );
  if( mri_has( ds, "images.dimensions" ) )
    {
      const char* dimstr= mri_get_string(ds,"images.dimensions");
      if (!strcmp(dimstr,"xyz")) {
	/* Everything is fine */
      }
      else if (!strcmp(dimstr,"vxyz")) {
	if (!mri_has(ds,"images.extent.v") || mri_get_int(ds,"images.extent.v")!=1)
	  Abort("%s: input dataset %s has invalid v extent, or v extent != 1.\n",
		progname,fname);
      }
      else if (!strcmp(dimstr,"vxyzt")) {
	if (!mri_has(ds,"images.extent.v") || mri_get_int(ds,"images.extent.v")!=1)
	  Abort("%s: input dataset %s has invalid v extent, or v extent != 1.\n",
		progname,fname);
	if (!mri_has(ds,"images.extent.t") || mri_get_int(ds,"images.extent.t")!=1)
	  Abort("%s: input dataset %s has invalid t extent, or t extent != 1.\n",
		progname,fname);
      }
      else if (!strcmp(dimstr,"xyzt")) {
	if (!mri_has(ds,"images.extent.t") || mri_get_int(ds,"images.extent.t")!=1)
	  Abort("%s: input dataset %s has invalid t extent, or t extent != 1.\n",
		progname,fname);
      }
      else 
	Abort("%s: input dataset %s must have dimensions (v)xyz(t)!\n",
	      progname,fname);
    }
  else
    Abort( "%s: %s does not have the images.dimensions key.", 
	   progname,fname);
  if (!mri_has(ds,"images.extent.x"))
    Abort("%s: %s is missing images.extent.x tag!\n",progname,fname);
  if (!mri_has(ds,"images.extent.y"))
    Abort("%s: %s is missing images.extent.y tag!\n",progname,fname);
  if (!mri_has(ds,"images.extent.z"))
    Abort("%s: %s is missing images.extent.z tag!\n",progname,fname);
}

static long* buildNeighborOffsetTable3D(int *nNbrs, long dx, long dy, long dz)
{
  long* nbrOffsets= NULL;
  long loop;
  long i,j,k;
  long offset;
  double* dummy;

  *nNbrs= (3*3*3)-1; /* because we are in 3D */

  if (!(nbrOffsets= (long*)malloc((*nNbrs)*sizeof(long))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,(*nNbrs)*sizeof(long));

  /* Initialize the neighbor table.  We just use dvol because it is convenient;
   * this doesn't actually access any memory locations. 
   */
  loop= 0;
  for (i=-1;i<=1;i++)
    for (j=-1;j<=1;j++)
      for (k=-1;k<=1;k++) {
	offset= &PLOC(dummy,i,j,k)-&PLOC(dummy,0,0,0);
	if (offset!=0) nbrOffsets[loop++]= offset;
      }

  if (debug_flag) {
    fprintf(stderr,"Offset table follows:\n");
    for (i=0; i<loop; i++) fprintf(stderr,"%ld: %ld\n",i,nbrOffsets[i]);
  }

  assert(loop==(*nNbrs));

  return nbrOffsets;
}

static double* allocateAndClearPaddedDouble(long dx, long dy, long dz)
{
  double* result= NULL;
  double* runner= NULL;
  double* lim= NULL;
  if (!(result= (double*)malloc((dx+2)*(dy+2)*(dz+2)*sizeof(double)))) 
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,(dx+2)*(dy+2)*(dz+2)*sizeof(double));
  runner= result;
  lim= result+((dx+2)*(dy+2)*(dz+2));
  while (runner<lim) *runner++= 0;
  return result;
}

static long* allocateAndClearPaddedLong(long dx, long dy, long dz)
{
  long* result= NULL;
  long* runner= NULL;
  long* lim= NULL;
  if (!(result= (long*)malloc((dx+2)*(dy+2)*(dz+2)*sizeof(long)))) 
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,(dx+2)*(dy+2)*(dz+2)*sizeof(long));
  runner= result;
  lim= result+((dx+2)*(dy+2)*(dz+2));
  while (runner<lim) *runner++= 0;
  return result;
}

static int compareLiveVoxels(const void* p1, const void* p2)
{
  /* Sort in descending order, breaking ties based on
   * proximity to the high corner of the grid.
   */
  LiveVoxel* v1= (LiveVoxel*)p1;
  LiveVoxel* v2= (LiveVoxel*)p2;
  if (v1->val<v2->val) return 1;
  if (v2->val<v1->val) return -1;
  return ((int)(v2->offset - v1->offset));
}

static long* allocOneLong(long offset)
{
  long* result= (long*)malloc(sizeof(long));
  if (!result) Abort("%s: unable to allocate 1 long!\n");
  *result= offset;
  return result;
}

static void floodFillConstantRegion(long* labelVol, const long* maskVol,
				    const double* dvol,
				    long seedOffset, long seedVal,
				    const long* nbrOffsets, int nNbrs, 
				    long fill)
{
  /* This is inefficient, but it should seldom be used */
  SList* stack= slist_create();
  long* oldOffsetPtr; 
  long count= 0;
  
  labelVol[seedOffset]= fill;
  count++;
  slist_push(stack,allocOneLong(seedOffset));
  while (!slist_empty(stack)) {
    int loop;
    oldOffsetPtr= (long*)slist_pop(stack);
    for (loop=0; loop<nNbrs; loop++) {
      long thisNbrOffset= nbrOffsets[loop];
      if (maskVol[*oldOffsetPtr + thisNbrOffset]
	  && labelVol[*oldOffsetPtr + thisNbrOffset]==LBL_EMPTY
	  && dvol[*oldOffsetPtr+thisNbrOffset]==seedVal) {
	labelVol[*oldOffsetPtr+thisNbrOffset]= fill;
	count++;
	slist_push(stack,allocOneLong(*oldOffsetPtr+thisNbrOffset));
      }
    }
    free(oldOffsetPtr);
  }
  slist_destroy(stack,NULL);
  if (debug_flag) 
    fprintf(stderr,"Flood filled %ld voxels having val %g with %ld\n",
	    count, seedVal, fill);
}

static long applyWatershed(long* labelVol, 
			   const double* dvol, const long* maskVol,
			   const long* nbrOffsets, long nNbrs,
			   long dx, long dy, long dz, const Algorithm* alg)
{
  long i, j, k;
  LiveVoxel* queue= NULL;
  LiveVoxel* queueRunner= NULL;
  long voxelsInQueue= 0;
  long regionCount= 0;
  long loop;

  /* Make a big buffer, insert all the points, and sort them. */
  if (!(queue= (LiveVoxel*)malloc(dx*dy*dz*sizeof(LiveVoxel))))
    Abort("%s: unable to allocate %d bytes!\n",progname,
	  dx*dy*dz*sizeof(LiveVoxel));
  queueRunner= queue;
  voxelsInQueue= 0;
  for (i=0; i<dx; i++)
    for (j=0; j<dy; j++)
      for (k=0; k<dz; k++) {
	if (PLOC(maskVol,i,j,k)) {
	  queueRunner->val= PLOC(dvol,i,j,k);
	  queueRunner->offset= &PLOC(dvol,i,j,k)-dvol;
	  queueRunner++;
	}
      }
  voxelsInQueue= queueRunner-queue;
  if (debug_flag) fprintf(stderr,"%ld voxels are within the mask\n",
			  voxelsInQueue);
  qsort(queue, voxelsInQueue, sizeof(LiveVoxel), compareLiveVoxels);

  regionCount= 0;
  for (queueRunner=queue; queueRunner<queue+voxelsInQueue; queueRunner++) {
    LiveVoxel thisVox= *queueRunner;
    if (labelVol[thisVox.offset]==LBL_EMPTY) {
      /* This is a maximum which hasn't yet been labeled.
       * There may be a constant region at this minimum which we must fill.
       */
      long i, j, k;
      Region* rgn= NULL;
      regionCount++;
      if (regionTableSize<regionCount) growRegionTable();
      floodFillConstantRegion(labelVol, maskVol, dvol, 
			      thisVox.offset, thisVox.val,
			      nbrOffsets, nNbrs, regionCount);
      indicesFromPOffset(thisVox.offset,&i,&j,&k,dx,dy,dz);
      rgn= &regionTable[regionCount-1];
      rgn->max_val= dvol[thisVox.offset];
      rgn->max_offset= thisVox.offset;
      rgn->i= i;
      rgn->j= j;
      rgn->k= k;
      rgn->nvoxels= 0;
    }
    /* bounds don't spread unless the algorithm explicitly says they do */
    if (alg->fatBoundaries || labelVol[thisVox.offset]!=LBL_BOUNDARY) { 
      for (loop=0; loop<nNbrs; loop++) {
	long thisNbrOffset= thisVox.offset-nbrOffsets[loop]; /*bckwrd offset!*/
	if (maskVol[thisNbrOffset]                    /* nbr is in mask */
	    && (dvol[thisNbrOffset]<=dvol[thisVox.offset])) /*nbr not uphill*/
	  {
	    if (labelVol[thisNbrOffset]==LBL_EMPTY) {
	      labelVol[thisNbrOffset]= labelVol[thisVox.offset];
	      if (alg->alwaysFlood)
		floodFillConstantRegion(labelVol, maskVol, dvol,
					thisNbrOffset, dvol[thisNbrOffset],
					nbrOffsets, nNbrs, 
					labelVol[thisNbrOffset]);
	    }
	    else if ((labelVol[thisNbrOffset] != labelVol[thisVox.offset])
		     && alg->markBoundaries)
	      labelVol[thisNbrOffset]= LBL_BOUNDARY;
	  }
      }
    }
  }

  free(queue);
  return regionCount;
}

static void parseAlgorithm( Algorithm* alg, const char* algstring )
{
  char* str= strdup(algstring);
  char* tok= NULL;
  char* here= NULL;

  tok= strtok_r(str, ",", &here);
  while (tok) {
    if (!strcasecmp(tok,"markboundaries")) alg->markBoundaries= 1;
    else if (!strcasecmp(tok,"nomarkboundaries")) alg->markBoundaries= 0;
    else if (!strcasecmp(tok,"fatboundaries")) alg->fatBoundaries= 1;
    else if (!strcasecmp(tok,"nofatboundaries")) alg->fatBoundaries= 0;
    else if (!strcasecmp(tok,"alwaysflood")) { 
      alg->alwaysFlood= 1;
      alg->breakTies= 0;
    }
    else if (!strcasecmp(tok,"floodmaxima")) {
      alg->alwaysFlood= 0;
      alg->breakTies= 1;
    }
    else Abort("%s: unrecognized token <%s> in algorithm string!\n",
	       progname, tok);
    tok= strtok_r(NULL, ",", &here);
  }

  free(str);
}

static void algorithmToString( const Algorithm* alg, char* buf, int bufsize )
{
  snprintf(buf,bufsize,"%s,%s,%s",
	   alg->markBoundaries?"markboundaries":"nomarkboundaries",
	   alg->fatBoundaries?"fatboundaries":"nofatboundaries",
	   (alg->alwaysFlood && !alg->breakTies)?"alwaysflood":"floodmaxima");
}

int main( argc, argv ) 
     int argc;
     char **argv;
{
  MRI_Dataset *Input = NULL, *Output = NULL, *Mask= NULL;
  char infile[512], outfile[512], maskfile[512], algstring[512];
  long dx, dy, dz;
  long *labelVol= NULL, *maskVol= NULL; /* padded */
  long *ioVol= NULL; /* not padded */
  double *dvol= NULL; /* padded */
  long i, j, k, loop;
  long* buf= NULL;
  double* dbuf= NULL;
  long* nbrOffsets= NULL;
  int nNbrs= 0;
  int use_mask;
  long regionCount= 0;
  Algorithm alg;

  progname= argv[0];

  /*** Initialize algorithm ***/
  parseAlgorithm( &alg, DEFAULT_ALG_STRING );

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

  /* Get filenames */
  use_mask= cl_get( "mask|m", "%option %s", maskfile );
  debug_flag= cl_present("debug");
  verbose_flag= cl_present("v|verbose");
  if (cl_get("algorithm|alg", "%option %s", algstring))
    parseAlgorithm( &alg, algstring );
  if (!cl_get("","%s",infile)) {
    fprintf(stderr,"%s: required input filename not given.\n",progname);
    Help("usage");
    exit(-1);
  }
  if (!cl_get("","%s",outfile)) {
    fprintf(stderr,"%s: required output filename not given.\n",progname);
    Help("usage");
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

  /* Print version number */
  if (verbose_flag) Message( "# %s\n", rcsid );

  /* Open input dataset */
  if ( !strcmp( infile, outfile ) || !strcmp( infile, maskfile )
      || !strcmp( outfile, maskfile ) )
    Abort( "Input, output, and mask files must be distinct." );
  Input = mri_open_dataset( infile, MRI_READ );
  if (use_mask) Mask  = mri_open_dataset( maskfile, MRI_READ );
  else Mask= NULL;

  /* Check that program will function on datasets */
  checkFormat(Input,infile);
  if (Mask) checkFormat(Mask,maskfile);

  dx= mri_get_int(Input,"images.extent.x");
  dy= mri_get_int(Input,"images.extent.y");
  dz= mri_get_int(Input,"images.extent.z");

  if (Mask) {
    if ( mri_get_int(Mask,"images.extent.x") != dx 
	 || mri_get_int(Mask,"images.extent.y") != dy 
	 || mri_get_int(Mask,"images.extent.z") != dz )
      Abort("%s: input and mask are not commensurate!\n",progname);
  }

  /* Set output dataset */
  Output = mri_copy_dataset( outfile, Input );
  hist_add_cl( Output, argc, argv );
  mri_set_string( Output, "images.datatype", "int32" );

  /* Initialize and clear padded volumes */
  labelVol= allocateAndClearPaddedLong(dx,dy,dz);
  maskVol= allocateAndClearPaddedLong(dx,dy,dz);
  dvol= allocateAndClearPaddedDouble(dx,dy,dz);

  /* Build the neighbor table */
  nbrOffsets= buildNeighborOffsetTable3D(&nNbrs, dx, dy, dz);

  /* Load the data.  Remember that the padding region is beyond the
   * support of the mask, so its values don't matter.
   */
  dbuf= mri_get_chunk(Input, "images", dx*dy*dz, 0, MRI_DOUBLE);
  for (k=0; k<dz; k++)
    for (j=0; j<dy; j++)
      for (i=0; i<dx; i++) {
	PLOC(dvol,i,j,k)= LOC(dbuf,i,j,k);
      }
  mri_close_dataset( Input );
  if (Mask) {
    buf= mri_get_chunk(Mask, "images", dx*dy*dz, 0, MRI_LONG);
    for (k=0; k<dz; k++)
      for (j=0; j<dy; j++)
	for (i=0; i<dx; i++) {
	  PLOC(maskVol,i,j,k)= LOC(buf,i,j,k);
	}
    mri_close_dataset( Mask );
  }
  else {
    /* No mask means a mask of all 1's */
    for (k=0; k<dz; k++)
      for (j=0; j<dy; j++)
	for (i=0; i<dx; i++) {
	  PLOC(maskVol,i,j,k)= 1;
	}
  }

  if (verbose_flag) {
    algorithmToString( &alg, algstring, sizeof(algstring) );
    Message("# Algorithm is <%s>\n",algstring);
  }
  regionCount= applyWatershed(labelVol, dvol, maskVol, nbrOffsets, nNbrs, 
			      dx, dy, dz, &alg);

  /* Write and close output, conveniently counting voxels as we go. */
  if (!(ioVol= (long*)malloc(dx*dy*dz*sizeof(long)))) 
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,dx*dy*dz*sizeof(long));
  for (loop=0; loop<regionCount; loop++)
    regionTable[loop].nvoxels= 0;
  for (k=0; k<dz; k++)
    for (j=0; j<dy; j++)
      for (i=0; i<dx; i++) {
	long val= PLOC(labelVol,i,j,k);
	if (val>0) regionTable[val-1].nvoxels++;
	LOC(ioVol,i,j,k)= val;
      }
  mri_set_chunk( Output, "images", dx*dy*dz, 0, MRI_LONG, ioVol );
  mri_close_dataset( Output );
  free(ioVol);

  /* Write tabular data */
  for (loop=0; loop<regionCount; loop++) {
    Region* rgn= regionTable+loop;
    fprintf(stdout,
	    "region:%ld peakval:%g peakloc:%ld,%ld,%ld nvoxels:%ld\n",
	    loop+1,rgn->max_val,rgn->i,rgn->j,rgn->k,rgn->nvoxels);
  }

  /* Clean up */
  free(nbrOffsets);
  free(dvol);
  free(maskVol);
  free(labelVol);

  if (verbose_flag)
    Message( "# Watershed classification complete (%ld regions).\n",
	     regionCount);
  return 0;
}

