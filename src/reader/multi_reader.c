/************************************************************
 *                                                          *
 *  multi_reader.c                                         *
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
 *  Original programming by Mark Fitzgerald  5-96           *
 *  Modified to exclusively use libmri calls for output,    *
 *       Greg Hood (PSC), 9-98                              *
 *  Modified to read header files, and to use LX2           *
 *       resampling stuff, Joel Welling (PSC/Stats), 5-1999 *
 ************************************************************/

#include <search.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
#include "mri.h"
#if (SGI64 || SGI5 || SGIMP)
#include <bstring.h>
#endif

#include "bio.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"
#include "smartreader.h"

static char rcsid[] = "$Id: multi_reader.c,v 1.37 2006/12/14 22:30:03 welling Exp $";

/* Notes-
 * -We could also miss cases where, for example, 2 subfiles have 
 *  different voxel sizes.
 */

/* Two unit vectors are considered to align if their dot product
 * is greater than or equal to this value.  It's used for checking to 
 * see if multiple slices make a volume.
 */
#define ALIGNED_DOT_LIMIT 0.99999
/* A step between image planes is considered to be "acceptable" if its
 * size is some integer times a fixed minimum step size, plus or minus
 * the value below.  Also, the corner coordinates must not differ by
 * more than this far in the slice plane.
 */
#define STEP_SIZE_MARGIN 0.00005
/* Two slices are considered to be in the same place if they are closer
 * together than this.
 */
#define SAME_SLICE_POSITION_CHANGE 0.0005

/* Handler-specific hook data */
typedef struct multi_data_struct {
  SList* kids;
  SList* kidsInFirstVolume;
  long long offset_shift;
  long long bytes_this_kid;
  long long skips_this_kid;
  int concatMode;
} MultiData;

/* A struct to hold info used in tree-sorting slices by offset */
typedef struct slice_tree_entry_struct {
  double disp; /* Must be first field! */
  int sliceIndex;
  int howManySeen;
} SliceTreeEntry;

static void destroySelf( FileHandler* self )
{
  MultiData* data= (MultiData*)(self->hook);
  FileHandler* kid= NULL;

  while (!slist_empty(data->kids)) {
    FileHandler* kid= (FileHandler*)slist_pop(data->kids);
    FH_DESTROYSELF(kid);
  }

  if (data->kidsInFirstVolume) slist_destroy(data->kidsInFirstVolume,NULL);
  slist_destroy(data->kids, NULL);
  baseDestroySelf(self);
}

static long long calcBytes( FileHandler* kid, KVHash* subInfo )
{
  /* This returns the total number of *useable* bytes in the kid */
  char* fname= kid->fileName;
  char* dimstr= kvGetString(subInfo,"dimstr");
  char* here;
  long long result;
  char buf[64];

  if (dimstr) {
    result= 1;
    here= dimstr;
    while (*here) {
      sprintf(buf,"d%c",*here);
      result *= kvGetInt(subInfo,buf);
      here++;
    }
    result *= srdrTypeSize[ kvGetInt(subInfo,"datatype_in") ];
    if (debug)
      fprintf(stderr,"calcBytes: final result for <%s> is %lld\n",
	      fname,result);
  }
  else {
    /* Apparently this is a raw file for which no dimensions have been set. */
    return kid->totalLengthBytes;
  }

  return result;
}

static long long calcTotalSkip( FileHandler* kid, KVHash* subInfo )
{
  /* This returns the total number of skipped bytes in the kid,
   * not including the kid's start_offset.
   */
  char* fname= kid->fileName;
  char* dimstr= kvGetString(subInfo,"dimstr");
  char* here;
  long long result;
  char buf[64];

  if (dimstr) {
    result= 0;
    here= dimstr;
    while (*here) {
      sprintf(buf,"d%c",*here);
      result *= kvGetInt(subInfo,buf);
      sprintf(buf,"skip.%c",*here);
      if (kvLookup(subInfo,buf))
	result += kvGetLong(subInfo,buf);
      here++;
    }
    if (debug)
      fprintf(stderr,"calcTotalSkip: final result for <%s> is %lld\n",
	      fname,result);
  }
  else {
    /* Apparently this is a raw file for which no dimensions have been set. */
    if (debug) 
      fprintf(stderr,"calcTotalSkip: guessing 0 for a naked raw file\n");
    return 0;
  }

  return result;
}

static int calcSliceOrientation( KVHash* info, double* tlc, 
				 double* x_norm, double* y_norm, 
				 double* z_norm )
{
  double trc[3];
  double brc[3];
  double delta[3];
  double scale;

  if (!getVec3(info,"slice_tlc",tlc)) return 0;
  if (!getVec3(info,"slice_trc",trc)) return 0;
  if (!getVec3(info,"slice_brc",brc)) return 0;

#ifdef never
  fprintf(stderr,"Corners: (%f %f %f) (%f %f %f) (%f %f %f)\n",
	  tlc[0],tlc[1],tlc[2],
	  trc[0],trc[1],trc[2],
	  brc[0],brc[1],brc[2]);
#endif

  /* Orient the slice in the standard Fiasco coordinate system. */
  subtractVec3(delta,trc,tlc);
  scale= normVec3(delta);
  if (scale!=0.0) multVec3( x_norm, delta, 1.0/scale );
  else return 0;

  subtractVec3(delta,trc,brc);
  scale= normVec3(delta);
  if (scale!=0.0) multVec3( y_norm, delta, 1.0/scale );
  else return 0;

  crossVec3( z_norm, x_norm, y_norm );

#ifdef never
  fprintf(stderr,"Norms: (%f %f %f) (%f %f %f) (%f %f %f)\n",
	  x_norm[0], x_norm[1], x_norm[2],
	  y_norm[0], y_norm[1], y_norm[2],
	  z_norm[0], z_norm[1], z_norm[2]);
#endif

  return 1;
}

static int checkSliceOrientation( KVHash* info, double* base_tlc, 
				  double* x_norm_in, double* y_norm_in, 
				  double* z_norm_in, double* disp )
{
  double tlc[3];
  double x_norm[3];
  double y_norm[3];
  double z_norm[3];
  double delta[3];
  double delta_norm[3];

  calcSliceOrientation( info, tlc, x_norm, y_norm, z_norm );
  if (dotVec3(x_norm, x_norm_in) < ALIGNED_DOT_LIMIT) return 0;
  if (dotVec3(y_norm, y_norm_in) < ALIGNED_DOT_LIMIT) return 0;
  if (dotVec3(z_norm, z_norm_in) < ALIGNED_DOT_LIMIT) return 0;

  subtractVec3(delta, tlc, base_tlc);
  *disp= dotVec3(delta, z_norm);
  copyVec3(delta_norm,delta);
  normalizeVec3(delta_norm);
  if (debug) {
    static int checkNum= 1; /* offset by 1 to account for 1st slice */
    fprintf(stderr,"%d: Norm dot products are %f and %f on disp %f\n",
	    checkNum++,dotVec3(delta_norm,x_norm),dotVec3(delta_norm,y_norm),
	    *disp);
    fprintf(stderr,"   delta_norm: %f %f %f\n",delta_norm[0],delta_norm[1],
	    delta_norm[2]);
    fprintf(stderr,"       x_norm: %f %f %f\n",x_norm[0],x_norm[1],x_norm[2]);
    fprintf(stderr,"       y_norm: %f %f %f\n",y_norm[0],y_norm[1],y_norm[2]);
    fprintf(stderr,"       z_norm: %f %f %f\n",z_norm[0],z_norm[1],z_norm[2]);
    fprintf(stderr,"          tlc: %f %f %f\n",tlc[0],tlc[1],tlc[2]);
    fprintf(stderr,"     base_tlc: %f %f %f\n",base_tlc[0],base_tlc[1],
	    base_tlc[2]);
  }
  if (*disp < STEP_SIZE_MARGIN ) {
    /* We've come back around to another copy of the first slice */
    return 1;
  }
  if (dotVec3(delta_norm,x_norm)>STEP_SIZE_MARGIN
      || dotVec3(delta_norm,y_norm)>STEP_SIZE_MARGIN) {
    if (verbose_flg) fprintf(stderr,"Slice stack is skewed.\n");
    return 0;
  }

  return 1;
}

static char pickDim( KVHash* info )
{
  char* dimstr= kvGetString(info,"dimstr");
  char mydim;

  if (!dimstr) mydim= 't';
  else if (!strcmp(dimstr,"xyz")) mydim= 't';
  else if (!strcmp(dimstr,"xyz")) mydim= 't';
  else if (!strcmp(dimstr,"vxyz")) mydim= 't';
  else if (!strcmp(dimstr,"xy")) mydim= 'z';
  else if (!strcmp(dimstr,"vxy")) mydim= 'z';  
  else if (!strcmp(dimstr,"x")) mydim= 'y';
  else if (!strcmp(dimstr,"vx")) mydim= 'y';  
  else {
    for (mydim='a'; mydim<='z'; mydim++) {
      if (!strchr(dimstr,mydim)) break;
    }
    if (mydim>'z')
      Abort("%s: multi_reader: can't find an unused dimension index!\n",
	    progname);
  }
  return mydim;
}

static int useLastDim( FileHandler* self, KVHash* info )
{
  MultiData* data= (MultiData*)(self->hook);
  char* dimstr= kvGetString(info,"dimstr");
  if (data->concatMode) return 1;
  else if (dimstr!=NULL) {
    char* last= dimstr + (strlen(dimstr)-1);
    char buf[64];
    char mydim;
    
    /* If the last extent is 1, we'll use it */
    sprintf(buf,"d%c",*last);
    return( kvGetInt(info,buf)==1 );
  }
  else return 0;
}

static SliceTreeEntry* createSliceTreeEntry( double disp, int sliceIndex )
{
  SliceTreeEntry* result= NULL;
  if (!(result=(SliceTreeEntry*)malloc(sizeof(SliceTreeEntry))))
    Abort("unable to allocate %d bytes!\n",sizeof(SliceTreeEntry));
  result->disp= disp;
  result->sliceIndex= sliceIndex;
  return result;
}

static void destroySliceTreeEntry( SliceTreeEntry* target )
{
  free(target);
}

static int compareSliceTreeEntries( const void* p1, const void* p2 )
{
  SliceTreeEntry* ste1= (SliceTreeEntry*)p1;
  SliceTreeEntry* ste2= (SliceTreeEntry*)p2;
  int result;
  double sep= fabs(ste1->disp - ste2->disp);
  if (sep<=SAME_SLICE_POSITION_CHANGE) result= 0;
  else if (ste1->disp<ste2->disp) result= -1;
  else if (ste1->disp>ste2->disp) result= 1;
  else result= 0;
  return result;
}

static int compareSliceTreeEntryPtrs( const void* p1, const void* p2 )
{
  SliceTreeEntry** steptr1= (SliceTreeEntry**)p1;
  SliceTreeEntry** steptr2= (SliceTreeEntry**)p2;
  return compareSliceTreeEntries( *steptr1, *steptr2 );
}

static int isOrderlyVolume( FileHandler* self, KVHash* info, int* d1, int* d2 )
{
  MultiData* data= (MultiData*)(self->hook);
  KVHash* multi= kvGetHash(info,"multi_hash");
  char* dimstr= kvGetString(info,"dimstr");
  char lastDim= dimstr[strlen(dimstr)-1];
  char buf[64];
  double firstTlc[3];
  double x_norm[3];
  double y_norm[3];
  double z_norm[3];
  int slice_range= 0;
  int slice_reps= 0;
  int mode= 0;
  int count= 0;
  int haveSeenRepeats= 0;
  int highestCountSeen= 0;
  void* treeRoot= NULL;
  int foldIsValid= 1;

  /* slice number information seems to be present in some form.
   * (we'll call it that, because it would be slices
   * if lastDim=='z').  
   *
   * We have a copy of the info for the first "slice" in
   * the top-level info hash, so we can look there to see
   * if there is any point in proceeding.
   */
  sprintf(buf,"%c",lastDim);
  if (testVec3(info,"slice_tlc") && testVec3(info,"slice_trc") 
      && testVec3(info,"slice_brc")) mode= 1;
  else if (kvLookup(info,buf)) {
    mode= 0; /* the old way of doing things */

    /* This tag is at best confusing, so delete the
     * copy we got from the first child. 
     */
    kvDelete(info,buf);
  }
  else return 0; /* we can't distinguish slices */

  /* Scan the slices to identify repeated slices */
  for (slist_totop(data->kids); !slist_atend(data->kids);
       slist_next(data->kids)) {
    FileHandler* kid= (FileHandler*)slist_get(data->kids);
    KVHash* subInfo= kvGetHash(multi,kid->fileName);
    SliceTreeEntry* thisSTE= NULL;
    SliceTreeEntry* matchingSTE= NULL;
    SliceTreeEntry** stePtrPtr= NULL;
    double dist= 0.0;
 
    if (slist_attop(data->kids)) {
      if (mode==0) {
	dist= 0.0;
	firstTlc[0]= firstTlc[1]= 0.0;
	firstTlc[2]= (double)kvGetInt(subInfo,buf);
	x_norm[0]= 1.0; x_norm[0]= 0.0; x_norm[0]= 0.0;
	y_norm[0]= 0.0; y_norm[0]= 1.0; y_norm[0]= 0.0;
	z_norm[0]= 0.0; z_norm[0]= 0.0; z_norm[0]= 1.0;
      }
      else {
	assert(mode==1);
	if (!calcSliceOrientation(subInfo, firstTlc, x_norm, y_norm, z_norm))
	  Abort("%s: unexpectedly failed to find orientation of first slice!\n",
		progname);
      }
    }
    else {
      if (mode==0) {
	dist= (double)kvGetInt(subInfo,buf) - firstTlc[2];
      }
      else {
	if (!(checkSliceOrientation(subInfo, firstTlc, 
				    x_norm, y_norm, z_norm, &dist)))
	  foldIsValid= 0;
      }
    }

    thisSTE= createSliceTreeEntry(dist, count);
    stePtrPtr= (SliceTreeEntry**)tsearch((void*)thisSTE, &(treeRoot),
					 compareSliceTreeEntries);
    if (!stePtrPtr)
      Abort("%s: unable to allocate tsearch tree node!\n",progname);
    matchingSTE= *stePtrPtr;
    if (matchingSTE->sliceIndex == thisSTE->sliceIndex) {
      if (haveSeenRepeats) {
	if (debug)
	  fprintf(stderr,"Slice %d disp %f is a stray!\n",count,thisSTE->disp);
	if (foldIsValid)
	  Warning(1,"%s: multi_reader: irregular slice pattern!\n",
		  progname);
	foldIsValid= 0;
      }
      matchingSTE->howManySeen= 1;
      if (!(data->kidsInFirstVolume))
	data->kidsInFirstVolume= slist_create();
      slist_append(data->kidsInFirstVolume,thisSTE);
      if (debug)
	fprintf(stderr,"Slice %d disp %f is unique\n",
		thisSTE->sliceIndex,thisSTE->disp);
    }
    else {
      destroySliceTreeEntry(thisSTE);
      matchingSTE->howManySeen++;
      if (!haveSeenRepeats) {
	slice_range= slist_count(data->kidsInFirstVolume);
	haveSeenRepeats= 1;
      }
      if ((count%slice_range != matchingSTE->sliceIndex)
	  || (matchingSTE->howManySeen<highestCountSeen)){
	if (foldIsValid)
	  Warning(1,
		  "%s: multi_reader: slice sequence is inconsistent at %s!\n",
		  progname,kid->fileName);
	foldIsValid= 0;
      }
      if (debug)
	fprintf(stderr,"Slice %d disp %f repeats slice %d disp %f (seen %d)\n",
		thisSTE->sliceIndex,thisSTE->disp,
		matchingSTE->sliceIndex,matchingSTE->disp,
		matchingSTE->howManySeen);
      if (matchingSTE->howManySeen>highestCountSeen) 
	highestCountSeen=matchingSTE->howManySeen;
    }
    count++;
  }
  if (!haveSeenRepeats) {
    slice_range= slist_count(data->kidsInFirstVolume);
    haveSeenRepeats= 1;
  }

  /* Destroy the search tree, checking for uniform number of slices
   * in the process.
   */
  for (slist_totop(data->kidsInFirstVolume);
       !slist_atend(data->kidsInFirstVolume);
       slist_next(data->kidsInFirstVolume)) {
    SliceTreeEntry* thisSTE= 
      (SliceTreeEntry*)slist_get(data->kidsInFirstVolume);
    int haveShownSliceCountMessage= 0;
    if (slice_reps==0) slice_reps= thisSTE->howManySeen;
    else {
      if (thisSTE->howManySeen != slice_reps) {
	if (!haveShownSliceCountMessage) {
	  Warning(1,
		  "%s: multi_reader: different counts for different slices!\n",
		  progname);
	  haveShownSliceCountMessage= 1;
	}
	foldIsValid= 0;
      }
    }
    tdelete(thisSTE, &treeRoot, compareSliceTreeEntries);
  }

  /* Sort the slice records; this sort will put them in positional
   * order.
   */
  slist_sort(data->kidsInFirstVolume, compareSliceTreeEntryPtrs);

  *d1= slice_range;
  *d2= slice_reps;
  if (debug)
    fprintf(stderr,
	    "foldIsValid= %d, slice_range= %d; reps= %d; total count %d\n",
	    foldIsValid,slice_range,*d2,slist_count(data->kids));
  return foldIsValid;
}

static int reorderSlicesToFormVolume( FileHandler* self, KVHash* info )
{
  MultiData* data= (MultiData*)(self->hook);
  KVHash* multi= kvGetHash(info,"multi_hash");
  KVHash* defs= kvGetHash(info,"definitions");
  KVHash* extNames= kvGetHash(info,"external_names");
  FileHandler** fh_table= NULL;
  char* dimstr= kvGetString(info,"dimstr");
  char buf[64];
  int dz;
  int dt;
  int z;
  int t;
  int totSlices;
  
  if (!(data->kidsInFirstVolume))
    Abort("%s: internal error: tried to form a volume from data which was never folded!\n",
	  progname);
  
  if (!(fh_table=(FileHandler**)malloc(slist_count(data->kids)
				       *sizeof(FileHandler*)))) 
    Abort("%s: unable to allocate %d bytes!\n",progname,
	  slist_count(data->kids)*sizeof(FileHandler*));
  
  /* This is a volume, so strlen(dimstr) is at least 3 */
  sprintf(buf,"d%c",dimstr[strlen(dimstr)-1]);
  dt= kvGetInt(info,buf);
  if (dt==slist_count(data->kids)) {
    /* There is only one instance of the volume */
    dz= dt;
    dt= 1;
  }
  else {
    sprintf(buf,"d%c",dimstr[strlen(dimstr)-2]); 
    dz= kvGetInt(info,buf);
  }
  totSlices= dz*dt;
  if (totSlices != slist_count(data->kids)) {
    if (verbose_flg) 
      fprintf(stderr,
	      "Cannot verify slice order; files are not one slice each!\n");
    return 0;
  }
  
  if (!(fh_table=(FileHandler**)malloc(totSlices
				       *sizeof(FileHandler*)))) 
    Abort("%s: unable to allocate %d bytes!\n",progname,
	  slist_count(data->kids)*sizeof(FileHandler*));
  
  slist_totop(data->kids);
  if (debug) {
    fprintf(stderr,"dz= %d, dt= %d, Length of kidlist: %d\n",dz,dt,slist_count(data->kids));
    fprintf(stderr,"Kids in first volume: %d\n",slist_count(data->kidsInFirstVolume));
  }
  for (t=0; t<dt; t++)
    for (z=0; z<dz; z++)
      fh_table[t*dz+z]= (FileHandler*)slist_pop(data->kids);
  for (t=0; t<dt; t++) {
    for (slist_totop(data->kidsInFirstVolume);
	 !slist_atend(data->kidsInFirstVolume);
	 slist_next(data->kidsInFirstVolume)) {
      SliceTreeEntry* thisSTE= 
	(SliceTreeEntry*)slist_get(data->kidsInFirstVolume);
      z= thisSTE->sliceIndex;
      slist_append(data->kids, fh_table[t*dz+z]);
    }
  }
  free(fh_table);
  
  return 1;

}

static int locateVolumeCorners( FileHandler* self, KVHash* info )
{
  MultiData* data= (MultiData*)(self->hook);
  KVHash* multi= kvGetHash(info,"multi_hash");
  KVHash* defs= kvGetHash(info,"definitions");
  KVHash* extNames= kvGetHash(info,"external_names");
  double tlc[3];
  double x_norm[3];
  double y_norm[3];
  double z_norm[3];
  double vol_tlf[3];
  double vol_tlb[3];
  double vol_trf[3];
  double vol_trb[3];
  double vol_blf[3];
  double vol_blb[3];
  double vol_brf[3];
  double vol_brb[3];
  double voxel_x;
  double voxel_y;
  double slice_sep;
  FileHandler* kid1= NULL;
  KVHash* subInfo1= NULL;
  FileHandler* kid2= NULL;
  KVHash* subInfo2= NULL;

#ifdef never
  /* Start by deleting image corners we may have inherited from
   * a child slice.
   */
  kvDeleteAll(info,"slice_tlc.0");
  kvDeleteAll(info,"slice_tlc.1");
  kvDeleteAll(info,"slice_tlc.2");
  kvDeleteAll(info,"slice_trc.0");
  kvDeleteAll(info,"slice_trc.1");
  kvDeleteAll(info,"slice_trc.2");
  kvDeleteAll(info,"slice_brc.0");
  kvDeleteAll(info,"slice_brc.1");
  kvDeleteAll(info,"slice_brc.2");
#endif

  slist_totop(data->kids);
  kid1= (FileHandler*)slist_get(data->kids);
  subInfo1= kvGetHash(multi,kid1->fileName);
  slist_next(data->kids);
  kid2= (FileHandler*)slist_get(data->kids);
  subInfo2= kvGetHash(multi,kid2->fileName);
  if (!kvLookup(subInfo1,"slice_tlc.0") || !kvLookup(subInfo2,"slice_tlc.0")) {
    /* Give up if these slices doesn't know their location */
    return 0;
  }
  voxel_x= kvGetDouble(subInfo1,"voxel_x");
  voxel_y= kvGetDouble(subInfo1,"voxel_y");
  if (!calcSliceOrientation( subInfo1, tlc, x_norm, y_norm, z_norm ))
    return 0;
  if (!checkSliceOrientation( subInfo2, tlc, x_norm, y_norm, z_norm, 
			      &slice_sep ))
    return 0;

  if (slice_sep<=0.0)
    Abort("%s: unexpectedly found inverted or coplanar slices!\n",progname);

  kvDefDouble(info,"slice_thickness", kvGetDouble(info,"voxel_z"));
  kvDefDouble(info,"slice_gap", 
	      slice_sep - kvGetDouble(info,"voxel_z"));
  kvDefDouble(info,"voxel_z", slice_sep);
  copyVec3(vol_blb, tlc);
  xplusbyVec3(vol_blf, vol_blb, y_norm, 
	      -kvGetInt(info,"dy")*kvGetDouble(info,"voxel_y"));
  xplusbyVec3(vol_brb, vol_blb, x_norm, 
	      kvGetInt(info,"dx")*kvGetDouble(info,"voxel_x"));
  xplusbyVec3(vol_brf, vol_blf, x_norm, 
	      kvGetInt(info,"dx")*kvGetDouble(info,"voxel_x"));
  xplusbyVec3(vol_tlf, vol_blf, z_norm, 
	      kvGetInt(info,"dz")*kvGetDouble(info,"voxel_z"));
  xplusbyVec3(vol_tlb, vol_blb, z_norm, 
	      kvGetInt(info,"dz")*kvGetDouble(info,"voxel_z"));
  xplusbyVec3(vol_trf, vol_brf, z_norm, 
	      kvGetInt(info,"dz")*kvGetDouble(info,"voxel_z"));
  xplusbyVec3(vol_trb, vol_brb, z_norm, 
	      kvGetInt(info,"dz")*kvGetDouble(info,"voxel_z"));
  defVec3(info, "tlf", vol_tlf);
  defVec3(info, "trf", vol_trf);
  defVec3(info, "tlb", vol_tlb);
  defVec3(info, "trb", vol_trb);
  defVec3(info, "blf", vol_blf);
  defVec3(info, "brf", vol_brf);
  defVec3(info, "blb", vol_blb);
  defVec3(info, "brb", vol_brb);

  return 1;
}

static int canConcatenate( FileHandler* f1, KVHash* info1,
			   FileHandler* f2, KVHash* info2 )
{
  const char* d1= kvGetString(info1,"dimstr");
  const char* d2= kvGetString(info2,"dimstr");
  if (strcmp(d1,d2))
    return 0;
  else {
    int i;
    for (i=0; i<strlen(d1)-1; i++) {
      char buf[8];
      snprintf(buf,sizeof(buf),"d%c",d1[i]);
      if (kvGetInt(info1,buf)!=kvGetInt(info2,buf))
	return 0;
    }
    return 1;
  }
}

static void scanForConsistency( FileHandler* self, KVHash* info )
{
  /* Scan for consistency over files.  If the files are not 
   * sufficiently similar, simply abort.
   */
  MultiData* data= (MultiData*)(self->hook);
  KVHash* multi= kvGetHash(info,"multi_hash");
  KVHash* defs= kvGetHash(info,"definitions");
  KVHash* extNames= kvGetHash(info,"external_names");
  KVHash* subInfo;
  KVHash* subInfo0;
  FileHandler* kid;
  FileHandler* kid0;
  int datatype= 0;
  int handler_datatype_out= 0;
  char buf[8];
  char* dimstr= NULL;
  int lastExtent= 0;
  
  /* Scan for consistency */
  slist_totop(data->kids);
  while (!(slist_atend(data->kids))) {
    kid= kid0= (FileHandler*)slist_get(data->kids);
    subInfo= subInfo0= kvGetHash(multi,kid->fileName);

    if (debug) fprintf(stderr,"Scanning <%s> for consistency\n",
		       kid->fileName);

    if (slist_attop(data->kids)) {
      /* First file - export potentially useful information. */
      KVHash* subDefs= kvGetHash(subInfo,"definitions");
      KVHash* subExtNames= kvGetHash(subInfo,"external_names");
      kvCopyUniqueExceptHashes(info, subInfo);
      kvCopyUniqueExceptHashes(defs, subDefs);
      kvCopyUniqueExceptHashes(extNames, subExtNames);

      /* We need to make sure other files are at least minimally consistent */
      datatype= kvGetInt(subInfo,"datatype_in");
      handler_datatype_out= kvGetInt(subInfo,"handler_datatype_out");
      dimstr= kvGetString(subInfo,"dimstr");
      snprintf(buf,sizeof(buf),"d%c",dimstr[strlen(dimstr)-1]);
      lastExtent= kvGetInt(subInfo,buf);
    }
    else {
      if (kvGetInt(subInfo,"datatype_in") != datatype)
	Abort("%s: multi_reader: files do not have the same datatype!\n",
	      progname);
      if (kvGetInt(subInfo,"handler_datatype_out") != handler_datatype_out)
	Abort("%s: multi_reader: files do not have the same datatype!\n",
	      progname);
      if (strcmp(dimstr,kvGetString(subInfo,"dimstr")))
	Abort("%s: multi_reader: files do not have the same dim string!\n",
	      progname);
      /* Buf still contains 'dx', where x is the last dimension of kid0 */
      if (kvGetInt(subInfo,buf)==lastExtent) {
	/* Leave data->concatMode off, which will result in the data
	 * being blocked rather than concatenated within the slowest dimension.
	 */
      }
      else {
	if (canConcatenate(kid0,subInfo0,kid,subInfo)) {
	  data->concatMode= 1;
	}
	else {
	  Abort("%s: multi_reader: file structure does not allow concatenation!\n",
		progname);
	}
      }
      
    }

    slist_next(data->kids);
  }
}

static void setReorderPattern( FileHandler* self, KVHash* info )
{
  MultiData* data= (MultiData*)(self->hook);
  KVHash* multi= kvGetHash(info,"multi_hash");
  KVHash* defs= kvGetHash(info,"definitions");
  KVHash* extNames= kvGetHash(info,"external_names");
  int dz= slist_count(data->kidsInFirstVolume);
  char* buf= NULL;
  int* indexTable= NULL;
  char wordbuf[64];
  const char* name;
  int i;

  if (!(data->kidsInFirstVolume))
    Abort("%s: internal error: tried to pick reorder pattern data which was never folded!\n",
	  progname);

  if (!(indexTable=(int*)malloc(dz*sizeof(int))))
    Abort("%s: unable to allocate %d bytes!\n",progname,dz*sizeof(int));

  i= 0;
  for (slist_totop(data->kidsInFirstVolume);
       !slist_atend(data->kidsInFirstVolume);
       slist_next(data->kidsInFirstVolume)) {
    SliceTreeEntry* thisSTE= 
      (SliceTreeEntry*)slist_get(data->kidsInFirstVolume);
    indexTable[thisSTE->sliceIndex]= i;
    i++;
  }

  if (name= slp_findSlicePatternNameFromTable(dz, indexTable)) {
    kvDefString(info,"reorder_pattern",name);
  }
  else {
    /* As a last resort, enumerate the reorder pattern with a string. */
    int bufsize= 4*dz;
    if (!(buf=(char*)malloc(bufsize*sizeof(char))))
      Abort("%s: unable to allocate %d bytes!\n",progname,
	    bufsize*sizeof(char));

    buf[0]= '\0';
    for (i=0; i<dz; i++) {
      sprintf(wordbuf,"%d",indexTable[i]);
      if (i!=0) strcat(buf,",");
      strncat(buf,wordbuf,bufsize-(strlen(buf)+2));
      if (strlen(buf)+strlen(wordbuf)>=bufsize)
	Abort("%s: ran out of space writing reorder_pattern string!\n",
	      progname);
    }
    kvDefString(info,"reorder_pattern",buf);
    free(buf);
  }

  free(indexTable);
}

static void deriveCollectiveInfo( FileHandler* self, KVHash* info )
{
  MultiData* data= (MultiData*)(self->hook);
  KVHash* multi= kvGetHash(info,"multi_hash");
  KVHash* defs= kvGetHash(info,"definitions");
  KVHash* extNames= kvGetHash(info,"external_names");
  char* dimstr;
  char buf[64];
  char buf2[64];
  char* here;
  char lastDim= 0;
  char appendDim= 0;
  char skipDim= 0;
  int lastExtent= 0;
  int appendExtent= 0;
  int needsVolumeConstruction= 0;
  
  slist_totop(data->kids);

  kvDefLong(info,"start_offset",0); /* we will hide kid's start offsets */

  if (kvLookup(info,"dimstr")) {
    if ( strlen(kvGetString(info,"dimstr")) > sizeof(buf)-3 )
      Abort("%s: multi_reader: dimension string is too long!\n", progname);
    strcpy(buf, (char*)kvGetString(info,"dimstr"));
    if (strlen(buf)>0) lastDim= *(buf + strlen(buf) - 1);
  }
  else buf[0]= '\0';
  
  if (isOrderlyVolume(self, info, &lastExtent, &appendExtent)) {

    if (lastDim=='x' || lastDim=='y' || lastDim=='z')
      needsVolumeConstruction= 1;
    else
      needsVolumeConstruction= 0;

    if (needsVolumeConstruction) {
      if (useLastDim(self,info)) {
	if (appendExtent>1) {
	  appendDim= pickDim( info );
	}
	else {
	  appendDim= 0;
	  appendExtent= 0;
	}
      }
      else {
	if (appendExtent>1) {
	  /* We need to append *2* dimensions */
	  lastDim= pickDim( info );
	  here= buf + strlen(buf);
	  *here++= lastDim;
	  *here++= '\0';
	  kvDefString(info,"dimstr",buf);
	  appendDim= pickDim( info );
	}
	else {
	  appendDim= pickDim( info );
	  appendExtent= lastExtent;
	  
	}
      }
    }
    else {
      /* This is an odd case- somehow an orderly volume, but not
       * made up of slices.  Let's just ignore the fact that it
       * might benefit from volume construction.
       */
      needsVolumeConstruction= 0;
      lastExtent= slist_count(data->kids);
      appendExtent= 0;
      if (!useLastDim(self,info)) {
	/* can't use last dim, so we need a new dim */
	/* append one, and treat it like it was there all along. */
	lastDim= pickDim( info );
	here= buf + strlen(buf);
	*here++= lastDim;
	*here++= '\0';
	kvDefString(info,"dimstr",buf);
      }
    }    
  }
  else {
    needsVolumeConstruction= 0;
    if (useLastDim(self,info)) {
      char buf3[8];
      snprintf(buf3,sizeof(buf3),"d%c",lastDim);
      if (data->concatMode) {
	lastExtent= 0;
	for (slist_totop(data->kids); !slist_atend(data->kids);
	     slist_next(data->kids)) {
	  FileHandler* kid= (FileHandler*)slist_get(data->kids);
	  KVHash* subInfo= kvGetHash(multi,kid->fileName);
	  lastExtent += kvGetInt(subInfo,buf3);
	}
      }
      else {
	lastExtent= slist_count(data->kids);
      }
      appendDim= 0;
      appendExtent= 0;
    }
    else {
      appendDim= pickDim(info);
      appendExtent= slist_count(data->kids);
    }
  }

  /* Make sure reading is broken at the slice boundaries */
  if (strlen(buf)>1) skipDim= *(strchr(buf,lastDim)-1);
  else skipDim= 0;

  if (appendExtent != 0) {
    here= buf + strlen(buf);
    *here++= appendDim;
    *here++= '\0';
    kvDefString(info,"dimstr",buf);
    sprintf(buf2,"d%c",appendDim);
    kvDefInt(info,buf2,appendExtent);
  }

  if (lastExtent != 0) {
    kvDefString(info,"dimstr",buf);
    sprintf(buf2,"d%c",lastDim);
    kvDefInt(info,buf2,lastExtent);
  }

  if (skipDim != 0) {
    /* Make sure this dim has an associated skip */
    sprintf(buf,"skip.%c",skipDim);
    if (!kvLookup(info,buf)) 
      kvDefLong(info, buf, 0);
  }
  
  dimstr= kvGetString(info,"dimstr");
  if ((strchr(dimstr,'x') != NULL) && (strchr(dimstr,'y') != NULL) 
      && (strchr(dimstr,'z') != NULL) && needsVolumeConstruction) {
    /* Well, it appears to be a volume of some sort. */
    if (reorderSlicesToFormVolume(self, info)
	&& locateVolumeCorners(self, info)) {
      sprintf(buf2,"description.%c",lastDim);
      kvDefString(info,buf2,"gridded image-space");
      kvDefBoolean(info,"reorder",0);
      setReorderPattern(self, info);
    }
  }
  slist_totop(data->kids);
}

static int fhCompare( const void* p1, const void* p2 )
{
  FileHandler** f1= (FileHandler**)p1;
  FileHandler** f2= (FileHandler**)p2;
  return FH_COMPARE((*f1),(*f2));
}

static void processHeader( FileHandler* self, KVHash* info, SList* chunkStack )
{
  MultiData* data= (MultiData*)(self->hook);
  KVHash* multi= kvFactory(KV_DEFAULT_SIZE);
  FileHandler* kid;

  slist_totop(data->kids);
  while ((kid= (FileHandler*)slist_next(data->kids)) != NULL) {
    KVHash* subInfo= kvCloneUnique(info);
    int oldDebug= debug;
    if (debug) {
      fprintf(stderr,"Multi: Processing header for <%s>\n",
	      kid->fileName);
      debug= 0; /* silence reading of kids; can debug that by reading
		 * them singly.
		 */
    }
    FH_PROCESSHEADER( kid, subInfo, chunkStack );
    FH_CLOSE( kid ); /* to avoid having too many open files */
    debug= oldDebug;
    kvDefHash(multi, kid->fileName, subInfo);
  }
  slist_totop(data->kids);

  /* We'll attach our own little world to the the outer universe */
  kvDefHash( info, "multi_hash", multi );

  /* Are all the kids sufficiently alike? If not, abort. */
  scanForConsistency( self, info ); /* doesn't return if tests fail */

  /* Sort the kids into the "appropriate" order */
  slist_sort(data->kids, fhCompare);

  /* Get set up to handle the first read, for which we'll use
   * the first available file.
   */
  slist_totop(data->kids);
  data->offset_shift= 0;
  kid= (FileHandler*)slist_get(data->kids);
  data->bytes_this_kid= calcBytes(kid, 
				  kvGetHash(multi,kid->fileName));
  data->skips_this_kid= calcTotalSkip(kid, 
				      kvGetHash(multi,kid->fileName));

  /* Now we need to make up a description of the whole collection */
  deriveCollectiveInfo( self, info );
}

static void multiRead( FileHandler* self, KVHash* info,
		       long long offset, long n,
		       SRDR_Datatype datatype,
		       void* obuf )
{
  MultiData* data= (MultiData*)(self->hook);
  FileHandler* kid= (FileHandler*)slist_get(data->kids);
  KVHash* multi= kvGetHash(info,"multi_hash");
  KVHash* kidInfo= kvGetHash(multi,kid->fileName);

  FH_REOPEN(kid);

  if (debug) 
    fprintf(stderr,
	    "Testing for read %d bytes from %s, rel offset %lld, end %lld\n",
	    n*srdrTypeSize[datatype], kid->fileName, offset - data->offset_shift, 
	    data->bytes_this_kid+data->skips_this_kid);

  while (offset >= (data->offset_shift
		    +data->bytes_this_kid+data->skips_this_kid)) {

    /* Skip to next file */
    if (debug) fprintf(stderr,"Moving past end of subfile %s\n",kid->fileName);
    data->offset_shift += data->bytes_this_kid+data->skips_this_kid;
    FH_CLOSE(kid);

    slist_next(data->kids);
    if (!(kid= (FileHandler*)slist_get(data->kids)))
      Abort("%s: multiRead: read at %lld is past end of last file!\n",
	    progname,offset);
    if (debug) fprintf(stderr,"new subfile is <%s>\n",kid->fileName);
    FH_REOPEN(kid);
    kidInfo= kvGetHash(multi,kid->fileName);
    data->bytes_this_kid= calcBytes(kid, kidInfo);
    data->skips_this_kid= calcTotalSkip(kid, kidInfo);
  }

  if ((offset + n*srdrTypeSize[datatype]) 
      > (data->offset_shift + data->bytes_this_kid + data->skips_this_kid)) {
    /* This read will cross a file end, which isn't allowed. */
    Abort("%s: multiRead: framing error reading %d %s's at %lld!\n",
	  progname,n,srdrTypeName[datatype],offset);
  }

  if (debug) 
    fprintf(stderr,
	    "Trying to read %d bytes from %s, rel offset %lld\n",
	    n*srdrTypeSize[datatype], kid->fileName, offset - data->offset_shift);

  FH_READ( kid, kidInfo, 
	   (offset - data->offset_shift + kvGetLong(kidInfo,"start_offset")), 
	   n, datatype, obuf );
}

FileHandler* multiFileHandlerFactory()
{
  FileHandler* result= baseFactory("NotARealFile");
  MultiData* data;
  
  result->typeName= strdup("Multi(***uninitialized***)");
  
  result->destroySelf= destroySelf;
  result->read= multiRead;
  result->processHeader= processHeader;
  
  if (!(data= (MultiData*)malloc(sizeof(MultiData))))
    Abort("%s: unable to allocate %d bytes!\n",sizeof(MultiData));
  
  data->kids= slist_create();
  data->kidsInFirstVolume= NULL;
  data->offset_shift= 0;
  data->bytes_this_kid= 0;
  data->concatMode= 0;

  result->hook= data;

  return result;
}

void multiFileHandlerAddFile(FileHandler* self, FileHandler* kid)
{
  MultiData* data= (MultiData*)(self->hook);
  
  if (slist_empty(data->kids)) {
    /* How nice! Our first child! */
    char* typeName;
    int typeNameLength;

    typeNameLength= strlen("Multi(...)") + strlen(kid->typeName)+8;
    if (!(typeName=(char*)malloc(typeNameLength)))
      Abort("%s: unable to allocate %d bytes!\n",progname,typeNameLength);
    sprintf(typeName,"Multi(%s...)",kid->typeName);
    free(self->typeName);
    self->typeName= typeName;
  }

  slist_append(data->kids, kid);
  self->totalLengthBytes += kid->totalLengthBytes;
}


