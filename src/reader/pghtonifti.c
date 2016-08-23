/************************************************************
 *                                                          *
 *  pghtonifti.c                                             *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 2007 Department of Statistics,         *
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
 *  Derived from smartreader, Joel Welling 6/07
 ************************************************************/

/* Notes-
 */

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include "mri.h"
#if (SGI64 || SGI5 || SGIMP)
#include <bstring.h>
#endif
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <assert.h>

#include "bio.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"
#include "smartreader.h"
#include "nifti1.h"

static char rcsid[] = "$Id: pghtonifti.c,v 1.9 2007/07/07 18:46:06 welling Exp $";

/* Notes-
 */

/* We need to place a limit on the maximum number of "things" in one IO op */
#define MAX_BLOCK (16*1024*1024)

int debug = 0;          /* Global debug flag                         */

int verbose_flg = 0;        /* Global verbosity value (0=off) */

char* progname= NULL; /* program name */

static void initStuff(KVHash* info)
{
  initInfoHash(info);

  /* Set up some default values specific to the first chunk */
  kvDefBoolean(info,"big_endian_input",1); /* fmri data usually bigendian */
  kvDefInt(info,"datatype_in",SRDR_INT16);
  kvDefDouble(info,"autoscale_range",1.0);
  kvDefInt(info,"start_offset",0); 
  kvDefString(info,"chunkname","images");
  kvDefString(info,"chunkfile",".dat");
}

static int reconcile(KVHash* info)
{
  /* For automatically-generated chunks there may be no predefined
   * datatype_out.  In this case, do what the handler likes best.
   */
  if (!kvLookup(info,"datatype_out")) {
    kvDefInt(info, "datatype_out", kvGetInt(info,"handler_datatype_out"));
    if (verbose_flg) 
      fprintf(stderr,"reconciling chunk %s to type %s\n",
	      kvGetString(info,"chunkname"), srdrTypeName[kvGetInt(info,"datatype_out")]);
  }

  return 1;
}

static void parse_command_line( KVHash* info, 
				char* inputFileName, char* outFileName,
				int argc, char* argv[] )
{
  char string[512];
  int itmp;
  double dtmp;
  KVHash* defs= kvGetHash(info,"definitions");

  cl_scan( argc, argv );

  /* Get filenames */
  if (!cl_get( "", "%s", inputFileName )) {
    fprintf(stderr,"%s: Input file name not given!\n",progname);
    Help("usage");
    exit(-1);
  }
  if (!cl_get( "", "%s", outFileName )) {
    fprintf(stderr,"%s: Output file name not given!\n",progname);
    Help("usage");
    exit(-1);
  }

  debug= cl_present( "debug" );
  verbose_flg= cl_present( "verbose|v" );
  if (cl_present("nii")) kvDefBoolean(info,"cl_nii_mode",1);

  kvDefString(info,"cl_input_file_name",inputFileName);

  if (cl_cleanup_check()) {
    int i;
    fprintf(stderr,"%s: invalid argument in command line:\n    ",argv[0]);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(-1);
  }

  kvDefInt(info,"handler_datatype_out", kvGetInt(info,"datatype_in"));
}

static int canTranslate(SList* chunkList)
{
  int found= 0;
  ChunkHandlerPair* imgPair= NULL;
  KVHash* info= NULL;
  const char* dimstr;

  slist_totop(chunkList);
  while (!slist_atend(chunkList)) {
    imgPair= (ChunkHandlerPair*)slist_next(chunkList);
    if (kvLookup(imgPair->info,"chunkname") 
	&& !strcmp(kvGetString(imgPair->info,"chunkname"),"images")) {
      found= 1;
      break;
    }
  }
  slist_totop(chunkList);

  if (!found) {
    Error("%s: input file has no 'images' chunk!\n",progname);
    return 0;
  }

  info= imgPair->info;
  dimstr= kvGetString(info,"dimstr");
  if (strncmp(dimstr,"vxyztabc",strlen(dimstr))
      && strncmp(dimstr,"xyztabc",strlen(dimstr))) {
    Error("%s: input file image dimensions <%s> do not translate!\n",
	  progname, dimstr);
    return 0;
  }

  if (*dimstr == 'v') {
    int dv= kvGetInt(info,"dv");
    if (dv<1 || dv>3) {
      Error("%s: input file image v dimension of %d does not translate!\n",
	    progname,dv);
      return 0;
    }
  }

  return 1;
}

static void openOutputFiles(const char* fname, KVHash* info,
			    FILE** head, FILE** brik)
{
  char* buf;
  char* here;

  if (!(buf=(char*)malloc(strlen(fname)+6)))
    Abort("%s: unable to allocate %d bytes!\n",progname,strlen(fname)+6);
  strcpy(buf,fname);

  here= strrchr(buf,'.');
  if (here && (!strcmp(here,".hdr") || !strcmp(here,".img")
	       || !strcmp(here,".nii"))) {
    *here= '\0';
  }
  else here= buf + strlen(buf);

  if (kvLookup(info,"cl_nii_mode") && kvGetBoolean(info,"cl_nii_mode")) {
    strcat(buf,".nii");
    if (!(*head=fopen(buf,"w")))
      Abort("%s: cannot open <%s> for writing!\n",progname,buf);
    *brik= *head; /* Write data and header to same file */
  }
  else {
    strcat(buf,".hdr");
    if (!(*head=fopen(buf,"w")))
      Abort("%s: cannot open <%s> for writing!\n",progname,buf);
    *here= '\0';
    strcat(buf,".img");
    if (!(*brik=fopen(buf,"w")))
      Abort("%s: cannot open <%s> for writing!\n",progname,buf);
  }
  free(buf);
}

static void closeOutputFiles(KVHash* info, FILE* head, FILE* brik)
{
  if (kvLookup(info,"cl_nii_mode") && kvGetBoolean(info,"cl_nii_mode")) {
    (void)fclose(head);
  }
  else {
    (void)fclose(head);
    (void)fclose(brik);
  }
}


static void kvDefGeneric( KVHash* kvh, const char* key, KVPair* p )
{
  switch (kvType(p)) {
  case KV_STRING: kvDefString(kvh, key, p->v.s);
    break;
  case KV_LONG: kvDefLong(kvh, key, p->v.l);
    break;
  case KV_DOUBLE: kvDefDouble(kvh, key, p->v.d);
    break;
  case KV_HASH: kvDefHash(kvh, key, kvCloneUnique(p->v.h));
    break;
  case KV_BOOLEAN: kvDefBoolean(kvh, key, (p->v.l != 0));
    break;
  case KV_INT: kvDefInt(kvh, key, (int)p->v.l);
    break;
  }
}		  

static KVHash* mergeTags( SList* chunkList )
{
  KVHash* result= kvFactory(KV_DEFAULT_SIZE);
  SList* otherChunks= slist_create();

  while (!slist_empty(chunkList)) {
    ChunkHandlerPair* chPair= (ChunkHandlerPair*)slist_pop(chunkList);
    KVHash* info= chPair->info;
    KVHash* defs= kvGetHash(info,"definitions");
    KVHash* extNames= kvGetHash(info,"external_names");

    if (!strcmp(kvGetString(info,"chunkname"), "images")) {
      KVIterator* kvi= kvUniqueIteratorFactory(info);
      kvDefHash(result,"definitions",
		kvCloneUnique(kvGetHash(info,"definitions")));
      kvDefHash(result,"external_names",
		kvCloneUnique(kvGetHash(info,"external_names")));
      while (kvIteratorHasMorePairs(kvi)) {
	KVPair* p= kvIteratorNextPair(kvi);

	kvDefGeneric(result, kvKey(p), p);
      }
      /* chPair->handler is still referenced, for use transferring
       * data in the main chunk.
       */
      kvDestroy(chPair->info);
      free(chPair);
    }
    else {
      /* We will snag history, "nifti_", and misc. chunk
       * information at this point.
       */
      const char* chunkname= kvGetString(info,"chunkname");
      KVIterator* kvi= kvUniqueIteratorFactory(info);
      while (kvIteratorHasMorePairs(kvi)) {
	KVPair* p= kvIteratorNextPair(kvi);
	if (!strncmp(kvKey(p), "history.",strlen("history."))
	    || !strncmp(kvKey(p), "nifti_",strlen("nifti_"))
	    || !strncmp(kvKey(p), "cl_", strlen("cl_")))
	  kvDefGeneric(result, kvKey(p), p);
      }
      slist_append(otherChunks, chPair);
    }
  }  

  /* Slosh the chunks other than "images" back into the now-empty
   * chunk list.
   */
  while (!slist_empty(otherChunks)) 
    slist_append(chunkList, slist_pop(otherChunks));
  slist_destroy(otherChunks,NULL);

  slist_totop(chunkList);

  return result;
}

static void transfer( FileHandler* handler, KVHash* info, FILE* niftiImg,
		      nifti_1_header* hdr )
{
  long long totalCount;
  long long offset= 0;
  long countThisBlock;
  void* buf= NULL;
  int datatype;
  const char* dimstr= kvGetString(info,"dimstr");
  const char* here= NULL;

  here= dimstr;
  totalCount= 1;
  while (*here) {
    char buf[8];
    snprintf(buf,sizeof(buf),"skip.%c",*here);
    if (kvLookup(info,buf)) 
      Abort("%s: internal error: the input filehandler requests skip.%c, which is not supported!\n",
	    progname,*here);
    snprintf(buf,sizeof(buf),"d%c",*here);
    totalCount *= kvGetInt(info,buf);
    here++;
  }

  datatype= kvGetInt(info,"handler_datatype_out");

  if (debug) fprintf(stderr,"Transferring %lld of type %s\n",
		     totalCount,srdrTypeName[datatype]);

  if (!(buf= (void*)malloc(MAX_BLOCK*srdrTypeSize[datatype])))
    Abort("%s: unable to allocate %d bytes!\n",progname,
	  MAX_BLOCK*srdrTypeSize[datatype]);

  bigfile_fseek( niftiImg, (long)(hdr->vox_offset), SEEK_SET );
  while (totalCount>0) {
    if (MAX_BLOCK<totalCount) countThisBlock= MAX_BLOCK;
    else countThisBlock= totalCount;

    FH_READ( handler, info, offset*srdrTypeSize[datatype], countThisBlock,
	     datatype, buf );
    if (fwrite(buf, srdrTypeSize[datatype], (size_t) countThisBlock, niftiImg)
	!= (size_t)countThisBlock) {
      perror("Error writing IMG");
      Abort("%s: error writing output img data!\n",progname);
    }

    totalCount -= countThisBlock;
    offset += countThisBlock;
  }

  free(buf);
}


static void addCLToHistory( KVHash* info, int argc, char* argv[] )
{
  int i;
  int totLength;
  char buf[64];
  char* clbuf= NULL;
  char* here;
  int thisEntry;

  /* Find the final existing history entry */
  i= 1;
  while (1) {
    snprintf(buf, sizeof(buf), "history.%d", i);
    if (!kvLookup(info,buf)) {
      snprintf(buf, sizeof(buf), "history.%03d", i);
      if (!kvLookup(info,buf)) break;
    }
    i++;
  } 
  thisEntry= i;

  totLength= 1;
  for (i=0; i<argc; i++) totLength += strlen(argv[i])+1;

  if (!(clbuf=(char*)malloc(totLength*sizeof(char))))
    Abort("%s: unable to allocate %d bytes!\n",progname,
	  totLength*sizeof(char));

  here= clbuf;
  sprintf(here,"%s",argv[0]);
  here+= strlen(argv[0]);
  for (i=1; i<argc; i++) {
    sprintf(here," %s",argv[i]);
    here += strlen(argv[i])+1;
  }

  kvDefString(info, buf, clbuf);

  free(clbuf);
}

static void writeNiftiHeader( FILE* f, nifti_1_header* hdr )
{
  unsigned char extensionCodes[4]= {0,0,0,0};

  /* This routine actually exports the header */
  /* libbio will abort on errors */
  FWrInt32(f,hdr->sizeof_hdr);
  FWrUInt8Array(f, &(hdr->data_type[0]), 10);
  FWrUInt8Array(f, &(hdr->db_name[0]), 18);
  FWrInt32(f,hdr->extents);
  FWrInt16(f,hdr->session_error);
  FWrUInt8(f,hdr->regular);
  FWrUInt8(f,hdr->dim_info);
  FWrInt16Array(f, &(hdr->dim[0]), 8);
  FWrFloat32Array(f, &(hdr->intent_p1), 3);
  FWrInt16Array(f, &(hdr->intent_code), 4);
  FWrFloat32Array(f, &(hdr->pixdim[0]), 11);
  FWrInt16(f,hdr->slice_end);
  FWrUInt8(f,hdr->slice_code);
  FWrUInt8(f,hdr->xyzt_units);
  FWrFloat32Array(f, &(hdr->cal_max), 4);
  FWrInt32Array(f, &(hdr->glmax), 2);
  FWrUInt8Array(f, &(hdr->descrip[0]), 80);
  FWrUInt8Array(f, &(hdr->aux_file[0]), 24);
  FWrInt16(f,hdr->qform_code);
  FWrInt16(f,hdr->sform_code);
  FWrFloat32Array(f, &(hdr->quatern_b),18);
  FWrUInt8Array(f, &(hdr->intent_name[0]), 16);
  FWrUInt8Array(f, &(hdr->magic[0]), 4);

  /* Write an empty 'extensions' code */
  FWrUInt8Array(f, extensionCodes, 4);
}

static void setLegacyOrientation( KVHash* info, nifti_1_header* hdr)
{
  int i;

  /* The best we can do here is to indicate the fact that *all* the 
   * Pgh MRI index directions run negative to the NIfTI coordinate 
   * directions.
   */
  hdr->qform_code= 0;
  for (i=0; i<3; i++)
    if (hdr->pixdim[i]!=0.0) hdr->pixdim[i] *= -1;
}

static int trySetQuaternionOrientation( KVHash* info, nifti_1_header* hdr)
{
  Transform T;
  Transform V;
  Transform L;
  Transform invL;
  Transform pghToNifti;
  Transform invVoxScale;
  Transform qfacTrans;
  Quat q;
  Vec4 tmp= {0.0, 0.0, 0.0, 1.0};
  Vec4 tmp2= {0.0, 0.0, 0.0, 1.0};
  double determinant;
  float qfac;
  int i;

  /* We need the transformation from index values to Pgh MRI coordinate
   * locations for the points associated with the blf, blb, brf, and tlf
   * vectors.  If that transform is T, we can write:
   * 
   *             V = TL
   *
   * where V is the matrix for which the coordinates of the corner points
   * form the column vectors and L is the matrix for which the index offsets
   * form the column vectors.  Given V and the inverse of L, we've got T!
   *
   * T is itself the product of several matrices:
   * 
   *     T = GRSQ
   *
   * where G is the transformation from NIfTI to Pgh coords, R
   * is the rotation given by the NIfTI quaternion, S is a
   * diagonal matrix which applies scaling by voxel size, and
   * Q is a matrix applying the 'qfac' transformation which NIfTI
   * uses to make the coordinate system right-handed.
   */

  /* Fill out the V and L transforms if possible; otherwise bail */
  if (kvLookup(info,"dz") && kvGetInt(info,"dz")>1
      && testVec3(info,"blf") && testVec3(info,"blb") && testVec3(info,"brf")
      && testVec3(info,"tlf")) {

    /* All of these '-0's would be '-1's if the voxel-center-vs-voxel-corner
     * issue were resolved properly.
     */
    /* coords and index offsets of blf */
    getVec3(info, "blf", tmp);
    for (i=0; i<4; i++) V[4*i]= tmp[i];
    L[0]= 0.0;
    L[4]= (double)kvGetInt(info,"dy")-0;
    L[8]= 0.0;
    L[12]= 1.0;
    
    /* coords and index offsets of blf */
    getVec3(info, "blb", tmp);
    for (i=0; i<4; i++) V[(4*i)+1]= tmp[i];
    L[1]= 0.0;
    L[5]= 0.0;
    L[9]= 0.0;
    L[13]= 1.0;
    
    /* coords and index offsets of brf */
    getVec3(info, "brf", tmp);
    for (i=0; i<4; i++) V[(4*i)+2]= tmp[i];
    L[2]= (double)kvGetInt(info,"dx")-0;
    L[6]= (double)kvGetInt(info,"dy")-0;
    L[10]= 0.0;
    L[14]= 1.0;
    
    /* coords and index offsets of tlf */
    getVec3(info, "tlf", tmp);
    for (i=0; i<4; i++) V[(4*i)+3]= tmp[i];
    L[3]= 0.0;
    L[7]= (double)kvGetInt(info,"dy")-0;
    L[11]= (double)kvGetInt(info,"dz")-0;
    L[15]= 1.0;

  }
  else if (testVec3(info,"slice_norm") && testVec3(info,"slice_tlc")
	   && testVec3(info,"slice_trc") && testVec3(info,"slice_blc")
	   && kvLookup(info,"voxel_x") && kvLookup(info,"voxel_y")) {

    /* Slice top left corner is first voxel in slice */
    getVec3(info,"slice_tlc",tmp);
    for (i=0; i<4; i++) V[4*i]= tmp[i];
    L[0]= L[4]= L[8]= 0.0; L[12]= 1.0;

    /* Slice top right corner is last voxel in first row */
    getVec3(info,"slice_trc",tmp);
    for (i=0; i<4; i++) V[(4*i)+1]= tmp[i];
    L[1]= (double)kvGetInt(info,"dx")-0; L[5]= L[9]= 0.0; L[13]= 1.0;

    /* Slice bottom left corner is last voxel in first column */
    getVec3(info,"slice_blc",tmp);
    for (i=0; i<4; i++) V[(4*i)+2]= tmp[i];
    L[6]= (double)kvGetInt(info,"dy")-0; L[2]= L[10]= 0.0; L[14]= 1.0;

    /* For one additional point out-of-plane, translate the slice_tlc
     * along the slice normal direction.
     */
    getVec3(info,"slice_norm",tmp);
    normalizeVec3(tmp); /* just in case */
    getVec3(info,"slice_tlc",tmp2);
    for (i=0; i<4; i++) V[(4*i)+3]= tmp[i]+tmp2[i];
    L[3]= L[7]= 0.0; L[11]= L[15]= 1.0;
  }
  else return 0; /* We don't have enough orientation info to proceed */

  if (kvLookup(info,"nifti_qform_code"))
    hdr->qform_code= (short)kvGetInt(info,"nifti_qform_code");
  else hdr->qform_code= NIFTI_XFORM_SCANNER_ANAT;

  if (debug) {
    fprintf(stderr,"Setting up for quaternion orientation\n");
    fprintf(stderr,"V matrix:\n");
    trans_dump(stderr,V);
    fprintf(stderr,"L matrix:\n");
    trans_dump(stderr,L);
  }

  if (!trans_inverse(invL, L))
    Abort("%s: inversion of offset matrix L failed!\n",progname);
  if (debug) {
    fprintf(stderr,"Inverse of L matrix:\n");
    trans_dump(stderr,invL);
  }

  trans_copy(T, V);
  trans_mult_right(T,invL);
  if (debug) {
    fprintf(stderr,"Offset transformation in Pgh MRI coord system:\n");
    trans_dump(stderr,T);
  }
  
  /* Conversion from Pgh MRI coords to NIfTI coordinates */
  /* Note that it is its own inverse, that is, it's idempotent */
  /* pghToNifti is the matrix inverse of G above */
  trans_identity( pghToNifti );
  pghToNifti[0]= -1.0; /* right-to-left vs. left-to-right */
  pghToNifti[10]= -1.0; /* superior-to-inferior vs. inferior-to-superior */

  trans_mult_left(pghToNifti,T);
  if (debug) {
    fprintf(stderr,"Offset transformation in NIfTI coordinate system:\n");
    trans_dump(stderr,T);
  }

  /* Set and remove the qoffset values */
  hdr->qoffset_x= (float)T[3];
  T[3]= 0.0;
  hdr->qoffset_y= (float)T[7];
  T[7]= 0.0;
  hdr->qoffset_z= (float)T[11];
  T[11]= 0.0;

  /* The determinant of the upper 3x3 submatrix of T now gives us qfac */
  /* qfacTrans is idempotent, so we can use it as its own inverse */
  determinant= (T[0]*T[5]*T[10] + T[1]*T[6]*T[8] + T[2]*T[4]*T[9])
    - (T[8]*T[5]*T[2] + T[4]*T[1]*T[10] + T[0]*T[9]*T[6]);
  if (determinant<0.0) qfac= -1.0;
  else qfac= 1.0;
  if (debug) fprintf(stderr,"determinant= %g -> qfac= %g\n",determinant,qfac);
  hdr->pixdim[0]= qfac;
  trans_identity(qfacTrans);
  qfacTrans[10] *= qfac;
  trans_mult_right(T,qfacTrans);
  if (debug) {
    fprintf(stderr,"After determinant correction:\n");
    trans_dump(stderr,T);
  }

  /* Scale out voxel size */
  trans_identity(invVoxScale);
  invVoxScale[0] = 1.0/kvGetDouble(info,"voxel_x");
  invVoxScale[5] = 1.0/kvGetDouble(info,"voxel_y");
  invVoxScale[10] = 1.0/kvGetDouble(info,"voxel_z");
  trans_mult_right(T,invVoxScale);
  if (debug) {
    fprintf(stderr,"Without offsets and corrected for voxel scale:\n");
    trans_dump(stderr,T);
  }

  trans_to_quat( &q, T );
  if (debug) fprintf(stderr,"Resulting quaternion: (%g, %g, %g, %g)\n",
		     q.x,q.y,q.z,q.w);
  hdr->quatern_b= (float)q.x;
  hdr->quatern_c= (float)q.y;
  hdr->quatern_d= (float)q.z;

  return 1;
}

static void setQuaternionOrientation2D( KVHash* info, nifti_1_header* hdr)
{
  Abort("Not implemented!\n");
}

static void fillHeader( KVHash* info, nifti_1_header* hdr )
{
  const char* dimstr= kvGetString(info,"dimstr");
  const char* allowedDims= "xyztabc"; /* leading v handled separately */
  int niftiDataType= DT_UNKNOWN;
  const char* here= NULL;
  int i;
  char buf[64];

  /* Fill out the header data structure */

  /* We'll start by zeroing out everything */
  bzero(hdr,sizeof(nifti_1_header));

  /* Work through the relevant fields one by one */
  hdr->sizeof_hdr= 348;
  
  if (kvLookup(info,"cl_nii_mode") && kvGetBoolean(info,"cl_nii_mode")) {
    snprintf(hdr->magic,sizeof(hdr->magic),"n+1");
    hdr->vox_offset= 352.0; /* must be a multiple of 16 >= 352 */
  }
  else {
    snprintf(hdr->magic,sizeof(hdr->magic),"ni1");
    hdr->vox_offset= 0.0; /* must be a multiple of 16 >= 0 */
  }
  
  if (kvLookup(info,"nifti_descrip"))
    snprintf(hdr->descrip,sizeof(hdr->descrip),"%s",
	     kvGetString(info,"nifti_descrip"));
  else if (kvLookup(info,"tag"))
    snprintf(hdr->descrip,sizeof(hdr->descrip),"%s",
	     kvGetString(info,"tag"));
  if (kvLookup(info,"nifti_aux_file"))
    snprintf(hdr->aux_file,sizeof(hdr->aux_file),"%s",
	     kvGetString(info,"nifti_aux_file"));

  /* Only a couple of 'v' data formats are supported by nifti */
  if (*dimstr=='v') {
    if (kvGetInt(info,"dv")==1) {
      /* Ignore it */
    }
    else if (kvGetInt(info,"dv")==2) {
      if (kvGetInt(info,"handler_datatype_out")==SRDR_FLOAT32) {
	niftiDataType= DT_COMPLEX64;
      }
      if (kvGetInt(info,"handler_datatype_out")==SRDR_FLOAT64) {
	niftiDataType= DT_COMPLEX128;
      }
    }
    else if (kvGetInt(info,"dv")==3 
	     && kvGetInt(info,"handler_datatype_out")==SRDR_UINT8) {
      niftiDataType= DT_RGB24;
    }
    else Abort("%s: Vectors of length %d and type %s are not supported by nifti!\n",
	       progname, kvGetInt(info,"dv"), 
	       srdrTypeName[kvGetInt(info,"handler_datatype_out")]);
    dimstr++; /* step past that pesky v */
  }

  if (niftiDataType==DT_UNKNOWN) {
    switch (kvGetInt(info,"handler_datatype_out")) {
    case SRDR_UINT8: niftiDataType= DT_UINT8; break;
    case SRDR_INT16: niftiDataType= DT_INT16; break;
    case SRDR_UINT16: 
      /* Since UINT16 is not a native type, all such data should have been
       * promoted to INT32 on reading.  The Pgh MRI format has no UINT16
       * file datatype.
       */
      Abort("%s: Internal error: no stored dataset should actually be UINT16!\n",
	    progname);
      break;
    case SRDR_INT32: niftiDataType= DT_INT32; break;
    case SRDR_FLOAT32: niftiDataType= DT_FLOAT32; break;
    case SRDR_FLOAT64: niftiDataType= DT_FLOAT64; break;
    case SRDR_INT64: niftiDataType= DT_INT64; break;
    default:
      Abort("%s: internal error: unimplemented datatype %s!\n",
	    progname,srdrTypeName[kvGetInt(info,"handler_datatype_out")]);
    }
  }
  hdr->datatype= niftiDataType;
  switch (niftiDataType) {
  case DT_BINARY:
    hdr->bitpix= 1; /* but we don't support that case! */
    break;
  case DT_UINT8:
  case DT_INT8:
    hdr->bitpix= 8;
    break;
  case DT_INT16:
  case DT_UINT16:
    hdr->bitpix= 16;
    break;
  case DT_FLOAT32:
  case DT_INT32:
  case DT_UINT32:
    hdr->bitpix= 32;
    break;
  case DT_FLOAT64:
  case DT_COMPLEX64:
  case DT_UINT64:
  case DT_INT64:
    hdr->bitpix= 64;
    break;
  case DT_FLOAT128:
  case DT_COMPLEX128:
    hdr->bitpix= 128;
    break;
  case DT_COMPLEX256:
    hdr->bitpix= 256;
    break;
  case DT_RGB24:
    hdr->bitpix= 24;
    break;
  default:
    hdr->bitpix= 0; /* the type was unknown- shouldn't happen */
    break;
  }

  if (strncmp(dimstr,allowedDims,strlen(dimstr)))
    Abort("%s: nifti only supports dimension strings in the pattern %s. Permute or remap?\n",
	  progname, allowedDims);
    
  /* Now we just slog through the dimensions that are present. */
  hdr->dim[0]= strlen(dimstr);
  for (here=dimstr; *here; here++) {
    int offset= (int)(here-dimstr) + 1;
    switch (*here) {
    case 'x':
      {
	hdr->dim[offset]= kvGetInt(info,"dx");
	if (kvLookup(info,"voxel_x")) {
	  hdr->pixdim[offset]= kvGetDouble(info,"voxel_x");
	  hdr->xyzt_units |= NIFTI_UNITS_MM;
	}
      }
      break;
    case 'y':
      {
	hdr->dim[offset]= kvGetInt(info,"dy");
	if (kvLookup(info,"voxel_y")) {
	  hdr->pixdim[offset]= kvGetDouble(info,"voxel_y");
	  hdr->xyzt_units |= NIFTI_UNITS_MM;
	}
      }
      break;
    case 'z':
      {
	hdr->dim[offset]= kvGetInt(info,"dz");
	if (kvLookup(info,"voxel_z")) {
	  hdr->pixdim[offset]= kvGetDouble(info,"voxel_z");
	  hdr->xyzt_units |= NIFTI_UNITS_MM;
	}
      }
      break;
    case 't':
      {
	hdr->dim[offset]= kvGetInt(info,"dt");
	if (kvLookup(info,"TR")) {
	  float sliceTime;
	  if (kvLookup(info,"dz")) 
	    sliceTime= (float)kvGetInt(info,"TR")/kvGetInt(info,"dz");
	  else sliceTime= (float)kvGetInt(info,"TR");
	  hdr->pixdim[offset]= sliceTime;
	  hdr->xyzt_units |= NIFTI_UNITS_USEC;
	}
      }
      break;
    case 'a':
      {
	hdr->dim[offset]= kvGetInt(info,"da");
      }
      break;
    case 'b':
      {
	hdr->dim[offset]= kvGetInt(info,"db");
      }
      break;
    case 'c':
      {
	hdr->dim[offset]= kvGetInt(info,"dc");
      }
      break;
    }
  }

  /* We won't do a scaling transformation on the data, but we should at
   * least recognize it if the user has done so for us.
   */
  if (kvLookup(info,"nifti_scl_slope") && kvLookup(info,"nifti_scl_inter")) {
    hdr->scl_slope= (float)kvGetDouble(info,"nifti_scl_slope");
    hdr->scl_inter= (float)kvGetDouble(info,"nifti_scl_inter");
  }

  /* Intent codes if we've got 'em */
  if (kvLookup(info,"nifti_intent_name")) 
    snprintf(hdr->intent_name,sizeof(hdr->intent_name),"%s",
	     kvGetString(info,"nifti_intent_name"));
  if (kvLookup(info,"nifti_intent_p1"))
    hdr->intent_p1= kvGetDouble(info,"nifti_intent_p1");
  if (kvLookup(info,"nifti_intent_p2"))
    hdr->intent_p2= kvGetDouble(info,"nifti_intent_p2");
  if (kvLookup(info,"nifti_intent_p3"))
    hdr->intent_p3= kvGetDouble(info,"nifti_intent_p3");

  /* Slice ordering information.  Pgh MRI always assumes Z is the
   * slice direction.  Maybe I should regret that?
   */
  /* Pgh MRI doesn't allow specification of slice direction; it's 
   * assumed to be in 'z'.  Also, one might end up reordering the
   * data, which would further confuse the issue.  But we'll do
   * the best we can with slice sequence information.
   *
   * NIfTI slice ordering is defined with respect to NIfTI dimensions,
   * which run inferior to superior.  Pgh MRI dimensions run
   * superior to inferior.
   */
  hdr->dim_info= FPS_INTO_DIM_INFO(0,0,3);
  if (kvLookup(info,"dz") && kvLookup(info,"TR")) {
    hdr->slice_duration= (float)kvGetInt(info,"TR")/(float)kvGetInt(info,"dz");
    if (kvLookup(info,"reorder_pattern")) {
      if (!strcmp(kvGetString(info,"reorder_pattern"),"reversed_sequential"))
	hdr->slice_code= NIFTI_SLICE_SEQ_INC;
      else if (!strcmp(kvGetString(info,"reorder_pattern"),"sequential"))
	hdr->slice_code= NIFTI_SLICE_SEQ_DEC;
      else if (!strcmp(kvGetString(info,"reorder_pattern"),
		       "reversed_even/odd"))
	hdr->slice_code= NIFTI_SLICE_ALT_INC;
      else if (!strcmp(kvGetString(info,"reorder_pattern"),
		       "even/odd"))
	hdr->slice_code= NIFTI_SLICE_ALT_DEC;
      else if (!strcmp(kvGetString(info,"reorder_pattern"),
		       "reversed_odd/even"))
	hdr->slice_code= NIFTI_SLICE_ALT_INC2;
      else if (!strcmp(kvGetString(info,"reorder_pattern"),
		       "odd/even"))
	hdr->slice_code= NIFTI_SLICE_ALT_DEC2;
    }
    hdr->slice_start= 0;
    hdr->slice_end= kvGetInt(info,"dz")-1;
  }

  /* Save the NIfTI affine transformation info if present. */
  if (kvLookup(info,"nifti_sform_code"))
    hdr->sform_code= kvGetInt(info,"nifti_sform_code");
  for (i=0; i<4; i++) {
    snprintf(buf,sizeof(buf),"nifti_srow_x_%d",i);
    if (kvLookup(info,buf))
      hdr->srow_x[i]= (float)kvGetDouble(info,buf);
    snprintf(buf,sizeof(buf),"nifti_srow_y_%d",i);
    if (kvLookup(info,buf))
      hdr->srow_y[i]= (float)kvGetDouble(info,buf);
    snprintf(buf,sizeof(buf),"nifti_srow_z_%d",i);
    if (kvLookup(info,buf))
      hdr->srow_z[i]= kvGetDouble(info,buf);
  }
  
  if (kvLookup(info,"dx") && kvLookup(info,"dy")) {
    /* If all the needed corner coords are present, set up a quaternion-type
     * orientation structure.  Otherwise set up a 'legacy' structure.
     */
    if (!trySetQuaternionOrientation(info,hdr))
      setLegacyOrientation(info,hdr);      
  }
}

int main(int argc, char* argv[])
{
  char inputFileName[512], outFileName[512];
  KVHash* info;
  KVHash* wroteThese= NULL;
  KVHash* goodInfo= NULL;
  KVHash* niftiNames= NULL;
  FILE* niftiHead= NULL;
  FILE* niftiImg= NULL;
  SList* chunkList= NULL;
  SList* processedChunkList= NULL;
  ChunkHandlerPair* chPair= NULL;
  FileHandler* fileHandler= NULL;
  nifti_1_header hdr;

  /* Print version number */
  Message( "# %s\n", rcsid );

  progname= argv[0];

  /* Verify data size assumptions used by the nifti1.h header.
   */
  assert(sizeof(int)==4);
  assert(sizeof(float)==4);
  assert(sizeof(short)==2);

  /* Check to see if help was requested */
  if( ( argc > 1 ) && !strcmp( argv[1], "-help" ) )
    {
      if( argc == 2 )
	Help( "selecttopic" );
      else
	Help( argv[2] );
    }

  /* We'll keep lists of chunks before and after header processing. */
  chunkList= slist_create();
  processedChunkList= slist_create();

  /* (almost) All knowledge will accumulate as key-value pairs */
  info= kvFactory(KV_DEFAULT_SIZE);
  wroteThese= kvFactory(KV_DEFAULT_SIZE);
  niftiNames= kvFactory(KV_DEFAULT_SIZE);
#ifdef never
  loadStringTransTable(niftiNames, niftiNameTransTable);
#endif

  /* Miscellaneous initialization tasks- definition and translation
   * tables, default values, opportunity for libbio to figure out
   * endian-ness of this machine.
   */
  InitBIO();
  initStuff(info);

  /* Parse command line */
  parse_command_line(info, inputFileName, outFileName, argc, argv);

  /* Is the input really a Pgh MRI file? */
  if (strcmp(inputFileName+strlen(inputFileName)-4, ".mri")) {
    if (strlen(inputFileName)>sizeof(inputFileName)-5)
      Abort("%s: Input file name <%s> is too long!\n",
	    progname,inputFileName);
    strcat(inputFileName, ".mri");
  }
  if (!pghmriTester(inputFileName))
    Abort("%s: %s is not a Pittsburgh MRI file!\n",
	  progname, inputFileName);

  /* The pghmri handler will be the first pair on our chunk list. */
  if (!(chPair= (ChunkHandlerPair*)malloc(sizeof(ChunkHandlerPair))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,sizeof(ChunkHandlerPair));
  chPair->handler= pghmriFactory(inputFileName, info);
  chPair->info= info;
  slist_append(chunkList,chPair);

  /* We can now loop over the chunks (perhaps discovering more
   * chunks in the process).  Make sure to snag the FileHandler
   * for the actual data.
   */
  while (!slist_empty(chunkList)) {
    chPair= (ChunkHandlerPair*)slist_pop(chunkList);
    FH_PROCESSHEADER(chPair->handler, chPair->info, chunkList);

    if (!smart_reconcile(chPair->info)) {
      smart_dump(chPair->info);
      Abort("Available information is irreconcilable; can't proceed!\n");
    }
    
    slist_append(processedChunkList, chPair);
    
    /* If we have to do type conversion, we'll put a wrapper around the
     * FileHandler to do so.
     */
    if (kvGetInt(chPair->info,"handler_datatype_out") 
	== kvGetInt(chPair->info,"datatype_out")) {
      /* No action required */
    }
    else {
      if (debug) fprintf(stderr,"Chunk %s needs a converter (%s to %s)\n",
			 kvGetString(chPair->info,"chunkname"),
			 srdrTypeName[kvGetInt(chPair->info,
					       "handler_datatype_out")],
			 srdrTypeName[kvGetInt(chPair->info,
					       "datatype_out")]);
      chPair->handler= 
	convertFactory( chPair->handler, 
			kvGetInt(chPair->info,"handler_datatype_out"),
			kvGetInt(chPair->info,"datatype_out") );
    }
    if (!strcmp(kvGetString(chPair->info,"chunkname"),"images"))
      fileHandler= chPair->handler;
  }
  slist_destroy(chunkList,NULL);

  if (!canTranslate(processedChunkList))
    Abort("%s: %s cannot be translated to NIFTI format.\n",
	  progname, inputFileName);

  /* Unfortunately, useful information is scattered among the
   * several chunks.  Gather it all in one place.  We keep the
   * chunks other than "images" in the chunk list, so relevant
   * info can be extracted from them.
   */
  goodInfo= mergeTags( processedChunkList );

  /* Add this command to the dataset history.  We do this last because
   * history information may have been added by various chunk handlers
   * as well.
   */
  addCLToHistory( goodInfo, argc, argv );

  /* This better be enough info to generate the header and
   * transfer the data!
   */
  if (verbose_flg) {
    smart_dump(goodInfo);
    slist_totop(processedChunkList);
    while (!slist_atend(processedChunkList)) {
      ChunkHandlerPair* chPair= 
	(ChunkHandlerPair*)slist_next(processedChunkList);
      fprintf(stderr,"Supplementary chunk: %s\n",
	      kvGetString(chPair->info, "chunkname"));
    }
    slist_totop(processedChunkList);
  }

  /* Open files- there may actually be only one */
  openOutputFiles(outFileName, goodInfo, &niftiHead, &niftiImg);

  /* Fill out the header info */
  fillHeader(goodInfo, &hdr);

  /* Emit the header to the appropriate file */
  writeNiftiHeader(niftiHead, &hdr);

  /* Move the data... */
  transfer( fileHandler, goodInfo, niftiImg, &hdr );

  /* Close output dataset */
  closeOutputFiles(goodInfo, niftiHead, niftiImg);

  /* Clean up */
  slist_destroy(processedChunkList,NULL);
  kvDestroy(goodInfo);

  fprintf(stderr,"#      Data converted to NIFTI format.\n");
  return 0;

}

