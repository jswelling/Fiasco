/************************************************************
 *                                                          *
 *  nifti_reader.c                                       *
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
 *  Original programming by Mark Fitzgerald  5-96           *
 *  Modified to exclusively use libmri calls for output,    *
 *       Greg Hood (PSC), 9-98                              *
 *  Modified to read header files, and to use LX2           *
 *       resampling stuff, Joel Welling (PSC/Stats), 5-1999 *
 ************************************************************/
/* Notes-
 * -I'm not handling cal_min and cal_max values.
 * -I need to do something better with intent codes.
 * -The NIfTI spec says that coords are for voxel centers, which
 *  implies that the volume should be (n-1)*voxelsize on a side.
 *  This doesn't work out correctly, though, especially when compared
 *  to DICOM.  Thus we're fudging things and putting the far corners
 *  of the volume n*voxelsize away.
 * -Does nifti say anything about slice *gap*?  Should I assume it's 
 *  zero?
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#if (SGI64 || SGI5 || SGIMP)
#include <bstring.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <math.h>
#include <assert.h>

#include "nifti1.h"
#include "mri.h"
#include "bio.h"
#include "fmri.h"
#include "array.h"
#include "stdcrg.h"
#include "misc.h"
#include "smartreader.h"
#include "nr_sub.h"
#include "rcn.h"

static char rcsid[] = "$Id: nifti_reader.c,v 1.9 2007/07/07 18:46:06 welling Exp $";

/* A reader for PNG files */
typedef struct png_data_struct {
  nifti_1_header hdr;
  char* datafileName;
} NiftiData;

static void readNiftiHeader( FILE* f, nifti_1_header* hdr )
{
  /* libbio will abort on errors */
  hdr->sizeof_hdr= FRdInt32(f);
  FRdUInt8Array(f, &(hdr->data_type[0]), 10);
  FRdUInt8Array(f, &(hdr->db_name[0]), 18);
  hdr->extents= FRdInt32(f);
  hdr->session_error= FRdInt16(f);
  hdr->regular= FRdUInt8(f);
  hdr->dim_info= FRdUInt8(f);
  FRdInt16Array(f, &(hdr->dim[0]), 8);
  FRdFloat32Array(f, &(hdr->intent_p1), 3);
  FRdInt16Array(f, &(hdr->intent_code), 4);
  FRdFloat32Array(f, &(hdr->pixdim[0]), 11);
  hdr->slice_end= FRdInt16(f);
  hdr->slice_code= FRdUInt8(f);
  hdr->xyzt_units= FRdUInt8(f);
  FRdFloat32Array(f, &(hdr->cal_max), 4);
  FRdInt32Array(f, &(hdr->glmax), 2);
  FRdUInt8Array(f, &(hdr->descrip[0]), 80);
  FRdUInt8Array(f, &(hdr->aux_file[0]), 24);
  hdr->qform_code= FRdInt16(f);
  hdr->sform_code= FRdInt16(f);
  FRdFloat32Array(f, &(hdr->quatern_b),18);
  FRdUInt8Array(f, &(hdr->intent_name[0]), 16);
  FRdUInt8Array(f, &(hdr->magic[0]), 4);
}

static float scaleDistToMM( float rawDist, char xyzt_units )
{
  /* Pgh MRI always stores distance in mm */
  if (xyzt_units & NIFTI_UNITS_MM) return rawDist;
  else if (xyzt_units & NIFTI_UNITS_METER) return 1000.0*rawDist;
  else if (xyzt_units & NIFTI_UNITS_MICRON) return 0.001*rawDist;
  else return rawDist; /* all we can do is guess */
}

static float scaleTimeToUSec( float rawTime, char xyzt_units )
{
  /* Pgh MRI always stores time as microseconds */
  if (xyzt_units & NIFTI_UNITS_SEC) return 1.0e6*rawTime;
  else if (xyzt_units & NIFTI_UNITS_MSEC) return 1.0e3*rawTime;
  else if (xyzt_units & NIFTI_UNITS_USEC) return rawTime;
  else if (xyzt_units & NIFTI_UNITS_HZ) {
    Warning(1,"%s: Warning: data is supposedly in frequency space, in HZ!\n",
	    progname);
    return(1.0e6*rawTime);
  }
  else if (xyzt_units & NIFTI_UNITS_PPM) {
    Warning(1,"%s: Warning: data is supposedly in frequency space, in PPM!\n",
	    progname);
    return(60.0*1.0e6*rawTime);
  }
  else if (xyzt_units & NIFTI_UNITS_HZ) {
    Warning(1,"%s: Warning: data is supposedly in frequency space, in rad/sec!\n",
	    progname);
    return(1.0e6*rawTime/(2.0*M_PI));
  }
  else return 1.0e6*rawTime; /* all we can do is guess */
}

static void transformAndSet( FileHandler* self, KVHash* info, Transform t,
			     char* name, 
			     double xpix, double ypix, double zpix )
{
  Vec4 v;
  v[0]= xpix;
  v[1]= ypix;
  v[2]= zpix;
  v[3]= 1.0;
  trans_vec_mult(t,v);
  defVec3(info, name, v);
}

static void orientByLegacyMethod( FileHandler* self, KVHash* info )
{
  NiftiData* data= (NiftiData*)self->hook;
  Vec4 v1;
  Vec4 v2;
  Vec4 v3;
  Transform niftiR;
  Transform niftiPreTrans;
  Transform niftiToPgh;
  Transform wholeTrans;
  int qfac;
  double tmp;

  kvDefString(info,"nifti_orient_mode","legacy");

  /* This pre-transformation does the rescaling before application of 
   * the rotation and offset.
   */
  trans_identity( niftiPreTrans );
  niftiPreTrans[0]= data->hdr.pixdim[1];
  niftiPreTrans[5]= data->hdr.pixdim[2];
  niftiPreTrans[10]= data->hdr.pixdim[3];

  /* Conversion from NIfTI coordinates to Pgh MRI coords */
  trans_identity( niftiToPgh );
  niftiToPgh[0]= -1.0; /* left-to-right vs. right-to-left */
  niftiToPgh[10]= -1.0; /* inferior-to-superior vs. superior-to-inferior */

  /* Assemble a transformation containing everything up to this point */
  trans_identity( wholeTrans );
  trans_mult_left( niftiPreTrans, wholeTrans );
  trans_mult_left( niftiToPgh, wholeTrans );

  /* Add a final translation to bring the center of the scan volume
   * to the origin, as is typically done with Pgh MRI files when
   * no other info is available.
   */
  v1[0]= kvGetInt(info,"dx")-1;
  v1[1]= kvGetInt(info,"dy")-1;
  v1[2]= kvGetInt(info,"dz")-1;
  v1[3]= 1.0;
  trans_vec_mult(wholeTrans,v1);
  wholeTrans[3] -= 0.5*v1[0];
  wholeTrans[7] -= 0.5*v1[1];
  wholeTrans[11] -= 0.5*v1[2];

    /* Some diagnostics */
  if (debug) {
    fprintf(stderr,"PreTrans:\n");
    trans_dump(stderr,niftiPreTrans);
    fprintf(stderr,"NiftoToPgh:\n");
    trans_dump(stderr,niftiToPgh);
    fprintf(stderr,"Net transformation:\n");
    trans_dump(stderr,wholeTrans);
  }

  /* Let's get the slice normal as a sanity check. */
  v1[0]= 0.0; /* this is our origin */
  v1[1]= 0.0;
  v1[2]= 0.0;
  v1[3]= 1.0;
  trans_vec_mult(wholeTrans,v1);
  v2[0]= 0.0; /* slice normal comes from edge cross prod at slice top left */
  v2[1]= 0.0;
  v2[2]= 1.0;
  v2[3]= 1.0;
  trans_vec_mult(wholeTrans,v2);
  subtractVec3(v3,v2,v1);
  normalizeVec3(v3);
  defVec3(info,"slice_norm",v3);
  
  /* If this is a 3D volume, we define volume corners.  Otherwise
   * we give some slice corner info.
   */
  if (kvLookup(info,"dz")) {
    /* bottom left front is at voxel (0, dy-1, 0) */
    transformAndSet(self,info,wholeTrans,"blf",
		    0.0, kvGetInt(info,"dy")-1, 0.0);
    
    /* bottom left back is at voxel (0, 0, 0) */
    transformAndSet(self, info, wholeTrans, "blb",
		    0.0, 0.0, 0.0);
    
    /* bottom right front is at voxel (dx-1, dy-1, 0) */
    transformAndSet(self, info, wholeTrans, "brf",
		    kvGetInt(info,"dx")-1, kvGetInt(info,"dy")-1, 0.0);
    
    /* bottom right back is at voxel (dx-1, 0, 0) */
    transformAndSet(self, info, wholeTrans, "brb",
		    kvGetInt(info,"dx")-1, 0.0, 0.0);
    
    /* top left front is at voxel (0, dy-1, dz-1) */
    transformAndSet(self, info, wholeTrans, "tlf",
		    0.0, kvGetInt(info,"dy")-1, kvGetInt(info,"dz")-1);
    
    /* top left back is at voxel (0, 0, dz-1) */
    transformAndSet(self, info, wholeTrans, "tlb",
		    0.0, 0.0, kvGetInt(info,"dz")-1);
    
    /* top right front is at voxel (dx-1, dy-1, dz-1) */
    transformAndSet(self, info, wholeTrans, "trf",
		    kvGetInt(info,"dx")-1, kvGetInt(info,"dy")-1, 
		    kvGetInt(info,"dz")-1);
    
    /* top right back is at voxel (dx-1, 0, dz-1) */
    transformAndSet(self,info,wholeTrans,"trb",
		    kvGetInt(info,"dx")-1, 0.0, kvGetInt(info,"dz")-1);
  }
  else {
    /* slice tlc is the first voxel in the array */
    transformAndSet(self,info,wholeTrans,"slice_tlc",0.0,0.0,0.0);

    /* slice trc is the last voxel in the first row */
    transformAndSet(self,info,wholeTrans,"slice_trc",
		    kvGetInt(info,"dx")-1, 0.0, 0.0);

    /* slice blc is the first voxel in the last row */
    transformAndSet(self,info,wholeTrans,"slice_blc",
		    0.0, kvGetInt(info,"dy")-1, 0.0);

    /* slice brc is the last voxel in the slice */
    transformAndSet(self,info,wholeTrans,"slice_brc",
		    kvGetInt(info,"dx")-1, kvGetInt(info,"dy")-1, 0.0);
  }

}

static void orientByQuaternionMethod( FileHandler* self, KVHash* info )
{
  NiftiData* data= (NiftiData*)self->hook;
  Quat q;
  Vec4 v1;
  Vec4 v2;
  Vec4 v3;
  Transform niftiR;
  Transform niftiPreTrans;
  Transform niftiToPgh;
  Transform invNiftiToPgh;
  Transform wholeTrans;
  int qfac;
  double tmp;

  kvDefString(info,"nifti_orient_mode","quaternion");
  kvDefInt(info,"nifti_qform_code",data->hdr.qform_code);

  /* Extract qfac, quaternion and qoffset from the header */
  if (data->hdr.pixdim[0]==0.0) qfac= 1;
  else qfac= data->hdr.pixdim[0];
  q.x= data->hdr.quatern_b;
  q.y= data->hdr.quatern_c;
  q.z= data->hdr.quatern_d;
  tmp= 1.0-(q.x*q.x + q.y*q.y + q.z*q.z);
  if (tmp<0.0) tmp= 0.0;
  q.w= sqrt(tmp);
  quat_to_trans(niftiR, &q, data->hdr.qoffset_x, data->hdr.qoffset_y, 
		data->hdr.qoffset_z);

  /* This pre-transformation does the rescaling before application of 
   * the rotation and offset.
   */
  trans_identity( niftiPreTrans );
  niftiPreTrans[0]= data->hdr.pixdim[1];
  niftiPreTrans[5]= data->hdr.pixdim[2];
  niftiPreTrans[10]= qfac*data->hdr.pixdim[3];

  /* Conversion from NIfTI coordinates to Pgh MRI coords */
  trans_identity( niftiToPgh );
  niftiToPgh[0]= -1.0; /* left-to-right vs. right-to-left */
  niftiToPgh[10]= -1.0; /* inferior-to-superior vs. superior-to-inferior */
  /* niftiToPgh happens to be idempotent */
  trans_copy(invNiftiToPgh, niftiToPgh);

  /* Assemble a transformation containing everything up to this point */
  trans_identity( wholeTrans );
  trans_mult_left( niftiPreTrans, wholeTrans );
  trans_mult_left( niftiR, wholeTrans );
  trans_mult_left( niftiToPgh, wholeTrans );

    /* Some diagnostics */
  if (debug) {
    fprintf(stderr,"PreTrans:\n");
    trans_dump(stderr,niftiPreTrans);
    fprintf(stderr,"Orienting quaternion: %g %g %g %g\n",q.x,q.y,q.z,q.w);
    fprintf(stderr,"Offset vector: %g %g %g\n",
	    data->hdr.qoffset_x, data->hdr.qoffset_y, data->hdr.qoffset_z);
    fprintf(stderr,"Corresponding transform:\n");
    trans_dump(stderr,niftiR);
    fprintf(stderr,"NiftoToPgh:\n");
    trans_dump(stderr,niftiToPgh);
    fprintf(stderr,"Net transformation:\n");
    trans_dump(stderr,wholeTrans);
  }

  /* Let's get the slice normal as a sanity check. */
  v1[0]= 0.0; /* this is our origin */
  v1[1]= 0.0;
  v1[2]= 0.0;
  v1[3]= 1.0;
  trans_vec_mult(wholeTrans,v1);
  v2[0]= 0.0; 
  v2[1]= 0.0;
  v2[2]= 1.0;
  v2[3]= 1.0;
  trans_vec_mult(wholeTrans,v2);
  subtractVec3(v3,v2,v1);
  normalizeVec3(v3);
  defVec3(info,"slice_norm",v3);
  
  /* If this is a 3D volume, we define volume corners.  Otherwise
   * we give some slice corner info.
   */
  if (kvLookup(info,"dz")) {
    /* bottom left front is at voxel (0, dy-1, 0) */
    transformAndSet(self,info,wholeTrans,"blf",
		    0.0, kvGetInt(info,"dy")-0, 0.0);
    
    /* bottom left back is at voxel (0, 0, 0) */
    transformAndSet(self, info, wholeTrans, "blb",
		    0.0, 0.0, 0.0);
    
    /* bottom right front is at voxel (dx-1, dy-1, 0) */
    transformAndSet(self, info, wholeTrans, "brf",
		    kvGetInt(info,"dx")-0, kvGetInt(info,"dy")-0, 0.0);
    
    /* bottom right back is at voxel (dx-1, 0, 0) */
    transformAndSet(self, info, wholeTrans, "brb",
		    kvGetInt(info,"dx")-0, 0.0, 0.0);
    
    /* top left front is at voxel (0, dy-1, dz-1) */
    transformAndSet(self, info, wholeTrans, "tlf",
		    0.0, kvGetInt(info,"dy")-0, kvGetInt(info,"dz")-0);
    
    /* top left back is at voxel (0, 0, dz-1) */
    transformAndSet(self, info, wholeTrans, "tlb",
		    0.0, 0.0, kvGetInt(info,"dz")-0);
    
    /* top right front is at voxel (dx-1, dy-1, dz-1) */
    transformAndSet(self, info, wholeTrans, "trf",
		    kvGetInt(info,"dx")-0, kvGetInt(info,"dy")-0, 
		    kvGetInt(info,"dz")-0);
    
    /* top right back is at voxel (dx-1, 0, dz-1) */
    transformAndSet(self,info,wholeTrans,"trb",
		    kvGetInt(info,"dx")-0, 0.0, kvGetInt(info,"dz")-0);
  }
  else {
    /* slice tlc is the first voxel in the array */
    transformAndSet(self,info,wholeTrans,"slice_tlc",0.0,0.0,0.0);

    /* slice trc is the last voxel in the first row */
    transformAndSet(self,info,wholeTrans,"slice_trc",
		    kvGetInt(info,"dx")-1, 0.0, 0.0);

    /* slice blc is the first voxel in the last row */
    transformAndSet(self,info,wholeTrans,"slice_trc",
		    0.0, kvGetInt(info,"dy")-1, 0.0);

    /* slice brc is the last voxel in the slice */
    transformAndSet(self,info,wholeTrans,"slice_trc",
		    kvGetInt(info,"dx")-1, kvGetInt(info,"dy")-1, 0.0);
  }

}

static void processHeader( FileHandler* self, KVHash* info, SList* cStack )
{
  NiftiData* data= (NiftiData*)self->hook;
  KVHash* defs;
  KVHash* extNames;
  char buf[512];      /* scratch space */
  char* dimstr;
  
  /* Call the base class method */
  baseProcessHeader( self, info, cStack );

  /* Definitions and extNames added by the base method */
  defs= kvGetHash(info,"definitions");
  extNames= kvGetHash(info,"external_names");
  
  /* We test to see if the data does have the expected endian order.
   */
  if (!(self->file = fopen(self->fileName,"r"))) {
    perror("Error opening header");
    Abort("nifti_reader: unable to read or parse header from <%s>!\n",
	  self->fileName);
  }
  readNiftiHeader(self->file, &(data->hdr));
  if (data->hdr.sizeof_hdr != 348) {
    /* Oops, try it the other way! */
    bio_big_endian_input= (bio_big_endian_input ? 0 : 1);
    if (bigfile_fseek(self->file, 0, SEEK_SET)) {
      perror("seek failed");
      Abort("%s: cannot seek on input file %s!\n", progname, self->fileName);
    }
    readNiftiHeader(self->file, &(data->hdr));
    if (data->hdr.sizeof_hdr != 348) {
      Abort("%s: internal error: nifti_reader unexpectedly can't find file length!\n",
	    progname);
    }
  }
  if (fclose(self->file)) perror("Error closing file");
  self->file= NULL; /* reminder to self that data file may be different */

  kvDefBoolean(info,"big_endian_input",bio_big_endian_input);
  if (debug) fprintf(stderr,"Header indicates %s input\n",
		     bio_big_endian_input ? "bigendian" : "littleendian");

  if (!strcmp(data->hdr.magic,"ni1")) {
    char* dot= rindex(self->fileName,'.');
    char* slash= rindex(self->fileName,'/');
    if (dot) {
      if (slash) {
	if (slash>dot) {
	  /* Some stupid directory name in the path has a dot in it! */
	  snprintf(buf,sizeof(buf),"%s.img",self->fileName);
	}
	else {
	  int nchars= dot-(self->fileName);
	  if (nchars>sizeof(buf)-5) nchars= sizeof(buf)-5;
	  strncpy(buf,self->fileName,nchars);
	  buf[nchars]= '\0';
	  strncat(buf,".img",sizeof(buf));
	}
      }
      else {
	int nchars= (dot-self->fileName);
	if (nchars>sizeof(buf)-1) nchars= sizeof(buf)-5;
	strncpy(buf,self->fileName,nchars);
	buf[nchars]= '\0';
	strncat(buf,".img",sizeof(buf));
      }
    }
    else {
      /* No extension */
      snprintf(buf,sizeof(buf),"%s.img",self->fileName);
    }
    kvDefString(info,"datafile",buf);
  }
  else {
    /* Data is stored in this file */
    kvDefString(info,"datafile",self->fileName);
  }
  kvDefLong(info,"start_offset",(int)rint(data->hdr.vox_offset));
  data->datafileName= strdup(kvGetString(info,"datafile"));

  /* Snag descriptive strings */
  if (strlen(data->hdr.descrip)>0)
    kvDefString(info,"nifti_descrip",data->hdr.descrip);
  if (strlen(data->hdr.aux_file)>0)
    kvDefString(info,"nifti_aux_file",data->hdr.aux_file);

  /* Handle data dimensionality, pixdim and intent_code info */
  /* How embarrassing that I can't see how to do this with one switch */
  switch (data->hdr.dim[0]) {
  case 1: 
    kvDefString(info,"dimstr","x");
    break;
  case 2: 
    kvDefString(info,"dimstr","xy");
    break;
  case 3: 
    kvDefString(info,"dimstr","xyz");
    break;
  case 4: 
    kvDefString(info,"dimstr","xyzt");
    break;
  case 5: 
    kvDefString(info,"dimstr","xyzta");
    break;
  case 6: 
    kvDefString(info,"dimstr","xyztab");
    break;
  case 7: 
    kvDefString(info,"dimstr","xyztabc");
    break;
  default: 
    Abort("%s: nifti_reader identified this as a NIfTI-1 file, but dim[0] is an un-allowed value (%d)\n",
	  progname, data->hdr.dim[0]);
  }

  switch (data->hdr.dim[0]) {
  case 7:
    {
      kvDefInt(info,"dc",data->hdr.dim[7]);
    }
  case 6:
    {
      kvDefInt(info,"db",data->hdr.dim[6]);
    }
  case 5:
    {
      kvDefInt(info,"da",data->hdr.dim[5]);
    }
  case 4:
    {
      kvDefInt(info,"dt",data->hdr.dim[4]);
      kvDefInt(info,"TR",
	       (long)(scaleTimeToUSec(data->hdr.pixdim[4], 
				     data->hdr.xyzt_units)
		      *data->hdr.dim[3]));
      kvDefString(defs,"TR","TR (us)");
    }
  case 3:
    {
      kvDefInt(info,"dz",data->hdr.dim[3]);
      kvDefDouble(info,"voxel_spacing.z",
		  scaleDistToMM(data->hdr.pixdim[3], data->hdr.xyzt_units));
    }
  case 2:
    {
      kvDefInt(info,"dy",data->hdr.dim[2]);
      kvDefDouble(info,"voxel_spacing.y",
		  scaleDistToMM(data->hdr.pixdim[2], data->hdr.xyzt_units));
      kvDefDouble(info,"voxel_y",kvGetDouble(info,"voxel_spacing.y"));
    }
  case 1:
    {
      kvDefInt(info,"dx",data->hdr.dim[1]);
      kvDefDouble(info,"voxel_spacing.x",
		  scaleDistToMM(data->hdr.pixdim[1], data->hdr.xyzt_units));
      kvDefDouble(info,"voxel_x",kvGetDouble(info,"voxel_spacing.x"));
    }
    break;
  }

  /* Data type. Scaled data is always floats (or doubles); otherwise
   * we'll use the stored data type.
   */
  switch (data->hdr.datatype) {
  case NIFTI_TYPE_UINT8:
    {
      kvDefInt(info,"datatype_in",SRDR_UINT8);
      kvDefInt(info,"handler_datatype_out",SRDR_UINT8);
    }
    break;
  case NIFTI_TYPE_INT16:
    {
      kvDefInt(info,"datatype_in",SRDR_INT16);
      kvDefInt(info,"handler_datatype_out",SRDR_INT16);
    }
    break;
  case NIFTI_TYPE_INT32:
    {
      kvDefInt(info,"datatype_in",SRDR_INT32);
      kvDefInt(info,"handler_datatype_out",SRDR_INT32);
    }
    break;
  case NIFTI_TYPE_FLOAT32:
    {
      kvDefInt(info,"datatype_in",SRDR_FLOAT32);
      kvDefInt(info,"handler_datatype_out",SRDR_FLOAT32);
    }
    break;
  case NIFTI_TYPE_COMPLEX64:
    {
      /* Must insert a fast dimension for real/complex */
      snprintf(buf,sizeof(buf),"v%s",kvGetString(info,"dimstr"));
      kvDefString(info,"dimstr",buf);
      kvDefInt(info,"dv",2);
      kvDefInt(info,"datatype_in",SRDR_FLOAT32);
      kvDefInt(info,"handler_datatype_out",SRDR_FLOAT32);
    }
    break;
  case NIFTI_TYPE_FLOAT64:
    {
      kvDefInt(info,"datatype_in",SRDR_FLOAT64);
      kvDefInt(info,"handler_datatype_out",SRDR_FLOAT64);
    }
    break;
  case NIFTI_TYPE_RGB24:
    {
      /* Must insert a fast dimension for RGB */
      snprintf(buf,sizeof(buf),"v%s",kvGetString(info,"dimstr"));
      kvDefString(info,"dimstr",buf);
      kvDefInt(info,"dv",3);
      kvDefInt(info,"datatype_in",SRDR_UINT8);
      kvDefInt(info,"handler_datatype_out",SRDR_UINT8);
    }
    break;
  case NIFTI_TYPE_INT8:
    {
      /* We're just going to treat them like bytes. */
      kvDefInt(info,"datatype_in",SRDR_UINT8);
      kvDefInt(info,"handler_datatype_out",SRDR_UINT8);
    }
    break;
  case NIFTI_TYPE_UINT16:
    {
      kvDefInt(info,"datatype_in",SRDR_UINT16);
      /* The owning package will wrap a 'converter' to run these into INT32s */
      kvDefInt(info,"handler_datatype_out",SRDR_UINT16);
    }
    break;
  case NIFTI_TYPE_INT64:
    {
      kvDefInt(info,"datatype_in",SRDR_INT64);
      kvDefInt(info,"handler_datatype_out",SRDR_INT64);
    }
    break;
  case NIFTI_TYPE_COMPLEX128:
    {
      /* Must insert a fast dimension for real/complex */
      snprintf(buf,sizeof(buf),"v%s",kvGetString(info,"dimstr"));
      kvDefString(info,"dimstr",buf);
      kvDefInt(info,"dv",2);
      kvDefInt(info,"datatype_in",SRDR_FLOAT64);
      kvDefInt(info,"handler_datatype_out",SRDR_FLOAT64);
    }
    break;
  case NIFTI_TYPE_UINT32:
    {
      Abort("%s: nifti reader: input data type UINT32 is not supported!\n",
	    progname);
    }
    break;
  case NIFTI_TYPE_UINT64:
    {
      Abort("%s: nifti reader: input data type UINT64 is not supported!\n",
	    progname);
    }
    break;
  case NIFTI_TYPE_FLOAT128:
    {
      Abort("%s: nifti reader: input data type FLOAT128 is not supported!\n",
	    progname);
    }
    break;
  case NIFTI_TYPE_COMPLEX256:
    {
      Abort("%s: nifti reader: input data type COMPLEX256 is not supported!\n",
	    progname);
    }
    break;
  default:
    {
      Abort("%s: nifti reader: unrecognized data type code %d!\n",
	    progname, data->hdr.datatype);
    }
    break;
  }
  if (data->hdr.scl_slope!=0.0) {
    if (kvGetInt(info,"datatype_in")==SRDR_FLOAT64
	|| kvGetInt(info,"datatype_in")==SRDR_INT64)
      kvDefInt(info,"handler_datatype_out",SRDR_FLOAT64);
    else
      kvDefInt(info,"handler_datatype_out",SRDR_FLOAT32);
  }

  /* It's hard to know what to do with intent codes, but we can at least
   * make that information available.
   */
  if (strlen(data->hdr.intent_name)>0)
    kvDefString(info,"nifti_intent_name",data->hdr.intent_name);
  if (data->hdr.intent_code != NIFTI_INTENT_NONE) {
    kvDefInt(info,"nifti_intent_code",data->hdr.intent_code);
    kvDefDouble(info,"nifti_intent_p1",data->hdr.intent_p1);
    kvDefDouble(info,"nifti_intent_p2",data->hdr.intent_p2);
    kvDefDouble(info,"nifti_intent_p3",data->hdr.intent_p3);
  }

  /* Pgh MRI doesn't allow specification of slice direction; it's 
   * assumed to be in 'z'.  Also, one might end up reordering the
   * data, which would further confuse the issue.  But we'll do
   * the best we can with slice sequence information.
   *
   * NIfTI slice ordering is defined with respect to NIfTI dimensions,
   * which run inferior to superior.  Pgh MRI dimensions run
   * superior to inferior.
   */
  if ((data->hdr.slice_code!=0) && (data->hdr.slice_duration>0)
      && (DIM_INFO_TO_SLICE_DIM(data->hdr.dim_info)!=0)) {
    if (DIM_INFO_TO_SLICE_DIM(data->hdr.dim_info)!=3) {
      Warning(1,"%s: nifti_reader: can't translate slice order info for non-axial slices!\n",
	      progname);
    }
    else if ((data->hdr.slice_start != 0) || 
	     ((data->hdr.slice_end != 0) &&
	      (data->hdr.slice_end!=kvGetInt(info,"dz")-1))) {
      Warning(1,"%s: nifti_reader: lost slice order info; padded slice sequences are untranslatable\n",
	      progname);
    }
    else {
      kvDefBoolean(info,"reorder",1);
      switch (data->hdr.slice_code) {
      case NIFTI_SLICE_UNKNOWN:
	/* Do nothing */
	break;
      case NIFTI_SLICE_SEQ_INC:
	kvDefString(info,"reorder_pattern","reversed_sequential");
	break;
      case NIFTI_SLICE_SEQ_DEC:
	kvDefString(info,"reorder_pattern","sequential");
	break;
      case NIFTI_SLICE_ALT_INC:
	kvDefString(info,"reorder_pattern","reversed_even/odd");
	break;
      case NIFTI_SLICE_ALT_DEC:
	kvDefString(info,"reorder_pattern","even/odd");
	break;
      case NIFTI_SLICE_ALT_INC2:
	kvDefString(info,"reorder_pattern","reversed_odd/even");
	break;
      case NIFTI_SLICE_ALT_DEC2:
	kvDefString(info,"reorder_pattern","odd/even");
	break;
      default:
	Abort("%s: nifti_reader: unknown slice sequence code %d!\n",
	      progname,(int)data->hdr.slice_code);
      }
    }
  }

  /* Save the NIfTI affine transformation info if present.  We can't
   * use it because the Pgh MRI coordinate system is always rectangular.
   */
  if (data->hdr.sform_code != 0) {
    int i;
    kvDefInt(info,"nifti_sform_code",data->hdr.sform_code);
    for (i=0; i<4; i++) { 
      snprintf(buf,sizeof(buf),"nifti_srow_x_%d",i);
      kvDefDouble(info,buf,data->hdr.srow_x[i]);
      snprintf(buf,sizeof(buf),"nifti_srow_y_%d",i);
      kvDefDouble(info,buf,data->hdr.srow_y[i]);
      snprintf(buf,sizeof(buf),"nifti_srow_z_%d",i);
      kvDefDouble(info,buf,data->hdr.srow_z[i]);
    }
  }

  /* Extract coordinate system orientation info */
  if (data->hdr.qform_code>0) {
    orientByQuaternionMethod(self,info);
  }
  else {
    orientByLegacyMethod(self,info);
  }
}

static void niftiReopen( FileHandler* self )
{
  NiftiData* data= (NiftiData*)(self->hook);
  /* We actually want to reopen the img file rather than the hdr. */
  if (!(self->file)) {
    if (!(self->file= fopen(data->datafileName,"r")))
      Abort("%s: unable to open file <%s> for reading!\n",
	    progname,data->datafileName);
  }
}

#define CONVERT_BACKWARDS(inType, outType, buf, n, slope, inter) \
{ \
  long i; \
  inType* tinBuf= (inType*)buf; \
  outType* toutBuf= (outType*)buf; \
  for (i=n; i>=0; i--) toutBuf[i]= slope*tinBuf[i] + inter; \
} 

static void niftiRead( FileHandler* self, KVHash* info,
			 long long offset, long n,
			 SRDR_Datatype datatype_out, void* obuf )
{
  NiftiData* data= (NiftiData*)(self->hook);

  FH_REOPEN(self);

  if (data->hdr.scl_slope == 0.0)
    baseRead( self, info, offset, n, datatype_out, obuf );
  else {
    long i;
    /* Data is no larger than obuf, because of type mapping rules.
     * Read it in, than rescale in place.
     */
    baseRead( self, info, offset, n, kvGetInt(info,"datatype_in"), obuf );
    if (kvGetInt(info,"handler_datatype_out")==SRDR_FLOAT32) {
      switch (kvGetInt(info,"datatype_in")) {
      case SRDR_UINT8:
	CONVERT_BACKWARDS(char, float, obuf, n, 
			  data->hdr.scl_slope, data->hdr.scl_inter);
	break;
      case SRDR_INT16:
	CONVERT_BACKWARDS(short, float, obuf, n, 
			  data->hdr.scl_slope, data->hdr.scl_inter);
	break;
      case SRDR_UINT16:
	{
	  /* Annoying- probably a common case, too */
	  long i;
	  int val;
	  short* shortBuf= (short*)obuf;
	  float* floatBuf= (float*)obuf;
	  for (i=n; i>=0; i--) {
	    if (shortBuf[i]>=0) val= shortBuf[i];
	    else val= 65536 + shortBuf[i];
	    floatBuf[i]= data->hdr.scl_slope*val + data->hdr.scl_inter;
	  }
	}
	break;
      case SRDR_INT32:
	CONVERT_BACKWARDS(int, float, obuf, n, 
			  data->hdr.scl_slope, data->hdr.scl_inter);
	break;
      case SRDR_FLOAT32:
	CONVERT_BACKWARDS(float, float, obuf, n, 
			  data->hdr.scl_slope, data->hdr.scl_inter);
	break;
	default:
	  Abort("%s: nifti_reader internal error: confused about scaling floats!\n",
		progname);
      }
    }
    else if (kvGetInt(info,"handler_datatype_out")==SRDR_FLOAT64) {
      switch (kvGetInt(info,"datatype_in")) {
      case SRDR_FLOAT64:
	CONVERT_BACKWARDS(double, double, obuf, n, 
			  data->hdr.scl_slope, data->hdr.scl_inter);
	break;
      case SRDR_INT64:
	CONVERT_BACKWARDS(double, double, obuf, n, 
			  data->hdr.scl_slope, data->hdr.scl_inter);
	break;
      default:
	Abort("%s: nifti_reader internal error: confused about scaling doubles!\n",
	      progname);
      }
    }
    else Abort("%s: nifti_reader internal error: confused about scaling data!\n",
	       progname);
  }
  
}

static void niftiDestroySelf( FileHandler* self )
{
  NiftiData* data= (NiftiData*)(self->hook);

  if (data->datafileName) free(data->datafileName);
  baseDestroySelf(self);
}

FileHandler* niftiFactory(char* fname, KVHash* info)
{
  FileHandler* result= baseFactory(fname);
  NiftiData* data;

  /* Verify data size assumptions used by the nifti1.h header.
   */
  assert(sizeof(int)==4);
  assert(sizeof(float)==4);
  assert(sizeof(short)==2);

  if (!(data=(NiftiData*)malloc(sizeof(NiftiData))))
    Abort("%s: unable to allocate %d bytes!\n",progname,sizeof(NiftiData));
  data->datafileName= NULL;

  result->hook= data;
  result->processHeader= processHeader;
  result->read= niftiRead;
  result->reopen= niftiReopen;
  result->typeName= strdup( "NIFTI Image" );
  result->destroySelf= niftiDestroySelf;
  return result;
}

static int niftiHeaderTest( const char* filename )
{
  FILE* f;
  char buf[4];
  int ierror= 0;
  nifti_1_header hdr;

  if ((f = fopen(filename,"r"))!=NULL)
    {
      readNiftiHeader(f, &hdr);
      if (fclose(f)) {
	perror("Error closing header");
	ierror=1;
      }
    }
  else {
    perror("Error opening header");
    ierror= 1;
  }
    
  if (ierror) return 0;
  else {
    if (debug) fprintf(stderr,"sizeof_hdr= %d\n",hdr.sizeof_hdr);
    /* This line tests both the correctly ordered code and the
     * wrong-endian version of the correctly ordered code.
     */
    if (hdr.sizeof_hdr != 348 && hdr.sizeof_hdr != 1543569408) return 0;
    if (debug) fprintf(stderr,"magic: <%s>\n",hdr.magic);
    if (strcmp(hdr.magic,"ni1") && strcmp(hdr.magic,"n+1"))
      return 0;
  }
  return 1;
}

int niftiTester(const char* filename)
{
  return niftiHeaderTest(filename);
}

