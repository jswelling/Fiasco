/************************************************************
 *                                                          *
 *  par_translate.c                                         *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 2000 Department of Statistics,         *
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
 *  Original programming by Joel Welling, 6/00              *
 ************************************************************/
/* This package implements smoothing methods. */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <time.h>
#include <math.h>
#include "parameters.h"
#include "fmri.h"
#include "misc.h"
#include "stdcrg.h"

static char rcsid[] = "$Id: par_translate.c,v 1.11 2007/03/21 23:59:16 welling Exp $";

/* Notes-
 * 
 */

#define RAD2DEG (180.0/M_PI)
#define DEG2RAD (M_PI/180.0)

static char* progname;

static Transform identity= {
  1.0, 0.0, 0.0, 0.0,
  0.0, 1.0, 0.0, 0.0,
  0.0, 0.0, 1.0, 0.0,
  0.0, 0.0, 0.0, 1.0
};

static Transform afni3dToEstireg3d= {
  0.0, 1.0, 0.0, 0.0,
  0.0, 0.0,-1.0, 0.0,
 -1.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 1.0  
}; 

static Transform airToEstireg3d= {
  1.0, 0.0, 0.0, 0.0,
  0.0,-1.0, 0.0, 0.0,
  0.0, 0.0, 1.0, 0.0,
  0.0, 0.0, 0.0, 1.0  
}; /* note that this is an inversion! */

double prm_access( PrmStruct* in, int z, int t, int p )
{
  if (z>=in->dz || t >=in->dt || p >= in->npar
      || z<0 || t<0 || p<0) {
    Abort("%s: out-of-range location (%d %d %d) in prm_access!\n",
	  progname, z, t, p);
  }
  return in->data[(t*in->dz + z)*in->npar + p];
}

static PrmStruct* prm_new()
{
  PrmStruct* result;

  if (!(result= (PrmStruct*)malloc(sizeof(PrmStruct))))
    Abort("%s: unable to allocate %d bytes!\n",progname,sizeof(PrmStruct));

  result->data= NULL;
  result->dim_x= result->dim_y= result->dim_z= 0;
  result->dt= result->dz= 0;
  result->npar= result->datasize= 0;
  result->vox_x= result->vox_y= result->vox_z= 0.0;
  result->type=PRM_NONE;

  return result;
}

void prm_free( PrmStruct* s )
{
  free(s->data);
  free(s);
}

static void prm_realloc( PrmStruct* in, int dz_new, int dt_new, int np_new )
{
  double* data_new;
  int new_datasize;
  int z;
  int t;
  int p;
  int i;

  if (dz_new<=in->dz && dt_new<=in->dt && np_new<=in->npar) return;

  new_datasize= dz_new*dt_new*np_new;
  if (!(data_new=(double*)malloc(new_datasize*sizeof(double)))) {
    Abort("%s: unable to allocate %d bytes!\n",new_datasize*sizeof(double));
  }
  for (i=0; i<new_datasize; i++) data_new[i]= 0.0;

  if (in->data != NULL) {
    for (z=0; z<in->dz; z++)
      for (t=0; t<in->dt; t++)
	for (p=0; p<in->npar; p++)
	  data_new[(t*dz_new + z)*np_new + p]= prm_access(in,z,t,p);
  }

  in->dz= dz_new;
  in->dt= dt_new;
  in->npar= np_new;
  in->datasize= new_datasize;
  if (in->data != NULL) free(in->data);
  in->data= data_new;
}

static void copy_struct_info( PrmStruct* out, PrmStruct* in )
{
  /* Copies geometry info but not parameter data */
  out->type= in->type;
  out->dim_x= in->dim_x;
  out->dim_y= in->dim_y;
  out->dim_z= in->dim_z;
  out->npar= in->npar;
  out->vox_x= in->vox_x;
  out->vox_y= in->vox_y;
  out->vox_z= in->vox_z;
  out->dz= out->dt= 0;
  out->datasize= 0;
}

PrmStruct* prm_copy( PrmStruct* in )
{
  PrmStruct* result;
  int i;

  if (!(result= (PrmStruct*)malloc(sizeof(PrmStruct))))
    Abort("%s: unable to allocate %d bytes!\n",progname,sizeof(PrmStruct));
  copy_struct_info(result, in);
  result->dz= in->dz;
  result->dt= in->dt;
  result->datasize= in->datasize;
  if (!(result->data= (double*)malloc(result->datasize*sizeof(double))))
    Abort("%s: unable to allocate %d doubles!\n",progname,
	  result->datasize);
  for (i=0; i<result->datasize; i++) result->data[i]= in->data[i];

  return result;
}

void prm_add( PrmStruct* in, int z, int t, int np, double* data_in )
{
  int i;

  prm_realloc(in, z+1, t+1, np);
  for (i=0; i<np; i++) 
    in->data[(t*in->dz + z)*in->npar + i]= data_in[i];
}

int prm_needs_dims( prm_type type )
{
  switch (type) {
  case PRM_AFNI3D: return 1;
  case PRM_ESTIREG3D: return 1;
  case PRM_AIR: return 0;
  default: return 0;
  }
}

int prm_needs_vox( prm_type type )
{
  switch (type) {
  case PRM_AFNI3D: return 1;
  case PRM_ESTIREG3D: return 1;
  case PRM_AIR: return 0;
  default: return 0;
  }
}

int prm_can_load( prm_type type )
{
  switch (type) {
  case PRM_AFNI3D: return 1;
  case PRM_ESTIREG3D: return 1;
  case PRM_AIR: return 1;
  default: return 0;
  }
}

PrmStruct* prm_load( prm_type type, int dz, int dt, 
		     int dim_x, int dim_y, int dim_z,
		     float vox_x, float vox_y, float vox_z,
		     char* fname )
{
  FILE* fp;
  PrmStruct* result= NULL;
  int linenum;
  int numread;
  double tdata[64]; /* some common scratch space */
  int idata[64]; /* ditto */
  double* parameters= NULL;
  int paramSize= 0;
  char scanline[512];

  result= prm_new();
  result->type= type;

  if (!(fp=fopen(fname,"r"))) {
    Abort("%s: unable to open %s for reading!\n",progname,fname);
  }
  
  switch (type) {

  case PRM_AIR:
    {
      int ready_to_read= 0;

      linenum = 0; /* in this case it's actually block number */
      while ( !feof(fp) && !ferror( fp ) ) 
	{
	  /* Scan a line, waiting for V= */
	  if (fgets(scanline, sizeof(scanline), fp)) {
	    if (ready_to_read) {
	      ready_to_read= 0;
	      if ((numread= sscanf( scanline, "[%lg%lg%lg%lg",
				    &tdata[0],&tdata[1],&tdata[2],&tdata[3]))
		  != 4) {
	        Abort( "%s: Block %ld of %s is corrupt!\n", 
		       progname, linenum+1,fname );
	      }
	      if (!fgets(scanline, sizeof(scanline), fp)) {
	        Abort( "%s: Block %ld of %s is corrupt!\n", 
		       progname,linenum+1,fname );
	      }
	      if ((numread= sscanf( scanline, "%lg%lg%lg%lg",
				    &tdata[4],&tdata[5],&tdata[6],&tdata[7]))
		  != 4) {
	        Abort( "%s: Block %ld of %s is corrupt!\n", 
		       progname,linenum+1,fname );
	      }
	      if (!fgets(scanline, sizeof(scanline), fp)) {
	        Abort( "%s: Block %ld of %s is corrupt!\n", 
		       progname,linenum+1,fname );
	      }
	      if ((numread= sscanf( scanline, "%lg%lg%lg%lg",
				    &tdata[8],&tdata[9],&tdata[10],&tdata[11]))
		  != 4) {
	        Abort( "%s: Block %ld of %s is corrupt!\n", 
		       progname, linenum+1,fname );
	      }
	      if (!fgets(scanline, sizeof(scanline), fp)) {
	        Abort( "%s: Block %ld of %s is corrupt!\n", 
		       progname,linenum+1,fname );
	      }
	      if ((numread= sscanf( scanline, "%lg%lg%lg%lg",
				    &tdata[12],&tdata[13],&tdata[14],
				    &tdata[15]))
		  != 4) {
	        Abort( "%s: Block %ld of %s is corrupt!\n", 
		       progname,linenum+1,fname );
	      }
	      prm_add(result, 0, linenum, 16, tdata);
	      linenum++;
	    }
	    else {
	      /* Wait for next block, grabbing data dims as they pass */
	      if (!strncmp(scanline,"file dimensions:",16)) {
		int dx, dy, dz;
		if (sscanf( scanline,
			    "file dimensions:%d by%d by%d pixels (x,y,z)",
			    &dx, &dy, &dz ) != 3) {
		  Abort( "%s: Block %ld of %s is corrupt!\n", 
			 progname,linenum+1,fname );
		}
		if (result->dim_x == 0) result->dim_x= dx;
		else if (result->dim_x != dx) {
		  Abort( "%s: data dimensions change at block %ld of %s!\n", 
			 progname,linenum+1,fname );
		}
		if (result->dim_y == 0) result->dim_y= dy;
		else if (result->dim_y != dy) {
		  Abort( "%s: data dimensions change at block %ld of %s!\n", 
			 progname,linenum+1,fname );
		}
		if (result->dim_z == 0) result->dim_z= dz;
		else if (result->dim_z != dz) {
		  Abort( "%s: data dimensions change at block %ld of %s!\n", 
			 progname,linenum+1,fname );
		}
	      }
	      if (!strncmp(scanline,"voxel dimensions:",17)) {
		float vx, vy, vz;
		if (sscanf( scanline,
			    "voxel dimensions:%g by%g by%g (x,y,z)",
			    &vx, &vy, &vz ) != 3) {
		  Abort( "%s: Block %ld of %s is corrupt!\n", 
			 progname,linenum+1,fname );
		}
		if (result->vox_x == 0.0) result->vox_x= vx;
		else if (result->vox_x != vx) {
		  Abort( "%s: voxel dimensions change at block %ld of %s!\n", 
			 progname,linenum+1,fname );
		}
		if (result->vox_y == 0.0) result->vox_y= vy;
		else if (result->vox_y != vy) {
		  Abort( "%s: voxel dimensions change at block %ld of %s!\n", 
			 progname,linenum+1,fname );
		}
		if (result->vox_z == 0.0) result->vox_z= vz;
		else if (result->vox_z != vz) {
		  Abort( "%s: voxel dimensions change at block %ld of %s!\n", 
			 progname,linenum+1,fname );
		}
	      }
	      else if (!strncmp(scanline,"V=",2)) {
		ready_to_read= 1;
	      }
	    }
	  }
	}

    }
    break;

  case PRM_AFNI3D:
    {
      linenum = -1;
      while( !feof( fp ) && !ferror( fp ) )
	{
	  linenum++;
	  
	  /* Scan a line */
	  if (fgets(scanline, sizeof(scanline), fp)) {
	    numread = sscanf( scanline, "%ld%lg%lg%lg%lg%lg%lg%lg%lg",
			      &idata[0],
			      &tdata[0], &tdata[1], &tdata[2], &tdata[3], 
			      &tdata[4], &tdata[5], &tdata[6], &tdata[7] );
	    if( numread < 9 )
	      {
		Warning( 1, "Line %ld of %s is too short (%ld) -- Ignoring.\n",
			 linenum, fname, numread );
		continue;
	      }
	    
	    prm_add(result, 0, idata[0], 8, tdata);
	  }
	}
      /* AFNI files give coordinates in mm.  Their center of rotation is
       * at the geometrical center, AFNI coordinates (dim_x-1)/2 etc.
       */
      if (dim_x==0 || dim_y==0 || dim_z==0
	  || vox_x==0.0 || vox_y==0.0 || vox_z==0.0) {
	Abort("%s: internal error: dimension or vox size not valid!\n",
	      progname);
      }
      result->dim_x= dim_x;
      result->dim_y= dim_y;
      result->dim_z= dim_z;
      result->vox_x= vox_x;
      result->vox_y= vox_y;
      result->vox_z= vox_z;
    }
    break;
  case PRM_ESTIREG3D:
    {
      int neg_sign= 0;

      linenum = -1;
      while( !feof( fp ) && !ferror( fp ) )
	{
	  linenum++;
	  
	  /* Scan a line, ignoring comments (which begin with '#') */
	  if (fgets(scanline, sizeof(scanline), fp)
	      && strlen(scanline)>0 && scanline[0] != '#') {
	    numread = sscanf( scanline, "%ld%lg%lg%lg%lg%lg%lg%lg%lg",
			      &idata[0],
			      &tdata[0], &tdata[1], &tdata[2], &tdata[3], 
			      &tdata[4], &tdata[5], &tdata[6], &tdata[7] );
	    if( numread < 9 )
	      {
		Warning( 1, "Line %ld of %s is too short (%ld) -- Ignoring.\n",
			 linenum, fname, numread );
		continue;
	      }
	    
	    /* Normalize quaternion.  Since abs(w) is near 1, we assume the
	     * values in x, y, and z are more accurate. */
	    neg_sign= ( tdata[3]<0.0 );
	    tdata[3]= 
	      sqrt( 1.0 - 
		    (tdata[0]*tdata[0]+tdata[1]*tdata[1]+tdata[2]*tdata[2]) );
	    if (neg_sign) tdata[3] *= -1.0;

	    prm_add(result, 0, idata[0], 8, tdata);
	  }
	}

      /* FIASCO files give coordinates in voxels.  Their center of rotation is
       * at the Fourier center, Fiasco coordinates dim_x/2 etc.
       */
      if (dim_x==0 || dim_y==0 || dim_z==0
	  || vox_x==0.0 || vox_y==0.0 || vox_z==0.0) {
	Abort("%s: internal error: dimension or vox size not valid!\n",
	      progname);
      }
      result->dim_x= dim_x;
      result->dim_y= dim_y;
      result->dim_z= dim_z;
      result->vox_x= vox_x;
      result->vox_y= vox_y;
      result->vox_z= vox_z;
    }
    break;
  case PRM_NONE:
    {
    }
    break;
  }

  (void)fclose(fp);
  return result;
}

int prm_can_save( prm_type type )
{
  switch (type) {
  case PRM_AFNI3D: return 1;
  case PRM_ESTIREG3D: return 1;
  default: return 0;
  }
}

void prm_save( PrmStruct* ps, char* fname ) 
{
  FILE* fp;

  if (!(fp=fopen(fname,"w"))) {
    Abort("%s: unable to open %s for writing!\n",progname,fname);
  }
  
  switch (ps->type) {
  case PRM_AFNI3D:
    {
      int t;
      for (t=0; t<ps->dt; t++) {
	fprintf(fp,
     "%d %11.5lg %11.5lg %11.5lg %11.5lg %11.5lg %11.5lg %11.5lg %11.5lg\n",
		t,
		prm_access(ps,0,t,0), prm_access(ps,0,t,1),
		prm_access(ps,0,t,2), prm_access(ps,0,t,3),
		prm_access(ps,0,t,4), prm_access(ps,0,t,5),
		prm_access(ps,0,t,6), prm_access(ps,0,t,7));
      }
    }
    break;
  case PRM_ESTIREG3D:
    {
      time_t tm;
      int t;
      
      tm= time(NULL);
      fprintf(fp,"##Format: order:index_t, type:raw\n");
      fprintf(fp,"##Format: names:(3d_q_x,3d_q_y,3d_q_z,3d_q_w,");
      fprintf(fp,"3d_deltax,3d_deltay,3d_deltaz,mse)\n");
      fprintf(fp,"# Alignment parameters translated %s",
	      asctime(localtime(&tm)));

      for (t=0; t<ps->dt; t++) {
	fprintf(fp,
     "%d %11.5lg %11.5lg %11.5lg %11.5lg %11.5lg %11.5lg %11.5lg %11.5lg\n",
		t,
		prm_access(ps,0,t,0), prm_access(ps,0,t,1),
		prm_access(ps,0,t,2), prm_access(ps,0,t,3),
		prm_access(ps,0,t,4), prm_access(ps,0,t,5),
		prm_access(ps,0,t,6), prm_access(ps,0,t,7));
      }
    }
    break;
  case PRM_NONE:
    {
    }
    break;
  }

  (void)fclose(fp);
}

static void get_normalizing_trans( prm_type type, Transform t )
{
  int i;

  switch (type) {
  case PRM_ESTIREG3D:
    for (i=0; i<16; i++) t[i]= identity[i];
    break;
  case PRM_AIR:
    for (i=0; i<16; i++) t[i]= airToEstireg3d[i];
    break;
  case PRM_AFNI3D:
    for (i=0; i<16; i++) t[i]= afni3dToEstireg3d[i];
    break;
  default:
    Abort("%s: internal error: no normalizing transformation for type %s!\n",
	  progname,prm_type_to_name(type));
  }
}

static void normalize_transform( prm_type type, Transform tInOut,
				 float vox_x, float vox_y, float vox_z )
{
  /* This routine is in charge of taking a transformation matrix in a
   * package's native coordinates and transforming it to Fiasco coordinates 
   * with displacements in mm.  The voxel sizes given are ordered for
   * package native coordinate axes.
   */
  int i;
  Transform t;
  Transform t2;

  for (i=0; i<16; i++) t[i]= tInOut[i]; /* scratch copy */

  get_normalizing_trans(type,t2);

  /* Scale translation as needed, and set up matrix */
  switch (type) {
  case PRM_ESTIREG3D:
    {
      /* displacements given in voxels; compensate for possible anisotropy */
      t[3] *= vox_x;
      t[7] *= vox_y;
      t[11] *= vox_z;
    }
    break;

  case PRM_AFNI3D:
    {
      /* no scaling required; already in mm */
    }
    break;

  case PRM_AIR:
    {
      Vec4 shift;
      /* displacements given in voxels; compensate for possible anisotropy */
      t[3] *= vox_x;
      t[7] *= vox_y;
      t[11] *= vox_z;
      /* AIR seems to give the *inverse* of the expected transform;
       * compensate for this. 
       */
      shift[0]= t[3];
      shift[1]= t[7];
      shift[2]= t[11];
      t[3]= t[7]= t[11]= 0.0;
      trans_transpose(t);
      t[3]= -shift[0];
      t[7]= -shift[1];
      t[11]= -shift[2];
    }
    break;
  default:
    Abort("%s: internal error: cannot normalize coordinates of type %s!\n",
	  progname,prm_type_to_name(type));
  }

  trans_mult_left(t2,t);
  trans_transpose(t2);
  trans_mult_right(t,t2);
  for (i=0; i<16; i++) tInOut[i]= t[i];
  
}

static void denormalize_transform( prm_type type, Transform tInOut, 
				   float vox_x, float vox_y, float vox_z )
{
  /* This routine is in charge of taking a transformation matrix in
   * Fiasco coordinates with displacements in mm and translating it to
   * a package's native coordinate system.  The voxel sizes are given
   * in package native coordinate ordering.
   */
  int i;
  Transform t;
  Transform t2;

  for (i=0; i<16; i++) t[i]= tInOut[i]; /* scratch copy */

  get_normalizing_trans(type, t2);
  trans_mult_right(t,t2);
  trans_transpose(t2);
  trans_mult_left(t2,t);

  /* Scale translation as needed, and set up matrix */
  switch (type) {
  case PRM_ESTIREG3D:
    {
      /* displacements given in voxels; compensate for possible anisotropy */
      t[3] /= vox_x;
      t[7] /= vox_y;
      t[11] /= vox_z;
    }
    break;

  case PRM_AFNI3D:
    {
      /* no scaling required; already in mm */
    }
    break;

  case PRM_AIR:
    {
      Vec4 shift;
      /* AIR seems to give the *inverse* of the expected transform;
       * compensate for this. 
       */
      shift[0]= t[3];
      shift[1]= t[7];
      shift[2]= t[11];
      t[3]= t[7]= t[11]= 0.0;
      trans_transpose(t);
      t[3]= -shift[0];
      t[7]= -shift[1];
      t[11]= -shift[2];
      /* displacements given in voxels; compensate for possible anisotropy */
      t[3] /= vox_x;
      t[7] /= vox_y;
      t[11] /= vox_z;
    }
    break;
  default:
    Abort("%s: internal error: cannot normalize coordinates of type %s!\n",
	  progname,prm_type_to_name(type));
  }

  for (i=0; i<16; i++) tInOut[i]= t[i];
  
}

static void transform_vox_info( prm_type typeIn, prm_type typeOut,
				int* dim_x, int* dim_y, int* dim_z,
				float* vox_x, float* vox_y, float* vox_z )
{
  /* This routine translates voxel info between coordinate models in place */
  Transform tIn;
  Transform tOut;
  Vec4 dimVec;
  Vec4 voxVec;

  get_normalizing_trans(typeIn,tIn);
  get_normalizing_trans(typeOut,tOut);
  trans_transpose(tOut);

  dimVec[0]= *dim_x;
  dimVec[1]= *dim_y;
  dimVec[2]= *dim_z;
  dimVec[3]= 1.0;
  voxVec[0]= *vox_x;
  voxVec[1]= *vox_y;
  voxVec[2]= *vox_z;
  voxVec[3]= 1.0;
  trans_vec_mult( tIn, dimVec );
  trans_vec_mult( tOut, dimVec );
  trans_vec_mult( tIn, voxVec );
  trans_vec_mult( tOut, voxVec );

  *dim_x= (int)floor(fabs(dimVec[0])+0.5);
  *dim_y= (int)floor(fabs(dimVec[1])+0.5);
  *dim_z= (int)floor(fabs(dimVec[2])+0.5);
  *vox_x= fabs(voxVec[0]);
  *vox_y= fabs(voxVec[1]);
  *vox_z= fabs(voxVec[2]);
}

static void get_center_shift( prm_type t, Vec4 shift,
			      int dim_x, int dim_y, int dim_z,
			      float vox_x, float vox_y, float vox_z )
{
  /* This routine returns the translation from the center of rotation in
   * the given native format type to the center of rotation in Fiasco
   * coordinates.  Distances in mm.
   */
  switch (t) {
  case PRM_ESTIREG3D:
    {
      shift[0]= shift[1]= shift[2]= 0.0;
      shift[3]= 1.0;
    }
    break;
    
  case PRM_AFNI3D:
    {
      /* AFNI centers on the geometrical center.  Remember reversal of
       * Fiasco's Y direction wrt. data order.
       */
      shift[0]= 0.5*vox_x;
      shift[1]= -0.5*vox_y;
      shift[2]= 0.5*vox_z;
    }
    break;
    
  case PRM_AIR:
    {
      /* AIR centers on the 0,0,0 voxel (first in data order).  Remember 
       * reversal of Fiasco's Y direction wrt. data order.
       */
      shift[0]= 0.5*dim_x*vox_x;
      shift[1]= -0.5*dim_y*vox_y;
      shift[2]= 0.5*dim_z*vox_z;
      shift[3]= 1.0;
    }
    break;
    
  default:
    Abort("%s: internal error: cannot get center shift of type %s!\n",
	  progname,prm_type_to_name(t));
  }
}

PrmStruct* prm_translate( PrmStruct* in, prm_type type )
{
  PrmStruct* result= NULL;

  if (type==in->type) result= prm_copy(in);
  else {
    int t;
    Transform tStd;
    Transform tShiftIn, tShiftOut;
    int dimXStd, dimYStd, dimZStd;
    float voxXStd, voxYStd, voxZStd;
    Vec4 ctrShiftIn, ctrShiftOut;

    result= prm_new();
    result->type= type;

    /* Translate the voxel info first to standard coords and then to output */
    dimXStd= in->dim_x;
    dimYStd= in->dim_y;
    dimZStd= in->dim_z;
    voxXStd= in->vox_x;
    voxYStd= in->vox_y;
    voxZStd= in->vox_z;
    transform_vox_info(in->type, PRM_ESTIREG3D, &dimXStd, &dimYStd, &dimZStd,
		       &voxXStd, &voxYStd, &voxZStd);
    
    result->dim_x= dimXStd;
    result->dim_y= dimYStd;
    result->dim_z= dimZStd;
    result->vox_x= voxXStd;
    result->vox_y= voxYStd;
    result->vox_z= voxZStd;
    transform_vox_info(PRM_ESTIREG3D, type, 
		       &(result->dim_x),&(result->dim_y),&(result->dim_z),
		       &(result->vox_x),&(result->vox_y),&(result->vox_z));

    /* Translate the lot of them */
    for (t= in->dt-1; t>=0; t--) { /* backwards to avoid some mallocs */
      int i;

      /* First, get them into transformations in input coordinates */
      switch (in->type) {
      case PRM_ESTIREG3D:
	{
	  Quat q;
	  q.x= prm_access(in,0,t,0);
	  q.y= prm_access(in,0,t,1);
	  q.z= prm_access(in,0,t,2);
	  q.w= prm_access(in,0,t,3);
	  quat_to_trans(tStd, &q,
			prm_access(in,0,t,4),
			prm_access(in,0,t,5),
			prm_access(in,0,t,6));
	}
	break;
      case PRM_AFNI3D:
	{
	  Quat q;
	  quat_from_euler_RzRyRx( &q, 
				  DEG2RAD*prm_access(in,0,t,0),
				  DEG2RAD*prm_access(in,0,t,1),
				  DEG2RAD*prm_access(in,0,t,2));
	  quat_to_trans(tStd, &q,
			prm_access(in,0,t,3),
			prm_access(in,0,t,4),
			prm_access(in,0,t,5));
	}
	break;
      case PRM_AIR:
	{
	  /* Check that we can decompose this matrix */
	  if ((prm_access(in,0,t,12) != 0.0)
	      || (prm_access(in,0,t,13) != 0.0) 
	      || (prm_access(in,0,t,14) != 0.0) 
	      || (prm_access(in,0,t,15) == 0.0)) {
	    Abort("%s: AIR matrices are not rigid body; cannot translate!\n",
		  progname);
	  }
	  for (i=0; i<16; i++) tStd[i]= prm_access(in,0,t,i);
	}
	break;
      default:
	Abort("%s: internal error: cannot handle 3D coordinates of type %s!\n",
	      progname,prm_type_to_name(in->type));
      }

      /* Normalize */
      normalize_transform(in->type, tStd, in->vox_x, in->vox_y, in->vox_z);

      /* Compensate for differing centers of rotation */
      get_center_shift( in->type, ctrShiftIn,
			dimXStd, dimYStd, dimZStd, voxXStd, voxYStd, voxZStd );
      get_center_shift( result->type, ctrShiftOut,
			dimXStd, dimYStd, dimZStd, voxXStd, voxYStd, voxZStd );
      for (i=0; i<3; i++) ctrShiftOut[i] *= -1.0;
      for (i=0; i<16; i++) tShiftIn[i]= identity[i];
      tShiftIn[3]= ctrShiftIn[0]-ctrShiftOut[0];
      tShiftIn[7]= ctrShiftIn[1]-ctrShiftOut[1];
      tShiftIn[11]= ctrShiftIn[2]-ctrShiftOut[2];
      for (i=0; i<16; i++) tShiftOut[i]= identity[i];
      tShiftOut[3]= ctrShiftOut[0]-ctrShiftIn[0];
      tShiftOut[7]= ctrShiftOut[1]-ctrShiftIn[1];
      tShiftOut[11]= ctrShiftOut[2]-ctrShiftIn[2];
      trans_mult_right(tStd,tShiftIn);

      trans_mult_left(tShiftOut,tStd);

      /* Denormalize */
      denormalize_transform(result->type, tStd, 
			    result->vox_x, result->vox_y, result->vox_z);

      /* Pack the results back up */
      switch (result->type) {
      case PRM_ESTIREG3D:
	{
	  Quat q;
	  double scratch[8];
	  trans_to_quat(&q,tStd);
	  scratch[0]= q.x;
	  scratch[1]= q.y;
	  scratch[2]= q.z;
	  scratch[3]= q.w;
	  scratch[4]= tStd[3];
	  scratch[5]= tStd[7];
	  scratch[6]= tStd[11];
	  if (in->type==PRM_ESTIREG3D) scratch[7]= prm_access(in,0,t,7);
	  else if (in->type==PRM_AFNI3D) scratch[7]= prm_access(in,0,t,7);
	  else scratch[7]= 0.0;
	  prm_add( result, 0, t, 8, scratch );
	}
	break;

      case PRM_AFNI3D:
	{
	  Quat q;
	  double scratch[8];
	  double rX, rY, rZ; /* Euler angles */
	  trans_to_quat(&q,tStd);
	  if (!quat_to_euler_RzRyRx(&q, &rX, &rY, &rZ)) {
	    Error("%s: euler decomposition failed for time %d!\n",
		  progname,t);
	    rX= rY= rZ= 0.0;
	  }
	  scratch[0]= RAD2DEG*rX;
	  scratch[1]= RAD2DEG*rY;
	  scratch[2]= RAD2DEG*rZ;
	  scratch[3]= tStd[3];
	  scratch[4]= tStd[7];
	  scratch[5]= tStd[11];
	  /* mse's before and after */
	  if (in->type==PRM_AFNI3D) {
	    scratch[6]= prm_access(in,0,t,6);
	    scratch[7]= prm_access(in,0,t,7);
	  }
	  else if (in->type==PRM_ESTIREG3D) {
	    scratch[6]= 0.0;
	    scratch[7]= prm_access(in,0,t,7);
	  }
	  else {
	    scratch[6]= 0.0;
	    scratch[7]= 0.0;
	  }
	  prm_add( result, 0, t, 8, scratch );	  
	}
	break;

      case PRM_AIR:
	{
	  int i;
	  prm_add( result, 0, t, 16, tStd );
	}
	break;

      default:
	Abort("%s: internal error: cannot handle 3D coordinates of type %s!\n",
	      progname,prm_type_to_name(result->type));
      }
    }
  }

  return result;
}

prm_type prm_name_to_type( char* name )
{
  if (!strcasecmp(name,"afni3d")) return PRM_AFNI3D;
  else if (!strcasecmp(name,"estireg3d")) return PRM_ESTIREG3D;
  else if (!strcasecmp(name,"air")) return PRM_AIR;
  else return PRM_NONE;
}

char* prm_type_to_name( prm_type type )
{
  switch (type) {
  case PRM_AFNI3D: return "afni3d";
  case PRM_ESTIREG3D: return "estireg3d";
  case PRM_AIR: return "air";
  case PRM_NONE: return "none";
  }
  return "none";
}

int main( int argc, char** argv )
{
  char infile[512], outfile[512], instring[512], outstring[512];
  prm_type informat, outformat;
  int dt= 0;
  int dt_known= 0;
  int dz= 0;
  int dz_known= 0;
  int dim_x= 0;
  int dim_y= 0;
  int dim_z= 0;
  int dims_known= 0;
  float vox_x= 0;
  float vox_y= 0;
  float vox_z= 0;
  int vox_known= 0;
  PrmStruct* inPrm= NULL;
  PrmStruct* outPrm= NULL;

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

  /* Get formats */
  if (!cl_get("i", "%option %s",&instring)) {
    fprintf(stderr,"%s: input format not given.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  if (!cl_get("o", "%option %s",&outstring)) {
    fprintf(stderr,"%s: output format not given.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  dt_known= cl_get("t","%option %d",&dt);
  dz_known= cl_get("z","%option %d",&dz);

  dims_known= ( cl_get("nx","%option %d",&dim_x)
		&& cl_get("ny","%option %d",&dim_y)
		&& cl_get("nz","%option %d",&dim_z) );

  vox_known= ( cl_get("vx","%option %f",&vox_x)
		&& cl_get("vy","%option %f",&vox_y)
		&& cl_get("vz","%option %f",&vox_z) );

  /* Get filenames */
  if (!cl_get("", "%s", infile)) {
    fprintf(stderr,"%s: Input file name not given.\n",argv[0]);
    Help( "usage" );
    exit(-1);
  }
  if (!cl_get("", "%s", outfile)) {
    fprintf(stderr,"%s: Output file name not given.\n",argv[0]);
    Help( "usage" );
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

  if (dims_known && (dim_x<=0 || dim_y<=0 || dim_z<=0)) {
    Abort("%s: invalid grid dimension given!\n",argv[0]);
  }

  if (vox_known && (vox_x<=0.0 || vox_y<=0.0 || vox_z<=0.0)) {
    Abort("%s: invalid voxel size given!\n",argv[0]);
  }

  if ((informat= prm_name_to_type(instring)) == PRM_NONE) {
    Abort("%s: unrecognized input format <%s>\n",argv[0],instring);
  }
  if ((outformat= prm_name_to_type(outstring)) == PRM_NONE) {
    Abort("%s: unrecognized output format <%s>\n",argv[0],outstring);
  }

  if (!prm_can_load(informat)) {
    Abort("%s: Reading of type <%s> is not supported.\n",argv[0],instring);
  }

  if (!prm_can_save(outformat)) {
    Abort("%s: Writing of type <%s> is not supported.\n",argv[0],outstring);
  }

  if (prm_needs_dims(informat) && !dims_known) {
    Abort("%s: This translation requires image dimensions.\n",argv[0]);
  }

  if (prm_needs_vox(informat) && !vox_known) {
    Abort("%s: This translation requires voxel sizes.\n",argv[0]);
  }

  /* Voxel dims and sizes are given in Fiasco coordinates.  We need
   * to translate them to the internal coordinates of the input
   * coordinate system.
   */
  transform_vox_info( PRM_ESTIREG3D, informat, &dim_x, &dim_y, &dim_z,
		      &vox_x, &vox_y, &vox_z );

  if (!(inPrm= prm_load(informat, dz, dt, dim_x, dim_y, dim_z,
			vox_x, vox_y, vox_z, infile))) {
    Abort("%s: format error reading input file <%s> as %s\n",
	  argv[0], infile, instring);
  }

  Message("# translating <%s> to <%s>\n",
	  prm_type_to_name(informat), prm_type_to_name(outformat));

  if (!(outPrm= prm_translate(inPrm, outformat))) {
    Abort("%s: input format %s cannot be translated to %s\n",
	  argv[0], instring, outstring);
  }
  
  prm_save(outPrm, outfile);
  prm_free(inPrm);
  prm_free(outPrm);
  
  return 0;
}
