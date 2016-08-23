/************************************************************
 *                                                          *
 *  displace3d.c                                            *
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
#include <time.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "misc.h"
#include "stdcrg.h"
#include "slist.h"

static char rcsid[] = "$Id: displace3d.c,v 1.8 2006/10/23 20:59:57 welling Exp $";

#define KEYBUF_SIZE 512

typedef struct regpar3d_struct {
  long t;
  Quat q;
  double x;
  double y;
  double z;
  double mse;
} RegPar3D;

static int debug_flag= 0;
static int dx, dy, dz;
static double xvoxel, yvoxel, zvoxel;
static char* progname;

static void safe_copy(char* str1, char* str2) {
  strncpy(str1, str2, KEYBUF_SIZE);
  str1[KEYBUF_SIZE-1]= '\0';
}

static void safe_concat(char* str1, char* str2) {
  strncat(str1, str2, (KEYBUF_SIZE-strlen(str1))-1);
}

static int safe_get_extent(MRI_Dataset* ds, char* chunk, char* dim)
{
  char key_buf[KEYBUF_SIZE];
  char dim_buf[4];
  dim_buf[0]= *dim;
  dim_buf[1]= '\0';
  safe_copy(key_buf,chunk);
  safe_concat(key_buf,".extent.");
  safe_concat(key_buf,dim_buf);
  if (mri_has(ds,key_buf)) return mri_get_int(ds,key_buf);
  else Abort("%s: input missing tag %s!\n",progname,key_buf);
  return 0; /* not reached */
}

static int weight_valid(MRI_Dataset* weightDS, 
			const int has_xdim,
			const int has_ydim,
			const int has_zdim)
{
  char key_buf[KEYBUF_SIZE];
  char* dimstr;

  if (!mri_has(weightDS,"images.dimensions")) {
    Error("%s: weight file has no images.dimensions tag!\n",progname);
    return 0;
  }
  dimstr= mri_get_string(weightDS,"images.dimensions");
  if (!strncmp(dimstr,"xyz",3)) {
    /* This will definitely work. */
  }
  else if (!strncmp(dimstr,"vxyz",4)) {
    if (safe_get_extent(weightDS,"images","v") != 1) {
      Abort("%s: weight dataset must have v extent 1!\n",progname);
      return 0;
    }
  }
  else {
    Error("%s: weight dataset must have dimensions (v)xyz(...)!\n",progname);
    return 0;
  }

  if (has_xdim && safe_get_extent(weightDS,"images","x")!=dx) {
    Error("%s: weight file x extent doesn't match command line value!\n",
	  progname);
    return 0;
  }
  if (has_ydim && safe_get_extent(weightDS,"images","y")!=dy) {
    Error("%s: weight file y extent doesn't match command line value!\n",
	  progname);
    return 0;
  }
  if (has_zdim && safe_get_extent(weightDS,"images","z")!=dz) {
    Error("%s: weight file y extent doesn't match command line value!\n",
	  progname);
    return 0;
  }

  return 1;
}

static FILE* initParFile(char* parfname, char* inparfname, int filtered_flag)
{
  FILE* result;
  time_t tm;
  if (!(result= fopen(parfname,"w"))) {
    Abort("%s: unable to open <%s> for writing!\n",progname,parfname);
  }

  tm= time(NULL);
  fprintf(result,"##Format: order:index_t, type:%s\n",
	  filtered_flag ? "filtered":"raw");
  fprintf(result,"##Format: names:(3d_rotate,3d_translate,3d_displace)\n");
  fprintf(result,"# Generated at %s",asctime(localtime(&tm)));
  fprintf(result,"# Input registration parameters: %s\n",inparfname);
  fprintf(result,"# dims %d, %d, %d; voxel sizes %g, %g %g\n",
	  dx, dy, dz, xvoxel, yvoxel, zvoxel);
  fflush(result);
  return result;
}

static void regpar3d_copy(RegPar3D* out, RegPar3D* in)
{
  out->t= in->t;
  quat_copy(&(out->q),&(in->q));
  out->x= in->x;
  out->y= in->y;
  out->z= in->z;
}

static void regpar3d_invert(RegPar3D* p)
{
  Transform t;
  Transform tInv;
  Vec4 vec;
  if (debug_flag)
    fprintf(stderr,"Inverting (%g %g %g %g) %g %g %g\n",
	    p->q.x, p->q.y, p->q.z, p->q.w, p->x, p->y, p->z);
  quat_to_trans(t,&(p->q),p->x,p->y,p->z);
  if (!trans_inverse(tInv,t)) 
    Abort("%s: unable to invert a rotation transform!\n",progname);
  trans_to_quat(&(p->q),tInv);
  p->x= tInv[3];
  p->y= tInv[7];
  p->z= tInv[11];
}

static void load_reg_params( char* parfile, SList* parList )
{
  FILE *fp = NULL;
  char scanline[512];
  long linenum, numread;
  int inverse_mode_set= 0;
  int inverse_mode= 0;
  int neg_sign;

  fp = efopen( parfile, "r" );
  linenum = -1;
  while( !feof( fp ) && !ferror( fp ) )
    {
      RegPar3D* thisPar= NULL;
      linenum++;

      /* Scan a line, ignoring comments (which begin with '#') */
      if (fgets(scanline, sizeof(scanline), fp)) {
	if (strlen(scanline)>0 && scanline[0] == '#') {
	  /* Comment- check for format info */
	  const char* formatLoc= strstr(scanline,"##Format:");
	  if (formatLoc) {
	    const char* nameLoc= strstr(formatLoc,"names:");
	    if (nameLoc) {
	      char* start= strchr(nameLoc,'(');
	      if (start) {
		char* end= NULL;
		start++; /* skip off the delimiter */
		end= strchr(start,')');
		if (end) {
		  char* range= strdup(start);
		  int hits= 0;
		  int unhits= 0;
		  range[(int)(end-start)]= '\0';
		  hits += (strstr(range,"3d_qbar_x") ? 1:0);
		  hits += (strstr(range,"3d_qbar_y") ? 1:0);
		  hits += (strstr(range,"3d_qbar_z") ? 1:0);
		  hits += (strstr(range,"3d_qbar_w") ? 1:0);
		  hits += (strstr(range,"3d_deltabarx") ? 1:0);
		  hits += (strstr(range,"3d_deltabary") ? 1:0);
		  hits += (strstr(range,"3d_deltabarz") ? 1:0);
		  unhits += (strstr(range,"3d_q_x") ? 1:0);
		  unhits += (strstr(range,"3d_q_y") ? 1:0);
		  unhits += (strstr(range,"3d_q_z") ? 1:0);
		  unhits += (strstr(range,"3d_q_w") ? 1:0);
		  unhits += (strstr(range,"3d_deltax") ? 1:0);
		  unhits += (strstr(range,"3d_deltay") ? 1:0);
		  unhits += (strstr(range,"3d_deltaz") ? 1:0);
		  free(range);
		  if (hits==7) {
		    inverse_mode= 1;
		    inverse_mode_set= 1;
		  }
		  else if (unhits==7) {
		    inverse_mode= 0;
		    inverse_mode_set= 1;
		  }
		  else Abort("%s: unrecognized field names in par file!\n",
			     progname);
		}
		else
		  Abort("%s: badly formatted Format:names: entry in par file!\n",
			progname);
		
	      }
	      else Abort("%s: badly formatted Format:names: entry in par file!\n",
			 progname);
	    }
	  }
	}
	else {
	  RegPar3D p;
	  if (!inverse_mode_set)
	    Abort("%s: this parameter file lacks needed format information!\n",
		  progname);
	  numread = sscanf( scanline, "%ld%lg%lg%lg%lg%lg%lg%lg%*g",
			    &(p.t), &(p.q.x), &(p.q.y), &(p.q.z), &(p.q.w),
			    &(p.x), &(p.y), &(p.z) );
	  if( numread < 8 )
	    {
	      Warning( 1, "Line %ld of %s is too short (%ld) -- Ignoring.\n",
		       linenum, parfile, numread );
	      continue;
	    }
	  
	  /* Normalize quaternion.  Since abs(w) is near 1, we assume the
	   * values in x, y, and z are more accurate. */
	  neg_sign= ( p.q.w<0.0 );
	  p.q.w= sqrt( 1.0 - (p.q.x*p.q.x + p.q.y*p.q.y + p.q.z*p.q.z) );
	  if (neg_sign) p.q.w *= -1.0;
	  
	  /* Scale the translation info from voxels to mm */
	  p.x *= xvoxel;
	  p.y *= yvoxel;
	  p.z *= zvoxel;

	  /* Invert the quaternion if necessary */
	  if (inverse_mode) regpar3d_invert(&p);

	  if (debug_flag)
	    fprintf(stderr,"%d: loaded (%g %g %g %g) %g/%g %g/%g %g/%g\n",
		    p.t,p.q.x, p.q.y, p.q.z, p.q.w, 
		    p.x, xvoxel, p.y, yvoxel, p.z, zvoxel);

	  /* Put parameters into appropriate storage */
	  if (!(thisPar=(RegPar3D*)malloc(sizeof(RegPar3D))))
	    Abort("%s: unable to allocate %d bytes!\n",progname,
		  sizeof(RegPar3D));
	  regpar3d_copy(thisPar,&p);
	  slist_append(parList,thisPar);
	}
      }
    }
  efclose( fp );
  if (slist_empty(parList))
    Abort("%s: no valid parameter lines read!\n",progname);
}

static void writeOutput( RegPar3D *par, float* weights,
			 FILE* ofp ) {
  double trans;
  double rot; /* in radians */
  double disp;
  double weight_total;
  double rx, ry, rz, theta;
  int count;
  int i, j, k;
  Quat tq;
  Quat qbar;
  double loc[3], rot_loc[3];
  
  /* Calculate the easy bits */
  trans= sqrt( par->x*par->x + par->y*par->y + par->z*par->z );
  quat_to_axis_angle( &(par->q), &rx, &ry, &rz, &rot );

  /* Calculate displacement statistic.  Since the model is that 
   * translations are applied after rotations, we can just add the
   * translation in at the end. 
   */
  quat_copy(&qbar,&(par->q));
  quat_conjugate(&qbar);
  disp= 0.0;
  weight_total= 0.0;
  for (k=0; k<dz; k++) 
    for (j=0; j<dy; j++)
      for (i=0; i<dx; i++) {
	double wt= weights[(k*dy+j)*dx+i];
	tq.x= loc[0] = (float) (xvoxel*( i - ( dx / 2 ) ));    
	tq.y= loc[1] = (float) (yvoxel*( j - ( dy / 2 ) ));    
	tq.z= loc[2] = (float) (zvoxel*( k - ( dz / 2 ) ));    
	tq.w= 0.0;
	quat_mult_left(&(par->q),&tq);
	quat_mult_right(&tq,&qbar);
	rot_loc[0]= tq.x;
	rot_loc[1]= tq.y;
	rot_loc[2]= tq.z;
	rot_loc[0] -= loc[0];
	rot_loc[1] -= loc[1];
	rot_loc[2] -= loc[2];
	disp += wt*sqrt( rot_loc[0]*rot_loc[0] + rot_loc[1]*rot_loc[1] +
			 rot_loc[2]*rot_loc[2] );
	weight_total += wt;
      }
  if (debug_flag) {
    fprintf(stderr,"t= %d: raw disp %lg over %d voxels, total weight %f\n",
	    par->t,disp,dx*dy*dz,weight_total);
    fprintf(stderr,"       trans= %g, rot= %g (radians)\n",trans,rot);
  }
  disp /= weight_total;
  disp += trans;
  
  /* Write 'em out */
  fprintf(ofp,"%d %11.5g %11.5g %11.5g\n",
	  par->t,(180.0/M_PI)*rot,trans,disp);
  
}

int main( int argc, char** argv ) 
{
  char infile[512], parfile[512], wtfile[512];
  FILE *ifp = NULL, *ofp = NULL;
  MRI_Dataset* weightDS= NULL;
  int x;
  int y;
  int z;
  int t;
  int filtered_flag= 0;
  SList* parList= slist_create();
  int recCount= 0;
  int has_weight= 0;
  int has_xdim= 0;
  int has_ydim= 0;
  int has_zdim= 0;
  float* weights= NULL;

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

  if (cl_present( "input|i" ))
    Abort ("Option input|i has been replaced by infile outfile format.  Please see help file.\n");
  if (cl_present( "parameters|p" ))
     Abort ("Option parameters|p has been replaced by estimates|est|e.  Please see help file.\n");
  if (cl_present( "x" ))
     Abort ("Option x has been replaced by xdimension|xdm.  Please see help file.\n");
  if (cl_present( "y" ))
     Abort ("Option y has been replaced by ydimension|ydm.  Please see help file.\n");
  if (cl_present( "z" ))
     Abort ("Option z has been replaced by zdimension|zdm.  Please see help file.\n");

  
  /* Get filenames */
  cl_get( "estimates|est|e", "%option %s[%]", "displace3d.par", parfile );
  has_xdim= cl_get( "xdimension|xdm", "%option %ld", &dx );
  has_ydim= cl_get( "ydimension|ydm", "%option %ld",  &dy );
  has_zdim= cl_get( "zdimension|zdm", "%option %ld", &dz );
  has_weight= cl_get( "weight|wgt|w", "%option %s", wtfile );
  if (!cl_get( "xvx|xvoxel", "%option %lf", &xvoxel )) {
    fprintf(stderr,"%s: required argument xvoxel omitted.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  if (!cl_get( "yvx|yvoxel", "%option %lf", &yvoxel )) {
    fprintf(stderr,"%s: required argument yvoxel omitted.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  if (!cl_get( "zvx|zvoxel", "%option %lf", &zvoxel )) {
  fprintf(stderr,"%s: required argument zvoxel omitted.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  filtered_flag= cl_present("filtered|ftd");
  debug_flag= cl_present("debug|deb");

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

  if ((has_xdim && ( dx <= 0 )) 
      || (has_ydim && ( dy <= 0 )) 
      || (has_zdim && (dz <= 0) )) {
    Abort( "Dimensions not valid (%ld,%ld,%dz).", dx, dy, dz );
  }

  if( ( xvoxel <= 0.0 ) || ( yvoxel <= 0.0 ) || (zvoxel <= 0.0) ) {
    Abort( "Voxel sizes not valid (%g, %g, %g).", xvoxel, yvoxel, zvoxel );
  }

  if (!((has_xdim && has_ydim && has_zdim) || has_weight))
    Abort( "Dimension info not provided (explicitly or by weight file).");
  
  if (has_weight) {
    weightDS= mri_open_dataset(wtfile,MRI_READ);
    if (!weight_valid(weightDS,has_xdim,has_ydim,has_zdim))
      Abort("%s: weight information is not useable.\n",argv[0]);
    if (!has_xdim) dx= safe_get_extent(weightDS,"images","x");
    if (!has_ydim) dy= safe_get_extent(weightDS,"images","y");
    if (!has_zdim) dz= safe_get_extent(weightDS,"images","z");
  }
  else {
    weightDS= NULL;
  }
  
  /* Open files */
  if (debug_flag) fprintf(stderr,"reading <%s>, writing <%s>\n",
			  infile,parfile);
  if (!(ifp= fopen(infile,"r"))) {
    Abort("%s: unable to open <%s> for reading!\n",progname,infile);
  }
  ofp= initParFile(parfile, infile, filtered_flag);

  /* Allocate and load (or fake) weights */
  if (!(weights=(float*)malloc(dx*dy*dz*sizeof(float))))
    Abort("%s: unable to allocate %d bytes!\n",progname);
  if (has_weight) mri_read_chunk(weightDS, "images", dx*dy*dz, 0,
				 MRI_FLOAT, weights);
  else {
    for (z=0; z<dz; z++)
      for (y=0; y<dy; y++)
	for (x=0; x<dx; x++) weights[(z*dy+y)*dx+x]= 1.0;
  }

  /* Read in all input lines */
  load_reg_params( infile, parList );

  /* Generate output, record by record */
  recCount= 0;
  slist_totop(parList);
  while (!slist_atend(parList)) {
    RegPar3D* thisPar= (RegPar3D*)slist_next(parList);
    writeOutput(thisPar,weights,ofp);
    recCount++;
  }

  /* Close files */
  if (fclose(ifp)) {
    perror("Error closing input:");
  }
  if (fclose(ofp)) {
    perror("Error closing output:");
  }
  if (has_weight) mri_close_dataset(weightDS);

  slist_destroy(parList,free);

  Message( "#      Displacement statistics calculated (%d records).\n",
    recCount);

  return 0;
}

