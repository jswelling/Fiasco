/************************************************************
 *                                                          *
 *  ireg3d.c                                                *
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
 *     2-99: 3D version, Joel Welling                       *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "stdcrg.h"
#include "misc.h"

static char rcsid[] = "$Id: ireg3d.c,v 1.25 2007/03/21 23:52:54 welling Exp $";

/* Access for 3D arrays */
#define MEM_BCK(matrix,nx,ny,nz,x,y,z) matrix[((((z)*ny)+(y))*nx)+(x)]
#define MEM(matrix,nx,ny,nz,x,y,z) matrix[((((x)*ny)+(y))*nz)+(z)]

typedef struct regpar3d_struct {
  Quat q;
  double x;
  double y;
  double z;
} RegPar3D;

static int debug_flag= 0;
static int verbose_flag= 0;
static int fourier_mode= 1; 
static int shear4_flag= 0;
static int shear7_flag= 0;
static int shear13_flag= 0;
static InterpolatorType interpolatorType= INTRP_LINEAR; /*if not fourier_mode*/
static Interpolator* interpolator= NULL;
static int real_flag= 0;
static const char* progname;
static FComplex *c_image = NULL, *c_corr_image = NULL;
static long complex_buffer_size= 0;
static char* check = NULL;
static long check_buffer_size= 0;

static void regpar3d_copy(RegPar3D* out, RegPar3D* in)
{
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

static void load_reg_params(char* parfile, RegPar3D* par, long dt)
{
  FILE *fp = NULL;
  char scanline[512];
  long linenum, numread, tt;
  int inverse_mode_set= 0;
  int inverse_mode= 0;
  int neg_sign;
  int* covered= NULL;
  int i;

  /* Make sure the input file covers all times */
  if (!(covered=(int*)malloc(dt*sizeof(int))))
    Abort("%s: unable to allocate %d bytes!\n",progname,dt*sizeof(int));
  for (i=0; i<dt; i++) covered[i]= 0;

  fp = efopen( parfile, "r" );
  linenum = -1;
  while( !feof( fp ) && !ferror( fp ) )
    {
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
			    &tt, &(p.q.x), &(p.q.y), &(p.q.z), &(p.q.w),
			    &(p.x), &(p.y), &(p.z) );
	  if( numread < 8 )
	    {
	      Warning( 1, "Line %ld of %s is too short (%ld) -- Ignoring.\n",
		       linenum, parfile, numread );
	      continue;
	    }
	  
	  /* Check to see if image and slice numbers are in bounds */
	  if( ( tt < 0 ) || ( tt >= dt ) )
	    {
	      Warning( 1,
		       "Image out of bounds (%ld) for line %ld of %s -- Ignoring.\n",
		       tt, linenum, parfile );
	      continue;
	    }
	  
	  /* Normalize quaternion.  Since abs(w) is near 1, we assume the
	   * values in x, y, and z are more accurate. */
	  neg_sign= ( p.q.w<0.0 );
	  p.q.w= sqrt( 1.0 - (p.q.x*p.q.x + p.q.y*p.q.y + p.q.z*p.q.z) );
	  if (neg_sign) p.q.w *= -1.0;
	  
	  if (inverse_mode) regpar3d_invert(&p);

	  if (debug_flag)
	    fprintf(stderr,"%d: loaded (%g %g %g %g) %g %g %g\n",
		    tt,p.q.x, p.q.y, p.q.z, p.q.w, p.x, p.y, p.z);

	  /* Put parameters into appropriate storage */
	  regpar3d_copy(&(par[tt]),&p);
	  covered[tt]= 1;
	}
      }
    }
  efclose( fp );
  for (i=0; i<dt; i++) {
    if (!covered[i]) 
      Abort("%s: no input parameters for some times (for example %d)!\n",
	    progname,i);
  }
  free(covered);
}

void checkFourierBufferSize( long dv, long dx, long dy, long dz )
{
  long bufsize= dx*dy*dz;
  long x, y, z;

  if (complex_buffer_size<bufsize) {
    if (c_image) free(c_image);
    if (c_corr_image) free(c_corr_image);
    
    if (!(c_image= (FComplex*)malloc(bufsize*sizeof(FComplex)))) 
      Abort("%s: unable to allocate %d bytes!\n",
	    progname,bufsize*sizeof(FComplex));
    if (!(c_corr_image= (FComplex*)malloc(bufsize*sizeof(FComplex)))) 
      Abort("%s: unable to allocate %d bytes!\n",
	    progname,bufsize*sizeof(FComplex));

    if( dv == 1 ) {
      for(x=0; x<dx; x++)
	for(y=0; y<dy; y++)
	  for (z = 0; z<dz; z++)
	    MEM(c_image,dx,dy,dz,x,y,z).imag= 0.0;
    }
    
    complex_buffer_size= bufsize;
  }
}

void checkCheckBufferSize( long dv, long dx, long dy, long dz )
{
  long bufsize= dx*dy*dz;
  if (check_buffer_size<bufsize) {
    if (check) free(check);
    if (!(check= (char*)malloc(bufsize*sizeof(char))))
      Abort("%s: unable to allocate %d bytes!\n",
	    progname,bufsize*sizeof(char));
    check_buffer_size= bufsize;
  }
}

void applyFourierTransReal( double* image_in, double* corr_image,
			    long dv, long dx, long dy, long dz,
			    double length_x, double length_y, 
			    double length_z,
			    RegPar3D* par)
{
  long x, y, z;
  checkFourierBufferSize( dv, dx, dy, dz );
  for (z=0; z<dz; z++) 
    for(y=0; y<dy; y++)
      for(x=0; x<dx; x++)
	MEM(c_image,dx,dy,dz,x,y,z).real= MEM_BCK(image_in,dx,dy,dz,x,y,z);

  /* Perform registration. */
  fourier_shift_rot3d( &(par->q), par->x, par->y, par->z,
		       c_image, c_corr_image, dx, dy, dz,
		       length_x, length_y, length_z, 1 );

  /* Permute the registered images.  This requires flipping the
   * data back to x-fastest order.
   */
  /* If i-space, input was i-space, so result must be real as well */
  for(x=0; x<dx; x++) 
    for(y=0; y<dy; y++)
      for (z=0; z<dz; z++) {
	MEM_BCK(corr_image,dx,dy,dz,x,y,z)= 
	  MEM(c_corr_image,dx,dy,dz,x,y,z).real;
      }
}

void applyFourierTransComplex( FComplex* c_image_in, FComplex* c_image_out,
			       long dv, long dx, long dy, long dz,
			       double length_x, double length_y, 
			       double length_z,
			       RegPar3D* par)
{
  long x, y, z;
  checkFourierBufferSize( dv, dx, dy, dz );
  for (z=0; z<dz; z++) 
    for(y=0; y<dy; y++)
      for(x=0; x<dx; x++) {
	MEM(c_image,dx,dy,dz,x,y,z).real= 
	  MEM_BCK(c_image_in,dx,dy,dz,x,y,z).real;
	MEM(c_image,dx,dy,dz,x,y,z).imag= 
	  MEM_BCK(c_image_in,dx,dy,dz,x,y,z).imag;
      }
  /* Perform registration.  If our input data was complex, we
   * make the assumption that the shift is to be implemented
   * by a future conversion to image space- only the phases
   * needed for the shift are set.
   */
  /* Must do FFT in Z, since data is only in k-space in XY */
  fft3d( c_image, dx, dy, dz, +1, "z" );
  fourier_shift_rot3d( &(par->q), 0.0, 0.0, 0.0,
		       c_image, c_corr_image, dx, dy, dz,
		       length_x, length_y, length_z, real_flag );
  fshrot3d_set_shift_phases( par->x, par->y, par->z,
			     c_corr_image, dx, dy, dz,
			     length_x, length_y, length_z, 
			     1, real_flag );
  /* Undo earlier Z FFT */
  fft3d( c_corr_image, dx, dy, dz, -1, "z" );

  /* Write out registered image.  This requires flipping the
   * data back to x-fastest order;  we use the other complex
   * buffer as space to do this in the complex case.
   */
  /* here we permute c_corr_image into c_image for output */
  for(x=0; x<dx; x++) 
    for(y=0; y<dy; y++)
      for (z=0; z<dz; z++) {
	MEM_BCK(c_image_out,dx,dy,dz,x,y,z).real=
	  MEM(c_corr_image,dx,dy,dz,x,y,z).real;
	MEM_BCK(c_image_out,dx,dy,dz,x,y,z).imag=
	  MEM(c_corr_image,dx,dy,dz,x,y,z).imag;
      }
}

void applyInterpolatorTransReal( double* image_in, double* image_out,
				 long dv, long dx, long dy, long dz,
				 double length_x, double length_y, 
				 double length_z,
				 RegPar3D* par)
{
  Transform trans;
  Transform invTrans;
  double xTrans, yTrans, zTrans;
  checkCheckBufferSize(dv,dx,dy,dz);
  /* Remember that displacements are in voxels- we must rescale
   * to mm, because that is what intrp_warpApply expects.
   */
  xTrans= par->x*(length_x/dx);
  yTrans= par->y*(length_y/dy);
  zTrans= par->z*(length_z/dz);

  quat_to_trans( trans, &(par->q), xTrans, yTrans, zTrans );
  if (!trans_inverse(invTrans,trans))
    Abort("%s: unable to take the inverse of a transformation matrix!\n",
	  progname);
#ifdef never
  fprintf(stderr,"Transform follows: \n");
  trans_dump(stderr,trans);
  fprintf(stderr,"Inverse transform follows: \n");
  trans_dump(stderr,invTrans);
#endif
  intrp_warpApply( interpolator, invTrans, image_in, image_out, check, 
		   dx, dy, dz, 1, length_x, length_y, length_z );
}

void applyInterpolatorTransComplex( FComplex* image_in, FComplex* corr_image,
				    long dv, long dx, long dy, long dz,
				    double length_x, double length_y, 
				    double length_z,
				    RegPar3D* par)
{
  Abort("%s: internal error: complex data is only supported in fourier mode!\n",
	progname);
}

static void printProgress(long t, long dt)
{
  if (dt>1) {
    /* Print progress report after completing an image */
    if( !t )
      Message( "      " );
    if( ( t == ( dt - 1 ) ) || !( ( t + 1 ) % 60 ) )
      {
	Message( "# %ld\n", (long) ( t + 1 ) );
	if( t != ( dt - 1 ) )
	  Message( "      " );
      }
    else
      Message( "#" );
  }
}

int main( int argc, char* argv[] ) 
{

  MRI_Dataset *Input = NULL, *Output = NULL;
  char infile[512], hdrfile[512], outfile[512], parfile[512];
  char modeString[512];
  char qual_string[512];
  long dv=1, dx=1, dy=1, dz=1, dt=1;
  long x, y, t, z;
  RegPar3D* par= NULL;
  Quat q;
  const char* dimstr= NULL;
  double length_x, length_y, length_z;
  double xvoxel, yvoxel, zvoxel;

  progname= argv[0];

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

  if (cl_present("trilinear"))
    Abort("Option -trilinear has been replaced by -interp linear.  Please see help file.\n");

  /* Get filenames */
  strcpy(outfile,".dat"); /* we no longer let the user set datafile name */
  cl_get( "headerout|h", "%option %s[%]", "ireg3d.mri", hdrfile );
  cl_get( "input|i", "%option %s[%]", "input.mri", infile );
  cl_get( "parameters|p", "%option %s[%]", "reg3d.par", parfile );
  if (!cl_get( "x|xvoxel", "%option %lf", &xvoxel )) {
    fprintf(stderr,"%s: required argument xvoxel omitted.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  if (!cl_get( "y|yvoxel", "%option %lf", &yvoxel )) {
    fprintf(stderr,"%s: required argument yvoxel omitted.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  if (!cl_get( "z|zvoxel", "%option %lf", &zvoxel )) {
    fprintf(stderr,"%s: required argument zvoxel omitted.\n",argv[0]);
    Help("usage");
    exit(-1);
  }
  debug_flag= cl_present("dbg|debug");
  verbose_flag= cl_present("v|verbose");
  real_flag= cl_present("real");
  shear4_flag= cl_present("shear4");
  shear7_flag= cl_present("shear7");
  shear13_flag= cl_present("shear13");
  if (cl_get("interp","%option %s",modeString)) {
    if (strcasecmp(modeString,"fourier")) {
      fourier_mode= 0;
      interpolatorType= intrp_typeFromName(modeString);
      if (interpolatorType==INTRP_UNKNOWN) {
	fprintf(stderr,"%s: unknown mode <%s>.\n",argv[0],modeString);
	Help("usage");
	exit(-1);
      }
    }
    else fourier_mode= 1;
  }
  if (fourier_mode) {
    if (shear4_flag + shear7_flag + shear13_flag > 1) {
      fprintf(stderr,"%s: incompatible flags were specified.\n",argv[0]);
      Help("usage");
      exit(-1);
    }
  }
  else {
    if (shear4_flag || shear7_flag || shear13_flag) {
      fprintf(stderr,"%s: shear pattern specified, but node is not fourier!\n",
	      argv[0]);
      Help("usage");
      exit(-1);
    }
  }
  
  cl_get("qualmeasure","%option %s[%]", "ssqr", qual_string);

  if (cl_cleanup_check()) {
    int i;
    fprintf(stderr,"%s: invalid argument in command line:\n    ",argv[0]);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(-1);
  }

  /*** End command-line parsing ***/

  /* Establish shear4 pattern as default */
  if (fourier_mode && (!shear4_flag && !shear7_flag && !shear13_flag)) {
    shear4_flag= 1;
  }

  /* Print version number */
  if (verbose_flag) Message( "# %s\n", rcsid );

  /* Open input dataset */
  if( !strcmp( infile, hdrfile ) )
    Abort( "Input and output files must be distinct." );
  Input = mri_open_dataset( infile, MRI_READ );

  /* Check that program will function on data-set */
  if( !mri_has( Input, "images" ) )
    Abort( "%s operates only on standard images.", argv[0] );
  if( mri_has( Input, "images.dimensions" ) )
    {
      dimstr= mri_get_string( Input, "images.dimensions" );
      if ( !strcmp( dimstr,"vxyzt") || !strcmp(dimstr,"vxyz") ) {
	if ( mri_has( Input, "images.extent.v" ) ) {
	  dv= mri_get_int(Input,"images.extent.v");
	  if( ( dv < 1 ) || ( dv > 2 ) )
	    Abort( "%s takes only reals or complex numbers of the form (v)xyz(t).",
		   progname );
	  if (dv==2 && !fourier_mode) {
	    Abort( "%s: only fourier mode is supported for complex inputs.",
		   progname );
	  }
	}
	else
	  Abort( "%s: %s does not have the images.extent.v key.",
		 progname, infile );
      }
      else if ( !strcmp( dimstr,"xyzt") || !strcmp(dimstr,"xyz") ) {
	dv= 1;
      }
      else {
	Abort( "%s takes only reals or complex numbers of the form (v)xyzt.",
	       argv[0] );
      }
    }
  else
    Abort( "%s does not have the images.dimensions key.", infile );

  /* Set parameters in local variables */
  if( !mri_has( Input, "images.extent.x" ) ||
      !mri_has( Input, "images.extent.y" ) ||
      !mri_has( Input, "images.extent.z" ) )
    Abort( "images.extent key(s) missing from header." );
  dx = mri_get_int( Input, "images.extent.x" );
  dy = mri_get_int( Input, "images.extent.y" );
  dz = mri_get_int( Input, "images.extent.z" );
  if ( !strcmp(dimstr,"vxyzt") || !strcmp(dimstr,"xyzt") ) {
    if (mri_has(Input, "images.extent.t"))
      dt = mri_get_int( Input, "images.extent.t" );
    else Abort("%s: %s is missing the images.extent.t key.",
	       argv[0],infile);
  }
  else {
    dt= 1;
  }
  if( ( dt <= 0 ) || ( dx <= 0 ) || ( dy <= 0 ) || ( dz <= 0 ) )
    Abort( "images.extent key(s) is non-positive." );

  /* Set output dataset */
  Output = mri_copy_dataset( hdrfile, Input );
  hist_add_cl( Output, argc, argv );
  mri_set_string( Output, "images.datatype", "float32" );
  mri_set_string( Output, "images.file", outfile );
  mri_set_string( Output, "images.dimensions", "vxyzt" );
  mri_set_int( Output, "images.extent.v", (int) dv );
  mri_set_int( Output, "images.extent.t", (int) dt );

  /* Allocate parameter storage */
  if (!(par= (RegPar3D*)malloc(dt*sizeof(RegPar3D))))
    Abort("%s: unable to allocate %d bytes!\n",
	  dt*sizeof(RegPar3D));
  
  /* Initialize registration parameters */
  /*   to zero in case they are missing */
  quat_identity(&q);
  for( t = 0; t < dt; t++ ) {
    quat_copy(&(par[t].q), &q);
    par[t].x= par[t].y= par[t].z= 0.0;
  }

  /* Read in estimated registration parameters */
  load_reg_params( parfile, par, dt );

  /* Set up voxel edge dimensions. */
  length_x= dx*xvoxel;
  length_y= dy*yvoxel;
  length_z= dz*zvoxel;

  /* Method-specific initialization */
  if (fourier_mode) {
    /* Turn on diagnostics if requested, and clear fshrot3d's
     * operation counters.
     */
    fshrot3d_set_debug(debug_flag);
    if (shear7_flag) fshrot3d_set(FR3D_SHEAR_PATTERN, FR3D_SHEAR_7);
    else if (shear13_flag) fshrot3d_set(FR3D_SHEAR_PATTERN, FR3D_SHEAR_13);
    else fshrot3d_set(FR3D_SHEAR_PATTERN, FR3D_SHEAR_4);
    if (!strcasecmp(qual_string,"cox")) 
      fshrot3d_set(FR3D_QUAL_MEASURE, FR3D_QUAL_COX);
    else if (!strcasecmp(qual_string,"sabs")) 
      fshrot3d_set(FR3D_QUAL_MEASURE, FR3D_QUAL_SUM_ABS);
    else if (!strcasecmp(qual_string,"ssqr")) 
      fshrot3d_set(FR3D_QUAL_MEASURE, FR3D_QUAL_SUM_SQR);
    else if (!strcasecmp(qual_string,"ucell")) 
      fshrot3d_set(FR3D_QUAL_MEASURE, FR3D_QUAL_UNIT_CELL);
    else {
      fprintf(stderr,"%s: invalid quality measure name given.\n",argv[0]);
      Help("usage");
      exit(-1);
    }
    fshrot3d_clear_shear_counts();
  }
  else {
    interpolator= intrp_createInterpolator3DByType(interpolatorType,
						   dx,dy,dz,1);
    intrp_warpClearCounts();
    if (debug_flag) {
      interpolator->setInt(interpolator,INTRP_OPT_DEBUG,1);
      interpolator->dumpSelf(interpolator,stderr);
    }
  }

  if (dv==1) {
    double *image_in = NULL, *image_out = NULL;
    long long block_offset= 0;
    long block_size= dx*dy*dz;
    if (!(image_out= (double*)malloc(block_size*sizeof(double))))
      Abort("%s: unable to allocate %d bytes!\n",progname,
	    block_size*sizeof(double));
    for (t=0; t<dt; t++) {
      image_in= (double*)mri_get_chunk(Input, "images", block_size,
				      block_offset, MRI_DOUBLE);
      if (fourier_mode) {
	applyFourierTransReal( image_in, image_out, 
			       dv, dx, dy, dz, 
			       length_x, length_y, length_z,
			       &(par[t]) );
      }
      else {
	applyInterpolatorTransReal( image_in, image_out, 
				    dv, dx, dy, dz, 
				    length_x, length_y, length_z,
				    &(par[t]) );
      }
      mri_set_chunk( Output, "images", block_size, block_offset, MRI_DOUBLE,
		     image_out );
      printProgress(t,dt);
      block_offset += block_size;
    }
    free(image_out);
  }
  else {
    FComplex *c_image_in= NULL, *c_image_out= NULL;
    long long block_offset= 0;
    long block_size= 2*dx*dy*dz;

    if (!(c_image_out= (FComplex*)malloc(dx*dy*dz*sizeof(FComplex))))
      Abort("%s: unable to allocate %d bytes!\n",
	    progname,dx*dy*dz*sizeof(FComplex));
    for (t=0; t<dt; t++) {
      c_image_in= (FComplex*)mri_get_chunk(Input, "images", block_size,
					   block_offset, MRI_FLOAT);
      if (fourier_mode) {
	applyFourierTransComplex( c_image_in, c_image_out,
				  dv, dx, dy, dz,
				  length_x, length_y, length_z,
				  &(par[t]) );
      }
      else {
	applyInterpolatorTransComplex( c_image_in, c_image_out,
				       dv, dx, dy, dz,
				       length_x, length_y, length_z,
				       &(par[t]) );
      }
      mri_set_chunk( Output, "images", block_size, block_offset, MRI_FLOAT,
		     (float*)c_image_out );
      printProgress(t,dt);
      block_offset += block_size;
    }
    free(c_image_out);
  }

  /* Write out op counts */
  if (verbose_flag) {
    if (fourier_mode) {
      int shear_count_x, shear_count_y, shear_count_z, phase_count;
      double shears_per_call;
      fshrot3d_get_shear_counts(&shear_count_x, &shear_count_y, &shear_count_z,
				&phase_count, &shears_per_call);
      Message("# Total shears in x, y, z, rephases: %d %d %d %d (%f per call)\n",
	      shear_count_x,shear_count_y, shear_count_z, phase_count,
	      shears_per_call );
    }
    else {
      long count;
      intrp_warpGetCounts(&count);
      Message("# Total transforms by %s interpolation: %ld\n",
	      intrp_nameFromType(interpolatorType),count );
    }
  }

  /* Write and close data-sets */
  mri_close_dataset( Input );
  mri_close_dataset( Output );
  
  /* Clean up */
  if (interpolator) interpolator->destroySelf(interpolator);
  free(par);

  if (verbose_flag) Message( "#      Image transformation complete.\n" );

  return 0;
}

