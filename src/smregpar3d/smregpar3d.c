/************************************************************
 *                                                          *
 *  smregpar3d.c                                            *
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
 *     4-00: modified for 3D by Joel Welling
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <time.h>
#include <math.h>
#include "mri.h"
#include "fmri.h"
#include "misc.h"
#include "stdcrg.h"

static char rcsid[] = "$Id: smregpar3d.c,v 1.13 2007/07/30 18:35:48 welling Exp $";

typedef struct regpar3d_struct {
  Quat q;
  double x;
  double y;
  double z;
  double mse;
} RegPar3D;

static char* progname;

static FILE* initParFile(char* parfname, char* infname, char* inparfname,
			 const char* parNameString)
{
  FILE* result;
  time_t tm;
  if (!(result= fopen(parfname,"w"))) {
    Abort("%s: unable to open <%s> for writing!\n",progname,parfname);
  }

  tm= time(NULL);
  fprintf(result,"##Format: order:index_t, type:filtered\n");
  fprintf(result,"##Format: names:(%s)\n",parNameString);
  fprintf(result,"# Alignment parameters smoothed %s",
	  asctime(localtime(&tm)));
  fprintf(result,"# Input file: %s\n",infname);
  fprintf(result,"# Unsmoothed parameters: %s\n",inparfname);
  fflush(result);
  return result;
}

static void saveRegParams(FILE* ofile, RegPar3D* p, long t) {
  fprintf(ofile, "%ld %11.5lg %11.5lg %11.5lg %11.5lg %11.5lg %11.5lg %11.5lg %11.5lg\n",
	  t, p->q.x, p->q.y, p->q.z, p->q.w, p->x, p->y, p->z, p->mse);
  fflush(ofile);
}

static void loadRegParams(char* parfile, RegPar3D* par, long dt, 
			  char** parNameString)
{
  FILE *fp = NULL;
  char scanline[512];
  long linenum, numread, tt;
  Quat q;
  double x_shift, y_shift, z_shift, mse;
  int neg_sign;

  fp = efopen( parfile, "r" );
  linenum = -1;
  while( !feof( fp ) && !ferror( fp ) )
    {
      linenum++;

      /* Scan a line, ignoring comments (which begin with '#') */
      if (fgets(scanline, sizeof(scanline), fp)
	  && strlen(scanline)>0) {
	if (scanline[0]=='#') { /* comment */
	  if (strstr(scanline,"##Format:")) {
	    char* here;
	    if (here=strstr(scanline,"names:")) {
	      char* there= NULL;
	      while (*here && *here!='(') here++;
	      if (*here != '(') 
		Abort("%s: incorrectly formatted format string!\n",progname);
	      here++; /* skip over paren */
	      there= here;
	      while (*there && *there!=')') there++;
	      if (*there != ')') 
		Abort("%s: incorrectly formatted format string!\n",progname);
	      *parNameString= strdup(here);
	      *(*parNameString+(int)(there-here))= '\0';
	      fprintf(stderr,"Got string <%s>\n",*parNameString);
	    }
	  }
	}
	else {
	  numread = sscanf( scanline, "%ld%lg%lg%lg%lg%lg%lg%lg%lg",
			    &tt, &(q.x), &(q.y), &(q.z), &(q.w),
			    &x_shift, &y_shift, &z_shift, &mse );
	  if( numread < 9 )
	    {
	      Warning( 1, "Line %ld of %s is too short (%ld) -- Ignoring.\n",
		       linenum, parfile, numread );
	      continue;
	    }
	  
	  /* Check to see if image number is in bounds */
	  if( ( tt < 0 ) || ( tt >= dt ) )
	    {
	      Warning( 1,
		       "Image out of bounds (%ld) for line %ld of %s -- Ignoring.\n",
		       tt, linenum, parfile );
	      continue;
	    }
	  
	  /* Normalize quaternion.  Since abs(w) is near 1, we assume the
	   * values in x, y, and z are more accurate. */
	  neg_sign= ( q.w<0.0 );
	  q.w= sqrt( 1.0 - (q.x*q.x + q.y*q.y + q.z*q.z) );
	  if (neg_sign) q.w *= -1.0;
	  
	  /* Put parameters into appropriate storage */
	  quat_copy(&(par[tt].q),&q);
	  par[tt].x= x_shift;
	  par[tt].y= y_shift;
	  par[tt].z= z_shift;
	  par[tt].mse= mse;
	}
      }
    }
  efclose( fp );
}

/* This function provides the metric used to check the threshold
 *  cutoff value for smoothing.
 */
static float thresh_value( Smoother* smoother, float** dtbl, int ndata,
			   int n1, int n2 )
{
  float result= 0.0;

  if (n1==n2) return 0.0;
  else {
    /* through the magic of quaternions, this is actually a measure of 
     * mean displacement integrated over the unit sphere.
     */
    float sqdiffx;
    float sqdiffy;
    float sqdiffz;
    float sqw;
    Quat q1, q2;

    sqdiffx= dtbl[4][n1] - dtbl[4][n2];
    sqdiffx *= sqdiffx;
    
    sqdiffy= dtbl[5][n1] - dtbl[5][n2];
    sqdiffy *= sqdiffy;
    
    sqdiffz= dtbl[6][n1] - dtbl[6][n2];
    sqdiffz *= sqdiffz;

    q1.x= dtbl[0][n1]; q1.y= dtbl[1][n1]; q1.z= dtbl[2][n1]; q1.w= dtbl[3][n1];
    q2.x= dtbl[0][n2]; q2.y= dtbl[1][n2]; q2.z= dtbl[2][n2]; q2.w= dtbl[3][n2];
    quat_mult_right(&q1,quat_conjugate(&q2));
    
    sqw= q1.w * q1.w;
    result= sqdiffx + sqdiffy + sqdiffz + 2.0*(1.0-sqw);
    result *= (3.0/(4.0*M_PI)); /* scale for unit sphere */
  }

  return result;
}

int main( int argc, char **argv ) 
{
  MRI_Dataset *Input = NULL;
  FILE *ofp = NULL;
  char infile[512], iparfile[512], oparfile[512], kernel[512];
  float bandwidth;
  float cutoff;
  float threshold;
  unsigned char **missing = NULL, **image_missing= NULL;
  char scanline[512];
  long linenum, numread;
  long dx, dy, dt, dz, t, z, tt, zz;
  int i;
  int whole_image_missing;
  RegPar3D *par = NULL, *corr_par = NULL;
  float *ordmse = NULL, medianmse, iqrmse, msecutoff;
  long msecount;
  long fixed_image, temp_fixed_image;
  float rotparadj, sqdiffx, sqdiffy, sqdiffrot, distance;
  float weight, totalweight;
  Smoother* smoother;
  sm_type smoother_type;
  float* smooth_ibuf= NULL;
  float* smooth_obuf= NULL;
  float* smooth_itbl[8];
  float* smooth_otbl[8];
  int this_image_missing;
  char* parNameString= NULL;

  progname= argv[0];

  /* Print version number */
  Message( "# %s\n", rcsid );

  sm_init();

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
  /* some systems seem to have trouble with setting the default
   * value for fixed_image, so we handle it as a special case.
   */
  cl_get( "headerinput|h", "%option %s[%]", "input.mri", infile );
  cl_get( "parameterin", "%option %s[%]", "reg3d.par", iparfile );
  cl_get( "parameterout", "%option %s[%]", "smreg3d.par", oparfile );
  cl_get( "cutoff|c", "%option %f[%]", 10.0, &cutoff );
  fixed_image= -1;
  cl_get( "fixed|f", "%option %ld", &fixed_image );
  cl_get( "bandwidth|b", "%option %f[%]", 3.0, &bandwidth );
  cl_get( "kernel|k", "%option %s[%]", "gaussian", kernel );
  cl_get( "threshold|t", "%option %f[%]", 0.01, &threshold );

  if( !strcasecmp( kernel, "gaussian" ) || !strcasecmp( kernel, "g" ) )
    {
      smoother_type= SM_GAUSSIAN;
    }
  else if( !strcasecmp( kernel, "triangular" ) || !strcasecmp( kernel, "t" ) )
    {
      smoother_type= SM_TRIANGULAR;
    }
  else if( !strcasecmp( kernel, "exponential" ) || !strcasecmp( kernel, "e" ) )
     {
      smoother_type= SM_POWER;
    }
  else {
    Abort( "Kernel unrecognized (%s).", kernel );
  }

  sm_set_params( smoother_type, bandwidth, 0.0, threshold, 
		 (SmootherThreshTest)thresh_value);
  sm_parse_cl_opts();
  sm_get_params( &smoother_type, &bandwidth, NULL, &threshold, NULL );

  if (cl_cleanup_check()) {
    int i;
    fprintf(stderr,"%s: invalid argument in command line:\n    ",argv[0]);
    for (i=0; i<argc; i++) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"\n");
    Help( "usage" );
    exit(-1);
  }

  /*** End command-line parsing ***/


  /* Open input dataset */
  Input = mri_open_dataset( infile, MRI_MODIFY );
  hist_add_cl(Input,argc,argv); /* since we will update the missing chunk */

  /* Check that program will function on data-set */
  if( !mri_has( Input, "images" ) )
    Abort( "%s operates only on standard images.", argv[0] );
  if( !mri_has( Input, "images.dimensions" ) ||
      ( strcmp( mri_get_string( Input, "images.dimensions" ), "vxyzt" ) &&
	strcmp( mri_get_string( Input, "images.dimensions" ), "xyzt" ) ) )
    Abort( "%s only works on standard images in (v)xyzt format.", infile );

  /* Read/Create missing image indicators */
  missing = get_missing( Input );

  /* Set parameters in local variables */
  if( !mri_has( Input, "images.extent.t" ) ||
      !mri_has( Input, "images.extent.x" ) ||
      !mri_has( Input, "images.extent.y" ) ||
      !mri_has( Input, "images.extent.z" ) )
    Abort( "images.extent key(s) missing from header." );
  dt = mri_get_int( Input, "images.extent.t" );
  dx = mri_get_int( Input, "images.extent.x" );
  dy = mri_get_int( Input, "images.extent.y" );
  dz = mri_get_int( Input, "images.extent.z" );
  if( ( dt <= 0 ) || ( dx <= 0 ) || ( dy <= 0 ) || ( dz <= 0 ) )
    Abort( "%s: images.extent key(s) is non-positive.\n", argv[0] );
  if (dt<2) Abort("%s: not enough data to smooth.\n", argv[0]);

  /* Allocate storage and initialize */
  if (!(par=(RegPar3D*)malloc(dt*sizeof(RegPar3D))))
    Abort("%s: unable to allocate %d bytes!\n",argv[0],dt*sizeof(RegPar3D));
  if (!(corr_par=(RegPar3D*)malloc(dt*sizeof(RegPar3D))))
    Abort("%s: unable to allocate %d bytes!\n",argv[0],dt*sizeof(RegPar3D));
  if (!(ordmse=(float*)malloc(dt*sizeof(float))))
    Abort("%s: unable to allocate %d bytes!\n",argv[0],dt*sizeof(float));
  image_missing= Matrix( dt, 1, unsigned char );

  for( t = 0; t < dt; t++ ) {
    par[t].q.x= par[t].q.y= par[t].q.z= 0.0; 
    par[t].q.w= 1.0;
    par[t].x= par[t].y= par[t].z= 0.0;
    par[t].mse= -1.0;
    corr_par[t]= par[t];
  }

  if (!(smooth_ibuf= (float*)malloc(8*dt*sizeof(float))))
    Abort("%s: unable to allocate %d floats!\n",argv[0],8*dt);
  for (i=0; i<8; i++) smooth_itbl[i]= smooth_ibuf+(i*dt);

  if (!(smooth_obuf= (float*)malloc(8*dt*sizeof(float))))
    Abort("%s: unable to allocate %d floats!\n",argv[0],8*dt);
  for (i=0; i<8; i++) smooth_otbl[i]= smooth_obuf+(i*dt);

  /* Read in estimated registration parameters */
  loadRegParams(iparfile, par, dt, &parNameString);

  /* Create smoother with specified params */
  smoother= sm_create_smoother();

  /* Build imagewise missing data */
  for (t=0; t<dt; t++) {
    this_image_missing= 1;
    for (z=0; z<dz; z++) {
      if (missing[t][z]==0) {
	this_image_missing= 0;
	break;
      }
    }
    image_missing[t][0]= this_image_missing;
  }

  /* Find fixed image */
  temp_fixed_image = fixed_image;
  if( ( fixed_image >= 0 ) && ( fixed_image < dt ) ) {
    while (1) {
      if (image_missing[temp_fixed_image][0]) {
	temp_fixed_image++;
	if( temp_fixed_image == dt ) temp_fixed_image = 0; /* wrap around */
	if( temp_fixed_image == fixed_image ) { /* give up */
	  temp_fixed_image = -1;
	  break;
	}
      }
      else break;
    }
  }

  /* Do a pass of smoothing including all fields, to
   * produce smoothed MSE values that will be used in
   * deciding which slices are marked missing.  We patch
   * over the zero value of mse at the fixed image to
   * avoid disturbing the smoothed mse curve.
   */
  for (t=0; t<dt; t++) {
    smooth_itbl[0][t]= par[t].q.x;
    smooth_itbl[1][t]= par[t].q.y;
    smooth_itbl[2][t]= par[t].q.z;
    smooth_itbl[3][t]= par[t].q.w;
    smooth_itbl[4][t]= par[t].x;
    smooth_itbl[5][t]= par[t].y;
    smooth_itbl[6][t]= par[t].z;
    smooth_itbl[7][t]= par[t].mse;
  }
  if (temp_fixed_image>=0) { /* patch over that image */
    if (temp_fixed_image>0) {
      if (temp_fixed_image<dt-1) {
	smooth_itbl[7][temp_fixed_image]= 
	  0.5*(par[temp_fixed_image-1].mse + par[temp_fixed_image+1].mse);
      }
      else smooth_itbl[7][temp_fixed_image]= par[temp_fixed_image-1].mse;
    }
    else smooth_itbl[7][temp_fixed_image]= par[1].mse;
  }

  SM_SMOOTH_GROUP( smoother, smooth_itbl, smooth_otbl, 8, dt, 
		   image_missing, 0 );
#ifdef never
  SM_SMOOTH_GROUP( smoother, smooth_itbl, smooth_otbl, 8, dt, 
		   NULL, 0 );
#endif

  /* Use this smoothing result as the smoothed mse */
  for (t=0; t<dt; t++) corr_par[t].mse= smooth_otbl[7][t];
  if (temp_fixed_image>=0)
    corr_par[temp_fixed_image].mse= par[temp_fixed_image].mse;

  /* MAKE "MISSING" EVALUATIONS BASED ON MSE'S */
          
  /* Sort MSE's ignoring fixed images and missing values */
  msecount = 0;
  for( t = 0; t < dt; t++ )
    if ( (t != temp_fixed_image) && (par[t].mse>=0) )
      ordmse[msecount++] = par[t].mse-corr_par[t].mse;
  fqsrt( ordmse, 0, (long) ( msecount - 1 ) );

  /* Calculate median and iqr for MSE's */
  medianmse = ordmse[(long) ( msecount / 2 )];
  iqrmse = ordmse[(long) ( 3 * msecount / 4 )] -
    ordmse[(long) ( msecount / 4 )];
  
  fprintf(stdout,"Median delta MSE= %g; IQR= %g\n",medianmse,iqrmse);
  
  /* Set maximum allowable value for uncontaminated MSE */
  msecutoff = medianmse + cutoff * iqrmse;

  /* If MSE is too far above smoothed MSE for an image, or is
   * zero or negative, set the image to missing.
   */
  for( t = 0; t < dt; t++ ) {
    if( (t != temp_fixed_image) 
	&& (( par[t].mse <= 0.0 ) 
	    || ( par[t].mse - corr_par[t].mse > msecutoff )) ) {
      for (z=0; z<dz; z++) missing[t][z] = (unsigned char) 1;
      image_missing[t][0] = (unsigned char) 1;
    }
  }
    
  /* END "MISSING" EVALUATIONS */     
  
  /* Smooth registration parameters again using new missing info,
   * without smoothing MSE.  Data for fixed image is unchanged.
   */

  SM_SMOOTH_GROUP( smoother, smooth_itbl, smooth_otbl, 7, dt, 
		   image_missing, 0 );
#ifdef never
  SM_SMOOTH_GROUP( smoother, smooth_itbl, smooth_otbl, 7, dt, 
		   NULL, 0 );
#endif

  for (t=0; t<dt; t++) {
    corr_par[t].q.x= smooth_otbl[0][t];
    corr_par[t].q.y= smooth_otbl[1][t];
    corr_par[t].q.z= smooth_otbl[2][t];
    corr_par[t].q.w= smooth_otbl[3][t];
    corr_par[t].x= smooth_otbl[4][t];
    corr_par[t].y= smooth_otbl[5][t];
    corr_par[t].z= smooth_otbl[6][t];
  }
  if (temp_fixed_image>=0) {
    corr_par[temp_fixed_image].q.x= smooth_otbl[0][temp_fixed_image];
    corr_par[temp_fixed_image].q.y= smooth_otbl[1][temp_fixed_image];
    corr_par[temp_fixed_image].q.z= smooth_otbl[2][temp_fixed_image];
    corr_par[temp_fixed_image].q.w= smooth_otbl[3][temp_fixed_image];
    corr_par[temp_fixed_image].x= smooth_otbl[4][temp_fixed_image];
    corr_par[temp_fixed_image].y= smooth_otbl[5][temp_fixed_image];
    corr_par[temp_fixed_image].z= smooth_otbl[6][temp_fixed_image];
  }
    
  /* END SMOOTHING */

  /* Normalize quaternions */
  for (t=0; t<dt; t++) {
    Quat* q= &(corr_par[t].q);
    if (q->x==0.0 && q->y==0.0 && q->z==0.0) q->w= 1.0;
    quat_normalize(q);
  }

  /* Write out missing chunk, since it may have changed */
  mri_set_chunk( Input, "missing", (long) ( dt * dz ), 0,
		 MRI_UNSIGNED_CHAR, *missing );  

  /* Close data set */
  mri_close_dataset( Input );
      
  /* Write out new parameters */
  ofp= initParFile( oparfile, infile, iparfile, parNameString );
  free(parNameString);
  for (t=0; t<dt; t++) 
    saveRegParams( ofp, corr_par+t, t );
  if (fclose(ofp)) {
    perror("Error closing output:");
  }

  /* Clean up a bit */
  sm_destroy(smoother);
  free( par );
  free( corr_par );
  free( ordmse );
  free( smooth_ibuf );
  free( smooth_obuf );
  FreeMatrix( image_missing );

  Message( "#      3D Registration parameter smoothing complete.\n" );

  return 0;
}



