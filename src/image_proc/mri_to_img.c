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
/* This module was derived from imtops by Joel Welling in 10/00. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <mri.h>
#include <fmri.h>
#include <stdcrg.h>
#include <misc.h>

#ifdef USE_PNG
#ifndef LINUX
#include <setjmp.h>
#endif
#include <png.h>
#endif

#ifdef CRAY
#include <fp.h>
#define finite( foo ) isfinite( foo )
#endif

#ifdef DARWIN
#define finite( foo ) isfinite( foo )
#endif

static char rcsid[] = "$Id: mri_to_img.c,v 1.18 2007/03/21 23:51:36 welling Exp $";

/* Notes...
 */

#define STR_SIZE 256 /* maximum size of strings */

#define DEFAULT_FORMAT IMG_PS
#define GRAY 1
#define COLOR 3
#define OPACITY 4
#define RED 0
#define GREEN 1
#define BLUE 2
#define JENN 1

typedef enum { IMG_PS, IMG_PNG } ImgFormat;

typedef struct opt_struct {
  ImgFormat format;
  void* hook; /* for things needed by the image writing routines */
  int b_flg, c_flg, d_flg, f_flg, g_flg;
  int h_flg, l_flg, n_flg, o_flg, p_flg;
  int r_flg, t_flg, w_flg, x_flg;
  int all_times_flg;
  int mosaic_flg;
  int columns;
  float window;
  float level;
  float gamma;
  float min;
  float max;
  float hist_scale;
  int subsa;
  int xorig;
  int yorig;
  int imht;
  int imwth;
} Options;

typedef struct img_info_struct {
  float* data;
  int vdim;
  int xdim;
  int ydim;
  float min;
  float max;
  int range_valid;
} ImageInfo;

static char* progname;

static void safe_copy(char* str1, const char* str2) {
  strncpy(str1, str2, STR_SIZE);
  str1[STR_SIZE-1]= '\0';
}

static void safe_concat(char* str1, const char* str2) {
  strncat(str1, str2, (STR_SIZE-strlen(str1))-1);
}

static int safe_get_extent(MRI_Dataset* ds, const char* chunk, char* dim)
{
  char key_buf[STR_SIZE];
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

static void calc_sizes(MRI_Dataset* ds, const char* this_chunk, 
		       const char* dimstr, char* selected_dim,
		       int* fast_blocksize_out, int* slow_blocksize_out ) {
  /* This routine will fail if selected_dim is not in dimstr! */
  int fast_blocksize;
  int slow_blocksize;
  char* dimstr_copy;
  char* this_dim;

  dimstr_copy= strdup(dimstr);

  fast_blocksize= 1;

  this_dim= dimstr_copy;
  while (*this_dim != *selected_dim) 
    fast_blocksize *= safe_get_extent(ds,this_chunk,this_dim++);

  this_dim++; /* step over selected dim */

  slow_blocksize= 1;
  while (*this_dim)
    slow_blocksize *= safe_get_extent(ds, this_chunk, this_dim++);

  *fast_blocksize_out= fast_blocksize;
  *slow_blocksize_out= slow_blocksize;

  free(dimstr_copy);
}

static int dataset_valid(MRI_Dataset* ds, const char* chunk, Options* opts)
{
  char keybuf[STR_SIZE];
  char* dimstr;
  int zdim;
  int xdim;
  int vdim;
  char* x_here;
  char* y_here;
  char* z_here;
  char* v_here;
  int fastblock;
  int slowblock;

  sprintf(keybuf,"%s.dimensions",chunk); /* chunk name size checked earlier */
  if (!mri_has(ds,keybuf)) {
    Error("%s: required tag <%s> missing!\n",progname,keybuf);
    return 0;
  }

  dimstr= mri_get_string(ds,keybuf);

  if (!(x_here=strchr(dimstr, 'x'))) {
    Error("%s: required x dimension is not present!\n",progname);
    return 0;
  }

  if (!(y_here=strchr(dimstr, 'y'))) {
    Error("%s: required y dimension is not present!\n",progname);
    return 0;
  }

  if (!(z_here=strchr(dimstr, 'z'))) {
    zdim= 1;
  }
  else zdim= safe_get_extent(ds,chunk,"z");

  if (!(v_here=strchr(dimstr, 'v'))) {
    vdim= 1;
  }
  else vdim= safe_get_extent(ds,chunk,"v");

  if (v_here && (x_here != (v_here+1))) {
    Error("%s: dimension v must be before dimension x!\n",progname);
    return 0;
  }

  if (y_here != (x_here+1)) {
    Error("%s: dimensions x and y are not sequential!\n",progname);
    return 0;
  }

  if (z_here && (z_here != (y_here+1))) {
    Error("%s: dimensions x, y, z are not sequential!\n",progname);
    return 0;
  }

  xdim= safe_get_extent(ds,chunk,"x");

  calc_sizes(ds, chunk, dimstr, "y", &fastblock, &slowblock);

  if (fastblock != vdim * xdim) {
    Error("%s: Need %d value(s) per pixel, found %d!\n",
	  progname, vdim, fastblock/xdim);
    return 0;
  }

  if (!opts->all_times_flg && (slowblock != zdim)) {
    Error("%s: too many images!\n",progname);
    return 0;
  }

  return 1;
}

static ImageInfo* load_slice(MRI_Dataset* ds, const char* chunk,
			     const int vdim, const int xdim, const int ydim, 
			     const int zdim, const int tdim, const int z, 
			     const int t)
{
  ImageInfo* result;
  int imgSize;
  int i;
  float min;
  float max;
  int bounds_set;
  float* inbuf;

  if (!(result= (ImageInfo*)malloc(sizeof(ImageInfo))))
    Abort("%s: unable to allocate %d bytes!\n",progname,sizeof(ImageInfo));

  imgSize= xdim*ydim*vdim;
  if (!(result->data= (float*)malloc(imgSize*sizeof(float))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,imgSize*sizeof(float));

  result->xdim= xdim;
  result->ydim= ydim;
  result->vdim= vdim;

  inbuf= mri_get_chunk(ds, chunk, imgSize, (t*zdim + z)*imgSize, MRI_FLOAT);

  bounds_set= 0;
  for (i=0; i<imgSize; i++) {
    result->data[i]= inbuf[i];
    if (finite(inbuf[i])) {
      if (bounds_set) {
	if (inbuf[i]<min) min= inbuf[i];
	if (inbuf[i]>max) max= inbuf[i];
      }
      else {
	min= max= inbuf[i];
	bounds_set= 1;
      }
    }
  }

  if (bounds_set) {
    result->min= min;
    result->max= max;
    result->range_valid= 1;
  }
  else {
    result->min= 0.0;
    result->max= 0.0;
    result->range_valid= 0;
  }

  return result;
}

static void free_image(ImageInfo* img)
{
  free(img->data);
  free(img);
}

static void adjust_bounds(ImageInfo* img, Options* opts)
{
  float min, max;
  int min_set, max_set;
  int tot;
  int i;

  tot= img->xdim * img->ydim;

  if (img->range_valid) {
    min= img->min;
    max= img->max;
    min_set= max_set= 1;
  }
  else {
    min= 0.0;
    max= 0.0;
    min_set= max_set= 0;
  }

  if (opts->w_flg) {
    min = opts->level - opts->window/2; 
    max = opts->level + opts->window/2; 
    min_set= max_set= 1;
  }
  else {
    if (opts->n_flg) {
      min= opts->min;
      min_set= 1;
    }
    if (opts->x_flg) {
      max= opts->max;
      max_set= 1;
    }
  }

  if (opts->b_flg)
    { max = min + (max - min)/2; }
  if (opts->d_flg)
    { min = min + (max - min)/2; }
  
  if (min_set && max_set && min!=max) {
    img->range_valid= 1;
    img->min= min;
    img->max= max;
  }
  else {
    Error("%s: can't find range! (min = %f, max= %f)\n",
	  progname,min,max);
    img->range_valid= 0;
  }
}

static ImageInfo* load_mosaic(MRI_Dataset* ds, const char* chunk,
			      const int vdim, const int xdim, 
			      const int ydim, const int zdim, 
			      const int tdim, const int t, const
Options* opts)
{
  ImageInfo* result;
  ImageInfo* small;
  int imgSize;
  int i;
  int j;
  int v;
  int ismall;
  int jsmall;
  int iblock;
  int jblock;
  int xdimBig;
  int ydimBig;
  int z;
  int columns, rows;
  float min;
  float max;
  int range_valid;
  float* inbuf;
  float fzdim, fcolumns, frows;

  /* Find the right decomposition */

  if (!opts->columns){
    columns= 1;
    while (columns*columns < zdim) { columns++;}
    rows= columns;
  }
  else {
    columns= opts->columns;
    rows = ceil((float)zdim/(float)columns);
  }
 
  if (!(result= (ImageInfo*)malloc(sizeof(ImageInfo))))
    Abort("%s: unable to allocate %d bytes!\n",progname,sizeof(ImageInfo));

  imgSize= columns*rows*xdim*ydim*vdim;

  if (!(result->data= (float*)malloc(imgSize*sizeof(float))))
    Abort("%s: unable to allocate %d bytes!\n",
	  progname,imgSize*sizeof(float));
  for (i=0; i<imgSize; i++) result->data[i]= 0.0;

  result->xdim= xdimBig= columns*xdim;
  result->ydim= ydimBig= rows*ydim;
  result->vdim= vdim;

  iblock= jblock= 0;
  range_valid= 0;
  min= max= 0.0;
  for (z=0; z<zdim; z++) {

    small= load_slice(ds, chunk, vdim, xdim, ydim, zdim, tdim, z, t);

    /* Merge in the new subimage's range */
    if (small->range_valid) {
      if (range_valid) {
	if (small->min < min) min= small->min;
	if (small->max > max) max= small->max;
      }
      else {
	range_valid= 1;
	min= small->min;
	max= small->max;
      }
    }

    /* Copy the new subimage into place */
    for (jsmall=0; jsmall<ydim; jsmall++) {
      j= jblock*ydim + jsmall;
      for (ismall=0; ismall<xdim; ismall++) {
	i= iblock*xdim + ismall;
	for (v=0; v<vdim; v++) {
	  result->data[(((j*xdimBig)+i)*vdim)+v]= small->data[(((jsmall*xdim) + ismall)*vdim)+v];
	}
      }
    }

      free_image(small);
      iblock++;
      if (iblock>=columns) {
	iblock= 0;
	jblock++;
      }
  }


  result->min= min;
  result->max= max;
  result->range_valid= range_valid;

  return result;
}

static int get_image_byte( ImageInfo* img, const Options* opts, 
			   const int i, const int j )
{
  float datum;
  int c;

  if (!(img->range_valid)) return 0;

  datum= img->data[((j*opts->subsa)*img->xdim) + (i*opts->subsa)];
  
  /* a way to keep nans out of histograms */
  if (!finite(datum))
    c= -1;
  else {
    float bounded_datum= 
      (datum>img->max) ? img->max:((datum<img->min) ? img->min : datum);
    c = (int)(((bounded_datum - img->min)/(img->max-img->min))  * 255.0) ;
  }
  if ((c<0) && (c!=-1)) c= 0;
  if (c>255) c= 255;

  if ((opts->r_flg) && (c!=-1)) c= 255-c;

  return c;
}

static int get_color_image_byte( ImageInfo* img, const Options* opts, 
				 const int i, const int j, const int color)
{
  float datum;
  int c;

  if (!(img->range_valid)) return 0;

  datum= img->data[(((((j*opts->subsa)*img->xdim) + (i*opts->subsa))*img->vdim) + color)];
  if (!finite(datum)) c= (opts->r_flg ? 255 : 0);
  else {
    c = (int)(((datum - img->min)/(img->max-img->min))  * 255.0) ;
  }
  if (c<0) c= 0;
  if (c>255) c= 255;

  if (opts->r_flg) c= 255-c;

  
#ifdef never
    fprintf(stderr, "min = %4.1f, max = %4.1f \n", img->min, img->max);
    fprintf(stderr, "datum = %4.1f, i = %d, j = % d, color = %d, c = %d \n", datum, i, j, color, c);
#endif   

  return c;
}

static FILE* init_out( char* fname, Options* opts, char* command_line )
{
  FILE* fp= NULL;

  switch (opts->format) {
  case IMG_PS:
    {
      char* hexcv;
      int i;

      if (fname==NULL || strlen(fname)==0) fp= stdout;
      else if (!(fp= fopen(fname,"w")))
	Abort("%s: unable to open %s for writing!\n",progname,fname);
      
      /* start pumping out the post script */
      fprintf(fp,"%%!PS-Adobe-1.0\n");
      fprintf(fp,"%%%%BoundingBox: 0 0 612 792\n");
      fprintf(fp,"%%%%Creator: %s\n", command_line );
      fprintf(fp,"%%%%Title: %s\n",fname);
      fprintf(fp,"/BodyF /Helvetica findfont 15 scalefont def \n" );
      fprintf(fp,"/StartPage { /SavedPage save def \n" );
      fprintf(fp,"  BodyF setfont 0 setgray } def \n" );
      fprintf(fp,"/EndPage   {showpage SavedPage restore } def \n" );
      fprintf(fp,"/outsidecircletext \n");
      fprintf(fp," { circtextdict begin \n");
      fprintf(fp,"  /radius exch def /centerangle exch def \n");
      fprintf(fp,"  /ptsize exch def /str exch def \n");
      fprintf(fp,"  /xradius radius ptsize 4 div add def \n");
      fprintf(fp,"  gsave \n");
      fprintf(fp,"   centerangle str findhalfangle add rotate \n");
      fprintf(fp,"   str \n");
      fprintf(fp,"    {/charcode exch def \n");
      fprintf(fp,"     (A) dup 0 charcode put outsideplacechar \n");
      fprintf(fp,"   } forall \n");
      fprintf(fp,"   grestore \n");
      fprintf(fp,"   end \n");
      fprintf(fp,"   } def \n");
      fprintf(fp,"/circtextdict 16 dict def \n");
      fprintf(fp,"circtextdict begin \n");
      fprintf(fp," /findhalfangle \n");
      fprintf(fp,"   {stringwidth pop 2 div \n");
      fprintf(fp,"    2 xradius mul pi mul div 360 mul \n");
      fprintf(fp,"   } def \n");
      fprintf(fp,"/outsideplacechar \n");
      fprintf(fp,"  {/char exch def \n");
      fprintf(fp,"   /halfangle char findhalfangle def \n");
      fprintf(fp,"   gsave \n");
      fprintf(fp,"    halfangle neg rotate \n");
      fprintf(fp,"    radius 0 translate \n");
      fprintf(fp,"    -90 rotate \n");
      fprintf(fp,"    char stringwidth pop 2 div neg 0 moveto \n");
      fprintf(fp,"    char show \n");
      fprintf(fp,"   grestore \n");
      fprintf(fp,"   halfangle 2 mul neg rotate \n");
      fprintf(fp,"  } def \n");
      fprintf(fp," /pi 3.1415923 def \n");
      fprintf(fp,"end \n");
      fprintf(fp,"%%%%EndProlog\n");

      if (!(hexcv= (char*)malloc(3*256))) 
	Abort("%s: unable to allocate %d bytes!\n",3*256);
      for (i=0; i<16; i++)   sprintf(hexcv+3*i,"0%1x",i);
      for (i=16; i<256; i++) sprintf(hexcv+3*i,"%2x",i);
      
      opts->hook= (void*)hexcv;
    }
    break;

  case IMG_PNG:
    {
      if (fname==NULL || strlen(fname)==0) fp= stdout;
      else if (!(fp= fopen(fname,"w")))
	Abort("%s: unable to open %s for writing!\n",progname,fname);
    }
  }

  return fp;
}

static void close_out( FILE* fp, Options* opts )
{
  switch (opts->format) {
  case IMG_PS:
    {
      if (fclose(fp)) perror("Error closing output");
      free( opts->hook );
    }
    break;
  case IMG_PNG:
    {
      if (fclose(fp)) perror("Error closing output");
    }
    break;
  }
}

static void drawPSLogo( FILE* fp, double xloc, double yloc, double scale )
{
  char* evar= NULL;
  fprintf(fp,"gsave \n");
  fprintf(fp,"%g %g translate %g %g scale \n",xloc,yloc,scale,scale);

  fprintf(fp,"0 0 moveto \n");
  fprintf(fp,"/Helvetica findfont 15 scalefont setfont \n");
  evar= getenv("FIASCO_VERSION");
  fprintf(fp,"(%s) show \n", (evar != NULL) ? evar : "???");
  fprintf(fp,"0 -2 moveto \n");
  
  fprintf(fp,"gsave \n");
  fprintf(fp,"/Helvetica findfont 20 scalefont setfont \n");
  fprintf(fp,"11 -2 translate \n");
  fprintf(fp,"(FIASCO) \n");
  fprintf(fp,"20 90 20 outsidecircletext \n");
  
  fprintf(fp,"(FIAT) \n");
  fprintf(fp,"20 270 20 outsidecircletext \n");
  fprintf(fp,"grestore \n");
  
  fprintf(fp,"-6 -2.5 moveto \n");
  fprintf(fp,"27 -2.5 lineto stroke \n");
  
  fprintf(fp,"gsave \n");
  fprintf(fp,"27 -5 moveto \n");
  fprintf(fp,"180 rotate \n");
  fprintf(fp,"(%s) show \n",FIAT_VERSION_STRING);
  fprintf(fp,"grestore \n");

  fprintf(fp,"grestore \n");
}

static void do_page( const char* creDate, const char* prtDate,
		     const char* pageLabel1, 
		     const char* pageLabel2, const char* pageLabel3,
		     ImageInfo* img, FILE* fp, Options* opts )
{
  static int pageNum= 0;
  int hist[255];
  int i;

  if (opts->h_flg)
    for (i=0; i<256; i++) hist[i]=0;

  pageNum++;

  switch (opts->format) {
  case IMG_PNG:
    {
#ifdef USE_PNG   
      png_structp png_ptr;
      png_infop info_ptr;
      png_text *notes;
      png_bytep row_ptr;
      int i;
      int j;
      int row_out = img->ydim/opts->subsa;
      int col_out = img->xdim/opts->subsa;
      int bar_width = 0;
      char minBuf[64];
      char maxBuf[64];
      int bytesPerPixel= 0;
      
      if (pageNum>1) {
	/* checks elsewhere should prevent this from happening */
	Abort("%s: Output for png must be single images!\n",progname);
      }

      png_ptr= png_create_write_struct
	(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
      if (!png_ptr)
	Abort("%s: unable to create png data structure! (1)\n",progname);
      
      info_ptr = png_create_info_struct(png_ptr);
      if (!info_ptr) {
	png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
	Abort("%s: unable to create png data structure! (2)\n",progname);
      }

      if (setjmp(png_ptr->jmpbuf)) {
	png_destroy_write_struct(&png_ptr, &info_ptr);
	Abort("%s: fatal error in png libraries!\n",progname);
      }

      png_init_io(png_ptr, fp);

      if (opts->g_flg) bar_width= col_out/16;
      else bar_width= 0;

      info_ptr->width= col_out + bar_width;
      info_ptr->height= row_out;

      switch (img->vdim) {
      case GRAY:
	info_ptr->color_type= PNG_COLOR_TYPE_GRAY;
	bytesPerPixel= 1;
	break;
      case OPACITY:
      case COLOR:
	info_ptr->color_type= PNG_COLOR_TYPE_RGB;
	bytesPerPixel= 3;
      }

      info_ptr->bit_depth= 8;

      info_ptr->gamma= opts->gamma;

      if (!(notes= (png_text*)malloc(10*sizeof(png_text))))
	Abort("%s: unable to allocate %d bytes!\n",10*sizeof(png_text));

      notes[0].key= "Title";
      notes[0].text= (png_charp)pageLabel1;
      notes[0].text_length= strlen(pageLabel1);
      notes[0].compression= -1;
      notes[1].key= "Label2";
      notes[1].text= (png_charp)pageLabel2;
      notes[1].text_length= strlen(pageLabel2);
      notes[1].compression= -1;
      notes[2].key= "Label3";
      notes[2].text= (png_charp)pageLabel3;
      notes[2].text_length= strlen(pageLabel3);
      notes[2].compression= -1;
      notes[3].key= "Software";
      notes[3].text= (png_charp)"mri_to_img";
      notes[3].text_length= strlen(notes[3].text);
      notes[3].compression= -1;
      notes[4].key= "Creation Time";
      notes[4].text= (png_charp)creDate;
      notes[4].text_length= strlen(creDate);
      notes[4].compression= -1;
      notes[5].key= "Translation Time";
      notes[5].text= (png_charp)prtDate;
      notes[5].text_length= strlen(prtDate);
      notes[5].compression= -1;
      sprintf(minBuf,"%g",img->min);
      notes[6].key= "Minimum";
      notes[6].text= minBuf;
      notes[6].text_length= strlen(minBuf);
      notes[6].compression= -1;
      sprintf(maxBuf,"%g",img->max);
      notes[7].key= "Maximum";
      notes[7].text= maxBuf;
      notes[7].text_length= strlen(maxBuf);
      notes[7].compression= -1;

      info_ptr->num_text= 8;
      info_ptr->max_text= 10;
      info_ptr->text= notes;

      png_write_info(png_ptr, info_ptr);

      if (!(row_ptr=(png_bytep)malloc(bytesPerPixel*(col_out+bar_width))))
	Abort("%s: unable to allocate %d bytes!\n",progname,
	      bytesPerPixel*(col_out+bar_width));
      for (j=0; j<row_out; j++) {

	switch(img->vdim) {

	case GRAY:
	  if (bar_width != 0) {
	    int barVal= (j*255)/(row_out-1);
	    if (!(opts->r_flg)) barVal= 255-barVal;
	    for (i=0; i<bar_width; i++) row_ptr[i]= barVal;
	  }
	  for (i=0; i<col_out; i++)
	    row_ptr[bar_width + i]= get_image_byte(img, opts, i, j);
	  break;

	case OPACITY:
	case COLOR:
	  if (bar_width != 0) {
	    int barVal= (j*255)/(row_out-1);
	    if (!(opts->r_flg)) barVal= 255-barVal;
	    for (i=0; i<bar_width; i++) {
	      row_ptr[bytesPerPixel*i]= barVal;
	      row_ptr[bytesPerPixel*i+1]= barVal;
	      row_ptr[bytesPerPixel*i+2]= barVal;
	    }
	  }
	  for (i=0; i<col_out; i++) {
	    row_ptr[bar_width+((bytesPerPixel * i) + RED)] = 
		   get_color_image_byte(img, opts, i, j, RED);
	    row_ptr[bar_width+((bytesPerPixel * i) + GREEN)]= 
		   get_color_image_byte(img, opts, i, j, GREEN);
	    row_ptr[bar_width+((bytesPerPixel * i) + BLUE)]= 
		   get_color_image_byte(img, opts, i, j, BLUE);
	  }
	}

	png_write_row(png_ptr, row_ptr);

      }
      free(row_ptr);

      png_write_end(png_ptr, info_ptr);
      png_destroy_write_struct(&png_ptr, &info_ptr); /* seems to free notes */
#endif
    }
    break;

  case IMG_PS:
    {
      int row_out = img->ydim/opts->subsa;
      int col_out = img->xdim/opts->subsa;
      int xorig=opts->xorig;
      int yorig;
      int imht= opts->imht;
      int imwth= opts->imwth;
      int i,j;
      char* hexcv;
      char outstr[80];
      char* cptr;
      float scale;
      int num_nan = 0;

      hexcv= (char*)opts->hook;

      /*      opts->h_flg > 0 ? yorig=288 : yorig = opts->yorig;*/
      if (opts->h_flg) yorig=206;
      else yorig=opts->yorig;

      /* Page setup */

      fprintf( fp, "%%%%Page: %ld %ld \n", pageNum, pageNum );
      fprintf( fp, "(%ld) (mri_to_img) StartPage \n", pageNum);      


      if(!opts->f_flg) {
	drawPSLogo(fp, 50.0, 52.0, 1.0);
      }

      fprintf(fp,"294 65 moveto \n");
      fprintf(fp,"/Helvetica findfont 15 scalefont setfont \n");
      if (opts->c_flg) {
	fprintf(fp,"(Created: %s) show \n",(creDate != NULL) ? creDate : " ");
      }
      else fprintf(fp,"%%%%CREATIONMSG%%%%\n");
      
      fprintf(fp,"300 50 moveto \n");
      fprintf(fp,"/Helvetica findfont 15 scalefont setfont \n");
      if(opts->p_flg) {
	fprintf(fp,"(Printed: %s) show \n",(prtDate != NULL) ? prtDate : " ");
      }
      else fprintf(fp,"%%%%PRINTMSG%%%%\n");

      if (opts->t_flg)
	{
	  fprintf(fp,"gsave\n");
	  fprintf(fp,"/Helvetica findfont 25 scalefont setfont\n");
	  if (!img->range_valid) {
	    fprintf(fp,"150 635 moveto\n");
	    fprintf(fp,"(Image invalid; bad range!) show\n");
	  }
	  fprintf(fp,"150 610 moveto\n");
	  fprintf(fp,"(%s) show\n",(pageLabel1 != NULL) ? pageLabel1 : " ");
	  fprintf(fp,"150 585 moveto\n");
	  fprintf(fp,"(%s) show\n",(pageLabel2 != NULL) ? pageLabel2 : " ");
	  fprintf(fp,"150 560 moveto\n");
	  fprintf(fp,"(%s) show\n",(pageLabel3 != NULL) ? pageLabel3 : " ");
	  fprintf(fp,"grestore\n");
	}

      /* Output the page proper */

      if (opts->gamma != 1.0)
	{
	  fprintf(fp,"/transarray [ %% def\n");
	  for (i=0; i< 256; i++)
	    fprintf(fp," %d",
		    (int) floor(100.0*pow(1.0*i/255,1.0/opts->gamma) + 0.5));
	  fprintf(fp,"\n] def %% transarray\n");
	  fprintf(fp," { \n 255 mul cvi %% multiply gray x 256 & make int\n");
	  fprintf(fp," transarray exch get    %% look up gray val in array\n");
	  fprintf(fp," 100 div    %% return gray setting to 0-1 range\n");
	  fprintf(fp," } settransfer\n");
	}
      
      if (opts->g_flg)
	{
	  fprintf(fp,"gsave\n");
	  fprintf(fp,"/Helvetica findfont 12 scalefont setfont\n");
	  fprintf(fp,"%d %d moveto\n",xorig-110, yorig+imht-10);
	  fprintf(fp,"(%6.3g) show\n", img->max);
	  fprintf(fp,"%d %d moveto\n",xorig-110, (int)(yorig+(imht-10)*.75));
	  fprintf(fp,"(%6.3g) show\n", .75*img->max +.25*img->min);
	  fprintf(fp,"%d %d moveto\n",xorig-110, (int)(yorig+(imht-10)*.5));
	  fprintf(fp,"(%6.3g) show\n", .5*img->max+ .5*img->min);
	  fprintf(fp,"%d %d moveto\n",xorig-110, (int)(yorig+(imht-10)*.25));
	  fprintf(fp,"(%6.3g) show\n", .25*img->max+.75*img->min);
	  fprintf(fp,"%d %d moveto\n",xorig-110, (int)(yorig+5));
	  fprintf(fp,"(%6.3g) show\n", img->min);
	  fprintf(fp,"%d %d add %d translate\n",xorig-25,-imwth/16,yorig);
	  fprintf(fp,"%d %d scale\n",imwth/16,imht);
	  if (opts->r_flg){
	    fprintf(fp,"1 %d 8 [1 0 0 %d 0 %d]\n",imht,- imht,imht);
	  }else{
	    fprintf(fp,"1 %d 8 [1 0 0 %d 0 0]\n",imht,imht);
	  }
	  fprintf(fp,"{<");
	  for (i=0; i < imht; i++)
	    fprintf(fp,"%02x",(unsigned int)((256*i)/imht));
	  fprintf(fp,">}\n");
	  fprintf(fp,"image\n");
	  fprintf(fp,"grestore\n");
	}

      /* #ifdef JENN */
      fprintf(fp,"gsave\n");
      fprintf(fp,"\n/picstr %d string def\n",col_out);
      fprintf(fp,"%d %d translate\n %d %d scale\n",xorig,yorig,imwth,imht) ;
      fprintf(fp,"%d %d 8 [%d 0 0 %d 0 %d]\n", 
	      col_out, row_out, col_out, -row_out, row_out) ;
      /* #endif */

      /* #ifdef JENN */
      switch (img->vdim) {

      case GRAY:
	{
	  fprintf(fp,"{currentfile picstr readhexstring pop} image\n\n");
     
	  cptr= outstr;
	  for (j=0; j<row_out; j++) {
	    for (i=0; i<col_out; i++) {
	      /* coords get flipped to satisfy postscript orientation */
	      int val= get_image_byte( img, opts, i, j );
       	      /* fprintf(stderr, "%d\n", val); */
	      if ( val == -1){
		 val = (opts->r_flg ? 255 : 0);
		 num_nan++;
      	         /* fprintf(stderr, "%d\n", num_nan); */
	      }
	      else if (opts->h_flg) hist[val]++;
	 
	  /* #ifdef JENN */
	      *cptr++= hexcv[3*val];
	      *cptr++= hexcv[3*val + 1];
	      if (cptr-outstr >= sizeof(outstr)-2) {
		*cptr++= '\n';
		*cptr++= '\0';
		fputs(outstr, fp);
		cptr= outstr;
	      }
	    }
	  }
	  if (cptr>outstr) {
	    *cptr++= '\n';
	    *cptr++= '\0';
	    fputs(outstr, fp);
	    cptr= outstr;
	  }
	  /* #endif */


	  /* #ifdef JENN */
	  if (opts->h_flg) {

	    /* fprintf(stderr, "number of pixels = %d\n", row_out*col_out);*/
	    /* fprintf(stderr, "number of nans = %d\n", num_nan);*/
	    fprintf(fp,"\n");
	    fprintf(fp,"grestore\n");
	    fprintf(fp,"/Helvetica findfont 12 scalefont setfont\n");
	    fprintf(fp,"%d %d moveto\n",xorig-55, yorig-40);
	    fprintf(fp,"(%6.3g) show\n", img->min);
	    fprintf(fp,"%d %d moveto\n",xorig-55, yorig-155);
	    fprintf(fp,"(%6.3g) show\n", img->max);
	    fprintf(fp,"%d %d moveto\n", xorig, yorig-30);
	    /*	    if (!opts->mosaic_flg) scale=(float)1;
		    else */
	    /*	    scale1=14256; */
	    scale= (float)opts->hist_scale/(row_out*col_out-num_nan);
	    /*	    fprintf(fp,"%f %f scale\n", (float)28512/(row_out*col_out), (float)144/256);
	     */
	    fprintf(fp,"%f %f scale\n", scale, (float)144/256);

	    /* fprintf(stderr, "row out=%d, col out=%d\n", row_out,col_out);*/
	    /* fprintf(stderr, "scale=%f\n", scale);*/

	    fprintf(fp,"-90 rotate\n");
	    for (i=0  ; i<256; i++) {
	      /* fprintf(stderr, "hist[%d]=%d\n",i,hist[i]); */
	      if (hist[i]==0) fprintf(fp,"%d %d rlineto\n",imwth/256,0);
	      else {
		fprintf(fp,"%d %d rlineto\n",0,hist[i]);
		fprintf(fp,"%d %d rlineto\n",imwth/256,0);
		fprintf(fp,"%d %d rlineto\n",0,-hist[i]);
	      /* fprintf(stderr, "hist[%d]=%d\n", i, hist[i]); */
	      }
	    }
	    fprintf(fp,"stroke\n");
	  }
	  /* #endif */


	}
	break;

      case OPACITY:
      case COLOR:
	{
	  int red, green, blue;
	  fprintf(fp,"{currentfile picstr readhexstring pop} false 3 colorimage\n\n");
      
	  cptr= outstr;
	  for (j=0; j<row_out; j++) {
	    for (i=0; i<col_out; i++) {
	      /* coords get flipped to satisfy postscript orientation */
	      red= get_color_image_byte( img, opts, i, j, RED );
	      *cptr++= hexcv[3*red];
	      *cptr++= hexcv[3*red + 1];
	      if (cptr-outstr >= sizeof(outstr)-2) {
		*cptr++= '\n';
		*cptr++= '\0';
		fputs(outstr, fp);
		cptr= outstr;
	      }
	      green= get_color_image_byte( img, opts, i, j, GREEN );
	      *cptr++= hexcv[3*green];
	      *cptr++= hexcv[3*green + 1];
	      if (cptr-outstr >= sizeof(outstr)-2) {
		*cptr++= '\n';
		*cptr++= '\0';
		fputs(outstr, fp);
		cptr= outstr;
	      }
	      blue= get_color_image_byte( img, opts, i, j, BLUE );
	      *cptr++= hexcv[3*blue];
	      *cptr++= hexcv[3*blue + 1];
	      if (cptr-outstr >= sizeof(outstr)-2) {
		*cptr++= '\n';
		*cptr++= '\0';
		fputs(outstr, fp);
		cptr= outstr;
	      }
	    }
	  }
	  if (cptr>outstr) {
	    *cptr++= '\n';
	    *cptr++= '\0';
	    fputs(outstr, fp);
	    cptr= outstr;
	  }
	}
      }

      /* #endif */

      /* Finish the page */

      fprintf(fp,"\n") ;
      fprintf(fp,"grestore\n");
      fprintf(fp,"EndPage\n");
    
    }
  }
}

int main(argc, argv)
int argc;
char *argv[];
{
  MRI_Dataset* ds;
  FILE* ofile;
  Options opts;
  
  char infile[STR_SIZE];
  char chunk[STR_SIZE];
  char keybuf[STR_SIZE];
  char outfile[STR_SIZE];
  char command[STR_SIZE];
  char **command_v;

  int xdim, ydim, zdim, vdim, tdim;
  int z;
  int t;
  int i; 
  int junk;
  int command_c;
  char *pageLabel1;
  char *pageLabel2;
  char *pageLabel3;
  char *creDate;
  char *prtDate;
  char creDateBuf[STR_SIZE];
  char prtDateBuf[STR_SIZE];
  int format_set_count= 0;
  
  opts.format= DEFAULT_FORMAT;
  opts.hook= NULL;

  progname= argv[0];  

  /* Check to see if help was requested */
  if (testHelp(&argc, argv)) exit(0);

  command_c=argc;
  command_v=Malloc(argc, char *);
  for (i=0;i<command_c; i++){
    command_v[i] = strdup(argv[i]);
    sprintf(command, "%s%s%s", command, command_v[i]," ");
  }
 
  /*** Parse command line and get input params ***/
  cl_scan( argc, argv );


  
  if (cl_present("png")) {
    opts.format= IMG_PNG;
    format_set_count++;
  }
  if (cl_present("ps")) {
    opts.format= IMG_PS;
    format_set_count++;
  }

  /* Deprecate old options */ 

   if (cl_present( "b" ))
     Abort ("Option b(brighten) has been replaced by brighten|bri.  Please see help file.\n");
   if (cl_present( "c" ))
     Abort ("Option c(created_date) has been replaced by date_created|dtc.  Please see help file.\n");
   if (cl_present( "d" ))
     Abort ("Option d(darken) has been replaced by darken|drk.  Please see help file.\n");
   if (cl_present( "g" ))
     Abort ("Option g(gamma) has been replaced by gray_scale|scg.  Please see help file.\n");
   if (cl_present( "p" ))
     Abort ("Option p(printed date) has been replaced by date_printed|dtp.  Please see help file.\n");
   if (cl_present( "t" ))
     Abort ("Option t(title) has been replaced by title|ttl.  Please see help file.\n");
   if (cl_present( "w" ))
     Abort ("Option w(window width) has been replaced by window_width|win.  Please see help file.\n");
   if (cl_present( "l" ))
     Abort ("Option l(level) has been replaced by level|lev.  Please see help file.\n");
   if (cl_present( "x" ))
     Abort ("Option x(maximum level) has been replaced by black_max|bmx.  Please see help file.\n");
   if (cl_present( "n" ))
     Abort ("Option n(minimum level) has been replaced by black_min|bmn.  Please see help file.\n");
   if (cl_present( "subsa" ))
     Abort ("Option subsa(subsampling value) has been replaced by subsample|sub.  Please see help file.\n");
   if (cl_present( "imht" ))
     Abort ("Option imht(image height) has been replaced by image_height|imh.  Please see help file.\n");
   if (cl_present( "imwt" ))
     Abort ("Option imwt(image weight) has been replaced by image_width|imw.  Please see help file.\n");
   if (cl_present( "o" ))
     Abort ("Option o(outfile) has been replaced by an input/output format.  Please see help file.\n");
   if (cl_present( "k" ))
     Abort ("Option k(color range broaden) no longer works.  Please see help file.\n");
  if (cl_present( "dcr" ))
     Abort ("Option dcr(date created) has been replaced by date_created|dtc.  Please see help file.\n");
  if (cl_present( "dpr" ))
     Abort ("Option dpr(date printed) has been replaced by date_printed|dtp.  Please see help file.\n");
  if (cl_present( "sgr" ))
     Abort ("Option sgr(gray scale) has been replaced by gray_scale|scg.  Please see help file.\n");
  if (cl_present( "tit" ))
     Abort ("Option tit(title) has been replaced by title|tit.  Please see help file.\n");
  if (cl_present( "wwt" ))
     Abort ("Option wwt(window width) has been replaced by window_width|win.  Please see help file.\n");


  opts.b_flg= cl_present("brighten|bri");
  opts.c_flg= cl_present("date_created|dtc");
  opts.d_flg= cl_present("darken|dar");
  opts.f_flg= cl_present("no_logo|nol");
  opts.g_flg= cl_present("gray_scale|scg");
  opts.h_flg= cl_present("histogram|hst");
  opts.n_flg= cl_get("black_min|bmn","%option %f", &(opts.min));
  opts.x_flg= cl_get("black_max|bmx","%option %f", &(opts.max));
  opts.p_flg= cl_present("date_printed|dtp");
  opts.r_flg= cl_present("reverse|rev|r");
  opts.t_flg= cl_present("title|ttl");
  opts.window= 1024.0;
  opts.w_flg= cl_get("window_width|win","%option %f", &(opts.window));
  opts.level= 1024.0;
  opts.l_flg= cl_get("level|lev","%option %f", &(opts.level));
  cl_get("columns|col","%option %d[0]", &(opts.columns));
  cl_get("gamma|gam","%option %f[1.0]", &(opts.gamma));
  cl_get("subsample|sub","%option %d[1]", &(opts.subsa));
  cl_get("xorig|xor","%option %d[144]", &(opts.xorig));
  cl_get("yorig|yor","%option %d[144]", &(opts.yorig));
  cl_get("image_height|imh","%option %d[400]", &(opts.imht));
  cl_get("image_width|imw","%option %d[400]", &(opts.imwth));
  cl_get("hist_scale|hsc", "%option %f[28512]", &(opts.hist_scale));
  opts.all_times_flg= cl_present("all_times|all");

  opts.mosaic_flg= (cl_present("mosaic|mos") || opts.columns>0);

  cl_get("chunk|chu|c","%option %s[%]","images",chunk);

  if(!cl_get("", "%s", infile)) {
    fprintf(stderr, "%s: Input file name not given.\n", argv[0]);
    exit(-1);
  }
  if(!cl_get("", "%s", outfile)) {
    fprintf(stderr, "%s: Output file name not given.\n", argv[0]);
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

  /* Too many formats? */
  if (format_set_count>1) {
    Abort("%s: multiple incompatible output formats specified.\n",argv[0]);
  }

  /* Diddle switches for specific formats */
  switch (opts.format) {
  case IMG_PS: break;
  case IMG_PNG: 
    {
#ifndef USE_PNG
      Abort("%s: this executable was compiled without PNG support!\n",argv[0]);
#endif
      if (!(opts.mosaic_flg==1))
	Warning(1, "%s: -png only supports mosaic, mosaic flag set to 1\n", argv[0]);
      opts.mosaic_flg= 1;
      if (opts.all_times_flg)
	Error("%s: -all flag incompatible with -png!\n",argv[0]);
      opts.all_times_flg= 0;
    }
    break;
  }

  /* Check filename for validity; add ".mri" if needed. */
  if (strlen(infile)==0) {
    fprintf(stderr,"%s: required filename not given\n",argv[0]);
    Help( "usage" );
    exit(-1);
  }
  if (strlen(infile)<4 || strcmp(infile+strlen(infile)-4,".mri"))
    safe_concat(infile,".mri");

  /* Check for wildly long chunk names, which can cause us to overrun
   * editing buffers.
   */
  if (strlen(chunk)+64 >= sizeof(keybuf)-1) 
    Abort("%s: chunk name too long!\n",argv[0]);

  ds= mri_open_dataset(infile,MRI_READ);
  if (!dataset_valid(ds,chunk,&opts))
    Abort("%s: Chunk <%s> in %s cannot be converted.\n",
	  argv[0],chunk,infile);
  xdim= safe_get_extent(ds,chunk,"x");
  ydim= safe_get_extent(ds,chunk,"y");
  sprintf(keybuf,"%s.extent.z",chunk);
  if (mri_has(ds,keybuf)) zdim= mri_get_int(ds,keybuf);
  else zdim= 1;
  sprintf(keybuf,"%s.extent.v",chunk);
  if (mri_has(ds,keybuf)) vdim= mri_get_int(ds,keybuf);
  else vdim = 1;
  if (!((vdim == 1) || (vdim == 3) || (vdim == 4)))
    Abort("%s: mri_to_image only works for v=1,3, and 4.\n", argv[0]);
  sprintf(keybuf,"%s.dimensions",chunk);
  calc_sizes(ds, chunk, 
	     mri_get_string( ds, keybuf ),
	     "y", &junk, &tdim);
  tdim /= zdim;

  /* Initialize output */
  ofile= init_out(outfile, &opts, command);

  /* Do the pages */
  pageLabel1= getenv("F_PS_TITLE1");
  if (!pageLabel1) pageLabel1= infile;
  creDate= getenv("F_CREDATE");
  if (!creDate) {
    struct stat st;
    stat(infile,&st);
    strncpy(creDateBuf, ctime(&(st.st_mtime)), STR_SIZE);
    creDate= creDateBuf;
  }
  prtDate= getenv("F_PRTDATE");
  if (!prtDate) {
    time_t t;
    t= time(NULL);
    strncpy(prtDateBuf, ctime(&t), STR_SIZE);
    prtDate= prtDateBuf;
  }
  
  for (t=0; t<tdim; t++) {
    char title2Buf[STR_SIZE];

    pageLabel2= getenv("F_PS_TITLE2");
    if (!pageLabel2) {
      sprintf(title2Buf,"Image %d\n",t);
      pageLabel2= title2Buf;
    }

    if (opts.mosaic_flg) {
      ImageInfo* img;
      char title3Buf[STR_SIZE];
      
	  pageLabel3= getenv("F_PS_TITLE3");
	  if (!pageLabel3) {
	    sprintf(title3Buf,"%d Slices\n",zdim);
	    pageLabel3= title3Buf;
	  }
	  img= load_mosaic(ds, chunk, vdim, xdim, ydim, zdim, tdim, t,
&opts);
	  adjust_bounds(img, &opts);
	  do_page(creDate, prtDate,
		  pageLabel1, pageLabel2, pageLabel3, img, ofile, &opts);
	  free_image(img);

    }

    else {
      for (z=0; z<zdim; z++) {
	ImageInfo* img;
	char title3Buf[STR_SIZE];
	pageLabel3= getenv("F_PS_TITLE3");
	if (!pageLabel3) {
	  sprintf(keybuf,"%s.scan.date",chunk);
	  if (mri_has(ds,keybuf)) {
	    safe_copy(title3Buf,"Acquired ");
	    safe_concat(title3Buf,mri_get_string(ds,keybuf));
	    sprintf(keybuf,"%s.scan.time",chunk);
	    if (mri_has(ds,keybuf)) {
	      safe_concat(title3Buf," ");
	      safe_concat(title3Buf,mri_get_string(ds,keybuf));
	    }
	    pageLabel3= title3Buf;
	  }
	  else {
	    sprintf(title3Buf,"Slice %d\n",z);
	    pageLabel3= title3Buf;
	  }
	}
	img= load_slice(ds, chunk, vdim, xdim, ydim, zdim, tdim, z, t);
	adjust_bounds(img, &opts);
	do_page(creDate, prtDate,
		pageLabel1, pageLabel2, pageLabel3, img, ofile, &opts);
	free_image(img);
      }
    }
  }

  /* Clean up */
  close_out(ofile, &opts);

  return 0;
}     /* End of main program */


