/************************************************************
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *     Copyright (c) 1995 Department of Statistics,         *
 *                        Carnegie Mellon University        *
# *                                                          *
# *  This program is distributed in the hope that it will    *
# *  be useful, but WITHOUT ANY WARRANTY; without even the   *
# *  implied warranty of MERCHANTABILITY or FITNESS FOR A    *
# *  PARTICULAR PURPOSE.  Neither Carnegie Mellon University *
# *  nor any of the authors assume any liability for         *
# *  damages, incidental or otherwise, caused by the         *
# *  installation or use of this software.                   *
# *                                                          *
# *  CLINICAL APPLICATIONS ARE NOT RECOMMENDED, AND THIS     *
# *  SOFTWARE HAS NOT BEEN EVALUATED BY THE UNITED STATES    *
# *  FDA FOR ANY CLINICAL USE.                               *
# *                                                          *
 *                                                          *
 *  Original programming by Audris Mockus                   *
 ************************************************************/

/* Heavy mods by Joel Welling (starting with Audris' read.c)
 * to support Pgh MRI format, 7/97.
 */

/*  $Id: read_mri.c,v 1.8 2007/03/21 23:44:49 welling Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include "mri.h"

static char rcsid[] = "$Id: read_mri.c,v 1.8 2007/03/21 23:44:49 welling Exp $";

#define NAME_LENGTH 1023
static char Name[NAME_LENGTH+1];
static char DataName[NAME_LENGTH+1];
static char valid= 0;

#define MAX_DIM 10
static int dim_array[MAX_DIM];
static char dim_string[MAX_DIM+1];

static int Type;
static int Ndim;
static int Veclen;
static struct stat StatBuf;
void* DataBuf= NULL;

int
readFieldHead (char ** name,
		 int * type,
		 int * ndim,
		 int * veclen,
		 int * dims)
{
  int i;
  char* typestr= NULL;
  char* dimstr= NULL;
  char* runner;
  MRI_Dataset *TheField;
  FILE* datafile;
  
  fprintf(stderr,"point 1; name is <%s>\n",*name);
  if (strlen(*name)>NAME_LENGTH) {
    fprintf(stderr,"readFieldHead: name too long!\n");
    return -1;
  }

  fprintf(stderr,"point 2\n");
  strcpy(Name,*name);
  fprintf(stderr,"point 2.5\n");
  TheField= mri_open_dataset(*name, MRI_READ);
  fprintf(stderr,"point 3\n");

  if (mri_has(TheField,"images.file")) {
    char* ext;
    char* data_fname;
    int length;
    
    data_fname= mri_get_string( TheField, "images.file" );
    length= strlen(data_fname);
    for (ext= data_fname + (length-1); ext>=data_fname; ext--)
      if (*ext == '.') break;
    if (ext==data_fname) {
      /* extension only is given */
      int basename_length;
fprintf(stderr,"point 3.1\n");
      strcpy(DataName,Name);
      basename_length= strlen(DataName);
      if (basename_length>4 && !strcmp(DataName+(basename_length-4),".mri"))
	*(DataName+(basename_length-4))= '\0';
      if (strlen(DataName) + length > NAME_LENGTH) {
	fprintf(stderr,"readFieldHead: generated data file name too long!\n");
	return -1;
      }
      strcat(DataName,data_fname);
    }
    else {
fprintf(stderr,"point 3.2\n");
      if (length>NAME_LENGTH) {
	fprintf(stderr,"readFieldHead: data file name too long!\n");
	return -1;
      }
      strcpy(DataName,data_fname);
    }
  }
  else {
    fprintf(stderr,
	    "readFieldHead: input data file must have separate data file.\n");
    return -1;
  }
  
  fprintf(stderr,"point 4: Name is <%s>, DataName is <%s>\n",Name,DataName);
  if( !mri_has( TheField, "images" ) ) {
    fprintf(stderr,"readFieldHead: input data is not standard images.\n");
    return -1;
  }
  
  if ( !(typestr= mri_get_string( TheField, "images.datatype" )) ) {
    fprintf(stderr,"readFieldHead: input data type info not present.\n");
    return -1;    
  }
  if (!strcmp(typestr,"uint8")) {
    *type= 0;
  } else if (!strcmp(typestr,"int16")) {
    *type= 4;
  } else if (!strcmp(typestr,"int32")) {
    *type= 1;
  } else if (!strcmp(typestr,"float32")) {
    *type= 2;
  } else if (!strcmp(typestr,"float64")) {
    *type= 3;
  } else {
    fprintf(stderr,"readFieldHead: unknown input data type %s.\n",typestr);
    return -1;    
  }
  Type= *type;
  fprintf(stderr,"readFieldHead: %s is type %d\n",*name,*type);

  if ( !(dimstr= mri_get_string( TheField, "images.dimensions" )) ) {
    fprintf(stderr,"readFieldHead: input data lacks dimension info.\n");
    return -1;    
  }
  fprintf(stderr,"readFieldHead: image dimensions <%s>\n",dimstr);
  if (strlen(dimstr)>MAX_DIM) {
    fprintf(stderr,"readFieldHead: input data dimensionality too high.\n");
    return -1;
  }
  strcpy(dim_string,dimstr);

  if (mri_has( TheField, "images.extent.v" )) {
    *ndim= strlen(dimstr)-1;
    *veclen= mri_get_int( TheField, "images.extent.v" );
  }
  else {
    *ndim= strlen(dimstr);
    *veclen= 1;
  }
  Ndim= *ndim;
  Veclen= *veclen;
  fprintf(stderr,"readFieldHead: dim string <%s>; ndim %d, veclen %d\n",
	  dimstr,*ndim,*veclen);

  /* walk dimstr, finding extent of each dimension */
  i= 0;
  for (runner= dimstr; *runner; runner++) {
    char attrstring[64];
    if (*runner != 'v') {
      sprintf(attrstring,"images.extent.%c",*runner);
      dim_array[i]= dims [i] = mri_get_int(TheField, attrstring);
      fprintf(stderr,"readFieldHead: dim %c, extent %d\n",*runner,dims[i]);
      i++;
    }
  }

fprintf(stderr,"readFieldHead: unmap\n");
  /* Unmap any leftovers */
  if (DataBuf) {
    munmap((void*)DataBuf,StatBuf.st_size);
    DataBuf= NULL;
  }

fprintf(stderr,"readFieldHead: open\n");
  /* Open the data file */
  if (!(datafile= fopen(DataName,"r"))) {
    perror("readFieldHead: cannot open input data");
    return -1;
  }

fprintf(stderr,"readFieldHead: stat\n");
  /* Get file characteristics */
  if (stat(DataName,&StatBuf)) {
    perror("readFieldHead: cannot stat data");
    return -1;
  }

fprintf(stderr,"readFieldHead: mmap\n");
  /* Memory map the thing */
#if ( SGI5 || SGI64 || SGIMP || SGI64MP || SGI5MP || SUN4SOL2 || DARWIN )
  DataBuf= (void*)mmap(0, StatBuf.st_size, 
		       PROT_READ, MAP_SHARED,
		       fileno(datafile), 0);
#else
  DataBuf= (void*)mmap(0, StatBuf.st_size, 
		       PROT_READ, MAP_FILE | MAP_SHARED,
		       fileno(datafile), 0);
#endif
  if (DataBuf==(void*)-1) {
    perror("readFieldHead: cannot mmap data");
    return -1;
  }

fprintf(stderr,"readFieldHead: clean up\n");
  valid= 1;
  mri_close_dataset(TheField);
  fclose(datafile);
fprintf(stderr,"returning from readFieldHead\n");

  return 0;
}

int
readField (int * dimsL, int * dimsU,
			  float * value)
{
  int sliceId;
  int timeId;
  int x_extent;
  int y_extent;
  int v_extent;
  int i;
  int j;
  int v;
  float* image;
  float* runner= value;
  MRI_Dataset *TheField;

fprintf(stderr,"mmap test follows\n");
fprintf(stderr,"%d %d %d\n",*(int*)DataBuf,*((int*)DataBuf+1),
	*((int*)DataBuf+2));

fprintf(stderr,"readField point 1\n");
  if (!Name){
    fprintf (stderr, "readField: No field is open\n");
    return -1;
  }

  /* At this point we have to reopen the field, because Splus won't
   * let us keep the field pointer (allocated with S_alloc) around
   * between calls.
   */
  TheField= mri_open_dataset(Name, MRI_READ);
  
fprintf(stderr,"readField point 2\n");
  if (mri_has(TheField,"images.extent.v")) {
fprintf(stderr,"readField point 3\n");
    v_extent= mri_get_int(TheField,"images.extent.v");
  }
  else v_extent= 1;
fprintf(stderr,"readField point 4\n");
  x_extent= mri_get_int(TheField,"images.extent.x");
  y_extent= mri_get_int(TheField,"images.extent.y");

  for (timeId= dimsL[3]; timeId<=dimsU[3]; timeId++) {
    for (sliceId= dimsL[2]; sliceId<=dimsL[2]; sliceId++) {
      image= mri_get_image(TheField,timeId,sliceId,MRI_FLOAT32);
      for (j= dimsL[1]; j<=dimsU[1]; j++) {
	for (i= dimsL[0]; i<=dimsU[0]; i++) {
	  for (v=0; v<v_extent; v++) {
	    *runner++ = (float)(image[ ((j*x_extent + i)*v_extent) + v ]);
	  }
	}
      }
    }
  }

  mri_close_dataset(TheField);
fprintf(stderr,"readField return\n");
  return 0;
}

int
writeField (char ** name,
		 int * type,
		 int * ndim,
		 int * veclen,
		 int * dims, float * val)
{
  MRI_Dataset *dset;
  char* dimstr= NULL;
  char* runner;
  int k, kmax;
  int t, tmax;

  dset= mri_open_dataset(*name,MRI_WRITE);

  switch (*type) {
  case 0:
    mri_set_string(dset,"images.datatype","uint8");
    break;
  case 1:
    mri_set_string(dset,"images.datatype","int32");
    break;
  case 2:
    mri_set_string(dset,"images.datatype","float32");
    break;
  case 3:
    mri_set_string(dset,"images.datatype","float64");
    break;
  case 4:
    mri_set_string(dset,"images.datatype","int16");
    break;
  default:
    fprintf(stderr,"writeField: unknown field type %d\n",*type);
    return -1;
  }

  if (*veclen>1) {
    kmax= *veclen;
    if (*ndim==2) {
      dimstr= "vxy";
      kmax= 1;
      tmax= 1;
    }
    else if (*ndim==3) {
      dimstr= "vxyz";
      kmax= dims[2];
      tmax= 1;
    }
    else if (*ndim==4) {
      dimstr= "vxyzt";
      kmax= dims[2];
      tmax= dims[3];
    }
  }
  else {
    if (*ndim==2) {
      dimstr= "xy";
      kmax= 1;
      tmax= 1;
    }
    else if (*ndim==3) {
      dimstr= "xyz";
      kmax= dims[2];
      tmax= 1;
    }
    else if (*ndim==4) {
      dimstr= "xyzt";
      kmax= dims[2];
      tmax= dims[3];
    }
  }
  if (!dimstr) {
    fprintf(stderr,"writeField: data dimensionality %d not supported.",
	    *ndim);
    return -1;
  }
fprintf(stderr,"writeField: ended up with dimstr <%s>\n",dimstr);

  mri_set_string(dset,"images.dimensions",dimstr);

  for (runner= dimstr; *runner; runner++) {
    switch (*runner) {
    case 'v': 
      mri_set_int(dset,"images.extent.v",*veclen);
      break;
    case 'x':
      mri_set_int(dset,"images.extent.x",dims[0]);
      break;
    case 'y':
      mri_set_int(dset,"images.extent.y",dims[1]);
      break;
    case 'z':
      mri_set_int(dset,"images.extent.z",dims[2]);
      break;
    case 't':
      mri_set_int(dset,"images.extent.t",dims[3]);
      break;
    }
  }

  switch (*type) {
  case 0:
    {
      char* frunner= (char*)val;
      for (t=0; t<tmax; t++)
	for (k=0; k<kmax; k++) {
	  mri_set_image(dset, t, k, MRI_UINT8, frunner);
	  frunner += *veclen * dims[0] * dims[1];
	}
    }
  case 1:
    {
      long* frunner= (long*) val;
      for (t=0; t<tmax; t++)
	for (k=0; k<kmax; k++) {
	  mri_set_image(dset, t, k, MRI_INT32, frunner);
	  frunner += *veclen * dims[0] * dims[1];
	}
    }
  case 2:
    {
      float* frunner= val;
      for (t=0; t<tmax; t++)
	for (k=0; k<kmax; k++) {
	  mri_set_image(dset, t, k, MRI_FLOAT32, frunner);
	  frunner += *veclen * dims[0] * dims[1];
	}
    }
  case 3:
    {
      double* frunner= (double*)val;
      for (t=0; t<tmax; t++)
	for (k=0; k<kmax; k++) {
	  mri_set_image(dset, t, k, MRI_FLOAT64, frunner);
	  frunner += *veclen * dims[0] * dims[1];
	}
    }
  case 4:
    {
      short* frunner= (short*)val;
      for (t=0; t<tmax; t++)
	for (k=0; k<kmax; k++) {
	  mri_set_image(dset, t, k, MRI_INT16, frunner);
	  frunner += *veclen * dims[0] * dims[1];
	}
    }
  break;
  }
      
  mri_close_dataset(dset);
  return 0;
}







