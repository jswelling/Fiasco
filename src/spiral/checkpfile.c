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

static char rcsid[] = "$Id: checkpfile.c,v 1.14 2007/03/22 00:09:18 welling Exp $";

/*
 *	File:		checkpfile.c
 *
 *	Desc.:	-d <dir>	dir is the directory containing the pfiles
 *		-c  		sanity check a set of pfiles:
 *			 	-- same date, same xyz dimension,  same xyz voxel 
 *				   size, same FOV, same #spirals, ref file OK
 *              -s          print out summary for all pfiles, human readable
 *              -x          print out summary for all pfiles, csh readable
 *              -l          list each pfile info, no check, no ref finding
 *              -g          print out geometry info from 1st pfile
 *		-b <n>      print out TLC TRC BRC in RAS coords for slice n, 0<=n<z
 *              -y <n>      set size of XY dim for recon matrix (defaults to 64)
 *              -h          print help menu
 *
 *			Criteria for finding ref file:
 *			 	-- ref file has earliest date/time stamp
 *				-- ref pfile is smallest size
 *
 *	Input:		[-d <dir> -c -s -x -l -h] arbitrarily ordered list of pfiles
 *
 *  Return:		0 = sanity checks all OK or not performed (just -l flag)
 *			1 = at least 1 sanity check failed, or system error
 *
 *  -x flag output for shell scripts:
 *  <x_dim> <y_dim> <z_dim> <x_size> <y_size> <z_size> <#func images> <ref name>
 *
 *
 *  Kate Fissell 4/23/97
 *  adapted from Greg Hood's tiny.c program (Pittsburgh Supercomputing Center)
 *
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#if (SGI64 || SGI5 || SGIMP) 
#include <bstring.h>
#endif
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>

/* Fiasco 3.9 */
#include "misc.h"
#include "bio.h"
#include "array.h"
#include "frozen_header_info.h"

#include "rttraj.h"
#include "spiral.h"
#include "checkpfile.h"

void print_pinfo(P_INFO *p);
int compare_name();
int compare_size();
int compare_time();
int check(P_INFO *p, P_INFO *tmp, int n);

float ***Alloc3DFloatArray();


#define ACQ_OUT_OFF	5		/* allowable diff between acq. and out matrix */
#define DEFAULT_XYDIM   64   		/* default for XY dim of recon matrix */


int main(int argc, char* argv[])
{

  extern int optind;
  extern char *optarg;
  extern int errno;
  struct stat filestat;


  unsigned char header[FRZ_RDB_HEADER_SIZE_BYTES];
  Context	ic;
  P_INFO *pfiles,*tmp;

  int gflag, cflag, sflag, xflag, lflag, hflag, bflag;
  int numfiles,totalimages;
  int ret,i,j,c,fildes,cret;
  int gslice= -1; /* neg value means uninitialized */
  int xydim;
  char *pathname,filename[256];
  float b1[3], b2[3], b3[3];
  float ctr[3], nor[3];

  numfiles = argc-1;
  cflag=FALSE;
  sflag=FALSE;
  xflag=FALSE;
  lflag=FALSE;
  hflag=FALSE;
  gflag=FALSE;
  bflag=FALSE;
  pathname=NULL;
  cret = -1;
  totalimages = 0;
  xydim = DEFAULT_XYDIM;


  /* get input flags */
  while ((c = getopt(argc, argv, "d:csxglhb:y:")) != EOF)
    switch(c) {
    case 'd':
      pathname = optarg;
      break;
    case 'c':
      cflag=TRUE;
      break;
    case 's':
      sflag=TRUE;
      break;
    case 'x':
      xflag=TRUE;
      break;
    case 'g':
      gflag=TRUE;
      break;
    case 'l':
      lflag=TRUE;
      break;
    case 'b':
      bflag=TRUE;
      gslice = atoi(optarg);
      break;
    case 'y':
      xydim = atoi(optarg);
      break;
    case 'h':
      hflag = TRUE;
      break;
    case '?':
    default:
      hflag = TRUE;
    }

  numfiles = argc-optind;
  if (numfiles < 1) 
    fprintf(stdout, "\nNo input files given.\n");

  if ((argc < 2) || !(cflag || sflag || xflag || gflag || lflag || bflag || hflag) || (numfiles<1))
    hflag = TRUE;

  if (hflag) {
    fprintf(stdout, "\n\nUsage: %s [-d <dir> -c -s -x -l -b <slice> -h] <pfiles...>",argv[0]);
    fprintf(stdout, "\n \
\t-d <dir> Directory where pfiles are located.\n \
\t-c \tRun consistency check on pfiles. Errors printed to stdout.\n \
\t-s \tPrint to stdout summary information about pfiles.\n \
\t-x \tPrint to stdout 1 line of info for input to shell scripts.\n \
\t-g \tPrint to stdout 1 line of geom info: TLC, TRC, BRC, TLZ, NOR, CTR, Z_SIZE, SPACING, Z_DIM \n \
\t-l \tPrint to stdout info about each file.\n \
\t\tHW params = risetime, fsgcm, ts, gts, densamp, chop  (ask Doug) \n \
\t-b <slice> \tPrint TLC TRC BRC for slice.  0 <= slice < z\n \
\t-y <xy dim> \tSet xy size for recon matrix (defaults to 64) \n \
\t-h \tPrint this help info.\n \
\t<pfiles...> \tArbitrarily ordered list of pfiles.\n\n \
\t\tAll flags are optional, but you gotta pick at least one.\n \
\n");
    exit(0);
  }


  /* allocate structs for pinfo */
  pfiles = (P_INFO *)malloc((size_t)(numfiles*sizeof(P_INFO)));
  tmp = (P_INFO *)malloc((size_t)(numfiles*sizeof(P_INFO)));
  if ((pfiles == NULL) || (tmp == NULL)) {
    fprintf(stderr, "\nmalloc error: cannot allocate pfile info structs.\n");
    exit(1);
  }

  /* loop to fill p_info structs */
  for (i=0; i<numfiles; i++) {
    bzero((void *)(pfiles+i), sizeof(P_INFO));
    if (pathname != NULL)
      sprintf(filename,"%s/%s",pathname,argv[optind+i]);
    else
      sprintf(filename,"%s",argv[optind+i]);


    /* check file, get file size */
    fildes = open(filename, O_RDONLY); 
    if (fildes == -1) {
      perror(filename);
      exit(1);
    }
    ret = fstat(fildes, &filestat);
    if (ret != 0) {
      perror(filename);
      exit(1);
    }
    close(fildes);

    /* check that file is a file, not a dir */
    if (! S_ISREG(filestat.st_mode)) {
      fprintf(stderr, "\n%s is not a regular file.\n",filename);
      exit(1);
    }

    /* get size */
    pfiles[i].filesize = (long) filestat.st_size;


    /* read header info */
    ic.big_endian_input = TRUE;
    ic.big_endian_output = bio_big_endian_machine;
    ReadFileHeader(filename, &ic, header);

    /* copy header info to pinfo struct */
    if (strrchr(argv[optind+i], '/') == NULL)
      strcpy(pfiles[i].filename, argv[optind+i]);
    else
      strcpy(pfiles[i].filename, strrchr(argv[optind+i],'/')+1);

    pfiles[i].x_dim = xydim;
    pfiles[i].y_dim = xydim;
    pfiles[i].z_dim = ic.nslices;

    pfiles[i].fov = ic.opfov;
    pfiles[i].x_size = pfiles[i].fov / pfiles[i].x_dim;
    pfiles[i].y_size = pfiles[i].fov / pfiles[i].y_dim;
    pfiles[i].z_size = ic.slthick;

    strncpy(pfiles[i].time, ic.time, 8);
    strncpy(pfiles[i].date, ic.date, 10);
    pfiles[i].nimages = ic.nimages;
    pfiles[i].rep = ic.nph1;
    pfiles[i].repmul = ic.nphmult;
    pfiles[i].nspirals = ic.npr;
    pfiles[i].minfilesize = 
      ic.ncoils*ic.nslices*ic.nph1*(ic.nphmult*ic.npr + 1)*ic.ndat*4 + FRZ_RDB_HEADER_SIZE_BYTES;
    pfiles[i].tr = ic.tr;
    pfiles[i].te = ic.te;
    pfiles[i].flip = ic.flip;
    pfiles[i].spacing = ic.spacing;
    pfiles[i].nex = ic.nex;
    pfiles[i].readout = ic.ndat;
    pfiles[i].coils = ic.ncoils;
    pfiles[i].coil_rec_len = ic.coil_record_length;

    ctr[0] = 
      BRdFloat32(&header[FRZ_RDB_MRIMAGEDATATYPE_OFF + FRZ_IMAGEHEAD_CTR_R_OFF]);
    ctr[1] =
      BRdFloat32(&header[FRZ_RDB_MRIMAGEDATATYPE_OFF + FRZ_IMAGEHEAD_CTR_A_OFF]);
    ctr[2] =
      BRdFloat32(&header[FRZ_RDB_MRIMAGEDATATYPE_OFF + FRZ_IMAGEHEAD_CTR_S_OFF]);
    nor[0] = 
      BRdFloat32(&header[FRZ_RDB_MRIMAGEDATATYPE_OFF + FRZ_IMAGEHEAD_NORM_R_OFF]);
    nor[1] = 
      BRdFloat32(&header[FRZ_RDB_MRIMAGEDATATYPE_OFF + FRZ_IMAGEHEAD_NORM_A_OFF]);
    nor[2] = 
      BRdFloat32(&header[FRZ_RDB_MRIMAGEDATATYPE_OFF + FRZ_IMAGEHEAD_NORM_S_OFF]);

    ReadSliceCoords(&ic, header, 0, pfiles[i].gw_first, b2, b3,
		    &(pfiles[i].i_offset));
    ReadSliceCoords(&ic, header, pfiles[0].z_dim-1, pfiles[i].tlz, 
		    b2, b3, NULL);

    for(j=0; j<=2; j++) {
      pfiles[i].ctr[j] = ctr[j];
      pfiles[i].nor[j] = nor[j];
      pfiles[i].tlc[j] = ic.tlc[j];
      pfiles[i].trc[j] = ic.trc[j];
      pfiles[i].brc[j] = ic.brc[j];
    }
    pfiles[i].rotation = 
      BRdInt16(&header[FRZ_RDBHEAD_RDB_HDR_ROTATION_OFF]);
    pfiles[i].transpose = BRdInt16(&header[FRZ_RDBHEAD_RDB_HDR_TRANSPOSE_OFF]);
    strncpy(pfiles[i].psd, ic.psd, 34);
    pfiles[i].hwparams[0] = ic.risetime;
    pfiles[i].hwparams[1] = ic.fsgcm;
    pfiles[i].hwparams[2] = ic.ts;
    pfiles[i].hwparams[3] = ic.gts;
    pfiles[i].hwparams[4] = ic.densamp;
    pfiles[i].hwparams[5] = (float)ic.chop;
    totalimages += pfiles[i].nimages;

    /* allocate t2k array, get res, free it */
    ic.t2k = Alloc3DFloatArray(2, ic.npr, ic.ndat);
#ifdef never
    pfiles[i].acq_mat = getrttrajvd(ic.ndat, ic.npr, ic.ts, ic.gts, 0.0, ic.fsgcm, ic.opfov, ic.risetime, ic.densamp, 1.0, 1.0, ic.t2k[0], ic.t2k[1]);
#endif
    pfiles[i].acq_mat = getrttrajghg(ic.opxres, ic.npr, ic.ts, ic.gts, 
				     ic.fsgcm, ic.opfov, ic.slewrate,
				     ic.gts*21000, ic.gtype, 
				     ic.t2k[0], ic.t2k[1]);
    Free3DFloatArray(ic.t2k);

  }

  /* copy pfile info to tmp buffer for sorting */
  bcopy((void *)pfiles, (void *)tmp, numfiles*sizeof(P_INFO));


  /* l flag, just list pfile info */
  if (lflag)
    for (i=0; i<numfiles; i++) 
      print_pinfo(pfiles+i);


  /* b flag, get slice coords */
  if (bflag) {
    if ((gslice<0) || (gslice>pfiles[0].z_dim-1)) {
      fprintf(stderr, "\nSlice out of range [0 %d]: %d.\n",pfiles[0].z_dim-1,gslice);
      exit(1);
    }

    ReadSliceCoords(&ic, header, gslice, b1, b2, b3, NULL);
    
    fprintf(stdout, "%.3f %.3f %.3f   %.3f %.3f %.3f   %.3f %.3f %.3f",
	    b1[0],b1[1],b1[2],b2[0],b2[1],b2[2],b3[0],b3[1],b3[2]);
  }
  
  /* c flag, check consistency */
  if (cflag) {
    cret = check(pfiles,tmp,numfiles);
    if (cret == 0 )
      fprintf(stdout, "\npfiles are consistent.\n");
  }
  
  /* s flag, print summary, (must have run consistency check) */
  if (sflag) {
    if (cret == -1)
      cret = check(pfiles,tmp,numfiles);
    if (cret != 0)
      fprintf(stdout, "\npfiles inconsistent, cannot print summary.\n");
    else {
      if (numfiles > 1) {
	/* sort to put ref file first, then by time */
	qsort((void *)pfiles, numfiles, sizeof(P_INFO), compare_size);
	qsort((void *)(pfiles+1), numfiles-1, sizeof(P_INFO), compare_time);
	fprintf(stdout,"\n");
	fprintf(stdout, "\nFiles: \t\t\t%d pfile(s), %s (%s) - %s (%s), %s", numfiles-1, pfiles[1].filename,pfiles[1].time,pfiles[numfiles-1].filename,pfiles[numfiles-1].time,pfiles[0].date);
	fprintf(stdout, "\nReference file: \t%s (%d bytes) (%s)",pfiles[0].filename,pfiles[0].filesize,pfiles[0].time);
	fprintf(stdout, "\nTotal images: \t\t%d reference %d functional",pfiles[0].nimages, totalimages-pfiles[0].nimages);
	fprintf(stdout, "\nx,y,z dim: \t\t%d %d %d",pfiles[0].x_dim,pfiles[0].y_dim,pfiles[0].z_dim);
	fprintf(stdout, "\nx,y,z size: \t\t%.3f %.3f %.3f",pfiles[0].x_size,pfiles[0].y_size,pfiles[0].z_size);
	fprintf(stdout, "\nSeq.,TR,TE,flip: \t%s %.3f %.3f %d",
		pfiles[0].psd, (float)pfiles[0].tr/1000.0, 
		(float)pfiles[0].te/1000.0,pfiles[0].flip); 
	fprintf(stdout, "\nAcq., readout: \t\t%d %d", pfiles[0].acq_mat, pfiles[0].readout);
	fprintf(stdout,"\n");
	fprintf(stdout,"\n");
      }
      else
	print_pinfo(pfiles);
    }
  }


  /* x flag, print 1 line for use in shell scripts */
  /* this assumes the user ran a consistency check in a previous invocation */
  if (xflag) {
    qsort((void *)pfiles, numfiles, sizeof(P_INFO), compare_size);
    fprintf(stdout, "%d %d %d %.3f %.3f %.3f %d %s",pfiles[0].x_dim, pfiles[0].y_dim, pfiles[0].z_dim, pfiles[0].x_size, pfiles[0].y_size, pfiles[0].z_size, totalimages-pfiles[0].nimages, pfiles[0].filename);

  }

  if (gflag) {
    /* print tlc, trc, brc, tlz, nor, ctr, z_size, spacing, z_dim */
    fprintf(stdout, "%.3f %.3f %.3f    %.3f %.3f %.3f    %.3f %.3f %.3f    %.3f %.3f %.3f   %.3f %.3f %.3f   %.3f %.3f %.3f   %.3f %.3f %d", pfiles[0].tlc[0],pfiles[0].tlc[1],pfiles[0].tlc[2], pfiles[0].trc[0],pfiles[0].trc[1],pfiles[0].trc[2], pfiles[0].brc[0],pfiles[0].brc[1],pfiles[0].brc[2], pfiles[0].tlz[0],pfiles[0].tlz[1],pfiles[0].tlz[2],  pfiles[0].nor[0],pfiles[0].nor[1],pfiles[0].nor[2], pfiles[0].ctr[0],pfiles[0].ctr[1],pfiles[0].ctr[2],  pfiles[0].z_size,  pfiles[0].spacing, pfiles[0].z_dim);
  }

  /* return results of consistency check if it was run */
  if (cret == -1)
    return(0);
  return(cret);
}

/****************************************************************************/
/*                                                                          */
/*                                                                          */
/*                                                                          */
/****************************************************************************/
int check(P_INFO *p, P_INFO *tmp, int n)
{
  int errors;
  int wrap,i;
  char smalltime[256];

  errors=0;


  /* check that each file is >= its min size */
  for (i=0; i<n; i++)
    if (p[i].filesize < p[i].minfilesize) {
      fprintf(stdout, "\nError: %s is too small: %d (%d)",p[i].filename,p[i].filesize,p[i].minfilesize);
      errors++;
    }

  /* check that acq matrix size is the same for all and not more than 5 off dim */
  for (i=1; i<n; i++)
    if (p[i].acq_mat != p[i-1].acq_mat) {
      fprintf(stdout, "\nError: %s and %s have different size acquisition matrices: %d %d.",p[i].filename,p[i-1].filename,p[i].acq_mat,p[i-1].acq_mat);
      errors++;
    }
  if ((p[0].acq_mat > p[0].x_dim+ACQ_OUT_OFF) || (p[0].acq_mat < p[0].x_dim-ACQ_OUT_OFF)) {
    fprintf(stdout, "\nWarning: Discrepancy between acquisition matrix (%d) and output matrix (%d) is too great.",p[0].acq_mat, p[0].x_dim);
  }


  /* if only 1 file, these checks are enough */
  if (n == 1)
    return(errors);


  /* check that xyz, dim, size is the same for set */
  for (i=1; i<n; i++) {
    if (p[i].x_size != p[i-1].x_size) {
      fprintf(stdout, "\nError: x_size changed: %s (%f) %s (%f)",p[i].filename, p[i].x_size, p[i-1].filename, p[i-1].x_size);
      errors++;
    }
    if (p[i].y_size != p[i-1].y_size) {
      fprintf(stdout, "\nError: y_size changed: %s (%f) %s (%f)",p[i].filename, p[i].y_size, p[i-1].filename, p[i-1].y_size);
      errors++;
    }
    if (p[i].z_size != p[i-1].z_size) {
      fprintf(stdout, "\nError: z_size changed: %s (%f) %s (%f)",p[i].filename, p[i].z_size, p[i-1].filename, p[i-1].z_size);
      errors++;
    }
    if (p[i].x_dim != p[i-1].x_dim) {
      fprintf(stdout, "\nError: x_dim changed: %s (%d) %s (%d)",p[i].filename, p[i].x_dim, p[i-1].filename, p[i-1].x_dim);
      errors++;
    }
    if (p[i].y_dim != p[i-1].y_dim) {
      fprintf(stdout, "\nError: y_dim changed: %s (%d) %s (%d)",p[i].filename, p[i].y_dim, p[i-1].filename, p[i-1].y_dim);
      errors++;
    }
    if (p[i].z_dim != p[i-1].z_dim) {
      fprintf(stdout, "\nError: z_dim changed: %s (%d) %s (%d)",p[i].filename, p[i].z_dim, p[i-1].filename, p[i-1].z_dim);
      errors++;
    }
  }


  /* check that date is the same for set */
  for (i=1; i<n; i++) 
    if (strncmp(p[i].date,p[i-1].date,10)) {
      fprintf(stdout, "\nError: date changed: %s (%s)  %s (%s)",p[i].filename, p[i].date, p[i-1].filename, p[i-1].date);
      errors++;
    }

  /* check for duplicate names */
  qsort((void *)tmp, n, sizeof(P_INFO), compare_name);
  for (i=1; i<n; i++) 
    if (! strcmp(tmp[i].filename,tmp[i-1].filename)) {
      fprintf(stdout, "\nError: duplicate filename: %s %s",tmp[i].filename, tmp[i-1].filename);
      errors++;
    }

  /* check that there is unique smallest file */
  qsort((void *)tmp, n, sizeof(P_INFO), compare_size);
  if (tmp[0].filesize == tmp[1].filesize) {
    fprintf(stdout, "\nError: pfile set does not include 1 unique smallest file; %s %s are both %d bytes.\n",tmp[0].filename,tmp[1].filename,tmp[1].filesize);
    errors++;
  }


  /* check that smallest size is also 1st time  */
  /* note: 2 pfiles can have the same time, since it is only precise to minutes */
  qsort((void *)tmp, n, sizeof(P_INFO), compare_time);
  strcpy(smalltime, tmp[0].filename);
  qsort((void *)tmp, n, sizeof(P_INFO), compare_size);

  /* if names don't match could be a problem */
  if (strcmp(smalltime, tmp[0].filename)) {

    /* if time matches 2nd, try that */
    if (!strcmp(tmp[0].time,tmp[1].time)) {
      /* try 2nd name */
      if (strcmp(smalltime, tmp[1].filename)) {
	fprintf(stdout, "\nError: Earliest file (%s) is not smallest file (%s)",smalltime,tmp[0].filename);
	errors++;
      }
    }
    /* times don't match, error */
    else {
      fprintf(stdout, "\nError: Earliest file (%s) is not smallest file (%s)",smalltime,tmp[0].filename);
      errors++;
    }
  }


  /* check that time ordering = name ordering */
  qsort((void *)tmp, n, sizeof(P_INFO), compare_size); /* get ref first */
  /* sort all but ref */
  qsort((void *)(tmp+1), n-1, sizeof(P_INFO), compare_time);
  wrap = 0;
  for (i=1; i<n; i++) 
    if (strcmp(tmp[i].filename,tmp[i-1].filename) < 1) {
      wrap++;
    }
  if (wrap > 1) {
    fprintf(stdout, "\nError: files break ascending order twice:\n");
    for (i=0; i<n; i++)
      fprintf(stdout, "%s (%s)  ", tmp[i].filename, tmp[i].time);
    fprintf(stdout, "\n");
    errors++;
  }
  if (wrap == 1)
    if (strcmp(tmp[n-1].filename,tmp[0].filename) >= 0) {
      fprintf(stdout, "\nError: filenames have wrapped, but last has overtaken first: %s %s\n",tmp[n-1].filename,tmp[0].filename);
      for (i=0; i<n; i++)
	fprintf(stdout, "%s (%s)  ", tmp[i].filename, tmp[i].time);
      fprintf(stdout, "\n");
      errors++;
    }


  fprintf(stdout, "\n");
  return(errors);
}

/****************************************************************************/
/*                                                                          */
/*                                                                          */
/*                                                                          */
/****************************************************************************/
int compare_name(p1,p2)
     P_INFO *p1, *p2;
{

  return(strcmp(p1->filename,p2->filename));


}


int compare_size(p1,p2)
     P_INFO *p1, *p2;
{

  if (p1->filesize < p2->filesize)
    return(-1);
  if (p1->filesize > p2->filesize)
    return(1);

  return(0);
}


int compare_time(p1,p2)
     P_INFO *p1, *p2;
{
  int h1,m1,h2,m2;
  char *t1, *t2;

  t1 = strdup(p1->time);
  t2 = strdup(p2->time);
  h1 = atoi(strtok(t1,":"));
  m1 = atoi(strtok(NULL,":"));
  h2 = atoi(strtok(t2,":"));
  m2 = atoi(strtok(NULL,":"));

  if (h1 < h2)
    return(-1);
  if (h1 > h2)
    return(1);
  if (m1 < m2)
    return(-1);
  if (m1 > m2)
    return(1);

  /* if time is identical, sort by name */
  return(strcmp(p1->filename,p2->filename));
}
/****************************************************************************/
/*                                                                          */
/*                                                                          */
/*                                                                          */
/****************************************************************************/
void print_pinfo(P_INFO *p)
{

  fprintf(stdout, "\nName: \t\t\t%s",p->filename);
  fprintf(stdout, "\nDate: \t\t\t%s  %s",p->date,p->time);
  fprintf(stdout, "\nx,y,z dim: \t\t%d %d %d",p->x_dim,p->y_dim,p->z_dim);
  fprintf(stdout, "\nrepeats*mul = #img: \t%d*%d = %d",p->rep, p->repmul,p->nimages);
  fprintf(stdout, "\nx,y,z size: \t\t%.3f %.3f %.3f",p->x_size,p->y_size,p->z_size);
  fprintf(stdout, "\nSlice spacing (mm): \t%.3f", p->spacing);
  fprintf(stdout, "\nFOV:  \t\t\t%.3f",p->fov);
  fprintf(stdout, "\nCenter coord: \t\t(%.3f, %.3f, %.3f)",p->ctr[0],p->ctr[1],p->ctr[2]);
  fprintf(stdout, "\nNormal coord: \t\t(%.3f, %.3f, %.3f)",p->nor[0],p->nor[1],p->nor[2]);
  fprintf(stdout, "\nTop left coord: \t(%.3f, %.3f, %.3f)",p->tlc[0],p->tlc[1],p->tlc[2]);
  fprintf(stdout, "\nTop right coord: \t(%.3f, %.3f, %.3f)",p->trc[0],p->trc[1],p->trc[2]);
  fprintf(stdout, "\nBottom right coord: \t(%.3f, %.3f, %.3f)",p->brc[0],p->brc[1],p->brc[2]);
  fprintf(stdout, "\nTLZ (from gw):   \t(%.3f, %.3f, %.3f)",p->tlz[0],p->tlz[1],p->tlz[2]);
  fprintf(stdout, "\ni_offset, rot, trans: \t%.3f  %d  %d",p->i_offset,p->rotation,p->transpose);
  fprintf(stdout, "\ngw_first: \t\t(%.3f, %.3f, %.3f)",p->gw_first[0],p->gw_first[1],p->gw_first[2]);
  fprintf(stdout, "\nPulse seq. #:\t\t%s",p->psd);
  fprintf(stdout, "\nTR (millisec): \t\t%.3f", (float)(p->tr/1000.0));
  fprintf(stdout, "\nTE (millisec): \t\t%.3f", (float)(p->te/1000.0));
  fprintf(stdout, "\nFlip angle (degrees): \t%d", p->flip);
  fprintf(stdout, "\nNumber of spirals: \t%d",p->nspirals);
  fprintf(stdout, "\nNEX: \t\t\t%.3f",p->nex);
  fprintf(stdout, "\nAcquisition matrix: \t%d",p->acq_mat);
  fprintf(stdout, "\nCoil #, rec. length: \t%d %d",p->coils,p->coil_rec_len);
  fprintf(stdout, "\nReadout: \t\t%d",p->readout);
  fprintf(stdout, "\nHW params: \t\t%.3f %.3f %.3f %.3f %.3f %d",p->hwparams[0],p->hwparams[1],p->hwparams[2]*1e6,p->hwparams[3]*1e6,p->hwparams[4],(int)(p->hwparams[5]));
  fprintf(stdout, "\nFile size [min size]: \t%d [%d]",p->filesize,p->minfilesize);
  fprintf(stdout,"\n");



  return;

}

/****************************************************************************/
/*                                                                          */
/*                                                                          */
/*                                                                          */
/****************************************************************************/
/****
int
Misc_Round (float f)
{
  int i;
 
  i = floor(f + (1.0/2.0));
  if ((((float) i) - f) == (1.0/2.0))
    i &= ~1;
  return(i);
}
***/

/****************************************************************************/
/*                                                                          */
/*                                                                          */
/*                                                                          */
/****************************************************************************/
/* Abort is used for unrecoverable errors which the user should be notified
   about */
/****
void
Misc_Abort (char *fmt, ...)
{
  va_list args;
  FILE *f;
  Filename fn;
  char hn[128];
 
  va_start(args, fmt);
  vfprintf(stdout, fmt, args);
  va_end(args);
 
  exit(1);
}
***/

/****************************************************************************/
/*                                                                          */
/*                                                                          */
/*                                                                          */
/****************************************************************************/
/* Error is used for recoverable errors which the user should be notified
   about */

/****
void
Misc_Error (char *fmt, ...)
{
  va_list args;
 
  va_start(args, fmt);
  vfprintf(stdout, fmt, args);
  va_end(args);
}
***/


