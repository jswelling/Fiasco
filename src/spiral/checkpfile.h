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

/*
 *	Definitions for checkpfiles.c
 *
 *
 */

#define SWAP(a,b) {tempr=(a);(a)=(b);(b)=tempr;}

#define  IM_sctime	(38916+22)    /* float - scan time (s) */
#define  IM_slthick	(38916+26)    /* float - Slice Thickness (mm) */
#define  IM_tr	    (38916+194)   /* 4 byte int, TR */
#define  IM_te	    (38916+202)   /* 4 byte int, TE */
#define  IM_flip    (38916+254)   /* 2 byte short flip angle */
#define  IM_spacing (38916+116)   /* spacing between scans */
#define  IM_nex 	(38916+218)   /* NEX*/

#define  IM_ctr_r  (38916+130)   /* R coord of CTR */
#define  IM_ctr_a  (38916+134)   /* A coord of CTR */
#define  IM_ctr_s  (38916+138)   /* S coord of CTR */

#define  IM_nor_r  (38916+142)   /* R coord of Norm */
#define  IM_nor_a  (38916+146)   /* A coord of Norm */
#define  IM_nor_s  (38916+150)   /* S coord of Norm */

#define  IM_tlc_r  (38916+154)   /* R coord of TLC */
#define  IM_tlc_a  (38916+158)   /* A coord of TLC */
#define  IM_tlc_s  (38916+162)   /* S coord of TLC */
#define  IM_trc_r  (38916+166)   /* R coord of TRC */
#define  IM_trc_a  (38916+170)   /* A coord of TRC */
#define  IM_trc_s  (38916+174)   /* S coord of TRC */
#define  IM_brc_r  (38916+178)   /* R coord of BRC */
#define  IM_brc_a  (38916+182)   /* A coord of BRC */
#define  IM_brc_s  (38916+186)   /* S coord of BRC */
#define  IM_psd    (38916+308)   /* pulse seq number */
#define  DAT_ACQ_TAB	10240	/* start of RDB_DATA_ACQ_TAB */
#define  GW1		10244	/* gw_point1[3] for 1st slice */

/* information about pfile needed to do recon/air and to do sanity checking */
typedef struct pfile_info {

  char	filename[256];
  long	filesize;
  long	minfilesize;
  char	date[10];
  char	time[8];
  int             x_dim;
  int             y_dim;
  int		z_dim;
  float	x_size;
  float	y_size;
  float	z_size;
  double	fov;
  int		nimages;
  int     rep;
  int     repmul;
  int		nspirals;
  int		acq_mat;
  int 	tr;
  int		te;
  int	sctime;
  short   flip;
  float   spacing;
  float 	nex;				
  float 	ctr[3];        
  float 	nor[3];        
  float 	tlc[3];        
  float 	trc[3];       
  float 	brc[3];       
  float 	tlz[3];         /* TLZ, calculated from tlc #slices, normal */
  float gw_first[3];	/* gw point 1 for 1st slice */
  short rotation;       /* rdb_hdr_rotation */
  short transpose;      /* rdb_hdr_transpose */
  float i_offset;       /* shift to isocenter after prescription */
  int	readout;
  int     coils;
  int     coil_rec_len;
  float	hwparams[6];    /* risetime, fsgcm, ts, gts, densamp, chop */
  char    psd[34];  
} P_INFO;


