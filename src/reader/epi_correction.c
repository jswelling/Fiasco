/* Revision 1.2 1997/12/2  Bryan J. Mock (GE Medical Systems)
 * Added some read() data functions for reading the rowflip.param file 
 * and also moved the reading in of the bandpass asymmetry, vrgf.dat, 
 * and ref.dat file to here from the main().  Started to add some baseline
 * processing for the fast receiver (rmv_baseln()). 
 *
 * Revision 1.1 - initial revision
 * Bryan J Mock (UW-Madison - 1/26/96): epi_correction.c contains:
 *
 *       -bp_correction(): routine to perform bandpass asymmetry 
 *        correction for receiver chain imperfections.  Uses bcrv*.dat
 *	  calibration files taken from the bin directory on the
 *        scanner. Code looks for EPIRECON_PATH environment variable
 *	  to tell it where to look for files.  Default is current
 *        working directory. Routine provided by GE.
 * 
 *      - epi_correction(): applies ahn correction using ref.dat
 *        obtained from GE scanner console in /usr/g/bin.  Assumes
 *	  ref.dat file is in the current working directory (not
 *	  necessarily where the raw data file is).  This routine does
 *        a 1D-fft to phase correct the data, fft's back to k-space,
 *	  then 2dfft result - allows for homodyne - slow but works
 *
 *      - phase_correct(): does rotation (complex multiplication) of
 *        fft-data by linear and const phase coeffcients obtained from
 *	  ref.dat
 *
 *      - spat_filter(): filter data by applying a Hamming window (-H) or
 *        fermi filter (-F) (using GE radius and width from header) if
 *        proper flag set from command line 
 *    
 *      - homodyne(): supports partial Ky epi scanning - provided by
 *        GE.  Modified to return image rather than k-space.    
 *
 *      - fovar_recon(): zero fills along x and picks proper image
 *        replication to account for fract. FOV (0.5 only) EPI
 *       
 *      - vrgf_filt(): applies filter to combine data sampled on ramps
 *                     Provided by GE.
 *
 */
static const char rcsid[] = "$Id: epi_correction.c,v 1.8 2007/03/22 00:03:47 welling Exp $";

#include<stdio.h>
/* #include "sunmath.h" */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "misc.h"
#include "bio.h"
#include "fmri.h"
#include "stdcrg.h"
#include "nr_sub.h"
#include "rcn.h"
#ifdef never
#include "io_signa_lx.h"
#endif
#define MAXS 4048			/*maximum number of samples*/
#define MAXV 512			/*maximum number of views*/
#define MAXC 4				/*maximum number of coils*/
#define MAXSHORT 32767
#ifndef M_PI
#define M_PI    3.14159265358979323846
#endif
#define UW_COMPILE 1
extern int debug;

/*************************************************************************/
/* complex rotation to correct epi phase error */
int phase_correct(fcomplex data[], int xres, float lp, float pc)
{
  fcomplex temp;  /* temporary storage of rotation values */
  fcomplex crot;  
  int i;

    for(i=0;i<xres;i++)
      {
        crot = nr_Complex(cos(lp*(i-xres/2.0)+pc),sin(lp*(i-xres/2.0)+pc));
	temp = data[i];
	data[i] = nr_Cmul(temp,crot);
      }
  return(0);
} /* END phase_correct() */

/*************************************************************************/
/* spatially filter each image based on FLAG setting, filter applied */
/* to raw k-space domplex data before reconstruction */
/* filters available: Hamming and None, Fermi still needs to be debugged */
int spat_filter(fcomplex data[], int xres, int yres, float fw, float
  fr, char *FLAG, int rep)
{ 

  int i, j, k;                           
  float hcoef1 = .54, hcoef2 = .46; /* hamming window coefficients */
  float ff;	                   /* fermi filter */
  float rad, rx, ry;	 	   /* distances for calc. fermi filter */		     
  float *filterx, *filtery;         /* filter val in each dimension*/
  int xind, yind;                   /* indices for filter application */
				    /*  must convert to radius for */
				    /*  fermi */
      		      
  /* allocate meory for filter values */
  if((filterx = (float *)calloc((size_t)xres,sizeof(float)))==NULL)
  {
      Abort("allocation failure for filter array.\n");
  }

  if((filtery = (float *)calloc((size_t)yres,sizeof(float)))==NULL)
  {
      Abort("allocation failure for filter array.\n");
  }

  /* determine filter to use and set up filter values */
  if(strcmp("-H",FLAG)==0) { /* hamming window */
     if (debug) Message("# Applying Hamming filter\n");
      for(i=0;i<xres;i++)
        filterx[i]=hcoef1-hcoef2*cos((float)(i+1)*(2*M_PI)/(float)xres);
      for(i=0;i<yres;i++)
        filtery[i]=hcoef1-hcoef2*cos((float)(i+1)*(2*M_PI)/(float)yres);
        /* printf("Hamming Filtered\n");*/
   }
      
/* determine Hamming indices and apply */
   if(strcmp("-H",FLAG)==0) {       /* hamming window */ 
      for(i=0; i<xres*yres; i++) {  /* set up indices */
        xind = (int)((floor)(i))%xres+1;
        yind = (floor)(i/(xres))+1;  
        /* apply Hamming filter */ 
	data[i].r = filterx[xind]*filtery[yind]*data[i].r;
        data[i].i = filterx[xind]*filtery[yind]*data[i].i;
     } 
    } /* end hamming if */
  
if(strcmp("-F",FLAG) == 0) { /* fermi filter*/
     if (debug) Message("# Applying Fermi filter\n");
     for (i=0;i<yres;i++) {
     ry=(i-yres/2+.5)*(i-yres/2+.5);
	for (j=0;j<xres;j++) {
	    k=i*xres+j;
            rx=(j-xres/2+.5)*(j-xres/2+.5);
	    rad=sqrt(rx+ry);
            ff=1.0/(1.0+exp((rad-fr)/fw));
	    data[k].r *= ff;
	    data[k].i *= ff;
	 }
     }    
} /* end fermi if */
  
  free(filterx);
  free(filtery);
  return(0);
} /* END spat_filter() */

/**************************************************************************/
/* Perform two consecutive 1D-FFT on data in complex array using Numerical */
/* Recipes routine, does complex rotation using linear and */
/* const. phase parameters from ref.dat file after transform along */
/* xdim, then does ydim transform */
fcomplex *epi_correction(fcomplex data[], int xres, int acq_xres, int
yres, int acq_yres, int scale, float *lp, float *poff)
{
  int i,j;
  fcomplex *ahn_data;    /* 1D storage array for processing */

  /* allocate vector for 1D transforms */
  if((ahn_data = (fcomplex *)calloc((size_t)xres,sizeof(fcomplex)))==NULL)
  {
      Abort("allocation failure for 1D storage array.\n");
    }

/* Transfer raw data into real and imaginary arrays one line at a time */
/* perfrom 1D transform, correct phase, transform back to k-space */

for(j=0; j<acq_yres; j++) {
  for(i=0; i<xres; i++)
      { /* store raw data in row format */
	ahn_data[i].r = (float)data[j*xres + i].r;
	ahn_data[i].i = (float)data[j*xres + i].i;
      }

        /* do transform */
        i1dfft(ahn_data,xres,sqrt(scale),-1); 

       /* apply phase correction for this line */
        phase_correct(ahn_data, xres, lp[j], poff[j]);

	/* transform back */
        i1dfft(ahn_data,xres,sqrt(scale),1); 

  /* fill data with result */
  for(i=0; i<xres; i++)    
      {   
	data[j*xres + i].r = (float)ahn_data[i].r;
	data[j*xres + i].i = (float)ahn_data[i].i;       
      }
}  /* end outer yres loop */

  /* clean up */
  free(ahn_data);

  /* pass back corrected data */
  return(data);

} /* end epi_correction() */


/*********************************************************************/
/* BJM add (1/26/96) */
/* homodyne recon for fractional Ky */  
int homodyne_recon(fcomplex *data,int xres,int yres, int acq_yres,int
scale, float hnw, int numover, int rep)
{
  int i,j;
  fcomplex *hold;      /* 1D holding array */
  fcomplex *homo_temp; /* 2D holding array */
  
/* allocate holding vector */
  hold  = (fcomplex *)calloc((size_t)xres,sizeof(fcomplex));
  homo_temp = (fcomplex *)calloc((size_t)xres*yres,sizeof(fcomplex));

/* setup trnasformed data */
for(j=0; j<yres; j++) {
  for(i=0; i<xres; i++)
      { /* store raw data in row format */
	hold[i].r = (float)data[j*xres + i].r;
	hold[i].i = (float)data[j*xres + i].i;
      }

   /* do transform along x dim */
   i1dfft(hold,xres,sqrt(scale),-1); 

/* fill homodyne array with result */
  for(i=0; i<xres; i++)    
      {   
	homo_temp[j*xres + i].r = (float)hold[i].r;
	homo_temp[j*xres + i].i = (float)hold[i].i;       
      }
}  /* end outer yres loop */

 /* pass fft result to homodyne function */
  homodyne(homo_temp, xres, acq_yres,scale, hnw, numover, rep);  

/* replace data with homodyne result */  
for (i=0; i<xres; i++) {
  for (j=0; j<yres; j++) {
    /* fill data array */
    data[i+yres*j].r = homo_temp[i+yres*j].r; 
    data[i+yres*j].i = homo_temp[i+yres*j].i; 
  }
}

  /* clean up */
  free(hold);
  free(homo_temp);

  return(0);
} /* END homodyne_recon */

/*************************************************************************/
/* Apply homodyne processing to partial ky data */
/*************************************************************************/

int homodyne(fcomplex hd_data[], int xres, int acq_yres, int scale,
  float hnw, int nover, int rep)
{ 

  int           j, l, i, alt = 1; 
  int           nrec, nview;         /* recon pts,  and acq_yres */
  float		j1, j2, x;           /* homodyne transition pts */
  float		*wh, *wl;            /* filters */
  fcomplex	*temp1, *temp2;      /* weighted fft arrays */

/* calculation of recon pts. and homdyne trans pts. */
  nview = acq_yres;
  nrec=log((float)nview)/log(2.0);
  nrec=2*pow(2.0,(float)nrec);
  nover=(int)(nview-nrec/2);
  j2=nview-2*hnw;
  j1=nrec-1-j2;

/* allocate memory for above arrays */
  wh = (float *)calloc((size_t)nrec,sizeof(float));
  wl = (float *)calloc((size_t)nrec,sizeof(float));
  temp1 = (fcomplex *)calloc((size_t)nrec,sizeof(fcomplex));
  temp2 = (fcomplex *)calloc((size_t)nrec,sizeof(fcomplex));

  Message("# Partial Nex (Homodyne) Reconstruction\n");
  Message("# Number of acquired views = %d\n",nview);
  Message("# All of k-space = %d views\n",nrec);
  Message("# Overscanned half of k-space by %d views\n",nover);
  Message("# Homodyne transition points and width = %5.1f %5.1f %5.1f\n",j1,j2,hnw); 
 /* printf("ky offset = %f\n",hshift);*/ 

/*calculate the partial nex filters*/
            for (j=0;j<nrec;j++) {
                wh[j]=2.0-1.0/(1.0+exp((j1-j)/hnw))-1.0/(1.0+exp((j2-j)/hnw));
                wl[j]=1.0/(1.0+exp((j-j2)/hnw))-1.0/(1.0+exp((j-j1)/hnw));
                }
            
/*calculate the column fft*/
	    for (i=0;i<xres;i++) {
	    	for (j=0;j<nrec;j++) {
		    l=j*xres+i;
/* weighted arrays temp1 & temp2 */
		    temp1[j].r = hd_data[l].r*wh[j];
		    temp1[j].i = hd_data[l].i*wh[j];
                    temp2[j].r = hd_data[l].r*wl[j];
                    temp2[j].i = hd_data[l].i*wl[j];
		    }
/* do fft's */
  	    	    i1dfft(temp1,nrec,sqrt(scale), -1);
	    	    i1dfft(temp2,nrec,sqrt(scale), -1);
                    
/* fill hd_data with homodyne result */
                for (j=0;j<nrec;j++) {
                    l=j*xres+i;
		    x= sqrt(temp2[j].r*temp2[j].r+temp2[j].i*temp2[j].i);
                    hd_data[l].r= alt*(temp2[j].r*temp1[j].r+temp2[j].i*temp1[j].i)/x;
                    hd_data[l].i = 0;
                    alt*= -1.0;
                    }
		}

     /*chop data along x */
     choprawx(hd_data,xres,nrec);

/* free allocated memory */
  free(wl);
  free(wh);
  free(temp1);
  free(temp2);

return(0);

} /* End fractional nex processing */

/************************ar*************************************************/
/* support fractional FOV in Ky */
/* coded at GE - modified by BJM 1/26/96 */
int fovar_recon(fcomplex *data,int rcxres,int interm_xres,int
  rcyres,int scale, float fovar)
{
  int i,j,k,l;
  int padw = 0;             /* pads data */
  fcomplex *interm_data;
  
/* allocate memory for intermediate data array */
  interm_data = (fcomplex *)calloc((size_t)rcxres*rcyres,sizeof(fcomplex));
 
/* fill interm_data with fractional FOV data */  
 for (j=0;j<rcyres;j++) 
   for (i=0;i<rcxres;i++) {
        interm_data[j*rcyres+i].r = data[j*rcyres+i].r;
        interm_data[j*rcyres+i].i = data[j*rcyres+i].i; 
      }

  /* pad along x - then chop resulting image */
  padx(&interm_data,interm_xres,rcxres,rcyres,0);

  /* transform result */
  i2dfft(interm_data,interm_xres,rcyres,scale,-1); 

  /* Load data array for display - pick off replication */ 
  padw = (rcxres - (int)((float)rcyres*fovar))/2.0;

  for (i=0;i<rcyres;i++) {
    for (j=0;j<rcxres;j++) {
      l = i*rcxres+j;
      k = (i-padw)*interm_xres+j;
      if (i>=padw && i<rcyres-padw) {
        data[l].r = interm_data[k].r;
        data[l].i = interm_data[k].i; 
      } else {
        data[l].r = 0.0;
        data[l].i = 0.0; 
	}
    }
  } 

/* clean up */
  free(interm_data);

  return(0);

} /* end fovar_recon() */

/***************************************************************************/
/* VRGF Support routine - combines ramp sampled data */
int vrgf_filt(fcomplex* data_out, fcomplex *data_in,
	      float **vrgf,int vrgfop,int vrgfip,int rcyres)
{

  int i,j, k, l;              
  float filt_tmp_real,filt_tmp_imag;
  fcomplex *vrgfobuf;

  /* allocate memory for vrgf coefficient data, output buffer, and */
  /* matrix passed back to main */
    vrgfobuf = (fcomplex *)calloc((size_t)vrgfop,sizeof(fcomplex)); 

  /* Apply VRGF to data */
  for (k = 0; k < rcyres; k++) {
    for (j = 0; j < vrgfop; j++) {
      filt_tmp_real = 0.0;
      filt_tmp_imag = 0.0;
     for (i = 0; i < vrgfip; i++) {
        l = k*vrgfip + i;
	filt_tmp_real += (float)data_in[l].r * vrgf[i][j];
	filt_tmp_imag += (float)data_in[l].i * vrgf[i][j];
        }
      vrgfobuf[j] = nr_Complex(filt_tmp_real,filt_tmp_imag); /*VRGF output pt j*/
      }

    for (i = 0; i < vrgfop; i++) {        /*VRGF result for this row*/
      l = k*vrgfop + i;
      data_out[l] = vrgfobuf[i];
      }
    }

/* clean up */
  free(vrgfobuf);
  return(0);

} /* end vrgf_filt() */
 


/***************************************************************************/
/*Bandpass Asymmetry Correction Routine Based on Kevin King's recon */
/*code. Adapted by DMW. 
Modified by BJM (UW - Madison 12/13/96)
   -added logprintf, pass rep number, slice number, and correction vectors
Note: would run faster if file opened was passed once instead of for
every rep.  May change this in the future. - NOW changed, bp file read
in function read_bp() below....
*/
int bp_corr(fcomplex *data,int xres,int yres,float bandw, float *amag,
  float *aphase, int nbpc, int rep, int slice)
{
  int i,j, k, l, m;
  float u,v,x,y,z;	/*temporary variables*/
/*  int fast_rec_flag;	flag indicating fast receiver */
  float dt;		/*a/d sample time (usec)*/
  int nshot;            /*number of views between filter direction changes*/
  int flip_flag;        /*flag to alternate filter direction*/ 
  int flip_filter;      /*flag to flip analog filter correction*/
/*  int nbpc;      */       /*number of filter correction file points*/
  int nsamp = xres;	/*number of samples per view*/
  int npts;		/*analog filter correction size*/
  int power;		/*fft power for analog filter correction*/
  int ncoil=1;		/*number of coils for multicoil*/
  int nview=yres;	/*number of views of acquired data*/
  int nrec=yres;        /*number of views reconstructed*/
  float hshift=0.0;	/*ky offset shift*/
/*  static float kr[MAXS*MAXV*(MAXC+1)]; raw data real part*/
/*  static float ki[MAXS*MAXV*(MAXC+1)];     raw data imag part*/
/*  static float amag[512]; */        /*filter correction magnitude*/
/*  static float aphase[512]; */      /*filter correction phase*/
  static float afr[2*MAXS];       /*filter correction real part*/
  static float afi[2*MAXS];       /*filter correction imag part*/
  static fcomplex temp1[MAXS];	/*temporary array for FFT*/
/*  FILE *f18;			pointer to aliased mode coefficients*/
/*  char message[80];		error message string*/
  
  dt = 1000.0/(2.0*bandw);
  nshot = 1;

/* BJM: only need to flip filter once - but flip flag is always 0 */
if((rep == 0) && (slice == 0)){
  flip_flag = 0;
  flip_filter = 1;
} else {
  flip_filter = 0;
  flip_flag = 0;
}

            if (flip_filter) {
                amag[0]=amag[nbpc-1];
                aphase[0]=aphase[nbpc-1];
                for (i=1;i<nbpc/2;i++) {
                    x=amag[i];
                    amag[i]=amag[nbpc-i];
                    amag[nbpc-i]=x;
                    x=aphase[i];
                    aphase[i]=aphase[nbpc-i];
                    aphase[nbpc-i]=x;
                    }
                }
/*filter correction*/
	    /* BJM 1/30/97: need log and exp to base 2, */
            /* log2() & exp2() dont exist in Solaris or SGI */
            /* math lib - only sunOS  */
            power=(log10((float)nsamp)/log10(2));
            npts=pow(2.0,(float)power);
            if (npts<nsamp) {
                npts*=2;
                power++;
                }
            if (debug) Message("# BP Filter Correction fft Size = %d\n",npts);
/* f18=fopen("fcor.as","w"); */
            for (i=0;i<npts;i++) {
                x=nbpc/2.+((float)(i-npts/2))/dt*(float)nbpc/(float)npts;
                j=x;
                x=x-j;
                y=(1.0-x)*amag[j]+x*amag[j+1];
                z=(1.0-x)*aphase[j]+x*aphase[j+1];
                afr[i]=y*cos(z*M_PI/180);
                afi[i]=y*sin(z*M_PI/180);
/*
 printf("i = %d areal = %f aimag = %f y = %f z = %f\n",i,afr[i],afi[i],y,z); 
 printf("mag = %f phase = %f x = %f\n\n", amag[i],aphase[i],x); */
                }
/* fclose(f18);
 f18=fopen("fft.as","w"); */
	    for (m=0;m<ncoil;m++) {
                for (i=0;i<nview;i++) {
                    for(j=0;j<npts;j++) {
                        temp1[j].r=0;
                        temp1[j].i=0;
                        }
                    for (j=0;j<nsamp;j++) {
			k=.5*(yres-nrec)*xres+.5*(xres-nsamp)+i*xres+m*xres*yres+j;
                        temp1[j].r=data[k].r;
                        temp1[j].i=data[k].i;
                        }
                    if ((flip_flag==1) && ((i/nshot)%2)) {
                        for (j=0;j<nsamp/2;j++) {
                            x=temp1[nsamp-1-j].r;
                            y=temp1[nsamp-1-j].i;
                            temp1[nsamp-1-j].r=temp1[j].r;
                            temp1[nsamp-1-j].i=temp1[j].i;
                            temp1[j].r=x;
                            temp1[j].i=y;
                            }
                        }
  	    	    i1dfft(temp1,nsamp,1, 1);

                    for (j=0;j<npts;j++) {
  /*printf("temp = %f \n",temp1[j]);*/
                        x=temp1[j].r*afr[j]-temp1[j].i*afi[j];
                        y=temp1[j].r*afi[j]+temp1[j].i*afr[j];
                        temp1[j].r=x;
                        temp1[j].i=y;
                        }

  	    	    i1dfft(temp1,nsamp,1, -1);

                    if ((flip_flag==1) && ((i/nshot)%2)) {
                        for (j=0;j<nsamp/2;j++) {
                            x=temp1[nsamp-1-j].r;
                            y=temp1[nsamp-1-j].i;
                            temp1[nsamp-1-j].r=temp1[j].r;
                            temp1[nsamp-1-j].i=temp1[j].i;
                            temp1[j].r=x;
                            temp1[j].i=y;
                            }
                        }
                    for (j=0;j<nsamp;j++) {
			k=.5*(yres-nrec)*xres+.5*(xres-nsamp)+i*xres+m*xres*yres+j;
                        data[k].r=temp1[j].r;
                        data[k].i=temp1[j].i;
                        }
                    }
		}
  /* fclose(f18); */
   
  return(0);

} /* end bp_corr() */
 
/******************** LOAD BANDPASS Asymmetry Data ********************/
/* BJM 1/3/97 - add this function to speed up recon when bp correction */
/* used */
int read_bp(float bandw, float *amag, float *aphase, int *nbpc, 
	    int *NBP_flag, char* bpassdir) 
{
  int index;
  float x;                /* temp var. */
  char cname[1024];	  /* analog filter correction file name */
  FILE *fbp;              /* file pointer */

  /* Determine bp correction file to read based on bandwidth */
  if (bandw > 62.5) /* fast receiver */
    if (bandw <= 100.0) sprintf(cname,"%s/bcrvf1.dat",bpassdir);
    else if (bandw <= 200.0) sprintf(cname,"%s/bcrvf2.dat",bpassdir);
    else if (bandw <= 300.0)  sprintf(cname,"%s/bcrvf3.dat",bpassdir);
    else if (bandw <= 400.0)  sprintf(cname,"%s/bcrvf4.dat",bpassdir);
    else sprintf(cname,"%s/bcrvf5.dat",bpassdir);
  else sprintf(cname,"%s/bcrvs1.dat",bpassdir); /* digital receiver */
  
  /*get analog filter correction values*/
  if (!(fbp=fopen(cname,"r"))) {
    Message("# Couldn't open BP Asymmetry correction file %s\n",cname);
    *NBP_flag = 0;
    return(0);
  }
  else Message("#      Bandpass Asymmetry Correction.\n"); 
  Message("#      BP Filename is %s\n",cname);
  index=0;
  while (fscanf(fbp,"%f %f %f",&x,amag+index,aphase+index)!=EOF) index++;
  *nbpc = index;
  Message("#      %d Points Read in from BP Filter Correction File\n",*nbpc);
  fclose(fbp);
  
  return(0);
}
/*********************End bandpass READ *******************************/

/* BJM 12/2/97: moved from main () to here */
/* LOAD VRGF Data for Ramped sampled data if necessary */
/* if rcxres < xres */  
int vrgf_read(float **vrgf,float *fbuf,int xres,int rcxres, int *VRGF_flag, 
	      char* fname)
{
  int i,j;
  FILE *fpvrgf;

  /* Read VRGF file generated by Signa  */
  if ((fpvrgf = fopen(fname,"r"))!=NULL){
     for(j=0;j<rcxres;j++) {
       FRdFloat32Array(fpvrgf,fbuf,xres);
        for(i=0;i<xres;i++) {
          vrgf[i][j] = fbuf[i];
        }
       *VRGF_flag = 1;
     }
     Message("#      Applying vrgf filter to epi data\n");
     fclose(fpvrgf);
     } else {
      /* set VRG_flag = 0 - don't even try to apply filter */
      *VRGF_flag =0;  
      Message("# Error Opening <%s> file\n",fname);
      Message("# Can't reconstruct vrgf data without this file!\n");
      Message("# Pull file from /usr/g/bin on signa and try again.\n");
    }
  return 0;
} /* end vrgf_read() */

/* BJM 12/2/97: moved from main () to here */
/* read reference data from local directory if available*/
int ref_dat_read(float *lp, float *poff,int xres,int yres,int *NO_REF_flag,
		 char* fname)

{
   int loffset;
   FILE *fpref;

   if ((fpref = fopen(fname,"rb"))!=NULL)
    {
      loffset = 512*sizeof(float); 
      fseek(fpref, (long) 0, (int) 0); 
      FRdFloat32Array(fpref,poff,yres);
      if (ferror(fpref)) { *NO_REF_flag= 0; return 0; }
      fseek(fpref, (long) loffset,(int) 0);
      FRdFloat32Array(fpref,lp,yres);
      if (ferror(fpref)) { *NO_REF_flag= 0; return 0; }
      fclose(fpref);
      Message("#      Reference data read from %s\n",fname);
    }
    else
    {
      *NO_REF_flag = 0;
    }

   return(0);
} /* end ref_dat_read() */

/* BJM 12/2/97: moved from main () to here */
/* read rowflip.param data from local directory if available*/
int rowflip_read(int *rwflp,int yres,int *ROWFLP_flag, int rep)
{
   int j;
   int header_lines = 4; /*skip header (3) + Ky = 0 (1) info */
   int max_char_per_line = 100;
   int kyline, flipval;
   char dummy[100];
   FILE *fpflp;

   if ((fpflp = fopen("rowflip.param","rb+"))!=NULL)
    {
      /* skip header of rowflip.param */
      for (j=0; j<header_lines; j++)
          fgets(dummy,max_char_per_line,fpflp);

      /* now read in ky line + flip info */
      for(j=0; j<yres; j++)  {
      	fscanf(fpflp,"%d %d\n",&kyline, &flipval);
	rwflp[j] = flipval;
      }
      fclose(fpflp);
      *ROWFLP_flag = 1;
      Message("# Reading rowflip.param file.\n");
    }
    else
    {
      Message("# Can't open rowflip.param file.\n"); 
      *ROWFLP_flag = 0;
    }

   return(0);

} /* end rowflip() */ 

/* BJM 2/27/98: avg baseline data */
fcomplex process_baseln(fcomplex *bl_data, int xres)
{
  int i,j;
  fcomplex avgIQ;

 /* add up I and Q channels */
 for(i=0;i<xres;i++)     
      {
       avgIQ.r += (float) bl_data[i].r;
       avgIQ.i += (float) bl_data[i].i;
      }

  return(avgIQ);

} 

/* BJM 12/2/97: uses conventional baseline processing to remove */
/*              dc artifact from fast receiver epi data */
int rmv_baseln(fcomplex *data, fcomplex bline, int xres, int yres)
{
  int i,j;
 
 /* subtract off baseline from every data point */

 for(j=0;j<yres;j++) {
    for(i=0;i<xres;i++)     
      {
        data[j*xres + i].r -= (float)bline.r;
        data[j*xres + i].i -= (float)bline.i;
      }
  }

  return(0);

} /* rmv_baseln() */

#ifdef never
/* BJM 12/2/97: get baseline data from pfile and then averge pos         */
/*              and neg views then subtract - this assumes the baselines */
/*              are not avg'd on the fly and shoved into the first view  */
/*              it's untested since this also requires a psd change      */
fcomplex *get_baseln(char filename[],int nbln, int xres, int yres, 
                     int psize, int headsize, int rep)
{
  int i,j,k;
  int offset;
  fcomplex *bview;
  scomplex *blraw;
  icomplex *iblraw;  /* in case of ext dyn range data */
  FILE *fpbl;

   /*Allocate memory for raw bl data */
   if (psize < sizeof(long)) {
       if ((blraw =(scomplex *)calloc((size_t)xres,sizeof(scomplex)))==NULL)
          {
            Abort("Allocation failure for raw (scomplex).\n");
          }   
   } 
   else 
   {
   /*Allocate memory for raw data - extended dynamic range */
    if ((iblraw =(icomplex *)calloc((size_t)xres,sizeof(icomplex)))==NULL)
        {
            Abort("Allocation failure for raw (icomplex).\n");
        }
  }

  /* Allocate memory for processed bl data */
   if ((bview =(fcomplex *)calloc((size_t)xres,sizeof(fcomplex)))==NULL)
    {
      Abort("allocation failure for bview data.\n");
    }
  
  /* offset into pfile by POOL_HEAR size - passed in from main*/
  offset = (int)headsize;     
  
  /* loop over the number of baselines collected */
  if(nbln > 0)
  for(k=0;k<nbln;k++) {

 /* open the pfile and get baseline data */
    if ((fpbl = fopen(filename,"rb+"))!=NULL)
      {
        fseek(fpbl, (long) offset, (int) 0);
          if(psize < sizeof(int))
	   FRdInt16Array(fpbl,blraw,2*xres);
          else
	   FRdInt32Array(iblraw,iblraw,2*xres);
        fclose(fpbl);
      }
    else
      return(NULL);

    /* avg. them together */
    for(i=0;i<xres;i++)     
      {
        bview[j*xres + i].r += (float)blraw[i].r;
        bview[j*xres + i].i += (float)blraw[i].i;
      }

    /* update offset by xres for next pass */
    offset += xres;
 } /* end loop over nbln */

  /* divide by numbr of blines/2 */
    for(i=0;i<xres;i++)     
      {
        bview[j*xres + i].r /= ((float)(nbln/2.0));
        bview[j*xres + i].i /= ((float)(nbln/2.0));
      }

if(nbln != 0)
 Message("# Subtracting %s baselines for fast receiver data correction\n", nbln);
else
  Message("# No baseline were collected.  Not doing baseline correction...\n");

  return(bview);

} /* get_baseln() */
#endif









