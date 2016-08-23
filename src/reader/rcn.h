#define RECON_FLAG
#include "control_data.h"
#include <stdio.h>

/* read & write routines */
#ifdef never
fcomplex *readraw(int,int,int,int,float,int,int,int,int,int,char *);
#endif
int readcl(int, char **, char *);
int writeimage(float *, int, int, char *, float,float,int,int);
int writeimage_sdt(float *, int, int, char *, float,float,int,int,int,int, FILE *,int);
void logprintf(int, char *, ...);  /* print to file and stdout */
void stimulate_header(FILE *,int,int,float,float,int,int,float,float,float, char *);
int myreadcl(int, char **, char *, char *, float *, int *, int *, int
  *, int *,int *, int *,char *);

/* chopping, flipping, rotating, ect */
int fliprowx(fcomplex *,int,int,int,int,int *,int);
int fliprawy(fcomplex *,int ,int);
int fliprawx(fcomplex *,int ,int);
int reorder_ileave(fcomplex *, int, int, int);
int choprawx(fcomplex *,int,int);
int choprawy(fcomplex *,int,int);
int choprawy_ileave(fcomplex *,int,int,int,int *,int); 
int padx(fcomplex **,int,int,int, int);
int pady(fcomplex **, int, int, int, int);
int xshift(int,fcomplex *,int,int);
int xshift2(int,fcomplex **,int,int);
int dorot(fcomplex *,int,int,float);
int rotate_image(float *,int,int,float);
int pixshiftx(fcomplex *,int,int,int);
int pixshifty(fcomplex *,int,int,int);

/* fft routines */
int i2dfft(fcomplex *,int,int,float,int);
int i1dfft(fcomplex *, int, float, int);

/* image processing stuff */
float *mag(fcomplex *, int, int);
float *lograw(fcomplex *,int,int);
fcomplex *avg_data(fcomplex *, fcomplex *, int, int);
fcomplex *parse_data(fcomplex *, fcomplex *, int, int);
float *cd(fcomplex *, fcomplex *,int ,int);
float *perp_cd(fcomplex *, fcomplex *, int, int);
float *bk_phase(fcomplex *, fcomplex *, int, int);
float *pd(fcomplex *,fcomplex *, int, int);
float *pamag(fcomplex **,int, int, int);
float *pacd(fcomplex **,fcomplex **,int,int,int);
float *papd(fcomplex **,fcomplex **,int, int,int);
int spat_filter(fcomplex *,int,int, float, float, char *,int);
float *avg(float **,int, int, int);
float Carg(fcomplex);

/* epi correction routines */
fcomplex *epi_correction(fcomplex *,int,int,int,int,int, float *,float *);
int homodyne_recon(fcomplex *, int, int, int, int, float, int, int);
int homodyne(fcomplex *,int,int,int,float,int,int);
int fovar_recon(fcomplex *, int, int,int,int,float);
int phase_correct(fcomplex *, int, float, float);
int rt_relphase_cor(fcomplex *,fcomplex *,int,int);
int vrgf_filt(fcomplex*,fcomplex *,float **,int,int,int);
int bp_corr(fcomplex *,int,int,float,float *,float *,int,int,int); 
int read_bp(float, float *, float *,int *, int *, char *);
int rowflip_read(int *,int,int *,int);
int ref_dat_read(float *,float *,int,int,int *,char*);
int vrgf_read(float **,float *,int,int,int *,char*);
fcomplex *get_baseln(char [],int,int,int,int,int,int);
int rmv_baseln(fcomplex *,fcomplex,int,int);
fcomplex process_baseln(fcomplex *,int);

/* Short int complex data structure to be compatible with complex.c */
/* Numerical Recipes data */
typedef struct SCOMPLEX 
{
  short r,i;
}
scomplex;

typedef struct ICOMPLEX 
{
  int r,i;
}
icomplex;










