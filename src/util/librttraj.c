#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "misc.h"
#include "rttraj.h"

/* Used by getrttrajghg */
#define MAX_PG_WAMP 32766
#define GRESMAX 21000
#define DTSMAX 16384

/* Used by getrttrajvd */
#define MAXDECRATIO 32 /* maximum allowed decimation of input ts */
#define MAX_PG_WAMP 32766
#define GAM 4257.0       

static char rcsid[] = "$Id: librttraj.c,v 1.3 2003/05/13 20:53:51 welling Exp $";

static float *fspace( int size );

/*
 * generates a k-space trajectory for a scan that 
 * does real-time spiral gradient generation.   
 * returns the number of samples
 *
 * does dyun algo using glover code
 * This bit put together by D. Noll, 6/99
*/

/*   multi- shot spiral design 
*    uses Duyn's approximate slewrate limited design 
*    augmented with archimedian gmax limit 
*    inputs (args) 
*        opfov = FOV, cm 
*        opxres = matrix size 
*	 Tmax = longest acquisition allowed, s 
*	 dts = output sample spacing, s 
*        gtype = trajectory type 
*    inputs (CVs) 
*        nl = number of interleaves 
*        gamp = design grad max, G/cm 
*        gslew = design slew rate, mT/m/ms 
*	 nramp = number of rampdown points 
*    outputs 
*        Gx, Gy 
*        grev 
*    time is in sec 
* 
*		rev 0 12/26/98	original 
*		rev 1 4/15/99	little better calc of ts 
*/ 

extern int getrttrajghg(int opxres, int nl, double dts, double gts, 
			double gamp, double opfov,
                        double gslew, double Tmax, int gtype, 
			float** kx, float**ky )
{
    
  float gamma = 2*M_PI*4.257e3; 
  float gambar = 4.257e3; 
  float dt; 
  float S0;
  int i, j, n; 
  int gres,indl;
  float Ts, T, ts, t, theta, dthdt, thetas, x, y; 
  float a1, a2, beta, Gmax, bx, by, c, s; 
  float kxt, kyt; 
  float gx[2*GRESMAX], gy[2*GRESMAX]; 
  int Gx[GRESMAX], Gy[GRESMAX]; 
  float Kx1[GRESMAX], Ky1[GRESMAX]; 
  float gabs, gmax; 
  short gxi, gyi; 
  float q; 
  float *r1, *r2;
#define GTYPE	0 
#define GTYPE_REVERSED 1
#define LIMIT(x, xmin, xmax)   ( (x<xmax)? ((x>xmin)? x:xmin):xmax ) 
 
  q = 5;		/* nice number  */ 
  S0 = gslew*100; 
  dt = gts*.5; 
 
  r1 = fspace(nl);
  r2 = fspace(nl);
  for (i = 0; i < nl; i++)
    {
      r1[i] = cos(2*M_PI*((double)i/nl+0.25));
      r2[i] = sin(2*M_PI*((double)i/nl+0.25));
      if (!kx[i]) kx[i] = fspace(DTSMAX);
      if (!ky[i]) ky[i] = fspace(DTSMAX);
    }
  
  /*
    fprintf(stderr,"entering genspiral  opxres, opfov, nl = %d  %6.2f cm  %d\n",
    opxres, opfov, nl); 
    fprintf(stderr,"dt = %f, dts = %f, gslew = %f, gamp = %f, Tmax = %f\n",
    dt*1e6, dts*1e6, gslew, gamp, Tmax);
  */
 
  if(gtype != GTYPE && gtype != GTYPE_REVERSED)  { 
    printf("wrong trajectory type %d\n",gtype); 
    exit(0);
  } 
 
/*  slew-rate limited approximation  */ 
 
  Ts = .666667/nl*sqrt(pow((double)(M_PI*opxres), (double)3.0)/(gamma*opfov*S0)); 
  if(Ts > Tmax)  { 
    printf("slew limited readout too long"); 
    exit(0);
  } 
  a2 = opxres*M_PI/(nl*pow(Ts, .666667)); 
  a1 = 1.5*S0/a2; 
  beta = S0*gamma*opfov/nl; 
  Gmax = a1*pow(Ts, .333333); 
  gmax = 0; 
  i = 0; 
  for (t = 0; t<=Ts; t+=dt)  { 
    x = pow(t, 1.333333); 
    theta = .5*beta*t*t/(q + .5*beta/a2*x); 
    y = q+.5*beta/a2*x; 
    dthdt = beta*t*(q+.166667*beta/a2*x)/(y*y); 
    c = cos(theta); 
    s = sin(theta); 
    gx[i] = nl/(opfov*gamma)*dthdt*(c - theta*s); 
    gy[i] = nl/(opfov*gamma)*dthdt*(s + theta*c); 
    gabs = hypot(gx[i], gy[i]); 
    if(gabs>=gamp)  { 
      if(gmax==0)  gmax = hypot(gamp/theta, gamp); 
      if(gabs>gmax) break; 
    } 
    ts = t; 
    thetas = theta; 
    i++; 
  } 
 
/*  gmax limited approximation  */ 
 
  if(Gmax > gamp)  { 
    T = ((M_PI*opxres/nl)*(M_PI*opxres/nl) - thetas*thetas)/(2*gamma*gamp*opfov/nl) + ts; 
    if(T > Tmax)  { 
      printf("gmax limited readout too long"); 
      exit(0);
    } 
    for (t=ts+dt; t<=T; t+=dt)  { 
      theta = sqrt(thetas*thetas + 2*gamma*gamp*opfov*(t-ts)/nl); 
      c = cos(theta); 
      s = sin(theta); 
      gx[i] = gamp*(c/theta - s); 
      gy[i++] = gamp*(s/theta + c); 
    } 
  } 
 
/*  decimate by 2 to get back to 4us sampling */ 
 
  n = 0; 
  Kx1[0] = Ky1[0] = 0.0;
  for (j=0; j<i; j+=2)  { 
    x = LIMIT(gx[j], -gamp, gamp); 
    gxi = x*MAX_PG_WAMP/gamp; 
    Gx[n] = 2*(gxi/2); 
    y = LIMIT(gy[j], -gamp, gamp); 
    gyi = y*MAX_PG_WAMP/gamp; 
    Gy[n] = 2*(gyi/2); 
    Kx1[n+1] = Kx1[n] + Gx[n]*(gamp*gts*opfov*gambar/MAX_PG_WAMP);
    Ky1[n+1] = Ky1[n] + Gy[n]*(gamp*gts*opfov*gambar/MAX_PG_WAMP);
    n++;
  } 
 
  gres = n - 1; 
   
  for (i = 0; i < nl; i++)
  {
    kx[i][0] = 0.0;
    ky[i][0] = 0.0;
  }
  n=0;
  for (ts=dts; ts<=t; ts+=dts)  { 
    indl = floor(ts/gts);
    if (indl < gres) {
      kxt = Kx1[indl]*(1- ts/gts + indl) + Kx1[indl+1]*(ts/gts - indl);
      kyt = Ky1[indl]*(1- ts/gts + indl) + Ky1[indl+1]*(ts/gts - indl);
      for (i = 0; i < nl; i++)
      {
	kx[i][n+1] = (kxt*r1[i] - kyt*r2[i]);
	ky[i][n+1] = (kxt*r2[i] + kyt*r1[i]);
	n++;
      }
    }
  }
  /*
  fprintf(stderr,"Tslew, Ttot (ms), tot points, nsamples =  %f  %f  %d %d\n", 
	 ts*1000, t*1000, gres, n);
  */

#ifdef never
  if (gtype == GTYPE_REVERSED) {
    /*reverse the k-space*/
    for(i = 0; i < gres; i++) {
      for (j = 0; j < nl; j++) {
	float tempx = kx[j][i];
	float tempy = ky[j][i];
	kx[j][gres-i-1] = tempx;
	ky[j][gres-i-1] = tempy;
      }
    }
  }
#endif

  return n;
}    

/*
 * generates a k-space trajectory for a scan that 
 * does real-time spiral gradient generation.   
 * returns the k-space matrix size.
 *
 * currently just does a slew-rate design.
 *
 * reads rawh as an external structure to get information about
 * gradient system from the rhuser variables.
 *
 * notes:
 *  (0) gts/ts and dr stuff is *not* debugged.  in fact, it known to be wrong.
 *  (1) ts should be an integral multiple of gts
 *
 * Copyright 1996, Craig Meyer, Leland Stanford Junior University.
 * Copyright 1996, Douglas Noll, University of Pittsburgh
*/

int getrttrajvd (int npts, int nl, double ts, double gts, double angle, 
		 double fsgcm, double opfov, int risetime, int densamp,
		 double kxscale, double kyscale, float **kx, float **ky,
		 int* res) 
{
    float *r1, *r2, kxt, kyt;
    int i, m, n, dnpts, loop, resolution;
    int dr;
    int dentrans, den1, kdenrad;
    double decratio;
    double g0, thetan_1, theta, deltheta, taun_1, taun, tauhat;
    double absg,gtilde,B,t1,t2,t3,ac,tgx,tgy;
    double tkx, tky, oldkx, oldky, *pres;
    double S, OM, s, om, A;
    double OMF, omf, denrad, scoffset, denoffset, scthat,fractrans, realn, ksv;
    double distance;
    short gx, gy;
    int samplecount= 0;

    dr = ts/gts;
    if (dr != ts/gts)
    {
	fprintf(stderr, "ts must be an integral multiple of gts\n");
	exit(1);
    }
    
    r1 = (float *) malloc(nl*sizeof(float));
    r2 = (float *) malloc(nl*sizeof(float));
    for (i = 0; i < nl; i++)
    {
	r1[i] = cos(2*M_PI*((double)i/nl+angle/360.0+0.25));
	r2[i] = sin(2*M_PI*((double)i/nl+angle/360.0+0.25));
	if (!kx[i]) kx[i] = fspace(npts);
	if (!ky[i]) ky[i] = fspace(npts);	
    }

    /* initialize kx and ky arrays to 0.0 */
    for (i = 0; i < nl; ++i)
      for (m = 0; m < npts; ++m)
	{
	  kx[i][m] = 0.0;
	  ky[i][m] = 0.0;
	}

    /* other rhuser cv's are read in calling program */
#ifdef never
    fprintf(stderr,"fsgcm = %.1f, risetime = %d, opfov = %.1f, densamp = %d\n",
	   fsgcm, risetime, opfov, densamp);
    fprintf(stderr,"kxscale = %f kyscale = %f\n", kxscale, kyscale); 
#endif
    
    A = MAX_PG_WAMP;
    S = (gts/1e-6)*A/risetime;
    OM = (2*M_PI)/(nl) * (opfov/10.0)/(1/(GAM*fsgcm*gts));
    /* distance for one k-space unit */
    distance = 1.0 / (opfov/10.0*GAM*fsgcm*gts/A); 

    OMF = OM*nl;
    dentrans = densamp/2;

    ac = A; loop = 1; decratio = 1;
    while (loop)
    {
	loop = 0; dnpts = dr*npts*decratio; om = OM/decratio; s = S/decratio;
	g0 = 0; gx = g0; gy = 0; absg = hypot(g0, 0.0);
	oldkx = 0; oldky = 0; tkx = gx; tky = gy; kxt = tkx; kyt = tky;
	thetan_1 = 0; taun = 0; n = 0;
	omf = OMF/decratio; den1 = 0;
	while (n < (dnpts-1))
	{ 
	    taun_1 = taun; taun = hypot(tkx,tky)/A; tauhat = taun;
	    realn = ((float) n)/((float) decratio);
	    if (realn > (float) densamp) {
	      if (den1 == 0)  {
		ksv = taun;
		den1 = 1;
	      }
	      if (realn > (float) (densamp+dentrans)) {
		scoffset = scthat;
		denoffset = taun_1;
	        scthat = scoffset + om*(tauhat - denoffset);
	      }
	      else {
		scoffset = scthat;
		denoffset = taun_1;
		fractrans = (realn - densamp)/((float) dentrans); 
		fractrans = 1 - ( (fractrans-1)*(fractrans-1));
	        scthat = scoffset + (omf + (om-omf)*fractrans)*(tauhat - denoffset);
	      }
	    }
	    else
	      scthat = omf*tauhat;
	    theta = atan2(scthat,1.0)+scthat;
	    if (absg < ac)
	    {
		deltheta = theta-thetan_1;
		B = 1.0/(1.0+tan(deltheta)*tan(deltheta));
		gtilde = absg;
		t1 = s*s;
		t2 = gtilde*gtilde*(1-B);
		if (t2 > t1)
		{
		    decratio *= 2.0;
		    if (decratio > MAXDECRATIO)
		    {
			fprintf(stderr, 
				"k-space calculation failed.");
			exit(1);
		    }
		    loop = 1;
		    break;
		}
		t3 = sqrt(t1-t2);
		absg = sqrt(B)*gtilde+t3;
		if (absg > ac)
		    absg = ac;
	    }
	    tgx = absg*cos(theta);
	    tgy = absg*sin(theta);
	    tkx += tgx;
	    tky += tgy;
	    thetan_1=theta;
	    if (!(n % ((int) Round(decratio))))
	    {
		m = n/Round(decratio);
		gx = ((int) Round(tkx-oldkx))/decratio;
		gx &= ~1;
		gy = ((int) Round(tky-oldky))/decratio;
		gy &= ~1;
		kxt += gx;
		kyt += gy;
		oldkx = tkx;
		oldky = tky;
		if (!(m % dr))
		{
		    m /= dr;
		    if ((m+1)>(npts-1))
			break;
		    for (i = 0; i < nl; i++)
		    {
                        kx[i][m+1] = kxscale*(kxt*r1[i] - kyt*r2[i])/distance;
			ky[i][m+1] = kyscale*(kxt*r2[i] + kyt*r1[i])/distance;
		    }
		    samplecount= m+2;
		}
	    }
	    n++;
	}
    }
    kdenrad = ceil(2*nl*om*ksv/(2*M_PI));
#ifdef never
    fprintf(stderr,"kdenread = %d\n",kdenrad);
#endif
    resolution = ceil(2*nl*om*taun/(2*M_PI)); 

    free(r1);
    free(r2);
    *res= resolution;
    return samplecount;
}    

static float *fspace(int size)
{
    float *buffer;
    if (!(buffer = (float *)calloc( size,sizeof (float) ) ))
    { fprintf(stderr,"calloc: space not assigned"); exit (0); }
    return(buffer);
}

