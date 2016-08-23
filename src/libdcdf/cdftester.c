#include <stdio.h>
#include <stdlib.h>
#include <dcdflib.h>

int main(int argc, char* argv[]) {
  double df;
  double p, pout;
  double q, qout;
  double x, xout;
  int retval= 0;
  int one= 1;
  int two= 2;
  double t1=0.0, t2=0.0, t3=0.0;
  int status;
  double bound;
  int i;

  if (argc<2) {
    printf("usage: %s x p df [t1 [t2 [t3]]]\n",argv[0]);
    return -1;
  }

  if (argc>=2) x= atof(argv[1]);
  if (argc>=3) p= atof(argv[2]);
  if (argc>=4) df= atof(argv[3]);
  if (argc>=5) t1= atof(argv[4]);
  if (argc>=6) t2= atof(argv[5]);
  if (argc>=7) t3= atof(argv[6]);

  q= 1.0-p;

#ifdef never
  /* Test normal distribution */
  printf("normal dist\n");

  cdf_nor(&one, &pout, &qout, &x, &t1, &t2, &status, &bound);
  printf("%ld %ld %g: ??? %f %f %f -> %f\n",retval,status,bound,x,t1,t2,pout);
   
  cdf_nor(&two, &p, &q, &xout, &t1, &t2, &status, &bound);
  printf("%ld %ld %g: %f ??? %f %f -> %f\n",retval,status,bound,p,t1,t2,xout);
#endif
   
#ifdef never
  /* Test chi squared distribution */
  printf("chi squared dist\n");

  cdf_chi(&one, &pout, &qout, &x, &df, &status, &bound);
  printf("%ld %ld %g: ??? %f %f -> %f\n",retval,status,bound,x,df,pout);
   
  cdf_chi(&two, &p, &q, &xout, &df, &status, &bound);
  printf("%ld %ld %g: %f ??? %f -> %f\n",retval,status,bound,p,df,xout);
#endif

#ifdef never
  /* Test F distribution */
  printf("F dist\n");

  cdf_f(&one, &pout, &qout, &x, &df, &t1, &status, &bound);
  printf("%ld %ld %g: ??? %f %f %f -> %f\n",retval,status,bound,x,df,t1,pout);
   
  cdf_f(&two, &p, &q, &xout, &df, &t1, &status, &bound);
  printf("%ld %ld %g: %f ??? %f %f -> %f\n",retval,status,bound,p,df,t1,xout);
#endif

#ifdef never
  /* Test T distribution */
  printf("T dist\n");

  cdf_t(&one, &pout, &qout, &x, &df, &status, &bound);
  printf("%ld %g: ??? %f %f -> %f\n",status,bound,x,df,pout);
   
  cdf_t(&two, &p, &q, &xout, &df, &status, &bound);
  printf("%ld %g: %f ??? %f -> %f\n",status,bound,p,df,xout);
#endif

#ifdef never
  /* Test Poisson distribution */
  printf("Poisson dist\n");

  cdf_poi(&one, &pout, &qout, &x, &df, &status, &bound);
  printf("%ld %ld %g: ??? %f %f -> %f\n",retval,status,bound,x,df,pout);
   
  cdf_poi(&two, &p, &q, &xout, &df, &status, &bound);
  printf("%ld %ld %g: %f ??? %f -> %f\n",retval,status,bound,p,df,xout);
#endif

#ifdef never
  /* Test binomial distribution */
  printf("Binomial dist\n");

  cdf_bin(&one, &pout, &qout, &x, &df, &t1, &status, &bound);
  printf("%ld %ld %g: ??? %f %f %f -> %f\n",retval,status,bound,x,df,t1,pout);
   
  cdf_bin(&two, &p, &q, &xout, &df, &t1, &status, &bound);
  printf("%ld %ld %g: %f ??? %f %f -> %f\n",retval,status,bound,p,df,t1,xout);
#endif

#ifdef never
  /* Test beta distribution */
  printf("Beta dist\n");

  cdf_bet(&one, &pout, &qout, &x, &df, &t1, &status, &bound);
  printf("%ld %ld %g: ??? %f %f %f -> %f\n",retval,status,bound,x,df,t1,pout);
   
  cdf_bet(&two, &p, &q, &xout, &df, &t1, &status, &bound);
  printf("%ld %ld %g: %f ??? %f %f -> %f\n",retval,status,bound,p,df,t1,xout);
#endif

#ifdef never
  /* Test gamma distribution */
  printf("Gamma dist\n");

  cdf_gam(&one, &pout, &qout, &x, &df, &t1, &status, &bound);
  printf("%ld %ld %g: ??? %f %f %f -> %f\n",retval,status,bound,x,df,t1,pout);
   
  cdf_gam(&two, &p, &q, &xout, &df, &t1, &status, &bound);
  printf("%ld %ld %g: %f ??? %f %f -> %f\n",retval,status,bound,p,df,t1,xout);
#endif

  for (i=1; i<=10; i++) {
    printf("cdf_ipmpar[%d]= %d\n",i,cdf_ipmpar(&i));
  }

  return 0;
}
