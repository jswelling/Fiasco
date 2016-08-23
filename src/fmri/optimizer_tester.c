/************************************************************
 *                                                          *
 *  optimizer_tester.c                                      *
 *                                                          *
 *  Permission is hereby granted to any individual or       *
 *  institution for use, copying, or redistribution of      *
 *  this code and associated documentation, provided        *
 *  that such code and documentation are not sold for       *
 *  profit and the following copyright notice is retained   *
 *  in the code and documentation:                          *
 *                                                          *
 *  Copyright (c) 1999 Department of Statistics             *
 *                     Carnegie Mellon University           *
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
 *  Original programming by Joel Welling 2/2004             *
 ************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fmri.h>

typedef struct Data {
  double A;
  double B;
  double C;
  double D;
} Data;

static double value( ScalarFunction* sf, const double* par, const int nPar )
{
  double result;
  Data* data= (Data*)(sf->data);
  if (nPar != 3)
    Abort("Wrong number of parameters in optimizer_tester:value()!\n");
  result= (par[0]-data->A)*(par[0]-data->A)
    + (par[1]-data->B)*(par[1]-data->B)
    + (par[2]-data->C)*(par[2]-data->C)
    + data->D;
  fprintf(stderr,"eval! %g %g %g -> %g\n",par[0],par[1],par[2],result);

  return result;
}

static void reset( ScalarFunction* sf, const double* par, const int nPar )
{
  Data* data= (Data*)(sf->data);
  if (nPar != 3)
    Abort("Wrong number of parameters in optimizer_tester:value()!\n");
  fprintf(stderr,"reset! %g %g %g\n",par[0],par[1],par[2]);
}

static void destroy( ScalarFunction* sf )
{
  free(sf);
}

int main( int argc, char* argv[] ) 
{
  ScalarFunction* sf= NULL;
  Optimizer* opt= NULL;
  Data data;
  char inbuf[81];
  double par[3];
  int i;
  double best;
  int code;

  if (!(sf=(ScalarFunction*)malloc(sizeof(ScalarFunction*))))
    Abort("%s: unable to allocate %d bytes!\n",argv[0]);
  sf->value= value;
  sf->reset= reset;
  sf->destroySelf= destroy;
  sf->data= &data;
  sf->nPar= 3;

  data.A= 3.0;
  data.B= 5.0;
  data.C= -1.0;
  data.D= 7.0;

  printf("Enter the character encoding of an optimizer\n");
  fgets(inbuf,80,stdin);
  inbuf[strlen(inbuf)-1]= '\0'; /* strip trailing newline */
  opt= optimizerFromStringRep( inbuf );
  fprintf(stderr,"Stringrep is <%s> vs. your input <%s>\n",
	  opt->getStringRep(opt),inbuf);
  opt->setDebugLevel(opt,4);
  fprintf(stderr,"debugLevel set to %d\n",opt->getDebugLevel(opt));
  for (i=0; i<sizeof(par)/sizeof(double); i++) par[i]= 0.0;
  code= opt->go(opt, sf, par, sizeof(par)/sizeof(double), &best);
  fprintf(stderr,"Result code: %d\n",code);
  for (i=0; i<sizeof(par)/sizeof(double); i++)
    fprintf(stderr,"%d: %f\n",i,par[i]);
  fprintf(stderr,"Minimum value is %f\n",best);
  opt->destroySelf(opt);
}
