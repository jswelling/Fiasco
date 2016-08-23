#include <stdio.h>
#include <stdlib.h>
#include <fmri.h>

static const char* knownSlicePatternNames[]= {
  "sequential",
  "reversed_sequential",
  "even/odd",
  "reversed_even/odd",
  "odd/even",
  "reversed_odd/even",
  "halves_low_first",
  "reversed_halves_low_first",
  "halves_high_first",
  "reversed_halves_high_first",
  NULL /* must be last */
};

int lengthsToTry[]= { 26,27,1,2 };
#define NLENGTHS (sizeof(lengthsToTry)/sizeof(int))

void shuffle( int len, int* array, int* scratch, const char* pattern )
{
  int* tbl= slp_generateSlicePatternTable( len, pattern );
  int i;
  printf("Applying shuffle %s\n",pattern);
  for (i=0; i<len; i++) scratch[i]= array[tbl[i]];
  for (i=0; i<len; i++) array[i]= scratch[i];
  free(tbl);
}

void unshuffle( int len, int* array, int* scratch, const char* pattern )
{
  int* tbl= slp_generateInvertedSlicePatternTable( len, pattern );
  int i;
  printf("Applying unshuffle %s\n",pattern);
  for (i=0; i<len; i++) scratch[i]= array[tbl[i]];
  for (i=0; i<len; i++) array[i]= scratch[i];
  free(tbl);
}

int* mkScratch(int len)
{
  int i;
  int* result= (int*)malloc(len*sizeof(int));
  if (!result) {
    fprintf(stderr,"Unable to allocate %d bytes!\n",len*sizeof(int));
    exit(-1);
  }
  return result;
}

int* mkArray(int len)
{
  int i;
  int* result= mkScratch(len);
  for (i=0;i<len;i++) result[i]= i;
  return result;
}

int main(int argc, char* argv[])
{
  int i;

  for (i=0; i<NLENGTHS; i++) {
    int len= lengthsToTry[i];
    int* array= mkArray(len);
    int* scratch= mkScratch(len);
    int j;
    printf("#### Trying length %d\n",len);
    for (j=0; knownSlicePatternNames[j]; j++)
      shuffle(len,array,scratch,knownSlicePatternNames[j]);
    j -= 1;
    for ( ; j>=0; j--)
      unshuffle(len,array,scratch,knownSlicePatternNames[j]);
#ifdef never
#endif
    for (j=0; j<len; j++) printf("%d: %d\n",j,array[j]);
    free(array);
    free(scratch);
  }
}
