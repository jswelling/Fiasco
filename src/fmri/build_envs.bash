CFLAGS="-D_BSD_SOURCE -D_XOPEN_SOURCE=500 -DFORTRAN_ADD_UNDERSCORE -Wimplicit -DLINUX -O -fPIC -DLINUXX86_64 -I/home/welling/git/Fiasco//include/LINUXX86_64 -DMPI -DFFTW3 -DNOAFS -I/home/welling/anaconda3/include/python3.6"
LFLAGS="-L/home/welling/git/Fiasco//lib/LINUXX86_64 -lfftw3 -lm -lfmri -lmri -lpar -lbio -lacct -lmisc -lcrg -llapack -lblas -lm"
