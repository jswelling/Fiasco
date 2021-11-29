CFLAGS="-D_BSD_SOURCE -D_XOPEN_SOURCE=500 -DFORTRAN_ADD_UNDERSCORE -Wimplicit -DLINUX -O -fPIC -DLINUXX86_64 -I/home/welling/git/Fiasco//include/LINUXX86_64 -DMPI -DFFTW3 -DUSE_PNG -DNOAFS -I/home/welling/anaconda3/envs/fiascoEnv/include/python3.9"
LFLAGS="-L/home/welling/git/Fiasco//lib/LINUXX86_64 -lpng -lfftw3 -lm -lfmri -lmri -lpar -lbio -lacct -lmisc -lcrg -llapack -lblas -lm"
