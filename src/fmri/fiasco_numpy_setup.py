import sys
import os
from distutils.core import setup, Extension

incDirList= []
libDirList= []
libList= []
defList= []
undefList= []
otherCompileFlags= []
otherLinkFlags= []

if 'CFLAGS' in os.environ:
    cflags= os.environ['CFLAGS']
else:
    cflags= ''

words= cflags.split()
for word in words:
    if word.startswith('-D'):
        word= word[2:]
        loc= word.find('=')
        if (loc>=0):
            defList.append((word[:loc-1],word[loc+1:]))
        else:
            defList.append((word,None))
    elif word.startswith('-U'):
        undefList.append(word[2:])
    elif word.startswith('-I'):
        incDirList.append(word[2:])
    else:
        otherCompileFlags.append(word)

if 'LFLAGS' in os.environ:
    lflags= os.environ['LFLAGS']
else:
    lflags= ''

words= lflags.split()
for word in words:
    if word.startswith('-L'):
        libDirList.append(word[2:])
    elif word.startswith('-l'):
        libList.append(word[2:])
    else:
        otherLinkFlags.append(word)

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    import numpy
    numpy_include = numpy.get_include()
except ImportError:
    sys.exit("The python package 'numpy' is not installed!")
except AttributeError:
    numpy_include = numpy.get_numpy_include()
incDirList.append(numpy_include)

setup(name="fiasco_numpy",
      version='1.0',
      py_modules=['fiasco_numpy'],
      ext_modules=[Extension('_fiasco_numpy',
                             ['glm.c','glm_irls.c',
                              'optimizer.c','interpolator.c',
                              'kalmanfilter.c',
                              'fiasco_numpy_wrap.c'],
                             define_macros=defList,
                             undef_macros=undefList,
                             include_dirs=incDirList,
                             extra_compile_args=otherCompileFlags,
                             library_dirs=libDirList,
                             libraries=libList,
                             extra_link_args=otherLinkFlags
                             )
                   ],
      description='Python numpy interface for Fiasco/FIAT libs',
      author='Joel Welling',
      author_email='welling@psc.edu',
      url='http://www.stat.cmu.edu/~fiasco',
      )
