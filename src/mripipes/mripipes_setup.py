from __future__ import print_function
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

if 'SRCFILES' not in os.environ:
    sys.exit("Required environment variable SRCFILES is not set!")
srcFileList= os.environ['SRCFILES'].split()

## print("incDirList: <%s>"%incDirList)
## print("libDirList: <%s>"%libDirList)
## print("libList: <%s>"%libList)
## print("defList: <%s>"%defList)
## print("undefList: <%s>"%undefList)
## print("otherCompileFlags: <%s>"%otherCompileFlags)
## print("otherLinkFlags: <%s>"%otherLinkFlags)
print("srcFileList: <%s>"%srcFileList)

setup(name="mripipes",
      version='1.0',
      py_modules=['mripipes'],
      ext_modules=[Extension('_mripipes',srcFileList,
                             define_macros=defList,
                             undef_macros=undefList,
                             include_dirs=incDirList,
                             extra_compile_args=otherCompileFlags,
                             library_dirs=libDirList,
                             libraries=libList,
                             extra_link_args=otherLinkFlags
                             )
                   ],
      description='Python mripipes functionality for Fiasco/FIAT',
      author='Joel Welling',
      author_email='welling@psc.edu',
      url='http://www.stat.cmu.edu/~fiasco',
      )
