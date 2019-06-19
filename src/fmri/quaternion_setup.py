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

# print('incDirList: ', incDirList)
# print('libDirList: ', libDirList)
# print('libList: ', libList)
# print('defList: ', defList)
# print('undefList: ', undefList)
# print('otherCompileFlags: ', otherCompileFlags)
# print('otherLinkFlags: ', otherLinkFlags)

setup(name="quaternion",
      version='1.0',
      py_modules=['quaternion'],
      ext_modules=[Extension('_quaternion',
                             ['quaternion.c','quaternion_wrap.c'],
                             define_macros=defList,
                             undef_macros=undefList,
                             include_dirs=incDirList,
                             extra_compile_args=otherCompileFlags,
                             library_dirs=libDirList,
                             libraries=libList,
                             extra_link_args=otherLinkFlags
                             )
                   ],
      description='Python quaternion functionality for Fiasco/FIAT',
      author='Joel Welling',
      author_email='welling@psc.edu',
      url='http://www.stat.cmu.edu/~fiasco',
      )
