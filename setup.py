import subprocess
from setuptools import setup, Extension
from distutils.command.build import build as DistutilsBuild
from distutils.command.install import install as DistutilsInstall
from pprint import pprint
from pathlib import Path
import itertools

# Thanks to https://stackoverflow.com/users/2650249/hoefling
# for https://stackoverflow.com/a/50370153/4821395
def partition(pred, iterable):
    t1, t2 = itertools.tee(iterable)
    return itertools.filterfalse(pred, t1), filter(pred, t2)


class MyInstall(DistutilsInstall):

    def run(self):
        pprint(dir(self))
        self.do_pre_install_stuff()
        DistutilsInstall.run(self)
        self.do_post_install_stuff()

    def do_pre_install_stuff(self):
        print(f'DO_PRE_INSTALL_STUFF IS HAPPENING in {self.build_base}')
        print('DO_PRE_INSTALL_STUFF FINISHED')
        
    def do_post_install_stuff(self):
        print('DO_POST_INSTALL_STUFF IS HAPPENING')
        
class MyBuild(DistutilsBuild):

    def finalize_options(self):
        """
        We need the calls to swig (triggered by build_ext) to happen
        before the addition of py_modules, because swig *creates* the
        python modules!
        """
        super().finalize_options()
        condition = lambda el: el[0] == 'build_ext'
        rest, sub_build_ext = partition(condition, self.sub_commands)
        self.sub_commands[:] = list(sub_build_ext) + list(rest)
    
    def run(self):
        pprint(dir(self))
        self.do_pre_build_stuff()
        DistutilsBuild.run(self)
        self.do_post_build_stuff()

    def do_pre_build_stuff(self):
        print(f'DO_PRE_BUILD_STUFF IS HAPPENING in {self.build_base}')
        top_dir = Path(__file__).parent
        subprocess.run([top_dir / 'configure'])
        # scan the newly created config.mk for what Fiasco's architecture name
        with open('config.mk') as f:
            for line in f:
                words = line.split()
                words = [word.strip() for word in words]
                if words[0].upper() == 'ARCH' and words[1] == '=':
                    fiasco_arch = words[2]
                    break
            else:
                raise RuntimeError('Could not find config parameter ARCH')
        subprocess.run(['make', '-f', top_dir / 'Makefile'])
        (top_dir / 'FiascoFiat').mkdir(exist_ok=True)
        include_dir = top_dir / 'FiascoFiat' / 'include'
        include_dir.mkdir(exist_ok=True)
        for path in (top_dir / 'include' / fiasco_arch).glob('**/*.h'):
            (include_dir / path.name).unlink(missing_ok=True)
            (include_dir / path.name).symlink_to(path)
        (include_dir / 'lapack.h').unlink(missing_ok=True)
        (include_dir / 'lapack.h').symlink_to(top_dir / 'src' / 'fmri' / 'lapack.h')
        quaternion_dir = top_dir / 'FiascoFiat' / 'quaternion'
        quaternion_dir.mkdir(exist_ok=True)
        for fname in ['quaternion.i', 'quaternion.c']:
            target = top_dir / 'src' / 'fmri' / fname
            (quaternion_dir / fname).unlink(missing_ok=True)
            (quaternion_dir / fname).symlink_to(target)
        
    
        print(f'DO_PRE_BUILD_STUFF FINISHED')
        
    def do_post_build_stuff(self):
        print('DO_POST_BUILD_STUFF IS HAPPENING')
        
setup(
    cmdclass={'install':MyInstall, 'build':MyBuild},
    name='FiascoFiat',
    version='5.3.1',
    #packages=['quaternion'],
    install_requires=['subprocess', 'pathlib'],
    ext_modules=[
        Extension('_quaternion',
                  ['FiascoFiat/quaternion/quaternion.i',
                   'FiascoFiat/quaternion/quaternion.c'],
                  include_dirs=['FiascoFiat/include'],
                  swig_opts=['-IFiascoFiat/include'],
                  libraries=['lapack', 'blas', 'misc'],
                  library_dirs=['lib/LINUXX86_64'],
                  define_macros=[('FORTRAN_ADD_UNDERSCORE', None)],
#                  extra_link_args=['-llapack', '-lblas'],
        )
        ],
    py_modules=['FiascoFiat.quaternion.quaternion'],
)
