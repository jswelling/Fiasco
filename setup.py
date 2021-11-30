import subprocess
import sys
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


def link_in(pathstr, fiasco_fiat_dir):
    top_dir = fiasco_fiat_dir.parent
    target = top_dir / pathstr
    (fiasco_fiat_dir / target.name).unlink(missing_ok=True)
    (fiasco_fiat_dir / target.name).symlink_to(target)
        
        
class MyInstall(DistutilsInstall):

    def run(self):
        self.do_pre_install_stuff()
        DistutilsInstall.run(self)
        self.do_post_install_stuff()

    def do_pre_install_stuff(self):
        pass
        
    def do_post_install_stuff(self):
        pass


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
        self.do_pre_build_stuff()
        DistutilsBuild.run(self)
        self.do_post_build_stuff()

    def do_pre_build_stuff(self):
        top_dir = Path(__file__).parent
        subprocess.run([top_dir / 'configure'])
        subprocess.run(['make', '-f', top_dir / 'Makefile'])
        fiasco_fiat_dir = top_dir / 'FiascoFiat'
        fiasco_fiat_dir.mkdir(exist_ok=True)
        (fiasco_fiat_dir / '__init__.py').touch()
        include_dir = fiasco_fiat_dir / 'include'
        include_dir.mkdir(exist_ok=True)
        fiasco_arch = read_config_dict(required=['ARCH'])['ARCH']
        for path in (top_dir / 'include' / fiasco_arch).glob('**/*.h'):
            (include_dir / path.name).unlink(missing_ok=True)
            (include_dir / path.name).symlink_to(path)
        (include_dir / 'lapack.h').unlink(missing_ok=True)
        (include_dir / 'lapack.h').symlink_to(top_dir / 'src' / 'fmri' / 'lapack.h')
        fiasco_fiat_dir.mkdir(exist_ok=True)
        fname_list = [
            'src/fmri/quaternion.i', 'src/fmri/quaternion.c',
            'src/csh/fiasco_utils.py',
            'src/fmri/fiasco_numpy.i', 'src/fmri/numpy.i',
            'src/mripipes/mripipes.i', 'src/mripipes/mripipes.c',
        ]
        fname_list.extend([str(path)
                           for path in
                           (top_dir/'src'/'mripipes').glob('*_tool.c')])
        for fname in fname_list:
            link_in(fname, fiasco_fiat_dir)        
    
    def do_post_build_stuff(self):
        pass

def read_config_dict(required=[]):
    """
    Return the contents of the config.mk file generated by the
    native configuration routine in the form of a dict
    """
    top_dir = Path(__file__).parent
    subprocess.run([top_dir / 'configure'])
    assert Path('config.mk').exists(), 'config.mk was not found'
    # scan the newly created config.mk for what Fiasco's architecture name
    rslt = {}
    with open('config.mk') as f:
        for line in f:
            words = line.split()
            words = [word.strip() for word in words]
            if len(words) >= 3 and words[1] == '=' and words[0][0] != '#':
                rslt[words[0].upper()] = ' '.join(words[2:])
    for req in required:
        if req.upper() not in rslt:
            raise RuntimeError(f'Could not find config parameter {req}')
    return rslt


def parse_link_args_for_libs(arg_str):
    """
    Given an arg_str like "-L/baz -L/blrfl -lfoo -lbar", return
    a tuple of lists like (["foo", "bar"], ["/baz", "/blrfl'])
    """
    lib_l = []
    libdir_l = []
    for word in arg_str.split():
        word = word.strip()
        if word.startswith('-l'):
            lib_l.append(word[2:])
        elif word.startswith('-L'):
            libdir_l.append(word[2:])
        else:
            raise RuntimeError(f'Unexpected word {word} in arg string {arg_str}')
    return lib_l, libdir_l


def get_extra_macros(config):
    rslt = []
    if config['ARCH'] not in ["HPPA", "HPPA20", "CRAY", "T3D", "T3E"]:
        rslt.append(('FORTRAN_ADD_UNDERSCORE', None))
    return rslt


def gen_quaternion_extension():
    config = read_config_dict(required=['ARCH', 'LAPACK_LIBS'])
    libraries = ['misc']
    library_dirs = [f"lib/{config['ARCH']}"]
    libs, dirs = parse_link_args_for_libs(config['LAPACK_LIBS'])
    libraries.extend(libs)
    library_dirs.extend(dirs)
    return Extension(
        '_quaternion',
        ['FiascoFiat/quaternion.i', 
         'FiascoFiat/quaternion.c'],
        include_dirs=['FiascoFiat/include'],
        swig_opts=['-IFiascoFiat/include'],
        libraries = libraries,
        library_dirs= library_dirs,
        define_macros=get_extra_macros(config)
    )


def gen_quaternion_pymodule():
    return 'FiascoFiat.quaternion'


def gen_fiasco_utils_pymodule():
    return 'FiascoFiat.fiasco_utils'


def gen_fiasco_numpy_extension():
    config = read_config_dict(required=['ARCH', 'LAPACK_LIBS'])
    libraries = ['fmri', 'misc']
    library_dirs = [f"lib/{config['ARCH']}"]
    libs, dirs = parse_link_args_for_libs(config['LAPACK_LIBS'])
    libraries.extend(libs)
    library_dirs.extend(dirs)
    # Obtain the numpy include directory.
    # This logic works across numpy versions.
    try:
        import numpy
        numpy_include = numpy.get_include()
    except ImportError:
        sys.exit("The python package 'numpy' is not installed!")
    except AttributeError:
        numpy_include = numpy.get_numpy_include()
    
    return Extension(
        '_fiasco_numpy',
        ['FiascoFiat/fiasco_numpy.i'],
        include_dirs=['FiascoFiat/include', numpy_include],
        swig_opts=['-IFiascoFiat/include', f'-I{numpy_include}'],
        libraries = libraries,
        library_dirs= library_dirs,
        define_macros=get_extra_macros(config)
    )
    

def gen_fiasco_numpy_pymodule():
    return 'FiascoFiat.fiasco_numpy'


def gen_mripipes_extension():
    config = read_config_dict(required=['ARCH', 'LAPACK_LIBS'])
    src_files = ['FiascoFiat/mripipes.i', 'FiascoFiat/mripipes.c']
    libraries = ['mripipestools', 'mri', 'fmri', 'dcdf', 'misc', 'bio', 'crg']
    library_dirs = [f"lib/{config['ARCH']}"]
    libs, dirs = parse_link_args_for_libs(config['LAPACK_LIBS'])
    libraries.extend(libs)
    library_dirs.extend(dirs)
    
    return Extension(
        '_mripipes',
        src_files,
        include_dirs=['FiascoFiat/include'],
        swig_opts=['-IFiascoFiat/include'],
        libraries = libraries,
        library_dirs= library_dirs,
        define_macros=get_extra_macros(config)
    )


def gen_mripipes_pymodule():
    return 'FiascoFiat.mripipes'


setup(
    cmdclass={'install':MyInstall, 'build':MyBuild},
    name='FiascoFiat',
    version='5.3.1',
    install_requires=['subprocess', 'pathlib', 'numpy'],
    ext_modules=[
        gen_quaternion_extension(),
        gen_fiasco_numpy_extension(),
        gen_mripipes_extension(),
    ],
    py_modules=[
        gen_quaternion_pymodule(),
        gen_fiasco_numpy_pymodule(),
        gen_fiasco_utils_pymodule(),
        gen_mripipes_pymodule()
    ],
)
