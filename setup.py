import os
from subprocess import check_output
import sys

# Numpy/mpi4py must be installed prior to installing TACS
import numpy
import mpi4py

# Import distutils
from setuptools import setup
from distutils.core import Extension as Ext
from Cython.Build import cythonize
from Cython.Compiler import Options

# Convert from local to absolute directories
def get_global_dir(files):
    tmr_root = os.path.abspath(os.path.dirname(__file__))
    new = []
    for f in files:
        new.append(os.path.join(tmr_root, f))
    return new

def get_mpi_flags():
    # Split the output from the mpicxx command
    args = check_output(['mpicxx', '-show']).decode('utf-8').split()

    # Determine whether the output is an include/link/lib command
    inc_dirs, lib_dirs, libs = [], [], []
    for flag in args:
        if flag[:2] == '-I':
            inc_dirs.append(flag[2:])
        elif flag[:2] == '-L':
            lib_dirs.append(flag[2:])
        elif flag[:2] == '-l':
            libs.append(flag[2:])

    return inc_dirs, lib_dirs, libs

inc_dirs, lib_dirs, libs = get_mpi_flags()

# Add the numpy/mpi4py/tacs/paropt include directories
inc_dirs.extend([numpy.get_include(), mpi4py.get_include()])

# Add stacs libraries
rel_inc_dirs = ['cpp']
rel_lib_dirs = ['cpp']
inc_dirs.extend(get_global_dir(rel_inc_dirs))
lib_dirs.extend(get_global_dir(rel_lib_dirs))
libs.extend(['stacs'])
runtime_lib_dirs = get_global_dir(['cpp'])

# Add the TACS libraries
import tacs
if 'tacs' in sys.modules:
    inc_dirs.extend(tacs.get_include())
    inc_dirs.extend(tacs.get_cython_include())
    tacs_lib_dirs, tacs_libs = tacs.get_libraries()
    lib_dirs.extend(tacs_lib_dirs)
    libs.extend(tacs_libs)
    runtime_lib_dirs.extend(tacs_lib_dirs)

# Add the PSPACE libraries
import pspace
if 'pspace' in sys.modules:
    inc_dirs.extend(pspace.get_include())
    inc_dirs.extend(pspace.get_cython_include())
    pspace_lib_dirs, pspace_libs = pspace.get_libraries()
    lib_dirs.extend(pspace_lib_dirs)
    libs.extend(pspace_libs)
    runtime_lib_dirs.extend(pspace_lib_dirs)

exts = []
mod = 'STACS'
exts.append(Ext('stacs.%s'%(mod), sources=['stacs/%s.pyx'%(mod)],
                include_dirs=inc_dirs, libraries=libs,
                library_dirs=lib_dirs, runtime_library_dirs=runtime_lib_dirs))
for e in exts:
    e.cython_directives = {'embedsignature': True,
                           'binding': True}
setup(name='stacs',
      version=1.0,
      description='Stochastic TACS : UQ/OUU plugin for TACS',
      author='Komahan Boopathy',
      author_email='komibuddy@gmail.com',
      ext_modules=cythonize(exts, language='c++',
                            include_path=inc_dirs)
      )
