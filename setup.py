import os
from subprocess import check_output
import sys

# Numpy/mpi4py must be installed prior to installing TACSUQ
import numpy
import mpi4py

# Import distutils
from setuptools import setup
from distutils.core import Extension as Ext
from Cython.Build import cythonize

# Convert from local to absolute directories
def get_global_dir(files):
    tacsuq_root = os.path.abspath(os.path.dirname(__file__))
    new = []
    for f in files:
        new.append(os.path.join(tacsuq_root, f))
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

# Add tacsuq-dev/lib as a runtime directory
runtime_lib_dirs = get_global_dir(['stacs'])

# Relative paths for the include/library directories
rel_inc_dirs = ['cpp']
rel_lib_dirs = ['cpp']
libs.extend(['stacs'])

# Convert from relative to absolute directories
inc_dirs.extend(get_global_dir(rel_inc_dirs))
lib_dirs.extend(get_global_dir(rel_lib_dirs))

# include TACS
inc_dirs.extend(get_global_dir(['../tacs/src/',
                                '../tacs/src/elements/',
                                '../tacs/src/functions/']))
# include PSPACE
inc_dirs.extend(get_global_dir(['../pspace/cpp/']))

# Add the numpy/mpi4py directories
inc_dirs.extend([numpy.get_include(), mpi4py.get_include()])

exts = []
for mod in ['cwrap']:
    exts.append(Ext('stacs.%s'%(mod), sources=['stacs/%s.pyx'%(mod)],
                    include_dirs=inc_dirs, libraries=libs, 
                    library_dirs=lib_dirs, runtime_library_dirs=runtime_lib_dirs,
                    cython_directives={"embedsignature": True, "binding": True}))

setup(name='stacs',
      version=1.0,
      description='Stochastic TACS : UQ/OUU plugin for TACS',
      author='Komahan Boopathy',
      author_email='komibuddy@gmail.com',
      ext_modules=cythonize(exts, include_path=inc_dirs))