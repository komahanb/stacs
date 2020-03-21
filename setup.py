import os
from subprocess import check_output
import sys

# Numpy/mpi4py must be installed prior to installing TACS
import numpy
import mpi4py
import pspace
import tacs
import tmr

# Import distutils
from setuptools import setup
from distutils.core import Extension as Ext
from Cython.Build import cythonize

# Convert from local to absolute directories
def get_global_dir(files):
    tacs_root = os.path.abspath(os.path.dirname(__file__))
    new = []
    for f in files:
        new.append(os.path.join(tacs_root, f))
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

# Add the numpy/mpi4py directories
inc_dirs.extend([numpy.get_include(), mpi4py.get_include()])

# STACS (local cpp files)
inc_dirs.extend(get_global_dir(['stacs']))
inc_dirs.extend(get_global_dir(['src/include']))
lib_dirs.extend(get_global_dir(['lib']))
libs.extend(['stacs'])

# TACS
inc_dirs.extend(tacs.get_include())
inc_dirs.extend(tacs.get_cython_include())
lib_dirs.extend(tacs.get_libraries()[0])
libs.extend(tacs.get_libraries()[1])

# TMR
inc_dirs.extend(tmr.get_include())
inc_dirs.extend(tmr.get_cython_include())
lib_dirs.extend(tmr.get_libraries()[0])
libs.extend(tmr.get_libraries()[1])

# PSPACE
inc_dirs.extend(pspace.get_include())
inc_dirs.extend(pspace.get_cython_include())
lib_dirs.extend(pspace.get_libraries()[0])
libs.extend(pspace.get_libraries()[1])

# The provide where the run time libraries are present
runtime_lib_dirs = [] 
runtime_lib_dirs.extend(lib_dirs)

exts = []
for mod in ['STACS']:
    exts.append(Ext('stacs.%s'%(mod), sources=['stacs/%s.pyx'%(mod)],
                    include_dirs=inc_dirs,
                    libraries=libs, 
                    library_dirs=lib_dirs,
                    runtime_library_dirs=runtime_lib_dirs,
                    cython_directives={"embedsignature": True, "binding": True}))

setup(name='stacs',
      version=1.0,
      description='Stochastic TACS -- Extension of TACS',
      author='Komahan Boopathy',
      author_email='komibuddy@gmail.com',
      ext_modules=cythonize(exts, include_path=inc_dirs))
