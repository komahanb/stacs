# distutils: language = c++
from __future__ import print_function, division

# Include the mpi4py header
cdef extern from "mpi-compat.h":
    pass

from mpi4py.libmpi cimport *
cimport mpi4py.MPI as MPI

# Import numpy
cimport numpy as np
import numpy as np
np.import_array()

#from TACS cimport *
from pspace.PSPACE cimport *
from tacs.TACS cimport *
from STACS cimport *

# Import TACS
include "TacsDefs.pxi"

# Import the definition required for const strings
from libc.string cimport const_char
from libc.stdlib cimport malloc, free

# Import C methods for python
from cpython cimport PyObject, Py_INCREF

## cdef class PyStochasticElement(elements.Element):
##     def __cinit__(self,
##                   elements.Element element, PyParameterContainer pcontainer,
##                   NULL):
##         self.ptr = None
##         self.ptr.incref()        
##         return

#cdef class PyGreeting:
#    def __cinit__(self, ):
#        return
