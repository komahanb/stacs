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
from tacs.elements cimport *
from STACS cimport *


# Import TACS
include "TacsDefs.pxi"

# Import the definition required for const strings
from libc.string cimport const_char
from libc.stdlib cimport malloc, free

# Import C methods for python
from cpython cimport PyObject, Py_INCREF

cdef class PyStochasticElement(Element):
    def __cinit__(self, Element elem):
        cdef TACSStochasticElement *sptr
        sptr = new TACSStochasticElement(elem.ptr, NULL, NULL)
        #self.ptr = sptr
        #self.ptr.incref()
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()
        return
