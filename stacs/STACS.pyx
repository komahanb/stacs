# distutils: language = c++
from __future__ import print_function, division

# Import TACS
# include "TacsDefs.pxi"

# Import C methods for python
from cpython cimport PyObject, Py_INCREF

## cdef class PyStochasticElement(Element):
##     cdef TACSStochasticElement *sptr
##     def __cinit__(self, Element elem):#, PyParameterContainer pc):
##         self.sptr = new TACSStochasticElement(elem.ptr, NULL, NULL)
##         self.ptr = self.sptr
##         self.ptr.incref()
##         return

##     def __dealloc__(self):
##         if self.ptr:
##             self.ptr.decref()
##         return
