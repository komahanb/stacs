# distutils: language = c++
from STACS cimport *

# Import numpy
cimport numpy as np
import numpy as np
np.import_array()

# Import C methods for python
from cpython.object cimport PyObject

include "TacsDefs.pxi"

cdef inplace_array_1d(int nptype, int dim1, void *data_ptr):
    '''Return a numpy version of the array'''
    cdef int size = 1
    cdef np.npy_intp shape[1]
    cdef np.ndarray ndarray
    shape[0] = <np.npy_intp>dim1
    ndarray = np.PyArray_SimpleNewFromData(size, shape, nptype, data_ptr)
    return ndarray

cdef class PyStochasticElement(Element):
    def __cinit__(self, Element elem,
                  PyParameterContainer pc,
                  update):
        self.sptr = new TACSStochasticElement(elem.ptr, pc.ptr, NULL)
        self.sptr.incref()        
        self.ptr = self.sptr
    def __dealloc__(self):        
        if self.sptr:
            self.sptr.decref()
    def getDeterministicElement(self):
        delem = Element()
        delem.ptr = self.sptr.getDeterministicElement() 
        delem.ptr.incref()
        return delem
    def updateElement(self, Element elem, np.ndarray[TacsScalar, ndim=1, mode='c'] vals):
        self.sptr.updateElement(elem.ptr, <TacsScalar*> vals.data)
    def setPythonCallback(self, cb):
        self.sptr.setPythonCallback(<PyObject*>cb)

cdef class PyMomentSpaceTimeIntegral(Function):
    def __cinit__(self,
                  Assembler assembler, Function func, PyParameterContainer pc,
                  int quantity_type, int moment_type):
        self.sptr = new TACSStochasticFunction( assembler.ptr, func.ptr, pc.ptr,
                                                quantity_type, moment_type )
        self.sptr.incref()        
        self.ptr = self.sptr
        return
    
    def __dealloc__(self):        
        if self.sptr:
            self.sptr.decref()
        return

    def getFunctionValue(self):
        return self.sptr.getFunctionValue()

cdef class PyMomentMaxSpaceTimeIntegral(Function):
    def __cinit__(self,
                  Assembler assembler, Function func, PyParameterContainer pc,
                  int quantity_type, int moment_type,
                  int ksweight):
        self.sptr = new TACSKSStochasticFunction( assembler.ptr, func.ptr, pc.ptr,
                                                  quantity_type, moment_type,
                                                  ksweight)
        self.sptr.incref()        
        self.ptr = self.sptr
        return
    
    def __dealloc__(self):        
        if self.sptr:
            self.sptr.decref()
        return

    def getFunctionValue(self):
        return self.sptr.getFunctionValue()
