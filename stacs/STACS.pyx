# distutils: language = c++
from STACS cimport *

# Import numpy
cimport numpy as np
import numpy as np
np.import_array()

# Import C methods for python
from cpython cimport PyObject, Py_INCREF, Py_DECREF

include "TacsDefs.pxi"

cdef inplace_array_1d(int nptype, int dim1, void *data_ptr):
    '''Return a numpy version of the array'''
    cdef int size = 1
    cdef np.npy_intp shape[1]
    cdef np.ndarray ndarray
    shape[0] = <np.npy_intp>dim1
    ndarray = np.PyArray_SimpleNewFromData(size, shape, nptype, data_ptr)
    return ndarray

cdef void updateCB(TACSElement *elem, TacsScalar *yvals, void *pyptr):
    _yvals = inplace_array_1d(TACS_NPY_SCALAR, 5, <void*> yvals)
    (<object>pyptr).update(_yvals)
    return

cdef class PyStochasticElement(Element):
    def __cinit__(self, Element elem,
                  PyParameterContainer pc,
                  update):
        self.sptr = new TACSStochasticElement(elem.ptr, pc.ptr, &updateCB)
        self.sptr.incref()        
        self.sptr.setPythonCallback(<PyObject*>update)
        self.ptr = self.sptr
        Py_INCREF(update)
        Py_INCREF(self)
    def __dealloc__(self):        
        if self.sptr:
            self.sptr.decref()
            Py_DECREF(self)
    def getDeterministicElement(self):
        delem = Element()
        delem.ptr = self.sptr.getDeterministicElement() 
        delem.ptr.incref()
        return delem
    def updateElement(self, Element elem, np.ndarray[TacsScalar, ndim=1, mode='c'] vals):
        self.sptr.updateElement(elem.ptr, <TacsScalar*> vals.data)
    def setPythonCallback(self, cb):
        self.sptr.setPythonCallback(<PyObject*>cb)

cdef class MutableElement3D(Element):
    def __cinit__(self, Element elem):
        self.sptr = new TACSMutableElement3D(elem.ptr)
        self.ptr = self.sptr
        self.ptr.incref()
        return
    def __dealloc__(self):        
        if self.ptr:
            self.ptr.decref()
    def setDensity(self, int rho):
        self.sptr.setDensity(rho)

cdef class PyMomentSpaceTimeIntegral(Function):
    def __cinit__(self,
                  Assembler assembler, Function func, PyParameterContainer pc,
                  int quantity_type, int moment_type):
        self.sptr = new TACSStochasticFunction( assembler.ptr, func.ptr, pc.ptr,
                                                quantity_type, moment_type )
        self.sptr.incref()        
        self.ptr = self.sptr
        Py_INCREF(self) #? do we need this?
        return
    
    def __dealloc__(self):        
        if self.sptr:
            self.sptr.decref()
            Py_DECREF(self)  #? do we need this?
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
        Py_INCREF(self) #? do we need this?
        return
    
    def __dealloc__(self):        
        if self.sptr:
            self.sptr.decref()
            Py_DECREF(self)  #? do we need this?
        return

    def getFunctionValue(self):
        return self.sptr.getFunctionValue()
