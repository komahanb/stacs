# Typdefs required for either real or complex mode
include "TacsTypedefs.pxi"

# Include the mpi4py header
cdef extern from "mpi-compat.h":
    pass

from pspace.PSPACE cimport *
from tacs.elements cimport *
from tacs.functions cimport *

# C++ class TACSStochasticElement
cdef extern from "TACSStochasticElement.h":
    cdef cppclass TACSStochasticElement(TACSElement):
        TACSStochasticElement( TACSElement *_delem,
                               ParameterContainer *_pc,
                               void (*_update)(TACSElement*, TacsScalar*, void*) except * )
        TACSElement* getDeterministicElement()
        void updateElement(TACSElement* elem, TacsScalar* vals)
        void setPythonCallback(PyObject *cbptr)

# Python wrapped class to C++ class TACSStochasticElement
cdef class PyStochasticElement(Element):
    cdef TACSStochasticElement *sptr

# C++ class     
cdef extern from "TACSStochasticFunction.h":
    cdef cppclass TACSStochasticFunction(TACSFunction):
        TACSStochasticFunction( TACSAssembler *tacs,
                                TACSFunction *dfunc,
                                ParameterContainer *pc,
                                int quantityType,
                                int moment_type )
        TacsScalar getFunctionValue()
        
# Python wrapped class
cdef class PyMomentSpaceTimeIntegral(Function):
    cdef TACSStochasticFunction *sptr

# C++ class 
cdef extern from "TACSKSStochasticFunction.h":
    cdef cppclass TACSKSStochasticFunction(TACSFunction):
        TACSKSStochasticFunction( TACSAssembler *tacs,
                                  TACSFunction *dfunc,
                                  ParameterContainer *pc,
                                  int quantityType,
                                  int moment_type,
                                  int ksweight)
        TacsScalar getFunctionValue()

# Python wrapped class
cdef class PyMomentMaxSpaceTimeIntegral(Function):
    cdef TACSKSStochasticFunction *sptr
