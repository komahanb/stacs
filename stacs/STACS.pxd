# distutils: language = c++

# Typdefs required for either real or complex mode
include "TacsTypedefs.pxi"

# Include the mpi4py header
cdef extern from "mpi-compat.h":
    pass

#from TACS cimport *
from tacs.elements cimport *
from pspace.PSPACE cimport *

cdef extern from "TACSStochasticElement.h":
    cdef cppclass TACSStochasticElement(TACSElement):
        TACSStochasticElement( TACSElement *_delem,
                               ParameterContainer *_pc,
                               void (*_update)(TACSElement*, TacsScalar*, void*) )
        TACSElement* getDeterministicElement()
        void updateElement(TACSElement* elem, TacsScalar* vals)
        void setPythonCallback(PyObject *cbptr)

# A simple test element for TACS
cdef extern from "smd.h":
    cdef cppclass SMD(TACSElement):
        SMD( TacsScalar, TacsScalar, TacsScalar)
        void setMass(TacsScalar)
        void setStiffness(TacsScalar)
        void setDamping(TacsScalar)
