# Include the mpi4py header
cdef extern from "mpi-compat.h":
    pass

# Typdefs required for either real or complex mode
#include "TacsTypedefs.pxi"

#from TACS cimport *
from tacs.elements cimport *
from pspace.PSPACE cimport *

cdef extern from "TACSStochasticElement.h":
    cdef cppclass TACSStochasticElement(TACSElement):
        TACSStochasticElement( TACSElement *_delem,
                               ParameterContainer *_pc,
                               void (*_update)(TACSElement*, TacsScalar*) )
