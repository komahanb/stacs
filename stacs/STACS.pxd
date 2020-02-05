# Typdefs required for either real or complex mode
include "TacsTypedefs.pxi"

cdef extern from "TACSElement.h":
    cdef cppclass TACSElement:
        pass

cdef extern from "ParameterContainer.h":
    cdef cppclass ParameterContainer:
        pass

cdef extern from "TACSStochasticElement.h":
    cdef cppclass TACSStochasticElement(TACSElement):
        TACSStochasticElement( TACSElement *_delem,
                               ParameterContainer *_pc,
                               void (*_update)(TACSElement*, TacsScalar*) )
