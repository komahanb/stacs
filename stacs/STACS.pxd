# distutils: language = c++

from pspace.PSPACE cimport *

include "TacsTypedefs.pxi"

# Include the mpi4py header
cdef extern from "mpi-compat.h":
    pass

from tacs.elements cimport *

cdef extern from "TACSStochasticElement.h":
    cdef cppclass TACSStochasticElement(TACSElement):
        TACSStochasticElement( TACSElement *_delem,
                               ParameterContainer *_pc,
                               void (*_update)(TACSElement*, TacsScalar*, void*) )
        TACSElement* getDeterministicElement()
        void updateElement(TACSElement* elem, TacsScalar* vals)
        void setPythonCallback(PyObject *cbptr)

cdef extern from "TACSMutableElement3D.h":
    cdef cppclass TACSMutableElement3D(TACSElement):
        TACSMutableElement3D( TACSElement *_elem )
        void setDensity( TacsScalar _rho )

# A simple test element for TACS
cdef extern from "smd.h":
    cdef cppclass SMD(TACSElement):
        SMD( TacsScalar, TacsScalar, TacsScalar, TacsScalar, TacsScalar)
        void setMass(TacsScalar)
        void setStiffness(TacsScalar)
        void setDamping(TacsScalar)
        void setInitPosition(TacsScalar)
        void setInitVelocity(TacsScalar)

from tacs.functions cimport *

cdef extern from "TACSStochasticFunction.h":
    cdef cppclass TACSStochasticFunction(TACSFunction):
        TACSStochasticFunction( TACSAssembler *tacs,
                                TACSFunction *dfunc,
                                ParameterContainer *pc,
                                int quantityType,
                                int moment_type )
        TacsScalar getFunctionValue()

cdef extern from "TACSKSStochasticFunction.h":
    cdef cppclass TACSKSStochasticFunction(TACSFunction):
        TACSKSStochasticFunction( TACSAssembler *tacs,
                                  TACSFunction *dfunc,
                                  ParameterContainer *pc,
                                  int quantityType,
                                  int moment_type,
                                  int ksweight)
        TacsScalar getFunctionValue()
