# distutils: language = c++

# Import numpy
cimport numpy as np
import numpy as np
np.import_array()

# Import TACS
include "TacsDefs.pxi"

#from tacs.TACS cimport *
#from tacs.constitutive cimport *
#from tacs.elements cimport *

# Import STACS
#from cwrap cimport *

#from pspace.cwrap cimport *

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
