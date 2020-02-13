import numpy as np
from mpi4py import MPI
from pspace import PSPACE
from tacs import TACS, elements

class SMDUpdate:
    def __init__(self, elem):
        self.element = elem
        return

    def update(self, vals):
        self.element.setMass(vals[0])
        self.element.setStiffness(vals[1])
        self.element.setDamping(vals[2])
        #self.element.m = vals[0]
        #self.element.c = vals[1] 
        #self.element.k = vals[2]
        return

# Define an element in TACS using the pyElement feature
class SpringMassDamper(elements.pyElement):
    def __init__(self, num_disps, num_nodes, m, c, k):
        self.m = m
        self.c = c
        self.k = k    

    def getInitConditions(self, index, X, v, dv, ddv):
        '''Define the initial conditions'''
        v[0] = -0.5
        dv[0] = 1.0
        return

    def addResidual(self, index, time, X, v, dv, ddv, res):
        '''Add the residual of the governing equations'''
        res[0] += self.m*ddv[0] + self.c*dv[0] + self.k*v[0]
        return    

    def addJacobian(self, index, time, alpha, beta, gamma, X, v, dv, ddv, res, mat):
        '''Add the Jacobian of the governing equations'''
        res[0] += self.m*ddv[0] + self.c*dv[0] + self.k*v[0]
        mat[0] += alpha*self.k + beta*self.c + gamma*self.m
        return

def createAssembler(m=1.0, c=0.5, k=5.0, pc=None):
    num_disps = 1
    num_nodes = 1
    #spr = SpringMassDamper(num_disps, num_nodes, m, c, k)
    spr = PSPACE.PySMD(m, c, k)
    elem = spr
    ndof_per_node = 1
    num_owned_nodes = 1
    num_elems = 1
    if pc is not None:
        cb = SMDUpdate(spr)
        elem = PSPACE.PyStochasticElement(spr, pc, cb)
        ndof_per_node = ndof_per_node*pc.getNumBasisTerms()
    
    # Add user-defined element to TACS
    comm = MPI.COMM_WORLD
    assembler = TACS.Assembler.create(comm, ndof_per_node, num_owned_nodes, num_elems)

    ptr = np.array([0, 1], dtype=np.intc)
    conn = np.array([0], dtype=np.intc)
    assembler.setElementConnectivity(ptr, conn)
    assembler.setElements([elem])
    assembler.initialize()

    return assembler

pfactory = PSPACE.PyParameterFactory()
y1 = pfactory.createNormalParameter(mu=1.0, sigma=0.1, dmax=2)
y2 = pfactory.createUniformParameter(a=1.0, b=0.1, dmax=2)
y3 = pfactory.createExponentialParameter(mu=1.0, beta=0.1, dmax=2)

basis_type = 1
pc = PSPACE.PyParameterContainer(basis_type)
pc.addParameter(y1)
pc.addParameter(y2)
pc.addParameter(y3)

pc.initialize()

# Create TACS
m = 1.0
c = 0.5
k = 5.0
tf = 10.0
assembler = createAssembler(m=m, c=c, k=k, pc=pc)

# Create Integrator
t0 = 0.0
tf = 1.0
num_steps = 100
order = 2
integrator = TACS.BDFIntegrator(assembler, t0, tf, num_steps, order)
integrator.setPrintLevel(1)
integrator.integrate()
