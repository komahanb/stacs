import numpy as np
from mpi4py import MPI
import matplotlib.pyplot as plt

from pspace import PSPACE
from tacs import TACS, elements

legend=True

dmax=2
if dmax > 0:
    legend = False
    
def getJacobian(pc):
    M = pc.getNumQuadraturePoints()
    N = pc.getNumBasisTerms()
    A = np.zeros((N, N))
    for q in range(M):
        wq, zq, yq = pc.quadrature(q)
        for i in range(N):
            psiziq = pc.basis(i, zq)
            for j in range(N):
                psizjq = pc.basis(j, zq)                
                A[i,j] += wq*psiziq*psizjq*(yq[0]+yq[1]+yq[2])
    return A

class SMDUpdate:
    def __init__(self, elem):
        self.element = elem
        return

    def update(self, vals):
        self.element.setMass(vals[0])
        self.element.setStiffness(vals[1])
        self.element.setInitVelocity(vals[2])
        return

class ForceUpdate:
    def __init__(self, elem):
        self.element = elem
        return

    def update(self, vals):
        self.element.amplitude = vals[3]
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

# Define an element in TACS using the pyElement feature
class ForcingElement(elements.pyElement):
    def __init__(self, num_disps, num_nodes, amplitude, omega):
        self.amplitude = amplitude
        self.omega = omega

    def getInitConditions(self, index, X, v, dv, ddv):
        '''Define the initial conditions'''
        return

    def addResidual(self, index, time, X, v, dv, ddv, res):
        '''Add the residual of the governing equations'''
        res[0] += self.amplitude*np.sin(self.omega*time)
        return    

    def addJacobian(self, index, time, alpha, beta, gamma, X, v, dv, ddv, res, mat):
        '''Add the Jacobian of the governing equations'''
        self.addResidual(index, time, X, v, dv, ddv, res)
        return
    
def createAssembler(m=5.0, c=0.5, k=5.0, u0=-0.5, udot0=1.0, pc=None):
    num_disps = 1
    num_nodes = 1

    # Spring element
    #spr = SpringMassDamper(num_disps, num_nodes, m, c, k)
    dspr  = PSPACE.PySMD(m, c, k, u0, udot0)
    sprcb = SMDUpdate(dspr)
    sspr  = PSPACE.PyStochasticElement(dspr, pc, sprcb)

    dforce = ForcingElement(num_disps, num_nodes, amplitude=1.0, omega=10.0)
    forcecb = ForceUpdate(dforce)
    sforce = PSPACE.PyStochasticElement(dforce, pc, forcecb)

    ndof_per_node = 1*pc.getNumBasisTerms()
    num_owned_nodes = 1
    num_elems = 1        
    
    # Add user-defined element to TACS
    comm = MPI.COMM_WORLD
    assembler = TACS.Assembler.create(comm, ndof_per_node, num_owned_nodes, num_elems)

    ptr = np.array([0, 1], dtype=np.intc)
    conn = np.array([0], dtype=np.intc)
    assembler.setElementConnectivity(ptr, conn)

    # Set elements
    assembler.setElements([sspr])

    # set Auxiliary elements
    aux = TACS.AuxElements()
    aux.addElement(0, sforce)
    assembler.setAuxElements(aux)
    
    assembler.initialize()

    return assembler

def sgmmoments(bdf, num_steps, nterms):
    umean = np.zeros((num_steps+1))
    udotmean = np.zeros((num_steps+1))
    uddotmean = np.zeros((num_steps+1))
    time = np.zeros((num_steps+1))

    uvar = np.full_like(umean, 0)
    udotvar = np.full_like(udotmean, 0)
    uddotvar = np.full_like(uddotmean, 0)

    # Compute mean and variance at each time step
    for k in range(0, num_steps+1):
        # Extract the state vector
        t, uvec, udotvec, uddotvec = bdf.getStates(k)
        u = uvec.getArray()
        udot = udotvec.getArray()
        uddot = uddotvec.getArray()
        
        # Compute moments
        time[k] = t
        umean[k] = u[0]
        udotmean[k] = udot[0]
        uddotmean[k] = uddot[0]
        for i in range(1,nterms):
            uvar[k] += u[i]**2 
            udotvar[k] += udot[i]**2
            uddotvar[k] += uddot[i]**2
        
    return time, umean, udotmean, uddotmean, uvar, udotvar, uddotvar

pfactory = PSPACE.PyParameterFactory()
y1 = pfactory.createExponentialParameter(mu=4.0, beta=1.0, dmax=dmax) # mass
y2 = pfactory.createNormalParameter(mu=5.0, sigma=0.5, dmax=dmax) # stiff
y3 = pfactory.createUniformParameter(a=0.50, b=1.50, dmax=dmax) # init velocity
y4 = pfactory.createNormalParameter(mu=1.0, sigma=0.2, dmax=dmax) # amplitude

basis_type=1
pc = PSPACE.PyParameterContainer(basis_type)
pc.addParameter(y1)
pc.addParameter(y2)
pc.addParameter(y3)
pc.addParameter(y4)

pc.initialize()

## print("nterms ", pc.getNumBasisTerms())
## pmax = np.array(([5,5,5,5]))
## print(pmax)
## pc.initializeQuadrature(pmax)
## nqpts = pc.getNumQuadraturePoints()
## wq, zq, yq = pc.quadrature(nqpts)
## for i in range(nqpts):
##     print(i, wq, zq, yq)
## stop

## stop

## from pspace.plotter import plot_jacobian
## A = getJacobian(pc)
## plot_jacobian(A, 'smd-sparsity.pdf')

## stop

# Create TACS
m = 1.0
c = 0.5
k = 5.0
u0 = -0.5
udot0 = 1.0
assembler = createAssembler(m=m, c=c, k=k, u0=u0, udot0=udot0, pc=pc)

# Create Integrator
t0 = 0.0
tf = 10.0
num_steps = 1000
order = 2
integrator = TACS.BDFIntegrator(assembler, t0, tf, num_steps, order)
integrator.setPrintLevel(1)
integrator.integrate()

nterms = pc.getNumBasisTerms()

time, umean, udotmean, uddotmean, uvar, udotvar, uddotvar = sgmmoments(integrator, num_steps, nterms)

    
# Compute moments

###################################################################
# plot results
###################################################################

# Configure 
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'out'

# Optionally set font to Computer Modern to avoid common missing
# font errors
params = {
  'axes.labelsize': 20,
  'legend.fontsize': 14,
  'xtick.labelsize': 20,
  'ytick.labelsize': 20,
  'text.usetex': True}
plt.rcParams.update(params)

# Latex math
plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath}']
#plt.rcParams['font.family'] = 'sans-serif'
#plt.rcParams['font.sans-serif'] = 'courier'
plt.rcParams['font.size'] = 18
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['lines.color'] = 'r'

# Make sure everything is within the frame
plt.rcParams.update({'figure.autolayout': True})

# Set marker size
markerSize = 7.0 #11.0
mew = 2.0

# bar chart settings
lalpha    = 0.9
rev_alpha = 0.9/1.5

# These are the "Tableau 20" colors as RGB.    
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (25, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]    
  
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.)

# Define colors
poscolor = tableau20[2]
posbandcolor = tableau20[3]

velcolor = tableau20[16]
velbandcolor = tableau20[17]

acccolor = tableau20[0]
accbandcolor = tableau20[1]

alpha_level = 0.35

plt.figure()
fig, ax = plt.subplots()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.plot(time, umean    , '-', label='$u(t)$'       , mew=mew, ms=markerSize, color=poscolor, mec='black')
plt.plot(time[1:], udotmean[1:] , '-', label='$\dot{u}(t)$' , mew=mew, ms=markerSize, color=velcolor, mec='black')
plt.plot(time[1:], uddotmean[1:], '-', label='$\ddot{u}(t)$', mew=mew, ms=markerSize, color=acccolor, mec='black')
plt.xlabel('time [s]')
plt.ylabel('expectation')
if legend is True:
    plt.legend(loc='right', frameon=False)
plt.xlim(left=0.0,right=10.0)
plt.ylim(top=2.0,bottom=-2.0)
plt.savefig('smd-galerkin-expectation-%d-complete.pdf' % dmax,
            bbox_inches='tight', pad_inches=0.05)

plt.figure()
fig, ax = plt.subplots()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.plot(time, uvar    , '-', label='$u(t)$'       , mew=mew, ms=markerSize, color=poscolor, mec='black')
plt.plot(time[1:], udotvar[1:] , '-', label='$\dot{u}(t)$' , mew=mew, ms=markerSize, color=velcolor, mec='black')
plt.plot(time[1:], uddotvar[1:], '-', label='$\ddot{u}(t)$', mew=mew, ms=markerSize, color=acccolor, mec='black')
plt.xlabel('time [s]')
plt.ylabel('variance')
if legend is True:
    plt.legend(loc='upper right', frameon=False)
plt.savefig('smd-galerkin-variance-%d-complete.pdf' % dmax,
            bbox_inches='tight', pad_inches=0.05)

plt.figure()
sigma = 1.0
fig, ax = plt.subplots()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

plt.plot(time, umean    , '-', label='${E}[{u}(t)]\pm{S}[{u}(t)]$', mew=mew, ms=markerSize, color=poscolor, mec='black')
plt.fill_between(time, umean, umean + sigma*np.sqrt(uvar), color=posbandcolor, alpha=alpha_level)
plt.fill_between(time, umean, umean - sigma*np.sqrt(uvar), color=posbandcolor, alpha=alpha_level)

plt.plot(time, udotmean , '-', label='${E}[\dot{u}(t)]\pm{S}[\dot{u}(t)]$', mew=mew, ms=markerSize, color=velcolor, mec='black')
plt.fill_between(time, udotmean, udotmean + sigma*np.sqrt(udotvar), color=velbandcolor, alpha=alpha_level)
plt.fill_between(time, udotmean, udotmean - sigma*np.sqrt(udotvar), color=velbandcolor, alpha=alpha_level)

plt.plot(time[1:], uddotmean[1:], '-', label='${E}[\ddot{u}(t)]\pm{S}[\ddot{u}(t)]$', mew=mew, ms=markerSize, color=acccolor, mec='black')
plt.fill_between(time[1:], uddotmean[1:], uddotmean[1:] + sigma*np.sqrt(uddotvar[1:]), color=accbandcolor, alpha=alpha_level)
plt.fill_between(time[1:], uddotmean[1:], uddotmean[1:] - sigma*np.sqrt(uddotvar[1:]), color=accbandcolor, alpha=alpha_level)

plt.xlim(left=0.0,right=10.0)
plt.ylim(top=2.0,bottom=-2.0)

plt.xlabel('time [s]')
#plt.ylabel('response')
if legend is True:
    plt.legend(loc='upper right', frameon=False)
plt.savefig('smd-galerkin-one-sigma-%d-complete.pdf' % dmax,
            bbox_inches='tight', pad_inches=0.05)

plt.figure()
sigma = 2.0
fig, ax = plt.subplots()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

plt.plot(time, umean    , '-', label='${E}[{u}(t)]\pm 2\cdot{S}[{u}(t)]$', mew=mew, ms=markerSize, color=poscolor, mec='black')
plt.fill_between(time, umean, umean + sigma*np.sqrt(uvar), color=posbandcolor, alpha=alpha_level)
plt.fill_between(time, umean, umean - sigma*np.sqrt(uvar), color=posbandcolor, alpha=alpha_level)

plt.plot(time, udotmean , '-', label='${E}[\dot{u}(t)]\pm 2\cdot{S}[\dot{u}(t)]$', mew=mew, ms=markerSize, color=velcolor, mec='black')
plt.fill_between(time, udotmean, udotmean + sigma*np.sqrt(udotvar), color=velbandcolor, alpha=alpha_level)
plt.fill_between(time, udotmean, udotmean - sigma*np.sqrt(udotvar), color=velbandcolor, alpha=alpha_level)

plt.plot(time[1:], uddotmean[1:], '-', label='${E}[\ddot{u}(t)]\pm 2\cdot{S}[\ddot{u}(t)]$', mew=mew, ms=markerSize, color=acccolor, mec='black')
plt.fill_between(time[1:], uddotmean[1:], uddotmean[1:] + sigma*np.sqrt(uddotvar[1:]), color=accbandcolor, alpha=alpha_level)
plt.fill_between(time[1:], uddotmean[1:], uddotmean[1:] - sigma*np.sqrt(uddotvar[1:]), color=accbandcolor, alpha=alpha_level)

plt.xlim(left=0.0,right=10.0)
plt.ylim(top=2.0,bottom=-2.0)

plt.xlabel('time [s]')
#plt.ylabel('response')
if legend is True:
    plt.legend(loc='upper right', frameon=False)
plt.savefig('smd-galerkin-two-sigma-%d-complete.pdf' % dmax,
            bbox_inches='tight', pad_inches=0.05)


plt.figure()
sigma = 3.0
fig, ax = plt.subplots()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

plt.plot(time, umean    , '-', label='${E}[{u}(t)]\pm 3\cdot{S}[{u}(t)]$', mew=mew, ms=markerSize, color=poscolor, mec='black')
plt.fill_between(time, umean, umean + sigma*np.sqrt(uvar), color=posbandcolor, alpha=alpha_level)
plt.fill_between(time, umean, umean - sigma*np.sqrt(uvar), color=posbandcolor, alpha=alpha_level)

plt.plot(time, udotmean , '-', label='${E}[\dot{u}(t)]\pm 3\cdot{S}[\dot{u}(t)]$', mew=mew, ms=markerSize, color=velcolor, mec='black')
plt.fill_between(time, udotmean, udotmean + sigma*np.sqrt(udotvar), color=velbandcolor, alpha=alpha_level)
plt.fill_between(time, udotmean, udotmean - sigma*np.sqrt(udotvar), color=velbandcolor, alpha=alpha_level)

plt.plot(time[1:], uddotmean[1:], '-', label='${E}[\ddot{u}(t)]\pm 3\cdot{S}[\ddot{u}(t)]$', mew=mew, ms=markerSize, color=acccolor, mec='black')
plt.fill_between(time[1:], uddotmean[1:], uddotmean[1:] + sigma*np.sqrt(uddotvar[1:]), color=accbandcolor, alpha=alpha_level)
plt.fill_between(time[1:], uddotmean[1:], uddotmean[1:] - sigma*np.sqrt(uddotvar[1:]), color=accbandcolor, alpha=alpha_level)

plt.xlim(left=0.0,right=10.0)
plt.ylim(top=2.0,bottom=-2.0)

plt.xlabel('time [s]')
#plt.ylabel('response')

if legend is True:
    plt.legend(loc='upper right', frameon=False)
plt.savefig('smd-galerkin-three-sigma-%d-complete.pdf' % dmax,
            bbox_inches='tight', pad_inches=0.05)
