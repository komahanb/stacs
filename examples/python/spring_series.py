import numpy as np
from mpi4py import MPI
import matplotlib.pyplot as plt

from pspace import PSPACE
from stacs import STACS
from tacs import TACS, elements

from ProjectionDriver import ProjectionDriver, SamplingDriver

legend=True

dmax=7
if dmax > 0:
    legend = True
    
def getJacobian(pc):
    M = pc.getNumQuadraturePoints()
    N = pc.getNumBasisTerms()
    A = np.zeros((N, N))
    for q in range(M):
        wq, zq, yq = pc.quadrature(q)
        for i in range(N):
            psiziq = pc.basis(i, z)
            for j in range(N):
                psizjq = pc.basis(j, zq)                
                A[i,j] += wq*psiziq*psizjq*(yq[0]+yq[1]+yq[2])
    return A

class SMDUpdate:
    def __init__(self, elem):
        self.element = elem
        return

    def update(self, vals):
        self.element.m[0,0] = vals[0]
        self.element.m[1,1] = vals[1]
        self.element.m[2,2] = vals[2]
        return

## class ForceUpdate:
##     def __init__(self, elem):
##         self.element = elem
##         return

##     def update(self, vals):
##         self.element.amplitude = vals[3]
##         return

# Define an element in TACS using the pyElement feature
class SpringMassDamper(elements.pyElement):
    def __init__(self, num_disps, num_nodes, m, k):
        self.m = m
        self.k = k
        return

    def getInitConditions(self, index, X, v, dv, ddv):
        '''Define the initial conditions'''
        v[0] = 0.0        
        v[1] = 0.0
        v[2] = 0.0
        dv[2] = 0.1    
        return

    def addResidual(self, index, time, X, v, dv, ddv, res):
        '''Add the residual of the governing equations'''
        res[:] += np.squeeze(np.array(np.dot(self.m,ddv))) + np.squeeze(np.array(np.dot(self.k,v)))[:]
        return    

    def addJacobian(self, index, time, alpha, beta, gamma, X, v, dv, ddv, res, mat):
        res[:] += np.squeeze(np.array(np.dot(self.m,ddv))) + np.squeeze(np.array(np.dot(self.k,v)))[:]
        mat[:] += np.squeeze(np.array(np.reshape(gamma*self.m + alpha*self.k, (9,))))[:]
        return

## # Define an element in TACS using the pyElement feature
## class ForcingElement(elements.pyElement):
##     def __init__(self, num_disps, num_nodes, amplitude, omega):
##         self.amplitude = amplitude
##         self.omega = omega

##     def getInitConditions(self, index, X, v, dv, ddv):
##         '''Define the initial conditions'''
##         return

##     def addResidual(self, index, time, X, v, dv, ddv, res):
##         '''Add the residual of the governing equations'''
##         res[0] += self.amplitude*np.sin(self.omega*time)
##         return    

##     def addJacobian(self, index, time, alpha, beta, gamma, X, v, dv, ddv, res, mat):
##         '''Add the Jacobian of the governing equations'''
##         self.addResidual(index, time, X, v, dv, ddv, res)
##         return

def createSolver(params, pc):
    '''
    Creating solver
    '''
    if pc is not None:
        pc.initialize()
        print("number of basis terms = ", pc.getNumBasisTerms())
        
    if params is not None:
        m1 = params[0]
        m2 = params[1]
        m3 = params[2]
    else:
        m1 = 1.0
        m2 = 10.0
        m3 = 100.0
    
    # Create TACS
    M = np.matrix([[m1, 0.0, 0.0],
                   [0.0, m2, 0.0],
                   [0.0, 0.0, m3]])
    
    k1 = 1.0; k2 = 10.0; k3 = 100.0; k4 = 1000.0;
    K = np.matrix([[k1 + k2, -k2, 0.0],
                   [-k2, k2 + k3, -k3],
                   [0.0, -k3, k3 + k4]])
    
    # Spring element
    num_disps = 3
    num_nodes = 1
    dspr = SpringMassDamper(num_disps, num_nodes, M, K)

    if pc is not None:
        sprcb = SMDUpdate(dspr)
        sspr  = STACS.PyStochasticElement(dspr, pc, sprcb)

    # dforce = ForcingElement(num_disps, num_nodes, amplitude=1.0, omega=10.0)
    # forcecb = ForceUpdate(dforce)
    # sforce = STACS.PyStochasticElement(dforce, pc, forcecb)

    if pc is not None:
        ndof_per_node = num_disps*pc.getNumBasisTerms()
    else:
        ndof_per_node = num_disps
        
    num_owned_nodes = 1
    num_elems = 1        
    
    # Add user-defined element to TACS
    comm = MPI.COMM_WORLD
    assembler = TACS.Assembler.create(comm, ndof_per_node, num_owned_nodes, num_elems)

    ptr = np.array([0, 1], dtype=np.intc)
    conn = np.array([0], dtype=np.intc)
    assembler.setElementConnectivity(ptr, conn)

    # Set elements
    if pc is not None:
        assembler.setElements([sspr])
    else:
        assembler.setElements([dspr])
        
    # set Auxiliary elements
    # aux = TACS.AuxElements()
    # aux.addElement(0, sforce)
    # assembler.setAuxElements(aux)
    
    assembler.initialize()

    # Create Integrator
    t0 = 0.0
    tf = 2.0
    num_steps = 100
    order = 2
    integrator = TACS.BDFIntegrator(assembler, t0, tf, num_steps, order)
    integrator.setPrintLevel(1)
    integrator.integrate()
   
    return integrator

def ssmmoments(num_steps, pc):
    #num_steps  = bdf.getNumTimeSteps()
    nqpts = np.array([10, 10, 10], dtype=np.intc)
    pc.initializeQuadrature(nqpts)
    nterms     = pc.getNumBasisTerms()    
    vmean      = np.zeros((3,num_steps+1))
    v2mean     = np.zeros((3,num_steps+1))
    vvar       = np.full_like(umean, 0)
    time       = np.zeros((num_steps+1))
    for q in range(pc.getNumQuadraturePoints()):       
        wq, zq, yq = pc.quadrature(q)        
        sbdf = createSolver(yq, None)
        num_steps = sbdf.getNumTimeSteps()
        for k in range(0, num_steps+1):
            t, uvec, udotvec, uddotvec = sbdf.getStates(k)
            time[k] = t
            u = uvec.getArray(); udot = udotvec.getArray(); uddot = uddotvec.getArray()
            vmean[:,k] += wq*u[:]
            v2mean[:,k] += wq*(u[:]**2)
        vvar = v2mean - vmean*vmean
    return time, vmean, vvar

def sgmmoments(bdf, pc):
    num_steps = bdf.getNumTimeSteps()
    nterms = pc.getNumBasisTerms()    
    umean = np.zeros((num_steps+1))
    udotmean = np.zeros((num_steps+1))
    uddotmean = np.zeros((num_steps+1))
    time = np.zeros((num_steps+1))
    uvar = np.full_like(umean, 0)
    udotvar = np.full_like(udotmean, 0)
    uddotvar = np.full_like(uddotmean, 0)

    dof = 0
    
    # Compute mean and variance at each time step
    for k in range(0, num_steps+1):
        # Extract the state vector
        t, uvec, udotvec, uddotvec = bdf.getStates(k)
        u = uvec.getArray()
        udot = udotvec.getArray()
        uddot = uddotvec.getArray()
        
        # Compute moments
        time[k] = t
        umean[k] = u[dof]
        udotmean[k] = u[dof+1]
        uddotmean[k] = u[dof+2]
        for i in range(1,nterms):
            uvar[k] += u[i*3]**2 
            udotvar[k] += u[i*3+1]**2
            uddotvar[k] += u[i*3+2]**2
        
    return time, umean, udotmean, uddotmean, uvar, udotvar, uddotvar

if __name__ == "__main__":
    
    pfactory = PSPACE.PyParameterFactory()
    y1 = pfactory.createUniformParameter(a=1.0, b=2.0, dmax=dmax)          # mass 1
    y2 = pfactory.createExponentialParameter(mu=10.0, beta=1.0, dmax=dmax) # mass 2
    y3 = pfactory.createNormalParameter(mu=100.0, sigma=10.0, dmax=dmax)   # mass 3
    
    pc = PSPACE.PyParameterContainer(basis_type=1)
    pc.addParameter(y1)
    pc.addParameter(y2)
    pc.addParameter(y3)
      

    bdf = createSolver(None, pc)
    time, umean, udotmean, uddotmean, uvar, udotvar, uddotvar = sgmmoments(bdf, pc)

    bdf = createSolver(None, None)
    time, vmean, vvar = ssmmoments(bdf.getNumTimeSteps(), pc)

    print(pc.getNumBasisTerms(), np.linalg.norm(vmean[0,:] - umean), np.linalg.norm(vvar[0,:] - uvar))

    ###################################################################
    # plot results
    ###################################################################
    
    # Configure 
    plt.rcParams['xtick.direction'] = 'out'
    plt.rcParams['ytick.direction'] = 'out'
    
    # Optionally set font to Computer Modern to avoid common missing
    # font errors
    params = {
      'axes.labelsize': 16,
      'legend.fontsize': 14,
      'xtick.labelsize': 16,
      'ytick.labelsize': 16,
      'text.usetex': True}
    plt.rcParams.update(params)
    
    # Latex math
    plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath}\usepackage{amsmath}\usepackage{amssymb}']
    #plt.rcParams['font.family'] = 'sans-serif'
    #plt.rcParams['font.sans-serif'] = 'courier'
    plt.rcParams['font.size'] = 16
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
    
    
    poscolor = 'r'
    posbandcolor = 'r'
    
    velcolor = 'g'
    velbandcolor = 'g'
    
    acccolor = 'b'
    accbandcolor = 'b'
    
    alpha_level = 0.3
    
    ## 
    ## plt.figure()
    ## fig, ax = plt.subplots()
    ## ax.spines['right'].set_visible(False)
    ## ax.spines['top'].set_visible(False)
    ## ax.xaxis.set_ticks_position('bottom')
    ## ax.yaxis.set_ticks_position('left')
    ## plt.plot(time, umean    , '-', label='$u(t)$'       , mew=mew, ms=markerSize, color=poscolor, mec='black')
    ## plt.plot(time[1:], udotmean[1:] , '-', label='$\dot{u}(t)$' , mew=mew, ms=markerSize, color=velcolor, mec='black')
    ## plt.plot(time[1:], uddotmean[1:], '-', label='$\ddot{u}(t)$', mew=mew, ms=markerSize, color=acccolor, mec='black')
    ## plt.xlabel('time [s]')
    ## plt.ylabel('expectation')
    ## if legend is True:
    ##     plt.legend(loc='right', frameon=False)
    ## #plt.xlim(left=0.0,right=10.0)
    ## #plt.ylim(top=2.0,bottom=-2.0)
    ## plt.savefig('smd-galerkin-expectation-%d-complete.pdf' % dmax,
    ##             bbox_inches='tight', pad_inches=0.05)
    ## 
    ## plt.figure()
    ## fig, ax = plt.subplots()
    ## ax.spines['right'].set_visible(False)
    ## ax.spines['top'].set_visible(False)
    ## ax.xaxis.set_ticks_position('bottom')
    ## ax.yaxis.set_ticks_position('left')
    ## plt.plot(time, uvar    , '-', label='$u(t)$'       , mew=mew, ms=markerSize, color=poscolor, mec='black')
    ## plt.plot(time[1:], udotvar[1:] , '-', label='$\dot{u}(t)$' , mew=mew, ms=markerSize, color=velcolor, mec='black')
    ## plt.plot(time[1:], uddotvar[1:], '-', label='$\ddot{u}(t)$', mew=mew, ms=markerSize, color=acccolor, mec='black')
    ## plt.xlabel('time [s]')
    ## plt.ylabel('variance')
    ## if legend is True:
    ##     plt.legend(loc='upper right', frameon=False)
    ## plt.savefig('smd-galerkin-variance-%d-complete.pdf' % dmax,
    ##             bbox_inches='tight', pad_inches=0.05)
    ## 
    ##

    

    plt.figure()
    sigma = 1.0
    fig, ax = plt.subplots()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    
    plt.plot(time, umean*100    , '-', label='$\mathrm{projection~u_1}$', mew=mew, ms=markerSize, color=poscolor, mec='black', markevery=5)
    plt.fill_between(time, umean*100, (umean + sigma*np.sqrt(uvar))*100, color=posbandcolor, alpha=alpha_level)
    plt.fill_between(time, umean*100, (umean - sigma*np.sqrt(uvar))*100, color=posbandcolor, alpha=alpha_level)

    plt.plot(time, udotmean*100    , '-', label='$\mathrm{projection~u_2}$', mew=mew, ms=markerSize, color=velcolor, mec='black', markevery=5)
    plt.fill_between(time, udotmean*100, (udotmean + sigma*np.sqrt(udotvar))*100, color=velbandcolor, alpha=alpha_level)
    plt.fill_between(time, udotmean*100, (udotmean - sigma*np.sqrt(udotvar))*100, color=velbandcolor, alpha=alpha_level)

    plt.plot(time, uddotmean*100    , '-', label='$\mathrm{projection~u_3}$', mew=mew, ms=markerSize, color=acccolor, mec='black', markevery=5)
    plt.fill_between(time, uddotmean*100, (uddotmean + sigma*np.sqrt(uddotvar))*100, color=accbandcolor, alpha=alpha_level)
    plt.fill_between(time, uddotmean*100, (uddotmean - sigma*np.sqrt(uddotvar))*100, color=accbandcolor, alpha=alpha_level)

    plt.plot(time, vmean[0,:]*100    , 'o', label='$\mathrm{sampling~u_1}$', mew=mew, ms=markerSize/1.5, color=poscolor, mec=poscolor, mfc='white', markevery=20, alpha=1)
    plt.plot(time, (vmean[0,:] + sigma*np.sqrt(vvar[0,:]))*100, 'o', mew=mew, ms=markerSize/2, color=poscolor, mec=posbandcolor, mfc='white', markevery=10, alpha=alpha_level)
    plt.plot(time, (vmean[0,:] - sigma*np.sqrt(vvar[0,:]))*100, 'o', mew=mew, ms=markerSize/2, color=poscolor, mec=posbandcolor, mfc='white', markevery=10, alpha=alpha_level)

    plt.plot(time, vmean[1,:]*100    , 'o', label='$\mathrm{sampling~u_2}$', mew=mew, ms=markerSize/1.5, color=velcolor, mec=velcolor, mfc='white', markevery=20, alpha=1)
    plt.plot(time, (vmean[1,:] + sigma*np.sqrt(vvar[1,:]))*100, 'o', mew=mew, ms=markerSize/2, color=velcolor, mec=velbandcolor, mfc='white', markevery=10, alpha=alpha_level)
    plt.plot(time, (vmean[1,:] - sigma*np.sqrt(vvar[1,:]))*100, 'o', mew=mew, ms=markerSize/2, color=velcolor, mec=velbandcolor, mfc='white', markevery=10, alpha=alpha_level)

    plt.plot(time, vmean[2,:]*100    , 'o', label='$\mathrm{sampling~u_3}$', mew=mew, ms=markerSize/1.5, color=acccolor, mec=acccolor, mfc='white', markevery=20, alpha=1)
    plt.plot(time, (vmean[2,:] + sigma*np.sqrt(vvar[2,:]))*100, 'o', mew=mew, ms=markerSize/2, color=acccolor, mec=accbandcolor, mfc='white', markevery=10, alpha=alpha_level)
    plt.plot(time, (vmean[2,:] - sigma*np.sqrt(vvar[2,:]))*100, 'o', mew=mew, ms=markerSize/2, color=acccolor, mec=accbandcolor, mfc='white', markevery=10, alpha=alpha_level)
    
    plt.ylim(top=8,bottom=-8)
    
    plt.xlabel('time')
    plt.ylabel('$\mathbb{E}[\mathrm{u(t)}]\pm\mathbb{S}[\mathrm{u(t)}]~\mathrm{x}~10^{-2}$')
    if legend is True:
        plt.legend(loc='lower left', frameon=False,ncol=2,handletextpad=0.2)           
    plt.savefig('smd-galerkin-one-sigma-complete-dmax%d-nterms%d.pdf' % (dmax,pc.getNumBasisTerms()),
                bbox_inches='tight', pad_inches=0.00)




    
    ## plt.figure()
    ## sigma = 2.0
    ## fig, ax = plt.subplots()
    ## ax.spines['right'].set_visible(False)
    ## ax.spines['top'].set_visible(False)
    ## ax.xaxis.set_ticks_position('bottom')
    ## ax.yaxis.set_ticks_position('left')
    
    ## plt.plot(time, umean    , '-', label='${E}[{u}(t)]\pm 2\cdot{S}[{u}(t)]$', mew=mew, ms=markerSize, color=poscolor, mec='black')
    ## plt.fill_between(time, umean, umean + sigma*np.sqrt(uvar), color=posbandcolor, alpha=alpha_level)
    ## plt.fill_between(time, umean, umean - sigma*np.sqrt(uvar), color=posbandcolor, alpha=alpha_level)
    
    ## plt.plot(time, udotmean , '-', label='${E}[\dot{u}(t)]\pm 2\cdot{S}[\dot{u}(t)]$', mew=mew, ms=markerSize, color=velcolor, mec='black')
    ## plt.fill_between(time, udotmean, udotmean + sigma*np.sqrt(udotvar), color=velbandcolor, alpha=alpha_level)
    ## plt.fill_between(time, udotmean, udotmean - sigma*np.sqrt(udotvar), color=velbandcolor, alpha=alpha_level)
    
    ## plt.plot(time[1:], uddotmean[1:], '-', label='${E}[\ddot{u}(t)]\pm 2\cdot{S}[\ddot{u}(t)]$', mew=mew, ms=markerSize, color=acccolor, mec='black')
    ## plt.fill_between(time[1:], uddotmean[1:], uddotmean[1:] + sigma*np.sqrt(uddotvar[1:]), color=accbandcolor, alpha=alpha_level)
    ## plt.fill_between(time[1:], uddotmean[1:], uddotmean[1:] - sigma*np.sqrt(uddotvar[1:]), color=accbandcolor, alpha=alpha_level)
    
    ## #plt.xlim(left=0.0,right=10.0)
    ## #plt.ylim(top=2.0,bottom=-2.0)
    
    ## plt.xlabel('time [s]')
    ## #plt.ylabel('response')
    ## if legend is True:
    ##     plt.legend(loc='upper right', frameon=False)
    ## plt.savefig('smd-galerkin-two-sigma-dmax%d-complete.pdf' % dmax,
    ##             bbox_inches='tight', pad_inches=0.05)
    
    
    ## plt.figure()
    ## sigma = 3.0
    ## fig, ax = plt.subplots()
    ## ax.spines['right'].set_visible(False)
    ## ax.spines['top'].set_visible(False)
    ## ax.xaxis.set_ticks_position('bottom')
    ## ax.yaxis.set_ticks_position('left')
    
    ## plt.plot(time, umean    , '-', label='${E}[{u}(t)]\pm 3\cdot{S}[{u}(t)]$', mew=mew, ms=markerSize, color=poscolor, mec='black')
    ## plt.fill_between(time, umean, umean + sigma*np.sqrt(uvar), color=posbandcolor, alpha=alpha_level)
    ## plt.fill_between(time, umean, umean - sigma*np.sqrt(uvar), color=posbandcolor, alpha=alpha_level)
    
    ## plt.plot(time, udotmean , '-', label='${E}[\dot{u}(t)]\pm 3\cdot{S}[\dot{u}(t)]$', mew=mew, ms=markerSize, color=velcolor, mec='black')
    ## plt.fill_between(time, udotmean, udotmean + sigma*np.sqrt(udotvar), color=velbandcolor, alpha=alpha_level)
    ## plt.fill_between(time, udotmean, udotmean - sigma*np.sqrt(udotvar), color=velbandcolor, alpha=alpha_level)
    
    ## plt.plot(time[1:], uddotmean[1:], '-', label='${E}[\ddot{u}(t)]\pm 3\cdot{S}[\ddot{u}(t)]$', mew=mew, ms=markerSize, color=acccolor, mec='black')
    ## plt.fill_between(time[1:], uddotmean[1:], uddotmean[1:] + sigma*np.sqrt(uddotvar[1:]), color=accbandcolor, alpha=alpha_level)
    ## plt.fill_between(time[1:], uddotmean[1:], uddotmean[1:] - sigma*np.sqrt(uddotvar[1:]), color=accbandcolor, alpha=alpha_level)
    
    ## #plt.xlim(left=0.0,right=10.0)
    ## #plt.ylim(top=2.0,bottom=-2.0)
    
    ## plt.xlabel('time [s]')
    ## #plt.ylabel('response')
    
    ## if legend is True:
    ##     plt.legend(loc='upper right', frameon=False)
    ## plt.savefig('smd-galerkin-three-sigma-%d-complete.pdf' % dmax,
    ##             bbox_inches='tight', pad_inches=0.05)
