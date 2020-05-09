import numpy as np

def moments(bdf, num_steps, nterms):
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



