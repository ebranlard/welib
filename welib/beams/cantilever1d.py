import numpy as np
from scipy.integrate import cumulative_trapezoid
#from welib.essentials import *


def cumtrapz_flipped(v, z, initial=0):
    return np.flip( cumulative_trapezoid(np.flip(v), np.flip(z), initial = initial) )



def deflection(p, z, EI, method='cumtrapz', Ftip=0, Mtip=0):
    """ 
    Compute deflection from a cantilever beam with stiffness EI(z) for a loading p(z)

    INPUTS:
     - p: array of size(n) of load per unit span
     - r : array of size(n) for position along the beam
     - EI: scalar or array of size(n) for position along the beam

    """
    n = len(z)
    EI = np.asarray(EI)
    
    # Shear and moment
    S = cumtrapz_flipped(p, z, initial=0) - Ftip
    M = cumtrapz_flipped(S, z, initial=0) + Mtip
    # Curvature
    kappa = M / EI
    # Slope and deflection
    theta = cumulative_trapezoid(kappa, z, initial=0)
    u     = cumulative_trapezoid(theta, z, initial=0)

    return u, theta, kappa, S, M



def deflection_all(pz, r, EI, method='cumtrapz', Ftip=0, Mtip=0):
    """ 
    Compute deflection from a cantilever beam with stiffness EI(z) for a loading p(z)

    TODO: NEED TO DECIDE ON A SIGN CONVENTION. FOR NOW, u, M are positive


    INPUTS:
     - pz: array of size(n) of load per unit span
     - r : array of size(n) for position along the beam
     - EI: scalar or array of size(n) for position along the beam

    """
    n = len(r)
    EI = np.asarray(EI)
    
    # Initialization including boundary condition at the tip
    My = np.zeros(n)
    Tz = np.zeros(n)
    
    # Loop to get T and M from p
    if method=='cumtrapz':
        Tz = cumtrapz_flipped(pz, r, initial = 0) - Ftip
        My = cumtrapz_flipped(Tz, r, initial = 0)
    elif method=='manual':
        Tz[-1] = -Ftip
        My[-1] = Mtip
        for i in range(n-1, 0, -1):
            Tz[i-1] = Tz[i] - 0.5 * (pz[i-1] + pz[i]) * (r[i] - r[i-1]) 
            My[i-1] = My[i] - Tz[i] * (r[i] - r[i-1]) - (1/6 * pz[i-1] + 1/3 * pz[i]) * (r[i] - r[i-1])**2
    else:
        raise NotImplementedError()

    # Going from principal directions to y and z
    ky =- My / EI
    
    # Initialization includes boundary condition at the hub
    thetay = np.zeros(n)
    uz = np.zeros(n)
    
    if method=='cumtrapz':
        thetay = cumulative_trapezoid(ky, r, initial = 0)
        uz     = - cumulative_trapezoid(thetay, r, initial = 0)
    elif method=='manual':
        for i in range(n-1):
            thetay[i+1] = thetay[i] + 0.5 * (ky[i+1] + ky[i]) * (r[i+1] - r[i])
            uz[i+1] = uz[i] - thetay[i] * (r[i+1] - r[i]) - (1/6 * ky[i+1] + 1/3 * ky[i]) * (r[i+1] - r[i])**2


    return uz, thetay, ky, Tz, My

def compute_modes(n_modes, r, EI, m, maxiter=500, abs_tol=0.01, method='cumtrapz'):
    """ """
    modes = np.zeros((n_modes, len(r)))
    freqs = np.zeros(n_modes)
    n = len(r)
    if not hasattr(EI, '__len__'):
        EI = np.asarray([EI]*n)
    if not hasattr(m, '__len__'):
        m = np.asarray([m]*n)
    
    for iMode in range(n_modes):
        pz = np.ones(n)
        omega2_prev = 0
        # --- Iteration loop
        for step in range(maxiter):
            uz , thetay, ky, Sy, Mz = deflection(pz, r, EI, method=method)
            # Orthogonalization against previous modes
            for iPrevMode in range(iMode):
                uz_p = modes[iPrevMode,:]
                const = np.trapezoid(uz_p * m * uz, r) / np.trapezoid(uz_p * m * uz_p, r)
                uz -=  const * uz_p
            # Normalization
            if uz[-1]==0:
                omega2, norm_factor = pz[-1] , 1   # Would unlikely happen
            else:
                omega2 = pz[-1] / (uz[-1] * m[-1])
                norm_factor = np.sqrt(uz[-1]**2)
            pz = omega2 * m * uz / norm_factor
            if abs(omega2 - omega2_prev) < abs_tol:
                #print(omega2, pz / (uz* m))
                break
            # Prepare for next iteration
            omega2_prev = omega2
        if step == maxiter - 1:
            raise Exception(f'Mode {iMode+1} did not converge')
        
        modes[iMode,:] = uz
        freqs[iMode]   = np.sqrt(omega2)/(2*np.pi)
        
    return freqs, modes

if __name__ == '__main__':
    from theory import UniformBeamDeflection
    from theory import UniformBeamBendingModes
    import matplotlib.pyplot as plt

    # --------------------------------------------------------------------------------}
    # --- Deflections 
    # --------------------------------------------------------------------------------{
    z = np.linspace(0, 1, 10)
    method ='default'
    method ='cumtrapz'

    params = {"F": 10, "EI": 2000, "p": 5, "p0": 8}
    types = []
#     types += ['point_load']
    types += ['uniform']
#     types += ['increasing'] # Shear is off
#     types += ['decreasing']

    for loading_type in types:
        u, S, M, M0, umax, p = UniformBeamDeflection('clamped-free', loading_type,  z, params)
        uz, thetay, ky, Sy, Mz = deflection(p, z, params['EI'], method=method)

        fig,axes = plt.subplots(3, 2, sharey=False, figsize=(6.4,8.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.30)
        ax = axes[0,0]
        ax.plot(z, u  , label='')
        ax.plot(z, uz , label='')
        ax.set_ylabel('Deflection')

        ax = axes[1,0]
    #     ax.plot(z, u  , label='')
        ax.plot(z, thetay , label='')
        ax.set_ylabel('Slope')

        ax = axes[2,0]
    #     ax.plot(z, u  , label='')
        ax.plot(z, ky , label='')
        ax.set_ylabel('Curvature')

        ax = axes[0,1]
        ax.plot(z, M  , label='')
        ax.plot(z, Mz , label='')
        ax.set_ylabel('Moment')

        ax = axes[1,1]
        ax.plot(z, S  , label='')
        ax.plot(z, Sy , label='')
        ax.set_ylabel('Shear')

        ax = axes[2,1]
        ax.plot(z, p  , label='')
    #     ax.plot(z, Sy , label='')
        ax.set_ylabel('Loading')
        ax.set_xlabel('')
#     ax.legend()


    # --------------------------------------------------------------------------------}
    # --- Modes 
    # --------------------------------------------------------------------------------{
    L  = 100    
    EI = 1.86e+12
    m  = 8.82e+03
    z = np.linspace(0, L, 100)
    nModes = 5
    rho, A = m , 1

    f_th, x, modes_th, ModesV, ModesK = UniformBeamBendingModes('unloaded-clamped-free', EI, rho, A, L, x=z, norm='tip', nModes=nModes)

    freq, modes = compute_modes(nModes, z, EI, m, maxiter=500, abs_tol=0.01)

    for i, (f,m) in enumerate(zip(freq, modes)):
        fr = f_th[i]
        print(f'Mode {i+1}: Frequency = {fr:.2f} {f:.2f} Hz')
        plt.figure()
        plt.plot(z, modes_th[i,:], 'k-', label='Num')
        plt.plot(z, m, '+', label='Num')

    plt.show()










    plt.show()

