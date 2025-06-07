import numpy as np
import pandas as pd
try:
    from numpy import trapezoid
except:
    from numpy import trapz as trapezoid

# import welib.hydro.wavekin as wk
from welib.hydro.wavekin import kinematics2d


def inline_load(u, du, D, CD, CM, rho):
    """
    Morison force per length on a structure defined along the z axis
    INPUTS:
     - u: velocity, array like of shape (nz x ...)
     - du: acceleration, array like of shape (nz x ...)
     - z: vertical distance along the structure (typically negative and 0 at sea bed) 
     - D: structure diameter, scalar or array of shape (nz x ...)
     - CD: drag coefficient, scalar or array of shape (nz x ...)
     - CM: inertia coefficient, scalar or array of shape (nz x ...)
    """
    p_drag    = .5 * rho * CD * D * u * np.abs(u)
    p_inertia = rho * CM * np.pi*(D**2)/4 * du 
#     dM_ref  = P * (z-zref)                # [Nm/m] 
#     F_total = trapezoid(P, z, axis=0)      # [N] integrated load
#     M_total = trapezoid(dM_ref, z, axis=0) # [Nm]

    return p_drag + p_inertia


def monopileHydroLoads1D(ai, fi, ki, epsi, h, t, z, x, D, CD, CM, rho, u_struct):  
    """ 
    Compute hydrodynamic loads (Morison equation) on monopile
    Wave kinematics are computed on the fly at undeflected monopile position

    INPUTS:
       fi  : (nf-array) frequencies, for wave kinematics
       ai  : (nf-array) amplitudes
       ki  : (nf-array) wave numbers
       epsi: (nf-array) phases
       h: water depth >0, the sea bed is at z=-h [m]
       t: time [s]
       z: (n-array) vertical positions defining monopile (and possibly tower) sections
    """
        
    # Limit calculations to below water and above sea bed
    bWet = np.logical_and(z<=0, z>=-h)  # TODO eta
    D    = D[bWet].flatten()
    CD   = CD[bWet].flatten()
    CM   = CM[bWet].flatten()
    zWet = z[bWet].flatten()
    u_struct = u_struct[bWet].flatten()
    
    # Wave kinematics
    u,du = kinematics2d(ai, fi, ki, epsi, h, t, zWet, x=x, Wheeler=False)
    u = u-u_struct
    p = inline_load(u, du, D, CD, CM, rho) # Morison inline force N/m

    # Return vector of length of input
    z_ref=np.min(zWet) 
    p_hydro = np.zeros(len(z))
    p_hydro[bWet] = p
    dM_hydro = p_hydro * (z-z_ref)

    # Integration of loads
    F_hydro    = trapezoid(p               , zWet) # [N]
    M_sb_hydro = trapezoid(p * (zWet-z_ref), zWet) # Sea bed moment [Nm]

    return u, du, p_hydro, dM_hydro, F_hydro, M_sb_hydro



if __name__ == '__main__':
    pass
