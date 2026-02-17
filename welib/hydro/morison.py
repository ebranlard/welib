import numpy as np
import pandas as pd
try:
    from numpy import trapezoid
except:
    from numpy import trapz as trapezoid

from welib.hydro.wavekin import kinematics2d


def inline_load(u_rel, a_wav, a_rel, D, rho, Cd, Cp=None, Ca=None, CM=None):
    """
    Morison force per length on a structure defined along the z axis
    INPUTS:
     - u_rel: relative velocity u_wav-u_str, array like of shape (nz x ...)
     - a_wav: acceleration, array like of shape (nz x ...)
     - a_rel: relative acceleration, a_wav-a_st array like of shape (nz x ...)
     - z: vertical distance along the structure (typically negative and 0 at sea bed) 
     - D: structure diameter, scalar or array of shape (nz x ...)
     - Cd: drag coefficient, scalar or array of shape (nz x ...)
     - Cp: Froude-Krylov, pressure coefficient, scalar or array of shape (nz x ...)
     - Ca: Added mass coefficient, scalar or array of shape (nz x ...)
     - CM: inertia coefficient (CM=Cp+Ca), scalar or array of shape (nz x ...)

    NOTES:
     - For linear waves, Cp 
     - Ca is added mass only (not including FK).
     - p = 0.5 ? Cd D u_rel|u_rel| + ? (1+Ca) A a_wave - ? Ca A a_struct
    """

    # Drag (relative velocity)
    p_drag = 0.5 * rho * Cd * D * u_rel * np.abs(u_rel)

    a_str = a_wav - a_rel

    A = np.pi * D**2 / 4
    if Ca is None and CM is not None:
        Cp = CM-1
        Ca = 1
    p_FK      = rho * Cp * A * a_wav # Froude-Krylov (pressure)
    p_AM      = rho * Ca * A * a_rel  # Added-mass (relative acceleration)
    p_AM_wav  = rho * Ca * A * a_wav  # Added-mass (relative acceleration)
    p_AM_str  =-rho * Ca * A * a_str  # Added-mass (relative acceleration)
    p_inertia =      p_FK + p_AM
    #p_inertia =  rho * CM * A * a_wav # NOTE: Not including Ca contrib here

    p_tot = p_drag + p_inertia
    return p_tot, p_drag, p_AM, p_FK, p_AM_wav, p_AM_str




def monopileHydroLoads1D(t, ai, fi, ki, epsi, h, z, x, D, rho, Cd, Cp, Ca, u_struct, a_struct=None, CM=None):  
    """ 
    Compute hydrodynamic loads (Morison equation) on monopile
    Wave kinematics are computed on the fly at undeflected monopile position

    Acceleration is optional as it may be accounted for elsewhere

    INPUTS:
     - t: time, scalar [s]
     - fi  : (nf-array) frequencies, for wave kinematics
     - ai  : (nf-array) amplitudes
     - ki  : (nf-array) wave numbers
     - epsi: (nf-array) phases
     - h: water depth >0, the sea bed is at z=-h [m]
     - z: (n-array) vertical positions defining monopile (and possibly tower) sections
    """
    if a_struct is None:
        a_struct = np.zeros_like(u_struct)
        
    # Limit calculations to below water and above sea bed
    bWet = np.logical_and(z<=0, z>=-h)  # TODO eta
    D    = D[bWet].flatten()
    Cd   = Cd[bWet].flatten()
    Cp   = Cp[bWet].flatten()
    Ca   = Ca[bWet].flatten()
    zWet = z[bWet].flatten()
    if CM is not None:
        CM   = CM[bWet].flatten()
    u_struct = u_struct[bWet].flatten()
    a_struct = a_struct[bWet].flatten()
    
    # Wave kinematics
    u_wave, a_wave = kinematics2d(ai, fi, ki, epsi, h, t, zWet, x=x, Wheeler=False)
    # Relative motion
    u_rel = u_wave - u_struct
    a_rel = a_wave - a_struct

    p, p_drag, p_AM, p_FK, p_AM_wav, p_AM_str = inline_load(u_rel, a_wave, a_rel, D, rho, Cd=Cd, Cp=Cp, Ca=Ca, CM=CM) # Morison inline force N/m

    # Return vector of length of input
    z_ref=np.min(zWet) 
    p_hydro = np.zeros(len(z)); p_hydro[bWet] = p

    P_drag   = np.zeros(len(z)); P_drag[bWet]   = p_drag
    P_AM     = np.zeros(len(z)); P_AM[bWet]     = p_AM
    P_FK     = np.zeros(len(z)); P_FK[bWet]     = p_FK
    P_AM_wav = np.zeros(len(z)); P_AM_wav[bWet] = p_AM_wav
    P_AM_str = np.zeros(len(z)); P_AM_str[bWet] = p_AM_str

    dM_hydro = p_hydro * (z-z_ref)

    # Integration of loads
    F_hydro    = trapezoid(p               , zWet) # [N]
    M_sb_hydro = trapezoid(p * (zWet-z_ref), zWet) # Sea bed moment [Nm]

    outH=dict()
    outH['u_wav'] = u_wave
    outH['a_wav'] = a_wave
    outH['bWet'] = bWet
    outH['u_rel'] = u_rel
    outH['a_rel'] = a_rel
    outH['p_hydro']  = p_hydro
    outH['p_drag']   = P_drag
    outH['p_AM']     = P_AM
    outH['p_FK']     = P_FK
    outH['p_AM_str'] = P_AM_str
    outH['p_AM_wav'] = P_AM_wav
    outH['dM_hydro'] = dM_hydro
    outH['F_hydro'] = F_hydro
    outH['M_sb_hydro'] = M_sb_hydro
    return outH




if __name__ == '__main__':
    pass
