import numpy as np
import pandas as pd

# import welib.hydro.wavekin as wk


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
#     F_total = np.trapz(P, z, axis=0)      # [N] integrated load
#     M_total = np.trapz(dM_ref, z, axis=0) # [Nm]

    return p_drag + p_inertia

if __name__ == '__main__':
    pass
