import numpy as np


def cylindervert_addedmass(R, z1, z2, rho, Ca=1, AxCa=1,
        m_f=0, z_f=0, m_mg=0, z_mg=0):
    """
     Hydrodynamic added mass matrix of a vertical cylinder 

     - R    :  cylinder radius (m)
     - z1   :  vertical position of point 1 (m) positive above water level
     - z2   :  vertical position of point 2 (m) positive above water level
     - rho  :  water density (kg/m^3)
     - Ca   :  added mass coefficient perpendicular to the cylinder
     - AxCa :  added mass coefficient normal to the cylinder (for heave)

     TODO TODO TODO Not finished (empirical equations)
     TODO might need modif if fully submerged
    """
    if z1<z2:
        raise Exception('z1 should be above z2')
    if z1<0:
        # Fully submerged
        ztop = z1
        A0=0
        nAx=2
    else:
        # Partially submerged
        ztop = 0 
        A0 = np.pi*R**2 # undisplaced waterplane area of platform (m^2)
        nAx=1

    h   = ztop-z2      # submerged height
    z_b = (ztop+z2)/2  # coordinates of the center of buoyancy of the undisplaced platform (m)
    V0  = np.pi*R**2*h # undisplaced volume of platform (m^3)

    M=np.zeros((6,6))
    M[0,0] =  Ca * rho*V0
    M[1,1] =  Ca * rho*V0
    M[2,2] =  nAx*AxCa * 2/3*rho*np.pi * R**3  # rho*V0* D/(3*h)
    M[4,0] =  M[0,0]*z_b # TODO empirical
    M[3,1] = -M[0,0]*z_b # TODO empirical
    M[0,4] =  M[0,0]*z_b # TODO empirical
    M[1,3] = -M[0,0]*z_b # TODO empirical

    T1 =Ca*rho*np.pi*R**2 * h**3 /3  # Ca * rho*V0 * h**2/3 # TODO a bit empirical
    M[3,3] = T1
    M[4,4] = T1
    return M
