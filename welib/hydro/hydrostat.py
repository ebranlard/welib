import numpy as np


def cylindervert_hydrostat(R, z1, z2, rho, g,
        m_f=0, z_f=0, m_mg=0, z_mg=0):
    """
     Hydrostatic restoring matrix of a vertical cylinder 

     - R    :  cylinder radius (m)
     - z1   :  vertical position of point 1 (m) positive above water level
     - z2   :  vertical position of point 2 (m) positive above water level
     - rho  :  water density (kg/m^3)
     - g    : gravity (m/s^2)

     OPTIONAL INPUTS:
     - m_mg:  total mass of marine growth (kg)
     - z_mg:  coordinates of the center of mass of the undisplaced marine growth mass (m)
     - m_f :  total mass of ballasting/flooding (kg)
     - z_f :  coordinates of the center of mass of the undisplaced filled fluid (flooding or ballasting) mass (m)

     - A_{0}: undisplaced waterplane area of platform (m^2)
     - V_{0}: undisplaced volume of platform (m^3)
    """
    if z1<z2:
        raise Exception('z1 should be above z2')
    if z1<0:
        # Fully submerged
        A0=0
        h=z1-z2 # submerged height
    else:
        h  = -z2        # submerged height
        A0 = np.pi*R**2 # undisplaced waterplane area of platform (m^2)

    V0 = np.pi*R**2*h         # undisplaced volume of platform (m^3)
    z_b = h/2         # coordinates of the center of buoyancy of the undisplaced platform (m)

    K=np.zeros((6,6))
    K[2,2] = rho * g * A0
    K[3,3] = rho * g * A0 * R**2/4 + rho*g*V0*z_b- m_mg * g * z_mg - m_f *g *z_f
    K[4,4] = K[3,3]
    return K
