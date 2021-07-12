"""
2D Lamb-Oseen vortex

See:
 [1] Chapter 1, p.58, Branlard - Wind turbine aerodynamics and vorticity based methods, Springer 2017


"""

import numpy as np
    
def lo_omega(X,Y,Gamma=1,t=1,nu=1, polarIn=False): 
    """ 
    Vorticity distribution for 2D Lamb-Oseen vortex
    """
    if polarIn:
        r = X
    else:
        r = np.sqrt(X ** 2 + Y ** 2)
    
    omega_z = Gamma/(4*np.pi*nu*t) * (np.exp(- r**2/(4*nu*t)))
    return omega_z

def lo_u(X,Y,Gamma=1,t=1,nu=1, polarIn=False, polarOut=False): 
    """ 
    Velocity for 2D Lamb-Oseen vortex
    """

    if polarIn:
        r = X
        theta = Y
    else:
        r = np.sqrt(X ** 2 + Y ** 2)
        theta = np.arctan2(Y,X)
    Ut = np.zeros(X.shape)
    # r>0
    bPos = r>1e-12
    r=r[bPos]
    Ut[bPos] = Gamma/(2*np.pi*r) * (1 - np.exp(-r**2/(4*nu*t)))
    if polarOut:
        U = Ut * 0
        V = Ut
    else:
        U = -np.sin(theta)*Ut
        V =  np.cos(theta)*Ut
    return U,V
