""" 
2D inviscid vortex patch (ivp)

omega = (1-r**2)**k

See:
 [1] Chapter 33, p.402, Branlard - Wind turbine aerodynamics and vorticity based methods, Springer 2017


"""

import numpy as np
    
def ivp_omega(X,Y,k=2,polarIn=False): 
    """ 
    Vorticity distribution for 2D inviscid vortex patch (ivp)
    """
    X=np.asarray(X)
    Y=np.asarray(Y)
    if polarIn:
        r = X
        theta = Y
    else:
        r = np.sqrt(X ** 2 + Y ** 2)
        theta = np.arctan2(Y,X)
    omega_z=np.zeros(r.shape)
    # r<=1
    bLow=r<=1
    rl = r[bLow]
    omega_z[bLow] = (1-rl**2)**k
    # r>1
    omega_z[np.abs(r)>1] = 0
    # r=0
    omega_z[np.abs(r)<1e-14] = 1
    return omega_z

def ivp_Gamma(r,k=2): 
    """ 
    Circulation distribution for 2D inviscid vortex patch (ivp)
    """
    r=np.asarray(r)
    Gamma=np.zeros(r.shape)
    # r<=1
    bLow=r<=1
    rl = r[bLow]
    Gamma[bLow] = np.pi/(k+1)* (1-(1-rl**2)**(1+k))
    # r>1
    Gamma[np.abs(r)>1] = np.pi/(k+1)
    return Gamma


def ivp_u(X, Y, k=2, polarIn=False, polarOut=False): 
    """
    Velocity distribution for 2D inviscid vortex patch (ivp)

    if polarIn is True:
       X=r
       Y=theta
    """
    if polarIn:
        r = X
        theta = Y
    else:
        r = np.sqrt(X ** 2 + Y ** 2)
        theta = np.arctan2(Y,X)
    Ut      = np.zeros(X.shape)
    # r<=1
    bLow= np.logical_and(r<=1, r>1e-14)
    rl = r[bLow]
    Ut[bLow]       = (1 - (1-rl**2)**k +  (rl**2 *(1-rl**2)**k))/(2*(1+k)*rl)
    # r>1
    bHigh=r>1
    Ut[bHigh] = 1.0/(2*(1+k)*r[bHigh])
    # r=0
    bZero=r<1e-14
    Ut[bZero]      = 0
    if polarOut:
        U = Ut * 0
        V = Ut
    else:
        U = -np.sin(theta)*Ut
        V =  np.cos(theta)*Ut
    return U,V
