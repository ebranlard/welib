"""
2D Gaussian vorticity patch (vpg), similar to LambOseen

Also Asymetric gaussain vortex patch (vpga) with sigmaX and sigmaY

See:
 [1] Chapter 1, p.58, Branlard - Wind turbine aerodynamics and vorticity based methods, Springer 2017


"""

import numpy as np
    
def vpg_omega(X,Y,Gamma=1, sigma=1, polarIn=False): 
    """ 
    Vorticity distribution for 2D Gaussian vortex patch
    """
    if polarIn:
        r = X
    else:
        r = np.sqrt(X ** 2 + Y ** 2)
    
    omega_z = Gamma/(np.pi*sigma) * (np.exp(- r**2/sigma**2))
    return omega_z


def vpg_u(X,Y,Gamma=1, sigma=1, polarIn=False, polarOut=False): 
    """ 
    Velocity for 2D Gaussian vortex patch
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
    Ut[bPos] = Gamma/(2*np.pi*r) * (1 - np.exp(-r**2/(sigma**2)))
    if polarOut:
        U = Ut * 0
        V = Ut
    else:
        U = -np.sin(theta)*Ut
        V =  np.cos(theta)*Ut
    return U,V

# --------------------------------------------------------------------------------}
# --- Asymmetric Gaussian 
# --------------------------------------------------------------------------------{
def vpga_omega(X,Y,Gamma=1, sigmaX=3, sigmaY=1, polarIn=False): 
    """ 
    Vorticity distribution for 2D Gaussian vortex patch asymmetric  
    """
    if polarIn:
        r     = X
        theta = Y
        X     = r*np*cos(theta)
        Y     = r*np*sin(theta)
    omega_z = Gamma/(np.pi*sigmaX*sigmaY) * (np.exp(- (X**2/sigmaX**2 + Y**2/sigmaY**2)))
    return omega_z


def vpga_u(X,Y,Gamma=1, sigmaX=3, sigmaY=1, polarIn=False, polarOut=False): 
    """ 
    Velocity for 2D Gaussian vortex patch asymmetric
    """
    print('>>>> TODO vpga is wrong for now')

    if polarIn:
        r     = X
        theta = Y
        X     = r*np*cos(theta)
        Y     = r*np*sin(theta)
    else:
        r = np.sqrt(X ** 2 + Y ** 2)
        theta = np.arctan2(Y,X)
    # TODO, conformal map?
    Ut = np.zeros(X.shape)

    # Try1
    Y2= Y*sigmaX/sigmaY
    X2= X*sigmaY/sigmaX
    U1,V1 = vpg_u(X,Y2,Gamma=Gamma, sigma=sigmaX, polarIn=False, polarOut=False)
    U2,V2 = vpg_u(X2,Y,Gamma=Gamma, sigma=sigmaY, polarIn=False, polarOut=False)
    U=U1*np.sqrt(sigmaX)
    V=V1*np.sqrt(sigmaX)
    #U=U1*(sigmaY/sigmaX)**2
    #V=V1*sigmaY/sigmaX

    # Try2
    #U1,V1 = vpg_u(X,Y,Gamma=Gamma, sigma=sigmaX, polarIn=False, polarOut=False)
    #k=1/sigmaY**2-1/sigmaX**2
    #U1*=(1-24*k)
    #V1*=(1-24*k)
    #U=U1
    #V=V1

    # r>0
    bPos = r>1e-12
    r=r[bPos]

    #Ut[bPos] = Gamma/(2*np.pi*r) * (1 - (np.exp(- (X**2/sigmaX**2 + Y**2/sigmaY**2))))
    #V=V*sigmaY/sigmaX*3 # <<<< TODO
    #U*=1.5
    ## r>0
    #bPos = r>1e-12
    #r=r[bPos]
    #Ut[bPos] = Gamma/(2*np.pi*r) * (1 - np.exp(-r**2/(sigma**2)))
    #if polarOut:
    #    U = Ut * 0
    #    V = Ut
    #else:
    #    #U = -sigmaX*np.sin(theta)*Ut
    #    #V =  sigmaY*np.cos(theta)*Ut
    #    U = -np.sin(theta)*Ut
    #    V =  np.cos(theta)*Ut
    return U,V


