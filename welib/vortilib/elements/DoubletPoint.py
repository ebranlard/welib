""" 
2D and 3D doublet point

Reference: 
    [1] Branlard (2017) Wind turbine aerodynamics and vorticity-based methods
"""
import numpy as np



# --------------------------------------------------------------------------------}
# --- 2D 
# --------------------------------------------------------------------------------{
def dp2d_u(X, Y, Pd, Mu=1, alpha=0, regParam=0, regMethod=None, grad=False): 
    """ 
    Velocity from one double at multiple control points with coordinates X and Y
    """
    DX = X - Pd[0]
    DY = Y - Pd[1]
    theta = np.arctan2(DY, DX)
    r2 = DX ** 2 + DY ** 2
    rX =  DX
    rY =  DY
    tX = -DY
    tY =  DX
    if regMethod is None:
        k = Mu/r2 
        U = k * (-rX*np.cos(alpha-theta) + tX*np.sin(alpha-theta))  # [1] Eq. 32.14
        V = k * (-rY*np.cos(alpha-theta) + tY*np.sin(alpha-theta))
    else:
        k = Mu/r2 * (1 - np.exp(- r2 / regParam ** 2))
        U = k * (-rX*np.cos(alpha-theta) + tX*np.sin(alpha-theta))
        V = k * (-rY*np.cos(alpha-theta) + tY*np.sin(alpha-theta))
    return U, V



def dp2d_psi(x, y, Pd, Mu=1, alpha=0 ):
    """ Streamfunction from one 2D doublet point """
    z = x-Pd[0] + (y-Pd[1]) * 1j
    fac = Mu*np.exp(1j*alpha)
    psi = np.imag(fac/z) # [1] Eq 32.15
    #psi = - Mu * y/(x**2+y**2)
    return psi

def dp2d_phi(x, y, Pd, Mu=1, alpha=0 ):
    """ Potential from one 2D doublet point """
    z = x-Pd[0] + (y-Pd[1]) * 1j
    fac = Mu*np.exp(1j*alpha)
    phi = np.real(fac/z) # [1] Eq 32.15
    #phi = Mu * x/(x**2+y**2)
    return phi

