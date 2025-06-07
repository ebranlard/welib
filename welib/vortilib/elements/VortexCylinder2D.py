""" 
Cyclic or acyclic - 2D vortex cylinder (for 3D see VortexCylinder.py)


Reference: 
    [1] Branlard (2017) Wind turbine aerodynamics and vorticity-based methods

"""
import numpy as np

def vc_coords(R, ne=100):
    theta = np.linspace(-np.pi,np.pi,ne)
    xc    = R * np.cos(theta)
    yc    = R * np.sin(theta)
    return xc, yc

def vc_coords_yx(xc, R, ne=100):
    yc  = R*np.sin(np.arccos(xc/R))
    return yc

def vc_Cptheta(theta, U0=1, R=1, Gamma=0, alpha=0):
    """
        Cp(theta) for a lifting cylinder
    For Gamma=0:
        Cp(theta) 1 - 4 * sin(theta)**2
    """
    Cp = 1- (2*np.sin(theta+alpha) - Gamma/(2*np.pi*R*U0))**2
    return Cp

def vc_Cpx(x, U0=1, R=1, Gamma=0, alpha=0):
    """
        Cp(x) for a lifting cylinder
    For Gamma=0:
        Cp(theta) 1 - 4 * sin(theta)**2
    """
    theta = np.arccos(x/R)
    return vc_Cptheta(theta, U0=U0, R=R, Gamma=Gamma, alpha=alpha)

def vc_Cp(X, Y, P=[0,0], U0=1, R=1, Gamma=0, alpha=0):
    U, V = vc_u(X, Y, P, U0=U0, R=R, Gamma=Gamma, alpha=alpha)
    U2 = U**2 + V**2 
    Cp = 1-U2/U0**2
    return Cp

def vc_stag(U0=1, R=1, Gamma=0, alpha=0):
    """ 
    return location of stagnation points
    """
    ratio = Gamma/(4*np.pi*U0*R)
    if abs(ratio)==1:
        # Only one point at the bottom
        x=np.array([R*np.cos(-np.pi/2 + alpha)])
        y=np.array([R*np.sin(-np.pi/2 + alpha)])
    elif abs(ratio)<1:
        # Two points
        theta1 = np.arcsin(ratio)
        theta2 = np.pi-np.arcsin(ratio)
        theta = np.array([theta1, theta2])
        x= R*np.cos(theta+alpha)
        y= R*np.sin(theta+alpha)
    elif abs(ratio)>1:
        rat = ratio*R
        r1 = rat + np.sqrt(rat**2 - R**2)
        r2 = rat - np.sqrt(rat**2 - R**2)
        r = np.array([r1, r2])
        x= r*np.cos(np.pi/2+alpha)
        y= r*np.sin(np.pi/2+alpha)
    return x, y


def vc_u(X, Y, P=[0,0], U0=1, R=1, Gamma=0, alpha=0):
    """
    2D Cylinder with possible free stream and circulation
    """
    r     = np.sqrt((X-P[0])**2 + (Y-P[1])**2) # Radius at mesh locations
    theta = np.arctan2(Y-P[1], X-P[0])          # 
    U = U0*np.cos(alpha) - (U0*((R/r)**2)*np.cos(2*theta - alpha))  - Gamma*np.sin(theta)/(2*np.pi*r)
    V = U0*np.sin(alpha) - (U0*((R/r)**2)*np.sin(2*theta - alpha))  + Gamma*np.cos(theta)/(2*np.pi*r)
    return U, V
# U(r<rc-10^-14)=0;
# V(r<rc-10^-14)=0;


def vc_psi(x, y, R=1, U0=1):
    Mu = U0*R**2
    psi = U0 * y - Mu * y/(x**2+y**2) - Gamma/(4*np.pi)*np.log(x**2+y**2)
    return psi


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    # --- Compare Rankine nose and cylinder
    U0=1
    R=2
    Gamma_ref =  -4*np.pi*U0*R
    Gamma=Gamma_ref*(1/(2*np.pi))
    alpha=0*np.pi/180
    Cl = -2 *Gamma /(R*U0) # = 8 * np.pi*np.sin(alpha)
    print('Cl', Cl)

    Sigma = 2*U0*R
    x0 = Sigma/(2*np.pi*U0)

    x_cyl  = np.linspace(-R, R, 200)
    Cp_cyl = vc_Cpx(x_cyl, U0 = U0, R = R, Gamma = Gamma, alpha = alpha)

    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(x_cyl, Cp_cyl, label='Cylinder')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.legend()

    plt.show()
