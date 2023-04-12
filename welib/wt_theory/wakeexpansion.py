""" 
Wake expansion models

Main variables used:
- xb: dimensionless downstream distance, xb=x/R
- CT: Thrust coefficient
- a : axial induction factor


References:
 [1] Branlard (2017) Wind turbine aerodynamics and vorticity based methods
 [2] Rathman et al. (2006) Turbine wake model for wind ressource software
 [3] Frandsen et al. (2006) Analytical modelling of wind speed deficit in large offshore wind farms
"""
import numpy as np

def wake_expansion(xb, CT=None, a=None, model='cylinder', **kwargs):
    """ 
    Switch between different wake expansion models
    """
    xb = np.asarray(xb)
    if model=='momentum':
        return wake_expansion_momentum(xb, CT=CT, a=a)
    elif model=='cylinder':
        return wake_expansion_cylinder(xb, CT=CT, a=a)
    elif model=='Rathmann':
        return wake_expansion_Rathmann(xb, CT=CT, a=a, **kwargs)
    elif model=='Frandsen':
        return wake_expansion_Frandsen(xb, CT=CT, a=a, **kwargs)
    else:
        raise NotImplementedError()

def wake_expansion_momentum(xb=None, CT=None, a=None):
    """ 
    Wake expansion based on momentum theory
        (Rw/R)^2 = (1-a)/(1-2a)     [1] Eq. 15.1 
    """
    if CT is not None:
        a    = 1/2*(1-np.sqrt(1-CT))
    r = np.sqrt((1-a)/(1-2*a))
    if xb is not None:
        xb = np.asarray(xb)
        r = xb*0 + r
    return r

def wake_expansion_cylinder(xb, CT=None, a=None):
    """ 
    Wake expansion based on vortex cylinder axis velocity

    See [1] Eq. 15.1 

    INPUTS: 
    """
    xb = np.asarray(xb)
    if CT is not None:
        a    = 1/2*(1-np.sqrt(1-CT))
    r2 = (1-a) /  (1-a * (1+ xb/np.sqrt(1+xb**2) ) )
    return np.sqrt(r2)

def wake_expansion_Rathmann(xb, CT=None, a=None, alpha=0.7):
    """ 
    [2]
    """
    xb = np.asarray(xb)
    assert(np.all(xb>=0))
    if CT is not None:
        a    = 1/2*(1-np.sqrt(1-CT))
    beta  = (1 - a / 2) / (1 - a)
    r =  (1/2*alpha * xb) ** (1/2)
    r = np.clip(r, beta**(1/2), np.inf )
    return r

def wake_expansion_Frandsen(xb, CT=None, a=None, alpha=0.7,k=2):
    """ 
    [3]
    """
    xb = np.asarray(xb)
    if CT is not None:
        a = (1 - np.sqrt(1 - CT))
    
    beta = (1 - a / 2) / (1 - a)
    r = (beta**(k/2) + 1 / 2 * alpha * xb) ** (1 / k)
    r[xb==0] = 1
    return r

def wake_expansion_kxm(xb, k=0.05, m=1):
    """ """
    xb = np.asarray(xb)
    return 1 + k*(xb/2)**m


def downstreamDistanceForGivenExpansion(CT, expansion, model='cylinder', method='interp',**kwargs):
    """ 
    Return downstream distance for a given expansion
    """
    if method=='interp':
        xb        = np.linspace(0,20,20000)
        expansions = wake_expansion(xb, CT  = CT, model = model, **kwargs)
        xReached = np.interp(expansion, expansions, xb)

    elif method=='analytical':
        if model=='cylinder':
            a = 1/2*(1-np.sqrt(1-CT))
            k = (expansion**2*(1-a) + a - 1 ) / (expansion**2 * a )
            xReached = k/np.sqrt(1-k**2)
            return xReached
        else:
            raise NotImplementedError
    return xReached




if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import welib.tools
    
    xb = np.linspace(0,20, 100) 
    rwC = wake_expansion_cylinder(xb=xb, CT=0.8)
    rwR = wake_expansion_Rathmann(xb=xb, CT=0.8)
    rwF = wake_expansion_Frandsen(xb=xb, CT=0.8)
    rwK = wake_expansion_kxm     (xb=xb, k=0.05, m=1)
    rw0 = wake_expansion_momentum(CT=0.8)

    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(xb, xb*0+rw0, 'k-', label='NW - Momentum theory')
    ax.plot(xb, rwC, 'k--'    , label='NW - Vortex Cylinder')
    ax.plot(xb, rwF, ':'      , label='FW - Frandsen')
    ax.plot(xb, rwR, ':'      , label='FW - Rathmann')
    ax.plot(xb, rwK, ':'      , label=r'FW $k x^m$')
    ax.set_xlabel('z/R [-]')
    ax.set_ylabel('r/R [-]')
    ax.legend()
    plt.show()



   
