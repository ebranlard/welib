import numpy as np
    
def fMeshProgressiveExtent(xstart = None,xend = None,dx = None,p = None): 
    # ]xstart xend]
    
    dx = np.abs(dx)
    s = sign(xend - xstart)
    dx = dx * s
    p = 1 + p
    x = []
    
    xlast = xstart
    xi = xlast + dx * p
    dx_last = dx * p
    while np.abs(xi) < np.abs(xend):

        x = np.array([x,xi])
        xlast = xi
        xi = xlast + dx_last * p
        dx_last = dx_last * p

    
    x = __builtint__.sorted(x)
    return x