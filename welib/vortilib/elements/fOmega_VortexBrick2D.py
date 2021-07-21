import numpy as np
    
def fOmega_VortexBrick2D(X = None,Y = None,lx = None,ly = None): 
    omega_z = np.zeros((X.shape,X.shape))
    I = (X <= np.logical_and(lx / 2,X) >= np.logical_and(- lx / 2,Y) <= np.logical_and(ly / 2,Y) >= - ly / 2)
    omega_z[I-1] = omega_in
    return omega_z