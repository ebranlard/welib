import numpy as np
    
def fOmega_VortexBrick3D(X = None,Y = None,Z = None,lx = None,ly = None,lz = None): 
    omega_z = np.zeros((X.shape,X.shape))
    I = (X <= np.logical_and(lx / 2,X) >= np.logical_and(- lx / 2,Y) <= np.logical_and(ly / 2,Y) >= np.logical_and(- ly / 2,Z) <= np.logical_and(lz / 2,Z) >= - lz / 2)
    omega_z[I-1] = omega_in
    return omega_z